#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"
#include <stdatomic.h>
#include <stdlib.h>

#define SIZE 30 * 1024            // 30720 double，30 KB.
__gsm__ char gx[6 * 1024 * 1024]; // 保证 16 字节对齐(双精度：786,432)

#ifndef BLOCK_SIZEN
#define BLOCK_SIZEN 780000
#endif

// __gsm__ double gsm_x[512 * 128 * 6]; // 3MB GSM for x / 393216
// __gsm__ double gsm_y[512 * 128 * 6]; // 3MB GSM for y / 393216

__global__ void kernel_null()
{
}

// =================================================================
// 双精度：矩阵向量乘（SPU 版本）
// =================================================================
__global__ void kernel_SpMV_SPU(int num_cols,
                                int *A_i,
                                int *A_j,
                                double *A_data,
                                double *x_data,
                                double *y_data,
                                int *row_bounds)
{
    int tid = get_thread_id();

    int startrow = row_bounds[tid];
    int stoprow = row_bounds[tid + 1];

    for (int i = startrow; i < stoprow; i++)
    {
        for (int jj = A_i[i]; jj < A_i[i + 1]; jj++)
        {
            int j = A_j[jj];
            if (j >= 0 && j < num_cols)
            {
                y_data[i] += A_data[jj] * x_data[j];
            }
        }
    }
}

// =============================================================
// 内积版本三个内核函数，x 主机端预处理，vector_load 加载
// 版本一(A)：VPU 端按行计算，寄存器规约
// 版本二(B)：CPU 端规约，VPU->CPU 顺序执行（OMP 并行情形下效果最优）
// 版本三(C)：CPU 端规约，VPU/CPU 异步并发
// =============================================================
__global__ void kernel_v2(int coreNums,
                          int num_tasks,
                          int *A_i,
                          double *A_data,
                          double *x_data,
                          double *y_data,
                          int *row_bounds)
{
    int tid = get_thread_id();

    // if (tid == 0)
    // {
    //     hthread_printf("kernel_v2 \n");
    // }

    int AM_size = 96000 / 3;
    int step = AM_size;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result_vector;

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {

        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];
        if (task_end_row <= task_start_row)
            continue;

        int start_nnz = A_i[task_start_row];
        int stop_nnz = A_i[task_end_row];
        int length = stop_nnz - start_nnz;
        if (length <= 0)
            continue;

        // 加载当前分块到 AM
        vector_load(&A_data[start_nnz], A_cache, length * sizeof(double));
        vector_load(&x_data[start_nnz], x_cache, length * sizeof(double));

        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);

            // 向量乘法计算
            result_vector = vec_svbcast(0.0);
            result_vector = vec_muli(A_value, x_value);

            // 将结果存回 AM 缓存
            vec_st(result_vector, offset, x_cache);
        }

        if (vec_iters > 0)
        {
            vector_store(x_cache, &y_data[start_nnz], vec_iters * 16 * sizeof(double));
        }

        // 尾部（长度不足16个元素）标量处理
        int tail_len = length % 16;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 16; // 尾部在当前块内的偏移
            for (int k = 0; k < tail_len; k++)
            {
                int global_idx = start_nnz + remain_start + k;
                double A_scalar = A_data[global_idx];
                double x_scalar = x_data[global_idx];
                y_data[global_idx] = A_scalar * x_scalar;
            }
        }
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// =====================================================
// 将向量 x 延长, SPMV 转为 DOT 计算
// 任务划分：    先 DDR 后 AM 的二级划分
// x 加载方式：  vector_load
// 规约方式：    CPU端规约
// 任务执行方式: 顺序执行
// =====================================================
__global__ void kernel_v3(int coreNums,
                          int num_tasks,     // 任务划分总总段数
                          int nbuf_per_core, // 单个队列缓冲区个数 N
                          int *A_i,
                          double *A_data,
                          double *x_data,
                          double *result_buffers,
                          int *row_bounds,
                          volatile long long *d_in,
                          volatile long long *d_out)
{
    int tid = get_thread_id();

    // 定义常量，异步并发版本需与主机端保持一致
    int AM_size = 96000 / 3;
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {
        // 1. 检查自己的队列是否已满
        while (d_in[tid] - d_out[tid] >= nbuf_per_core)
        {
            // 队列已满，原地等待消费者处理
        }

        // 2. 计算当前任务在 result_buffers 中的写入位置
        int local_idx = d_in[tid] % nbuf_per_core;
        int global_idx = tid * nbuf_per_core + local_idx;
        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;

        // 3. 获取当前任务的计算范围
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];
        int start_nnz = A_i[start_row];
        int stop_nnz = A_i[end_row];
        int length = stop_nnz - start_nnz;

        vector_load(&A_data[start_nnz], A_cache, length * sizeof(double));
        vector_load(&x_data[start_nnz], x_cache, length * sizeof(double));

        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, offset, x_cache);
        }

        if (vec_iters > 0)
        {
            vector_store(x_cache, current_buffer_ptr, vec_iters * 16 * sizeof(double));
        }

        int tail_len = length % 16;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 16;
            for (int k = 0; k < tail_len; k++)
            {
                int data_idx = start_nnz + remain_start + k;
                int buffer_idx = remain_start + k;
                current_buffer_ptr[buffer_idx] = A_data[data_idx] * x_data[data_idx];
            }
        }

        current_buffer_ptr[AM_size] = (double)start_row;
        current_buffer_ptr[AM_size + 1] = (double)end_row;

        // __sync_synchronize();
        d_in[tid]++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
}

// =====================================================
// 将向量 x 延长, SPMV 转为 DOT 计算
// 任务划分：    先 DDR 后 AM 的二级划分
// x 加载方式：  vector_load
// 规约方式：    CPU端规约
// 任务执行方式: 顺序执行，与 v3 相比，引入了四个单独的标记
// =====================================================
__global__ void kernel_v3B(int coreNums,
                           int num_tasks,     // 任务划分总总段数
                           int nbuf_per_core, // 单个队列缓冲区个数 N
                           int *A_i,
                           double *A_data,
                           double *x_data,
                           double *result_buffers,
                           int *row_bounds,
                           volatile long *d0_in,
                           volatile long *d0_out,
                           volatile long *d1_in,
                           volatile long *d1_out,
                           volatile long *d2_in,
                           volatile long *d2_out,
                           volatile long *d3_in,
                           volatile long *d3_out)
{
    int tid = get_thread_id();

    volatile long *d_in_p[4] = {d0_in, d1_in, d2_in, d3_in};
    volatile long *d_out_p[4] = {d0_out, d1_out, d2_out, d3_out};

    // 获取当前线程应该操作的 in/out 指针
    volatile long *my_d_in = d_in_p[tid];
    volatile long *my_d_out = d_out_p[tid];

    // 定义常量，异步并发版本需与主机端保持一致
    int AM_size = 96000 / 3;
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {
        // 1. 检查自己的队列是否已满
        while ((*my_d_in) - (*my_d_out) >= nbuf_per_core)
        {
            // 队列已满，原地等待消费者处理
        }

        // 2. 计算当前任务在 result_buffers 中的写入位置
        int local_idx = (*my_d_in) % nbuf_per_core;
        int global_idx = tid * nbuf_per_core + local_idx;
        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;

        // 3. 获取当前任务的计算范围
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];
        int start_nnz = A_i[start_row];
        int stop_nnz = A_i[end_row];
        int length = stop_nnz - start_nnz;
        vector_load(&A_data[start_nnz], A_cache, length * sizeof(double));
        vector_load(&x_data[start_nnz], x_cache, length * sizeof(double));

        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, offset, x_cache);
        }

        if (vec_iters > 0)
        {
            vector_store(x_cache, current_buffer_ptr, vec_iters * 16 * sizeof(double));
        }
        int tail_len = length % 16;

        if (tail_len > 0)
        {
            int remain_start = vec_iters * 16;
            for (int k = 0; k < tail_len; k++)
            {
                int data_idx = start_nnz + remain_start + k;
                int buffer_idx = remain_start + k;
                current_buffer_ptr[buffer_idx] = A_data[data_idx] * x_data[data_idx];
            }
        }

        current_buffer_ptr[AM_size] = (double)start_row;
        current_buffer_ptr[AM_size + 1] = (double)end_row;
        (*my_d_in)++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
}

// 增加对零行的判断 =======================
__global__ void spmv_baseline(int *A_i,
                              int *A_j,
                              double *A_data,
                              double *x_data,
                              double *y_data,
                              int *row_bounds)
{
    int tid = get_thread_id();

    // 初始化本地缓存和向量寄存器
    int AM_size = 96000 / 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result;

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    int warned_long_row = 0;

    // 按“行组”(row group) 为单位，处理分配给本线程的所有行
    while (current_row < tid_end)
    {
        int current_row_nnz = A_i[current_row + 1] - A_i[current_row];
        if (current_row_nnz > AM_size && !warned_long_row)
        {
            hthread_printf("[spmv_baseline warning] tid = %d, row = %d has %d nonzeros, "
                           "which exceeds AM_size = %d. This row may overflow A_cache/x_cache.\n",
                           tid, current_row, current_row_nnz, AM_size);
            warned_long_row = 1;
        }

        // ====================================================================
        // 步骤 1: 确定本次要处理的行组的边界
        // ====================================================================
        int group_end_row = current_row;
        int nnz_in_group = 0;
        while (group_end_row < tid_end)
        {
            int next_row_nnz = A_i[group_end_row + 1] - A_i[group_end_row];
            if (nnz_in_group + next_row_nnz > AM_size)
            {
                break;
            }
            nnz_in_group += next_row_nnz;
            group_end_row++;
        }
        if (group_end_row == current_row)
            group_end_row++; // 确保至少处理一行

        // ====================================================================
        // 【新增步骤】: 扫描非空行边界，并预处理空行
        // ====================================================================
        int group_first_nnz_row = -1;
        int group_last_nnz_row = -1;

        // 从前向后找到第一个非空行
        for (int i = current_row; i < group_end_row; i++)
        {
            if (A_i[i + 1] > A_i[i])
            {
                group_first_nnz_row = i;
                break;
            }
        }

        // 如果 group_first_nnz_row 未变，说明整个行组都是空的
        if (group_first_nnz_row == -1)
        {
            for (int i = current_row; i < group_end_row; i++)
            {
                y_data[i] = 0.0;
            }
        }
        else
        {
            // 如果存在非空行，则从后向前找到最后一个非空行
            for (int i = group_end_row - 1; i >= current_row; i--)
            {
                if (A_i[i + 1] > A_i[i])
                {
                    group_last_nnz_row = i;
                    break;
                }
            }

            // 预处理前导空行
            for (int i = current_row; i < group_first_nnz_row; i++)
            {
                y_data[i] = 0.0;
            }

            // 预处理尾随空行
            for (int i = group_last_nnz_row + 1; i < group_end_row; i++)
            {
                y_data[i] = 0.0;
            }

            // 定义收缩后的、纯净的计算区间
            int proc_start_row = group_first_nnz_row;
            int proc_end_row = group_last_nnz_row + 1; // 区间为 [start, end)

            // ====================================================================
            // 步骤 2: 在【收缩后的边界】上处理不对齐元素
            // ====================================================================
            int group_start_nnz_idx = A_i[proc_start_row];
            int group_end_nnz_idx = A_i[proc_end_row];

            double head_sum = 0.0; // 用于存储块头部不对齐元素的计算结果
            double tail_sum = 0.0; // 用于存储块尾部不对齐元素的计算结果

            int dma_start_idx = group_start_nnz_idx;
            int dma_end_idx = group_end_nnz_idx;

            if (dma_start_idx % 2 != 0)
            {
                double A_val = A_data[dma_start_idx];
                int column_index = A_j[dma_start_idx] / 8;
                double x_val = x_data[column_index];
                head_sum = A_val * x_val;
                dma_start_idx++;
            }

            if (dma_end_idx % 2 != 0)
            {
                dma_end_idx--;
                double A_val = A_data[dma_end_idx];
                int column_index = A_j[dma_end_idx] / 8;
                double x_val = x_data[column_index];
                tail_sum = A_val * x_val;
            }

            // ====================================================================
            // 步骤 3: 在【收缩后的边界】上执行DMA加载
            // ====================================================================
            int dma_load_count = dma_end_idx - dma_start_idx;
            if (dma_load_count > 0)
            {
                vector_load(&A_data[dma_start_idx], A_cache, dma_load_count * sizeof(double));
                unsigned int ch = dma_sg(x_data,
                                         &A_j[dma_start_idx],
                                         1, dma_load_count * sizeof(double), 0,
                                         x_cache,
                                         1, dma_load_count * sizeof(double), 0);
                dma_wait(ch);
            }

            // ====================================================================
            // 步骤 4: 【行级别】在行组内逐行进行计算
            // ====================================================================
            for (int i = proc_start_row; i < proc_end_row; i++)
            {
                result = vec_svbcast(0.0);
                double sum = 0.0;

                if (i == proc_start_row)
                {
                    sum += head_sum;
                }
                if (i == proc_end_row - 1)
                {
                    sum += tail_sum;
                }

                int row_nnz_start = A_i[i];
                int row_nnz_end = A_i[i + 1];

                if (row_nnz_end - row_nnz_start == 0)
                {
                    y_data[i] = sum;
                    continue;
                }

                if (row_nnz_start < dma_start_idx)
                    row_nnz_start = dma_start_idx;
                if (row_nnz_end > dma_end_idx)
                    row_nnz_end = dma_end_idx;

                int vec_length = row_nnz_end - row_nnz_start;
                if (vec_length > 0)
                {
                    int cache_offset = row_nnz_start - dma_start_idx;
                    int vec_iter = vec_length / 16;
                    int tail_len = vec_length % 16;

                    // --- 向量化主循环 ---
                    for (int j = 0; j < vec_iter; j++)
                    {
                        int element_offset = cache_offset + j * 16;
                        lvector double A_value = vec_ld(element_offset, A_cache);
                        lvector double x_value = vec_ld(element_offset, x_cache);
                        result = vec_mula(A_value, x_value, result);
                    }

                    if (tail_len > 0)
                    {
                        lvector double temp = vec_svbcast(0.0);
                        lvector double tail_result = vec_svbcast(0.0);
                        int tail_offset = cache_offset + vec_iter * 16;
                        lvector double A_value = vec_ld(tail_offset, A_cache);
                        lvector double x_value = vec_ld(tail_offset, x_cache);

                        tail_result = vec_muli(A_value, x_value);
                        vec_st(temp, tail_len, &tail_result);
                        result += tail_result;
                    }
                }

                y_data[i] = reduce_sum_fp64(result) + sum;
            }
        }

        // 更新下一组的起始行
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// CPU/DSP 并发执行版本 ========================================
__global__ void kernel_dotSpMV_FP64C(int num_tasks, // 任务划分总总段数
                                     int num_cores, // 调用 DSP 数量
                                     int N,         // 单个队列缓冲区个数 N
                                     int *A_i,
                                     double *A_data,
                                     double *x_data,
                                     int *row_bounds,
                                     int *per_core_tasks,
                                     double *result_buffers,
                                     volatile long long *d_in,
                                     volatile long long *d_out)
{
    int tid = get_thread_id();

    // 定义常量，异步并发版本需与主机端保持一致
    int AM_size = 96000 / 3;
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    int my_num_tasks = per_core_tasks[tid]; // 从主机端获取每个核心的任务数
    int task_count = 0;

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        while (1)
        {
            long long in = d_in[tid];
            long long out = d_out[tid];
            if (in - out < N)
                break;

            // int i = 0;
            // while (i < 24)
            // {
            //     i++; // 只自增 i，什么都不干
            // }

            volatile int dummy = 0;
            for (int i = 0; i < 36; i++)
            {
                dummy = dummy; // 强制读写，防优化
            }
        }

        // 2. 计算将要使用的缓冲区全局索引 (使用参数 N)
        int local_idx = d_in[tid] % N;
        int global_idx = tid * N + local_idx;

        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];

        if (start_row >= end_row)
        {
            __sync_synchronize();
            d_in[tid]++;
            continue;
        }

        int start_nnz = A_i[start_row];
        int stop_nnz = A_i[end_row];
        int length = stop_nnz - start_nnz;

        vector_load(&A_data[start_nnz], A_cache, length * sizeof(double));
        vector_load(&x_data[start_nnz], x_cache, length * sizeof(double));

        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, offset, x_cache);
        }

        if (vec_iters > 0)
        {
            vector_store(x_cache, current_buffer_ptr, vec_iters * 16 * sizeof(double));
        }

        int tail_len = length % 16;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 16;
            for (int k = 0; k < tail_len; k++)
            {
                int data_idx = start_nnz + remain_start + k;
                int buffer_idx = remain_start + k;
                current_buffer_ptr[buffer_idx] = A_data[data_idx] * x_data[data_idx];
            }
        }

        current_buffer_ptr[AM_size] = (double)start_row;
        current_buffer_ptr[AM_size + 1] = (double)end_row;

        __sync_synchronize();
        d_in[tid]++;
        task_count++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
}

// 2025.12.18 ============================
// GSM 版本 baseline  ====================
__global__ void spmv_baseline_C(int num_rows,
                                int n_cores,
                                int barrier_id,
                                int *A_i,
                                int *A_j,
                                double *A_data,
                                double *x_data,
                                double *y_data,
                                int *row_bounds)
{
    int tid = get_thread_id();

    // 初始化本地缓存和向量寄存器
    int AM_size = 96000 / 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result;

    // 1. 获取 GSM 上的 x 向量起始地址
    double *gsm_x = (double *)gx;

    // 计算每个核心负责搬运的搬运块大小（x 向量的均匀切分）
    int x_chunk_size = ((num_rows + n_cores - 1) / n_cores + 1) & ~1; // 强制偶数
    int x_start_idx = tid * x_chunk_size;
    int x_width = (x_start_idx + x_chunk_size <= num_rows) ? x_chunk_size : (num_rows - x_start_idx);

    if (num_rows <= 786432)
    {
        // 只有当起始索引未越界时才执行搬运
        if (x_start_idx < num_rows && x_width > 0)
        {
            // 使用 dma_p2p 从主存 x_data 搬运到共享内存 gsm_x
            unsigned int c = dma_p2p(&x_data[x_start_idx], 1, x_width * sizeof(double), 0,
                                     &gsm_x[x_start_idx], 1, x_width * sizeof(double), 0,
                                     false, 0);
            dma_wait(c);
        }

        // 所有核心在此同步，确保全量 x 已经在 GSM 准备就绪
        core_barrier(barrier_id, n_cores);
    }
    else
    {
        if (tid == 0)
        {
            hthread_printf("Warning: num_rows (%d) exceeds GSM capacity (786432). \n", num_rows);
        }
        return;
    }

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int warned_long_row = 0;

    // 按“行组”(row group) 为单位，处理分配给本线程的所有行
    while (current_row < tid_end)
    {
        int current_row_nnz = A_i[current_row + 1] - A_i[current_row];
        if (current_row_nnz > AM_size && !warned_long_row)
        {
            hthread_printf("[spmv_baseline warning] tid = %d, row = %d has %d nonzeros, "
                           "which exceeds AM_size = %d. This row may overflow A_cache/x_cache.\n",
                           tid, current_row, current_row_nnz, AM_size);
            warned_long_row = 1;
        }

        // ====================================================================
        // 步骤 1: 确定本次要处理的行组的边界
        // ====================================================================
        int group_end_row = current_row;
        int nnz_in_group = 0;
        while (group_end_row < tid_end)
        {
            int next_row_nnz = A_i[group_end_row + 1] - A_i[group_end_row];
            if (nnz_in_group + next_row_nnz > AM_size)
            {
                break;
            }
            nnz_in_group += next_row_nnz;
            group_end_row++;
        }
        if (group_end_row == current_row)
            group_end_row++; // 确保至少处理一行

        // ====================================================================
        // 【新增步骤】: 扫描非空行边界，并预处理空行
        // ====================================================================
        int group_first_nnz_row = -1;
        int group_last_nnz_row = -1;

        // 从前向后找到第一个非空行
        for (int i = current_row; i < group_end_row; i++)
        {
            if (A_i[i + 1] > A_i[i])
            {
                group_first_nnz_row = i;
                break;
            }
        }

        // 如果 group_first_nnz_row 未变，说明整个行组都是空的
        if (group_first_nnz_row == -1)
        {
            for (int i = current_row; i < group_end_row; i++)
            {
                // 假设 y=Ax，空行的结果就是0。
                // 如果是 y=alpha*Ax+beta*y，这里应为 y[i]*=beta
                y_data[i] = 0.0;
            }
        }
        else
        {
            // 如果存在非空行，则从后向前找到最后一个非空行
            for (int i = group_end_row - 1; i >= current_row; i--)
            {
                if (A_i[i + 1] > A_i[i])
                {
                    group_last_nnz_row = i;
                    break;
                }
            }

            // 预处理前导空行
            for (int i = current_row; i < group_first_nnz_row; i++)
            {
                y_data[i] = 0.0;
            }

            // 预处理尾随空行
            for (int i = group_last_nnz_row + 1; i < group_end_row; i++)
            {
                y_data[i] = 0.0;
            }

            // 定义收缩后的、纯净的计算区间
            int proc_start_row = group_first_nnz_row;
            int proc_end_row = group_last_nnz_row + 1; // 区间为 [start, end)

            // ====================================================================
            // 步骤 2: 在【收缩后的边界】上处理不对齐元素
            // ====================================================================
            int group_start_nnz_idx = A_i[proc_start_row];
            int group_end_nnz_idx = A_i[proc_end_row];

            double head_sum = 0.0; // 用于存储块头部不对齐元素的计算结果
            double tail_sum = 0.0; // 用于存储块尾部不对齐元素的计算结果

            int dma_start_idx = group_start_nnz_idx;
            int dma_end_idx = group_end_nnz_idx;

            if (dma_start_idx % 2 != 0)
            {
                double A_val = A_data[dma_start_idx];
                int column_index = A_j[dma_start_idx] / 8;
                double x_val = x_data[column_index];
                head_sum = A_val * x_val;
                dma_start_idx++;
            }

            if (dma_end_idx % 2 != 0)
            {
                dma_end_idx--;
                double A_val = A_data[dma_end_idx];
                int column_index = A_j[dma_end_idx] / 8;
                double x_val = x_data[column_index];
                tail_sum = A_val * x_val;
            }

            // ====================================================================
            // 步骤 3: 执行 dma_sg 从 GSM 加载到 AM (x_cache)
            // ====================================================================
            int dma_load_count = dma_end_idx - dma_start_idx;
            if (dma_load_count > 0)
            {
                // 1. 连续搬运 A_data 到 A_cache
                vector_load(&A_data[dma_start_idx], A_cache, dma_load_count * sizeof(double));

                // 2. 从 GSM 离散搬运 x 到 x_cache
                unsigned int ch = dma_sg(gsm_x,                           // 源基地址: GSM 缓冲区
                                         &A_j[dma_start_idx],             // 离散的字节偏移量 (已存 *8)
                                         1,                               // 每次读取 1 个元素
                                         dma_load_count * sizeof(double), // 总字节数
                                         0,                               // 标志位
                                         x_cache,                         // 目的地址：片上 AM
                                         1,                               // 连续存放
                                         dma_load_count * sizeof(double), // 存放总字节数
                                         0);                              // 标志位
                dma_wait(ch);
            }

            // ====================================================================
            // 步骤 4: 【行级别】在行组内逐行进行计算
            // ====================================================================
            for (int i = proc_start_row; i < proc_end_row; i++)
            {
                result = vec_svbcast(0.0);
                double sum = 0.0;

                // 因为边界已经收缩，这里的逻辑现在是正确的
                if (i == proc_start_row)
                {
                    sum += head_sum;
                }
                if (i == proc_end_row - 1)
                {
                    sum += tail_sum;
                }

                int row_nnz_start = A_i[i];
                int row_nnz_end = A_i[i + 1];

                // 中间的空行在这里会被自然跳过 (vec_length=0)
                if (row_nnz_end - row_nnz_start == 0)
                {
                    y_data[i] = sum; // 空行只应用边界值(如果恰好是边界)
                    continue;
                }

                if (row_nnz_start < dma_start_idx)
                    row_nnz_start = dma_start_idx;
                if (row_nnz_end > dma_end_idx)
                    row_nnz_end = dma_end_idx;

                int vec_length = row_nnz_end - row_nnz_start;
                if (vec_length > 0)
                {
                    int cache_offset = row_nnz_start - dma_start_idx;
                    int vec_iter = vec_length / 16;
                    int tail_len = vec_length % 16;

                    // --- 向量化主循环 ---
                    for (int j = 0; j < vec_iter; j++)
                    {
                        int element_offset = cache_offset + j * 16;
                        lvector double A_value = vec_ld(element_offset, A_cache);
                        lvector double x_value = vec_ld(element_offset, x_cache);
                        result = vec_mula(A_value, x_value, result);
                    }

                    if (tail_len > 0)
                    {
                        lvector double temp = vec_svbcast(0.0);
                        lvector double tail_result = vec_svbcast(0.0);
                        int tail_offset = cache_offset + vec_iter * 16;
                        lvector double A_value = vec_ld(tail_offset, A_cache);
                        lvector double x_value = vec_ld(tail_offset, x_cache);

                        tail_result = vec_muli(A_value, x_value);
                        vec_st(temp, tail_len, &tail_result);
                        result += tail_result;
                    }
                }

                y_data[i] = reduce_sum_fp64(result) + sum;
                // y_temp = reduce_sum_fp64(result) + sum;
                // y_data[i] = beta* y_data[i] + alpha*y_temp
            }
        }

        // 更新下一组的起始行
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// =====================================================
// =====================================================

// =====================================================
// x 预加载到 GSM
// 任务划分：    基于 AM 的细粒度划分 （v1 VS v2）
// x 加载方式：  vector_load
// 规约方式：    CPU端规约
// 任务执行方式: 顺序执行
// 论文中采用版本
// =====================================================
__global__ void SpMV_GSM_FP64_V2(int num_cores,
                                 int num_rows,
                                 int num_tasks,
                                 int barrier_id,
                                 int *A_i,
                                 int *A_j,
                                 double *A_data,
                                 double *x_data,
                                 double *y_data,
                                 int *row_bounds)
{
    int tid = get_thread_id();

    int AM_size = 96000 / 3;
    int step = AM_size;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *result_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result_vector = vec_svbcast(0.0);

    /* 预加载 x 到 GSM ------------------------------------------------------ */
    double *gsm_x = (double *)gx;
    int base = num_rows / num_cores;   // 每个核心基础份量
    int remain = num_rows % num_cores; // 余数，前 remain 个核心多拿 1 个
    int my_start = tid * base + (tid < remain ? tid : remain);
    int my_width = base + (tid < remain ? 1 : 0);
    if (my_start >= num_rows)
    {
        my_width = 0;
    }
    else if (my_start + my_width > num_rows)
    {
        my_width = num_rows - my_start;
    }

    if (my_width > 0 && num_rows <= 786432)
    {
        unsigned int ch = dma_p2p(&x_data[my_start],
                                  1, my_width * sizeof(double), 0,
                                  &gsm_x[my_start],
                                  1, my_width * sizeof(double), 0,
                                  false, 0);
        dma_wait(ch);
    }
    else if (num_rows > 786432 && tid == 0)
    {
        hthread_printf("Warning: num_rows (%d) exceeds GSM capacity (786432). \n", num_rows);
    }

    core_barrier(barrier_id, num_cores);
    /* ============== x 搬运完成 ============== */

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];
        int start_nnz = A_i[task_start_row];
        int end_nnz = A_i[task_end_row];
        int length = end_nnz - start_nnz;

        if (length <= 0)
            continue;

        int dma_load_start = start_nnz;
        int dma_load_end = end_nnz;
        // 处理头部不对齐
        if (dma_load_start % 2 != 0)
        {
            double A_val = A_data[dma_load_start];
            int col = A_j[dma_load_start] / 8;
            double x_val = gsm_x[col];              // ← 用预加载的 GSM
            y_data[dma_load_start] = A_val * x_val; // ← 结果直接写回 y_tmp
            dma_load_start++;
        }
        // 处理尾部不对齐
        if (dma_load_end % 2 != 0)
        {
            dma_load_end--;
            double A_val = A_data[dma_load_end];
            int col = A_j[dma_load_end] / 8;
            double x_val = gsm_x[col];
            y_data[dma_load_end] = A_val * x_val;
        }
        int dma_load_count = dma_load_end - dma_load_start;

        if (dma_load_count > 0)
        {
            // 连续搬运 A_data 到 A_cache
            vector_load(&A_data[dma_load_start], A_cache, dma_load_count * sizeof(double));

            // gather 对应的 x 值到连续缓冲
            unsigned int ch = dma_sg(gsm_x,                           // 源基地址: GSM 缓冲区
                                     &A_j[dma_load_start],            // 字节偏移量 (已 *8)
                                     1,                               // 每次读取 1 个元素
                                     dma_load_count * sizeof(double), // 总字节数
                                     0,                               // 标志位
                                     x_cache,                         // 目的地址：片上 AM
                                     1,                               // 连续存放
                                     dma_load_count * sizeof(double), // 存放总字节数
                                     0);                              // 标志位
            dma_wait(ch);
        }

        // 向量计算
        int vec_len = 16;
        int vec_iters = dma_load_count / vec_len;

        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * vec_len;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);

            // 向量乘法计算
            // result_vector = vec_svbcast(0.0);
            result_vector = vec_muli(A_value, x_value);

            // 将结果存回 AM 缓存
            vec_st(result_vector, offset, result_cache);
        }

        // 批量写回 DDR （注意：使用正确的 nnz 起始位置！）
        if (vec_iters > 0)
        {
            vector_store(result_cache,
                         &y_data[dma_load_start], // ← 用 dma_load_start（对齐后的 nnz 起始）
                         vec_iters * vec_len * sizeof(double));
        }

        // 处理剩余尾部
        int remain_count = dma_load_count % vec_len;
        if (remain_count > 0)
        {
            int remain_start = vec_iters * vec_len;
            for (int k = 0; k < remain_count; k++)
            {
                int global_idx = dma_load_start + remain_start + k;
                double A_scalar = A_data[global_idx];
                int col = A_j[global_idx] / 8;
                double x_scalar = gsm_x[col];
                y_data[global_idx] = A_scalar * x_scalar;
            }
        }
    }

    // 释放内存
    vector_free(A_cache);
    vector_free(x_cache);
    vector_free(result_cache);
}

__global__ void SpMV_GSM_FP64_V3(int num_tasks, // 任务划分总总段数
                                 int num_cores, // 调用 DSP 数量
                                 int N,         // 单个队列缓冲区个数 N
                                 int num_rows,
                                 int barrier_id,
                                 int *A_i,
                                 int *A_j,
                                 double *A_data,
                                 double *x_data,
                                 int *row_bounds,
                                 int *per_core_tasks,
                                 double *result_buffers,
                                 volatile long long *d_in,
                                 volatile long long *d_out)
{
    int tid = get_thread_id();

    // 定义常量，异步并发版本需与主机端保持一致
    int AM_size = 96000 / 3;
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *result_cache = vector_malloc(AM_size * sizeof(double));

    int my_num_tasks = per_core_tasks[tid]; // 从主机端获取每个核心的任务数
    int task_count = 0;

    /* 预加载 x 到 GSM ------------------------------------------------------ */
    double *gsm_x = (double *)gx;
    int base = num_rows / num_cores;   // 每个核心基础份量
    int remain = num_rows % num_cores; // 余数，前 remain 个核心多拿 1 个
    int my_start = tid * base + (tid < remain ? tid : remain);
    int my_width = base + (tid < remain ? 1 : 0);
    if (my_start >= num_rows)
    {
        my_width = 0;
    }
    else if (my_start + my_width > num_rows)
    {
        my_width = num_rows - my_start;
    }

    if (my_width > 0 && num_rows <= 786432)
    {
        unsigned int ch = dma_p2p(&x_data[my_start],
                                  1, my_width * sizeof(double), 0,
                                  &gsm_x[my_start],
                                  1, my_width * sizeof(double), 0,
                                  false, 0);
        dma_wait(ch);
    }
    else if (num_rows > 786432 && tid == 0)
    {
        hthread_printf("Warning: num_rows (%d) exceeds GSM capacity (786432). \n", num_rows);
    }

    core_barrier(barrier_id, num_cores);
    /* ========================== x 搬运完成 ========================== */

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        while (1)
        {
            long long in = d_in[tid];
            long long out = d_out[tid];
            if (in - out < N)
                break;

            volatile int dummy = 0;
            for (int i = 0; i < 36; i++)
            {
                dummy = dummy; // 强制读写，防优化
            }
        }

        // 2. 计算将要使用的缓冲区全局索引
        int local_idx = d_in[tid] % N;
        int global_idx = tid * N + local_idx;
        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];

        if (start_row >= end_row)
        {
            __sync_synchronize();
            d_in[tid]++;
            continue;
        }

        int start_nnz = A_i[start_row];
        int end_nnz = A_i[end_row];
        int length = end_nnz - start_nnz;
        if (length <= 0)
            continue;

        int write_offset = 0;
        int dma_load_start = start_nnz;
        int dma_load_end = end_nnz;

        // 处理头部不对齐
        if (dma_load_start % 2 != 0)
        {
            double A_val = A_data[dma_load_start];
            int col = A_j[dma_load_start] / 8;
            double x_val = gsm_x[col]; // ← 用预加载的 GSM
            current_buffer_ptr[0] = A_val * x_val;

            write_offset++;
            dma_load_start++;
        }
        // 处理尾部不对齐
        if (dma_load_end % 2 != 0)
        {
            dma_load_end--;
            double A_val = A_data[dma_load_end];
            int col = A_j[dma_load_end] / 8;
            double x_val = gsm_x[col];
            current_buffer_ptr[length - 1] = A_val * x_val;
        }
        int dma_load_count = dma_load_end - dma_load_start;

        if (dma_load_count > 0)
        {
            // 连续搬运 A_data 到 A_cache
            vector_load(&A_data[dma_load_start], A_cache, dma_load_count * sizeof(double));

            // gather 对应的 x 值到连续缓冲
            unsigned int ch = dma_sg(gsm_x,                           // 源基地址: GSM 缓冲区
                                     &A_j[dma_load_start],            // 字节偏移量 (已 *8)
                                     1,                               // 每次读取 1 个元素
                                     dma_load_count * sizeof(double), // 总字节数
                                     0,                               // 标志位
                                     x_cache,                         // 目的地址：片上 AM
                                     1,                               // 连续存放
                                     dma_load_count * sizeof(double), // 存放总字节数
                                     0);                              // 标志位
            dma_wait(ch);
        }

        // 向量计算
        int vec_len = 16;
        int vec_iters = dma_load_count / vec_len;

        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * vec_len;
            lvector double A_value = vec_ld(offset, A_cache);
            lvector double x_value = vec_ld(offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, offset, result_cache);
        }

        if (vec_iters > 0)
        {
            vector_store(result_cache,
                         &current_buffer_ptr[write_offset],
                         vec_iters * vec_len * sizeof(double));

            write_offset += (vec_iters * vec_len); // 更新局部偏移
        }

        // 处理剩余尾部
        int remain_count = dma_load_count % vec_len;
        if (remain_count > 0)
        {
            int remain_start = vec_iters * vec_len;
            for (int k = 0; k < remain_count; k++)
            {
                int global_idx = dma_load_start + remain_start + k;
                double A_scalar = A_data[global_idx];
                int col = A_j[global_idx] / 8;
                double x_scalar = gsm_x[col];

                current_buffer_ptr[write_offset] = A_scalar * x_scalar;
                write_offset++;
            }
        }

        current_buffer_ptr[AM_size] = (double)start_row;
        current_buffer_ptr[AM_size + 1] = (double)end_row;

        __sync_synchronize();
        d_in[tid]++;
        task_count++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
    vector_free(result_cache);
}