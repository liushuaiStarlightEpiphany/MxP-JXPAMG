#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"
#include <stdatomic.h>
#include <stdlib.h>

// =============================================================
// 内积版本三个内核函数，x 主机端预处理，vector_load 加载
// 版本一(A)：VPU 端按行计算，寄存器规约
// 版本二(B)：CPU 端规约，VPU->CPU 顺序执行（OMP 并行情形下效果最优）
// 版本三(C)：CPU 端规约，VPU/CPU 异步并发
// =============================================================
// 内积版本 SPMV ===> 寄存器规约
__global__ void kernel_dotSpMV_FP64A(int *A_i,       // CSR row pointers
                                     double *A_data, // CSR non-zero values
                                     double *x_data, // 向量 x 的数据 (已预重排)
                                     double *y_data, // 输出向量 y
                                     int *row_bounds)
{
    int tid = get_thread_id();
    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    int AM_size = 96000 / 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result, tail_result, temp;

    // 外层循环：以“行组”为单位进行迭代
    int current_row = tid_start;
    while (current_row < tid_end)
    {
        // 1. 确定下一个行组的边界 (group_end_row)
        int group_end_row = current_row;
        int nnz_in_group = 0;
        while (group_end_row < tid_end)
        {
            // 查看下一行的非零元个数
            int next_row_nnz = A_i[group_end_row + 1] - A_i[group_end_row];
            if (nnz_in_group + next_row_nnz > AM_size)
            {
                break;
            }
            nnz_in_group += next_row_nnz;
            group_end_row++;
        }

        // 2. 加载整个行组的数据
        int group_start_nnz_idx = A_i[current_row];
        int nnz_to_load = A_i[group_end_row] - group_start_nnz_idx;

        if (nnz_to_load > 0)
        {
            vector_load(&A_data[group_start_nnz_idx], A_cache, nnz_to_load * sizeof(double));
            vector_load(&x_data[group_start_nnz_idx], x_cache, nnz_to_load * sizeof(double));
        }

        // 3. 内层循环：计算行组内的每一行
        for (int i = current_row; i < group_end_row; i++)
        {
            result = vec_svbcast(0.0);
            int start = A_i[i];
            int length = A_i[i + 1] - start;

            // 计算当前行数据在 cache 中的偏移量
            int cache_offset = start - group_start_nnz_idx;
            int vec_iter = length / 16;
            int tail_len = length % 16;

            for (int j = 0; j < vec_iter; j++)
            {
                int element_offset = cache_offset + j * 16;
                lvector double A_value = vec_ld(element_offset, A_cache);
                lvector double x_value = vec_ld(element_offset, x_cache);
                result = vec_mula(A_value, x_value, result);
            }
            if (tail_len > 0)
            {
                temp = vec_svbcast(0.0);
                tail_result = vec_svbcast(0.0);
                int tail_offset = cache_offset + vec_iter * 16;
                lvector double A_value = vec_ld(tail_offset, A_cache);
                lvector double x_value = vec_ld(tail_offset, x_cache);

                tail_result = vec_muli(A_value, x_value);
                vec_st(temp, tail_len, &tail_result);
                result += tail_result;
            }

            y_data[i] = reduce_sum_fp64(result);
        }

        current_row = group_end_row;
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
__global__ void kernel_v1_mix(int *A_i,
                              double *A_data,
                              double *x_data,
                              double *y_data,
                              int *row_bounds)
{
    int tid = get_thread_id();

    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int start = A_i[tid_start];
    int stop = A_i[tid_end];

    // double 情形
    // int AM_size = 96000 / 2;
    // float 情形
    int AM_size = 96000;
    int step = AM_size;
    lvector float *A_cache = vector_malloc(AM_size * sizeof(float));
    lvector float *x_cache = vector_malloc(AM_size * sizeof(float));
    lvector float result_vector;

    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&A_data[i], A_cache, length * sizeof(double));
        vector_load(&x_data[i], x_cache, length * sizeof(double));

        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;

            // 从 AM 缓存中取数
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
            vector_store(x_cache, &y_data[i], vec_iters * 16 * sizeof(double));
        }

        // 尾部（长度不足16个元素）标量处理
        int tail_len = length % 16;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 16; // 尾部在当前块内的偏移
            for (int k = 0; k < tail_len; k++)
            {
                int global_idx = i + remain_start + k;

                // 获取每个标量元素
                double A_scalar = A_data[global_idx];
                double x_scalar = x_data[global_idx];

                // 计算标量乘法
                double result_scalar = A_scalar * x_scalar;

                // 直接写入全局内存
                y_data[global_idx] = result_scalar;
            }
        }
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// =====================================================
// 将向量 x 延长, SPMV 转为 DOT 计算
// 任务划分：    基于 AM 的细粒度划分 （v1 VS v2）
// x 加载方式：  vector_load
// 规约方式：    CPU端规约
// 任务执行方式: 顺序执行
// =====================================================
__global__ void kernel_v2(int coreNums,
                          int num_tasks,
                          int *A_i,
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
    lvector double result_vector;

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {

        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];
        int start_nnz = A_i[task_start_row];
        int stop_nnz = A_i[task_end_row];
        int length = stop_nnz - start_nnz;

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

// ============================================================================
// 分隔符 ======================================================================

// base line 版本 ========================
// A: vector load
// x: dma_sg
// SVR 寄存器规约
// 考虑了矩阵中存在零行的情况，增加对零行的判断，确保没对齐的元素能够正确累加
__global__ void spmv_baseline(int *A_i,
                              int *A_j,
                              double *A_data,
                              double *x_data,
                              double *y_data,
                              int *row_bounds)
{
    int tid = get_thread_id();

    int AM_size = 96000 / 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result;

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    // 按“行组”(row group) 为单位，处理分配给本线程的所有行
    while (current_row < tid_end)
    {
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
                // 假设 y=Ax，空行的结果就是 0
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
            }
        }

        // 更新下一组的起始行
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// 与 baseline 的区别：规约方式
// baseline：二叉树形式，先按照 16 进行向量加，再 reduce
// baseline_B：16 个为批次进行 reduce，然后进行累加
// 论文中采取这个版本 *********************************
__global__ void spmv_baseline_B(int *A_i,
                                int *A_j,
                                double *A_data,
                                double *x_data,
                                double *y_data,
                                int *row_bounds)
{
    int tid = get_thread_id();

    // if (tid == 0)
    // {
    //     hthread_printf("spmv_baseline_B \n");
    // }

    // 初始化本地缓存和向量寄存器
    int AM_size = 96000 / 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));
    lvector double result;

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    // 按“行组”(row group) 为单位，处理分配给本线程的所有行
    while (current_row < tid_end)
    {
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
                        result = vec_muli(A_value, x_value);
                        y_data[i] += reduce_sum_fp64(result);
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
                        y_data[i] += reduce_sum_fp64(tail_result);
                    }
                }

                y_data[i] = y_data[i] + sum;
            }
        }

        // 更新下一组的起始行
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// ====================================================================
// 内积版本 SPMV, 无需对 x 做预处理
// DSP 端计算，CPU 端规约，DSP/CPU 顺序执行
// A: DDR->AM(vector_load)    x: DDR->AM(dma_sg)
// ====================================================================
__global__ void spmv_fp64_v1(int *A_i,
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

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    // 按“行组”(row group) 为单位，处理分配给本线程的所有行
    while (current_row < tid_end)
    {
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
        // 步骤 2: 【块级别】处理整个行组的边界不对齐元素
        // ====================================================================

        int group_start_nnz_idx = A_i[current_row];
        int group_end_nnz_idx = A_i[group_end_row];
        int dma_start_idx = group_start_nnz_idx;
        int dma_end_idx = group_end_nnz_idx;

        // --- 处理块的头部 (仅当起始地址为奇数时) ---
        if (dma_start_idx % 2 != 0)
        {
            double A_val = A_data[dma_start_idx];
            int column_index = A_j[dma_start_idx] / 8;
            double x_val = x_data[column_index];
            y_data[dma_start_idx] = A_val * x_val;
            dma_start_idx++;
        }

        // --- 处理块的尾部 (仅当结束地址为奇数时) ---
        if (dma_end_idx % 2 != 0)
        {
            dma_end_idx--; // DMA加载区域向前移动一位
            double A_val = A_data[dma_end_idx];
            int column_index = A_j[dma_end_idx] / 8;
            double x_val = x_data[column_index];
            y_data[dma_end_idx] = A_val * x_val;
        }

        // ====================================================================
        // 步骤 3: 【块级别】执行对齐、偶数长度的DMA加载
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
        int vec_iters = dma_load_count / 16;
        int tail_len = dma_load_count % 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int element_offset = j * 16;
            lvector double A_value = vec_ld(element_offset, A_cache);
            lvector double x_value = vec_ld(element_offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, element_offset, x_cache);
        }

        // ====================================================================
        // 步骤 5: 【写回】将计算好的向量化结果一次性写回全局内存 AM->DDR
        // ====================================================================
        if (vec_iters > 0)
        {
            vector_store(x_cache, &y_data[dma_start_idx], vec_iters * 16 * sizeof(double));
        }

        if (tail_len > 0)
        {
            int scalar_start_offset = vec_iters * 16;
            for (int k = scalar_start_offset; k < dma_load_count; k++)
            {
                int global_idx = dma_start_idx + k;
                double A_val = A_data[global_idx];
                int column_index = A_j[global_idx] / 8;
                double x_val = x_data[column_index];
                y_data[global_idx] = A_val * x_val;
            }
        }

        // 更新下一组的起始行
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

// ====================================================================
// 内积版本 SPMV, 无需对 x 做预处理
// DSP 端计算，CPU 端规约，DSP/CPU 异步并发
// A: DDR->AM(vector_load)    x: DDR->AM(dma_sg)
// ====================================================================
__global__ void spmv_fp64_v2(int num_tasks,            // 任务划分总总段数
                             int num_cores,            // 调用 DSP 数量
                             int N,                    // 单个队列缓冲区个数 N
                             int *A_i,                 // CSR A_i
                             int *A_j,                 // CSR A_j
                             double *A_data,           // CSR A_data
                             double *x_data,           // 向量 x_data
                             double *result_buffers,   // 结果 y_data
                             int *row_bounds,          // 任务划分数组
                             volatile long long *d_in, // 同步指针
                             volatile long long *d_out)
{
    int tid = get_thread_id();
    if (tid == 1)
    {
        hthread_printf("spmv_fp64_v2 \n");
    }

    // 初始化本地缓存和向量寄存器
    int AM_size = 96000 / 3;
    // 单个缓冲区大小
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        // 1. 检查自己的队列是否已满
        while (d_in[tid] - d_out[tid] >= N)
        {
            // 队列已满，原地等待消费者处理
        }

        // 2. 计算将要使用的缓冲区全局索引 (使用参数 N)
        int local_idx = d_in[tid] % N;
        int global_idx = tid * N + local_idx;

        // 3. 获取当前任务的计算范围
        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];

        // 获取当前任务的非零值起始全局索引，用于计算相对偏移
        int start_nnz_of_task = A_i[start_row];
        int stop_nnz_of_task = A_i[end_row];

        int dma_start_idx = start_nnz_of_task;
        int dma_end_idx = stop_nnz_of_task;

        // --- 处理块的头部 (仅当起始地址为奇数时) ---
        if (dma_start_idx % 2 != 0)
        {
            double A_val = A_data[dma_start_idx];
            int column_index = A_j[dma_start_idx] / 8;
            double x_val = x_data[column_index];

            int offset_in_buffer = dma_start_idx - start_nnz_of_task;
            current_buffer_ptr[offset_in_buffer] = A_val * x_val;

            dma_start_idx++;
        }

        // --- 处理块的尾部 (仅当结束地址为奇数时) ---
        if (dma_end_idx % 2 != 0)
        {
            dma_end_idx--; // DMA加载区域向前移动一位
            double A_val = A_data[dma_end_idx];
            int column_index = A_j[dma_end_idx] / 8;
            double x_val = x_data[column_index];

            int offset_in_buffer = dma_end_idx - start_nnz_of_task;
            current_buffer_ptr[offset_in_buffer] = A_val * x_val;
        }
        // ====================================================================
        // 步骤 3: 【块级别】执行对齐、偶数长度的DMA加载
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
        int vec_iters = dma_load_count / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int element_offset = j * 16;
            lvector double A_value = vec_ld(element_offset, A_cache);
            lvector double x_value = vec_ld(element_offset, x_cache);
            lvector double result_vector = vec_muli(A_value, x_value);
            vec_st(result_vector, element_offset, x_cache);
        }

        // ====================================================================
        // 步骤 5: 【写回】将计算好的向量化结果一次性写回全局内存 AM->DDR
        // ====================================================================
        if (vec_iters > 0)
        {
            int offset_in_buffer = dma_start_idx - start_nnz_of_task;
            vector_store(x_cache, current_buffer_ptr + offset_in_buffer, vec_iters * 16 * sizeof(double));
        }

        // ====================================================================
        // 步骤 6: 【尾部】行快内尾部不足 16 的元素标量处理
        // ====================================================================
        int tail_len = dma_load_count % 16;
        if (tail_len > 0)
        {
            int scalar_start_offset = vec_iters * 16;
            for (int k = scalar_start_offset; k < dma_load_count; k++)
            {
                int global_idx = dma_start_idx + k;
                double A_val = A_data[global_idx];
                int column_index = A_j[global_idx] / 8;
                double x_val = x_data[column_index];
                int offset_in_buffer = global_idx - start_nnz_of_task;
                current_buffer_ptr[offset_in_buffer] = A_val * x_val;
            }
        }
        // 写入元数据，供CPU消费者使用
        current_buffer_ptr[AM_size] = (double)start_row;
        current_buffer_ptr[AM_size + 1] = (double)end_row;

        __sync_synchronize();
        d_in[tid]++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
}

// CPU/DSP 并发执行版本 ========================================
// 论文中采用的版本 ********************************************
__global__ void kernel_dotSpMV_FP64C(int num_tasks, // 任务划分总总段数
                                     int num_cores, // 调用 DSP 数量
                                     int N,         // 单个队列缓冲区个数 N
                                     int *A_i,
                                     double *A_data,
                                     double *x_data,
                                     int *row_bounds,
                                     int *tasks_per_core,
                                     double *result_buffers,
                                     volatile long long *d_in,
                                     volatile long long *d_out)
{
    int tid = get_thread_id();

    // hthread_printf("DSP Core %d reporting for duty! (Kernel has started)\n", tid);
    // if (tid == 0)
    // {
    //     hthread_printf("kernel_dotSpMV_FP64E \n");
    // }

    // 定义常量，异步并发版本需与主机端保持一致
    int AM_size = 96000 / 3;
    int single_buffer_stride = AM_size + 2;
    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    int my_num_tasks = tasks_per_core[tid]; // 从主机端获取每个核心的任务数
    int task_count = 0;

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        // 1. 检查自己的队列是否已满
        while (d_in[tid] - d_out[tid] >= N)
        {
            // 队列已满，原地等待消费者处理
        }

        // 2. 计算将要使用的缓冲区全局索引 (使用参数 N)
        int local_idx = d_in[tid] % N;
        int global_idx = tid * N + local_idx;

        double *current_buffer_ptr = result_buffers + global_idx * single_buffer_stride;
        int start_row = row_bounds[task_id];
        int end_row = row_bounds[task_id + 1];

        if (start_row >= end_row)
        {
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

        // __sync_synchronize();
        d_in[tid]++;
        task_count++;
    }

    vector_free(A_cache);
    vector_free(x_cache);
}