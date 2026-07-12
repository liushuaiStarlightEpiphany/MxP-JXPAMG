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

__global__ void spmv_baseline(int num_rows,
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

    // if (tid == 0)
    //         hthread_printf("spmv_baseline (new)\n");

    int AM_size = 96000 / 3;

    lvector double *A_cache = vector_malloc(AM_size * sizeof(double));
    lvector double *x_cache = vector_malloc(AM_size * sizeof(double));

    if (A_cache == NULL || x_cache == NULL)
    {
        if (tid == 0)
            hthread_printf("ERROR vector_malloc failed\n");
        return;
    }

    lvector double result;

    double *gsm_x = (double *)gx;

    int x_chunk_size = ((num_rows + n_cores - 1) / n_cores + 1) & ~1;
    int x_start_idx = tid * x_chunk_size;
    int x_width = 0;

    if (x_start_idx < num_rows)
    {
        x_width = x_chunk_size;
        if (x_start_idx + x_width > num_rows)
            x_width = num_rows - x_start_idx;
    }

    if (num_rows <= 786432)
    {
        if (x_width > 0)
        {
            unsigned int c = dma_p2p(&x_data[x_start_idx], 1, x_width * sizeof(double), 0,
                                     &gsm_x[x_start_idx], 1, x_width * sizeof(double), 0,
                                     false, 0);
            dma_wait(c);
        }

        core_barrier(barrier_id, n_cores);
    }
    else
    {
        if (tid == 0)
            hthread_printf("Warning: num_rows (%d) exceeds GSM capacity (786432).\n", num_rows);

        vector_free(x_cache);
        vector_free(A_cache);
        return;
    }

    int current_row = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];

    if (current_row < 0 || tid_end < current_row || tid_end > num_rows)
    {
        hthread_printf("ERROR bad row_bounds tid=%d current_row=%d tid_end=%d num_rows=%d\n",
                       tid, current_row, tid_end, num_rows);
        vector_free(x_cache);
        vector_free(A_cache);
        return;
    }

    while (current_row < tid_end)
    {
        int group_end_row = current_row;
        int nnz_in_group = 0;

        while (group_end_row < tid_end)
        {
            int next_row_nnz = A_i[group_end_row + 1] - A_i[group_end_row];

            if (next_row_nnz < 0)
            {
                hthread_printf("ERROR bad A_i tid=%d row=%d A_i[row]=%d A_i[row+1]=%d\n",
                               tid, group_end_row, A_i[group_end_row], A_i[group_end_row + 1]);
                vector_free(x_cache);
                vector_free(A_cache);
                return;
            }

            if (nnz_in_group + next_row_nnz > AM_size)
                break;

            nnz_in_group += next_row_nnz;
            group_end_row++;
        }

        if (group_end_row == current_row)
            group_end_row++;

        if (group_end_row <= current_row || group_end_row > tid_end)
        {
            hthread_printf("ERROR no progress tid=%d current_row=%d group_end_row=%d tid_end=%d\n",
                           tid, current_row, group_end_row, tid_end);
            vector_free(x_cache);
            vector_free(A_cache);
            return;
        }

        int group_first_nnz_row = -1;
        int group_last_nnz_row = -1;

        for (int i = current_row; i < group_end_row; i++)
        {
            if (A_i[i + 1] > A_i[i])
            {
                group_first_nnz_row = i;
                break;
            }
        }

        if (group_first_nnz_row == -1)
        {
            for (int i = current_row; i < group_end_row; i++)
                y_data[i] = 0.0;
        }
        else
        {
            for (int i = group_end_row - 1; i >= current_row; i--)
            {
                if (A_i[i + 1] > A_i[i])
                {
                    group_last_nnz_row = i;
                    break;
                }
            }

            for (int i = current_row; i < group_first_nnz_row; i++)
                y_data[i] = 0.0;

            for (int i = group_last_nnz_row + 1; i < group_end_row; i++)
                y_data[i] = 0.0;

            int proc_start_row = group_first_nnz_row;
            int proc_end_row = group_last_nnz_row + 1;

            int group_start_nnz_idx = A_i[proc_start_row];
            int group_end_nnz_idx = A_i[proc_end_row];

            double head_sum = 0.0;
            double tail_sum = 0.0;

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

            int dma_load_count = dma_end_idx - dma_start_idx;

            if (dma_load_count < 0)
            {
                hthread_printf("ERROR negative dma_load_count tid=%d current_row=%d count=%d\n",
                               tid, current_row, dma_load_count);
                vector_free(x_cache);
                vector_free(A_cache);
                return;
            }

            if (dma_load_count > AM_size)
            {
                for (int i = proc_start_row; i < proc_end_row; i++)
                {
                    double sum = 0.0;

                    for (int jj = A_i[i]; jj < A_i[i + 1]; jj++)
                    {
                        int col = A_j[jj] / 8;
                        sum += A_data[jj] * x_data[col];
                    }

                    y_data[i] = sum;
                }

                current_row = group_end_row;
                continue;
            }

            if (dma_load_count > 0)
            {
                vector_load(&A_data[dma_start_idx], A_cache, dma_load_count * sizeof(double));

                unsigned int ch = dma_sg(gsm_x,
                                         &A_j[dma_start_idx],
                                         1,
                                         dma_load_count * sizeof(double),
                                         0,
                                         x_cache,
                                         1,
                                         dma_load_count * sizeof(double),
                                         0);
                dma_wait(ch);
            }

            for (int i = proc_start_row; i < proc_end_row; i++)
            {
                result = vec_svbcast(0.0);
                double sum = 0.0;

                if (i == proc_start_row)
                    sum += head_sum;

                if (i == proc_end_row - 1)
                    sum += tail_sum;

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

                    for (int j = 0; j < vec_iter; j++)
                    {
                        int element_offset = cache_offset + j * 16;

                        lvector double A_value = vec_ld(element_offset, A_cache);
                        lvector double x_value = vec_ld(element_offset, x_cache);

                        result = vec_mula(A_value, x_value, result);
                    }

                    if (tail_len > 0)
                    {
                        int tail_offset = cache_offset + vec_iter * 16;

                        lvector double A_value = vec_ld(tail_offset, A_cache);
                        lvector double x_value = vec_ld(tail_offset, x_cache);

                        lvector double temp = vec_svbcast(0.0);
                        lvector double tail_result = vec_svbcast(0.0);

                        tail_result = vec_muli(A_value, x_value);
                        vec_st(temp, tail_len, &tail_result);

                        result += tail_result;
                    }
                }

                y_data[i] = reduce_sum_fp64(result) + sum;
            }
        }

        current_row = group_end_row;
    }

    vector_free(x_cache);
    vector_free(A_cache);
}

// =====================================================
__global__ void SpMV_GSM_FP64(int num_rows,
                              int num_cores,
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

__global__ void SpMV_DOT_FP64(int coreNums,
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
    //     hthread_printf("SpMV_DOT_FP64 \n");
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