#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"

extern __gsm__ char gx[6 * 1024 * 1024]; // 6MB = 6 * 1024 * 1024 字节
const int FP32_block_size = 2;           // 1x2 BSR

// =================================================================
// 单精度：矩阵向量乘（SPU 版本）
// =================================================================
__global__ void kernel_SpMV_SPU_FP32(int num_rows,
                                     int num_cols,
                                     int rslice,
                                     int num_vectors,
                                     int *A_i,
                                     int *A_j,
                                     float *A_data,
                                     float *x_data,
                                     float *y_data,
                                     float *y_temp)
{
    int tid = get_thread_id();

    int startrow = tid * rslice < num_rows ? tid * rslice : num_rows;          // 线程起始行，限制不超过num_rows
    int stoprow = startrow + rslice < num_rows ? startrow + rslice : num_rows; // 线程结束行，限制不超过num_rows

    if (num_vectors == 1)
    {
        // 遍历行
        for (int i = startrow; i < stoprow; i++)
        {
            for (int jj = A_i[i]; jj < A_i[i + 1]; jj++)
            {
                int j = A_j[jj];
                if (j >= 0 && j < num_cols)
                {
                    y_temp[i] += A_data[jj] * x_data[j];
                }
            }
            y_data[i] += y_temp[i];
        }
    }
    else
    {
        hthread_printf("num_vectors > 1 \n");
    }
}

// 内核函数规约版本 ==============================
__global__ void kernel_dotSpMV_FP32(int *A_i,
                                    float *A_data,
                                    float *x_data,
                                    float *y_data,
                                    int *row_bounds)
{
    int tid = get_thread_id();
    // if (tid == 0)
    // {
    //     hthread_printf("kernel_dotSpMV_FP32 \n");
    // }

    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int AM_size = 96000;

    lvector float *A_cache = vector_malloc(AM_size * sizeof(float));
    lvector float *x_cache = vector_malloc(AM_size * sizeof(float));
    lvector double FP64_temp = vec_svbcast(1.0);
    lvector float result, tail_result, temp;
    float last_result = 0;

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
                // 如果加上下一行就超出缓存，则当前行组到此为
                // 注意：如果第一行就超缓存，需要特殊处理，保证至少处理一行
                // if (group_end_row == current_row)

                // {
                //    group_end_row++;
                // }
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
            vector_load(&A_data[group_start_nnz_idx], A_cache, nnz_to_load * sizeof(float));
            vector_load(&x_data[group_start_nnz_idx], x_cache, nnz_to_load * sizeof(float));
        }

        // 3. 内层循环：计算行组内的每一行 (不再需要 if-else)
        for (int i = current_row; i < group_end_row; i++)
        {
            last_result = 0.0f;
            result = vec_vsvbcast_sdf2(0.0);
            int start = A_i[i];
            int length = A_i[i + 1] - start;

            // 计算当前行数据在 cache 中的偏移量
            int cache_offset = start - group_start_nnz_idx;
            int vec_iter = length / 32;
            int tail_len = length % 32;
            for (int j = 0; j < vec_iter; j++)
            {
                int element_offset = (cache_offset + j * 32) / 2;
                lvector float A_value = vec_ld(element_offset, A_cache);
                lvector float x_value = vec_ld(element_offset, x_cache);
                result = vec_mula(A_value, x_value, result);
                // if (tid == 0 && i == 0)
                // {
                //     hthread_printf("result2 = %.12e ---\n", reduce_sum_float(result));
                // }
            }
            if (tail_len > 0)
            {
                int tail_offset = (cache_offset + vec_iter * 32) / 2;
                // if (tid == 0 && i == 0)
                // {
                //     hthread_printf("tail_len %d\n", tail_len);
                // }
                lvector float A_value = vec_ld(tail_offset, A_cache);
                lvector float x_value = vec_ld(tail_offset, x_cache);
                temp = vec_vsvbcast_sdf2(0.0);
                tail_result = vec_vsvbcast_sdf2(0.0);
                tail_result = vec_muli(A_value, x_value);
                vec_st(temp, (tail_len / 2), &tail_result);
                result += tail_result;
            }

            // y_data[i] = (float)reduce_sum_fp32(result); // float 规约 ==> double 累加
            // y_data[i] = reduce_sum_fp32(result);

            lvector double r0 = vec_fstdl(result);
            lvector double r1 = vec_fstdh(result);
            lvector double r = vec_mula(r0, FP64_temp, r1);
            y_data[i] = reduce_sum_fp64(r);
        }

        // 4. 更新 current_row，准备处理下一个行组
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void kernel_SpMV_FP32A(int *A_i,      // CSR row pointers
                                  float *A_data, // CSR non-zero values
                                  float *x_data, // 向量 x 的数据 (已预重排)
                                  float *y_data, // 输出向量 y
                                  int *row_bounds)
{
    int tid = get_thread_id();

    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int start = A_i[tid_start];
    int stop = A_i[tid_end];

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
        vector_load(&A_data[i], A_cache, length * sizeof(float));
        vector_load(&x_data[i], x_cache, length * sizeof(float));

        int vec_iters = length / 32;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;

            // 从 AM 缓存中取数
            lvector float A_value = vec_ld(offset, A_cache);
            lvector float x_value = vec_ld(offset, x_cache);

            // 向量乘法计算
            result_vector = vec_vsvbcast_sdf2(0.0);
            result_vector = vec_muli(A_value, x_value);

            // 将结果存回缓存
            vec_st(result_vector, offset, x_cache);
        }

        // 尾部标量处理
        int tail_len = length % 32;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 32; // 尾部在当前块内的偏移
            for (int k = 0; k < tail_len; k++)
            {
                int global_idx = i + remain_start + k;

                // 获取每个标量元素
                float A_scalar = A_data[global_idx];
                float x_scalar = x_data[global_idx];

                // 计算标量乘法
                float result_scalar = A_scalar * x_scalar;

                // 直接写入全局内存
                y_data[global_idx] = result_scalar;
            }
        }

        // vector_store(源地址, 目的地址, 需要写回的向量大小);
        if (vec_iters > 0)
        {
            vector_store(x_cache, &y_data[i], vec_iters * 32 * sizeof(float));
        }
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void spmv_baseline_FP32(int num_rows,
                                   int x_ext_size,
                                   int n_cores,
                                   int barrier_id,
                                   int *A_i,
                                   int *A_j,
                                   float *A_data,
                                   float *x_ext,
                                   float *y_data,
                                   int *core_row_bounds)
{
    int tid = get_thread_id();

    int AM_size = (96000 / 3) * 2;

    lvector float *A_cache = vector_malloc(AM_size * sizeof(float));
    lvector float *x_cache = vector_malloc(AM_size * sizeof(float));

    if (A_cache == NULL || x_cache == NULL)
    {
        if (tid == 0)
            hthread_printf("ERROR vector_malloc failed\n");
        return;
    }

    lvector float result;
    float *gsm_x = (float *)gx;

    /* 搬运 x_ext 到 GSM。x_ext[2*i] 是有效值，x_ext[2*i+1] 是 padding。 */
    int x_chunk_size = ((x_ext_size + n_cores - 1) / n_cores + 1) & ~1;
    int x_start_idx = tid * x_chunk_size;
    int x_width = 0;

    if (x_start_idx < x_ext_size)
    {
        x_width = x_chunk_size;
        if (x_start_idx + x_width > x_ext_size)
            x_width = x_ext_size - x_start_idx;
    }

    /* 如果 786432 是 double 槽位容量，则 FP32 x_ext 上限是 2 * 786432。 */
    if (x_ext_size <= 786432 * 2)
    {
        if (x_width > 0)
        {
            unsigned int c = dma_p2p(&x_ext[x_start_idx], 1, x_width * sizeof(float), 0,
                                     &gsm_x[x_start_idx], 1, x_width * sizeof(float), 0,
                                     false, 0);
            dma_wait(c);
        }

        core_barrier(barrier_id, n_cores);
    }
    else
    {
        if (tid == 0)
            hthread_printf("Warning: x_ext_size (%d) exceeds GSM capacity.\n", x_ext_size);

        vector_free(x_cache);
        vector_free(A_cache);
        return;
    }

    int current_row = core_row_bounds[tid];
    int tid_end = core_row_bounds[tid + 1];

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
                y_data[i] = 0.0f;
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
                y_data[i] = 0.0f;

            for (int i = group_last_nnz_row + 1; i < group_end_row; i++)
                y_data[i] = 0.0f;

            int proc_start_row = group_first_nnz_row;
            int proc_end_row = group_last_nnz_row + 1;

            int group_start_nnz_idx = A_i[proc_start_row];
            int group_end_nnz_idx = A_i[proc_end_row];

            float head_sum = 0.0f;
            float tail_sum = 0.0f;

            int dma_start_idx = group_start_nnz_idx;
            int dma_end_idx = group_end_nnz_idx;

            /*
             * FP32 使用 x_ext + 8B offset:
             * A_j[k] = original_col * 8.
             * gsm_x[A_j[k] / 4] == gsm_x[2 * original_col].
             */
            if (dma_start_idx % 2 != 0)
            {
                float A_val = A_data[dma_start_idx];
                int x_ext_index = A_j[dma_start_idx] / 4;
                float x_val = gsm_x[x_ext_index];

                head_sum = A_val * x_val;
                dma_start_idx++;
            }

            if (dma_end_idx % 2 != 0)
            {
                dma_end_idx--;

                float A_val = A_data[dma_end_idx];
                int x_ext_index = A_j[dma_end_idx] / 4;
                float x_val = gsm_x[x_ext_index];

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
                    float sum = 0.0f;

                    for (int jj = A_i[i]; jj < A_i[i + 1]; jj++)
                    {
                        int x_ext_index = A_j[jj] / 4;
                        sum += A_data[jj] * gsm_x[x_ext_index];
                    }

                    y_data[i] = sum;
                }

                current_row = group_end_row;
                continue;
            }

            if (dma_load_count > 0)
            {
                vector_load(&A_data[dma_start_idx],
                            A_cache,
                            dma_load_count * sizeof(float));

                /*
                 * A_j is byte offset in x_ext.
                 * dma_sg reads sizeof(float) bytes from each 8B slot start.
                 * x_cache becomes dense FP32 values.
                 */
                unsigned int ch = dma_sg(gsm_x,
                                         &A_j[dma_start_idx],
                                         1,
                                         dma_load_count * sizeof(float),
                                         0,
                                         x_cache,
                                         1,
                                         dma_load_count * sizeof(float),
                                         0);
                dma_wait(ch);
            }
            for (int i = proc_start_row; i < proc_end_row; i++)
            {
                result = vec_vsvbcast_sdf2(0.0f);
                float sum = 0.0f;

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

                int row_start = row_nnz_start;
                int row_end = row_nnz_end;

                /* FP32 vec_ld offset uses 8-byte units, so align row start/end to 2 floats. */
                if (row_start < row_end && (row_start % 2 != 0))
                {
                    int x_ext_index = A_j[row_start] / 4;
                    sum += A_data[row_start] * gsm_x[x_ext_index];
                    row_start++;
                }

                if (row_start < row_end && (row_end % 2 != 0))
                {
                    row_end--;
                    int x_ext_index = A_j[row_end] / 4;
                    sum += A_data[row_end] * gsm_x[x_ext_index];
                }

                int aligned_len = row_end - row_start;

                if (aligned_len > 0)
                {
                    int cache_offset = row_start - dma_start_idx;
                    int cache_offset_8b = cache_offset / 2;

                    int vec_len = 32;
                    int vec_iter = aligned_len / vec_len;
                    int tail_len = aligned_len % vec_len;

                    for (int j = 0; j < vec_iter; j++)
                    {
                        int offset = cache_offset_8b + j * 16;

                        lvector float A_value = vec_ld(offset, A_cache);
                        lvector float x_value = vec_ld(offset, x_cache);

                        result = vec_mula(A_value, x_value, result);
                    }

                    if (tail_len > 0)
                    {
                        int tail_offset_float = cache_offset + vec_iter * vec_len;
                        int tail_offset_8b = tail_offset_float / 2;

                        lvector float A_value = vec_ld(tail_offset_8b, A_cache);
                        lvector float x_value = vec_ld(tail_offset_8b, x_cache);

                        lvector float temp = vec_vsvbcast_sdf2(0.0f);
                        lvector float tail_result = vec_vsvbcast_sdf2(0.0f);

                        tail_result = vec_muli(A_value, x_value);
                        vec_st(temp, tail_len, &tail_result);

                        result += tail_result;
                    }
                }

                y_data[i] = reduce_sum_fp32(result) + sum;
            }
        }

        current_row = group_end_row;
    }

    vector_free(x_cache);
    vector_free(A_cache);
}

// ========================================================
__global__ void kernel_SpMV_FP32B(int coreNums,
                                  int num_tasks,
                                  int *A_i,
                                  float *A_data,
                                  float *x_data,
                                  int *row_bounds)
{
    int tid = get_thread_id();

    // AM 上限（cluster共享）要保守：这里沿用你FP16B风格
    int AM_elems = (96000 / 3) * 2; // 你原来的设定：约 64000 floats
    const int step = AM_elems;

    lvector float *A_cache = vector_malloc((size_t)AM_elems * sizeof(float));
    lvector float *x_cache = vector_malloc((size_t)AM_elems * sizeof(float));
    lvector float result_vector;

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];
        if (task_end_row <= task_start_row)
            continue;

        int start_nnz = A_i[task_start_row];
        int stop_nnz = A_i[task_end_row];
        if (stop_nnz <= start_nnz)
            continue;

        // 分段处理，确保每次 load/store 不超过 AM cache 容量
        for (int base = start_nnz; base < stop_nnz; base += step)
        {
            int end = (base + step < stop_nnz) ? (base + step) : stop_nnz;
            int length = end - base;
            if (length <= 0)
                continue;

            // // FP32: 16B 对齐 = 4 floats
            // int aligned = length & ~3;   // 向下取整到 4 的倍数
            // int tail = length - aligned; // 0..3

            // FP32: 8B 对齐 = 2 floats
            int aligned = length & ~1;   // 向下取整到 2 的倍数
            int tail = length - aligned; // 0..1

            // 只对齐部分走 vector_load/vector_store（避免 DMA/AM 越界）
            if (aligned > 0)
            {
                vector_load(&A_data[base], A_cache, (size_t)aligned * sizeof(float));
                vector_load(&x_data[base], x_cache, (size_t)aligned * sizeof(float));

                // 你原来是一组 32 个 float 做一次 vec_ld/vec_st
                int vec_iters = aligned / 32;

                for (int j = 0; j < vec_iters; j++)
                {
                    int offset = j * 16; // offset 单位是 16B block（你们平台约定）
                    lvector float A_value = vec_ld(offset, A_cache);
                    lvector float x_value = vec_ld(offset, x_cache);

                    result_vector = vec_vsvbcast_sdf2(0.0);
                    result_vector = vec_muli(A_value, x_value);

                    vec_st(result_vector, offset, x_cache);
                }

                if (vec_iters > 0)
                {
                    vector_store(x_cache, &x_data[base], (size_t)vec_iters * 32 * sizeof(float));
                }

                // 处理 aligned 里剩下但未满 32 的那一段（aligned % 32）
                int rem = aligned - vec_iters * 32; // 0..31
                for (int k = 0; k < rem; k++)
                {
                    int g = base + vec_iters * 32 + k;
                    x_data[g] = A_data[g] * x_data[g];
                }
            }

            // 尾巴 tail (<4) 直接标量
            for (int k = 0; k < tail; k++)
            {
                int g = base + aligned + k;
                x_data[g] = A_data[g] * x_data[g];
            }
        }
    }

    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void SpMV_GSM_FP32_BLOCK2(int num_rows,
                                     int x_pair_size,
                                     int num_cores,
                                     int num_tasks,
                                     int barrier_id,
                                     int *A_i_block,
                                     int *A_j_block,
                                     float *A_data_block,
                                     float *x_pair,
                                     float *dy_data,
                                     int *row_bounds)
{
    int tid = get_thread_id();

    int AM_size = (96000 / 3) * 2;
    int max_blocks_per_load = (AM_size / 2) & ~1;

    float *gsm_x = (float *)gx;

    int pair_count = x_pair_size / 2;
    int base_pair = pair_count / num_cores;
    int remain_pair = pair_count % num_cores;

    int my_pair_start = tid * base_pair + (tid < remain_pair ? tid : remain_pair);
    int my_pair_width = base_pair + (tid < remain_pair ? 1 : 0);

    int my_start = 2 * my_pair_start;
    int my_width = 2 * my_pair_width;

    if (my_start >= x_pair_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pair_size)
    {
        my_width = x_pair_size - my_start;
        my_width &= ~1;
    }

    if (x_pair_size <= 786432 * 2)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pair[my_start],
                                      1, my_width * sizeof(float), 0,
                                      &gsm_x[my_start],
                                      1, my_width * sizeof(float), 0,
                                      false, 0);
            dma_wait(ch);
        }

        core_barrier(barrier_id, num_cores);
    }
    else
    {
        if (tid == 0)
            hthread_printf("Warning: x_pair_size (%d) exceeds GSM capacity.\n", x_pair_size);
        return;
    }

    lvector float *A_cache = vector_malloc(AM_size * sizeof(float));
    lvector float *x_cache = vector_malloc(AM_size * sizeof(float));
    lvector float *result_cache = vector_malloc(AM_size * sizeof(float));

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];

        int start = A_i_block[task_start_row];
        int end = A_i_block[task_end_row];

        if (end <= start)
            continue;

        if (start & 1)
        {
            int b = start;
            int x_base = A_j_block[b] / 4;

            dy_data[2 * b] =
                A_data_block[2 * b] * gsm_x[x_base];

            dy_data[2 * b + 1] =
                A_data_block[2 * b + 1] * gsm_x[x_base + 1];

            start++;
        }

        if (start < end && (end & 1))
        {
            int b = end - 1;
            int x_base = A_j_block[b] / 4;

            dy_data[2 * b] =
                A_data_block[2 * b] * gsm_x[x_base];

            dy_data[2 * b + 1] =
                A_data_block[2 * b + 1] * gsm_x[x_base + 1];

            end--;
        }

        if (end <= start)
            continue;

        for (int base = start; base < end; base += max_blocks_per_load)
        {
            int stop = base + max_blocks_per_load;
            if (stop > end)
                stop = end;

            int num_block = stop - base;
            if (num_block <= 0)
                continue;

            int float_count = 2 * num_block;

            vector_load(&A_data_block[2 * base],
                        A_cache,
                        float_count * sizeof(float));

            unsigned int ch = dma_sg(gsm_x,
                                     &A_j_block[base],
                                     1,
                                     float_count * sizeof(float),
                                     0,
                                     x_cache,
                                     1,
                                     float_count * sizeof(float),
                                     0);
            dma_wait(ch);

            int vec_len = 32;
            int vec_iters = float_count / vec_len;

            for (int j = 0; j < vec_iters; j++)
            {
                int offset = j * 16;

                lvector float A_value = vec_ld(offset, A_cache);
                lvector float x_value = vec_ld(offset, x_cache);

                lvector float result_vector = vec_muli(A_value, x_value);

                vec_st(result_vector, offset, result_cache);
            }

            if (vec_iters > 0)
            {
                vector_store(result_cache,
                             &dy_data[2 * base],
                             vec_iters * vec_len * sizeof(float));
            }

            int remain_float = float_count - vec_iters * vec_len;
            if (remain_float > 0)
            {
                int start_float = vec_iters * vec_len;

                for (int k = 0; k < remain_float; k++)
                {
                    int global_float = 2 * base + start_float + k;
                    int block_id = global_float / 2;
                    int lane = global_float & 1;
                    int x_base = A_j_block[block_id] / 4;

                    dy_data[global_float] =
                        A_data_block[global_float] * gsm_x[x_base + lane];
                }
            }
        }
    }

    vector_free(result_cache);
    vector_free(x_cache);
    vector_free(A_cache);
}