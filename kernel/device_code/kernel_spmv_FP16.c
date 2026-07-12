#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"

extern __gsm__ char gx[6 * 1024 * 1024]; // 6MB = 6 * 1024 * 1024 字节

// ==========================================================
// ==========================================================
// 待补充
// ==========================================================

// ==========================================================
// 基于 dot 思想的 spmv ======================================
// ==========================================================
__global__ void kernel_dotSpMV_FP16(int *A_i,
                                    __fp16 *A_data,
                                    __fp16 *x_data,
                                    __fp16 *y_data,
                                    int *row_bounds)
{
    int tid = get_thread_id();
    // if (tid == 0)
    // {
    //     hthread_printf("kernel_dotSpMV_FP16 \n");
    // }

    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int AM_size = 96000 * 2;

    lvector __fp16 *A_cache = vector_malloc(AM_size * sizeof(__fp16));
    lvector __fp16 *x_cache = vector_malloc(AM_size * sizeof(__fp16));
    lvector double FP64_temp = vec_svbcast(0.0);
    lvector double FP64_one = vec_svbcast(1.0);
    lvector __fp16 result, tail_result, temp;
    __fp16 last_result = 0;

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
            vector_load(&A_data[group_start_nnz_idx], A_cache, nnz_to_load * sizeof(__fp16));
            vector_load(&x_data[group_start_nnz_idx], x_cache, nnz_to_load * sizeof(__fp16));
        }

        // 3. 内层循环：计算行组内的每一行 (不再需要 if-else)
        for (int i = current_row; i < group_end_row; i++)
        {
            last_result = 0.0f; // 重置 last_result
            result = vec_abs((lvector __fp16)FP64_temp);
            int start = A_i[i];
            int length = A_i[i + 1] - start;

            // 计算当前行数据在 cache 中的偏移量
            int cache_offset = start - group_start_nnz_idx;
            int vec_iter = length / 64;
            int tail_len = length % 64;

            for (int j = 0; j < vec_iter; j++)
            {
                int element_offset = (cache_offset + j * 64) / 4;
                lvector __fp16 A_value = vec_ld(element_offset, A_cache);
                lvector __fp16 x_value = vec_ld(element_offset, x_cache);
                result = vec_mula(A_value, x_value, result);
            }
            if (tail_len > 0)
            {
                int tail_offset = (cache_offset + vec_iter * 64) / 4;
                lvector __fp16 A_value = vec_ld(tail_offset, A_cache);
                lvector __fp16 x_value = vec_ld(tail_offset, x_cache);
                temp = vec_abs((lvector __fp16)FP64_temp);
                tail_result = vec_abs((lvector __fp16)FP64_temp);
                tail_result = vec_muli(A_value, x_value);
                vec_st(temp, (tail_len / 4), &tail_result);
                result += tail_result;
            }
            // fp16 规约 ==> double 累加
            // y_data[i] = (float)reduce_sum_fp16(result);

            lvector float f0 = vec_fstdl(result);
            lvector double r0 = vec_fstdl(f0);
            lvector double r1 = vec_fstdh(f0);

            lvector float f1 = vec_fstdh(result);
            lvector double r2 = vec_fstdl(f1);
            lvector double r3 = vec_fstdh(f1);

            lvector double s0 = vec_mula(r0, FP64_one, r1);      // s0 = r0 + r1
            lvector double s1 = vec_mula(r2, FP64_one, r3);      // s1 = r2 + r3
            lvector double r_total = vec_mula(s0, FP64_one, s1); // r_total = (r0+r1) + (r2+r3)

            y_data[i] = reduce_sum_fp64(r_total);
        }
        // 4. 更新 current_row，准备处理下一个行组
        current_row = group_end_row;
    }

    // 释放内存
    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void kernel_SpMV_FP16A(int *A_i,       // CSR row pointers
                                  __fp16 *A_data, // CSR non-zero values
                                  __fp16 *x_data, // 向量 x 的数据 (已预重排)
                                  __fp16 *y_data, // 输出向量 y
                                  int *row_bounds)
{
    int tid = get_thread_id();

    int tid_start = row_bounds[tid];
    int tid_end = row_bounds[tid + 1];
    int start = A_i[tid_start];
    int stop = A_i[tid_end];

    int AM_size = 96000 * 2;
    int step = AM_size;
    lvector __fp16 *A_cache = vector_malloc(AM_size * sizeof(__fp16));
    lvector __fp16 *x_cache = vector_malloc(AM_size * sizeof(__fp16));
    lvector __fp16 result_vector;
    lvector double FP64_temp = vec_svbcast(0.0);

    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&A_data[i], A_cache, length * sizeof(__fp16));
        vector_load(&x_data[i], x_cache, length * sizeof(__fp16));

        int vec_iters = length / 64;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;

            // 从 AM 缓存中取数
            lvector __fp16 A_value = vec_ld(offset, A_cache);
            lvector __fp16 x_value = vec_ld(offset, x_cache);

            // 向量乘法计算
            result_vector = vec_abs((lvector __fp16)FP64_temp);
            result_vector = vec_muli(A_value, x_value);

            // 将结果存回缓存
            vec_st(result_vector, offset, x_cache);
        }

        // vector_store(源地址, 目的地址, 需要写回的向量大小);
        if (vec_iters > 0)
        {
            vector_store(x_cache, &y_data[i], vec_iters * 64 * sizeof(__fp16));
        }

        // 尾部标量处理
        int tail_len = length % 64;
        if (tail_len > 0)
        {
            int remain_start = vec_iters * 64; // 尾部在当前块内的偏移
            for (int k = 0; k < tail_len; k++)
            {
                int global_idx = i + remain_start + k;

                // 获取每个标量元素
                __fp16 A_scalar = A_data[global_idx];
                __fp16 x_scalar = x_data[global_idx];

                // 计算标量乘法
                __fp16 result_scalar = A_scalar * x_scalar;

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
// FP16 版本：与 kernel_SpMV_FP32B 对应
// - 任务划分：基于 AM 的细粒度划分（row_bounds/task）
// - 加载方式：vector_load 到 AM
// - 计算：AM 内 vec_muli(A, x) 写回 AM
// - 写回：vector_store 回 x_data（或 y_data）
// - 规约：CPU 端按 row_sum 做
// =====================================================
// __global__ void kernel_SpMV_FP16B(int coreNums,
//                                   int num_tasks,
//                                   int *A_i,
//                                   __fp16 *A_data,
//                                   __fp16 *x_data,
//                                   int *row_bounds)
// {
//     int tid = get_thread_id();

//     if (tid == 0)
//     {
//         hthread_printf("kernel_SpMV_FP16B \n");
//     }

//     // 注意：AM 上限 768KB（cluster 共享），这里的 AM_size 必须非常保守
//     int AM_elems = (96000 / 3) * 4;
//     const int step = AM_elems;

//     lvector __fp16 *A_cache = vector_malloc(AM_elems * sizeof(__fp16));
//     lvector __fp16 *x_cache = vector_malloc(AM_elems * sizeof(__fp16));
//     lvector __fp16 result_vector;

//     for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
//     {
//         int task_start_row = row_bounds[task_id];
//         int task_end_row = row_bounds[task_id + 1];

//         // 防御：空任务
//         if (task_end_row <= task_start_row)
//             continue;

//         int start_nnz = A_i[task_start_row];
//         int stop_nnz = A_i[task_end_row];

//         // 分段处理，确保每次 load/store 不超过 AM cache 容量
//         for (int base = start_nnz; base < stop_nnz; base += step)
//         {
//             int end = (base + step < stop_nnz) ? (base + step) : stop_nnz;
//             int length = end - base;
//             if (length <= 0)
//                 continue;

//             if (tid == 0 && task_id == 0)
//             {
//                 size_t bytes = (size_t)length * sizeof(__fp16);
//                 hthread_printf("base=%d len=%d bytes=%lu bytes/16=%lu\n",
//                                base, length, (unsigned long)bytes, (unsigned long)(bytes >> 4));
//             }

//             // load 到 AM
//             vector_load(&A_data[base], A_cache, length * sizeof(__fp16));
//             vector_load(&x_data[base], x_cache, length * sizeof(__fp16));

//             // // ===== 向量乘法 =====
//             // int vec_iters = length / 64;
//             // for (int j = 0; j < vec_iters; j++)
//             // {
//             //     int offset = j * 16;

//             //     lvector __fp16 A_value = vec_ld(offset, A_cache);
//             //     lvector __fp16 x_value = vec_ld(offset, x_cache);

//             //     // 乘法
//             //     result_vector = vec_muli(A_value, x_value);

//             //     // 写回 AM
//             //     vec_st(result_vector, offset, x_cache);
//             // }

//             // // store 回全局（只写回完整向量部分）
//             // if (vec_iters > 0)
//             // {
//             //     vector_store(x_cache, &x_data[base], vec_iters * 64 * sizeof(__fp16));
//             // }

//             // // ===== 尾部标量处理 =====
//             // int tail_len = length % 64;
//             // if (tail_len > 0)
//             // {
//             //     int remain_start = vec_iters * 64;
//             //     for (int k = 0; k < tail_len; k++)
//             //     {
//             //         int g = base + remain_start + k;
//             //         __fp16 a = A_data[g];
//             //         __fp16 x = x_data[g];
//             //         x_data[g] = a * x;
//             //     }
//             // }
//         }
//     }

//     vector_free(x_cache);
//     vector_free(A_cache);
// }

// __global__ void kernel_SpMV_FP16B(int coreNums,
//                                   int num_tasks,
//                                   int *A_i,
//                                   __fp16 *A_data,
//                                   __fp16 *x_data,
//                                   int *row_bounds)
// {
//     int tid = get_thread_id();

//     // if (tid == 0)
//     // {
//     //     hthread_printf("kernel_SpMV_FP16B \n");
//     // }

//     const int AM_elems = (96000 / 3) * 4; // 64000, 4的倍数
//     const int step = AM_elems;

//     lvector __fp16 *A_cache = vector_malloc(AM_elems * sizeof(__fp16));
//     lvector __fp16 *x_cache = vector_malloc(AM_elems * sizeof(__fp16));
//     lvector __fp16 result_vector;

//     if (!A_cache || !x_cache)
//     {
//         if (A_cache)
//             vector_free(A_cache);
//         if (x_cache)
//             vector_free(x_cache);
//         return;
//     }

//     for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
//     {
//         int task_start_row = row_bounds[task_id];
//         int task_end_row = row_bounds[task_id + 1];
//         if (task_end_row <= task_start_row)
//             continue;

//         int start_nnz = A_i[task_start_row];
//         int stop_nnz = A_i[task_end_row];
//         if (stop_nnz <= start_nnz)
//             continue;

//         // ---- 关键1：先把起点推到 4 对齐（保证FP16地址8B对齐）----
//         int base = start_nnz;
//         while (base < stop_nnz && (base & 3) != 0)
//         {
//             x_data[base] = A_data[base] * x_data[base];
//             base++;
//         }

//         for (; base < stop_nnz; base += step)
//         {
//             int end = (base + step < stop_nnz) ? (base + step) : stop_nnz;
//             int length = end - base;
//             if (length <= 0)
//                 continue;

//             // ---- 关键2：DMA长度也按 4 对齐（保证字节数8B倍数）----
//             int body = (length / 4) * 4; // body个fp16 => body*2字节为8B倍数
//             if (body > 0)
//             {
//                 vector_load(&A_data[base], A_cache, body * sizeof(__fp16));
//                 vector_load(&x_data[base], x_cache, body * sizeof(__fp16));

//                 // 向量部分：按64元素做
//                 int vec_elems = (body / 64) * 64;
//                 int vec_iters = vec_elems / 64;

//                 for (int j = 0; j < vec_iters; j++)
//                 {
//                     int offset = j * 16;
//                     lvector __fp16 A_value = vec_ld(offset, A_cache);
//                     lvector __fp16 x_value = vec_ld(offset, x_cache);
//                     result_vector = vec_muli(A_value, x_value);
//                     vec_st(result_vector, offset, x_cache);
//                 }

//                 // vec_elems 之后（但仍在 body 内）的部分，标量算（避免漏算）
//                 for (int g = base + vec_elems; g < base + body; g++)
//                 {
//                     x_data[g] = A_data[g] * x_data[g];
//                 }

//                 vector_store(x_cache, &x_data[base], vec_elems * sizeof(__fp16));
//                 // 注意：上面只store了vec_elems部分，因为只有这部分写回到x_cache了
//                 // 如果你想把body全store，需要把标量部分也写回x_cache（不划算）
//             }

//             // 尾部标量（body之外）
//             for (int g = base + body; g < end; g++)
//             {
//                 x_data[g] = A_data[g] * x_data[g];
//             }
//         }
//     }

//     vector_free(x_cache);
//     vector_free(A_cache);
// }

__global__ void kernel_SpMV_FP16B(int coreNums,
                                  int num_tasks,
                                  int *A_i,
                                  __fp16 *A_data,
                                  __fp16 *x_data,
                                  int *row_bounds)
{
    int tid = get_thread_id();

    // if (tid == 0)
    // {
    //     hthread_printf("kernel_SpMV_FP16B \n");
    // }

    // 注意：AM 上限 768KB（cluster 共享），这里的 AM_size 必须非常保守
    int AM_elems = (96000 / 3) * 4;
    const int step = AM_elems;

    lvector __fp16 *A_cache = vector_malloc(AM_elems * sizeof(__fp16));
    lvector __fp16 *x_cache = vector_malloc(AM_elems * sizeof(__fp16));
    lvector __fp16 result_vector;

    for (int task_id = tid; task_id < num_tasks; task_id += coreNums)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];
        if (task_end_row <= task_start_row)
            continue;

        int start_nnz = A_i[task_start_row];
        int stop_nnz = A_i[task_end_row];

        // 分段处理，确保每次 load/store 不超过 AM cache 容量
        for (int base = start_nnz; base < stop_nnz; base += step)
        {
            int end = (base + step < stop_nnz) ? (base + step) : stop_nnz;
            int length = end - base;
            if (length <= 0)
                continue;

            int aligned = length & ~7; // 向下取整到 8 的倍数（8个fp16=16B）
            int tail = length - aligned;

            // if (tid == 0 && task_id == 0)
            // {
            //     size_t bytes = (size_t)length * sizeof(__fp16);
            //     hthread_printf("base=%d len=%d bytes=%lu bytes/16=%lu\n",
            //                    base, length, (unsigned long)bytes, (unsigned long)(bytes >> 4));
            // }

            // 只对齐部分走 vector_load
            if (aligned > 0)
            {
                vector_load(&A_data[base], A_cache, aligned * sizeof(__fp16));
                vector_load(&x_data[base], x_cache, aligned * sizeof(__fp16));

                int vec_iters = aligned / 64; // 你后面按64个fp16一组计算
                for (int j = 0; j < vec_iters; j++)
                {
                    int offset = j * 16; // 你之前确认 offset 单位是 16B block
                    lvector __fp16 A_value = vec_ld(offset, A_cache);
                    lvector __fp16 x_value = vec_ld(offset, x_cache);
                    lvector __fp16 r = vec_muli(A_value, x_value);
                    vec_st(r, offset, x_cache);
                }

                if (vec_iters > 0)
                {
                    vector_store(x_cache, &x_data[base], vec_iters * 64 * sizeof(__fp16));
                }

                // 处理 aligned 里剩下但未满 64 的那一段（aligned%64），建议也用标量或再做一次小向量
                int rem = aligned - vec_iters * 64;
                for (int k = 0; k < rem; k++)
                {
                    int g = base + vec_iters * 64 + k;
                    x_data[g] = A_data[g] * x_data[g];
                }
            }

            // 尾巴 tail (<8) 直接标量
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

__global__ void SpMV_GSM_FP16_BLOCK4(int num_rows,
                                     int x_pack_size,
                                     int num_cores,
                                     int num_tasks,
                                     int barrier_id,
                                     int *A_i_block,
                                     int *A_j_block,
                                     __fp16 *A_data_block,
                                     __fp16 *x_pack,
                                     __fp16 *dy_data,
                                     int *row_bounds)
{
    int tid = get_thread_id();

    int AM_size = (96000 / 3) * 4;
    int max_blocks_per_load = (AM_size / 4) & ~1;

    __fp16 *gsm_x = (__fp16 *)gx;

    int pack_count = x_pack_size / 4;
    int base_pack = pack_count / num_cores;
    int remain_pack = pack_count % num_cores;

    int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
    int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

    int my_start = 4 * my_pack_start;
    int my_width = 4 * my_pack_width;

    if (my_start >= x_pack_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pack_size)
    {
        my_width = x_pack_size - my_start;
        my_width &= ~3;
    }

    if (x_pack_size <= 786432 * 4)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pack[my_start],
                                      1, my_width * sizeof(__fp16), 0,
                                      &gsm_x[my_start],
                                      1, my_width * sizeof(__fp16), 0,
                                      false, 0);
            dma_wait(ch);
        }

        core_barrier(barrier_id, num_cores);
    }
    else
    {
        if (tid == 0)
            hthread_printf("Warning: x_pack_size (%d) exceeds GSM capacity.\n", x_pack_size);
        return;
    }

    lvector __fp16 *A_cache = vector_malloc(AM_size * sizeof(__fp16));
    lvector __fp16 *x_cache = vector_malloc(AM_size * sizeof(__fp16));
    // lvector __fp16 *result_cache = vector_malloc(AM_size * sizeof(__fp16));

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
            int x_base = A_j_block[b] / 2;

            dy_data[4 * b] =
                A_data_block[4 * b] * gsm_x[x_base];

            dy_data[4 * b + 1] =
                A_data_block[4 * b + 1] * gsm_x[x_base + 1];

            dy_data[4 * b + 2] =
                A_data_block[4 * b + 2] * gsm_x[x_base + 2];

            dy_data[4 * b + 3] =
                A_data_block[4 * b + 3] * gsm_x[x_base + 3];

            start++;
        }

        if (start < end && (end & 1))
        {
            int b = end - 1;
            int x_base = A_j_block[b] / 2;

            dy_data[4 * b] =
                A_data_block[4 * b] * gsm_x[x_base];

            dy_data[4 * b + 1] =
                A_data_block[4 * b + 1] * gsm_x[x_base + 1];

            dy_data[4 * b + 2] =
                A_data_block[4 * b + 2] * gsm_x[x_base + 2];

            dy_data[4 * b + 3] =
                A_data_block[4 * b + 3] * gsm_x[x_base + 3];

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

            int half_count = 4 * num_block;

            vector_load(&A_data_block[4 * base],
                        A_cache,
                        half_count * sizeof(__fp16));

            unsigned int ch = dma_sg(gsm_x,
                                     &A_j_block[base],
                                     1,
                                     half_count * sizeof(__fp16),
                                     0,
                                     x_cache,
                                     1,
                                     half_count * sizeof(__fp16),
                                     0);
            dma_wait(ch);

            int vec_len = 64;
            int vec_iters = half_count / vec_len;

            for (int j = 0; j < vec_iters; j++)
            {
                int offset = j * 16;

                lvector __fp16 A_value = vec_ld(offset, A_cache);
                lvector __fp16 x_value = vec_ld(offset, x_cache);

                lvector __fp16 result_vector = vec_muli(A_value, x_value);

                // vec_st(result_vector, offset, result_cache);
                vec_st(result_vector, offset, x_cache);
            }

            if (vec_iters > 0)
            {
                // vector_store(result_cache, &dy_data[4 * base], vec_iters * vec_len * sizeof(__fp16));
                vector_store(x_cache, &dy_data[4 * base], vec_iters * vec_len * sizeof(__fp16));
            }

            int remain_half = half_count - vec_iters * vec_len;
            if (remain_half > 0)
            {
                int start_half = vec_iters * vec_len;

                for (int k = 0; k < remain_half; k++)
                {
                    int global_half = 4 * base + start_half + k;
                    int block_id = global_half / 4;
                    int lane = global_half & 3;
                    int x_base = A_j_block[block_id] / 2;

                    dy_data[global_half] = A_data_block[global_half] * gsm_x[x_base + lane];
                }
            }
        }
    }

    // vector_free(result_cache);
    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void SpMV_GSM_FP16_BLOCK4_F32PROD_PAIRSUM_VEC(int num_rows,
                                                         int x_pack_size,
                                                         int num_cores,
                                                         int num_tasks,
                                                         int barrier_id,
                                                         int *A_i_block,
                                                         int *A_j_block,
                                                         __fp16 *A_data_block,
                                                         __fp16 *x_pack,
                                                         float *block_pair_data,
                                                         int *row_bounds)
{
    int tid = get_thread_id();

    int AM_half_size = (96000 / 3) * 4;
    int max_blocks_per_load = (AM_half_size / 4) & ~15;

    __fp16 *gsm_x = (__fp16 *)gx;

    int pack_count = x_pack_size / 4;
    int base_pack = pack_count / num_cores;
    int remain_pack = pack_count % num_cores;

    int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
    int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

    int my_start = 4 * my_pack_start;
    int my_width = 4 * my_pack_width;

    if (my_start >= x_pack_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pack_size)
    {
        my_width = x_pack_size - my_start;
        my_width &= ~3;
    }

    if (x_pack_size <= 786432 * 4)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pack[my_start],
                                      1, my_width * sizeof(__fp16), 0,
                                      &gsm_x[my_start],
                                      1, my_width * sizeof(__fp16), 0,
                                      false, 0);
            dma_wait(ch);
        }

        core_barrier(barrier_id, num_cores);
    }
    else
    {
        if (tid == 0)
            hthread_printf("Warning: x_pack_size (%d) exceeds GSM capacity.\n", x_pack_size);
        return;
    }

    lvector __fp16 *A_cache = vector_malloc(AM_half_size * sizeof(__fp16));
    lvector __fp16 *x_cache = vector_malloc(AM_half_size * sizeof(__fp16));
    lvector float *pair_cache = vector_malloc(AM_half_size * sizeof(float));

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
            int x_base = A_j_block[b] / 2;

            float s0 =
                (float)A_data_block[4 * b] * (float)gsm_x[x_base] +
                (float)A_data_block[4 * b + 2] * (float)gsm_x[x_base + 2];

            float s1 =
                (float)A_data_block[4 * b + 1] * (float)gsm_x[x_base + 1] +
                (float)A_data_block[4 * b + 3] * (float)gsm_x[x_base + 3];

            block_pair_data[2 * b] = s0;
            block_pair_data[2 * b + 1] = s1;

            start++;
        }

        // if (start < end && (end & 1))
        // {
        //     int b = end - 1;
        //     int x_base = A_j_block[b] / 2;

        //     float s0 =
        //         (float)A_data_block[4 * b] * (float)gsm_x[x_base] +
        //         (float)A_data_block[4 * b + 2] * (float)gsm_x[x_base + 2];

        //     float s1 =
        //         (float)A_data_block[4 * b + 1] * (float)gsm_x[x_base + 1] +
        //         (float)A_data_block[4 * b + 3] * (float)gsm_x[x_base + 3];

        //     block_pair_data[2 * b] = s0;
        //     block_pair_data[2 * b + 1] = s1;

        //     end--;
        // }

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

            int vec_blocks = num_block & ~15;

            if (vec_blocks > 0)
            {
                int half_count = 4 * vec_blocks;
                int float_pair_count = 2 * vec_blocks;

                vector_load(&A_data_block[4 * base],
                            A_cache,
                            half_count * sizeof(__fp16));

                unsigned int ch = dma_sg(gsm_x,
                                         &A_j_block[base],
                                         1,
                                         half_count * sizeof(__fp16),
                                         0,
                                         x_cache,
                                         1,
                                         half_count * sizeof(__fp16),
                                         0);
                dma_wait(ch);

                int vec_iters = half_count / 64;

                for (int j = 0; j < vec_iters; j++)
                {
                    int h_offset = j * 16;
                    int f_offset = j * 16;

                    lvector __fp16 A_h = vec_ld(h_offset, A_cache);
                    lvector __fp16 x_h = vec_ld(h_offset, x_cache);

                    lvector float A_lo = vec_fstdl(A_h);
                    lvector float A_hi = vec_fstdh(A_h);

                    lvector float x_lo = vec_fstdl(x_h);
                    lvector float x_hi = vec_fstdh(x_h);

                    lvector float r_lo = vec_muli(A_lo, x_lo);
                    lvector float r_hi = vec_muli(A_hi, x_hi);

                    lvector float pair_sum = r_lo + r_hi;
                    vec_st(pair_sum, f_offset, pair_cache);
                }

                vector_store(pair_cache,
                             &block_pair_data[2 * base],
                             float_pair_count * sizeof(float));
            }

            for (int b = base + vec_blocks; b < stop; b++)
            {
                int x_base = A_j_block[b] / 2;

                float s0 =
                    (float)A_data_block[4 * b] * (float)gsm_x[x_base] +
                    (float)A_data_block[4 * b + 2] * (float)gsm_x[x_base + 2];

                float s1 =
                    (float)A_data_block[4 * b + 1] * (float)gsm_x[x_base + 1] +
                    (float)A_data_block[4 * b + 3] * (float)gsm_x[x_base + 3];

                block_pair_data[2 * b] = s0;
                block_pair_data[2 * b + 1] = s1;
            }
        }
    }

    vector_free(pair_cache);
    vector_free(x_cache);
    vector_free(A_cache);
}

__global__ void SpMV_GSM_FP16_BLOCK4_XF32_VEC_F32PROD_PAIRSUM(int num_rows,
                                                              int x_pack_size,
                                                              int num_cores,
                                                              int num_tasks,
                                                              int barrier_id,
                                                              int *A_i_block,
                                                              int *A_j_lo_f32,
                                                              int *A_j_hi_f32,
                                                              __fp16 *A_data_block,
                                                              float *x_pack_f32,
                                                              float *block_pair_data,
                                                              int *row_bounds)
{
    int tid = get_thread_id();

    // if (tid == 0)
    // {
    //     hthread_printf("SpMV_GSM_FP16_BLOCK4_XF32_VEC_F32PROD_PAIRSUM \n");
    // }

    const int XF32_VEC_BLOCKS = 4096;
    const int XF32_HALF_COUNT = 4 * XF32_VEC_BLOCKS;
    const int XF32_PAIR_COUNT = 2 * XF32_VEC_BLOCKS;

    float *gsm_x = (float *)gx;

    int pack_count = x_pack_size / 4;
    int base_pack = pack_count / num_cores;
    int remain_pack = pack_count % num_cores;

    int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
    int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

    int my_start = 4 * my_pack_start;
    int my_width = 4 * my_pack_width;

    if (my_start >= x_pack_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pack_size)
    {
        my_width = x_pack_size - my_start;
        my_width &= ~3;
    }

    if (x_pack_size <= 786432)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pack_f32[my_start],
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
        {
            hthread_printf("Warning: x_pack_size (%d) exceeds GSM FP32 capacity.\n",
                           x_pack_size);
        }
        return;
    }

    lvector __fp16 *A_cache =
        (lvector __fp16 *)vector_malloc(XF32_HALF_COUNT * sizeof(__fp16));

    lvector float *x_lo_cache =
        (lvector float *)vector_malloc(XF32_PAIR_COUNT * sizeof(float));

    lvector float *x_hi_cache =
        (lvector float *)vector_malloc(XF32_PAIR_COUNT * sizeof(float));

    lvector float *pair_cache =
        (lvector float *)vector_malloc(XF32_PAIR_COUNT * sizeof(float));

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];

        int start = A_i_block[task_start_row];
        int end = A_i_block[task_end_row];

        if (end <= start)
            continue;

        int base = start;

        for (; base < end; base += XF32_VEC_BLOCKS)
        {
            int stop = base + XF32_VEC_BLOCKS;
            if (stop > end)
                stop = end;

            int num_block = stop - base;
            if (num_block <= 0)
                continue;

            /*
             * Each vector iteration processes 16 blocks:
             *
             *   16 blocks * 4 half/block = 64 half
             *
             * This matches one lvector __fp16.
             */
            int vec_blocks = num_block & ~15;

            if (vec_blocks > 0)
            {
                int half_count = 4 * vec_blocks;
                int pair_count = 2 * vec_blocks;

                vector_load(&A_data_block[4 * base],
                            A_cache,
                            half_count * sizeof(__fp16));

                /*
                 * A_j_lo_f32[base] gathers x0,x1 for each block.
                 * A_j_hi_f32[base] gathers x2,x3 for each block.
                 *
                 * Each offset points to two contiguous FP32 values.
                 */
                unsigned int ch_lo = dma_sg(gsm_x,
                                            &A_j_lo_f32[base],
                                            1,
                                            pair_count * sizeof(float),
                                            0,
                                            x_lo_cache,
                                            1,
                                            pair_count * sizeof(float),
                                            0);

                unsigned int ch_hi = dma_sg(gsm_x,
                                            &A_j_hi_f32[base],
                                            1,
                                            pair_count * sizeof(float),
                                            0,
                                            x_hi_cache,
                                            1,
                                            pair_count * sizeof(float),
                                            0);

                dma_wait(ch_lo);
                dma_wait(ch_hi);

                /*
                 * One vector iteration consumes 64 half values from A_cache,
                 * corresponding to 16 blocks.
                 *
                 * For each 64-half vector:
                 *   vec_fstdl extracts A indices:
                 *      0,1, 4,5, ..., 60,61
                 *   vec_fstdh extracts A indices:
                 *      2,3, 6,7, ..., 62,63
                 *
                 * Therefore:
                 *   A_lo matches (a0,a1) of 16 blocks
                 *   A_hi matches (a2,a3) of 16 blocks
                 *
                 * x_lo_cache stores (x0,x1) of 16 blocks,
                 * x_hi_cache stores (x2,x3) of 16 blocks.
                 */
                int vec_iters = vec_blocks / 16;

                for (int j = 0; j < vec_iters; j++)
                {
                    int h_offset = j * 16;
                    int f_offset = j * 16;

                    lvector __fp16 A_h = vec_ld(h_offset, A_cache);

                    lvector float A_lo = vec_fstdl(A_h);
                    lvector float A_hi = vec_fstdh(A_h);

                    lvector float x_lo = vec_ld(f_offset, x_lo_cache);
                    lvector float x_hi = vec_ld(f_offset, x_hi_cache);

                    lvector float r_lo = vec_muli(A_lo, x_lo);
                    lvector float r_hi = vec_muli(A_hi, x_hi);

                    lvector float pair_sum = r_lo + r_hi;

                    vec_st(pair_sum, f_offset, pair_cache);
                }

                vector_store(pair_cache,
                             &block_pair_data[2 * base],
                             pair_count * sizeof(float));
            }

            /*
             * Scalar tail inside this tile.
             * Only the remaining num_block % 16 blocks go here.
             */
            for (int b = base + vec_blocks; b < stop; b++)
            {
                int x0_idx = A_j_lo_f32[b] / sizeof(float);
                int x2_idx = A_j_hi_f32[b] / sizeof(float);

                float a0 = (float)A_data_block[4 * b];
                float a1 = (float)A_data_block[4 * b + 1];
                float a2 = (float)A_data_block[4 * b + 2];
                float a3 = (float)A_data_block[4 * b + 3];

                float x0 = gsm_x[x0_idx];
                float x1 = gsm_x[x0_idx + 1];
                float x2 = gsm_x[x2_idx];
                float x3 = gsm_x[x2_idx + 1];

                block_pair_data[2 * b] =
                    a0 * x0 + a2 * x2;

                block_pair_data[2 * b + 1] =
                    a1 * x1 + a3 * x3;
            }
        }
    }

    vector_free(pair_cache);
    vector_free(x_hi_cache);
    vector_free(x_lo_cache);
    vector_free(A_cache);
}

// ==============================
// __global__ void
// SpMV_GSM_FP16_BLOCK2T_XF32_VEC_F32PROD_SUM(int num_rows,
//                                            int x_pack_size,
//                                            int num_cores,
//                                            int num_tasks,
//                                            int barrier_id,
//                                            int *A_i_block2,
//                                            int *A_j_x0_f32,
//                                            int *A_j_x1_f32,
//                                            __fp16 *A_data_block2_t,
//                                            float *x_pack_f32,
//                                            float *dy_block2_f32,
//                                            int *row_bounds)
// {
//     int tid = get_thread_id();

//     /*
//      * One vector iteration processes 32 BLOCK2 entries:
//      *
//      *   32 blocks * 2 half/block = 64 half
//      *
//      * vec_fstdl extracts lane-0 coefficients of 32 blocks.
//      * vec_fstdh extracts lane-1 coefficients of 32 blocks.
//      */
//     const int B2_VEC_BLOCKS = 4096;
//     const int B2_HALF_COUNT = 2 * B2_VEC_BLOCKS;
//     const int B2_FLOAT_COUNT = B2_VEC_BLOCKS;

//     float *gsm_x = (float *)gx;

//     int pack_count = x_pack_size / 4;
//     int base_pack = pack_count / num_cores;
//     int remain_pack = pack_count % num_cores;

//     int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
//     int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

//     int my_start = 4 * my_pack_start;
//     int my_width = 4 * my_pack_width;

//     if (my_start >= x_pack_size)
//     {
//         my_width = 0;
//     }
//     else if (my_start + my_width > x_pack_size)
//     {
//         my_width = x_pack_size - my_start;
//         my_width &= ~3;
//     }

//     if (x_pack_size <= 786432)
//     {
//         if (my_width > 0)
//         {
//             unsigned int ch = dma_p2p(&x_pack_f32[my_start],
//                                       1, my_width * sizeof(float), 0,
//                                       &gsm_x[my_start],
//                                       1, my_width * sizeof(float), 0,
//                                       false, 0);
//             dma_wait(ch);
//         }

//         core_barrier(barrier_id, num_cores);
//     }
//     else
//     {
//         if (tid == 0)
//         {
//             hthread_printf("Warning: x_pack_size (%d) exceeds GSM FP32 capacity in BLOCK2T.\n",
//                            x_pack_size);
//         }
//         return;
//     }

//     lvector __fp16 *A_cache =
//         (lvector __fp16 *)vector_malloc(B2_HALF_COUNT * sizeof(__fp16));

//     lvector float *x0_cache =
//         (lvector float *)vector_malloc(B2_FLOAT_COUNT * sizeof(float));

//     lvector float *x1_cache =
//         (lvector float *)vector_malloc(B2_FLOAT_COUNT * sizeof(float));

//     lvector float *sum_cache =
//         (lvector float *)vector_malloc(B2_FLOAT_COUNT * sizeof(float));

//     if (A_cache == 0 || x0_cache == 0 || x1_cache == 0 || sum_cache == 0)
//     {
//         if (tid == 0)
//             hthread_printf("Error: vector_malloc failed in BLOCK2T-XF32 kernel.\n");

//         if (sum_cache)
//             vector_free(sum_cache);
//         if (x1_cache)
//             vector_free(x1_cache);
//         if (x0_cache)
//             vector_free(x0_cache);
//         if (A_cache)
//             vector_free(A_cache);

//         return;
//     }

//     for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
//     {
//         int task_start_row = row_bounds[task_id];
//         int task_end_row = row_bounds[task_id + 1];

//         int start = A_i_block2[task_start_row];
//         int end = A_i_block2[task_end_row];

//         if (end <= start)
//             continue;

//         for (int base = start; base < end; base += B2_VEC_BLOCKS)
//         {
//             int stop = base + B2_VEC_BLOCKS;
//             if (stop > end)
//                 stop = end;

//             int cur = base;

//             /*
//              * The transposed BLOCK2 layout groups two consecutive blocks.
//              * If cur is odd, process one block scalar first to make cur even.
//              */
//             if (cur < stop && (cur & 1))
//             {
//                 int b = cur;
//                 int g = b >> 1;

//                 float a0 = (float)A_data_block2_t[4 * g + 1];
//                 float a1 = (float)A_data_block2_t[4 * g + 3];

//                 int x0_idx = A_j_x0_f32[b] / sizeof(float);
//                 int x1_idx = A_j_x1_f32[b] / sizeof(float);

//                 float x0 = gsm_x[x0_idx];
//                 float x1 = gsm_x[x1_idx];

//                 dy_block2_f32[b] = a0 * x0 + a1 * x1;

//                 cur++;
//             }

//             int num_block = stop - cur;
//             if (num_block <= 0)
//                 continue;

//             /*
//              * Each vector iteration handles 32 BLOCK2 entries.
//              */
//             // int vec_blocks = num_block & ~31;
//             int vec_blocks = 0;

//             if (vec_blocks > 0)
//             {
//                 int half_count = 2 * vec_blocks;
//                 int float_count = vec_blocks;

//                 /*
//                  * Since cur is even:
//                  *
//                  *   group index = cur / 2
//                  *   half offset = 4 * (cur / 2) = 2 * cur
//                  */
//                 vector_load(&A_data_block2_t[2 * cur],
//                             A_cache,
//                             half_count * sizeof(__fp16));

//                 /*
//                  * x0 and x1 are gathered separately.
//                  * Each offset points to one FP32 value.
//                  */
//                 unsigned int ch0 = dma_sg(gsm_x,
//                                           &A_j_x0_f32[cur],
//                                           1,
//                                           float_count * sizeof(float),
//                                           0,
//                                           x0_cache,
//                                           1,
//                                           float_count * sizeof(float),
//                                           0);

//                 unsigned int ch1 = dma_sg(gsm_x,
//                                           &A_j_x1_f32[cur],
//                                           1,
//                                           float_count * sizeof(float),
//                                           0,
//                                           x1_cache,
//                                           1,
//                                           float_count * sizeof(float),
//                                           0);

//                 dma_wait(ch0);
//                 dma_wait(ch1);

//                 int vec_iters = vec_blocks / 32;

//                 for (int j = 0; j < vec_iters; j++)
//                 {
//                     int h_offset = j * 16;
//                     int f_offset = j * 16;

//                     lvector __fp16 A_h = vec_ld(h_offset, A_cache);

//                     /*
//                      * Because of transposed BLOCK2 storage:
//                      *
//                      * vec_fstdl(A_h) -> a0 of 32 blocks
//                      * vec_fstdh(A_h) -> a1 of 32 blocks
//                      */
//                     lvector float A0 = vec_fstdl(A_h);
//                     lvector float A1 = vec_fstdh(A_h);

//                     lvector float x0 = vec_ld(f_offset, x0_cache);
//                     lvector float x1 = vec_ld(f_offset, x1_cache);

//                     lvector float r0 = vec_muli(A0, x0);
//                     lvector float r1 = vec_muli(A1, x1);

//                     lvector float sum = r0 + r1;

//                     vec_st(sum, f_offset, sum_cache);
//                 }

//                 vector_store(sum_cache,
//                              &dy_block2_f32[cur],
//                              float_count * sizeof(float));
//             }

//             /*
//              * Scalar tail.
//              */
//             for (int b = cur + vec_blocks; b < stop; b++)
//             {
//                 int g = b >> 1;

//                 float a0, a1;

//                 if ((b & 1) == 0)
//                 {
//                     a0 = (float)A_data_block2_t[4 * g + 0];
//                     a1 = (float)A_data_block2_t[4 * g + 2];
//                 }
//                 else
//                 {
//                     a0 = (float)A_data_block2_t[4 * g + 1];
//                     a1 = (float)A_data_block2_t[4 * g + 3];
//                 }

//                 int x0_idx = A_j_x0_f32[b] / sizeof(float);
//                 int x1_idx = A_j_x1_f32[b] / sizeof(float);

//                 float x0 = gsm_x[x0_idx];
//                 float x1 = gsm_x[x1_idx];

//                 dy_block2_f32[b] = a0 * x0 + a1 * x1;
//             }
//         }
//     }

//     vector_free(sum_cache);
//     vector_free(x1_cache);
//     vector_free(x0_cache);
//     vector_free(A_cache);
// }

__global__ void
SpMV_GSM_FP16_BLOCK2G_XF32_VEC_F32PROD_PAIR(int num_rows,
                                            int x_pack_size,
                                            int num_cores,
                                            int num_tasks,
                                            int barrier_id,
                                            int *A_i_block2,
                                            int *A_j_lo16_f32,
                                            int *A_j_hi16_f32,
                                            __fp16 *A_data_block2_g,
                                            float *x_pack_f32,
                                            float *dy_block2_pair,
                                            int *row_bounds)
{
    int tid = get_thread_id();

    const int B2_VEC_BLOCKS = 32768; /* must be multiple of 32 */
    const int B2_HALF_COUNT = 2 * B2_VEC_BLOCKS;
    const int B2_X_FLOAT_COUNT = B2_VEC_BLOCKS;

    float *gsm_x = (float *)gx;

    int pack_count = x_pack_size / 4;
    int base_pack = pack_count / num_cores;
    int remain_pack = pack_count % num_cores;

    int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
    int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

    int my_start = 4 * my_pack_start;
    int my_width = 4 * my_pack_width;

    if (my_start >= x_pack_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pack_size)
    {
        my_width = x_pack_size - my_start;
        my_width &= ~3;
    }

    if (x_pack_size <= 786432)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pack_f32[my_start],
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
            hthread_printf("Warning: x_pack_size exceeds GSM FP32 capacity in BLOCK2G.\n");
        return;
    }

    lvector __fp16 *A_cache =
        (lvector __fp16 *)vector_malloc(B2_HALF_COUNT * sizeof(__fp16));

    lvector float *x_lo_cache =
        (lvector float *)vector_malloc(B2_X_FLOAT_COUNT * sizeof(float));

    lvector float *x_hi_cache =
        (lvector float *)vector_malloc(B2_X_FLOAT_COUNT * sizeof(float));

    lvector float *prod_cache =
        (lvector float *)vector_malloc(32 * sizeof(float));

    // if (A_cache == 0 || x_lo_cache == 0 || x_hi_cache == 0 || prod_cache == 0)
    // {
    //     if (tid == 0)
    //         hthread_printf("Error: vector_malloc failed in BLOCK2G-XF32 kernel.\n");

    //     if (prod_cache)
    //         vector_free(prod_cache);
    //     if (x_hi_cache)
    //         vector_free(x_hi_cache);
    //     if (x_lo_cache)
    //         vector_free(x_lo_cache);
    //     if (A_cache)
    //         vector_free(A_cache);

    //     return;
    // }

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];

        int start = A_i_block2[task_start_row];
        int end = A_i_block2[task_end_row];

        if (end <= start)
            continue;

        for (int base = start; base < end; base += B2_VEC_BLOCKS)
        {
            int stop = base + B2_VEC_BLOCKS;
            if (stop > end)
                stop = end;

            int cur = base;

            /*
             * Scalar prefix until cur is aligned to a 32-block group.
             */
            for (; cur < stop && (cur & 31); cur++)
            {
                int g = cur >> 5;
                int r = cur & 31;
                int lane_base;
                float a0, a1;
                int off;
                int x0_idx;

                if (r < 16)
                    lane_base = 64 * g + 4 * r;
                else
                    lane_base = 64 * g + 4 * (r - 16) + 2;

                a0 = (float)A_data_block2_g[lane_base + 0];
                a1 = (float)A_data_block2_g[lane_base + 1];

                if (r < 16)
                    off = A_j_lo16_f32[16 * g + r];
                else
                    off = A_j_hi16_f32[16 * g + (r - 16)];

                x0_idx = off / sizeof(float);

                dy_block2_pair[2 * cur] =
                    a0 * gsm_x[x0_idx];

                dy_block2_pair[2 * cur + 1] =
                    a1 * gsm_x[x0_idx + 1];
            }

            int num_block = stop - cur;
            int vec_blocks = num_block & ~31;

            if (vec_blocks > 0)
            {
                int half_count = 2 * vec_blocks;
                int pair_count = vec_blocks / 2; /* number of 8-byte x-pair gathers */
                int x_float_count = vec_blocks;  /* pair_count * 2 */

                /*
                 * cur is 32-block aligned.
                 * A_data layout uses 64 half per 32 blocks.
                 * Therefore half offset is 2 * cur.
                 */
                vector_load(&A_data_block2_g[2 * cur],
                            A_cache,
                            half_count * sizeof(__fp16));

                int off16 = (cur >> 5) * 16;

                /*
                 * Each offset gathers two contiguous FP32 values, x0 and x1.
                 */
                unsigned int ch_lo = dma_sg(gsm_x,
                                            &A_j_lo16_f32[off16],
                                            1,
                                            x_float_count * sizeof(float),
                                            0,
                                            x_lo_cache,
                                            1,
                                            x_float_count * sizeof(float),
                                            0);

                unsigned int ch_hi = dma_sg(gsm_x,
                                            &A_j_hi16_f32[off16],
                                            1,
                                            x_float_count * sizeof(float),
                                            0,
                                            x_hi_cache,
                                            1,
                                            x_float_count * sizeof(float),
                                            0);

                dma_wait(ch_lo);
                dma_wait(ch_hi);

                int vec_iters = vec_blocks / 32;

                for (int j = 0; j < vec_iters; j++)
                {
                    int h_offset = j * 16;
                    int f_offset = j * 16;
                    int out_block = cur + 32 * j;

                    lvector __fp16 A_h = vec_ld(h_offset, A_cache);

                    /*
                     * A_lo: 16 blocks, natural [a0,a1] order.
                     * A_hi: next 16 blocks, natural [a0,a1] order.
                     */
                    lvector float A_lo = vec_fstdl(A_h);
                    lvector float A_hi = vec_fstdh(A_h);

                    lvector float x_lo = vec_ld(f_offset, x_lo_cache);
                    lvector float x_hi = vec_ld(f_offset, x_hi_cache);

                    lvector float r_lo = vec_muli(A_lo, x_lo);
                    lvector float r_hi = vec_muli(A_hi, x_hi);

                    vec_st(r_lo, 0, prod_cache);
                    vector_store(prod_cache,
                                 &dy_block2_pair[2 * out_block],
                                 32 * sizeof(float));

                    vec_st(r_hi, 0, prod_cache);
                    vector_store(prod_cache,
                                 &dy_block2_pair[2 * (out_block + 16)],
                                 32 * sizeof(float));
                }
            }

            /*
             * Scalar tail.
             */
            for (int b = cur + vec_blocks; b < stop; b++)
            {
                int g = b >> 5;
                int r = b & 31;
                int lane_base;
                int off;
                int x0_idx;
                float a0, a1;

                if (r < 16)
                {
                    lane_base = 64 * g + 4 * r;
                    off = A_j_lo16_f32[16 * g + r];
                }
                else
                {
                    lane_base = 64 * g + 4 * (r - 16) + 2;
                    off = A_j_hi16_f32[16 * g + (r - 16)];
                }

                a0 = (float)A_data_block2_g[lane_base + 0];
                a1 = (float)A_data_block2_g[lane_base + 1];

                x0_idx = off / sizeof(float);

                dy_block2_pair[2 * b] =
                    a0 * gsm_x[x0_idx];

                dy_block2_pair[2 * b + 1] =
                    a1 * gsm_x[x0_idx + 1];
            }
        }
    }

    vector_free(prod_cache);
    vector_free(x_hi_cache);
    vector_free(x_lo_cache);
    vector_free(A_cache);
}

__global__ void
SpMV_GSM_FP16_BLOCK2G_XF32_VEC_F32PROD_PAIR_1SG(int num_rows,
                                                int x_pack_size,
                                                int num_cores,
                                                int num_tasks,
                                                int barrier_id,
                                                int *A_i_block2,
                                                int *A_j_block_f32,
                                                __fp16 *A_data_block2_g,
                                                float *x_pack_f32,
                                                float *dy_block2_pair,
                                                int *row_bounds)
{
    int tid = get_thread_id();

    const int B2_VEC_BLOCKS = 32768; /* must be multiple of 32 */
    const int B2_HALF_COUNT = 2 * B2_VEC_BLOCKS;
    const int B2_X_FLOAT_COUNT = 2 * B2_VEC_BLOCKS;

    float *gsm_x = (float *)gx;

    int pack_count = x_pack_size / 4;
    int base_pack = pack_count / num_cores;
    int remain_pack = pack_count % num_cores;

    int my_pack_start = tid * base_pack + (tid < remain_pack ? tid : remain_pack);
    int my_pack_width = base_pack + (tid < remain_pack ? 1 : 0);

    int my_start = 4 * my_pack_start;
    int my_width = 4 * my_pack_width;

    if (my_start >= x_pack_size)
    {
        my_width = 0;
    }
    else if (my_start + my_width > x_pack_size)
    {
        my_width = x_pack_size - my_start;
        my_width &= ~3;
    }

    if (x_pack_size <= 786432)
    {
        if (my_width > 0)
        {
            unsigned int ch = dma_p2p(&x_pack_f32[my_start],
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
            hthread_printf("Warning: x_pack_size exceeds GSM FP32 capacity in BLOCK2G-1SG.\n");
        return;
    }

    lvector __fp16 *A_cache =
        (lvector __fp16 *)vector_malloc(B2_HALF_COUNT * sizeof(__fp16));

    /*
     * One cache for all x pairs:
     *   32 BLOCK2 -> 64 FP32 values.
     */
    lvector float *x_pair_cache =
        (lvector float *)vector_malloc(B2_X_FLOAT_COUNT * sizeof(float));

    lvector float *prod_cache =
        (lvector float *)vector_malloc(32 * sizeof(float));

    if (A_cache == 0 || x_pair_cache == 0 || prod_cache == 0)
    {
        if (tid == 0)
            hthread_printf("Error: vector_malloc failed in BLOCK2G-XF32-1SG kernel.\n");

        if (prod_cache)
            vector_free(prod_cache);
        if (x_pair_cache)
            vector_free(x_pair_cache);
        if (A_cache)
            vector_free(A_cache);

        return;
    }

    for (int task_id = tid; task_id < num_tasks; task_id += num_cores)
    {
        int task_start_row = row_bounds[task_id];
        int task_end_row = row_bounds[task_id + 1];

        int start = A_i_block2[task_start_row];
        int end = A_i_block2[task_end_row];

        if (end <= start)
            continue;

        for (int base = start; base < end; base += B2_VEC_BLOCKS)
        {
            int stop = base + B2_VEC_BLOCKS;
            if (stop > end)
                stop = end;

            int cur = base;

            /*
             * Scalar prefix until cur is aligned to a 32-block group.
             */
            for (; cur < stop && (cur & 31); cur++)
            {
                int g = cur >> 5;
                int r = cur & 31;
                int lane_base;
                int off;
                int x0_idx;
                float a0, a1;

                if (r < 16)
                    lane_base = 64 * g + 4 * r;
                else
                    lane_base = 64 * g + 4 * (r - 16) + 2;

                a0 = (float)A_data_block2_g[lane_base + 0];
                a1 = (float)A_data_block2_g[lane_base + 1];

                off = A_j_block_f32[cur];
                x0_idx = off / sizeof(float);

                dy_block2_pair[2 * cur] =
                    a0 * gsm_x[x0_idx];

                dy_block2_pair[2 * cur + 1] =
                    a1 * gsm_x[x0_idx + 1];
            }

            int num_block = stop - cur;
            int vec_blocks = num_block & ~31;

            if (vec_blocks > 0)
            {
                int half_count = 2 * vec_blocks;
                int x_float_count = 2 * vec_blocks;

                /*
                 * cur is 32-block aligned.
                 * A_data layout uses 64 half per 32 blocks.
                 * Therefore half offset is 2 * cur.
                 */
                vector_load(&A_data_block2_g[2 * cur],
                            A_cache,
                            half_count * sizeof(__fp16));

                /*
                 * Single dma_sg:
                 *   each A_j_block_f32 entry gathers one FP32 pair:
                 *      x[2*k], x[2*k+1]
                 *
                 * x_pair_cache layout:
                 *   block 0 : x0, x1
                 *   block 1 : x0, x1
                 *   ...
                 */
                unsigned int ch = dma_sg(gsm_x,
                                         &A_j_block_f32[cur],
                                         1,
                                         x_float_count * sizeof(float),
                                         0,
                                         x_pair_cache,
                                         1,
                                         x_float_count * sizeof(float),
                                         0);
                dma_wait(ch);

                int vec_iters = vec_blocks / 32;

                for (int j = 0; j < vec_iters; j++)
                {
                    /*
                     * A:
                     *   one lvector __fp16 = 64 half = 32 BLOCK2 coefficients
                     *
                     * x_pair_cache:
                     *   32 BLOCK2 pairs = 64 float
                     *   first 16 BLOCK2 pairs  -> first lvector float
                     *   second 16 BLOCK2 pairs -> second lvector float
                     */
                    int h_offset = j * 16;

                    /*
                     * vec_ld offset is in 8-byte units.
                     * 32 float = 128 bytes = 16 units.
                     * 64 float = 256 bytes = 32 units.
                     */
                    int x_offset_lo = j * 32;
                    int x_offset_hi = j * 32 + 16;

                    int out_block = cur + 32 * j;

                    lvector __fp16 A_h = vec_ld(h_offset, A_cache);

                    lvector float A_lo = vec_fstdl(A_h);
                    lvector float A_hi = vec_fstdh(A_h);

                    lvector float x_lo = vec_ld(x_offset_lo, x_pair_cache);
                    lvector float x_hi = vec_ld(x_offset_hi, x_pair_cache);

                    lvector float r_lo = vec_muli(A_lo, x_lo);
                    lvector float r_hi = vec_muli(A_hi, x_hi);

                    vec_st(r_lo, 0, prod_cache);
                    vector_store(prod_cache,
                                 &dy_block2_pair[2 * out_block],
                                 32 * sizeof(float));

                    vec_st(r_hi, 0, prod_cache);
                    vector_store(prod_cache,
                                 &dy_block2_pair[2 * (out_block + 16)],
                                 32 * sizeof(float));
                }
            }

            /*
             * Scalar tail.
             */
            for (int b = cur + vec_blocks; b < stop; b++)
            {
                int g = b >> 5;
                int r = b & 31;
                int lane_base;
                int off;
                int x0_idx;
                float a0, a1;

                if (r < 16)
                    lane_base = 64 * g + 4 * r;
                else
                    lane_base = 64 * g + 4 * (r - 16) + 2;

                a0 = (float)A_data_block2_g[lane_base + 0];
                a1 = (float)A_data_block2_g[lane_base + 1];

                off = A_j_block_f32[b];
                x0_idx = off / sizeof(float);

                dy_block2_pair[2 * b] =
                    a0 * gsm_x[x0_idx];

                dy_block2_pair[2 * b + 1] =
                    a1 * gsm_x[x0_idx + 1];
            }
        }
    }

    vector_free(prod_cache);
    vector_free(x_pair_cache);
    vector_free(A_cache);
}
