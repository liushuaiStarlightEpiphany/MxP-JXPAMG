#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

extern __gsm__ char gx[1024 * 128 * 6];

// __gsm__ double gsm_x[512 * 128 * 6]; // 3MB GSM for x / 393216
// __gsm__ double gsm_y[512 * 128 * 6]; // 3MB GSM for y / 393216

// SPU 版本 DOT ==============================
__global__ void kernel_dot_SPU(int n,
                               double *x,
                               double *y,
                               int *split,
                               double *result)
{
    int tid = get_thread_id();
    int start = split[tid];
    int stop = split[tid + 1];

    double local_sum = 0.0;

    for (int i = start; i < stop; i++)
    {
        local_sum += x[i] * y[i];
    }

    result[tid] = local_sum;
}

// VPU 版本 DOT (冯老板——补零为 16 的倍数) ================
__global__ void kernel_dot_v1(int new_n,
                              int thread_seg,
                              double *dx,
                              double *dy,
                              double *tid_result)
{
    int tid = get_thread_id();

    // if (tid == 0)
    // {
    //     printf("kernel_dot_v1 \n");
    // }

    int tid_start = 16 * thread_seg * tid;
    int tid_end = 16 * thread_seg * (tid + 1) < new_n ? 16 * thread_seg * (tid + 1) : new_n;
    int tid_length = tid_end - tid_start;

    int am_cache = 32000;
    lvector double *x_cache = (lvector double *)vector_malloc(sizeof(double) * am_cache);
    lvector double *y_cache = (lvector double *)vector_malloc(sizeof(double) * am_cache);
    lvector double *result_cache = (lvector double *)vector_malloc(sizeof(double) * am_cache);

    int seg = tid_length % am_cache == 0 ? (tid_length / am_cache) : (tid_length / am_cache + 1);
    lvector double x_vector, y_vector, result_vector;
    int offset;
    double sum = 0;

    for (int i = 0; i < seg; i++)
    {
        vector_load(dx + tid_start + i * am_cache, x_cache, sizeof(double) * am_cache);
        vector_load(dy + tid_start + i * am_cache, y_cache, sizeof(double) * am_cache);

        int seg_start = tid_start + i * am_cache;
        int seg_stop = seg_start + am_cache < tid_end ? seg_start + am_cache : tid_end;
        int seg_length = seg_stop - seg_start;

        int index = seg_length % 16 == 0 ? seg_length / 16 : (seg_length / 16) + 1;

        for (int j = 0; j < index; j++)
        {
            offset = j * 16;
            x_vector = vec_ld(offset, x_cache);
            y_vector = vec_ld(offset, y_cache);
            result_vector = vec_muli(x_vector, y_vector);
            sum += reduce_sum_fp64(result_vector);
        }
    }

    tid_result[tid] = sum;

    vector_free(x_cache);
    vector_free(y_cache);
    vector_free(result_cache);
}

// VPU 版本 DOT ====================================
// 数据流：DDR -> GSM -> AM, 大于小于 780000 均适用
// 不延长均匀划分版本，单缓冲区    ===================
// __global__ void kernel_dot_v5(int n,
//                               int coreNum,
//                               int barrier_id,
//                               double *x,
//                               double *y,
//                               double *result)
// {
//     int tid = get_thread_id();
//     const int BLOCK_SIZE = 390000;

//     // int total = 1024 * 128 * 6; // = 786432
//     // int half = total / 2;       // = 393216
//     // double *gsm_x = &gx[0];     // GSM 缓冲区前半段
//     // double *gsm_y = &gx[half];  // GSM 缓冲区后半段

//     // 总 GSM 空间大小（单位：字节）
//     const int total_bytes = 6 * 1024 * 1024;     // 6MB = 6×1024×1024
//     const int half_bytes = total_bytes / 2;      // 3MB = 前后各用一半
//     double *gsm_x = (double *)(gx);              // 起始地址
//     double *gsm_y = (double *)(gx + half_bytes); // 偏移 3MB（按字节）

//     // 初始化向量化缓冲区
//     int am_cache = 32000;
//     int step = min(n, am_cache);
//     lvector double *x_cache = vector_malloc(am_cache * sizeof(double));
//     lvector double *y_cache = vector_malloc(am_cache * sizeof(double));
//     lvector double x_vector, y_vector, result_vector;
//     double local_sum = 0.0;

//     for (int offset = 0; offset < n; offset += BLOCK_SIZE)
//     {
//         // 当前批次大小
//         int chunk_size = (n - offset) < BLOCK_SIZE ? (n - offset) : BLOCK_SIZE;
//         if (tid == 0)
//         {
//             unsigned int cx = dma_p2p(&x[offset], 1, chunk_size * sizeof(double), 0,
//                                       &gsm_x[0], 1, chunk_size * sizeof(double), 0, false, 0);
//             unsigned int cy = dma_p2p(&y[offset], 1, chunk_size * sizeof(double), 0,
//                                       &gsm_y[0], 1, chunk_size * sizeof(double), 0, false, 0);
//             dma_wait(cx);
//             dma_wait(cy);
//         }
//         core_barrier(barrier_id, coreNum);

//         // 线程任务划分
//         int p = chunk_size / coreNum;
//         int r = chunk_size % coreNum;
//         int row_start = offset + p * tid + min(tid, r);
//         int row_stop = row_start + p + (tid < r ? 1 : 0);
//         int trswidth = row_stop - row_start;
//         int local_gsm_base = row_start - offset; // 在 GSM 中的偏移

//         // 向量化计算内积
//         for (int ii = 0; ii < trswidth; ii += step)
//         {
//             int stride = min(ii + step, trswidth);
//             int length = stride - ii;

//             if (length != 0)
//             {
//                 vector_load(&gsm_x[local_gsm_base + ii], x_cache, length * sizeof(double));
//                 vector_load(&gsm_y[local_gsm_base + ii], y_cache, length * sizeof(double));

//                 result_vector = vec_svbcast(0.0);
//                 int vec_iters = length / 16;
//                 for (int j = 0; j < vec_iters; j++)
//                 {
//                     int j_offset = j * 16;
//                     x_vector = vec_ld(j_offset, x_cache);
//                     y_vector = vec_ld(j_offset, y_cache);
//                     result_vector = vec_mula(x_vector, y_vector, result_vector);
//                 }
//                 local_sum += reduce_sum_fp64(result_vector);

//                 int rem_start = vec_iters * 16;
//                 for (int j = 0; j < length - rem_start; j++)
//                 {
//                     local_sum += gsm_x[local_gsm_base + ii + rem_start + j] *
//                                  gsm_y[local_gsm_base + ii + rem_start + j];
//                 }
//             }
//         }
//         core_barrier(barrier_id, coreNum);
//     }

//     // 局部结果写出
//     result[tid] = local_sum;
//     // hthread_printf("dot v5 tid = %d, local_sum = %lf\n", tid, local_sum);

//     // core_barrier(barrier_id, coreNum);

//     // 释放向量缓冲区
//     vector_free(x_cache);
//     vector_free(y_cache);
// }

// VPU 版本 DOT ====================================
// 数据流：DDR -> GSM -> AM, 大于小于 780000 均适用
// 不延长均匀划分版本，双缓冲区    ===================
__global__ void kernel_dot_v6(int n,
                              int coreNum,
                              int barrier_id,
                              double *x,
                              double *y,
                              double *result)
{
    int tid = get_thread_id();
    const int BLOCK_SIZE = 195000;
    int total = 1024 * 128 * 6; // = 786432
    int half = total / 2;       // = 393216

    // GSM四缓冲区：x0/x1，y0/y1
    double *gsm_x0 = (double *)&gx[0];
    double *gsm_x1 = (double *)&gx[half / 2];
    double *gsm_y0 = (double *)&gx[half];
    double *gsm_y1 = (double *)&gx[half + half / 2];

    // 向量化缓冲区
    int am_cache = 32000;
    int step = min(n, am_cache);
    lvector double *x_cache = vector_malloc(am_cache * sizeof(double));
    lvector double *y_cache = vector_malloc(am_cache * sizeof(double));
    lvector double x_vector, y_vector, result_vector;
    double local_sum = 0.0;

    int buffer_flag = 0;
    int offset = 0;

    // 预加载第一个块（tid=0）
    if (tid == 0 && offset < n)
    {
        int chunk_size = min(BLOCK_SIZE, n - offset);
        unsigned int cx = dma_p2p(&x[offset], 1, chunk_size * sizeof(double), 0,
                                  gsm_x0, 1, chunk_size * sizeof(double), 0, false, 0);
        unsigned int cy = dma_p2p(&y[offset], 1, chunk_size * sizeof(double), 0,
                                  gsm_y0, 1, chunk_size * sizeof(double), 0, false, 0);
        dma_wait(cx);
        dma_wait(cy);
    }
    core_barrier(barrier_id, coreNum);

    for (offset = 0; offset < n; offset += BLOCK_SIZE)
    {
        int chunk_size = min(BLOCK_SIZE, n - offset);

        // 预加载下一块
        if (tid == 0 && offset + BLOCK_SIZE < n)
        {
            int next_chunk = min(BLOCK_SIZE, n - offset - BLOCK_SIZE);
            double *next_x = buffer_flag ? gsm_x0 : gsm_x1;
            double *next_y = buffer_flag ? gsm_y0 : gsm_y1;
            unsigned int cx = dma_p2p(&x[offset + BLOCK_SIZE], 1, next_chunk * sizeof(double), 0,
                                      next_x, 1, next_chunk * sizeof(double), 0, false, 0);
            unsigned int cy = dma_p2p(&y[offset + BLOCK_SIZE], 1, next_chunk * sizeof(double), 0,
                                      next_y, 1, next_chunk * sizeof(double), 0, false, 0);
            dma_wait(cx);
            dma_wait(cy);
        }

        core_barrier(barrier_id, coreNum);

        // 当前缓冲区选择
        double *cur_gsm_x = buffer_flag ? gsm_x1 : gsm_x0;
        double *cur_gsm_y = buffer_flag ? gsm_y1 : gsm_y0;

        // 线程任务划分
        int p = chunk_size / coreNum;
        int r = chunk_size % coreNum;
        int row_start = p * tid + (tid < r ? tid : r);
        int row_stop = row_start + p + (tid < r ? 1 : 0);
        int trswidth = row_stop - row_start;
        int local_base = row_start - offset; // 当前块内偏移

        // 向量化内积计算
        for (int ii = 0; ii < trswidth; ii += step)
        {
            int stride = min(ii + step, trswidth);
            int length = stride - ii;

            if (length > 0)
            {
                vector_load(&cur_gsm_x[local_base + ii], x_cache, length * sizeof(double));
                vector_load(&cur_gsm_y[local_base + ii], y_cache, length * sizeof(double));

                result_vector = vec_svbcast(0.0);
                int vec_iters = length / 16;
                for (int j = 0; j < vec_iters; j++)
                {
                    int j_offset = j * 16;
                    x_vector = vec_ld(j_offset, x_cache);
                    y_vector = vec_ld(j_offset, y_cache);
                    result_vector = vec_mula(x_vector, y_vector, result_vector);
                }
                local_sum += reduce_sum_fp64(result_vector);

                int rem_start = vec_iters * 16;
                for (int j = 0; j < length - rem_start; j++)
                {
                    local_sum += cur_gsm_x[local_base + ii + rem_start + j] *
                                 cur_gsm_y[local_base + ii + rem_start + j];
                }
            }
        }

        core_barrier(barrier_id, coreNum);
        buffer_flag = 1 - buffer_flag;
    }

    result[tid] = local_sum;
    vector_free(x_cache);
    vector_free(y_cache);
}

// xingxin =================================
__global__ void kernel_dot_FP64(double *x,
                                double *y,
                                double *result,
                                int *row_bounds)
{
    int tid = get_thread_id();
    int start = row_bounds[tid];
    int stop = row_bounds[tid + 1];
    int vec_iters;

    // 向量缓存和累加器初始化
    int step = 48000;
    lvector double *x_cache = vector_malloc(step * sizeof(double));
    lvector double *y_cache = vector_malloc(step * sizeof(double));
    lvector double result_vector;

    double local_sum = 0.0;
    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&x[i], x_cache, length * sizeof(double));
        vector_load(&y[i], y_cache, length * sizeof(double));

        result_vector = vec_svbcast(0.0);
        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            // hthread_printf("向量处理");
            int offset = j * 16;
            lvector double x_vector = vec_ld(offset, x_cache);
            lvector double y_vector = vec_ld(offset, y_cache);
            result_vector = vec_mula(x_vector, y_vector, result_vector);
        }
        local_sum += reduce_sum_fp64(result_vector);
    }

    result[tid] = local_sum;
    // hthread_printf("dot v3 tid = %d, local_sum = %lf\n", tid, local_sum);

    vector_free(x_cache);
    vector_free(y_cache);
}

// xingxin =================================
__global__ void kernel_dot_v3(int chunk_size,
                              double *x,
                              double *y,
                              double *result)
{
    int tid = get_thread_id();
    int start = tid * chunk_size;
    int stop = start + chunk_size;
    int vec_iters;

    // 向量缓存和累加器初始化
    int step = 32000;
    lvector double *x_cache = vector_malloc(step * sizeof(double));
    lvector double *y_cache = vector_malloc(step * sizeof(double));
    lvector double result_vector;

    double local_sum = 0.0;
    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&x[i], x_cache, length * sizeof(double));
        vector_load(&y[i], y_cache, length * sizeof(double));

        result_vector = vec_svbcast(0.0);
        int vec_iters = length / 16;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector double x_vector = vec_ld(offset, x_cache);
            lvector double y_vector = vec_ld(offset, y_cache);
            result_vector = vec_mula(x_vector, y_vector, result_vector);
        }
        local_sum += reduce_sum_fp64(result_vector);
    }

    result[tid] = local_sum;
    // hthread_printf("dot v3 tid = %d, local_sum = %lf\n", tid, local_sum);

    vector_free(x_cache);
    vector_free(y_cache);
}
