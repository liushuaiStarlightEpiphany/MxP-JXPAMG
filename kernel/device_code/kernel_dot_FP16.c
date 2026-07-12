#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"

extern __gsm__ char gx[6 * 1024 * 1024]; // 6MB = 6 * 1024 * 1024 字节

#define min(a, b) ((a) < (b) ? (a) : (b))

// SPU 版本 DOT ==============================
__global__ void kernel_dot_SPU_FP16(int n,
                                    __fp16 *x,
                                    __fp16 *y,
                                    int *split,
                                    __fp16 *result)
{
    int tid = get_thread_id();
    int start = split[tid];
    int stop = split[tid + 1];

    __fp16 local_sum = 0.0;

    for (int i = start; i < stop; i++)
    {
        local_sum += x[i] * y[i];
    }

    result[tid] = local_sum;
}

__global__ void kernel_dot_FP16(int chunk_size,
                                __fp16 *x,
                                __fp16 *y,
                                __fp16 *result)
{
    int tid = get_thread_id();
    if (tid == 0)
    {
        hthread_printf("\n kernel_dot_FP16 \n");
    }
    int start = tid * chunk_size;
    int stop = start + chunk_size;
    int vec_iters;

    // 向量缓存和累加器初始化
    int step = 2 * 2 * 32000;
    lvector __fp16 *x_cache = vector_malloc(step * sizeof(__fp16));
    lvector __fp16 *y_cache = vector_malloc(step * sizeof(__fp16));
    __fp16 local_sum = 0.0;

    lvector double temp;
    lvector __fp16 result_vector;
    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&x[i], x_cache, length * sizeof(__fp16));
        vector_load(&y[i], y_cache, length * sizeof(__fp16));

        temp = vec_svbcast(0.0);
        result_vector = vec_abs((lvector __fp16)temp);
        int vec_iters = length / 64;

        // 主循环 === 64 批次处理
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector __fp16 x_vector = vec_ld(offset, x_cache);
            lvector __fp16 y_vector = vec_ld(offset, y_cache);
            result_vector = vec_mula(x_vector, y_vector, result_vector);
        }
        local_sum += reduce_sum_fp16(result_vector);
    }

    result[tid] = local_sum;
    // hthread_printf("dot v3 tid = %d, local_sum = %lf\n", tid, local_sum);

    vector_free(x_cache);
    vector_free(y_cache);
}
