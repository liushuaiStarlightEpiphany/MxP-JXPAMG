#include <compiler/m3000.h>
#include "hthread_device.h"
#include "../vec_sum.h"

// extern __gsm__ char gx[1024 * 128 * 6];

extern __gsm__ char gx[6 * 1024 * 1024]; // 6MB = 6 * 1024 * 1024 字节

#define min(a, b) ((a) < (b) ? (a) : (b))

// SPU 版本 DOT ==============================
__global__ void kernel_dot_SPU_FP32(int n,
                                    float *x,
                                    float *y,
                                    int *split,
                                    float *result)
{
    int tid = get_thread_id();
    int start = split[tid];
    int stop = split[tid + 1];

    float local_sum = 0.0f;

    for (int i = start; i < stop; i++)
    {
        local_sum += x[i] * y[i];
    }

    result[tid] = local_sum;
}

// VPU + CPU 版本 DOT ===============================================
__global__ void kernel_dot_FP32(int chunk_size,
                                float *x,
                                float *y,
                                float *result)
{
    int tid = get_thread_id();
    if (tid == 0)
    {
        hthread_printf("\n kernel_dot_FP32 \n");
    }
    int start = tid * chunk_size;
    int stop = start + chunk_size;
    int vec_iters;

    // 向量缓存和累加器初始化
    int step = 2 * 32000;
    lvector float *x_cache = vector_malloc(step * sizeof(float));
    lvector float *y_cache = vector_malloc(step * sizeof(float));
    lvector float result_vector;

    // float local_sum = 0.0;
    double local_sum = 0.0;
    for (int i = start; i < stop; i += step)
    {
        int stride = (i + step < stop) ? (i + step) : stop;
        int length = stride - i;

        // 加载当前分块到 AM
        vector_load(&x[i], x_cache, length * sizeof(float));
        vector_load(&y[i], y_cache, length * sizeof(float));

        result_vector = vec_vsvbcast_sdf2(0.0);
        int vec_iters = length / 32;
        for (int j = 0; j < vec_iters; j++)
        {
            int offset = j * 16;
            lvector float x_vector = vec_ld(offset, x_cache);
            lvector float y_vector = vec_ld(offset, y_cache);
            result_vector = vec_mula(x_vector, y_vector, result_vector);
        }
        local_sum += reduce_sum_fp32(result_vector);
    }

    result[tid] = (float)local_sum;
    // hthread_printf("FP32 dot v3 tid = %d, local_sum = %lf\n", tid, local_sum);

    vector_free(x_cache);
    vector_free(y_cache);
}
