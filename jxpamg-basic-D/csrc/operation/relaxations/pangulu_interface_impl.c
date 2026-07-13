#ifndef CALCULATE_TYPE_R64
#define CALCULATE_TYPE_R64
#endif

#include "pangulu_interface.h"
#include "pangulu.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <unistd.h>  // for getenv
#include <pthread.h>  // for pthread_self
#include <sys/syscall.h>  // for syscall
#include <time.h>  // for clock_gettime
#include <sched.h>  // for CPU affinity

static inline size_t round_up_size(size_t value, size_t align)
{
    return (value + align - 1u) & ~(align - 1u);
}

static inline long long round_up_ll(long long value, long long multiple)
{
    return (value + multiple - 1LL) / multiple * multiple;
}

static void* malloc_page(size_t nbytes)
{
    void *p = NULL;
    /* Use 64-byte (512-bit) alignment to satisfy DMA requirements */
    /* 64-byte alignment is required for PanguLU on some architectures */
    size_t align = 64u;
    size_t size = round_up_size(nbytes, align);
    if (posix_memalign(&p, align, size) != 0) {
        return NULL;
    }
    /* Verify alignment is at least 8 bytes (64 bits) */
    if (((uintptr_t)p) % 8 != 0) {
        free(p);
        return NULL;
    }
    return p;
}

static inline int pad16(int n)
{
    int r = n % 16;
    return r ? (n + 16 - r) : n;
}

static inline int pick_nb(int n, int nb_hint)
{
    (void)n;
    if (nb_hint > 0) {
        return nb_hint;
    }
    return 32;
}

// 获取当前线程 ID（系统级）
static inline long get_thread_id(void)
{
    return (long)syscall(SYS_gettid);
}

// 获取当前 CPU ID
static inline int get_cpu_id(void)
{
    return sched_getcpu();
}

// 获取当前时间戳（纳秒）
static inline long long get_timestamp_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long long)ts.tv_sec * 1000000000LL + (long long)ts.tv_nsec;
}

// 打印线程和 CPU 信息
static void print_thread_info(const char *prefix)
{
    long tid = get_thread_id();
    int cpu = get_cpu_id();
    pthread_t ptid = pthread_self();
    long long ts = get_timestamp_ns();
    fprintf(stderr, "[JX PGLU] %s: tid=%ld ptid=%lu cpu=%d timestamp=%lld ns\n",
            prefix, tid, (unsigned long)ptid, cpu, ts);
    fflush(stderr);
}

static uint64_t hash64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

static uint64_t csr_fingerprint(int n, long long nnz,
                                const sparse_pointer_t *rowptr,
                                const sparse_index_t   *colidx)
{
    uint64_t h = 1469598103934665603ULL;
    h ^= (uint64_t)n;   h *= 1099511628211ULL;
    h ^= (uint64_t)nnz; h *= 1099511628211ULL;

    int step_r = (n + 1) > 64 ? (n + 1) / 64 : 1;
    for (int i = 0; i <= n; i += step_r) {
        h ^= (uint64_t)rowptr[i];
        h *= 1099511628211ULL;
    }

    long long step_c = nnz > 4096 ? nnz / 4096 : 1;
    for (long long k = 0; k < nnz; k += step_c) {
        h ^= (uint64_t)colidx[k];
        h *= 1099511628211ULL;
    }

    return hash64(h);
}

typedef struct {
    uint64_t          key;
    int               n;
    int               n_pad;
    int               nb;
    long long         nnz;
    long long         nnz_pad;

    sparse_pointer_t *rowptr_pad;
    sparse_index_t   *colidx_pad;
    sparse_value_t   *aval_pad;

    void             *handle;

    int               initialized;
} PGLU_Cache;

static PGLU_Cache g_cache = {0};

static void cache_clear(void)
{
    if (g_cache.handle) {
        pangulu_finalize(&g_cache.handle);
        g_cache.handle = NULL;
    }

    free(g_cache.rowptr_pad); g_cache.rowptr_pad = NULL;
    free(g_cache.colidx_pad); g_cache.colidx_pad = NULL;
    free(g_cache.aval_pad);   g_cache.aval_pad   = NULL;

    memset(&g_cache, 0, sizeof(g_cache));
}

int jx_pglu_ensure_factorized(
    int                       n,
    long long                 nnz,
    const sparse_pointer_t   *rowptr,
    const sparse_index_t     *colidx,
    const sparse_value_t     *aval,
    int                       nb_hint,
    int                       nthread_hint,
    int                      *n_pad_out)
{
    if (n <= 0 || nnz < 0 || !rowptr || !colidx || !aval) {
        return -1;
    }

    size_t row_bytes = (size_t)(n + 1) * sizeof(*rowptr);
    size_t col_bytes = (size_t)nnz     * sizeof(*colidx);

    fprintf(stderr, "[JX PGLU] ensure_factorized begin: n=%d nnz=%lld nb_hint=%d nthread_hint=%d\n",
            n, nnz, nb_hint, nthread_hint);

    sparse_pointer_t *row_copy = (sparse_pointer_t *)malloc_page(row_bytes);
    sparse_index_t   *col_copy = (sparse_index_t   *)malloc_page(col_bytes);
    if (!row_copy || !col_copy) {
        fprintf(stderr, "[JX PGLU] alloc row/col copy failed\n");
        free(row_copy);
        free(col_copy);
        return -2;
    }

    // printf("[DEBUG] rowptr (size %d):\n", n + 1);
    // for (int i = 0; i <= n; i++) {
    //     printf("%llu ", (unsigned long long)rowptr[i]);
    // }
    // printf("\n\n\n");
    
    // // 打印 colidx
    // printf("[DEBUG] colidx (size %d):\n", nnz);
    // for (int i = 0; i < nnz; i++) {
    //     printf("%u ", colidx[i]);
    // }
    // printf("\n\n\n");
    
    // // 打印 value
    // printf("[DEBUG] value (size %d):\n", nnz);
    // for (int i = 0; i < nnz; i++) {
    //     printf("%.6lf ", aval[i]);
    // }
    // printf("\n\n\n");
    
    memcpy(row_copy, rowptr, row_bytes);
    memcpy(col_copy, colidx, col_bytes);

    // 验证列索引是否在有效范围内 [0, n-1]
    int colidx_invalid = 0;
    int error_count = 0;
    for (long long k = 0; k < nnz; ++k) {
        if (col_copy[k] < 0 || col_copy[k] >= n) {
            if (error_count < 10) {
                fprintf(stderr, "[JX PGLU] ERROR: invalid colidx[%lld]=%u (n=%d)\n",
                        k, col_copy[k], n);
            }
            colidx_invalid = 1;
            error_count++;
        }
    }
    if (colidx_invalid) {
        fprintf(stderr, "[JX PGLU] colidx validation failed: %d invalid entries (total nnz=%lld), aborting\n",
                error_count, nnz);
        free(row_copy);
        free(col_copy);
        return -5;
    }

    uint64_t key_now = csr_fingerprint(n, nnz, row_copy, col_copy);
    int nb_use = pick_nb(n, nb_hint);
    int n_pad = pad16(n);

    fprintf(stderr, "[JX PGLU] fingerprint=0x%016llx n_pad=%d nb_use=%d\n",
            (unsigned long long)key_now, n_pad, nb_use);

    int need_build = 1;
    if (g_cache.initialized) {
        if (g_cache.key   == key_now &&
            g_cache.n     == n      &&
            g_cache.nnz   == nnz    &&
            g_cache.n_pad == n_pad  &&
            g_cache.nb    == nb_use) {
            need_build = 0;
            fprintf(stderr, "[JX PGLU] cache hit -> reuse existing factorization\n");
        } else {
            fprintf(stderr, "[JX PGLU] cache miss -> rebuild (old key=0x%016llx)\n",
                    (unsigned long long)g_cache.key);
            cache_clear();
        }
    }

    if (!need_build) {
        free(row_copy);
        free(col_copy);
        if (n_pad_out) {
            *n_pad_out = g_cache.n_pad;
        }
        fprintf(stderr, "[JX PGLU] ensure_factorized end (reused)\n");
        return 0;
    }

    int pad_rows = n_pad - n;
    long long nnz_pad = nnz + (long long)pad_rows;

    fprintf(stderr, "[JX PGLU] allocating padded buffers: n_pad=%d nnz_pad=%lld\n",
            n_pad, nnz_pad);

    sparse_pointer_t *rowpad = (sparse_pointer_t *)malloc_page((size_t)(n_pad + 1) * sizeof(*rowpad));
    sparse_index_t   *colpad = (sparse_index_t   *)malloc_page((size_t)nnz_pad     * sizeof(*colpad));
    sparse_value_t   *avalpad= (sparse_value_t   *)malloc_page((size_t)nnz_pad     * sizeof(*avalpad));
    if (!rowpad || !colpad || !avalpad) {
        fprintf(stderr, "[JX PGLU] padded allocation failed\n");
        free(rowpad); free(colpad); free(avalpad);
        free(row_copy); free(col_copy);
        return -3;
    }
    
    /* Immediately check alignment after allocation */
    uintptr_t rowpad_addr_check = (uintptr_t)rowpad;
    uintptr_t colpad_addr_check = (uintptr_t)colpad;
    uintptr_t avalpad_addr_check = (uintptr_t)avalpad;
    fprintf(stderr, "[JX PGLU] after malloc_page: rowpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)rowpad, rowpad_addr_check % 8, rowpad_addr_check % 64);
    fprintf(stderr, "[JX PGLU] after malloc_page: colpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)colpad, colpad_addr_check % 8, colpad_addr_check % 64);
    fprintf(stderr, "[JX PGLU] after malloc_page: avalpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)avalpad, avalpad_addr_check % 8, avalpad_addr_check % 64);
    fflush(stderr);

    memcpy(rowpad, row_copy, (size_t)(n + 1) * sizeof(*rowpad));
    memcpy(colpad, col_copy, (size_t)nnz * sizeof(*colpad));
    memcpy(avalpad, aval,    (size_t)nnz * sizeof(*avalpad));

    free(row_copy);
    free(col_copy);

    long long dest = nnz;
    for (int r = n; r < n_pad; ++r) {
        rowpad[r] = (sparse_pointer_t)dest;
        colpad[dest] = (sparse_index_t)r;
        avalpad[dest] = 1.0;
        ++dest;
    }
    rowpad[n_pad] = (sparse_pointer_t)dest;

    /* Print all debug info BEFORE calling pangulu_init */
    fprintf(stderr, "[JX PGLU] ========== BEFORE pangulu_init ==========\n");
    fprintf(stderr, "[JX PGLU] calling pangulu_init (nb=%d)\n", nb_use);
    fprintf(stderr, "[JX PGLU] pangulu_init params: n_pad=%d, nnz_pad=%lld\n", n_pad, nnz_pad);
    fprintf(stderr, "[JX PGLU] pangulu_init params: rowpad=%p, colpad=%p, avalpad=%p\n",
            (void*)rowpad, (void*)colpad, (void*)avalpad);
    fprintf(stderr, "[JX PGLU] pangulu_init params: sizeof(sparse_pointer_t)=%zu, sizeof(sparse_index_t)=%zu, sizeof(sparse_value_t)=%zu\n",
            sizeof(sparse_pointer_t), sizeof(sparse_index_t), sizeof(sparse_value_t));
    fprintf(stderr, "[JX PGLU] rowpad alignment: %zu mod 8, %zu mod 64\n", 
            (uintptr_t)rowpad % 8, (uintptr_t)rowpad % 64);
    fprintf(stderr, "[JX PGLU] colpad alignment: %zu mod 8, %zu mod 64\n", 
            (uintptr_t)colpad % 8, (uintptr_t)colpad % 64);
    fprintf(stderr, "[JX PGLU] avalpad alignment: %zu mod 8, %zu mod 64\n", 
            (uintptr_t)avalpad % 8, (uintptr_t)avalpad % 64);
    fflush(stderr);

    pangulu_init_options init_opts;
    memset(&init_opts, 0, sizeof(init_opts));
    init_opts.nb = nb_use;
    // 设置多线程数：如果指定了 nthread_hint 则使用，否则使用默认值
    // 默认值：如果环境变量 OMP_NUM_THREADS 存在则使用，否则使用 32
    if (nthread_hint > 0) {
        init_opts.nthread = nthread_hint;
    } else {
        const char *omp_threads = getenv("OMP_NUM_THREADS");
        if (omp_threads) {
            init_opts.nthread = atoi(omp_threads);
        } else {
            init_opts.nthread = 32;  // 默认值
        }
    }
    fprintf(stderr, "[JX PGLU] pangulu_init options: nb=%d, nthread=%d\n", 
            init_opts.nb, init_opts.nthread);
    
    // 打印调用 pangulu_init 前的线程信息
    print_thread_info("BEFORE pangulu_init call");
    fprintf(stderr, "[JX PGLU] Environment check: OMP_NUM_THREADS=%s, PGLU_NUM_THREADS=%s\n",
            getenv("OMP_NUM_THREADS") ? getenv("OMP_NUM_THREADS") : "(unset)",
            getenv("PGLU_NUM_THREADS") ? getenv("PGLU_NUM_THREADS") : "(unset)");
    fflush(stderr);

    void *handle = NULL;
    long long t_start = get_timestamp_ns();
    pangulu_init((sparse_index_t)n_pad, (sparse_pointer_t)nnz_pad,
                 rowpad, colpad, avalpad, &init_opts, &handle);
    long long t_end = get_timestamp_ns();
    long long t_elapsed = t_end - t_start;
    
    // 打印调用 pangulu_init 后的线程信息
    print_thread_info("AFTER pangulu_init call");
    fprintf(stderr, "[JX PGLU] pangulu_init elapsed time: %lld ns (%.6f ms)\n",
            t_elapsed, t_elapsed / 1000000.0);
    fflush(stderr);
    if (!handle) {
        fprintf(stderr, "[JX PGLU] pangulu_init returned NULL\n");
        free(rowpad); free(colpad); free(avalpad);
        return -4;
    }

    fprintf(stderr, "[JX PGLU] ========== AFTER pangulu_init ==========\n");
    fprintf(stderr, "[JX PGLU] pangulu_init ok handle=%p\n", handle);
    fflush(stderr);

    /* Verify alignment of arrays passed to PanguLU */
    uintptr_t rowpad_addr = (uintptr_t)rowpad;
    uintptr_t colpad_addr = (uintptr_t)colpad;
    uintptr_t avalpad_addr = (uintptr_t)avalpad;
    fprintf(stderr, "[JX PGLU] ========== ALIGNMENT CHECK ==========\n");
    fprintf(stderr, "[JX PGLU] alignment check: rowpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)rowpad, rowpad_addr % 8, rowpad_addr % 64);
    fprintf(stderr, "[JX PGLU] alignment check: colpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)colpad, colpad_addr % 8, colpad_addr % 64);
    fprintf(stderr, "[JX PGLU] alignment check: avalpad=%p (8-byte align=%zu, 64-byte align=%zu)\n",
            (void*)avalpad, avalpad_addr % 8, avalpad_addr % 64);
    fflush(stderr);
    if (rowpad_addr % 8 != 0 || colpad_addr % 8 != 0 || avalpad_addr % 8 != 0) {
        fprintf(stderr, "[JX PGLU] ERROR: arrays are not 8-byte aligned!\n");
        fflush(stderr);
        free(rowpad); free(colpad); free(avalpad);
        pangulu_finalize(&handle);
        return -6;
    }

    pangulu_gstrf_options gstrf_opts;
    memset(&gstrf_opts, 0, sizeof(gstrf_opts));
    fprintf(stderr, "[JX PGLU] ========== BEFORE pangulu_gstrf ==========\n");
    fprintf(stderr, "[JX PGLU] calling pangulu_gstrf (n_pad=%d, nnz_pad=%lld)\n", n_pad, nnz_pad);
    fprintf(stderr, "[JX PGLU] array sizes: rowpad=%zu bytes, colpad=%zu bytes, avalpad=%zu bytes\n",
            (size_t)(n_pad + 1) * sizeof(*rowpad),
            (size_t)nnz_pad * sizeof(*colpad),
            (size_t)nnz_pad * sizeof(*avalpad));
    fflush(stderr);
    
    // 打印调用 pangulu_gstrf 前的线程信息
    print_thread_info("BEFORE pangulu_gstrf call");
    long long t_gstrf_start = get_timestamp_ns();
    
    pangulu_gstrf(&gstrf_opts, &handle);
    
    long long t_gstrf_end = get_timestamp_ns();
    long long t_gstrf_elapsed = t_gstrf_end - t_gstrf_start;
    
    // 打印调用 pangulu_gstrf 后的线程信息
    print_thread_info("AFTER pangulu_gstrf call");
    fprintf(stderr, "[JX PGLU] pangulu_gstrf elapsed time: %lld ns (%.6f ms)\n",
            t_gstrf_elapsed, t_gstrf_elapsed / 1000000.0);
    fprintf(stderr, "[JX PGLU] pangulu_gstrf done\n");
    fflush(stderr);

    g_cache.key        = key_now;
    g_cache.n          = n;
    g_cache.nnz        = nnz;
    g_cache.n_pad      = n_pad;
    g_cache.nnz_pad    = nnz_pad;
    g_cache.nb         = nb_use;
    g_cache.rowptr_pad = rowpad;
    g_cache.colidx_pad = colpad;
    g_cache.aval_pad   = avalpad;
    g_cache.handle     = handle;
    g_cache.initialized = 1;

    if (n_pad_out) {
        *n_pad_out = n_pad;
    }

    fprintf(stderr, "[JX PGLU] ensure_factorized end (rebuild)\n");

    return 0;
}

int jx_pglu_solve(const double *rhs_pad, double *sol_pad)
{
    if (!g_cache.initialized || !g_cache.handle || !rhs_pad || !sol_pad) {
        fprintf(stderr, "[JX PGLU] solve: invalid state (initialized=%d handle=%p rhs=%p sol=%p)\n",
                g_cache.initialized, g_cache.handle, (const void*)rhs_pad, (void*)sol_pad);
        return -1;
    }

    fprintf(stderr, "[JX PGLU] solve begin: n_pad=%d\n", g_cache.n_pad);
    
    // 打印调用 solve 前的线程信息
    print_thread_info("BEFORE solve call");

    sparse_value_t *rhs = (sparse_value_t *)malloc_page((size_t)g_cache.n_pad * sizeof(*rhs));
    sparse_value_t *sol = (sparse_value_t *)malloc_page((size_t)g_cache.n_pad * sizeof(*sol));
    if (!rhs || !sol) {
        fprintf(stderr, "[JX PGLU] solve: staging alloc failed\n");
        free(rhs);
        free(sol);
        return -2;
    }

    memcpy(rhs, rhs_pad, (size_t)g_cache.n_pad * sizeof(*rhs));
    memcpy(sol, rhs,      (size_t)g_cache.n_pad * sizeof(*sol));

    pangulu_gstrs_options gstrs_opts;
    memset(&gstrs_opts, 0, sizeof(gstrs_opts));
    fprintf(stderr, "[JX PGLU] calling pangulu_gstrs\n");
    fflush(stderr);
    
    // 打印调用 pangulu_gstrs 前的线程信息
    print_thread_info("BEFORE pangulu_gstrs call");
    long long t_gstrs_start = get_timestamp_ns();
    
    pangulu_gstrs(sol, &gstrs_opts, &g_cache.handle);
    
    long long t_gstrs_end = get_timestamp_ns();
    long long t_gstrs_elapsed = t_gstrs_end - t_gstrs_start;
    
    // 打印调用 pangulu_gstrs 后的线程信息
    print_thread_info("AFTER pangulu_gstrs call");
    fprintf(stderr, "[JX PGLU] pangulu_gstrs elapsed time: %lld ns (%.6f ms)\n",
            t_gstrs_elapsed, t_gstrs_elapsed / 1000000.0);
    fprintf(stderr, "[JX PGLU] pangulu_gstrs done\n");
    fflush(stderr);

    memcpy(sol_pad, sol, (size_t)g_cache.n_pad * sizeof(*sol));

    free(rhs);
    free(sol);

    fprintf(stderr, "[JX PGLU] solve end\n");

    return 0;
}

void jx_pglu_release(void)
{
    fprintf(stderr, "[JX PGLU] release cache\n");
    cache_clear();
}

