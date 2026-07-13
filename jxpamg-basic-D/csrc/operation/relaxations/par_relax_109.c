//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*
 *  par_relax_109.c  -- RELAX_TYPE=9: Direct coarse solve via PanguLU (cached)
 */

#include "jx_hpcsr.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#ifdef _WIN32
#include <direct.h>
#include <io.h>
#include <sys/stat.h>
#include <windows.h>  // for GetCurrentThreadId, GetTickCount
#define mkdir(path, mode) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>  // for sysconf
#include <pthread.h>  // for pthread_self
#include <sys/syscall.h>  // for syscall
#include <time.h>  // for clock_gettime
#include <sched.h>  // for CPU affinity
#endif

#ifndef CALCULATE_TYPE_R64
#define CALCULATE_TYPE_R64
#endif

#include "pangulu_interface.h"

// 获取当前线程 ID（系统级）
static inline long get_thread_id(void)
{
#ifdef _WIN32
    return (long)GetCurrentThreadId();
#else
    return (long)syscall(SYS_gettid);
#endif
}

// 获取当前 CPU ID
static inline int get_cpu_id(void)
{
#ifdef _WIN32
    return -1;  // Windows 不支持
#else
    return sched_getcpu();
#endif
}

// 获取当前时间戳（纳秒）
static inline long long get_timestamp_ns(void)
{
#ifdef _WIN32
    // Windows 实现（简化版）
    return (long long)GetTickCount() * 1000000LL;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long long)ts.tv_sec * 1000000000LL + (long long)ts.tv_nsec;
#endif
}

// 打印线程和 CPU 信息
static void print_thread_info(const char *prefix)
{
#ifdef _WIN32
    long tid = get_thread_id();
    long long ts = get_timestamp_ns();
    fprintf(stderr, "[par_relax_109] %s: tid=%ld timestamp=%lld ns\n",
            prefix, tid, ts);
#else
    long tid = get_thread_id();
    int cpu = get_cpu_id();
    pthread_t ptid = pthread_self();
    long long ts = get_timestamp_ns();
    fprintf(stderr, "[par_relax_109] %s: tid=%ld ptid=%lu cpu=%d timestamp=%lld ns\n",
            prefix, tid, (unsigned long)ptid, cpu, ts);
#endif
    fflush(stderr);
}

static void* aligned_malloc64(size_t nbytes)
{
    void *p = NULL;
    size_t size = (nbytes + 4095u) & ~4095u;
    if (posix_memalign(&p, 4096u, size) != 0) {
        return NULL;
    }
    return p;
}

static int export_global_csr_from_par(
    jx_ParCSRMatrix *par_matrix,
    int *n_out, long long *nnz_out,
    sparse_pointer_t **rowptr_out,
    sparse_index_t   **colidx_out,
    double           **aval_out,
    jx_CSRMatrix    **A_CSR_out)
{
    jx_CSRMatrix *A = jx_ParCSRMatrixToCSRMatrixAll(par_matrix);
    if (!A) {
        fprintf(stderr, "[par_relax_109] to-CSR failed\n");
        return -1;
    }

    JX_Int  *Ai = jx_CSRMatrixI(A);
    JX_Int  *Aj = jx_CSRMatrixJ(A);
    JX_Real *Ax = jx_CSRMatrixData(A);

    int n   = (int)jx_CSRMatrixNumRows(A);
    long long nnz = (long long)Ai[n];

    sparse_pointer_t *rp = (sparse_pointer_t*)aligned_malloc64((size_t)(n + 1) * sizeof(*rp));
    sparse_index_t   *ci = (sparse_index_t*)  aligned_malloc64((size_t)nnz     * sizeof(*ci));
    double           *av = (double*)          aligned_malloc64((size_t)nnz     * sizeof(*av));
    if (!rp || !ci || !av) {
        free(rp); free(ci); free(av);
        jx_CSRMatrixDestroy(A);
        fprintf(stderr, "[par_relax_109] CSR alloc failed\n");
        return -2;
    }

    for (int i = 0; i <= n; ++i) {
        rp[i] = (sparse_pointer_t)Ai[i];
    }
    for (long long k = 0; k < nnz; ++k) {
        ci[k] = (sparse_index_t)Aj[k];
        av[k] = (double)Ax[k];
    }

    *n_out = n;
    *nnz_out = nnz;
    *rowptr_out = rp;
    *colidx_out = ci;
    *aval_out = av;
    *A_CSR_out = A;
    return 0;
}

static int export_rhs_pad(jx_ParVector *rhs, int n_pad, double **rhs_out)
{
    jx_Vector *v = rhs->local_vector;
    if (!v) {
        return -1;
    }

    double *buf = (double*)aligned_malloc64((size_t)n_pad * sizeof(double));
    if (!buf) {
        return -2;
    }

    int copy_n = (int)((v->size < n_pad) ? v->size : n_pad);
    for (int i = 0; i < copy_n; ++i) {
        buf[i] = (double)v->data[i];
    }
    for (int i = copy_n; i < n_pad; ++i) {
        buf[i] = 0.0;
    }

    *rhs_out = buf;
    return 0;
}

static int export_global_rhs_original(jx_ParVector *par_rhs, int n, double **rhs_out)
{
    if (!par_rhs) {
        fprintf(stderr, "[par_relax_109] export_global_rhs_original: par_rhs is NULL\n");
        return -1;
    }
    
    if (n <= 0) {
        fprintf(stderr, "[par_relax_109] export_global_rhs_original: invalid n=%d\n", n);
        return -1;
    }
    
    jx_Vector *v_global = jx_ParVectorToVectorAll(par_rhs);
    if (!v_global) {
        fprintf(stderr, "[par_relax_109] export_global_rhs_original: jx_ParVectorToVectorAll returned NULL\n");
        return -1;
    }

    JX_Real *v_data = jx_VectorData(v_global);
    if (!v_data) {
        fprintf(stderr, "[par_relax_109] export_global_rhs_original: jx_VectorData returned NULL\n");
        jx_SeqVectorDestroy(v_global);
        return -1;
    }
    
    int v_size = (int)jx_VectorSize(v_global);
    fprintf(stderr, "[par_relax_109] export_global_rhs_original: v_size=%d, n=%d\n", v_size, n);

    double *buf = (double*)aligned_malloc64((size_t)n * sizeof(double));
    if (!buf) {
        fprintf(stderr, "[par_relax_109] export_global_rhs_original: aligned_malloc64 failed for size %zu\n", (size_t)n * sizeof(double));
        jx_SeqVectorDestroy(v_global);
        return -2;
    }

    int copy_n = (v_size < n) ? v_size : n;
    for (int i = 0; i < copy_n; ++i) {
        buf[i] = (double)v_data[i];
    }
    for (int i = copy_n; i < n; ++i) {
        buf[i] = 0.0;
    }

    jx_SeqVectorDestroy(v_global);
    *rhs_out = buf;
    fprintf(stderr, "[par_relax_109] export_global_rhs_original: success, copied %d values\n", copy_n);
    return 0;
}

static void import_solution_unpad(const double *sol_pad, int n, jx_ParVector *par_app)
{
    jx_Vector *x = par_app->local_vector;
    for (int i = 0; i < n; ++i) {
        x->data[i] = (JX_Real)sol_pad[i];
    }
}

static int ensure_directory_exists(const char *dir_path)
{
#ifdef _WIN32
    struct _stat st = {0};
    if (_stat(dir_path, &st) == -1) {
        if (_mkdir(dir_path) != 0) {
            return -1;
        }
    }
#else
    struct stat st = {0};
    if (stat(dir_path, &st) == -1) {
        char tmp[512];
        char *p = NULL;
        size_t len;
        
        snprintf(tmp, sizeof(tmp), "%s", dir_path);
        len = strlen(tmp);
        if (tmp[len - 1] == '/') {
            tmp[len - 1] = 0;
        }
        for (p = tmp + 1; *p; p++) {
            if (*p == '/') {
                *p = 0;
                if (stat(tmp, &st) == -1) {
                    if (mkdir(tmp, 0755) != 0) {
                        return -1;
                    }
                }
                *p = '/';
            }
        }
        if (mkdir(tmp, 0755) != 0) {
            return -1;
        }
    }
#endif
    return 0;
}

static int save_csr_matrix_to_file(
    int n, long long nnz,
    sparse_pointer_t *rowptr,
    sparse_index_t   *colidx,
    double           *aval,
    const char *base_dir,
    int call_count)
{
    if (ensure_directory_exists(base_dir) != 0) {
        fprintf(stderr, "[par_relax_109] Failed to create directory: %s\n", base_dir);
        return -1;
    }

    char filename[512];
#ifdef _WIN32
    snprintf(filename, sizeof(filename), "%s\\matrix_%04d.mtx", base_dir, call_count);
#else
    snprintf(filename, sizeof(filename), "%s/matrix_%04d.mtx", base_dir, call_count);
#endif

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "[par_relax_109] Failed to open file for writing: %s\n", filename);
        return -2;
    }

    fprintf(fp, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%d %d %lld\n", n, n, nnz);

    for (int i = 0; i < n; ++i) {
        long long start = (long long)rowptr[i];
        long long end   = (long long)rowptr[i + 1];
        for (long long k = start; k < end; ++k) {
            int j = (int)colidx[k];
            double val = aval[k];
            fprintf(fp, "%d %d %.16e\n", i + 1, j + 1, val);
        }
    }

    fclose(fp);
    fprintf(stderr, "[par_relax_109] Saved CSR matrix to: %s (n=%d, nnz=%lld)\n", filename, n, nnz);
    return 0;
}


JX_Int jx_PAMGRelax109(jx_ParCSRMatrix *par_matrix,
                     jx_ParVector    *par_rhs,
                     JX_Int          *cf_marker,
                     JX_Int           relax_points,
                     JX_Real          relax_weight,
                     JX_Real          omega,
                     jx_ParVector    *par_app,
                     jx_ParVector    *par_tmp)
{
    (void)cf_marker;
    (void)relax_points;
    (void)relax_weight;
    (void)omega;
    (void)par_tmp;

    static int call_count = 0;
    call_count++;
    
    const char *mtxdata_dir_env = getenv("PGMTXDATA_DIR");
    const char *mtxdata_dir = (mtxdata_dir_env != NULL) ? mtxdata_dir_env : "./mtxdata";
    
    MPI_Comm comm = jx_ParCSRMatrixComm(par_matrix);
    JX_Int my_rank_int = 0;
    jx_MPI_Comm_rank(comm, &my_rank_int);
    int my_rank = (int)my_rank_int;

    int n = 0;
    long long nnz = 0;
    sparse_pointer_t *rowptr = NULL;
    sparse_index_t   *colidx = NULL;
    double           *aval   = NULL;
    jx_CSRMatrix     *A_CSR  = NULL;

    if (export_global_csr_from_par(par_matrix, &n, &nnz,
                                   &rowptr, &colidx, &aval, &A_CSR) != 0) {
        return -1;
    }

    fprintf(stderr, "[par_relax_109] ensure_factorized request: n=%d nnz=%lld\n", n, nnz);

    // 在调用 jx_pglu_ensure_factorized 之前保存矩阵和右端项
    // 这样可以确保即使 pangulu_init 卡住，数据也已经保存
    if (my_rank == 0) {
        save_csr_matrix_to_file(n, nnz, rowptr, colidx, aval, mtxdata_dir, call_count);
        
        fprintf(stderr, "[par_relax_109] Starting to export global RHS vector (n=%d)\n", n);
        
        double *rhs_global = NULL;
        if (export_global_rhs_original(par_rhs, n, &rhs_global) == 0) {
            fprintf(stderr, "[par_relax_109] Successfully exported global RHS (n=%d)\n", n);
            
            char filename_rhs[512];
#ifdef _WIN32
            snprintf(filename_rhs, sizeof(filename_rhs), "%s\\rhs_%04d.txt", mtxdata_dir, call_count);
#else
            snprintf(filename_rhs, sizeof(filename_rhs), "%s/rhs_%04d.txt", mtxdata_dir, call_count);
#endif
            
            FILE *fp_rhs = fopen(filename_rhs, "w");
            if (fp_rhs) {
                fprintf(fp_rhs, "%d\n", n);
                for (int i = 0; i < n; ++i) {
                    fprintf(fp_rhs, "%.16e\n", rhs_global[i]);
                }
                if (fclose(fp_rhs) == 0) {
                    fprintf(stderr, "[par_relax_109] Successfully saved global RHS vector to: %s (n=%d)\n", filename_rhs, n);
                } else {
                    fprintf(stderr, "[par_relax_109] WARNING: Error closing file: %s\n", filename_rhs);
                }
            } else {
                fprintf(stderr, "[par_relax_109] ERROR: Failed to open file for writing: %s (errno=%d)\n", filename_rhs, errno);
            }
            free(rhs_global);
        } else {
            fprintf(stderr, "[par_relax_109] ERROR: Failed to export global RHS\n");
        }
    }

    int n_pad = 0;
    const int nb_hint = 256;
    // 获取线程数：优先使用环境变量，否则使用默认值
    const char *omp_threads = getenv("OMP_NUM_THREADS");
    int nthread_hint = (omp_threads != NULL) ? atoi(omp_threads) : 32;
    // 如果环境变量未设置，尝试从系统获取 CPU 核心数
    if (nthread_hint <= 0) {
        long nproc = sysconf(_SC_NPROCESSORS_ONLN);
        if (nproc > 0) {
            nthread_hint = (int)nproc;
        } else {
            nthread_hint = 32;  // 默认值
        }
    }
    fprintf(stderr, "[par_relax_109] Using nthread=%d for PanguLU parallel computation\n", nthread_hint);
    
    // 打印调用 ensure_factorized 前的线程信息
    print_thread_info("BEFORE ensure_factorized call");
    fprintf(stderr, "[par_relax_109] Matrix info: n=%d nnz=%lld\n", n, nnz);
    fprintf(stderr, "[par_relax_109] Environment: OMP_NUM_THREADS=%s, PGLU_NUM_THREADS=%s\n",
            getenv("OMP_NUM_THREADS") ? getenv("OMP_NUM_THREADS") : "(unset)",
            getenv("PGLU_NUM_THREADS") ? getenv("PGLU_NUM_THREADS") : "(unset)");
    fflush(stderr);
    
    long long t_ensure_start = get_timestamp_ns();
    int prep = jx_pglu_ensure_factorized(n, nnz, rowptr, colidx, aval, nb_hint, nthread_hint, &n_pad);
    long long t_ensure_end = get_timestamp_ns();
    long long t_ensure_elapsed = t_ensure_end - t_ensure_start;
    
    // 打印调用 ensure_factorized 后的线程信息
    print_thread_info("AFTER ensure_factorized call");
    fprintf(stderr, "[par_relax_109] ensure_factorized result: code=%d n_pad=%d\n", prep, n_pad);
    fprintf(stderr, "[par_relax_109] ensure_factorized elapsed time: %lld ns (%.6f ms)\n",
            t_ensure_elapsed, t_ensure_elapsed / 1000000.0);
    fflush(stderr);

    free(rowptr);
    free(colidx);
    free(aval);
    jx_CSRMatrixDestroy(A_CSR);

    if (prep != 0) {
        fprintf(stderr, "[par_relax_109] ensure_factorized failed -> abort\n");
        return -2;
    }

    double *rhs_pad = NULL;
    if (export_rhs_pad(par_rhs, n_pad, &rhs_pad) != 0) {
        fprintf(stderr, "[par_relax_109] export_rhs_pad failed\n");
        return -3;
    }

    fprintf(stderr, "[par_relax_109] solve begin (n_pad=%d)\n", n_pad);
    
    // 打印调用 solve 前的线程信息
    print_thread_info("BEFORE solve call");

    double *sol_pad = (double*)aligned_malloc64((size_t)n_pad * sizeof(double));
    if (!sol_pad) {
        fprintf(stderr, "[par_relax_109] alloc sol_pad failed\n");
        free(rhs_pad);
        return -4;
    }

    long long t_solve_start = get_timestamp_ns();
    int solve_ret = jx_pglu_solve(rhs_pad, sol_pad);
    long long t_solve_end = get_timestamp_ns();
    long long t_solve_elapsed = t_solve_end - t_solve_start;
    
    // 打印调用 solve 后的线程信息
    print_thread_info("AFTER solve call");
    fprintf(stderr, "[par_relax_109] solve end code=%d\n", solve_ret);
    fprintf(stderr, "[par_relax_109] solve elapsed time: %lld ns (%.6f ms)\n",
            t_solve_elapsed, t_solve_elapsed / 1000000.0);
    fflush(stderr);
    if (solve_ret != 0) {
        fprintf(stderr, "[par_relax_109] jx_pglu_solve error -> abort\n");
        free(rhs_pad);
        free(sol_pad);
        return -5;
    }

    import_solution_unpad(sol_pad, n, par_app);

    free(rhs_pad);
    free(sol_pad);

    fprintf(stderr, "[par_relax_109] iteration done\n");

    return 0;
}

JX_Int jx_hpPAMGRelax109(jx_hpCSRMatrix *hp_matrix,
                       jx_ParVector *par_rhs,
                       JX_Int *cf_marker,
                       JX_Int relax_points,
                       JX_Real relax_weight,
                       JX_Real omega,
                       jx_ParVector *par_app,
                       jx_ParVector *Vtemp)
{
    return jx_PAMGRelax109(jx_hpCSRMatrixPar(hp_matrix), par_rhs,
                         cf_marker, relax_points, relax_weight,
                         omega, par_app, Vtemp);
}

void jx_PAMGRelax109_release_cache(void)
{
    jx_pglu_release();
}
