#include "beidoublas_spgemm.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <hthread_host.h>

#define CLUSTERID 1

/* block 大小来自 thmkl.h 中的定义：
   blc_len_m=4, blc_len_n=4, blc_size=16 */

/*---------- 前置声明 ----------*/

/* DSP SpGEMM 实现在 dsp_spgemm.h 中（头文件函数），其依赖的 thmkl_d_csrtombsr
   也在 util.h（被 dsp_spgemm.h 包含）中定义 */
#ifdef SPGEMM_ENABLE_DSP
#include "dsp_spgemm.h"
#else
/* 未启用 DSP 时的适配：在代码中引用会链接失败 */
#endif

/*---------- 算法选择 ----------*/

spgemm_algo_t spgemm_algo_from_env(void)
{
    const char *env = getenv(SPGEMM_ALGO_ENV_VAR);
    if (!env) return SPGEMM_ALGO_HASH;
    if (strcasecmp(env, SPGEMM_ALGO_MERGE_STR) == 0) return SPGEMM_ALGO_MERGE;
    if (strcasecmp(env, SPGEMM_ALGO_MIX_STR)   == 0) return SPGEMM_ALGO_MIX;
    if (strcasecmp(env, SPGEMM_ALGO_ORG_STR)   == 0) return SPGEMM_ALGO_ORG;
    if (strcasecmp(env, SPGEMM_ALGO_SPA_STR)   == 0) return SPGEMM_ALGO_SPA;
    if (strcasecmp(env, SPGEMM_ALGO_DSP_STR)   == 0) return SPGEMM_ALGO_DSP;
    return SPGEMM_ALGO_HASH;
}

const char* spgemm_algo_name(spgemm_algo_t algo)
{
    switch (algo) {
        case SPGEMM_ALGO_HASH:  return SPGEMM_ALGO_HASH_STR;
        case SPGEMM_ALGO_MERGE: return SPGEMM_ALGO_MERGE_STR;
        case SPGEMM_ALGO_MIX:   return SPGEMM_ALGO_MIX_STR;
        case SPGEMM_ALGO_ORG:   return SPGEMM_ALGO_ORG_STR;
        case SPGEMM_ALGO_SPA:   return SPGEMM_ALGO_SPA_STR;
        case SPGEMM_ALGO_DSP:   return SPGEMM_ALGO_DSP_STR;
        default:                return "unknown";
    }
}

int spgemm_algo_available(spgemm_algo_t algo)
{
    switch (algo) {
        case SPGEMM_ALGO_HASH:
        case SPGEMM_ALGO_MERGE:
        case SPGEMM_ALGO_MIX:
        case SPGEMM_ALGO_ORG:
        case SPGEMM_ALGO_SPA:
            return 1;   /* 库导出，理论上可用 */
        case SPGEMM_ALGO_DSP:
#ifdef SPGEMM_ENABLE_DSP
            return 1;
#else
            return 0;
#endif
        default:
            return 0;
    }
}

/*---------- 内部工具：获取 SpGEMM 函数指针 ----------*/

typedef sparse_status_t (*spgemm_func_t)(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);

static spgemm_func_t spgemm_func_by_algo(spgemm_algo_t algo)
{
    switch (algo) {
        case SPGEMM_ALGO_HASH:  return thmkl_sparse_spmm_hash;
        case SPGEMM_ALGO_MERGE: return thmkl_sparse_spmm_merge;
        case SPGEMM_ALGO_MIX:   return thmkl_sparse_spmm_mix;
        case SPGEMM_ALGO_ORG:   return thmkl_sparse_spmm_org;
        case SPGEMM_ALGO_SPA:   return thmkl_sparse_spmm_spa;
        default:                return thmkl_sparse_spmm_hash;
    }
}

/*---------- MT3000 设备初始化 ----------*/

static int beidoublas_dev_open = 0;

static int beidoublas_init(void)
{
    if (beidoublas_dev_open) return 0;
    if (hthread_dev_open(CLUSTERID) != 0) {
        fprintf(stderr, "[beidoublas] hthread_dev_open failed\n");
        return -1;
    }
    beidoublas_dev_open = 1;
    return 0;
}

#ifdef SPGEMM_ENABLE_DSP
/* DSP 内核加载标志 */
static int dsp_kernel_loaded = 0;

/* DSP SpGEMM 调用后关闭设备，释放线程资源 */
static void beidoublas_shutdown(void)
{
    if (!beidoublas_dev_open) return;
    hthread_dev_close(CLUSTERID);
    beidoublas_dev_open = 0;
    dsp_kernel_loaded = 0;
}
#endif

/*---------- 格式转换 ----------*/

/* jx_CSRMatrix → sparse_matrix_t（拷贝所有数组，BeidouBLAS 接管所有权） */
static int jx_to_beidou(jx_CSRMatrix *mat, sparse_matrix_t *out)
{
    THMKL2_INT rows = mat->num_rows;
    THMKL2_INT cols = mat->num_cols;
    THMKL2_INT nnz  = mat->num_nonzeros;

    THMKL2_INT *rows_start = (THMKL2_INT*)malloc((rows + 1) * sizeof(THMKL2_INT));
    THMKL2_INT *rows_end   = (THMKL2_INT*)malloc(rows * sizeof(THMKL2_INT));
    THMKL2_INT *col_idx    = (THMKL2_INT*)malloc(nnz * sizeof(THMKL2_INT));
    double      *values    = (double*)malloc(nnz * sizeof(double));

    if (!rows_start || !rows_end || !col_idx || !values) {
        free(rows_start); free(rows_end); free(col_idx); free(values);
        return -1;
    }

    for (THMKL2_INT i = 0; i <= rows; i++)
        rows_start[i] = mat->i[i];
    for (THMKL2_INT i = 0; i < rows; i++)
        rows_end[i] = mat->i[i + 1];
    for (THMKL2_INT k = 0; k < nnz; k++) {
        col_idx[k] = mat->j[k];
        values[k] = (double)mat->data[k];
    }

    sparse_status_t st = thmkl_sparse_d_create_csr(
        out, SPARSE_INDEX_BASE_ZERO, rows, cols,
        rows_start, rows_end, col_idx, values);

    if (st != 0) {
        free(rows_start); free(rows_end); free(col_idx); free(values);
    }
    return (st == 0) ? 0 : -1;
}

/* sparse_matrix_t → jx_CSRMatrix（拷贝数据到新分配的内存） */
static jx_CSRMatrix* beidou_to_jx(sparse_matrix_t mat)
{
    th_sparse_matrix_t sm = (th_sparse_matrix_t)mat;
    jx_CSRMatrix *out = (jx_CSRMatrix*)malloc(sizeof(jx_CSRMatrix));
    if (!out) return NULL;

    out->num_rows     = sm->m;
    out->num_cols     = sm->n;
    out->num_nonzeros = sm->nnz;
    out->rownnz       = NULL;
    out->num_rownnz   = 0;
    out->owns_data    = 1;

    out->i = (JX_Int*)malloc((sm->m + 1) * sizeof(JX_Int));
    out->j = (JX_Int*)malloc(sm->nnz * sizeof(JX_Int));
    out->data = (JX_Real*)malloc(sm->nnz * sizeof(JX_Real));
    if (!out->i || !out->j || !out->data) {
        free(out->i); free(out->j); free(out->data); free(out);
        return NULL;
    }

    for (int k = 0; k <= sm->m; k++)
        out->i[k] = (JX_Int)sm->ia[k];
    for (int k = 0; k < sm->nnz; k++) {
        out->j[k] = (JX_Int)sm->ja[k];
        out->data[k] = (JX_Real)sm->val[k];
    }
    return out;
}

#ifdef SPGEMM_ENABLE_DSP
/* 将 BSR 格式的输出矩阵（block size = blc_len_m × blc_len_n）展开为 CSR */
/* （仅在 DSP 启用时需要） */
static jx_CSRMatrix* bsr_to_jx(sparse_matrix_t mat,
                               int expected_rows, int expected_cols)
{
    th_sparse_matrix_t sm = (th_sparse_matrix_t)mat;
    int bsr_m = sm->bsr_m;
    int bsr_n = sm->bsr_n;
    int bsr_nnb = sm->bsr_nnb;
    int *bsria = sm->bsria;
    int *bsrja = sm->bsrja;
    double *bsrval = sm->bsrval;
    unsigned short *bit_map = sm->bit_map;

    int rows_pad = bsr_m * blc_len_m;
    int cols_pad = bsr_n * blc_len_n;
    int rows = (expected_rows > 0 && expected_rows < rows_pad) ? expected_rows : rows_pad;
    int cols = (expected_cols > 0 && expected_cols < cols_pad) ? expected_cols : cols_pad;

    /* 第一遍：统计实际非零元数（按 bit_map 筛选，排除 padding 区域） */
    int nnz = 0;
    for (int bi = 0; bi < bsr_m; bi++) {
        for (int bj = bsria[bi]; bj < bsria[bi + 1]; bj++) {
            int bcol = bsrja[bj];
            unsigned short bm = bit_map[bj];
            double *blk = &bsrval[bj * blc_size];
            for (int mi = 0; mi < blc_len_m; mi++) {
                int real_row = bi * blc_len_m + mi;
                if (real_row >= rows) continue;
                for (int mj = 0; mj < blc_len_n; mj++) {
                    int real_col = bcol * blc_len_n + mj;
                    if (real_col >= cols) continue;
                    if (bm & (1 << (mi * blc_len_n + mj)))
                        nnz++;
                }
            }
        }
    }

    jx_CSRMatrix *out = jx_CSRMatrixCreate(rows, cols, nnz);
    jx_CSRMatrixInitialize(out);
    JX_Int *ir = (JX_Int*)jx_CSRMatrixI(out);
    JX_Int *jc = (JX_Int*)jx_CSRMatrixJ(out);
    JX_Real *vd = (JX_Real*)jx_CSRMatrixData(out);

    /* 第二遍：展开 BSR 到 CSR（COO 格式中间存储，跳过 padding 区域） */
    JX_Int *coo_row = (JX_Int*)malloc(nnz * sizeof(JX_Int));
    JX_Int *coo_col = (JX_Int*)malloc(nnz * sizeof(JX_Int));
    JX_Real *coo_val = (JX_Real*)malloc(nnz * sizeof(JX_Real));
    if (!coo_row || !coo_col || !coo_val) {
        free(coo_row); free(coo_col); free(coo_val);
        jx_CSRMatrixDestroy(out);
        return NULL;
    }

    int pos = 0;
    for (int bi = 0; bi < bsr_m; bi++) {
        for (int bj = bsria[bi]; bj < bsria[bi + 1]; bj++) {
            int bcol = bsrja[bj];
            unsigned short bm = bit_map[bj];
            double *blk = &bsrval[bj * blc_size];
            for (int mi = 0; mi < blc_len_m; mi++) {
                int real_row = bi * blc_len_m + mi;
                if (real_row >= rows) continue;
                for (int mj = 0; mj < blc_len_n; mj++) {
                    int real_col = bcol * blc_len_n + mj;
                    if (real_col >= cols) continue;
                    if (bm & (1 << (mi * blc_len_n + mj))) {
                        coo_row[pos] = real_row;
                        coo_col[pos] = real_col;
                        coo_val[pos] = (JX_Real)blk[mi * blc_len_n + mj];
                        pos++;
                    }
                }
            }
        }
    }

    /* 按行排序转为 CSR */
    int *row_nnz = (int*)calloc(rows, sizeof(int));
    for (int k = 0; k < nnz; k++)
        row_nnz[coo_row[k]]++;
    int accum = 0;
    for (int i = 0; i < rows; i++) {
        ir[i] = accum;
        accum += row_nnz[i];
    }
    ir[rows] = accum;

    int *row_pos = (int*)calloc(rows, sizeof(int));
    for (int k = 0; k < nnz; k++) {
        int r = coo_row[k];
        int p = ir[r] + row_pos[r];
        jc[p] = coo_col[k];
        vd[p] = coo_val[k];
        row_pos[r]++;
    }

    free(coo_row); free(coo_col); free(coo_val);
    free(row_nnz); free(row_pos);
    return out;
}
#endif /* SPGEMM_ENABLE_DSP */

#ifdef SPGEMM_ENABLE_DSP
/* DSP 内核加载（一次） */
static int load_dsp_kernel(void)
{
    if (dsp_kernel_loaded) return 0;
    const char *kernel_path = getenv("SPGEMM_KERNEL");
    if (!kernel_path)
        kernel_path = "/vol8/home/xtu_pcy/l_s/BeidouBLAS_unified/dsp/spgemm_kernel.dat";
    if (hthread_dat_load(CLUSTERID, kernel_path) != 0) {
        fprintf(stderr, "[beidoublas] hthread_dat_load(%s) failed\n", kernel_path);
        return -1;
    }
    dsp_kernel_loaded = 1;
    return 0;
}
#endif

/* 释放 BeidouBLAS 矩阵（CPU CSR 版本） */
static void free_beidou(sparse_matrix_t mat)
{
    if (!mat) return;
    th_sparse_matrix_t sm = (th_sparse_matrix_t)mat;
    free(sm->ia);
    free(sm->ja);
    free(sm->val);
    hthread_free(sm);
}

#ifdef SPGEMM_ENABLE_DSP
/* 释放 DSP BSR 输出矩阵（仅 BSR 字段有效） */
static void free_beidou_dsp(sparse_matrix_t mat)
{
    if (!mat) return;
    th_sparse_matrix_t sm = (th_sparse_matrix_t)mat;
    if (sm->bsria) hthread_free(sm->bsria);
    if (sm->bsrja) hthread_free(sm->bsrja);
    if (sm->bsrval) hthread_free(sm->bsrval);
    if (sm->bit_map) hthread_free(sm->bit_map);
    hthread_free(sm);
}

/* 释放经 thmkl_d_csrtombsr 转换后的输入矩阵（CSR + BSR 字段均需释放） */
static void free_beidou_bsr_input(sparse_matrix_t mat)
{
    if (!mat) return;
    th_sparse_matrix_t sm = (th_sparse_matrix_t)mat;
    if (sm->bsria) hthread_free(sm->bsria);
    if (sm->bsrja) hthread_free(sm->bsrja);
    if (sm->bsrval) hthread_free(sm->bsrval);
    if (sm->bit_map) hthread_free(sm->bit_map);
    free(sm->ia);
    free(sm->ja);
    free(sm->val);
    hthread_free(sm);
}
#endif

/*---------- 公共 API 实现 ----------*/

int spgemm_adapter_multiply_ex(
    jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C,
    spgemm_algo_t algo)
{
    sparse_matrix_t beidou_A = NULL;
    sparse_matrix_t beidou_B = NULL;
    sparse_matrix_t beidou_C = NULL;
    int ret = -1;

    if (beidoublas_init() != 0) return -1;

    if (!spgemm_algo_available(algo)) {
        fprintf(stderr, "[beidoublas] 算法 %s 不可用，退回 hash\n", spgemm_algo_name(algo));
        algo = SPGEMM_ALGO_HASH;
    }

    spgemm_func_t spgemm = spgemm_func_by_algo(algo);
    const char *algo_name = spgemm_algo_name(algo);

    fprintf(stderr, "[beidoublas] 算法=%s A: %d x %d nnz=%d  B: %d x %d nnz=%d\n",
            algo_name, A->num_rows, A->num_cols, A->num_nonzeros,
            B->num_rows, B->num_cols, B->num_nonzeros);

    if (jx_to_beidou(A, &beidou_A) != 0) goto cleanup;
    if (jx_to_beidou(B, &beidou_B) != 0) goto cleanup;

    int st = -1;
    if (algo == SPGEMM_ALGO_DSP) {
#if defined(SPGEMM_ENABLE_DSP)
        /* DSP 路径：CSR → BSR → DSP SpGEMM → BSR → CSR */
        if (load_dsp_kernel() != 0) {
            fprintf(stderr, "[beidoublas] DSP 内核加载失败\n");
            goto cleanup;
        }
        /* 保存预期输出维度（BSR 转换后丢失精确 CSR 维度） */
        int expect_rows = A->num_rows;
        int expect_cols = B->num_cols;
        /* CSR → BSR 转换（就地） */
        thmkl_d_csrtombsr(&beidou_A);
        thmkl_d_csrtombsr(&beidou_B);

        double dsp_time = 0.0;
        st = thmkl_sparse_spmm_dsp(
            SPARSE_OPERATION_NON_TRANSPOSE, beidou_A, beidou_B, &beidou_C, &dsp_time);

        fprintf(stderr, "[beidoublas] DSP time=%.2f ms\n", dsp_time);

        /* thmkl_sparse_spmm_dsp 的 info 变量未初始化，返回值不可靠，
           因此只要 beidou_C 非空且内容合理即视为成功 */
        if (beidou_C) {
            th_sparse_matrix_t sm = (th_sparse_matrix_t)beidou_C;
            if (sm->bsr_m > 0 && sm->bsr_nnb >= 0) {
                *C = bsr_to_jx(beidou_C, expect_rows, expect_cols);
                if (*C) {
                    fprintf(stderr, "[beidoublas] C: %d x %d nnz=%d\n",
                            (*C)->num_rows, (*C)->num_cols, (*C)->num_nonzeros);
                    ret = 0;
                }
            }
        }
        fprintf(stderr, "[beidoublas] DSP %s\n", ret == 0 ? "成功" : "失败");
#else
        fprintf(stderr, "[beidoublas] DSP 未启用（需编译时加 -DSPGEMM_ENABLE_DSP）\n");
#endif
        goto cleanup;
    } else {
        st = spgemm(
            SPARSE_OPERATION_NON_TRANSPOSE, beidou_A, beidou_B, &beidou_C);

        /* thmkl_sparse_spmm_org 因库代码未初始化返回值，状态码不可靠，
           改以输出矩阵是否成功创建来判断。其他算法走正常 st==0 检查 */
        if (algo == SPGEMM_ALGO_ORG) {
            if (!beidou_C) {
                fprintf(stderr, "[beidoublas] %s 失败（输出为空）\n", algo_name);
                goto cleanup;
            }
        } else if (st != 0) {
            fprintf(stderr, "[beidoublas] %s 失败 返回=%d\n", algo_name, (int)st);
            goto cleanup;
        }

        {
            th_sparse_matrix_t sm = (th_sparse_matrix_t)beidou_C;
            fprintf(stderr, "[beidoublas] C: %d x %d nnz=%d\n", sm->m, sm->n, sm->nnz);
        }

        *C = beidou_to_jx(beidou_C);
        if (!*C) goto cleanup;

        ret = 0;
    }
cleanup:
#if defined(SPGEMM_ENABLE_DSP)
    if (algo == SPGEMM_ALGO_DSP) {
        free_beidou_bsr_input(beidou_A);
        free_beidou_bsr_input(beidou_B);
        free_beidou_dsp(beidou_C);
        beidoublas_shutdown();
    } else {
        free_beidou(beidou_A);
        free_beidou(beidou_B);
        free_beidou(beidou_C);
    }
#else
    free_beidou(beidou_A);
    free_beidou(beidou_B);
    free_beidou(beidou_C);
#endif
    return ret;
}

/* 使用 BeidouBLAS 计算数值，填入预分配好的符号结构 C（i, j 已由 JXPAMG 计算好）。
 * A, B 为输入；C 需有正确的 i, j 数组，此函数仅填充 C->data。
 * 对于 DSP 算法：先 DSP 尝试，结构不匹配则退回 JXPAMG 数值。
 * 对于 CPU 算法：直接调 BeidouBLAS 库函数，结构不匹配则退回 JXPAMG。
 * 返回 0 成功，非 0 失败。 */
int spgemm_adapter_numeric_multiply(
    jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix *C)
{
    spgemm_algo_t algo = spgemm_algo_from_env();
    const char *algo_name = spgemm_algo_name(algo);

    fprintf(stderr, "[bdblas-num] algo=%s A:%dx%d nnz=%d B:%dx%d nnz=%d C:%dx%d nnz=%d\n",
            algo_name, A->num_rows, A->num_cols, A->num_nonzeros,
            B->num_rows, B->num_cols, B->num_nonzeros,
            C->num_rows, C->num_cols, C->num_nonzeros);

#if defined(SPGEMM_ENABLE_DSP)
    if (algo == SPGEMM_ALGO_DSP) {
        int dsp_ok = 0;
        sparse_matrix_t beidou_A = NULL;
        sparse_matrix_t beidou_B = NULL;
        sparse_matrix_t beidou_C = NULL;

        if (beidoublas_init() == 0 && load_dsp_kernel() == 0) {
            int expect_rows = A->num_rows;
            int expect_cols = B->num_cols;

            if (jx_to_beidou(A, &beidou_A) == 0 && jx_to_beidou(B, &beidou_B) == 0) {
                thmkl_d_csrtombsr(&beidou_A);
                thmkl_d_csrtombsr(&beidou_B);

                double dsp_time = 0.0;
                thmkl_sparse_spmm_dsp(
                    SPARSE_OPERATION_NON_TRANSPOSE, beidou_A, beidou_B, &beidou_C, &dsp_time);
                fprintf(stderr, "[bdblas-num] DSP time=%.2f ms\n", dsp_time);

                if (beidou_C) {
                    th_sparse_matrix_t sm = (th_sparse_matrix_t)beidou_C;
                    if (sm->bsr_m > 0 && sm->bsr_nnb >= 0) {
                        jx_CSRMatrix *Cdsp = bsr_to_jx(beidou_C, expect_rows, expect_cols);
                        if (Cdsp) {
                            int match = (Cdsp->num_rows == C->num_rows &&
                                         Cdsp->num_cols == C->num_cols &&
                                         Cdsp->num_nonzeros == C->num_nonzeros);
                            if (match) {
                                /* 按列索引匹配数值：对每行，在 DSP 中找对应列的值填入 C */
                                for (int i = 0; i < C->num_rows && match; i++) {
                                    int len = C->i[i+1] - C->i[i];
                                    if (len <= 0) continue;
                                    for (int k = 0; k < len; k++) {
                                        JX_Int col = C->j[C->i[i] + k];
                                        int found = 0;
                                        for (int k2 = 0; k2 < len; k2++) {
                                            if (Cdsp->j[Cdsp->i[i] + k2] == col) {
                                                C->data[C->i[i] + k] = Cdsp->data[Cdsp->i[i] + k2];
                                                found = 1;
                                                break;
                                            }
                                        }
                                        if (!found) { match = 0; break; }
                                    }
                                }
                            }
                            if (match) {
                                dsp_ok = 1;
                                fprintf(stderr, "[bdblas-num] DSP 数值填充成功 ✓\n");
                            } else {
                                fprintf(stderr, "[bdblas-num] DSP 列集合不匹配（某列在 JXPAMG 中但 DSP 中没有）\n");
                            }
                            jx_CSRMatrixDestroy(Cdsp);
                        }
                    }
                }
                free_beidou_bsr_input(beidou_A);
                free_beidou_bsr_input(beidou_B);
                free_beidou_dsp(beidou_C);
            } else {
                free_beidou(beidou_A);
                free_beidou(beidou_B);
            }
            beidoublas_shutdown();
        }
        if (dsp_ok) return 0;
        goto jxpamg_fallback;
    }
#endif /* SPGEMM_ENABLE_DSP */

    /* --- CPU 算法路径 (hash/merge/mix/org/spa) --- */
    {
        sparse_matrix_t beidou_A = NULL;
        sparse_matrix_t beidou_B = NULL;
        sparse_matrix_t beidou_C = NULL;

        if (beidoublas_init() != 0) goto jxpamg_fallback;

        if (jx_to_beidou(A, &beidou_A) != 0) { goto jxpamg_fallback; }
        if (jx_to_beidou(B, &beidou_B) != 0) { free_beidou(beidou_A); goto jxpamg_fallback; }

        spgemm_func_t spgemm = spgemm_func_by_algo(algo);
        int st = spgemm(SPARSE_OPERATION_NON_TRANSPOSE, beidou_A, beidou_B, &beidou_C);

        if (st != 0 || !beidou_C) {
            fprintf(stderr, "[bdblas-num] CPU %s 失败\n", algo_name);
            free_beidou(beidou_A);
            free_beidou(beidou_B);
            free_beidou(beidou_C);
            goto jxpamg_fallback;
        }

        th_sparse_matrix_t sm = (th_sparse_matrix_t)beidou_C;
        int match = (sm->m == C->num_rows && sm->n == C->num_cols &&
                     sm->nnz == C->num_nonzeros);
        if (match) {
            for (int i = 0; i < C->num_rows && match; i++) {
                int len = C->i[i+1] - C->i[i];
                for (int k = 0; k < len; k++) {
                    JX_Int col = C->j[C->i[i] + k];
                    int found = 0;
                    for (int k2 = 0; k2 < len; k2++) {
                        if ((JX_Int)sm->ja[sm->ia[i] + k2] == col) {
                            C->data[C->i[i] + k] = (JX_Real)sm->val[sm->ia[i] + k2];
                            found = 1;
                            break;
                        }
                    }
                    if (!found) { match = 0; break; }
                }
            }
        }
        if (match) {
            free_beidou(beidou_A);
            free_beidou(beidou_B);
            free_beidou(beidou_C);
            fprintf(stderr, "[bdblas-num] CPU %s 数值填充成功 ✓\n", algo_name);
            return 0;
        }

        fprintf(stderr, "[bdblas-num] CPU %s 结构不匹配，退回 JXPAMG\n", algo_name);
        free_beidou(beidou_A);
        free_beidou(beidou_B);
        free_beidou(beidou_C);
    }

jxpamg_fallback:
    {
        jx_CSRMatrix *Cfull = jx_CSRMatrixMultiply(A, B);
        if (!Cfull) {
            fprintf(stderr, "[bdblas-num] JXPAMG 全量乘法失败\n");
            return -1;
        }
        if (Cfull->num_nonzeros != C->num_nonzeros) {
            fprintf(stderr, "[bdblas-num] JXPAMG nnz 不匹配 %d vs %d\n",
                    Cfull->num_nonzeros, C->num_nonzeros);
            jx_CSRMatrixDestroy(Cfull);
            return -1;
        }
        for (int k = 0; k < C->num_nonzeros; k++)
            C->data[k] = Cfull->data[k];
        jx_CSRMatrixDestroy(Cfull);
        fprintf(stderr, "[bdblas-num] JXPAMG 数值填充成功 ✓\n");
        return 0;
    }
}
