#ifndef BEIDOUBLAS_SPGEMM_H
#define BEIDOUBLAS_SPGEMM_H

#include "jx_mv.h"
#include "thmkl.h"

/* 可用的 BeidouBLAS SpGEMM 算法枚举 */
typedef enum {
    SPGEMM_ALGO_HASH  = 0,   /* 哈希表法（默认，已验证） */
    SPGEMM_ALGO_MERGE = 1,   /* 归并法 */
    SPGEMM_ALGO_MIX   = 2,   /* 混合法 */
    SPGEMM_ALGO_ORG   = 3,   /* 原始法（库函数返回值未初始化，适配器已做绕过处理） */
    SPGEMM_ALGO_SPA   = 4,   /* SPA（稀疏累加器）法 */
    SPGEMM_ALGO_DSP   = 5,   /* DSP 硬件加速（实验性，需编译时加 -DSPGEMM_ENABLE_DSP；大矩阵可能崩溃） */
} spgemm_algo_t;

/* 算法名称字符串，用于从环境变量读取 */
#define SPGEMM_ALGO_HASH_STR  "hash"
#define SPGEMM_ALGO_MERGE_STR "merge"
#define SPGEMM_ALGO_MIX_STR   "mix"
#define SPGEMM_ALGO_ORG_STR   "org"
#define SPGEMM_ALGO_SPA_STR   "spa"
#define SPGEMM_ALGO_DSP_STR   "dsp"
#define SPGEMM_ALGO_ENV_VAR   "SPGEMM_ALGO"

/*========== 库函数前置声明 ==========*/
/* 声明均在 thmkl.h 中有同名 API，但库仅导出 _bsr 和内部变体 */

/* 库导出的 CSR SpGEMM（4 参数，CSR 直入） */
sparse_status_t thmkl_sparse_spmm_hash(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);
sparse_status_t thmkl_sparse_spmm_merge(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);
sparse_status_t thmkl_sparse_spmm_mix(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);
sparse_status_t thmkl_sparse_spmm_org(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);
sparse_status_t thmkl_sparse_spmm_spa(
    sparse_operation_t, sparse_matrix_t, sparse_matrix_t, sparse_matrix_t*);

/* CSR → BSR 就地转换（库导出，声明在 thmkl.h 但未暴露） */
void thmkl_d_csrtobsr(sparse_matrix_t *A);

/*========== 公共适配器 API ==========*/

/* 从环境变量 SPGEMM_ALGO 读取算法，失败返回默认值 */
spgemm_algo_t spgemm_algo_from_env(void);

/* 将算法枚举转为名称字符串 */
const char* spgemm_algo_name(spgemm_algo_t algo);

/* 原始的 jx_CSRMatrix → jx_CSRMatrix SpGEMM（使用 env 选择的算法）；
 * 返回 0 成功，非 0 失败。
 * 输入 A、B 不会被修改；
 * 输出 *C 是新分配的 jx_CSRMatrix，调用者负责 jx_CSRMatrixDestroy。 */
int spgemm_adapter_multiply(
    jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C);

/* 指定算法的变体 */
int spgemm_adapter_multiply_ex(
    jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C,
    spgemm_algo_t algo);

/* 判断算法是否可用（DSP 算法需要硬件和内核文件） */
int spgemm_algo_available(spgemm_algo_t algo);

/* 数值填充：C 的 i, j 已由 JXPAMG 符号分解预计算好，
 * 此函数使用 BeidouBLAS 计算数值填入 C->data。
 * 返回 0 成功，非 0 失败。 */
int spgemm_adapter_numeric_multiply(
    jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix *C);

#endif /* BEIDOUBLAS_SPGEMM_H */
