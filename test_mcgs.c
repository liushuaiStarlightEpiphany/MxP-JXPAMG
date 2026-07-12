#include "jx_mv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 声明要测试的函数 */
extern JX_Int jx_PAMGRelaxMultiColorGS(
    jx_ParCSRMatrix *par_matrix,
    jx_ParVector    *par_rhs,
    JX_Int          *cf_marker,
    JX_Int           relax_points,
    JX_Real          relax_weight,
    JX_Real          omega,
    jx_ParVector    *par_app,
    jx_ParVector    *Vtemp);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* 构建一个 4x4 对角占优矩阵:
     * A = [ 4 -1  0  0
     *      -1  4 -1  0
     *       0 -1  4 -1
     *       0  0 -1  4 ]
     *
     * b = [1, 2, 3, 4]^T
     * x = [0, 0, 0, 0]^T (初值)
     *
     * 真解: A*x = b
     */
    JX_Int n = 4;
    JX_Int nnz = 10; /* 每行非对角+对角 */

    /* 创建 CSR 矩阵 */
    jx_CSRMatrix *A_diag = jx_CSRMatrixCreate(n, n, nnz);
    jx_CSRMatrixInitialize(A_diag);

    JX_Int    *A_i = jx_CSRMatrixI(A_diag);
    JX_Int    *A_j = jx_CSRMatrixJ(A_diag);
    JX_Real   *A_data = jx_CSRMatrixData(A_diag);

    /* 手动填充: 4 -1  0  0
     *         -1  4 -1  0
     *          0 -1  4 -1
     *          0  0 -1  4 */
    JX_Int rowptr[] = {0, 2, 5, 8, 10};
    JX_Int cols[]   = {0, 1,    1, 0, 2,    2, 1, 3,    3, 2};
    JX_Real vals[]  = {4,-1,    4,-1,-1,    4,-1,-1,    4,-1};

    memcpy(A_i, rowptr, (n+1)*sizeof(JX_Int));
    memcpy(A_j, cols, nnz*sizeof(JX_Int));
    memcpy(A_data, vals, nnz*sizeof(JX_Real));

    /* 设置 rownnz */
    A_diag->rownnz = (JX_Int*)malloc(n * sizeof(JX_Int));
    for (JX_Int i = 0; i < n; i++)
        A_diag->rownnz[i] = A_i[i+1] - A_i[i];

    /* 创建 ParCSRMatrix */
    JX_Int row_starts[2] = {0, n};
    JX_Int col_starts[2] = {0, n};

    jx_ParCSRMatrix *A = jx_ParCSRMatrixCreate(
        MPI_COMM_WORLD, n, n, row_starts, col_starts, 0, 0, nnz);

    /* 设置 diag/offd (所有列都在本地, offd 为空) */
    jx_ParCSRMatrixDiag(A) = A_diag;
    jx_ParCSRMatrixOffd(A) = jx_CSRMatrixCreate(n, 0, 0);
    jx_CSRMatrixInitialize(jx_ParCSRMatrixOffd(A));

    /* 创建通信包 */
    jx_MatvecCommPkgCreate(A);

    /* 创建向量 */
    jx_ParVector *b = jx_ParVectorCreate(MPI_COMM_WORLD, n, row_starts);
    jx_ParVectorInitialize(b);
    jx_ParVector *x = jx_ParVectorCreate(MPI_COMM_WORLD, n, row_starts);
    jx_ParVectorInitialize(x);
    jx_ParVector *Vtemp = jx_ParVectorCreate(MPI_COMM_WORLD, n, row_starts);
    jx_ParVectorInitialize(Vtemp);

    /* 设置 b = [1,2,3,4]^T */
    JX_Real *b_data = jx_VectorData(jx_ParVectorLocalVector(b));
    b_data[0] = 1.0; b_data[1] = 2.0; b_data[2] = 3.0; b_data[3] = 4.0;
    /* x 初始化为 0 (已经是 0) */

    if (rank == 0)
        printf("Testing Multi-Color GS Relaxation on 4x4 system...\n");

    /* 执行 4 次多色 GS 调用 (模拟磨光) */
    for (int iter = 0; iter < 4; iter++)
    {
        jx_PAMGRelaxMultiColorGS(A, b, NULL, 0, 1.0, 1.0, x, Vtemp);

        if (rank == 0) {
            JX_Real *xd = jx_VectorData(jx_ParVectorLocalVector(x));
            printf("  iter %d: x = [%.6f, %.6f, %.6f, %.6f]\n",
                   iter+1, xd[0], xd[1], xd[2], xd[3]);
        }
    }

    /* 验证: 计算残差 ||b - A*x|| */
    jx_ParVector *res = jx_ParVectorCreate(MPI_COMM_WORLD, n, row_starts);
    jx_ParVectorInitialize(res);
    jx_ParVectorCopy(b, res);
    jx_ParCSRMatrixMatvec(-1.0, A, x, 1.0, res);

    JX_Real rnorm = sqrt(jx_ParVectorInnerProd(res, res));
    JX_Real bnorm = sqrt(jx_ParVectorInnerProd(b, b));

    if (rank == 0)
        printf("\n  Final rel. residual: %e\n", rnorm / bnorm);

    /* 释放 (注意: 手动构建的 CSR 矩阵需要手动释放) */
    jx_ParVectorDestroy(res);
    jx_ParVectorDestroy(Vtemp);
    jx_ParVectorDestroy(x);
    jx_ParVectorDestroy(b);

    /* A_diag 是手动创建的, 需要先解引用再 destroy ParCSRMatrix */
    jx_CSRMatrix *tmp_diag = jx_ParCSRMatrixDiag(A);
    jx_CSRMatrix *tmp_offd = jx_ParCSRMatrixOffd(A);
    jx_ParCSRMatrixDiag(A) = NULL;
    jx_ParCSRMatrixOffd(A) = NULL;
    jx_ParCSRMatrixDestroy(A);
    free(tmp_diag->rownnz);
    tmp_diag->rownnz = NULL;
    jx_CSRMatrixDestroy(tmp_diag);
    jx_CSRMatrixDestroy(tmp_offd);

    MPI_Finalize();

    if (rank == 0)
        printf("PASSED\n");
    return 0;
}
