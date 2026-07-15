// test_bsr_cprgmres_mixed.c - CPR-BiCGSTAB混合精度测试程序
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

// 双精度头文件
#include "jx_mv.h"
#include "jx_bsr_mv.h"
#include "jx_parbsr_mv.h"
#include "jx_cpr.h"
#include "jx_bsr_decoup.h"
#include "jx_pamg.h"
#include "jx_krylov.h"

// 单精度头文件（假设头文件前缀是jxf）
#include "jxf_mv.h"
#include "jxf_bsr_mv.h"
#include "jxf_parbsr_mv.h"
#include "jxf_cpr.h"
#include "jxf_bsr_decoup.h"
#include "jxf_pamg.h"
#include "jxf_krylov.h"

// 定义精度类型
typedef double Real_double;  // 双精度
typedef float  Real_float;   // 单精度

typedef int Int;

// 函数声明
jx_BSRMatrix* read_bsr_matrix(const char* filename, int binary);
jx_Vector* read_rhs_file(const char* filename, int expected_size);
jx_ParBSRMatrix* distribute_bsr_matrix(MPI_Comm comm, jx_BSRMatrix* A_bsr, 
                                       JX_BigInt* row_part, JX_BigInt* col_part);

// 读取右端项文件函数
jx_Vector* read_rhs_file(const char* filename, int expected_size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening RHS file: %s\n", filename);
        return NULL;
    }
    
    // 读取向量维度
    int vector_size;
    if (fscanf(file, "%d", &vector_size) != 1) {
        printf("Error reading vector size from RHS file\n");
        fclose(file);
        return NULL;
    }
    
    // 检查维度是否匹配
    if (expected_size != -1 && vector_size != expected_size) {
        printf("Warning: RHS vector size (%d) does not match expected size (%d)\n", 
               vector_size, expected_size);
    }
    
    // 创建向量
    jx_Vector* rhs = jx_SeqVectorCreate(vector_size);
    jx_SeqVectorInitialize(rhs);
    JX_Real* rhs_data = jx_VectorData(rhs);
    
    // 初始化向量为0
    for (int i = 0; i < vector_size; i++) {
        // rhs_data[i] = 0.0;
        fscanf(file, "%lf", &rhs_data[i]);
    }
    
    // 读取索引和数值
    // int index;
    // double value;
    // int count = 0;
    // while (fscanf(file, "%d %lf", &index, &value) == 2) {
    //     if (index >= 0 && index < vector_size) {
    //         rhs_data[index] = value;
    //         count++;
    //     } else {
    //         printf("Warning: Index %d out of range [0, %d)\n", index, vector_size);
    //     }
    // }
    
    // printf("Read %d non-zero entries from RHS file (total size: %d)\n", 
    //        count, vector_size);
    
    fclose(file);
    return rhs;
}

// // =============================================================
// // 精度转换函数（修正版）
// // =============================================================

// // 1. 双精度BSR矩阵转换为单精度BSR矩阵
// jxf_BSRMatrix* convert_BSR_double_to_float(jx_BSRMatrix* A_double)
// {
//     if (!A_double) return NULL;
    
//     int nb = jx_BSRMatrixBlockSize(A_double);
//     int nrows = jx_BSRMatrixNumRows(A_double);
//     int ncols = jx_BSRMatrixNumCols(A_double);
//     int nnz = jx_BSRMatrixNumNonzeros(A_double);
    
//     // 创建单精度矩阵
//     jxf_BSRMatrix* A_float = jxf_BSRMatrixCreate(nb, nrows, ncols, nnz);
//     if (!A_float) return NULL;
    
//     // 复制索引信息（这些是整数，不需要转换）
//     int* I_double = jx_BSRMatrixI(A_double);
//     int* J_double = jx_BSRMatrixJ(A_double);
//     double* data_double = jx_BSRMatrixData(A_double);
    
//     int* I_float = jxf_BSRMatrixI(A_float);
//     int* J_float = jxf_BSRMatrixJ(A_float);
//     float* data_float = jxf_BSRMatrixData(A_float);
    
//     // 复制索引
//     for (int i = 0; i <= nrows; i++) {
//         I_float[i] = I_double[i];
//     }
//     for (int i = 0; i < nnz; i++) {
//         J_float[i] = J_double[i];
//     }
    
//     // 转换数据：double → float
//     int data_size = nnz * nb * nb;
//     for (int i = 0; i < data_size; i++) {
//         data_float[i] = (float)data_double[i];
//     }
    
//     return A_float;
// }

// // 2. 单精度BSR矩阵转换为双精度BSR矩阵
// jx_BSRMatrix* convert_BSR_float_to_double(jxf_BSRMatrix* A_float)
// {
//     if (!A_float) return NULL;
    
//     int nb = jxf_BSRMatrixBlockSize(A_float);
//     int nrows = jxf_BSRMatrixNumRows(A_float);
//     int ncols = jxf_BSRMatrixNumCols(A_float);
//     int nnz = jxf_BSRMatrixNumNonzeros(A_float);
    
//     // 创建双精度矩阵
//     jx_BSRMatrix* A_double = jx_BSRMatrixCreate(nb, nrows, ncols, nnz);
//     if (!A_double) return NULL;
    
//     // 复制索引信息
//     int* I_float = jxf_BSRMatrixI(A_float);
//     int* J_float = jxf_BSRMatrixJ(A_float);
//     float* data_float = jxf_BSRMatrixData(A_float);
    
//     int* I_double = jx_BSRMatrixI(A_double);
//     int* J_double = jx_BSRMatrixJ(A_double);
//     double* data_double = jx_BSRMatrixData(A_double);
    
//     // 复制索引
//     for (int i = 0; i <= nrows; i++) {
//         I_double[i] = I_float[i];
//     }
//     for (int i = 0; i < nnz; i++) {
//         J_double[i] = J_float[i];
//     }
    
//     // 转换数据：float → double
//     int data_size = nnz * nb * nb;
//     for (int i = 0; i < data_size; i++) {
//         data_double[i] = (double)data_float[i];
//     }
    
//     return A_double;
// }

// // 3. 双精度向量转换为单精度向量
// jxf_Vector* convert_Vector_double_to_float(jx_Vector* v_double)
// {
//     if (!v_double) return NULL;
    
//     int n = jx_VectorSize(v_double);
//     double* data_double = jx_VectorData(v_double);
    
//     // 创建单精度向量
//     jxf_Vector* v_float = jxf_SeqVectorCreate(n);
//     jxf_SeqVectorInitialize(v_float);
//     float* data_float = jxf_VectorData(v_float);
    
//     // 转换数据
//     for (int i = 0; i < n; i++) {
//         data_float[i] = (float)data_double[i];
//     }
    
//     return v_float;
// }

// // 4. 单精度向量转换为双精度向量
// jx_Vector* convert_Vector_float_to_double(jxf_Vector* v_float)
// {
//     if (!v_float) return NULL;
    
//     int n = jxf_VectorSize(v_float);
//     float* data_float = jxf_VectorData(v_float);
    
//     // 创建双精度向量
//     jx_Vector* v_double = jx_SeqVectorCreate(n);
//     jx_SeqVectorInitialize(v_double);
//     double* data_double = jx_VectorData(v_double);
    
//     // 转换数据
//     for (int i = 0; i < n; i++) {
//         data_double[i] = (double)data_float[i];
//     }
    
//     return v_double;
// }

// // 5. 双精度并行BSR矩阵转换为单精度并行BSR矩阵
// jxf_ParBSRMatrix* convert_ParBSR_double_to_float(jx_ParBSRMatrix* A_double)
// {
//     if (!A_double) return NULL;

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     MPI_Comm comm = jx_ParBSRMatrixComm(A_double);
    
//     // 转换对角部分
//     jxf_BSRMatrix* diag_float = convert_BSR_double_to_float(jx_ParBSRMatrixDiag(A_double));
//     if (!diag_float) return NULL;

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 转换非对角部分
//     jxf_BSRMatrix* offd_float = NULL;
//     jx_BSRMatrix* offd_double = jx_ParBSRMatrixOffd(A_double);
//     if (offd_double) {
//         offd_float = convert_BSR_double_to_float(offd_double);
//         if (!offd_float) {
//             jxf_BSRMatrixDestroy(diag_float);
//             return NULL;
//         }
//     }

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 创建并行矩阵结构
//     jxf_ParBSRMatrix* A_float = (jxf_ParBSRMatrix*)malloc(sizeof(jxf_ParBSRMatrix));
//     if (!A_float) {
//         jxf_BSRMatrixDestroy(diag_float);
//         if (offd_float) jxf_BSRMatrixDestroy(offd_float);
//         return NULL;
//     }

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 初始化结构体
//     memset(A_float, 0, sizeof(jxf_ParBSRMatrix));

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 复制基本信息
//     A_float->comm = comm;
//     A_float->global_num_rows = jx_ParBSRMatrixGlobalNumRows(A_double);
//     A_float->global_num_cols = jx_ParBSRMatrixGlobalNumCols(A_double);
//     A_float->first_row_index = jx_ParBSRMatrixFirstRowIndex(A_double);
//     A_float->first_col_diag = jx_ParBSRMatrixFirstColDiag(A_double);
//     A_float->last_row_index = jx_ParBSRMatrixLastRowIndex(A_double);
//     A_float->last_col_diag = jx_ParBSRMatrixLastColDiag(A_double);
  
//     // 设置矩阵块
//     A_float->diag = diag_float;
//     A_float->offd = offd_float;

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 复制列映射（注意：jx_ParBSRMatrixColMapOffd 返回的是 JX_BigInt*，需要转换为 JXF_BigInt*）
//     JX_BigInt* col_map_double = jx_ParBSRMatrixColMapOffd(A_double);
//     int col_map_size = offd_double ? jx_BSRMatrixNumCols(offd_double) : 0;
//     if (col_map_size > 0 && col_map_double) {
//         A_float->col_map_offd = (JXF_BigInt*)malloc(col_map_size * sizeof(JXF_BigInt));
//         if (A_float->col_map_offd) {
//             for (int i = 0; i < col_map_size; i++) {
//                 A_float->col_map_offd[i] = (JXF_BigInt)col_map_double[i];
//             }
//         }
//     } else {
//         A_float->col_map_offd = NULL;
//     }

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);

//     // 复制分区信息（行划分）
//     JX_Int* row_starts_double = jx_ParBSRMatrixRowStarts(A_double);
//     if (row_starts_double) {
//         // 假设行划分数组有2个元素（根据注释）
//         A_float->row_starts = (JXF_Int*)malloc(2 * sizeof(JXF_Int));
//         if (A_float->row_starts) {
//             A_float->row_starts[0] = (JXF_Int)row_starts_double[0];
//             A_float->row_starts[1] = (JXF_Int)row_starts_double[1];
//             A_float->owns_row_starts = 1;  // 单精度矩阵拥有所有权
//         }
//     }

//     printf(" [%s:%d]  \n", __FUNCTION__, __LINE__);
//     fflush(stdout);
        
//     // 复制分区信息（列划分）
//     JX_Int* col_starts_double = jx_ParBSRMatrixColStarts(A_double);
//     if (col_starts_double) {
//         // 假设列划分数组有2个元素（根据注释）
//         A_float->col_starts = (JXF_Int*)malloc(2 * sizeof(JXF_Int));
//         if (A_float->col_starts) {
//             A_float->col_starts[0] = (JXF_Int)col_starts_double[0];
//             A_float->col_starts[1] = (JXF_Int)col_starts_double[1];
//             A_float->owns_col_starts = 1;  // 单精度矩阵拥有所有权
//         }
//     }
    
//     // 通信包需要重新创建，不直接复制
//     A_float->comm_pkg = NULL;
//     A_float->comm_pkgT = NULL;
    
//     // 所有权标记
//     A_float->owns_data = 1;  // 单精度矩阵拥有数据所有权
    
//     // 其他字段
//     A_float->num_nonzeros = jx_ParBSRMatrixNumNonzeros(A_double);
//     A_float->d_num_nonzeros = (JXF_Real)jx_ParBSRMatrixDNumNonzeros(A_double);
    
//     // 缓冲区字段
//     A_float->rowindices = NULL;
//     A_float->rowvalues = NULL;
//     A_float->getrowactive = 0;
    
//     return A_float;
// }

// // 6. 双精度并行向量转换为单精度并行向量
// jxf_ParVector* convert_ParVector_double_to_float(jx_ParVector* v_double)
// {
//     if (!v_double) return NULL;
    
//     MPI_Comm comm = jx_ParVectorComm(v_double);
//     JX_Int global_size = jx_ParVectorGlobalSize(v_double);
    
//     // 获取分区信息
//     JXF_Int* partitioning = NULL;
//     JX_Int* partitioning_double = jx_ParVectorPartitioning(v_double);
//     if (partitioning_double) {
//         // 假设分区数组有2个元素（根据注释）
//         partitioning = (JXF_Int*)malloc(2 * sizeof(JXF_Int));
//         if (partitioning) {
//             partitioning[0] = (JXF_Int)partitioning_double[0];
//             partitioning[1] = (JXF_Int)partitioning_double[1];
//         }
//     }
    
//     // 创建单精度并行向量
//     jxf_ParVector* v_float = jxf_ParVectorCreate(comm, global_size, partitioning);
//     if (!v_float) {
//         if (partitioning) free(partitioning);
//         return NULL;
//     }
    
//     jxf_ParVectorInitialize(v_float);
    
//     // 转换本地数据
//     jx_Vector* local_double = jx_ParVectorLocalVector(v_double);
//     jxf_Vector* local_float = jxf_ParVectorLocalVector(v_float);
    
//     double* data_double = jx_VectorData(local_double);
//     float* data_float = jxf_VectorData(local_float);
//     int n = jx_VectorSize(local_double);
    
//     for (int i = 0; i < n; i++) {
//         data_float[i] = (float)data_double[i];
//     }
    
//     if (partitioning) free(partitioning);
//     return v_float;
// }

// // 7. 单精度并行向量转换为双精度并行向量
// jx_ParVector* convert_ParVector_float_to_double(jxf_ParVector* v_float)
// {
//     if (!v_float) return NULL;
    
//     MPI_Comm comm = jxf_ParVectorComm(v_float);
//     JXF_Int global_size = jxf_ParVectorGlobalSize(v_float);
    
//     // 获取分区信息
//     JX_Int* partitioning = NULL;
//     JXF_Int* partitioning_float = jxf_ParVectorPartitioning(v_float);
//     if (partitioning_float) {
//         // 假设分区数组有2个元素（根据注释）
//         partitioning = (JX_Int*)malloc(2 * sizeof(JX_Int));
//         if (partitioning) {
//             partitioning[0] = (JX_Int)partitioning_float[0];
//             partitioning[1] = (JX_Int)partitioning_float[1];
//         }
//     }
    
//     // 创建双精度并行向量
//     jx_ParVector* v_double = jx_ParVectorCreate(comm, global_size, partitioning);
//     if (!v_double) {
//         if (partitioning) free(partitioning);
//         return NULL;
//     }
    
//     jx_ParVectorInitialize(v_double);
    
//     // 转换本地数据
//     jxf_Vector* local_float = jxf_ParVectorLocalVector(v_float);
//     jx_Vector* local_double = jx_ParVectorLocalVector(v_double);
    
//     float* data_float = jxf_VectorData(local_float);
//     double* data_double = jx_VectorData(local_double);
//     int n = jxf_VectorSize(local_float);
    
//     for (int i = 0; i < n; i++) {
//         data_double[i] = (double)data_float[i];
//     }
    
//     if (partitioning) free(partitioning);
//     return v_double;
// }

// // 8. 复制双精度向量到单精度向量（原地转换）
// int copy_Vector_double_to_float(jx_Vector* src_double, jxf_Vector* dst_float)
// {
//     if (!src_double || !dst_float) return 0;
    
//     int n_src = jx_VectorSize(src_double);
//     int n_dst = jxf_VectorSize(dst_float);
    
//     if (n_src != n_dst) return 0;
    
//     double* data_src = jx_VectorData(src_double);
//     float* data_dst = jxf_VectorData(dst_float);
    
//     for (int i = 0; i < n_src; i++) {
//         data_dst[i] = (float)data_src[i];
//     }
    
//     return 1;
// }

// // 9. 复制单精度向量到双精度向量（原地转换）
// int copy_Vector_float_to_double(jxf_Vector* src_float, jx_Vector* dst_double)
// {
//     if (!src_float || !dst_double) return 0;
    
//     int n_src = jxf_VectorSize(src_float);
//     int n_dst = jx_VectorSize(dst_double);
    
//     if (n_src != n_dst) return 0;
    
//     float* data_src = jxf_VectorData(src_float);
//     double* data_dst = jx_VectorData(dst_double);
    
//     for (int i = 0; i < n_src; i++) {
//         data_dst[i] = (double)data_src[i];
//     }
    
//     return 1;
// }

// // 10. 并行向量的复制版本
// int copy_ParVector_double_to_float(jx_ParVector* src_double, jxf_ParVector* dst_float)
// {
//     if (!src_double || !dst_float) return 0;
    
//     jx_Vector* local_src = jx_ParVectorLocalVector(src_double);
//     jxf_Vector* local_dst = jxf_ParVectorLocalVector(dst_float);
    
//     return copy_Vector_double_to_float(local_src, local_dst);
// }

// int copy_ParVector_float_to_double(jxf_ParVector* src_float, jx_ParVector* dst_double)
// {
//     if (!src_float || !dst_double) return 0;
    
//     jxf_Vector* local_src = jxf_ParVectorLocalVector(src_float);
//     jx_Vector* local_dst = jx_ParVectorLocalVector(dst_double);
    
//     return copy_Vector_float_to_double(local_src, local_dst);
// }


// =============================================================
// BSR精度转换函数（极简版）
// =============================================================

// 1. 双精度BSR矩阵转换为单精度BSR矩阵
JXF_Int
jxmp_BSRMatrixDtoF(jx_BSRMatrix *A, jxf_BSRMatrix *B, JXF_Int copy_data)
{
    if (!A || !B) return -1;
    
    JXF_Int ierr = 0;
    JXF_Int num_rows = jx_BSRMatrixNumRows(A);
    JXF_Int nnz = jx_BSRMatrixNumNonzeros(A);
    JXF_Int blk_size = jx_BSRMatrixBlockSize(A);
    
    JX_Int *A_i = jx_BSRMatrixI(A);
    JX_Int *A_j = jx_BSRMatrixJ(A);
    JX_Real *A_data = jx_BSRMatrixData(A);
    
    JXF_Int *B_i = jxf_BSRMatrixI(B);
    JXF_Int *B_j = jxf_BSRMatrixJ(B);
    JXF_Real *B_data = jxf_BSRMatrixData(B);
    
    // 只需要复制索引（整数，可能类型相同）
    for (JXF_Int i = 0; i <= num_rows; i++) {
        B_i[i] = (JXF_Int)A_i[i];
    }
    
    for (JXF_Int i = 0; i < nnz; i++) {
        B_j[i] = (JXF_Int)A_j[i];
    }
    
    // 如果需要，转换浮点数数据
    if (copy_data && A_data && B_data) {
        JXF_Int data_size = nnz * blk_size * blk_size;
        #pragma omp parallel for
        for (JXF_Int i = 0; i < data_size; i++) {
            B_data[i] = (JXF_Real)A_data[i];  // double -> float
        }
    }
    
    return ierr;
}

// 2. 双精度并行BSR矩阵转换为单精度并行BSR矩阵
JXF_Int
jxmp_ParBSRMatrixDtoF(jx_ParBSRMatrix *A, jxf_ParBSRMatrix *B, JXF_Int copy_data)
{
    if (!A || !B) return -1;
    
    JXF_Int ierr = 0;
    JXF_Int i;
    
    // 转换对角部分
    if (A->diag && B->diag) {
        ierr = jxmp_BSRMatrixDtoF(A->diag, B->diag, copy_data);
        if (ierr != 0) return ierr;
    }
    
    // 转换非对角部分
    if (A->offd && B->offd) {
        ierr = jxmp_BSRMatrixDtoF(A->offd, B->offd, copy_data);
        if (ierr != 0) return ierr;
    }
    
    // 复制列映射（整数数组）
    if (A->col_map_offd && B->col_map_offd) {
        JXF_Int num_cols_offd = A->offd ? jx_BSRMatrixNumCols(A->offd) : 0;
        for (i = 0; i < num_cols_offd; i++) {
            B->col_map_offd[i] = (JXF_BigInt)A->col_map_offd[i];
        }
    }
    
    // 注意：row_starts和col_starts已经在jxf_ParBSRMatrixCreate中设置好了
    // 不需要在这里再次复制，除非需要类型转换
    
    // 但为了完整性，如果需要整型类型转换：
    // if (A->row_starts && B->row_starts) {
    //     // 假设长度相同（都是nprocs+1）
    //     // 这里我们不知道长度，所以这是一个潜在的问题
    //     // 实际上，在Create时已经复制了，所以可能不需要这里
    // }
    
    return ierr;
}

// 这些函数已经存在于solver_poisson_F.c中，可以直接使用

// 双精度并行向量 -> 单精度并行向量
JXF_Int
jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y)
{
    if (!x || !y) return -1;
    
    jx_Vector *x_local = jx_ParVectorLocalVector(x);
    jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
    
    return jxmp_SeqVectorDtoF(x_local, y_local);
}

// 单精度并行向量 -> 双精度并行向量  
JXF_Int
jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y)
{
    if (!x || !y) return -1;
    
    jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
    jx_Vector *y_local = jx_ParVectorLocalVector(y);
    
    return jxmp_SeqVectorFtoD(x_local, y_local);
}

// 序列向量的转换（CSR版本已有，可直接使用）
JXF_Int
jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y)
{
    if (!x || !y) return -1;
    
    JX_Real *x_data = jx_VectorData(x);
    JXF_Real *y_data = jxf_VectorData(y);
    JX_Int size = jx_VectorSize(x);
    JX_Int i;

    size *= jx_VectorNumVectors(x);

    for (i = 0; i < size; i++) {
        y_data[i] = (JXF_Real)x_data[i];
    }

    return 0;
}

JXF_Int
jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y)
{
    if (!x || !y) return -1;
    
    JXF_Real *x_data = jxf_VectorData(x);
    JX_Real *y_data = jx_VectorData(y);
    JX_Int size = jxf_VectorSize(x);
    JX_Int i = 0;

    size *= jxf_VectorNumVectors(x);

    for (i = 0; i < size; i++) {
        y_data[i] = (JX_Real)x_data[i];
    }

    return 0;
}

// 全局变量：预分配的双精度工作向量（用于混合精度CPR预条件器的Stage2）
static jx_ParVector* g_double_r = NULL;
static jx_ParVector* g_double_x = NULL;
static jx_ParVector* g_double_w = NULL;
static jx_ParBSRMatrix* g_double_A = NULL;



// =============================================================
// 混合精度CPR预条件器（单精度接口，双精度Stage2磨光）
// stage2_type=4: 双精度BGS；stage2_type=5: 双精度HSGS
// =============================================================

JXF_Int JXF_CPRPrecond_MxP(jxf_CPRPrecond *cpr,
                          jxf_ParBSRMatrix* par_matrix,
                          jxf_ParVector *par_rhs,
                          jxf_ParVector *par_sol)
{
    if (cpr == NULL || !cpr->is_initialized) return jxf_error_flag;
    if (par_rhs == NULL || par_sol == NULL) return jxf_error_flag;

    MPI_Comm comm = cpr->comm;
    int myid;
    MPI_Comm_rank(comm, &myid);
    JXF_Real t1, t2;
    JXF_Int i;

    JXF_Int block_size = cpr->block_size;
    JXF_Int pressure_index = cpr->pressure_index;

    // Stage1: 单精度AMG求解压力子系统
    if (jxf_RestrictPressureVector(par_rhs, cpr->rp, block_size, pressure_index) != JXF_SUCCESS)
        return jxf_error_flag;
    jxf_ParVectorSetConstantValues(cpr->xp, 0.0);
    t1 = jxf_MPI_Wtime();
    for (i = 0; i < cpr->stage1_maxit; i++) {
        JXF_PAMGSolve(cpr->stage1_solver,
                     (JXF_ParCSRMatrix)cpr->A_pressure,
                     (JXF_Vector)cpr->rp,
                     (JXF_Vector)cpr->xp);
    }
    t2 = jxf_MPI_Wtime();
    cpr->stage1_solve_time += t2 - t1;
    jxf_ProlongatePressureVector(cpr->xp, par_sol, block_size, pressure_index, 0);

    // Stage2: 双精度块磨光
    t1 = jxf_MPI_Wtime();

    if (g_double_A != NULL && g_double_r != NULL && g_double_x != NULL && g_double_w != NULL) {
        // 复制解和残差到双精度向量
        jxmp_ParVectorFtoD(par_sol, g_double_x);
        jxmp_ParVectorFtoD(par_rhs, g_double_r);

        JX_Real relax_weight = 1.0;
        JX_Real omega = 0.8;
        JX_Int forward_or_backward = 1;

        for (i = 0; i < cpr->stage2_maxit; i++) {
            if (cpr->stage2_solver_type == 4)
                jx_ParBSRHGSRelax(g_double_A, g_double_x, g_double_r, g_double_w,
                                 relax_weight, omega, forward_or_backward);
            else if (cpr->stage2_solver_type == 5)
                jx_ParBSRHSGSRelax(g_double_A, g_double_x, g_double_r, g_double_w,
                                  relax_weight, omega, forward_or_backward);
        }

        // 转回单精度
        jxmp_ParVectorDtoF(g_double_x, par_sol);
    }

    t2 = jxf_MPI_Wtime();
    cpr->stage2_solve_time += t2 - t1;

    return JXF_SUCCESS;
}

// =============================================================
// 迭代精化(IR)主函数
// =============================================================

int solve_with_ir_mixed_precision(
    jx_ParBSRMatrix* A_double,     // 双精度矩阵
    jx_ParVector* b_double,        // 双精度右端项
    jx_ParVector* x_double,        // 双精度解（输入/输出）
    int max_ir_iterations,         // 最大IR迭代次数
    Real_double ir_tolerance,      // IR收敛容差
    int inner_max_iterations,      // 内部单精度求解器最大迭代
    Real_float inner_tolerance,    // 内部求解器容差
    int stage2_type                // Stage2磨光类型: 2=BGS(单精), 4=双精BGS, 5=双精HSGS
)
{

    int myid;
    MPI_Comm_rank(A_double->comm, &myid);
    
        // 获取矩阵信息
    int blk_size = jx_ParBSRMatrixBlockSize(A_double);
    int local_nrows = jx_ParBSRMatrixNumRows(A_double);
    int global_nrows = jx_ParBSRMatrixGlobalNumRows(A_double);
    int local_scalar_size = local_nrows * blk_size;
    int global_scalar_size = global_nrows * blk_size;

    // 1. 将矩阵转换为单精度（仅需一次）
    // jxf_ParBSRMatrix* A_float = convert_ParBSR_double_to_float(A_double);
    // if (!A_float) {
    //     if (myid == 0) printf("Error converting matrix to single precision\n");
    //     return -1;
    // }

    // 1. 创建单精度矩阵（使用双精度矩阵的结构信息）
    double IR_cprgmres_start = MPI_Wtime();
    jxf_ParBSRMatrix* A_float = jxf_ParBSRMatrixCreate(
        A_double->comm,
        jx_ParBSRMatrixBlockSize(A_double),
        (JXF_BigInt)jx_ParBSRMatrixGlobalNumRows(A_double),
        (JXF_BigInt)jx_ParBSRMatrixGlobalNumCols(A_double),
        A_double->row_starts,  // 直接使用相同的分区
        A_double->col_starts,
        A_double->offd ? jx_BSRMatrixNumCols(A_double->offd) : 0,
        A_double->diag ? jx_BSRMatrixNumNonzeros(A_double->diag) : 0,
        A_double->offd ? jx_BSRMatrixNumNonzeros(A_double->offd) : 0
    );

    if (!A_float) {
        if (myid == 0) printf("Error converting matrix to single precision\n");
        return -1;
    }

    // 初始化矩阵（分配内部数组）
    jxf_ParBSRMatrixInitialize(A_float);

    // 2. 只转换数值数据
    if (jxmp_ParBSRMatrixDtoF(A_double, A_float, 1) != 0) {
        if (myid == 0) printf("Error converting matrix to single precision\n");
        return -1;
    }


    // 2. 创建单精度工作向量
    // JX_BigInt global_size = jxf_ParVectorGlobalSize(
    //     convert_ParVector_double_to_float(b_double));
    JX_Int  *partitioning;
    jxf_ParBSRMatrixGetRowPartitioning(A_float, &partitioning);
    // 创建时使用相同的分区
    jxf_ParVector* r_float = jxf_ParVectorCreate(
        A_float->comm,
        global_scalar_size,
        partitioning  // 使用与双精度向量相同的分区
    );

    jxf_ParVector* x_float = jxf_ParVectorCreate(
        A_double->comm,
        global_scalar_size,
        partitioning  // 使用与双精度向量相同的分区
    );

    jxf_ParVectorInitialize(r_float);
    jxf_ParVectorInitialize(x_float);


    // jxf_ParVector* b_float = jxf_ParVectorCreate(comm, global_scalar_size, NULL);
    // jxf_ParVector* x_float = jxf_ParVectorCreate(comm, global_scalar_size, NULL);
    // jxf_ParVector* r_float = jxf_ParVectorCreate(comm, global_scalar_size, NULL);
    
    if (!x_float || !r_float) {
        if (myid == 0) printf("Error creating single precision vectors\n");
        // 清理内存...
        return -1;
    }
    
    // jxf_ParVectorInitialize(b_float);
    // jxf_ParVectorInitialize(x_float);
    // jxf_ParVectorInitialize(r_float);
    
    // 开始计时
        // BiCGSTAB设置阶段
    // double IR_cprgmres_start = MPI_Wtime(); // moved earlier


    // 在IR主循环外作setup？

     // ------------------------------------------------------------
        // 3. 创建CPR预条件器
        // ------------------------------------------------------------
        if (myid == 0) {
            printf("\nCreating CPR preconditioner...\n");
        }
        
        jxf_CPRPrecond* cpr = JXF_CPRCreate(MPI_COMM_WORLD);
        
        // 设置CPR参数
        JXF_CPRSetParameter(cpr, "pressure_index", 0);    // 压力变量索引
        JXF_CPRSetParameter(cpr, "block_size", blk_size);
        JXF_CPRSetParameter(cpr, "stage1_maxit", 1);      // 阶段1迭代次数
        JXF_CPRSetParameter(cpr, "stage2_maxit", 2);      // 阶段2迭代次数
        JXF_CPRSetParameter(cpr, "stage1_solver_type", 1); // AMG求解器
        JXF_CPRSetParameter(cpr, "stage2_solver_type", stage2_type);
        JXF_CPRSetParameter(cpr, "print_level", 1); // 打印等级
        
        JXF_CPRSetRealParameter(cpr, "threshold", 1e-15); // 提取压力矩阵阈值
        
        // 设置CPR
        double setup_start = MPI_Wtime();
        JXF_Int setup_result = JXF_CPRSetup(cpr, A_float);
        double setup_end = MPI_Wtime();
        
        if (setup_result != JXF_SUCCESS) {
            if (myid == 0) printf("Error setting up CPR preconditioner\n");
            // _exit(0);
            return 1;
        }
           
        if (myid == 0) {
            printf("CPR setup completed in %.6f seconds\n", setup_end - setup_start);
        }


        // ------------------------------------------------------------
        // 4. 创建BiCGSTAB求解器并设置CPR为预条件器
        // ------------------------------------------------------------
        if (myid == 0) {
            printf("\nCreating BiCGSTAB solver with CPR preconditioner...\n");
        }
        fflush(stdout); 
        
        // BiCGSTAB参数
                int max_iter = inner_max_iterations;      // 最大迭代次数
        float tol = inner_tolerance;          // 收敛容差
        int print_level = 1;      // 打印级别
        int is_check_restarted = 1; // 检查重启

        // 创建BiCGSTAB求解器
        // BSR BiCGSTAB
        JXF_Solver bicgstab_solver;
        jxf_BiCGSTABFunctions *fg = jxf_BiCGSTABFunctionsCreate(
            jxf_ParBSRKrylovCreateVector, jxf_ParBSRKrylovDestroyVector,
            jxf_ParBSRKrylovMatvecCreate, jxf_ParBSRKrylovMatvec, jxf_ParBSRKrylovMatvecDestroy,
            jxf_ParBSRKrylovInnerProd, jxf_ParBSRKrylovCopyVector,
            jxf_ParBSRKrylovClearVector, jxf_ParBSRKrylovScaleVector, jxf_ParBSRKrylovAxpy,
            jxf_ParBSRKrylovCommInfo,
            JXF_CPRSetup, JXF_CPRPrecond);
        bicgstab_solver = (JXF_Solver)jxf_BiCGSTABCreate(fg);

        // 设置BiCGSTAB参数
                        JXF_BiCGSTABSetMaxIter(bicgstab_solver, max_iter);
        JXF_BiCGSTABSetTol(bicgstab_solver, tol);
        JXF_BiCGSTABSetLogging(bicgstab_solver, 1);
        JXF_BiCGSTABSetPrintLevel(bicgstab_solver, print_level);

        // 设置预条件器（根据stage2_type选择）
        if (stage2_type == 4 || stage2_type == 5) {
            // 混合精度预条件器：单精AMG + 双精BGS
            JXF_CPRSetParameter(cpr, "stage2_solver_type", stage2_type);
            JXF_BiCGSTABSetPrecond(bicgstab_solver,
                              (JXF_PtrToSolverFcn)JXF_CPRPrecond_MxP,
                              (JXF_PtrToSolverFcn)JXF_CPRSetup,
                              (JXF_Solver)cpr);
        } else {
            // 原始单精度预条件器
            JXF_CPRSetParameter(cpr, "stage2_solver_type", stage2_type);
            JXF_BiCGSTABSetPrecond(bicgstab_solver,
                              (JXF_PtrToSolverFcn)JXF_CPRPrecond,
                              (JXF_PtrToSolverFcn)JXF_CPRSetup,
                              (JXF_Solver)cpr);
        }
        
        // BiCGSTAB设置阶段
        double gmres_setup_start = MPI_Wtime();
        JXF_BiCGSTABSetup(bicgstab_solver, (JXF_Matrix)A_float, 
                      (JXF_Vector)r_float, (JXF_Vector)x_float);
        double gmres_setup_end = MPI_Wtime();
        
        if (myid == 0) {
            printf("BiCGSTAB setup completed in %.6f seconds\n", 
                   gmres_setup_end - gmres_setup_start);
            fflush(stdout);  
        }

    // 3. 创建双精度残差向量和工作向量
    jx_ParVector* r_double = jx_ParVectorCreate(A_double->comm, global_scalar_size, partitioning);
    jx_ParVector* dx_double = jx_ParVectorCreate(A_double->comm, global_scalar_size, partitioning);
    if (!r_double || !dx_double) {
        if (myid == 0) printf("Error creating double precision work vectors\n");
        return -1;
    }
    jx_ParVectorInitialize(r_double);
    jx_ParVectorInitialize(dx_double);

    // 预分配全局双精度向量（用于混合精度预条件器stage2）
    g_double_r = jx_ParVectorCreate(A_double->comm, global_scalar_size, partitioning);
    g_double_x = jx_ParVectorCreate(A_double->comm, global_scalar_size, partitioning);
    g_double_w = jx_ParVectorCreate(A_double->comm, global_scalar_size, partitioning);
    g_double_A = A_double;
    if (!g_double_r || !g_double_x || !g_double_w) {
        if (myid == 0) printf("Error creating global double work vectors\n");
        return -1;
    }
    jx_ParVectorInitialize(g_double_r);
    jx_ParVectorInitialize(g_double_x);
    jx_ParVectorInitialize(g_double_w);

        
    // 4. IR主循环
    Real_double initial_residual = 0.0;
    Real_double current_residual = 0.0;
    
    for (int ir_iter = 0; ir_iter < max_ir_iterations; ir_iter++) {
        
        // 4.1 计算双精度残差：r = b - A*x
        jx_ParVectorCopy(b_double, r_double);
        jx_ParBSRMatrixMatvec(-1.0, A_double, x_double, 1.0, r_double);
        
        // 4.2 计算残差范数
        current_residual = jx_ParVectorNorm2(r_double);
        if (ir_iter == 0) {
            initial_residual = current_residual;
            if (myid == 0) {
                printf("IR Iteration %2d: Initial residual = %.6e\n", 
                       ir_iter, initial_residual);
            }
        }
        
        // 4.3 检查收敛
        if (current_residual / initial_residual < ir_tolerance) {
            if (myid == 0) {
                printf("IR converged after %d iterations\n", ir_iter);
                printf("Final relative residual: %.6e\n", 
                       current_residual / initial_residual);
            }
            break;
        }
        
        if (myid == 0) {
            printf("IR Iteration %2d: Residual = %.6e (rel = %.6e)\n",
                   ir_iter, current_residual, current_residual / initial_residual);
        }
        
        // 4.4 将残差转换为单精度
        // copy_ParVector_double_to_float(r_double, r_float);

        // 在IR循环中转换残差：
        // 只需要复制数值，结构已经相同
        if (jxmp_ParVectorDtoF(r_double, r_float) != 0) {
            // 错误处理
        }

        
        // 4.5 使用单精度求解器求解：A_float * dx_float = r_float
        jxf_ParVectorSetConstantValues(x_float, 0.0);  // 清零
        
        // 这里调用您的单精度CPR-BiCGSTAB求解器
        // // 您需要实现一个类似的单精度求解函数
        // int inner_result = solve_cpr_gmres_single_precision(
        //     A_float, r_float, x_float, 
        //     inner_max_iterations, inner_tolerance);
        

        
        // ------------------------------------------------------------
        // 5. 使用CPR-BiCGSTAB求解线性系统
        // ------------------------------------------------------------
        // if (myid == 0) {
        //     printf("\nSolving linear system with CPR-BiCGSTAB...\n");
        //     printf("BiCGSTAB parameters:\n");
        //     printf("  Krylov subspace dimension: %d\n", k_dim);
        //     printf("  Maximum iterations: %d\n", max_iter);
        //     printf("  Tolerance: %.1e\n", tol);
        //     printf("\n");
        //     fflush(stdout);  
        // }
        
        // // 初始残差
        // jxf_ParVector* residual = jxf_ParVectorCreate(
        //     MPI_COMM_WORLD, global_scalar_size, NULL);
        // jxf_ParVectorInitialize(residual);

        // jxf_ParVectorCopy(par_rhs, residual);
        // jxf_ParBSRMatrixMatvec(-1.0, A_parbsr, par_sol, 1.0, residual);

        // Real init_res_norm = jxf_ParVectorNorm2(residual);
        // if (myid == 0) {
        //     printf("Initial residual norm: %.6e\n", init_res_norm);
        // }

        // 求解
        // double solve_start = MPI_Wtime();
        JXF_BiCGSTABSolve(bicgstab_solver, (JXF_Matrix)A_float,
                      (JXF_Matrix)A_float, (JXF_Vector)r_float, (JXF_Vector)x_float);
        // double solve_end = MPI_Wtime();

        // // 最终残差
        // jxf_ParVectorCopy(par_rhs, residual);
        // jxf_ParBSRMatrixMatvec(-1.0, A_parbsr, par_sol, 1.0, residual);
        // Real final_res_norm = jxf_ParVectorNorm2(residual);
        
        // // 获取求解器统计信息
        // JXF_Int num_iterations = 0;
        // JXF_Real final_rel_res_norm = 0.0;
        
        // JXF_BiCGSTABGetNumIterations(bicgstab_solver, &num_iterations);
        // JXF_BiCGSTABGetFinalRelativeResidualNorm(bicgstab_solver, &final_rel_res_norm);


        // if (inner_result != 0) {
        //     if (myid == 0) printf("Inner solver failed at IR iteration %d\n", ir_iter);
        //     break;
        // }
        
        // 4.6 将单精度解转换为双精度并更新解
        if (jxmp_ParVectorFtoD(x_float, dx_double) != 0) {
            if (myid == 0) printf("Error converting correction to double precision\n");
            break;
        }
        
        // 更新解：x_double = x_double + dx_double
        jx_ParVectorAxpy(1.0, dx_double, x_double);
    }
    
    double IR_cprgmres_end = MPI_Wtime();
    
    if (myid == 0) {
        printf("IR_cprgmres completed in %.6f seconds\n", 
                IR_cprgmres_end - IR_cprgmres_start);
        fflush(stdout);  
    }


    // 5. 清理内存
    jxf_ParBSRMatrixDestroy(A_float);
    jxf_ParVectorDestroy(x_float);
    jxf_ParVectorDestroy(r_float);
    jx_ParVectorDestroy(r_double);
    jx_ParVectorDestroy(dx_double);
    if (g_double_r) { jx_ParVectorDestroy(g_double_r); g_double_r = NULL; }
    if (g_double_x) { jx_ParVectorDestroy(g_double_x); g_double_x = NULL; }
    if (g_double_w) { jx_ParVectorDestroy(g_double_w); g_double_w = NULL; }
    g_double_A = NULL;
    JXF_BiCGSTABDestroy(bicgstab_solver);
    JXF_CPRDestroy(&cpr);
    
    return 0;
}

// =============================================================
// 单精度CPR-BiCGSTAB求解函数（需要您根据现有代码适配）
// =============================================================

int solve_cpr_gmres_single_precision(
    jxf_ParBSRMatrix* A_float,
    jxf_ParVector* b_float,
    jxf_ParVector* x_float,
    int max_iterations,
    Real_float tolerance)
{
    // 这里是您的单精度CPR-BiCGSTAB求解器实现
    // 您可以复制现有的test_bsr_cprgmres.c代码，但：
    // 1. 将所有jx_前缀改为jxf_
    // 2. 将所有Real/double改为float
    // 3. 简化参数设置，因为这是内部求解器
    
    // 由于时间关系，我先提供一个框架
    int myid;
    MPI_Comm_rank(A_float->comm, &myid);
    
    if (myid == 0) {
        printf("  [Inner single-precision solver] Starting...\n");
    }
    
    // TODO: 实现单精度CPR-BiCGSTAB
    // 1. 创建单精度CPR预条件器
    // 2. 创建单精度BiCGSTAB求解器
    // 3. 设置参数
    // 4. 求解
    
    // 简化：假设成功
    jxf_ParVectorCopy(b_float, x_float);  // 占位，实际应调用求解器
    
    return 0;
}

// =============================================================
// 主函数（混合精度版本）
// =============================================================

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int myid, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    jx_Vector       *b_Ser;
    jx_Vector       *x_Ser;

    jx_ParVector     *par_rhs;
    jx_ParVector     *par_sol;
    
    if (argc < 3) {
        if (myid == 0) {
            printf("Usage: %s <matrix_file> <matrix_type> [binary] [rhs_file]\n", argv[0]);
            printf("  matrix_type: 0=CSR, 1=BSR\n");
            printf("  binary: 0=text, 1=binary (default: 0)\n");
            printf("  rhs_file: right-hand side file (optional)\n");
            printf("  stage2_type: 2=BGS(单精), 4=双精BGS, 5=双精HSGS (default: 2)\n");
        }
        _exit(0);
        return 1;
    }
    
    char* filename = argv[1];
    int matrix_type = atoi(argv[2]);
    int binary = (argc > 3) ? atoi(argv[3]) : 0;
    char* rhs_file = (argc > 4) ? argv[4] : NULL;
    int stage2_type = (argc > 5) ? atoi(argv[5]) : 2;

    // 测试BSR矩阵的CPR-BiCGSTAB
    if (myid == 0) {
        printf("\n=================================================\n");
        printf("Testing CPR-BiCGSTAB for BSR Matrix\n");
        printf("Matrix file: %s\n", filename);
        printf("Matrix type: BSR\n");
        printf("Binary format: %s\n", binary ? "yes" : "no");
        if (rhs_file) printf("RHS file: %s\n", rhs_file);
        printf("=================================================\n\n");
    }
    
    // ------------------------------------------------------------
    // 1. 读取并分发BSR矩阵
    // ------------------------------------------------------------
    jx_BSRMatrix* A_bsr = NULL;
    if (myid == 0) {
        if (binary) {
            A_bsr = jx_BSRMatrixRead_Binary(filename);
        } else {
            A_bsr = jx_BSRMatrixRead(filename);
        }
        
        if (!A_bsr) {
            printf("Error reading BSR matrix from %s\n", filename);
            _exit(1);
        }
        
        printf("BSR Matrix Info (myid 0):\n");
        printf("  Block size: %d\n", jx_BSRMatrixBlockSize(A_bsr));
        printf("  Number of block rows: %d\n", jx_BSRMatrixNumRows(A_bsr));
        printf("  Number of block cols: %d\n", jx_BSRMatrixNumCols(A_bsr));
        printf("  Number of block nonzeros: %d\n", jx_BSRMatrixNumNonzeros(A_bsr));
        
        int blk_size = jx_BSRMatrixBlockSize(A_bsr);
        int nrows = jx_BSRMatrixNumRows(A_bsr);
        int ncols = jx_BSRMatrixNumCols(A_bsr);
        int Scalar_nrows = nrows * blk_size;
        printf("  Scalar dimensions: %d x %d\n", nrows * blk_size, ncols * blk_size);

        // 将矩阵保存成二进制格式
        // jx_BSRMatrixPrint_Binary(A_bsr,filename_output);
        
        if (rhs_file) {
            // 从文件读取右端项
            printf("Reading RHS from file: %s\n", rhs_file);
            b_Ser = read_rhs_file(rhs_file, Scalar_nrows);
            if (!b_Ser) {
                printf("Error reading RHS file, using default RHS (all ones)\n");
                // 如果读取失败，使用默认的1向量
                b_Ser = jx_SeqVectorCreate(Scalar_nrows);
                jx_SeqVectorInitialize(b_Ser);
                JX_Real *b_data = jx_VectorData(b_Ser);
                for (int i = 0; i < Scalar_nrows; i++) {
                    b_data[i] = 1.0;
                }
            }
        } else {
            // 使用默认的1向量
            b_Ser = jx_SeqVectorCreate(Scalar_nrows);
            jx_SeqVectorInitialize(b_Ser);
            JX_Real *b_data = jx_VectorData(b_Ser);
            for (int i = 0; i < Scalar_nrows; i++) {
                b_data[i] = 1.0;
            }
            printf("Using default RHS (all ones)\n");
        }
        
        // 创建初始解向量（全零）
        x_Ser = jx_SeqVectorCreate(Scalar_nrows);
        jx_SeqVectorInitialize(x_Ser);
        JX_Real *x_data = jx_VectorData(x_Ser);
        for (int i = 0; i < Scalar_nrows; i++) {
            x_data[i] = 0.0;
        }
    }
    
        // ------------------------------------------------------------
        // 2. 在串行数据上进行解耦
        // ------------------------------------------------------------
    if (myid == 0){

        printf("\nApplying decoupling to serial BSR system...\n");
        
        // 选择解耦方法
        JX_DecoupType decoup_type = JX_DECOUP_TIMPES;  // 推荐使用ABF
        // JX_DecoupType decoup_type = JX_DECOUP_ANL;  // 推荐使用ABF
        JX_Int is_thermal = 0;  // 假设非热力模型
        
        // printf("DEBUG: decoup_type value = %d\n", decoup_type);
        // printf("DEBUG: JX_DECOUP_TIMPES macro value = %d\n", JX_DECOUP_TIMPES);

        printf("  Decoupling method: %s\n", jx_DecoupTypeToString(decoup_type));
        printf("  Is thermal model: %s\n", is_thermal ? "yes" : "no");
        
        // 备份原始矩阵（如果需要保留原始数据）
        // jx_BSRMatrix* A_bsr_orig = jx_BSRMatrixDuplicate(A_bsr);
        // jx_Vector* b_Ser_orig = jx_SeqVectorDuplicate(b_Ser);
        
        // 执行解耦（同时修改矩阵和右端项）
        double decoup_start = MPI_Wtime();
        JX_Int decoup_result = jx_BSRMatrixDecouple(A_bsr, b_Ser, decoup_type, is_thermal);
        double decoup_end = MPI_Wtime();
        
        if (decoup_result != JX_SUCCESS) {
            printf("Error applying decoupling to serial system\n");
            // cleanup error, but solver may have completed
        _exit(1);
        }
        
        printf("  Decoupling completed in %.6f seconds\n", decoup_end - decoup_start);
        printf("  Matrix and RHS have been modified by decoupling\n");
        
        // 可选：验证解耦效果
        // 可以检查解耦后的对角块是否接近单位阵
        if (1) {
            int nb = jx_BSRMatrixBlockSize(A_bsr);
            JX_Real* data = jx_BSRMatrixData(A_bsr);
            JX_Int* rpt = jx_BSRMatrixI(A_bsr);
            int nblocks = jx_BSRMatrixNumRows(A_bsr);
            
            // 检查第一个对角块
            if (nblocks > 0) {
                JX_Real* first_diag_block = data;
                printf("  First diagonal block after decoupling:\n");
                for (int i = 0; i < nb; i++) {
                    for (int j = 0; j < nb; j++) {
                        printf("    %12.6e", first_diag_block[i*nb + j]);
                    }
                    printf("\n");
                }
            }
        }
    
        // // 广播是否成功读取和解耦的标志
        // int read_success = (myid == 0) ? 1 : 0;
        // MPI_Bcast(&read_success, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // if (!read_success) {
        //     if (myid == 0) printf("Error during matrix reading/decoupling\n");
        //     _exit(0);
        //     return 1;
        // }
    }

    // 生成分区
    JX_BigInt* global_row_start = NULL;
    JX_BigInt* global_col_start = NULL;
    if (myid == 0) {
        jx_GeneratePartitioning(jx_BSRMatrixNumRows(A_bsr), size, &global_row_start);
        jx_GeneratePartitioning(jx_BSRMatrixNumCols(A_bsr), size, &global_col_start);
    }
    
    // 分发矩阵到各进程
    jx_ParBSRMatrix* A_parbsr = jx_BSRMatrixToParBSRMatrix(
        MPI_COMM_WORLD, A_bsr, global_row_start, global_col_start);
    
    // jx_ParBSRMatrixReorder(A_parbsr);
    // jx_ParBSRMatrixSetNumNonzeros(A_parbsr);
    
    // 清理串行矩阵内存
    // if (myid == 0) {
    //     jx_BSRMatrixDestroy(A_bsr);
    // }
    
    if (!A_parbsr) {
        if (myid == 0) printf("Error creating parallel BSR matrix\n");
        _exit(0);
        return 1;
    }
    

    
    // if (myid == 0) {
    //     printf("\nParallel BSR Matrix Info:\n");
    //     printf("  Block size: %d\n", blk_size);
    //     printf("  Global block rows: %d\n", global_nrows);
    //     printf("  Global scalar rows: %d\n", global_scalar_size);
    //     printf("\n");
    // }
    
    // ------------------------------------------------------------
    // 2. 创建右端项和初始解向量
    // ------------------------------------------------------------
    // 创建分区数组（标量行分区）
    JX_Int  *partitioning;
    jx_ParBSRMatrixGetRowPartitioning(A_parbsr, &partitioning);
    

    // 创建并行向量
    par_rhs = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), b_Ser, partitioning);
    par_sol = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), x_Ser, partitioning);
    
    if (!par_rhs || !par_sol) {
        printf("Rank %d: ERROR: Failed to create parallel vectors!\n", myid);
        // cleanup error, but solver may have completed
        _exit(1);
    }




    
    // 然后调用混合精度求解器
    if (myid == 0) {
        printf("\n=================================================\n");
        printf("Mixed Precision Iterative Refinement (IR) Solver\n");
        printf("Stage1=单精AMG, Stage2=%s\n",
               stage2_type == 4 ? "双精BGS" :
               stage2_type == 5 ? "双精HSGS" : "单精BGS");
        printf("Stage2 type: %d\n", stage2_type);
        printf("=================================================\n\n");
    }
    
    // // 假设已经有了双精度的A_double, b_double, x_double
    // jx_ParBSRMatrix* A_double = ...;  // 您的现有代码
    // jx_ParVector* b_double = ...;     // 您的现有代码
    // jx_ParVector* x_double = ...;     // 您的现有代码
    
    // IR参数
    int max_ir_iterations = 100;
    Real_double ir_tolerance = 1e-4;
    int inner_max_iterations = 3;
    Real_float inner_tolerance = 1e-2;
    
    // 调用混合精度求解
    int result = solve_with_ir_mixed_precision(
        A_parbsr, par_rhs, par_sol,
        max_ir_iterations, ir_tolerance,
        inner_max_iterations, inner_tolerance,
        stage2_type);
    
    _exit(0);
    return result;
}