//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

//========================================================================//
//  JXPAMG BSR Matrix-Vector Structures                                   //
//========================================================================//

/*!
 *  jx_bsr_mv.h -- head files for bsr matrix-vector operation
 *  移植自FASPXX的 _bsr_mv_struct.h
 *  Date: 2025/10/08
 *  Created by zlj
 */

#include "jx_mv.h" 
#ifndef JX_BSR_MV_HEADER
#define JX_BSR_MV_HEADER

#ifndef JX_UTIL_HEADER 
#include "jx_util.h"
#endif

/*----------------------------------------------------------------*
 *                   BSR Matrix Structure                        *
 *----------------------------------------------------------------*/

/*!
 * \struct jx_BSRMatrix
 * \brief BSR矩阵结构体
 * 
 * 注意：此矩阵使用基于0的索引
 * 注意：所有维度均以块为单位（即num_rows是块行数）
 * 注意：每个子块的存储方式为行主序
 */
typedef struct
{
   JX_Real  *data;          /* 数据数组，存储块矩阵元素 */
   JX_Int   *i;            /* 行指针数组（类似CSR的ia） */
   JX_Int   *j;            /* 列索引数组（类似CSR的ja） */
   JX_BigInt   *big_j;
   JX_Int    block_size;   /* 块大小（每个块的行数和列数） */
   JX_Int    num_rows;     /* 块行数 */
   JX_Int    num_cols;     /* 块列数 */
   JX_Int    num_nonzeros; /* 块非零元个数 */
   JX_Int    owns_data;    /* 标记：是否拥有数据所有权 */

} jx_BSRMatrix;

/*----------------------------------------------------------------*
 *             BSR矩阵结构体访问宏定义                             *
 *----------------------------------------------------------------*/

#define jx_BSRMatrixData(matrix)            ((matrix)->data)
#define jx_BSRMatrixI(matrix)               ((matrix)->i)
#define jx_BSRMatrixJ(matrix)               ((matrix)->j)
#define jx_BSRMatrixBigJ(matrix)            ((matrix)->big_j)
#define jx_BSRMatrixBlockSize(matrix)       ((matrix)->block_size)
#define jx_BSRMatrixNumRows(matrix)         ((matrix)->num_rows)
#define jx_BSRMatrixNumCols(matrix)         ((matrix)->num_cols)
#define jx_BSRMatrixNumNonzeros(matrix)     ((matrix)->num_nonzeros)
#define jx_BSRMatrixOwnsData(matrix)        ((matrix)->owns_data)

/*----------------------------------------------------------------*
 *             BSR矩阵基本操作函数声明                             *
 *----------------------------------------------------------------*/

/* 创建和销毁函数 */
jx_BSRMatrix *jx_BSRMatrixCreate(JX_Int num_rows, JX_Int num_cols, 
                                  JX_Int num_nonzeros, JX_Int block_size);
JX_Int jx_BSRMatrixInitialize(jx_BSRMatrix *matrix);
JX_Int jx_BSRMatrixDestroy(jx_BSRMatrix *matrix);

/* 矩阵-向量运算 */
JX_Int jx_BSRMatrixMatvec(JX_Real alpha, jx_BSRMatrix *A, 
                           jx_Vector *x, JX_Real beta, jx_Vector *y);

/* 格式转换函数 */
jx_CSRMatrix *jx_BSRMatrixToCSRMatrix(jx_BSRMatrix *bsr_mat);
jx_BSRMatrix *jx_CSRMatrixToBSRMatrix(jx_CSRMatrix *csr_mat, JX_Int block_size);
jx_CSRMatrix* jx_BSRMatrixGetSubmatrix(jx_BSRMatrix* matrix, JX_Int subposition);

/* 工具函数 */
JX_Int jx_BSRMatrixSetDataOwner(jx_BSRMatrix *matrix, JX_Int owns_data);
JX_Int jx_BSRMatrixPrint(jx_BSRMatrix *matrix, const char *file_name);
JX_Int jx_BSRMatrixCopy(jx_BSRMatrix *src, jx_BSRMatrix *dst, JX_Int copy_data);

jx_BSRMatrix* jx_BSRMatrixRead(const char* file_name);
jx_BSRMatrix* jx_BSRMatrixRead_Binary(const char* file_name);

#endif /* JX_BSR_MV_HEADER */