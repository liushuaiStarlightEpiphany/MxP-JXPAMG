//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

//========================================================================//
//  JXFPAMG BSR Matrix-Vector Structures                                   //
//========================================================================//

/*!
 *  jxf_bsr_mv.h -- head files for bsr matrix-vector operation
 *  移植自FASPXX的 _bsr_mv_struct.h
 *  Date: 2025/10/08
 *  Created by zlj
 */

#include "jxf_mv.h" 
#ifndef JXF_BSR_MV_HEADER
#define JXF_BSR_MV_HEADER

#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

/*----------------------------------------------------------------*
 *                   BSR Matrix Structure                        *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_BSRMatrix
 * \brief BSR矩阵结构体
 * 
 * 注意：此矩阵使用基于0的索引
 * 注意：所有维度均以块为单位（即num_rows是块行数）
 * 注意：每个子块的存储方式为行主序
 */
typedef struct
{
   JXF_Real  *data;          /* 数据数组，存储块矩阵元素 */
   JXF_Int   *i;            /* 行指针数组（类似CSR的ia） */
   JXF_Int   *j;            /* 列索引数组（类似CSR的ja） */
   JXF_BigInt   *big_j;
   JXF_Int    block_size;   /* 块大小（每个块的行数和列数） */
   JXF_Int    num_rows;     /* 块行数 */
   JXF_Int    num_cols;     /* 块列数 */
   JXF_Int    num_nonzeros; /* 块非零元个数 */
   JXF_Int    owns_data;    /* 标记：是否拥有数据所有权 */

} jxf_BSRMatrix;

/*----------------------------------------------------------------*
 *             BSR矩阵结构体访问宏定义                             *
 *----------------------------------------------------------------*/

#define jxf_BSRMatrixData(matrix)            ((matrix)->data)
#define jxf_BSRMatrixI(matrix)               ((matrix)->i)
#define jxf_BSRMatrixJ(matrix)               ((matrix)->j)
#define jxf_BSRMatrixBigJ(matrix)            ((matrix)->big_j)
#define jxf_BSRMatrixBlockSize(matrix)       ((matrix)->block_size)
#define jxf_BSRMatrixNumRows(matrix)         ((matrix)->num_rows)
#define jxf_BSRMatrixNumCols(matrix)         ((matrix)->num_cols)
#define jxf_BSRMatrixNumNonzeros(matrix)     ((matrix)->num_nonzeros)
#define jxf_BSRMatrixOwnsData(matrix)        ((matrix)->owns_data)

/*----------------------------------------------------------------*
 *             BSR矩阵基本操作函数声明                             *
 *----------------------------------------------------------------*/

/* 创建和销毁函数 */
jxf_BSRMatrix *jxf_BSRMatrixCreate(JXF_Int num_rows, JXF_Int num_cols, 
                                  JXF_Int num_nonzeros, JXF_Int block_size);
JXF_Int jxf_BSRMatrixInitialize(jxf_BSRMatrix *matrix);
JXF_Int jxf_BSRMatrixDestroy(jxf_BSRMatrix *matrix);

/* 矩阵-向量运算 */
JXF_Int jxf_BSRMatrixMatvec(JXF_Real alpha, jxf_BSRMatrix *A, 
                           jxf_Vector *x, JXF_Real beta, jxf_Vector *y);

/* 格式转换函数 */
jxf_CSRMatrix *jxf_BSRMatrixToCSRMatrix(jxf_BSRMatrix *bsr_mat);
jxf_BSRMatrix *jxf_CSRMatrixToBSRMatrix(jxf_CSRMatrix *csr_mat, JXF_Int block_size);
jxf_CSRMatrix* jxf_BSRMatrixGetSubmatrix(jxf_BSRMatrix* matrix, JXF_Int subposition);

/* 工具函数 */
JXF_Int jxf_BSRMatrixSetDataOwner(jxf_BSRMatrix *matrix, JXF_Int owns_data);
JXF_Int jxf_BSRMatrixPrint(jxf_BSRMatrix *matrix, const char *file_name);
JXF_Int jxf_BSRMatrixCopy(jxf_BSRMatrix *src, jxf_BSRMatrix *dst, JXF_Int copy_data);

jxf_BSRMatrix* jxf_BSRMatrixRead(const char* file_name);
jxf_BSRMatrix* jxf_BSRMatrixRead_Binary(const char* file_name);

#endif /* JXF_BSR_MV_HEADER */