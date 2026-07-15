//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

//========================================================================//
//  JXPAMG Parallel BSR Matrix-Vector Structures                          //
//========================================================================//

/*!
 *  jx_parbsr_mv.h -- head files for Parallel bsr matrix-vector operation
 *  移植自FASPXX的 _parbsr_mv_struct.h
 *  Date: 2025/10/08
 *  Created by zlj
 */

#ifndef JX_PARBSR_MV_HEADER
#define JX_PARBSR_MV_HEADER

#ifndef JX_BSR_MV_HEADER
#include "jx_bsr_mv.h"
#endif

#ifndef JX_PARCSR_MV_HEADER
#include "jx_mv.h"  /* 包含并行通信相关结构体 */
#endif

/*----------------------------------------------------------------*
 *             Parallel BSR Matrix Structure                     *
 *----------------------------------------------------------------*/

/*!
 * \struct jx_ParBSRMatrix
 * \brief 并行BSR矩阵结构体
 */
typedef struct
{
   MPI_Comm          comm;               /* MPI通信域 */
   JX_BigInt            global_num_rows;    /* 全局块行数 */
   JX_BigInt            global_num_cols;    /* 全局块列数 */
   JX_BigInt            first_row_index;    /* 本地第一行的全局索引 */
   JX_BigInt            first_col_diag;     /* 本地对角块第一列的全局索引 */
   JX_BigInt            last_row_index;     /* 本地最后一行的全局索引 */
   JX_BigInt            last_col_diag;      /* 本地对角块最后一列的全局索引 */
   
   jx_BSRMatrix     *diag;               /* 指向本地对角块矩阵 */
   jx_BSRMatrix     *offd;               /* 指向本地非对角块矩阵 */
   JX_BigInt           *col_map_offd;       /* 将offd的列映射到全局列索引 */
   
   JX_Int           *row_starts;         /* 行划分数组 [start, end) */
   JX_Int           *col_starts;         /* 列划分数组 [start, end) */
   
    // JX_BigInt row_starts[2];
    // /* row_starts[0] is start of local rows
    //     row_starts[1] is start of next processor's rows */
    // JX_BigInt col_starts[2];
    // /* col_starts[0] is start of local columns
    //     col_starts[1] is start of next processor's columns */
   jx_ParCSRCommPkg *comm_pkg;           /* 通信包（用于矩阵向量乘） */
   jx_ParCSRCommPkg *comm_pkgT;          /* 转置通信包 */
   
   JX_Int            owns_data;          /* 标记：是否拥有数据所有权 */
   JX_Int            owns_row_starts;    /* 标记：是否拥有行划分所有权 */
   JX_Int            owns_col_starts;    /* 标记：是否拥有列划分所有权 */
   
   JX_BigInt            num_nonzeros;       /* 全局非零元总数 */
   JX_Real           d_num_nonzeros;     /* 全局非零元总数（实数形式） */
   
   /* 用于GetRow函数的缓冲区 */
   JX_Int           *rowindices;
   JX_Real          *rowvalues;
   JX_Int            getrowactive;
   
   jx_IJAssumedPart *assumed_partition;  /* 假设分区信息 */

} jx_ParBSRMatrix;

/*----------------------------------------------------------------*
 *         并行BSR矩阵结构体访问宏定义                            *
 *----------------------------------------------------------------*/

#define jx_ParBSRMatrixComm(matrix)             ((matrix)->comm)
#define jx_ParBSRMatrixGlobalNumRows(matrix)    ((matrix)->global_num_rows)
#define jx_ParBSRMatrixGlobalNumCols(matrix)    ((matrix)->global_num_cols)
#define jx_ParBSRMatrixFirstRowIndex(matrix)    ((matrix)->first_row_index)
#define jx_ParBSRMatrixFirstColDiag(matrix)     ((matrix)->first_col_diag)
#define jx_ParBSRMatrixLastRowIndex(matrix)     ((matrix)->last_row_index)
#define jx_ParBSRMatrixLastColDiag(matrix)      ((matrix)->last_col_diag)
#define jx_ParBSRMatrixBlockSize(matrix)        ((matrix)->diag->block_size)
#define jx_ParBSRMatrixDiag(matrix)             ((matrix)->diag)
#define jx_ParBSRMatrixOffd(matrix)             ((matrix)->offd)
#define jx_ParBSRMatrixColMapOffd(matrix)       ((matrix)->col_map_offd)
#define jx_ParBSRMatrixRowStarts(matrix)        ((matrix)->row_starts)
#define jx_ParBSRMatrixColStarts(matrix)        ((matrix)->col_starts)
#define jx_ParBSRMatrixCommPkg(matrix)          ((matrix)->comm_pkg)
#define jx_ParBSRMatrixCommPkgT(matrix)         ((matrix)->comm_pkgT)
#define jx_ParBSRMatrixOwnsData(matrix)         ((matrix)->owns_data)
#define jx_ParBSRMatrixOwnsRowStarts(matrix)    ((matrix)->owns_row_starts)
#define jx_ParBSRMatrixOwnsColStarts(matrix)    ((matrix)->owns_col_starts)

#define jx_ParBSRMatrixNumRows(matrix)          \
    jx_BSRMatrixNumRows(jx_ParBSRMatrixDiag(matrix))
#define jx_ParBSRMatrixNumCols(matrix)          \
    jx_BSRMatrixNumCols(jx_ParBSRMatrixDiag(matrix))
    
#define jx_ParBSRMatrixNumNonzeros(matrix)      ((matrix)->num_nonzeros)
#define jx_ParBSRMatrixDNumNonzeros(matrix)     ((matrix)->d_num_nonzeros)
#define jx_ParBSRMatrixRowindices(matrix)       ((matrix)->rowindices)
#define jx_ParBSRMatrixRowvalues(matrix)        ((matrix)->rowvalues)
#define jx_ParBSRMatrixGetrowactive(matrix)     ((matrix)->getrowactive)
#define jx_ParBSRMatrixAssumedPartition(matrix) ((matrix)->assumed_partition)

/*----------------------------------------------------------------*
 *         并行BSR矩阵操作函数声明                                *
 *----------------------------------------------------------------*/

/* 创建和销毁函数 */
jx_ParBSRMatrix *jx_ParBSRMatrixCreate(MPI_Comm comm,
                                        JX_Int block_size, 
                                        JX_BigInt global_num_rows,
                                        JX_BigInt global_num_cols,
                                        JX_BigInt *row_starts,
                                        JX_BigInt *col_starts,
                                        JX_Int num_cols_offd,
                                        JX_Int num_nonzeros_diag,
                                        JX_Int num_nonzeros_offd);

JX_Int jx_ParBSRMatrixInitialize(jx_ParBSRMatrix *matrix);
JX_Int jx_ParBSRMatrixDestroy(jx_ParBSRMatrix *matrix);

/* 格式转换函数 */
jx_ParCSRMatrix *jx_ParBSRMatrixToParCSRMatrix(jx_ParBSRMatrix *par_bsr);
jx_ParBSRMatrix *jx_ParCSRMatrixToParBSRMatrix(jx_ParCSRMatrix *par_csr, 
                                                JX_Int block_size);
jx_ParCSRMatrix* jx_ParBSRMatrixGetSubmatrix(jx_ParBSRMatrix* matrix, JX_Int subposition, JX_Real threshold);

/* 矩阵-向量运算 */
JX_Int jx_ParBSRMatrixMatvec(JX_Real alpha, jx_ParBSRMatrix *A,
                              jx_ParVector *x, JX_Real beta, jx_ParVector *y);

/* 工具函数 */
JX_Int jx_ParBSRMatrixGetLocalRange(jx_ParBSRMatrix *matrix,
                                     JX_Int *row_start, JX_Int *row_end,
                                     JX_Int *col_start, JX_Int *col_end);

JX_Int jx_ParBSRMatrixSetDataOwner(jx_ParBSRMatrix *matrix, JX_Int owns_data);
JX_Int jx_ParBSRMatrixSetRowStartsOwner(jx_ParBSRMatrix *matrix, 
                                         JX_Int owns_row_starts);
JX_Int jx_ParBSRMatrixSetColStartsOwner(jx_ParBSRMatrix *matrix, 
                                         JX_Int owns_col_starts);

/* 通信相关 */
JX_Int jx_ParBSRMatrixCreateCommPkg(jx_ParBSRMatrix *matrix);

jx_ParBSRMatrix* jx_BSRMatrixToParBSRMatrix(MPI_Comm comm, jx_BSRMatrix* A, JX_BigInt* global_row_starts, JX_BigInt* global_col_starts);
jx_ParCSRCommHandle *jx_ParBSRCommHandleCreate( JX_Int job, JX_Int bnnz, jx_ParCSRCommPkg *comm_pkg, void *send_data, void  *recv_data );
JX_Int jx_ParBSRMatrixGetRowPartitioning( jx_ParBSRMatrix *matrix, JX_Int **row_partitioning_ptr );
JX_Int jx_ParBSRMatrixGetColPartitioning( jx_ParBSRMatrix *matrix, JX_Int **row_partitioning_ptr );
#endif /* JX_PARBSR_MV_HEADER */