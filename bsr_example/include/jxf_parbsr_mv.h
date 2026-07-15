//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

//========================================================================//
//  JXFPAMG Parallel BSR Matrix-Vector Structures                          //
//========================================================================//

/*!
 *  jxf_parbsr_mv.h -- head files for Parallel bsr matrix-vector operation
 *  移植自FASPXX的 _parbsr_mv_struct.h
 *  Date: 2025/10/08
 *  Created by zlj
 */

#ifndef JXF_PARBSR_MV_HEADER
#define JXF_PARBSR_MV_HEADER

#ifndef JXF_BSR_MV_HEADER
#include "jxf_bsr_mv.h"
#endif

#ifndef JXF_PARCSR_MV_HEADER
#include "jxf_mv.h"  /* 包含并行通信相关结构体 */
#endif

/*----------------------------------------------------------------*
 *             Parallel BSR Matrix Structure                     *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_ParBSRMatrix
 * \brief 并行BSR矩阵结构体
 */
typedef struct
{
   MPI_Comm          comm;               /* MPI通信域 */
   JXF_BigInt            global_num_rows;    /* 全局块行数 */
   JXF_BigInt            global_num_cols;    /* 全局块列数 */
   JXF_BigInt            first_row_index;    /* 本地第一行的全局索引 */
   JXF_BigInt            first_col_diag;     /* 本地对角块第一列的全局索引 */
   JXF_BigInt            last_row_index;     /* 本地最后一行的全局索引 */
   JXF_BigInt            last_col_diag;      /* 本地对角块最后一列的全局索引 */
   
   jxf_BSRMatrix     *diag;               /* 指向本地对角块矩阵 */
   jxf_BSRMatrix     *offd;               /* 指向本地非对角块矩阵 */
   JXF_BigInt           *col_map_offd;       /* 将offd的列映射到全局列索引 */
   
   JXF_Int           *row_starts;         /* 行划分数组 [start, end) */
   JXF_Int           *col_starts;         /* 列划分数组 [start, end) */
   
    // JXF_BigInt row_starts[2];
    // /* row_starts[0] is start of local rows
    //     row_starts[1] is start of next processor's rows */
    // JXF_BigInt col_starts[2];
    // /* col_starts[0] is start of local columns
    //     col_starts[1] is start of next processor's columns */
   jxf_ParCSRCommPkg *comm_pkg;           /* 通信包（用于矩阵向量乘） */
   jxf_ParCSRCommPkg *comm_pkgT;          /* 转置通信包 */
   
   JXF_Int            owns_data;          /* 标记：是否拥有数据所有权 */
   JXF_Int            owns_row_starts;    /* 标记：是否拥有行划分所有权 */
   JXF_Int            owns_col_starts;    /* 标记：是否拥有列划分所有权 */
   
   JXF_BigInt            num_nonzeros;       /* 全局非零元总数 */
   JXF_Real           d_num_nonzeros;     /* 全局非零元总数（实数形式） */
   
   /* 用于GetRow函数的缓冲区 */
   JXF_Int           *rowindices;
   JXF_Real          *rowvalues;
   JXF_Int            getrowactive;
   
   jxf_IJAssumedPart *assumed_partition;  /* 假设分区信息 */

} jxf_ParBSRMatrix;

/*----------------------------------------------------------------*
 *         并行BSR矩阵结构体访问宏定义                            *
 *----------------------------------------------------------------*/

#define jxf_ParBSRMatrixComm(matrix)             ((matrix)->comm)
#define jxf_ParBSRMatrixGlobalNumRows(matrix)    ((matrix)->global_num_rows)
#define jxf_ParBSRMatrixGlobalNumCols(matrix)    ((matrix)->global_num_cols)
#define jxf_ParBSRMatrixFirstRowIndex(matrix)    ((matrix)->first_row_index)
#define jxf_ParBSRMatrixFirstColDiag(matrix)     ((matrix)->first_col_diag)
#define jxf_ParBSRMatrixLastRowIndex(matrix)     ((matrix)->last_row_index)
#define jxf_ParBSRMatrixLastColDiag(matrix)      ((matrix)->last_col_diag)
#define jxf_ParBSRMatrixBlockSize(matrix)        ((matrix)->diag->block_size)
#define jxf_ParBSRMatrixDiag(matrix)             ((matrix)->diag)
#define jxf_ParBSRMatrixOffd(matrix)             ((matrix)->offd)
#define jxf_ParBSRMatrixColMapOffd(matrix)       ((matrix)->col_map_offd)
#define jxf_ParBSRMatrixRowStarts(matrix)        ((matrix)->row_starts)
#define jxf_ParBSRMatrixColStarts(matrix)        ((matrix)->col_starts)
#define jxf_ParBSRMatrixCommPkg(matrix)          ((matrix)->comm_pkg)
#define jxf_ParBSRMatrixCommPkgT(matrix)         ((matrix)->comm_pkgT)
#define jxf_ParBSRMatrixOwnsData(matrix)         ((matrix)->owns_data)
#define jxf_ParBSRMatrixOwnsRowStarts(matrix)    ((matrix)->owns_row_starts)
#define jxf_ParBSRMatrixOwnsColStarts(matrix)    ((matrix)->owns_col_starts)

#define jxf_ParBSRMatrixNumRows(matrix)          \
    jxf_BSRMatrixNumRows(jxf_ParBSRMatrixDiag(matrix))
#define jxf_ParBSRMatrixNumCols(matrix)          \
    jxf_BSRMatrixNumCols(jxf_ParBSRMatrixDiag(matrix))
    
#define jxf_ParBSRMatrixNumNonzeros(matrix)      ((matrix)->num_nonzeros)
#define jxf_ParBSRMatrixDNumNonzeros(matrix)     ((matrix)->d_num_nonzeros)
#define jxf_ParBSRMatrixRowindices(matrix)       ((matrix)->rowindices)
#define jxf_ParBSRMatrixRowvalues(matrix)        ((matrix)->rowvalues)
#define jxf_ParBSRMatrixGetrowactive(matrix)     ((matrix)->getrowactive)
#define jxf_ParBSRMatrixAssumedPartition(matrix) ((matrix)->assumed_partition)

/*----------------------------------------------------------------*
 *         并行BSR矩阵操作函数声明                                *
 *----------------------------------------------------------------*/

/* 创建和销毁函数 */
jxf_ParBSRMatrix *jxf_ParBSRMatrixCreate(MPI_Comm comm,
                                        JXF_Int block_size, 
                                        JXF_BigInt global_num_rows,
                                        JXF_BigInt global_num_cols,
                                        JXF_BigInt *row_starts,
                                        JXF_BigInt *col_starts,
                                        JXF_Int num_cols_offd,
                                        JXF_Int num_nonzeros_diag,
                                        JXF_Int num_nonzeros_offd);

JXF_Int jxf_ParBSRMatrixInitialize(jxf_ParBSRMatrix *matrix);
JXF_Int jxf_ParBSRMatrixDestroy(jxf_ParBSRMatrix *matrix);

/* 格式转换函数 */
jxf_ParCSRMatrix *jxf_ParBSRMatrixToParCSRMatrix(jxf_ParBSRMatrix *par_bsr);
jxf_ParBSRMatrix *jxf_ParCSRMatrixToParBSRMatrix(jxf_ParCSRMatrix *par_csr, 
                                                JXF_Int block_size);
jxf_ParCSRMatrix* jxf_ParBSRMatrixGetSubmatrix(jxf_ParBSRMatrix* matrix, JXF_Int subposition, JXF_Real threshold);

/* 矩阵-向量运算 */
JXF_Int jxf_ParBSRMatrixMatvec(JXF_Real alpha, jxf_ParBSRMatrix *A,
                              jxf_ParVector *x, JXF_Real beta, jxf_ParVector *y);

/* 工具函数 */
JXF_Int jxf_ParBSRMatrixGetLocalRange(jxf_ParBSRMatrix *matrix,
                                     JXF_Int *row_start, JXF_Int *row_end,
                                     JXF_Int *col_start, JXF_Int *col_end);

JXF_Int jxf_ParBSRMatrixSetDataOwner(jxf_ParBSRMatrix *matrix, JXF_Int owns_data);
JXF_Int jxf_ParBSRMatrixSetRowStartsOwner(jxf_ParBSRMatrix *matrix, 
                                         JXF_Int owns_row_starts);
JXF_Int jxf_ParBSRMatrixSetColStartsOwner(jxf_ParBSRMatrix *matrix, 
                                         JXF_Int owns_col_starts);

/* 通信相关 */
JXF_Int jxf_ParBSRMatrixCreateCommPkg(jxf_ParBSRMatrix *matrix);

jxf_ParBSRMatrix* jxf_BSRMatrixToParBSRMatrix(MPI_Comm comm, jxf_BSRMatrix* A, JXF_BigInt* global_row_starts, JXF_BigInt* global_col_starts);
jxf_ParCSRCommHandle *jxf_ParBSRCommHandleCreate( JXF_Int job, JXF_Int bnnz, jxf_ParCSRCommPkg *comm_pkg, void *send_data, void  *recv_data );
JXF_Int jxf_ParBSRMatrixGetRowPartitioning( jxf_ParBSRMatrix *matrix, JXF_Int **row_partitioning_ptr );
JXF_Int jxf_ParBSRMatrixGetColPartitioning( jxf_ParBSRMatrix *matrix, JXF_Int **row_partitioning_ptr );
#endif /* JXF_PARBSR_MV_HEADER */