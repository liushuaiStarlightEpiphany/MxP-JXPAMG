//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_mv.h -- head files for matrix-vector operation
 *  Date: 2011/09/08
 *  Modified Date: 2012/10/22
 *
 *  Created by peghoty
 *  Modified by Yue Xiaoqiang
 */

#ifndef JXF_MV_HEADER
#define JXF_MV_HEADER

#ifndef JXF_UTIL_HEADER
#include "jxf_util.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

typedef struct jxf_SpMVPrecondFP32Data_struct jxf_SpMVPrecondFP32Data;
typedef struct jxh_SpMVPrecondFP16Data_struct jxh_SpMVPrecondFP16Data;

/*!
 * \struct jxf_Vector
 */
typedef struct
{
   JXF_Real *data;
   JXF_Int size;
   JXF_Int owns_data;
   JXF_Int num_vectors;
   JXF_Int multivec_storage_method;
   JXF_Int vecstride, idxstride;

} jxf_Vector;

#define jxf_VectorData(vector) ((vector)->data)
#define jxf_VectorSize(vector) ((vector)->size)
#define jxf_VectorOwnsData(vector) ((vector)->owns_data)
#define jxf_VectorNumVectors(vector) ((vector)->num_vectors)
#define jxf_VectorMultiVecStorageMethod(vector) ((vector)->multivec_storage_method)
#define jxf_VectorVectorStride(vector) ((vector)->vecstride)
#define jxf_VectorIndexStride(vector) ((vector)->idxstride)

/*!
 * \struct jxf_CSRMatrix
 */
typedef struct
{
   JXF_Real *data;
   JXF_Int *i;
   JXF_Int *j;
   JXF_Int num_rows;
   JXF_Int num_cols;
   JXF_Int num_nonzeros;
   JXF_Int *rownnz;
   JXF_Int num_rownnz;
   JXF_Int owns_data;

   jxf_SpMVPrecondFP32Data *spmv_precond_fp32;

} jxf_CSRMatrix;

#define jxf_CSRMatrixData(matrix) ((matrix)->data)
#define jxf_CSRMatrixI(matrix) ((matrix)->i)
#define jxf_CSRMatrixJ(matrix) ((matrix)->j)
#define jxf_CSRMatrixNumRows(matrix) ((matrix)->num_rows)
#define jxf_CSRMatrixNumCols(matrix) ((matrix)->num_cols)
#define jxf_CSRMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxf_CSRMatrixRownnz(matrix) ((matrix)->rownnz)
#define jxf_CSRMatrixNumRownnz(matrix) ((matrix)->num_rownnz)
#define jxf_CSRMatrixOwnsData(matrix) ((matrix)->owns_data)
#define jxf_CSRMatrixSpMVPrecondFP32(matrix) ((matrix)->spmv_precond_fp32)

/*!
 * \struct jxh_CSRMatrix
 * \brief Half-precision CSR matrix with scaling factor.
 */
typedef struct
{
   JXH_Real *data;
   JXF_Int *i;
   JXF_Int *j;
   JXF_Int num_rows;
   JXF_Int num_cols;
   JXF_Int num_nonzeros;
   JXF_Int *rownnz;
   JXF_Int num_rownnz;
   JXF_Int owns_data;

   /*
    * Scaling factor for half-precision storage.
    *
    * F -> H:
    *     data_h[j] = (JXH_Real)(data_f[j] / scale)
    *
    * H -> F:
    *     data_f[j] = (JXF_Real)data_h[j] * scale
    */
   JXF_Real scale;

   /*
    * HmulFacc path:
    *   A:  FP16
    *   x:  FP16 packed
    *   dy: FP16 product buffer
    */
   jxh_SpMVPrecondFP16Data *spmv_precond_fp16;

} jxh_CSRMatrix;
#define jxh_CSRMatrixData(matrix) ((matrix)->data)
#define jxh_CSRMatrixI(matrix) ((matrix)->i)
#define jxh_CSRMatrixJ(matrix) ((matrix)->j)
#define jxh_CSRMatrixNumRows(matrix) ((matrix)->num_rows)
#define jxh_CSRMatrixNumCols(matrix) ((matrix)->num_cols)
#define jxh_CSRMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxh_CSRMatrixRownnz(matrix) ((matrix)->rownnz)
#define jxh_CSRMatrixNumRownnz(matrix) ((matrix)->num_rownnz)
#define jxh_CSRMatrixOwnsData(matrix) ((matrix)->owns_data)
#define jxh_CSRMatrixScale(matrix) ((matrix)->scale)
#define jxh_CSRMatrixSpMVPrecondFP16(matrix) ((matrix)->spmv_precond_fp16)

#define JXF_USING_PERSISTENT_COMM
#ifdef JXF_USING_PERSISTENT_COMM
typedef enum Jxf_CommPkgJobType
{
   JXF_COMM_PKG_JOB_COMPLEX = 0,
   JXF_COMM_PKG_JOB_COMPLEX_TRANSPOSE,
   JXF_COMM_PKG_JOB_INT,
   JXF_COMM_PKG_JOB_INT_TRANSPOSE,
   JXF_NUM_OF_COMM_PKG_JOB_TYPE
} Jxf_CommPkgJobType;

typedef struct
{
   void *send_data;
   void *recv_data;

   JXF_Int num_requests;
   MPI_Request *requests;

   JXF_Int own_send_data, own_recv_data;

} jxf_ParCSRPersistentCommHandle;
#endif

/*!
 * \struct jxf_ParCSRCommPkg
 * \brief Structure containing information for doing communications
 */
typedef struct
{
   MPI_Comm comm;
   JXF_Int num_sends;
   JXF_Int *send_procs;
   JXF_Int *send_map_starts;
   JXF_Int *send_map_elmts;

   JXF_Int num_recvs;
   JXF_Int *recv_procs;
   JXF_Int *recv_vec_starts;

   /* remote communication information */
   MPI_Datatype *send_mpi_types;
   MPI_Datatype *recv_mpi_types;

   // if reoder MatvecCommPkg
#if JXF_REODER_SEND_RECV
   JXF_Int if_reorder_MatvecCommPkg;
   JXF_Int *send_reorder_procs;
   JXF_Int *recv_reorder_procs;
   JXF_Int *send_reorder_map;
   JXF_Int *recv_reorder_map;
#endif
#ifdef JXF_USING_PERSISTENT_COMM
   jxf_ParCSRPersistentCommHandle *persistent_comm_handles[JXF_NUM_OF_COMM_PKG_JOB_TYPE];
#endif

} jxf_ParCSRCommPkg;

#define jxf_ParCSRCommPkgComm(comm_pkg) (comm_pkg->comm)
#define jxf_ParCSRCommPkgNumSends(comm_pkg) (comm_pkg->num_sends)
#define jxf_ParCSRCommPkgSendMapStarts(comm_pkg) (comm_pkg->send_map_starts)
#define jxf_ParCSRCommPkgSendMapElmts(comm_pkg) (comm_pkg->send_map_elmts)
#define jxf_ParCSRCommPkgSendMapElmt(comm_pkg, i) (comm_pkg->send_map_elmts[i])
#define jxf_ParCSRCommPkgNumRecvs(comm_pkg) (comm_pkg->num_recvs)
#define jxf_ParCSRCommPkgRecvVecStarts(comm_pkg) (comm_pkg->recv_vec_starts)
#define jxf_ParCSRCommPkgSendMPITypes(comm_pkg) (comm_pkg->send_mpi_types)
#define jxf_ParCSRCommPkgSendMPIType(comm_pkg, i) (comm_pkg->send_mpi_types[i])
#define jxf_ParCSRCommPkgRecvMPITypes(comm_pkg) (comm_pkg->recv_mpi_types)
#define jxf_ParCSRCommPkgRecvMPIType(comm_pkg, i) (comm_pkg->recv_mpi_types[i])

// add by mrz
#if JXF_REODER_SEND_RECV
#define jxf_hpIfReorderMatvecCommPkg(comm_pkg) (comm_pkg->if_reorder_MatvecCommPkg)
#define jxf_hpSendReorderMap(comm_pkg) (comm_pkg->send_reorder_map)
#define jxf_hpRecvReorderMap(comm_pkg) (comm_pkg->recv_reorder_map)
#define jxf_hpSendReorderMapIndex(comm_pkg, i) (comm_pkg->send_reorder_map[i])
#define jxf_hpRecvReorderMapIndex(comm_pkg, i) (comm_pkg->recv_reorder_map[i])
#define jxf_hpSendReorderProcs(comm_pkg) (comm_pkg->send_reorder_procs)
#define jxf_hpRecvReorderProcs(comm_pkg) (comm_pkg->recv_reorder_procs)
#define jxf_hpSendReorderProc(comm_pkg, i) (comm_pkg->send_reorder_procs[i])
#define jxf_hpRecvReorderProc(comm_pkg, i) (comm_pkg->recv_reorder_procs[i])
#define jxf_hpRecvVecStartReorder(comm_pkg, i, k) (comm_pkg->recv_vec_starts[jxf_hpRecvReorderMapIndex(comm_pkg, i) + k])
#define jxf_hpSendMapStartReorder(comm_pkg, i, k) (comm_pkg->send_map_starts[jxf_hpSendReorderMapIndex(comm_pkg, i) + k])
#endif
#define jxf_ParCSRCommPkgRecvVecStart(comm_pkg, i) (comm_pkg->recv_vec_starts[i])
#define jxf_ParCSRCommPkgSendMapStart(comm_pkg, i) (comm_pkg->send_map_starts[i])
#define jxf_ParCSRCommPkgSendProcs(comm_pkg) (comm_pkg->send_procs)
#define jxf_ParCSRCommPkgSendProc(comm_pkg, i) (comm_pkg->send_procs[i])
#define jxf_ParCSRCommPkgRecvProcs(comm_pkg) (comm_pkg->recv_procs)
#define jxf_ParCSRCommPkgRecvProc(comm_pkg, i) (comm_pkg->recv_procs[i])

/*!
 * \struct jxf_ParCSRCommHandle
 */
typedef struct
{
   jxf_ParCSRCommPkg *comm_pkg;
   void *send_data;
   void *recv_data;

   JXF_Int num_requests;
   MPI_Request *requests;

} jxf_ParCSRCommHandle;

#define jxf_ParCSRCommHandleCommPkg(comm_handle) (comm_handle->comm_pkg)
#define jxf_ParCSRCommHandleSendData(comm_handle) (comm_handle->send_data)
#define jxf_ParCSRCommHandleRecvData(comm_handle) (comm_handle->recv_data)
#define jxf_ParCSRCommHandleNumRequests(comm_handle) (comm_handle->num_requests)
#define jxf_ParCSRCommHandleRequests(comm_handle) (comm_handle->requests)
#define jxf_ParCSRCommHandleRequest(comm_handle, i) (comm_handle->requests[i])

/*!
 * \struct jxf_AuxParCSRMatrix
 */
typedef struct
{
   JXF_Int local_num_rows; /* defines number of rows on this processors */
   JXF_Int local_num_cols; /* defines number of cols of diag */

   JXF_Int need_aux; /* if need_aux = 1, aux_j, aux_data are used to
    generate the parcsr matrix (default),
    for need_aux = 0, data is put directly into
    parcsr structure (requires the knowledge of
    offd_i and diag_i ) */

   JXF_Int *row_length; /* row_length_diag[i] contains number of stored
        elements in i-th row */
   JXF_Int *row_space;  /* row_space_diag[i] contains space allocated to
         i-th row */
   JXF_Int **aux_j;     /* contains collected column indices */
   JXF_Real **aux_data; /* contains collected data */

   JXF_Int *indx_diag;         /* indx_diag[i] points to first empty space of portion
              in diag_j , diag_data assigned to row i */
   JXF_Int *indx_offd;         /* indx_offd[i] points to first empty space of portion
              in offd_j , offd_data assigned to row i */
   JXF_Int max_off_proc_elmts; /* length of off processor stash set for
          SetValues and AddTOValues */
   JXF_Int current_num_elmts;  /* current no. of elements stored in stash */
   JXF_Int off_proc_i_indx;    /* pointer to first empty space in
          set_off_proc_i_set */
   JXF_Int *off_proc_i;        /* length 2*num_off_procs_elmts, contains info pairs
            (code, no. of elmts) where code contains global
            row no. if  SetValues, and (-global row no. -1)
            if  AddToValues*/
   JXF_Int *off_proc_j;        /* contains column indices */
   JXF_Real *off_proc_data;    /* contains corresponding data */
   JXF_Int cancel_indx;        /* number of elements that have to be deleted due
              to setting values from another processor */
} jxf_AuxParCSRMatrix;

#define jxf_AuxParCSRMatrixLocalNumRows(matrix) ((matrix)->local_num_rows)
#define jxf_AuxParCSRMatrixLocalNumCols(matrix) ((matrix)->local_num_cols)
#define jxf_AuxParCSRMatrixNeedAux(matrix) ((matrix)->need_aux)
#define jxf_AuxParCSRMatrixRowLength(matrix) ((matrix)->row_length)
#define jxf_AuxParCSRMatrixRowSpace(matrix) ((matrix)->row_space)
#define jxf_AuxParCSRMatrixAuxJ(matrix) ((matrix)->aux_j)
#define jxf_AuxParCSRMatrixAuxData(matrix) ((matrix)->aux_data)
#define jxf_AuxParCSRMatrixIndxDiag(matrix) ((matrix)->indx_diag)
#define jxf_AuxParCSRMatrixIndxOffd(matrix) ((matrix)->indx_offd)
#define jxf_AuxParCSRMatrixMaxOffProcElmts(matrix) ((matrix)->max_off_proc_elmts)
#define jxf_AuxParCSRMatrixCurrentNumElmts(matrix) ((matrix)->current_num_elmts)
#define jxf_AuxParCSRMatrixOffProcIIndx(matrix) ((matrix)->off_proc_i_indx)
#define jxf_AuxParCSRMatrixOffProcI(matrix) ((matrix)->off_proc_i)
#define jxf_AuxParCSRMatrixOffProcJ(matrix) ((matrix)->off_proc_j)
#define jxf_AuxParCSRMatrixOffProcData(matrix) ((matrix)->off_proc_data)
#define jxf_AuxParCSRMatrixCancelIndx(matrix) ((matrix)->cancel_indx)

/*!
 * \struct jxf_AuxParVector
 */
typedef struct
{
   JXF_Int max_off_proc_elmts; /* length of off processor stash for
                                      SetValues and AddToValues*/
   JXF_Int current_num_elmts;  /* current no. of elements stored in stash */
   JXF_Int *off_proc_i;        /* contains column indices */
   JXF_Real *off_proc_data;    /* contains corresponding data */
   JXF_Int cancel_indx;        /* number of elements that have to be deleted due
                               to setting values from another processor */
} jxf_AuxParVector;

#define jxf_AuxParVectorMaxOffProcElmts(matrix) ((matrix)->max_off_proc_elmts)
#define jxf_AuxParVectorCurrentNumElmts(matrix) ((matrix)->current_num_elmts)
#define jxf_AuxParVectorOffProcI(matrix) ((matrix)->off_proc_i)
#define jxf_AuxParVectorOffProcData(matrix) ((matrix)->off_proc_data)
#define jxf_AuxParVectorCancelIndx(matrix) ((matrix)->cancel_indx)

/*!
 * \struct jxf_IJMatrix
 */
typedef struct jxf_IJMatrix_struct
{
   MPI_Comm comm;

   JXF_Int *row_partitioning; /* distribution of rows across processors */
   JXF_Int *col_partitioning; /* distribution of columns */

   JXF_Int object_type;   /* Indicates the type of "object" */
   void *object;          /* Structure for storing local portion */
   void *translator;      /* optional storage_type specfic structure
                             for holding additional local info */
   void *assumed_part;    /* IJMatrix assumed partition */
   JXF_Int assemble_flag; /* indicates whether matrix has been assembled */

   JXF_Int global_first_row; /* these for data items are necessary */
   JXF_Int global_first_col; /* to be able to avoind using the global */
   JXF_Int global_num_rows;  /* global partition */
   JXF_Int global_num_cols;
   JXF_Int omp_flag;
   JXF_Int print_level;

} jxf_IJMatrix;

#define jxf_IJMatrixComm(matrix) ((matrix)->comm)
#define jxf_IJMatrixRowPartitioning(matrix) ((matrix)->row_partitioning)
#define jxf_IJMatrixColPartitioning(matrix) ((matrix)->col_partitioning)
#define jxf_IJMatrixObjectType(matrix) ((matrix)->object_type)
#define jxf_IJMatrixObject(matrix) ((matrix)->object)
#define jxf_IJMatrixTranslator(matrix) ((matrix)->translator)
#define jxf_IJMatrixAssumedPart(matrix) ((matrix)->assumed_part)
#define jxf_IJMatrixAssembleFlag(matrix) ((matrix)->assemble_flag)
#define jxf_IJMatrixGlobalFirstRow(matrix) ((matrix)->global_first_row)
#define jxf_IJMatrixGlobalFirstCol(matrix) ((matrix)->global_first_col)
#define jxf_IJMatrixGlobalNumRows(matrix) ((matrix)->global_num_rows)
#define jxf_IJMatrixGlobalNumCols(matrix) ((matrix)->global_num_cols)
#define jxf_IJMatrixOMPFlag(matrix) ((matrix)->omp_flag)
#define jxf_IJMatrixPrintLevel(matrix) ((matrix)->print_level)

/*!
 * \struct jxf_IJVector
 */
typedef struct jxf_IJVector_struct
{
   MPI_Comm comm;
   JXF_Int *partitioning;    /* Indicates partitioning over tasks */
   JXF_Int object_type;      /* Indicates the type of "local storage" */
   void *object;             /* Structure for storing local portion */
   void *translator;         /* Structure for storing off processor information */
   void *assumed_part;       /* IJ Vector assumed partition */
   JXF_Int global_first_row; /* these for data items are necessary */
   JXF_Int global_num_rows;  /* to be able to avoind using the global */
                             /* global partition */
   JXF_Int print_level;
} jxf_IJVector;

#define jxf_IJVectorComm(vector) ((vector)->comm)
#define jxf_IJVectorPartitioning(vector) ((vector)->partitioning)
#define jxf_IJVectorObjectType(vector) ((vector)->object_type)
#define jxf_IJVectorObject(vector) ((vector)->object)
#define jxf_IJVectorTranslator(vector) ((vector)->translator)
#define jxf_IJVectorAssumedPart(vector) ((vector)->assumed_part)
#define jxf_IJVectorGlobalFirstRow(vector) ((vector)->global_first_row)
#define jxf_IJVectorGlobalNumRows(vector) ((vector)->global_num_rows)
#define jxf_IJVectorPrintLevel(vector) ((vector)->print_level)

/*!
 * \struct jxf_IJAssumedPart
 */
typedef struct
{
   JXF_Int length;
   JXF_Int row_start;
   JXF_Int row_end;
   JXF_Int storage_length;
   JXF_Int *proc_list;
   JXF_Int *row_start_list;
   JXF_Int *row_end_list;
   JXF_Int *sort_index;

} jxf_IJAssumedPart;

/*!
 * \struct jxf_ParCSRMatrix
 */
typedef struct
{
   MPI_Comm comm;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int first_row_index;
   JXF_Int first_col_diag;
   JXF_Int last_row_index;
   JXF_Int last_col_diag;
   jxf_CSRMatrix *diag;
   jxf_CSRMatrix *offd;
   jxf_CSRMatrix *diagT;
   jxf_CSRMatrix *offdT;
   JXF_Int *col_map_offd;
   JXF_Int *row_starts;
   JXF_Int *col_starts;
   jxf_ParCSRCommPkg *comm_pkg;
   jxf_ParCSRCommPkg *comm_pkgT;
   JXF_Int owns_data;
   JXF_Int owns_row_starts;
   JXF_Int owns_col_starts;
   JXF_Int num_nonzeros;
   JXF_Real d_num_nonzeros;
   JXF_Int *rowindices;
   JXF_Real *rowvalues;
   JXF_Int getrowactive;
   jxf_IJAssumedPart *assumed_partition; /* only populated if no_global_partition option
                                           is used (compile-time option) */
} jxf_ParCSRMatrix;

#define jxf_ParCSRMatrixComm(matrix) ((matrix)->comm)
#define jxf_ParCSRMatrixGlobalNumRows(matrix) ((matrix)->global_num_rows)
#define jxf_ParCSRMatrixGlobalNumCols(matrix) ((matrix)->global_num_cols)
#define jxf_ParCSRMatrixFirstRowIndex(matrix) ((matrix)->first_row_index)
#define jxf_ParCSRMatrixFirstColDiag(matrix) ((matrix)->first_col_diag)
#define jxf_ParCSRMatrixLastRowIndex(matrix) ((matrix)->last_row_index)
#define jxf_ParCSRMatrixLastColDiag(matrix) ((matrix)->last_col_diag)
#define jxf_ParCSRMatrixDiag(matrix) ((matrix)->diag)
#define jxf_ParCSRMatrixOffd(matrix) ((matrix)->offd)
#define jxf_ParCSRMatrixDiagT(matrix) ((matrix)->diagT)
#define jxf_ParCSRMatrixOffdT(matrix) ((matrix)->offdT)
#define jxf_ParCSRMatrixColMapOffd(matrix) ((matrix)->col_map_offd)
#define jxf_ParCSRMatrixRowStarts(matrix) ((matrix)->row_starts)
#define jxf_ParCSRMatrixColStarts(matrix) ((matrix)->col_starts)
#define jxf_ParCSRMatrixCommPkg(matrix) ((matrix)->comm_pkg)
#define jxf_ParCSRMatrixCommPkgT(matrix) ((matrix)->comm_pkgT)
#define jxf_ParCSRMatrixOwnsData(matrix) ((matrix)->owns_data)
#define jxf_ParCSRMatrixOwnsRowStarts(matrix) ((matrix)->owns_row_starts)
#define jxf_ParCSRMatrixOwnsColStarts(matrix) ((matrix)->owns_col_starts)
#define jxf_ParCSRMatrixNumRows(matrix) jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(matrix))
#define jxf_ParCSRMatrixNumCols(matrix) jxf_CSRMatrixNumCols(jxf_ParCSRMatrixDiag(matrix))
#define jxf_ParCSRMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxf_ParCSRMatrixDNumNonzeros(matrix) ((matrix)->d_num_nonzeros)
#define jxf_ParCSRMatrixRowindices(matrix) ((matrix)->rowindices)
#define jxf_ParCSRMatrixRowvalues(matrix) ((matrix)->rowvalues)
#define jxf_ParCSRMatrixGetrowactive(matrix) ((matrix)->getrowactive)
#define jxf_ParCSRMatrixAssumedPartition(matrix) ((matrix)->assumed_partition)

/*!
 * \struct jxh_ParCSRMatrix
 * \brief Half-precision ParCSR matrix.
 */
typedef struct
{
   MPI_Comm comm;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int first_row_index;
   JXF_Int first_col_diag;
   JXF_Int last_row_index;
   JXF_Int last_col_diag;

   jxh_CSRMatrix *diag;
   jxh_CSRMatrix *offd;
   jxh_CSRMatrix *diagT;
   jxh_CSRMatrix *offdT;

   JXF_Int *col_map_offd;
   JXF_Int *row_starts;
   JXF_Int *col_starts;

   /* 通信包通常不依赖浮点精度，可以沿用 jxf_ParCSRCommPkg */
   jxf_ParCSRCommPkg *comm_pkg;
   jxf_ParCSRCommPkg *comm_pkgT;

   JXF_Int owns_data;
   JXF_Int owns_row_starts;
   JXF_Int owns_col_starts;
   JXF_Int num_nonzeros;
   JXF_Real d_num_nonzeros;
   JXF_Int *rowindices;
   JXH_Real *rowvalues;
   JXF_Int getrowactive;
   jxf_IJAssumedPart *assumed_partition;

} jxh_ParCSRMatrix;
#define jxh_ParCSRMatrixComm(matrix) ((matrix)->comm)
#define jxh_ParCSRMatrixGlobalNumRows(matrix) ((matrix)->global_num_rows)
#define jxh_ParCSRMatrixGlobalNumCols(matrix) ((matrix)->global_num_cols)
#define jxh_ParCSRMatrixFirstRowIndex(matrix) ((matrix)->first_row_index)
#define jxh_ParCSRMatrixFirstColDiag(matrix) ((matrix)->first_col_diag)
#define jxh_ParCSRMatrixLastRowIndex(matrix) ((matrix)->last_row_index)
#define jxh_ParCSRMatrixLastColDiag(matrix) ((matrix)->last_col_diag)
#define jxh_ParCSRMatrixDiag(matrix) ((matrix)->diag)
#define jxh_ParCSRMatrixOffd(matrix) ((matrix)->offd)
#define jxh_ParCSRMatrixDiagT(matrix) ((matrix)->diagT)
#define jxh_ParCSRMatrixOffdT(matrix) ((matrix)->offdT)
#define jxh_ParCSRMatrixColMapOffd(matrix) ((matrix)->col_map_offd)
#define jxh_ParCSRMatrixRowStarts(matrix) ((matrix)->row_starts)
#define jxh_ParCSRMatrixColStarts(matrix) ((matrix)->col_starts)
#define jxh_ParCSRMatrixCommPkg(matrix) ((matrix)->comm_pkg)
#define jxh_ParCSRMatrixCommPkgT(matrix) ((matrix)->comm_pkgT)
#define jxh_ParCSRMatrixOwnsData(matrix) ((matrix)->owns_data)
#define jxh_ParCSRMatrixOwnsRowStarts(matrix) ((matrix)->owns_row_starts)
#define jxh_ParCSRMatrixOwnsColStarts(matrix) ((matrix)->owns_col_starts)
#define jxh_ParCSRMatrixNumRows(matrix) jxh_CSRMatrixNumRows(jxh_ParCSRMatrixDiag(matrix))
#define jxh_ParCSRMatrixNumCols(matrix) jxh_CSRMatrixNumCols(jxh_ParCSRMatrixDiag(matrix))
#define jxh_ParCSRMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxh_ParCSRMatrixDNumNonzeros(matrix) ((matrix)->d_num_nonzeros)
#define jxh_ParCSRMatrixRowindices(matrix) ((matrix)->rowindices)
#define jxh_ParCSRMatrixRowvalues(matrix) ((matrix)->rowvalues)
#define jxh_ParCSRMatrixGetrowactive(matrix) ((matrix)->getrowactive)
#define jxh_ParCSRMatrixAssumedPartition(matrix) ((matrix)->assumed_partition)

/*!
 * \struct jxf_ParVector
 */
typedef struct
{
   MPI_Comm comm;

   JXF_Int global_size; // global size of the vector
   JXF_Int first_index; // first index(in global sense) of the local vector
   JXF_Int last_index;  // last  index(in global sense) of the local vector
   JXF_Int *partitioning;
   jxf_Vector *local_vector;
   JXF_Int owns_data;
   JXF_Int owns_partitioning;

   jxf_IJAssumedPart *assumed_partition; /* only populated if no_global_partition option
                                           is used (compile-time option) AND this partition
                                           needed(for setting off-proc elements, for example) */
} jxf_ParVector;

#define jxf_ParVectorComm(vector) ((vector)->comm)
#define jxf_ParVectorGlobalSize(vector) ((vector)->global_size)
#define jxf_ParVectorFirstIndex(vector) ((vector)->first_index)
#define jxf_ParVectorLastIndex(vector) ((vector)->last_index)
#define jxf_ParVectorPartitioning(vector) ((vector)->partitioning)
#define jxf_ParVectorLocalVector(vector) ((vector)->local_vector)
#define jxf_ParVectorOwnsData(vector) ((vector)->owns_data)
#define jxf_ParVectorOwnsPartitioning(vector) ((vector)->owns_partitioning)
#define jxf_ParVectorNumVectors(vector) (jxf_VectorNumVectors(jxf_ParVectorLocalVector(vector)))
#define jxf_ParVectorAssumedPartition(vector) ((vector)->assumed_partition)

/*!
 * \struct jxf_ParVecCommPkg
 * \brief Communication Package for the swith between two
 *        parallel vectors with different partitionings.
 * \author peghoty
 * \date 2011/09/07
 */
typedef struct
{
   MPI_Comm comm;

   JXF_Int num_sends;
   JXF_Int *send_procs;
   JXF_Int *send_starts;

   JXF_Int num_recvs;
   JXF_Int *recv_procs;
   JXF_Int *recv_starts;

} jxf_ParVecCommPkg;

#define jxf_ParVecCommPkgComm(comm_pkg) (comm_pkg->comm)
#define jxf_ParVecCommPkgNumSends(comm_pkg) (comm_pkg->num_sends)
#define jxf_ParVecCommPkgSendProcs(comm_pkg) (comm_pkg->send_procs)
#define jxf_ParVecCommPkgSendProc(comm_pkg, i) (comm_pkg->send_procs[i])
#define jxf_ParVecCommPkgSendStarts(comm_pkg) (comm_pkg->send_starts)
#define jxf_ParVecCommPkgSendStart(comm_pkg, i) (comm_pkg->send_starts[i])
#define jxf_ParVecCommPkgNumRecvs(comm_pkg) (comm_pkg->num_recvs)
#define jxf_ParVecCommPkgRecvProcs(comm_pkg) (comm_pkg->recv_procs)
#define jxf_ParVecCommPkgRecvProc(comm_pkg, i) (comm_pkg->recv_procs[i])
#define jxf_ParVecCommPkgRecvStarts(comm_pkg) (comm_pkg->recv_starts)
#define jxf_ParVecCommPkgRecvStart(comm_pkg, i) (comm_pkg->recv_starts[i])

/*!
 * \struct jxf_ParVecCommHandle
 * \brief Communication handle for the swith between two
 *        parallel vectors with different partitionings.
 * \author peghoty
 * \date 2011/09/07
 */
typedef struct
{
   jxf_ParVecCommPkg *comm_pkg;
   void *send_data;
   void *recv_data;

   JXF_Int num_requests;
   MPI_Request *requests;

} jxf_ParVecCommHandle;

#define jxf_ParVecCommHandleCommPkg(comm_handle) (comm_handle->comm_pkg)
#define jxf_ParVecCommHandleSendData(comm_handle) (comm_handle->send_data)
#define jxf_ParVecCommHandleRecvData(comm_handle) (comm_handle->recv_data)
#define jxf_ParVecCommHandleNumRequests(comm_handle) (comm_handle->num_requests)
#define jxf_ParVecCommHandleRequests(comm_handle) (comm_handle->requests)
#define jxf_ParVecCommHandleRequest(comm_handle, i) (comm_handle->requests[i])

/*!
 * \struct jxf_ProcListElements
 */
typedef struct
{
   JXF_Int length;
   JXF_Int storage_length;
   JXF_Int *id;
   JXF_Int *vec_starts;
   JXF_Int element_storage_length;
   JXF_Int *elements;
   JXF_Real *d_elements; // ?
   void *v_elements;

} jxf_ProcListElements;

/*!
 * \struct jxf_CSRBlockMatrix
 */
typedef struct
{
   JXF_Real *data;
   JXF_Int *i;
   JXF_Int *j;
   JXF_Int block_size;
   JXF_Int num_rows;
   JXF_Int num_cols;
   JXF_Int num_nonzeros;
   JXF_Int owns_data;

} jxf_CSRBlockMatrix;

#define jxf_CSRBlockMatrixData(matrix) ((matrix)->data)
#define jxf_CSRBlockMatrixI(matrix) ((matrix)->i)
#define jxf_CSRBlockMatrixJ(matrix) ((matrix)->j)
#define jxf_CSRBlockMatrixBlockSize(matrix) ((matrix)->block_size)
#define jxf_CSRBlockMatrixNumRows(matrix) ((matrix)->num_rows)
#define jxf_CSRBlockMatrixNumCols(matrix) ((matrix)->num_cols)
#define jxf_CSRBlockMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxf_CSRBlockMatrixOwnsData(matrix) ((matrix)->owns_data)

/*!
 * \struct jxf_ParCSRBlockMatrix
 */
typedef struct
{
   MPI_Comm comm;

   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int first_row_index;
   JXF_Int first_col_diag;
   JXF_Int last_row_index;
   JXF_Int last_col_diag;

   jxf_CSRBlockMatrix *diag;
   jxf_CSRBlockMatrix *offd;

   JXF_Int *col_map_offd;
   JXF_Int *row_starts;
   JXF_Int *col_starts;

   jxf_ParCSRCommPkg *comm_pkg;
   jxf_ParCSRCommPkg *comm_pkgT;

   JXF_Int owns_data;
   JXF_Int owns_row_starts;
   JXF_Int owns_col_starts;
   JXF_Int num_nonzeros;
   JXF_Real d_num_nonzeros;
   JXF_Int *rowindices;
   JXF_Real *rowvalues;
   JXF_Int getrowactive;

   jxf_IJAssumedPart *assumed_partition; /* only populated if no_global_partition option
                                           is used (compile-time option) */
} jxf_ParCSRBlockMatrix;

#define jxf_ParCSRBlockMatrixComm(matrix) ((matrix)->comm)
#define jxf_ParCSRBlockMatrixGlobalNumRows(matrix) ((matrix)->global_num_rows)
#define jxf_ParCSRBlockMatrixGlobalNumCols(matrix) ((matrix)->global_num_cols)
#define jxf_ParCSRBlockMatrixFirstRowIndex(matrix) ((matrix)->first_row_index)
#define jxf_ParCSRBlockMatrixFirstColDiag(matrix) ((matrix)->first_col_diag)
#define jxf_ParCSRBlockMatrixLastRowIndex(matrix) ((matrix)->last_row_index)
#define jxf_ParCSRBlockMatrixLastColDiag(matrix) ((matrix)->last_col_diag)
#define jxf_ParCSRBlockMatrixBlockSize(matrix) ((matrix)->diag->block_size)
#define jxf_ParCSRBlockMatrixDiag(matrix) ((matrix)->diag)
#define jxf_ParCSRBlockMatrixOffd(matrix) ((matrix)->offd)
#define jxf_ParCSRBlockMatrixColMapOffd(matrix) ((matrix)->col_map_offd)
#define jxf_ParCSRBlockMatrixRowStarts(matrix) ((matrix)->row_starts)
#define jxf_ParCSRBlockMatrixColStarts(matrix) ((matrix)->col_starts)
#define jxf_ParCSRBlockMatrixCommPkg(matrix) ((matrix)->comm_pkg)
#define jxf_ParCSRBlockMatrixCommPkgT(matrix) ((matrix)->comm_pkgT)
#define jxf_ParCSRBlockMatrixOwnsData(matrix) ((matrix)->owns_data)
#define jxf_ParCSRBlockMatrixOwnsRowStarts(matrix) ((matrix)->owns_row_starts)
#define jxf_ParCSRBlockMatrixOwnsColStarts(matrix) ((matrix)->owns_col_starts)
#define jxf_ParCSRBlockMatrixNumRows(matrix) jxf_CSRBlockMatrixNumRows(jxf_ParCSRBlockMatrixDiag(matrix))
#define jxf_ParCSRBlockMatrixNumCols(matrix) jxf_CSRBlockMatrixNumCols(jxf_ParCSRBlockMatrixDiag(matrix))
#define jxf_ParCSRBlockMatrixNumNonzeros(matrix) ((matrix)->num_nonzeros)
#define jxf_ParCSRBlockMatrixDNumNonzeros(matrix) ((matrix)->d_num_nonzeros)
#define jxf_ParCSRBlockMatrixRowindices(matrix) ((matrix)->rowindices)
#define jxf_ParCSRBlockMatrixRowvalues(matrix) ((matrix)->rowvalues)
#define jxf_ParCSRBlockMatrixGetrowactive(matrix) ((matrix)->getrowactive)
#define jxf_ParCSRBlockMatrixAssumedPartition(matrix) ((matrix)->assumed_partition)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/mat_vec_data/seq_csr_matrix.c */
jxf_CSRMatrix *jxf_CSRMatrixCreate(JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros);
jxh_CSRMatrix *jxh_CSRMatrixCreate(JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros);
JXF_Int jxf_CSRMatrixInitialize(jxf_CSRMatrix *matrix);
JXF_Int jxh_CSRMatrixInitialize(jxh_CSRMatrix *matrix);
JXF_Int jxf_CSRMatrixDestroy(jxf_CSRMatrix *matrix);
JXF_Int jxh_CSRMatrixDestroy(jxh_CSRMatrix *matrix);
JXF_Int jxf_CSRMatrixPrint(jxf_CSRMatrix *matrix, char *file_name);
jxf_CSRMatrix *jxf_CSRMatrixRead(char *file_name, JXF_Int file_base);
jxf_CSRMatrix *jxf_CSRMatrixRead3(char *file_name);
jxf_CSRMatrix *jxf_CSRMatrixBinaryRead(char *file_name);
jxf_CSRMatrix *jxf_CSRMatrixBinaryRead_FASP(char *file_name);
jxf_CSRMatrix *jxf_CSRMatrixRead2(char *file_name, JXF_Int file_base);
JXF_Int jxf_CSRMatrixCopy(jxf_CSRMatrix *A, jxf_CSRMatrix *B, JXF_Int copy_data);
JXF_Int jxf_CSRMatrixGetBandWidth(jxf_CSRMatrix *A, JXF_Int *nbl_ptr, JXF_Int *nbr_ptr);
JXF_Int jxf_CSRMatrixSetRownnz(jxf_CSRMatrix *matrix);
JXF_Int jxf_CSRMatrixGetLoadBalancedPartitionBoundary(jxf_CSRMatrix *A, JXF_Int idx);
JXF_Int jxf_CSRMatrixGetLoadBalancedPartitionBegin(jxf_CSRMatrix *A);
JXF_Int jxf_CSRMatrixGetLoadBalancedPartitionEnd(jxf_CSRMatrix *A);

/* csrc/mat_vec_data/seq_csr_matop.c */
JXF_Int jxf_CSRMatrixReorder(jxf_CSRMatrix *A);
JXF_Int jxf_CSRMatrixTranspose(jxf_CSRMatrix *A, jxf_CSRMatrix **AT, JXF_Int data);
void jxf_CSRMatrixReorderColumnNumber12(jxf_CSRMatrix *A);
void jxf_CSRMatrixReorderColumnNumberAll(jxf_CSRMatrix *A);
jxf_CSRMatrix *jxf_CSRMatrixAdd(jxf_CSRMatrix *A, jxf_CSRMatrix *B);
jxf_CSRMatrix *jxf_CSRMatrixMultiply(jxf_CSRMatrix *A, jxf_CSRMatrix *B);
void jxf_CSRMatrixScale(jxf_CSRMatrix *A, JXF_Real alpha);
jxf_CSRMatrix *jxf_CSRMatrixSymmetrization(jxf_CSRMatrix *A);
jxf_CSRMatrix *jxf_CSRMatrixMGReorderByVariables(jxf_CSRMatrix *A, JXF_Int num_groups);

/* csrc/mat_vec_data/seq_vector.c */
jxf_Vector *jxf_SeqVectorCreate(JXF_Int size);
JXF_Int jxf_SeqVectorInitialize(jxf_Vector *vector);
JXF_Int jxf_SeqVectorSetDataOwner(jxf_Vector *vector, JXF_Int owns_data);
JXF_Int jxf_SeqVectorDestroy(jxf_Vector *vector);
JXF_Int jxf_SeqVectorCopy(jxf_Vector *x, jxf_Vector *y);
JXF_Int jxf_SeqVectorSetConstantValues(jxf_Vector *x, JXF_Real value);
JXF_Real jxf_SeqVectorInnerProd(jxf_Vector *x, jxf_Vector *y);
JXF_Int jxf_SeqVectorPrint(jxf_Vector *vector, char *file_name);
JXF_Int jxf_SeqVectorAxpy(JXF_Real alpha, jxf_Vector *x, jxf_Vector *y);
JXF_Int jxf_SeqVectorMassDotpTwo(jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Int unroll, JXF_Real *result_x, JXF_Real *result_y);
JXF_Int jxf_SeqVectorMassDotpTwo8(jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Real *result_x, JXF_Real *result_y);
JXF_Int jxf_SeqVectorMassDotpTwo4(jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Real *result_x, JXF_Real *result_y);
JXF_Int jxf_SeqVectorMassAxpy(JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k, JXF_Int unroll);
JXF_Int jxf_SeqVectorMassAxpy8(JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k);
JXF_Int jxf_SeqVectorMassAxpy4(JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k);
JXF_Int jxf_SeqVectorMassInnerProd(jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Int unroll, JXF_Real *result);
JXF_Int jxf_SeqVectorMassInnerProd8(jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Real *result);
JXF_Int jxf_SeqVectorMassInnerProd4(jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Real *result);
jxf_Vector *jxf_SeqVectorMGReorderByVariables(jxf_Vector *x, JXF_Int num_groups);
jxf_Vector *jxf_SeqVectorCloneDeep(jxf_Vector *x);
jxf_Vector *jxf_SeqVectorRead(char *file_name);
jxf_Vector *jxf_SeqVectorBinaryRead(char *file_name);
jxf_Vector *jxf_SeqVectorRead2(char *file_name);
jxf_Vector *jxf_SeqMultiVectorCreate(JXF_Int size, JXF_Int num_vectors);
JXF_Int jxf_SeqVectorScale(JXF_Real alpha, jxf_Vector *y);
JXF_Real jxf_SeqVectorNorm1(jxf_Vector *x);
JXF_Real jxf_SeqVectorPointWiseRelNorm1(jxf_Vector *x, jxf_Vector *y);
JXF_Int jxf_SeqVecZXY(jxf_Vector *z, JXF_Real alpha, jxf_Vector *x, jxf_Vector *y);
JXF_Int jxf_SeqVecMul(jxf_Vector *x, jxf_Vector *y, jxf_Vector *z);

/* csrc/mat_vec_data/seq_csr_matvec.c */

JXF_Int
jxf_CSRMatrixMatvec(JXF_Real alpha,
                    jxf_CSRMatrix *A,
                    jxf_Vector *x,
                    JXF_Real beta,
                    jxf_Vector *y,
                    JXF_Int myid);

JXF_Int
jxf_CSRMatrixMatvec_origin(JXF_Real alpha,
                           jxf_CSRMatrix *A,
                           jxf_Vector *x,
                           JXF_Real beta,
                           jxf_Vector *y);

JXF_Int
jxf_CSRMatrixMatvec_v1(JXF_Real alpha,
                       jxf_CSRMatrix *A,
                       jxf_Vector *x,
                       JXF_Real beta,
                       jxf_Vector *y,
                       JXF_Int myid);

JXF_Int
jxf_CSRMatrixMatvec_v2(JXF_Real alpha,
                       jxf_CSRMatrix *A,
                       jxf_Vector *x,
                       JXF_Real beta,
                       jxf_Vector *y,
                       JXF_Int myid);

JXF_Int
jxf_CSRMatrixMatvec_v3(JXF_Real alpha,
                       jxf_CSRMatrix *A,
                       jxf_Vector *x,
                       JXF_Real beta,
                       jxf_Vector *y,
                       JXF_Int myid);

JXF_Int
jxh_CSRMatrixMatvec_HF32(JXF_Real alpha,
                         jxh_CSRMatrix *A,
                         jxf_Vector *x,
                         JXF_Real beta,
                         jxf_Vector *y,
                         JXF_Int myid);
JXF_Int
jxh_CSRMatrixMatvec_HF32_origin(JXF_Real alpha,
                                jxh_CSRMatrix *A,
                                jxf_Vector *x,
                                JXF_Real beta,
                                jxf_Vector *y);
JXF_Int
jxh_CSRMatrixMatvec_HF32_v1(JXF_Real alpha,
                            jxh_CSRMatrix *A,
                            jxf_Vector *x,
                            JXF_Real beta,
                            jxf_Vector *y,
                            JXF_Int myid);

JXF_Int
jxh_CSRMatrixMatvec_HF32_XF32_v1(JXF_Real alpha,
                                 jxh_CSRMatrix *A,
                                 jxf_Vector *x,
                                 JXF_Real beta,
                                 jxf_Vector *y,
                                 JXF_Int myid);

JXF_Int
jxh_CSRMatrixMatvec_HF32_XF32_VEC_v1(JXF_Real alpha,
                                     jxh_CSRMatrix *A,
                                     jxf_Vector *x,
                                     JXF_Real beta,
                                     jxf_Vector *y,
                                     JXF_Int myid);
JXF_Int
jxh_CSRMatrixMatvec_HmulFacc(JXF_Real alpha,
                             jxh_CSRMatrix *A,
                             jxf_Vector *x,
                             JXF_Real beta,
                             jxf_Vector *y,
                             JXF_Int myid);
JXF_Int
jxh_CSRMatrixMatvec_HmulFacc_origin(JXF_Real alpha,
                                    jxh_CSRMatrix *A,
                                    jxf_Vector *x,
                                    JXF_Real beta,
                                    jxf_Vector *y);
JXF_Int
jxh_CSRMatrixMatvec_HmulFacc_v1(JXF_Real alpha,
                                jxh_CSRMatrix *A,
                                jxf_Vector *x,
                                JXF_Real beta,
                                jxf_Vector *y,
                                JXF_Int myid);

JXF_Int
jxf_CSRMatrixMatvecOutOfPlace(JXF_Real alpha,
                              jxf_CSRMatrix *A,
                              jxf_Vector *x,
                              JXF_Real beta,
                              jxf_Vector *b,
                              jxf_Vector *y,
                              JXF_Int offset);
JXF_Int
jxf_CSRMatrixMatvecT(JXF_Real alpha,
                     jxf_CSRMatrix *A,
                     jxf_Vector *x,
                     JXF_Real beta,
                     jxf_Vector *y);
JXF_Int
jxh_CSRMatrixMatvecT(JXF_Real alpha,
                     jxh_CSRMatrix *A,
                     jxf_Vector *x,
                     JXF_Real beta,
                     jxf_Vector *y);
JXF_Int
jxh_CSRMatrixMatvecT_HmulFacc(JXF_Real alpha,
                              jxh_CSRMatrix *A,
                              jxf_Vector *x,
                              JXF_Real beta,
                              jxf_Vector *y);

/* csrc/mat_vec_data/seq_csr_ls.c */
void jxf_ReadLinearSystemData(char *MatFile,
                              char *RhsFile,
                              JXF_Int **IA_ptr,
                              JXF_Int **JA_ptr,
                              JXF_Real **AA_ptr,
                              JXF_Real **F_ptr,
                              JXF_Int *n_ptr);
void jxf_InitializeLinearSystemStruct(jxf_CSRMatrix **matrix_ptr,
                                      jxf_Vector **vectorF_ptr,
                                      jxf_Vector **vectorX_ptr);
void jxf_FillData2Struct(JXF_Int n,
                         JXF_Int *ia,
                         JXF_Int *ja,
                         JXF_Real *aa,
                         JXF_Real *f,
                         JXF_Real *x,
                         jxf_CSRMatrix *matrix,
                         jxf_Vector *vectorF,
                         jxf_Vector *vectorX);

/* csrc/mat_vec_data/par_csr_matop.c */
JXF_Int jxf_ParCSRMatrixTranspose(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix **AT_ptr, JXF_Int data);
jxf_CSRMatrix *jxf_ParCSRMatrixExtractBExt(jxf_ParCSRMatrix *B, jxf_ParCSRMatrix *A, JXF_Int data);
void jxf_ParCSRMatrixExtractBExt_Arrays(JXF_Int **pB_ext_i,
                                        JXF_Int **pB_ext_j,
                                        JXF_Real **pB_ext_data,
                                        JXF_Int **pB_ext_row_map,
                                        JXF_Int *num_nonzeros,
                                        JXF_Int data,
                                        JXF_Int find_row_map,
                                        MPI_Comm comm,
                                        jxf_ParCSRCommPkg *comm_pkg,
                                        JXF_Int num_cols_B,
                                        JXF_Int num_recvs,
                                        JXF_Int num_sends,
                                        JXF_Int first_col_diag,
                                        JXF_Int first_row_index,
                                        JXF_Int *recv_vec_starts,
                                        JXF_Int *send_map_starts,
                                        JXF_Int *send_map_elmts,
                                        JXF_Int *diag_i,
                                        JXF_Int *diag_j,
                                        JXF_Int *offd_i,
                                        JXF_Int *offd_j,
                                        JXF_Int *col_map_offd,
                                        JXF_Real *diag_data,
                                        JXF_Real *offd_data);
jxf_ParCSRMatrix *
jxf_ParMatmul(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix *B);
jxf_ParCSRMatrix *
jxf_ParTMatmul(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix *B);
void jxf_ParMatmul_RowSizes(JXF_Int **C_diag_i,
                            JXF_Int **C_offd_i,
                            JXF_Int *A_diag_i,
                            JXF_Int *A_diag_j,
                            JXF_Int *A_offd_i,
                            JXF_Int *A_offd_j,
                            JXF_Int *B_diag_i,
                            JXF_Int *B_diag_j,
                            JXF_Int *B_offd_i,
                            JXF_Int *B_offd_j,
                            JXF_Int *B_ext_diag_i,
                            JXF_Int *B_ext_diag_j,
                            JXF_Int *B_ext_offd_i,
                            JXF_Int *B_ext_offd_j,
                            JXF_Int *map_B_to_C,
                            JXF_Int *C_diag_size,
                            JXF_Int *C_offd_size,
                            JXF_Int num_rows_diag_A,
                            JXF_Int num_cols_offd_A,
                            JXF_Int allsquare,
                            JXF_Int num_cols_diag_B,
                            JXF_Int num_cols_offd_B,
                            JXF_Int num_cols_offd_C);
JXF_Int
jxf_ParcsrAdd(JXF_Real alpha,
              jxf_ParCSRMatrix *A,
              JXF_Real beta,
              jxf_ParCSRMatrix *B,
              jxf_ParCSRMatrix **Cout);
JXF_Int jxf_ParCSRMatrixDropSmallEntries(jxf_ParCSRMatrix *A, JXF_Real tol, JXF_Int type);

/* csrc/mat_vec_data/par_csr_matrix.c */
jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile(char *filename, JXF_Int num_functions, JXF_Int file_base);
jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile5(char *filename, JXF_Int num_functions, JXF_Int file_base);
jxf_ParCSRMatrix *
jxf_ParCSRMatrixCreate(MPI_Comm comm,
                       JXF_Int global_num_rows,
                       JXF_Int global_num_cols,
                       JXF_Int *row_starts,
                       JXF_Int *col_starts,
                       JXF_Int num_cols_offd,
                       JXF_Int num_nonzeros_diag,
                       JXF_Int num_nonzeros_offd);
jxh_ParCSRMatrix *
jxh_ParCSRMatrixCreate(MPI_Comm comm,
                       JXF_Int global_num_rows,
                       JXF_Int global_num_cols,
                       JXF_Int *row_starts,
                       JXF_Int *col_starts,
                       JXF_Int num_cols_offd,
                       JXF_Int num_nonzeros_diag,
                       JXF_Int num_nonzeros_offd);
JXF_Int jxf_ParCSRMatrixInitialize(jxf_ParCSRMatrix *matrix);
JXF_Int jxh_ParCSRMatrixInitialize(jxh_ParCSRMatrix *matrix);
jxf_CSRMatrix *jxf_ParCSRMatrixToCSRMatrixAll(jxf_ParCSRMatrix *par_matrix);
jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix(MPI_Comm comm, jxf_CSRMatrix *A, JXF_Int *row_starts, JXF_Int *col_starts);
jxf_ParCSRMatrix *
jxf_CSRMatrixToParCSRMatrix_FromGivenPro(MPI_Comm comm,
                                         JXF_Int srcid,
                                         jxf_CSRMatrix *A,
                                         JXF_Int *row_starts,
                                         JXF_Int *col_starts);
JXF_Int jxf_ParCSRMatrixDestroy_sp(jxf_ParCSRMatrix *matrix);
jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix_sp(MPI_Comm comm, jxf_CSRMatrix *A);
JXF_Int jxf_GenerateDiagAndOffd(jxf_CSRMatrix *A, jxf_ParCSRMatrix *matrix, JXF_Int first_col_diag, JXF_Int last_col_diag);
jxf_CSRMatrix *jxf_MergeDiagAndOffd(jxf_ParCSRMatrix *par_matrix);
jxf_ParCSRMatrix *jxf_ParCSRMatrixRead(MPI_Comm comm, const char *file_name, JXF_Int file_base);
JXF_Int jxf_ParCSRMatrixPrint(jxf_ParCSRMatrix *matrix, const char *file_name);
JXF_Int jxf_ParCSRMatrixCopy(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JXF_Int copy_data);

JXF_Int jxf_ParCSRMatrixDestroy(jxf_ParCSRMatrix *matrix);
JXF_Int
jxf_ParCSRMatrixGetLocalRange(jxf_ParCSRMatrix *matrix,
                              JXF_Int *row_start,
                              JXF_Int *row_end,
                              JXF_Int *col_start,
                              JXF_Int *col_end);
JXF_Int
JXF_ParCSRMatrixGetRow(JXF_ParCSRMatrix matrix,
                       JXF_Int row,
                       JXF_Int *size,
                       JXF_Int **col_ind,
                       JXF_Real **values);
JXF_Int
jxf_ParCSRMatrixGetRow(jxf_ParCSRMatrix *matrix,
                       JXF_Int row,
                       JXF_Int *size,
                       JXF_Int **col_ind,
                       JXF_Real **values);
JXF_Int
JXF_ParCSRMatrixRestoreRow(JXF_ParCSRMatrix matrix,
                           JXF_Int row,
                           JXF_Int *size,
                           JXF_Int **col_ind,
                           JXF_Real **values);
JXF_Int
jxf_ParCSRMatrixRestoreRow(jxf_ParCSRMatrix *matrix,
                           JXF_Int row,
                           JXF_Int *size,
                           JXF_Int **col_ind,
                           JXF_Real **values);
JXF_Int jxf_ParCSRMatrixCreateAssumedPartition(jxf_ParCSRMatrix *matrix);
JXF_Int
jxf_LocateAssummedPartition(MPI_Comm comm,
                            JXF_Int row_start,
                            JXF_Int row_end,
                            JXF_Int global_first_row,
                            JXF_Int global_num_rows,
                            jxf_IJAssumedPart *part,
                            JXF_Int myid);
JXF_Int
jxf_GetAssumedPartitionProcFromRow(MPI_Comm comm,
                                   JXF_Int row,
                                   JXF_Int global_first_row,
                                   JXF_Int global_num_rows,
                                   JXF_Int *proc_id);
JXF_Int
jxf_FillResponseParToCSRMatrix(void *p_recv_contact_buf,
                               JXF_Int contact_size,
                               JXF_Int contact_proc,
                               void *ro,
                               MPI_Comm comm,
                               void **p_send_response_buf,
                               JXF_Int *response_message_size);
JXF_Int jxf_ParCSRMatrixDestroyAssumedPartition(jxf_ParCSRMatrix *matrix);
JXF_Int jxf_ParCSRMatrixGetRowPartitioning(jxf_ParCSRMatrix *matrix, JXF_Int **row_partitioning_ptr);
JXF_Int jxf_ParCSRMatrixGetColPartitioning(jxf_ParCSRMatrix *matrix, JXF_Int **col_partitioning_ptr);
JXF_Int jxf_ParCSRMatrixSetNumNonzeros(jxf_ParCSRMatrix *matrix);
JXF_Int jxh_ParCSRMatrixSetNumNonzeros(jxh_ParCSRMatrix *matrix);
JXF_Int jxf_ParCSRMatrixSetDNumNonzeros(jxf_ParCSRMatrix *matrix);
JXF_Int jxh_ParCSRMatrixSetDNumNonzeros(jxh_ParCSRMatrix *matrix);
JXF_Int jxf_ParCSRMatrixSetDataOwner(jxf_ParCSRMatrix *matrix, JXF_Int owns_data);
JXF_Int jxf_ParCSRMatrixSetRowStartsOwner(jxf_ParCSRMatrix *matrix, JXF_Int owns_row_starts);
JXF_Int jxf_ParCSRMatrixSetColStartsOwner(jxf_ParCSRMatrix *matrix, JXF_Int owns_col_starts);
JXF_Int
jxf_BuildCSRMatrixMPIDataType(JXF_Int num_nonzeros,
                              JXF_Int num_rows,
                              JXF_Real *a_data,
                              JXF_Int *a_i,
                              JXF_Int *a_j,
                              MPI_Datatype *csr_matrix_datatype);
JXF_Int
jxf_Build_Par_LinearSystem(MPI_Comm comm,
                           char *filenamemat,
                           char *filenamerhs,
                           jxf_ParCSRMatrix **par_matrix_ptr,
                           jxf_ParVector **par_rhs_ptr);
JXF_Int jxf_PutDiagFirst(JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n);
JXF_Int jxf_ParCSRMatrixGetDims(jxf_ParCSRMatrix *A, JXF_Int *M, JXF_Int *N);
JXF_Int jxf_ParCSRMatrixGetComm(jxf_ParCSRMatrix *matrix, MPI_Comm *comm);
jxf_CSRMatrix *jxf_MergeDiagAndOffdDropSmall(jxf_ParCSRMatrix *par_matrix, JXF_Real drop_tol);

/* csrc/mat_vec_data/par_vector.c */
jxf_ParVector *jxf_BuildRhsParFromOneFile(char *filename, jxf_ParCSRMatrix *A, JXF_Int file_type);
jxf_ParVector *jxf_BuildRhsParFromOneFile3(char *filename, jxf_ParCSRMatrix *A);
jxf_ParVector *jxf_ParVectorCreate(MPI_Comm comm, JXF_Int global_size, JXF_Int *partitioning);
JXF_Int jxf_ParVectorInitialize(jxf_ParVector *vector);
JXF_Int jxf_ParVectorCopy(jxf_ParVector *x, jxf_ParVector *y);
JXF_Int jxf_ParVectorDestroy(jxf_ParVector *vector);
JXF_Int jxf_ParVectorDestroyAssumedPartition(jxf_ParVector *vector);
JXF_Int jxf_ParVectorSetConstantValues(jxf_ParVector *v, JXF_Real value);
JXF_Real jxf_ParVectorInnerProd(jxf_ParVector *x, jxf_ParVector *y);
JXF_Int jxf_ParVectorMassInnerProd(jxf_ParVector *x, jxf_ParVector **y, JXF_Int k, JXF_Int unroll, JXF_Real *result);
JXF_Real jxf_ParVectorNorm1(jxf_ParVector *x);
JXF_Real jxf_ParVectorNorm2(jxf_ParVector *x);
JXF_Real jxf_ParVectorPointWiseRelNorm1(jxf_ParVector *x, jxf_ParVector *y);
JXF_Int jxf_ParVectorAxpy(JXF_Real alpha, jxf_ParVector *x, jxf_ParVector *y);
JXF_Int jxf_ParVectorMassDotpTwo(jxf_ParVector *x, jxf_ParVector *y, jxf_ParVector **z, JXF_Int k, JXF_Int unroll, JXF_Real *result_x, JXF_Real *result_y);
JXF_Int jxf_ParVectorMassAxpy(JXF_Real *alpha, jxf_ParVector **x, jxf_ParVector *y, JXF_Int k, JXF_Int unroll);
JXF_Int jxf_ParVectorScale(JXF_Real alpha, jxf_ParVector *y);
jxf_Vector *jxf_ParVectorToVectorAll(jxf_ParVector *par_v);
JXF_Int jxf_ParVectorToVector_Alloctaed(jxf_ParVector *par_v, jxf_Vector *vector);
JXF_Int jxf_ParVectorPrint(jxf_ParVector *vector, const char *file_name);
jxf_ParVector *jxf_ParVectorRead(MPI_Comm comm, const char *file_name);
JXF_Int jxf_GeneratePartitioning(JXF_Int length, JXF_Int num_procs, JXF_Int **part_ptr);
JXF_Int jxf_GenerateLocalPartitioning(JXF_Int length, JXF_Int num_procs, JXF_Int myid, JXF_Int **part_ptr);
JXF_Int
jxf_FillResponseParToVectorAll(void *p_recv_contact_buf,
                               JXF_Int contact_size,
                               JXF_Int contact_proc,
                               void *ro,
                               MPI_Comm comm,
                               void **p_send_response_buf,
                               JXF_Int *response_message_size);
jxf_ParVector *jxf_VectorToParVector(MPI_Comm comm, jxf_Vector *v, JXF_Int *vec_starts);
jxf_ParVector *jxf_VectorToParVector_FromGivenPro(MPI_Comm comm, JXF_Int srcid, jxf_Vector *v, JXF_Int *vec_starts);
jxf_ParVector *jxf_VectorToParVector_sp(MPI_Comm comm, jxf_Vector *vector);
JXF_Int
jxf_VectorToParVector_Allocated(MPI_Comm comm,
                                jxf_Vector *v,
                                JXF_Int *vec_starts,
                                jxf_ParVector *par_vector);
JXF_Int
jxf_VectorToParVector_Allocated2(MPI_Comm comm,
                                 jxf_Vector *v,
                                 JXF_Int *vec_starts,
                                 jxf_ParVector *par_vector);
JXF_Int
jxf_VectorToParVector_Allocated_FromGivenPro(MPI_Comm comm,
                                             JXF_Int id,
                                             jxf_Vector *v,
                                             JXF_Int *vec_starts,
                                             jxf_ParVector *par_vector);
jxf_ParVector *jxf_ParMultiVectorCreate(MPI_Comm comm, JXF_Int global_size, JXF_Int *partitioning, JXF_Int num_vectors);
JXF_Int jxf_ParVectorSetDataOwner(jxf_ParVector *vector, JXF_Int owns_data);
JXF_Int jxf_ParVectorSetPartitioningOwner(jxf_ParVector *vector, JXF_Int owns_partitioning);
jxf_ParVecCommPkg *jxf_ParVecCommPkgCreate(MPI_Comm comm, JXF_Int *px, JXF_Int *py);
jxf_ParVecCommHandle *
jxf_ParVecCommHandleCreate(JXF_Int job,
                           jxf_ParVecCommPkg *comm_pkg,
                           void *send_data,
                           void *recv_data);
JXF_Int jxf_ParVecCommHandleDestroy(jxf_ParVecCommHandle *comm_handle);
JXF_Int jxf_ParVecCommPkgDestroy(jxf_ParVecCommPkg *comm_pkg);
JXF_Int
jxf_ParVecSwitch(jxf_ParVecCommPkg *comm_pkg,
                 JXF_Int job,
                 jxf_ParVector *x,
                 jxf_ParVector *y);
JXF_Int jxf_ParVecMul(jxf_ParVector *x, jxf_ParVector *y, jxf_ParVector *z);
JXF_Int jxf_ParVecZXY(jxf_ParVector *z, JXF_Real alpha, jxf_ParVector *x, jxf_ParVector *y);

/* csrc/mat_vec_data/par_csr_matvec.c */
JXF_Int
jxf_ParCSRMatrixMatvec(JXF_Real alpha,
                       jxf_ParCSRMatrix *A,
                       jxf_ParVector *x,
                       JXF_Real beta,
                       jxf_ParVector *y);
JXF_Int
jxh_ParCSRMatrixMatvec(JXF_Real alpha,
                       jxh_ParCSRMatrix *A,
                       jxf_ParVector *x,
                       JXF_Real beta,
                       jxf_ParVector *y);
JXF_Int
jxf_ParCSRMatrixMatvecOutOfPlace(JXF_Real alpha,
                                 jxf_ParCSRMatrix *A,
                                 jxf_ParVector *x,
                                 JXF_Real beta,
                                 jxf_ParVector *b,
                                 jxf_ParVector *y);
JXF_Int
jxf_ParCSRMatrixMatvecT(JXF_Real alpha,
                        jxf_ParCSRMatrix *A,
                        jxf_ParVector *x,
                        JXF_Real beta,
                        jxf_ParVector *y);
JXF_Int
jxh_ParCSRMatrixMatvecT(JXF_Real alpha,
                        jxh_ParCSRMatrix *A,
                        jxf_ParVector *x,
                        JXF_Real beta,
                        jxf_ParVector *y);

/* csrc/mat_vec_data/par_csr_com.c */
JXF_Int jxf_MatvecCommPkgCreate(jxf_ParCSRMatrix *par_A);
JXF_Int jxh_MatvecCommPkgCreate(jxh_ParCSRMatrix *par_A);
jxf_ParCSRCommHandle *
jxf_ParCSRCommHandleCreate(JXF_Int job,
                           jxf_ParCSRCommPkg *comm_pkg,
                           void *send_data,
                           void *recv_data);
JXF_Int jxf_ParCSRCommHandleDestroy(jxf_ParCSRCommHandle *comm_handle);
void jxf_MatvecCommPkgCreate_core(MPI_Comm comm,
                                  JXF_Int *col_map_offd,
                                  JXF_Int first_col_diag,
                                  JXF_Int *col_starts,
                                  JXF_Int num_cols_diag,
                                  JXF_Int num_cols_offd,
                                  JXF_Int firstColDiag,
                                  JXF_Int *colMapOffd,
                                  JXF_Int data,
                                  JXF_Int *p_num_recvs,
                                  JXF_Int **p_recv_procs,
                                  JXF_Int **p_recv_vec_starts,
                                  JXF_Int *p_num_sends,
                                  JXF_Int **p_send_procs,
                                  JXF_Int **p_send_map_starts,
                                  JXF_Int **p_send_map_elmts);
JXF_Int jxf_MatvecCommPkgDestroy(jxf_ParCSRCommPkg *comm_pkg);

/* csrc/mat_vec_data/new_commpkg.c */
JXF_Int
jxf_NewCommPkgCreate_core(MPI_Comm comm,
                          JXF_Int *col_map_off_d,
                          JXF_Int first_col_diag,
                          JXF_Int col_start,
                          JXF_Int col_end,
                          JXF_Int num_cols_off_d,
                          JXF_Int global_num_cols,
                          JXF_Int *p_num_recvs,
                          JXF_Int **p_recv_procs,
                          JXF_Int **p_recv_vec_starts,
                          JXF_Int *p_num_sends,
                          JXF_Int **p_send_procs,
                          JXF_Int **p_send_map_starts,
                          JXF_Int **p_send_map_elements,
                          jxf_IJAssumedPart *apart);
JXF_Int
jxf_RangeFillResponseIJDetermineRecvProcs(void *p_recv_contact_buf,
                                          JXF_Int contact_size,
                                          JXF_Int contact_proc,
                                          void *ro,
                                          MPI_Comm comm,
                                          void **p_send_response_buf,
                                          JXF_Int *response_message_size);
JXF_Int
jxf_FillResponseIJDetermineSendProcs(void *p_recv_contact_buf,
                                     JXF_Int contact_size,
                                     JXF_Int contact_proc,
                                     void *ro,
                                     MPI_Comm comm,
                                     void **p_send_response_buf,
                                     JXF_Int *response_message_size);
JXF_Int
jxf_GetAssumedPartitionRowRange(MPI_Comm comm,
                                JXF_Int proc_id,
                                JXF_Int global_first_row,
                                JXF_Int global_num_rows,
                                JXF_Int *row_start,
                                JXF_Int *row_end);

/* csrc/mat_vec_data/ij_matrix.c */
JXF_Int
JXF_IJMatrixCreate(MPI_Comm comm,
                   JXF_Int ilower,
                   JXF_Int iupper,
                   JXF_Int jlower,
                   JXF_Int jupper,
                   JXF_IJMatrix *matrix);
JXF_Int
jxf_IJMatrixCreateParCSR(jxf_IJMatrix *matrix);
JXF_Int
jxf_AuxParCSRMatrixCreate(jxf_AuxParCSRMatrix **aux_matrix,
                          JXF_Int local_num_rows,
                          JXF_Int local_num_cols,
                          JXF_Int *sizes);
JXF_Int
JXF_IJMatrixInitialize(JXF_IJMatrix matrix);
JXF_Int
jxf_IJMatrixInitializehpCSR(jxf_IJMatrix *matrix);
JXF_Int
jxf_AuxParCSRMatrixInitialize(jxf_AuxParCSRMatrix *matrix);
JXF_Int
JXF_IJMatrixDestroy(JXF_IJMatrix matrix);
JXF_Int
jxf_IJMatrixDestroyhpCSR(jxf_IJMatrix *matrix);
JXF_Int
jxf_AuxParCSRMatrixDestroy(jxf_AuxParCSRMatrix *matrix);
JXF_Int
jxf_AssumedPartitionDestroy(jxf_IJAssumedPart *apart);
JXF_Int
JXF_IJMatrixSetObjectType(JXF_IJMatrix matrix, JXF_Int type);
JXF_Int
JXF_IJMatrixSetValues(JXF_IJMatrix matrix,
                      JXF_Int nrows,
                      JXF_Int *ncols,
                      const JXF_Int *rows,
                      const JXF_Int *cols,
                      const JXF_Real *values);
JXF_Int
jxf_IJMatrixSetValuesOMPhpCSR(jxf_IJMatrix *matrix,
                              JXF_Int nrows,
                              JXF_Int *ncols,
                              const JXF_Int *rows,
                              const JXF_Int *cols,
                              const JXF_Real *values);
JXF_Int
jxf_IJMatrixSetValueshpCSR(jxf_IJMatrix *matrix,
                           JXF_Int nrows,
                           JXF_Int *ncols,
                           const JXF_Int *rows,
                           const JXF_Int *cols,
                           const JXF_Real *values);
JXF_Int
JXF_IJMatrixAssemble(JXF_IJMatrix matrix);
JXF_Int
jxf_IJMatrixAssemblehpCSR(jxf_IJMatrix *matrix);
JXF_Int
jxf_IJMatrixAssembleOffProcValshpCSR(jxf_IJMatrix *matrix,
                                     JXF_Int off_proc_i_indx,
                                     JXF_Int max_off_proc_elmts,
                                     JXF_Int current_num_elmts,
                                     JXF_Int *off_proc_i,
                                     JXF_Int *off_proc_j,
                                     JXF_Real *off_proc_data);
JXF_Int
jxf_IJMatrixCreateAssumedPartition(jxf_IJMatrix *matrix);
JXF_Int
jxf_FillResponseIJOffProcVals(void *p_recv_contact_buf,
                              JXF_Int contact_size,
                              JXF_Int contact_proc,
                              void *ro,
                              MPI_Comm comm,
                              void **p_send_response_buf,
                              JXF_Int *response_message_size);
JXF_Int
jxf_IJMatrixAddToValueshpCSR(jxf_IJMatrix *matrix,
                             JXF_Int nrows,
                             JXF_Int *ncols,
                             const JXF_Int *rows,
                             const JXF_Int *cols,
                             const JXF_Real *values);
JXF_Int
JXF_IJMatrixAddToValues(JXF_IJMatrix matrix,
                        JXF_Int nrows,
                        JXF_Int *ncols,
                        const JXF_Int *rows,
                        const JXF_Int *cols,
                        const JXF_Real *values);
JXF_Int
JXF_IJMatrixRead(const char *filename,
                 MPI_Comm comm,
                 JXF_Int type,
                 JXF_IJMatrix *matrix_ptr);
JXF_Int
JXF_IJMatrixPrint(JXF_IJMatrix matrix, const char *filename);
JXF_Int
jxf_IJMatrixAddToValuesOMPhpCSR(jxf_IJMatrix *matrix,
                                JXF_Int nrows,
                                JXF_Int *ncols,
                                const JXF_Int *rows,
                                const JXF_Int *cols,
                                const JXF_Real *values);
JXF_Int
JXF_IJMatrixGetObject(JXF_IJMatrix matrix, void **object);
JXF_Int
JXF_IJMatrixGetLocalRange(JXF_IJMatrix matrix,
                          JXF_Int *ilower,
                          JXF_Int *iupper,
                          JXF_Int *jlower,
                          JXF_Int *jupper);

/* csrc/mat_vec_data/ij_vector.c */
JXF_Int
JXF_IJVectorCreate(MPI_Comm comm,
                   JXF_Int jlower,
                   JXF_Int jupper,
                   JXF_IJVector *vector);
JXF_Int
jxf_IJVectorCreatePar(jxf_IJVector *vector, JXF_Int *IJpartitioning);
JXF_Int
jxf_AuxParVectorCreate(jxf_AuxParVector **aux_vector);
JXF_Int
JXF_IJVectorInitialize(JXF_IJVector vector);
JXF_Int
jxf_IJVectorInitializePar(jxf_IJVector *vector);
JXF_Int
jxf_AuxParVectorInitialize(jxf_AuxParVector *vector);
JXF_Int
JXF_IJVectorDestroy(JXF_IJVector vector);
JXF_Int
jxf_IJVectorDestroyPar(jxf_IJVector *vector);
JXF_Int
jxf_AuxParVectorDestroy(jxf_AuxParVector *vector);
JXF_Int
JXF_IJVectorSetObjectType(JXF_IJVector vector, JXF_Int type);
JXF_Int
JXF_IJVectorSetValues(JXF_IJVector vector,
                      JXF_Int nvalues,
                      const JXF_Int *indices,
                      const JXF_Real *values);
JXF_Int
jxf_IJVectorSetValuesPar(jxf_IJVector *vector,
                         JXF_Int num_values,
                         const JXF_Int *indices,
                         const JXF_Real *values);
JXF_Int
JXF_IJVectorAssemble(JXF_IJVector vector);
JXF_Int
jxf_IJVectorAssemblePar(jxf_IJVector *vector);
JXF_Int
jxf_IJVectorAssembleOffProcValsPar(jxf_IJVector *vector,
                                   JXF_Int max_off_proc_elmts,
                                   JXF_Int current_num_elmts,
                                   JXF_Int *off_proc_i,
                                   JXF_Real *off_proc_data);
JXF_Int
jxf_IJVectorCreateAssumedPartition(jxf_IJVector *vector);
JXF_Int
jxf_FindProc(JXF_Int *list, JXF_Int value, JXF_Int list_length);
JXF_Int
JXF_IJVectorGetValues(JXF_IJVector vector,
                      JXF_Int nvalues,
                      const JXF_Int *indices,
                      JXF_Real *values);
JXF_Int
jxf_IJVectorGetValuesPar(jxf_IJVector *vector,
                         JXF_Int num_values,
                         const JXF_Int *indices,
                         JXF_Real *values);
JXF_Int
JXF_IJVectorAddToValues(JXF_IJVector vector,
                        JXF_Int nvalues,
                        const JXF_Int *indices,
                        const JXF_Real *values);
JXF_Int
JXF_IJVectorRead(const char *filename,
                 MPI_Comm comm,
                 JXF_Int type,
                 JXF_IJVector *vector_ptr);
JXF_Int
JXF_IJVectorPrint(JXF_IJVector vector, const char *filename);
JXF_Int
jxf_IJVectorAddToValuesPar(jxf_IJVector *vector,
                           JXF_Int num_values,
                           const JXF_Int *indices,
                           const JXF_Real *values);

JXF_Int
JXF_IJVectorGetObject(JXF_IJVector vector,
                      void **object);

/* csrc/operation/matops/par_gs.c */
JXF_Int
jxf_gselim(JXF_Real *A_matrix, JXF_Real *x_vector, JXF_Int size);
JXF_Int
jxh_gselim(JXF_Real *A_matrix, JXF_Real *x_vector, JXF_Int size);
JXF_Int
jxh_gselim_HmulFacc(__fp16 *A_matrix_h, JXF_Real *x_vector, JXF_Int size);

/* csrc/operation/matops/par_gs_piv.c */
JXF_Int
jxf_gselim_piv(JXF_Real *A_matrix, JXF_Real *x_vector, JXF_Int size);

/* csrc/operation/matops/par_rap.c */
jxf_CSRMatrix *
jxf_ExchangeRAPData(jxf_CSRMatrix *RAP_int, jxf_ParCSRCommPkg *comm_pkg_RT);
JXF_Int
jxf_PAMGBuildCoarseOperator(jxf_ParCSRMatrix *RT,
                            jxf_ParCSRMatrix *par_A,
                            jxf_ParCSRMatrix *par_P,
                            jxf_ParCSRMatrix **RAP_ptr);
JXF_Int
jxf_PAMGBuildCoarseOperatorKT(jxf_ParCSRMatrix *RT,
                              jxf_ParCSRMatrix *par_A,
                              jxf_ParCSRMatrix *par_P,
                              JXF_Int keepTranspose,
                              jxf_ParCSRMatrix **RAP_ptr);
JXF_Int
jxf_PAMGBuildCoarseOperatorKT_SplitRAP(jxf_ParCSRMatrix *RT,
                                       jxf_ParCSRMatrix *par_A,
                                       jxf_ParCSRMatrix *par_P,
                                       JXF_Int keepTranspose,
                                       jxf_ParCSRMatrix **RAP_ptr);
JXF_Int
jxf_PAMGBuildCoarseOperatorKT_origin(jxf_ParCSRMatrix *RT,
                                     jxf_ParCSRMatrix *par_A,
                                     jxf_ParCSRMatrix *par_P,
                                     JXF_Int keepTranspose,
                                     jxf_ParCSRMatrix **RAP_ptr);

/* csrc/operation/matops/par_rap_omp.c */
JXF_Int
jxf_PAMGBuildCoarseOperatorOMP(jxf_ParCSRMatrix *RT,
                               jxf_ParCSRMatrix *par_A,
                               jxf_ParCSRMatrix *par_P,
                               jxf_ParCSRMatrix **RAP_ptr,
                               JXF_Int *icor_yoo);

/* csrc/operation/relaxations/par_relax_0.c */
JXF_Int
jxf_PAMGRelax0(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI0(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_1.c */
JXF_Int
jxf_PAMGRelax1(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI1(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_2.c */
JXF_Int
jxf_PAMGRelax2(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI2(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_3.c */
JXF_Int
jxf_PAMGRelax3(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxh_PAMGRelax3(jxh_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxh_PAMGRelax3_HmulFacc(jxh_ParCSRMatrix *par_matrix,
                        jxf_ParVector *par_rhs,
                        JXF_Int *cf_marker,
                        JXF_Int relax_points,
                        JXF_Real relax_weight,
                        JXF_Real omega,
                        jxf_ParVector *par_app,
                        jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI3(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_4.c */
JXF_Int
jxf_PAMGRelax4(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI4(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_5.c */
JXF_Int
jxf_PAMGRelax5(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI5(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_6.c */
JXF_Int
jxf_PAMGRelax6(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI6(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_7.c */
JXF_Int
jxf_PAMGRelax7(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxAI7(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);

/* csrc/operation/relaxations/par_relax_8.c */
JXF_Int
jxf_PAMGRelax8(jxf_ParCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp);

/* csrc/operation/relaxations/par_relax_9.c */
JXF_Int
jxf_PAMGRelax9(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxh_PAMGRelax9(jxh_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxh_CoarseSolveFP16_Switch(void *amg_vdata,
                           jxf_ParVector *par_rhs,
                           jxf_ParVector *par_app);
JXF_Int
jxh_CoarseSolve(void *amg_vdata,
                jxf_ParVector *par_rhs,
                jxf_ParVector *par_app);
JXF_Int
jxh_CoarseSolve_HmulFacc(void *amg_vdata,
                         jxf_ParVector *par_rhs,
                         jxf_ParVector *par_app);

/* csrc/operation/relaxations/par_relax_13.c */
JXF_Int
jxf_PAMGRelax13(jxf_ParCSRMatrix *A,
                jxf_ParVector *f,
                JXF_Int *cf_marker,
                JXF_Int relax_points,
                JXF_Real relax_weight,
                JXF_Real omega,
                JXF_Real *l1_norms,
                jxf_ParVector *u,
                jxf_ParVector *Vtemp,
                jxf_ParVector *Ztemp);

/* csrc/operation/relaxations/par_relax_14.c */
JXF_Int
jxf_PAMGRelax14(jxf_ParCSRMatrix *A,
                jxf_ParVector *f,
                JXF_Int *cf_marker,
                JXF_Int relax_points,
                JXF_Real relax_weight,
                JXF_Real omega,
                JXF_Real *l1_norms,
                jxf_ParVector *u,
                jxf_ParVector *Vtemp,
                jxf_ParVector *Ztemp);

/* csrc/hx/l1norm.c */
JXF_Int jxf_ParCSRComputeL1Norms(jxf_ParCSRMatrix *A, JXF_Int option, JXF_Int *cf_marker, JXF_Real **l1_norm_ptr);
JXF_Int jxf_ParCSRComputeL1NormsThreads(jxf_ParCSRMatrix *A, JXF_Int option, JXF_Int num_threads, JXF_Int *cf_marker, JXF_Real **l1_norm_ptr);

/* 新增 天河字段 */
JXF_Int jxf_CSRMatrixSpMVPrecondFP32Create(jxf_CSRMatrix *A, JXF_Int x_size, JXF_Int myid);
void jxf_CSRMatrixSpMVPrecondFP32Destroy(jxf_CSRMatrix *A);
JXF_Int jxf_ParCSRMatrixSpMVPrecondFP32Create(jxf_ParCSRMatrix *A, JXF_Int myid);

/* 新增 天河字段 */
JXF_Int jxh_CSRMatrixSpMVPrecondFP16Create(jxh_CSRMatrix *A, JXF_Int x_size, JXF_Int myid);
void jxh_CSRMatrixSpMVPrecondFP16Destroy(jxh_CSRMatrix *A);
JXF_Int jxh_ParCSRMatrixSpMVPrecondFP16Create(jxh_ParCSRMatrix *A, JXF_Int myid);

#endif
