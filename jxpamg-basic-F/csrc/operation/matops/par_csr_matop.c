//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_matop.c -- basic operations for parallel matrix.
 *  Date: 2011/09/05
 */ 

#include "jxf_mv.h"

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

/*!
 * \fn JXF_Int jxf_ParCSRMatrixCopy
 * \brief Copies A to B.
 * \note If copy_data = 0, only the structure of A is copied to B.
 *       This routine does not check whether the dimensions of
 *       A and B are compatible.
 * \date 2011/09/08
 */
JXF_Int 
jxf_ParCSRMatrixCopy( jxf_ParCSRMatrix *A, 
                     jxf_ParCSRMatrix *B, 
                     JXF_Int              copy_data )
{
   jxf_CSRMatrix *A_diag;
   jxf_CSRMatrix *A_offd;
   JXF_Int *col_map_offd_A;
   jxf_CSRMatrix *B_diag;
   jxf_CSRMatrix *B_offd;
   JXF_Int *col_map_offd_B;
   JXF_Int num_cols_offd;
   JXF_Int i;

   if (!A)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (!B)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   A_diag = jxf_ParCSRMatrixDiag(A);
   A_offd = jxf_ParCSRMatrixOffd(A);
   col_map_offd_A = jxf_ParCSRMatrixColMapOffd(A);
   B_diag = jxf_ParCSRMatrixDiag(B);
   B_offd = jxf_ParCSRMatrixOffd(B);
   col_map_offd_B = jxf_ParCSRMatrixColMapOffd(B);
   num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   jxf_CSRMatrixCopy(A_diag, B_diag, copy_data);
   jxf_CSRMatrixCopy(A_offd, B_offd, copy_data);
   if (num_cols_offd && col_map_offd_B == NULL)
   {
      col_map_offd_B = jxf_CTAlloc(JXF_Int,num_cols_offd);
      jxf_ParCSRMatrixColMapOffd(B) = col_map_offd_B;
   }
   for (i = 0; i < num_cols_offd; i ++)
   {
      col_map_offd_B[i] = col_map_offd_A[i];
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixTranspose
 * \brief Transpose a parallel CSR matrix.
 * \param *A pointer to the matrix to be transposed.
 * \param **AT_ptr pointer to pointer to the transposed matrix.
 * \param data flag to indicate whether transpose the data part
 *        data = 0: Only transpose the sparse pattern;
 *        data = 1: Transpose both the sparse pattern and the data part.
 * \date 2011/09/05
 */
JXF_Int
jxf_ParCSRMatrixTranspose( jxf_ParCSRMatrix  *A,
                          jxf_ParCSRMatrix **AT_ptr,
                          JXF_Int               data ) 
{
   jxf_ParCSRCommHandle *comm_handle = NULL;
   MPI_Comm             comm = jxf_ParCSRMatrixComm(A);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_CSRMatrix        *A_diag   = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix        *A_offd   = jxf_ParCSRMatrixOffd(A);
   JXF_Int  num_cols = jxf_ParCSRMatrixNumCols(A);
   JXF_Int  first_row_index = jxf_ParCSRMatrixFirstRowIndex(A);
   JXF_Int *row_starts = jxf_ParCSRMatrixRowStarts(A);
   JXF_Int *col_starts = jxf_ParCSRMatrixColStarts(A);

   JXF_Int	      num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int        ierr = 0;
   JXF_Int	      num_sends = 0, num_recvs = 0, num_cols_offd_AT; 
   JXF_Int	      i, j, k, index, counter, j_row;
   JXF_Int        value;

   jxf_ParCSRMatrix *AT;
   jxf_CSRMatrix    *AT_diag;
   jxf_CSRMatrix    *AT_offd;
   jxf_CSRMatrix    *AT_tmp;
   
   JXF_Int first_row_index_AT, first_col_diag_AT;
   JXF_Int local_num_rows_AT, local_num_cols_AT;
   
   JXF_Int    *AT_tmp_i    = NULL;
   JXF_Int    *AT_tmp_j    = NULL;
   JXF_Real *AT_tmp_data = NULL;

   JXF_Int    *AT_buf_i    = NULL;
   JXF_Int    *AT_buf_j    = NULL;
   JXF_Real *AT_buf_data = NULL;

   JXF_Int    *AT_offd_i;
   JXF_Int    *AT_offd_j;
   JXF_Real *AT_offd_data;
   JXF_Int    *col_map_offd_AT;
   JXF_Int    *row_starts_AT;
   JXF_Int    *col_starts_AT;

   JXF_Int num_procs, my_id;

   JXF_Int *recv_procs      = NULL;
   JXF_Int *send_procs      = NULL;
   JXF_Int *recv_vec_starts = NULL;
   JXF_Int *send_map_starts = NULL;
   JXF_Int *send_map_elmts  = NULL;
   JXF_Int *tmp_recv_vec_starts;
   JXF_Int *tmp_send_map_starts;
   
   jxf_ParCSRCommPkg *tmp_comm_pkg;

   jxf_MPI_Comm_size(comm, &num_procs);   
   jxf_MPI_Comm_rank(comm, &my_id);
  
   num_cols_offd_AT = 0;
   counter = 0;
   AT_offd_j = NULL;
   AT_offd_data = NULL;
   col_map_offd_AT = NULL;
 
   /*---------------------------------------------------------------------
    * If there exists no CommPkg for A, a CommPkg is generated using
    * equally load balanced partitionings
    *--------------------------------------------------------------------*/
    
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A); 
   }

   if (num_procs > 1)
   {
      jxf_CSRMatrixTranspose (A_offd, &AT_tmp, data);

      AT_tmp_i = jxf_CSRMatrixI(AT_tmp);
      AT_tmp_j = jxf_CSRMatrixJ(AT_tmp);
      if (data) AT_tmp_data = jxf_CSRMatrixData(AT_tmp);

      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
      recv_procs = jxf_ParCSRCommPkgRecvProcs(comm_pkg);
      send_procs = jxf_ParCSRCommPkgSendProcs(comm_pkg);
      recv_vec_starts = jxf_ParCSRCommPkgRecvVecStarts(comm_pkg);
      send_map_starts = jxf_ParCSRCommPkgSendMapStarts(comm_pkg);
      send_map_elmts = jxf_ParCSRCommPkgSendMapElmts(comm_pkg);

      AT_buf_i = jxf_CTAlloc(JXF_Int,send_map_starts[num_sends]); 

      for (i = 0; i < AT_tmp_i[num_cols_offd]; i ++)
      {
         AT_tmp_j[i] += first_row_index;
      }

      for (i = 0; i < num_cols_offd; i ++)
      {
         AT_tmp_i[i] = AT_tmp_i[i+1] - AT_tmp_i[i];
      }
	
      comm_handle = jxf_ParCSRCommHandleCreate(12, comm_pkg, AT_tmp_i, AT_buf_i);
   }

   jxf_CSRMatrixTranspose(A_diag, &AT_diag, data);

   AT_offd_i = jxf_CTAlloc(JXF_Int, num_cols+1);

   if (num_procs > 1)
   {   
      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;

      tmp_send_map_starts = jxf_CTAlloc(JXF_Int, num_sends+1);
      tmp_recv_vec_starts = jxf_CTAlloc(JXF_Int, num_recvs+1);

      tmp_send_map_starts[0] = send_map_starts[0];
      for (i = 0; i < num_sends; i ++)
      {
	 tmp_send_map_starts[i+1] = tmp_send_map_starts[i];
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
 	 {
	    tmp_send_map_starts[i+1] += AT_buf_i[j];
	    AT_offd_i[send_map_elmts[j]+1] += AT_buf_i[j];
	 }
      }
      for (i = 0; i < num_cols; i ++)
      {
         AT_offd_i[i+1] += AT_offd_i[i];
      }

      tmp_recv_vec_starts[0] = recv_vec_starts[0];
      for (i = 0; i < num_recvs; i ++)
      {
	 tmp_recv_vec_starts[i+1] = tmp_recv_vec_starts[i];
         for (j = recv_vec_starts[i]; j < recv_vec_starts[i+1]; j ++)
         {
            tmp_recv_vec_starts[i+1] += AT_tmp_i[j];
         }
      }

      tmp_comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
      jxf_ParCSRCommPkgComm(tmp_comm_pkg) = comm;
      jxf_ParCSRCommPkgNumSends(tmp_comm_pkg) = num_sends;
      jxf_ParCSRCommPkgNumRecvs(tmp_comm_pkg) = num_recvs;
      jxf_ParCSRCommPkgRecvProcs(tmp_comm_pkg) = recv_procs;
      jxf_ParCSRCommPkgSendProcs(tmp_comm_pkg) = send_procs;
      jxf_ParCSRCommPkgRecvVecStarts(tmp_comm_pkg) = tmp_recv_vec_starts;
      jxf_ParCSRCommPkgSendMapStarts(tmp_comm_pkg) = tmp_send_map_starts;
#if JXF_REODER_SEND_RECV
jxf_hpIfReorderMatvecCommPkg(tmp_comm_pkg) = 0;
  /*-----------------------------------------------------------
   *  reoderMatvecCommPkg
   *----------------------------------------------------------*/
   jxf_hpCSRCommPkgReorder(tmp_comm_pkg);
#endif
      AT_buf_j = jxf_CTAlloc(JXF_Int,tmp_send_map_starts[num_sends]);
      comm_handle = jxf_ParCSRCommHandleCreate(12, tmp_comm_pkg, AT_tmp_j, AT_buf_j);
      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;

      if (data)
      {
         AT_buf_data = jxf_CTAlloc(JXF_Real,tmp_send_map_starts[num_sends]);
         comm_handle = jxf_ParCSRCommHandleCreate(2,tmp_comm_pkg,AT_tmp_data,AT_buf_data);
         jxf_ParCSRCommHandleDestroy(comm_handle);
         comm_handle = NULL;
      }

      jxf_TFree(tmp_recv_vec_starts);
      jxf_TFree(tmp_send_map_starts);
#if JXF_REODER_SEND_RECV
   if(jxf_hpSendReorderMap(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderMap(tmp_comm_pkg));
   }
   if(jxf_hpRecvReorderMap(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderMap(tmp_comm_pkg));
   }
   if(jxf_hpSendReorderProcs(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderProcs(tmp_comm_pkg));
   }
   if(jxf_hpRecvReorderProcs(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderProcs(tmp_comm_pkg));
   }
#endif
      jxf_TFree(tmp_comm_pkg);
      jxf_CSRMatrixDestroy(AT_tmp);

      if (AT_offd_i[num_cols])
      {
         AT_offd_j = jxf_CTAlloc(JXF_Int, AT_offd_i[num_cols]);
         if (data) AT_offd_data = jxf_CTAlloc(JXF_Real, AT_offd_i[num_cols]);
      }
      else
      {
         AT_offd_j = NULL;
         AT_offd_data = NULL;
      }
	 
      counter = 0;
      for (i = 0; i < num_sends; i ++)
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
	 {
	    j_row = send_map_elmts[j];
	    index = AT_offd_i[j_row];
	    for (k = 0; k < AT_buf_i[j]; k ++)
	    {
	       if (data) AT_offd_data[index] = AT_buf_data[counter];
	       AT_offd_j[index++] = AT_buf_j[counter++];
	    }
	    AT_offd_i[j_row] = index;
	 }
      }
      for (i = num_cols; i > 0; i --)
      {
         AT_offd_i[i] = AT_offd_i[i-1];
      }
      AT_offd_i[0] = 0;

      if (counter)
      {
         jxf_qsort0(AT_buf_j, 0, counter - 1);
         num_cols_offd_AT = 1;
	 value = AT_buf_j[0];
         for (i = 1; i < counter; i ++)
	 {
	    if (value < AT_buf_j[i])
	    {
	       AT_buf_j[num_cols_offd_AT++] = AT_buf_j[i];
	       value = AT_buf_j[i];
	    }
	 }
      }

      if (num_cols_offd_AT)
      {
         col_map_offd_AT = jxf_CTAlloc(JXF_Int, num_cols_offd_AT);
      }
      else
      {
         col_map_offd_AT = NULL;
      }
      
      for (i = 0; i < num_cols_offd_AT; i ++)
      {
	 col_map_offd_AT[i] = AT_buf_j[i];
      }

      jxf_TFree(AT_buf_i);
      jxf_TFree(AT_buf_j);
      if (data) jxf_TFree(AT_buf_data);

      for (i = 0; i < counter; i ++)
      {
	 AT_offd_j[i] = jxf_BinarySearch(col_map_offd_AT, AT_offd_j[i], num_cols_offd_AT);
      }
   }

   AT_offd = jxf_CSRMatrixCreate(num_cols,num_cols_offd_AT,counter);
   jxf_CSRMatrixI(AT_offd) = AT_offd_i;
   jxf_CSRMatrixJ(AT_offd) = AT_offd_j;
   jxf_CSRMatrixData(AT_offd) = AT_offd_data;
   

#ifdef JXF_NO_GLOBAL_PARTITION

   row_starts_AT = jxf_CTAlloc(JXF_Int, 2);
   for (i = 0; i < 2; i ++)
   {
      row_starts_AT[i] = col_starts[i];
   }   

   if (row_starts != col_starts)
   {
      col_starts_AT = jxf_CTAlloc(JXF_Int, 2);
      for (i = 0; i < 2; i++)
      {
         col_starts_AT[i] = row_starts[i];
      }
   }
   else
   {
      col_starts_AT = row_starts_AT;
   }

   first_row_index_AT =  row_starts_AT[0];
   first_col_diag_AT  =  col_starts_AT[0];

   local_num_rows_AT = row_starts_AT[1] - first_row_index_AT ;
   local_num_cols_AT = col_starts_AT[1] - first_col_diag_AT;

#else

   row_starts_AT = jxf_CTAlloc(JXF_Int, num_procs+1);
   for (i = 0; i < num_procs + 1; i ++)
   {
      row_starts_AT[i] = col_starts[i];
   }

   if (row_starts != col_starts)
   {
      col_starts_AT = jxf_CTAlloc(JXF_Int, num_procs + 1);
      for (i = 0; i < num_procs+1; i ++)
      {
         col_starts_AT[i] = row_starts[i];
      }
   }
   else
   {
      col_starts_AT = row_starts_AT;
   }
   first_row_index_AT = row_starts_AT[my_id];
   first_col_diag_AT  = col_starts_AT[my_id];

   local_num_rows_AT = row_starts_AT[my_id+1] - first_row_index_AT ;
   local_num_cols_AT = col_starts_AT[my_id+1] - first_col_diag_AT;

#endif


   AT = jxf_CTAlloc(jxf_ParCSRMatrix,1);
   jxf_ParCSRMatrixComm(AT) = comm;
   jxf_ParCSRMatrixDiag(AT) = AT_diag;
   jxf_ParCSRMatrixOffd(AT) = AT_offd;
   jxf_ParCSRMatrixGlobalNumRows(AT) = jxf_ParCSRMatrixGlobalNumCols(A);
   jxf_ParCSRMatrixGlobalNumCols(AT) = jxf_ParCSRMatrixGlobalNumRows(A);
   jxf_ParCSRMatrixRowStarts(AT)  = row_starts_AT;
   jxf_ParCSRMatrixColStarts(AT)  = col_starts_AT;
   jxf_ParCSRMatrixColMapOffd(AT) = col_map_offd_AT;
 
   jxf_ParCSRMatrixFirstRowIndex(AT) = first_row_index_AT;
   jxf_ParCSRMatrixFirstColDiag(AT)  = first_col_diag_AT;

   jxf_ParCSRMatrixLastRowIndex(AT) = first_row_index_AT + local_num_rows_AT - 1;
   jxf_ParCSRMatrixLastColDiag(AT)  = first_col_diag_AT  + local_num_cols_AT - 1;

   jxf_ParCSRMatrixOwnsData(AT) = 1;
   jxf_ParCSRMatrixOwnsRowStarts(AT) = 1;
   jxf_ParCSRMatrixOwnsColStarts(AT) = 1;
   if (row_starts_AT == col_starts_AT)
   {
      jxf_ParCSRMatrixOwnsColStarts(AT) = 0;
   }

   jxf_ParCSRMatrixCommPkg(AT) = NULL;
   jxf_ParCSRMatrixCommPkgT(AT) = NULL;

   jxf_ParCSRMatrixRowindices(AT) = NULL;
   jxf_ParCSRMatrixRowvalues(AT)  = NULL;
   jxf_ParCSRMatrixGetrowactive(AT) = 0;

   *AT_ptr = AT;
  
   return ierr;
}

/*!
 * \fn jxf_CSRMatrix *jxf_ParCSRMatrixExtractBExt
 * \brief Extracts rows from B which are located on other processors
 *        and needed for multiplication with A locally. The rows
 *        are returned as CSRMatrix.
 * \date 2011/09/05
 */ 
jxf_CSRMatrix * 
jxf_ParCSRMatrixExtractBExt( jxf_ParCSRMatrix *B, 
                            jxf_ParCSRMatrix *A, 
                            JXF_Int              data )
{
   MPI_Comm comm = jxf_ParCSRMatrixComm(B);
   JXF_Int first_col_diag = jxf_ParCSRMatrixFirstColDiag(B);
   JXF_Int first_row_index = jxf_ParCSRMatrixFirstRowIndex(B);
   JXF_Int *col_map_offd = jxf_ParCSRMatrixColMapOffd(B);

   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   JXF_Int  num_recvs;
   JXF_Int *recv_vec_starts;
   JXF_Int  num_sends;
   JXF_Int *send_map_starts;
   JXF_Int *send_map_elmts;
 
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(B);

   JXF_Int    *diag_i    = jxf_CSRMatrixI(diag);
   JXF_Int    *diag_j    = jxf_CSRMatrixJ(diag);
   JXF_Real *diag_data = jxf_CSRMatrixData(diag);

   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(B);

   JXF_Int    *offd_i    = jxf_CSRMatrixI(offd);
   JXF_Int    *offd_j    = jxf_CSRMatrixJ(offd);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);

   JXF_Int num_cols_B, num_nonzeros;
   JXF_Int num_rows_B_ext;

   jxf_CSRMatrix *B_ext;

   JXF_Int    *B_ext_i;
   JXF_Int    *B_ext_j;
   JXF_Real *B_ext_data;
   JXF_Int    *idummy;

  /*---------------------------------------------------------------------
   * If there exists no CommPkg for A, a CommPkg is generated using
   * equally load balanced partitionings 
   *--------------------------------------------------------------------*/
   
   if (!jxf_ParCSRMatrixCommPkg(A))
   {
      jxf_MatvecCommPkgCreate(A);
   }
    
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
   recv_vec_starts = jxf_ParCSRCommPkgRecvVecStarts(comm_pkg);
   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   send_map_starts = jxf_ParCSRCommPkgSendMapStarts(comm_pkg);
   send_map_elmts = jxf_ParCSRCommPkgSendMapElmts(comm_pkg);
 
   num_cols_B = jxf_ParCSRMatrixGlobalNumCols(B);
   num_rows_B_ext = recv_vec_starts[num_recvs];

   jxf_ParCSRMatrixExtractBExt_Arrays( &B_ext_i, &B_ext_j, &B_ext_data, &idummy,
                                      &num_nonzeros,
                                      data, 0, comm, comm_pkg,
                                      num_cols_B, num_recvs, num_sends,
                                      first_col_diag, first_row_index,
                                      recv_vec_starts, send_map_starts, send_map_elmts,
                                      diag_i, diag_j, offd_i, offd_j, col_map_offd,
                                      diag_data, offd_data);

   B_ext = jxf_CSRMatrixCreate(num_rows_B_ext,num_cols_B,num_nonzeros);
   jxf_CSRMatrixI(B_ext) = B_ext_i;
   jxf_CSRMatrixJ(B_ext) = B_ext_j;
   if (data) jxf_CSRMatrixData(B_ext) = B_ext_data;

   return B_ext;
}

/*!
 * \fn void jxf_ParCSRMatrixExtractBExt_Arrays
 * \brief This function was formerly part of jxf_ParCSRMatrixExtractBExt
 *        but the code was removed so it can be used for a corresponding function
 *        for Boolean matrices.
 * \date 2011/09/05
 */ 
void 
jxf_ParCSRMatrixExtractBExt_Arrays( JXF_Int              **pB_ext_i, 
                                   JXF_Int              **pB_ext_j, 
                                   JXF_Real           **pB_ext_data, 
                                   JXF_Int              **pB_ext_row_map,
                                   JXF_Int               *num_nonzeros,
                                   JXF_Int                data, 
                                   JXF_Int                find_row_map, 
                                   MPI_Comm           comm, 
                                   jxf_ParCSRCommPkg  *comm_pkg,
                                   JXF_Int                num_cols_B, 
                                   JXF_Int                num_recvs, 
                                   JXF_Int                num_sends,
                                   JXF_Int                first_col_diag, 
                                   JXF_Int                first_row_index,
                                   JXF_Int               *recv_vec_starts, 
                                   JXF_Int               *send_map_starts, 
                                   JXF_Int               *send_map_elmts,
                                   JXF_Int               *diag_i, 
                                   JXF_Int               *diag_j, 
                                   JXF_Int               *offd_i, 
                                   JXF_Int               *offd_j, 
                                   JXF_Int               *col_map_offd,
                                   JXF_Real            *diag_data, 
                                   JXF_Real            *offd_data)
{
   jxf_ParCSRCommHandle *comm_handle;
   jxf_ParCSRCommPkg    *tmp_comm_pkg;
   JXF_Int *B_int_i;
   JXF_Int *B_int_j;
   JXF_Int *B_ext_i;
   JXF_Int *B_ext_j;
   JXF_Real *B_ext_data = NULL;
   JXF_Real *B_int_data = NULL;
   JXF_Int *B_int_row_map = NULL;
   JXF_Int *B_ext_row_map = NULL;
   JXF_Int num_procs, my_id;
   JXF_Int *jdata_recv_vec_starts;
   JXF_Int *jdata_send_map_starts;
 
   JXF_Int i, j, k, counter;
   JXF_Int start_index;
   JXF_Int j_cnt, j_cnt_rm, jrow;
   JXF_Int num_rows_B_ext;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   num_rows_B_ext = recv_vec_starts[num_recvs];
   
   if (num_rows_B_ext < 0)  /* no B_ext, no communication */
   {  
      *pB_ext_i = NULL;
      *pB_ext_j = NULL;
      if ( data ) *pB_ext_data = NULL;
      if ( find_row_map ) *pB_ext_row_map = NULL;
      *num_nonzeros = 0;
      return;
   }
   
   B_int_i = jxf_CTAlloc(JXF_Int, send_map_starts[num_sends]+1);
   B_ext_i = jxf_CTAlloc(JXF_Int, num_rows_B_ext+1);
   *pB_ext_i = B_ext_i;
   
   if (find_row_map) 
   {
      B_int_row_map = jxf_CTAlloc(JXF_Int, send_map_starts[num_sends]+1);
      B_ext_row_map = jxf_CTAlloc(JXF_Int, num_rows_B_ext+1);
      *pB_ext_row_map = B_ext_row_map;
   }

  /*--------------------------------------------------------------------------
   * generate B_int_i through adding number of row-elements of offd and diag
   * for corresponding rows. B_int_i[j+1] contains the number of elements of
   * a row j (which is determined through send_map_elmts) 
   *--------------------------------------------------------------------------*/
   
   B_int_i[0] = 0;
   j_cnt = 0;
   j_cnt_rm = 0;
   *num_nonzeros = 0;
   for (i = 0; i < num_sends; i ++)
   {
      for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
      {
         jrow = send_map_elmts[j];
         B_int_i[++j_cnt] = offd_i[jrow+1] - offd_i[jrow] + diag_i[jrow+1] - diag_i[jrow];
	 *num_nonzeros += B_int_i[j_cnt];
      }
      if (find_row_map) 
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++) 
         {
            jrow = send_map_elmts[j];
            B_int_row_map[j_cnt_rm++] = jrow + first_row_index;
         }
      }
   }

  /*---------------------------------------------------------------------------------
   * initialize communication 
   *--------------------------------------------------------------------------------*/
   comm_handle = jxf_ParCSRCommHandleCreate( 11,comm_pkg, &B_int_i[1],&(B_ext_i[1]) );
   if ( find_row_map ) 
   {
      /* scatter/gather B_int row numbers to form array of B_ext row numbers */
      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, B_int_row_map, B_ext_row_map);
   }

   B_int_j = jxf_CTAlloc(JXF_Int, *num_nonzeros);
   if (data) B_int_data = jxf_CTAlloc(JXF_Real, *num_nonzeros);

   jdata_send_map_starts = jxf_CTAlloc(JXF_Int, num_sends+1);
   jdata_recv_vec_starts = jxf_CTAlloc(JXF_Int, num_recvs+1);
   start_index = B_int_i[0];
   jdata_send_map_starts[0] = start_index;
   counter = 0;
   for (i = 0; i < num_sends; i ++)
   {
	*num_nonzeros = counter;
	for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
	{
	    jrow = send_map_elmts[j];
	    for (k = diag_i[jrow]; k < diag_i[jrow+1]; k ++) 
	    {
		B_int_j[counter] = diag_j[k] + first_col_diag;
		if (data) B_int_data[counter] = diag_data[k];
		counter ++;
  	    }
	    for (k = offd_i[jrow]; k < offd_i[jrow+1]; k ++) 
	    {
		B_int_j[counter] = col_map_offd[offd_j[k]];
		if (data) B_int_data[counter] = offd_data[k];
		counter ++;
  	    }
	   
	}
	*num_nonzeros = counter - *num_nonzeros;

	start_index += *num_nonzeros;
        jdata_send_map_starts[i+1] = start_index;
   }

   tmp_comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
   jxf_ParCSRCommPkgComm(tmp_comm_pkg) = comm;
   jxf_ParCSRCommPkgNumSends(tmp_comm_pkg) = num_sends;
   jxf_ParCSRCommPkgNumRecvs(tmp_comm_pkg) = num_recvs;
   jxf_ParCSRCommPkgSendProcs(tmp_comm_pkg) = jxf_ParCSRCommPkgSendProcs(comm_pkg);
   jxf_ParCSRCommPkgRecvProcs(tmp_comm_pkg) = jxf_ParCSRCommPkgRecvProcs(comm_pkg);
   jxf_ParCSRCommPkgSendMapStarts(tmp_comm_pkg) = jdata_send_map_starts; 
#if JXF_REODER_SEND_RECV
  /*-----------------------------------------------------------
   *  reoderMatvecCommPkg
   *----------------------------------------------------------*/
   jxf_hpIfReorderMatvecCommPkg(tmp_comm_pkg) = 0;
   jxf_hpCSRCommPkgReorder(tmp_comm_pkg);
#endif
   jxf_ParCSRCommHandleDestroy(comm_handle);
   comm_handle = NULL;

  /*--------------------------------------------------------------------------
   * after communication exchange B_ext_i[j+1] contains the number of elements
   * of a row j ! 
   * evaluate B_ext_i and compute *num_nonzeros for B_ext 
   *--------------------------------------------------------------------------*/

   for (i = 0; i < num_recvs; i ++)
   {
      for (j = recv_vec_starts[i]; j < recv_vec_starts[i+1]; j ++)
      {
         B_ext_i[j+1] += B_ext_i[j];
      }
   }
   *num_nonzeros = B_ext_i[num_rows_B_ext];

   *pB_ext_j = jxf_CTAlloc(JXF_Int, *num_nonzeros);
   B_ext_j = *pB_ext_j;
   if (data) 
   {
      *pB_ext_data = jxf_CTAlloc(JXF_Real, *num_nonzeros);
      B_ext_data = *pB_ext_data;
   }

   for (i = 0; i < num_recvs; i ++)
   {
	start_index = B_ext_i[recv_vec_starts[i]];
	*num_nonzeros = B_ext_i[recv_vec_starts[i+1]] - start_index;
	jdata_recv_vec_starts[i+1] = B_ext_i[recv_vec_starts[i+1]];
   }

   jxf_ParCSRCommPkgRecvVecStarts(tmp_comm_pkg) = jdata_recv_vec_starts;

   comm_handle = jxf_ParCSRCommHandleCreate(11, tmp_comm_pkg, B_int_j, B_ext_j);
   jxf_ParCSRCommHandleDestroy(comm_handle);
   comm_handle = NULL;

   if (data)
   {
      comm_handle = jxf_ParCSRCommHandleCreate(1,tmp_comm_pkg,B_int_data,B_ext_data);
      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;
   }

   jxf_TFree(jdata_send_map_starts);
   jxf_TFree(jdata_recv_vec_starts);

#if JXF_REODER_SEND_RECV
   if(jxf_hpSendReorderMap(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderMap(tmp_comm_pkg));
   }
   if(jxf_hpRecvReorderMap(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderMap(tmp_comm_pkg));
   }
   if(jxf_hpSendReorderProcs(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderProcs(tmp_comm_pkg));
   }
   if(jxf_hpRecvReorderProcs(tmp_comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderProcs(tmp_comm_pkg));
   }
#endif
   jxf_TFree(tmp_comm_pkg);
   jxf_TFree(B_int_i);
   jxf_TFree(B_int_j);
   if (data) jxf_TFree(B_int_data);
   if (find_row_map) jxf_TFree(B_int_row_map);
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_ParMatmul
 * \brief multiplies two ParCSRMatrices A and B and returns the product in ParCSRMatrix C
 * \note C does not own the partitionings since its row_starts is owned by A and col_starts by B
 * \date 2017/02/25
 */
jxf_ParCSRMatrix *
jxf_ParMatmul( jxf_ParCSRMatrix *A, jxf_ParCSRMatrix *B )
{
   MPI_Comm      comm = jxf_ParCSRMatrixComm(A);

   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(A);
   
   JXF_Real       *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int          *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int          *A_diag_j = jxf_CSRMatrixJ(A_diag);

   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(A);
   
   JXF_Real       *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int          *A_offd_i = jxf_CSRMatrixI(A_offd);
   JXF_Int          *A_offd_j = jxf_CSRMatrixJ(A_offd);

   JXF_Int          *row_starts_A = jxf_ParCSRMatrixRowStarts(A);
   JXF_Int           num_rows_diag_A = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int           num_cols_diag_A = jxf_CSRMatrixNumCols(A_diag);
   JXF_Int           num_cols_offd_A = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_CSRMatrix *B_diag = jxf_ParCSRMatrixDiag(B);
   
   JXF_Real    *B_diag_data = jxf_CSRMatrixData(B_diag);
   JXF_Int       *B_diag_i = jxf_CSRMatrixI(B_diag);
   JXF_Int       *B_diag_j = jxf_CSRMatrixJ(B_diag);

   jxf_CSRMatrix *B_offd = jxf_ParCSRMatrixOffd(B);
   JXF_Int          *col_map_offd_B = jxf_ParCSRMatrixColMapOffd(B);
   
   JXF_Real    *B_offd_data = jxf_CSRMatrixData(B_offd);
   JXF_Int       *B_offd_i = jxf_CSRMatrixI(B_offd);
   JXF_Int       *B_offd_j = jxf_CSRMatrixJ(B_offd);

   JXF_Int        first_col_diag_B = jxf_ParCSRMatrixFirstColDiag(B);
   JXF_Int        last_col_diag_B;
   JXF_Int       *col_starts_B = jxf_ParCSRMatrixColStarts(B);
   JXF_Int        num_rows_diag_B = jxf_CSRMatrixNumRows(B_diag);
   JXF_Int        num_cols_diag_B = jxf_CSRMatrixNumCols(B_diag);
   JXF_Int        num_cols_offd_B = jxf_CSRMatrixNumCols(B_offd);

   jxf_ParCSRMatrix *C;
   JXF_Int          *col_map_offd_C = NULL;
   JXF_Int          *map_B_to_C = NULL;

   jxf_CSRMatrix *C_diag;

   JXF_Real   *C_diag_data;
   JXF_Int       *C_diag_i;
   JXF_Int       *C_diag_j;

   jxf_CSRMatrix *C_offd;

   JXF_Real   *C_offd_data=NULL;
   JXF_Int       *C_offd_i=NULL;
   JXF_Int       *C_offd_j=NULL;

   JXF_Int        C_diag_size;
   JXF_Int        C_offd_size;
   JXF_Int        num_cols_offd_C = 0;
   
   jxf_CSRMatrix *Bs_ext = NULL;
   
   JXF_Real    *Bs_ext_data = NULL;
   JXF_Int       *Bs_ext_i = NULL;
   JXF_Int       *Bs_ext_j = NULL;

   JXF_Real    *B_ext_diag_data = NULL;
   JXF_Int       *B_ext_diag_i;
   JXF_Int       *B_ext_diag_j = NULL;
   JXF_Int        B_ext_diag_size;

   JXF_Real    *B_ext_offd_data = NULL;
   JXF_Int       *B_ext_offd_i;
   JXF_Int       *B_ext_offd_j = NULL;
   JXF_Int        B_ext_offd_size;

   JXF_Int        n_rows_A, n_cols_A;
   JXF_Int        n_rows_B, n_cols_B;
   JXF_Int        allsquare = 0;
   JXF_Int        num_procs;
   JXF_Int       *my_diag_array;
   JXF_Int       *my_offd_array;
   JXF_Int        max_num_threads;

   JXF_Real    zero = 0.0;

   n_rows_A = jxf_ParCSRMatrixGlobalNumRows(A);
   n_cols_A = jxf_ParCSRMatrixGlobalNumCols(A);
   n_rows_B = jxf_ParCSRMatrixGlobalNumRows(B);
   n_cols_B = jxf_ParCSRMatrixGlobalNumCols(B);

   max_num_threads = jxf_NumThreads();
   my_diag_array = jxf_CTAlloc(JXF_Int, max_num_threads);
   my_offd_array = jxf_CTAlloc(JXF_Int, max_num_threads);

   if (n_cols_A != n_rows_B || num_cols_diag_A != num_rows_diag_B)
   {
      jxf_error_w_msg(JXF_ERROR_GENERIC," Error! Incompatible matrix dimensions!\n");
      return NULL;
   }
   if ( num_rows_diag_A==num_cols_diag_B) allsquare = 1;

   /*-----------------------------------------------------------------------
    *  Extract B_ext, i.e. portion of B that is stored on neighbor procs
    *  and needed locally for matrix matrix product 
    *-----------------------------------------------------------------------*/

   jxf_MPI_Comm_size(comm, &num_procs);

   if (num_procs > 1)
   {
      /*---------------------------------------------------------------------
       * If there exists no CommPkg for A, a CommPkg is generated using
       * equally load balanced partitionings within 
       * jxf_ParCSRMatrixExtractBExt
       *--------------------------------------------------------------------*/
      Bs_ext = jxf_ParCSRMatrixExtractBExt(B,A,1);
      Bs_ext_data = jxf_CSRMatrixData(Bs_ext);
      Bs_ext_i    = jxf_CSRMatrixI(Bs_ext);
      Bs_ext_j    = jxf_CSRMatrixJ(Bs_ext);
   }
   B_ext_diag_i = jxf_CTAlloc(JXF_Int, num_cols_offd_A+1);
   B_ext_offd_i = jxf_CTAlloc(JXF_Int, num_cols_offd_A+1);
   B_ext_diag_size = 0;
   B_ext_offd_size = 0;
   last_col_diag_B = first_col_diag_B + num_cols_diag_B -1;

   JXF_Int *temp = NULL;
#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
   {
     JXF_Int size, rest, ii;
     JXF_Int ns, ne;
     JXF_Int i1, i, j;
     JXF_Int my_offd_size, my_diag_size;
     JXF_Int cnt_offd, cnt_diag;

     JXF_Int num_threads = jxf_NumActiveThreads();

     size = num_cols_offd_A/num_threads;
     rest = num_cols_offd_A - size*num_threads;
     ii = jxf_GetThreadNum();
     if (ii < rest)
     {
       ns = ii*size+ii;
       ne = (ii+1)*size+ii+1;
     }
     else
     {
       ns = ii*size+rest;
       ne = (ii+1)*size+rest;
     }

     my_diag_size = 0;
     my_offd_size = 0;
     for (i=ns; i < ne; i++)
     {
       B_ext_diag_i[i] = my_diag_size;
       B_ext_offd_i[i] = my_offd_size;
       for (j=Bs_ext_i[i]; j < Bs_ext_i[i+1]; j++)
         if (Bs_ext_j[j] < first_col_diag_B || Bs_ext_j[j] > last_col_diag_B)
            my_offd_size++;
         else
            my_diag_size++;
     }
     my_diag_array[ii] = my_diag_size;
     my_offd_array[ii] = my_offd_size;

#if JXF_USING_OPENMP
#pragma omp barrier
#endif

     if (ii)
     {
       my_diag_size = my_diag_array[0];
       my_offd_size = my_offd_array[0];
       for (i1 = 1; i1 < ii; i1++)
       {
          my_diag_size += my_diag_array[i1];
          my_offd_size += my_offd_array[i1];
       }

       for (i1 = ns; i1 < ne; i1++)
       {
          B_ext_diag_i[i1] += my_diag_size;
          B_ext_offd_i[i1] += my_offd_size;
       }
     }
     else
     {
       B_ext_diag_size = 0;
       B_ext_offd_size = 0;
       for (i1 = 0; i1 < num_threads; i1++)
       {
          B_ext_diag_size += my_diag_array[i1];
          B_ext_offd_size += my_offd_array[i1];
       }
       B_ext_diag_i[num_cols_offd_A] = B_ext_diag_size;
       B_ext_offd_i[num_cols_offd_A] = B_ext_offd_size;

       if (B_ext_diag_size)
       {
          B_ext_diag_j = jxf_CTAlloc(JXF_Int, B_ext_diag_size);
          B_ext_diag_data = jxf_CTAlloc(JXF_Real, B_ext_diag_size);
       }
       if (B_ext_offd_size)
       {
          B_ext_offd_j = jxf_CTAlloc(JXF_Int, B_ext_offd_size);
          B_ext_offd_data = jxf_CTAlloc(JXF_Real, B_ext_offd_size);
       }
       if (B_ext_offd_size || num_cols_offd_B)
          temp = jxf_CTAlloc(JXF_Int, B_ext_offd_size+num_cols_offd_B);
     }

#if JXF_USING_OPENMP
#pragma omp barrier
#endif

     cnt_offd = B_ext_offd_i[ns];
     cnt_diag = B_ext_diag_i[ns];
     for (i=ns; i < ne; i++)
     {
       for (j=Bs_ext_i[i]; j < Bs_ext_i[i+1]; j++)
         if (Bs_ext_j[j] < first_col_diag_B || Bs_ext_j[j] > last_col_diag_B)
         {
            temp[cnt_offd] = Bs_ext_j[j];
            B_ext_offd_j[cnt_offd] = Bs_ext_j[j];
            B_ext_offd_data[cnt_offd++] = Bs_ext_data[j];
         }
         else
         {
            B_ext_diag_j[cnt_diag] = Bs_ext_j[j] - first_col_diag_B;
            B_ext_diag_data[cnt_diag++] = Bs_ext_data[j];
         }
     }

#if JXF_USING_OPENMP
#pragma omp barrier
#endif

     if (ii == 0)
     {
      JXF_Int        cnt;

      if (num_procs > 1)
      {
         jxf_CSRMatrixDestroy(Bs_ext);
         Bs_ext = NULL;
      }

      cnt = 0;
      if (B_ext_offd_size || num_cols_offd_B)
      {
         cnt = B_ext_offd_size;
         for (i=0; i < num_cols_offd_B; i++)
            temp[cnt++] = col_map_offd_B[i];
         if (cnt)
         {
            JXF_Int        value;
            jxf_qsort0(temp, 0, cnt-1);
            num_cols_offd_C = 1;
            value = temp[0];
            for (i=1; i < cnt; i++)
            {
               if (temp[i] > value)
               {
                  value = temp[i];
                  temp[num_cols_offd_C++] = value;
               }
            }
         }

         if (num_cols_offd_C)
            col_map_offd_C = jxf_CTAlloc(JXF_Int,num_cols_offd_C);

         for (i=0; i < num_cols_offd_C; i++)
            col_map_offd_C[i] = temp[i];

         jxf_TFree(temp);
      }
     }

#if JXF_USING_OPENMP
#pragma omp barrier
#endif

     for (i=ns; i < ne; i++)
        for (j=B_ext_offd_i[i]; j < B_ext_offd_i[i+1]; j++)
            B_ext_offd_j[j] = jxf_BinarySearch(col_map_offd_C, B_ext_offd_j[j], num_cols_offd_C);

    } /* end parallel region */

    jxf_TFree(my_diag_array);
    jxf_TFree(my_offd_array);

     if (num_cols_offd_B)
     {
         JXF_Int i, cnt;
         map_B_to_C = jxf_CTAlloc(JXF_Int,num_cols_offd_B);

         cnt = 0;
         for (i=0; i < num_cols_offd_C; i++)
            if (col_map_offd_C[i] == col_map_offd_B[cnt])
            {
               map_B_to_C[cnt++] = i;
               if (cnt == num_cols_offd_B) break;
            }
      }

   jxf_ParMatmul_RowSizes(
      /*&C_diag_i, &C_offd_i, &B_marker,*/
      &C_diag_i, &C_offd_i, 
      A_diag_i, A_diag_j, A_offd_i, A_offd_j,
      B_diag_i, B_diag_j, B_offd_i, B_offd_j,
      B_ext_diag_i, B_ext_diag_j, B_ext_offd_i, B_ext_offd_j,
      map_B_to_C,
      &C_diag_size, &C_offd_size,
      num_rows_diag_A, num_cols_offd_A, allsquare,
      num_cols_diag_B, num_cols_offd_B,
      num_cols_offd_C
      );

   /*-----------------------------------------------------------------------
    *  Allocate C_diag_data and C_diag_j arrays.
    *  Allocate C_offd_data and C_offd_j arrays.
    *-----------------------------------------------------------------------*/
 
   last_col_diag_B = first_col_diag_B + num_cols_diag_B - 1;
   C_diag_data = jxf_CTAlloc(JXF_Real, C_diag_size);
   C_diag_j    = jxf_CTAlloc(JXF_Int, C_diag_size);
   if (C_offd_size)
   { 
      C_offd_data = jxf_CTAlloc(JXF_Real, C_offd_size);
      C_offd_j    = jxf_CTAlloc(JXF_Int, C_offd_size);
   } 

   /*-----------------------------------------------------------------------
    *  Second Pass: Fill in C_diag_data and C_diag_j.
    *  Second Pass: Fill in C_offd_data and C_offd_j.
    *-----------------------------------------------------------------------*/

   /*-----------------------------------------------------------------------
    *  Initialize some stuff.
    *-----------------------------------------------------------------------*/
#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
   {
    JXF_Int *B_marker = NULL;
    JXF_Int ns, ne, size, rest, ii;
    JXF_Int i1, i2, i3, jj2, jj3;
    JXF_Int jj_row_begin_diag, jj_count_diag;
    JXF_Int jj_row_begin_offd, jj_count_offd;
    JXF_Int num_threads;
    JXF_Real a_entry; /*, a_b_product;*/

    ii = jxf_GetThreadNum();
    num_threads = jxf_NumActiveThreads();
    size = num_rows_diag_A/num_threads;
    rest = num_rows_diag_A - size*num_threads;
    if (ii < rest)
    {
       ns = ii*size+ii;
       ne = (ii+1)*size+ii+1;
    }
    else
    {
       ns = ii*size+rest;
       ne = (ii+1)*size+rest;
    }
    jj_count_diag = C_diag_i[ns];
    jj_count_offd = C_offd_i[ns];
    if (num_cols_diag_B || num_cols_offd_C)
    B_marker = jxf_CTAlloc(JXF_Int, num_cols_diag_B+num_cols_offd_C);
    for (i1 = 0; i1 < num_cols_diag_B+num_cols_offd_C; i1++)
      B_marker[i1] = -1;

    /*-----------------------------------------------------------------------
     *  Loop over interior c-points.
     *-----------------------------------------------------------------------*/

    for (i1 = ns; i1 < ne; i1++)
    {

      /*--------------------------------------------------------------------
       *  Create diagonal entry, C_{i1,i1} 
       *--------------------------------------------------------------------*/

      jj_row_begin_diag = jj_count_diag;
      jj_row_begin_offd = jj_count_offd;
      if ( allsquare ) 
      {
         B_marker[i1] = jj_count_diag;
         C_diag_data[jj_count_diag] = zero;
         C_diag_j[jj_count_diag] = i1;
         jj_count_diag++;
      }

      /*-----------------------------------------------------------------
       *  Loop over entries in row i1 of A_offd.
       *-----------------------------------------------------------------*/
         
      if (num_cols_offd_A)
      {
         for (jj2 = A_offd_i[i1]; jj2 < A_offd_i[i1+1]; jj2++)
         {
            i2 = A_offd_j[jj2];
            a_entry = A_offd_data[jj2];
            
            /*-----------------------------------------------------------
             *  Loop over entries in row i2 of B_ext.
             *-----------------------------------------------------------*/

            for (jj3 = B_ext_offd_i[i2]; jj3 < B_ext_offd_i[i2+1]; jj3++)
            {
               i3 = num_cols_diag_B+B_ext_offd_j[jj3];
                  
               /*--------------------------------------------------------
                *  Check B_marker to see that C_{i1,i3} has not already
                *  been accounted for. If it has not, create a new entry.
                *  If it has, add new contribution.
                *--------------------------------------------------------*/

               if (B_marker[i3] < jj_row_begin_offd)
               {
                  B_marker[i3] = jj_count_offd;
                  C_offd_data[jj_count_offd] = a_entry*B_ext_offd_data[jj3];
                  C_offd_j[jj_count_offd] = i3-num_cols_diag_B;
                  jj_count_offd++;
               }
               else
                  C_offd_data[B_marker[i3]] += a_entry*B_ext_offd_data[jj3];
            }
            for (jj3 = B_ext_diag_i[i2]; jj3 < B_ext_diag_i[i2+1]; jj3++)
            {
               i3 = B_ext_diag_j[jj3];
               if (B_marker[i3] < jj_row_begin_diag)
               {
                  B_marker[i3] = jj_count_diag;
                  C_diag_data[jj_count_diag] = a_entry*B_ext_diag_data[jj3];
                  C_diag_j[jj_count_diag] = i3;
                  jj_count_diag++;
               }
               else
                  C_diag_data[B_marker[i3]] += a_entry*B_ext_diag_data[jj3];
            }
         }
      }

      /*-----------------------------------------------------------------
       *  Loop over entries in row i1 of A_diag.
       *-----------------------------------------------------------------*/

      for (jj2 = A_diag_i[i1]; jj2 < A_diag_i[i1+1]; jj2++)
      {
         i2 = A_diag_j[jj2];
         a_entry = A_diag_data[jj2];
            
         /*-----------------------------------------------------------
          *  Loop over entries in row i2 of B_diag.
          *-----------------------------------------------------------*/

         for (jj3 = B_diag_i[i2]; jj3 < B_diag_i[i2+1]; jj3++)
         {
            i3 = B_diag_j[jj3];
                  
            /*--------------------------------------------------------
             *  Check B_marker to see that C_{i1,i3} has not already
             *  been accounted for. If it has not, create a new entry.
             *  If it has, add new contribution.
             *--------------------------------------------------------*/

            if (B_marker[i3] < jj_row_begin_diag)
            {
               B_marker[i3] = jj_count_diag;
               C_diag_data[jj_count_diag] = a_entry*B_diag_data[jj3];
               C_diag_j[jj_count_diag] = i3;
               jj_count_diag++;
            }
            else
            {
               C_diag_data[B_marker[i3]] += a_entry*B_diag_data[jj3];
            }
         }
         if (num_cols_offd_B)
         {
            for (jj3 = B_offd_i[i2]; jj3 < B_offd_i[i2+1]; jj3++)
            {
               i3 = num_cols_diag_B+map_B_to_C[B_offd_j[jj3]];
                  
               /*--------------------------------------------------------
                *  Check B_marker to see that C_{i1,i3} has not already
                *  been accounted for. If it has not, create a new entry.
                *  If it has, add new contribution.
                *--------------------------------------------------------*/

               if (B_marker[i3] < jj_row_begin_offd)
               {
                  B_marker[i3] = jj_count_offd;
                  C_offd_data[jj_count_offd] = a_entry*B_offd_data[jj3];
                  C_offd_j[jj_count_offd] = i3-num_cols_diag_B;
                  jj_count_offd++;
               }
               else
               {
                  C_offd_data[B_marker[i3]] += a_entry*B_offd_data[jj3];
               }
            }
         }
      }
    }
    jxf_TFree(B_marker);
   } /*end parallel region */

   C = jxf_ParCSRMatrixCreate(comm, n_rows_A, n_cols_B, row_starts_A,
                 col_starts_B, num_cols_offd_C, C_diag_size, C_offd_size);

   /* Note that C does not own the partitionings */
   jxf_ParCSRMatrixSetRowStartsOwner(C,0);
   jxf_ParCSRMatrixSetColStartsOwner(C,0);

   C_diag = jxf_ParCSRMatrixDiag(C);
   jxf_CSRMatrixData(C_diag) = C_diag_data; 
   jxf_CSRMatrixI(C_diag) = C_diag_i; 
   jxf_CSRMatrixJ(C_diag) = C_diag_j; 

   C_offd = jxf_ParCSRMatrixOffd(C);
   jxf_CSRMatrixI(C_offd) = C_offd_i; 
   jxf_ParCSRMatrixOffd(C) = C_offd;

   if (num_cols_offd_C)
   {
      jxf_CSRMatrixData(C_offd) = C_offd_data; 
      jxf_CSRMatrixJ(C_offd) = C_offd_j; 
      jxf_ParCSRMatrixColMapOffd(C) = col_map_offd_C;

   }

   /*-----------------------------------------------------------------------
    *  Free various arrays
    *-----------------------------------------------------------------------*/

   jxf_TFree(B_ext_diag_i);
   if (B_ext_diag_size)
   {
      jxf_TFree(B_ext_diag_j);
      jxf_TFree(B_ext_diag_data);
   }
   jxf_TFree(B_ext_offd_i);
   if (B_ext_offd_size)
   {
      jxf_TFree(B_ext_offd_j);
      jxf_TFree(B_ext_offd_data);
   }
   if (num_cols_offd_B) jxf_TFree(map_B_to_C);

   return C;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_ParTMatmul
 * \brief multiplies two ParCSRMatrices transpose(A) and B and returns the product in ParCSRMatrix C
 * \note C does not own the partitionings since its row_starts is owned by A and col_starts by B.
 * \date 2017/02/25
 */
jxf_ParCSRMatrix *
jxf_ParTMatmul( jxf_ParCSRMatrix *A, jxf_ParCSRMatrix *B )
{
   MPI_Comm 	     comm = jxf_ParCSRMatrixComm(A);
   jxf_ParCSRCommPkg *comm_pkg_A = jxf_ParCSRMatrixCommPkg(A);

   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix *AT_diag = NULL;

   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(A);
   jxf_CSRMatrix *AT_offd = NULL;
   
   JXF_Int	num_rows_diag_A = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int	num_cols_diag_A = jxf_CSRMatrixNumCols(A_diag);
   
   jxf_CSRMatrix *B_diag = jxf_ParCSRMatrixDiag(B);
   
   jxf_CSRMatrix *B_offd = jxf_ParCSRMatrixOffd(B);
   JXF_Int	        *col_map_offd_B = jxf_ParCSRMatrixColMapOffd(B);
   
   JXF_Int	first_col_diag_B = jxf_ParCSRMatrixFirstColDiag(B);
   JXF_Int *col_starts_A = jxf_ParCSRMatrixColStarts(A);
   JXF_Int *col_starts_B = jxf_ParCSRMatrixColStarts(B);
   JXF_Int	num_rows_diag_B = jxf_CSRMatrixNumRows(B_diag);
   JXF_Int	num_cols_diag_B = jxf_CSRMatrixNumCols(B_diag);
   JXF_Int	num_cols_offd_B = jxf_CSRMatrixNumCols(B_offd);

   jxf_ParCSRMatrix *C;
   JXF_Int	           *col_map_offd_C = NULL;
   JXF_Int	           *map_B_to_C = NULL;

   jxf_CSRMatrix *C_diag = NULL;
   jxf_CSRMatrix *C_tmp_diag = NULL;

   JXF_Real    *C_diag_data = NULL;
   JXF_Int       *C_diag_i = NULL;
   JXF_Int       *C_diag_j = NULL;
   JXF_Int	first_col_diag_C;
   JXF_Int	last_col_diag_C;

   jxf_CSRMatrix *C_offd = NULL;
   jxf_CSRMatrix *C_tmp_offd = NULL;
   jxf_CSRMatrix *C_int = NULL;
   jxf_CSRMatrix *C_ext = NULL;
   JXF_Int   *C_ext_i = NULL;
   JXF_Int   *C_ext_j = NULL;
   JXF_Real   *C_ext_data = NULL;
   JXF_Int   *C_ext_diag_i = NULL;
   JXF_Int   *C_ext_diag_j = NULL;
   JXF_Real   *C_ext_diag_data = NULL;
   JXF_Int   *C_ext_offd_i = NULL;
   JXF_Int   *C_ext_offd_j = NULL;
   JXF_Real   *C_ext_offd_data = NULL;
   JXF_Int	C_ext_size = 0;
   JXF_Int	C_ext_diag_size = 0;
   JXF_Int	C_ext_offd_size = 0;

   JXF_Int   *C_tmp_diag_i;
   JXF_Int   *C_tmp_diag_j = NULL;
   JXF_Real   *C_tmp_diag_data = NULL;
   JXF_Int   *C_tmp_offd_i = NULL;
   JXF_Int   *C_tmp_offd_j = NULL;
   JXF_Real   *C_tmp_offd_data = NULL;

   JXF_Real          *C_offd_data=NULL;
   JXF_Int       *C_offd_i=NULL;
   JXF_Int       *C_offd_j=NULL;

   JXF_Int       *temp;
   JXF_Int       *send_map_starts_A = NULL;
   JXF_Int       *send_map_elmts_A = NULL;
   JXF_Int        num_sends_A = 0;

   JXF_Int		    num_cols_offd_C = 0;
   
   JXF_Int		   *P_marker;

   JXF_Int              i, j;
   JXF_Int              i1, j_indx;
   
   JXF_Int		    n_rows_A, n_cols_A;
   JXF_Int		    n_rows_B, n_cols_B;
   /*JXF_Int              allsquare = 0;*/
   JXF_Int              cnt, cnt_offd, cnt_diag;
   JXF_Int 		    value;
   JXF_Int 		    num_procs, my_id;
   JXF_Int                max_num_threads;
   JXF_Int               *C_diag_array = NULL;
   JXF_Int               *C_offd_array = NULL;

   JXF_Int first_row_index, first_col_diag;
   JXF_Int local_num_rows, local_num_cols;

   n_rows_A = jxf_ParCSRMatrixGlobalNumRows(A);
   n_cols_A = jxf_ParCSRMatrixGlobalNumCols(A);
   n_rows_B = jxf_ParCSRMatrixGlobalNumRows(B);
   n_cols_B = jxf_ParCSRMatrixGlobalNumCols(B);

   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   max_num_threads = jxf_NumThreads();

   if (n_rows_A != n_rows_B || num_rows_diag_A != num_rows_diag_B)
   {
        jxf_error_w_msg(JXF_ERROR_GENERIC," Error! Incompatible matrix dimensions!\n");
	return NULL;
   }

   /*if (num_cols_diag_A == num_cols_diag_B) allsquare = 1;*/

   jxf_CSRMatrixTranspose(A_diag, &AT_diag, 1);
   jxf_CSRMatrixTranspose(A_offd, &AT_offd, 1);

   C_tmp_diag = jxf_CSRMatrixMultiply(AT_diag, B_diag);
   C_ext_size = 0;
   if (num_procs > 1) 
   {
      jxf_CSRMatrix *C_int_diag;
      jxf_CSRMatrix *C_int_offd;
      C_tmp_offd = jxf_CSRMatrixMultiply(AT_diag, B_offd);
      C_int_diag = jxf_CSRMatrixMultiply(AT_offd, B_diag);
      C_int_offd = jxf_CSRMatrixMultiply(AT_offd, B_offd);
      jxf_ParCSRMatrixDiag(B) = C_int_diag;
      jxf_ParCSRMatrixOffd(B) = C_int_offd;
      C_int = jxf_MergeDiagAndOffd(B);
      jxf_ParCSRMatrixDiag(B) = B_diag;
      jxf_ParCSRMatrixOffd(B) = B_offd;
      C_ext = jxf_ExchangeRAPData(C_int, comm_pkg_A);
      C_ext_i = jxf_CSRMatrixI(C_ext);
      C_ext_j = jxf_CSRMatrixJ(C_ext);
      C_ext_data = jxf_CSRMatrixData(C_ext);
      C_ext_size = C_ext_i[jxf_CSRMatrixNumRows(C_ext)];
      jxf_CSRMatrixDestroy(C_int);
      jxf_CSRMatrixDestroy(C_int_diag);
      jxf_CSRMatrixDestroy(C_int_offd);
   }
   else
   {
     C_tmp_offd = jxf_CSRMatrixCreate(num_cols_diag_A, 0, 0);
     jxf_CSRMatrixInitialize(C_tmp_offd);
   }
   jxf_CSRMatrixDestroy(AT_diag);
   jxf_CSRMatrixDestroy(AT_offd);

   /*-----------------------------------------------------------------------
    *  Add contents of C_ext to C_tmp_diag and C_tmp_offd
    *  to obtain C_diag and C_offd
    *-----------------------------------------------------------------------*/

   /* check for new nonzero columns in C_offd generated through C_ext */

   first_col_diag_C = first_col_diag_B;
   last_col_diag_C = first_col_diag_B + num_cols_diag_B - 1;

   C_tmp_diag_i = jxf_CSRMatrixI(C_tmp_diag);
   if (C_ext_size || num_cols_offd_B)
   {
      JXF_Int C_ext_num_rows;
      num_sends_A = jxf_ParCSRCommPkgNumSends(comm_pkg_A);
      send_map_starts_A = jxf_ParCSRCommPkgSendMapStarts(comm_pkg_A);
      send_map_elmts_A = jxf_ParCSRCommPkgSendMapElmts(comm_pkg_A);
      C_ext_num_rows =  send_map_starts_A[num_sends_A];
     
      C_ext_diag_i = jxf_CTAlloc(JXF_Int, C_ext_num_rows+1);
      C_ext_offd_i = jxf_CTAlloc(JXF_Int, C_ext_num_rows+1);
      temp = jxf_CTAlloc(JXF_Int, C_ext_size+num_cols_offd_B);
      C_ext_diag_size = 0;
      C_ext_offd_size = 0;
      for (i=0; i < C_ext_num_rows; i++)
      {
         for (j=C_ext_i[i]; j < C_ext_i[i+1]; j++)
            if (C_ext_j[j] < first_col_diag_C || C_ext_j[j] > last_col_diag_C)
	       temp[C_ext_offd_size++] = C_ext_j[j];
            else
               C_ext_diag_size++;
         C_ext_diag_i[i+1] = C_ext_diag_size;
         C_ext_offd_i[i+1] = C_ext_offd_size;
      }
      cnt = C_ext_offd_size;
      for (i=0; i < num_cols_offd_B; i++)
         temp[cnt++] = col_map_offd_B[i];

      if (cnt)
      {
	  jxf_qsort0(temp,0,cnt-1);
          value = temp[0];
          num_cols_offd_C = 1;
          for (i=1; i < cnt; i++)
          {
 	     if (temp[i] > value)
	     {
		value = temp[i];
		temp[num_cols_offd_C++] = value;
	     }
  	  }
       }

       if (num_cols_offd_C)
	  col_map_offd_C = jxf_CTAlloc(JXF_Int, num_cols_offd_C);
       for (i=0; i < num_cols_offd_C; i++)
	  col_map_offd_C[i] = temp[i];

       jxf_TFree(temp);
   
      if (C_ext_diag_size)
      {
         C_ext_diag_j = jxf_CTAlloc(JXF_Int, C_ext_diag_size);
         C_ext_diag_data = jxf_CTAlloc(JXF_Real, C_ext_diag_size);
      }
      if (C_ext_offd_size)
      {
         C_ext_offd_j = jxf_CTAlloc(JXF_Int, C_ext_offd_size);
         C_ext_offd_data = jxf_CTAlloc(JXF_Real, C_ext_offd_size);
      }

      C_tmp_diag_j = jxf_CSRMatrixJ(C_tmp_diag);
      C_tmp_diag_data = jxf_CSRMatrixData(C_tmp_diag);

      C_tmp_offd_i = jxf_CSRMatrixI(C_tmp_offd);
      C_tmp_offd_j = jxf_CSRMatrixJ(C_tmp_offd);
      C_tmp_offd_data = jxf_CSRMatrixData(C_tmp_offd);

      cnt_offd = 0;
      cnt_diag = 0;
      for (i=0; i < C_ext_num_rows; i++)
      {
         for (j=C_ext_i[i]; j < C_ext_i[i+1]; j++)
            if (C_ext_j[j] < first_col_diag_C || C_ext_j[j] > last_col_diag_C)
            {
               C_ext_offd_j[cnt_offd] = jxf_BinarySearch(col_map_offd_C,
                                           C_ext_j[j],
                                           num_cols_offd_C);
               C_ext_offd_data[cnt_offd++] = C_ext_data[j];
            }
            else
            {
               C_ext_diag_j[cnt_diag] = C_ext_j[j] - first_col_diag_C;
               C_ext_diag_data[cnt_diag++] = C_ext_data[j];
            }
      }
   }

   if (C_ext)
   {
      jxf_CSRMatrixDestroy(C_ext);
      C_ext = NULL;
   }

   if (num_cols_offd_B)
   {
      map_B_to_C = jxf_CTAlloc(JXF_Int,num_cols_offd_B);

      cnt = 0;
      for (i=0; i < num_cols_offd_C; i++)
         if (col_map_offd_C[i] == col_map_offd_B[cnt])
         {
            map_B_to_C[cnt++] = i;
            if (cnt == num_cols_offd_B) break;
         }
      for (i=0; 
	i < jxf_CSRMatrixI(C_tmp_offd)[jxf_CSRMatrixNumRows(C_tmp_offd)]; i++)
      {
         j_indx = C_tmp_offd_j[i];
         C_tmp_offd_j[i] = map_B_to_C[j_indx];
      }
   }

   /*-----------------------------------------------------------------------
    *  Need to compute C_diag = C_tmp_diag + C_ext_diag
    *  and  C_offd = C_tmp_offd + C_ext_offd   !!!!
    *  First generate structure
    *-----------------------------------------------------------------------*/

   if (C_ext_size || num_cols_offd_B)
   {
     C_diag_i = jxf_CTAlloc(JXF_Int, num_cols_diag_A+1);
     C_offd_i = jxf_CTAlloc(JXF_Int, num_cols_diag_A+1);

     C_diag_array = jxf_CTAlloc(JXF_Int, max_num_threads);
     C_offd_array = jxf_CTAlloc(JXF_Int, max_num_threads);

#if JXF_USING_OPENMP
#pragma omp parallel
#endif
     {
        JXF_Int *B_marker = NULL;
        JXF_Int *B_marker_offd = NULL;
        JXF_Int ik, jk, j1, j2, jcol;
        JXF_Int ns, ne, ii, nnz_d, nnz_o;
        JXF_Int rest, size;
        JXF_Int num_threads = jxf_NumActiveThreads();

        size = num_cols_diag_A/num_threads;
        rest = num_cols_diag_A - size*num_threads;
        ii = jxf_GetThreadNum();
        if (ii < rest)
        {
           ns = ii*size+ii;
           ne = (ii+1)*size+ii+1;
        }
        else
        {
           ns = ii*size+rest;
           ne = (ii+1)*size+rest;
        }

        B_marker = jxf_CTAlloc(JXF_Int, num_cols_diag_B);
        B_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd_C);

        for (ik = 0; ik < num_cols_diag_B; ik++)
           B_marker[ik] = -1;

        for (ik = 0; ik < num_cols_offd_C; ik++)
           B_marker_offd[ik] = -1;

        nnz_d = 0;
        nnz_o = 0;
        for (ik = ns; ik < ne; ik++)
        {
          for (jk = C_tmp_diag_i[ik]; jk < C_tmp_diag_i[ik+1]; jk++)
          {
             jcol = C_tmp_diag_j[jk];
             B_marker[jcol] = ik;
	     nnz_d++;
          }
          for (jk = C_tmp_offd_i[ik]; jk < C_tmp_offd_i[ik+1]; jk++)
          {
             jcol = C_tmp_offd_j[jk];
             B_marker_offd[jcol] = ik;
	     nnz_o++;
          }
          for (jk = 0; jk < num_sends_A; jk++)
            for (j1 = send_map_starts_A[jk]; j1 < send_map_starts_A[jk+1]; j1++)
             if (send_map_elmts_A[j1] == ik)
             {
                for (j2 = C_ext_diag_i[j1]; j2 < C_ext_diag_i[j1+1]; j2++)
                {
                    jcol = C_ext_diag_j[j2];
 	            if (B_marker[jcol] < ik)
                    {
                       B_marker[jcol] = ik;
	               nnz_d++;
                    }
                }
                for (j2 = C_ext_offd_i[j1]; j2 < C_ext_offd_i[j1+1]; j2++)
                {
                    jcol = C_ext_offd_j[j2];
 	            if (B_marker_offd[jcol] < ik)
                    {
                       B_marker_offd[jcol] = ik;
	               nnz_o++;
                    }
                }
                break;
             }
               C_diag_array[ii] = nnz_d;
               C_offd_array[ii] = nnz_o;
        }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
       if (ii == 0)
       {
          nnz_d = 0;
          nnz_o = 0;
          for (ik = 0; ik < num_threads-1; ik++)
          {
             C_diag_array[ik+1] += C_diag_array[ik];
             C_offd_array[ik+1] += C_offd_array[ik];
          }
          nnz_d = C_diag_array[num_threads-1];
          nnz_o = C_offd_array[num_threads-1];
          C_diag_i[num_cols_diag_A] = nnz_d;
          C_offd_i[num_cols_diag_A] = nnz_o;

          C_diag = jxf_CSRMatrixCreate(num_cols_diag_A, num_cols_diag_A, nnz_d);
          C_offd = jxf_CSRMatrixCreate(num_cols_diag_A, num_cols_offd_C, nnz_o);
          jxf_CSRMatrixI(C_diag) = C_diag_i;
          jxf_CSRMatrixInitialize(C_diag);
          C_diag_j = jxf_CSRMatrixJ(C_diag);
          C_diag_data = jxf_CSRMatrixData(C_diag);
          jxf_CSRMatrixI(C_offd) = C_offd_i;
          jxf_CSRMatrixInitialize(C_offd);
          C_offd_j = jxf_CSRMatrixJ(C_offd);
          C_offd_data = jxf_CSRMatrixData(C_offd); 
       }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif

   /*-----------------------------------------------------------------------
    *  Need to compute C_diag = C_tmp_diag + C_ext_diag
    *  and  C_offd = C_tmp_offd + C_ext_offd   !!!!
    *  Now fill in values
    *-----------------------------------------------------------------------*/
   
     for (ik = 0; ik < num_cols_diag_B; ik++)
        B_marker[ik] = -1;

     for (ik = 0; ik < num_cols_offd_C; ik++)
        B_marker_offd[ik] = -1;

   /*-----------------------------------------------------------------------
    *  Populate matrices
    *-----------------------------------------------------------------------*/

      nnz_d = 0;
      nnz_o = 0;
        nnz_o = 0;
        if (ii)
        {
           nnz_d = C_diag_array[ii-1];
           nnz_o = C_offd_array[ii-1];
        }
        for (ik = ns; ik < ne; ik++)
        {
           C_diag_i[ik] = nnz_d;
           C_offd_i[ik] = nnz_o;
           for (jk = C_tmp_diag_i[ik]; jk < C_tmp_diag_i[ik+1]; jk++)
           {
              jcol = C_tmp_diag_j[jk];
              C_diag_j[nnz_d] = jcol;
              C_diag_data[nnz_d] = C_tmp_diag_data[jk];
              B_marker[jcol] = nnz_d;
              nnz_d++;
           }
           for (jk = C_tmp_offd_i[ik]; jk < C_tmp_offd_i[ik+1]; jk++)
           {
              jcol = C_tmp_offd_j[jk];
              C_offd_j[nnz_o] = jcol;
              C_offd_data[nnz_o] = C_tmp_offd_data[jk];
              B_marker_offd[jcol] = nnz_o;
              nnz_o++;
           }
           for (jk = 0; jk < num_sends_A; jk++)
              for (j1 = send_map_starts_A[jk]; j1 < send_map_starts_A[jk+1]; j1++)
                 if (send_map_elmts_A[j1] == ik)
                 {
                    for (j2 = C_ext_diag_i[j1]; j2 < C_ext_diag_i[j1+1]; j2++)
                    {
                       jcol = C_ext_diag_j[j2];
                       if (B_marker[jcol] < C_diag_i[ik])
                       {
                          C_diag_j[nnz_d] = jcol;
                          C_diag_data[nnz_d] = C_ext_diag_data[j2];
                          B_marker[jcol] = nnz_d;
                          nnz_d++;
                       }
                       else
                          C_diag_data[B_marker[jcol]] += C_ext_diag_data[j2];
                    }
                    for (j2 = C_ext_offd_i[j1]; j2 < C_ext_offd_i[j1+1]; j2++)
                    {
                       jcol = C_ext_offd_j[j2];
                       if (B_marker_offd[jcol] < C_offd_i[ik])
                       {
                          C_offd_j[nnz_o] = jcol;
                          C_offd_data[nnz_o] = C_ext_offd_data[j2];
                          B_marker_offd[jcol] = nnz_o;
                          nnz_o++;
                       }
                       else
                          C_offd_data[B_marker_offd[jcol]] += C_ext_offd_data[j2];
                    }
                    break;
                 }
        }
        jxf_TFree(B_marker);
        jxf_TFree(B_marker_offd);
     } /*end parallel region */
     jxf_TFree(C_diag_array);
     jxf_TFree(C_offd_array);
   }

   /*C = jxf_ParCSRMatrixCreate(comm, n_cols_A, n_cols_B, col_starts_A,
	col_starts_B, num_cols_offd_C, nnz_diag, nnz_offd);

   jxf_CSRMatrixDestroy(jxf_ParCSRMatrixDiag(C));
   jxf_CSRMatrixDestroy(jxf_ParCSRMatrixOffd(C)); */
#ifdef JXF_NO_GLOBAL_PARTITION
   /* row_starts[0] is start of local rows.  row_starts[1] is start of next 
      processor's rows */
   first_row_index = col_starts_A[0];
   local_num_rows = col_starts_A[1]-first_row_index ;
   first_col_diag = col_starts_B[0];
   local_num_cols = col_starts_B[1]-first_col_diag;
#else
   first_row_index = col_starts_A[my_id];
   local_num_rows = col_starts_A[my_id+1]-first_row_index;
   first_col_diag = col_starts_B[my_id];
   local_num_cols = col_starts_B[my_id+1]-first_col_diag;
#endif

   C = jxf_CTAlloc(jxf_ParCSRMatrix, 1);
   jxf_ParCSRMatrixComm(C) = comm;
   jxf_ParCSRMatrixGlobalNumRows(C) = n_cols_A;
   jxf_ParCSRMatrixGlobalNumCols(C) = n_cols_B;
   jxf_ParCSRMatrixFirstRowIndex(C) = first_row_index;
   jxf_ParCSRMatrixFirstColDiag(C) = first_col_diag;
   jxf_ParCSRMatrixLastRowIndex(C) = first_row_index + local_num_rows - 1;
   jxf_ParCSRMatrixLastColDiag(C) = first_col_diag + local_num_cols - 1;

   jxf_ParCSRMatrixColMapOffd(C) = NULL;

   jxf_ParCSRMatrixAssumedPartition(C) = NULL;

   jxf_ParCSRMatrixRowStarts(C) = col_starts_A;
   jxf_ParCSRMatrixColStarts(C) = col_starts_B;

   jxf_ParCSRMatrixCommPkg(C) = NULL;
   jxf_ParCSRMatrixCommPkgT(C) = NULL;

   /* set defaults */
   jxf_ParCSRMatrixOwnsData(C) = 1;
   jxf_ParCSRMatrixRowindices(C) = NULL;
   jxf_ParCSRMatrixRowvalues(C) = NULL;
   jxf_ParCSRMatrixGetrowactive(C) = 0;

/* Note that C does not own the partitionings */
   jxf_ParCSRMatrixSetRowStartsOwner(C,0);
   jxf_ParCSRMatrixSetColStartsOwner(C,0);

   if (C_diag) jxf_ParCSRMatrixDiag(C) = C_diag;
   else jxf_ParCSRMatrixDiag(C) = C_tmp_diag;
   if (C_offd) jxf_ParCSRMatrixOffd(C) = C_offd;
   else jxf_ParCSRMatrixOffd(C) = C_tmp_offd;

   if (num_cols_offd_C)
   {
      JXF_Int jj_count_offd, nnz_offd;
      JXF_Int *new_col_map_offd_C = NULL;

      P_marker = jxf_CTAlloc(JXF_Int,num_cols_offd_C);
      for (i=0; i < num_cols_offd_C; i++)
         P_marker[i] = -1;

      jj_count_offd = 0;
      nnz_offd = C_offd_i[num_cols_diag_A];
      for (i=0; i < nnz_offd; i++)
      {
         i1 = C_offd_j[i];
         if (P_marker[i1])
         {
            P_marker[i1] = 0;
            jj_count_offd++;
         }
      }

      if (jj_count_offd < num_cols_offd_C)
      {
         new_col_map_offd_C = jxf_CTAlloc(JXF_Int,jj_count_offd);
         jj_count_offd = 0;
         for (i=0; i < num_cols_offd_C; i++)
            if (!P_marker[i])
            {
               P_marker[i] = jj_count_offd;
               new_col_map_offd_C[jj_count_offd++] = col_map_offd_C[i];
            }

         for (i=0; i < nnz_offd; i++)
         {
            i1 = C_offd_j[i];
            C_offd_j[i] = P_marker[i1];
         }

         num_cols_offd_C = jj_count_offd;
         jxf_TFree(col_map_offd_C);
         col_map_offd_C = new_col_map_offd_C;
         jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(C)) = num_cols_offd_C;
      }
      jxf_TFree(P_marker);
   }
   jxf_ParCSRMatrixColMapOffd(C) = col_map_offd_C;
   /*-----------------------------------------------------------------------
    *  Free various arrays
    *-----------------------------------------------------------------------*/

   if (C_ext_size || num_cols_offd_B)
   {
      jxf_TFree(C_ext_diag_i);
      jxf_TFree(C_ext_offd_i);
   }
   if (C_ext_diag_size)
   {
      jxf_TFree(C_ext_diag_j);
      jxf_TFree(C_ext_diag_data);
   }
   if (C_ext_offd_size)
   {
      jxf_TFree(C_ext_offd_j);
      jxf_TFree(C_ext_offd_data);
   }
   if (num_cols_offd_B) jxf_TFree(map_B_to_C);

   if (C_diag) jxf_CSRMatrixDestroy(C_tmp_diag);
   if (C_offd) jxf_CSRMatrixDestroy(C_tmp_offd);

   return C;
}

/*!
 * \fn void jxf_ParMatmul_RowSizes
 * \brief The following function was formerly part of jxf_ParMatmul
 *        but was removed so it can also be used for
 *        multiplication of Boolean matrices
 * \date 2017/02/25
 */
void
jxf_ParMatmul_RowSizes( JXF_Int **C_diag_i,
                       JXF_Int **C_offd_i,
                       JXF_Int  *A_diag_i,
                       JXF_Int  *A_diag_j,
                       JXF_Int  *A_offd_i,
                       JXF_Int  *A_offd_j,
                       JXF_Int  *B_diag_i,
                       JXF_Int  *B_diag_j,
                       JXF_Int  *B_offd_i,
                       JXF_Int  *B_offd_j,
                       JXF_Int  *B_ext_diag_i,
                       JXF_Int  *B_ext_diag_j, 
                       JXF_Int  *B_ext_offd_i,
                       JXF_Int  *B_ext_offd_j,
                       JXF_Int  *map_B_to_C,
                       JXF_Int  *C_diag_size,
                       JXF_Int  *C_offd_size,
                       JXF_Int   num_rows_diag_A,
                       JXF_Int   num_cols_offd_A,
                       JXF_Int   allsquare,
                       JXF_Int   num_cols_diag_B,
                       JXF_Int   num_cols_offd_B,
                       JXF_Int   num_cols_offd_C )
{
   JXF_Int i1, i2, i3, jj2, jj3;
   JXF_Int jj_count_diag, jj_count_offd, jj_row_begin_diag, jj_row_begin_offd;
   JXF_Int start_indexing = 0; /* start indexing for C_data at 0 */
   JXF_Int num_threads = jxf_NumThreads();
   JXF_Int *jj_count_diag_array;
   JXF_Int *jj_count_offd_array;
   JXF_Int ii, size, rest;

   *C_diag_i = jxf_CTAlloc(JXF_Int, num_rows_diag_A+1);
   *C_offd_i = jxf_CTAlloc(JXF_Int, num_rows_diag_A+1);
   jj_count_diag_array = jxf_CTAlloc(JXF_Int, num_threads);
   jj_count_offd_array = jxf_CTAlloc(JXF_Int, num_threads);
   /*-----------------------------------------------------------------------
    *  Loop over rows of A
    *-----------------------------------------------------------------------*/
   size = num_rows_diag_A/num_threads;
   rest = num_rows_diag_A - size*num_threads;

#if JXF_USING_OPENMP
#pragma omp parallel private(ii, i1, jj_row_begin_diag, jj_row_begin_offd, jj_count_diag, jj_count_offd, jj2, i2, jj3, i3) 
#endif
   /*for (ii=0; ii < num_threads; ii++)*/
   {
    JXF_Int *B_marker = NULL;
    JXF_Int ns, ne;
    ii = jxf_GetThreadNum();
    if (ii < rest)
    {
       ns = ii*size+ii;
       ne = (ii+1)*size+ii+1;
    }
    else
    {
       ns = ii*size+rest;
       ne = (ii+1)*size+rest;
    }
    jj_count_diag = start_indexing;
    jj_count_offd = start_indexing;

    if (num_cols_diag_B || num_cols_offd_C)
    B_marker = jxf_CTAlloc(JXF_Int, num_cols_diag_B+num_cols_offd_C);
    for (i1 = 0; i1 < num_cols_diag_B+num_cols_offd_C; i1++)
      B_marker[i1] = -1;

    for (i1 = ns; i1 < ne; i1++)
    {
      /*--------------------------------------------------------------------
       *  Set marker for diagonal entry, C_{i1,i1} (for square matrices). 
       *--------------------------------------------------------------------*/
 
      jj_row_begin_diag = jj_count_diag;
      jj_row_begin_offd = jj_count_offd;
      if ( allsquare ) {
         B_marker[i1] = jj_count_diag;
         jj_count_diag++;
      }

      /*-----------------------------------------------------------------
       *  Loop over entries in row i1 of A_offd.
       *-----------------------------------------------------------------*/
         
      if (num_cols_offd_A)
      {
         for (jj2 = A_offd_i[i1]; jj2 < A_offd_i[i1+1]; jj2++)
         {
            i2 = A_offd_j[jj2];
 
            /*-----------------------------------------------------------
             *  Loop over entries in row i2 of B_ext.
             *-----------------------------------------------------------*/
 
            for (jj3 = B_ext_offd_i[i2]; jj3 < B_ext_offd_i[i2+1]; jj3++)
            {
               i3 = num_cols_diag_B+B_ext_offd_j[jj3];
                  
               /*--------------------------------------------------------
                *  Check B_marker to see that C_{i1,i3} has not already
                *  been accounted for. If it has not, mark it and increment
                *  counter.
                *--------------------------------------------------------*/

               if (B_marker[i3] < jj_row_begin_offd)
               {
                  B_marker[i3] = jj_count_offd;
                  jj_count_offd++;
               } 
            }
            for (jj3 = B_ext_diag_i[i2]; jj3 < B_ext_diag_i[i2+1]; jj3++)
            {
               i3 = B_ext_diag_j[jj3];
                  
               if (B_marker[i3] < jj_row_begin_diag)
               {
                  B_marker[i3] = jj_count_diag;
                  jj_count_diag++;
               } 
            }
         }
      }

      /*-----------------------------------------------------------------
       *  Loop over entries in row i1 of A_diag.
       *-----------------------------------------------------------------*/
         
      for (jj2 = A_diag_i[i1]; jj2 < A_diag_i[i1+1]; jj2++)
      {
         i2 = A_diag_j[jj2];
 
         /*-----------------------------------------------------------
          *  Loop over entries in row i2 of B_diag.
          *-----------------------------------------------------------*/
 
         for (jj3 = B_diag_i[i2]; jj3 < B_diag_i[i2+1]; jj3++)
         {
            i3 = B_diag_j[jj3];
                  
            /*--------------------------------------------------------
             *  Check B_marker to see that C_{i1,i3} has not already
             *  been accounted for. If it has not, mark it and increment
             *  counter.
             *--------------------------------------------------------*/
 
            if (B_marker[i3] < jj_row_begin_diag)
            {
               B_marker[i3] = jj_count_diag;
               jj_count_diag++;
            }
         }

         /*-----------------------------------------------------------
          *  Loop over entries in row i2 of B_offd.
          *-----------------------------------------------------------*/

         if (num_cols_offd_B)
         { 
            for (jj3 = B_offd_i[i2]; jj3 < B_offd_i[i2+1]; jj3++)
            {
               i3 = num_cols_diag_B+map_B_to_C[B_offd_j[jj3]];
                  
               /*--------------------------------------------------------
                *  Check B_marker to see that C_{i1,i3} has not already
                *  been accounted for. If it has not, mark it and increment
                *  counter.
                *--------------------------------------------------------*/
 
               if (B_marker[i3] < jj_row_begin_offd)
               {
                  B_marker[i3] = jj_count_offd;
                  jj_count_offd++;
               }
            }
         }
      }
            
      /*--------------------------------------------------------------------
       * Set C_diag_i and C_offd_i for this row.
       *--------------------------------------------------------------------*/
 
      (*C_diag_i)[i1] = jj_row_begin_diag;
      (*C_offd_i)[i1] = jj_row_begin_offd;
      
    }
    jj_count_diag_array[ii] = jj_count_diag;
    jj_count_offd_array[ii] = jj_count_offd;

    jxf_TFree(B_marker);
#if JXF_USING_OPENMP
#pragma omp barrier
#endif

    if (ii)
    {
       jj_count_diag = jj_count_diag_array[0];
       jj_count_offd = jj_count_offd_array[0];
       for (i1 = 1; i1 < ii; i1++)
       {
          jj_count_diag += jj_count_diag_array[i1];
          jj_count_offd += jj_count_offd_array[i1];
       }

       for (i1 = ns; i1 < ne; i1++)
       {
          (*C_diag_i)[i1] += jj_count_diag;
          (*C_offd_i)[i1] += jj_count_offd;
       }
    }
    else
    {
       (*C_diag_i)[num_rows_diag_A] = 0;
       (*C_offd_i)[num_rows_diag_A] = 0;
       for (i1 = 0; i1 < num_threads; i1++)
       {
          (*C_diag_i)[num_rows_diag_A] += jj_count_diag_array[i1];
          (*C_offd_i)[num_rows_diag_A] += jj_count_offd_array[i1];
       }
    }
   } /* end parallel loop */
 
   /*-----------------------------------------------------------------------
    *  Allocate C_diag_data and C_diag_j arrays.
    *  Allocate C_offd_data and C_offd_j arrays.
    *-----------------------------------------------------------------------*/

   *C_diag_size = (*C_diag_i)[num_rows_diag_A];
   *C_offd_size = (*C_offd_i)[num_rows_diag_A];

   jxf_TFree(jj_count_diag_array);
   jxf_TFree(jj_count_offd_array);
}

/* C = alpha * A + beta * B
 * A and B are assumed to have the same row and column partitionings */
JXF_Int
jxf_ParcsrAdd( JXF_Real alpha,
              jxf_ParCSRMatrix *A,
              JXF_Real beta,
              jxf_ParCSRMatrix *B,
              jxf_ParCSRMatrix **Cout )
{
   MPI_Comm         comm     = jxf_ParCSRMatrixComm(A);
   JXF_Int        num_procs, my_id;
   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   JXF_Int i, j;

   /* diag part of A */
   jxf_CSRMatrix *A_diag   = jxf_ParCSRMatrixDiag(A);
   JXF_Real   *A_diag_a = jxf_CSRMatrixData(A_diag);
   JXF_Int       *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int       *A_diag_j = jxf_CSRMatrixJ(A_diag);
   /* off-diag part of A */
   jxf_CSRMatrix *A_offd   = jxf_ParCSRMatrixOffd(A);
   JXF_Real   *A_offd_a = jxf_CSRMatrixData(A_offd);
   JXF_Int       *A_offd_i = jxf_CSRMatrixI(A_offd);
   JXF_Int       *A_offd_j = jxf_CSRMatrixJ(A_offd);

   JXF_Int        num_cols_A_offd = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int    *col_map_offd_A  = jxf_ParCSRMatrixColMapOffd(A);
   JXF_Int       *A2C_offd = jxf_TAlloc(JXF_Int, num_cols_A_offd);

   JXF_Int     nrow_global = jxf_ParCSRMatrixGlobalNumRows(A);
   JXF_Int     ncol_global = jxf_ParCSRMatrixGlobalNumCols(A);
   JXF_Int        nrow_local = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int        ncol_local = jxf_CSRMatrixNumCols(A_diag);
   JXF_Int        nnz_diag_A = A_diag_i[nrow_local];
   JXF_Int        nnz_offd_A = A_offd_i[nrow_local];

   /* diag part of B */
   jxf_CSRMatrix *B_diag   = jxf_ParCSRMatrixDiag(B);
   JXF_Real   *B_diag_a = jxf_CSRMatrixData(B_diag);
   JXF_Int       *B_diag_i = jxf_CSRMatrixI(B_diag);
   JXF_Int       *B_diag_j = jxf_CSRMatrixJ(B_diag);
   /* off-diag part of B */
   jxf_CSRMatrix *B_offd   = jxf_ParCSRMatrixOffd(B);
   JXF_Real   *B_offd_a = jxf_CSRMatrixData(B_offd);
   JXF_Int       *B_offd_i = jxf_CSRMatrixI(B_offd);
   JXF_Int       *B_offd_j = jxf_CSRMatrixJ(B_offd);

   JXF_Int        num_cols_B_offd = jxf_CSRMatrixNumCols(B_offd);
   JXF_Int    *col_map_offd_B  = jxf_ParCSRMatrixColMapOffd(B);
   JXF_Int       *B2C_offd = jxf_TAlloc(JXF_Int, num_cols_B_offd);

   jxf_assert(nrow_global == jxf_ParCSRMatrixGlobalNumRows(B));
   jxf_assert(ncol_global == jxf_ParCSRMatrixGlobalNumCols(B));
   jxf_assert(nrow_local == jxf_CSRMatrixNumRows(B_diag));
   jxf_assert(ncol_local == jxf_CSRMatrixNumCols(B_diag));
   JXF_Int        nnz_diag_B = B_diag_i[nrow_local];
   JXF_Int        nnz_offd_B = B_offd_i[nrow_local];

   /* C */
   jxf_ParCSRMatrix *C;
   JXF_Int       *row_starts_C, *col_starts_C;
   jxf_CSRMatrix    *C_diag;
   jxf_CSRMatrix    *C_offd;

   JXF_Int        num_cols_C_offd = num_cols_A_offd + num_cols_B_offd;
   JXF_Int    *col_map_offd_C = jxf_TAlloc(JXF_Int, num_cols_C_offd);

   JXF_Int        nnz_diag_C_alloc = nnz_diag_A + nnz_diag_B;
   JXF_Int        nnz_offd_C_alloc = nnz_offd_A + nnz_offd_B;
   JXF_Int        nnz_diag_C = 0, nnz_offd_C = 0;

   JXF_Int     *C_diag_i = jxf_CTAlloc(JXF_Int,     nrow_local + 1);
   JXF_Int     *C_diag_j = jxf_CTAlloc(JXF_Int,     nnz_diag_C_alloc);
   JXF_Real *C_diag_a = jxf_CTAlloc(JXF_Real, nnz_diag_C_alloc);
   JXF_Int     *C_offd_i = jxf_CTAlloc(JXF_Int,     nrow_local + 1);
   JXF_Int     *C_offd_j = jxf_CTAlloc(JXF_Int,     nnz_offd_C_alloc);
   JXF_Real *C_offd_a = jxf_CTAlloc(JXF_Real, nnz_offd_C_alloc);

   jxf_union2( num_cols_A_offd, col_map_offd_A, num_cols_B_offd, col_map_offd_B,
                 &num_cols_C_offd, col_map_offd_C, A2C_offd, B2C_offd );

   JXF_Int     *marker_diag = jxf_TAlloc(JXF_Int, ncol_local);
   JXF_Int     *marker_offd = jxf_TAlloc(JXF_Int, num_cols_C_offd);

   for (i = 0; i < ncol_local; i++)
   {
      marker_diag[i] = -1;
   }
   for (i = 0; i < num_cols_C_offd; i++)
   {
      marker_offd[i] = -1;
   }

   /* main loop for each row i */
   for (i = 0; i < nrow_local; i++)
   {
      JXF_Int diag_i_start = nnz_diag_C;
      JXF_Int offd_i_start = nnz_offd_C;

      for (j = A_diag_i[i]; j < A_diag_i[i+1]; j++)
      {
         JXF_Int     col = A_diag_j[j];
         JXF_Real val = A_diag_a[j];
         if (marker_diag[col] < diag_i_start)
         {
            /* this col has not been seen before, create new entry */
            marker_diag[col] = nnz_diag_C;
            C_diag_j[nnz_diag_C] = col;
            C_diag_a[nnz_diag_C] = alpha * val;
            nnz_diag_C ++;
         }
         else
         {
            /* this should not happen */
            jxf_printf("jxfpamg warning: invalid ParCSR matrix %s %s %d\n",
                         __FILE__, __func__, __LINE__);
         }
      }

      for (j = B_diag_i[i]; j < B_diag_i[i+1]; j++)
      {
         JXF_Int     col = B_diag_j[j];
         JXF_Real val = B_diag_a[j];
         if (marker_diag[col] < diag_i_start /*&& jxf_abs(val) > 0.0*/)
         {
            /* this col has not been seen before, create new entry */
            marker_diag[col] = nnz_diag_C;
            C_diag_j[nnz_diag_C] = col;
            C_diag_a[nnz_diag_C] = beta * val;
            nnz_diag_C ++;
         }
         else
         {
            /* existing entry, update */
            JXF_Int p = marker_diag[col];

            jxf_assert(C_diag_j[p] == col);

            C_diag_a[p] += beta * val;
         }
      }

      C_diag_i[i+1] = nnz_diag_C;

      if (num_procs <= 1)
      {
         continue;
      }

      for (j = A_offd_i[i]; j < A_offd_i[i+1]; j++)
      {
         JXF_Int     colA = A_offd_j[j];
         JXF_Int     colC = A2C_offd[colA];
         JXF_Real val  = A_offd_a[j];
         if (marker_offd[colC] < offd_i_start)
         {
            /* this col has not been seen before, create new entry */
            marker_offd[colC] = nnz_offd_C;
            C_offd_j[nnz_offd_C] = colC;
            C_offd_a[nnz_offd_C] = alpha * val;
            nnz_offd_C ++;
         }
         else
         {
            /* this should not happen */
            jxf_printf("jxfpamg warning: invalid ParCSR matrix %s %s %d\n",
                         __FILE__, __func__, __LINE__);
         }
      }

      for (j = B_offd_i[i]; j < B_offd_i[i+1]; j++)
      {
         JXF_Int     colB = B_offd_j[j];
         JXF_Int     colC = B2C_offd[colB];
         JXF_Real val  = B_offd_a[j];
         if (marker_offd[colC] < offd_i_start /*&& jxf_abs(val) > 0.0*/)
         {
            /* this col has not been seen before, create new entry */
            marker_offd[colC] = nnz_offd_C;
            C_offd_j[nnz_offd_C] = colC;
            C_offd_a[nnz_offd_C] = beta * val;
            nnz_offd_C ++;
         }
         else
         {
            /* existing entry, update */
            JXF_Int p = marker_offd[colC];

            jxf_assert(C_offd_j[p] == colC);

            C_offd_a[p] += beta * val;
         }
      }

      C_offd_i[i+1] = nnz_offd_C;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   j = 2;
#else
   j = num_procs + 1;
#endif

   row_starts_C = jxf_TAlloc(JXF_Int, j);
   col_starts_C = jxf_TAlloc(JXF_Int, j);
   memcpy(row_starts_C, jxf_ParCSRMatrixRowStarts(A), j*sizeof(JXF_Int));
   memcpy(col_starts_C, jxf_ParCSRMatrixColStarts(A), j*sizeof(JXF_Int));

   /* Now, we should have everything of Parcsr matrix C */
   C = jxf_ParCSRMatrixCreate(comm,
                                nrow_global,
                                ncol_global,
                                row_starts_C,
                                col_starts_C,
                                num_cols_C_offd,
                                nnz_diag_C,
                                nnz_offd_C);

   C_diag = jxf_ParCSRMatrixDiag(C);
   jxf_CSRMatrixData(C_diag) = C_diag_a;
   jxf_CSRMatrixI(C_diag)    = C_diag_i;
   jxf_CSRMatrixJ(C_diag)    = C_diag_j;

   C_offd = jxf_ParCSRMatrixOffd(C);
   jxf_CSRMatrixData(C_offd) = C_offd_a;
   jxf_CSRMatrixI(C_offd)    = C_offd_i;
   jxf_CSRMatrixJ(C_offd)    = C_offd_j;

   jxf_ParCSRMatrixColMapOffd(C) = col_map_offd_C;

   jxf_ParCSRMatrixSetNumNonzeros(C);
   jxf_ParCSRMatrixDNumNonzeros(C) = (JXF_Real) jxf_ParCSRMatrixNumNonzeros(C);

   /* create CommPkg of C */
   jxf_MatvecCommPkgCreate(C);

   *Cout = C;

   /* done */
   jxf_TFree(A2C_offd);
   jxf_TFree(B2C_offd);
   jxf_TFree(marker_diag);
   jxf_TFree(marker_offd);

   return jxf_error_flag;
}

/* drop the entries that are not on the diagonal and smaller than
 * its row norm: type 1: 1-norm, 2: 2-norm, -1: infinity norm */
JXF_Int
jxf_ParCSRMatrixDropSmallEntries( jxf_ParCSRMatrix *A, JXF_Real tol, JXF_Int type )
{
   JXF_Int i, j, k, nnz_diag, nnz_offd, A_diag_i_i, A_offd_i_i;

   MPI_Comm         comm     = jxf_ParCSRMatrixComm(A);
   /* diag part of A */
   jxf_CSRMatrix *A_diag   = jxf_ParCSRMatrixDiag(A);
   JXF_Real      *A_diag_a = jxf_CSRMatrixData(A_diag);
   JXF_Int       *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int       *A_diag_j = jxf_CSRMatrixJ(A_diag);
   /* off-diag part of A */
   jxf_CSRMatrix *A_offd   = jxf_ParCSRMatrixOffd(A);
   JXF_Real      *A_offd_a = jxf_CSRMatrixData(A_offd);
   JXF_Int       *A_offd_i = jxf_CSRMatrixI(A_offd);
   JXF_Int       *A_offd_j = jxf_CSRMatrixJ(A_offd);

   JXF_Int  num_cols_A_offd = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int *col_map_offd_A  = jxf_ParCSRMatrixColMapOffd(A);
   JXF_Int *marker_offd = NULL;

   JXF_Int first_row  = jxf_ParCSRMatrixFirstRowIndex(A);
   JXF_Int nrow_local = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int my_id, num_procs;
   /* MPI size and rank*/
   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (tol <= 0.0)
   {
      return jxf_error_flag;
   }

   marker_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);

   nnz_diag = nnz_offd = A_diag_i_i = A_offd_i_i = 0;
   for (i = 0; i < nrow_local; i++)
   {
      /* compute row norm */
      JXF_Real row_nrm = 0.0;
      for (j = A_diag_i_i; j < A_diag_i[i+1]; j++)
      {
         JXF_Real v = A_diag_a[j];
         if (type == 1)
         {
            row_nrm += fabs(v);
         }
         else if (type == 2)
         {
            row_nrm += v*v;
         }
         else
         {
            row_nrm = jxf_max(row_nrm, fabs(v));
         }
      }
      if (num_procs > 1)
      {
         for (j = A_offd_i_i; j < A_offd_i[i+1]; j++)
         {
            JXF_Real v = A_offd_a[j];
            if (type == 1)
            {
               row_nrm += fabs(v);
            }
            else if (type == 2)
            {
               row_nrm += v*v;
            }
            else
            {
               row_nrm = jxf_max(row_nrm, fabs(v));
            }
         }
      }

      if (type == 2)
      {
         row_nrm = sqrt(row_nrm);
      }

      /* drop small entries based on tol and row norm */
      for (j = A_diag_i_i; j < A_diag_i[i+1]; j++)
      {
         JXF_Int     col = A_diag_j[j];
         JXF_Real val = A_diag_a[j];
         if (i == col || fabs(val) >= tol * row_nrm)
         {
            A_diag_j[nnz_diag] = col;
            A_diag_a[nnz_diag] = val;
            nnz_diag ++;
         }
      }
      if (num_procs > 1)
      {
         for (j = A_offd_i_i; j < A_offd_i[i+1]; j++)
         {
            JXF_Int     col = A_offd_j[j];
            JXF_Real val = A_offd_a[j];
            /* in normal cases: diagonal entry should not
             * appear in A_offd (but this can still be possible) */
            if (i + first_row == col_map_offd_A[col] || fabs(val) >= tol * row_nrm)
            {
               if (0 == marker_offd[col])
               {
                  marker_offd[col] = 1;
               }
               A_offd_j[nnz_offd] = col;
               A_offd_a[nnz_offd] = val;
               nnz_offd ++;
            }
         }
      }
      A_diag_i_i = A_diag_i[i+1];
      A_offd_i_i = A_offd_i[i+1];
      A_diag_i[i+1] = nnz_diag;
      A_offd_i[i+1] = nnz_offd;
   }

   jxf_CSRMatrixNumNonzeros(A_diag) = nnz_diag;
   jxf_CSRMatrixNumNonzeros(A_offd) = nnz_offd;
   jxf_ParCSRMatrixSetNumNonzeros(A);
   jxf_ParCSRMatrixDNumNonzeros(A) = (JXF_Real) jxf_ParCSRMatrixNumNonzeros(A);

   for (i = 0, k = 0; i < num_cols_A_offd; i++)
   {
      if (marker_offd[i])
      {
         col_map_offd_A[k] = col_map_offd_A[i];
         marker_offd[i] = k++;
      }
   }
   /* num_cols_A_offd = k; */
   jxf_CSRMatrixNumCols(A_offd) = k;
   for (i = 0; i < nnz_offd; i++)
   {
      A_offd_j[i] = marker_offd[A_offd_j[i]];
   }

   if ( jxf_ParCSRMatrixCommPkg(A) )
   {
      jxf_MatvecCommPkgDestroy( jxf_ParCSRMatrixCommPkg(A) );
   }
   jxf_MatvecCommPkgCreate(A);

   jxf_TFree(marker_offd);

   return jxf_error_flag;
}
