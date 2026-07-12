//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_matrix.c -- basic operations for parallel matrices.
 *  Date: 2011/09/07
 */

#include "jxf_mv.h"

/*!
 * \fn jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile
 * \brief Build a parallel matrix by reading data from a given file.
 * \param filename pointer to the matrix file.
 * \param num_functions number of unknown functions
 * \note The storage format of the matrix is expected to be CSR, i.e.
 *         the 1st row:                           n                          // JXF_Int
 *         the 2 to (n+2)-th row:                ia[k]     k=0,1,...,n       // JXF_Int
 *         the (n+3) to (n+nz+2)-th row:         ja[k]     k=0,1,...,nz-1    // JXF_Int
 *         the (n+nz+3) to ((n+2nz+2))-th row:    a[k]     k=0,1,...,nz-1    // JXF_Real
 *       where n is the order of the matrix, nz is the number of non-zeros of the matrix.
 * \author peghoty, Yue Xiaoqiang
 * \date 2009/07/10, 2017/10/31
 */
jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile(char *filename, JXF_Int num_functions, JXF_Int file_base)
{
   jxf_ParCSRMatrix *A = NULL;
   jxf_CSRMatrix *A_CSR = NULL;

   JXF_Int my_id, num_procs;
   JXF_Int i, rest, size, num_nodes, num_dofs;

   JXF_Int *row_part = NULL;
   JXF_Int *col_part = NULL;

   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
   jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   if (my_id == 0)
   {
      if (file_base == 2)
      {
         A_CSR = jxf_CSRMatrixBinaryRead(filename);
      }
      else if (file_base == 12)
      {
         A_CSR = jxf_CSRMatrixBinaryRead_FASP(filename);
      }
      else if (file_base == 5)
      {
         A_CSR = jxf_CSRMatrixRead3(filename);
      }
      else
      {
         A_CSR = jxf_CSRMatrixRead(filename, file_base);
      }
      jxf_CSRMatrixReorder(A_CSR);
   }

   if ((my_id == 0) && (num_functions > 1))
   {
      num_dofs = jxf_CSRMatrixNumRows(A_CSR);
      num_nodes = num_dofs / num_functions;
      if (num_dofs != num_functions * num_nodes)
      {
         row_part = NULL;
         col_part = NULL;
      }
      else
      {
         row_part = jxf_CTAlloc(JXF_Int, num_procs + 1);
         row_part[0] = 0;
         size = num_nodes / num_procs;
         rest = num_nodes - size * num_procs;
         for (i = 0; i < num_procs; i++)
         {
            row_part[i + 1] = row_part[i] + size * num_functions;
            if (i < rest)
               row_part[i + 1] += num_functions;
         }
         col_part = row_part;
      }
   }

   /* ser -> par transfering */
   A = jxf_CSRMatrixToParCSRMatrix(MPI_COMM_WORLD, A_CSR, row_part, col_part);

   if (my_id == 0)
   {
      jxf_CSRMatrixDestroy(A_CSR);
   }

   return (A);
}

jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile5(char *filename, JXF_Int num_functions, JXF_Int file_base)
{
   jxf_ParCSRMatrix *A = NULL;
   jxf_CSRMatrix *A_CSR = NULL;

   JXF_Int my_id, num_procs;
   JXF_Int i, rest, size, num_nodes, num_dofs;

   JXF_Int *row_part = NULL;
   JXF_Int *col_part = NULL;

   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
   jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   if (my_id == 0)
   {
      A_CSR = jxf_CSRMatrixRead2(filename, file_base);
      jxf_CSRMatrixReorder(A_CSR);
   }

   if ((my_id == 0) && (num_functions > 1))
   {
      num_dofs = jxf_CSRMatrixNumRows(A_CSR);
      num_nodes = num_dofs / num_functions;
      if (num_dofs != num_functions * num_nodes)
      {
         row_part = NULL;
         col_part = NULL;
      }
      else
      {
         row_part = jxf_CTAlloc(JXF_Int, num_procs + 1);
         row_part[0] = 0;
         size = num_nodes / num_procs;
         rest = num_nodes - size * num_procs;
         for (i = 0; i < num_procs; i++)
         {
            row_part[i + 1] = row_part[i] + size * num_functions;
            if (i < rest)
               row_part[i + 1] += num_functions;
         }
         col_part = row_part;
      }
   }

   /* ser -> par transfering */
   A = jxf_CSRMatrixToParCSRMatrix(MPI_COMM_WORLD, A_CSR, row_part, col_part);

   if (my_id == 0)
   {
      jxf_CSRMatrixDestroy(A_CSR);
   }

   return (A);
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_ParCSRMatrixCreate
 * \brief Create a parallel CSR Matrix.
 * \note If create is called for JXF_NO_GLOBAL_PARTITION and row_starts and col_starts are NOT
 *       null, then it is assumed that they are array of length 2 containing the start row of
 *       the calling processor followed by the start row of the next processor. - AHB 6/05
 * \date 2009/07/10
 */
jxf_ParCSRMatrix *
jxf_ParCSRMatrixCreate(MPI_Comm comm,
                       JXF_Int global_num_rows,
                       JXF_Int global_num_cols,
                       JXF_Int *row_starts,
                       JXF_Int *col_starts,
                       JXF_Int num_cols_offd,
                       JXF_Int num_nonzeros_diag,
                       JXF_Int num_nonzeros_offd)
{
   jxf_ParCSRMatrix *matrix;

   JXF_Int num_procs, my_id;
   JXF_Int local_num_rows, local_num_cols;
   JXF_Int first_row_index, first_col_diag;

   matrix = jxf_CTAlloc(jxf_ParCSRMatrix, 1);

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   if (!row_starts)
   {
#ifdef JXF_NO_GLOBAL_PARTITION
      jxf_GenerateLocalPartitioning(global_num_rows, num_procs, my_id, &row_starts);
#else
      jxf_GeneratePartitioning(global_num_rows, num_procs, &row_starts);
#endif
   }

   if (!col_starts)
   {
      if (global_num_rows == global_num_cols)
      {
         col_starts = row_starts;
      }
      else
      {
#ifdef JXF_NO_GLOBAL_PARTITION
         jxf_GenerateLocalPartitioning(global_num_cols, num_procs, my_id, &col_starts);
#else
         jxf_GeneratePartitioning(global_num_cols, num_procs, &col_starts);
#endif
      }
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   /* row_starts[0] is start of local rows.  row_starts[1] is start of next processor's rows */
   first_row_index = row_starts[0];
   local_num_rows = row_starts[1] - first_row_index;
   first_col_diag = col_starts[0];
   local_num_cols = col_starts[1] - first_col_diag;
#else
   first_row_index = row_starts[my_id];
   local_num_rows = row_starts[my_id + 1] - first_row_index;
   first_col_diag = col_starts[my_id];
   local_num_cols = col_starts[my_id + 1] - first_col_diag;
#endif

   jxf_ParCSRMatrixComm(matrix) = comm;
   jxf_ParCSRMatrixDiag(matrix) = jxf_CSRMatrixCreate(local_num_rows, local_num_cols, num_nonzeros_diag);
   jxf_ParCSRMatrixOffd(matrix) = jxf_CSRMatrixCreate(local_num_rows, num_cols_offd, num_nonzeros_offd);
   jxf_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   jxf_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;
   jxf_ParCSRMatrixFirstRowIndex(matrix) = first_row_index;
   jxf_ParCSRMatrixFirstColDiag(matrix) = first_col_diag;

   jxf_ParCSRMatrixLastRowIndex(matrix) = first_row_index + local_num_rows - 1;
   jxf_ParCSRMatrixLastColDiag(matrix) = first_col_diag + local_num_cols - 1;

   jxf_ParCSRMatrixColMapOffd(matrix) = NULL;
   jxf_ParCSRMatrixAssumedPartition(matrix) = NULL;

   /*-------------------------------------------------------------------------
    *   When NO_GLOBAL_PARTITION is set we could make these null, instead
    * of leaving the range.  If that change is made, then when this create
    * is called from functions like the matrix-matrix multiply, be careful
    * not to generate a new partition.
    *-------------------------------------------------------------------------*/

   jxf_ParCSRMatrixRowStarts(matrix) = row_starts;
   jxf_ParCSRMatrixColStarts(matrix) = col_starts;

   jxf_ParCSRMatrixDiagT(matrix) = NULL;
   jxf_ParCSRMatrixOffdT(matrix) = NULL;

   jxf_ParCSRMatrixCommPkg(matrix) = NULL;
   jxf_ParCSRMatrixCommPkgT(matrix) = NULL;

   /* set defaults */
   jxf_ParCSRMatrixOwnsData(matrix) = 1;
   jxf_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   jxf_ParCSRMatrixOwnsColStarts(matrix) = 1;
   if (row_starts == col_starts)
   {
      jxf_ParCSRMatrixOwnsColStarts(matrix) = 0;
   }
   jxf_ParCSRMatrixRowindices(matrix) = NULL;
   jxf_ParCSRMatrixRowvalues(matrix) = NULL;
   jxf_ParCSRMatrixGetrowactive(matrix) = 0;

   return matrix;
}

/*!
 * \fn jxh_ParCSRMatrix *jxh_ParCSRMatrixCreate
 * \brief Create a half-precision parallel CSR Matrix.
 */
jxh_ParCSRMatrix *
jxh_ParCSRMatrixCreate(MPI_Comm comm,
                       JXF_Int global_num_rows,
                       JXF_Int global_num_cols,
                       JXF_Int *row_starts,
                       JXF_Int *col_starts,
                       JXF_Int num_cols_offd,
                       JXF_Int num_nonzeros_diag,
                       JXF_Int num_nonzeros_offd)
{
   jxh_ParCSRMatrix *matrix;

   JXF_Int num_procs, my_id;
   JXF_Int local_num_rows, local_num_cols;
   JXF_Int first_row_index, first_col_diag;

   matrix = jxf_CTAlloc(jxh_ParCSRMatrix, 1);

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   if (!row_starts)
   {
#ifdef JXF_NO_GLOBAL_PARTITION
      jxf_GenerateLocalPartitioning(global_num_rows, num_procs, my_id, &row_starts);
#else
      jxf_GeneratePartitioning(global_num_rows, num_procs, &row_starts);
#endif
   }

   if (!col_starts)
   {
      if (global_num_rows == global_num_cols)
      {
         col_starts = row_starts;
      }
      else
      {
#ifdef JXF_NO_GLOBAL_PARTITION
         jxf_GenerateLocalPartitioning(global_num_cols, num_procs, my_id, &col_starts);
#else
         jxf_GeneratePartitioning(global_num_cols, num_procs, &col_starts);
#endif
      }
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   first_row_index = row_starts[0];
   local_num_rows = row_starts[1] - first_row_index;
   first_col_diag = col_starts[0];
   local_num_cols = col_starts[1] - first_col_diag;
#else
   first_row_index = row_starts[my_id];
   local_num_rows = row_starts[my_id + 1] - first_row_index;
   first_col_diag = col_starts[my_id];
   local_num_cols = col_starts[my_id + 1] - first_col_diag;
#endif

   jxh_ParCSRMatrixComm(matrix) = comm;

   jxh_ParCSRMatrixDiag(matrix) = jxh_CSRMatrixCreate(local_num_rows, local_num_cols, num_nonzeros_diag);

   jxh_ParCSRMatrixOffd(matrix) = jxh_CSRMatrixCreate(local_num_rows, num_cols_offd, num_nonzeros_offd);

   jxh_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   jxh_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;

   jxh_ParCSRMatrixFirstRowIndex(matrix) = first_row_index;
   jxh_ParCSRMatrixFirstColDiag(matrix) = first_col_diag;

   jxh_ParCSRMatrixLastRowIndex(matrix) = first_row_index + local_num_rows - 1;
   jxh_ParCSRMatrixLastColDiag(matrix) = first_col_diag + local_num_cols - 1;

   jxh_ParCSRMatrixColMapOffd(matrix) = NULL;
   jxh_ParCSRMatrixAssumedPartition(matrix) = NULL;

   jxh_ParCSRMatrixRowStarts(matrix) = row_starts;
   jxh_ParCSRMatrixColStarts(matrix) = col_starts;

   jxh_ParCSRMatrixDiagT(matrix) = NULL;
   jxh_ParCSRMatrixOffdT(matrix) = NULL;

   jxh_ParCSRMatrixCommPkg(matrix) = NULL;
   jxh_ParCSRMatrixCommPkgT(matrix) = NULL;

   /* set defaults */
   jxh_ParCSRMatrixOwnsData(matrix) = 1;
   jxh_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   jxh_ParCSRMatrixOwnsColStarts(matrix) = 1;

   if (row_starts == col_starts)
   {
      jxh_ParCSRMatrixOwnsColStarts(matrix) = 0;
   }
   jxh_ParCSRMatrixRowindices(matrix) = NULL;
   jxh_ParCSRMatrixRowvalues(matrix) = NULL;
   jxh_ParCSRMatrixGetrowactive(matrix) = 0;

   return matrix;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixInitialize
 * \brief Initialize a parallel CSR Matrix.
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixInitialize(jxf_ParCSRMatrix *matrix)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_CSRMatrixInitialize(jxf_ParCSRMatrixDiag(matrix));
   jxf_CSRMatrixInitialize(jxf_ParCSRMatrixOffd(matrix));
   jxf_ParCSRMatrixColMapOffd(matrix) = jxf_CTAlloc(JXF_Int, jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(matrix)));

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxh_ParCSRMatrixInitialize
 * \brief Initialize a half-precision parallel CSR Matrix.
 */
JXF_Int
jxh_ParCSRMatrixInitialize(jxh_ParCSRMatrix *matrix)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxh_CSRMatrixInitialize(jxh_ParCSRMatrixDiag(matrix));
   jxh_CSRMatrixInitialize(jxh_ParCSRMatrixOffd(matrix));

   jxh_ParCSRMatrixColMapOffd(matrix) = 
      jxf_CTAlloc(JXF_Int,                  jxh_CSRMatrixNumCols(jxh_ParCSRMatrixOffd(matrix)));

   return jxf_error_flag;
}

/*!
 * \fn jxf_CSRMatrix *jxf_ParCSRMatrixToCSRMatrixAll
 * \brief Generates a CSRMatrix from a ParCSRMatrix on all
 *        processors that have parts of the ParCSRMatrix.
 * \date 2009/07/10
 */
jxf_CSRMatrix *
jxf_ParCSRMatrixToCSRMatrixAll(jxf_ParCSRMatrix *par_matrix)
{
   MPI_Comm comm = jxf_ParCSRMatrixComm(par_matrix);
   jxf_CSRMatrix *matrix;
   jxf_CSRMatrix *local_matrix;
   JXF_Int num_rows = jxf_ParCSRMatrixGlobalNumRows(par_matrix);
   JXF_Int num_cols = jxf_ParCSRMatrixGlobalNumCols(par_matrix);
#ifndef JXF_NO_GLOBAL_PARTITION
   JXF_Int *row_starts = jxf_ParCSRMatrixRowStarts(par_matrix);
#endif
   JXF_Int *matrix_i;
   JXF_Int *matrix_j;
   JXF_Real *matrix_data;

   JXF_Int *local_matrix_i;
   JXF_Int *local_matrix_j;
   JXF_Real *local_matrix_data;

   JXF_Int i, j;
   JXF_Int local_num_rows;
   JXF_Int local_num_nonzeros;
   JXF_Int num_nonzeros;
   JXF_Int num_data;
   JXF_Int num_requests;
   JXF_Int vec_len, offset;
   JXF_Int start_index;
   JXF_Int proc_id;
   JXF_Int num_procs, my_id;
   JXF_Int num_types;
   JXF_Int *used_procs;

   MPI_Request *requests;
   MPI_Status *status;

#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int *new_vec_starts;
   JXF_Int num_contacts;
   JXF_Int contact_proc_list[1];
   JXF_Int contact_send_buf[1];
   JXF_Int contact_send_buf_starts[2];
   JXF_Int max_response_size;
   JXF_Int *response_recv_buf = NULL;
   JXF_Int *response_recv_buf_starts = NULL;

   jxf_DataExchangeResponse response_obj;
   jxf_ProcListElements send_proc_obj;

   JXF_Int *send_info = NULL;
   MPI_Status status1;
   JXF_Int count, tag1 = 11112, tag2 = 22223, tag3 = 33334;
   JXF_Int start;
#endif

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION

   local_num_rows = jxf_ParCSRMatrixLastRowIndex(par_matrix) -
                    jxf_ParCSRMatrixFirstRowIndex(par_matrix) + 1;

   local_matrix = jxf_MergeDiagAndOffd(par_matrix); /* creates matrix */
   local_matrix_i = jxf_CSRMatrixI(local_matrix);
   local_matrix_j = jxf_CSRMatrixJ(local_matrix);
   local_matrix_data = jxf_CSRMatrixData(local_matrix);

   /*-----------------------------------------------------------------------------
    * Determine procs that have vector data and store their ids in used_procs
    * we need to do an exchange data for this. If I own row then I will contact
    * processor 0 with the endpoint of my local range。
    *-----------------------------------------------------------------------------*/

   if (local_num_rows > 0)
   {
      num_contacts = 1;
      contact_proc_list[0] = 0;
      contact_send_buf[0] = jxf_ParCSRMatrixLastRowIndex(par_matrix);
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 1;
   }
   else
   {
      num_contacts = 0;
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 0;
   }

   /*--------------------------------------------------------------
    * build the response object
    * send_proc_obj will  be for saving info from contacts
    *------------------------------------------------------------*/
   send_proc_obj.length = 0;
   send_proc_obj.storage_length = 10;
   send_proc_obj.id = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length);
   send_proc_obj.vec_starts = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length + 1);
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = 10;
   send_proc_obj.elements = jxf_CTAlloc(JXF_Int, send_proc_obj.element_storage_length);

   max_response_size = 0; /* each response is null */
   response_obj.fill_response = jxf_FillResponseParToCSRMatrix;
   response_obj.data1 = NULL;
   response_obj.data2 = &send_proc_obj; /* this is where we keep info from contacts */

   jxf_DataExchangeList(num_contacts,
                        contact_proc_list, contact_send_buf,
                        contact_send_buf_starts, sizeof(JXF_Int),
                        sizeof(JXF_Int), &response_obj,
                        max_response_size, 1,
                        comm, (void **)&response_recv_buf,
                        &response_recv_buf_starts);

   /*-------------------------------------------------------------------------------
    * now processor 0 should have a list of ranges for processors that have rows -
    * these are in send_proc_obj - it needs to create the new list of processors
    * and also an array of vec starts - and send to those who own row.
    *-------------------------------------------------------------------------------*/

   if (my_id)
   {
      if (local_num_rows)
      {
         /* look for a message from processor 0 */
         jxf_MPI_Probe(0, tag1, comm, &status1);
         jxf_MPI_Get_count(&status1, JXF_MPI_INT, &count);

         send_info = jxf_CTAlloc(JXF_Int, count);
         jxf_MPI_Recv(send_info, count, JXF_MPI_INT, 0, tag1, comm, &status1);

         /* now unpack */
         num_types = send_info[0];
         used_procs = jxf_CTAlloc(JXF_Int, num_types);
         new_vec_starts = jxf_CTAlloc(JXF_Int, num_types + 1);

         for (i = 1; i <= num_types; i++)
         {
            used_procs[i - 1] = send_info[i];
         }
         for (i = num_types + 1; i < count; i++)
         {
            new_vec_starts[i - num_types - 1] = send_info[i];
         }
      }
      else /* clean up and exit */
      {
         jxf_TFree(send_proc_obj.vec_starts);
         jxf_TFree(send_proc_obj.id);
         jxf_TFree(send_proc_obj.elements);
         if (response_recv_buf)
            jxf_TFree(response_recv_buf);
         if (response_recv_buf_starts)
            jxf_TFree(response_recv_buf_starts);

         if (jxf_CSRMatrixOwnsData(local_matrix))
         {
            jxf_CSRMatrixDestroy(local_matrix);
         }
         else
         {
            jxf_TFree(local_matrix);
         }
         return NULL;
      }
   }
   else /* my_id ==0 */
   {
      num_types = send_proc_obj.length;
      used_procs = jxf_CTAlloc(JXF_Int, num_types);
      new_vec_starts = jxf_CTAlloc(JXF_Int, num_types + 1);

      new_vec_starts[0] = 0;
      for (i = 0; i < num_types; i++)
      {
         used_procs[i] = send_proc_obj.id[i];
         new_vec_starts[i + 1] = send_proc_obj.elements[i] + 1;
      }
      jxf_qsort0(used_procs, 0, num_types - 1);
      jxf_qsort0(new_vec_starts, 0, num_types);

      /*now we need to put into an array to send */
      count = 2 * num_types + 2;
      send_info = jxf_CTAlloc(JXF_Int, count);
      send_info[0] = num_types;
      for (i = 1; i <= num_types; i++)
      {
         send_info[i] = used_procs[i - 1];
      }
      for (i = num_types + 1; i < count; i++)
      {
         send_info[i] = new_vec_starts[i - num_types - 1];
      }
      requests = jxf_CTAlloc(MPI_Request, num_types);
      status = jxf_CTAlloc(MPI_Status, num_types);

      /* don't send to myself - these are sorted so my id would be first */
      start = 0;
      if (used_procs[0] == 0)
      {
         start = 1;
      }

      for (i = start; i < num_types; i++)
      {
         jxf_MPI_Isend(send_info, count, JXF_MPI_INT, used_procs[i], tag1, comm, &requests[i - start]);
      }
      jxf_MPI_Waitall(num_types - start, requests, status);

      jxf_TFree(status);
      jxf_TFree(requests);
   }

   /* clean up */
   jxf_TFree(send_proc_obj.vec_starts);
   jxf_TFree(send_proc_obj.id);
   jxf_TFree(send_proc_obj.elements);
   jxf_TFree(send_info);
   if (response_recv_buf)
      jxf_TFree(response_recv_buf);
   if (response_recv_buf_starts)
      jxf_TFree(response_recv_buf_starts);

   /* now proc 0 can exit if it has no rows */
   if (!local_num_rows)
   {
      if (jxf_CSRMatrixOwnsData(local_matrix))
      {
         jxf_CSRMatrixDestroy(local_matrix);
      }
      else
      {
         jxf_TFree(local_matrix);
      }

      jxf_TFree(new_vec_starts);
      jxf_TFree(used_procs);

      return NULL;
   }

   /* everyone left has rows and knows: new_vec_starts, num_types, and used_procs */

   /* this matrix should be rather small */
   matrix_i = jxf_CTAlloc(JXF_Int, num_rows + 1);

   num_requests = 4 * num_types;
   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   /* exchange contents of local_matrix_i - here we are sending to ourself also */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = new_vec_starts[i + 1] - new_vec_starts[i];
      jxf_MPI_Irecv(&matrix_i[new_vec_starts[i] + 1], vec_len, JXF_MPI_INT,
                    proc_id, tag2, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      jxf_MPI_Isend(&local_matrix_i[1], local_num_rows, JXF_MPI_INT,
                    proc_id, tag2, comm, &requests[j++]);
   }

   jxf_MPI_Waitall(j, requests, status);

   /* generate matrix_i from received data */
   /* global numbering? */
   offset = matrix_i[new_vec_starts[1]];
   for (i = 1; i < num_types; i++)
   {
      for (j = new_vec_starts[i]; j < new_vec_starts[i + 1]; j++)
      {
         matrix_i[j + 1] += offset;
      }
      offset = matrix_i[new_vec_starts[i + 1]];
   }

   num_nonzeros = matrix_i[num_rows];

   matrix = jxf_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   jxf_CSRMatrixI(matrix) = matrix_i;
   jxf_CSRMatrixInitialize(matrix);
   matrix_j = jxf_CSRMatrixJ(matrix);
   matrix_data = jxf_CSRMatrixData(matrix);

   /*--------------------------------------------------------------
    * Generate datatypes for further data exchange and exchange
    * remaining data, i.e. column info and actual data.
    *-------------------------------------------------------------*/

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      start_index = matrix_i[new_vec_starts[i]];
      num_data = matrix_i[new_vec_starts[i + 1]] - start_index;
      jxf_MPI_Irecv(&matrix_data[start_index], num_data, JXF_MPI_REAL,
                    used_procs[i], tag1, comm, &requests[j++]);
      jxf_MPI_Irecv(&matrix_j[start_index], num_data, JXF_MPI_INT,
                    used_procs[i], tag3, comm, &requests[j++]);
   }
   local_num_nonzeros = local_matrix_i[local_num_rows];
   for (i = 0; i < num_types; i++)
   {
      jxf_MPI_Isend(local_matrix_data, local_num_nonzeros, JXF_MPI_REAL,
                    used_procs[i], tag1, comm, &requests[j++]);
      jxf_MPI_Isend(local_matrix_j, local_num_nonzeros, JXF_MPI_INT,
                    used_procs[i], tag3, comm, &requests[j++]);
   }

   jxf_MPI_Waitall(num_requests, requests, status);

   jxf_TFree(new_vec_starts);

#else

   local_num_rows = row_starts[my_id + 1] - row_starts[my_id];

   /* if my_id contains no data, return NULL */

   if (!local_num_rows)
   {
      return NULL;
   }

   local_matrix = jxf_MergeDiagAndOffd(par_matrix);
   local_matrix_i = jxf_CSRMatrixI(local_matrix);
   local_matrix_j = jxf_CSRMatrixJ(local_matrix);
   local_matrix_data = jxf_CSRMatrixData(local_matrix);

   matrix_i = jxf_CTAlloc(JXF_Int, num_rows + 1);

   /* determine procs that have vector data and store their ids in used_procs */

   num_types = 0;
   for (i = 0; i < num_procs; i++)
   {
      if (row_starts[i + 1] - row_starts[i] && i - my_id)
      {
         num_types++;
      }
   }
   num_requests = 4 * num_types;

   used_procs = jxf_CTAlloc(JXF_Int, num_types);
   j = 0;
   for (i = 0; i < num_procs; i++)
   {
      if (row_starts[i + 1] - row_starts[i] && i - my_id)
      {
         used_procs[j++] = i;
      }
   }

   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   /* exchange contents of local_matrix_i */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = row_starts[proc_id + 1] - row_starts[proc_id];
      jxf_MPI_Irecv(&matrix_i[row_starts[proc_id] + 1], vec_len, JXF_MPI_INT,
                    proc_id, 0, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      jxf_MPI_Isend(&local_matrix_i[1], local_num_rows, JXF_MPI_INT,
                    proc_id, 0, comm, &requests[j++]);
   }

   vec_len = row_starts[my_id + 1] - row_starts[my_id];
   for (i = 1; i <= vec_len; i++)
   {
      matrix_i[row_starts[my_id] + i] = local_matrix_i[i];
   }

   jxf_MPI_Waitall(j, requests, status);

   /* generate matrix_i from received data */

   offset = matrix_i[row_starts[1]];
   for (i = 1; i < num_procs; i++)
   {
      for (j = row_starts[i]; j < row_starts[i + 1]; j++)
      {
         matrix_i[j + 1] += offset;
      }
      offset = matrix_i[row_starts[i + 1]];
   }

   num_nonzeros = matrix_i[num_rows];

   matrix = jxf_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   jxf_CSRMatrixI(matrix) = matrix_i;
   jxf_CSRMatrixInitialize(matrix);
   matrix_j = jxf_CSRMatrixJ(matrix);
   matrix_data = jxf_CSRMatrixData(matrix);

   /*-------------------------------------------------------------
    * Generate datatypes for further data exchange and exchange
    * remaining data, i.e. column info and actual data
    *------------------------------------------------------------*/

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      start_index = matrix_i[row_starts[proc_id]];
      num_data = matrix_i[row_starts[proc_id + 1]] - start_index;
      jxf_MPI_Irecv(&matrix_data[start_index], num_data, JXF_MPI_REAL,
                    used_procs[i], 0, comm, &requests[j++]);
      jxf_MPI_Irecv(&matrix_j[start_index], num_data, JXF_MPI_INT,
                    used_procs[i], 0, comm, &requests[j++]);
   }
   local_num_nonzeros = local_matrix_i[local_num_rows];
   for (i = 0; i < num_types; i++)
   {
      jxf_MPI_Isend(local_matrix_data, local_num_nonzeros, JXF_MPI_REAL,
                    used_procs[i], 0, comm, &requests[j++]);
      jxf_MPI_Isend(local_matrix_j, local_num_nonzeros, JXF_MPI_INT,
                    used_procs[i], 0, comm, &requests[j++]);
   }

   start_index = matrix_i[row_starts[my_id]];
   for (i = 0; i < local_num_nonzeros; i++)
   {
      matrix_j[start_index + i] = local_matrix_j[i];
      matrix_data[start_index + i] = local_matrix_data[i];
   }
   jxf_MPI_Waitall(num_requests, requests, status);

   start_index = matrix_i[row_starts[my_id]];
   for (i = 0; i < local_num_nonzeros; i++)
   {
      matrix_j[start_index + i] = local_matrix_j[i];
      matrix_data[start_index + i] = local_matrix_data[i];
   }
   jxf_MPI_Waitall(num_requests, requests, status);

#endif

   if (jxf_CSRMatrixOwnsData(local_matrix))
   {
      jxf_CSRMatrixDestroy(local_matrix);
   }
   else
   {
      jxf_TFree(local_matrix);
   }

   if (num_requests)
   {
      jxf_TFree(requests);
      jxf_TFree(status);
      jxf_TFree(used_procs);
   }

   return matrix;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix
 * \brief Generates a ParCSRMatrix distributed across the
 *        processors in comm from a CSRMatrix on proc 0.
 * \note This shouldn't be used with the JXF_NO_GLOBAL_PARTITON option!
 * \date 2009/07/10
 */
jxf_ParCSRMatrix *
jxf_CSRMatrixToParCSRMatrix(MPI_Comm comm,
                            jxf_CSRMatrix *A,
                            JXF_Int *row_starts,
                            JXF_Int *col_starts)
{
   JXF_Int *global_data;
   JXF_Int global_size;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int *local_num_rows;

   JXF_Int num_procs, my_id;
   JXF_Int *local_num_nonzeros = NULL;
   JXF_Int num_nonzeros;

   JXF_Real *a_data = NULL;
   JXF_Int *a_i = NULL;
   JXF_Int *a_j = NULL;

   jxf_CSRMatrix *local_A;

   MPI_Request *requests;
   MPI_Status *status, status0;
   MPI_Datatype *csr_matrix_datatypes;

   jxf_ParCSRMatrix *par_matrix;

   JXF_Int first_col_diag;
   JXF_Int last_col_diag;

   JXF_Int i, j, ind;

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   global_data = jxf_CTAlloc(JXF_Int, 2 * num_procs + 6);
   if (my_id == 0)
   {
      global_size = 3;
      if (row_starts)
      {
         if (col_starts)
         {
            if (col_starts != row_starts)
            {
               /*----------------------------------------------------------
                * contains code for what to expect,
                * if 0: row_starts = col_starts, only row_starts given
                * if 1: only row_starts given, col_starts = NULL
                * if 2: both row_starts and col_starts given
                * if 3: only col_starts given, row_starts = NULL
                *---------------------------------------------------------*/
               global_data[3] = 2;
               global_size = 2 * num_procs + 6;
               for (i = 0; i < num_procs + 1; i++)
               {
                  global_data[i + 4] = row_starts[i];
               }
               for (i = 0; i < num_procs + 1; i++)
               {
                  global_data[i + num_procs + 5] = col_starts[i];
               }
            }
            else
            {
               global_data[3] = 0;
               global_size = num_procs + 5;
               for (i = 0; i < num_procs + 1; i++)
               {
                  global_data[i + 4] = row_starts[i];
               }
            }
         }
         else
         {
            global_data[3] = 1;
            global_size = num_procs + 5;
            for (i = 0; i < num_procs + 1; i++)
            {
               global_data[i + 4] = row_starts[i];
            }
         }
      }
      else
      {
         if (col_starts)
         {
            global_data[3] = 3;
            global_size = num_procs + 5;
            for (i = 0; i < num_procs + 1; i++)
            {
               global_data[i + 4] = col_starts[i];
            }
         }
      }
      global_data[0] = jxf_CSRMatrixNumRows(A);
      global_data[1] = jxf_CSRMatrixNumCols(A);
      global_data[2] = global_size;
      a_data = jxf_CSRMatrixData(A);
      a_i = jxf_CSRMatrixI(A);
      a_j = jxf_CSRMatrixJ(A);
   }
   jxf_MPI_Bcast(global_data, 3, JXF_MPI_INT, 0, comm);
   global_num_rows = global_data[0];
   global_num_cols = global_data[1];

   global_size = global_data[2];
   if (global_size > 3)
   {
      jxf_MPI_Bcast(&global_data[3], global_size - 3, JXF_MPI_INT, 0, comm);
      if (my_id > 0)
      {
         if (global_data[3] < 3)
         {
            row_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
            for (i = 0; i < num_procs + 1; i++)
            {
               row_starts[i] = global_data[i + 4];
            }
            if (global_data[3] == 0)
            {
               col_starts = row_starts;
            }
            if (global_data[3] == 2)
            {
               col_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
               for (i = 0; i < num_procs + 1; i++)
               {
                  col_starts[i] = global_data[i + num_procs + 5];
               }
            }
         }
         else
         {
            col_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
            for (i = 0; i < num_procs + 1; i++)
            {
               col_starts[i] = global_data[i + 4];
            }
         }
      }
   }
   jxf_TFree(global_data);

   local_num_rows = jxf_CTAlloc(JXF_Int, num_procs);
   csr_matrix_datatypes = jxf_CTAlloc(MPI_Datatype, num_procs);

   par_matrix = jxf_ParCSRMatrixCreate(comm, global_num_rows, global_num_cols, row_starts, col_starts, 0, 0, 0);

   row_starts = jxf_ParCSRMatrixRowStarts(par_matrix);
   col_starts = jxf_ParCSRMatrixColStarts(par_matrix);

   for (i = 0; i < num_procs; i++)
   {
      local_num_rows[i] = row_starts[i + 1] - row_starts[i];
   }

   if (my_id == 0)
   {
      local_num_nonzeros = jxf_CTAlloc(JXF_Int, num_procs);
      for (i = 0; i < num_procs - 1; i++)
      {
         local_num_nonzeros[i] = a_i[row_starts[i + 1]] - a_i[row_starts[i]];
      }
      local_num_nonzeros[num_procs - 1] = a_i[global_num_rows] - a_i[row_starts[num_procs - 1]];
   }
   jxf_MPI_Scatter(local_num_nonzeros, 1, JXF_MPI_INT, &num_nonzeros, 1, JXF_MPI_INT, 0, comm);

   if (my_id == 0)
   {
      num_nonzeros = local_num_nonzeros[0];
   }

   local_A = jxf_CSRMatrixCreate(local_num_rows[my_id], global_num_cols, num_nonzeros);

   if (my_id == 0)
   {
      requests = jxf_CTAlloc(MPI_Request, num_procs - 1);
      status = jxf_CTAlloc(MPI_Status, num_procs - 1);
      j = 0;
      for (i = 1; i < num_procs; i++)
      {
         ind = a_i[row_starts[i]];
         jxf_BuildCSRMatrixMPIDataType(local_num_nonzeros[i],
                                       local_num_rows[i],
                                       &a_data[ind],
                                       &a_i[row_starts[i]],
                                       &a_j[ind],
                                       &csr_matrix_datatypes[i]);
         jxf_MPI_Isend(MPI_BOTTOM, 1, csr_matrix_datatypes[i], i, 0, comm, &requests[j++]);
         jxf_MPI_Type_free(&csr_matrix_datatypes[i]);
      }
      jxf_CSRMatrixData(local_A) = a_data;
      jxf_CSRMatrixI(local_A) = a_i;
      jxf_CSRMatrixJ(local_A) = a_j;
      jxf_CSRMatrixOwnsData(local_A) = 0;
      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
      jxf_TFree(local_num_nonzeros);
   }
   else
   {
      jxf_CSRMatrixInitialize(local_A);
      jxf_BuildCSRMatrixMPIDataType(num_nonzeros,
                                    local_num_rows[my_id],
                                    jxf_CSRMatrixData(local_A),
                                    jxf_CSRMatrixI(local_A),
                                    jxf_CSRMatrixJ(local_A),
                                    csr_matrix_datatypes);
      jxf_MPI_Recv(MPI_BOTTOM, 1, csr_matrix_datatypes[0], 0, 0, comm, &status0);
      jxf_MPI_Type_free(csr_matrix_datatypes);
   }

   first_col_diag = col_starts[my_id];
   last_col_diag = col_starts[my_id + 1] - 1;

   jxf_GenerateDiagAndOffd(local_A, par_matrix, first_col_diag, last_col_diag);

   /* set pointers back to NULL before destroying */
   if (my_id == 0)
   {
      jxf_CSRMatrixData(local_A) = NULL;
      jxf_CSRMatrixI(local_A) = NULL;
      jxf_CSRMatrixJ(local_A) = NULL;
   }
   jxf_CSRMatrixDestroy(local_A);
   jxf_TFree(local_num_rows);
   jxf_TFree(csr_matrix_datatypes);

   return par_matrix;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix_FromGivenPro
 * \brief Generates a ParCSRMatrix distributed across the
 *        processors in comm from a CSRMatrix on proc srcid.
 * \note This shouldn't be used with the JXF_NO_GLOBAL_PARTITON option!
 * \note 在jxf_CSRMatrixToParCSRMatrix的基础上改写的程序还有问题，这里重新写.
 * \author peghoty
 * \date 2012/03/03
 */
jxf_ParCSRMatrix *
jxf_CSRMatrixToParCSRMatrix_FromGivenPro(MPI_Comm comm,
                                         JXF_Int srcid,
                                         jxf_CSRMatrix *A,
                                         JXF_Int *row_starts,
                                         JXF_Int *col_starts)
{
   JXF_Int *global_data = NULL;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int *local_num_rows = NULL;
   JXF_Int *local_num_nonzeros = NULL;

   JXF_Int num_procs, my_id;
   JXF_Int num_rows;
   JXF_Int num_nonzeros;

   JXF_Real *a_data = NULL;
   JXF_Int *a_i = NULL;
   JXF_Int *a_j = NULL;

   jxf_CSRMatrix *local_A = NULL;

   MPI_Request *requests;
   MPI_Status *status;
   MPI_Status status0;

   jxf_ParCSRMatrix *par_matrix = NULL;

   JXF_Int first_col_diag;
   JXF_Int last_col_diag;

   JXF_Int i, j, comsize;

   JXF_Int start_ia, start_ja;
   JXF_Int tag_ia = 111;
   JXF_Int tag_ja = 222;
   JXF_Int tag_aa = 333;

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   comsize = 3 * (num_procs - 1);

   global_data = jxf_CTAlloc(JXF_Int, 2);
   if (my_id == srcid)
   {
      global_data[0] = jxf_CSRMatrixNumRows(A);
      global_data[1] = jxf_CSRMatrixNumCols(A);
      a_data = jxf_CSRMatrixData(A);
      a_i = jxf_CSRMatrixI(A);
      a_j = jxf_CSRMatrixJ(A);
   }
   jxf_MPI_Bcast(global_data, 2, JXF_MPI_INT, srcid, comm);
   global_num_rows = global_data[0];
   global_num_cols = global_data[1];
   jxf_TFree(global_data);

   local_num_rows = jxf_CTAlloc(JXF_Int, num_procs);
   local_num_nonzeros = jxf_CTAlloc(JXF_Int, num_procs);

   par_matrix = jxf_ParCSRMatrixCreate(comm, global_num_rows, global_num_cols, row_starts, col_starts, 0, 0, 0);

   row_starts = jxf_ParCSRMatrixRowStarts(par_matrix);
   col_starts = jxf_ParCSRMatrixColStarts(par_matrix);

   for (i = 0; i < num_procs; i++)
   {
      local_num_rows[i] = row_starts[i + 1] - row_starts[i];
   }

   if (my_id == srcid)
   {
      for (i = 0; i < num_procs; i++)
      {
         local_num_nonzeros[i] = a_i[row_starts[i + 1]] - a_i[row_starts[i]];
      }
   }
   jxf_MPI_Bcast(local_num_nonzeros, num_procs, JXF_MPI_INT, srcid, comm);

   num_nonzeros = local_num_nonzeros[my_id];
   num_rows = local_num_rows[my_id];

   local_A = jxf_CSRMatrixCreate(num_rows, global_num_cols, num_nonzeros);
   jxf_CSRMatrixInitialize(local_A);

   if (my_id == srcid)
   {
      requests = jxf_CTAlloc(MPI_Request, comsize);
      status = jxf_CTAlloc(MPI_Status, comsize);
      j = 0;
      for (i = 0; i < num_procs; i++)
      {
         if (i != srcid)
         {
            start_ia = row_starts[i];
            start_ja = a_i[start_ia];
            jxf_MPI_Isend(&a_i[start_ia], local_num_rows[i] + 1, JXF_MPI_INT, i, tag_ia, comm, &requests[j++]);
            jxf_MPI_Isend(&a_j[start_ja], local_num_nonzeros[i], JXF_MPI_INT, i, tag_ja, comm, &requests[j++]);
            jxf_MPI_Isend(&a_data[start_ja], local_num_nonzeros[i], JXF_MPI_REAL, i, tag_aa, comm, &requests[j++]);
         }
      }

      start_ia = row_starts[my_id];
      start_ja = a_i[start_ia];
      for (i = 0; i <= num_rows; i++)
      {
         jxf_CSRMatrixI(local_A)[i] = a_i[start_ia + i];
      }
      for (i = 0; i < num_nonzeros; i++)
      {
         jxf_CSRMatrixData(local_A)[i] = a_data[start_ja + i];
         jxf_CSRMatrixJ(local_A)[i] = a_j[start_ja + i];
      }

      jxf_MPI_Waitall(comsize, requests, status);

      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      jxf_MPI_Recv(jxf_CSRMatrixI(local_A), num_rows + 1, JXF_MPI_INT, srcid, tag_ia, comm, &status0);
      jxf_MPI_Recv(jxf_CSRMatrixJ(local_A), num_nonzeros, JXF_MPI_INT, srcid, tag_ja, comm, &status0);
      jxf_MPI_Recv(jxf_CSRMatrixData(local_A), num_nonzeros, JXF_MPI_REAL, srcid, tag_aa, comm, &status0);
   }

   first_col_diag = col_starts[my_id];
   last_col_diag = col_starts[my_id + 1] - 1;
   jxf_GenerateDiagAndOffd(local_A, par_matrix, first_col_diag, last_col_diag);

   jxf_CSRMatrixDestroy(local_A);
   jxf_TFree(local_num_rows);
   jxf_TFree(local_num_nonzeros);

   return (par_matrix);
}

/*!
 * \fn JXF_Int jxf_CSRMatrixToParCSRMatrix_sp
 * \brief Transfer a serial CSR matrix to a parallel one (only for single-processor case).
 * \author peghoty
 * \date 2011/09/27
 */
jxf_ParCSRMatrix *
jxf_CSRMatrixToParCSRMatrix_sp(MPI_Comm comm, jxf_CSRMatrix *A)
{
   jxf_ParCSRMatrix *par_matrix = NULL;

   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   JXF_Int num_nonzeros = jxf_CSRMatrixNumNonzeros(A);

   /* Create a paralle matrix */
   par_matrix = jxf_ParCSRMatrixCreate(comm, num_rows, num_cols, NULL, NULL, 0, num_nonzeros, 0);

   /* deal with the diag */
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(par_matrix);
   jxf_CSRMatrixData(diag) = jxf_CSRMatrixData(A);
   jxf_CSRMatrixI(diag) = jxf_CSRMatrixI(A);
   jxf_CSRMatrixJ(diag) = jxf_CSRMatrixJ(A);

   /* deal with the offd */
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(par_matrix);
   JXF_Int *offd_i = jxf_CTAlloc(JXF_Int, num_rows + 1);
   JXF_Int i;
   for (i = 0; i < num_rows + 1; i++)
   {
      offd_i[i] = 0;
   }
   jxf_CSRMatrixNumCols(offd) = 0;
   jxf_CSRMatrixI(offd) = offd_i;

   return (par_matrix);
}

/*!
 * \fn JXF_Int jxf_GenerateDiagAndOffd
 * \brief Generate the Diag and Offd for a parallel CSR matrix.
 * \date 2011/09/05
 */
JXF_Int
jxf_GenerateDiagAndOffd(jxf_CSRMatrix *A,
                        jxf_ParCSRMatrix *matrix,
                        JXF_Int first_col_diag,
                        JXF_Int last_col_diag)
{
   JXF_Int i, j;
   JXF_Int jo, jd;
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   JXF_Real *a_data = jxf_CSRMatrixData(A);
   JXF_Int *a_i = jxf_CSRMatrixI(A);
   JXF_Int *a_j = jxf_CSRMatrixJ(A);

   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(matrix);
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(matrix);

   JXF_Int *col_map_offd;

   JXF_Real *diag_data, *offd_data;
   JXF_Int *diag_i, *offd_i;
   JXF_Int *diag_j, *offd_j;
   JXF_Int *marker;
   JXF_Int num_cols_diag, num_cols_offd;
   JXF_Int first_elmt = a_i[0];
   JXF_Int num_nonzeros = a_i[num_rows] - first_elmt;
   JXF_Int counter;

   num_cols_diag = last_col_diag - first_col_diag + 1;
   num_cols_offd = 0;

   if (num_cols - num_cols_diag)
   {
      jxf_CSRMatrixInitialize(diag);
      diag_i = jxf_CSRMatrixI(diag);

      jxf_CSRMatrixInitialize(offd);
      offd_i = jxf_CSRMatrixI(offd);
      marker = jxf_CTAlloc(JXF_Int, num_cols);

      for (i = 0; i < num_cols; i++)
      {
         marker[i] = 0;
      }

      jo = 0;
      jd = 0;
      for (i = 0; i < num_rows; i++)
      {
         offd_i[i] = jo;
         diag_i[i] = jd;

         for (j = a_i[i] - first_elmt; j < a_i[i + 1] - first_elmt; j++)
         {
            if (a_j[j] < first_col_diag || a_j[j] > last_col_diag)
            {
               if (!marker[a_j[j]])
               {
                  marker[a_j[j]] = 1;
                  num_cols_offd++;
               }
               jo++;
            }
            else
            {
               jd++;
            }
         }
      }
      offd_i[num_rows] = jo;
      diag_i[num_rows] = jd;

      jxf_ParCSRMatrixColMapOffd(matrix) = jxf_CTAlloc(JXF_Int, num_cols_offd);
      col_map_offd = jxf_ParCSRMatrixColMapOffd(matrix);

      counter = 0;
      for (i = 0; i < num_cols; i++)
      {
         if (marker[i])
         {
            col_map_offd[counter] = i;
            marker[i] = counter;
            counter++;
         }
      }
      jxf_CSRMatrixNumNonzeros(diag) = jd;
      jxf_CSRMatrixInitialize(diag);
      diag_data = jxf_CSRMatrixData(diag);
      diag_j = jxf_CSRMatrixJ(diag);

      jxf_CSRMatrixNumNonzeros(offd) = jo;
      jxf_CSRMatrixNumCols(offd) = num_cols_offd;
      jxf_CSRMatrixInitialize(offd);
      offd_data = jxf_CSRMatrixData(offd);
      offd_j = jxf_CSRMatrixJ(offd);

      jo = 0;
      jd = 0;
      for (i = 0; i < num_rows; i++)
      {
         for (j = a_i[i] - first_elmt; j < a_i[i + 1] - first_elmt; j++)
         {
            if (a_j[j] < first_col_diag || a_j[j] > last_col_diag)
            {
               offd_data[jo] = a_data[j];
               offd_j[jo++] = marker[a_j[j]];
            }
            else
            {
               diag_data[jd] = a_data[j];
               diag_j[jd++] = a_j[j] - first_col_diag;
            }
         }
      }
      jxf_TFree(marker);
   }
   else
   {
      jxf_CSRMatrixNumNonzeros(diag) = num_nonzeros;
      jxf_CSRMatrixInitialize(diag);
      diag_data = jxf_CSRMatrixData(diag);
      diag_i = jxf_CSRMatrixI(diag);
      diag_j = jxf_CSRMatrixJ(diag);

      for (i = 0; i < num_nonzeros; i++)
      {
         diag_data[i] = a_data[i];
         diag_j[i] = a_j[i];
      }
      offd_i = jxf_CTAlloc(JXF_Int, num_rows + 1);

      for (i = 0; i < num_rows + 1; i++)
      {
         diag_i[i] = a_i[i];
         offd_i[i] = 0;
      }

      jxf_CSRMatrixNumCols(offd) = 0;
      jxf_CSRMatrixI(offd) = offd_i;
   }

   return jxf_error_flag;
}

/*!
 * \fn jxf_CSRMatrix *jxf_MergeDiagAndOffd
 * \brief Merge the Diag and Offd for a parallel CSR matrix.
 * \date 2011/09/08
 */
jxf_CSRMatrix *
jxf_MergeDiagAndOffd(jxf_ParCSRMatrix *par_matrix)
{
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(par_matrix);
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(par_matrix);
   jxf_CSRMatrix *matrix;

   JXF_Int num_cols = jxf_ParCSRMatrixGlobalNumCols(par_matrix);
   JXF_Int first_col_diag = jxf_ParCSRMatrixFirstColDiag(par_matrix);
   JXF_Int *col_map_offd = jxf_ParCSRMatrixColMapOffd(par_matrix);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(diag);

   JXF_Int *diag_i = jxf_CSRMatrixI(diag);
   JXF_Int *diag_j = jxf_CSRMatrixJ(diag);
   JXF_Real *diag_data = jxf_CSRMatrixData(diag);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);

   JXF_Int *matrix_i;
   JXF_Int *matrix_j;
   JXF_Real *matrix_data;

   JXF_Int num_nonzeros, i, j;
   JXF_Int count;

   num_nonzeros = diag_i[num_rows] + offd_i[num_rows];

   matrix = jxf_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   jxf_CSRMatrixInitialize(matrix);

   matrix_i = jxf_CSRMatrixI(matrix);
   matrix_j = jxf_CSRMatrixJ(matrix);
   matrix_data = jxf_CSRMatrixData(matrix);

   count = 0;
   matrix_i[0] = 0;
   for (i = 0; i < num_rows; i++)
   {
      for (j = diag_i[i]; j < diag_i[i + 1]; j++)
      {
         matrix_data[count] = diag_data[j];
         matrix_j[count++] = diag_j[j] + first_col_diag;
      }
      for (j = offd_i[i]; j < offd_i[i + 1]; j++)
      {
         matrix_data[count] = offd_data[j];
         matrix_j[count++] = col_map_offd[offd_j[j]];
      }
      matrix_i[i + 1] = count;
   }

   return matrix;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_ParCSRMatrixRead
 * \brief Read a parallel matrix from given file(s).
 * \date 2011/09/08
 */
jxf_ParCSRMatrix *
jxf_ParCSRMatrixRead(MPI_Comm comm, const char *file_name, JXF_Int file_base)
{
   jxf_ParCSRMatrix *matrix;
   jxf_CSRMatrix *diag;
   jxf_CSRMatrix *offd;
   JXF_Int my_id, i, num_procs;
   char new_file_d[80], new_file_o[80], new_file_info[80];
   JXF_Int global_num_rows, global_num_cols, num_cols_offd;
   JXF_Int local_num_rows;
   JXF_Int *row_starts;
   JXF_Int *col_starts;
   JXF_Int *col_map_offd;
   FILE *fp;
   JXF_Int equal = 1;

#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int row_s, row_e, col_s, col_e;
#endif

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

#ifdef JXF_NO_GLOBAL_PARTITION
   row_starts = jxf_CTAlloc(JXF_Int, 2);
   col_starts = jxf_CTAlloc(JXF_Int, 2);
#else
   row_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
   col_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
#endif

   jxf_sprintf(new_file_d, "%s.D.%d", file_name, my_id);
   jxf_sprintf(new_file_o, "%s.O.%d", file_name, my_id);
   jxf_sprintf(new_file_info, "%s.INFO.%d", file_name, my_id);
   fp = fopen(new_file_info, "r");
   jxf_fscanf(fp, "%d", &global_num_rows);
   jxf_fscanf(fp, "%d", &global_num_cols);
   jxf_fscanf(fp, "%d", &num_cols_offd);

#ifdef JXF_NO_GLOBAL_PARTITION
   /* the bgl input file should only contain the EXACT range for local processor */
   jxf_fscanf(fp, "%d %d %d %d", &row_s, &row_e, &col_s, &col_e);
   row_starts[0] = row_s;
   row_starts[1] = row_e;
   col_starts[0] = col_s;
   col_starts[1] = col_e;
#else
   for (i = 0; i < num_procs; i++)
   {
      jxf_fscanf(fp, "%d %d", &row_starts[i], &col_starts[i]);
   }
   row_starts[num_procs] = global_num_rows;
   col_starts[num_procs] = global_num_cols;
#endif

   col_map_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);

   for (i = 0; i < num_cols_offd; i++)
   {
      jxf_fscanf(fp, "%d", &col_map_offd[i]);
   }

   fclose(fp);

#ifdef JXF_NO_GLOBAL_PARTITION
   for (i = 1; i >= 0; i--)
      if (row_starts[i] != col_starts[i])
      {
         equal = 0;
         break;
      }
#else
   for (i = num_procs; i >= 0; i--)
   {
      if (row_starts[i] != col_starts[i])
      {
         equal = 0;
         break;
      }
   }
#endif
   if (equal)
   {
      jxf_TFree(col_starts);
      col_starts = row_starts;
   }

   diag = jxf_CSRMatrixRead(new_file_d, file_base);
   local_num_rows = jxf_CSRMatrixNumRows(diag);

   if (num_cols_offd)
   {
      offd = jxf_CSRMatrixRead(new_file_o, file_base);
   }
   else
   {
      offd = jxf_CSRMatrixCreate(local_num_rows, 0, 0);
      jxf_CSRMatrixInitialize(offd);
   }

   matrix = jxf_CTAlloc(jxf_ParCSRMatrix, 1);

   jxf_ParCSRMatrixComm(matrix) = comm;
   jxf_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   jxf_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;

#ifdef JXF_NO_GLOBAL_PARTITION
   jxf_ParCSRMatrixFirstRowIndex(matrix) = row_s;
   jxf_ParCSRMatrixFirstColDiag(matrix) = col_s;
   jxf_ParCSRMatrixLastRowIndex(matrix) = row_e - 1;
   jxf_ParCSRMatrixLastColDiag(matrix) = col_e - 1;
#else
   jxf_ParCSRMatrixFirstRowIndex(matrix) = row_starts[my_id];
   jxf_ParCSRMatrixFirstColDiag(matrix) = col_starts[my_id];
   jxf_ParCSRMatrixLastRowIndex(matrix) = row_starts[my_id + 1] - 1;
   jxf_ParCSRMatrixLastColDiag(matrix) = col_starts[my_id + 1] - 1;
#endif

   jxf_ParCSRMatrixRowStarts(matrix) = row_starts;
   jxf_ParCSRMatrixColStarts(matrix) = col_starts;
   jxf_ParCSRMatrixCommPkg(matrix) = NULL;

   /* set defaults */
   jxf_ParCSRMatrixOwnsData(matrix) = 1;
   jxf_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   jxf_ParCSRMatrixOwnsColStarts(matrix) = 1;
   if (row_starts == col_starts)
   {
      jxf_ParCSRMatrixOwnsColStarts(matrix) = 0;
   }

   jxf_ParCSRMatrixDiag(matrix) = diag;
   jxf_ParCSRMatrixOffd(matrix) = offd;
   if (num_cols_offd)
   {
      jxf_ParCSRMatrixColMapOffd(matrix) = col_map_offd;
   }
   else
   {
      jxf_ParCSRMatrixColMapOffd(matrix) = NULL;
   }

   return matrix;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixPrint
 * \brief Print a parallel matrix into given file(s).
 * \date 2011/09/08
 */
JXF_Int
jxf_ParCSRMatrixPrint(jxf_ParCSRMatrix *matrix, const char *file_name)
{
   MPI_Comm comm;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int *col_map_offd;
#ifndef JXF_NO_GLOBAL_PARTITION
   JXF_Int *row_starts;
   JXF_Int *col_starts;
#endif
   JXF_Int my_id, i, num_procs;
   char new_file_d[80], new_file_o[80], new_file_info[80];
   FILE *fp;
   JXF_Int num_cols_offd = 0;
#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int row_s, row_e, col_s, col_e;
#endif
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   comm = jxf_ParCSRMatrixComm(matrix);
   global_num_rows = jxf_ParCSRMatrixGlobalNumRows(matrix);
   global_num_cols = jxf_ParCSRMatrixGlobalNumCols(matrix);
   col_map_offd = jxf_ParCSRMatrixColMapOffd(matrix);

#ifndef JXF_NO_GLOBAL_PARTITION
   row_starts = jxf_ParCSRMatrixRowStarts(matrix);
   col_starts = jxf_ParCSRMatrixColStarts(matrix);
#endif

   if (jxf_ParCSRMatrixOffd(matrix))
   {
      num_cols_offd = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(matrix));
   }

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   jxf_sprintf(new_file_d, "%s.D.%d", file_name, my_id);
   jxf_sprintf(new_file_o, "%s.O.%d", file_name, my_id);
   jxf_sprintf(new_file_info, "%s.INFO.%d", file_name, my_id);
   jxf_CSRMatrixPrint(jxf_ParCSRMatrixDiag(matrix), new_file_d);
   if (num_cols_offd != 0)
   {
      jxf_CSRMatrixPrint(jxf_ParCSRMatrixOffd(matrix), new_file_o);
   }

   fp = fopen(new_file_info, "w");
   jxf_fprintf(fp, "%d\n", global_num_rows);
   jxf_fprintf(fp, "%d\n", global_num_cols);
   jxf_fprintf(fp, "%d\n", num_cols_offd);

#ifdef JXF_NO_GLOBAL_PARTITION
   row_s = jxf_ParCSRMatrixFirstRowIndex(matrix);
   row_e = jxf_ParCSRMatrixLastRowIndex(matrix);
   col_s = jxf_ParCSRMatrixFirstColDiag(matrix);
   col_e = jxf_ParCSRMatrixLastColDiag(matrix);
   /* add 1 to the ends because this is a starts partition */
   jxf_fprintf(fp, "%d %d %d %d\n", row_s, row_e + 1, col_s, col_e + 1);
#else
   for (i = 0; i < num_procs; i++)
   {
      jxf_fprintf(fp, "%d %d\n", row_starts[i], col_starts[i]);
   }
#endif

   for (i = 0; i < num_cols_offd; i++)
   {
      jxf_fprintf(fp, "%d\n", col_map_offd[i]);
   }
   fclose(fp);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixDestroy
 * \brief Destroy a parallel CSR matrix.
 * \date 2011/09/08
 */
JXF_Int
jxf_ParCSRMatrixDestroy(jxf_ParCSRMatrix *matrix)
{
   if (matrix)
   {
      if (jxf_ParCSRMatrixOwnsData(matrix))
      {
         jxf_CSRMatrixDestroy(jxf_ParCSRMatrixDiag(matrix));
         jxf_CSRMatrixDestroy(jxf_ParCSRMatrixOffd(matrix));
         if (jxf_ParCSRMatrixDiagT(matrix))
         {
            jxf_CSRMatrixDestroy(jxf_ParCSRMatrixDiagT(matrix));
         }
         if (jxf_ParCSRMatrixOffdT(matrix))
         {
            jxf_CSRMatrixDestroy(jxf_ParCSRMatrixOffdT(matrix));
         }
         if (jxf_ParCSRMatrixColMapOffd(matrix))
         {
            jxf_TFree(jxf_ParCSRMatrixColMapOffd(matrix));
         }
         if (jxf_ParCSRMatrixCommPkg(matrix))
         {
            jxf_MatvecCommPkgDestroy(jxf_ParCSRMatrixCommPkg(matrix));
         }
         if (jxf_ParCSRMatrixCommPkgT(matrix))
         {
            jxf_MatvecCommPkgDestroy(jxf_ParCSRMatrixCommPkgT(matrix));
         }
      }
      if (jxf_ParCSRMatrixOwnsRowStarts(matrix))
      {
         jxf_TFree(jxf_ParCSRMatrixRowStarts(matrix));
      }
      if (jxf_ParCSRMatrixOwnsColStarts(matrix))
      {
         jxf_TFree(jxf_ParCSRMatrixColStarts(matrix));
      }

      jxf_TFree(jxf_ParCSRMatrixRowindices(matrix));
      jxf_TFree(jxf_ParCSRMatrixRowvalues(matrix));

      if (jxf_ParCSRMatrixAssumedPartition(matrix))
      {
         jxf_ParCSRMatrixDestroyAssumedPartition(matrix);
      }

      jxf_TFree(matrix);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetLocalRange
 * \brief Returns the row numbers of the rows stored on this processor.
 *       "End" is actually the row number of the last row on this processor.
 * \date 2011/09/08
 */
JXF_Int
jxf_ParCSRMatrixGetLocalRange(jxf_ParCSRMatrix *matrix,
                              JXF_Int *row_start,
                              JXF_Int *row_end,
                              JXF_Int *col_start,
                              JXF_Int *col_end)
{
   JXF_Int my_id;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_MPI_Comm_rank(jxf_ParCSRMatrixComm(matrix), &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION
   *row_start = jxf_ParCSRMatrixFirstRowIndex(matrix);
   *row_end = jxf_ParCSRMatrixLastRowIndex(matrix);
   *col_start = jxf_ParCSRMatrixFirstColDiag(matrix);
   *col_end = jxf_ParCSRMatrixLastColDiag(matrix);
#else
   *row_start = jxf_ParCSRMatrixRowStarts(matrix)[my_id];
   *row_end = jxf_ParCSRMatrixRowStarts(matrix)[my_id + 1] - 1;
   *col_start = jxf_ParCSRMatrixColStarts(matrix)[my_id];
   *col_end = jxf_ParCSRMatrixColStarts(matrix)[my_id + 1] - 1;
#endif

   return jxf_error_flag;
}

JXF_Int
JXF_ParCSRMatrixGetRow(JXF_ParCSRMatrix matrix,
                       JXF_Int row,
                       JXF_Int *size,
                       JXF_Int **col_ind,
                       JXF_Real **values)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParCSRMatrixGetRow((jxf_ParCSRMatrix *)matrix, row, size, col_ind, values);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetRow
 * \brief Returns global column indices and/or values for a given row in the global
 * matrix. Global row number is used, but the row must be stored locally or
 * an error is returned. This implementation copies from the two matrices that
 * store the local data, storing them in the hypre_ParCSRMatrix structure.
 * Only a single row can be accessed via this function at any one time; the
 * corresponding RestoreRow function must be called, to avoid bleeding memory,
 * and to be able to look at another row.
 * Either one of col_ind and values can be left null, and those values will
 * not be returned.
 * All indices are returned in 0-based indexing, no matter what is used under
 * the hood. EXCEPTION: currently this only works if the local CSR matrices
 * use 0-based indexing.
 * \date 2013/01/23
 */
JXF_Int
jxf_ParCSRMatrixGetRow(jxf_ParCSRMatrix *matrix,
                       JXF_Int row,
                       JXF_Int *size,
                       JXF_Int **col_ind,
                       JXF_Real **values)
{
   JXF_Int my_id;
   JXF_Int row_start, row_end;
   jxf_CSRMatrix *Aa;
   jxf_CSRMatrix *Ba;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   Aa = (jxf_CSRMatrix *)jxf_ParCSRMatrixDiag(matrix);
   Ba = (jxf_CSRMatrix *)jxf_ParCSRMatrixOffd(matrix);

   if (jxf_ParCSRMatrixGetrowactive(matrix))
      return -1;

   jxf_MPI_Comm_rank(jxf_ParCSRMatrixComm(matrix), &my_id);

   jxf_ParCSRMatrixGetrowactive(matrix) = 1;
#ifdef JXF_NO_GLOBAL_PARTITION
   row_start = jxf_ParCSRMatrixFirstRowIndex(matrix);
   row_end = jxf_ParCSRMatrixLastRowIndex(matrix) + 1;
#else
   row_end = jxf_ParCSRMatrixRowStarts(matrix)[my_id + 1];
   row_start = jxf_ParCSRMatrixRowStarts(matrix)[my_id];
#endif
   if (row < row_start || row >= row_end)
      return -1;

   if (!jxf_ParCSRMatrixRowvalues(matrix) && (col_ind || values))
   {
      JXF_Int max = 1, tmp;
      JXF_Int i;
      JXF_Int m = row_end - row_start;

      for (i = 0; i < m; i++)
      {
         tmp = jxf_CSRMatrixI(Aa)[i + 1] - jxf_CSRMatrixI(Aa)[i] + jxf_CSRMatrixI(Ba)[i + 1] - jxf_CSRMatrixI(Ba)[i];
         if (max < tmp)
         {
            max = tmp;
         }
      }
      jxf_ParCSRMatrixRowvalues(matrix) = (JXF_Real *)jxf_CTAlloc(JXF_Real, max);
      jxf_ParCSRMatrixRowindices(matrix) = (JXF_Int *)jxf_CTAlloc(JXF_Int, max);
   }

   JXF_Real *vworkA, *vworkB, *v_p;
   JXF_Int i, *cworkA, *cworkB, cstart = jxf_ParCSRMatrixFirstColDiag(matrix);
   JXF_Int nztot, nzA, nzB, lrow = row - row_start;
   JXF_Int *cmap, *idx_p;

   nzA = jxf_CSRMatrixI(Aa)[lrow + 1] - jxf_CSRMatrixI(Aa)[lrow];
   cworkA = &(jxf_CSRMatrixJ(Aa)[jxf_CSRMatrixI(Aa)[lrow]]);
   vworkA = &(jxf_CSRMatrixData(Aa)[jxf_CSRMatrixI(Aa)[lrow]]);

   nzB = jxf_CSRMatrixI(Ba)[lrow + 1] - jxf_CSRMatrixI(Ba)[lrow];
   cworkB = &(jxf_CSRMatrixJ(Ba)[jxf_CSRMatrixI(Ba)[lrow]]);
   vworkB = &(jxf_CSRMatrixData(Ba)[jxf_CSRMatrixI(Ba)[lrow]]);

   nztot = nzA + nzB;

   cmap = jxf_ParCSRMatrixColMapOffd(matrix);

   if (values || col_ind)
   {
      if (nztot)
      {
         JXF_Int imark = -1;

         if (values)
         {
            *values = v_p = jxf_ParCSRMatrixRowvalues(matrix);
            for (i = 0; i < nzB; i++)
            {
               if (cmap[cworkB[i]] < cstart)
                  v_p[i] = vworkB[i];
               else
                  break;
            }
            imark = i;
            for (i = 0; i < nzA; i++)
               v_p[imark + i] = vworkA[i];
            for (i = imark; i < nzB; i++)
               v_p[nzA + i] = vworkB[i];
         }
         if (col_ind)
         {
            *col_ind = idx_p = jxf_ParCSRMatrixRowindices(matrix);
            if (imark > -1)
            {
               for (i = 0; i < imark; i++)
               {
                  idx_p[i] = cmap[cworkB[i]];
               }
            }
            else
            {
               for (i = 0; i < nzB; i++)
               {
                  if (cmap[cworkB[i]] < cstart)
                     idx_p[i] = cmap[cworkB[i]];
                  else
                     break;
               }
               imark = i;
            }
            for (i = 0; i < nzA; i++)
               idx_p[imark + i] = cstart + cworkA[i];
            for (i = imark; i < nzB; i++)
               idx_p[nzA + i] = cmap[cworkB[i]];
         }
      }
      else
      {
         if (col_ind)
            *col_ind = 0;
         if (values)
            *values = 0;
      }
   }
   *size = nztot;

   return jxf_error_flag;
}

JXF_Int
JXF_ParCSRMatrixRestoreRow(JXF_ParCSRMatrix matrix,
                           JXF_Int row,
                           JXF_Int *size,
                           JXF_Int **col_ind,
                           JXF_Real **values)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParCSRMatrixRestoreRow((jxf_ParCSRMatrix *)matrix, row, size, col_ind, values);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixRestoreRow
 * \brief Restore Row
 * \date 2013/01/23
 */
JXF_Int
jxf_ParCSRMatrixRestoreRow(jxf_ParCSRMatrix *matrix,
                           JXF_Int row,
                           JXF_Int *size,
                           JXF_Int **col_ind,
                           JXF_Real **values)
{
   if (!jxf_ParCSRMatrixGetrowactive(matrix))
   {
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   jxf_ParCSRMatrixGetrowactive(matrix) = 0;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixCreateAssumedPartition
 * \brief  Each proc gets it own range. Then each needs to
 *         reconcile its actual range with its assumed
 *         range - the result is essentila a partition of its
 *         assumed range - this is the assumed partition.
 * \date 2011/09/08
 */
JXF_Int
jxf_ParCSRMatrixCreateAssumedPartition(jxf_ParCSRMatrix *matrix)
{
   JXF_Int global_num_cols;
   JXF_Int myid;
   JXF_Int row_start = 0, row_end = 0, col_start = 0, col_end = 0;

   MPI_Comm comm;

   jxf_IJAssumedPart *apart;

   global_num_cols = jxf_ParCSRMatrixGlobalNumCols(matrix);
   comm = jxf_ParCSRMatrixComm(matrix);

   /* find out my actualy range of rows and columns */
   jxf_ParCSRMatrixGetLocalRange(matrix, &row_start, &row_end, &col_start, &col_end);
   jxf_MPI_Comm_rank(comm, &myid);

   /* allocate space */
   apart = jxf_CTAlloc(jxf_IJAssumedPart, 1);

   /* get my assumed partitioning - we want partitioning of the vector that the
      matrix multiplies - so we use the col start and end */
   jxf_GetAssumedPartitionRowRange(comm, myid, 0, global_num_cols, &(apart->row_start), &(apart->row_end));

   /* allocate some space for the partition of the assumed partition */
   apart->length = 0;
   /* room for 10 owners of the assumed partition */
   apart->storage_length = 10; /* need to be >=1 */
   apart->proc_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_start_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_end_list = jxf_TAlloc(JXF_Int, apart->storage_length);

   /* now we want to reconcile our actual partition with the assumed partition */
   jxf_LocateAssummedPartition(comm, col_start, col_end, 0, global_num_cols, apart, myid);

   /* this partition will be saved in the matrix data structure until the matrix is destroyed */
   jxf_ParCSRMatrixAssumedPartition(matrix) = apart;

   return jxf_error_flag;
}

#define CONTACT(a, b) (contact_list[(a) * 3 + (b)])

/*!
 * \fn JXF_Int jxf_LocateAssummedPartition
 * \brief Reconcile assumed partition with actual partition. Essentially
 *        each processor ends of with a partition of its assumed partition.
 * \date 2011/09/08
 */
JXF_Int
jxf_LocateAssummedPartition(MPI_Comm comm,
                            JXF_Int row_start,
                            JXF_Int row_end,
                            JXF_Int global_first_row,
                            JXF_Int global_num_rows,
                            jxf_IJAssumedPart *part,
                            JXF_Int myid)
{
   JXF_Int i;

   JXF_Int *contact_list;
   JXF_Int contact_list_length, contact_list_storage;

   JXF_Int contact_row_start[2], contact_row_end[2], contact_ranges;
   JXF_Int owner_start, owner_end;
   JXF_Int tmp_row_start, tmp_row_end, complete;

   //   JXF_Int        locate_row_start[2], locate_ranges;
   JXF_Int locate_ranges;

   JXF_Int locate_row_count, rows_found;

   JXF_Int tmp_range[2];
   JXF_Int *si, *sortme;

   const JXF_Int flag1 = 17;

   MPI_Request *requests;
   MPI_Status status0, *statuses;

   /*------------------------------------------------------------------
    *  Contact ranges - which rows do I have that others  are assumed
    *  responsible for? (at most two ranges - maybe none)
    *----------------------------------------------------------------*/

   contact_row_start[0] = 0;
   contact_row_end[0] = 0;
   contact_row_start[1] = 0;
   contact_row_end[1] = 0;
   contact_ranges = 0;

   if (row_start <= row_end) /* must own at least one row */
   {
      /* no overlap - so all of my rows and only one range */
      if (part->row_end < row_start || row_end < part->row_start)
      {
         contact_row_start[0] = row_start;
         contact_row_end[0] = row_end;
         contact_ranges++;
      }
      else /* the two regions overlap - so one or two ranges */
      {
         /* check for contact rows on the low end of the local range */
         if (row_start < part->row_start)
         {
            contact_row_start[0] = row_start;
            contact_row_end[0] = part->row_start - 1;
            contact_ranges++;
         }
         if (part->row_end < row_end) /* check the high end */
         {
            if (contact_ranges) /* already found one range */
            {
               contact_row_start[1] = part->row_end + 1;
               contact_row_end[1] = row_end;
            }
            else
            {
               contact_row_start[0] = part->row_end + 1;
               contact_row_end[0] = row_end;
            }
            contact_ranges++;
         }
      }
   }

   /*-----------------------------------------------------------
    *  Contact: find out who is assumed responsible for
    *  these ranges of contact rows and contact them
    *-----------------------------------------------------------*/
   contact_list_length = 0;
   contact_list_storage = 5;
   contact_list = jxf_TAlloc(JXF_Int, contact_list_storage * 3); /* each contact needs 3 ints */

   for (i = 0; i < contact_ranges; i++)
   {
      /* get start and end row owners */
      jxf_GetAssumedPartitionProcFromRow(comm, contact_row_start[i], global_first_row,
                                         global_num_rows, &owner_start);
      jxf_GetAssumedPartitionProcFromRow(comm, contact_row_end[i], global_first_row,
                                         global_num_rows, &owner_end);

      if (owner_start == owner_end) /* same processor owns the whole range */
      {

         if (contact_list_length == contact_list_storage)
         {
            /*allocate more space*/
            contact_list_storage += 5;
            contact_list = jxf_TReAlloc(contact_list, JXF_Int, (contact_list_storage * 3));
         }
         CONTACT(contact_list_length, 0) = owner_start;          /* proc # */
         CONTACT(contact_list_length, 1) = contact_row_start[i]; /* start row */
         CONTACT(contact_list_length, 2) = contact_row_end[i];   /* end row */
         contact_list_length++;
      }
      else
      {
         complete = 0;
         while (!complete)
         {
            jxf_GetAssumedPartitionRowRange(comm, owner_start, global_first_row,
                                            global_num_rows, &tmp_row_start, &tmp_row_end);
            if (tmp_row_end >= contact_row_end[i])
            {
               tmp_row_end = contact_row_end[i];
               complete = 1;
            }
            if (tmp_row_start < contact_row_start[i])
            {
               tmp_row_start = contact_row_start[i];
            }

            if (contact_list_length == contact_list_storage)
            {
               /*allocate more space*/
               contact_list_storage += 5;
               contact_list = jxf_TReAlloc(contact_list, JXF_Int, (contact_list_storage * 3));
            }

            CONTACT(contact_list_length, 0) = owner_start;   /* proc # */
            CONTACT(contact_list_length, 1) = tmp_row_start; /* start row */
            CONTACT(contact_list_length, 2) = tmp_row_end;   /* end row */
            contact_list_length++;
            owner_start++; /* processors are seqential */
         }
      }
   }

   requests = jxf_CTAlloc(MPI_Request, contact_list_length);
   statuses = jxf_CTAlloc(MPI_Status, contact_list_length);

   /* send out messages */
   for (i = 0; i < contact_list_length; i++)
   {
      jxf_MPI_Isend(&CONTACT(i, 1), 2, JXF_MPI_INT, CONTACT(i, 0), flag1, comm, &requests[i]);
   }

   /*-----------------------------------------------------------
    *  Locate ranges -
    *  which rows in my assumed range do I not own
    *  (at most two ranges - maybe none)
    *  locate_row_count = total number of rows I must locate
    *-----------------------------------------------------------*/

   locate_row_count = 0;
   // locate_row_start[0] = 0;
   // locate_row_start[1] = 0;

   locate_ranges = 0;

   /* no overlap - so all of my assumed rows */
   if (part->row_end < row_start || row_end < part->row_start)
   {
      // locate_row_start[0] = part->row_start;
      locate_ranges++;
      locate_row_count += part->row_end - part->row_start + 1;
   }
   else /* the two regions overlap */
   {
      /* check for locate rows on the low end of the local range */
      if (part->row_start < row_start)
      {
         // locate_row_start[0] = part->row_start;
         locate_ranges++;
         locate_row_count += (row_start - 1) - part->row_start + 1;
      }
      if (row_end < part->row_end) /* check the high end */
      {
         if (locate_ranges) /* already have one range */
         {
            // locate_row_start[1] = row_end + 1;
         }
         else
         {
            // locate_row_start[0] = row_end + 1;
         }
         locate_ranges++;
         locate_row_count += part->row_end - (row_end + 1) + 1;
      }
   }

   /*-----------------------------------------------------------
    * Receive messages from other procs telling us where
    * all our  locate rows actually reside
    *-----------------------------------------------------------*/

   /* we will keep a partition of our assumed partition - list
      ourselves first.  We will sort later with an additional index.
      In practice, this should only contain a few processors */

   /* which part do I own? */
   tmp_row_start = jxf_max(part->row_start, row_start);
   tmp_row_end = jxf_min(row_end, part->row_end);

   if (tmp_row_start <= tmp_row_end)
   {
      part->proc_list[0] = myid;
      part->row_start_list[0] = tmp_row_start;
      part->row_end_list[0] = tmp_row_end;
      part->length++;
   }

   /* now look for messages that tell us which processor has our locate rows */
   /* these will be blocking receives as we know how many to expect and they
      should be waiting (and we don't want to continue on without them) */

   rows_found = 0;

   while (rows_found != locate_row_count)
   {
      jxf_MPI_Recv(tmp_range, 2, JXF_MPI_INT, MPI_ANY_SOURCE, flag1, comm, &status0);

      if (part->length == part->storage_length)
      {
         part->storage_length += 10;
         part->proc_list = jxf_TReAlloc(part->proc_list, JXF_Int, part->storage_length);
         part->row_start_list = jxf_TReAlloc(part->row_start_list, JXF_Int, part->storage_length);
         part->row_end_list = jxf_TReAlloc(part->row_end_list, JXF_Int, part->storage_length);
      }
      part->row_start_list[part->length] = tmp_range[0];
      part->row_end_list[part->length] = tmp_range[1];

      part->proc_list[part->length] = status0.MPI_SOURCE;
      rows_found += tmp_range[1] - tmp_range[0] + 1;

      part->length++;
   }

   /* In case the partition of the assumed partition is longish,
      we would like to know the sorted order */
   si = jxf_CTAlloc(JXF_Int, part->length);
   sortme = jxf_CTAlloc(JXF_Int, part->length);

   for (i = 0; i < part->length; i++)
   {
      si[i] = i;
      sortme[i] = part->row_start_list[i];
   }
   jxf_qsort2i(sortme, si, 0, (part->length) - 1);
   part->sort_index = si;

   /* free the requests */
   jxf_MPI_Waitall(contact_list_length, requests, statuses);

   jxf_TFree(statuses);
   jxf_TFree(requests);

   jxf_TFree(sortme);
   jxf_TFree(contact_list);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GetAssumedPartitionProcFromRow
 * \brief Assumed partition for IJ case. Given a particular row j,
 *        return the processor that is assumed to own that row.
 * \date 2011/09/08
 */
JXF_Int
jxf_GetAssumedPartitionProcFromRow(MPI_Comm comm,
                                   JXF_Int row,
                                   JXF_Int global_first_row,
                                   JXF_Int global_num_rows,
                                   JXF_Int *proc_id)
{
   JXF_Int num_procs;
   JXF_Int size, switch_row, extra;

   jxf_MPI_Comm_size(comm, &num_procs);

   /* this looks a bit odd, but we have to be very careful that
      this function and the next are inverses - and rounding
      errors make this difficult!!!!! */

   size = global_num_rows / num_procs;
   extra = global_num_rows - size * num_procs;
   switch_row = global_first_row + (size + 1) * extra;

   if (row >= switch_row)
   {
      *proc_id = extra + (row - switch_row) / size;
   }
   else
   {
      *proc_id = (row - global_first_row) / (size + 1);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_FillResponseParToCSRMatrix
 * \brief Fill response function for determining the send
 *        processors data exchange.
 * \date 2011/09/08
 */
JXF_Int
jxf_FillResponseParToCSRMatrix(void *p_recv_contact_buf,
                               JXF_Int contact_size,
                               JXF_Int contact_proc,
                               void *ro,
                               MPI_Comm comm,
                               void **p_send_response_buf,
                               JXF_Int *response_message_size)
{
   JXF_Int myid, i, index, count, elength;
   JXF_Int *recv_contact_buf = (JXF_Int *)p_recv_contact_buf;

   jxf_DataExchangeResponse *response_obj = ro;
   jxf_ProcListElements *send_proc_obj = response_obj->data2;

   jxf_MPI_Comm_rank(comm, &myid);

   /* check to see if we need to allocate more space in send_proc_obj for ids */
   if (send_proc_obj->length == send_proc_obj->storage_length)
   {
      send_proc_obj->storage_length += 10; /* add space for 10 more processors */
      send_proc_obj->id = jxf_TReAlloc(send_proc_obj->id, JXF_Int, send_proc_obj->storage_length);
      send_proc_obj->vec_starts = jxf_TReAlloc(send_proc_obj->vec_starts, JXF_Int, send_proc_obj->storage_length + 1);
   }

   /* initialize */
   count = send_proc_obj->length;
   index = send_proc_obj->vec_starts[count]; /* this is the number of elements */

   /* send proc */
   send_proc_obj->id[count] = contact_proc;

   /* do we need more storage for the elements? */
   if (send_proc_obj->element_storage_length < index + contact_size)
   {
      elength = jxf_max(contact_size, 10);
      elength += index;
      send_proc_obj->elements = jxf_TReAlloc(send_proc_obj->elements, JXF_Int, elength);
      send_proc_obj->element_storage_length = elength;
   }

   /* populate send_proc_obj */
   for (i = 0; i < contact_size; i++)
   {
      send_proc_obj->elements[index++] = recv_contact_buf[i];
   }
   send_proc_obj->vec_starts[count + 1] = index;
   send_proc_obj->length++;

   /*output - no message to return (confirmation) */
   *response_message_size = 0;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixDestroyAssumedPartition
 * \date 2011/09/08
 */
JXF_Int
jxf_ParCSRMatrixDestroyAssumedPartition(jxf_ParCSRMatrix *matrix)
{
   jxf_IJAssumedPart *apart;
   apart = jxf_ParCSRMatrixAssumedPartition(matrix);

   if (apart->storage_length > 0)
   {
      jxf_TFree(apart->proc_list);
      jxf_TFree(apart->row_start_list);
      jxf_TFree(apart->row_end_list);
      jxf_TFree(apart->sort_index);
   }

   jxf_TFree(apart);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetRowPartitioning
 * \author peghoty
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixGetRowPartitioning(jxf_ParCSRMatrix *matrix, JXF_Int **row_partitioning_ptr)
{
   JXF_Int *row_partitioning, *row_starts;
   JXF_Int num_procs, i;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_MPI_Comm_size(jxf_ParCSRMatrixComm(matrix), &num_procs);
   row_starts = jxf_ParCSRMatrixRowStarts(matrix);
   if (!row_starts)
      return -1;
   row_partitioning = jxf_CTAlloc(JXF_Int, num_procs + 1);
   for (i = 0; i < num_procs + 1; i++)
   {
      row_partitioning[i] = row_starts[i];
   }

   *row_partitioning_ptr = row_partitioning;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetColPartitioning
 * \author peghoty
 * \date 2011/08/27
 */
JXF_Int
jxf_ParCSRMatrixGetColPartitioning(jxf_ParCSRMatrix *matrix, JXF_Int **col_partitioning_ptr)
{
   JXF_Int *col_partitioning, *col_starts;
   JXF_Int num_procs, i;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_MPI_Comm_size(jxf_ParCSRMatrixComm(matrix), &num_procs);
   col_starts = jxf_ParCSRMatrixColStarts(matrix);
   if (!col_starts)
      return -1;
   col_partitioning = jxf_CTAlloc(JXF_Int, num_procs + 1);
   for (i = 0; i < num_procs + 1; i++)
   {
      col_partitioning[i] = col_starts[i];
   }

   *col_partitioning_ptr = col_partitioning;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetNumNonzeros
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixSetNumNonzeros(jxf_ParCSRMatrix *matrix)
{
   MPI_Comm comm;
   jxf_CSRMatrix *diag;
   JXF_Int *diag_i;
   jxf_CSRMatrix *offd;
   JXF_Int *offd_i;
   JXF_Int local_num_rows;
   JXF_Int total_num_nonzeros;
   JXF_Int local_num_nonzeros;
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   comm = jxf_ParCSRMatrixComm(matrix);
   diag = jxf_ParCSRMatrixDiag(matrix);
   diag_i = jxf_CSRMatrixI(diag);
   offd = jxf_ParCSRMatrixOffd(matrix);
   offd_i = jxf_CSRMatrixI(offd);
   local_num_rows = jxf_CSRMatrixNumRows(diag);

   local_num_nonzeros = diag_i[local_num_rows] + offd_i[local_num_rows];
   jxf_MPI_Allreduce(&local_num_nonzeros, &total_num_nonzeros, 1, JXF_MPI_INT, MPI_SUM, comm);
   jxf_ParCSRMatrixNumNonzeros(matrix) = total_num_nonzeros;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxh_ParCSRMatrixSetNumNonzeros
 * \brief Set the global number of nonzeros for a half-precision ParCSR matrix.
 */
JXF_Int
jxh_ParCSRMatrixSetNumNonzeros(jxh_ParCSRMatrix *matrix)
{
   MPI_Comm comm;
   jxh_CSRMatrix *diag;
   JXF_Int *diag_i;
   jxh_CSRMatrix *offd;
   JXF_Int *offd_i;
   JXF_Int local_num_rows;
   JXF_Int total_num_nonzeros;
   JXF_Int local_num_nonzeros;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   comm = jxh_ParCSRMatrixComm(matrix);
   diag = jxh_ParCSRMatrixDiag(matrix);
   diag_i = jxh_CSRMatrixI(diag);
   offd = jxh_ParCSRMatrixOffd(matrix);
   offd_i = jxh_CSRMatrixI(offd);

   local_num_rows = jxh_CSRMatrixNumRows(diag);
   local_num_nonzeros = diag_i[local_num_rows] + offd_i[local_num_rows];

   jxf_MPI_Allreduce(&local_num_nonzeros,
                     &total_num_nonzeros,
                     1,
                     JXF_MPI_INT,
                     MPI_SUM,
                     comm);

   jxh_ParCSRMatrixNumNonzeros(matrix) = total_num_nonzeros;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetDNumNonzeros
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixSetDNumNonzeros(jxf_ParCSRMatrix *matrix)
{
   MPI_Comm comm;
   jxf_CSRMatrix *diag;
   JXF_Int *diag_i;
   jxf_CSRMatrix *offd;
   JXF_Int *offd_i;
   JXF_Int local_num_rows;
   JXF_Real total_num_nonzeros;
   JXF_Real local_num_nonzeros;
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   comm = jxf_ParCSRMatrixComm(matrix);
   diag = jxf_ParCSRMatrixDiag(matrix);
   diag_i = jxf_CSRMatrixI(diag);
   offd = jxf_ParCSRMatrixOffd(matrix);
   offd_i = jxf_CSRMatrixI(offd);
   local_num_rows = jxf_CSRMatrixNumRows(diag);
   local_num_nonzeros = (JXF_Real)diag_i[local_num_rows] + (JXF_Real)offd_i[local_num_rows];
   jxf_MPI_Allreduce(&local_num_nonzeros, &total_num_nonzeros, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matrix) = total_num_nonzeros;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxh_ParCSRMatrixSetDNumNonzeros
 * \brief Set global number of nonzeros as JXF_Real for a half ParCSR matrix.
 */
JXF_Int
jxh_ParCSRMatrixSetDNumNonzeros(jxh_ParCSRMatrix *matrix)
{
   MPI_Comm comm;
   jxh_CSRMatrix *diag;
   JXF_Int *diag_i;
   jxh_CSRMatrix *offd;
   JXF_Int *offd_i;
   JXF_Int local_num_rows;
   JXF_Real total_num_nonzeros;
   JXF_Real local_num_nonzeros;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   comm = jxh_ParCSRMatrixComm(matrix);

   diag = jxh_ParCSRMatrixDiag(matrix);
   diag_i = jxh_CSRMatrixI(diag);

   offd = jxh_ParCSRMatrixOffd(matrix);
   offd_i = jxh_CSRMatrixI(offd);

   local_num_rows = jxh_CSRMatrixNumRows(diag);

   local_num_nonzeros =
       (JXF_Real)diag_i[local_num_rows] +
       (JXF_Real)offd_i[local_num_rows];

   jxf_MPI_Allreduce(&local_num_nonzeros,
                     &total_num_nonzeros,
                     1,
                     JXF_MPI_REAL,
                     MPI_SUM,
                     comm);

   jxh_ParCSRMatrixDNumNonzeros(matrix) = total_num_nonzeros;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetDataOwner
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixSetDataOwner(jxf_ParCSRMatrix *matrix,
                             JXF_Int owns_data)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_ParCSRMatrixOwnsData(matrix) = owns_data;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetRowStartsOwner
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixSetRowStartsOwner(jxf_ParCSRMatrix *matrix, JXF_Int owns_row_starts)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_ParCSRMatrixOwnsRowStarts(matrix) = owns_row_starts;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetColStartsOwner
 * \date 2009/07/10
 */
JXF_Int
jxf_ParCSRMatrixSetColStartsOwner(jxf_ParCSRMatrix *matrix, JXF_Int owns_col_starts)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_ParCSRMatrixOwnsColStarts(matrix) = owns_col_starts;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BuildCSRMatrixMPIDataType
 * \date 2009/07/10
 */
JXF_Int
jxf_BuildCSRMatrixMPIDataType(JXF_Int num_nonzeros,
                              JXF_Int num_rows,
                              JXF_Real *a_data,
                              JXF_Int *a_i,
                              JXF_Int *a_j,
                              MPI_Datatype *csr_matrix_datatype)
{
   JXF_Int block_lens[3];
   MPI_Aint displ[3];
   MPI_Datatype types[3];

   block_lens[0] = num_nonzeros;
   block_lens[1] = num_rows + 1;
   block_lens[2] = num_nonzeros;

   types[0] = JXF_MPI_REAL;
   types[1] = JXF_MPI_INT;
   types[2] = JXF_MPI_INT;

   jxf_MPI_Address(a_data, &displ[0]);
   jxf_MPI_Address(a_i, &displ[1]);
   jxf_MPI_Address(a_j, &displ[2]);
   jxf_MPI_Type_struct(3, block_lens, displ, types, csr_matrix_datatype);
   jxf_MPI_Type_commit(csr_matrix_datatype);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_Build_Par_LinearSystem
 * \brief Build a parallel linear system by reading data
 *        from given files.
 * \author peghoty
 * \date 2009/08/28
 */
JXF_Int
jxf_Build_Par_LinearSystem(MPI_Comm comm,
                           char *filenamemat,
                           char *filenamerhs,
                           jxf_ParCSRMatrix **par_matrix_ptr,
                           jxf_ParVector **par_rhs_ptr)
{
   FILE *fp = NULL;
   JXF_Int nrows, nonzeros;
   JXF_Int i;
   JXF_Int *ia, *ja;
   JXF_Real *a, *f;

   JXF_Int *row_part, *col_part; /* 并行矩阵的行列分划 */
   JXF_Int *partitioning;        /* 并行向量的分划，必须与并行矩阵的行分划一致 */

   jxf_CSRMatrix *matrix;
   jxf_Vector *rhs;

   jxf_ParCSRMatrix *par_matrix;
   jxf_ParVector *par_rhs;

   /*---------------------------------------------------
    *     read the matrix from file
    *--------------------------------------------------*/

   if ((fp = fopen(filenamemat, "r")) == NULL)
   {
      jxf_printf("Can't open the file %s !\n", filenamemat);
      exit(0);
   }

   jxf_fscanf(fp, "%d\n", &nrows);
   ia = (JXF_Int *)malloc((nrows + 1) * sizeof(JXF_Int));

   for (i = 0; i < nrows + 1; i++)
   {
      jxf_fscanf(fp, "%d\n", &ia[i]);
      ia[i]--;
   }

   nonzeros = ia[nrows];
   ja = (JXF_Int *)malloc(nonzeros * sizeof(JXF_Int));
   a = (JXF_Real *)malloc(nonzeros * sizeof(JXF_Real));

   for (i = 0; i < nonzeros; i++)
   {
      jxf_fscanf(fp, "%d\n", &ja[i]);
      ja[i]--;
   }

   for (i = 0; i < nonzeros; i++)
   {
      jxf_fscanf(fp, "%le\n", &a[i]);
   }

   fclose(fp);

   /*---------------------------------------------------
    *     Read the rhs from file
    *--------------------------------------------------*/

   f = (JXF_Real *)malloc(nrows * sizeof(JXF_Real));
   if ((fp = fopen(filenamerhs, "r")) == NULL)
   {
      jxf_printf(" Can't open the file %s!\n", filenamerhs);
      exit(0);
   }
   jxf_fscanf(fp, "%d\n", &nrows);
   for (i = 0; i < nrows; i++)
   {
      jxf_fscanf(fp, "%lf\n", &f[i]);
   }
   fclose(fp);

   /*--------------------------------------------------------------------------
    * build a parallel matrix and rhs
    *--------------------------------------------------------------------------*/

   /* 1. Make the diag element to be the first of each row */
   jxf_PutDiagFirst(ia, ja, a, nrows);

   /* 2. build a parallel matrix */
   matrix = jxf_CSRMatrixCreate(nrows, nrows, nonzeros);
   jxf_CSRMatrixI(matrix) = ia;
   jxf_CSRMatrixJ(matrix) = ja;
   jxf_CSRMatrixData(matrix) = a;
   row_part = NULL;
   col_part = NULL;
   par_matrix = jxf_CSRMatrixToParCSRMatrix(comm, matrix, row_part, col_part);

   /* 3. build a parallel rhs */
   rhs = jxf_SeqVectorCreate(nrows);
   jxf_VectorData(rhs) = f;
   jxf_ParCSRMatrixGetRowPartitioning(par_matrix, &partitioning);
   par_rhs = jxf_VectorToParVector(comm, rhs, partitioning);

   *par_matrix_ptr = par_matrix;
   *par_rhs_ptr = par_rhs;

   jxf_TFree(matrix);
   jxf_TFree(rhs);

   return 0;
}

/*!
 * \fn JXF_Int jxf_PutDiagFirst
 * \brief Reorder the matrix in order to make the diagonal entries
 *        are stored firstly in each row.
 * \author peghoty
 * \date 2009/06/17
 */
JXF_Int
jxf_PutDiagFirst(JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n)
{
   JXF_Int i, j;
   JXF_Int col, tmpcol;
   JXF_Real tmpa;

   for (i = 0; i < n; i++)
   {
      for (j = ia[i]; j < ia[i + 1]; j++)
      {
         col = ja[j];
         if (col == i && j != ia[i])
         {
            tmpcol = ja[ia[i]];
            ja[ia[i]] = col;
            ja[j] = tmpcol;

            tmpa = a[ia[i]];
            a[ia[i]] = a[j];
            a[j] = tmpa;
         }
      }
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetDims
 * \brief Get num_rows and num_cols.
 * \author Yue Xiaoqiang
 * \date 2013/01/23
 */
JXF_Int
jxf_ParCSRMatrixGetDims(jxf_ParCSRMatrix *A, JXF_Int *M, JXF_Int *N)
{
   if (!A)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   *M = jxf_ParCSRMatrixGlobalNumRows(A);
   *N = jxf_ParCSRMatrixGlobalNumCols(A);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixGetComm
 * \brief Get comm.
 * \author Yue Xiaoqiang
 * \date 2013/02/27
 */
JXF_Int
jxf_ParCSRMatrixGetComm(jxf_ParCSRMatrix *matrix, MPI_Comm *comm)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   *comm = jxf_ParCSRMatrixComm(matrix);

   return jxf_error_flag;
}
