//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_matrix.c -- operations for parallel matrices.
 *  Date: 2014/07/05
 */

#include "jxf_mv.h"
#include "jxf_ilu.h"

/*!
 * \fn jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile2
 * \brief Build a parallel matrix by reading data from a given file,
 *        including reorder by variables and along y-axis
 * \param filename pointer to the matrix file.
 * \param row_part pointer to the row-start partitioning.
 * \param row_part pointer to the col-start partitioning.
 * \note The storage format of the matrix is expected to be CSR, i.e.
 *         the 1st row:                           n                          // JXF_Int
 *         the 2 to (n+2)-th row:                ia[k]     k=0,1,...,n       // JXF_Int
 *         the (n+3) to (n+nz+2)-th row:         ja[k]     k=0,1,...,nz-1    // JXF_Int
 *         the (n+nz+3) to ((n+2nz+2))-th row:    a[k]     k=0,1,...,nz-1    // JXF_Real 
 *       where n is the order of the matrix, nz is the number of non-zeros of the matrix.
 * \author Yue Xiaoqiang
 * \date 2014/04/10
 */
jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile2( char *filename, 
                            JXF_Int   file_base,
                            JXF_Int  *row_part, 
                            JXF_Int  *col_part,
                            JXF_Int   num_equns,
                            JXF_Int   nx,
                            JXF_Int   ny )
{
   jxf_ParCSRMatrix  *A     = NULL;
   jxf_CSRMatrix     *A_CSR = NULL;
   jxf_CSRMatrix     *B_CSR = NULL;

   JXF_Int my_id;
   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0)
   {
      A_CSR = jxf_CSRMatrixRead(filename, file_base);
//      jxf_CSRMatrixReorder(A_CSR);
      B_CSR = jxf_CSRMatrixTDMGReorderByNodes(A_CSR, num_equns);
//      jxf_CSRMatrixReorder(B_CSR);
      jxf_CSRMatrixDestroy(A_CSR);
      A_CSR = jxf_CSRMatrixTDMGReorderAlongY(B_CSR, num_equns, nx, ny);
      jxf_CSRMatrixDestroy(B_CSR);
      jxf_CSRMatrixReorder(A_CSR);
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
 * \fn jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile3
 * \brief Build a parallel matrix by reading data from a given file,
 *        including reorder by variables and along y-axis, block-partitioned
 * \param filename pointer to the matrix file.
 * \param num_equns number of equations.
 * \param nx number of segments in x-direction.
 * \param ny number of segments in y-direction.
 * \note The storage format of the matrix is expected to be CSR, i.e.
 *         the 1st row:                           n                          // JXF_Int
 *         the 2 to (n+2)-th row:                ia[k]     k=0,1,...,n       // JXF_Int
 *         the (n+3) to (n+nz+2)-th row:         ja[k]     k=0,1,...,nz-1    // JXF_Int
 *         the (n+nz+3) to ((n+2nz+2))-th row:    a[k]     k=0,1,...,nz-1    // JXF_Real 
 *       where n is the order of the matrix, nz is the number of non-zeros of the matrix.
 * \author Yue Xiaoqiang
 * \date 2014/04/13
 */
jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile3( char *filename, JXF_Int file_base, JXF_Int num_equns, JXF_Int nx, JXF_Int ny )
{
   jxf_ParCSRMatrix  *A     = NULL;
   jxf_CSRMatrix     *A_CSR = NULL;
   jxf_CSRMatrix     *B_CSR = NULL;

   JXF_Int my_id;
   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0)
   {
      A_CSR = jxf_CSRMatrixRead(filename, file_base);
//      jxf_CSRMatrixReorder(A_CSR);
      B_CSR = jxf_CSRMatrixTDMGReorderByNodes(A_CSR, num_equns);
//      jxf_CSRMatrixReorder(B_CSR);
      jxf_CSRMatrixDestroy(A_CSR);
      A_CSR = jxf_CSRMatrixTDMGReorderAlongY(B_CSR, num_equns, nx, ny);
      jxf_CSRMatrixDestroy(B_CSR);
      jxf_CSRMatrixReorder(A_CSR);
   }
   
   /* ser -> par transfering */
   A = jxf_CSRMatrixToParCSRMatrix2(MPI_COMM_WORLD, A_CSR, num_equns, nx, ny);

   if (my_id == 0) 
   {
      jxf_CSRMatrixDestroy(A_CSR);
   }

   return (A);
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_ParCSRMatrixCreate2
 * \brief Create a parallel CSR Matrix.
 * \author Yue Xiaoqiang
 * \date 2014/04/13
 */
jxf_ParCSRMatrix *
jxf_ParCSRMatrixCreate2( MPI_Comm   comm,
                        JXF_Int        global_num_rows,
                        JXF_Int        global_num_cols,
                        JXF_Int       *row_starts,
                        JXF_Int       *col_starts,
                        JXF_Int        num_cols_offd,
                        JXF_Int        num_nonzeros_diag,
                        JXF_Int        num_nonzeros_offd,
                        JXF_Int        num_equns,
                        JXF_Int        nx,
                        JXF_Int        ny )
{
   jxf_ParCSRMatrix  *matrix;
   
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
      jxf_GeneratePartitioning2(num_procs, &row_starts, num_equns, nx, ny);
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
         jxf_GeneratePartitioning2(num_procs, &col_starts, num_equns, nx, ny);
#endif
      }
   }


#ifdef JXF_NO_GLOBAL_PARTITION
   /* row_starts[0] is start of local rows.  row_starts[1] is start of next processor's rows */
   first_row_index = row_starts[0];
   local_num_rows  = row_starts[1] - first_row_index ;
   first_col_diag  = col_starts[0];
   local_num_cols  = col_starts[1] - first_col_diag;
#else
   first_row_index = row_starts[my_id];
   local_num_rows  = row_starts[my_id+1] - first_row_index;
   first_col_diag  = col_starts[my_id];
   local_num_cols  = col_starts[my_id+1] - first_col_diag;
#endif


   jxf_ParCSRMatrixComm(matrix) = comm;
   jxf_ParCSRMatrixDiag(matrix) = jxf_CSRMatrixCreate(local_num_rows, local_num_cols, num_nonzeros_diag);
   jxf_ParCSRMatrixOffd(matrix) = jxf_CSRMatrixCreate(local_num_rows, num_cols_offd, num_nonzeros_offd);
   jxf_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   jxf_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;
   jxf_ParCSRMatrixFirstRowIndex(matrix) = first_row_index;
   jxf_ParCSRMatrixFirstColDiag(matrix)  = first_col_diag;
 
   jxf_ParCSRMatrixLastRowIndex(matrix) = first_row_index + local_num_rows - 1;
   jxf_ParCSRMatrixLastColDiag(matrix)  = first_col_diag  + local_num_cols - 1;

   jxf_ParCSRMatrixColMapOffd(matrix)       = NULL;
   jxf_ParCSRMatrixAssumedPartition(matrix) = NULL;


  /*------------------------------------------------------------------------- 
   *   When NO_GLOBAL_PARTITION is set we could make these null, instead
   * of leaving the range.  If that change is made, then when this create
   * is called from functions like the matrix-matrix multiply, be careful
   * not to generate a new partition.
   *-------------------------------------------------------------------------*/

   jxf_ParCSRMatrixRowStarts(matrix) = row_starts;
   jxf_ParCSRMatrixColStarts(matrix) = col_starts;

   jxf_ParCSRMatrixCommPkg(matrix)  = NULL;
   jxf_ParCSRMatrixCommPkgT(matrix) = NULL;

   /* set defaults */
   jxf_ParCSRMatrixOwnsData(matrix) = 1;
   jxf_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   jxf_ParCSRMatrixOwnsColStarts(matrix) = 1;
   if (row_starts == col_starts)
   {
      jxf_ParCSRMatrixOwnsColStarts(matrix) = 0;
   }
   jxf_ParCSRMatrixRowindices(matrix)   = NULL;
   jxf_ParCSRMatrixRowvalues(matrix)    = NULL;
   jxf_ParCSRMatrixGetrowactive(matrix) = 0;

   return matrix;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix2
 * \brief Generates a ParCSRMatrix distributed across the 
 *        processors in comm from a CSRMatrix on proc 0, block-partitioned
 * \note This shouldn't be used with the JXF_NO_GLOBAL_PARTITON option
 * \author Yue Xiaoqiang
 * \date 2014/04/13
 */
jxf_ParCSRMatrix *
jxf_CSRMatrixToParCSRMatrix2( MPI_Comm comm, jxf_CSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny )
{
   JXF_Int         *global_data;
   JXF_Int          global_size;
   JXF_Int          global_num_rows;
   JXF_Int          global_num_cols;
   JXF_Int         *local_num_rows;
   JXF_Int         *row_starts = NULL;
   JXF_Int         *col_starts = NULL;

   JXF_Int          num_procs, my_id;
   JXF_Int         *local_num_nonzeros = NULL;
   JXF_Int          num_nonzeros;
   
   JXF_Real       *a_data = NULL;
   JXF_Int          *a_i    = NULL;
   JXF_Int          *a_j    = NULL;
  
   jxf_CSRMatrix *local_A;

   MPI_Request  *requests;
   MPI_Status   *status, status0;
   MPI_Datatype *csr_matrix_datatypes;

   jxf_ParCSRMatrix *par_matrix;

   JXF_Int first_col_diag;
   JXF_Int last_col_diag;
 
   JXF_Int i, j, ind;

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   global_data = jxf_CTAlloc(JXF_Int, 2*num_procs + 6);
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
               global_size = 2*num_procs+6;
               for (i = 0; i < num_procs + 1; i ++)
               {
                  global_data[i+4] = row_starts[i];
               }
               for (i = 0; i < num_procs + 1; i ++)
               {
                  global_data[i+num_procs+5] = col_starts[i];
               }
	    }
	    else
	    {
               global_data[3] = 0;
               global_size = num_procs+5;
               for (i = 0; i < num_procs + 1; i ++)
               {
                  global_data[i+4] = row_starts[i];
               }
	    }
	 }
	 else
	 {
            global_data[3] = 1;
            global_size = num_procs+5;
	    for (i = 0; i < num_procs + 1; i ++)
	    {
	       global_data[i+4] = row_starts[i];
	    }
	 }
      }
      else 
      {
	 if (col_starts)
	 {
            global_data[3] = 3;
	    global_size = num_procs + 5;
	    for (i = 0; i < num_procs + 1; i ++)
	    {
	       global_data[i+4] = col_starts[i];
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
	    for (i = 0; i < num_procs + 1; i ++)
	    {
	       row_starts[i] = global_data[i+4];
	    }
	    if (global_data[3] == 0)
	    {
	       col_starts = row_starts;
	    }
	    if (global_data[3] == 2)
	    {
	       col_starts = jxf_CTAlloc(JXF_Int, num_procs+1);
	       for (i = 0; i < num_procs + 1; i ++)
	       {
	          col_starts[i] = global_data[i+num_procs+5];
	       }
	    }
	 }
	 else
	 {
	    col_starts = jxf_CTAlloc(JXF_Int, num_procs + 1);
	    for (i = 0; i < num_procs + 1; i ++)
	    {
	       col_starts[i] = global_data[i+4];
	    }
	 }
      }
   }
   jxf_TFree(global_data);

   local_num_rows = jxf_CTAlloc(JXF_Int, num_procs);
   csr_matrix_datatypes = jxf_CTAlloc(MPI_Datatype, num_procs);

   par_matrix = jxf_ParCSRMatrixCreate2(comm, global_num_rows, global_num_cols,
                                       row_starts, col_starts, 0, 0, 0, num_equns, nx, ny);

   row_starts = jxf_ParCSRMatrixRowStarts(par_matrix);
   col_starts = jxf_ParCSRMatrixColStarts(par_matrix);

   for (i = 0; i < num_procs; i ++)
   {
      local_num_rows[i] = row_starts[i+1] - row_starts[i];
   }

   if (my_id == 0)
   {
      local_num_nonzeros = jxf_CTAlloc(JXF_Int, num_procs);
      for (i = 0; i < num_procs - 1; i ++)
      {
         local_num_nonzeros[i] = a_i[row_starts[i+1]] - a_i[row_starts[i]];
      }
      local_num_nonzeros[num_procs-1] = a_i[global_num_rows] - a_i[row_starts[num_procs-1]];
   }
   jxf_MPI_Scatter(local_num_nonzeros, 1, JXF_MPI_INT, &num_nonzeros, 1, JXF_MPI_INT, 0, comm);

   if (my_id == 0) 
   {
      num_nonzeros = local_num_nonzeros[0];
   }

   local_A = jxf_CSRMatrixCreate(local_num_rows[my_id], global_num_cols, num_nonzeros);

   if (my_id == 0)
   {
        requests = jxf_CTAlloc (MPI_Request, num_procs - 1);
        status = jxf_CTAlloc(MPI_Status, num_procs - 1);
        j = 0;
        for ( i = 1; i < num_procs; i ++)
        {
                ind = a_i[row_starts[i]];
                jxf_BuildCSRMatrixMPIDataType( local_num_nonzeros[i], 
                                              local_num_rows[i],
                                              &a_data[ind],
                                              &a_i[row_starts[i]],
                                              &a_j[ind],
                                              &csr_matrix_datatypes[i] );
                jxf_MPI_Isend( MPI_BOTTOM, 1, csr_matrix_datatypes[i], i, 0, comm, &requests[j++] );
                jxf_MPI_Type_free(&csr_matrix_datatypes[i]);
        }
        jxf_CSRMatrixData(local_A) = a_data;
        jxf_CSRMatrixI(local_A) = a_i;
        jxf_CSRMatrixJ(local_A) = a_j;
        jxf_CSRMatrixOwnsData(local_A) = 0;
        jxf_MPI_Waitall(num_procs-1, requests, status);
        jxf_TFree(requests);
        jxf_TFree(status);
        jxf_TFree(local_num_nonzeros);
   }
   else
   {
        jxf_CSRMatrixInitialize(local_A);
        jxf_BuildCSRMatrixMPIDataType( num_nonzeros, 
                                      local_num_rows[my_id],
                                      jxf_CSRMatrixData(local_A),
                                      jxf_CSRMatrixI(local_A),
                                      jxf_CSRMatrixJ(local_A),
                                      csr_matrix_datatypes );
        jxf_MPI_Recv( MPI_BOTTOM, 1, csr_matrix_datatypes[0], 0, 0, comm, &status0 );
        jxf_MPI_Type_free(csr_matrix_datatypes);
   }

   first_col_diag = col_starts[my_id];
   last_col_diag  = col_starts[my_id+1] - 1;

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
 * \fn jxf_CSRMatrix *jxf_MergeDiagAndOffdDropSmall
 * \brief Merge the Diag and Offd for a parallel CSR matrix, except small items
 * \author Yue Xiaoqiang
 * \date 2014/04/24
 */
jxf_CSRMatrix *
jxf_MergeDiagAndOffdDropSmall( jxf_ParCSRMatrix *par_matrix, JXF_Real drop_tol )
{
   jxf_CSRMatrix  *diag = jxf_ParCSRMatrixDiag(par_matrix);
   jxf_CSRMatrix  *offd = jxf_ParCSRMatrixOffd(par_matrix);
   jxf_CSRMatrix  *matrix;

   JXF_Int           num_cols = jxf_ParCSRMatrixGlobalNumCols(par_matrix);
   JXF_Int           first_col_diag = jxf_ParCSRMatrixFirstColDiag(par_matrix);
   JXF_Int          *col_map_offd = jxf_ParCSRMatrixColMapOffd(par_matrix);
   JXF_Int           num_rows = jxf_CSRMatrixNumRows(diag);

   JXF_Int          *diag_i = jxf_CSRMatrixI(diag);
   JXF_Int          *diag_j = jxf_CSRMatrixJ(diag);
   JXF_Real       *diag_data = jxf_CSRMatrixData(diag);
   JXF_Int          *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int          *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Real       *offd_data = jxf_CSRMatrixData(offd);

   JXF_Int          *matrix_i;
   JXF_Int          *matrix_j;
   JXF_Real       *matrix_data;

   JXF_Int          num_nonzeros, i, j;
   JXF_Int          count;
   
   JXF_Real       tmpval;

   num_nonzeros = diag_i[num_rows] + offd_i[num_rows];

   matrix = jxf_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   jxf_CSRMatrixInitialize(matrix);

   matrix_i = jxf_CSRMatrixI(matrix);
   matrix_j = jxf_CSRMatrixJ(matrix);
   matrix_data = jxf_CSRMatrixData(matrix);

   count = 0;
   matrix_i[0] = 0;
   for (i = 0; i < num_rows; i ++)
   {
      tmpval = drop_tol * fabs(diag_data[diag_i[i]]); // Treat the diagonal as the maximal
      for (j = diag_i[i]; j < diag_i[i+1]; j ++)
      {
         if (fabs(diag_data[j]) > tmpval)
         {
            matrix_data[count] = diag_data[j];
            matrix_j[count++] = diag_j[j] + first_col_diag;
         }
      }
      for (j = offd_i[i]; j < offd_i[i+1]; j ++)
      {
         if (fabs(offd_data[j]) > tmpval)
         {
            matrix_data[count] = offd_data[j];
            matrix_j[count++] = col_map_offd[offd_j[j]];
         }
      }
      matrix_i[i+1] = count;
   }
   jxf_CSRMatrixNumNonzeros(matrix) = count; /* Reset num_nonzeros, while mat_j and mat_data needn't */

   return matrix;
}

/*!
 * \fn jxf_CSRMatrix *jxf_CSRMatrixMergeReorderIntUpRtDnLtQuasiBdy
 * \brief Merge and reorder in JXF_Int-up-rgt-dwn-lft
 * \author Yue Xiaoqiang
 * \date 2014/08/25
 */
jxf_CSRMatrix *
jxf_CSRMatrixMergeReorderIntUpRtDnLtQuasiBdy( MPI_Comm              comm,
                                             JXF_Int                  *row_starts,
                                             JXF_Int                  *permute,
                                             jxf_CSRMatrix         *ser_B,
                                             JXF_Int                   nz_srt,
                                             JXF_Int                   ng_pt,
                                             jxf_GridPartitionData *grid_data,
                                             JXF_Int                  *ex_len,
                                             JXF_Int                  *ey_len,
                                             JXF_Int                  *dx_len,
                                             JXF_Int                  *kx_len,
                                             JXF_Int                  *dy_len,
                                             JXF_Int                  *lx_len,
                                             JXF_Int                  *ly_len,
                                             JXF_Int                  *postn_a,
                                             JXF_Int                  *postn_b,
                                             JXF_Int                  *postn_c )
{
    JXF_Int num_rows = jxf_CSRMatrixNumRows(ser_B);
    JXF_Int num_cols = jxf_CSRMatrixNumCols(ser_B);
    JXF_Int *IB = jxf_CSRMatrixI(ser_B);
    JXF_Int *JB = jxf_CSRMatrixJ(ser_B);
    JXF_Real *BB = jxf_CSRMatrixData(ser_B);
    
    JXF_Int x_part_len = grid_data->x_part_len;
    JXF_Int y_part_len = grid_data->y_part_len;
    JXF_Int num_sideprocs = grid_data->num_nocrossside; // YUE: Needn't to consider the cross ones in 5-point stencil
    JXF_Int num_smallside = grid_data->num_smallside;
    JXF_Int num_largeside = grid_data->num_largeside;
    JXF_Int *xlo_array = grid_data->xlo_array;
    JXF_Int *xup_array = grid_data->xup_array;
    JXF_Int *ylo_array = grid_data->ylo_array;
    JXF_Int *yup_array = grid_data->yup_array;
    JXF_Int *smallprocs = grid_data->sideprocs;
    JXF_Int *largeprocs = smallprocs + num_smallside;
    JXF_Int *smallprcpos = grid_data->sideprcpos;
    JXF_Int *largeprcpos = smallprcpos + num_smallside;
    
    jxf_CSRMatrix *ser_A = NULL;
    JXF_Int *IA = NULL;
    JXF_Int *JA = NULL;
    JXF_Real *AA = NULL;
    
    JXF_Int ng_pt_mx = ng_pt * x_part_len;
    JXF_Int ng_pt_my = ng_pt * y_part_len;
    JXF_Int ng_pt_m2x = 2 * ng_pt_mx;
    JXF_Int ng_pt_m_px = ng_pt_mx + ng_pt;
    JXF_Int ng_pt_m_xo = ng_pt_mx - ng_pt;
    
    JXF_Int num_procs, my_id, recv_rows, rnt, cnt, ctt, i, j, k, l, pid, rele, col;
    JXF_Int plle, pile, pkle, pkke, plme, nloc, nvar, xloc, yloc, row_end, pyln, pxln;
    JXF_Int erow, jsrt, jend, pke, ple, psrt, pend, povr, psse, puue, ptte, pjje;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    if (grid_data->y_lower_idx == 0) // without DOWN
    {
        if (grid_data->x_lower_idx == 0) // without LEFT
        {
            recv_rows = 0;
           *ex_len = *ey_len = recv_rows;
            /* Attention the offsets */
            //jxf_assert(nz_srt == 0);
            ser_A = jxf_CSRMatrixCreate(num_rows+recv_rows, num_cols, jxf_CSRMatrixNumNonzeros(ser_B)+nz_srt);
            jxf_CSRMatrixInitialize(ser_A);
            IA = jxf_CSRMatrixI(ser_A) + recv_rows;
            JA = jxf_CSRMatrixJ(ser_A) + nz_srt;
            AA = jxf_CSRMatrixData(ser_A) + nz_srt;
            IA[0] = nz_srt;
            rnt = cnt = ctt = 0;
            plle = row_starts[my_id];
            pkle = row_starts[my_id+1];
            for (j = 0; j < num_rows; j ++) // No connection with LEFT or DOWN
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= plle) && (col < pkle)) // its interior, up and right points
                    {
                        JA[ctt] = col;
                    }
                    else // connections with UP and RIGHT only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_largeside; l ++)
                        {
                            pid = largeprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (largeprcpos[l] == 1) // my_id as pid's DOWN
                                {
                                    if (xloc == 0) // left-down corner
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pxln - xloc + 1;
                                        if (xlo_array[pid] != 0) // LEFT exists
                                        {
                                            erow += (pyln - 1);
                                        }
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                                else if (largeprcpos[l] == 3) // my_id as pid's LEFT
                                {
                                    if (yloc == 0) // left-down corner
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pyln - yloc + 1;
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
           *postn_a = *postn_b = *postn_c = cnt;
            //jxf_assert(cnt == num_rows);
           *dx_len = ng_pt_mx + ng_pt_m_xo; // senting start to RIGHT for TYPE 3
           *kx_len = ng_pt_mx; // senting interval to RIGHT for TYPE 3
           *dy_len = num_rows - ng_pt_m_xo; // senting start to UP for TYPE 2
           *lx_len = ng_pt_m_xo; // senting start to RIGHT for TYPE 1
           *ly_len = num_rows - ng_pt_mx; // senting start to UP for TYPE 1
        }
        else // with LEFT
        {
            recv_rows = ng_pt_my;
           *ex_len = *ey_len = recv_rows;
            /* Attention the offsets */
            ser_A = jxf_CSRMatrixCreate(num_rows+recv_rows, num_cols, jxf_CSRMatrixNumNonzeros(ser_B)+nz_srt);
            jxf_CSRMatrixInitialize(ser_A);
            IA = jxf_CSRMatrixI(ser_A) + recv_rows;
            JA = jxf_CSRMatrixJ(ser_A) + nz_srt;
            AA = jxf_CSRMatrixData(ser_A) + nz_srt;
            IA[0] = nz_srt;
            rnt = cnt = ctt = 0;
            pke = row_starts[my_id];
            ple = num_rows - recv_rows - ng_pt;
            ptte = pke + ng_pt;
            plme = num_rows - ng_pt;
            pkle = row_starts[my_id+1];
            psrt = 0;
            for (i = 0; i < y_part_len; i ++)
            {
                jsrt = i * ng_pt_mx;
                jend = jsrt + ng_pt_mx;
                psrt += ng_pt;
                pend = psrt + ng_pt;
                povr = psrt - ng_pt;
                pile = ple - i * ng_pt_m_xo;
                pkke = pke + jsrt;
                plle = pkke + ng_pt;
                psse = pkke + ng_pt_mx;
                for (j = jsrt+ng_pt; j < jend; j ++) // No connection with LEFT or DOWN
                {
                    permute[rnt++] = j;
                    IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                    cnt ++;
                    row_end = IB[j+1];
                    for (k = IB[j]; k < row_end; k ++)
                    {
                        col = JB[k]; // globally
                        if ((col >= pke) && (col < ptte)) // sent its left-down corner points to back, firstly
                        {
                            JA[ctt] = col + plme;
                        }
                        else if ((col >= pkke) && (col < plle)) // sent its left points to back, secondly
                        {
                            JA[ctt] = col + pile;
                        }
                        else if ((col >= plle) && (col < psse)) // its interior points
                        {
                            JA[ctt] = col - psrt;
                        }
                        else if ((col >= psse) && (col < pkle)) // its up points
                        {
                            JA[ctt] = col - pend;
                        }
                        else if ((col >= ptte) && (col < pkke)) // its down points
                        {
                            JA[ctt] = col - povr;
                        }
                        else // connections with UP and RIGHT only, if any. Yue Xiaoqiang on 2014/08/15
                        {
                            // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                            for (l = 0; l < num_largeside; l ++)
                            {
                                pid = largeprocs[l];
                                if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                                {
                                    pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                    pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                    rele = col - row_starts[pid]; // local index
                                    nloc = rele / ng_pt; // local grid-block index
                                    nvar = rele - nloc * ng_pt; // local variable index
                                    yloc = nloc / pxln;
                                    xloc = nloc - yloc * pxln;
                                    if (largeprcpos[l] == 1) // my_id as pid's DOWN
                                    {
                                        if (xloc == 0) // left-down corner
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pxln - xloc + 1;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                    else if (largeprcpos[l] == 3) // my_id as pid's LEFT
                                    {
                                        if (yloc == 0)
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pyln - yloc + 1;
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                }
                            }
                        }
                        AA[ctt++] = BB[k];
                    }
                }
            }
           *postn_a = *postn_b = cnt;
            psrt = ng_pt;
            for (i = 1; i < y_part_len; i ++)
            {
                jsrt = i * ng_pt_mx;
                jend = jsrt + ng_pt;
                psrt += ng_pt;
                pile = ple - i * ng_pt_m_xo;
                pend = pile - ng_pt_m_xo;
                povr = pile + ng_pt_m_xo;
                pkke = pke + jsrt;
                plle = pkke + ng_pt;
                psse = pkke + ng_pt_mx;
                for (j = jsrt; j < jend; j ++) // Connections with LEFT and UP
                {
                    permute[rnt++] = j;
                    IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                    cnt ++;
                    row_end = IB[j+1];
                    for (k = IB[j]; k < row_end; k ++)
                    {
                        col = JB[k]; // globally
                        if ((col >= pke) && (col < ptte)) // its left-down corner points
                        {
                            JA[ctt] = col + plme;
                        }
                        else if ((col >= pkke) && (col < plle)) // its interior points
                        {
                            JA[ctt] = col + pile;
                        }
                        else if ((col >= plle) && (col < psse)) // its right points
                        {
                            JA[ctt] = col - psrt;
                        }
                        else if ((col >= psse) && (col < pkle)) // its up points
                        {
                            JA[ctt] = col + pend;
                        }
                        else if ((col >= ptte) && (col < pkke)) // its down points
                        {
                            JA[ctt] = col + povr;
                        }
                        else // connections with LEFT and UP only, if any. Yue Xiaoqiang on 2014/08/15
                        {
                            // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                            for (l = 0; l < num_sideprocs; l ++)
                            {
                                pid = smallprocs[l];
                                if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                                {
                                    pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                    pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                    rele = col - row_starts[pid]; // local index
                                    nloc = rele / ng_pt; // local grid-block index
                                    nvar = rele - nloc * ng_pt; // local variable index
                                    yloc = nloc / pxln;
                                    xloc = nloc - yloc * pxln;
                                    if (smallprcpos[l] == 1) // my_id as pid's DOWN
                                    {
                                        if (xloc == 0) // left-down corner
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pxln - xloc + 1;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                    else if (smallprcpos[l] == 4) // my_id as pid's RIGHT
                                    {
                                        if (ylo_array[pid] == 0) // no DOWN
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                        else // DOWN exists
                                        {
                                            if (yloc == 0) // bottom-block of grid
                                            {
                                                erow = 2;
                                                if (xlo_array[pid] != 0) // LEFT exists
                                                {
                                                    erow += (pyln - 1);
                                                }
                                                JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                                break;
                                            }
                                            else
                                            {
                                                erow = nloc - pxln;
                                                if (xlo_array[pid] != 0) // LEFT exists
                                                {
                                                    erow -= yloc;
                                                }
                                                JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        AA[ctt++] = BB[k];
                    }
                }
            }
           *postn_c = cnt;
            pend = ple - ng_pt_m_xo;
            psse = pke + ng_pt_mx;
            for (j = 0; j < ng_pt; j ++) // Connections with LEFT only
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= pke) && (col < ptte)) // its interior points
                    {
                        JA[ctt] = col + plme;
                    }
                    else if ((col >= ptte) && (col < psse)) // its right points
                    {
                        JA[ctt] = col - ng_pt;
                    }
                    else if ((col >= psse) && (col < pkle)) // its up points
                    {
                        JA[ctt] = col + pend;
                    }
                    else // connections with LEFT only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_smallside; l ++)
                        {
                            pid = smallprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (smallprcpos[l] == 4) // my_id as pid's RIGHT
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        erow = nloc;
                                        if (xlo_array[pid] != 0) // LEFT exists
                                        {
                                            erow -= (yloc + 1);
                                        }
                                        JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                        break;
                                    }
                                    else // DOWN exists
                                    {
                                        if (yloc == 0) // bottom-block of grid
                                        {
                                            erow = 2;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                            JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                            break;
                                        }
                                        else
                                        {
                                            erow = nloc - pxln;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
            //jxf_assert(cnt == num_rows);
           *dx_len = 2 * ng_pt_m_xo - ng_pt; // senting start to RIGHT for TYPE 3
           *kx_len = ng_pt_m_xo; // senting interval to RIGHT for TYPE 3
           *dy_len = num_rows - ng_pt_m_xo - ng_pt_my; // senting start to UP for TYPE 2
           *lx_len = ng_pt_m_xo - ng_pt; // senting start to RIGHT for TYPE 1
           *ly_len = num_rows - 2 * ng_pt; // senting start to UP for TYPE 1
        }
    }
    else // with DOWN
    {
        if (grid_data->x_lower_idx == 0) // without LEFT
        {
            recv_rows = ng_pt_mx;
           *ex_len = *ey_len = recv_rows;
            /* Attention the offsets */
            ser_A = jxf_CSRMatrixCreate(num_rows+recv_rows, num_cols, jxf_CSRMatrixNumNonzeros(ser_B)+nz_srt);
            jxf_CSRMatrixInitialize(ser_A);
            IA = jxf_CSRMatrixI(ser_A) + recv_rows;
            JA = jxf_CSRMatrixJ(ser_A) + nz_srt;
            AA = jxf_CSRMatrixData(ser_A) + nz_srt;
            IA[0] = nz_srt;
            rnt = cnt = ctt = 0;
            pkke = row_starts[my_id];
            pjje = pkke + ng_pt;
            plle = pkke + recv_rows;
            ptte = num_rows - ng_pt;
            pile = ptte - recv_rows;
            pkle = row_starts[my_id+1];
            for (j = recv_rows; j < num_rows; j ++) // No connection with DOWN
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= pkke) && (col < pjje)) // sent its left-down corner points to back
                    {
                        JA[ctt] = col + ptte;
                    }
                    else if ((col >= pjje) && (col < plle)) // sent its down points to back
                    {
                        JA[ctt] = col + pile;
                    }
                    else if ((col >= plle) && (col < pkle)) // its interior, up and right points
                    {
                        JA[ctt] = col - recv_rows;
                    }
                    else // connections with UP and RIGHT only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_largeside; l ++)
                        {
                            pid = largeprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (largeprcpos[l] == 1) // my_id as pid's DOWN
                                {
                                    if (xloc == 0) // left-down corner
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pxln - xloc + 1;
                                        if (xlo_array[pid] != 0) // LEFT exists
                                        {
                                            erow += (pyln - 1);
                                        }
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                                else if (largeprcpos[l] == 3) // my_id as pid's LEFT
                                {
                                    if (yloc == 0)
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pyln - yloc + 1;
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
           *postn_a = cnt;
            for (j = ng_pt; j < recv_rows; j ++) // Connections with DOWN and RIGHT
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= pkke) && (col < pjje)) // its left-down corner points
                    {
                        JA[ctt] = col + ptte;
                    }
                    else if ((col >= pjje) && (col < plle)) // its interior points
                    {
                        JA[ctt] = col + pile;
                    }
                    else if ((col >= plle) && (col < pkle)) // its up points
                    {
                        JA[ctt] = col - recv_rows;
                    }
                    else // connections with RIGHT and DOWN only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_sideprocs; l ++)
                        {
                            pid = smallprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (smallprcpos[l] == 3) // my_id as pid's LEFT
                                {
                                    if (yloc == 0)
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pyln - yloc + 1;
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                                else if (smallprcpos[l] == 2) // my_id as pid's UP
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_mx + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                    else // DOWN exists
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_m2x + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc - ng_pt_mx;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
           *postn_b = *postn_c = cnt;
            for (j = 0; j < ng_pt; j ++) // Connections with DOWN only
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= pkke) && (col < pjje)) // its interior points
                    {
                        JA[ctt] = col + ptte;
                    }
                    else if ((col >= pjje) && (col < plle)) // its right points
                    {
                        JA[ctt] = col + pile;
                    }
                    else if ((col >= plle) && (col < pkle)) // its up points
                    {
                        JA[ctt] = col - recv_rows;
                    }
                    else // connections with DOWN only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_smallside; l ++)
                        {
                            pid = smallprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (smallprcpos[l] == 2) // my_id as pid's UP
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_mx + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                    else // DOWN exists
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_m2x + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc - ng_pt_mx;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
            //jxf_assert(cnt == num_rows);
           *dx_len = ng_pt_m_xo; // senting start to RIGHT for TYPE 3
           *kx_len = ng_pt_mx; // senting interval to RIGHT for TYPE 3
           *dy_len = num_rows - ng_pt_m_xo - ng_pt_mx; // senting start to UP for TYPE 2
           *lx_len = num_rows - 2 * ng_pt; // senting start to RIGHT for TYPE 1
           *ly_len = num_rows - 2 * ng_pt_mx; // senting start to UP for TYPE 1
        }
        else // with LEFT
        {
            recv_rows = ng_pt_mx + ng_pt_my;
           *ex_len = recv_rows;
           *ey_len = recv_rows - ng_pt; // YUE: exclude the reduplicate block
            /* Attention the offsets */
            ser_A = jxf_CSRMatrixCreate(num_rows+recv_rows, num_cols, jxf_CSRMatrixNumNonzeros(ser_B)+nz_srt);
            jxf_CSRMatrixInitialize(ser_A);
            IA = jxf_CSRMatrixI(ser_A) + recv_rows;
            JA = jxf_CSRMatrixJ(ser_A) + nz_srt;
            AA = jxf_CSRMatrixData(ser_A) + nz_srt;
            IA[0] = nz_srt;
            rnt = cnt = ctt = 0;
            pke = row_starts[my_id];
            ple = num_rows - recv_rows;
            pile = ple + ng_pt_m_xo;
            plme = num_rows - ng_pt;
            pkle = row_starts[my_id+1];
            psrt = ng_pt_mx;
            puue = pke + ng_pt;
            ptte = pke + ng_pt_mx;
            for (i = 1; i < y_part_len; i ++)
            {
                jsrt = i * ng_pt_mx;
                jend = jsrt + ng_pt_mx;
                psrt += ng_pt;
                pend = psrt + ng_pt;
                povr = psrt - ng_pt;
                pile -= ng_pt_m_xo;
                pkke = pke + jsrt;
                plle = pkke + ng_pt;
                psse = pkke + ng_pt_mx;
                for (j = jsrt+ng_pt; j < jend; j ++) // No connection with DOWN or LEFT
                {
                    permute[rnt++] = j;
                    IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                    cnt ++;
                    row_end = IB[j+1];
                    for (k = IB[j]; k < row_end; k ++)
                    {
                        col = JB[k]; // globally
                        if ((col >= pkke) && (col < plle)) // sent its left points to back, without LEFT
                        {
                            JA[ctt] = col + pile;
                        }
                        else if ((col >= plle) && (col < psse)) // its interior points
                        {
                            JA[ctt] = col - psrt;
                        }
                        else if ((col >= psse) && (col < pkle)) // its up points
                        {
                            JA[ctt] = col - pend;
                        }
                        else if ((col >= puue) && (col < ptte)) // its downest points
                        {
                            JA[ctt] = col + ple;
                        }
                        else if ((col >= ptte) && (col < pkke)) // its down points
                        {
                            JA[ctt] = col - povr;
                        }
                        else // connections with UP and RIGHT only, if any. Yue Xiaoqiang on 2014/08/15
                        {
                            // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                            for (l = 0; l < num_largeside; l ++)
                            {
                                pid = largeprocs[l];
                                if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                                {
                                    pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                    pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                    rele = col - row_starts[pid]; // local index
                                    nloc = rele / ng_pt; // local grid-block index
                                    nvar = rele - nloc * ng_pt; // local variable index
                                    yloc = nloc / pxln;
                                    xloc = nloc - yloc * pxln;
                                    if (largeprcpos[l] == 1) // my_id as pid's DOWN
                                    {
                                        if (xloc == 0) // left-down corner
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pxln - xloc + 1;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                    else if (largeprcpos[l] == 3) // my_id as pid's LEFT
                                    {
                                        if (yloc == 0)
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pyln - yloc + 1;
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                }
                            }
                        }
                        AA[ctt++] = BB[k];
                    }
                }
            }
           *postn_a = cnt;
            for (j = ng_pt; j < ng_pt_mx; j ++) // Connections with DOWN and RIGHT
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= puue) && (col < ptte)) // its interior points
                    {
                        JA[ctt] = col + ple;
                    }
                    else if ((col >= pke) && (col < puue)) // its leftest points
                    {
                        JA[ctt] = col + plme;
                    }
                    else if ((col >= ptte) && (col < pkle)) // its up points
                    {
                        JA[ctt] = col - ng_pt_m_px;
                    }
                    else // connections with DOWN and RIGHT only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_sideprocs; l ++)
                        {
                            pid = smallprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (smallprcpos[l] == 3) // my_id as pid's LEFT
                                {
                                    if (yloc == 0)
                                    {
                                        erow = 1;
                                    }
                                    else
                                    {
                                        erow = pyln - yloc + 1;
                                    }
                                    JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                    break;
                                }
                                else if (smallprcpos[l] == 2) // my_id as pid's UP
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_mx + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                    else // DOWN exists
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_m2x + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc - ng_pt_mx;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
           *postn_b = cnt;
            pile = ple + ng_pt_m_xo;
            psrt = ng_pt_mx;
            for (i = 1; i < y_part_len; i ++)
            {
                jsrt = i * ng_pt_mx;
                jend = jsrt + ng_pt;
                psrt += ng_pt;
                pile -= ng_pt_m_xo;
                pend = pile - ng_pt_m_xo;
                povr = pile + ng_pt_m_xo;
                pkke = pke + jsrt;
                plle = pkke + ng_pt;
                psse = pkke + ng_pt_mx;
                for (j = jsrt; j < jend; j ++) // Connections with UP and LEFT
                {
                    permute[rnt++] = j;
                    IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                    cnt ++;
                    row_end = IB[j+1];
                    for (k = IB[j]; k < row_end; k ++)
                    {
                        col = JB[k]; // globally
                        if ((col >= pkke) && (col < plle)) // its interior points
                        {
                            JA[ctt] = col + pile;
                        }
                        else if ((col >= plle) && (col < psse)) // its right points
                        {
                            JA[ctt] = col - psrt;
                        }
                        else if ((col >= psse) && (col < pkle)) // its up points
                        {
                            JA[ctt] = col + pend;
                        }
                        else if ((col >= pke) && (col < puue)) // its downest points
                        {
                            JA[ctt] = col + plme;
                        }
                        else if ((col >= puue) && (col < pkke)) // its down points
                        {
                            JA[ctt] = col + povr;
                        }
                        else // connections with UP and LEFT only, if any. Yue Xiaoqiang on 2014/08/15
                        {
                            // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                            for (l = 0; l < num_sideprocs; l ++)
                            {
                                pid = smallprocs[l];
                                if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                                {
                                    pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                    pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                    rele = col - row_starts[pid]; // local index
                                    nloc = rele / ng_pt; // local grid-block index
                                    nvar = rele - nloc * ng_pt; // local variable index
                                    yloc = nloc / pxln;
                                    xloc = nloc - yloc * pxln;
                                    if (smallprcpos[l] == 1) // my_id as pid's DOWN
                                    {
                                        if (xloc == 0) // left-down corner
                                        {
                                            erow = 1;
                                        }
                                        else
                                        {
                                            erow = pxln - xloc + 1;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                        }
                                        JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                        break;
                                    }
                                    else if (smallprcpos[l] == 4) // my_id as pid's RIGHT
                                    {
                                        if (ylo_array[pid] == 0) // no DOWN
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                        else // DOWN exists
                                        {
                                            if (yloc == 0) // bottom-block of grid
                                            {
                                                erow = 2;
                                                if (xlo_array[pid] != 0) // LEFT exists
                                                {
                                                    erow += (pyln - 1);
                                                }
                                                JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                                break;
                                            }
                                            else
                                            {
                                                erow = nloc - pxln;
                                                if (xlo_array[pid] != 0) // LEFT exists
                                                {
                                                    erow -= yloc;
                                                }
                                                JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        AA[ctt++] = BB[k];
                    }
                }
            }
           *postn_c = cnt;
            for (j = 0; j < ng_pt; j ++) // Connections with DOWN and LEFT
            {
                permute[rnt++] = j;
                IA[cnt+1] = IA[cnt] + (IB[j+1] - IB[j]);
                cnt ++;
                row_end = IB[j+1];
                for (k = IB[j]; k < row_end; k ++)
                {
                    col = JB[k]; // globally
                    if ((col >= puue) && (col < ptte)) // its right points
                    {
                        JA[ctt] = col + ple;
                    }
                    else if ((col >= pke) && (col < puue)) // its interior points
                    {
                        JA[ctt] = col + plme;
                    }
                    else if ((col >= ptte) && (col < pkle)) // its up points
                    {
                        JA[ctt] = col + ple;
                    }
                    else // connections with DOWN and LEFT only, if any. Yue Xiaoqiang on 2014/08/15
                    {
                        // Loop for the wanted PID, Yue Xiaoqiang on 2014/08/15
                        for (l = 0; l < num_smallside; l ++)
                        {
                            pid = smallprocs[l];
                            if ((col >= row_starts[pid]) && (col < row_starts[pid+1]))
                            {
                                pxln = xup_array[pid] - xlo_array[pid]; // local number of x points
                                pyln = yup_array[pid] - ylo_array[pid]; // local number of y points
                                rele = col - row_starts[pid]; // local index
                                nloc = rele / ng_pt; // local grid-block index
                                nvar = rele - nloc * ng_pt; // local variable index
                                yloc = nloc / pxln;
                                xloc = nloc - yloc * pxln;
                                if (smallprcpos[l] == 2) // my_id as pid's UP
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_mx + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= (yloc + 1);
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                    else // DOWN exists
                                    {
                                        if (xloc == 0) // leftest-block
                                        {
                                            if (xlo_array[pid] == 0) // no LEFT
                                            {
                                                JA[ctt] = row_starts[pid+1] - ng_pt_m2x + nvar;
                                                break;
                                            }
                                            else // LEFT exists
                                            {
                                                JA[ctt] = row_starts[pid+1] - 2 * ng_pt + nvar;
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            erow = nloc - ng_pt_mx;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                                else if (smallprcpos[l] == 4) // my_id as pid's RIGHT
                                {
                                    if (ylo_array[pid] == 0) // no DOWN
                                    {
                                        erow = nloc;
                                        if (xlo_array[pid] != 0) // LEFT exists
                                        {
                                            erow -= (yloc + 1);
                                        }
                                        JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                        break;
                                    }
                                    else // DOWN exists
                                    {
                                        if (yloc == 0) // bottom-block of grid
                                        {
                                            erow = 2;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow += (pyln - 1);
                                            }
                                            JA[ctt] = row_starts[pid+1] - ng_pt * erow + nvar;
                                            break;
                                        }
                                        else
                                        {
                                            erow = nloc - pxln;
                                            if (xlo_array[pid] != 0) // LEFT exists
                                            {
                                                erow -= yloc;
                                            }
                                            JA[ctt] = row_starts[pid] + ng_pt * erow + nvar;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    AA[ctt++] = BB[k];
                }
            }
            //jxf_assert(cnt == num_rows);
           *dx_len = ng_pt_m_xo - ng_pt; // senting start to RIGHT for TYPE 3
           *kx_len = ng_pt_m_xo; // senting interval to RIGHT for TYPE 3
           *dy_len = num_rows - 2 * ng_pt_m_xo - ng_pt_my; // senting start to UP for TYPE 2
           *lx_len = num_rows - ng_pt_m_px; // senting start to RIGHT for TYPE 1
           *ly_len = num_rows - 2 * ng_pt; // senting start to UP for TYPE 1
        }
    }
    
    return ser_A;
}
