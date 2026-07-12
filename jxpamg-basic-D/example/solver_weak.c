//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * solver_weak.c -- This is a program to test the parallel strategy for 
 * solving multi linear systems with the same sparse pattern.
 *
 * One of the following solvers can be used by
 * assigning parameter 'solver_id'
 *  solver_id = 0 : AMG
 *  solver_id = 11: CG
 *  solver_id = 12: AMG-CG
 *  solver_id = 13: DS-CG
 *  solver_id = 21: GMRES
 *  solver_id = 22: AMG-GMRES
 *  solver_id = 23: DS-GMRES
 *  solver_id = 31: BICGSTAB
 *  solver_id = 32: AMG-BICGSTAB
 *  solver_id = 33: DS-BICGSTAB
 *
 * The input data for the linear system to be solved should be
 * provided in the form as follows (the indices are of Fortran
 * style, i.e, starting from 1 not 0):
 *  1. matrix file
 *----------------------------------------
 *     n
 *     ia(i), i = 1(1)n+1
 *     ja(i), i = 1(1)nz
 *      a(i), i = 1(1)nz
 *----------------------------------------
 *  where nz = ia(n+1)-1.
 *  2. rhs file
 *----------------------------------------
 *     n
 *     f(i),i = 1(1)n
 *----------------------------------------
 *
 *  Created by peghoty  2010/12/02
 *  Xiangtan University
 *  peghoty@163.com
 *
 */

#include "jx_pamg.h"
#include "jx_krylov.h"
#include "jx_multils.h"

JX_Int
jx_ParMultiLSSolve ( JX_Int               global_num_ls,
                     jx_ParCSRMatrix **par_mat_array,
                     jx_ParVector    **par_rhs_array,
                     jx_ParVector    **par_app_array,
                     fsls_SolverParam *solver        );
fsls_ParLSData *
jx_ParMatVec2LSData( JX_Int               global_num_ls,
                     jx_ParCSRMatrix **par_mat_array,
                     jx_ParVector    **par_rhs_array,
                     jx_ParVector    **par_app_array );
fsls_CSRMatrix *
fsls_MatrixTransfer( jx_CSRMatrix *A );
fsls_Vector *
fsls_VectorTransfer( jx_Vector *x );
fsls_ParGLSData *
jx_ParMatVec2GLSData( JX_Int               global_num_ls,
                      jx_ParCSRMatrix **par_mat_array,
                      jx_ParVector    **par_rhs_array,
                      jx_ParVector    **par_app_array );
fsls_CSRMatrix *
jx_CombineDiagOffd( MPI_Comm      comm,
                    JX_Int           num_rows,
                    JX_Int           num_cols,
                    JX_Int           first_row_index,
                    jx_CSRMatrix *diag,
                    jx_CSRMatrix *offd,
                    JX_Int          *col_map_offd );                      
JX_Int
fsls_ParAppBack( fsls_ParLSData   *ls_data, 
                 jx_ParVector    **par_app_array );

int
main( int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   JX_Int myid, nprocs;
   JX_Real starttime, endtime;

   JX_Int arg_index = 0;
   JX_Int print_usage = 0;

   char **MatFile_array = NULL;
   char **RhsFile_array = NULL;

   jx_ParCSRMatrix **par_mat_array = NULL;
   jx_ParVector **par_rhs_array = NULL;
   jx_ParVector **par_app_array = NULL;

   fsls_SolverParam *solver = NULL;

   JX_Int i, LL;
   JX_Int *partitioning = NULL;

   JX_Int problem_id;
   JX_Int global_num_ls;
   JX_Int solver_id;      // solver id
   JX_Real tol;         // the tolerance
   JX_Int max_iter;       // the maximal number of iterations
   JX_Int min_iter;       // the minimal number of iterations
   JX_Int two_norm;       // L^2 norm(1) or preconditioner-norm(0)? used for PCG
   JX_Int k_dim;          // the dimension of krylov subspace, used for PGMRES
   JX_Int print_level;    // how much information to be printed
   JX_Int file_base;

   jx_MPI_Init(&argc, &argv);
   fsls_MPICommInformation(comm, &myid, &nprocs);

   global_num_ls = 9;
   problem_id    = 21;
   solver_id     = 12;
   tol           = 1.0e-8;
   max_iter      = 2000;
   min_iter      = 0;
   two_norm      = 0;
   k_dim         = 10;
   print_level   = 0;
   file_base     = 1;

   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-sid") == 0 )
      {
         arg_index ++;
         solver_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ls") == 0 )
      {
         arg_index ++;
         global_num_ls = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-pid") == 0 )
      {
         arg_index ++;
         problem_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage = 1;
         break;
      }
      else
      {
         arg_index ++;
      }
   }

   if (print_usage)
   {
      jx_printf("\n");
      jx_printf("  Usage: %s [<options>]\n", argv[0]);
      jx_printf("\n");
      jx_printf("  -ls      <val> : global_num_ls\n");
      jx_printf("  -sid     <val> : solver id\n");
      jx_printf("  -pid     <val> : problem id\n");
      jx_printf("  -help          : using help message\n\n");
      exit(1);
   }

   MatFile_array = fsls_CTAlloc(char *, global_num_ls);
   RhsFile_array = fsls_CTAlloc(char *, global_num_ls);

   //-----------------------------------------------------------
   //  Step 1: Provide the data file of all the linear systems
   //-----------------------------------------------------------
   switch (problem_id)
   {
      case 21:
      {
         for (i = 0; i < global_num_ls; i ++)
         {
            MatFile_array[i] = "./input/mat_20X20";
            RhsFile_array[i] = "./input/rhs_20X20";
         }
      }
      break;
   }

   //------------------------------------------------------------------------------
   //  Step 2: Build parallel matrices and rhs vectors by reading data from files
   //------------------------------------------------------------------------------
   par_mat_array = jx_CTAlloc(jx_ParCSRMatrix *, global_num_ls);
   par_rhs_array = jx_CTAlloc(jx_ParVector *, global_num_ls);
   par_app_array = jx_CTAlloc(jx_ParVector *, global_num_ls);
   for (i = 0; i < global_num_ls; i ++)
   {
      par_mat_array[i] = jx_BuildMatParFromOneFile(MatFile_array[i], 1, file_base);
      par_rhs_array[i] = jx_BuildRhsParFromOneFile(RhsFile_array[i], par_mat_array[i]);
   }

   //-------------------------------
   //  Step 3: Set initial guesses
   //-------------------------------
   LL = jx_ParVectorGlobalSize(par_rhs_array[0]);
   partitioning = jx_ParVectorPartitioning(par_rhs_array[0]);
   for (i = 0; i < global_num_ls; i ++)
   {
      par_app_array[i] = jx_ParVectorCreate(comm, LL, partitioning);
      jx_ParVectorSetPartitioningOwner(par_app_array[i], 0);
      jx_ParVectorInitialize(par_app_array[i]);
      jx_ParVectorSetConstantValues(par_app_array[i], 0.0);
   }

   //------------------------------------------------------------
   //  Step 4: Solve the global_num_ls linear system parallelly
   //------------------------------------------------------------
   starttime = jx_MPI_Wtime();

   /* create a solver parameter object */
   solver = fsls_SolverParamCreate();
   fsls_SolverParamSetSolverID(solver, solver_id);
   fsls_SolverParamSetTolerance(solver, tol);
   fsls_SolverParamSetMaxIter(solver, max_iter);
   fsls_SolverParamSetMinIter(solver, min_iter);
   fsls_SolverParamSetTwoNorm(solver, two_norm);
   fsls_SolverParamSetKDim(solver, k_dim);
   fsls_SolverParamSetPrintLevel(solver, print_level);

   jx_ParMultiLSSolve(global_num_ls, par_mat_array, par_rhs_array, par_app_array, solver);

   endtime = jx_MPI_Wtime();
   fsls_GetWallTime(comm, "JX_ParMultiLSSolve", starttime, endtime, 0, DECIMALPlACE);

   //---------------------------
   //  Step 5: Free some staff
   //---------------------------
   fsls_TFree(solver);
   fsls_TFree(MatFile_array);
   fsls_TFree(RhsFile_array);
   for (i = 0; i < global_num_ls; i ++)
   {
      jx_ParCSRMatrixDestroy(par_mat_array[i]);
      jx_ParVectorDestroy(par_rhs_array[i]);
      jx_ParVectorDestroy(par_app_array[i]);
   }

   jx_MPI_Finalize();

   return 0;
}

JX_Int
jx_ParMultiLSSolve( JX_Int               global_num_ls,
                    jx_ParCSRMatrix **par_mat_array,
                    jx_ParVector    **par_rhs_array,
                    jx_ParVector    **par_app_array,
                    fsls_SolverParam *solver        )
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_mat_array[0]);

   fsls_ParLSData   *ls_data  = NULL;
   fsls_ParGLSData  *gls_data = NULL;

   JX_Real starttime, endtime;

   JX_Int myid, nprocs;

   fsls_MPICommInformation(comm, &myid, &nprocs);

   if (global_num_ls < nprocs)
   {
      if (myid == 0)
      {
         jx_printf("\n >>> Warning: global_num_ls should be greater than or equal to nprocs!!\n\n");
         return -1;
      }
   }
   else
   {
      if (1)
      {
         //===========================
         //  Step 1: par -> gls_data
         //===========================
         starttime = jx_MPI_Wtime();
         gls_data = jx_ParMatVec2GLSData(global_num_ls, par_mat_array,
                                         par_rhs_array, par_app_array);
         endtime = jx_MPI_Wtime();
         fsls_GetWallTime(comm, "jx_ParMatVec2GLSData", starttime, endtime, 0, DECIMALPlACE);

         //===============================
         //  Step 2: gls_data -> ls_data
         //===============================
         starttime = jx_MPI_Wtime();
         ls_data = fsls_ParGLS2LSData(gls_data);
         endtime = jx_MPI_Wtime();
         fsls_GetWallTime(comm, "fsls_ParGLS2LSData", starttime, endtime, 0, DECIMALPlACE);

         fsls_ParGLSDataDestroy(gls_data);
      }
      else // this is just for comparison
      {
         ls_data = jx_ParMatVec2LSData(global_num_ls, par_mat_array, 
                                       par_rhs_array, par_app_array);
      }

      //=======================
      //  Step 3: Solve phase
      //=======================
      starttime = jx_MPI_Wtime();
      fsls_ParallelSolve(ls_data, solver);
      endtime = jx_MPI_Wtime();
      fsls_GetWallTime(comm, "fsls_ParallelSolve", starttime, endtime, 0, DECIMALPlACE);

      //=========================================
      //  Step 4: Return back the approximation
      //=========================================
      starttime = jx_MPI_Wtime();
      fsls_ParAppBack(ls_data, par_app_array);
      endtime = jx_MPI_Wtime();
      fsls_GetWallTime(comm, "fsls_ParAppBack", starttime, endtime, 0, DECIMALPlACE);

      //===========================
      //  Step 5: Free some stuff
      //===========================
      fsls_ParLSDataDestroy(ls_data);
   }

   return 0;
}

JX_Int
fsls_ParAppBack( fsls_ParLSData   *ls_data,
                 jx_ParVector    **par_app_array )
{
   MPI_Comm comm = fsls_ParLSDataComm(ls_data);
   JX_Int local_num_ls = fsls_ParLSDataLocalNumLS(ls_data);
   JX_Int global_num_ls = fsls_ParLSDataGlobalNumLS(ls_data);
   JX_Int global_size = fsls_ParLSDataNumRows(ls_data);
   JX_Int local_size = jx_VectorSize(jx_ParVectorLocalVector(par_app_array[0]));
   JX_Int *ls_partition = fsls_ParLSDataLSPartition(ls_data);
   JX_Int *num_ls_procs = fsls_ParLSDataNumLSProcs(ls_data);
   fsls_Vector **x_array = fsls_ParLSDataXArray(ls_data);
   JX_Int *dof_partition = jx_ParVectorPartitioning(par_app_array[0]);

   JX_Int app_recv_size;
   JX_Int app_send_size;

   JX_Real *app_recvbuf = NULL;
   JX_Real *app_sendbuf = NULL;
   JX_Real *app_commbuf = NULL;
   JX_Real *my_app_data = NULL;
   JX_Real *app_data    = NULL;

   MPI_Status *Status = NULL;
   MPI_Request *Requests = NULL;

   JX_Int i, j, k, m;
   JX_Int begin;
   JX_Int length;
   JX_Int num_requests;
   JX_Int myid, nprocs;

   fsls_MPICommInformation(comm, &myid, &nprocs);

   //-----------------------------------------
   // Allocate memory for send- and recv-buf
   //-----------------------------------------
   num_requests = 2*(nprocs-1);
   Requests = fsls_CTAlloc(MPI_Request, num_requests);
   Status   = fsls_CTAlloc(MPI_Status, num_requests);
   app_send_size = (global_size - local_size) * local_num_ls;
   app_recv_size = (global_num_ls - local_num_ls) * local_size;
   app_commbuf = fsls_CTAlloc(JX_Real, app_send_size + app_recv_size);
   app_recvbuf = app_commbuf;
   app_sendbuf = app_recvbuf + app_recv_size;

   //------------------------
   //  Fill the app_sendbuf
   //------------------------
   for (i = 0, k = 0; i < nprocs; i ++)
   {
      if (i != myid) // send data to other processors
      {
         begin  = dof_partition[i];
         length = dof_partition[i+1] - begin;
         for (j = 0; j < num_ls_procs[myid]; j ++)
         {
            app_data = fsls_VectorData(x_array[j]);
            memcpy(&app_sendbuf[k], &app_data[begin], length*sizeof(JX_Real));
            k += length;
         }
      }
   }

   //----------------
   //  Receive data
   //----------------
   for (i = 0, k = 0, m = 0; i < nprocs; i ++)
   {
      if (i != myid) // receive data from the i-th processors
      {
         length = local_size*(ls_partition[i+1] - ls_partition[i]);
         jx_MPI_Irecv(&app_recvbuf[k], length, MPI_DOUBLE, i, myid*321, comm, &Requests[m++]);
         k += length;
      }
   }

   //-------------
   //  Send data
   //-------------
   for (i = 0, k = 0; i < nprocs; i ++)
   {
      if (i != myid)  // send data to the i-th processor
      {
         length = local_num_ls*(dof_partition[i+1] - dof_partition[i]);
         jx_MPI_Isend(&app_sendbuf[k], length, MPI_DOUBLE, i, i*321, comm, &Requests[m++]);
         k += length;
      }
   }

   jx_MPI_Waitall(num_requests, Requests, Status);

   //--------------------------------------------------
   //  Copy the sub-app data in the current processor
   //--------------------------------------------------
   begin = dof_partition[myid];
   for (i = 0; i < local_num_ls; i ++)
   {
      begin = dof_partition[myid];
      length = dof_partition[myid+1] - begin;

      my_app_data = fsls_VectorData(x_array[i]);
      app_data    = jx_VectorData(jx_ParVectorLocalVector(par_app_array[ls_partition[myid]+i]));

      memcpy(app_data, &my_app_data[begin], length*sizeof(JX_Real));
   }

   //-----------------------------------------------
   //  Fill the sub-app data from other processors
   //-----------------------------------------------
   for (i = 0, k = 0; i < nprocs; i ++)
   {
      if (i != myid)
      {
         for (j = 0; j < num_ls_procs[i]; j ++)
         {
            app_data = jx_VectorData(jx_ParVectorLocalVector(par_app_array[ls_partition[i]+j]));
            memcpy(app_data, &app_recvbuf[k], local_size*sizeof(JX_Real));
            k += local_size;
         }
      }
   }

   //-------------------
   //  Free some staff
   //-------------------
   fsls_TFree(Requests);
   fsls_TFree(Status);
   fsls_TFree(app_commbuf);

   return 0;
}

fsls_ParLSData *
jx_ParMatVec2LSData( JX_Int               global_num_ls,
                     jx_ParCSRMatrix **par_mat_array,
                     jx_ParVector    **par_rhs_array,
                     jx_ParVector    **par_app_array )
{
   fsls_ParLSData *ls_data = NULL;

   MPI_Comm comm = jx_ParCSRMatrixComm(par_mat_array[0]);

   jx_CSRMatrix **AA_array = NULL;
   jx_Vector **bb_array = NULL;
   jx_Vector **xx_array = NULL;

   fsls_CSRMatrix **A_array = NULL;
   fsls_Vector **b_array = NULL;
   fsls_Vector **x_array = NULL;

   JX_Int  local_num_ls;
   JX_Int *ls_partition = NULL;

   JX_Int i;
   JX_Int num_rows, num_cols, num_nonzeros;
   JX_Int myid, nprocs;

   fsls_MPICommInformation(comm, &myid, &nprocs);

   AA_array = jx_CTAlloc(jx_CSRMatrix *, global_num_ls);
   bb_array = jx_CTAlloc(jx_Vector *, global_num_ls);
   xx_array = jx_CTAlloc(jx_Vector *, global_num_ls);

   A_array = jx_CTAlloc(fsls_CSRMatrix *, global_num_ls);
   b_array = jx_CTAlloc(fsls_Vector *, global_num_ls);
   x_array = jx_CTAlloc(fsls_Vector *, global_num_ls);

   for (i = 0; i < global_num_ls; i ++)
   {
      AA_array[i] = jx_ParCSRMatrixToCSRMatrixAll(par_mat_array[i]);
      bb_array[i] = jx_ParVectorToVectorAll(par_rhs_array[i]);
      xx_array[i] = jx_ParVectorToVectorAll(par_app_array[i]);

      A_array[i]  = fsls_MatrixTransfer(AA_array[i]);
      b_array[i]  = fsls_VectorTransfer(bb_array[i]);
      x_array[i]  = fsls_VectorTransfer(xx_array[i]);
   }

   num_rows = jx_CSRMatrixNumRows(AA_array[0]);
   num_cols = jx_CSRMatrixNumCols(AA_array[0]);
   num_nonzeros = jx_CSRMatrixNumNonzeros(AA_array[0]);

   //=======================================================================
   //  Generate 'ls_data' in each processor.
   //======================================================================= 
   ls_data = fsls_ParLSDataCreate(comm, num_rows, num_cols, num_nonzeros, global_num_ls, NULL, NULL);
   local_num_ls = fsls_ParLSDataLocalNumLS(ls_data);
   ls_partition = fsls_ParLSDataLSPartition(ls_data);
   fsls_ParLSDataAArray(ls_data) = fsls_CTAlloc(fsls_CSRMatrix *, local_num_ls);
   fsls_ParLSDataBArray(ls_data) = fsls_CTAlloc(fsls_Vector *, local_num_ls);
   fsls_ParLSDataXArray(ls_data) = fsls_CTAlloc(fsls_Vector *, local_num_ls);

   for (i = 0; i < local_num_ls; i ++)
   {
      fsls_ParLSDataAArray(ls_data)[i] = A_array[ls_partition[myid]+i];
      fsls_ParLSDataBArray(ls_data)[i] = b_array[ls_partition[myid]+i];
      fsls_ParLSDataXArray(ls_data)[i] = x_array[ls_partition[myid]+i];
   }

   return ls_data;
}

fsls_Vector *
fsls_VectorTransfer( jx_Vector *x )
{
   fsls_Vector *y = NULL;

   JX_Real *x_data = jx_VectorData(x);
   JX_Int     size   = jx_VectorSize(x);

   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);

   jx_SeqVectorDestroy(x);
   memcpy(fsls_VectorData(y), x_data, size*sizeof(JX_Real));

   return y;
}

fsls_CSRMatrix *
fsls_MatrixTransfer( jx_CSRMatrix *A )
{
   fsls_CSRMatrix *B = NULL;

   JX_Real *a = jx_CSRMatrixData(A);
   JX_Int *ia = jx_CSRMatrixI(A);
   JX_Int *ja = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   JX_Int num_nonzeros = jx_CSRMatrixNumNonzeros(A);

   B = fsls_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   fsls_CSRMatrixInitialize(B);

   memcpy(fsls_CSRMatrixI(B), ia, (num_rows+1)*sizeof(JX_Int));
   memcpy(fsls_CSRMatrixJ(B), ja, num_nonzeros*sizeof(JX_Int));
   memcpy(fsls_CSRMatrixData(B), a, num_nonzeros*sizeof(JX_Real));

   jx_CSRMatrixDestroy(A);

   return B;
}

fsls_ParGLSData *
jx_ParMatVec2GLSData( JX_Int               global_num_ls,
                      jx_ParCSRMatrix **par_mat_array,
                      jx_ParVector    **par_rhs_array,
                      jx_ParVector    **par_app_array )
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_mat_array[0]);
   JX_Int  global_num_dof = jx_ParVectorGlobalSize(par_rhs_array[0]);
   JX_Int *vec_partition = jx_ParVectorPartitioning(par_rhs_array[0]);
   JX_Int *dof_partition = NULL;
   JX_Int *num_dof_procs = NULL;

   JX_Int *col_map_offd  = NULL;

   jx_CSRMatrix	*diag = NULL;
   jx_CSRMatrix	*offd = NULL;

   fsls_CSRMatrix **subA_array = NULL;
   fsls_Vector **subb_array = NULL;
   fsls_Vector **subx_array = NULL;

   fsls_ParGLSData *gls_data = NULL;

   JX_Int num_rows;
   JX_Int num_cols;
   JX_Int first_row_index;

   JX_Int i;
   JX_Int myid,nprocs;

   fsls_MPICommInformation(comm, &myid, &nprocs);

   dof_partition = jx_CTAlloc(JX_Int, nprocs+1);
   num_dof_procs = jx_CTAlloc(JX_Int, nprocs);

   memcpy(dof_partition, vec_partition, (nprocs+1)*sizeof(JX_Int));

   for (i = 0; i < nprocs; i ++)
   {
      num_dof_procs[i] = dof_partition[i+1] - dof_partition[i];
   }

   /* create a fsls_ParGLSData structure */
   gls_data = fsls_ParGLSDataCreate(comm, global_num_ls, global_num_dof, dof_partition, num_dof_procs);

   fsls_ParGLSDataInitialize(gls_data);
   subA_array = fsls_ParGLSDataSubAArray(gls_data);
   subb_array = fsls_ParGLSDataSubBArray(gls_data);
   subx_array = fsls_ParGLSDataSubXArray(gls_data);

   num_rows = num_dof_procs[myid];
   num_cols = jx_ParCSRMatrixGlobalNumCols(par_mat_array[0]);

   for (i = 0; i < global_num_ls; i ++)
   {
      diag = jx_ParCSRMatrixDiag(par_mat_array[i]);
      offd = jx_ParCSRMatrixOffd(par_mat_array[i]);
      col_map_offd = jx_ParCSRMatrixColMapOffd(par_mat_array[i]);
      first_row_index = jx_ParCSRMatrixFirstRowIndex(par_mat_array[i]);
      subA_array[i] = jx_CombineDiagOffd(comm, num_rows, num_cols, first_row_index, diag, offd, col_map_offd);
   }

   for (i = 0; i < global_num_ls; i ++)
   {
      subb_array[i] = fsls_SeqVectorCreate(num_rows);
      fsls_SeqVectorInitialize(subb_array[i]);
      memcpy(fsls_VectorData(subb_array[i]),
             jx_VectorData(jx_ParVectorLocalVector(par_rhs_array[i])),
             num_rows*sizeof(JX_Real));
      subx_array[i] = fsls_SeqVectorCreate(num_rows);
      fsls_SeqVectorInitialize(subx_array[i]);
      memcpy(fsls_VectorData(subx_array[i]),
             jx_VectorData( jx_ParVectorLocalVector(par_app_array[i])),
             num_rows*sizeof(JX_Real));
   }

   return gls_data;
}

fsls_CSRMatrix *
jx_CombineDiagOffd( MPI_Comm      comm,
                    JX_Int           num_rows,
                    JX_Int           num_cols,
                    JX_Int           first_row_index,
                    jx_CSRMatrix *diag,
                    jx_CSRMatrix *offd,
                    JX_Int          *col_map_offd )
{
   fsls_CSRMatrix *A = NULL;

   JX_Real *a  = NULL;
   JX_Int *ia = NULL;
   JX_Int *ja = NULL;
   JX_Int num_nonzeros;

   JX_Real *a_diag = jx_CSRMatrixData(diag);
   JX_Int *ia_diag = jx_CSRMatrixI(diag);
   JX_Int *ja_diag = jx_CSRMatrixJ(diag);
   JX_Int nz_diag = jx_CSRMatrixNumNonzeros(diag);

   JX_Real *a_offd = jx_CSRMatrixData(offd);
   JX_Int *ia_offd = jx_CSRMatrixI(offd);
   JX_Int *ja_offd = jx_CSRMatrixJ(offd);
   JX_Int nz_offd = jx_CSRMatrixNumNonzeros(offd);

   JX_Int i, j, cnt;
   JX_Int myid, nprocs;

   fsls_MPICommInformation(comm, &myid, &nprocs);
   num_nonzeros = nz_diag + nz_offd;
   A = fsls_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   fsls_CSRMatrixInitialize(A);
   a  = fsls_CSRMatrixData(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);
   ia[0] = 0; cnt = 0;
   for (i = 0; i < num_rows; i ++)
   {
      for (j = ia_diag[i]; j < ia_diag[i+1]; j ++)
      {
         a[cnt] = a_diag[j];
         ja[cnt] = ja_diag[j] + first_row_index;
         cnt ++;
      }
      for (j = ia_offd[i]; j < ia_offd[i+1]; j ++)
      {
         a[cnt] = a_offd[j];
         ja[cnt] = col_map_offd[ja_offd[j]];
         cnt ++;
      }
      ia[i+1] = cnt;
   }

   return A;
}
