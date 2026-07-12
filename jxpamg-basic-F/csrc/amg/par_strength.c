//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_strength.c -- Generate a strengthen matrix of a parallel CSR matrix.
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGCreateS
 * \brief Generate a strengthen matrix of a given parallel CSR matrix.
 * \date 2011/09/03
 */ 
JXF_Int
jxf_PAMGCreateS( jxf_ParCSRMatrix    *par_A,
                JXF_Real              strength_threshold,
                JXF_Real              max_row_sum,
                JXF_Int                 num_functions,
                JXF_Int                *dof_func,
                jxf_ParCSRMatrix   **S_ptr )
{
   MPI_Comm              comm        = jxf_ParCSRMatrixComm(par_A);
   jxf_ParCSRCommPkg     *comm_pkg    = jxf_ParCSRMatrixCommPkg(par_A);
   jxf_ParCSRCommHandle  *comm_handle;
   jxf_CSRMatrix         *A_diag      = jxf_ParCSRMatrixDiag(par_A);
   JXF_Int                  *A_diag_i    = jxf_CSRMatrixI(A_diag);
   JXF_Real               *A_diag_data = jxf_CSRMatrixData(A_diag);


   jxf_CSRMatrix       *A_offd          = jxf_ParCSRMatrixOffd(par_A);
   JXF_Int                *A_offd_i        = jxf_CSRMatrixI(A_offd);
   JXF_Real             *A_offd_data     = NULL;
   JXF_Int                *A_diag_j        = jxf_CSRMatrixJ(A_diag);
   JXF_Int                *A_offd_j        = jxf_CSRMatrixJ(A_offd);

   JXF_Int 		      *row_starts      = jxf_ParCSRMatrixRowStarts(par_A);
   JXF_Int                 num_variables   = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int                 global_num_vars = jxf_ParCSRMatrixGlobalNumRows(par_A);
   JXF_Int 		       num_nonzeros_diag;
   JXF_Int 		       num_nonzeros_offd = 0;
   JXF_Int 		       num_cols_offd     = 0;

   jxf_ParCSRMatrix    *par_S;
   jxf_CSRMatrix       *S_diag;
   JXF_Int                *S_diag_i;
   JXF_Int                *S_diag_j;
   jxf_CSRMatrix       *S_offd;
   JXF_Int                *S_offd_i;
   JXF_Int                *S_offd_j = NULL;

   JXF_Real              diag, row_scale, row_sum;
   JXF_Int                 i, jA, jS;
                      
   JXF_Int                 ierr = 0;

   JXF_Int                *dof_func_offd;
   JXF_Int                 num_sends;
   JXF_Int                *int_buf_data;
   JXF_Int                 index, start, j;
   
   JXF_Int                 num_dd_strong = 0;
   JXF_Int                 num_not_dd = 0;
   JXF_Int                 num_not_dd_strong = 0;
   JXF_Int                 print_ai = 0;

  /*-------------------------------------------------------------------
   * Compute a  ParCSR strength matrix, par_S.
   *
   * For now, the "strength" of dependence/influence is defined in
   * the following way: i depends on j if
   *     aij > jxf_max(k != i) aik,    aii < 0
   * or
   *     aij < jxf_min(k != i) aik,    aii >= 0
   * Then S_ij = 1, else S_ij = 0.
   *
   * NOTE: the entries are negative initially, corresponding
   * to "unaccounted-for" dependence.
   *-----------------------------------------------------------------*/

   num_nonzeros_diag = A_diag_i[num_variables];
   num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   A_offd_i = jxf_CSRMatrixI(A_offd);
   num_nonzeros_offd = A_offd_i[num_variables];

   par_S = jxf_ParCSRMatrixCreate(comm, global_num_vars, global_num_vars,
                                 row_starts, row_starts, num_cols_offd,
                                 num_nonzeros_diag, num_nonzeros_offd);
   
   
   /* row_starts is owned by A, col_starts = row_starts */
   jxf_ParCSRMatrixSetRowStartsOwner(par_S, 0);
   S_diag = jxf_ParCSRMatrixDiag(par_S);
   jxf_CSRMatrixI(S_diag) = jxf_CTAlloc(JXF_Int, num_variables + 1);
   jxf_CSRMatrixJ(S_diag) = jxf_CTAlloc(JXF_Int, num_nonzeros_diag);
   S_offd = jxf_ParCSRMatrixOffd(par_S);
   jxf_CSRMatrixI(S_offd) = jxf_CTAlloc(JXF_Int, num_variables + 1);

   S_diag_i = jxf_CSRMatrixI(S_diag);
   S_diag_j = jxf_CSRMatrixJ(S_diag);
   S_offd_i = jxf_CSRMatrixI(S_offd);

   dof_func_offd = NULL;

   if (num_cols_offd)
   {
      A_offd_data = jxf_CSRMatrixData(A_offd);
      jxf_CSRMatrixJ(S_offd) = jxf_CTAlloc(JXF_Int, num_nonzeros_offd);
      S_offd_j = jxf_CSRMatrixJ(S_offd);
      jxf_ParCSRMatrixColMapOffd(par_S) = jxf_CTAlloc(JXF_Int, num_cols_offd);
      if (num_functions > 1)
      {
         dof_func_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      }
   }


  /*-------------------------------------------------------------------
   * Get the dof_func data for the off-processor columns
   *------------------------------------------------------------------*/

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   if (num_functions > 1)
   {
      int_buf_data = jxf_CTAlloc(JXF_Int,jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
	 start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	 for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
         {
            int_buf_data[index++] = dof_func[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
         }
      }
	
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, dof_func_offd);

      jxf_ParCSRCommHandleDestroy(comm_handle);   
      jxf_TFree(int_buf_data);
   }

   /* give par_S same nonzero structure as par_A */
   jxf_ParCSRMatrixCopy(par_A, par_S, 0);

#define JXF_SMP_PRIVATE i,diag,row_scale,row_sum,jA
#include "../../include/jxf_smp_forloop.h"
   for (i = 0; i < num_variables; i ++)
   {
      diag = A_diag_data[A_diag_i[i]];

      /* compute scaling factor and row sum */
      row_scale = 0.0;
      row_sum = diag;
      if (num_functions > 1)
      {
         if (diag < 0)
	 {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
            {
               if (dof_func[i] == dof_func[A_diag_j[jA]])
               {
                  row_scale = jxf_max(row_scale, A_diag_data[jA]);
                  row_sum += A_diag_data[jA];
               }
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
            {
               if (dof_func[i] == dof_func_offd[A_offd_j[jA]])
               {
                  row_scale = jxf_max(row_scale, A_offd_data[jA]);
                  row_sum += A_offd_data[jA];
               }
            }
         }
         else
         {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
            {
               if (dof_func[i] == dof_func[A_diag_j[jA]])
               {
                  row_scale = jxf_min(row_scale, A_diag_data[jA]);
                  row_sum += A_diag_data[jA];
               }
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
            {
               if (dof_func[i] == dof_func_offd[A_offd_j[jA]])
               {
                  row_scale = jxf_min(row_scale, A_offd_data[jA]);
                  row_sum += A_offd_data[jA];
               }
            }
         }
      }
      else
      {
         if (diag < 0)
	 {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
            {
               row_scale = jxf_max(row_scale, A_diag_data[jA]);
               row_sum += A_diag_data[jA];
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
            {
               row_scale = jxf_max(row_scale, A_offd_data[jA]);
               row_sum += A_offd_data[jA];
            }
         }
         else
         {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
            {
               row_scale = jxf_min(row_scale, A_diag_data[jA]);
               row_sum += A_diag_data[jA];
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
            {
               row_scale = jxf_min(row_scale, A_offd_data[jA]);
               row_sum += A_offd_data[jA];
            }
         }
      }

      if (diag > 0.0) {
         if (row_sum < 0.0) {
            //jxf_printf("i = %d, row_sum =%f, diag = %f\n ",i, row_sum, diag);
            num_not_dd++;
         } else if (row_sum > max_row_sum*diag) {
            //jxf_printf("i = %d, row_sum =%f, diag = %f\n ",i, row_sum, diag);
            num_dd_strong++;
         }
      } else if (diag < 0.0) {
         if (row_sum > 0.0) {
            //jxf_printf("i = %d, row_sum =%f, diag = %f\n ",i, row_sum, diag);
            num_not_dd++;
         } else if (fabs(row_sum) > max_row_sum*fabs(diag) ) {
            num_dd_strong++;
         }
      }
      if (row_sum / diag < -0.1 ) num_not_dd_strong++;

      row_sum = fabs( row_sum / diag );

      /* compute row entries of par_S */
      S_diag_j[A_diag_i[i]] = -1;
      if ( (row_sum > max_row_sum) && (max_row_sum < 1.0) )
      {
         /* make all dependencies weak */
         for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
         {
            S_diag_j[jA] = -1;
         }
         for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
         {
            S_offd_j[jA] = -1;
         }
      }
      else
      {
         if (num_functions > 1)
         { 
            if (diag < 0) 
            { 
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
               {
                  if (A_diag_data[jA] <= strength_threshold*row_scale || dof_func[i] != dof_func[A_diag_j[jA]])
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
               {
                  if (A_offd_data[jA] <= strength_threshold*row_scale || dof_func[i] != dof_func_offd[A_offd_j[jA]])
                  {
                     S_offd_j[jA] = -1;
                  }
               }
            }
            else
            {
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
               {
                  if (A_diag_data[jA] >= strength_threshold*row_scale || dof_func[i] != dof_func[A_diag_j[jA]])
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
               {
                  if (A_offd_data[jA] >= strength_threshold*row_scale || dof_func[i] != dof_func_offd[A_offd_j[jA]])
                  {
                     S_offd_j[jA] = -1;
                  }
               }
            }
         }
         else
         {
            if (diag < 0) 
            { 
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
               {
                  if (A_diag_data[jA] <= strength_threshold*row_scale)
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
               {
                  if (A_offd_data[jA] <= strength_threshold*row_scale)
                  {
                     S_offd_j[jA] = -1;
                  }
               }
            }
            else
            {
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA ++)
               {
                  if (A_diag_data[jA] >= strength_threshold*row_scale)
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
               {
                  if (A_offd_data[jA] >= strength_threshold*row_scale)
                  {
                     S_offd_j[jA] = -1;
                  }
               }
            }
         }
      }
   }

   if (print_ai == 1) {
   JXF_Real num_strong = num_dd_strong;
   JXF_Real num_not = num_not_dd;
   JXF_Real num_not_strong = num_not_dd_strong;
   JXF_Real num_var = num_variables;
   jxf_printf(" - num_dd_strong = %d, ratio_dd_strong = %f \n", num_dd_strong, num_strong/num_var);
   jxf_printf(" - num_not_dd = %d, ratio_not_dd = %f \n", num_not_dd, num_not/num_var);
   jxf_printf(" - num_not_dd_strong = %d, ratio_not_dd_strong = %f \n", num_not_dd_strong, num_not_strong/num_var);
   }

  /*--------------------------------------------------------------
   * "Compress" the strength matrix.
   *
   * NOTE: S has *NO DIAGONAL ELEMENT* on any row.  Caveat Emptor!
   *
   * NOTE: This "compression" section of code may be removed, and
   * coarsening will still be done correctly.  However, the routine
   * that builds interpolation would have to be modified first.
   *----------------------------------------------------------------*/

   /* RDF: not sure if able to thread this loop */
   jS = 0;
   for (i = 0; i < num_variables; i ++)
   {
      S_diag_i[i] = jS;
      for (jA = A_diag_i[i]; jA < A_diag_i[i+1]; jA ++)
      {
         if (S_diag_j[jA] > -1)
         {
            S_diag_j[jS] = S_diag_j[jA];
            jS ++;
         }
      }
   }
   S_diag_i[num_variables] = jS;
   jxf_CSRMatrixNumNonzeros(S_diag) = jS;

   /* RDF: not sure if able to thread this loop */
   jS = 0;
   for (i = 0; i < num_variables; i ++)
   {
      S_offd_i[i] = jS;
      for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA ++)
      {
         if (S_offd_j[jA] > -1)
         {
            S_offd_j[jS] = S_offd_j[jA];
            jS ++;
         }
      }
   }
   S_offd_i[num_variables] = jS;
   jxf_CSRMatrixNumNonzeros(S_offd) = jS;
   jxf_ParCSRMatrixCommPkg(par_S) = NULL;

   *S_ptr = par_S;

   jxf_TFree(dof_func_offd);

   return (ierr);
}



/*!
 * \fn JXF_Int jxf_PAMGCreateSabs
 * \brief Generates strength matrix
 * \date 2018/05/10
 */
JXF_Int
jxf_PAMGCreateSabs( jxf_ParCSRMatrix  *A,
                   JXF_Real           strength_threshold,
                   JXF_Real           max_row_sum,
                   JXF_Int            num_functions,
                   JXF_Int           *dof_func,
                   jxf_ParCSRMatrix **S_ptr)
{
   MPI_Comm 	        comm     = jxf_ParCSRMatrixComm(A);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_ParCSRCommHandle *comm_handle;
   jxf_CSRMatrix    *A_diag       = jxf_ParCSRMatrixDiag(A);
   JXF_Int          *A_diag_i     = jxf_CSRMatrixI(A_diag);
   JXF_Real         *A_diag_data  = jxf_CSRMatrixData(A_diag);


   jxf_CSRMatrix    *A_offd       = jxf_ParCSRMatrixOffd(A);
   JXF_Int          *A_offd_i     = jxf_CSRMatrixI(A_offd);
   JXF_Real         *A_offd_data  = NULL;
   JXF_Int          *A_diag_j     = jxf_CSRMatrixJ(A_diag);
   JXF_Int          *A_offd_j     = jxf_CSRMatrixJ(A_offd);

   JXF_Int 		 *row_starts      = jxf_ParCSRMatrixRowStarts(A);
   JXF_Int                 num_variables   = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int                 global_num_vars = jxf_ParCSRMatrixGlobalNumRows(A);
   JXF_Int 		  num_nonzeros_diag;
   JXF_Int 		  num_nonzeros_offd = 0;
   JXF_Int 		  num_cols_offd = 0;

   jxf_ParCSRMatrix *S;
   jxf_CSRMatrix    *S_diag;
   JXF_Int          *S_diag_i;
   JXF_Int          *S_diag_j;
   jxf_CSRMatrix    *S_offd;
   JXF_Int          *S_offd_i = NULL;
   JXF_Int          *S_offd_j = NULL;

   JXF_Real          diag, row_scale, row_sum;
   JXF_Int           i, jA, jS;
                      
   JXF_Int           ierr = 0;

   JXF_Int          *dof_func_offd;
   JXF_Int	    num_sends;
   JXF_Int	   *int_buf_data;
   JXF_Int	    index, start, j;
   
   /*--------------------------------------------------------------
    * Compute a  ParCSR strength matrix, S.
    *
    * For now, the "strength" of dependence/influence is defined in
    * the following way: i depends on j if
    *     aij > hypre_max (k != i) aik,    aii < 0
    * or
    *     aij < hypre_min (k != i) aik,    aii >= 0
    * Then S_ij = 1, else S_ij = 0.
    *
    * NOTE: the entries are negative initially, corresponding
    * to "unaccounted-for" dependence.
    *----------------------------------------------------------------*/

   num_nonzeros_diag = A_diag_i[num_variables];
   num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   A_offd_i = jxf_CSRMatrixI(A_offd);
   num_nonzeros_offd = A_offd_i[num_variables];

   S = jxf_ParCSRMatrixCreate(comm, global_num_vars, global_num_vars,
			row_starts, row_starts,
			num_cols_offd, num_nonzeros_diag, num_nonzeros_offd);


   /* row_starts is owned by A, col_starts = row_starts */
   jxf_ParCSRMatrixSetRowStartsOwner(S,0);
   S_diag = jxf_ParCSRMatrixDiag(S);
   jxf_CSRMatrixI(S_diag) = jxf_CTAlloc(JXF_Int,  num_variables+1);
   jxf_CSRMatrixJ(S_diag) = jxf_CTAlloc(JXF_Int,  num_nonzeros_diag);
   S_offd = jxf_ParCSRMatrixOffd(S);
   jxf_CSRMatrixI(S_offd) = jxf_CTAlloc(JXF_Int,  num_variables+1);

   S_diag_i = jxf_CSRMatrixI(S_diag);
   S_diag_j = jxf_CSRMatrixJ(S_diag);
   S_offd_i = jxf_CSRMatrixI(S_offd);

   dof_func_offd = NULL;

   if (num_cols_offd)
   {
        A_offd_data = jxf_CSRMatrixData(A_offd);
        jxf_CSRMatrixJ(S_offd) = jxf_CTAlloc(JXF_Int,  num_nonzeros_offd);
        S_offd_j = jxf_CSRMatrixJ(S_offd);
        jxf_ParCSRMatrixColMapOffd(S) = jxf_CTAlloc(JXF_Int,  num_cols_offd);
        if (num_functions > 1)
	{
	   dof_func_offd = jxf_CTAlloc(JXF_Int,  num_cols_offd);
	}
   }


  /*-------------------------------------------------------------------
    * Get the dof_func data for the off-processor columns
    *-------------------------------------------------------------------*/

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A); 
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   if (num_functions > 1)
   {
      int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i++)
      {
	 start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	 for (j=start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j++)
	 {
	    int_buf_data[index++] = dof_func[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
	 }
      }
	
      comm_handle = jxf_ParCSRCommHandleCreate( 11, comm_pkg, int_buf_data, dof_func_offd);

      jxf_ParCSRCommHandleDestroy(comm_handle);   
      jxf_TFree(int_buf_data);
   }

   /* give S same nonzero structure as A */
   jxf_ParCSRMatrixCopy(A,S,0);

#define JXF_SMP_PRIVATE i,diag,row_scale,row_sum,jA
#include "../../include/jxf_smp_forloop.h"
   for (i = 0; i < num_variables; i++)
   {
      diag = A_diag_data[A_diag_i[i]];

      /* compute scaling factor and row sum */
      row_scale = 0.0;
      row_sum = diag;
      if (num_functions > 1)
      {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA++)
            {
               if (dof_func[i] == dof_func[A_diag_j[jA]])
               {
                  row_scale = jxf_max(row_scale, fabs(A_diag_data[jA]));
                  row_sum += fabs(A_diag_data[jA]);
               }
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
            {
               if (dof_func[i] == dof_func_offd[A_offd_j[jA]])
               {
                  row_scale = jxf_max(row_scale, fabs(A_offd_data[jA]));
                  row_sum += fabs(A_offd_data[jA]);
               }
            }
      }
      else
      {
            for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA++)
            {
               row_scale = jxf_max(row_scale, fabs(A_diag_data[jA]));
               row_sum += fabs(A_diag_data[jA]);
            }
            for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
            {
               row_scale = jxf_max(row_scale, fabs(A_offd_data[jA]));
               row_sum += fabs(A_offd_data[jA]);
            }
      }

      /* compute row entries of S */
      S_diag_j[A_diag_i[i]] = -1;
      if ((fabs(row_sum) > fabs(diag)*max_row_sum) && (max_row_sum < 1.0))
      {
         /* make all dependencies weak */
         for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA++)
         {
            S_diag_j[jA] = -1;
         }
         for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
         {
            S_offd_j[jA] = -1;
         }
      }
      else
      {
         if (num_functions > 1)
         { 
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA++)
               {
                  if (fabs(A_diag_data[jA]) <= strength_threshold * row_scale || dof_func[i] != dof_func[A_diag_j[jA]])
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
               {
                  if (fabs(A_offd_data[jA]) <= strength_threshold * row_scale || dof_func[i] != dof_func_offd[A_offd_j[jA]])
                  {
                     S_offd_j[jA] = -1;
                  }
               }
         }
         else
         {
               for (jA = A_diag_i[i]+1; jA < A_diag_i[i+1]; jA++)
               {
                  if (fabs(A_diag_data[jA]) <= strength_threshold * row_scale)
                  {
                     S_diag_j[jA] = -1;
                  }
               }
               for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
               {
                  if (fabs(A_offd_data[jA]) <= strength_threshold * row_scale)
                  {
                     S_offd_j[jA] = -1;
                  }
               }
         }
      }
   }

   /*--------------------------------------------------------------
    * "Compress" the strength matrix.
    *
    * NOTE: S has *NO DIAGONAL ELEMENT* on any row.  Caveat Emptor!
    *
    * NOTE: This "compression" section of code may be removed, and
    * coarsening will still be done correctly.  However, the routine
    * that builds interpolation would have to be modified first.
    *----------------------------------------------------------------*/

   /* RDF: not sure if able to thread this loop */
   jS = 0;
   for (i = 0; i < num_variables; i++)
   {
      S_diag_i[i] = jS;
      for (jA = A_diag_i[i]; jA < A_diag_i[i+1]; jA++)
      {
         if (S_diag_j[jA] > -1)
         {
            S_diag_j[jS] = S_diag_j[jA];
            jS ++;
         }
      }
   }
   S_diag_i[num_variables] = jS;
   jxf_CSRMatrixNumNonzeros(S_diag) = jS;

   /* RDF: not sure if able to thread this loop */
   jS = 0;
   for (i = 0; i < num_variables; i++)
   {
      S_offd_i[i] = jS;
      for (jA = A_offd_i[i]; jA < A_offd_i[i+1]; jA++)
      {
         if (S_offd_j[jA] > -1)
         {
            S_offd_j[jS] = S_offd_j[jA];
            jS ++;
         }
      }
   }
   S_offd_i[num_variables] = jS;
   jxf_CSRMatrixNumNonzeros(S_offd) = jS;
   jxf_ParCSRMatrixCommPkg(S) = NULL;

   *S_ptr = S;

   jxf_TFree(dof_func_offd);

   return (ierr);
}



/*!
 * \fn JXF_Int jxf_PAMGCreateSCommPkg
 * \brief Create a communication package for S.
 * \date 2011/09/03
 */  
JXF_Int
jxf_PAMGCreateSCommPkg( jxf_ParCSRMatrix   *par_A, 
                       jxf_ParCSRMatrix   *par_S,
                       JXF_Int              **col_offd_S_to_A_ptr )
{
   MPI_Comm 	         comm = jxf_ParCSRMatrixComm(par_A);
   MPI_Status	        *status;
   MPI_Request	        *requests;
   jxf_ParCSRCommPkg     *comm_pkg_A = jxf_ParCSRMatrixCommPkg(par_A);
   jxf_ParCSRCommPkg     *comm_pkg_S;
   jxf_ParCSRCommHandle  *comm_handle;
   jxf_CSRMatrix         *A_offd = jxf_ParCSRMatrixOffd(par_A);
   JXF_Int                  *col_map_offd_A = jxf_ParCSRMatrixColMapOffd(par_A);

   jxf_CSRMatrix  *S_diag = jxf_ParCSRMatrixDiag(par_S);
   jxf_CSRMatrix  *S_offd = jxf_ParCSRMatrixOffd(par_S);
   JXF_Int           *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int           *S_offd_j = jxf_CSRMatrixJ(S_offd);
   JXF_Int           *col_map_offd_S = jxf_ParCSRMatrixColMapOffd(par_S);

   JXF_Int  *recv_procs_A = jxf_ParCSRCommPkgRecvProcs(comm_pkg_A);
   JXF_Int  *recv_vec_starts_A = jxf_ParCSRCommPkgRecvVecStarts(comm_pkg_A);
   JXF_Int  *send_procs_A = jxf_ParCSRCommPkgSendProcs(comm_pkg_A);
   JXF_Int  *send_map_starts_A = jxf_ParCSRCommPkgSendMapStarts(comm_pkg_A);
   JXF_Int  *recv_procs_S;
   JXF_Int  *recv_vec_starts_S;
   JXF_Int  *send_procs_S;
   JXF_Int  *send_map_starts_S;
   JXF_Int  *send_map_elmts_S;
   JXF_Int  *col_offd_S_to_A;

   JXF_Int  *S_marker;
   JXF_Int  *send_change;
   JXF_Int  *recv_change;

   JXF_Int num_variables   = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int num_cols_offd_A = jxf_CSRMatrixNumCols(A_offd);                 
   JXF_Int num_cols_offd_S;
   JXF_Int i, j, jcol;
   JXF_Int proc, cnt, proc_cnt, total_nz;
   JXF_Int first_row;
                      
   JXF_Int ierr = 0;

   JXF_Int num_sends_A = jxf_ParCSRCommPkgNumSends(comm_pkg_A);
   JXF_Int num_recvs_A = jxf_ParCSRCommPkgNumRecvs(comm_pkg_A);
   JXF_Int num_sends_S;
   JXF_Int num_recvs_S;
   JXF_Int num_nonzeros;

   num_nonzeros = S_offd_i[num_variables];

   S_marker = NULL;
   if (num_cols_offd_A)
   {
      S_marker = jxf_CTAlloc(JXF_Int,num_cols_offd_A);
   }

   for (i = 0; i < num_cols_offd_A; i ++)
   {
      S_marker[i] = -1;
   }

   for (i = 0; i < num_nonzeros; i ++)
   {
      jcol = S_offd_j[i];
      S_marker[jcol] = 0;
   }

   proc = 0;
   proc_cnt = 0;
   cnt = 0;
   num_recvs_S = 0;
   for (i = 0; i < num_recvs_A; i ++)
   {
      for (j = recv_vec_starts_A[i]; j < recv_vec_starts_A[i+1]; j ++)
      {
         if (!S_marker[j])
         {
            S_marker[j] = cnt;
            cnt ++;
            proc = 1;
         }
      }
      if (proc) 
      {
         num_recvs_S ++; 
         proc = 0;
      }
   }

   num_cols_offd_S = cnt;  
   recv_change = NULL;
   recv_procs_S = NULL;
   send_change = NULL;
   if (col_map_offd_S) 
   {
      jxf_TFree(col_map_offd_S);
   }
   col_map_offd_S = NULL;
   col_offd_S_to_A = NULL;
   
   if (num_recvs_A) 
   {
      recv_change = jxf_CTAlloc(JXF_Int, num_recvs_A);
   }
   
   if (num_sends_A) 
   {
      send_change = jxf_CTAlloc(JXF_Int, num_sends_A);
   }
   
   if (num_recvs_S) 
   {
      recv_procs_S = jxf_CTAlloc(JXF_Int, num_recvs_S);
   }
   
   recv_vec_starts_S = jxf_CTAlloc(JXF_Int, num_recvs_S + 1);
   
   if (num_cols_offd_S)
   {
      col_map_offd_S = jxf_CTAlloc(JXF_Int,num_cols_offd_S);
      col_offd_S_to_A = jxf_CTAlloc(JXF_Int,num_cols_offd_S);
   }
   
   if (num_cols_offd_S < num_cols_offd_A)
   {
      for (i = 0; i < num_nonzeros; i ++)
      {
         jcol = S_offd_j[i];
         S_offd_j[i] = S_marker[jcol];
      }

      proc = 0;
      proc_cnt = 0;
      cnt = 0;
      recv_vec_starts_S[0] = 0;
      for (i = 0; i < num_recvs_A; i ++)
      {
         for (j = recv_vec_starts_A[i]; j < recv_vec_starts_A[i+1]; j ++)
         {
            if (S_marker[j] != -1)
            {
               col_map_offd_S[cnt] = col_map_offd_A[j];
               col_offd_S_to_A[cnt++] = j;
               proc = 1;
            }
         }
         recv_change[i] = j - cnt - recv_vec_starts_A[i] + recv_vec_starts_S[proc_cnt];
         if (proc)
         {
            recv_procs_S[proc_cnt++] = recv_procs_A[i];
            recv_vec_starts_S[proc_cnt] = cnt;
            proc = 0;
         }
      }
   }
   else
   {
      for (i = 0; i < num_recvs_A; i ++)
      {
         for (j = recv_vec_starts_A[i]; j < recv_vec_starts_A[i+1]; j ++)
         {
            col_map_offd_S[j]  = col_map_offd_A[j];
            col_offd_S_to_A[j] = j;
         }
         recv_procs_S[i] = recv_procs_A[i];
         recv_vec_starts_S[i] = recv_vec_starts_A[i];
      }
      recv_vec_starts_S[num_recvs_A] = recv_vec_starts_A[num_recvs_A];
   } 

   requests = jxf_CTAlloc(MPI_Request, num_sends_A + num_recvs_A);
   
   j = 0;
   for (i = 0; i < num_sends_A; i ++)
   {
      jxf_MPI_Irecv(&send_change[i], 1, JXF_MPI_INT, send_procs_A[i], 0, comm, &requests[j++]);
   }

   for (i = 0; i < num_recvs_A; i ++)
   {
      jxf_MPI_Isend(&recv_change[i], 1, JXF_MPI_INT, recv_procs_A[i], 0, comm, &requests[j++]);
   }

   status = jxf_CTAlloc(MPI_Status, j);
   jxf_MPI_Waitall(j, requests, status);
   jxf_TFree(status);
   jxf_TFree(requests);

   num_sends_S = 0;
   total_nz = send_map_starts_A[num_sends_A];
   for (i = 0; i < num_sends_A; i ++)
   {
      if (send_change[i])
      {
	 if ((send_map_starts_A[i+1] - send_map_starts_A[i]) > send_change[i])
	 {
	    num_sends_S ++;
	 }
      }
      else
      {
	 num_sends_S ++;
      }
      total_nz -= send_change[i];
   }

   send_procs_S = NULL;
   if (num_sends_S)
   {
      send_procs_S = jxf_CTAlloc(JXF_Int, num_sends_S);
   }
   send_map_starts_S = jxf_CTAlloc(JXF_Int, num_sends_S + 1);
   send_map_elmts_S = NULL;
   
   if (total_nz)
   {
      send_map_elmts_S = jxf_CTAlloc(JXF_Int, total_nz);
   }

   proc = 0;
   proc_cnt = 0;
   for (i = 0; i < num_sends_A; i ++)
   {
      cnt = send_map_starts_A[i+1] - send_map_starts_A[i] - send_change[i];
      if (cnt)
      {
	 send_procs_S[proc_cnt++] = send_procs_A[i];
         send_map_starts_S[proc_cnt] = send_map_starts_S[proc_cnt-1] + cnt;
      }
   }

   comm_pkg_S = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
   jxf_ParCSRCommPkgComm(comm_pkg_S) = comm;

   jxf_ParCSRCommPkgNumRecvs(comm_pkg_S) = num_recvs_S;
   jxf_ParCSRCommPkgRecvProcs(comm_pkg_S) = recv_procs_S;
   jxf_ParCSRCommPkgRecvVecStarts(comm_pkg_S) = recv_vec_starts_S;
   jxf_ParCSRCommPkgNumSends(comm_pkg_S) = num_sends_S;
   jxf_ParCSRCommPkgSendProcs(comm_pkg_S) = send_procs_S;
   jxf_ParCSRCommPkgSendMapStarts(comm_pkg_S) = send_map_starts_S;
#if JXF_REODER_SEND_RECV
jxf_hpIfReorderMatvecCommPkg(comm_pkg_S) = 0;
  /*-----------------------------------------------------------
   *  reoderMatvecCommPkg
   *----------------------------------------------------------*/
jxf_hpCSRCommPkgReorder(comm_pkg_S);
#endif
   comm_handle = jxf_ParCSRCommHandleCreate(12, comm_pkg_S, col_map_offd_S, send_map_elmts_S);
   jxf_ParCSRCommHandleDestroy(comm_handle);

   first_row = jxf_ParCSRMatrixFirstRowIndex(par_A);
   if (first_row)
   {
      for (i = 0; i < send_map_starts_S[num_sends_S]; i ++)
      {
         send_map_elmts_S[i] -= first_row;
      }
   }

   jxf_ParCSRCommPkgSendMapElmts(comm_pkg_S) = send_map_elmts_S;
  
   jxf_ParCSRMatrixCommPkg(par_S) = comm_pkg_S;
   jxf_ParCSRMatrixColMapOffd(par_S) = col_map_offd_S;
   jxf_CSRMatrixNumCols(S_offd) = num_cols_offd_S;

   jxf_TFree(S_marker);
   jxf_TFree(send_change);
   jxf_TFree(recv_change);

   *col_offd_S_to_A_ptr = col_offd_S_to_A;

   return ierr;
}


/*!
 * \fn JXF_Int jxf_PAMGCreate2ndS
 * \brief Creates strength matrix on coarse points for
          second coarsening pass in aggressive coarsening (S*S+2S)
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGCreate2ndS( jxf_ParCSRMatrix  *S,
                   JXF_Int              *CF_marker,
                   JXF_Int               num_paths,
                   JXF_Int              *coarse_row_starts,
                   jxf_ParCSRMatrix **C_ptr )
{
   MPI_Comm 	   comm = jxf_ParCSRMatrixComm(S);
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(S);
   jxf_ParCSRCommPkg *tmp_comm_pkg;
   jxf_ParCSRCommHandle *comm_handle;

   jxf_CSRMatrix *S_diag = jxf_ParCSRMatrixDiag(S);

   JXF_Int *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int *S_diag_j = jxf_CSRMatrixJ(S_diag);

   jxf_CSRMatrix *S_offd = jxf_ParCSRMatrixOffd(S);

   JXF_Int *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int *S_offd_j = jxf_CSRMatrixJ(S_offd);

   JXF_Int num_cols_diag_S = jxf_CSRMatrixNumCols(S_diag);
   JXF_Int num_cols_offd_S = jxf_CSRMatrixNumCols(S_offd);

   jxf_ParCSRMatrix *S2;
   JXF_Int *col_map_offd_C = NULL;

   jxf_CSRMatrix *C_diag;

   JXF_Int *C_diag_data = NULL;
   JXF_Int *C_diag_i;
   JXF_Int *C_diag_j = NULL;

   jxf_CSRMatrix *C_offd;

   JXF_Int *C_offd_data = NULL;
   JXF_Int *C_offd_i;
   JXF_Int *C_offd_j = NULL;

   JXF_Int C_diag_size;
   JXF_Int C_offd_size;
   JXF_Int num_cols_offd_C = 0;

   JXF_Int *S_ext_diag_i = NULL;
   JXF_Int *S_ext_diag_j = NULL;
   JXF_Int S_ext_diag_size = 0;

   JXF_Int *S_ext_offd_i = NULL;
   JXF_Int *S_ext_offd_j = NULL;
   JXF_Int S_ext_offd_size = 0;

   JXF_Int *CF_marker_offd = NULL;

   JXF_Int *S_marker = NULL;
   JXF_Int *S_marker_offd = NULL;
   JXF_Int *temp = NULL;

   JXF_Int *fine_to_coarse = NULL;
   JXF_Int *fine_to_coarse_offd = NULL;
   JXF_Int *map_S_to_C = NULL;

   JXF_Int num_sends = 0;
   JXF_Int num_recvs = 0;
   JXF_Int *send_map_starts;
   JXF_Int *tmp_send_map_starts = NULL;
   JXF_Int *send_map_elmts;
   JXF_Int *recv_vec_starts;
   JXF_Int *tmp_recv_vec_starts = NULL;
   JXF_Int *int_buf_data = NULL;

   JXF_Int i, j, k;
   JXF_Int i1, i2, i3;
   JXF_Int jj1, jj2, jcol, jrow, j_cnt;

   JXF_Int jj_count_diag, jj_count_offd;
   JXF_Int jj_row_begin_diag, jj_row_begin_offd;
   JXF_Int cnt, cnt_offd, cnt_diag;
   JXF_Int num_procs, my_id;
   JXF_Int value, index;
   JXF_Int num_coarse;
   JXF_Int num_coarse_offd;
   JXF_Int num_nonzeros;
   JXF_Int num_nonzeros_diag;
   JXF_Int num_nonzeros_offd;
   JXF_Int global_num_coarse;
   JXF_Int my_first_cpt, my_last_cpt;

   JXF_Int *S_int_i = NULL;
   JXF_Int *S_int_j = NULL;
   JXF_Int *S_ext_i = NULL;
   JXF_Int *S_ext_j = NULL;

   /*-----------------------------------------------------------------------
    *  Extract S_ext, i.e. portion of B that is stored on neighbor procs
    *  and needed locally for matrix matrix product 
    *-----------------------------------------------------------------------*/

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION
   my_first_cpt = coarse_row_starts[0];
   my_last_cpt = coarse_row_starts[1]-1;
   if (my_id == (num_procs -1)) global_num_coarse = coarse_row_starts[1];
   jxf_MPI_Bcast(&global_num_coarse, 1, JXF_MPI_INT, num_procs-1, comm);
#else
   my_first_cpt = coarse_row_starts[my_id];
   my_last_cpt = coarse_row_starts[my_id+1]-1;
   global_num_coarse = coarse_row_starts[num_procs];
#endif

   if (num_cols_offd_S)
   {
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd_S);
      fine_to_coarse_offd = jxf_CTAlloc(JXF_Int, num_cols_offd_S);
   }

   if (num_cols_diag_S) fine_to_coarse = jxf_CTAlloc(JXF_Int, num_cols_diag_S);

   num_coarse = 0;
   for (i = 0; i < num_cols_diag_S; i ++)
   {
      if (CF_marker[i] > 0) 
      {
         fine_to_coarse[i] = num_coarse + my_first_cpt;
         num_coarse ++;
      }
      else
      {
         fine_to_coarse[i] = -1;
      }
   }

   if (num_procs > 1)
   {
      if (!comm_pkg)
      {
         jxf_MatvecCommPkgCreate(S);
         comm_pkg = jxf_ParCSRMatrixCommPkg(S);
      }
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      send_map_starts = jxf_ParCSRCommPkgSendMapStarts(comm_pkg);
      send_map_elmts = jxf_ParCSRCommPkgSendMapElmts(comm_pkg);
      num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
      recv_vec_starts = jxf_ParCSRCommPkgRecvVecStarts(comm_pkg);
      int_buf_data = jxf_CTAlloc(JXF_Int, send_map_starts[num_sends]);

      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
            int_buf_data[index++] = fine_to_coarse[send_map_elmts[j]];
      }

      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, fine_to_coarse_offd);

      for (i = 0; i < num_cols_diag_S; i ++)
         if (CF_marker[i] > 0) fine_to_coarse[i] -= my_first_cpt;

      jxf_ParCSRCommHandleDestroy(comm_handle);

      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
         {
            int_buf_data[index++] = CF_marker[send_map_elmts[j]];
         }
      }

      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);

      jxf_ParCSRCommHandleDestroy(comm_handle);
      jxf_TFree(int_buf_data);

      S_int_i = jxf_CTAlloc(JXF_Int, send_map_starts[num_sends]+1);
      S_ext_i = jxf_CTAlloc(JXF_Int, recv_vec_starts[num_recvs]+1);

     /*--------------------------------------------------------------------------
      * generate S_int_i through adding number of coarse row-elements of offd and diag
      * for corresponding rows. S_int_i[j+1] contains the number of coarse elements of
      * a row j (which is determined through send_map_elmts)
      *--------------------------------------------------------------------------*/
      S_int_i[0] = 0;
      j_cnt = 0;
      num_nonzeros = 0;
      for (i = 0; i < num_sends; i ++)
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
         {
            jrow = send_map_elmts[j];
            index = 0;
            for (k = S_diag_i[jrow]; k < S_diag_i[jrow+1]; k ++)
            {
	       if (CF_marker[S_diag_j[k]] > 0) index ++;
	    }
            for (k = S_offd_i[jrow]; k < S_offd_i[jrow+1]; k ++)
            {
	       if (CF_marker_offd[S_offd_j[k]] > 0) index ++;
	    }
            S_int_i[++j_cnt] = index;
            num_nonzeros += S_int_i[j_cnt];
         }
      }

     /*--------------------------------------------------------------------------
      * initialize communication
      *--------------------------------------------------------------------------*/
      if (num_procs > 1)
         comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, &S_int_i[1], &S_ext_i[1]);

      if (num_nonzeros) S_int_j = jxf_CTAlloc(JXF_Int, num_nonzeros);

      tmp_send_map_starts = jxf_CTAlloc(JXF_Int, num_sends+1);
      tmp_recv_vec_starts = jxf_CTAlloc(JXF_Int, num_recvs+1);

      tmp_send_map_starts[0] = 0;
      j_cnt = 0;
      for (i = 0; i < num_sends; i ++)
      {
         for (j = send_map_starts[i]; j < send_map_starts[i+1]; j ++)
         {
            jrow = send_map_elmts[j];
            for (k = S_diag_i[jrow]; k < S_diag_i[jrow+1]; k ++)
            {
               if (CF_marker[S_diag_j[k]] > 0)
		  S_int_j[j_cnt++] = fine_to_coarse[S_diag_j[k]]+my_first_cpt;
            }
            for (k = S_offd_i[jrow]; k < S_offd_i[jrow+1]; k ++)
            {
               if (CF_marker_offd[S_offd_j[k]] > 0)
                  S_int_j[j_cnt++] = fine_to_coarse_offd[S_offd_j[k]];
            }
         }
         tmp_send_map_starts[i+1] = j_cnt;
      }

      tmp_comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
      jxf_ParCSRCommPkgComm(tmp_comm_pkg) = comm;
      jxf_ParCSRCommPkgNumSends(tmp_comm_pkg) = num_sends;
      jxf_ParCSRCommPkgNumRecvs(tmp_comm_pkg) = num_recvs;
      jxf_ParCSRCommPkgSendProcs(tmp_comm_pkg) = jxf_ParCSRCommPkgSendProcs(comm_pkg);
      jxf_ParCSRCommPkgRecvProcs(tmp_comm_pkg) = jxf_ParCSRCommPkgRecvProcs(comm_pkg);
      jxf_ParCSRCommPkgSendMapStarts(tmp_comm_pkg) = tmp_send_map_starts;

      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;
     /*--------------------------------------------------------------------------
      * after communication exchange S_ext_i[j+1] contains the number of coarse elements
      * of a row j ! evaluate S_ext_i and compute num_nonzeros for S_ext
      *--------------------------------------------------------------------------*/
      for (i = 0; i < recv_vec_starts[num_recvs]; i ++)
         S_ext_i[i+1] += S_ext_i[i];

      num_nonzeros = S_ext_i[recv_vec_starts[num_recvs]];

      if (num_nonzeros) S_ext_j = jxf_CTAlloc(JXF_Int, num_nonzeros);

      tmp_recv_vec_starts[0] = 0;
      for (i = 0; i < num_recvs; i ++)
         tmp_recv_vec_starts[i+1] = S_ext_i[recv_vec_starts[i+1]];

      jxf_ParCSRCommPkgRecvVecStarts(tmp_comm_pkg) = tmp_recv_vec_starts;

   #if JXF_REODER_SEND_RECV
   jxf_hpIfReorderMatvecCommPkg(tmp_comm_pkg)   = 0;
   jxf_hpCSRCommPkgReorder(tmp_comm_pkg);
   #endif

      comm_handle = jxf_ParCSRCommHandleCreate(11, tmp_comm_pkg, S_int_j, S_ext_j);
      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;

      jxf_TFree(tmp_send_map_starts);
      jxf_TFree(tmp_recv_vec_starts);

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

      jxf_TFree(S_int_i);
      jxf_TFree(S_int_j);

      S_ext_diag_size = 0;
      S_ext_offd_size = 0;

      for (i = 0; i < num_cols_offd_S; i ++)
      {
         for (j = S_ext_i[i]; j < S_ext_i[i+1]; j ++)
         {
            if (S_ext_j[j] < my_first_cpt || S_ext_j[j] > my_last_cpt)
               S_ext_offd_size ++;
            else
               S_ext_diag_size ++;
         }
      }
      S_ext_diag_i = jxf_CTAlloc(JXF_Int, num_cols_offd_S+1);
      S_ext_offd_i = jxf_CTAlloc(JXF_Int, num_cols_offd_S+1);

      if (S_ext_diag_size)
      {
         S_ext_diag_j = jxf_CTAlloc(JXF_Int, S_ext_diag_size);
      }
      if (S_ext_offd_size)
      {
         S_ext_offd_j = jxf_CTAlloc(JXF_Int, S_ext_offd_size);
      }

      cnt_offd = 0;
      cnt_diag = 0;
      cnt = 0;
      num_coarse_offd = 0;
      for (i = 0; i < num_cols_offd_S; i ++)
      {
         if (CF_marker_offd[i] > 0) num_coarse_offd ++;

         for (j = S_ext_i[i]; j < S_ext_i[i+1]; j ++)
         {
            i1 = S_ext_j[j];
            if (i1 < my_first_cpt || i1 > my_last_cpt)
               S_ext_offd_j[cnt_offd++] = i1;
            else
               S_ext_diag_j[cnt_diag++] = i1 - my_first_cpt;
         }
         S_ext_diag_i[++cnt] = cnt_diag;
         S_ext_offd_i[cnt] = cnt_offd;
      }

      jxf_TFree(S_ext_i);
      jxf_TFree(S_ext_j);

      cnt = 0;
      if (S_ext_offd_size || num_coarse_offd)
      {
         temp = jxf_CTAlloc(JXF_Int, S_ext_offd_size+num_coarse_offd);
         for (i = 0; i < S_ext_offd_size; i ++)
            temp[i] = S_ext_offd_j[i];
         cnt = S_ext_offd_size;
         for (i = 0; i < num_cols_offd_S; i ++)
            if (CF_marker_offd[i] > 0) temp[cnt++] = fine_to_coarse_offd[i];
      }
      if (cnt)
      {
         jxf_qsort0(temp, 0, cnt-1);

         num_cols_offd_C = 1;
         value = temp[0];
         for (i = 1; i < cnt; i ++)
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

      for (i = 0; i < num_cols_offd_C; i ++)
         col_map_offd_C[i] = temp[i];

      if (S_ext_offd_size || num_coarse_offd)
         jxf_TFree(temp);

      for (i = 0 ; i < S_ext_offd_size; i ++)
         S_ext_offd_j[i] = jxf_BinarySearch(col_map_offd_C, S_ext_offd_j[i], num_cols_offd_C);
      if (num_cols_offd_S)
      {
         map_S_to_C = jxf_CTAlloc(JXF_Int, num_cols_offd_S);

         cnt = 0;
         for (i = 0; i < num_cols_offd_S; i ++)
         {
            if (CF_marker_offd[i] > 0)
            {
               while (fine_to_coarse_offd[i] > col_map_offd_C[cnt])
               {
                  cnt ++;
               }
               map_S_to_C[i] = cnt ++;
            }
            else
            {
               map_S_to_C[i] = -1;
            }
         }
      }
   }

   /*-----------------------------------------------------------------------
    *  Allocate and initialize some stuff.
    *-----------------------------------------------------------------------*/

   if (num_coarse) S_marker = jxf_CTAlloc(JXF_Int, num_coarse);

   for (i1 = 0; i1 < num_coarse; i1 ++)
      S_marker[i1] = -1;

   S_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd_C);

   for (i1 = 0; i1 < num_cols_offd_C; i1 ++)
      S_marker_offd[i1] = -1;

   C_diag_i = jxf_CTAlloc(JXF_Int, num_coarse+1);
   C_offd_i = jxf_CTAlloc(JXF_Int, num_coarse+1);

   /*-----------------------------------------------------------------------
    *  Loop over rows of S
    *-----------------------------------------------------------------------*/

   cnt = 0; 
   num_nonzeros_diag = 0; 
   num_nonzeros_offd = 0; 
   for (i1 = 0; i1 < num_cols_diag_S; i1 ++)
   {
      
      if (CF_marker[i1] > 0)
      {
         for (jj1 = S_diag_i[i1]; jj1 < S_diag_i[i1+1]; jj1 ++)
         {
             jcol = S_diag_j[jj1];
	     if (CF_marker[jcol] > 0)
	     {
	        S_marker[fine_to_coarse[jcol]] = i1;
	        num_nonzeros_diag ++;
	     }
         }
         for (jj1 = S_offd_i[i1]; jj1 < S_offd_i[i1+1]; jj1 ++)
         {
             jcol = S_offd_j[jj1];
	     if (CF_marker_offd[jcol] > 0)
	     {
	        S_marker_offd[map_S_to_C[jcol]] = i1;
	        num_nonzeros_offd ++;
	     }
         }
         for (jj1 = S_diag_i[i1]; jj1 < S_diag_i[i1+1]; jj1 ++)
         {
             i2 = S_diag_j[jj1];
             for (jj2 = S_diag_i[i2]; jj2 < S_diag_i[i2+1]; jj2 ++)
	     {
	        i3 = S_diag_j[jj2];
	        if (CF_marker[i3] > 0 && S_marker[fine_to_coarse[i3]] != i1)
	        {
                   S_marker[fine_to_coarse[i3]] = i1;
	           num_nonzeros_diag ++;
                }
             }
             for (jj2 = S_offd_i[i2]; jj2 < S_offd_i[i2+1]; jj2 ++)
	     {
	        i3 = S_offd_j[jj2];
	        if (CF_marker_offd[i3] > 0 && S_marker_offd[map_S_to_C[i3]] != i1)
	        {
                   S_marker_offd[map_S_to_C[i3]] = i1;
	           num_nonzeros_offd ++;
                }
             }
         }
         for (jj1 = S_offd_i[i1]; jj1 < S_offd_i[i1+1]; jj1 ++)
         {
             i2 = S_offd_j[jj1];
             for (jj2 = S_ext_diag_i[i2]; jj2 < S_ext_diag_i[i2+1]; jj2 ++)
	     {
	        i3 = S_ext_diag_j[jj2];
	        if (S_marker[i3] != i1)
	        {
                   S_marker[i3] = i1;
	           num_nonzeros_diag ++;
                }
             }
             for (jj2 = S_ext_offd_i[i2]; jj2 < S_ext_offd_i[i2+1]; jj2 ++)
	     {
	        i3 = S_ext_offd_j[jj2];
	        if (S_marker_offd[i3] != i1)
	        {
                   S_marker_offd[i3] = i1;
	           num_nonzeros_offd ++;
                }
             }
         }
         C_diag_i[++cnt] = num_nonzeros_diag;
         C_offd_i[cnt] = num_nonzeros_offd;
      }
   }

   if (num_nonzeros_diag)
   {
      C_diag_j = jxf_CTAlloc(JXF_Int, num_nonzeros_diag);
      C_diag_data = jxf_CTAlloc(JXF_Int, num_nonzeros_diag);
   }
   if (num_nonzeros_offd)
   {
      C_offd_j = jxf_CTAlloc(JXF_Int, num_nonzeros_offd);
      C_offd_data = jxf_CTAlloc(JXF_Int, num_nonzeros_offd);
   }

   for (i1 = 0; i1 < num_coarse; i1 ++)
      S_marker[i1] = -1;

   for (i1 = 0; i1 < num_cols_offd_C; i1 ++)
      S_marker_offd[i1] = -1;

   jj_count_diag = 0;
   jj_count_offd = 0;

   for (i1 = 0; i1 < num_cols_diag_S; i1 ++)
   {
      
      /*--------------------------------------------------------------------
       *  Set marker for diagonal entry, C_{i1,i1} (for square matrices). 
       *--------------------------------------------------------------------*/
 
      jj_row_begin_diag = jj_count_diag;
      jj_row_begin_offd = jj_count_offd;

      if (CF_marker[i1] > 0)
      {
         for (jj1 = S_diag_i[i1]; jj1 < S_diag_i[i1+1]; jj1 ++)
         {
             jcol = S_diag_j[jj1];
	     if (CF_marker[jcol] > 0)
	     {
	        S_marker[fine_to_coarse[jcol]] = jj_count_diag;
	        C_diag_j[jj_count_diag] = fine_to_coarse[jcol];
	        C_diag_data[jj_count_diag] = 2;
	        jj_count_diag ++;
	     }
         }
         for (jj1 = S_offd_i[i1]; jj1 < S_offd_i[i1+1]; jj1 ++)
         {
             jcol = S_offd_j[jj1];
	     if (CF_marker_offd[jcol] > 0)
	     {
	        index = map_S_to_C[jcol];
	        S_marker_offd[index] = jj_count_offd;
	        C_offd_j[jj_count_offd] = index;
	        C_offd_data[jj_count_offd] = 2;
	        jj_count_offd ++;
	     }
         }
         for (jj1 = S_diag_i[i1]; jj1 < S_diag_i[i1+1]; jj1 ++)
         {
             i2 = S_diag_j[jj1];
             for (jj2 = S_diag_i[i2]; jj2 < S_diag_i[i2+1]; jj2 ++)
	     {
	        i3 = S_diag_j[jj2];
	        if (CF_marker[i3] > 0)
	        {
		   if (S_marker[fine_to_coarse[i3]] < jj_row_begin_diag)
	           {
                      S_marker[fine_to_coarse[i3]] = jj_count_diag;
	              C_diag_j[jj_count_diag] = fine_to_coarse[i3];
	              C_diag_data[jj_count_diag] ++;
	              jj_count_diag ++;
                   }
	           else
	           {
	              C_diag_data[S_marker[fine_to_coarse[i3]]] ++;
	           }
                }
             }
             for (jj2 = S_offd_i[i2]; jj2 < S_offd_i[i2+1]; jj2 ++)
	     {
	        i3 = S_offd_j[jj2];
	        if (CF_marker_offd[i3] > 0)
                {
	           index = map_S_to_C[i3];
		   if (S_marker_offd[index] < jj_row_begin_offd)
	           {
                      S_marker_offd[index] = jj_count_offd;
	              C_offd_j[jj_count_offd] = index;
	              C_offd_data[jj_count_offd] ++;
	              jj_count_offd ++;
                   }
                   else
                   {
	              C_offd_data[S_marker_offd[index]] ++;
                   }
                }
             }
         }
         for (jj1 = S_offd_i[i1]; jj1 < S_offd_i[i1+1]; jj1 ++)
         {
             i2 = S_offd_j[jj1];
             for (jj2 = S_ext_diag_i[i2]; jj2 < S_ext_diag_i[i2+1]; jj2 ++)
	     {
	        i3 = S_ext_diag_j[jj2];
	        if (S_marker[i3] < jj_row_begin_diag)
	        {
                   S_marker[i3] = jj_count_diag;
	           C_diag_j[jj_count_diag] = i3;
	           C_diag_data[jj_count_diag] ++;
	           jj_count_diag ++;
                }
	        else
	        {
	           C_diag_data[S_marker[i3]] ++;
	        }
             }
             for (jj2 = S_ext_offd_i[i2]; jj2 < S_ext_offd_i[i2+1]; jj2 ++)
	     {
	        i3 = S_ext_offd_j[jj2];
	        if (S_marker_offd[i3] < jj_row_begin_offd)
	        {
                   S_marker_offd[i3] = jj_count_offd;
	           C_offd_j[jj_count_offd] = i3;
	           C_offd_data[jj_count_offd] ++;
	           jj_count_offd ++;
                }
                else
                {
	           C_offd_data[S_marker_offd[i3]] ++;
                }
             }
         }
      }
   }

   cnt = 0;

   for (i = 0; i < num_coarse; i ++)
   {
      for (j = C_diag_i[i]; j < C_diag_i[i+1]; j ++)
      {
         jcol = C_diag_j[j];
         if (C_diag_data[j] >= num_paths && jcol != i)
            C_diag_j[cnt++] = jcol;
      }
      C_diag_i[i] = cnt;
   }

   if (num_nonzeros_diag) jxf_TFree(C_diag_data);
   for (i = num_coarse; i > 0; i --)
      C_diag_i[i] = C_diag_i[i-1];

   C_diag_i[0] = 0;

   cnt = 0;
   for (i = 0; i < num_coarse; i ++)
   {
      for (j = C_offd_i[i]; j < C_offd_i[i+1]; j ++)
      {
         jcol = C_offd_j[j];
         if (C_offd_data[j] >= num_paths)
            C_offd_j[cnt++] = jcol;
      }
      C_offd_i[i] = cnt;
   }

   if (num_nonzeros_offd) jxf_TFree(C_offd_data);

   for (i = num_coarse; i > 0; i --)
      C_offd_i[i] = C_offd_i[i-1];

   C_offd_i[0] = 0;

   cnt = 0;
   for (i = 0; i < num_cols_diag_S; i ++)
   {
      if (CF_marker[i] > 0)
      {
         if (!(C_diag_i[cnt+1]-C_diag_i[cnt]) && !(C_offd_i[cnt+1]-C_offd_i[cnt]))
            CF_marker[i] = 2;
         cnt ++;
      }
   }

   C_diag_size = C_diag_i[num_coarse];
   C_offd_size = C_offd_i[num_coarse];

   S2 = jxf_ParCSRMatrixCreate(comm, global_num_coarse,
                              global_num_coarse, coarse_row_starts,
                              coarse_row_starts, num_cols_offd_C, C_diag_size, C_offd_size);


   jxf_ParCSRMatrixOwnsRowStarts(S2) = 0;

   C_diag = jxf_ParCSRMatrixDiag(S2);
   jxf_CSRMatrixI(C_diag) = C_diag_i;
   if (num_nonzeros_diag) jxf_CSRMatrixJ(C_diag) = C_diag_j;

   C_offd = jxf_ParCSRMatrixOffd(S2);
   jxf_CSRMatrixI(C_offd) = C_offd_i;
   jxf_ParCSRMatrixOffd(S2) = C_offd;

   if (num_cols_offd_C)
   {
      if (num_nonzeros_offd) jxf_CSRMatrixJ(C_offd) = C_offd_j;
      jxf_ParCSRMatrixColMapOffd(S2) = col_map_offd_C;
   }

   /*-----------------------------------------------------------------------
    *  Free various arrays
    *-----------------------------------------------------------------------*/

   jxf_TFree(S_marker);   
   jxf_TFree(S_marker_offd);   
   jxf_TFree(S_ext_diag_i);
   jxf_TFree(fine_to_coarse);
   if (S_ext_diag_size)
   {
      jxf_TFree(S_ext_diag_j);
   }
   jxf_TFree(S_ext_offd_i);
   if (S_ext_offd_size)
   {
      jxf_TFree(S_ext_offd_j);
   }
   if (num_cols_offd_S) 
   {
      jxf_TFree(map_S_to_C);
      jxf_TFree(CF_marker_offd);
      jxf_TFree(fine_to_coarse_offd);
   }

   *C_ptr = S2;

   return 0;
}


/*!
 * \fn JXF_Int jxf_PAMGCorrectCFMarker
 * \brief Corrects CF_marker after aggr. coarsening
 * \date 2015/08/14
 */ 
JXF_Int
jxf_PAMGCorrectCFMarker( JXF_Int *CF_marker, JXF_Int num_var, JXF_Int *new_CF_marker )
{
   JXF_Int i, cnt;

   cnt = 0;
   for (i = 0; i < num_var; i ++)
   {
      if (CF_marker[i] > 0)
      {
         if (CF_marker[i] == 1) CF_marker[i] = new_CF_marker[cnt++];
         else
         {
            CF_marker[i] = 1;
            cnt ++;
         }
      }
   }

   return 0;
}
