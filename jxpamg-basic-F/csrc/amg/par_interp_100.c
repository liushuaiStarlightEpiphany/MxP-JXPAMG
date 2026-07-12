//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_interp_100.c
 *  Date: 2018/05/22
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGBuildInterpOnePnt
 * \date 2018/05/22
 */
JXF_Int
jxf_PAMGBuildInterpOnePnt( jxf_ParCSRMatrix  *A,
                          JXF_Int           *CF_marker,
                          jxf_ParCSRMatrix  *S,
                          JXF_Int           *num_cpts_global,
                          JXF_Int            num_functions,
                          JXF_Int           *dof_func,
                          JXF_Int            debug_flag,
                          JXF_Int           *col_offd_S_to_A,
                          jxf_ParCSRMatrix **P_ptr )
{
   MPI_Comm 	         comm     = jxf_ParCSRMatrixComm(A);   
   jxf_ParCSRCommPkg     *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_ParCSRCommHandle  *comm_handle;

   jxf_CSRMatrix         *A_diag      = jxf_ParCSRMatrixDiag(A);
   JXF_Real              *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int               *A_diag_i    = jxf_CSRMatrixI(A_diag);
   JXF_Int               *A_diag_j    = jxf_CSRMatrixJ(A_diag);

   jxf_CSRMatrix         *A_offd      = jxf_ParCSRMatrixOffd(A);   
   JXF_Real              *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int               *A_offd_i    = jxf_CSRMatrixI(A_offd);
   JXF_Int               *A_offd_j    = jxf_CSRMatrixJ(A_offd);

   JXF_Int                num_cols_A_offd = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_CSRMatrix         *S_diag   = jxf_ParCSRMatrixDiag(S);
   JXF_Int               *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int               *S_diag_j = jxf_CSRMatrixJ(S_diag);

   jxf_CSRMatrix         *S_offd   = jxf_ParCSRMatrixOffd(S);   
   JXF_Int               *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int               *S_offd_j = jxf_CSRMatrixJ(S_offd);

   /* Interpolation matrix P */
   jxf_ParCSRMatrix      *P;
   /* csr's */
   jxf_CSRMatrix    *P_diag;
   jxf_CSRMatrix    *P_offd;
   /* arrays */
   JXF_Real         *P_diag_data;
   JXF_Int          *P_diag_i;
   JXF_Int          *P_diag_j;
   JXF_Real         *P_offd_data;
   JXF_Int          *P_offd_i;
   JXF_Int          *P_offd_j;
   JXF_Int           num_cols_offd_P;
   JXF_Int          *col_map_offd_P;
   /* CF marker off-diag part */
   JXF_Int          *CF_marker_offd = NULL;
   /* func type off-diag part */
   JXF_Int          *dof_func_offd  = NULL;
   /* nnz */
   JXF_Int           nnz_diag, nnz_offd, cnt_diag, cnt_offd;
   JXF_Int          *marker_diag, *marker_offd;
   /* local size */
   JXF_Int           n_fine = jxf_CSRMatrixNumRows(A_diag);
   /* number of C-pts */
   JXF_Int           n_cpts = 0;
   /* fine to coarse mapping: diag part and offd part */
   JXF_Int          *fine_to_coarse;
   JXF_Int          *fine_to_coarse_offd;
   JXF_Int           total_global_cpts, my_first_cpt;
   JXF_Int           my_id, num_procs;
   JXF_Int           num_sends, *int_buf_data;

   JXF_Int           i, j, i1, j1, k1, index, start;
   JXF_Int          *max_abs_cij;
   char               *max_abs_diag_offd;
   JXF_Real          max_abs_aij, vv;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION
   my_first_cpt = num_cpts_global[0];
   if (my_id == (num_procs -1)) total_global_cpts = num_cpts_global[1];
   jxf_MPI_Bcast(&total_global_cpts, 1, JXF_MPI_INT, num_procs-1, comm);
#else
   my_first_cpt = num_cpts_global[my_id];
   total_global_cpts = num_cpts_global[num_procs];
#endif

   /*-------------------------------------------------------------------
    * Get the CF_marker data for the off-processor columns
    *-------------------------------------------------------------------*/
   /* CF marker for the off-diag columns */
   if (num_cols_A_offd)
   {
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);
   }
   /* function type indicator for the off-diag columns */
   if (num_functions > 1 && num_cols_A_offd)
   {
      dof_func_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);
   }
   /* if CommPkg of A is not present, create it */
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A); 
   }
   /* number of sends to do (number of procs) */
   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   /* send buffer, of size send_map_starts[num_sends]), 
    * i.e., number of entries to send */
   int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   
   /* copy CF markers of elements to send to buffer 
    * RL: why copy them with two for loops? Why not just loop through all in one */
   index = 0;
   for (i = 0; i < num_sends; i++)
   {
      /* start pos of elements sent to send_proc[i] */
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      /* loop through all elems to send_proc[i] */
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j++)
      {
         /* CF marker of send_map_elemts[j] */
         int_buf_data[index++] = CF_marker[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
      }
   }
   /* create a handle to start communication. 11: for integer */
   comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);
   /* destroy the handle to finish communication */
   jxf_ParCSRCommHandleDestroy(comm_handle);

   /* do a similar communication for dof_func */
   if (num_functions > 1)
   {
      index = 0;
      for (i = 0; i < num_sends; i++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j++)
         {
            int_buf_data[index++] = dof_func[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
         }
      }
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, dof_func_offd);
      jxf_ParCSRCommHandleDestroy(comm_handle);   
   }

   /*-----------------------------------------------------------------------
    *  First Pass: Determine size of P and fill in fine_to_coarse mapping,
    *  and find the most strongly influencing C-pt for each F-pt
    *-----------------------------------------------------------------------*/
   /* nnz in diag and offd parts */
   cnt_diag = 0;
   cnt_offd = 0;
   max_abs_cij       = jxf_CTAlloc(JXF_Int, n_fine);
   max_abs_diag_offd = jxf_CTAlloc(char, n_fine);
   fine_to_coarse    = jxf_CTAlloc(JXF_Int, n_fine);

   /* markers initialized as zeros */
   marker_diag = jxf_CTAlloc(JXF_Int, n_fine);
   marker_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);

   for (i = 0; i < n_fine; i++)
   {
      /*--------------------------------------------------------------------
       *  If i is a C-point, interpolation is the identity. Also set up
       *  mapping vector.
       *--------------------------------------------------------------------*/
      if (CF_marker[i] >= 0)
      {
         fine_to_coarse[i] = my_first_cpt + n_cpts;
         n_cpts++;
         continue;
      }

      /* mark all the strong connections: in S */
      JXF_Int MARK = i + 1;
      /* loop through row i of S, diag part  */
      for (j = S_diag_i[i]; j < S_diag_i[i+1]; j++)
      {
         marker_diag[S_diag_j[j]] = MARK;
      }
      /* loop through row i of S, offd part  */
      if (num_procs > 1)
      {
          for (j = S_offd_i[i]; j < S_offd_i[i+1]; j++)
          {
             j1 = col_offd_S_to_A ? col_offd_S_to_A[S_offd_j[j]] : S_offd_j[j];
             marker_offd[j1] = MARK;
          }
      }

      fine_to_coarse[i] = -1;
      /*---------------------------------------------------------------------------
       *  If i is an F-pt, interpolation is from the most strongly influencing C-pt
       *  Find this C-pt and save it
       *--------------------------------------------------------------------------*/
      /* if we failed to find any strong C-pt, mark this point as an 'n' */
      char marker = 'n';
      /* max abs val */
      max_abs_aij = -1.0;
      /* loop through row i of A, diag part  */
      for (j = A_diag_i[i]; j < A_diag_i[i+1]; j++)
      {
         i1 = A_diag_j[j];
         vv = fabs(A_diag_data[j]);
#if 0
         /* !!! this is a hack just for code verification purpose !!!
            it basically says:
            1. if we see |a_ij| < 1e-14, force it to be 1e-14
            2. if we see |a_ij| == the max(|a_ij|) so far exactly, 
               replace it if the j idx is smaller
            Reasons:
            1. numerical round-off for eps-level values
            2. entries in CSR rows may be listed in different orders
         */
         vv = vv < 1e-14 ? 1e-14 : vv;
         if (CF_marker[i1] >= 0 && marker_diag[i1] == MARK && 
             vv == max_abs_aij && i1 < max_abs_cij[i])
         {
            /* mark it as a 'd' */
            marker         = 'd';
            max_abs_cij[i] = i1;
            max_abs_aij    = vv;
            continue;
         }
#endif
         /* it is a strong C-pt and has abs val larger than what have seen */
         if (CF_marker[i1] >= 0 && marker_diag[i1] == MARK && vv > max_abs_aij)
         {
            /* mark it as a 'd' */
            marker         = 'd';
            max_abs_cij[i] = i1;
            max_abs_aij    = vv;
         }
      }
      /* offd part */
      if (num_procs > 1)
      {
         for (j = A_offd_i[i]; j < A_offd_i[i+1]; j++)
         {
            i1 = A_offd_j[j];
            vv = fabs(A_offd_data[j]);
            if (CF_marker_offd[i1] >= 0 && marker_offd[i1] == MARK && vv > max_abs_aij)
            {
               /* mark it as an 'o' */
               marker         = 'o';
               max_abs_cij[i] = i1;
               max_abs_aij    = vv;
            }
         }
      }

      max_abs_diag_offd[i] = marker;
      
      if (marker == 'd')
      {
         cnt_diag ++;
      }
      else if (marker == 'o')
      {
         cnt_offd ++;
      }
   }

   nnz_diag = cnt_diag + n_cpts;
   nnz_offd = cnt_offd;
   
   /*------------- allocate arrays */
   P_diag_i    = jxf_CTAlloc(JXF_Int,  n_fine+1);
   P_diag_j    = jxf_CTAlloc(JXF_Int,  nnz_diag);
   P_diag_data = jxf_CTAlloc(JXF_Real, nnz_diag);

   /* not in ``if num_procs > 1'', 
    * allocation needed even for empty CSR */
   P_offd_i    = jxf_CTAlloc(JXF_Int,  n_fine+1);
   P_offd_j    = jxf_CTAlloc(JXF_Int,  nnz_offd);
   P_offd_data = jxf_CTAlloc(JXF_Real, nnz_offd);

   /* redundant */
   P_diag_i[0] = 0;
   P_offd_i[0] = 0;

   /* reset counters */
   cnt_diag = 0;
   cnt_offd = 0;

   /*-----------------------------------------------------------------------
    *  Send and receive fine_to_coarse info.
    *-----------------------------------------------------------------------*/ 
   fine_to_coarse_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd); 
   index = 0;
   for (i = 0; i < num_sends; i++)
   {
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j++)
      {
         int_buf_data[index++] = fine_to_coarse[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
      }
   }
   comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, fine_to_coarse_offd);  
   jxf_ParCSRCommHandleDestroy(comm_handle);   

   /*-----------------------------------------------------------------------
    *  Second Pass: Populate P
    *-----------------------------------------------------------------------*/
   for (i = 0; i < n_fine; i++)
   {
      if (CF_marker[i] >= 0)
      {
         /*--------------------------------------------------------------------
          *  If i is a C-point, interpolation is the identity. 
          *--------------------------------------------------------------------*/
         P_diag_j[cnt_diag] = fine_to_coarse[i] - my_first_cpt;
         P_diag_data[cnt_diag++] = 1.0;
      }
      else
      {
         /*---------------------------------------------------------------------------
          *  If i is an F-pt, interpolation is from the most strongly influencing C-pt
          *--------------------------------------------------------------------------*/
         if (max_abs_diag_offd[i] == 'd')
         {
            /* on diag part of P */
            j = max_abs_cij[i];
            P_diag_j[cnt_diag] = fine_to_coarse[j] - my_first_cpt;
            P_diag_data[cnt_diag++] = 1.0;
         }
         else if (max_abs_diag_offd[i] == 'o')
         {
            /* on offd part of P */
            j = max_abs_cij[i];
            P_offd_j[cnt_offd] = j;
            P_offd_data[cnt_offd++] = 1.0;
         }
      }

      P_diag_i[i+1] = cnt_diag;
      P_offd_i[i+1] = cnt_offd;
   }

   jxf_assert(cnt_diag == nnz_diag);
   jxf_assert(cnt_offd == nnz_offd);
   
   /* num of cols in the offd part of P */
   num_cols_offd_P = 0;
   
   /* marker_offd: all -1 */
   for (i = 0; i < num_cols_A_offd; i++)
   {
      marker_offd[i] = -1;
   }
   for (i = 0; i < nnz_offd; i++)
   {
      i1 = P_offd_j[i];
      if (marker_offd[i1] == -1)
      {
         num_cols_offd_P++;
         marker_offd[i1] = 1;
      }
   }

   /* col_map_offd_P: the col indices of the offd of P
    * we first keep them be the offd-idx of A */
   col_map_offd_P = jxf_CTAlloc(JXF_Int, num_cols_offd_P);
   for (i = 0, i1 = 0; i < num_cols_A_offd; i++)
   {
      if (marker_offd[i] == 1)
      {
         col_map_offd_P[i1++] = i;
      }
   }
   jxf_assert(i1 == num_cols_offd_P);

   /* now, adjust P_offd_j to local idx w.r.t col_map_offd_R
    * by searching */
   for (i = 0; i < nnz_offd; i++)
   {
      i1 = P_offd_j[i];
      k1 = jxf_BinarySearch(col_map_offd_P, i1, num_cols_offd_P);
      /* search must succeed */
      jxf_assert(k1 >= 0 && k1 < num_cols_offd_P);
      P_offd_j[i] = k1;
   }

   /* change col_map_offd_P to global coarse ids */
   for (i = 0; i < num_cols_offd_P; i++)
   {
      col_map_offd_P[i] = fine_to_coarse_offd[col_map_offd_P[i]];
   }

   /* Now, we should have everything of Parcsr matrix P */
   P = jxf_ParCSRMatrixCreate(comm,
                                jxf_ParCSRMatrixGlobalNumCols(A), /* global num of rows */
                                total_global_cpts, /* global num of cols */
                                jxf_ParCSRMatrixColStarts(A), /* row_starts */
                                num_cpts_global, /* col_starts */
                                num_cols_offd_P, /* num cols offd */
                                nnz_diag,
                                nnz_offd);

   P_diag = jxf_ParCSRMatrixDiag(P);
   jxf_CSRMatrixData(P_diag) = P_diag_data;
   jxf_CSRMatrixI(P_diag)    = P_diag_i;
   jxf_CSRMatrixJ(P_diag)    = P_diag_j;

   P_offd = jxf_ParCSRMatrixOffd(P);
   jxf_CSRMatrixData(P_offd) = P_offd_data;
   jxf_CSRMatrixI(P_offd)    = P_offd_i;
   jxf_CSRMatrixJ(P_offd)    = P_offd_j;
   /* P does not own ColStarts, since A does */
   jxf_ParCSRMatrixOwnsRowStarts(P) = 0;
   
   jxf_ParCSRMatrixColMapOffd(P) = col_map_offd_P;

   /* create CommPkg of P */
   jxf_MatvecCommPkgCreate(P);

   *P_ptr = P;

   /* free workspace */
   jxf_TFree(CF_marker_offd);
   jxf_TFree(dof_func_offd);
   jxf_TFree(int_buf_data);
   jxf_TFree(fine_to_coarse);
   jxf_TFree(fine_to_coarse_offd);
   jxf_TFree(marker_diag);
   jxf_TFree(marker_offd);
   jxf_TFree(max_abs_cij);
   jxf_TFree(max_abs_diag_offd);

   return 0;
}
