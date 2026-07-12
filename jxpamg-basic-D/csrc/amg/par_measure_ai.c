//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_splitting_ai.c
 *  Date: 2013/01/21
 */

#include "jx_pamg.h"

#define C_PT 1
#define F_PT -1
#define SF_PT -3
#define COMMON_C_PT 2
#define Z_PT -2

/*!
 * \fn JX_Int jx_PAMGMeasureAI
 * \brief AI Splitting routine.
 * measure_ai for point i:
 *
 *    > 0.0: ai-points
 *   =< 0.0: not ai-points
 *      (1) = 0.0: not ai-points, normal points.
 *      (2) =-1.0: not ai-points, parasitic points.
 *      (3) =-2.0: not ai-points, diagonal-dominant pionts.
 *                 including source points and isolated points.
 *      (4) other: not ai-points.
 * \date 2013/03/10
 */
JX_Int
jx_PAMGMeasureAI(jx_ParCSRMatrix *par_S,
                 jx_ParCSRMatrix *par_A,
                 JX_Int debug_flag,
                 JX_Real **measure_ai_ptr)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_S);
   jx_ParCSRCommPkg *comm_pkg = jx_ParCSRMatrixCommPkg(par_S);
   jx_ParCSRCommHandle *comm_handle = NULL;

   jx_CSRMatrix *S_diag = jx_ParCSRMatrixDiag(par_S);
   JX_Int *S_diag_i = jx_CSRMatrixI(S_diag);
   JX_Int *S_diag_j = jx_CSRMatrixJ(S_diag);

   jx_CSRMatrix *S_offd = jx_ParCSRMatrixOffd(par_S);
   JX_Int *S_offd_i = jx_CSRMatrixI(S_offd);
   JX_Int *S_offd_j = NULL;

   JX_Int num_variables = jx_CSRMatrixNumRows(S_diag);
   JX_Int num_cols_offd = 0;

   JX_Int num_sends = 0;
   JX_Int *int_buf_data;
   JX_Real *buf_data;

   JX_Real *measure_ai;
   JX_Real *measure_inf;

   JX_Int i, j;
   JX_Int index, start, my_id, num_procs;

   JX_Int ierr = 0;

   jx_ParCSRMatrix *par_ST;
   jx_CSRMatrix *ST_diag;
   JX_Int *ST_diag_i;
   JX_Int *ST_diag_j;
   jx_CSRMatrix *ST_offd;
   JX_Int *ST_offd_i;
   JX_Int *ST_offd_j;

   JX_Int nzi, nzit, jt;
   JX_Int nodep;
   JX_Int measure_dep, measure_dep_1;
   JX_Int *measure_dep_2; // number of neighbors for strong dependence each other.
                          //   JX_Real mai_1, mai_2;

   JX_Int ii;
   JX_Real buf_data_tmp;

   /***********************************************************************
    *  BEFORE THE INDEPENDENT SET COARSENING LOOP:
    *   measure_inf: calculate the measures, and communicate them
    *     (this array contains measures for both local and external nodes)
    *   CF_marker, CF_marker_offd: initialize CF_marker
    *     (separate arrays for local and external;
    *      0=unassigned, negative=F point, positive=C point)
    ***********************************************************************/

   /*-------------------------------------------------------------------------
    * Use the ParCSR strength matrix, par_S.
    *
    * For now, the "strength" of dependence/influence is defined in
    * the following way: i depends on j if
    *     aij > jx_max (k != i) aik,    aii < 0
    * or
    *     aij < jx_min (k != i) aik,    aii >= 0
    * Then S_ij = 1, else S_ij = 0.
    *
    * NOTE: S_data is not used; in stead, only strong columns are retained
    *       in S_j, which can then be used like S_data
    *-----------------------------------------------------------------------*/

   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);

   if (!comm_pkg)
   {
      comm_pkg = jx_ParCSRMatrixCommPkg(par_A);
   }

   if (!comm_pkg)
   {
      jx_MatvecCommPkgCreate(par_A);
      comm_pkg = jx_ParCSRMatrixCommPkg(par_A);
   }

   num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);

   int_buf_data = jx_CTAlloc(JX_Int, jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   buf_data = jx_CTAlloc(JX_Real, jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

   num_cols_offd = jx_CSRMatrixNumCols(S_offd);

   S_diag_j = jx_CSRMatrixJ(S_diag);

   if (num_cols_offd)
   {
      S_offd_j = jx_CSRMatrixJ(S_offd);
   }

   jx_ParCSRMatrixTranspose(par_S,
                            &par_ST,
                            0);

   ST_diag = jx_ParCSRMatrixDiag(par_ST);
   ST_diag_i = jx_CSRMatrixI(ST_diag);
   ST_diag_j = jx_CSRMatrixJ(ST_diag);

   ST_offd = jx_ParCSRMatrixOffd(par_ST);
   ST_offd_i = jx_CSRMatrixI(ST_offd);
   ST_offd_j = jx_CSRMatrixJ(ST_offd);

   /*---------------------------------------------------------------
    * Compute the measures
    *
    * The measures are currently given by the column sums of par_S.
    * Hence, measure_inf[i] is the number of influences
    * of variable i.
    *
    * The measures are augmented by a random number
    * between 0 and 1.
    *-------------------------------------------------------------*/

   measure_dep_2 = jx_CTAlloc(JX_Int, num_variables);
   measure_ai = jx_CTAlloc(JX_Real, num_variables + num_cols_offd);
   measure_inf = jx_CTAlloc(JX_Real, num_variables + num_cols_offd);

   /* first calculate the local part of the sums for the external nodes */
   for (i = 0; i < S_offd_i[num_variables]; i++)
   {
      measure_ai[num_variables + S_offd_j[i]] += 1.0;
      measure_inf[num_variables + S_offd_j[i]] += 1.0;
   }

   /* now send those locally calculated values for the external nodes to the neighboring processors */
   /* NOTE BY XU: for measure_ai, also should be considered parallel implementation. */
   if (num_procs > 1)
   {
      comm_handle = jx_ParCSRCommHandleCreate(2, comm_pkg, &measure_inf[num_variables], buf_data);
   }

   /* calculate the local part for the local nodes */
   for (i = 0; i < S_diag_i[num_variables]; i++)
   {
      measure_inf[S_diag_j[i]] += 1.0;
   }

   for (i = 0; i < num_variables; i++)
   {
      measure_ai[i] = 0.0;
   }

#if 0 
   /* calculate the local part for the local nodes */
   for (i = 0; i < num_variables; i ++) {
      for (nzi = S_diag_i[i]; nzi < S_diag_i[i+1]; nzi ++) {
          j = S_diag_j[nzi];
          nodep = 1;
          /* determine if the j depend on i */
          for (nzj = S_diag_i[j]; nzj < S_diag_i[j+1]; nzj ++) {
              if (S_diag_j[nzj] == i) {
                 nodep = 0;
                 break;
              }
          }
          if (nodep == 1) {
              measure_ai[j] += 1.0;
          }
      }
   }
#endif

#if 1
   /* calculate the local and neighbor parts using par_ST. */
   for (i = 0; i < num_variables; i++)
   {
      /* diagonal part. */
      for (nzi = S_diag_i[i]; nzi < S_diag_i[i + 1]; nzi++)
      {
         j = S_diag_j[nzi];
         nodep = 1;
         /* determine if the j depend on i */
         for (nzit = ST_diag_i[i]; nzit < ST_diag_i[i + 1]; nzit++)
         {
            jt = ST_diag_j[nzit];
            if (jt == j)
            {
               nodep = 0;
               measure_dep_2[i]++;
               // break;
            }
         }
         if (nodep == 1)
         {
            measure_ai[j] += 1.0;
         }
      }

      /* off-diagonal part. */
      for (nzi = S_offd_i[i]; nzi < S_offd_i[i + 1]; nzi++)
      {
         j = S_offd_j[nzi];
         nodep = 1;
         /* determine if the j depend on i */
         for (nzit = ST_offd_i[i]; nzit < ST_offd_i[i + 1]; nzit++)
         {
            jt = ST_offd_j[nzit];
            if (jt == j)
            {
               nodep = 0;
               measure_dep_2[i]++;
               // break;
            }
         }
         if (nodep == 1)
         {
            measure_ai[j] += 1.0;
         }
      }
   }
#endif

   /* finish the communication */
   if (num_procs > 1)
   {
      jx_ParCSRCommHandleDestroy(comm_handle);
   }

   /* need checked if num_procs. */
   /* now add the externally calculated part of the local nodes to the local nodes */
   index = 0;
   for (i = 0; i < num_sends; i++)
   {
      start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jx_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
      {
         ii = jx_ParCSRCommPkgSendMapElmt(comm_pkg, j);
         buf_data_tmp = buf_data[index++];
         measure_inf[ii] += buf_data_tmp;
      }
   }

   /* set the measures of the external nodes to zero */
   for (i = num_variables; i < num_variables + num_cols_offd; i++)
   {
      measure_ai[i] = 0.0;
      measure_inf[i] = 0.0;
   }

   /* normalizing measure_ai */
   if (debug_flag == 4)
   {
      jx_printf(" measure_ai: \n");
   }
   for (i = 0; i < num_variables; i++)
   {
      measure_dep = S_diag_i[i + 1] - S_diag_i[i];
      measure_dep = measure_dep + S_offd_i[i + 1] - S_offd_i[i];

      measure_dep_1 = measure_dep - measure_dep_2[i];
      if (measure_dep == 0)
      {
         if (measure_dep_1 != 0 || measure_dep_2[i] != 0)
         {
            jx_printf("for point %d : \n", i);
            jx_printf("measure_dep == 0, however, measure_dep_1 != 0 || measure_dep_2 !=0 \n");
            exit(0);
         }
      }

      if (measure_ai[i] > 0.0 && measure_dep > 0)
      {
#if 1
         if (measure_dep_1 == 0)
         {
            measure_ai[i] = measure_ai[i] / (measure_dep + measure_ai[i]) + 1.0;
         }
         else
         {
            measure_ai[i] = measure_ai[i] / (measure_dep + measure_ai[i]);
         }
#endif

         // measure_ai[i] = measure_ai[i]/(measure_dep+measure_ai[i]);

#if 0         
         mai_1 = measure_ai[i]/(measure_ai[i]+measure_dep);
         mai_2 = measure_ai[i]/(measure_ai[i]+measure_dep_1);
         measure_ai[i] = mai_1 + mai_2;
#endif
      }
      else if (measure_dep > 0)
      {
         measure_ai[i] = (measure_inf[i] - measure_dep) / measure_dep;
      }
      else // diagonal-dominant points: source point and isolate points.
      {
         measure_ai[i] = -2.0;
      }

      if (debug_flag == 4)
      {
         // if ( measure_ai[i] > 1.0e-10 || measure_ai[i] < -1.0e-10)
         //{
         jx_printf(" i= %d, %f \n", i, measure_ai[i]);
         //}
      }
   }

   jx_TFree(int_buf_data);
   jx_TFree(buf_data);
   jx_TFree(measure_dep_2);
   jx_TFree(measure_inf);
   jx_ParCSRMatrixDestroy(par_ST);

   *measure_ai_ptr = measure_ai;

   return (ierr);
}

// =========================================================
/* 1. 结构体比较函数：用于 qsort */
static int compare_ghost_pairs(const void *a, const void *b)
{
   typedef struct
   {
      JX_Int global_idx;
      JX_Int local_pos;
   } jx_GhostPair;

   const jx_GhostPair *pa = (const jx_GhostPair *)a;
   const jx_GhostPair *pb = (const jx_GhostPair *)b;

   if (pa->global_idx < pb->global_idx)
      return -1;
   if (pa->global_idx > pb->global_idx)
      return 1;
   return 0;
}

/* 2. 二分查找函数 */
static JX_Int binary_search_ghost_map(void *pairs_void, JX_Int size, JX_Int target)
{
   typedef struct
   {
      JX_Int global_idx;
      JX_Int local_pos;
   } jx_GhostPair;

   jx_GhostPair *pairs = (jx_GhostPair *)pairs_void;
   JX_Int low = 0, high = size - 1;

   while (low <= high)
   {
      JX_Int mid = low + (high - low) / 2;
      if (pairs[mid].global_idx == target)
         return mid;
      if (pairs[mid].global_idx < target)
         low = mid + 1;
      else
         high = mid - 1;
   }
   return -1;
}

JX_Int
jx_PAMGMeasureGAI_Base(jx_ParCSRMatrix *par_S,
                       jx_ParCSRMatrix *par_A,
                       JX_Int level,
                       JX_Int debug_flag,
                       JX_Real gai_threshold,
                       JX_Real **measure_gai_ptr,
                       JX_Int **marker_gai_ptr,
                       JX_Int **side_gai_ptr,
                       JX_Real *tau_out)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_S);
   JX_Int my_id, num_procs;
   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &num_procs);

   /* ---------------- threshold ---------------- */
   if (gai_threshold <= 0.0)
   {
      if (my_id == 0)
      {
         printf("[ERROR] jx_PAMGMeasureGAI_Base requires gai_threshold > 0.\n");
      }
      if (measure_gai_ptr)
         *measure_gai_ptr = NULL;
      if (marker_gai_ptr)
         *marker_gai_ptr = NULL;
      if (side_gai_ptr)
         *side_gai_ptr = NULL;
      if (tau_out)
         *tau_out = 0.0;
      return -1;
   }

   const JX_Real tau = gai_threshold;

   if (my_id == 0 && debug_flag)
   {
      printf("[INFO] jx_PAMGMeasureGAI_Base: fixed threshold = %.4f\n", (double)tau);
   }

   /* ---------------- S structure ---------------- */
   jx_CSRMatrix *S_diag = jx_ParCSRMatrixDiag(par_S);
   jx_CSRMatrix *S_offd = jx_ParCSRMatrixOffd(par_S);
   JX_Int *S_diag_i = jx_CSRMatrixI(S_diag);
   JX_Int *S_diag_j = jx_CSRMatrixJ(S_diag);
   JX_Int *S_offd_i = jx_CSRMatrixI(S_offd);
   JX_Int *S_offd_j = jx_CSRMatrixJ(S_offd);

   JX_Int num_variables = jx_CSRMatrixNumRows(S_diag);
   JX_Int num_cols_offd = jx_CSRMatrixNumCols(S_offd);
   JX_Int total_width = num_variables + num_cols_offd;

   /* ---------------- A values ---------------- */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(par_A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(par_A);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);

   /* ---------------- ST = S^T ---------------- */
   jx_ParCSRMatrix *par_ST = NULL;
   jx_ParCSRMatrixTranspose(par_S, &par_ST, 0);

   jx_CSRMatrix *ST_diag = jx_ParCSRMatrixDiag(par_ST);
   jx_CSRMatrix *ST_offd = jx_ParCSRMatrixOffd(par_ST);
   JX_Int *ST_diag_i = jx_CSRMatrixI(ST_diag);
   JX_Int *ST_diag_j = jx_CSRMatrixJ(ST_diag);
   JX_Int *ST_offd_i = jx_CSRMatrixI(ST_offd);
   JX_Int *ST_offd_j = jx_CSRMatrixJ(ST_offd);

   /* ---------------- ghost mapping ---------------- */
   JX_Int *ST_to_S_ghost_map = NULL;
   if (num_cols_offd > 0)
   {
      JX_Int num_ST_offd = jx_CSRMatrixNumCols(ST_offd);
      ST_to_S_ghost_map = jx_CTAlloc(JX_Int, num_ST_offd);

      JX_Int *mapS = jx_ParCSRMatrixColMapOffd(par_S);
      JX_Int *mapST = jx_ParCSRMatrixColMapOffd(par_ST);

      typedef struct
      {
         JX_Int global_idx;
         JX_Int local_pos;
      } jx_GhostPair;

      jx_GhostPair *pairs = jx_TAlloc(jx_GhostPair, num_cols_offd);

      for (JX_Int m = 0; m < num_cols_offd; m++)
      {
         pairs[m].global_idx = mapS[m];
         pairs[m].local_pos = m;
      }

      qsort(pairs, num_cols_offd, sizeof(jx_GhostPair), compare_ghost_pairs);

      for (JX_Int k = 0; k < num_ST_offd; k++)
      {
         JX_Int mid = binary_search_ghost_map(pairs, num_cols_offd, mapST[k]);
         ST_to_S_ghost_map[k] = (mid >= 0) ? pairs[mid].local_pos : -1;
      }

      jx_TFree(pairs);
   }

   /* ---------------- parameters ---------------- */
   const JX_Real GAI_ALPHA = 0.5;
   const JX_Real EPS_DEN = 1e-15;
   const JX_Real DELTA_SIDE = 0.05;

   /* ---------------- outputs ---------------- */
   JX_Real *measure_gai = jx_CTAlloc(JX_Real, num_variables);
   JX_Int *marker_gai = jx_CTAlloc(JX_Int, num_variables);
   JX_Int *side_gai = jx_CTAlloc(JX_Int, num_variables);

   /* ---------------- max row nnz ---------------- */
   JX_Int max_row_nnz = 0;
   for (JX_Int i = 0; i < num_variables; i++)
   {
      JX_Int nnz_i = (A_diag_i[i + 1] - A_diag_i[i]) +
                     (A_offd_i[i + 1] - A_offd_i[i]);
      if (nnz_i > max_row_nnz)
      {
         max_row_nnz = nnz_i;
      }
   }
   if (max_row_nnz < 16)
   {
      max_row_nnz = 16;
   }

#pragma omp parallel
   {
      JX_Real *a_row_vals = jx_CTAlloc(JX_Real, total_width);
      JX_Int *st_marker = jx_CTAlloc(JX_Int, total_width);
      JX_Int *s_marker = jx_CTAlloc(JX_Int, total_width);
      JX_Int *touched = jx_CTAlloc(JX_Int, max_row_nnz);

      for (JX_Int k = 0; k < total_width; k++)
      {
         a_row_vals[k] = 0.0;
         st_marker[k] = 0;
         s_marker[k] = 0;
      }

#pragma omp for schedule(static)
      for (JX_Int i = 0; i < num_variables; i++)
      {
         const JX_Int stamp = i + 1;
         JX_Int ntouch = 0;
         JX_Real row_sum_A = 0.0;

         /* ---------- cache A row values ---------- */
         for (JX_Int kk = A_diag_i[i]; kk < A_diag_i[i + 1]; kk++)
         {
            JX_Int col = A_diag_j[kk];
            if (col == i)
               continue;

            JX_Real v = fabs(A_diag_data[kk]);
            if (v != 0.0)
            {
               a_row_vals[col] = v;
               if (ntouch < max_row_nnz)
               {
                  touched[ntouch++] = col;
               }
               row_sum_A += v;
            }
         }

         for (JX_Int kk = A_offd_i[i]; kk < A_offd_i[i + 1]; kk++)
         {
            JX_Int col = A_offd_j[kk] + num_variables;
            JX_Real v = fabs(A_offd_data[kk]);
            if (v != 0.0)
            {
               a_row_vals[col] = v;
               if (ntouch < max_row_nnz)
               {
                  touched[ntouch++] = col;
               }
               row_sum_A += v;
            }
         }

         /* ---------- mark ST neighbors ---------- */
         for (JX_Int kk = ST_diag_i[i]; kk < ST_diag_i[i + 1]; kk++)
         {
            st_marker[ST_diag_j[kk]] = stamp;
         }

         for (JX_Int kk = ST_offd_i[i]; kk < ST_offd_i[i + 1]; kk++)
         {
            JX_Int st_col = ST_offd_j[kk];
            JX_Int s_pos = (ST_to_S_ghost_map ? ST_to_S_ghost_map[st_col] : -1);
            if (s_pos != -1)
            {
               st_marker[s_pos + num_variables] = stamp;
            }
         }

         /* ---------- mark S neighbors ---------- */
         for (JX_Int kk = S_diag_i[i]; kk < S_diag_i[i + 1]; kk++)
         {
            s_marker[S_diag_j[kk]] = stamp;
         }

         for (JX_Int kk = S_offd_i[i]; kk < S_offd_i[i + 1]; kk++)
         {
            s_marker[S_offd_j[kk] + num_variables] = stamp;
         }

         /* ---------- compute energies ---------- */
         JX_Real E_strong = 0.0;
         JX_Real E_dep = 0.0;
         JX_Real E_inf = 0.0;

         /* strong connections from row i */
         for (JX_Int kk = S_diag_i[i]; kk < S_diag_i[i + 1]; kk++)
         {
            JX_Int col = S_diag_j[kk];
            JX_Real v = a_row_vals[col];
            E_strong += v;
            if (st_marker[col] != stamp)
            {
               E_dep += v;
            }
         }

         for (JX_Int kk = S_offd_i[i]; kk < S_offd_i[i + 1]; kk++)
         {
            JX_Int col = S_offd_j[kk] + num_variables;
            JX_Real v = a_row_vals[col];
            E_strong += v;
            if (st_marker[col] != stamp)
            {
               E_dep += v;
            }
         }

         /* strong incoming connections from ST row i */
         for (JX_Int kk = ST_diag_i[i]; kk < ST_diag_i[i + 1]; kk++)
         {
            JX_Int col = ST_diag_j[kk];
            if (s_marker[col] != stamp)
            {
               E_inf += a_row_vals[col];
            }
         }

         for (JX_Int kk = ST_offd_i[i]; kk < ST_offd_i[i + 1]; kk++)
         {
            JX_Int st_col = ST_offd_j[kk];
            JX_Int s_pos = (ST_to_S_ghost_map ? ST_to_S_ghost_map[st_col] : -1);
            if (s_pos != -1)
            {
               JX_Int col = s_pos + num_variables;
               if (s_marker[col] != stamp)
               {
                  E_inf += a_row_vals[col];
               }
            }
         }

         /* ---------- normalized GAI measure ---------- */
         {
            JX_Real E_weak = row_sum_A - E_strong;
            if (E_weak < 0.0)
               E_weak = 0.0;

            JX_Real denom = E_strong + GAI_ALPHA * E_weak;

            JX_Real mu_dep = (denom > EPS_DEN) ? sqrt(E_dep / denom) : 0.0;
            JX_Real mu_inf = (denom > EPS_DEN) ? sqrt(E_inf / denom) : 0.0;
            JX_Real mu_gai = (mu_dep > mu_inf) ? mu_dep : mu_inf;

            JX_Int side_i = 0;
            if (mu_dep > (1.0 + DELTA_SIDE) * mu_inf)
            {
               side_i = -1;
            }
            else if (mu_inf > (1.0 + DELTA_SIDE) * mu_dep)
            {
               side_i = +1;
            }

            measure_gai[i] = mu_gai;
            marker_gai[i] = (mu_gai >= tau) ? 1 : 0;
            side_gai[i] = side_i;
         }

         /* ---------- clear touched ---------- */
         for (JX_Int t = 0; t < ntouch; t++)
         {
            a_row_vals[touched[t]] = 0.0;
         }
      } /* end omp for */

      jx_TFree(a_row_vals);
      jx_TFree(st_marker);
      jx_TFree(s_marker);
      jx_TFree(touched);
   } /* end omp parallel */

   /* ---------------- summary ---------------- */
   if (debug_flag)
   // if (1)
   {
      JX_Int local_gai = 0, global_gai = 0;
      for (JX_Int i = 0; i < num_variables; i++)
      {
         if (marker_gai[i])
            local_gai++;
      }

      jx_MPI_Allreduce(&local_gai, &global_gai, 1, JX_MPI_INT, MPI_SUM, comm);

      if (my_id == 0)
      {
         printf("\n--- PAMG GAI-Base (Fixed-Threshold Version) ---\n");
         printf("  [Level %d] Tau = %.4f\n", (int)level, (double)tau);
         printf("  [GAI-Point] %d\n", (int)global_gai);
         printf("------------------------------------------------\n");
      }
   }

   /* ---------------- outputs ---------------- */
   if (tau_out)
      *tau_out = tau;

   if (ST_to_S_ghost_map)
      jx_TFree(ST_to_S_ghost_map);
   if (par_ST)
      jx_ParCSRMatrixDestroy(par_ST);

   *measure_gai_ptr = measure_gai;
   *marker_gai_ptr = marker_gai;
   *side_gai_ptr = side_gai;

   return 0;
}