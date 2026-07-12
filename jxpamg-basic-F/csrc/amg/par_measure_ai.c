//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_splitting_ai.c
 *  Date: 2013/01/21
 */

#include "jxf_pamg.h"

#define C_PT 1
#define F_PT -1
#define SF_PT -3
#define COMMON_C_PT 2
#define Z_PT -2

/*!
 * \fn JXF_Int jxf_PAMGMeasureAI
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
JXF_Int
jxf_PAMGMeasureAI(jxf_ParCSRMatrix *par_S,
                 jxf_ParCSRMatrix *par_A,
                 JXF_Int debug_flag,
                 JXF_Real **measure_ai_ptr)
{
   MPI_Comm comm = jxf_ParCSRMatrixComm(par_S);
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(par_S);
   jxf_ParCSRCommHandle *comm_handle = NULL;

   jxf_CSRMatrix *S_diag = jxf_ParCSRMatrixDiag(par_S);
   JXF_Int *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int *S_diag_j = jxf_CSRMatrixJ(S_diag);

   jxf_CSRMatrix *S_offd = jxf_ParCSRMatrixOffd(par_S);
   JXF_Int *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int *S_offd_j = NULL;

   JXF_Int num_variables = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int num_cols_offd = 0;

   JXF_Int num_sends = 0;
   JXF_Int *int_buf_data;
   JXF_Real *buf_data;

   JXF_Real *measure_ai;
   JXF_Real *measure_inf;

   JXF_Int i, j;
   JXF_Int index, start, my_id, num_procs;

   JXF_Int ierr = 0;

   jxf_ParCSRMatrix *par_ST;
   jxf_CSRMatrix *ST_diag;
   JXF_Int *ST_diag_i;
   JXF_Int *ST_diag_j;
   jxf_CSRMatrix *ST_offd;
   JXF_Int *ST_offd_i;
   JXF_Int *ST_offd_j;

   JXF_Int nzi, nzit, jt;
   JXF_Int nodep;
   JXF_Int measure_dep, measure_dep_1;
   JXF_Int *measure_dep_2; // number of neighbors for strong dependence each other.
                          //   JXF_Real mai_1, mai_2;

   JXF_Int ii;
   JXF_Real buf_data_tmp;

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
    *     aij > jxf_max (k != i) aik,    aii < 0
    * or
    *     aij < jxf_min (k != i) aik,    aii >= 0
    * Then S_ij = 1, else S_ij = 0.
    *
    * NOTE: S_data is not used; in stead, only strong columns are retained
    *       in S_j, which can then be used like S_data
    *-----------------------------------------------------------------------*/

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (!comm_pkg)
   {
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A);
   }

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A);
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);

   int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

   num_cols_offd = jxf_CSRMatrixNumCols(S_offd);

   S_diag_j = jxf_CSRMatrixJ(S_diag);

   if (num_cols_offd)
   {
      S_offd_j = jxf_CSRMatrixJ(S_offd);
   }

   jxf_ParCSRMatrixTranspose(par_S,
                            &par_ST,
                            0);

   ST_diag = jxf_ParCSRMatrixDiag(par_ST);
   ST_diag_i = jxf_CSRMatrixI(ST_diag);
   ST_diag_j = jxf_CSRMatrixJ(ST_diag);

   ST_offd = jxf_ParCSRMatrixOffd(par_ST);
   ST_offd_i = jxf_CSRMatrixI(ST_offd);
   ST_offd_j = jxf_CSRMatrixJ(ST_offd);

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

   measure_dep_2 = jxf_CTAlloc(JXF_Int, num_variables);
   measure_ai = jxf_CTAlloc(JXF_Real, num_variables + num_cols_offd);
   measure_inf = jxf_CTAlloc(JXF_Real, num_variables + num_cols_offd);

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
      comm_handle = jxf_ParCSRCommHandleCreate(2, comm_pkg, &measure_inf[num_variables], buf_data);
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
      jxf_ParCSRCommHandleDestroy(comm_handle);
   }

   /* need checked if num_procs. */
   /* now add the externally calculated part of the local nodes to the local nodes */
   index = 0;
   for (i = 0; i < num_sends; i++)
   {
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
      {
         ii = jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j);
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
      jxf_printf(" measure_ai: \n");
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
            jxf_printf("for point %d : \n", i);
            jxf_printf("measure_dep == 0, however, measure_dep_1 != 0 || measure_dep_2 !=0 \n");
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
         jxf_printf(" i= %d, %f \n", i, measure_ai[i]);
         //}
      }
   }

   jxf_TFree(int_buf_data);
   jxf_TFree(buf_data);
   jxf_TFree(measure_dep_2);
   jxf_TFree(measure_inf);
   jxf_ParCSRMatrixDestroy(par_ST);

   *measure_ai_ptr = measure_ai;

   return (ierr);
}

// =========================================================
/* 1. 结构体比较函数：用于 qsort */
static int compare_ghost_pairs(const void *a, const void *b)
{
   typedef struct
   {
      JXF_Int global_idx;
      JXF_Int local_pos;
   } jxf_GhostPair;

   const jxf_GhostPair *pa = (const jxf_GhostPair *)a;
   const jxf_GhostPair *pb = (const jxf_GhostPair *)b;

   if (pa->global_idx < pb->global_idx)
      return -1;
   if (pa->global_idx > pb->global_idx)
      return 1;
   return 0;
}

/* 2. 二分查找函数 */
static JXF_Int binary_search_ghost_map(void *pairs_void, JXF_Int size, JXF_Int target)
{
   typedef struct
   {
      JXF_Int global_idx;
      JXF_Int local_pos;
   } jxf_GhostPair;

   jxf_GhostPair *pairs = (jxf_GhostPair *)pairs_void;
   JXF_Int low = 0, high = size - 1;

   while (low <= high)
   {
      JXF_Int mid = low + (high - low) / 2;
      if (pairs[mid].global_idx == target)
         return mid;
      if (pairs[mid].global_idx < target)
         low = mid + 1;
      else
         high = mid - 1;
   }
   return -1;
}

JXF_Int
jxf_PAMGMeasureGAI_Base(jxf_ParCSRMatrix *par_S,
                       jxf_ParCSRMatrix *par_A,
                       JXF_Int level,
                       JXF_Int debug_flag,
                       JXF_Real gai_threshold,
                       JXF_Real **measure_gai_ptr,
                       JXF_Int **marker_gai_ptr,
                       JXF_Int **side_gai_ptr,
                       JXF_Real *tau_out)
{
   MPI_Comm comm = jxf_ParCSRMatrixComm(par_S);
   JXF_Int my_id, num_procs;
   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   /* ---------------- threshold ---------------- */
   if (gai_threshold <= 0.0)
   {
      if (my_id == 0)
      {
         printf("[ERROR] jxf_PAMGMeasureGAI_Base requires gai_threshold > 0.\n");
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

   const JXF_Real tau = gai_threshold;

   if (my_id == 0 && debug_flag)
   {
      printf("[INFO] jxf_PAMGMeasureGAI_Base: fixed threshold = %.4f\n", (double)tau);
   }

   /* ---------------- S structure ---------------- */
   jxf_CSRMatrix *S_diag = jxf_ParCSRMatrixDiag(par_S);
   jxf_CSRMatrix *S_offd = jxf_ParCSRMatrixOffd(par_S);
   JXF_Int *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int *S_diag_j = jxf_CSRMatrixJ(S_diag);
   JXF_Int *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int *S_offd_j = jxf_CSRMatrixJ(S_offd);

   JXF_Int num_variables = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(S_offd);
   JXF_Int total_width = num_variables + num_cols_offd;

   /* ---------------- A values ---------------- */
   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(par_A);
   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(par_A);
   JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int *A_diag_j = jxf_CSRMatrixJ(A_diag);
   JXF_Real *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int *A_offd_i = jxf_CSRMatrixI(A_offd);
   JXF_Int *A_offd_j = jxf_CSRMatrixJ(A_offd);

   /* ---------------- ST = S^T ---------------- */
   jxf_ParCSRMatrix *par_ST = NULL;
   jxf_ParCSRMatrixTranspose(par_S, &par_ST, 0);

   jxf_CSRMatrix *ST_diag = jxf_ParCSRMatrixDiag(par_ST);
   jxf_CSRMatrix *ST_offd = jxf_ParCSRMatrixOffd(par_ST);
   JXF_Int *ST_diag_i = jxf_CSRMatrixI(ST_diag);
   JXF_Int *ST_diag_j = jxf_CSRMatrixJ(ST_diag);
   JXF_Int *ST_offd_i = jxf_CSRMatrixI(ST_offd);
   JXF_Int *ST_offd_j = jxf_CSRMatrixJ(ST_offd);

   /* ---------------- ghost mapping ---------------- */
   JXF_Int *ST_to_S_ghost_map = NULL;
   if (num_cols_offd > 0)
   {
      JXF_Int num_ST_offd = jxf_CSRMatrixNumCols(ST_offd);
      ST_to_S_ghost_map = jxf_CTAlloc(JXF_Int, num_ST_offd);

      JXF_Int *mapS = jxf_ParCSRMatrixColMapOffd(par_S);
      JXF_Int *mapST = jxf_ParCSRMatrixColMapOffd(par_ST);

      typedef struct
      {
         JXF_Int global_idx;
         JXF_Int local_pos;
      } jxf_GhostPair;

      jxf_GhostPair *pairs = jxf_TAlloc(jxf_GhostPair, num_cols_offd);

      for (JXF_Int m = 0; m < num_cols_offd; m++)
      {
         pairs[m].global_idx = mapS[m];
         pairs[m].local_pos = m;
      }

      qsort(pairs, num_cols_offd, sizeof(jxf_GhostPair), compare_ghost_pairs);

      for (JXF_Int k = 0; k < num_ST_offd; k++)
      {
         JXF_Int mid = binary_search_ghost_map(pairs, num_cols_offd, mapST[k]);
         ST_to_S_ghost_map[k] = (mid >= 0) ? pairs[mid].local_pos : -1;
      }

      jxf_TFree(pairs);
   }

   /* ---------------- parameters ---------------- */
   const JXF_Real GAI_ALPHA = 0.5;
   const JXF_Real EPS_DEN = 1e-15;
   const JXF_Real DELTA_SIDE = 0.05;

   /* ---------------- outputs ---------------- */
   JXF_Real *measure_gai = jxf_CTAlloc(JXF_Real, num_variables);
   JXF_Int *marker_gai = jxf_CTAlloc(JXF_Int, num_variables);
   JXF_Int *side_gai = jxf_CTAlloc(JXF_Int, num_variables);

   /* ---------------- max row nnz ---------------- */
   JXF_Int max_row_nnz = 0;
   for (JXF_Int i = 0; i < num_variables; i++)
   {
      JXF_Int nnz_i = (A_diag_i[i + 1] - A_diag_i[i]) +
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
      JXF_Real *a_row_vals = jxf_CTAlloc(JXF_Real, total_width);
      JXF_Int *st_marker = jxf_CTAlloc(JXF_Int, total_width);
      JXF_Int *s_marker = jxf_CTAlloc(JXF_Int, total_width);
      JXF_Int *touched = jxf_CTAlloc(JXF_Int, max_row_nnz);

      for (JXF_Int k = 0; k < total_width; k++)
      {
         a_row_vals[k] = 0.0;
         st_marker[k] = 0;
         s_marker[k] = 0;
      }

#pragma omp for schedule(static)
      for (JXF_Int i = 0; i < num_variables; i++)
      {
         const JXF_Int stamp = i + 1;
         JXF_Int ntouch = 0;
         JXF_Real row_sum_A = 0.0;

         /* ---------- cache A row values ---------- */
         for (JXF_Int kk = A_diag_i[i]; kk < A_diag_i[i + 1]; kk++)
         {
            JXF_Int col = A_diag_j[kk];
            if (col == i)
               continue;

            JXF_Real v = fabs(A_diag_data[kk]);
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

         for (JXF_Int kk = A_offd_i[i]; kk < A_offd_i[i + 1]; kk++)
         {
            JXF_Int col = A_offd_j[kk] + num_variables;
            JXF_Real v = fabs(A_offd_data[kk]);
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
         for (JXF_Int kk = ST_diag_i[i]; kk < ST_diag_i[i + 1]; kk++)
         {
            st_marker[ST_diag_j[kk]] = stamp;
         }

         for (JXF_Int kk = ST_offd_i[i]; kk < ST_offd_i[i + 1]; kk++)
         {
            JXF_Int st_col = ST_offd_j[kk];
            JXF_Int s_pos = (ST_to_S_ghost_map ? ST_to_S_ghost_map[st_col] : -1);
            if (s_pos != -1)
            {
               st_marker[s_pos + num_variables] = stamp;
            }
         }

         /* ---------- mark S neighbors ---------- */
         for (JXF_Int kk = S_diag_i[i]; kk < S_diag_i[i + 1]; kk++)
         {
            s_marker[S_diag_j[kk]] = stamp;
         }

         for (JXF_Int kk = S_offd_i[i]; kk < S_offd_i[i + 1]; kk++)
         {
            s_marker[S_offd_j[kk] + num_variables] = stamp;
         }

         /* ---------- compute energies ---------- */
         JXF_Real E_strong = 0.0;
         JXF_Real E_dep = 0.0;
         JXF_Real E_inf = 0.0;

         /* strong connections from row i */
         for (JXF_Int kk = S_diag_i[i]; kk < S_diag_i[i + 1]; kk++)
         {
            JXF_Int col = S_diag_j[kk];
            JXF_Real v = a_row_vals[col];
            E_strong += v;
            if (st_marker[col] != stamp)
            {
               E_dep += v;
            }
         }

         for (JXF_Int kk = S_offd_i[i]; kk < S_offd_i[i + 1]; kk++)
         {
            JXF_Int col = S_offd_j[kk] + num_variables;
            JXF_Real v = a_row_vals[col];
            E_strong += v;
            if (st_marker[col] != stamp)
            {
               E_dep += v;
            }
         }

         /* strong incoming connections from ST row i */
         for (JXF_Int kk = ST_diag_i[i]; kk < ST_diag_i[i + 1]; kk++)
         {
            JXF_Int col = ST_diag_j[kk];
            if (s_marker[col] != stamp)
            {
               E_inf += a_row_vals[col];
            }
         }

         for (JXF_Int kk = ST_offd_i[i]; kk < ST_offd_i[i + 1]; kk++)
         {
            JXF_Int st_col = ST_offd_j[kk];
            JXF_Int s_pos = (ST_to_S_ghost_map ? ST_to_S_ghost_map[st_col] : -1);
            if (s_pos != -1)
            {
               JXF_Int col = s_pos + num_variables;
               if (s_marker[col] != stamp)
               {
                  E_inf += a_row_vals[col];
               }
            }
         }

         /* ---------- normalized GAI measure ---------- */
         {
            JXF_Real E_weak = row_sum_A - E_strong;
            if (E_weak < 0.0)
               E_weak = 0.0;

            JXF_Real denom = E_strong + GAI_ALPHA * E_weak;

            JXF_Real mu_dep = (denom > EPS_DEN) ? sqrt(E_dep / denom) : 0.0;
            JXF_Real mu_inf = (denom > EPS_DEN) ? sqrt(E_inf / denom) : 0.0;
            JXF_Real mu_gai = (mu_dep > mu_inf) ? mu_dep : mu_inf;

            JXF_Int side_i = 0;
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
         for (JXF_Int t = 0; t < ntouch; t++)
         {
            a_row_vals[touched[t]] = 0.0;
         }
      } /* end omp for */

      jxf_TFree(a_row_vals);
      jxf_TFree(st_marker);
      jxf_TFree(s_marker);
      jxf_TFree(touched);
   } /* end omp parallel */

   /* ---------------- summary ---------------- */
   // if (debug_flag)
   if (1)
   {
      JXF_Int local_gai = 0, global_gai = 0;
      for (JXF_Int i = 0; i < num_variables; i++)
      {
         if (marker_gai[i])
            local_gai++;
      }

      jxf_MPI_Allreduce(&local_gai, &global_gai, 1, JXF_MPI_INT, MPI_SUM, comm);

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
      jxf_TFree(ST_to_S_ghost_map);
   if (par_ST)
      jxf_ParCSRMatrixDestroy(par_ST);

   *measure_gai_ptr = measure_gai;
   *marker_gai_ptr = marker_gai;
   *side_gai_ptr = side_gai;

   return 0;
}