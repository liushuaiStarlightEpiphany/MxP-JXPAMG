//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu_setup.c
 *
 */

#include "jx_ilu.h"
#include "jx_krylov.h"

typedef struct
{
   JX_Int index;
   JX_Real score;
} jx_ILULRCScore;

static int
jx_ILULRCScoreCompare(const void *a,
                      const void *b)
{
   const jx_ILULRCScore *sa = (const jx_ILULRCScore *)a;
   const jx_ILULRCScore *sb = (const jx_ILULRCScore *)b;

   if (sa->score < sb->score)
   {
      return 1;
   }
   if (sa->score > sb->score)
   {
      return -1;
   }
   return (sa->index > sb->index) - (sa->index < sb->index);
}

static int
jx_ILULRCIntCompare(const void *a,
                    const void *b)
{
   const JX_Int ia = *(const JX_Int *)a;
   const JX_Int ib = *(const JX_Int *)b;

   return (ia > ib) - (ia < ib);
}

static JX_Int
jx_ILUBuildLocalResidualCorrectionMask(jx_ParCSRMatrix *A,
                                       JX_Int *perm,
                                       JX_Int nLU,
                                       jx_ParCSRMatrix *L,
                                       JX_Real *D,
                                       jx_ParCSRMatrix *U,
                                       JX_Real threshold,
                                       JX_Int print_level,
                                       JX_Int **marker_ptr,
                                       JX_Int *num_selected_ptr)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   jx_ParVector *probe_rhs = NULL;
   jx_ParVector *probe_sol = NULL;
   jx_ParVector *probe_res = NULL;
   jx_ParVector *probe_work = NULL;
   jx_Vector *res_local;
   JX_Real *res_data;
   JX_Int *marker = NULL;
   JX_Real local_max = 0.0;
   JX_Real global_max = 0.0;
   JX_Real cutoff;
   JX_Int local_selected = 0;
   JX_Int global_selected = 0;
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int i, my_id;

   if (!marker_ptr || !num_selected_ptr)
   {
      return jx_error_flag;
   }

   *marker_ptr = NULL;
   *num_selected_ptr = 0;

   probe_rhs = jx_ParVectorCreate(comm,
                                  jx_ParCSRMatrixGlobalNumRows(A),
                                  jx_ParCSRMatrixRowStarts(A));
   probe_sol = jx_ParVectorCreate(comm,
                                  jx_ParCSRMatrixGlobalNumRows(A),
                                  jx_ParCSRMatrixRowStarts(A));
   probe_res = jx_ParVectorCreate(comm,
                                  jx_ParCSRMatrixGlobalNumRows(A),
                                  jx_ParCSRMatrixRowStarts(A));
   probe_work = jx_ParVectorCreate(comm,
                                   jx_ParCSRMatrixGlobalNumRows(A),
                                   jx_ParCSRMatrixRowStarts(A));

   jx_ParVectorInitialize(probe_rhs);
   jx_ParVectorInitialize(probe_sol);
   jx_ParVectorInitialize(probe_res);
   jx_ParVectorInitialize(probe_work);

   jx_ParVectorSetConstantValues(probe_rhs, 1.0);
   jx_ParVectorSetConstantValues(probe_sol, 0.0);
   jx_ParVectorSetConstantValues(probe_res, 0.0);
   jx_ParVectorSetConstantValues(probe_work, 0.0);

   jx_ILUSolveLU(A, probe_rhs, probe_sol, perm, nLU, L, D, U,
                 probe_res, probe_work);
   jx_ParCSRMatrixMatvecOutOfPlace(-1.0, A, probe_sol, 1.0,
                                   probe_rhs, probe_res);

   res_local = jx_ParVectorLocalVector(probe_res);
   res_data = jx_VectorData(res_local);
   for (i = 0; i < n; i++)
   {
      if (jx_abs(res_data[i]) > local_max)
      {
         local_max = jx_abs(res_data[i]);
      }
   }

   jx_MPI_Allreduce(&local_max, &global_max, 1, JX_MPI_REAL, MPI_MAX, comm);

   marker = jx_CTAlloc(JX_Int, n);
   if (global_max > JX_REAL_EPSILON)
   {
      if (threshold < 0.0)
      {
         threshold = 0.0;
      }
      cutoff = threshold * global_max;
      for (i = 0; i < n; i++)
      {
         if (jx_abs(res_data[i]) >= cutoff)
         {
            marker[i] = 1;
            local_selected++;
         }
      }
   }

   jx_MPI_Allreduce(&local_selected, &global_selected, 1, JX_MPI_INT, MPI_SUM, comm);
   jx_MPI_Comm_rank(comm, &my_id);
   if ((my_id == 0) && (print_level > 0))
   {
      jx_printf("ILU LRC fixed mask: threshold = %e, probe max residual = %e, selected rows = %d\n",
                threshold, global_max, global_selected);
   }

   *marker_ptr = marker;
   *num_selected_ptr = global_selected;

   jx_ParVectorDestroy(probe_rhs);
   jx_ParVectorDestroy(probe_sol);
   jx_ParVectorDestroy(probe_res);
   jx_ParVectorDestroy(probe_work);

   return jx_error_flag;
}

static JX_Int
jx_ILUBuildTopRatioLocalResidualCorrection(jx_ParCSRMatrix *A,
                                           jx_ParVector *probe_rhs,
                                           JX_Int *perm,
                                           JX_Int nLU,
                                           jx_ParCSRMatrix *L,
                                           JX_Real *D,
                                           jx_ParCSRMatrix *U,
                                           JX_Int base_max_iter,
                                           JX_Real ratio,
                                           JX_Int select_type,
                                           JX_Real energy_target,
                                           JX_Int print_level,
                                           JX_Int **marker_ptr,
                                           JX_Int **core_indices_ptr,
                                           JX_Int *core_size_ptr,
                                           JX_Int **local_indices_ptr,
                                           JX_Int **local_map_ptr,
                                           jx_CSRMatrix **local_mat_ptr,
                                           JX_Int **local_indexD_ptr,
                                           JX_Int **local_indexLU_ptr,
                                           JX_Real **local_valueLU_ptr,
                                           JX_Real **local_rhs_ptr,
                                           JX_Real **local_sol_ptr,
                                           JX_Real **local_work_ptr,
                                           JX_Int *num_selected_ptr,
                                           JX_Real *global_core_size_ptr,
                                           JX_Real *global_expanded_size_ptr,
                                           JX_Real *global_local_mat_nnz_ptr,
                                           JX_Real *global_local_ilu_nnz_ptr,
                                           JX_Real *probe_core_energy_ratio_ptr,
                                           JX_Real *probe_expanded_energy_ratio_ptr,
                                           JX_Real *offd_coupling_ratio_ptr)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   jx_ParVector *probe_sol = NULL;
   jx_ParVector *probe_res = NULL;
   jx_ParVector *probe_work = NULL;
   jx_Vector *res_local;
   JX_Real *res_data;
   jx_ILULRCScore *scores = NULL;
   JX_Int *marker = NULL;
   JX_Int *core_marker = NULL;
   JX_Int *core_indices = NULL;
   JX_Int *local_indices = NULL;
   JX_Int *local_map = NULL;
   jx_CSRMatrix *local_mat = NULL;
   JX_Int *B_i = NULL;
   JX_Int *B_j = NULL;
   JX_Real *B_data = NULL;
   JX_Int *indexD = NULL;
   JX_Int *indexLU = NULL;
   JX_Real *valueLU = NULL;
   JX_Real *local_rhs = NULL;
   JX_Real *local_sol = NULL;
   JX_Real *local_work = NULL;
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int target, expanded, local_nnz;
   JX_Int global_core_int, global_expanded_int;
   JX_Int i, j, k, old_row, new_col, diag_pos, my_id;
   JX_Int probe_iter;
   JX_Int local_bad_pivots, global_bad_pivots;
   JX_Int curve_id, curve_target;
   JX_Int *curve_marker = NULL;
   JX_Real local_nnz_real, local_ilu_nnz_real;
   JX_Real global_local_mat_nnz, global_local_ilu_nnz;
   JX_Real diag_value;
   JX_Real inv_diag_abs, pivot_abs, local_min_pivot, global_min_pivot;
   JX_Real local_total_energy, local_core_energy, local_expanded_energy;
   JX_Real global_total_energy, global_core_energy, global_expanded_energy;
   JX_Real local_diag_abs, local_offd_abs, global_diag_abs, global_offd_abs;
   JX_Real running_energy;
   JX_Real curve_ratios[5] = {0.05, 0.10, 0.20, 0.30, 0.40};
   JX_Real local_curve_core[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
   JX_Real local_curve_expanded[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
   JX_Real global_curve_core[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
   JX_Real global_curve_expanded[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

   *marker_ptr = NULL;
   *core_indices_ptr = NULL;
   *core_size_ptr = 0;
   *local_indices_ptr = NULL;
   *local_map_ptr = NULL;
   *local_mat_ptr = NULL;
   *local_indexD_ptr = NULL;
   *local_indexLU_ptr = NULL;
   *local_valueLU_ptr = NULL;
   *local_rhs_ptr = NULL;
   *local_sol_ptr = NULL;
   *local_work_ptr = NULL;
   *num_selected_ptr = 0;
   *global_core_size_ptr = 0.0;
   *global_expanded_size_ptr = 0.0;
   *global_local_mat_nnz_ptr = 0.0;
   *global_local_ilu_nnz_ptr = 0.0;
   *probe_core_energy_ratio_ptr = 0.0;
   *probe_expanded_energy_ratio_ptr = 0.0;
   *offd_coupling_ratio_ptr = 0.0;

   if (ratio < 0.0)
   {
      ratio = 0.0;
   }
   if (ratio > 1.0)
   {
      ratio = 1.0;
   }
   if (select_type != 1)
   {
      select_type = 0;
   }
   if (energy_target < 0.0)
   {
      energy_target = 0.0;
   }
   if (energy_target > 1.0)
   {
      energy_target = 1.0;
   }
   if (base_max_iter < 1)
   {
      base_max_iter = 1;
   }

   marker = jx_CTAlloc(JX_Int, n);
   core_marker = jx_CTAlloc(JX_Int, n);
   local_map = jx_TAlloc(JX_Int, n);
   for (i = 0; i < n; i++)
   {
      local_map[i] = -1;
   }

   if (target == 0)
   {
      *marker_ptr = marker;
      *local_map_ptr = local_map;
      jx_TFree(core_marker);
      return jx_error_flag;
   }

   probe_sol = jx_ParVectorCreate(comm,
                                  jx_ParCSRMatrixGlobalNumRows(A),
                                  jx_ParCSRMatrixRowStarts(A));
   probe_res = jx_ParVectorCreate(comm,
                                  jx_ParCSRMatrixGlobalNumRows(A),
                                  jx_ParCSRMatrixRowStarts(A));
   probe_work = jx_ParVectorCreate(comm,
                                   jx_ParCSRMatrixGlobalNumRows(A),
                                   jx_ParCSRMatrixRowStarts(A));
   jx_ParVectorInitialize(probe_sol);
   jx_ParVectorInitialize(probe_res);
   jx_ParVectorInitialize(probe_work);
   jx_ParVectorSetConstantValues(probe_sol, 0.0);
   jx_ParVectorSetConstantValues(probe_res, 0.0);
   jx_ParVectorSetConstantValues(probe_work, 0.0);

   for (probe_iter = 0; probe_iter < base_max_iter; probe_iter++)
   {
      jx_ILUSolveLU(A, probe_rhs, probe_sol, perm, nLU, L, D, U,
                    probe_res, probe_work);
   }
   jx_ParCSRMatrixMatvecOutOfPlace(-1.0, A, probe_sol, 1.0,
                                   probe_rhs, probe_res);

   res_local = jx_ParVectorLocalVector(probe_res);
   res_data = jx_VectorData(res_local);
   scores = jx_TAlloc(jx_ILULRCScore, n);
   for (i = 0; i < n; i++)
   {
      scores[i].index = i;
      scores[i].score = jx_abs(res_data[i]);
   }
   qsort(scores, (size_t)n, sizeof(jx_ILULRCScore), jx_ILULRCScoreCompare);

   local_total_energy = 0.0;
   for (i = 0; i < n; i++)
   {
      local_total_energy += res_data[i] * res_data[i];
   }
   jx_MPI_Allreduce(&local_total_energy, &global_total_energy, 1,
                    JX_MPI_REAL, MPI_SUM, comm);

   if (select_type == 1)
   {
      target = 0;
      running_energy = 0.0;
      if (local_total_energy > 0.0)
      {
         while (target < n && running_energy / local_total_energy < energy_target)
         {
            old_row = scores[target].index;
            running_energy += res_data[old_row] * res_data[old_row];
            target++;
         }
      }
      if (energy_target > 0.0 && target < 1 && n > 0)
      {
         target = 1;
      }
   }
   else
   {
      target = (JX_Int)(ratio * (JX_Real)n + 0.5);
      if (ratio > 0.0 && target < 1 && n > 0)
      {
         target = 1;
      }
   }
   if (target > n)
   {
      target = n;
   }
   if (target == 0)
   {
      *marker_ptr = marker;
      *local_map_ptr = local_map;
      jx_TFree(scores);
      jx_TFree(core_marker);
      jx_ParVectorDestroy(probe_sol);
      jx_ParVectorDestroy(probe_res);
      jx_ParVectorDestroy(probe_work);
      return jx_error_flag;
   }

   curve_marker = jx_CTAlloc(JX_Int, n);
   for (curve_id = 0; curve_id < 5; curve_id++)
   {
      curve_target = (JX_Int)(curve_ratios[curve_id] * (JX_Real)n + 0.5);
      if (curve_target < 1 && n > 0)
      {
         curve_target = 1;
      }
      if (curve_target > n)
      {
         curve_target = n;
      }
      for (i = 0; i < n; i++)
      {
         curve_marker[i] = 0;
      }
      for (i = 0; i < curve_target; i++)
      {
         old_row = scores[i].index;
         local_curve_core[curve_id] += res_data[old_row] * res_data[old_row];
         curve_marker[old_row] = 1;
      }
      for (i = 0; i < curve_target; i++)
      {
         old_row = scores[i].index;
         for (j = A_i[old_row]; j < A_i[old_row + 1]; j++)
         {
            new_col = A_j[j];
            if (new_col >= 0 && new_col < n)
            {
               curve_marker[new_col] = 1;
            }
         }
      }
      for (i = 0; i < n; i++)
      {
         if (curve_marker[i])
         {
            local_curve_expanded[curve_id] += res_data[i] * res_data[i];
         }
      }
   }
   jx_MPI_Allreduce(local_curve_core, global_curve_core, 5,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_MPI_Allreduce(local_curve_expanded, global_curve_expanded, 5,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_TFree(curve_marker);

   core_indices = jx_TAlloc(JX_Int, target);
   for (i = 0; i < target; i++)
   {
      old_row = scores[i].index;
      core_marker[old_row] = 1;
      marker[old_row] = 1;
      core_indices[i] = old_row;
   }
   qsort(core_indices, (size_t)target, sizeof(JX_Int), jx_ILULRCIntCompare);

   for (i = 0; i < target; i++)
   {
      old_row = core_indices[i];
      for (j = A_i[old_row]; j < A_i[old_row + 1]; j++)
      {
         new_col = A_j[j];
         if (new_col >= 0 && new_col < n)
         {
            marker[new_col] = 1;
         }
      }
   }

   local_core_energy = 0.0;
   local_expanded_energy = 0.0;
   for (i = 0; i < n; i++)
   {
      JX_Real e = res_data[i] * res_data[i];
      if (core_marker[i])
      {
         local_core_energy += e;
      }
      if (marker[i])
      {
         local_expanded_energy += e;
      }
   }
   jx_MPI_Allreduce(&local_core_energy, &global_core_energy, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_MPI_Allreduce(&local_expanded_energy, &global_expanded_energy, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   if (global_total_energy > 0.0)
   {
      *probe_core_energy_ratio_ptr = global_core_energy / global_total_energy;
      *probe_expanded_energy_ratio_ptr = global_expanded_energy / global_total_energy;
   }

   expanded = 0;
   for (i = 0; i < n; i++)
   {
      if (marker[i])
      {
         expanded++;
      }
   }

   local_indices = jx_TAlloc(JX_Int, expanded);
   k = 0;
   for (i = 0; i < n; i++)
   {
      if (marker[i])
      {
         local_map[i] = k;
         local_indices[k++] = i;
      }
   }

   local_diag_abs = 0.0;
   local_offd_abs = 0.0;
   for (i = 0; i < n; i++)
   {
      if (marker[i])
      {
         for (j = A_i[i]; j < A_i[i + 1]; j++)
         {
            local_diag_abs += jx_abs(A_data[j]);
         }
         for (j = A_offd_i[i]; j < A_offd_i[i + 1]; j++)
         {
            local_offd_abs += jx_abs(A_offd_data[j]);
         }
      }
   }
   jx_MPI_Allreduce(&local_diag_abs, &global_diag_abs, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_MPI_Allreduce(&local_offd_abs, &global_offd_abs, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   if (global_diag_abs > 0.0)
   {
      *offd_coupling_ratio_ptr = global_offd_abs / global_diag_abs;
   }

   local_nnz = 0;
   for (i = 0; i < expanded; i++)
   {
      old_row = local_indices[i];
      local_nnz++;
      for (j = A_i[old_row]; j < A_i[old_row + 1]; j++)
      {
         new_col = local_map[A_j[j]];
         if (new_col >= 0 && A_j[j] != old_row)
         {
            local_nnz++;
         }
      }
   }

   local_mat = jx_CSRMatrixCreate(expanded, expanded, local_nnz);
   jx_CSRMatrixInitialize(local_mat);
   B_i = jx_CSRMatrixI(local_mat);
   B_j = jx_CSRMatrixJ(local_mat);
   B_data = jx_CSRMatrixData(local_mat);

   k = 0;
   B_i[0] = 0;
   for (i = 0; i < expanded; i++)
   {
      old_row = local_indices[i];
      diag_value = 0.0;
      diag_pos = -1;
      for (j = A_i[old_row]; j < A_i[old_row + 1]; j++)
      {
         if (A_j[j] == old_row)
         {
            diag_pos = j;
            diag_value = A_data[j];
            break;
         }
      }
      if (diag_pos < 0 || jx_abs(diag_value) < MAT_TOL)
      {
         diag_value = 1.0e-06;
      }
      B_j[k] = i;
      B_data[k++] = diag_value;

      for (j = A_i[old_row]; j < A_i[old_row + 1]; j++)
      {
         new_col = local_map[A_j[j]];
         if (new_col >= 0 && A_j[j] != old_row)
         {
            B_j[k] = new_col;
            B_data[k++] = A_data[j];
         }
      }
      B_i[i + 1] = k;
   }

   jx_ILUZeroDecompositionA(local_mat, &indexD, &indexLU, &valueLU);
   local_bad_pivots = 0;
   local_min_pivot = 1.0e100;
   for (i = 0; i < expanded; i++)
   {
      inv_diag_abs = jx_abs(valueLU[indexD[i]]);
      if (inv_diag_abs > 0.0 && inv_diag_abs == inv_diag_abs)
      {
         pivot_abs = 1.0 / inv_diag_abs;
      }
      else
      {
         pivot_abs = 0.0;
      }
      if (pivot_abs < local_min_pivot)
      {
         local_min_pivot = pivot_abs;
      }
      if (pivot_abs < 1.0e-10)
      {
         local_bad_pivots++;
      }
   }
   local_rhs = jx_CTAlloc(JX_Real, expanded);
   local_sol = jx_CTAlloc(JX_Real, expanded);
   local_work = jx_CTAlloc(JX_Real, expanded);

   jx_MPI_Allreduce(&target, &global_core_int, 1, JX_MPI_INT, MPI_SUM, comm);
   jx_MPI_Allreduce(&expanded, &global_expanded_int, 1, JX_MPI_INT, MPI_SUM, comm);
   local_nnz_real = (JX_Real)local_nnz;
   local_ilu_nnz_real = (JX_Real)jx_CSRMatrixNumNonzeros(local_mat);
   jx_MPI_Allreduce(&local_nnz_real, &global_local_mat_nnz, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_MPI_Allreduce(&local_ilu_nnz_real, &global_local_ilu_nnz, 1,
                    JX_MPI_REAL, MPI_SUM, comm);
   jx_MPI_Allreduce(&local_bad_pivots, &global_bad_pivots, 1,
                    JX_MPI_INT, MPI_SUM, comm);
   jx_MPI_Allreduce(&local_min_pivot, &global_min_pivot, 1,
                    JX_MPI_REAL, MPI_MIN, comm);
   jx_MPI_Comm_rank(comm, &my_id);
   if ((my_id == 0) && (print_level > 0))
   {
      jx_printf("ILU LRC Top-ratio local ILU(0)+one-ring: base sweeps = %d, ratio = %e, select_type = %d, energy target = %e, core rows = %d, expanded rows = %d, expansion ratio = %.3f\n",
                base_max_iter, ratio, select_type, energy_target,
                global_core_int, global_expanded_int,
                global_core_int > 0 ? (JX_Real)global_expanded_int / (JX_Real)global_core_int : 0.0);
      jx_printf("ILU LRC local storage: A_SeSe nnz = %.0f, local ILU nnz = %.0f\n",
                global_local_mat_nnz, global_local_ilu_nnz);
      jx_printf("ILU LRC diagnostic: probe core energy = %.6e, expanded energy = %.6e, offd/diag coupling = %.6e\n",
                *probe_core_energy_ratio_ptr, *probe_expanded_energy_ratio_ptr,
                *offd_coupling_ratio_ptr);
      if (global_total_energy > 0.0)
      {
         jx_printf("ILU LRC energy curve: top05 core = %.6e expanded = %.6e\n",
                   global_curve_core[0] / global_total_energy,
                   global_curve_expanded[0] / global_total_energy);
         jx_printf("ILU LRC energy curve: top10 core = %.6e expanded = %.6e\n",
                   global_curve_core[1] / global_total_energy,
                   global_curve_expanded[1] / global_total_energy);
         jx_printf("ILU LRC energy curve: top20 core = %.6e expanded = %.6e\n",
                   global_curve_core[2] / global_total_energy,
                   global_curve_expanded[2] / global_total_energy);
         jx_printf("ILU LRC energy curve: top30 core = %.6e expanded = %.6e\n",
                   global_curve_core[3] / global_total_energy,
                   global_curve_expanded[3] / global_total_energy);
         jx_printf("ILU LRC energy curve: top40 core = %.6e expanded = %.6e\n",
                   global_curve_core[4] / global_total_energy,
                   global_curve_expanded[4] / global_total_energy);
      }
      if (global_bad_pivots > 0)
      {
         jx_printf("ILU LRC WARNING: local ILU has %d pivots with |pivot| < 1e-10, min |pivot| = %.3e\n",
                   global_bad_pivots, global_min_pivot);
      }
   }

   *marker_ptr = marker;
   *core_indices_ptr = core_indices;
   *core_size_ptr = target;
   *local_indices_ptr = local_indices;
   *local_map_ptr = local_map;
   *local_mat_ptr = local_mat;
   *local_indexD_ptr = indexD;
   *local_indexLU_ptr = indexLU;
   *local_valueLU_ptr = valueLU;
   *local_rhs_ptr = local_rhs;
   *local_sol_ptr = local_sol;
   *local_work_ptr = local_work;
   *num_selected_ptr = global_core_int;
   *global_core_size_ptr = (JX_Real)global_core_int;
   *global_expanded_size_ptr = (JX_Real)global_expanded_int;
   *global_local_mat_nnz_ptr = global_local_mat_nnz;
   *global_local_ilu_nnz_ptr = global_local_ilu_nnz;

   jx_TFree(scores);
   jx_TFree(core_marker);
   jx_ParVectorDestroy(probe_sol);
   jx_ParVectorDestroy(probe_res);
   jx_ParVectorDestroy(probe_work);

   return jx_error_flag;
}

/* Setup ILU data */
JX_Int
jx_ILUSetup(void *ilu_vdata,
            jx_ParCSRMatrix *A,
            jx_ParVector *f,
            jx_ParVector *u)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;
   // jx_ParILUData *schur_precond_ilu;
   // jx_ParNSHData *schur_solver_nsh;

   /* Pointers to ilu data */
   JX_Int logging = jx_ParILUDataLogging(ilu_data);
   JX_Int print_level = jx_ParILUDataPrintLevel(ilu_data);
   JX_Int ilu_type = jx_ParILUDataIluType(ilu_data);
   JX_Int nLU = jx_ParILUDataNLU(ilu_data);
   JX_Int nI = jx_ParILUDataNI(ilu_data);
   JX_Int fill_level = jx_ParILUDataLfil(ilu_data);
   JX_Int max_row_elmts = jx_ParILUDataMaxRowNnz(ilu_data);
   JX_Int max_iter = jx_ParILUDataMaxIter(ilu_data);
   JX_Int lrc_type = jx_ParILUDataLocalResidualCorrectionType(ilu_data);
   JX_Real lrc_threshold = jx_ParILUDataLocalResidualCorrectionThreshold(ilu_data);
   JX_Int lrc_select_type = jx_ParILUDataLocalResidualCorrectionSelectType(ilu_data);
   JX_Real lrc_energy_target = jx_ParILUDataLocalResidualCorrectionEnergyTarget(ilu_data);
   JX_Real *droptol = jx_ParILUDataDroptol(ilu_data);
   JX_Int *CF_marker_array = jx_ParILUDataCFMarkerArray(ilu_data);
   JX_Int *perm = jx_ParILUDataPerm(ilu_data);
   JX_Int *qperm = jx_ParILUDataQPerm(ilu_data);
   JX_Real tol_ddPQ = jx_ParILUDataTolDDPQ(ilu_data);

   jx_ParCSRMatrix *matA = jx_ParILUDataMatA(ilu_data);
   jx_ParCSRMatrix *matL = jx_ParILUDataMatL(ilu_data);
   JX_Real *matD = jx_ParILUDataMatD(ilu_data);
   jx_ParCSRMatrix *matU = jx_ParILUDataMatU(ilu_data);
   jx_ParCSRMatrix *matmL = jx_ParILUDataMatLModified(ilu_data);
   JX_Real *matmD = jx_ParILUDataMatDModified(ilu_data);
   jx_ParCSRMatrix *matmU = jx_ParILUDataMatUModified(ilu_data);
   jx_ParCSRMatrix *matS = jx_ParILUDataMatS(ilu_data);
   JX_Int n = jx_CSRMatrixNumRows(jx_ParCSRMatrixDiag(A));
   JX_Int reordering_type = jx_ParILUDataReorderingType(ilu_data);
   JX_Real nnzS; /* Total nnz in S */
   JX_Real nnzS_offd_local;
   JX_Real nnzS_offd;
   JX_Int size_C /* Total size of coarse grid */;

   jx_ParVector *Utemp = NULL;
   jx_ParVector *Ftemp = NULL;
   jx_ParVector *Xtemp = NULL;
   jx_ParVector *Ytemp = NULL;
   jx_ParVector *Ztemp = NULL;
   JX_Real *uext = NULL;
   JX_Real *fext = NULL;
   JX_Int *lrc_marker = NULL;
   JX_Int lrc_num_selected = 0;
   JX_Int lrc_core_size = 0;
   JX_Int *lrc_core_indices = NULL;
   JX_Int *lrc_local_indices = NULL;
   JX_Int *lrc_local_map = NULL;
   jx_CSRMatrix *lrc_local_mat = NULL;
   JX_Int *lrc_local_indexD = NULL;
   JX_Int *lrc_local_indexLU = NULL;
   JX_Real *lrc_local_valueLU = NULL;
   JX_Real *lrc_local_rhs = NULL;
   JX_Real *lrc_local_sol = NULL;
   JX_Real *lrc_local_work = NULL;
   jx_ParVector *rhs = NULL;
   jx_ParVector *x = NULL;

   /* TODO (VPM): Change F_array and U_array variable names */
   jx_ParVector *F_array = jx_ParILUDataF(ilu_data);
   jx_ParVector *U_array = jx_ParILUDataU(ilu_data);
   jx_ParVector *residual = jx_ParILUDataResidual(ilu_data);
   JX_Real *rel_res_norms = jx_ParILUDataRelResNorms(ilu_data);

   /* might need for Schur Complement */
   JX_Int *u_end = NULL;
   JX_Solver schur_solver = NULL;
   JX_Solver schur_precond = NULL;
   JX_Solver schur_precond_gotten = NULL;

   /* Whether or not to use exact (direct) triangular solves */
   JX_Int tri_solve = jx_ParILUDataTriSolve(ilu_data);

   /* help to build external */
   jx_ParCSRCommPkg *comm_pkg;
   JX_Int buffer_size;
   JX_Int num_sends;
   JX_Int send_size;
   JX_Int recv_size;
   JX_Int num_procs, my_id;

   /* ----- begin -----*/
   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);

   /* Free previously allocated data, if any not destroyed */
   jx_ParCSRMatrixDestroy(matL);
   matL = NULL;
   jx_ParCSRMatrixDestroy(matU);
   matU = NULL;
   jx_ParCSRMatrixDestroy(matmL);
   matmL = NULL;
   jx_ParCSRMatrixDestroy(matmU);
   matmU = NULL;
   jx_ParCSRMatrixDestroy(matS);
   matS = NULL;

   jx_TFree(matD);
   jx_TFree(matmD);
   jx_TFree(CF_marker_array);
   jx_TFree(jx_ParILUDataLocalResidualCorrectionMarker(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionCoreIndices(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalIndices(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalMap(ilu_data));
   jx_CSRMatrixDestroy(jx_ParILUDataLocalResidualCorrectionLocalMat(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalIndexD(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalIndexLU(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalValueLU(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalRhs(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalSol(ilu_data));
   jx_TFree(jx_ParILUDataLocalResidualCorrectionLocalWork(ilu_data));
   jx_ParILUDataLocalResidualCorrectionMarker(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionNumSelected(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionCoreSize(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionCoreIndices(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalIndices(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalMap(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalMat(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalIndexD(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalIndexLU(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalValueLU(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalRhs(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalSol(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionLocalWork(ilu_data) = NULL;
   jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionGlobalExpandedSize(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionGlobalLocalMatNnz(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionGlobalLocalIluNnz(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionLocalMatComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionLocalIluComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionExtraComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionSetupWorkComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionEffectiveComplexity(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionApplyCount(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionApplyTime(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionMatvecTime(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionLocalSolveTime(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionProbeCoreEnergyRatio(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionProbeExpandedEnergyRatio(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionOffdCouplingRatio(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionDeltaNormRatioSum(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionResidualRatioSum(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionResidualRatioMin(ilu_data) = 1.0e100;
   jx_ParILUDataLocalResidualCorrectionResidualRatioMax(ilu_data) = 0.0;
   jx_ParILUDataLocalResidualCorrectionAcceptLikeCount(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionRejectLikeCount(ilu_data) = 0;

   /* clear old l1_norm data, if created */
   jx_TFree(jx_ParILUDataL1Norms(ilu_data));

   /* setup temporary storage
    * first check is they've already here
    */
   jx_ParVectorDestroy(jx_ParILUDataUTemp(ilu_data));
   jx_ParVectorDestroy(jx_ParILUDataFTemp(ilu_data));
   jx_ParVectorDestroy(jx_ParILUDataRhs(ilu_data));
   jx_ParVectorDestroy(jx_ParILUDataX(ilu_data));
   jx_ParVectorDestroy(jx_ParILUDataResidual(ilu_data));
   jx_TFree(jx_ParILUDataUExt(ilu_data));
   jx_TFree(jx_ParILUDataFExt(ilu_data));
   jx_TFree(jx_ParILUDataUEnd(ilu_data));
   jx_TFree(jx_ParILUDataRelResNorms(ilu_data));

   jx_ParILUDataUTemp(ilu_data) = NULL;
   jx_ParILUDataFTemp(ilu_data) = NULL;
   jx_ParILUDataRhs(ilu_data) = NULL;
   jx_ParILUDataX(ilu_data) = NULL;
   jx_ParILUDataResidual(ilu_data) = NULL;

   /* Current ILU(0)-only version does not create Schur solver/preconditioner. */
   jx_ParILUDataSchurSolver(ilu_data) = NULL;
   jx_ParILUDataSchurPrecond(ilu_data) = NULL;

   /* Create work vectors */
   Utemp = jx_ParVectorCreate(jx_ParCSRMatrixComm(A),
                              jx_ParCSRMatrixGlobalNumRows(A),
                              jx_ParCSRMatrixRowStarts(A));
   jx_ParVectorInitialize(Utemp);
   jx_ParILUDataUTemp(ilu_data) = Utemp;

   Ftemp = jx_ParVectorCreate(jx_ParCSRMatrixComm(A),
                              jx_ParCSRMatrixGlobalNumRows(A),
                              jx_ParCSRMatrixRowStarts(A));
   jx_ParVectorInitialize(Ftemp);
   jx_ParILUDataFTemp(ilu_data) = Ftemp;

   /* set matrix, solution and rhs pointers */
   matA = A;
   F_array = f;
   U_array = u;

   /* Create perm array if necessary */
   if (!perm)
   {
      jx_ILUGetLocalPerm(matA, &perm, &nLU, reordering_type);
   }

   /* Factorization */
   switch (ilu_type)
   {
   case 0: /* BJ + ILU(0) */
      nLU = n;
      nI = n;

      jx_ILUSetupILU0(matA, perm, perm, nLU, nI,
                      &matL, &matD, &matU, &matS, &u_end);
      break;

   case 1: /* BJ + ILU(k) */
      nLU = n;
      nI = n;

      if (fill_level < 0)
      {
         fill_level = 0;
      }

      jx_ILUSetupILUK(matA, fill_level, perm, perm, nLU, nI,
                      &matL, &matD, &matU, &matS, &u_end);
      break;

   case 2: /* BJ + ILUT */
      nLU = n;
      nI = n;

      jx_ILUSetupILUT(matA, max_row_elmts, droptol, perm, perm, nLU, nI,
                      &matL, &matD, &matU, &matS, &u_end);
      break;

   default:
      if (my_id == 0)
      {
         jx_printf("[JX-ILU] ERROR: unsupported ilu_type = %d. "
                   "Current simplified version only supports "
                   "0: ILU(0), 1: ILU(k), 2: ILUT.\n",
                   ilu_type);
      }

      jx_error_w_msg(JX_ERROR_ARG,
                     "Unsupported ILU type in current simplified ILU setup.\n");

      return jx_error_flag;
   }

   /* Create additional temporary vector for iterative triangular solve */
   if (!tri_solve)
   {
      Ztemp = jx_ParVectorCreate(jx_ParCSRMatrixComm(A),
                                 jx_ParCSRMatrixGlobalNumRows(A),
                                 jx_ParCSRMatrixRowStarts(A));
      jx_ParVectorInitialize(Ztemp);
   }

   /* setup Schur solver - TODO (VPM): merge host and device paths below */
   /* No Schur solver is needed in the current host-side BJ ILU(0) version. */

   /* Set pointers to ilu data */
   jx_ParILUDataMatA(ilu_data) = matA;
   jx_ParILUDataXTemp(ilu_data) = Xtemp;
   jx_ParILUDataYTemp(ilu_data) = Ytemp;
   jx_ParILUDataZTemp(ilu_data) = Ztemp;
   jx_ParILUDataF(ilu_data) = F_array;
   jx_ParILUDataU(ilu_data) = U_array;
   jx_ParILUDataMatL(ilu_data) = matL;
   jx_ParILUDataMatD(ilu_data) = matD;
   jx_ParILUDataMatU(ilu_data) = matU;
   jx_ParILUDataMatLModified(ilu_data) = matmL;
   jx_ParILUDataMatDModified(ilu_data) = matmD;
   jx_ParILUDataMatUModified(ilu_data) = matmU;
   jx_ParILUDataMatS(ilu_data) = matS;
   jx_ParILUDataCFMarkerArray(ilu_data) = CF_marker_array;
   jx_ParILUDataPerm(ilu_data) = perm;
   jx_ParILUDataQPerm(ilu_data) = qperm;
   jx_ParILUDataNLU(ilu_data) = nLU;
   jx_ParILUDataNI(ilu_data) = nI;
   jx_ParILUDataUEnd(ilu_data) = u_end;
   jx_ParILUDataUExt(ilu_data) = uext;
   jx_ParILUDataFExt(ilu_data) = fext;

   if (lrc_type == 2)
   {
      jx_ILUBuildLocalResidualCorrectionMask(matA, perm, nLU,
                                             matL, matD, matU,
                                             lrc_threshold, print_level,
                                             &lrc_marker, &lrc_num_selected);

      jx_ParILUDataLocalResidualCorrectionMarker(ilu_data) = lrc_marker;
      jx_ParILUDataLocalResidualCorrectionNumSelected(ilu_data) = lrc_num_selected;
   }
   else if (lrc_type == 3 || lrc_type == 4)
   {
      jx_ILUBuildTopRatioLocalResidualCorrection(matA, F_array, perm, nLU,
                                                 matL, matD, matU,
                                                 max_iter,
                                                 lrc_threshold,
                                                 lrc_select_type,
                                                 lrc_energy_target,
                                                 print_level,
                                                 &lrc_marker,
                                                 &lrc_core_indices,
                                                 &lrc_core_size,
                                                 &lrc_local_indices,
                                                 &lrc_local_map,
                                                 &lrc_local_mat,
                                                 &lrc_local_indexD,
                                                 &lrc_local_indexLU,
                                                 &lrc_local_valueLU,
                                                 &lrc_local_rhs,
                                                 &lrc_local_sol,
                                                 &lrc_local_work,
                                                 &lrc_num_selected,
                                                 &jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionGlobalExpandedSize(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionGlobalLocalMatNnz(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionGlobalLocalIluNnz(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionProbeCoreEnergyRatio(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionProbeExpandedEnergyRatio(ilu_data),
                                                 &jx_ParILUDataLocalResidualCorrectionOffdCouplingRatio(ilu_data));

      jx_ParILUDataLocalResidualCorrectionMarker(ilu_data) = lrc_marker;
      jx_ParILUDataLocalResidualCorrectionCoreSize(ilu_data) = lrc_core_size;
      jx_ParILUDataLocalResidualCorrectionCoreIndices(ilu_data) = lrc_core_indices;
      jx_ParILUDataLocalResidualCorrectionLocalIndices(ilu_data) = lrc_local_indices;
      jx_ParILUDataLocalResidualCorrectionLocalMap(ilu_data) = lrc_local_map;
      jx_ParILUDataLocalResidualCorrectionLocalMat(ilu_data) = lrc_local_mat;
      jx_ParILUDataLocalResidualCorrectionLocalIndexD(ilu_data) = lrc_local_indexD;
      jx_ParILUDataLocalResidualCorrectionLocalIndexLU(ilu_data) = lrc_local_indexLU;
      jx_ParILUDataLocalResidualCorrectionLocalValueLU(ilu_data) = lrc_local_valueLU;
      jx_ParILUDataLocalResidualCorrectionLocalRhs(ilu_data) = lrc_local_rhs;
      jx_ParILUDataLocalResidualCorrectionLocalSol(ilu_data) = lrc_local_sol;
      jx_ParILUDataLocalResidualCorrectionLocalWork(ilu_data) = lrc_local_work;
      jx_ParILUDataLocalResidualCorrectionNumSelected(ilu_data) = lrc_num_selected;
   }

   /* compute operator complexity */
   jx_ParCSRMatrixSetDNumNonzeros(matA);
   nnzS = 0.0;

   /* size_C is the size of global coarse grid, upper left part */
   size_C = jx_ParCSRMatrixGlobalNumRows(matA);

   /* switch to compute complexity */
   jx_ParILUDataOperatorComplexity(ilu_data) = ((JX_Real)size_C + nnzS +
                                                jx_ParCSRMatrixDNumNonzeros(matL) +
                                                jx_ParCSRMatrixDNumNonzeros(matU)) /
                                               jx_ParCSRMatrixDNumNonzeros(matA);

   if ((lrc_type == 3 || lrc_type == 4) && jx_ParCSRMatrixDNumNonzeros(matA) > 0)
   {
      JX_Real annz = (JX_Real)jx_ParCSRMatrixDNumNonzeros(matA);
      jx_ParILUDataLocalResidualCorrectionLocalMatComplexity(ilu_data) =
         jx_ParILUDataLocalResidualCorrectionGlobalLocalMatNnz(ilu_data) / annz;
      jx_ParILUDataLocalResidualCorrectionLocalIluComplexity(ilu_data) =
         jx_ParILUDataLocalResidualCorrectionGlobalLocalIluNnz(ilu_data) / annz;
      jx_ParILUDataLocalResidualCorrectionExtraComplexity(ilu_data) =
         jx_ParILUDataLocalResidualCorrectionLocalMatComplexity(ilu_data) +
         jx_ParILUDataLocalResidualCorrectionLocalIluComplexity(ilu_data);
      jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data) =
         jx_ParILUDataOperatorComplexity(ilu_data) +
         jx_ParILUDataLocalResidualCorrectionExtraComplexity(ilu_data);
      jx_ParILUDataLocalResidualCorrectionSetupWorkComplexity(ilu_data) =
         jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data);
      jx_ParILUDataLocalResidualCorrectionEffectiveComplexity(ilu_data) =
         jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data);
   }

   /* TODO (VPM): Move ILU statistics printout to its own function */
   if ((my_id == 0) && (print_level > 0))
   {
      jx_printf("ILU SETUP: operator complexity = %f  \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      if (lrc_type == 3 || lrc_type == 4)
      {
         jx_printf("ILU LRC SETUP: core = %.0f, expanded = %.0f, expansion ratio = %.3f\n",
                   jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionGlobalExpandedSize(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data) > 0.0 ?
                      jx_ParILUDataLocalResidualCorrectionGlobalExpandedSize(ilu_data) /
                         jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data) : 0.0);
         jx_printf("ILU LRC SETUP: A_SeSe storage complexity = %f, local ILU storage complexity = %f, extra storage complexity = %f\n",
                   jx_ParILUDataLocalResidualCorrectionLocalMatComplexity(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionLocalIluComplexity(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionExtraComplexity(ilu_data));
         jx_printf("ILU LRC SETUP: storage complexity = %f, setup work complexity = %f, total effective complexity = %f\n",
                   jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionSetupWorkComplexity(ilu_data),
                   jx_ParILUDataLocalResidualCorrectionEffectiveComplexity(ilu_data));
      }
      if (jx_ParILUDataTriSolve(ilu_data))
      {
         jx_printf("ILU SOLVE: using direct triangular solves\n");
      }
      else
      {
         jx_printf("ILU SOLVE: using iterative triangular solves\n");
      }
   }

   if (logging > 1)
   {
      residual = jx_ParVectorCreate(jx_ParCSRMatrixComm(matA),
                                    jx_ParCSRMatrixGlobalNumRows(matA),
                                    jx_ParCSRMatrixRowStarts(matA));
      jx_ParVectorInitialize(residual);
      jx_ParILUDataResidual(ilu_data) = residual;
   }
   else
   {
      jx_ParILUDataResidual(ilu_data) = NULL;
   }

   rel_res_norms = jx_CTAlloc(JX_Real, jx_ParILUDataMaxIter(ilu_data));
   jx_ParILUDataRelResNorms(ilu_data) = rel_res_norms;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILU0
 *
 * Setup ILU(0)
 *
 * A = input matrix
 * perm = permutation array indicating ordering of rows.
 *        Perm could come from a CF_marker array or a reordering routine.
 *         When set to NULL, identity permutation is used.
 * qperm = permutation array indicating ordering of columns.
 *         When set to NULL, identity permutation is used.
 * nI = number of interial unknowns
 * nLU = size of incomplete factorization, nLU should obey nLU <= nI.
 *       Schur complement is formed if nLU < n
 * Lptr, Dptr, Uptr, Sptr = L, D, U, S factors.
 * will form global Schur Matrix if nLU < n
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetupILU0(jx_ParCSRMatrix *A,
                JX_Int *perm,
                JX_Int *qperm,
                JX_Int nLU,
                JX_Int nI,
                jx_ParCSRMatrix **Lptr,
                JX_Real **Dptr,
                jx_ParCSRMatrix **Uptr,
                jx_ParCSRMatrix **Sptr,
                JX_Int **u_end)
{
   return jx_ILUSetupMILU0(A, perm, qperm, nLU, nI, Lptr, Dptr, Uptr, Sptr, u_end, 0);
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILU0
 *
 * Setup modified ILU(0)
 *
 * A = input matrix
 * perm = permutation array indicating ordering of rows.
 *        Perm could come from a CF_marker array or a reordering routine.
 *        When set to NULL, indentity permutation is used.
 * qperm = permutation array indicating ordering of columns.
 *         When set to NULL, identity permutation is used.
 * nI = number of interior unknowns
 * nLU = size of incomplete factorization, nLU should obey nLU <= nI.
 *       Schur complement is formed if nLU < n
 * Lptr, Dptr, Uptr, Sptr = L, D, U, S factors.
 * modified set to 0 to use classical ILU
 * will form global Schur Matrix if nLU < n
 *--------------------------------------------------------------------------*/
JX_Int
jx_ILUSetupMILU0(jx_ParCSRMatrix *A,
                 JX_Int *permp,
                 JX_Int *qpermp,
                 JX_Int nLU,
                 JX_Int nI,
                 jx_ParCSRMatrix **Lptr,
                 JX_Real **Dptr,
                 jx_ParCSRMatrix **Uptr,
                 jx_ParCSRMatrix **Sptr,
                 JX_Int **u_end,
                 JX_Int modified)
{
   JX_Int i, ii, j, k, k1, k2, k3, ctrU, ctrL, ctrS;
   JX_Int lenl, lenu, jpiv, col, jpos;
   JX_Int *iw, *iL, *iU;
   JX_Real dd, t, dpiv, lxu, *wU, *wL;
   JX_Real drop;

   /* communication stuffs for S */
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   JX_Int S_offd_nnz, S_offd_ncols;
   jx_ParCSRCommPkg *comm_pkg;
   jx_ParCSRCommHandle *comm_handle;
   JX_Int num_sends, begin, end;
   JX_Int *send_buf = NULL;
   JX_Int num_procs, my_id;

   /* data objects for A */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);

   /* size of problem and schur system */
   JX_Int n = jx_CSRMatrixNumRows(A_diag);
   JX_Int m = n - nLU;
   JX_Int e = nI - nLU;
   JX_Int m_e = n - nI;
   JX_Real local_nnz, total_nnz;
   JX_Int *u_end_array;

   /* data objects for L, D, U */
   jx_ParCSRMatrix *matL;
   jx_ParCSRMatrix *matU;
   jx_CSRMatrix *L_diag;
   jx_CSRMatrix *U_diag;
   JX_Real *D_data;
   JX_Real *L_diag_data;
   JX_Int *L_diag_i;
   JX_Int *L_diag_j;
   JX_Real *U_diag_data;
   JX_Int *U_diag_i;
   JX_Int *U_diag_j;

   /* data objects for S */
   jx_ParCSRMatrix *matS = NULL;
   jx_CSRMatrix *S_diag;
   jx_CSRMatrix *S_offd;
   JX_Real *S_diag_data = NULL;
   JX_Int *S_diag_i = NULL;
   JX_Int *S_diag_j = NULL;
   JX_Int *S_offd_i = NULL;
   JX_Int *S_offd_j = NULL;
   JX_Int *S_offd_colmap = NULL;
   JX_Real *S_offd_data;
   JX_Int col_starts[2];
   JX_Int total_rows;

   /* memory management */
   JX_Int initial_alloc = 0;
   JX_Int capacity_L;
   JX_Int capacity_U;
   JX_Int capacity_S = 0;
   JX_Int nnz_A = A_diag_i[n];

   /* reverse permutation array */
   JX_Int *rperm;
   JX_Int *perm, *qperm;

   /* start setup
    * get communication stuffs first
    */
   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);
   comm_pkg = jx_ParCSRMatrixCommPkg(A);

   /* setup if not yet built */
   if (!comm_pkg)
   {
      jx_MatvecCommPkgCreate(A);
      comm_pkg = jx_ParCSRMatrixCommPkg(A);
   }

   /* check for correctness */
   if (nLU < 0 || nLU > n)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU out of range.\n");
   }
   if (e < 0)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU should not exceed nI.\n");
   }

   /* Allocate memory for u_end array */
   u_end_array = jx_TAlloc(JX_Int, nLU);

   /* Allocate memory for L,D,U,S factors */
   if (n > 0)
   {
      initial_alloc = (JX_Int)(nLU + ceil((nnz_A / 2.0) * nLU / n));
      capacity_S = (JX_Int)(m + ceil((nnz_A / 2.0) * m / n));
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   D_data = jx_TAlloc(JX_Real, n);
   L_diag_i = jx_TAlloc(JX_Int, n + 1);
   L_diag_j = jx_TAlloc(JX_Int, capacity_L);
   L_diag_data = jx_TAlloc(JX_Real, capacity_L);
   U_diag_i = jx_TAlloc(JX_Int, n + 1);
   U_diag_j = jx_TAlloc(JX_Int, capacity_U);
   U_diag_data = jx_TAlloc(JX_Real, capacity_U);
   S_diag_i = jx_TAlloc(JX_Int, m + 1);
   S_diag_j = jx_TAlloc(JX_Int, capacity_S);
   S_diag_data = jx_TAlloc(JX_Real, capacity_S);

   /* allocate working arrays */
   iw = jx_TAlloc(JX_Int, 3 * n);
   iL = iw + n;
   rperm = iw + 2 * n;
   wL = jx_TAlloc(JX_Real, n);

   ctrU = ctrL = ctrS = 0;
   L_diag_i[0] = U_diag_i[0] = S_diag_i[0] = 0;
   /* set marker array iw to -1 */
   for (i = 0; i < n; i++)
   {
      iw[i] = -1;
   }

   /* get reverse permutation (rperm).
    * create permutation if they are null
    * rperm holds the reordered indexes.
    * rperm only used for column
    */

   if (!permp)
   {
      perm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         perm[i] = i;
      }
   }
   else
   {
      perm = permp;
   }

   if (!qpermp)
   {
      qperm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         qperm[i] = i;
      }
   }
   else
   {
      qperm = qpermp;
   }

   for (i = 0; i < n; i++)
   {
      rperm[qperm[i]] = i;
   }

   /*---------  Begin Factorization. Work in permuted space  ----*/
   for (ii = 0; ii < nLU; ii++)
   {
      // get row i
      i = perm[ii];
      // get extents of row i
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      // track the drop
      drop = 0.0;

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL + ii;
      wU = wL + ii;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl = lenu = 0;
      iw[ii] = ii;
      /*-------------------- scan & unwrap column */
      for (j = k1; j < k2; j++)
      {
         col = rperm[A_diag_j[j]];
         t = A_diag_data[j];
         if (col < ii)
         {
            iw[col] = lenl;
            iL[lenl] = col;
            wL[lenl++] = t;
         }
         else if (col > ii)
         {
            iw[col] = lenu;
            iU[lenu] = col;
            wU[lenu++] = t;
         }
         else
         {
            dd = t;
         }
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jx_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jx_qsort3ir(iL, wL, iw, 0, (lenl - 1));
      for (j = 0; j < lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for (k = U_diag_i[jpiv]; k < U_diag_i[jpiv + 1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if (jpos < 0)
            {
               drop = drop - U_diag_data[k] * dpiv;
               continue;
            }

            lxu = -U_diag_data[k] * dpiv;
            if (col < ii)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if (col > ii)
            {
               /* dealing with U part */
               wU[jpos] += lxu;
            }
            else
            {
               /* diagonal update */
               dd += lxu;
            }
         }
      }
      /* modify when necessary */
      if (modified)
      {
         dd = dd + drop;
      }

      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for (j = 0; j < lenu; j++)
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      if (lenl > 0)
      {
         while ((ctrL + lenl) > capacity_L)
         {
            JX_Int tmp = capacity_L;
            capacity_L = (JX_Int)(capacity_L * EXPAND_FACT + 1);
            L_diag_j = jx_TReAlloc_v2(L_diag_j, JX_Int, tmp, JX_Int,
                                      capacity_L);
            L_diag_data = jx_TReAlloc_v2(L_diag_data, JX_Real, tmp, JX_Real,
                                         capacity_L);
         }
         jx_TMemcpy(&L_diag_j[ctrL], iL, JX_Int, lenl);
         jx_TMemcpy(&L_diag_data[ctrL], wL, JX_Real, lenl);
      }
      L_diag_i[ii + 1] = (ctrL += lenl);

      /* diagonal part (we store the inverse) */
      if (jx_abs(dd) < MAT_TOL)
      {
         dd = 1.0e-6;
      }
      D_data[ii] = 1. / dd;

      /* U part */
      /* Check that memory is sufficient */
      if (lenu > 0)
      {
         while ((ctrU + lenu) > capacity_U)
         {
            JX_Int tmp = capacity_U;
            capacity_U = (JX_Int)(capacity_U * EXPAND_FACT + 1);
            U_diag_j = jx_TReAlloc_v2(U_diag_j, JX_Int, tmp, JX_Int,
                                      capacity_U);
            U_diag_data = jx_TReAlloc_v2(U_diag_data, JX_Real, tmp, JX_Real,
                                         capacity_U);
         }
         jx_TMemcpy(&U_diag_j[ctrU], iU, JX_Int, lenu);
         jx_TMemcpy(&U_diag_data[ctrU], wU, JX_Real, lenu);
      }
      U_diag_i[ii + 1] = (ctrU += lenu);

      /* check and build u_end array */
      if (m > 0)
      {
         jx_qsort1(U_diag_j, U_diag_data, U_diag_i[ii], U_diag_i[ii + 1] - 1);
         jx_BinarySearch2(U_diag_j, nLU, U_diag_i[ii], U_diag_i[ii + 1] - 1, u_end_array + ii);
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }
   }

   // /*---------  Begin Factorization in Schur Complement part  ----*/
   // for (ii = nLU; ii < n; ii++)
   // {
   //    // get row i
   //    i = perm[ii];
   //    // get extents of row i
   //    k1 = A_diag_i[i];
   //    k2 = A_diag_i[i + 1];
   //    drop = 0.0;

   //    /*-------------------- unpack L & U-parts of row of A in arrays w */
   //    iU = iL + nLU + 1;
   //    wU = wL + nLU + 1;
   //    /*--------------------  diagonal entry */
   //    dd = 0.0;
   //    lenl = lenu = 0;
   //    iw[ii] = nLU;
   //    /*-------------------- scan & unwrap column */
   //    for (j = k1; j < k2; j++)
   //    {
   //       col = rperm[A_diag_j[j]];
   //       t = A_diag_data[j];
   //       if (col < nLU)
   //       {
   //          iw[col] = lenl;
   //          iL[lenl] = col;
   //          wL[lenl++] = t;
   //       }
   //       else if (col != ii)
   //       {
   //          iw[col] = lenu;
   //          iU[lenu] = col;
   //          wU[lenu++] = t;
   //       }
   //       else
   //       {
   //          dd = t;
   //       }
   //    }

   //    /* eliminate row */
   //    /*-------------------------------------------------------------------------
   //     *  In order to do the elimination in the correct order we must select the
   //     *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
   //     *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
   //     *  entering the elimination loop.
   //     *-----------------------------------------------------------------------*/
   //    //      jx_quickSortIR(iL, wL, iw, 0, (lenl-1));
   //    jx_qsort3ir(iL, wL, iw, 0, (lenl - 1));
   //    for (j = 0; j < lenl; j++)
   //    {
   //       jpiv = iL[j];
   //       /* get factor/ pivot element */
   //       dpiv = wL[j] * D_data[jpiv];
   //       /* store entry in L */
   //       wL[j] = dpiv;

   //       /* zero out element - reset pivot */
   //       iw[jpiv] = -1;
   //       /* combine current row and pivot row */
   //       for (k = U_diag_i[jpiv]; k < U_diag_i[jpiv + 1]; k++)
   //       {
   //          col = U_diag_j[k];
   //          jpos = iw[col];

   //          /* Only fill-in nonzero pattern (jpos != 0) */
   //          if (jpos < 0)
   //          {
   //             drop = drop - U_diag_data[k] * dpiv;
   //             continue;
   //          }

   //          lxu = -U_diag_data[k] * dpiv;
   //          if (col < nLU)
   //          {
   //             /* dealing with L part */
   //             wL[jpos] += lxu;
   //          }
   //          else if (col != ii)
   //          {
   //             /* dealing with U part */
   //             wU[jpos] += lxu;
   //          }
   //          else
   //          {
   //             /* diagonal update */
   //             dd += lxu;
   //          }
   //       }
   //    }
   //    if (modified)
   //    {
   //       dd = dd + drop;
   //    }
   //    /* restore iw (only need to restore diagonal and U part */
   //    iw[ii] = -1;
   //    for (j = 0; j < lenu; j++)
   //    {
   //       iw[iU[j]] = -1;
   //    }

   //    /* Update LDU factors */
   //    /* L part */
   //    /* Check that memory is sufficient */
   //    if (lenl > 0)
   //    {
   //       while ((ctrL + lenl) > capacity_L)
   //       {
   //          JX_Int tmp = capacity_L;
   //          capacity_L = (JX_Int)(capacity_L * EXPAND_FACT + 1);
   //          L_diag_j = jx_TReAlloc_v2(L_diag_j, JX_Int, tmp, JX_Int,
   //                                    capacity_L);
   //          L_diag_data = jx_TReAlloc_v2(L_diag_data, JX_Real, tmp, JX_Real,
   //                                       capacity_L);
   //       }
   //       jx_TMemcpy(&L_diag_j[ctrL], iL, JX_Int, lenl);
   //       jx_TMemcpy(&L_diag_data[ctrL], wL, JX_Real, lenl);
   //    }
   //    L_diag_i[ii + 1] = (ctrL += lenl);

   //    /* S part */
   //    /* Check that memory is sufficient */
   //    while ((ctrS + lenu + 1) > capacity_S)
   //    {
   //       JX_Int tmp = capacity_S;
   //       capacity_S = (JX_Int)(capacity_S * EXPAND_FACT + 1);
   //       S_diag_j = jx_TReAlloc_v2(S_diag_j, JX_Int, tmp, JX_Int,
   //                                 capacity_S);
   //       S_diag_data = jx_TReAlloc_v2(S_diag_data, JX_Real, tmp, JX_Real,
   //                                    capacity_S);
   //    }
   //    /* remember S in under a new index system! */
   //    S_diag_j[ctrS] = ii - nLU;
   //    S_diag_data[ctrS] = dd;
   //    for (j = 0; j < lenu; j++)
   //    {
   //       S_diag_j[ctrS + 1 + j] = iU[j] - nLU;
   //    }
   //    // jx_TMemcpy(S_diag_data+ctrS+1, wU, JX_Real, lenu);
   //    jx_TMemcpy(S_diag_data + ctrS + 1, wU, JX_Real, lenu);
   //    S_diag_i[ii - nLU + 1] = ctrS += (lenu + 1);
   // }

   /* Assemble LDUS matrices */
   /* zero out unfactored rows for U and D */
   for (k = nLU; k < n; k++)
   {
      U_diag_i[k + 1] = ctrU;
      D_data[k] = 1.;
   }

   // /* First create Schur complement if necessary
   //  * Check if we need to create Schur complement
   //  */
   // JX_Int big_m = (JX_Int)m;
   // jx_MPI_Allreduce(&big_m, &total_rows, 1, JX_MPI_BIG_INT, MPI_SUM, comm);

   // /* only form when total_rows > 0 */
   // if (total_rows > 0)
   // {
   //    /* now create S */
   //    /* need to get new column start */
   //    {
   //       JX_Int global_start;
   //       jx_MPI_Scan(&big_m, &global_start, 1, JX_MPI_BIG_INT, MPI_SUM, comm);
   //       col_starts[0] = global_start - m;
   //       col_starts[1] = global_start;
   //    }

   //    /* We did nothing to A_offd, so all the data kept, just reorder them
   //     * The create function takes comm, global num rows/cols,
   //     *    row/col start, num cols offd, nnz diag, nnz offd
   //     */
   //    S_offd_nnz = jx_CSRMatrixNumNonzeros(A_offd);
   //    S_offd_ncols = jx_CSRMatrixNumCols(A_offd);

   //    matS = jx_ParCSRMatrixCreate(comm,
   //                                 total_rows,
   //                                 total_rows,
   //                                 col_starts,
   //                                 col_starts,
   //                                 S_offd_ncols,
   //                                 ctrS,
   //                                 S_offd_nnz);

   //    /* first put diagonal data in */
   //    S_diag = jx_ParCSRMatrixDiag(matS);

   //    jx_CSRMatrixI(S_diag) = S_diag_i;
   //    jx_CSRMatrixData(S_diag) = S_diag_data;
   //    jx_CSRMatrixJ(S_diag) = S_diag_j;

   //    /* now start to construct offdiag of S */
   //    S_offd = jx_ParCSRMatrixOffd(matS);
   //    S_offd_i = jx_TAlloc(JX_Int, m + 1);
   //    S_offd_j = jx_TAlloc(JX_Int, S_offd_nnz);
   //    S_offd_data = jx_TAlloc(JX_Real, S_offd_nnz);
   //    S_offd_colmap = jx_CTAlloc(JX_Int, S_offd_ncols);

   //    /* simply use a loop to copy data from A_offd */
   //    S_offd_i[0] = 0;
   //    k3 = 0;
   //    for (i = 1; i <= e; i++)
   //    {
   //       S_offd_i[i] = k3;
   //    }
   //    for (i = 0; i < m_e; i++)
   //    {
   //       col = perm[i + nI];
   //       k1 = A_offd_i[col];
   //       k2 = A_offd_i[col + 1];
   //       for (j = k1; j < k2; j++)
   //       {
   //          S_offd_j[k3] = A_offd_j[j];
   //          S_offd_data[k3++] = A_offd_data[j];
   //       }
   //       S_offd_i[i + 1 + e] = k3;
   //    }

   //    /* give I, J, DATA to S_offd */
   //    jx_CSRMatrixI(S_offd) = S_offd_i;
   //    jx_CSRMatrixJ(S_offd) = S_offd_j;
   //    jx_CSRMatrixData(S_offd) = S_offd_data;

   //    /* now we need to update S_offd_colmap */

   //    /* get total num of send */
   //    num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
   //    begin = jx_ParCSRCommPkgSendMapStart(comm_pkg, 0);
   //    end = jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
   //    send_buf = jx_TAlloc(JX_Int, end - begin);
   //    /* copy new index into send_buf */
   //    for (i = begin; i < end; i++)
   //    {
   //       send_buf[i - begin] = rperm[jx_ParCSRCommPkgSendMapElmt(comm_pkg, i)] -
   //                             nLU + col_starts[0];
   //    }
   //    /* main communication */
   //    comm_handle = jx_ParCSRCommHandleCreate(21, comm_pkg, send_buf, S_offd_colmap);
   //    jx_ParCSRCommHandleDestroy(comm_handle);

   //    /* setup index */
   //    jx_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

   //    jx_ILUSortOffdColmap(matS);

   //    /* free */
   //    jx_TFree(send_buf);
   // } /* end of forming S */
   // /* create S finished */

   matL = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0,
                                ctrL,
                                0);

   L_diag = jx_ParCSRMatrixDiag(matL);
   jx_CSRMatrixI(L_diag) = L_diag_i;
   if (ctrL)
   {
      jx_CSRMatrixData(L_diag) = L_diag_data;
      jx_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jx_TFree(L_diag_j);
      jx_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)ctrL;
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0,
                                ctrU,
                                0);

   U_diag = jx_ParCSRMatrixDiag(matU);
   jx_CSRMatrixI(U_diag) = U_diag_i;
   if (ctrU)
   {
      jx_CSRMatrixData(U_diag) = U_diag_data;
      jx_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jx_TFree(U_diag_j);
      jx_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)ctrU;
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matU) = total_nnz;
   /* free memory */
   jx_TFree(wL);
   jx_TFree(iw);
   if (!matS)
   {
      /* we allocate some memory for S, need to free if unused */
      jx_TFree(S_diag_i);
   }

   if (!permp)
   {
      jx_TFree(perm);
   }
   if (!qpermp)
   {
      jx_TFree(qperm);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;
   *u_end = u_end_array;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUKSymbolic
 *
 * Setup ILU(k) symbolic factorization
 *
 * n = total rows of input
 * lfil = level of fill-in, the k in ILU(k)
 * perm = permutation array indicating ordering of factorization.
 * rperm = reverse permutation array, used here to avoid duplicate memory allocation
 * iw = working array, used here to avoid duplicate memory allocation
 * nLU = size of computed LDU factorization.
 * A/L/U/S_diag_i = the I slot of A, L, U and S
 * A/L/U/S_diag_j = the J slot of A, L, U and S
 *
 * Will form global Schur Matrix if nLU < n
 *--------------------------------------------------------------------------*/
JX_Int
jx_ILUSetupILUKSymbolic(JX_Int n,
                        JX_Int *A_diag_i,
                        JX_Int *A_diag_j,
                        JX_Int lfil,
                        JX_Int *perm,
                        JX_Int *rperm,
                        JX_Int *iw,
                        JX_Int nLU,
                        JX_Int *L_diag_i,
                        JX_Int *U_diag_i,
                        JX_Int *S_diag_i,
                        JX_Int **L_diag_j,
                        JX_Int **U_diag_j,
                        JX_Int **S_diag_j,
                        JX_Int **u_end)
{
   /*
    * 1: Setup and create buffers
    * A_diag_*: tempory pointer for the diagonal matrix of A and its '*' slot
    * ii: outer loop from 0 to nLU - 1
    * i: the real col number in diag inside the outer loop
    * iw:  working array store the reverse of active col number
    * iL: working array store the active col number
    * iLev: working array store the active level of current row
    * lenl/u: current position in iw and so
    * ctrL/U/S: global position in J
    */

   JX_Int *temp_L_diag_j, *temp_U_diag_j, *temp_S_diag_j = NULL, *u_levels;
   JX_Int *iL, *iLev;
   JX_Int ii, i, j, k, ku, lena, lenl, lenu, lenh, ilev, lev, col, icol;
   JX_Int m = n - nLU;
   JX_Int *u_end_array;

   /* memory management */
   JX_Int ctrL;
   JX_Int ctrU;
   JX_Int ctrS;
   JX_Int capacity_L;
   JX_Int capacity_U;
   JX_Int capacity_S = 0;
   JX_Int initial_alloc = 0;
   JX_Int nnz_A;

   /* set iL and iLev to right place in iw array */
   iL = iw + n;
   iLev = iw + 2 * n;

   /* setup initial memory used */
   nnz_A = A_diag_i[n];
   if (n > 0)
   {
      initial_alloc = (JX_Int)(nLU + ceil((nnz_A / 2.0) * nLU / n));
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   /* allocate other memory for L and U struct */
   temp_L_diag_j = jx_CTAlloc(JX_Int, capacity_L);
   temp_U_diag_j = jx_CTAlloc(JX_Int, capacity_U);

   if (m > 0)
   {
      capacity_S = (JX_Int)(m + ceil(nnz_A / 2.0 * m / n));
      temp_S_diag_j = jx_CTAlloc(JX_Int, capacity_S);
   }

   u_end_array = jx_TAlloc(JX_Int, nLU);
   u_levels = jx_CTAlloc(JX_Int, capacity_U);
   ctrL = ctrU = ctrS = 0;

   /* set initial value for working array */
   for (ii = 0; ii < n; ii++)
   {
      iw[ii] = -1;
   }

   /*
    * 2: Start of main loop
    * those in iL are NEW col index (after permutation)
    */
   for (ii = 0; ii < nLU; ii++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0; /* this is the current length of heap */
      lenu = ii;
      lena = A_diag_i[i + 1];
      /* put those already inside original pattern, and set their level to 0 */
      for (j = A_diag_i[i]; j < lena; j++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if (col < ii)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jx_ILUMinHeapAddIIIi(iL, iLev, iw, lenh);
         }
         else if (col > ii)
         {
            /* this is an entry in U */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      } /* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while (lenh > 0)
      {
         /*
          * k is now the new col index after permutation
          * the first element of the heap is the smallest
          */
         k = iL[0];
         ilev = iLev[0];
         /*
          * we now need to maintain the heap structure
          */
         jx_ILUMinHeapRemoveIIIi(iL, iLev, iw, lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k] = -1;
         jx_swap2i(iL, iLev, ii - lenl, lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k + 1];
         for (j = U_diag_i[k]; j < ku; j++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if (lev > lfil)
            {
               continue;
            }
            if (icol < 0)
            {
               /* not yet in */
               if (col < ii)
               {
                  /*
                   * if we add to the left L, we need to maintian the
                   *    heap structure
                   */
                  iL[lenh] = col;
                  iLev[lenh] = lev;
                  iw[col] = lenh++;
                  /*swap it with the element right after the heap*/

                  /* maintain the heap */
                  jx_ILUMinHeapAddIIIi(iL, iLev, iw, lenh);
               }
               else if (col > ii)
               {
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jx_min(lev, iLev[icol]);
            }
         } /* end of loop j for level update */
      } /* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii + 1] = L_diag_i[ii] + lenl;
      if (lenl > 0)
      {
         /* check if memory is enough */
         while (ctrL + lenl > capacity_L)
         {
            JX_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jx_TReAlloc_v2(temp_L_diag_j, JX_Int, tmp, JX_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for (j = 0; j < lenl; j++)
         {
            temp_L_diag_j[ctrL + j] = iL[ii - j - 1];
         }
         ctrL += lenl;
      }
      k = lenu - ii;
      U_diag_i[ii + 1] = U_diag_i[ii] + k;
      if (k > 0)
      {
         /* check if memory is enough */
         while (ctrU + k > capacity_U)
         {
            JX_Int tmp = capacity_U;
            capacity_U = (JX_Int)(capacity_U * EXPAND_FACT + 1);
            temp_U_diag_j = jx_TReAlloc_v2(temp_U_diag_j, JX_Int, tmp, JX_Int, capacity_U);
            u_levels = jx_TReAlloc_v2(u_levels, JX_Int, tmp, JX_Int, capacity_U);
         }
         // jx_TMemcpy(temp_U_diag_j+ctrU,iL+ii,JX_Int,k);
         jx_TMemcpy(temp_U_diag_j + ctrU, iL + ii, JX_Int, k);
         jx_TMemcpy(u_levels + ctrU, iLev + ii, JX_Int, k);
         ctrU += k;
      }
      if (m > 0)
      {
         jx_qsort2i(temp_U_diag_j, u_levels, U_diag_i[ii], U_diag_i[ii + 1] - 1);
         jx_BinarySearch2(temp_U_diag_j, nLU, U_diag_i[ii], U_diag_i[ii + 1] - 1, u_end_array + ii);
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }

      /* reset iw */
      for (j = ii; j < lenu; j++)
      {
         iw[iL[j]] = -1;
      }

   } /* end of main loop ii from 0 to nLU-1 */

   /* another loop to set EU^-1 and Schur complement */
   for (ii = nLU; ii < n; ii++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0;   /* this is the current length of heap */
      lenu = nLU; /* now this stores S, start from nLU */
      lena = A_diag_i[i + 1];
      /* put those already inside original pattern, and set their level to 0 */
      for (j = A_diag_i[i]; j < lena; j++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if (col < nLU)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jx_ILUMinHeapAddIIIi(iL, iLev, iw, lenh);
         }
         else if (col != ii) /* we for sure to add ii, avoid duplicate */
         {
            /* this is an entry in S */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      } /* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while (lenh > 0)
      {
         /*
          * k is now the new col index after permutation
          * the first element of the heap is the smallest
          */
         k = iL[0];
         ilev = iLev[0];
         /*
          * we now need to maintain the heap structure
          */
         jx_ILUMinHeapRemoveIIIi(iL, iLev, iw, lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k] = -1;
         jx_swap2i(iL, iLev, nLU - lenl, lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k + 1];
         for (j = U_diag_i[k]; j < ku; j++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if (lev > lfil)
            {
               continue;
            }
            if (icol < 0)
            {
               /* not yet in */
               if (col < nLU)
               {
                  /*
                   * if we add to the left L, we need to maintian the
                   *    heap structure
                   */
                  iL[lenh] = col;
                  iLev[lenh] = lev;
                  iw[col] = lenh++;
                  /*swap it with the element right after the heap*/

                  /* maintain the heap */
                  jx_ILUMinHeapAddIIIi(iL, iLev, iw, lenh);
               }
               else if (col != ii)
               {
                  /* S part */
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jx_min(lev, iLev[icol]);
            }
         } /* end of loop j for level update */
      } /* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii + 1] = L_diag_i[ii] + lenl;
      if (lenl > 0)
      {
         /* check if memory is enough */
         while (ctrL + lenl > capacity_L)
         {
            JX_Int tmp = capacity_L;
            capacity_L = (JX_Int)(capacity_L * EXPAND_FACT + 1);
            temp_L_diag_j = jx_TReAlloc_v2(temp_L_diag_j, JX_Int, tmp, JX_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for (j = 0; j < lenl; j++)
         {
            temp_L_diag_j[ctrL + j] = iL[nLU - j - 1];
         }
         ctrL += lenl;
      }
      k = lenu - nLU + 1;
      /* check if memory is enough */
      while (ctrS + k > capacity_S)
      {
         JX_Int tmp = capacity_S;
         capacity_S = (JX_Int)(capacity_S * EXPAND_FACT + 1);
         temp_S_diag_j = jx_TReAlloc_v2(temp_S_diag_j, JX_Int, tmp, JX_Int, capacity_S);
      }
      temp_S_diag_j[ctrS] = ii; /* must have diagonal */
      // jx_TMemcpy(temp_S_diag_j+ctrS+1,iL+nLU,JX_Int,k-1);
      jx_TMemcpy(temp_S_diag_j + ctrS + 1, iL + nLU, JX_Int, k - 1);
      ctrS += k;
      S_diag_i[ii - nLU + 1] = ctrS;

      /* reset iw */
      for (j = nLU; j < lenu; j++)
      {
         iw[iL[j]] = -1;
      }

   } /* end of main loop ii from nLU to n-1 */

   /*
    * 3: Update the struct for L, U and S
    */
   for (k = nLU; k < n; k++)
   {
      U_diag_i[k + 1] = U_diag_i[nLU];
   }
   /*
    * 4: Finishing up and free memory
    */
   jx_TFree(u_levels);

   *L_diag_j = temp_L_diag_j;
   *U_diag_j = temp_U_diag_j;
   *S_diag_j = temp_S_diag_j;
   *u_end = u_end_array;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUK
 *
 * Setup ILU(k) numeric factorization
 *
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * permp: permutation array indicating ordering of factorization.
 *        Perm could come from a CF_marker array or a reordering routine.
 * qpermp: column permutation array.
 * nLU: size of computed LDU factorization.
 * nI: number of interial unknowns, nI should obey nI >= nLU
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sprt: Schur Complement, if no Schur Complement, it will be set to NULL
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetupILUK(jx_ParCSRMatrix *A,
                JX_Int lfil,
                JX_Int *permp,
                JX_Int *qpermp,
                JX_Int nLU,
                JX_Int nI,
                jx_ParCSRMatrix **Lptr,
                JX_Real **Dptr,
                jx_ParCSRMatrix **Uptr,
                jx_ParCSRMatrix **Sptr,
                JX_Int **u_end)
{
   /*
    * 1: Setup and create buffers
    * matL/U: the ParCSR matrix for L and U
    * L/U_diag: the diagonal csr matrix of matL/U
    * A_diag_*: tempory pointer for the diagonal matrix of A and its '*' slot
    * ii = outer loop from 0 to nLU - 1
    * i = the real col number in diag inside the outer loop
    * iw =  working array store the reverse of active col number
    * iL = working array store the active col number
    */

   /* call ILU0 if lfil is 0 */
   if (lfil == 0)
   {
      return jx_ILUSetupILU0(A, permp, qpermp, nLU, nI, Lptr, Dptr, Uptr, Sptr, u_end);
   }

   JX_Real local_nnz, total_nnz;
   JX_Int i, ii, j, k, k1, k2, k3, kl, ku, jpiv, col, icol;
   JX_Int *iw;
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   JX_Int num_procs, my_id;

   /* data objects for A */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);

   /* data objects for L, D, U */
   jx_ParCSRMatrix *matL;
   jx_ParCSRMatrix *matU;
   jx_CSRMatrix *L_diag;
   jx_CSRMatrix *U_diag;
   JX_Real *D_data;
   JX_Real *L_diag_data = NULL;
   JX_Int *L_diag_i;
   JX_Int *L_diag_j = NULL;
   JX_Real *U_diag_data = NULL;
   JX_Int *U_diag_i;
   JX_Int *U_diag_j = NULL;

   /* data objects for S */
   jx_ParCSRMatrix *matS = NULL;
   jx_CSRMatrix *S_diag;
   jx_CSRMatrix *S_offd;
   JX_Real *S_diag_data = NULL;
   JX_Int *S_diag_i = NULL;
   JX_Int *S_diag_j = NULL;
   JX_Int *S_offd_i = NULL;
   JX_Int *S_offd_j = NULL;
   JX_Int *S_offd_colmap = NULL;
   JX_Real *S_offd_data;
   JX_Int S_offd_nnz, S_offd_ncols;
   JX_Int col_starts[2];
   JX_Int total_rows;

   /* communication */
   jx_ParCSRCommPkg *comm_pkg;
   jx_ParCSRCommHandle *comm_handle;
   JX_Int *send_buf = NULL;

   /* problem size */
   JX_Int n;
   JX_Int m;
   JX_Int e;
   JX_Int m_e;

   /* reverse permutation array */
   JX_Int *rperm;
   JX_Int *perm, *qperm;

   /* start setup */
   /* check input and get problem size */
   n = jx_CSRMatrixNumRows(A_diag);
   if (nLU < 0 || nLU > n)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU out of range.\n");
   }
   m = n - nLU;
   e = nI - nLU;
   m_e = n - nI;
   if (e < 0)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU should not exceed nI.\n");
   }

   /* Init I array anyway. S's might be freed later */
   D_data = jx_CTAlloc(JX_Real, n);
   L_diag_i = jx_CTAlloc(JX_Int, (n + 1));
   U_diag_i = jx_CTAlloc(JX_Int, (n + 1));
   S_diag_i = jx_CTAlloc(JX_Int, (m + 1));

   /* set Comm_Pkg if not yet built */
   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);
   comm_pkg = jx_ParCSRMatrixCommPkg(A);
   if (!comm_pkg)
   {
      jx_MatvecCommPkgCreate(A);
      comm_pkg = jx_ParCSRMatrixCommPkg(A);
   }

   /*
    * 2: Symbolic factorization
    * setup iw and rperm first
    */
   /* allocate work arrays */
   iw = jx_CTAlloc(JX_Int, 4 * n);
   rperm = iw + 3 * n;
   L_diag_i[0] = U_diag_i[0] = S_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    */

   if (!permp)
   {
      perm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         perm[i] = i;
      }
   }
   else
   {
      perm = permp;
   }

   if (!qpermp)
   {
      qperm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         qperm[i] = i;
      }
   }
   else
   {
      qperm = qpermp;
   }

   for (i = 0; i < n; i++)
   {
      rperm[qperm[i]] = i;
   }

   /* do symbolic factorization */
   jx_ILUSetupILUKSymbolic(n, A_diag_i, A_diag_j, lfil, perm, rperm, iw,
                           nLU, L_diag_i, U_diag_i, S_diag_i, &L_diag_j, &U_diag_j, &S_diag_j, u_end);

   /*
    * after this, we have our I,J for L, U and S ready, and L sorted
    * iw are still -1 after symbolic factorization
    * now setup helper array here
    */
   if (L_diag_i[n])
   {
      L_diag_data = jx_CTAlloc(JX_Real, L_diag_i[n]);
   }
   if (U_diag_i[n])
   {
      U_diag_data = jx_CTAlloc(JX_Real, U_diag_i[n]);
   }
   if (S_diag_i[m])
   {
      S_diag_data = jx_CTAlloc(JX_Real, S_diag_i[m]);
   }

   /*
    * 3: Begin real factorization
    * we already have L and U structure ready, so no extra working array needed
    */
   /* first loop for upper part */
   for (ii = 0; ii < nLU; ii++)
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii + 1];
      ku = U_diag_i[ii + 1];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      /* set up working arrays */
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      D_data[ii] = 0.0;
      iw[ii] = ii;
      for (j = U_diag_i[ii]; j < ku; j++)
      {
         col = U_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, D and U */
      for (j = k1; j < k2; j++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if (col < ii)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else if (col == ii)
         {
            D_data[ii] = A_diag_data[j];
         }
         else
         {
            U_diag_data[icol] = A_diag_data[j];
         }
      }
      /* elimination */
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv + 1];

         for (k = U_diag_i[jpiv]; k < ku; k++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if (icol < 0)
            {
               /* not in partern */
               continue;
            }
            if (col < ii)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j] * U_diag_data[k];
            }
            else if (col == ii)
            {
               /* diag part */
               D_data[icol] -= L_diag_data[j] * U_diag_data[k];
            }
            else
            {
               /* U part */
               U_diag_data[icol] -= L_diag_data[j] * U_diag_data[k];
            }
         }
      }
      /* reset working array */
      ku = U_diag_i[ii + 1];
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      iw[ii] = -1;
      for (j = U_diag_i[ii]; j < ku; j++)
      {
         col = U_diag_j[j];
         iw[col] = -1;
      }

      /* diagonal part (we store the inverse) */
      if (jx_abs(D_data[ii]) < MAT_TOL)
      {
         D_data[ii] = 1.0e-06;
      }
      D_data[ii] = 1. / D_data[ii];
   }

   /* Now lower part for Schur complement */
   for (ii = nLU; ii < n; ii++)
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii + 1];
      ku = S_diag_i[ii - nLU + 1];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      /* set up working arrays */
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      for (j = S_diag_i[ii - nLU]; j < ku; j++)
      {
         col = S_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, and S */
      for (j = k1; j < k2; j++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if (col < nLU)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else
         {
            S_diag_data[icol] = A_diag_data[j];
         }
      }
      /* elimination */
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv + 1];
         for (k = U_diag_i[jpiv]; k < ku; k++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if (icol < 0)
            {
               /* not in partern */
               continue;
            }
            if (col < nLU)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j] * U_diag_data[k];
            }
            else
            {
               /* S part */
               S_diag_data[icol] -= L_diag_data[j] * U_diag_data[k];
            }
         }
      }
      /* reset working array */
      for (j = L_diag_i[ii]; j < kl; j++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      ku = S_diag_i[ii - nLU + 1];
      for (j = S_diag_i[ii - nLU]; j < ku; j++)
      {
         col = S_diag_j[j];
         iw[col] = -1;
         /* remember to update index, S is smaller! */
         S_diag_j[j] -= nLU;
      }
   }

   /*
    * 4: Finishing up and free
    */

   /* First create Schur complement if necessary
    * Check if we need to create Schur complement
    */
   JX_Int big_m = (JX_Int)m;
   // jx_MPI_Allreduce(&big_m, &total_rows, 1, JX_MPI_BIG_INT, MPI_SUM, comm);
   jx_MPI_Allreduce(&big_m, &total_rows, 1, JX_MPI_INT, MPI_SUM, comm);
   /* only form when total_rows > 0 */
   if (total_rows > 0)
   {
      /* now create S */
      /* need to get new column start */
      {
         JX_Int global_start;
         // jx_MPI_Scan(&big_m, &global_start, 1, JX_MPI_BIG_INT, MPI_SUM, comm);
         jx_MPI_Scan(&big_m, &global_start, 1, JX_MPI_INT, MPI_SUM, comm);
         col_starts[0] = global_start - m;
         col_starts[1] = global_start;
      }

      /* We did nothing to A_offd, so all the data kept, just reorder them
       * The create function takes comm, global num rows/cols,
       *    row/col start, num cols offd, nnz diag, nnz offd
       */
      S_offd_nnz = jx_CSRMatrixNumNonzeros(A_offd);
      S_offd_ncols = jx_CSRMatrixNumCols(A_offd);

      matS = jx_ParCSRMatrixCreate(comm,
                                   total_rows,
                                   total_rows,
                                   col_starts,
                                   col_starts,
                                   S_offd_ncols,
                                   S_diag_i[m],
                                   S_offd_nnz);

      /* first put diagonal data in */
      S_diag = jx_ParCSRMatrixDiag(matS);

      jx_CSRMatrixI(S_diag) = S_diag_i;
      jx_CSRMatrixData(S_diag) = S_diag_data;
      jx_CSRMatrixJ(S_diag) = S_diag_j;

      /* now start to construct offdiag of S */
      S_offd = jx_ParCSRMatrixOffd(matS);
      S_offd_i = jx_TAlloc(JX_Int, m + 1);
      S_offd_j = jx_TAlloc(JX_Int, S_offd_nnz);
      S_offd_data = jx_TAlloc(JX_Real, S_offd_nnz);
      S_offd_colmap = jx_CTAlloc(JX_Int, S_offd_ncols);

      /* simply use a loop to copy data from A_offd */
      S_offd_i[0] = 0;
      k3 = 0;
      for (i = 1; i <= e; i++)
      {
         S_offd_i[i + 1] = k3;
      }
      for (i = 0; i < m_e; i++)
      {
         col = perm[i + nI];
         k1 = A_offd_i[col];
         k2 = A_offd_i[col + 1];
         for (j = k1; j < k2; j++)
         {
            S_offd_j[k3] = A_offd_j[j];
            S_offd_data[k3++] = A_offd_data[j];
         }
         S_offd_i[i + e + 1] = k3;
      }

      /* give I, J, DATA to S_offd */
      jx_CSRMatrixI(S_offd) = S_offd_i;
      jx_CSRMatrixJ(S_offd) = S_offd_j;
      jx_CSRMatrixData(S_offd) = S_offd_data;

      /* now we need to update S_offd_colmap */

      /* get total num of send */
      JX_Int num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
      JX_Int begin = jx_ParCSRCommPkgSendMapStart(comm_pkg, 0);
      JX_Int end = jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
      send_buf = jx_TAlloc(JX_Int, end - begin);

      /* copy new index into send_buf */
      for (i = begin; i < end; i++)
      {
         send_buf[i - begin] = rperm[jx_ParCSRCommPkgSendMapElmt(comm_pkg, i)] - nLU + col_starts[0];
      }

      /* main communication */
      comm_handle = jx_ParCSRCommHandleCreate(21, comm_pkg, send_buf, S_offd_colmap);
      jx_ParCSRCommHandleDestroy(comm_handle);

      /* setup index */
      jx_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

      jx_ILUSortOffdColmap(matS);

      /* free */
      jx_TFree(send_buf);
   } /* end of forming S */

   /* Assemble LDU matrices */
   /* zero out unfactored rows */
   for (k = nLU; k < n; k++)
   {
      D_data[k] = 1.;
   }

   matL = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0 /* num_cols_offd */,
                                L_diag_i[n],
                                0 /* num_nonzeros_offd */);

   L_diag = jx_ParCSRMatrixDiag(matL);
   jx_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[n] > 0)
   {
      jx_CSRMatrixData(L_diag) = L_diag_data;
      jx_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jx_TFree(L_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)(L_diag_i[n]);
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0,
                                U_diag_i[n],
                                0);

   U_diag = jx_ParCSRMatrixDiag(matU);
   jx_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[n] > 0)
   {
      jx_CSRMatrixData(U_diag) = U_diag_data;
      jx_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jx_TFree(U_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)(U_diag_i[n]);
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free */
   jx_TFree(iw);
   if (!matS)
   {
      /* we allocate some memory for S, need to free if unused */
      jx_TFree(S_diag_i);
   }

   if (!permp)
   {
      jx_TFree(perm);
   }

   if (!qpermp)
   {
      jx_TFree(qperm);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUT
 *
 * Setup ILU(t) numeric factorization
 *
 * A: input matrix
 * lfil: maximum nnz per row in L and U
 * tol: droptol array in ILUT
 *    tol[0]: matrix B
 *    tol[1]: matrix E and F
 *    tol[2]: matrix S
 * perm: permutation array indicating ordering of factorization.
 *       Perm could come from a CF_marker array or a reordering routine.
 * qperm: permutation array for column
 * nLU: size of computed LDU factorization.
 *      If nLU < n, Schur complement will be formed
 * nI: number of interial unknowns. nLU should obey nLU <= nI.
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sptr: Schur complement
 *
 * Keep the largest lfil entries that is greater than some tol relative
 *    to the input tol and the norm of that row in both L and U
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetupILUT(jx_ParCSRMatrix *A,
                JX_Int lfil,
                JX_Real *tol,
                JX_Int *permp,
                JX_Int *qpermp,
                JX_Int nLU,
                JX_Int nI,
                jx_ParCSRMatrix **Lptr,
                JX_Real **Dptr,
                jx_ParCSRMatrix **Uptr,
                jx_ParCSRMatrix **Sptr,
                JX_Int **u_end)
{
   /*
    * 1: Setup and create buffers
    * matL/U: the ParCSR matrix for L and U
    * L/U_diag: the diagonal csr matrix of matL/U
    * A_diag_*: tempory pointer for the diagonal matrix of A and its '*' slot
    * ii = outer loop from 0 to nLU - 1
    * i = the real col number in diag inside the outer loop
    * iw =  working array store the reverse of active col number
    * iL = working array store the active col number
    */
   JX_Real local_nnz, total_nnz;
   JX_Int i, ii, j, k, k1, k2, k3, kl, ku, col, icol, lenl, lenu, lenhu, lenhlr,
       lenhll, jpos, jrow;
   JX_Real inorm, itolb, itolef, itols, dpiv, lxu;
   JX_Int *iw, *iL;
   JX_Real *w;

   /* memory management */
   JX_Int ctrL;
   JX_Int ctrU;
   JX_Int initial_alloc = 0;
   JX_Int capacity_L;
   JX_Int capacity_U;
   JX_Int ctrS;
   JX_Int capacity_S = 0;
   JX_Int nnz_A;

   /* communication stuffs for S */
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   JX_Int S_offd_nnz, S_offd_ncols;
   jx_ParCSRCommPkg *comm_pkg;
   jx_ParCSRCommHandle *comm_handle;
   JX_Int num_procs, my_id;
   JX_Int col_starts[2];
   JX_Int total_rows;
   JX_Int num_sends;
   JX_Int begin, end;

   /* data objects for A */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);

   /* data objects for L, D, U */
   jx_ParCSRMatrix *matL;
   jx_ParCSRMatrix *matU;
   jx_CSRMatrix *L_diag;
   jx_CSRMatrix *U_diag;
   JX_Real *D_data;
   JX_Real *L_diag_data = NULL;
   JX_Int *L_diag_i;
   JX_Int *L_diag_j = NULL;
   JX_Real *U_diag_data = NULL;
   JX_Int *U_diag_i;
   JX_Int *U_diag_j = NULL;

   /* data objects for S */
   jx_ParCSRMatrix *matS = NULL;
   jx_CSRMatrix *S_diag;
   jx_CSRMatrix *S_offd;
   JX_Real *S_diag_data = NULL;
   JX_Int *S_diag_i = NULL;
   JX_Int *S_diag_j = NULL;
   JX_Int *S_offd_i = NULL;
   JX_Int *S_offd_j = NULL;
   JX_Int *S_offd_colmap = NULL;
   JX_Real *S_offd_data;
   JX_Int *send_buf = NULL;
   JX_Int *u_end_array;

   /* reverse permutation */
   JX_Int *rperm;
   JX_Int *perm, *qperm;

   /* problem size
    * m is n - nLU, num of rows of local Schur system
    * m_e is the size of interface nodes
    * e is the number of interial rows in local Schur Complement
    */
   JX_Int n;
   JX_Int m;
   JX_Int e;
   JX_Int m_e;

   /* start setup
    * check input first
    */
   n = jx_CSRMatrixNumRows(A_diag);
   if (nLU < 0 || nLU > n)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU out of range.\n");
   }
   m = n - nLU;
   e = nI - nLU;
   m_e = n - nI;
   if (e < 0)
   {
      jx_error_w_msg(JX_ERROR_ARG, "WARNING: nLU should not exceed nI.\n");
   }

   u_end_array = jx_TAlloc(JX_Int, nLU);

   /* start set up
    * setup communication stuffs first
    */
   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);
   comm_pkg = jx_ParCSRMatrixCommPkg(A);
   /* create if not yet built */
   if (!comm_pkg)
   {
      jx_MatvecCommPkgCreate(A);
      comm_pkg = jx_ParCSRMatrixCommPkg(A);
   }

   /* setup initial memory, in ILUT, just guess with max nnz per row */
   nnz_A = A_diag_i[nLU];
   if (n > 0)
   {
      initial_alloc = (JX_Int)(jx_min(nLU + ceil((nnz_A / 2.0) * nLU / n),
                                      nLU * lfil));
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   D_data = jx_CTAlloc(JX_Real, n);
   L_diag_i = jx_CTAlloc(JX_Int, (n + 1));
   U_diag_i = jx_CTAlloc(JX_Int, (n + 1));

   L_diag_j = jx_CTAlloc(JX_Int, capacity_L);
   U_diag_j = jx_CTAlloc(JX_Int, capacity_U);
   L_diag_data = jx_CTAlloc(JX_Real, capacity_L);
   U_diag_data = jx_CTAlloc(JX_Real, capacity_U);

   ctrL = ctrU = 0;

   ctrS = 0;
   S_diag_i = jx_CTAlloc(JX_Int, (m + 1));
   S_diag_i[0] = 0;

   /* only setup S part when n > nLU */
   if (m > 0)
   {
      capacity_S = (JX_Int)(jx_min(m + ceil((nnz_A / 2.0) * m / n), m * lfil));
      S_diag_j = jx_CTAlloc(JX_Int, capacity_S);
      S_diag_data = jx_CTAlloc(JX_Real, capacity_S);
   }

   /* setting up working array */
   iw = jx_CTAlloc(JX_Int, 3 * n);
   iL = iw + n;
   w = jx_CTAlloc(JX_Real, n);
   for (i = 0; i < n; i++)
   {
      iw[i] = -1;
   }
   L_diag_i[0] = U_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    * rperm[old] -> new
    * perm[new]  -> old
    */
   rperm = iw + 2 * n;

   if (!permp)
   {
      perm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         perm[i] = i;
      }
   }
   else
   {
      perm = permp;
   }

   if (!qpermp)
   {
      qperm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         qperm[i] = i;
      }
   }
   else
   {
      qperm = qpermp;
   }

   for (i = 0; i < n; i++)
   {
      rperm[perm[i]] = i;
   }
   /*
    * 2: Main loop of elimination
    * maintain two heaps
    * |----->*********<-----|-----*********|
    * |col heap***value heap|value in U****|
    */

   /* main outer loop for upper part */
   for (ii = 0; ii < nLU; ii++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      kl = ii - 1;
      /* reset row norm of ith row */
      inorm = .0;
      for (j = k1; j < k2; j++)
      {
         inorm += jx_abs(A_diag_data[j]);
      }
      if (inorm == .0)
      {
         jx_error_w_msg(JX_ERROR_ARG, "WARNING: ILUT with zero row.\n");
      }
      inorm /= (JX_Real)(k2 - k1);
      /* set the scaled tol for that row */
      itolb = tol[0] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      w[ii] = 0.0;
      iw[ii] = ii;
      /* copy in data from A */
      for (j = k1; j < k2; j++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if (col < ii)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jx_ILUMinHeapAddIRIi(iL, w, iw, lenhll);
         }
         else if (col == ii)
         {
            w[ii] = A_diag_data[j];
         }
         else
         {
            lenu++;
            jpos = lenu + ii;
            iL[jpos] = col;
            w[jpos] = A_diag_data[j];
            iw[col] = jpos;
         }
      }

      /*
       * main elimination
       * need to maintain 2 heaps for L, one heap for col and one heaps for value
       * maintian an array for U, and do qsplit with quick sort after that
       * while the heap of col is greater than zero
       */
      while (lenhll > 0)
      {

         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jx_ILUMinHeapRemoveIRIi(iL, w, iw, lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jx_swap2(iL, w, lenhll, kl - lenhlr);
         lenhlr++;
         jx_ILUMaxrHeapAddRabsI(w + kl, iL + kl, lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow + 1];
         for (j = U_diag_i[jrow]; j < ku; j++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = -dpiv * U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if ((icol == -1) &&
                ((col < nLU && jx_abs(lxu) < itolb) || (col >= nLU && jx_abs(lxu) < itolef)))
            {
               continue;
            }
            if (icol == -1)
            {
               if (col < ii)
               {
                  /* L part
                   * not already in L part
                   * put it to the end of heap
                   * might overwrite some small entries, no issue
                   */
                  iL[lenhll] = col;
                  w[lenhll] = lxu;
                  iw[col] = lenhll++;
                  /* add to heap, by col number */
                  jx_ILUMinHeapAddIRIi(iL, w, iw, lenhll);
               }
               else if (col == ii)
               {
                  w[ii] += lxu;
               }
               else
               {
                  /*
                   * not already in U part
                   * put is to the end of heap
                   */
                  lenu++;
                  jpos = lenu + ii;
                  iL[jpos] = col;
                  w[jpos] = lxu;
                  iw[col] = jpos;
               }
            }
            else
            {
               w[icol] += lxu;
            }
         }
      } /* while loop for the elimination of current row */

      if (jx_abs(w[ii]) < MAT_TOL)
      {
         w[ii] = 1.0e-06;
      }
      D_data[ii] = 1. / w[ii];
      iw[ii] = -1;

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii + 1] = L_diag_i[ii] + lenl;
      if (lenl > 0)
      {
         /* test if memory is enough */
         while (ctrL + lenl > capacity_L)
         {
            JX_Int tmp = capacity_L;
            capacity_L = (JX_Int)(capacity_L * EXPAND_FACT + 1);
            L_diag_j = jx_TReAlloc_v2(L_diag_j, JX_Int, tmp, JX_Int,
                                      capacity_L);
            L_diag_data = jx_TReAlloc_v2(L_diag_data, JX_Real, tmp, JX_Real,
                                         capacity_L);
         }
         ctrL += lenl;

         /* copy large data in */
         for (j = L_diag_i[ii]; j < ctrL; j++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jx_ILUMaxrHeapRemoveRabsI(w + kl, iL + kl, lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only U part
       */
      ku = lenu + ii;
      for (j = ii + 1; j <= ku; j++)
      {
         iw[iL[j]] = -1;
      }

      if (lenu < lfil)
      {
         /* we simply keep all of the data, no need to sort */
         lenhu = lenu;
      }
      else
      {
         /* need to sort the first small(hopefully) part of it */
         lenhu = lfil;
         /* quick split, only sort the first small part of the array */
         jx_ILUMaxQSplitRabsI(w, iL, ii + 1, ii + lenhu, ii + lenu);
      }

      U_diag_i[ii + 1] = U_diag_i[ii] + lenhu;
      if (lenhu > 0)
      {
         /* test if memory is enough */
         while (ctrU + lenhu > capacity_U)
         {
            JX_Int tmp = capacity_U;
            capacity_U = (JX_Int)(capacity_U * EXPAND_FACT + 1);
            U_diag_j = jx_TReAlloc_v2(U_diag_j, JX_Int, tmp, JX_Int,
                                      capacity_U);
            U_diag_data = jx_TReAlloc_v2(U_diag_data, JX_Real, tmp, JX_Real,
                                         capacity_U);
         }
         ctrU += lenhu;
         /* copy large data in */
         for (j = U_diag_i[ii]; j < ctrU; j++)
         {
            jpos = ii + 1 + j - U_diag_i[ii];
            U_diag_j[j] = iL[jpos];
            U_diag_data[j] = w[jpos];
         }
      }
      /* check and build u_end array */
      if (m > 0)
      {
         jx_qsort1(U_diag_j, U_diag_data, U_diag_i[ii], U_diag_i[ii + 1] - 1);
         jx_BinarySearch2(U_diag_j, nLU, U_diag_i[ii], U_diag_i[ii + 1] - 1, u_end_array + ii);
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }
   } /* end of ii loop from 0 to nLU-1 */

   /* now main loop for Schur comlement part */
   for (ii = nLU; ii < n; ii++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      kl = nLU - 1;
      /* reset row norm of ith row */
      inorm = .0;
      for (j = k1; j < k2; j++)
      {
         inorm += jx_abs(A_diag_data[j]);
      }
      if (inorm == .0)
      {
         jx_error_w_msg(JX_ERROR_ARG, "WARNING: ILUT with zero row.\n");
      }
      inorm /= (JX_Real)(k2 - k1);
      /* set the scaled tol for that row */
      itols = tol[2] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      /* copy in data from A */
      for (j = k1; j < k2; j++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if (col < nLU)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jx_ILUMinHeapAddIRIi(iL, w, iw, lenhll);
         }
         else if (col == ii)
         {
            /* the diagonla entry of S */
            iL[nLU] = col;
            w[nLU] = A_diag_data[j];
            iw[col] = nLU;
         }
         else
         {
            /* S part of it */
            lenu++;
            jpos = lenu + nLU;
            iL[jpos] = col;
            w[jpos] = A_diag_data[j];
            iw[col] = jpos;
         }
      }

      /*
       * main elimination
       * need to maintain 2 heaps for L, one heap for col and one heaps for value
       * maintian an array for S, and do qsplit with quick sort after that
       * while the heap of col is greater than zero
       */
      while (lenhll > 0)
      {
         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jx_ILUMinHeapRemoveIRIi(iL, w, iw, lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jx_swap2(iL, w, lenhll, kl - lenhlr);
         lenhlr++;
         jx_ILUMaxrHeapAddRabsI(w + kl, iL + kl, lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow + 1];
         for (j = U_diag_i[jrow]; j < ku; j++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = -dpiv * U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if ((icol == -1) &&
                ((col < nLU && jx_abs(lxu) < itolef) ||
                 (col >= nLU && jx_abs(lxu) < itols)))
            {
               continue;
            }
            if (icol == -1)
            {
               if (col < nLU)
               {
                  /* L part
                   * not already in L part
                   * put it to the end of heap
                   * might overwrite some small entries, no issue
                   */
                  iL[lenhll] = col;
                  w[lenhll] = lxu;
                  iw[col] = lenhll++;
                  /* add to heap, by col number */
                  jx_ILUMinHeapAddIRIi(iL, w, iw, lenhll);
               }
               else if (col == ii)
               {
                  /* the diagonla entry of S */
                  iL[nLU] = col;
                  w[nLU] = A_diag_data[j];
                  iw[col] = nLU;
               }
               else
               {
                  /*
                   * not already in S part
                   * put is to the end of heap
                   */
                  lenu++;
                  jpos = lenu + nLU;
                  iL[jpos] = col;
                  w[jpos] = lxu;
                  iw[col] = jpos;
               }
            }
            else
            {
               w[icol] += lxu;
            }
         }
      } /* while loop for the elimination of current row */

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii + 1] = L_diag_i[ii] + lenl;
      if (lenl > 0)
      {
         /* test if memory is enough */
         while (ctrL + lenl > capacity_L)
         {
            JX_Int tmp = capacity_L;
            capacity_L = (JX_Int)(capacity_L * EXPAND_FACT + 1);
            L_diag_j = jx_TReAlloc_v2(L_diag_j, JX_Int, tmp, JX_Int,
                                      capacity_L);
            L_diag_data = jx_TReAlloc_v2(L_diag_data, JX_Real, tmp, JX_Real,
                                         capacity_L);
         }
         ctrL += lenl;

         /* copy large data in */
         for (j = L_diag_i[ii]; j < ctrL; j++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jx_ILUMaxrHeapRemoveRabsI(w + kl, iL + kl, lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only S part
       */
      ku = lenu + nLU;
      for (j = nLU; j <= ku; j++)
      {
         iw[iL[j]] = -1;
      }

      /* no dropping at this point of time for S */
      // lenhu = lenu < lfil ? lenu : lfil;
      lenhu = lenu;
      /* quick split, only sort the first small part of the array */
      jx_ILUMaxQSplitRabsI(w, iL, nLU + 1, nLU + lenhu, nLU + lenu);
      /* we have diagonal in S anyway */
      /* test if memory is enough */
      while (ctrS + lenhu + 1 > capacity_S)
      {
         JX_Int tmp = capacity_S;
         capacity_S = (JX_Int)(capacity_S * EXPAND_FACT + 1);
         S_diag_j = jx_TReAlloc_v2(S_diag_j, JX_Int, tmp,
                                   JX_Int, capacity_S);
         S_diag_data = jx_TReAlloc_v2(S_diag_data, JX_Real, tmp,
                                      JX_Real, capacity_S);
      }

      ctrS += (lenhu + 1);
      S_diag_i[ii - nLU + 1] = ctrS;

      /* copy large data in, diagonal first */
      S_diag_j[S_diag_i[ii - nLU]] = iL[nLU] - nLU;
      S_diag_data[S_diag_i[ii - nLU]] = w[nLU];
      for (j = S_diag_i[ii - nLU] + 1; j < ctrS; j++)
      {
         jpos = nLU + j - S_diag_i[ii - nLU];
         S_diag_j[j] = iL[jpos] - nLU;
         S_diag_data[j] = w[jpos];
      }
   } /* end of ii loop from nLU to n-1 */

   /*
    * 3: Finishing up and free
    */

   /* First create Schur complement if necessary
    * Check if we need to create Schur complement
    */
   JX_Int big_m = (JX_Int)m;
   // jx_MPI_Allreduce(&big_m, &total_rows, 1, JX_MPI_BIG_INT, MPI_SUM, comm);
   jx_MPI_Allreduce(&big_m, &total_rows, 1, JX_MPI_INT, MPI_SUM, comm);

   /* only form when total_rows > 0 */
   if (total_rows > 0)
   {
      /* now create S */
      /* need to get new column start */
      {
         JX_Int global_start;
         // jx_MPI_Scan(&big_m, &global_start, 1, JX_MPI_BIG_INT, MPI_SUM, comm);
         jx_MPI_Scan(&big_m, &global_start, 1, JX_MPI_INT, MPI_SUM, comm);
         col_starts[0] = global_start - m;
         col_starts[1] = global_start;
      }
      /* We did nothing to A_offd, so all the data kept, just reorder them
       * The create function takes comm, global num rows/cols,
       *    row/col start, num cols offd, nnz diag, nnz offd
       */
      S_offd_nnz = jx_CSRMatrixNumNonzeros(A_offd);
      S_offd_ncols = jx_CSRMatrixNumCols(A_offd);

      matS = jx_ParCSRMatrixCreate(comm,
                                   total_rows,
                                   total_rows,
                                   col_starts,
                                   col_starts,
                                   S_offd_ncols,
                                   S_diag_i[m],
                                   S_offd_nnz);

      /* first put diagonal data in */
      S_diag = jx_ParCSRMatrixDiag(matS);

      jx_CSRMatrixI(S_diag) = S_diag_i;
      jx_CSRMatrixData(S_diag) = S_diag_data;
      jx_CSRMatrixJ(S_diag) = S_diag_j;

      /* now start to construct offdiag of S */
      S_offd = jx_ParCSRMatrixOffd(matS);
      S_offd_i = jx_TAlloc(JX_Int, m + 1);
      S_offd_j = jx_TAlloc(JX_Int, S_offd_nnz);
      S_offd_data = jx_TAlloc(JX_Real, S_offd_nnz);
      S_offd_colmap = jx_CTAlloc(JX_Int, S_offd_ncols);

      /* simply use a loop to copy data from A_offd */
      S_offd_i[0] = 0;
      k3 = 0;
      for (i = 1; i <= e; i++)
      {
         S_offd_i[i] = k3;
      }
      for (i = 0; i < m_e; i++)
      {
         col = perm[i + nI];
         k1 = A_offd_i[col];
         k2 = A_offd_i[col + 1];
         for (j = k1; j < k2; j++)
         {
            S_offd_j[k3] = A_offd_j[j];
            S_offd_data[k3++] = A_offd_data[j];
         }
         S_offd_i[i + e + 1] = k3;
      }

      /* give I, J, DATA to S_offd */
      jx_CSRMatrixI(S_offd) = S_offd_i;
      jx_CSRMatrixJ(S_offd) = S_offd_j;
      jx_CSRMatrixData(S_offd) = S_offd_data;

      /* now we need to update S_offd_colmap */

      /* get total num of send */
      num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
      begin = jx_ParCSRCommPkgSendMapStart(comm_pkg, 0);
      end = jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends);
      send_buf = jx_TAlloc(JX_Int, end - begin);
      /* copy new index into send_buf */
      for (i = begin; i < end; i++)
      {
         send_buf[i - begin] = rperm[jx_ParCSRCommPkgSendMapElmt(comm_pkg, i)] - nLU + col_starts[0];
      }

      /* main communication */
      comm_handle = jx_ParCSRCommHandleCreate(21, comm_pkg, send_buf, S_offd_colmap);
      /* need this to synchronize, Isend & Irecv used in above functions */
      jx_ParCSRCommHandleDestroy(comm_handle);

      /* setup index */
      jx_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

      jx_ILUSortOffdColmap(matS);

      /* free */
      jx_TFree(send_buf);
   } /* end of forming S */

   /* now start to construct L and U */
   for (k = nLU; k < n; k++)
   {
      /* set U after nLU to be 0, and diag to be one */
      U_diag_i[k + 1] = U_diag_i[nLU];
      D_data[k] = 1.;
   }

   /* create parcsr matrix */
   matL = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0,
                                L_diag_i[n],
                                0);

   L_diag = jx_ParCSRMatrixDiag(matL);
   jx_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[n] > 0)
   {
      jx_CSRMatrixData(L_diag) = L_diag_data;
      jx_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jx_TFree(L_diag_j);
      jx_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)(L_diag_i[n]);
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jx_ParCSRMatrixCreate(comm,
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixGlobalNumRows(A),
                                jx_ParCSRMatrixRowStarts(A),
                                jx_ParCSRMatrixColStarts(A),
                                0,
                                U_diag_i[n],
                                0);

   U_diag = jx_ParCSRMatrixDiag(matU);
   jx_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[n] > 0)
   {
      jx_CSRMatrixData(U_diag) = U_diag_data;
      jx_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jx_TFree(U_diag_j);
      jx_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JX_Real)(U_diag_i[n]);
   jx_MPI_Allreduce(&local_nnz, &total_nnz, 1, JX_MPI_REAL, MPI_SUM, comm);
   jx_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free working array */
   jx_TFree(iw);
   jx_TFree(w);

   if (!matS)
   {
      jx_TFree(S_diag_i);
   }

   if (!permp)
   {
      jx_TFree(perm);
   }

   if (!qpermp)
   {
      jx_TFree(qperm);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;
   *u_end = u_end_array;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_NSHSetup
 *--------------------------------------------------------------------------*/
JX_Int
jx_NSHSetup(void *nsh_vdata,
            jx_ParCSRMatrix *A,
            jx_ParVector *f,
            jx_ParVector *u)
{
   jx_printf("[JX-ILU] WARNING: jx_NSHSetup is called. \n");

   return jx_error_flag;
}

JX_Int
jx_ILUSetupILU0RAS(jx_ParCSRMatrix *A,
                   JX_Int *perm,
                   JX_Int nLU,
                   jx_ParCSRMatrix **Lptr,
                   JX_Real **Dptr,
                   jx_ParCSRMatrix **Uptr)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUSetupILU0RAS is called. \n");

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUKRASSymbolic
 *
 * ILU(k) symbolic factorization for RAS
 *
 * n = total rows of input
 * lfil = level of fill-in, the k in ILU(k)
 * perm = permutation array indicating ordering of factorization.
 * rperm = reverse permutation array, used here to avoid duplicate memory allocation
 * iw = working array, used here to avoid duplicate memory allocation
 * nLU = size of computed LDU factorization.
 * A/L/U/E_i = the I slot of A, L, U and E
 * A/L/U/E_j = the J slot of A, L, U and E
 *
 * Will form global Schur Matrix if nLU < n
 *--------------------------------------------------------------------------*/
JX_Int
jx_ILUSetupILUKRASSymbolic(JX_Int n,
                           JX_Int *A_diag_i,
                           JX_Int *A_diag_j,
                           JX_Int *A_offd_i,
                           JX_Int *A_offd_j,
                           JX_Int *E_i,
                           JX_Int *E_j,
                           JX_Int ext,
                           JX_Int lfil,
                           JX_Int *perm,
                           JX_Int *rperm,
                           JX_Int *iw,
                           JX_Int nLU,
                           JX_Int *L_diag_i,
                           JX_Int *U_diag_i,
                           JX_Int **L_diag_j,
                           JX_Int **U_diag_j)
{
   jx_printf("[JX-ILU] ERROR: jx_ILUSetupILUKRASSymbolic is called. \n");

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUKRAS
 *
 * ILU(k) numeric factorization for RAS
 *
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * perm: permutation array indicating ordering of factorization.
 *       Perm could come from a CF_marker array or a reordering routine.
 * nLU: size of computed LDU factorization.
 * Lptr, Dptr, Uptr: L, D, U factors.
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetupILUKRAS(jx_ParCSRMatrix *A,
                   JX_Int lfil,
                   JX_Int *perm,
                   JX_Int nLU,
                   jx_ParCSRMatrix **Lptr,
                   JX_Real **Dptr,
                   jx_ParCSRMatrix **Uptr)
{
   jx_printf("[JX-ILU] ERROR: jx_ILUSetupILUKRAS is called. \n");

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetupILUTRAS
 *
 * ILUT for RAS
 *
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * tol: droptol array in ILUT
 *    tol[0]: matrix B
 *    tol[1]: matrix E and F
 *    tol[2]: matrix S
 * perm: permutation array indicating ordering of factorization.
 *       Perm could come from a CF_marker: array or a reordering routine.
 * nLU: size of computed LDU factorization. If nLU < n, Schur compelemnt will be formed
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sptr: Schur complement
 *
 * Keep the largest lfil entries that is greater than some tol relative
 *    to the input tol and the norm of that row in both L and U
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetupILUTRAS(jx_ParCSRMatrix *A,
                   JX_Int lfil,
                   JX_Real *tol,
                   JX_Int *perm,
                   JX_Int nLU,
                   jx_ParCSRMatrix **Lptr,
                   JX_Real **Dptr,
                   jx_ParCSRMatrix **Uptr)
{
   jx_printf("[JX-ILU] ERROR: jx_ILUSetupILUTRAS is called. \n");
}
