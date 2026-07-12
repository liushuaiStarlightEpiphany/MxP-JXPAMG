#include "jx_pamg.h"

JX_Int
jx_PAMGBuildClassicalFix(jx_ParCSRMatrix *par_A,
                         JX_Int *CF_marker,
                         jx_ParCSRMatrix *par_S,
                         JX_Real *eta,
                         JX_Real eta_threshold,
                         JX_Int *num_cpts_global,
                         JX_Int num_functions,
                         JX_Int *dof_func,
                         JX_Int debug_flag,
                         JX_Real trunc_factor,
                         JX_Int max_elmts,
                         JX_Int *col_offd_S_to_A,
                         jx_ParCSRMatrix **P_ptr,
                         JX_Int *num_fpts_global_ptr,
                         JX_Int *num_ei_fpts_global_ptr)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_A);
   JX_Int my_id, num_procs;
   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &num_procs);
   JX_Int n_fine = jx_CSRMatrixNumRows(jx_ParCSRMatrixDiag(par_A));
   JX_Int local_f = 0, local_ei = 0;

   for (JX_Int i = 0; i < n_fine; i++)
      if (CF_marker[i] < 0 && CF_marker[i] != -3)
      { local_f++; if (eta && eta[i] >= eta_threshold) local_ei++; }

   if (my_id == 0)
      jx_printf("[ClsFix] F=%d, EI=%d, threshold=%.4f\n", local_f, local_ei, eta_threshold);

   /* Stage 1: classical for ALL F-points */
   jx_PAMGBuildInterp(par_A, CF_marker, par_S, num_cpts_global,
                       num_functions, dof_func, debug_flag,
                       trunc_factor, max_elmts, col_offd_S_to_A, P_ptr);

   /* Stage 2: for EI points, replace with extended interpolation
      Build extended P, copy EI rows, destroy extended P */
   if (local_ei > 0)
   {
      jx_ParCSRMatrix *P_ext = NULL;
      jx_PAMGBuildExtPIInterp(par_A, CF_marker, par_S, num_cpts_global,
                               num_functions, dof_func, debug_flag,
                               trunc_factor, max_elmts, col_offd_S_to_A, &P_ext);

      jx_ParCSRMatrix *P = *P_ptr;
      jx_CSRMatrix *cd = jx_ParCSRMatrixDiag(P);
      jx_CSRMatrix *co = jx_ParCSRMatrixOffd(P);
      jx_CSRMatrix *ed = jx_ParCSRMatrixDiag(P_ext);
      jx_CSRMatrix *eo = jx_ParCSRMatrixOffd(P_ext);
      JX_Int *cdi = jx_CSRMatrixI(cd), *edi = jx_CSRMatrixI(ed);
      JX_Int *cdj = jx_CSRMatrixJ(cd), *edj = jx_CSRMatrixJ(ed);
      JX_Real *cda = jx_CSRMatrixData(cd), *eda = jx_CSRMatrixData(ed);
      JX_Int *coi = co ? jx_CSRMatrixI(co) : NULL, *eoi = eo ? jx_CSRMatrixI(eo) : NULL;
      JX_Int *coj = co ? jx_CSRMatrixJ(co) : NULL, *eoj = eo ? jx_CSRMatrixJ(eo) : NULL;
      JX_Real *cod = co ? jx_CSRMatrixData(co) : NULL, *eod = eo ? jx_CSRMatrixData(eo) : NULL;

      jx_ParCSRMatrixColMapOffd(P_ext) = NULL;

      for (JX_Int i = 0; i < n_fine; i++)
      {
         if (CF_marker[i] >= 0 || CF_marker[i] == -3) continue;
         if (!eta || eta[i] < eta_threshold) continue;

         JX_Int cds = cdi[i], cde = cdi[i+1], cdn = cde - cds;
         JX_Int eds = edi[i], ede = edi[i+1], edn = ede - eds;
         JX_Int cons = 0, eons = 0, conn = 0, eonn = 0;
         if (coi && eoi) { cons = coi[i]; conn = coi[i+1] - cons; eons = eoi[i]; eonn = eoi[i+1] - eons; }

         JX_Int ediff = edn - cdn, odiff = eonn - conn;
         if (ediff > 0 || odiff > 0)
         {
            JX_Int nr = jx_CSRMatrixNumRows(cd);
            if (ediff > 0)
            {
               JX_Int nz = jx_CSRMatrixNumNonzeros(cd);
               jx_CSRMatrixJ(cd) = jx_TReAlloc(jx_CSRMatrixJ(cd), JX_Int, nz + ediff);
               jx_CSRMatrixData(cd) = jx_TReAlloc(jx_CSRMatrixData(cd), JX_Real, nz + ediff);
               cdj = jx_CSRMatrixJ(cd); cda = jx_CSRMatrixData(cd);
               for (JX_Int k = nz - 1; k >= cde; k--)
               { cdj[k+ediff] = cdj[k]; cda[k+ediff] = cda[k]; }
               for (JX_Int k = i+1; k <= nr; k++) cdi[k] += ediff;
            }
            if (odiff > 0 && coi)
            {
               JX_Int nz = jx_CSRMatrixNumNonzeros(co);
               jx_CSRMatrixJ(co) = jx_TReAlloc(jx_CSRMatrixJ(co), JX_Int, nz + odiff);
               jx_CSRMatrixData(co) = jx_TReAlloc(jx_CSRMatrixData(co), JX_Real, nz + odiff);
               coj = jx_CSRMatrixJ(co); cod = jx_CSRMatrixData(co);
               JX_Int coe = coi[i+1];
               for (JX_Int k = nz - 1; k >= coe; k--)
               { coj[k+odiff] = coj[k]; cod[k+odiff] = cod[k]; }
               for (JX_Int k = i+1; k <= jx_CSRMatrixNumRows(co); k++) coi[k] += odiff;
            }
         }

         for (JX_Int k = 0; k < edn; k++)
         { cdj[cdi[i] + k] = edj[eds + k]; cda[cdi[i] + k] = eda[eds + k]; }
         if (coj && eoj) for (JX_Int k = 0; k < eonn; k++)
         { coj[coi[i] + k] = eoj[eons + k]; cod[coi[i] + k] = eod[eons + k]; }
      }
      jx_ParCSRMatrixDestroy(P_ext);
   }

   if (num_fpts_global_ptr)
   { JX_Int g=0; jx_MPI_Allreduce(&local_f, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_fpts_global_ptr = g; }
   if (num_ei_fpts_global_ptr)
   { JX_Int g=0; jx_MPI_Allreduce(&local_ei, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_ei_fpts_global_ptr = g; }
   return 0;
}
