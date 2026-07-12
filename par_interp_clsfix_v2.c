//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//========================================================================//

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
   JX_Int local_num_fpts = 0, local_num_ei_fpts = 0;

   for (JX_Int i = 0; i < n_fine; i++)
   {
      if (CF_marker[i] >= 0 || CF_marker[i] == -3) continue;
      local_num_fpts++;
      if (eta && eta[i] >= eta_threshold) local_num_ei_fpts++;
   }

   /* Stage 1: build classical P */
   jx_ParCSRMatrix *P_cls = NULL;
   jx_PAMGBuildInterp(par_A, CF_marker, par_S,
                       num_cpts_global, num_functions,
                       dof_func, debug_flag,
                       trunc_factor, max_elmts,
                       col_offd_S_to_A, &P_cls);

   /* Stage 2: for EI points, replace rows with extended interpolation */
   jx_ParCSRMatrix *P_ext = NULL;
   jx_PAMGBuildExtPIInterp(par_A, CF_marker, par_S,
                            num_cpts_global, num_functions,
                            dof_func, debug_flag,
                            trunc_factor, max_elmts,
                            col_offd_S_to_A, &P_ext);

   /* Copy EI rows from extended to classical */
   jx_CSRMatrix *cl_d = jx_ParCSRMatrixDiag(P_cls);
   jx_CSRMatrix *cl_o = jx_ParCSRMatrixOffd(P_cls);
   jx_CSRMatrix *ex_d = jx_ParCSRMatrixDiag(P_ext);
   jx_CSRMatrix *ex_o = jx_ParCSRMatrixOffd(P_ext);
   JX_Int *cl_di = jx_CSRMatrixI(cl_d);
   JX_Int *cl_dj = jx_CSRMatrixJ(cl_d);
   JX_Real *cl_da = jx_CSRMatrixData(cl_d);
   JX_Int *cl_oi = cl_o ? jx_CSRMatrixI(cl_o) : NULL;
   JX_Int *cl_oj = cl_o ? jx_CSRMatrixJ(cl_o) : NULL;
   JX_Real *cl_od = cl_o ? jx_CSRMatrixData(cl_o) : NULL;
   JX_Int *ex_di = jx_CSRMatrixI(ex_d);
   JX_Int *ex_dj = jx_CSRMatrixJ(ex_d);
   JX_Real *ex_da = jx_CSRMatrixData(ex_d);
   JX_Int *ex_oi = ex_o ? jx_CSRMatrixI(ex_o) : NULL;
   JX_Int *ex_oj = ex_o ? jx_CSRMatrixJ(ex_o) : NULL;
   JX_Real *ex_od = ex_o ? jx_CSRMatrixData(ex_o) : NULL;

   jx_ParCSRMatrixColMapOffd(P_ext) = NULL;

   for (JX_Int i = 0; i < n_fine; i++)
   {
      if (CF_marker[i] >= 0 || CF_marker[i] == -3) continue;
      if (!eta || eta[i] < eta_threshold) continue;

      JX_Int cl_ds = cl_di[i], cl_de = cl_di[i + 1], cl_dn = cl_de - cl_ds;
      JX_Int ex_ds = ex_di[i], ex_de = ex_di[i + 1], ex_dn = ex_de - ex_ds;
      JX_Int cl_on = 0, ex_on = 0, cl_os = 0, ex_os = 0;

      if (cl_oi && ex_oi)
      {
         cl_os = cl_oi[i]; cl_on = cl_oi[i + 1] - cl_os;
         ex_os = ex_oi[i]; ex_on = ex_oi[i + 1] - ex_os;
      }

      JX_Int ed = ex_dn - cl_dn, eo = ex_on - cl_on;

      if (ed > 0 || eo > 0)
      {
         /* Expand row in classical P */
         JX_Int rn = jx_CSRMatrixNumRows(cl_d);
         if (ed > 0)
         {
            JX_Int sz = jx_CSRMatrixNumNonzeros(cl_d);
            jx_CSRMatrixJ(cl_d) = jx_TReAlloc(jx_CSRMatrixJ(cl_d), JX_Int, sz + ed);
            jx_CSRMatrixData(cl_d) = jx_TReAlloc(jx_CSRMatrixData(cl_d), JX_Real, sz + ed);
            cl_dj = jx_CSRMatrixJ(cl_d);
            cl_da = jx_CSRMatrixData(cl_d);
            for (JX_Int k = sz - 1; k >= cl_de; k--)
            { cl_dj[k + ed] = cl_dj[k]; cl_da[k + ed] = cl_da[k]; }
            for (JX_Int k = i + 1; k <= rn; k++) cl_di[k] += ed;
            cl_de += ed; cl_dn = ex_dn;
         }
         if (eo > 0 && cl_oi)
         {
            JX_Int sz = jx_CSRMatrixNumNonzeros(cl_o);
            jx_CSRMatrixJ(cl_o) = jx_TReAlloc(jx_CSRMatrixJ(cl_o), JX_Int, sz + eo);
            jx_CSRMatrixData(cl_o) = jx_TReAlloc(jx_CSRMatrixData(cl_o), JX_Real, sz + eo);
            cl_oj = jx_CSRMatrixJ(cl_o);
            cl_od = jx_CSRMatrixData(cl_o);
            JX_Int rno = jx_CSRMatrixNumRows(cl_o);
            JX_Int cl_oe = cl_oi[i + 1];
            for (JX_Int k = sz - 1; k >= cl_oe; k--)
            { cl_oj[k + eo] = cl_oj[k]; cl_od[k + eo] = cl_od[k]; }
            for (JX_Int k = i + 1; k <= rno; k++) cl_oi[k] += eo;
         }
      }

      for (JX_Int k = 0; k < ex_dn; k++)
      { cl_dj[cl_ds + k] = ex_dj[ex_ds + k]; cl_da[cl_ds + k] = ex_da[ex_ds + k]; }
      if (cl_oj && ex_oj) for (JX_Int k = 0; k < ex_on; k++)
      { cl_oj[cl_os + k] = ex_oj[ex_os + k]; cl_od[cl_os + k] = ex_od[ex_os + k]; }
   }

   jx_ParCSRMatrixDestroy(P_ext);

   *P_ptr = P_cls;

   if (num_fpts_global_ptr)
   { JX_Int g = 0; jx_MPI_Allreduce(&local_num_fpts, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_fpts_global_ptr = g; }
   if (num_ei_fpts_global_ptr)
   { JX_Int g = 0; jx_MPI_Allreduce(&local_num_ei_fpts, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_ei_fpts_global_ptr = g; }

   return 0;
}
