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

   /* 经典插值 */
   jx_PAMGBuildInterp(par_A, CF_marker, par_S, num_cpts_global,
                       num_functions, dof_func, debug_flag,
                       trunc_factor, max_elmts, col_offd_S_to_A, P_ptr);

   if (my_id == 0)
      jx_printf("[ClsFix] F=%d, EI=%d, threshold=%.4f\n", local_f, local_ei, eta_threshold);

   if (num_fpts_global_ptr)
   { JX_Int g = 0; jx_MPI_Allreduce(&local_f, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_fpts_global_ptr = g; }
   if (num_ei_fpts_global_ptr)
   { JX_Int g = 0; jx_MPI_Allreduce(&local_ei, &g, 1, JX_MPI_INT, MPI_SUM, comm); *num_ei_fpts_global_ptr = g; }
   return 0;
}
