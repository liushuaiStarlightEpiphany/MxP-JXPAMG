//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu_solve.c
 *
 */

#include "jx_ilu.h"
#include "jx_krylov.h"

static JX_Int
jx_ILUApplyLocalArrayILU0(jx_CSRMatrix *A,
                          JX_Int *indexD,
                          JX_Int *indexLU,
                          JX_Real *valueLU,
                          JX_Real *rhs,
                          JX_Real *sol,
                          JX_Real *work)
{
   JX_Int n = jx_CSRMatrixNumRows(A);
   JX_Int *IA = jx_CSRMatrixI(A);
   JX_Int i, j, row_end;
   JX_Real tmp;

   if (n <= 0)
   {
      return jx_error_flag;
   }

   work[0] = rhs[0];
   for (i = 1; i < n; i++)
   {
      tmp = rhs[i];
      row_end = indexD[i];
      for (j = IA[i]; j < row_end; j++)
      {
         tmp -= valueLU[j] * work[indexLU[j]];
      }
      work[i] = tmp;
   }

   row_end = indexD[n - 1];
   sol[n - 1] = work[n - 1] * valueLU[row_end];
   for (i = n - 2; i >= 0; i--)
   {
      tmp = work[i];
      row_end = IA[i + 1];
      for (j = indexD[i] + 1; j < row_end; j++)
      {
         tmp -= valueLU[j] * sol[indexLU[j]];
      }
      sol[i] = tmp * valueLU[indexD[i]];
   }

   return jx_error_flag;
}

static JX_Int
jx_ILUApplyLocalResidualCorrection(jx_ParCSRMatrix *A,
                                   jx_ParVector *f,
                                   jx_ParVector *u,
                                   jx_ParVector *ftemp,
                                   JX_Int *marker,
                                   JX_Int *core_indices,
                                   JX_Int core_size,
                                   JX_Int *local_indices,
                                   JX_Int *local_map,
                                   jx_CSRMatrix *local_mat,
                                   JX_Int *local_indexD,
                                   JX_Int *local_indexLU,
                                   JX_Real *local_valueLU,
                                   JX_Real *local_rhs,
                                   JX_Real *local_sol,
                                   JX_Real *local_work,
                                   JX_Int correction_type,
                                   JX_Int max_iter,
                                   JX_Real threshold,
                                   JX_Real omega,
                                   JX_Int print_level,
                                   JX_Int *apply_count,
                                   JX_Real *apply_time,
                                   JX_Real *matvec_time,
                                   JX_Real *local_solve_time,
                                   JX_Real *delta_norm_ratio_sum,
                                   JX_Real *residual_ratio_sum,
                                   JX_Real *residual_ratio_min,
                                   JX_Real *residual_ratio_max,
                                   JX_Int *accept_like_count,
                                   JX_Int *reject_like_count)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   jx_Vector *u_local = jx_ParVectorLocalVector(u);
   jx_Vector *r_local = jx_ParVectorLocalVector(ftemp);
   JX_Real *u_data = jx_VectorData(u_local);
   JX_Real *r_data = jx_VectorData(r_local);
   JX_Real alpha = -1.0;
   JX_Real beta = 1.0;
   JX_Real diag;
   JX_Int local_selected;
   JX_Int global_selected;
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int local_n = 0;
   JX_Int iter, i, j, my_id;
   JX_Int pos;
   JX_Real t0, t1, t2, t3, t4;
   JX_Real before_local_sq, before_global_sq, before_norm;
   JX_Real after_local_sq, after_global_sq, after_norm;
   JX_Real delta_local_sq, delta_global_sq, delta_norm;
   JX_Real sol_local_sq, sol_global_sq, sol_norm;
   JX_Real ratio;

   if (correction_type == 0 || max_iter <= 0 || omega == 0.0)
   {
      return jx_error_flag;
   }

   if (correction_type != 1 && correction_type != 2 &&
       correction_type != 3 && correction_type != 4)
   {
      return jx_error_flag;
   }

   (void)threshold;

   if (correction_type == 2 && !marker)
   {
      return jx_error_flag;
   }

   if (correction_type == 3 || correction_type == 4)
   {
      if (!core_indices || !local_indices || !local_map || !local_mat || !local_indexD || !local_indexLU ||
          !local_valueLU || !local_rhs || !local_sol || !local_work)
      {
         return jx_error_flag;
      }
      local_n = jx_CSRMatrixNumRows(local_mat);
      if (local_n <= 0)
      {
         return jx_error_flag;
      }
      if (core_size <= 0)
      {
         return jx_error_flag;
      }
   }

   jx_MPI_Comm_rank(comm, &my_id);

   for (iter = 0; iter < max_iter; iter++)
   {
      t0 = jx_MPI_Wtime();
      jx_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);
      t1 = jx_MPI_Wtime();
      if ((correction_type == 3 || correction_type == 4) && matvec_time)
      {
         *matvec_time += t1 - t0;
      }

      if (correction_type == 3 || correction_type == 4)
      {
         before_local_sq = 0.0;
         for (i = 0; i < n; i++)
         {
            before_local_sq += r_data[i] * r_data[i];
         }
         jx_MPI_Allreduce(&before_local_sq, &before_global_sq, 1,
                          JX_MPI_REAL, MPI_SUM, comm);
         before_norm = sqrt(before_global_sq);

         for (i = 0; i < local_n; i++)
         {
            local_rhs[i] = r_data[local_indices[i]];
            local_sol[i] = 0.0;
            local_work[i] = 0.0;
         }

         jx_ILUApplyLocalArrayILU0(local_mat, local_indexD, local_indexLU,
                                   local_valueLU, local_rhs, local_sol,
                                   local_work);
         t2 = jx_MPI_Wtime();

         delta_local_sq = 0.0;
         for (i = 0; i < core_size; i++)
         {
            pos = local_map[core_indices[i]];
            if (pos >= 0)
            {
               JX_Real delta_i = omega * local_sol[pos];
               delta_local_sq += delta_i * delta_i;
               u_data[core_indices[i]] += omega * local_sol[pos];
            }
         }
         jx_MPI_Allreduce(&delta_local_sq, &delta_global_sq, 1,
                          JX_MPI_REAL, MPI_SUM, comm);
         delta_norm = sqrt(delta_global_sq);

         sol_local_sq = 0.0;
         for (i = 0; i < n; i++)
         {
            sol_local_sq += u_data[i] * u_data[i];
         }
         jx_MPI_Allreduce(&sol_local_sq, &sol_global_sq, 1,
                          JX_MPI_REAL, MPI_SUM, comm);
         sol_norm = sqrt(sol_global_sq);

         t3 = jx_MPI_Wtime();
         jx_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);
         t4 = jx_MPI_Wtime();
         if (matvec_time)
         {
            *matvec_time += t4 - t3;
         }
         after_local_sq = 0.0;
         for (i = 0; i < n; i++)
         {
            after_local_sq += r_data[i] * r_data[i];
         }
         jx_MPI_Allreduce(&after_local_sq, &after_global_sq, 1,
                          JX_MPI_REAL, MPI_SUM, comm);
         after_norm = sqrt(after_global_sq);
         ratio = before_norm > 0.0 ? after_norm / before_norm : 1.0;

         if (delta_norm_ratio_sum)
         {
            *delta_norm_ratio_sum += sol_norm > 0.0 ? delta_norm / sol_norm : 0.0;
         }
         if (residual_ratio_sum)
         {
            *residual_ratio_sum += ratio;
         }
         if (residual_ratio_min && ratio < *residual_ratio_min)
         {
            *residual_ratio_min = ratio;
         }
         if (residual_ratio_max && ratio > *residual_ratio_max)
         {
            *residual_ratio_max = ratio;
         }
         if (ratio < 1.0)
         {
            if (accept_like_count)
            {
               (*accept_like_count)++;
            }
         }
         else if (reject_like_count)
         {
            (*reject_like_count)++;
         }

         if (apply_count)
         {
            (*apply_count)++;
         }
         if (apply_time)
         {
            *apply_time += t4 - t0;
         }
         if (local_solve_time)
         {
            *local_solve_time += t2 - t1;
         }

         if (print_level > 2)
         {
            local_selected = core_size;
            jx_MPI_Allreduce(&local_selected, &global_selected, 1, JX_MPI_INT, MPI_SUM, comm);
            if (my_id == 0)
            {
               jx_printf("    ILU local subspace residual correction %2d: solved rows = %d, corrected core rows = %d\n",
                         iter + 1, local_n, global_selected);
            }
         }

         continue;
      }

      local_selected = 0;
      for (i = 0; i < n; i++)
      {
         diag = 0.0;
         for (j = A_diag_i[i]; j < A_diag_i[i + 1]; j++)
         {
            if (A_diag_j[j] == i)
            {
               diag = A_diag_data[j];
               break;
            }
         }

         if (correction_type == 2 && marker[i] == 0)
         {
            continue;
         }

         if (jx_abs(diag) > DIVIDE_TOL)
         {
            u_data[i] += omega * r_data[i] / diag;
            local_selected++;
         }
      }

      if (print_level > 2)
      {
         jx_MPI_Allreduce(&local_selected, &global_selected, 1, JX_MPI_INT, MPI_SUM, comm);
         if (my_id == 0)
         {
            jx_printf("    ILU linear residual correction %2d: corrected rows = %d\n",
                      iter + 1, global_selected);
         }
      }
   }

   return jx_error_flag;
}

JX_Int
jx_ILUSolve(void *ilu_vdata,
            jx_ParCSRMatrix *A,
            jx_ParVector *f,
            jx_ParVector *u)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   /* Matrices */
   jx_ParCSRMatrix *matmL = jx_ParILUDataMatLModified(ilu_data);
   jx_ParCSRMatrix *matmU = jx_ParILUDataMatUModified(ilu_data);
   jx_ParCSRMatrix *matA = jx_ParILUDataMatA(ilu_data);
   jx_ParCSRMatrix *matL = jx_ParILUDataMatL(ilu_data);
   jx_ParCSRMatrix *matU = jx_ParILUDataMatU(ilu_data);
   jx_ParCSRMatrix *matS = jx_ParILUDataMatS(ilu_data);
   JX_Real *matD = jx_ParILUDataMatD(ilu_data);
   JX_Real *matmD = jx_ParILUDataMatDModified(ilu_data);

   /* Vectors */
   JX_Int ilu_type = jx_ParILUDataIluType(ilu_data);
   JX_Int *perm = jx_ParILUDataPerm(ilu_data);
   JX_Int *qperm = jx_ParILUDataQPerm(ilu_data);
   jx_ParVector *F_array = jx_ParILUDataF(ilu_data);
   jx_ParVector *U_array = jx_ParILUDataU(ilu_data);
   // jx_ParVector *D_array = jx_ParILUDataD(ilu_data);

   /* Solver settings */
   JX_Real tol = jx_ParILUDataTol(ilu_data);
   JX_Int logging = jx_ParILUDataLogging(ilu_data);
   JX_Int print_level = jx_ParILUDataPrintLevel(ilu_data);
   JX_Int max_iter = jx_ParILUDataMaxIter(ilu_data);
   JX_Int tri_solve = jx_ParILUDataTriSolve(ilu_data);
   JX_Int lower_jacobi_iters = jx_ParILUDataLowerJacobiIters(ilu_data);
   JX_Int upper_jacobi_iters = jx_ParILUDataUpperJacobiIters(ilu_data);
   JX_Int lrc_type = jx_ParILUDataLocalResidualCorrectionType(ilu_data);
   JX_Int lrc_max_iter = jx_ParILUDataLocalResidualCorrectionMaxIter(ilu_data);
   JX_Real lrc_threshold = jx_ParILUDataLocalResidualCorrectionThreshold(ilu_data);
   JX_Real lrc_omega = jx_ParILUDataLocalResidualCorrectionOmega(ilu_data);
   JX_Int *lrc_marker = jx_ParILUDataLocalResidualCorrectionMarker(ilu_data);
   JX_Int *lrc_core_indices = jx_ParILUDataLocalResidualCorrectionCoreIndices(ilu_data);
   JX_Int lrc_core_size = jx_ParILUDataLocalResidualCorrectionCoreSize(ilu_data);
   JX_Int *lrc_local_indices = jx_ParILUDataLocalResidualCorrectionLocalIndices(ilu_data);
   JX_Int *lrc_local_map = jx_ParILUDataLocalResidualCorrectionLocalMap(ilu_data);
   jx_CSRMatrix *lrc_local_mat = jx_ParILUDataLocalResidualCorrectionLocalMat(ilu_data);
   JX_Int *lrc_local_indexD = jx_ParILUDataLocalResidualCorrectionLocalIndexD(ilu_data);
   JX_Int *lrc_local_indexLU = jx_ParILUDataLocalResidualCorrectionLocalIndexLU(ilu_data);
   JX_Real *lrc_local_valueLU = jx_ParILUDataLocalResidualCorrectionLocalValueLU(ilu_data);
   JX_Real *lrc_local_rhs = jx_ParILUDataLocalResidualCorrectionLocalRhs(ilu_data);
   JX_Real *lrc_local_sol = jx_ParILUDataLocalResidualCorrectionLocalSol(ilu_data);
   JX_Real *lrc_local_work = jx_ParILUDataLocalResidualCorrectionLocalWork(ilu_data);
   // JX_Int IR_iters = jx_ParILUDataIRIters(ilu_data);
   JX_Real *norms = jx_ParILUDataRelResNorms(ilu_data);
   jx_ParVector *Ftemp = jx_ParILUDataFTemp(ilu_data);
   jx_ParVector *Utemp = jx_ParILUDataUTemp(ilu_data);
   jx_ParVector *Xtemp = jx_ParILUDataXTemp(ilu_data);
   jx_ParVector *Ytemp = jx_ParILUDataYTemp(ilu_data);
   JX_Real *fext = jx_ParILUDataFExt(ilu_data);
   JX_Real *uext = jx_ParILUDataUExt(ilu_data);
   jx_ParVector *residual = NULL;
   JX_Real alpha = -1.0;
   JX_Real beta = 1.0;
   JX_Real conv_factor = 0.0;
   JX_Real resnorm = 1.0;
   JX_Real init_resnorm = 0.0;
   JX_Real rel_resnorm;
   JX_Real rhs_norm = 0.0;
   JX_Real old_resnorm;
   JX_Real ieee_check = 0.0;
   JX_Real operat_cmplxty = jx_ParILUDataOperatorComplexity(ilu_data);
   JX_Int Solve_err_flag;
   JX_Int iter, num_procs, my_id;

   /* problem size */
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int nLU = n; // jx_ParILUDataNLU(ilu_data);
   JX_Int *u_end = jx_ParILUDataUEnd(ilu_data);

   /* Schur system solve */
   // JX_Solver schur_solver = jx_ParILUDataSchurSolver(ilu_data);
   // JX_Solver schur_precond = jx_ParILUDataSchurPrecond(ilu_data);
   jx_ParVector *rhs = jx_ParILUDataRhs(ilu_data);
   jx_ParVector *x = jx_ParILUDataX(ilu_data);

   if (logging > 1)
   {
      residual = jx_ParILUDataResidual(ilu_data);
   }

   jx_ParILUDataNumIterations(ilu_data) = 0;

   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);

   /*-----------------------------------------------------------------------
    *    Write the solver parameters 打印 ILU 求解器参数
    *-----------------------------------------------------------------------*/
   // if (my_id == 0 && print_level > 1)
   // {
   //    jx_ILUWriteSolverParams(ilu_data);
   // }

   /*-----------------------------------------------------------------------
    *    Initialize the solver error flag
    *-----------------------------------------------------------------------*/

   Solve_err_flag = 0;

   /*-----------------------------------------------------------------------
    *     write some initial info
    *-----------------------------------------------------------------------*/

   if (my_id == 0 && print_level > 1 && tol > 0.)
   {
      jx_printf("\n\n ILU SOLVER SOLUTION INFO:\n");
   }

   /*-----------------------------------------------------------------------
    *    Compute initial residual and print
    *-----------------------------------------------------------------------*/

   if (print_level > 1 || logging > 1 || tol > 0.)
   {
      if (logging > 1)
      {
         jx_ParVectorCopy(f, residual);
         if (tol > 0.0)
         {
            jx_ParCSRMatrixMatvec(alpha, A, u, beta, residual);
         }
         resnorm = sqrt(jx_ParVectorInnerProd(residual, residual));
      }
      else
      {
         jx_ParVectorCopy(f, Ftemp);
         if (tol > 0.0)
         {
            jx_ParCSRMatrixMatvec(alpha, A, u, beta, Ftemp);
         }
         resnorm = sqrt(jx_ParVectorInnerProd(Ftemp, Ftemp));
      }

      /* Since it does not diminish performance, attempt to return an error flag
         and notify users when they supply bad input. */
      if (resnorm != 0.)
      {
         ieee_check = resnorm / resnorm; /* INF -> NaN conversion */
      }
      if (ieee_check != ieee_check)
      {
         /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
            for ieee_check self-equality works on all IEEE-compliant compilers/
            machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
            by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
            found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
         if (print_level > 0)
         {
            jx_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
            jx_printf("ERROR -- jx_ILUSolve: INFs and/or NaNs detected in input.\n");
            jx_printf("User probably placed non-numerics in supplied A, x_0, or b.\n");
            jx_printf("ERROR detected by Hypre ...  END\n\n\n");
         }
         jx_error(JX_ERROR_GENERIC);

         return jx_error_flag;
      }

      init_resnorm = resnorm;
      rhs_norm = sqrt(jx_ParVectorInnerProd(f, f));
      if (rhs_norm > JX_REAL_EPSILON)
      {
         rel_resnorm = init_resnorm / rhs_norm;
      }
      else
      {
         /* rhs is zero, return a zero solution */
         jx_ParVectorSetConstantValues(U_array, 0.0);
         if (logging > 0)
         {
            rel_resnorm = 0.0;
            jx_ParILUDataFinalRelResidualNorm(ilu_data) = rel_resnorm;
         }

         return jx_error_flag;
      }
   }
   else
   {
      rel_resnorm = 1.;
   }

   if (my_id == 0 && print_level > 1)
   {
      jx_printf("                                            relative\n");
      jx_printf("               residual        factor       residual\n");
      jx_printf("               --------        ------       --------\n");
      jx_printf("    Initial    %e                 %e\n", init_resnorm,
                rel_resnorm);
   }

   matA = A;
   U_array = u;
   F_array = f;

   /************** Main Solver Loop - always do 1 iteration ************/
   iter = 0;

   /* 只要没有收敛，并且没有超过最大迭代次数，就继续做 ILU solve */
   while ((rel_resnorm >= tol || iter < 1) && iter < max_iter)
   {
      /* Do one solve on LUe=r */
      switch (ilu_type)
      {
      case 0:
      case 1:
         /* Basic Block Jacobi ILU */
         if (tri_solve == 1)
         {
            /* Direct triangular solve */
            jx_ILUSolveLU(matA, F_array, U_array, perm, n,
                          matL, matD, matU, Utemp, Ftemp);
         }
         else
         {
            /* Iterative triangular solve */
            jx_ILUSolveLUIter(matA, F_array, U_array, perm, n,
                              matL, matD, matU, Utemp, Ftemp,
                              lower_jacobi_iters, upper_jacobi_iters);
         }
         break;
      // case 10:
      // case 11:
      //    jx_ILUSolveSchurGMRES(matA, F_array, U_array, perm, perm, nLU, matL, matD, matU, matS,
      //                          Utemp, Ftemp, schur_solver, schur_precond, rhs, x, u_end); // GMRES
      //    break;
      // case 20:
      // case 21:
      //    jx_ILUSolveSchurNSH(matA, F_array, U_array, perm, nLU, matL, matD, matU, matS,
      //                        Utemp, Ftemp, schur_solver, rhs, x, u_end); // MR+NSH
      //    break;
      // case 30:
      // case 31:
      //    jx_ILUSolveLURAS(matA, F_array, U_array, perm, matL, matD, matU, Utemp, Utemp, fext, uext); // RAS
      //    break;
      // case 40:
      // case 41:
      //    jx_ILUSolveSchurGMRES(matA, F_array, U_array, perm, qperm, nLU, matL, matD, matU, matS,
      //                          Utemp, Ftemp, schur_solver, schur_precond, rhs, x, u_end); // GMRES
      //    break;
      default:
         jx_ILUSolveLU(matA, F_array, U_array, perm, n, matL, matD, matU, Utemp, Ftemp); // BJ
         break;
      }

      /* Type 4 is interleaved: after each BJ-ILU residual-correction sweep,
       * recompute the current residual and apply the fixed local ILU update.
       * Type 3 remains a postprocess correction after the base loop. */
      if (lrc_type != 3)
      {
         jx_ILUApplyLocalResidualCorrection(matA, F_array, U_array, Ftemp,
                                            lrc_marker,
                                            lrc_core_indices,
                                            lrc_core_size,
                                            lrc_local_indices,
                                            lrc_local_map,
                                            lrc_local_mat,
                                            lrc_local_indexD,
                                            lrc_local_indexLU,
                                            lrc_local_valueLU,
                                            lrc_local_rhs,
                                            lrc_local_sol,
                                            lrc_local_work,
                                            lrc_type, lrc_max_iter,
                                            lrc_threshold, lrc_omega,
                                            print_level,
                                            &jx_ParILUDataLocalResidualCorrectionApplyCount(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionApplyTime(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionMatvecTime(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionLocalSolveTime(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionDeltaNormRatioSum(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionResidualRatioSum(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionResidualRatioMin(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionResidualRatioMax(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionAcceptLikeCount(ilu_data),
                                            &jx_ParILUDataLocalResidualCorrectionRejectLikeCount(ilu_data));
      }

      /*---------------------------------------------------------------
       *    Compute residual and residual norm
       *----------------------------------------------------------------*/

      if (print_level > 1 || logging > 1 || tol > 0.)
      {
         old_resnorm = resnorm;

         if (logging > 1)
         {
            jx_ParVectorCopy(F_array, residual);
            jx_ParCSRMatrixMatvec(alpha, matA, U_array, beta, residual);
            resnorm = sqrt(jx_ParVectorInnerProd(residual, residual));
         }
         else
         {
            jx_ParVectorCopy(F_array, Ftemp);
            jx_ParCSRMatrixMatvec(alpha, matA, U_array, beta, Ftemp);
            resnorm = sqrt(jx_ParVectorInnerProd(Ftemp, Ftemp));
         }

         if (old_resnorm)
         {
            conv_factor = resnorm / old_resnorm;
         }
         else
         {
            conv_factor = resnorm;
         }

         if (rhs_norm > JX_REAL_EPSILON)
         {
            rel_resnorm = resnorm / rhs_norm;
         }
         else
         {
            rel_resnorm = resnorm;
         }

         norms[iter] = rel_resnorm;
      }

      ++iter;
      jx_ParILUDataNumIterations(ilu_data) = iter;
      jx_ParILUDataFinalRelResidualNorm(ilu_data) = rel_resnorm;

      if (my_id == 0 && print_level > 1)
      {
         jx_printf("    ILUSolve %2d   %e    %f     %e \n", iter,
                   resnorm, conv_factor, rel_resnorm);
      }
   }

   if (lrc_type == 3)
   {
      jx_ILUApplyLocalResidualCorrection(matA, F_array, U_array, Ftemp,
                                         lrc_marker,
                                         lrc_core_indices,
                                         lrc_core_size,
                                         lrc_local_indices,
                                         lrc_local_map,
                                         lrc_local_mat,
                                         lrc_local_indexD,
                                         lrc_local_indexLU,
                                         lrc_local_valueLU,
                                         lrc_local_rhs,
                                         lrc_local_sol,
                                         lrc_local_work,
                                         lrc_type, lrc_max_iter,
                                         lrc_threshold, lrc_omega,
                                         print_level,
                                         &jx_ParILUDataLocalResidualCorrectionApplyCount(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionApplyTime(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionMatvecTime(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionLocalSolveTime(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionDeltaNormRatioSum(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionResidualRatioSum(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionResidualRatioMin(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionResidualRatioMax(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionAcceptLikeCount(ilu_data),
                                         &jx_ParILUDataLocalResidualCorrectionRejectLikeCount(ilu_data));

      if (print_level > 1 || logging > 1 || tol > 0.)
      {
         if (logging > 1)
         {
            jx_ParVectorCopy(F_array, residual);
            jx_ParCSRMatrixMatvec(alpha, matA, U_array, beta, residual);
            resnorm = sqrt(jx_ParVectorInnerProd(residual, residual));
         }
         else
         {
            jx_ParVectorCopy(F_array, Ftemp);
            jx_ParCSRMatrixMatvec(alpha, matA, U_array, beta, Ftemp);
            resnorm = sqrt(jx_ParVectorInnerProd(Ftemp, Ftemp));
         }

         if (rhs_norm > JX_REAL_EPSILON)
         {
            rel_resnorm = resnorm / rhs_norm;
         }
         else
         {
            rel_resnorm = resnorm;
         }
         jx_ParILUDataFinalRelResidualNorm(ilu_data) = rel_resnorm;
      }
   }

   /* check convergence within max_iter */
   if (iter == max_iter && tol > 0. && rel_resnorm >= tol)
   {
      Solve_err_flag = 1;
      jx_error(JX_ERROR_CONV);
   }

   /*-----------------------------------------------------------------------
    *    Print closing statistics
    *    Add operator and grid complexity stats
    *-----------------------------------------------------------------------*/

   if (iter > 0 && init_resnorm)
   {
      conv_factor = pow((resnorm / init_resnorm), (1.0 / (JX_Real)iter));
   }
   else
   {
      conv_factor = 1.;
   }

   if (print_level > 1)
   {
      /*** compute operator and grid complexity (fill factor) here ?? ***/
      if (my_id == 0)
      {
         if (Solve_err_flag == 1)
         {
            jx_printf("\n\n==============================================");
            jx_printf("\n NOTE: Convergence tolerance was not achieved\n");
            jx_printf("      within the allowed %d iterations\n", max_iter);
            jx_printf("==============================================");
         }
         jx_printf("\n\n Average Convergence Factor = %f \n", conv_factor);
         jx_printf("                operator = %f\n", operat_cmplxty);
      }
   }
   return jx_error_flag;
}

/*--------------------------------------------------------------------
 * jx_ILUSolveSchurGMRES
 *
 * Schur Complement solve with GMRES on schur complement
 *
 * ParCSRMatrix S is already built in ilu data sturcture, here directly
 * use S, L, D and U factors only have local scope (no off-diag terms)
 * so apart from the residual calculation (which uses A), the solves
 * with the L and U factors are local.
 *
 * S is the global Schur complement
 * schur_solver is a GMRES solver
 * schur_precond is the ILU preconditioner for GMRES
 * rhs and x are helper vector for solving Schur system
 *--------------------------------------------------------------------*/
JX_Int
jx_ILUSolveSchurGMRES(jx_ParCSRMatrix *A,
                      jx_ParVector *f,
                      jx_ParVector *u,
                      JX_Int *perm,
                      JX_Int *qperm,
                      JX_Int nLU,
                      jx_ParCSRMatrix *L,
                      JX_Real *D,
                      jx_ParCSRMatrix *U,
                      jx_ParCSRMatrix *S,
                      jx_ParVector *ftemp,
                      jx_ParVector *utemp,
                      JX_Solver schur_solver,
                      JX_Solver schur_precond,
                      jx_ParVector *rhs,
                      jx_ParVector *x,
                      JX_Int *u_end)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUSolveSchurGMRES is called. \n");

   return jx_error_flag;
}

/*--------------------------------------------------------------------
 * jx_ILUSolveSchurNSH
 *
 * Newton-Schulz-Hotelling solve
 *
 * ParCSRMatrix S is already built in ilu data sturcture
 *
 * S here is the INVERSE of Schur Complement
 * L, D and U factors only have local scope (no off-diag terms)
 *  so apart from the residual calculation (which uses A), the solves
 *  with the L and U factors are local.
 * S is the inverse global Schur complement
 * rhs and x are helper vector for solving Schur system
 *--------------------------------------------------------------------*/
JX_Int
jx_ILUSolveSchurNSH(jx_ParCSRMatrix *A,
                    jx_ParVector *f,
                    jx_ParVector *u,
                    JX_Int *perm,
                    JX_Int nLU,
                    jx_ParCSRMatrix *L,
                    JX_Real *D,
                    jx_ParCSRMatrix *U,
                    jx_ParCSRMatrix *S,
                    jx_ParVector *ftemp,
                    jx_ParVector *utemp,
                    JX_Solver schur_solver,
                    jx_ParVector *rhs,
                    jx_ParVector *x,
                    JX_Int *u_end)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUSolveSchurNSH is called. \n");
}

/*--------------------------------------------------------------------
 * jx_ILUSolveLU
 *
 * Incomplete LU solve
 *
 * L, D and U factors only have local scope (no off-diagterms)
 *  so apart from the residual calculation (which uses A),
 *  the solves with the L and U factors are local.
 *
 * Note: perm contains the permutation of indexes corresponding to
 * user-prescribed reordering strategy. In the block Jacobi case, perm
 * may be NULL if no reordering is done (for performance, (perm == NULL)
 * assumes identity mapping of indexes). Hence we need to check the local
 * solves for this case and avoid segfaults. - DOK
 *--------------------------------------------------------------------*/
/* 串行三角求解 */
JX_Int
jx_ILUSolveLU(jx_ParCSRMatrix *A,
              jx_ParVector *f,
              jx_ParVector *u,
              JX_Int *perm,
              JX_Int nLU,
              jx_ParCSRMatrix *L,
              JX_Real *D,
              jx_ParCSRMatrix *U,
              jx_ParVector *ftemp,
              jx_ParVector *utemp)
{
   /* data objects for L and U */
   jx_CSRMatrix *L_diag = jx_ParCSRMatrixDiag(L);
   JX_Real *L_diag_data = jx_CSRMatrixData(L_diag);
   JX_Int *L_diag_i = jx_CSRMatrixI(L_diag);
   JX_Int *L_diag_j = jx_CSRMatrixJ(L_diag);
   jx_CSRMatrix *U_diag = jx_ParCSRMatrixDiag(U);
   JX_Real *U_diag_data = jx_CSRMatrixData(U_diag);
   JX_Int *U_diag_i = jx_CSRMatrixI(U_diag);
   JX_Int *U_diag_j = jx_CSRMatrixJ(U_diag);

   /* Vectors */
   jx_Vector *utemp_local = jx_ParVectorLocalVector(utemp);
   JX_Real *utemp_data = jx_VectorData(utemp_local);
   jx_Vector *ftemp_local = jx_ParVectorLocalVector(ftemp);
   JX_Real *ftemp_data = jx_VectorData(ftemp_local);
   JX_Real alpha = -1.0;
   JX_Real beta = 1.0;
   JX_Int i, j, k1, k2;

   /* Initialize Utemp to zero.
    * This is necessary for correctness, when we use optimized
    * vector operations in the case where sizeof(L, D or U) < sizeof(A)
    */
   // jx_ParVectorSetConstantValues( utemp, 0.);
   /* compute residual */
   jx_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* L solve - Forward solve */
   /* copy rhs to account for diagonal of L (which is identity) */
   if (perm)
   {
      for (i = 0; i < nLU; i++)
      {
         utemp_data[perm[i]] = ftemp_data[perm[i]];
      }
   }
   else
   {
      for (i = 0; i < nLU; i++)
      {
         utemp_data[i] = ftemp_data[i];
      }
   }

   /* Update with remaining (off-diagonal) entries of L */
   if (perm)
   {
      for (i = 0; i < nLU; i++)
      {
         k1 = L_diag_i[i];
         k2 = L_diag_i[i + 1];
         for (j = k1; j < k2; j++)
         {
            utemp_data[perm[i]] -= L_diag_data[j] * utemp_data[perm[L_diag_j[j]]];
         }
      }
   }
   else
   {
      for (i = 0; i < nLU; i++)
      {
         k1 = L_diag_i[i];
         k2 = L_diag_i[i + 1];
         for (j = k1; j < k2; j++)
         {
            utemp_data[i] -= L_diag_data[j] * utemp_data[L_diag_j[j]];
         }
      }
   }
   /*-------------------- U solve - Backward substitution */
   if (perm)
   {
      for (i = nLU - 1; i >= 0; i--)
      {
         /* first update with the remaining (off-diagonal) entries of U */
         k1 = U_diag_i[i];
         k2 = U_diag_i[i + 1];
         for (j = k1; j < k2; j++)
         {
            utemp_data[perm[i]] -= U_diag_data[j] * utemp_data[perm[U_diag_j[j]]];
         }

         /* diagonal scaling (contribution from D. Note: D is stored as its inverse) */
         utemp_data[perm[i]] *= D[i];
      }
   }
   else
   {
      for (i = nLU - 1; i >= 0; i--)
      {
         /* first update with the remaining (off-diagonal) entries of U */
         k1 = U_diag_i[i];
         k2 = U_diag_i[i + 1];
         for (j = k1; j < k2; j++)
         {
            utemp_data[i] -= U_diag_data[j] * utemp_data[U_diag_j[j]];
         }

         /* diagonal scaling (contribution from D. Note: D is stored as its inverse) */
         utemp_data[i] *= D[i];
      }
   }
   /* Update solution */
   jx_ParVectorAxpy(beta, utemp, u);

   return jx_error_flag;
}

/*--------------------------------------------------------------------
 * jx_ILUSolveLUIter
 *
 * Iterative incomplete LU solve
 *
 * L, D and U factors only have local scope (no off-diag terms)
 *  so apart from the residual calculation (which uses A), the solves
 *  with the L and U factors are local.
 *
 * Note: perm contains the permutation of indexes corresponding to
 * user-prescribed reordering strategy. In the block Jacobi case, perm
 * may be NULL if no reordering is done (for performance, (perm == NULL)
 * assumes identity mapping of indexes). Hence we need to check the local
 * solves for this case and avoid segfaults. - DOK
 *--------------------------------------------------------------------*/

JX_Int
jx_ILUSolveLUIter(jx_ParCSRMatrix *A,
                  jx_ParVector *f,
                  jx_ParVector *u,
                  JX_Int *perm,
                  JX_Int nLU,
                  jx_ParCSRMatrix *L,
                  JX_Real *D,
                  jx_ParCSRMatrix *U,
                  jx_ParVector *ftemp,
                  jx_ParVector *utemp,
                  JX_Int lower_jacobi_iters,
                  JX_Int upper_jacobi_iters)
{
   /* Data objects for L and U */
   jx_CSRMatrix *L_diag = jx_ParCSRMatrixDiag(L);
   JX_Real *L_diag_data = jx_CSRMatrixData(L_diag);
   JX_Int *L_diag_i = jx_CSRMatrixI(L_diag);
   JX_Int *L_diag_j = jx_CSRMatrixJ(L_diag);
   jx_CSRMatrix *U_diag = jx_ParCSRMatrixDiag(U);
   JX_Real *U_diag_data = jx_CSRMatrixData(U_diag);
   JX_Int *U_diag_i = jx_CSRMatrixI(U_diag);
   JX_Int *U_diag_j = jx_CSRMatrixJ(U_diag);

   /* Vectors */
   jx_Vector *utemp_local = jx_ParVectorLocalVector(utemp);
   JX_Real *utemp_data = jx_VectorData(utemp_local);
   jx_Vector *ftemp_local = jx_ParVectorLocalVector(ftemp);
   JX_Real *ftemp_data = jx_VectorData(ftemp_local);

   /* Local variables */
   JX_Real alpha = -1.0;
   JX_Real beta = 1.0;
   JX_Real sum;
   JX_Int i, j, k1, k2, kk;

   /* Initialize Utemp to zero.
    * This is necessary for correctness, when we use optimized
    * vector operations in the case where sizeof(L, D or U) < sizeof(A)
    */
   // jx_ParVectorSetConstantValues( utemp, 0.);
   /* compute residual */
   jx_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* L solve - Forward solve */
   /* copy rhs to account for diagonal of L (which is identity) */

   /* Initialize iteration to 0 */
   if (perm)
   {
      for (i = 0; i < nLU; i++)
      {
         utemp_data[perm[i]] = 0.0;
      }
   }
   else
   {
      for (i = 0; i < nLU; i++)
      {
         utemp_data[i] = 0.0;
      }
   }
   /* Jacobi iteration loop */
   for (kk = 0; kk < lower_jacobi_iters; kk++)
   {
      /* u^{k+1} = f - Lu^k */

      /* Do a SpMV with L and save the results in xtemp */
      if (perm)
      {
         for (i = nLU - 1; i >= 0; i--)
         {
            sum = 0.0;
            k1 = L_diag_i[i];
            k2 = L_diag_i[i + 1];
            for (j = k1; j < k2; j++)
            {
               sum += L_diag_data[j] * utemp_data[perm[L_diag_j[j]]];
            }
            utemp_data[perm[i]] = ftemp_data[perm[i]] - sum;
         }
      }
      else
      {
         for (i = nLU - 1; i >= 0; i--)
         {
            sum = 0.0;
            k1 = L_diag_i[i];
            k2 = L_diag_i[i + 1];
            for (j = k1; j < k2; j++)
            {
               sum += L_diag_data[j] * utemp_data[L_diag_j[j]];
            }
            utemp_data[i] = ftemp_data[i] - sum;
         }
      }
   } /* end jacobi loop */

   /* Initialize iteration to 0 */
   if (perm)
   {
      for (i = 0; i < nLU; i++)
      {
         ftemp_data[perm[i]] = 0.0;
      }
   }
   else
   {
      for (i = 0; i < nLU; i++)
      {
         ftemp_data[i] = 0.0;
      }
   }

   /* Jacobi iteration loop */
   for (kk = 0; kk < upper_jacobi_iters; kk++)
   {
      /* u^{k+1} = f - Uu^k */

      /* Do a SpMV with U and save the results in xtemp */
      if (perm)
      {
         for (i = 0; i < nLU; ++i)
         {
            sum = 0.0;
            k1 = U_diag_i[i];
            k2 = U_diag_i[i + 1];
            for (j = k1; j < k2; j++)
            {
               sum += U_diag_data[j] * ftemp_data[perm[U_diag_j[j]]];
            }
            ftemp_data[perm[i]] = D[i] * (utemp_data[perm[i]] - sum);
         }
      }
      else
      {
         for (i = 0; i < nLU; ++i)
         {
            sum = 0.0;
            k1 = U_diag_i[i];
            k2 = U_diag_i[i + 1];
            for (j = k1; j < k2; j++)
            {
               sum += U_diag_data[j] * ftemp_data[U_diag_j[j]];
            }
            ftemp_data[i] = D[i] * (utemp_data[i] - sum);
         }
      }
   } /* end jacobi loop */

   /* Update solution */
   jx_ParVectorAxpy(beta, ftemp, u);

   return jx_error_flag;
}

/*--------------------------------------------------------------------
 * jx_ILUSolveLURAS
 *
 * Incomplete LU solve RAS
 *
 * L, D and U factors only have local scope (no off-diag terms)
 *  so apart from the residual calculation (which uses A), the solves
 *  with the L and U factors are local.
 * fext and uext are tempory arrays for external data
 *--------------------------------------------------------------------*/
JX_Int
jx_ILUSolveLURAS(jx_ParCSRMatrix *A,
                 jx_ParVector *f,
                 jx_ParVector *u,
                 JX_Int *perm,
                 jx_ParCSRMatrix *L,
                 JX_Real *D,
                 jx_ParCSRMatrix *U,
                 jx_ParVector *ftemp,
                 jx_ParVector *utemp,
                 JX_Real *fext,
                 JX_Real *uext)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUSolveLURAS is called. \n");

   return jx_error_flag;
}

/******************************************************************************
 *
 * NSH functions.
 *
 *****************************************************************************/

/*--------------------------------------------------------------------
 * jx_NSHSolve
 *--------------------------------------------------------------------*/

JX_Int
jx_NSHSolve(void *nsh_vdata,
            jx_ParCSRMatrix *A,
            jx_ParVector *f,
            jx_ParVector *u)
{
   jx_printf("[JX-ILU] WARNING: jx_NSHSolve is called. \n");

   return jx_error_flag;
}

/*--------------------------------------------------------------------
 * jx_NSHSolveInverse
 *
 * Simply a matvec on residual with approximate inverse
 *
 * A: original matrix
 * f: rhs
 * u: solution
 * M: approximate inverse
 * ftemp, utemp: working vectors
 *--------------------------------------------------------------------*/

JX_Int
jx_NSHSolveInverse(jx_ParCSRMatrix *A,
                   jx_ParVector *f,
                   jx_ParVector *u,
                   jx_ParCSRMatrix *M,
                   jx_ParVector *ftemp,
                   jx_ParVector *utemp)
{
   jx_printf("[JX-ILU] WARNING: jx_NSHSolveInverse is called. \n");
   
   return jx_error_flag;
}
