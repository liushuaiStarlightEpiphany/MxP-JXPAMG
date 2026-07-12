//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu_solve.c
 *
 */

#include "jxf_ilu.h"
#include "jxf_krylov.h"

JXF_Int
jxf_ILUSolve( void *ilu_vdata, jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u )
{
   MPI_Comm             comm           = jxf_ParCSRMatrixComm(A);
   //   JXF_Int            i;

   jxf_ParILUData     *ilu_data      = (jxf_ParILUData*) ilu_vdata;

   /* get matrices */
   JXF_Int            ilu_type       = jxf_ParILUDataIluType(ilu_data);
   JXF_Int            *perm          = jxf_ParILUDataPerm(ilu_data);
   JXF_Int            *qperm         = jxf_ParILUDataQPerm(ilu_data);
   jxf_ParCSRMatrix   *matA          = jxf_ParILUDataMatA(ilu_data);
   jxf_ParCSRMatrix   *matL          = jxf_ParILUDataMatL(ilu_data);
   JXF_Real           *matD          = jxf_ParILUDataMatD(ilu_data);
   jxf_ParCSRMatrix   *matU          = jxf_ParILUDataMatU(ilu_data);
   jxf_ParCSRMatrix   *matS          = jxf_ParILUDataMatS(ilu_data);

   JXF_Int            iter, num_procs,  my_id;

   jxf_ParVector      *F_array       = jxf_ParILUDataF(ilu_data);
   jxf_ParVector      *U_array       = jxf_ParILUDataU(ilu_data);

   /* get settings */
   JXF_Real           tol            = jxf_ParILUDataTol(ilu_data);
   JXF_Int            logging        = jxf_ParILUDataLogging(ilu_data);
   JXF_Int            print_level    = jxf_ParILUDataPrintLevel(ilu_data);
   JXF_Int            max_iter       = jxf_ParILUDataMaxIter(ilu_data);
   JXF_Real           *norms         = jxf_ParILUDataRelResNorms(ilu_data);
   jxf_ParVector      *Ftemp         = jxf_ParILUDataFTemp(ilu_data);
   jxf_ParVector      *Utemp         = jxf_ParILUDataUTemp(ilu_data);
   JXF_Real           *fext          = jxf_ParILUDataFExt(ilu_data);
   JXF_Real           *uext          = jxf_ParILUDataUExt(ilu_data);
   jxf_ParVector      *residual = NULL;

   JXF_Real           alpha          = -1;
   JXF_Real           beta           = 1;
   JXF_Real           conv_factor    = 0.0;
   JXF_Real           resnorm        = 1.0;
   JXF_Real           init_resnorm   = 0.0;
   JXF_Real           rel_resnorm;
   JXF_Real           rhs_norm       = 0.0;
   JXF_Real           old_resnorm;
   JXF_Real           ieee_check     = 0.0;
   JXF_Real           operat_cmplxty = (ilu_data -> operator_complexity);

   JXF_Int            Solve_err_flag;

   /* problem size */
   JXF_Int            n              = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A));
   JXF_Int            nLU            = jxf_ParILUDataNLU(ilu_data);
   JXF_Int            *u_end         = jxf_ParILUDataUEnd(ilu_data);

   /* Schur system solve */
   JXF_Solver         schur_solver   = jxf_ParILUDataSchurSolver(ilu_data);
   JXF_Solver         schur_precond  = jxf_ParILUDataSchurPrecond(ilu_data);
   jxf_ParVector      *rhs           = jxf_ParILUDataRhs(ilu_data);
   jxf_ParVector      *x             = jxf_ParILUDataX(ilu_data);

   /* begin */

   if(logging > 1)
   {
      residual = (ilu_data -> residual);
   }

   (ilu_data -> num_iterations) = 0;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);

   /*-----------------------------------------------------------------------
    *    Write the solver parameters
    *-----------------------------------------------------------------------*/
   if (my_id == 0 && print_level > 1)
   {
      jxf_ILUWriteSolverParams(ilu_data);
   }

   /*-----------------------------------------------------------------------
    *    Initialize the solver error flag
    *-----------------------------------------------------------------------*/

   Solve_err_flag = 0;
   /*-----------------------------------------------------------------------
    *     write some initial info
    *-----------------------------------------------------------------------*/

   if (my_id == 0 && print_level > 1 && tol > 0.)
   {
      jxf_printf("\n\n ILU SOLVER SOLUTION INFO:\n");
   }


   /*-----------------------------------------------------------------------
    *    Compute initial residual and print
    *-----------------------------------------------------------------------*/
   if (print_level > 1 || logging > 1 || tol > 0.)
   {
      if ( logging > 1 )
      {
         jxf_ParVectorCopy(f, residual );
         if (tol > 0.0)
         {
            jxf_ParCSRMatrixMatvec(alpha, A, u, beta, residual );
         }
         resnorm = sqrt(jxf_ParVectorInnerProd( residual, residual ));
      }
      else
      {
         jxf_ParVectorCopy(f, Ftemp);
         if (tol > 0.0)
         {
            jxf_ParCSRMatrixMatvec(alpha, A, u, beta, Ftemp);
         }
         resnorm = sqrt(jxf_ParVectorInnerProd(Ftemp, Ftemp));
      }

      /* Since it is does not diminish performance, attempt to return an error flag
         and notify users when they supply bad input. */
      if (resnorm != 0.)
      {
         ieee_check = resnorm/resnorm; /* INF -> NaN conversion */
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
            jxf_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
            jxf_printf("ERROR -- jxf_ILUSolve: INFs and/or NaNs detected in input.\n");
            jxf_printf("User probably placed non-numerics in supplied A, x_0, or b.\n");
            jxf_printf("ERROR detected by Hypre ...  END\n\n\n");
         }
         jxf_error(JXF_ERROR_GENERIC);
         return jxf_error_flag;
      }

      init_resnorm = resnorm;
      rhs_norm = sqrt(jxf_ParVectorInnerProd(f, f));
      if (rhs_norm > JXF_REAL_EPSILON)
      {
         rel_resnorm = init_resnorm / rhs_norm;
      }
      else
      {
         /* rhs is zero, return a zero solution */
         jxf_ParVectorSetConstantValues(U_array, 0.0);
         if(logging > 0)
         {
            rel_resnorm = 0.0;
            (ilu_data -> final_rel_residual_norm) = rel_resnorm;
         }
         return jxf_error_flag;
      }
   }
   else
   {
      rel_resnorm = 1.;
   }

   if (my_id == 0 && print_level > 1)
   {
      jxf_printf("                                            relative\n");
      jxf_printf("               residual        factor       residual\n");
      jxf_printf("               --------        ------       --------\n");
      jxf_printf("    Initial    %e                 %e\n",init_resnorm,
            rel_resnorm);
   }

   matA = A;
   U_array = u;
   F_array = f;

   /************** Main Solver Loop - always do 1 iteration ************/
   iter = 0;

   while ((rel_resnorm >= tol || iter < 1)
         && iter < max_iter)
   {

      /* Do one solve on LUe=r */
      switch(ilu_type){
         case 0: case 1:
            jxf_ILUSolveLU(matA, F_array, U_array, perm, n, matL, matD, matU, Utemp, Ftemp); //BJ
            break;
         case 10: case 11:
            jxf_ILUSolveSchurGMRES(matA, F_array, U_array, perm, perm, nLU, matL, matD, matU, matS,
                  Utemp, Ftemp, schur_solver, schur_precond, rhs, x, u_end); //GMRES
            break;
         case 20: case 21:
            jxf_ILUSolveSchurNSH(matA, F_array, U_array, perm, nLU, matL, matD, matU, matS,
                  Utemp, Ftemp, schur_solver, rhs, x, u_end); //MR+NSH
            break;
         case 30: case 31:
            jxf_ILUSolveLURAS(matA, F_array, U_array, perm, matL, matD, matU, Utemp, Utemp, fext, uext); //RAS
            break;
         case 40: case 41:
            jxf_ILUSolveSchurGMRES(matA, F_array, U_array, perm, qperm, nLU, matL, matD, matU, matS,
                  Utemp, Ftemp, schur_solver, schur_precond, rhs, x, u_end); //GMRES
            break;
         default:
            jxf_ILUSolveLU(matA, F_array, U_array, perm, n, matL, matD, matU, Utemp, Ftemp); //BJ
            break;

      }

      /*---------------------------------------------------------------
       *    Compute residual and residual norm
       *----------------------------------------------------------------*/

      if (print_level > 1 || logging > 1 || tol > 0.)
      {
         old_resnorm = resnorm;

         if ( logging > 1 ) {
            jxf_ParVectorCopy(F_array, residual);
            jxf_ParCSRMatrixMatvec(alpha, matA, U_array, beta, residual );
            resnorm = sqrt(jxf_ParVectorInnerProd( residual, residual ));
         }
         else {
            jxf_ParVectorCopy(F_array, Ftemp);
            jxf_ParCSRMatrixMatvec(alpha, matA, U_array, beta, Ftemp);
            resnorm = sqrt(jxf_ParVectorInnerProd(Ftemp, Ftemp));
         }

         if (old_resnorm) conv_factor = resnorm / old_resnorm;
         else conv_factor = resnorm;
         if (rhs_norm > JXF_REAL_EPSILON)
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
      (ilu_data -> num_iterations) = iter;
      (ilu_data -> final_rel_residual_norm) = rel_resnorm;

      if (my_id == 0 && print_level > 1)
      {
         jxf_printf("    ILUSolve %2d   %e    %f     %e \n", iter,
               resnorm, conv_factor, rel_resnorm);
      }
   }

   /* check convergence within max_iter */
   if (iter == max_iter && tol > 0.)
   {
      Solve_err_flag = 1;
      jxf_error(JXF_ERROR_CONV);
   }

   /*-----------------------------------------------------------------------
    *    Print closing statistics
    *    Add operator and grid complexity stats
    *-----------------------------------------------------------------------*/

   if (iter > 0 && init_resnorm)
   {
      conv_factor = pow((resnorm/init_resnorm),(1.0/(JXF_Real) iter));
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
            jxf_printf("\n\n==============================================");
            jxf_printf("\n NOTE: Convergence tolerance was not achieved\n");
            jxf_printf("      within the allowed %d iterations\n",max_iter);
            jxf_printf("==============================================");
         }
         jxf_printf("\n\n Average Convergence Factor = %f \n",conv_factor);
         jxf_printf("                operator = %f\n",operat_cmplxty);
      }
   }
   return jxf_error_flag;
}

/* Schur Complement solve with GMRES on schur complement
 * ParCSRMatrix S is already built in ilu data sturcture, here directly use S
 * L, D and U factors only have local scope (no off-diagonal processor terms)
 * so apart from the residual calculation (which uses A), the solves with the
 * L and U factors are local.
 * S is the global Schur complement
 * schur_solver is a GMRES solver
 * schur_precond is the ILU preconditioner for GMRES
 * rhs and x are helper vector for solving Schur system
*/

JXF_Int
jxf_ILUSolveSchurGMRES(jxf_ParCSRMatrix *A, jxf_ParVector    *f,
                  jxf_ParVector    *u, JXF_Int *perm, JXF_Int *qperm,
                  JXF_Int nLU, jxf_ParCSRMatrix *L,
                  JXF_Real* D, jxf_ParCSRMatrix *U,
                  jxf_ParCSRMatrix *S,
                  jxf_ParVector *ftemp, jxf_ParVector *utemp,
                  JXF_Solver schur_solver, JXF_Solver schur_precond,
                  jxf_ParVector *rhs, jxf_ParVector *x, JXF_Int *u_end)
{
   /* data objects for communication */
   //   MPI_Comm          comm = jxf_ParCSRMatrixComm(A);

   /* data objects for L and U */
   jxf_CSRMatrix   *L_diag = jxf_ParCSRMatrixDiag(L);
   JXF_Real        *L_diag_data = jxf_CSRMatrixData(L_diag);
   JXF_Int         *L_diag_i = jxf_CSRMatrixI(L_diag);
   JXF_Int         *L_diag_j = jxf_CSRMatrixJ(L_diag);
   jxf_CSRMatrix   *U_diag = jxf_ParCSRMatrixDiag(U);
   JXF_Real        *U_diag_data = jxf_CSRMatrixData(U_diag);
   JXF_Int         *U_diag_i = jxf_CSRMatrixI(U_diag);
   JXF_Int         *U_diag_j = jxf_CSRMatrixJ(U_diag);
   jxf_Vector      *utemp_local = jxf_ParVectorLocalVector(utemp);
   JXF_Real        *utemp_data  = jxf_VectorData(utemp_local);
   jxf_Vector      *ftemp_local = jxf_ParVectorLocalVector(ftemp);
   JXF_Real        *ftemp_data  = jxf_VectorData(ftemp_local);

   JXF_Real        alpha;
   JXF_Real        beta;
   JXF_Int         i, j, k1, k2, col;

   /* problem size */
   JXF_Int         n = jxf_CSRMatrixNumRows(L_diag);
   //   JXF_Int         m = n - nLU;

   /* other data objects for computation */
   //   jxf_Vector      *f_local;
   //   JXF_Real        *f_data;
   jxf_Vector      *rhs_local;
   JXF_Real        *rhs_data;
   jxf_Vector      *x_local;
   JXF_Real        *x_data;

   /* begin */
   beta = 1.0;
   alpha = -1.0;

   /* compute residual */
   jxf_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* 1st need to solve LBi*xi = fi
    * L solve, solve xi put in u_temp upper
    */
   //   f_local = jxf_ParVectorLocalVector(f);
   //   f_data = jxf_VectorData(f_local);
   /* now update with L to solve */
   for(i = 0 ; i < nLU ; i ++)
   {
      utemp_data[qperm[i]] = ftemp_data[perm[i]];
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         utemp_data[qperm[i]] -= L_diag_data[j] * utemp_data[qperm[L_diag_j[j]]];
      }
   }

   /* 2nd need to compute g'i = gi - Ei*UBi^-1*xi
    * now put g'i into the f_temp lower
    */
   for(i = nLU ; i < n ; i ++)
   {
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         col = L_diag_j[j];
         ftemp_data[perm[i]] -= L_diag_data[j] * utemp_data[qperm[col]];
      }
   }

   /* 3rd need to solve global Schur Complement Sy = g'
    * for now only solve the local system
    * solve y put in u_temp lower
    * only solve whe S is not NULL
    */
   if(S)
   {
      /*initialize solution to zero for residual equation */
      jxf_ParVectorSetConstantValues(x, 0.0);
      /* setup vectors for solve */
      rhs_local   = jxf_ParVectorLocalVector(rhs);
      rhs_data    = jxf_VectorData(rhs_local);
      x_local     = jxf_ParVectorLocalVector(x);
      x_data      = jxf_VectorData(x_local);

      /* set rhs value */
      for(i = nLU ; i < n ; i ++)
      {
         rhs_data[i-nLU] = ftemp_data[perm[i]];
      }

      /* solve */
      JXF_GMRESSolve(schur_solver,(JXF_Matrix)S,(JXF_Matrix)S,(JXF_Vector)rhs,(JXF_Vector)x);

      /* copy value back to original */
      for(i = nLU ; i < n ; i ++)
      {
         utemp_data[qperm[i]] = x_data[i-nLU];
      }
   }

   /* 4th need to compute zi = xi - LBi^-1*yi
    * put zi in f_temp upper
    * only do this computation when nLU < n
    * U is unsorted, search is expensive when unnecessary
    */
   if(nLU < n)
   {
      for(i = 0 ; i < nLU ; i ++)
      {
         ftemp_data[perm[i]] = utemp_data[qperm[i]];
         k1 = u_end[i] ; k2 = U_diag_i[i+1];
         for(j = k1 ; j < k2 ; j ++)
         {
            col = U_diag_j[j];
            ftemp_data[perm[i]] -= U_diag_data[j] * utemp_data[qperm[col]];
         }
      }
      for(i = 0 ; i < nLU ; i ++)
      {
         utemp_data[qperm[i]] = ftemp_data[perm[i]];
      }
   }

   /* 5th need to solve UBi*ui = zi */
   /* put result in u_temp upper */
   for(i = nLU-1 ; i >= 0 ; i --)
   {
      k1 = U_diag_i[i] ; k2 = u_end[i];
      for(j = k1 ; j < k2 ; j ++)
      {
         col = U_diag_j[j];
         utemp_data[qperm[i]] -= U_diag_data[j] * utemp_data[qperm[col]];
      }
      utemp_data[qperm[i]] *= D[i];
   }

   /* done, now everything are in u_temp, update solution */
   jxf_ParVectorAxpy(beta, utemp, u);

   return jxf_error_flag;
}

/* Newton-Schulz-Hotelling solve
 * ParCSRMatrix S is already built in ilu data sturcture
 * S here is the INVERSE of Schur Complement
 * L, D and U factors only have local scope (no off-diagonal processor terms)
 * so apart from the residual calculation (which uses A), the solves with the
 * L and U factors are local.
 * S is the inverse global Schur complement
 * rhs and x are helper vector for solving Schur system
*/

JXF_Int
jxf_ILUSolveSchurNSH(jxf_ParCSRMatrix *A, jxf_ParVector    *f,
                  jxf_ParVector    *u, JXF_Int *perm,
                  JXF_Int nLU, jxf_ParCSRMatrix *L,
                  JXF_Real* D, jxf_ParCSRMatrix *U,
                  jxf_ParCSRMatrix *S,
                  jxf_ParVector *ftemp, jxf_ParVector *utemp,
                  JXF_Solver schur_solver,
                  jxf_ParVector *rhs, jxf_ParVector *x, JXF_Int *u_end)
{
   /* data objects for communication */
   //   MPI_Comm          comm = jxf_ParCSRMatrixComm(A);

   /* data objects for L and U */
   jxf_CSRMatrix   *L_diag = jxf_ParCSRMatrixDiag(L);
   JXF_Real        *L_diag_data = jxf_CSRMatrixData(L_diag);
   JXF_Int         *L_diag_i = jxf_CSRMatrixI(L_diag);
   JXF_Int         *L_diag_j = jxf_CSRMatrixJ(L_diag);
   jxf_CSRMatrix   *U_diag = jxf_ParCSRMatrixDiag(U);
   JXF_Real        *U_diag_data = jxf_CSRMatrixData(U_diag);
   JXF_Int         *U_diag_i = jxf_CSRMatrixI(U_diag);
   JXF_Int         *U_diag_j = jxf_CSRMatrixJ(U_diag);
   jxf_Vector      *utemp_local = jxf_ParVectorLocalVector(utemp);
   JXF_Real        *utemp_data  = jxf_VectorData(utemp_local);
   jxf_Vector      *ftemp_local = jxf_ParVectorLocalVector(ftemp);
   JXF_Real        *ftemp_data  = jxf_VectorData(ftemp_local);

   JXF_Real        alpha;
   JXF_Real        beta;
   JXF_Int         i, j, k1, k2, col;

   /* problem size */
   JXF_Int         n = jxf_CSRMatrixNumRows(L_diag);
   //   JXF_Int         m = n - nLU;

   /* other data objects for computation */
   //   jxf_Vector      *f_local;
   //   JXF_Real        *f_data;
   jxf_Vector      *rhs_local;
   JXF_Real        *rhs_data;
   jxf_Vector      *x_local;
   JXF_Real        *x_data;

   /* begin */
   beta = 1.0;
   alpha = -1.0;

   /* compute residual */
   jxf_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* 1st need to solve LBi*xi = fi
    * L solve, solve xi put in u_temp upper
    */
   //   f_local = jxf_ParVectorLocalVector(f);
   //   f_data = jxf_VectorData(f_local);
   /* now update with L to solve */
   for(i = 0 ; i < nLU ; i ++)
   {
      utemp_data[perm[i]] = ftemp_data[perm[i]];
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         utemp_data[perm[i]] -= L_diag_data[j] * utemp_data[perm[L_diag_j[j]]];
      }
   }

   /* 2nd need to compute g'i = gi - Ei*UBi^-1*xi
    * now put g'i into the f_temp lower
    */
   for(i = nLU ; i < n ; i ++)
   {
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         col = L_diag_j[j];
         ftemp_data[perm[i]] -= L_diag_data[j] * utemp_data[perm[col]];
      }
   }

   /* 3rd need to solve global Schur Complement Sy = g'
    * for now only solve the local system
    * solve y put in u_temp lower
    * only solve when S is not NULL
    */
   if(S)
   {
      /*initialize solution to zero for residual equation */
      jxf_ParVectorSetConstantValues(x, 0.0);
      /* setup vectors for solve */
      rhs_local   = jxf_ParVectorLocalVector(rhs);
      rhs_data    = jxf_VectorData(rhs_local);
      x_local     = jxf_ParVectorLocalVector(x);
      x_data      = jxf_VectorData(x_local);

      /* set rhs value */
      for(i = nLU ; i < n ; i ++)
      {
         rhs_data[i-nLU] = ftemp_data[perm[i]];
      }

      /* Solve Schur system with approx inverse
       * x = S*rhs
       */
      jxf_NSHSolve(schur_solver,S,rhs,x);

      /* copy value back to original */
      for(i = nLU ; i < n ; i ++)
      {
         utemp_data[perm[i]] = x_data[i-nLU];
      }
   }

   /* 4th need to compute zi = xi - LBi^-1*yi
    * put zi in f_temp upper
    * only do this computation when nLU < n
    * U is unsorted, search is expensive when unnecessary
    */
   if(nLU < n)
   {
      for(i = 0 ; i < nLU ; i ++)
      {
         ftemp_data[perm[i]] = utemp_data[perm[i]];
         k1 = u_end[i] ; k2 = U_diag_i[i+1];
         for(j = k1 ; j < k2 ; j ++)
         {
            col = U_diag_j[j];
            ftemp_data[perm[i]] -= U_diag_data[j] * utemp_data[perm[col]];
         }
      }
      for(i = 0 ; i < nLU ; i ++)
      {
         utemp_data[perm[i]] = ftemp_data[perm[i]];
      }
   }

   /* 5th need to solve UBi*ui = zi */
   /* put result in u_temp upper */
   for(i = nLU-1 ; i >= 0 ; i --)
   {
      k1 = U_diag_i[i] ; k2 = u_end[i];
      for(j = k1 ; j < k2 ; j ++)
      {
         col = U_diag_j[j];
         utemp_data[perm[i]] -= U_diag_data[j] * utemp_data[perm[col]];
      }
      utemp_data[perm[i]] *= D[i];
   }

   /* done, now everything are in u_temp, update solution */
   jxf_ParVectorAxpy(beta, utemp, u);

   return jxf_error_flag;
}

/* Incomplete LU solve
 * L, D and U factors only have local scope (no off-diagonal processor terms)
 * so apart from the residual calculation (which uses A), the solves with the
 * L and U factors are local.
*/

JXF_Int
jxf_ILUSolveLU(jxf_ParCSRMatrix *A, jxf_ParVector    *f,
                  jxf_ParVector    *u, JXF_Int *perm,
                  JXF_Int nLU, jxf_ParCSRMatrix *L,
                  JXF_Real* D, jxf_ParCSRMatrix *U,
                  jxf_ParVector *ftemp, jxf_ParVector *utemp)
{
   jxf_CSRMatrix *L_diag = jxf_ParCSRMatrixDiag(L);
   JXF_Real      *L_diag_data = jxf_CSRMatrixData(L_diag);
   JXF_Int       *L_diag_i = jxf_CSRMatrixI(L_diag);
   JXF_Int       *L_diag_j = jxf_CSRMatrixJ(L_diag);

   jxf_CSRMatrix *U_diag = jxf_ParCSRMatrixDiag(U);
   JXF_Real      *U_diag_data = jxf_CSRMatrixData(U_diag);
   JXF_Int       *U_diag_i = jxf_CSRMatrixI(U_diag);
   JXF_Int       *U_diag_j = jxf_CSRMatrixJ(U_diag);

   jxf_Vector    *utemp_local = jxf_ParVectorLocalVector(utemp);
   JXF_Real      *utemp_data  = jxf_VectorData(utemp_local);

   jxf_Vector    *ftemp_local = jxf_ParVectorLocalVector(ftemp);
   JXF_Real      *ftemp_data  = jxf_VectorData(ftemp_local);

   JXF_Real      alpha;
   JXF_Real      beta;
   JXF_Int       i, j, k1, k2;

   /* begin */
   alpha = -1.0;
   beta = 1.0;

   /* Initialize Utemp to zero.
    * This is necessary for correctness, when we use optimized
    * vector operations in the case where sizeof(L, D or U) < sizeof(A)
    */
   //jxf_ParVectorSetConstantValues( utemp, 0.);
   /* compute residual */
   jxf_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* L solve - Forward solve */
   /* copy rhs to account for diagonal of L (which is identity) */
   for( i = 0; i < nLU; i++ )
   {
      utemp_data[perm[i]] = ftemp_data[perm[i]];
   }
   /* update with remaining (off-diagonal) entries of L */
   for( i = 0; i < nLU; i++ )
   {
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      for(j=k1; j <k2; j++)
      {
         utemp_data[perm[i]] -= L_diag_data[j] * utemp_data[perm[L_diag_j[j]]];
      }
   }
   /*-------------------- U solve - Backward substitution */
   for( i = nLU-1; i >= 0; i-- )
   {
      /* first update with the remaining (off-diagonal) entries of U */
      k1 = U_diag_i[i] ; k2 = U_diag_i[i+1];
      for(j=k1; j <k2; j++)
      {
         utemp_data[perm[i]] -= U_diag_data[j] * utemp_data[perm[U_diag_j[j]]];
      }
      /* diagonal scaling (contribution from D. Note: D is stored as its inverse) */
      utemp_data[perm[i]] *= D[i];
   }

   /* Update solution */
   jxf_ParVectorAxpy(beta, utemp, u);


   return jxf_error_flag;
}


/* Incomplete LU solve RAS
 * L, D and U factors only have local scope (no off-diagonal processor terms)
 * so apart from the residual calculation (which uses A), the solves with the
 * L and U factors are local.
 * fext and uext are tempory arrays for external data
*/

JXF_Int
jxf_ILUSolveLURAS(jxf_ParCSRMatrix *A, jxf_ParVector    *f,
                  jxf_ParVector    *u, JXF_Int *perm,
                  jxf_ParCSRMatrix *L,
                  JXF_Real* D, jxf_ParCSRMatrix *U,
                  jxf_ParVector *ftemp, jxf_ParVector *utemp,
                  JXF_Real *fext, JXF_Real *uext)
{

   jxf_ParCSRCommPkg        *comm_pkg;
   jxf_ParCSRCommHandle     *comm_handle;
   JXF_Int                  num_sends, begin, end;

   jxf_CSRMatrix            *L_diag = jxf_ParCSRMatrixDiag(L);
   JXF_Real                 *L_diag_data = jxf_CSRMatrixData(L_diag);
   JXF_Int                  *L_diag_i = jxf_CSRMatrixI(L_diag);
   JXF_Int                  *L_diag_j = jxf_CSRMatrixJ(L_diag);

   jxf_CSRMatrix            *U_diag = jxf_ParCSRMatrixDiag(U);
   JXF_Real                 *U_diag_data = jxf_CSRMatrixData(U_diag);
   JXF_Int                  *U_diag_i = jxf_CSRMatrixI(U_diag);
   JXF_Int                  *U_diag_j = jxf_CSRMatrixJ(U_diag);

   JXF_Int                  n = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixDiag(A));
   JXF_Int                  m = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(A));
   //   JXF_Int                  buffer_size;
   JXF_Int                  n_total = m + n;

   JXF_Int                  idx;
   JXF_Int                  jcol;
   JXF_Int                  col;

   jxf_Vector               *utemp_local = jxf_ParVectorLocalVector(utemp);
   JXF_Real                 *utemp_data  = jxf_VectorData(utemp_local);

   jxf_Vector               *ftemp_local = jxf_ParVectorLocalVector(ftemp);
   JXF_Real                 *ftemp_data  = jxf_VectorData(ftemp_local);

   JXF_Real                 alpha;
   JXF_Real                 beta;
   JXF_Int                  i, j, k1, k2;

   /* begin */
   alpha = -1.0;
   beta = 1.0;

   /* prepare for communication */
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   /* setup if not yet built */
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* Initialize Utemp to zero.
    * This is necessary for correctness, when we use optimized
    * vector operations in the case where sizeof(L, D or U) < sizeof(A)
    */
   //jxf_ParVectorSetConstantValues( utemp, 0.);
   /* compute residual */
   jxf_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);

   /* communication to get external data */

   /* get total num of send */
   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   begin = jxf_ParCSRCommPkgSendMapStart(comm_pkg,0);
   end = jxf_ParCSRCommPkgSendMapStart(comm_pkg,num_sends);

   /* copy new index into send_buf */
   for(i = begin ; i < end ; i ++)
   {
      /* all we need is just send out data, we don't need to worry about the
       *    permutation of offd part, actually we don't need to worry about
       *    permutation at all
       * borrow uext as send buffer .
       */
      uext[i-begin] = ftemp_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)];
   }

   /* main communication */
   comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, uext, fext);
   jxf_ParCSRCommHandleDestroy(comm_handle);

   /* L solve - Forward solve */
   for( i = 0 ; i < n_total ; i ++)
   {
      k1 = L_diag_i[i] ; k2 = L_diag_i[i+1];
      if( i < n )
      {
         /* diag part */
         utemp_data[perm[i]] = ftemp_data[perm[i]];
         for(j=k1; j <k2; j++)
         {
            col = L_diag_j[j];
            if( col < n )
            {
               utemp_data[perm[i]] -= L_diag_data[j] * utemp_data[perm[col]];
            }
            else
            {
               jcol = col - n;
               utemp_data[perm[i]] -= L_diag_data[j] * uext[jcol];
            }
         }
      }
      else
      {
         /* offd part */
         idx = i - n;
         uext[idx] = fext[idx];
         for(j=k1; j <k2; j++)
         {
            col = L_diag_j[j];
            if(col < n)
            {
               uext[idx] -= L_diag_data[j] * utemp_data[perm[col]];
            }
            else
            {
               jcol = col - n;
               uext[idx] -= L_diag_data[j] * uext[jcol];
            }
         }
      }
   }

   /*-------------------- U solve - Backward substitution */
   for( i = n_total-1; i >= 0; i-- )
   {
      /* first update with the remaining (off-diagonal) entries of U */
      k1 = U_diag_i[i] ; k2 = U_diag_i[i+1];
      if( i < n )
      {
         /* diag part */
         for(j=k1; j <k2; j++)
         {
            col = U_diag_j[j];
            if( col < n )
            {
               utemp_data[perm[i]] -= U_diag_data[j] * utemp_data[perm[col]];
            }
            else
            {
               jcol = col - n;
               utemp_data[perm[i]] -= U_diag_data[j] * uext[jcol];
            }
         }
         /* diagonal scaling (contribution from D. Note: D is stored as its inverse) */
         utemp_data[perm[i]] *= D[i];
      }
      else
      {
         /* 2nd part of offd */
         idx = i - n;
         for(j=k1; j <k2; j++)
         {
            col = U_diag_j[j];
            if( col < n )
            {
               uext[idx] -= U_diag_data[j] * utemp_data[perm[col]];
            }
            else
            {
               jcol = col - n;
               uext[idx] -= U_diag_data[j] * uext[jcol];
            }
         }
         /* diagonal scaling (contribution from D. Note: D is stored as its inverse) */
         uext[idx] *= D[i];
      }

   }
   /* Update solution */
   jxf_ParVectorAxpy(beta, utemp, u);

   return jxf_error_flag;
}


/* solve functions for NSH */

JXF_Int
jxf_NSHSolve( void               *nsh_vdata,
                  jxf_ParCSRMatrix *A,
                  jxf_ParVector    *f,
                  jxf_ParVector    *u )
{
   MPI_Comm             comm = jxf_ParCSRMatrixComm(A);
   //   JXF_Int            i;

   jxf_ParNSHData     *nsh_data = (jxf_ParNSHData*) nsh_vdata;

   /* get matrices */
   jxf_ParCSRMatrix   *matA = jxf_ParNSHDataMatA(nsh_data);
   jxf_ParCSRMatrix   *matM = jxf_ParNSHDataMatM(nsh_data);

   JXF_Int            iter, num_procs,  my_id;

   jxf_ParVector      *F_array = jxf_ParNSHDataF(nsh_data);
   jxf_ParVector      *U_array = jxf_ParNSHDataU(nsh_data);

   /* get settings */
   JXF_Real           tol = jxf_ParNSHDataTol(nsh_data);
   JXF_Int            logging = jxf_ParNSHDataLogging(nsh_data);
   JXF_Int            print_level = jxf_ParNSHDataPrintLevel(nsh_data);
   JXF_Int            max_iter = jxf_ParNSHDataMaxIter(nsh_data);
   JXF_Real           *norms = jxf_ParNSHDataRelResNorms(nsh_data);
   jxf_ParVector      *Ftemp = jxf_ParNSHDataFTemp(nsh_data);
   jxf_ParVector      *Utemp = jxf_ParNSHDataUTemp(nsh_data);
   jxf_ParVector      *residual = NULL;

   JXF_Real           alpha = -1.0;
   JXF_Real           beta = 1.0;
   JXF_Real           conv_factor = 0.0;
   JXF_Real           resnorm = 1.0;
   JXF_Real           init_resnorm = 0.0;
   JXF_Real           rel_resnorm;
   JXF_Real           rhs_norm = 0.0;
   JXF_Real           old_resnorm;
   JXF_Real           ieee_check = 0.;
   JXF_Real           operat_cmplxty = jxf_ParNSHDataOperatorComplexity(nsh_data);

   JXF_Int            Solve_err_flag;

   /* problem size */
   //   JXF_Int            n = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A));

   /* begin */
   if(logging > 1)
   {
      residual = jxf_ParNSHDataResidual(nsh_data);
   }

   jxf_ParNSHDataNumIterations(nsh_data) = 0;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);

   /*-----------------------------------------------------------------------
    *    Write the solver parameters
    *-----------------------------------------------------------------------*/
   if (my_id == 0 && print_level > 1)
   {
      jxf_NSHWriteSolverParams(nsh_data);
   }

   /*-----------------------------------------------------------------------
    *    Initialize the solver error flag
    *-----------------------------------------------------------------------*/

   Solve_err_flag = 0;
   /*-----------------------------------------------------------------------
    *     write some initial info
    *-----------------------------------------------------------------------*/

   if (my_id == 0 && print_level > 1 && tol > 0.)
   {
      jxf_printf("\n\n Newton–Schulz–Hotelling SOLVER SOLUTION INFO:\n");
   }


   /*-----------------------------------------------------------------------
    *    Compute initial residual and print
    *-----------------------------------------------------------------------*/
   if (print_level > 1 || logging > 1 || tol > 0.)
   {
      if ( logging > 1 )
      {
         jxf_ParVectorCopy(f, residual );
         if (tol > 0.0)
         {
            jxf_ParCSRMatrixMatvec(alpha, A, u, beta, residual );
         }
         resnorm = sqrt(jxf_ParVectorInnerProd( residual, residual ));
      }
      else
      {
         jxf_ParVectorCopy(f, Ftemp);
         if (tol > 0.0)
         {
            jxf_ParCSRMatrixMatvec(alpha, A, u, beta, Ftemp);
         }
         resnorm = sqrt(jxf_ParVectorInnerProd(Ftemp, Ftemp));
      }

      /* Since it is does not diminish performance, attempt to return an error flag
         and notify users when they supply bad input. */
      if (resnorm != 0.)
      {
         ieee_check = resnorm/resnorm; /* INF -> NaN conversion */
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
            jxf_printf("\n\nERROR detected by Hypre ...  BEGIN\n");
            jxf_printf("ERROR -- jxf_NSHSolve: INFs and/or NaNs detected in input.\n");
            jxf_printf("User probably placed non-numerics in supplied A, x_0, or b.\n");
            jxf_printf("ERROR detected by Hypre ...  END\n\n\n");
         }
         jxf_error(JXF_ERROR_GENERIC);
         return jxf_error_flag;
      }

      init_resnorm = resnorm;
      rhs_norm = sqrt(jxf_ParVectorInnerProd(f, f));
      if (rhs_norm > JXF_REAL_EPSILON)
      {
         rel_resnorm = init_resnorm / rhs_norm;
      }
      else
      {
         /* rhs is zero, return a zero solution */
         jxf_ParVectorSetConstantValues(U_array, 0.0);
         if(logging > 0)
         {
            rel_resnorm = 0.0;
            jxf_ParNSHDataFinalRelResidualNorm(nsh_data) = rel_resnorm;
         }
         return jxf_error_flag;
      }
   }
   else
   {
      rel_resnorm = 1.;
   }

   if (my_id == 0 && print_level > 1)
   {
      jxf_printf("                                            relative\n");
      jxf_printf("               residual        factor       residual\n");
      jxf_printf("               --------        ------       --------\n");
      jxf_printf("    Initial    %e                 %e\n",init_resnorm,
            rel_resnorm);
   }

   matA = A;
   U_array = u;
   F_array = f;

   /************** Main Solver Loop - always do 1 iteration ************/
   iter = 0;

   while ((rel_resnorm >= tol || iter < 1)
         && iter < max_iter)
   {

      /* Do one solve on e = Mr */
      jxf_NSHSolveInverse(matA, f, u, matM, Utemp, Ftemp);

      /*---------------------------------------------------------------
       *    Compute residual and residual norm
       *----------------------------------------------------------------*/

      if (print_level > 1 || logging > 1 || tol > 0.)
      {
         old_resnorm = resnorm;

         if ( logging > 1 ) {
            jxf_ParVectorCopy(F_array, residual);
            jxf_ParCSRMatrixMatvec(alpha, matA, U_array, beta, residual );
            resnorm = sqrt(jxf_ParVectorInnerProd( residual, residual ));
         }
         else {
            jxf_ParVectorCopy(F_array, Ftemp);
            jxf_ParCSRMatrixMatvec(alpha, matA, U_array, beta, Ftemp);
            resnorm = sqrt(jxf_ParVectorInnerProd(Ftemp, Ftemp));
         }

         if (old_resnorm) conv_factor = resnorm / old_resnorm;
         else conv_factor = resnorm;
         if (rhs_norm > JXF_REAL_EPSILON)
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
      jxf_ParNSHDataNumIterations(nsh_data) = iter;
      jxf_ParNSHDataFinalRelResidualNorm(nsh_data) = rel_resnorm;

      if (my_id == 0 && print_level > 1)
      {
         jxf_printf("    NSHSolve %2d   %e    %f     %e \n", iter,
               resnorm, conv_factor, rel_resnorm);
      }
   }

   /* check convergence within max_iter */
   if (iter == max_iter && tol > 0.)
   {
      Solve_err_flag = 1;
      jxf_error(JXF_ERROR_CONV);
   }

   /*-----------------------------------------------------------------------
    *    Print closing statistics
    *    Add operator and grid complexity stats
    *-----------------------------------------------------------------------*/

   if (iter > 0 && init_resnorm)
   {
      conv_factor = pow((resnorm/init_resnorm),(1.0/(JXF_Real) iter));
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
            jxf_printf("\n\n==============================================");
            jxf_printf("\n NOTE: Convergence tolerance was not achieved\n");
            jxf_printf("      within the allowed %d iterations\n",max_iter);
            jxf_printf("==============================================");
         }
         jxf_printf("\n\n Average Convergence Factor = %f \n",conv_factor);
         jxf_printf("                operator = %f\n",operat_cmplxty);
      }
   }
   return jxf_error_flag;
}

/* NSH solve
 * Simply a matvec on residual with approximate inverse
 * A: original matrix
 * f: rhs
 * u: solution
 * M: approximate inverse
 * ftemp, utemp: working vectors
*/
JXF_Int
jxf_NSHSolveInverse(jxf_ParCSRMatrix *A, jxf_ParVector *f,
                        jxf_ParVector *u, jxf_ParCSRMatrix *M,
                        jxf_ParVector *ftemp, jxf_ParVector *utemp)
{
   JXF_Real      alpha;
   JXF_Real      beta;

   /* begin */
   alpha = -1.0;
   beta = 1.0;
   /* r = f-Au */
   jxf_ParCSRMatrixMatvecOutOfPlace(alpha, A, u, beta, f, ftemp);
   /* e = Mr */
   jxf_ParCSRMatrixMatvec(1.0, M, ftemp, 0.0, utemp);
   /* u = u + e */
   jxf_ParVectorAxpy(beta, utemp, u);
   return jxf_error_flag;
}
