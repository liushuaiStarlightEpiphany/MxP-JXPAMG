//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu_setup.c
 *
 */

#include "jxf_ilu.h"
#include "jxf_krylov.h"

/* Setup ILU data */
JXF_Int
jxf_ILUSetup( void *ilu_vdata, jxf_ParCSRMatrix *A )
{
   MPI_Comm           comm     = jxf_ParCSRMatrixComm(A);
   jxf_ParILUData     *ilu_data = (jxf_ParILUData*) ilu_vdata;
   jxf_ParILUData     *schur_precond_ilu;
   jxf_ParNSHData     *schur_solver_nsh;

   JXF_Int            i;
   // JXF_Int            num_threads;
   // JXF_Int            debug_flag           = 0;

   /* pointers to ilu data */
   JXF_Int            logging              = jxf_ParILUDataLogging(ilu_data);
   JXF_Int            print_level          = jxf_ParILUDataPrintLevel(ilu_data);
   JXF_Int            ilu_type             = jxf_ParILUDataIluType(ilu_data);
   JXF_Int            nLU                  = jxf_ParILUDataNLU(ilu_data);
   JXF_Int            nI                   = jxf_ParILUDataNI(ilu_data);
   JXF_Int            fill_level           = jxf_ParILUDataLfil(ilu_data);
   JXF_Int            max_row_elmts        = jxf_ParILUDataMaxRowNnz(ilu_data);
   JXF_Real           *droptol             = jxf_ParILUDataDroptol(ilu_data);
   JXF_Int            *CF_marker_array     = jxf_ParILUDataCFMarkerArray(ilu_data);
   JXF_Int            *perm                = jxf_ParILUDataPerm(ilu_data);
   JXF_Int            *qperm               = jxf_ParILUDataQPerm(ilu_data);
   JXF_Real           tol_ddPQ             = jxf_ParILUDataTolDDPQ(ilu_data);

   jxf_ParCSRMatrix   *matA                = jxf_ParILUDataMatA(ilu_data);
   jxf_ParCSRMatrix   *matL                = jxf_ParILUDataMatL(ilu_data);
   JXF_Real           *matD                = jxf_ParILUDataMatD(ilu_data);
   jxf_ParCSRMatrix   *matU                = jxf_ParILUDataMatU(ilu_data);
   jxf_ParCSRMatrix   *matS                = jxf_ParILUDataMatS(ilu_data);
   //   jxf_ParCSRMatrix   *matM                = NULL;
   JXF_Real           nnzS/* total nnz in S */;
   JXF_Int            nnzS_offd;
   JXF_Int            size_C/* total size of coarse grid */;

   JXF_Int            n                    = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A));
   //   JXF_Int            m;/* m = n-LU */
   /* reordering option */
   JXF_Int            reordering_type = jxf_ParILUDataReorderingType(ilu_data);
   JXF_Int            num_procs,  my_id;

   jxf_ParVector      *Utemp               = NULL;
   jxf_ParVector      *Ftemp               = NULL;
   JXF_Real           *uext                = NULL;
   JXF_Real           *fext                = NULL;
   jxf_ParVector      *rhs                 = NULL;
   jxf_ParVector      *x                   = NULL;
   jxf_ParVector      *residual            = jxf_ParILUDataResidual(ilu_data);
   JXF_Real           *rel_res_norms       = jxf_ParILUDataRelResNorms(ilu_data);

   /* might need for Schur Complement */
   JXF_Int            *u_end               = NULL;
   JXF_Solver         schur_solver         = NULL;
   JXF_Solver         schur_precond        = NULL;
   JXF_Solver         schur_precond_gotten = NULL;

   /* help to build external */
   jxf_ParCSRCommPkg  *comm_pkg;
   JXF_Int            buffer_size;
   JXF_Int            send_size;
   JXF_Int            recv_size;

   /* ----- begin -----*/

   //num_threads = jxf_NumThreads();

   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);

   /* Free Previously allocated data, if any not destroyed */
   if(matL)
   {
      jxf_ParCSRMatrixDestroy(matL);
      matL = NULL;
   }
   if(matU)
   {
      jxf_ParCSRMatrixDestroy(matU);
      matU = NULL;
   }
   if(matS)
   {
      jxf_ParCSRMatrixDestroy(matS);
      matS = NULL;
   }
   if(matD)
   {
      jxf_TFree(matD);
      matD = NULL;
   }
   if(CF_marker_array)
   {
      jxf_TFree(CF_marker_array);
      CF_marker_array = NULL;
   }


   /* clear old l1_norm data, if created */
   if(jxf_ParILUDataL1Norms(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataL1Norms(ilu_data));
      jxf_ParILUDataL1Norms(ilu_data) = NULL;
   }

   /* setup temporary storage
    * first check is they've already here
    */
   if (jxf_ParILUDataUTemp(ilu_data))
   {
      jxf_ParVectorDestroy(jxf_ParILUDataUTemp(ilu_data));
      jxf_ParILUDataUTemp(ilu_data) = NULL;
   }
   if (jxf_ParILUDataFTemp(ilu_data))
   {
      jxf_ParVectorDestroy(jxf_ParILUDataFTemp(ilu_data));
      jxf_ParILUDataFTemp(ilu_data) = NULL;
   }
   if (jxf_ParILUDataUExt(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataUExt(ilu_data));
      jxf_ParILUDataUExt(ilu_data) = NULL;
   }
   if ( jxf_ParILUDataFExt(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataFExt(ilu_data));
      jxf_ParILUDataFExt(ilu_data) = NULL;
   }
   if ( jxf_ParILUDataUEnd(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataUEnd(ilu_data));
      jxf_ParILUDataUEnd(ilu_data) = NULL;
   }
   if (jxf_ParILUDataRhs(ilu_data))
   {
      jxf_ParVectorDestroy(jxf_ParILUDataRhs(ilu_data));
      jxf_ParILUDataRhs(ilu_data) = NULL;
   }
   if (jxf_ParILUDataX(ilu_data))
   {
      jxf_ParVectorDestroy(jxf_ParILUDataX(ilu_data));
      jxf_ParILUDataX(ilu_data) = NULL;
   }
   if (jxf_ParILUDataResidual(ilu_data))
   {
      jxf_ParVectorDestroy(jxf_ParILUDataResidual(ilu_data));
      jxf_ParILUDataResidual(ilu_data) = NULL;
   }
   if (jxf_ParILUDataRelResNorms(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataRelResNorms(ilu_data));
      jxf_ParILUDataRelResNorms(ilu_data) = NULL;
   }
   if (jxf_ParILUDataSchurSolver(ilu_data))
   {
      switch(ilu_type){
         case 10: case 11:
            JXF_GMRESDestroy(jxf_ParILUDataSchurSolver(ilu_data)); //GMRES for Schur
            break;
         case 20: case 21:
            jxf_NSHDestroy(jxf_ParILUDataSchurSolver(ilu_data)); //NSH for Schur
            break;
         default:
            break;
      }
      (jxf_ParILUDataSchurSolver(ilu_data)) = NULL;
   }
   if(jxf_ParILUDataSchurPrecond(ilu_data))
   {
      switch(ilu_type){
         case 10: case 11:
            JXF_ILUDestroy(jxf_ParILUDataSchurPrecond(ilu_data)); //ILU as precond for Schur
            break;
         default:
            break;
      }
      (jxf_ParILUDataSchurPrecond(ilu_data)) = NULL;
   }
   /* start to create working vectors */
   Utemp = jxf_ParVectorCreate(jxf_ParCSRMatrixComm(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A));
   jxf_ParVectorInitialize(Utemp);
   jxf_ParVectorSetPartitioningOwner(Utemp,0);
   jxf_ParILUDataUTemp(ilu_data) = Utemp;

   Ftemp = jxf_ParVectorCreate(jxf_ParCSRMatrixComm(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A));
   jxf_ParVectorInitialize(Ftemp);
   jxf_ParVectorSetPartitioningOwner(Ftemp,0);
   jxf_ParILUDataFTemp(ilu_data) = Ftemp;
   /* set matrix, solution and rhs pointers */
   matA = A;

   // create perm arary if necessary
   if(perm == NULL)
   {
      switch(ilu_type)
      {
         case 10: case 11: case 20: case 21: case 30: case 31: /* symmetric */
            jxf_ILUGetInteriorExteriorPerm(matA, &perm, &nLU, reordering_type);
            break;
         case 40: case 41:/* ddPQ */
            jxf_ILUGetPermddPQ(matA, &perm, &qperm, tol_ddPQ, &nLU, &nI, reordering_type);
            break;
         case 0: case 1:
            jxf_ILUGetLocalPerm(matA, &perm, &nLU, reordering_type);
            break;
         default:
            jxf_ILUGetLocalPerm(matA, &perm, &nLU, reordering_type);
            break;
      }
   }
   //   m = n - nLU;
   /* factorization */
   switch(ilu_type)
   {
      case 0:  jxf_ILUSetupILUK(matA, fill_level, perm, perm, n, n, &matL, &matD, &matU, &matS, &u_end); //BJ + jxf_iluk()
               break;
      case 1:  jxf_ILUSetupILUT(matA, max_row_elmts, droptol, perm, perm, n, n, &matL, &matD, &matU, &matS, &u_end); //BJ + jxf_ilut()
               break;
      case 10: jxf_ILUSetupILUK(matA, fill_level, perm, perm, nLU, nLU, &matL, &matD, &matU, &matS, &u_end); //GMRES + jxf_iluk()
               break;
      case 11: jxf_ILUSetupILUT(matA, max_row_elmts, droptol, perm, perm, nLU, nLU, &matL, &matD, &matU, &matS, &u_end); //GMRES + jxf_ilut()
               break;
      case 20: jxf_ILUSetupILUK(matA, fill_level, perm, perm, nLU, nLU, &matL, &matD, &matU, &matS, &u_end); //Newton–Schulz–Hotelling + jxf_iluk()
               break;
      case 21: jxf_ILUSetupILUT(matA, max_row_elmts, droptol, perm, perm, nLU, nLU, &matL, &matD, &matU, &matS, &u_end); //Newton–Schulz–Hotelling + jxf_ilut()
               break;
      case 30: jxf_ILUSetupILUKRAS(matA, fill_level, perm, nLU, &matL, &matD, &matU); //RAS + jxf_iluk()
               break;
      case 31: jxf_ILUSetupILUTRAS(matA, max_row_elmts, droptol, perm, nLU, &matL, &matD, &matU); //RAS + jxf_ilut()
               break;
      case 40: jxf_ILUSetupILUK(matA, fill_level, perm, qperm, nLU, nI, &matL, &matD, &matU, &matS, &u_end); //ddPQ + GMRES + jxf_iluk()
               break;
      case 41: jxf_ILUSetupILUT(matA, max_row_elmts, droptol, perm, qperm, nLU, nI, &matL, &matD, &matU, &matS, &u_end); //ddPQ + GMRES + jxf_ilut()
               break;
      default: jxf_ILUSetupILU0(matA, perm, perm, n, n, &matL, &matD, &matU, &matS, &u_end);//BJ + jxf_ilu0()
               break;
   }
   /* setup Schur solver */
   switch(ilu_type)
   {
      case 10: case 11: case 40: case 41:
         if(matS)
         {
            /* setup GMRES parameters */
            JXF_GMRESCreate(comm, &schur_solver);

            JXF_GMRESSetKDim            (schur_solver, jxf_ParILUDataSchurGMRESKDim(ilu_data));
            JXF_GMRESSetMaxIter         (schur_solver, jxf_ParILUDataSchurGMRESMaxIter(ilu_data));/* we don't need that many solves */
            JXF_GMRESSetTol             (schur_solver, (ilu_data -> ss_tol));
            //JXF_GMRESSetTol             (schur_solver, 0.); /* set tol for schur solve to zero. Avoids triggering hypre error for non convergence -DOK*/
            JXF_GMRESSetAbsoluteTol     (schur_solver, (ilu_data -> ss_absolute_tol));
            JXF_GMRESSetLogging         (schur_solver, (ilu_data -> ss_logging));
            JXF_GMRESSetPrintLevel      (schur_solver, (ilu_data -> ss_print_level));/* set to zero now, don't print */
            JXF_GMRESSetRelChange       (schur_solver, (ilu_data -> ss_rel_change));

            /* setup preconditioner parameters */
            /* create precond, the default is ILU0 */
            JXF_ILUCreate               (&schur_precond);
            JXF_ILUSetType              (schur_precond, (ilu_data -> sp_ilu_type));
            JXF_ILUSetLevelOfFill       (schur_precond, (ilu_data -> sp_ilu_lfil));
            JXF_ILUSetMaxNnzPerRow      (schur_precond, (ilu_data -> sp_ilu_max_row_nnz));
            JXF_ILUSetDropThresholdArray(schur_precond, (ilu_data -> sp_ilu_droptol));
            jxf_ILUSetOwnDropThreshold  (schur_precond, 0);/* using exist droptol */
            JXF_ILUSetPrintLevel        (schur_precond, (ilu_data -> sp_print_level));
            JXF_ILUSetMaxIter           (schur_precond, (ilu_data -> sp_max_iter));
            //JXF_ILUSetTol               (schur_precond, (ilu_data -> sp_tol));
            JXF_ILUSetTol               (schur_precond, 0.); /* set tol for preconditioner to zero. Avoids triggering hypre error for non convergence -DOK*/

            /* add preconditioner to solver */
            JXF_GMRESSetPrecond(schur_solver,
                  (JXF_PtrToSolverFcn) JXF_ILUSolve,
                  (JXF_PtrToSolverFcn) JXF_ILUSetup,
                  schur_precond);
            JXF_GMRESGetPrecond(schur_solver, &schur_precond_gotten);
            if (schur_precond_gotten != (schur_precond))
            {
               jxf_printf("Schur complement got bad precond\n");
               return(-1);
            }

            /* need to create working vector rhs and x for Schur System */
            rhs = jxf_ParVectorCreate(comm,
                  jxf_ParCSRMatrixGlobalNumRows(matS),
                  jxf_ParCSRMatrixRowStarts(matS));
            jxf_ParVectorInitialize(rhs);
            jxf_ParVectorSetPartitioningOwner(rhs,0);
            x = jxf_ParVectorCreate(comm,
                  jxf_ParCSRMatrixGlobalNumRows(matS),
                  jxf_ParCSRMatrixRowStarts(matS));
            jxf_ParVectorInitialize(x);
            jxf_ParVectorSetPartitioningOwner(x,0);

            /* setup solver */
            JXF_GMRESSetup(schur_solver,(JXF_Matrix)matS,(JXF_Vector)rhs,(JXF_Vector)x);

            /* update ilu_data */
            jxf_ParILUDataSchurSolver   (ilu_data) = schur_solver;
            jxf_ParILUDataSchurPrecond  (ilu_data) = schur_precond;
            jxf_ParILUDataRhs           (ilu_data) = rhs;
            jxf_ParILUDataX             (ilu_data) = x;
         }
         break;
      case 20: case 21:
         if(matS)
         {
            /* approximate inverse preconditioner */
            schur_solver = (JXF_Solver)jxf_NSHCreate();

            /* set NSH parameters */
            jxf_NSHSetMaxIter           (schur_solver, jxf_ParILUDataSchurNSHSolveMaxIter(ilu_data));
            jxf_NSHSetTol               (schur_solver, jxf_ParILUDataSchurNSHSolveTol(ilu_data));
            jxf_NSHSetLogging           (schur_solver, jxf_ParILUDataSchurSolverLogging(ilu_data));
            jxf_NSHSetPrintLevel        (schur_solver, jxf_ParILUDataSchurSolverPrintLevel(ilu_data));
            jxf_NSHSetDropThresholdArray(schur_solver, jxf_ParILUDataSchurNSHDroptol(ilu_data));

            jxf_NSHSetNSHMaxIter        (schur_solver, jxf_ParILUDataSchurNSHMaxNumIter(ilu_data));
            jxf_NSHSetNSHMaxRowNnz      (schur_solver, jxf_ParILUDataSchurNSHMaxRowNnz(ilu_data));
            jxf_NSHSetNSHTol            (schur_solver, jxf_ParILUDataSchurNSHTol(ilu_data));

            jxf_NSHSetMRMaxIter         (schur_solver, jxf_ParILUDataSchurMRMaxIter(ilu_data));
            jxf_NSHSetMRMaxRowNnz       (schur_solver, jxf_ParILUDataSchurMRMaxRowNnz(ilu_data));
            jxf_NSHSetMRTol             (schur_solver, jxf_ParILUDataSchurMRTol(ilu_data));
            jxf_NSHSetColVersion        (schur_solver, jxf_ParILUDataSchurMRColVersion(ilu_data));

            /* need to create working vector rhs and x for Schur System */
            rhs = jxf_ParVectorCreate(comm,
                  jxf_ParCSRMatrixGlobalNumRows(matS),
                  jxf_ParCSRMatrixRowStarts(matS));
            jxf_ParVectorInitialize(rhs);
            jxf_ParVectorSetPartitioningOwner(rhs,0);
            x = jxf_ParVectorCreate(comm,
                  jxf_ParCSRMatrixGlobalNumRows(matS),
                  jxf_ParCSRMatrixRowStarts(matS));
            jxf_ParVectorInitialize(x);
            jxf_ParVectorSetPartitioningOwner(x,0);

            /* setup solver */
            jxf_NSHSetup(schur_solver,matS,rhs,x);

            jxf_ParILUDataSchurSolver(ilu_data) = schur_solver;
            jxf_ParILUDataRhs        (ilu_data) = rhs;
            jxf_ParILUDataX          (ilu_data) = x;
         }
         break;
      case 30 : case 31:
         /* now check communication package */
         comm_pkg = jxf_ParCSRMatrixCommPkg(matA);
         /* create if not yet built */
         if(!comm_pkg)
         {
            jxf_MatvecCommPkgCreate(matA);
            comm_pkg = jxf_ParCSRMatrixCommPkg(matA);
         }
         /* create uext and fext */
         send_size =  jxf_ParCSRCommPkgSendMapStart(comm_pkg,jxf_ParCSRCommPkgNumSends(comm_pkg))
            - jxf_ParCSRCommPkgSendMapStart(comm_pkg,0);
         recv_size = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(matA));
         buffer_size = send_size > recv_size ? send_size : recv_size;
         fext = jxf_TAlloc(JXF_Real,buffer_size);
         uext = jxf_TAlloc(JXF_Real,buffer_size);
         break;
      default:
         break;
   }
   /* set pointers to ilu data */
   jxf_ParILUDataMatA(ilu_data)            = matA;
   jxf_ParILUDataMatL(ilu_data)            = matL;
   jxf_ParILUDataMatD(ilu_data)            = matD;
   jxf_ParILUDataMatU(ilu_data)            = matU;
   jxf_ParILUDataMatS(ilu_data)            = matS;
   jxf_ParILUDataCFMarkerArray(ilu_data)   = CF_marker_array;
   jxf_ParILUDataPerm(ilu_data)            = perm;
   jxf_ParILUDataQPerm(ilu_data)           = qperm;
   jxf_ParILUDataNLU(ilu_data)             = nLU;
   jxf_ParILUDataNI(ilu_data)              = nI;
   jxf_ParILUDataUEnd(ilu_data)            = u_end;
   jxf_ParILUDataUExt(ilu_data)            = uext;
   jxf_ParILUDataFExt(ilu_data)            = fext;

   /* compute operator complexity */
   jxf_ParCSRMatrixSetDNumNonzeros(matA);
   nnzS = 0.0;
   /* size_C is the size of global coarse grid, upper left part */
   size_C = jxf_ParCSRMatrixGlobalNumRows(matA);
   /* switch to compute complexity */
   if(matS)
   {
      jxf_ParCSRMatrixSetDNumNonzeros(matS);
      nnzS = jxf_ParCSRMatrixDNumNonzeros(matS);
      /* if we have Schur system need to reduce it from size_C */
      size_C -= jxf_ParCSRMatrixGlobalNumRows(matS);
      switch(ilu_type)
      {
         case 10: case 11: case 40: case 41:
            /* now we need to compute the preoconditioner */
            schur_precond_ilu = (jxf_ParILUData*) (ilu_data -> schur_precond);
            /* borrow i for local nnz of S */
            i = jxf_CSRMatrixNumNonzeros(jxf_ParCSRMatrixOffd(matS));
            jxf_MPI_Allreduce(&i, &nnzS_offd, 1, JXF_MPI_INT, MPI_SUM, comm);
            nnzS = nnzS * (schur_precond_ilu -> operator_complexity) +nnzS_offd;
            break;
         case 20: case 21:
            schur_solver_nsh = (jxf_ParNSHData*) jxf_ParILUDataSchurSolver(ilu_data);
            nnzS = nnzS * (jxf_ParNSHDataOperatorComplexity(schur_solver_nsh));
            break;
         default:
            break;
      }
   }

   (ilu_data -> operator_complexity) =  ((JXF_Real)size_C + nnzS +
         jxf_ParCSRMatrixDNumNonzeros(matL) +
         jxf_ParCSRMatrixDNumNonzeros(matU)) /
      jxf_ParCSRMatrixDNumNonzeros(matA);
   if ((my_id == 0) && (print_level > 0))
   {
      jxf_printf("ILU SETUP: operator complexity = %f  \n", ilu_data -> operator_complexity);
   }

   if ( logging > 1 ) {
      residual =
         jxf_ParVectorCreate(jxf_ParCSRMatrixComm(matA),
               jxf_ParCSRMatrixGlobalNumRows(matA),
               jxf_ParCSRMatrixRowStarts(matA) );
      jxf_ParVectorInitialize(residual);
      jxf_ParVectorSetPartitioningOwner(residual,0);
      (ilu_data -> residual) = residual;
   }
   else{
      (ilu_data -> residual) = NULL;
   }
   rel_res_norms = jxf_CTAlloc(JXF_Real, (ilu_data -> max_iter));
   (ilu_data -> rel_res_norms) = rel_res_norms;

   return jxf_error_flag;
}

/* ILU(0)
 * A = input matrix
 * perm = permutation array indicating ordering of rows. Perm could come from a
 *    CF_marker array or a reordering routine.
 * qperm = permutation array indicating ordering of columns
 * nI = number of interial unknowns
 * nLU = size of incomplete factorization, nLU should obey nLU <= nI.
 *    Schur complement is formed if nLU < n
 * Lptr, Dptr, Uptr, Sptr = L, D, U, S factors.
 * will form global Schur Matrix if nLU < n
 */
JXF_Int
jxf_ILUSetupILU0(jxf_ParCSRMatrix *A, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI,
      jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end)
{
   JXF_Int                i, ii, j, k, k1, k2, k3, ctrU, ctrL, ctrS, lenl, lenu, jpiv, col, jpos;
   JXF_Int                *iw, *iL, *iU;
   JXF_Real               dd, t, dpiv, lxu, *wU, *wL;

   /* communication stuffs for S */
   MPI_Comm                 comm             = jxf_ParCSRMatrixComm(A);
   JXF_Int                S_offd_nnz, S_offd_ncols;
   jxf_ParCSRCommPkg      *comm_pkg;
   jxf_ParCSRCommHandle   *comm_handle;
   JXF_Int                num_sends, begin, end;
   JXF_Int                *send_buf        = NULL;
   JXF_Int                num_procs, my_id;

   /* data objects for A */
   jxf_CSRMatrix          *A_diag          = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix          *A_offd          = jxf_ParCSRMatrixOffd(A);
   JXF_Real               *A_diag_data     = jxf_CSRMatrixData(A_diag);
   JXF_Int                *A_diag_i        = jxf_CSRMatrixI(A_diag);
   JXF_Int                *A_diag_j        = jxf_CSRMatrixJ(A_diag);
   JXF_Real               *A_offd_data     = jxf_CSRMatrixData(A_offd);
   JXF_Int                *A_offd_i        = jxf_CSRMatrixI(A_offd);
   JXF_Int                *A_offd_j        = jxf_CSRMatrixJ(A_offd);

   /* size of problem and schur system */
   JXF_Int                n                =  jxf_CSRMatrixNumRows(A_diag);
   JXF_Int                m                = n - nLU;
   JXF_Int                e                = nI - nLU;
   JXF_Int                m_e              = n - nI;
   JXF_Real               local_nnz, total_nnz;
   JXF_Int                *u_end_array;
   JXF_Int                u_end_location;

   /* data objects for L, D, U */
   jxf_ParCSRMatrix       *matL;
   jxf_ParCSRMatrix       *matU;
   jxf_CSRMatrix          *L_diag;
   jxf_CSRMatrix          *U_diag;
   JXF_Real               *D_data;
   JXF_Real               *L_diag_data;
   JXF_Int                *L_diag_i;
   JXF_Int                *L_diag_j;
   JXF_Real               *U_diag_data;
   JXF_Int                *U_diag_i;
   JXF_Int                *U_diag_j;

   /* data objects for S */
   jxf_ParCSRMatrix       *matS = NULL;
   jxf_CSRMatrix          *S_diag;
   jxf_CSRMatrix          *S_offd;
   JXF_Real               *S_diag_data     = NULL;
   JXF_Int                *S_diag_i        = NULL;
   JXF_Int                *S_diag_j        = NULL;
   JXF_Int                *S_offd_i        = NULL;
   JXF_Int                *S_offd_j        = NULL;
   JXF_Int                *S_offd_colmap   = NULL;
   JXF_Real               *S_offd_data;
   JXF_Int                *col_starts;
   JXF_Int                total_rows;

   /* memory management */
   JXF_Int                initial_alloc    = 0;
   JXF_Int                capacity_L;
   JXF_Int                capacity_U;
   JXF_Int                capacity_S       = 0;
   JXF_Int                nnz_A            = A_diag_i[n];

   /* reverse permutation array */
   JXF_Int                *rperm;

   /* start setup
    * get communication stuffs first
    */
   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   /* setup if not yet built */
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* check for correctness */
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }
   if(e < 0)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU should not exceed nI.\n");
   }

   /* Allocate memory for u_end array */
   u_end_array    = jxf_TAlloc(JXF_Int, nLU);

   /* Allocate memory for L,D,U,S factors */
   if(n > 0)
   {
      initial_alloc  = nLU + ceil((nnz_A / 2.0)*nLU/n);
      capacity_S     = m + ceil((nnz_A / 2.0)*m/n);
   }
   capacity_L     = initial_alloc;
   capacity_U     = initial_alloc;

   D_data         = jxf_TAlloc(JXF_Real, n);
   L_diag_i       = jxf_TAlloc(JXF_Int, n+1);
   L_diag_j       = jxf_TAlloc(JXF_Int, capacity_L);
   L_diag_data    = jxf_TAlloc(JXF_Real, capacity_L);
   U_diag_i       = jxf_TAlloc(JXF_Int, n+1);
   U_diag_j       = jxf_TAlloc(JXF_Int, capacity_U);
   U_diag_data    = jxf_TAlloc(JXF_Real, capacity_U);
   S_diag_i       = jxf_TAlloc(JXF_Int, m+1);
   S_diag_j       = jxf_TAlloc(JXF_Int, capacity_S);
   S_diag_data    = jxf_TAlloc(JXF_Real, capacity_S);

   /* allocate working arrays */
   iw             = jxf_TAlloc(JXF_Int, 3*n);
   iL             = iw+n;
   rperm          = iw + 2*n;
   wL             = jxf_TAlloc(JXF_Real, n);

   ctrU        = ctrL        = ctrS        = 0;
   L_diag_i[0] = U_diag_i[0] = S_diag_i[0] = 0;
   /* set marker array iw to -1 */
   for( i = 0; i < n; i++ )
   {
      iw[i] = -1;
   }

   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    * rperm only used for column
    */
   for(i=0; i<n; i++)
   {
      rperm[qperm[i]] = i;
   }
   /*---------  Begin Factorization. Work in permuted space  ----*/
   for( ii = 0; ii < nLU; ii++ )
   {

      // get row i
      i = perm[ii];
      // get extents of row i
      k1=A_diag_i[i];
      k2=A_diag_i[i+1];

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL+ii;
      wU = wL+ii;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl  = lenu = 0;
      iw[ii] = ii;
      /*-------------------- scan & unwrap column */
      for(j=k1; j < k2; j++)
      {
         col = rperm[A_diag_j[j]];
         t = A_diag_data[j];
         if( col < ii )
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
            dd=t;
         }
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jxf_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jxf_qsort3ir(iL, wL, iw, 0, (lenl-1));
      for(j=0; j<lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for(k=U_diag_i[jpiv]; k<U_diag_i[jpiv+1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if(jpos < 0)
            {
               continue;
            }

            lxu = - U_diag_data[k] * dpiv;
            if(col < ii)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if(col > ii)
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
      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for( j = 0; j < lenu; j++ )
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      if(lenl > 0)
      {
         while((ctrL+lenl) > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         jxf_TMemcpy(&(L_diag_j)[ctrL], iL, JXF_Int, lenl);
         jxf_TMemcpy(&(L_diag_data)[ctrL], wL, JXF_Real, lenl);
      }
      L_diag_i[ii+1] = (ctrL+=lenl);

      /* diagonal part (we store the inverse) */
      if(fabs(dd) < MAT_TOL)
      {
         dd = 1.0e-6;
      }
      D_data[ii] = 1./dd;

      /* U part */
      /* Check that memory is sufficient */
      if(lenu > 0)
      {
         while((ctrU+lenu) > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
         }
         jxf_TMemcpy(&(U_diag_j)[ctrU], iU, JXF_Int, lenu);
         jxf_TMemcpy(&(U_diag_data)[ctrU], wU, JXF_Real, lenu);
      }
      U_diag_i[ii+1] = (ctrU+=lenu);

      /* check and build u_end array */
      if(m > 0)
      {
         jxf_qsort1(U_diag_j,U_diag_data,U_diag_i[ii],U_diag_i[ii+1]-1);
         u_end_location = jxf_BinarySearch2(U_diag_j,nLU,U_diag_i[ii],U_diag_i[ii+1]-1,u_end_array + ii);
         if(u_end_location >= 0)
         {
            u_end_array[ii] = u_end_location + 1;
         }
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }

   }
   /*---------  Begin Factorization in Schur Complement part  ----*/
   for( ii = nLU; ii < n; ii++ )
   {
      // get row i
      i = perm[ii];
      // get extents of row i
      k1=A_diag_i[i];
      k2=A_diag_i[i+1];

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL+nLU + 1;
      wU = wL+nLU + 1;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl  = lenu = 0;
      iw[ii] = nLU;
      /*-------------------- scan & unwrap column */
      for(j=k1; j < k2; j++)
      {
         col = rperm[A_diag_j[j]];
         t = A_diag_data[j];
         if( col < nLU )
         {
            iw[col] = lenl;
            iL[lenl] = col;
            wL[lenl++] = t;
         }
         else if (col != ii)
         {
            iw[col] = lenu;
            iU[lenu] = col;
            wU[lenu++] = t;
         }
         else
         {
            dd=t;
         }
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jxf_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jxf_qsort3ir(iL, wL, iw, 0, (lenl-1));
      for(j=0; j<lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for(k=U_diag_i[jpiv]; k<U_diag_i[jpiv+1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if(jpos < 0)
            {
               continue;
            }

            lxu = - U_diag_data[k] * dpiv;
            if(col < nLU)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if(col != ii)
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
      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for( j = 0; j < lenu; j++ )
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      if(lenl > 0)
      {
         while((ctrL+lenl) > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         jxf_TMemcpy(&(L_diag_j)[ctrL], iL, JXF_Int, lenl);
         jxf_TMemcpy(&(L_diag_data)[ctrL], wL, JXF_Real, lenl);
      }
      L_diag_i[ii+1] = (ctrL+=lenl);

      /* S part */
      /* Check that memory is sufficient */
      while((ctrS+lenu+1) > capacity_S)
      {
         JXF_Int tmp = capacity_S;
         capacity_S = capacity_S * EXPAND_FACT + 1;
         S_diag_j = jxf_TReAlloc_v2(S_diag_j, JXF_Int, tmp, JXF_Int, capacity_S);
         S_diag_data = jxf_TReAlloc_v2(S_diag_data, JXF_Real, tmp, JXF_Real, capacity_S);
      }
      /* remember S in under a new index system! */
      S_diag_j[ctrS] = ii - nLU;
      S_diag_data[ctrS] = dd;
      for(j = 0 ; j< lenu ; j ++)
      {
         S_diag_j[ctrS+1+j] = iU[j] - nLU;
      }
      //jxf_TMemcpy(S_diag_j+ctrS+1, iU, JXF_Int, lenu);
      jxf_TMemcpy(S_diag_data+ctrS+1, wU, JXF_Real, lenu);
      S_diag_i[ii-nLU+1] = ctrS+=(lenu+1);
   }
   /* Assemble LDUS matrices */
   /* zero out unfactored rows for U and D */
   for(k=nLU; k<n; k++)
   {
      U_diag_i[k+1] = ctrU;
      D_data[k] = 1.;
   }

   /* First create Schur complement if necessary
    * Check if we need to create Schur complement
    */
   JXF_Int big_m = (JXF_Int)m;
   jxf_MPI_Allreduce(&big_m, &total_rows, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* only form when total_rows > 0 */
   if( total_rows > 0 )
   {
      /* now create S */
      /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
      {
         JXF_Int global_start;
         col_starts = jxf_CTAlloc(JXF_Int,2);
         jxf_MPI_Scan( &big_m, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
         col_starts[0] = global_start - m;
         col_starts[1] = global_start;
      }
#else
      col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

      jxf_MPI_Allgather(&big_m,1,JXF_MPI_INT,&col_starts[1],
            1,JXF_MPI_INT,comm);

      for (i=2; i < num_procs+1; i++)
         col_starts[i] += col_starts[i-1];
#endif

      /* We did nothing to A_offd, so all the data kept, just reorder them
       * The create function takes comm, global num rows/cols,
       *    row/col start, num cols offd, nnz diag, nnz offd
       */
      S_offd_nnz = jxf_CSRMatrixNumNonzeros(A_offd);
      S_offd_ncols = jxf_CSRMatrixNumCols(A_offd);

      matS = jxf_ParCSRMatrixCreate( comm,
            total_rows,
            total_rows,
            col_starts,
            col_starts,
            S_offd_ncols,
            ctrS,
            S_offd_nnz);

      /* S owns different start/end */
      jxf_ParCSRMatrixSetColStartsOwner(matS,1);
      jxf_ParCSRMatrixSetRowStartsOwner(matS,0);/* square matrix, use same row and col start */

      /* first put diagonal data in */
      S_diag = jxf_ParCSRMatrixDiag(matS);

      jxf_CSRMatrixI(S_diag) = S_diag_i;
      jxf_CSRMatrixData(S_diag) = S_diag_data;
      jxf_CSRMatrixJ(S_diag) = S_diag_j;

      /* now start to construct offdiag of S */
      S_offd = jxf_ParCSRMatrixOffd(matS);
      S_offd_i = jxf_TAlloc(JXF_Int, m+1);
      S_offd_j = jxf_TAlloc(JXF_Int, S_offd_nnz);
      S_offd_data = jxf_TAlloc(JXF_Real, S_offd_nnz);
      S_offd_colmap = jxf_CTAlloc(JXF_Int, S_offd_ncols);

      /* simply use a loop to copy data from A_offd */
      S_offd_i[0] = 0;
      k3 = 0;
      for(i = 1 ; i <= e ; i ++)
      {
         S_offd_i[i] = k3;
      }
      for(i = 0 ; i < m_e ; i ++)
      {
         col = perm[i + nI];
         k1 = A_offd_i[col];
         k2 = A_offd_i[col+1];
         for(j = k1 ; j < k2 ; j ++)
         {
            S_offd_j[k3] = A_offd_j[j];
            S_offd_data[k3++] = A_offd_data[j];
         }
         S_offd_i[i+1+e] = k3;
      }

      /* give I, J, DATA to S_offd */
      jxf_CSRMatrixI(S_offd) = S_offd_i;
      jxf_CSRMatrixJ(S_offd) = S_offd_j;
      jxf_CSRMatrixData(S_offd) = S_offd_data;

      /* now we need to update S_offd_colmap */

      /* get total num of send */
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      begin = jxf_ParCSRCommPkgSendMapStart(comm_pkg,0);
      end = jxf_ParCSRCommPkgSendMapStart(comm_pkg,num_sends);
      send_buf = jxf_TAlloc(JXF_Int, end - begin);
      /* copy new index into send_buf */
#ifdef JXF_NO_GLOBAL_PARTITION
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[0];
      }
#else
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[my_id];
      }
#endif
      /* main communication */
      comm_handle = jxf_ParCSRCommHandleCreate(22, comm_pkg, send_buf, S_offd_colmap);
      jxf_ParCSRCommHandleDestroy(comm_handle);

      /* setup index */
      jxf_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

      jxf_ILUSortOffdColmap(matS);

      /* free */
      jxf_TFree(send_buf);
   }/* end of forming S */

   /* create S finished */

   matL = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         ctrL,
         0 );

   /* Have A own row/col partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (ctrL)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jxf_TFree(L_diag_j);
      jxf_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) ctrL;
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         ctrU,
         0 );

   /* Have A own row/col partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);
   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (ctrU)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jxf_TFree(U_diag_j);
      jxf_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) ctrU;
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;
   /* free memory */
   jxf_TFree(wL);
   jxf_TFree(iw);
   if(!matS)
   {
      /* we allocate some memory for S, need to free if unused */
      jxf_TFree(S_diag_i);
      //      jxf_TFree(col_starts);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;
   *u_end = u_end_array;

   return jxf_error_flag;
}

/* ILU(k) symbolic factorization
 * n = total rows of input
 * lfil = level of fill-in, the k in ILU(k)
 * perm = permutation array indicating ordering of factorization. Perm could come from a
 * rperm = reverse permutation array, used here to avoid duplicate memory allocation
 * iw = working array, used here to avoid duplicate memory allocation
 * nLU = size of computed LDU factorization.
 * A/L/U/S_diag_i = the I slot of A, L, U and S
 * A/L/U/S_diag_j = the J slot of A, L, U and S
 * will form global Schur Matrix if nLU < n
 */
JXF_Int
jxf_ILUSetupILUKSymbolic(JXF_Int n, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Int lfil, JXF_Int *perm,
      JXF_Int *rperm,   JXF_Int *iw,   JXF_Int nLU, JXF_Int *L_diag_i, JXF_Int *U_diag_i,
      JXF_Int *S_diag_i, JXF_Int **L_diag_j, JXF_Int **U_diag_j, JXF_Int **S_diag_j, JXF_Int **u_end)
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

   JXF_Int         *temp_L_diag_j, *temp_U_diag_j, *temp_S_diag_j = NULL, *u_levels;
   JXF_Int         *iL, *iLev;
   JXF_Int         ii, i, j, k, ku, lena, lenl, lenu, lenh, ilev, lev, col, icol;
   JXF_Int         m = n - nLU;
   JXF_Int         *u_end_array;
   JXF_Int         u_end_location;

   /* memory management */
   JXF_Int         ctrL;
   JXF_Int         ctrU;
   JXF_Int         ctrS;
   JXF_Int         capacity_L;
   JXF_Int         capacity_U;
   JXF_Int         capacity_S = 0;
   JXF_Int         initial_alloc = 0;
   JXF_Int         nnz_A;

   /* set iL and iLev to right place in iw array */
   iL                = iw + n;
   iLev              = iw + 2*n;

   /* setup initial memory used */
   nnz_A             = A_diag_i[n];
   if(n > 0)
   {
      initial_alloc     = nLU + ceil((nnz_A / 2.0) * nLU / n);
   }
   capacity_L        = initial_alloc;
   capacity_U        = initial_alloc;

   /* allocate other memory for L and U struct */
   temp_L_diag_j     = jxf_CTAlloc(JXF_Int, capacity_L);
   temp_U_diag_j     = jxf_CTAlloc(JXF_Int, capacity_U);

   if(m > 0)
   {
      capacity_S     = m + ceil(nnz_A / 2.0 * m / n);
      temp_S_diag_j  = jxf_CTAlloc(JXF_Int, capacity_S);
   }

   u_end_array       = jxf_TAlloc(JXF_Int, nLU);
   u_levels          = jxf_CTAlloc(JXF_Int, capacity_U);
   ctrL = ctrU = ctrS = 0;

   /* set initial value for working array */
   for(ii = 0 ; ii < n ; ii ++)
   {
      iw[ii] = -1;
   }

   /*
    * 2: Start of main loop
    * those in iL are NEW col index (after permutation)
    */
   for(ii = 0 ; ii < nLU ; ii ++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0;/* this is the current length of heap */
      lenu = ii;
      lena = A_diag_i[i+1];
      /* put those already inside original pattern, and set their level to 0 */
      for(j = A_diag_i[i] ; j < lena ; j ++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
         }
         else if(col > ii)
         {
            /* this is an entry in U */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      }/* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while(lenh > 0)
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
         jxf_ILUMinHeapRemoveIIIi(iL,iLev,iw,lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k]=-1;
         jxf_swap2i(iL,iLev,ii-lenl,lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k+1];
         for(j = U_diag_i[k] ; j < ku ; j ++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if(lev > lfil)
            {
               continue;
            }
            if(icol < 0)
            {
               /* not yet in */
               if(col < ii)
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
                  jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
               }
               else if(col > ii)
               {
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jxf_min(lev, iLev[icol]);
            }
         }/* end of loop j for level update */
      }/* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* check if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jxf_TReAlloc_v2(temp_L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for(j = 0 ; j < lenl ; j ++)
         {
            temp_L_diag_j[ctrL+j] = iL[ii-j-1];
         }
         ctrL += lenl;
      }
      k = lenu - ii;
      U_diag_i[ii+1] = U_diag_i[ii] + k;
      if(k > 0)
      {
         /* check if memory is enough */
         while(ctrU + k > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            temp_U_diag_j = jxf_TReAlloc_v2(temp_U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            u_levels = jxf_TReAlloc_v2(u_levels, JXF_Int, tmp, JXF_Int, capacity_U);
         }
         jxf_TMemcpy(temp_U_diag_j+ctrU,iL+ii,JXF_Int,k);
         jxf_TMemcpy(u_levels+ctrU,iLev+ii,JXF_Int,k);
         ctrU += k;
      }
      if(m > 0)
      {
         jxf_qsort2i(temp_U_diag_j,u_levels,U_diag_i[ii],U_diag_i[ii+1]-1);
         u_end_location = jxf_BinarySearch2(temp_U_diag_j,nLU,U_diag_i[ii],U_diag_i[ii+1]-1,u_end_array + ii);
         if(u_end_location >= 0)
         {
            u_end_array[ii] = u_end_location + 1;
         }
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }

      /* reset iw */
      for(j = ii ; j < lenu ; j ++)
      {
         iw[iL[j]] = -1;
      }

   }/* end of main loop ii from 0 to nLU-1 */

   /* another loop to set EU^-1 and Schur complement */
   for(ii = nLU ; ii < n ; ii ++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0;/* this is the current length of heap */
      lenu = nLU;/* now this stores S, start from nLU */
      lena = A_diag_i[i+1];
      /* put those already inside original pattern, and set their level to 0 */
      for(j = A_diag_i[i] ; j < lena ; j ++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if(col < nLU)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
         }
         else if(col != ii) /* we for sure to add ii, avoid duplicate */
         {
            /* this is an entry in S */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      }/* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while(lenh > 0)
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
         jxf_ILUMinHeapRemoveIIIi(iL,iLev,iw,lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k]=-1;
         jxf_swap2i(iL,iLev,nLU-lenl,lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k+1];
         for(j = U_diag_i[k] ; j < ku ; j ++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if(lev > lfil)
            {
               continue;
            }
            if(icol < 0)
            {
               /* not yet in */
               if(col < nLU)
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
                  jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
               }
               else if(col != ii)
               {
                  /* S part */
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jxf_min(lev, iLev[icol]);
            }
         }/* end of loop j for level update */
      }/* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* check if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jxf_TReAlloc_v2(temp_L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for(j = 0 ; j < lenl ; j ++)
         {
            temp_L_diag_j[ctrL+j] = iL[nLU-j-1];
         }
         ctrL += lenl;
      }
      k = lenu - nLU + 1;
      /* check if memory is enough */
      while(ctrS + k > capacity_S)
      {
         JXF_Int tmp = capacity_S;
         capacity_S = capacity_S * EXPAND_FACT + 1;
         temp_S_diag_j = jxf_TReAlloc_v2(temp_S_diag_j, JXF_Int, tmp, JXF_Int, capacity_S);
      }
      temp_S_diag_j[ctrS] = ii;/* must have diagonal */
      jxf_TMemcpy(temp_S_diag_j+ctrS+1,iL+nLU,JXF_Int,k-1);
      ctrS += k;
      S_diag_i[ii-nLU+1] = ctrS;

      /* reset iw */
      for(j = nLU ; j < lenu ; j ++)
      {
         iw[iL[j]] = -1;
      }

   }/* end of main loop ii from nLU to n-1 */

   /*
    * 3: Update the struct for L, U and S
    */
   for(k = nLU ; k < n ; k ++)
   {
      U_diag_i[k+1] = U_diag_i[nLU];
   }
   /*
    * 4: Finishing up and free memory
    */
   jxf_TFree(u_levels);

   *L_diag_j = temp_L_diag_j;
   *U_diag_j = temp_U_diag_j;
   *S_diag_j = temp_S_diag_j;
   *u_end = u_end_array;

   return jxf_error_flag;
}

/* ILU(k)
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * perm: permutation array indicating ordering of factorization. Perm could come from a
 *    CF_marker: array or a reordering routine.
 * qperm: column permutation array.
 * nLU: size of computed LDU factorization.
 * nI: number of interial unknowns, nI should obey nI >= nLU
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sprt: Schur Complement, if no Schur Complement is needed it will be set to NULL
 */
JXF_Int
jxf_ILUSetupILUK(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI,
      jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end)
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
   if(lfil == 0)
   {
      return jxf_ILUSetupILU0(A,perm,qperm,nLU,nI,Lptr,Dptr,Uptr,Sptr,u_end);
   }
   JXF_Real              local_nnz, total_nnz;
   JXF_Int               i, ii, j, k, k1, k2, k3, kl, ku, jpiv, col, icol;
   JXF_Int               *iw;
   MPI_Comm                comm = jxf_ParCSRMatrixComm(A);
   JXF_Int            num_procs,  my_id;

   /* data objects for A */
   jxf_CSRMatrix         *A_diag        = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix         *A_offd        = jxf_ParCSRMatrixOffd(A);
   JXF_Real              *A_diag_data   = jxf_CSRMatrixData(A_diag);
   JXF_Int               *A_diag_i      = jxf_CSRMatrixI(A_diag);
   JXF_Int               *A_diag_j      = jxf_CSRMatrixJ(A_diag);
   JXF_Real              *A_offd_data   = jxf_CSRMatrixData(A_offd);
   JXF_Int               *A_offd_i      = jxf_CSRMatrixI(A_offd);
   JXF_Int               *A_offd_j      = jxf_CSRMatrixJ(A_offd);

   /* data objects for L, D, U */
   jxf_ParCSRMatrix      *matL;
   jxf_ParCSRMatrix      *matU;
   jxf_CSRMatrix         *L_diag;
   jxf_CSRMatrix         *U_diag;
   JXF_Real              *D_data;
   JXF_Real              *L_diag_data   = NULL;
   JXF_Int               *L_diag_i;
   JXF_Int               *L_diag_j      = NULL;
   JXF_Real              *U_diag_data   = NULL;
   JXF_Int               *U_diag_i;
   JXF_Int               *U_diag_j      = NULL;

   /* data objects for S */
   jxf_ParCSRMatrix      *matS          = NULL;
   jxf_CSRMatrix         *S_diag;
   jxf_CSRMatrix         *S_offd;
   JXF_Real              *S_diag_data   = NULL;
   JXF_Int               *S_diag_i      = NULL;
   JXF_Int               *S_diag_j      = NULL;
   JXF_Int               *S_offd_i      = NULL;
   JXF_Int               *S_offd_j      = NULL;
   JXF_Int               *S_offd_colmap = NULL;
   JXF_Real              *S_offd_data;
   JXF_Int               S_offd_nnz, S_offd_ncols;
   JXF_Int               *col_starts;
   JXF_Int               total_rows;
   /* communication */
   jxf_ParCSRCommPkg     *comm_pkg;
   jxf_ParCSRCommHandle  *comm_handle;
   JXF_Int               *send_buf      = NULL;

   /* problem size */
   JXF_Int               n;
   JXF_Int               m;
   JXF_Int               e;
   JXF_Int               m_e;
   /* reverse permutation array */
   JXF_Int               *rperm;

   /* start setup */
   /* check input and get problem size */
   n =  jxf_CSRMatrixNumRows(A_diag);
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }
   m = n - nLU;
   e = nI - nLU;
   m_e = n - nI;
   if(e < 0)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU should not exceed nI.\n");
   }

   /* Init I array anyway. S's might be freed later */
   D_data = jxf_CTAlloc(JXF_Real, n);
   L_diag_i = jxf_CTAlloc(JXF_Int, (n+1));
   U_diag_i = jxf_CTAlloc(JXF_Int, (n+1));
   S_diag_i = jxf_CTAlloc(JXF_Int, (m+1));

   /* set Comm_Pkg if not yet built */
   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /*
    * 2: Symbolic factorization
    * setup iw and rperm first
    */
   /* allocate work arrays */
   iw = jxf_CTAlloc(JXF_Int, 4*n);
   rperm = iw + 3*n;
   L_diag_i[0] = U_diag_i[0] = S_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    */
   for(i=0; i<n; i++)
   {
      rperm[qperm[i]] = i;
   }

   /* do symbolic factorization */
   jxf_ILUSetupILUKSymbolic(n, A_diag_i, A_diag_j, lfil, perm, rperm, iw,
         nLU, L_diag_i, U_diag_i, S_diag_i, &L_diag_j, &U_diag_j, &S_diag_j, u_end);

   /*
    * after this, we have our I,J for L, U and S ready, and L sorted
    * iw are still -1 after symbolic factorization
    * now setup helper array here
    */
   if(L_diag_i[n])
   {
      L_diag_data = jxf_CTAlloc(JXF_Real, L_diag_i[n]);
   }
   if(U_diag_i[n])
   {
      U_diag_data = jxf_CTAlloc(JXF_Real, U_diag_i[n]);
   }
   if(S_diag_i[m])
   {
      S_diag_data = jxf_CTAlloc(JXF_Real, S_diag_i[m]);
   }

   /*
    * 3: Begin real factorization
    * we already have L and U structure ready, so no extra working array needed
    */
   /* first loop for upper part */
   for( ii = 0; ii < nLU; ii++ )
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii+1];
      ku = U_diag_i[ii+1];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      /* set up working arrays */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      D_data[ii] = 0.0;
      iw[ii] = ii;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, D and U */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if(col < ii)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else if(col == ii)
         {
            D_data[ii] = A_diag_data[j];
         }
         else
         {
            U_diag_data[icol] = A_diag_data[j];
         }
      }
      /* elimination */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv+1];

         for(k = U_diag_i[jpiv] ; k < ku ; k ++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if(icol < 0)
            {
               /* not in partern */
               continue;
            }
            if(col < ii)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else if(col == ii)
            {
               /* diag part */
               D_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else
            {
               /* U part */
               U_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
         }
      }
      /* reset working array */
      ku = U_diag_i[ii+1];
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      iw[ii] = -1;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = -1;
      }

      /* diagonal part (we store the inverse) */
      if(fabs(D_data[ii]) < MAT_TOL)
      {
         D_data[ii] = 1e-06;
      }
      D_data[ii] = 1./ D_data[ii];

   }

   /* Now lower part for Schur complement */
   for( ii = nLU; ii < n; ii++ )
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii+1];
      ku = S_diag_i[ii - nLU +1];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      /* set up working arrays */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      for(j = S_diag_i[ii - nLU] ; j < ku ; j ++)
      {
         col = S_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, and S */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if(col < nLU)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else
         {
            S_diag_data[icol] = A_diag_data[j];
         }
      }
      /* elimination */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv+1];
         for(k = U_diag_i[jpiv] ; k < ku ; k ++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if(icol < 0)
            {
               /* not in partern */
               continue;
            }
            if(col < nLU)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else
            {
               /* S part */
               S_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
         }
      }
      /* reset working array */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      ku = S_diag_i[ii-nLU+1];
      for(j = S_diag_i[ii-nLU] ; j < ku ; j ++)
      {
         col = S_diag_j[j];
         iw[col] = -1;
         /* remember to update index, S is smaller! */
         S_diag_j[j]-=nLU;
      }
   }

   /*
    * 4: Finishing up and free
    */

   /* First create Schur complement if necessary
    * Check if we need to create Schur complement
    */
   JXF_Int big_m = (JXF_Int)m;
   jxf_MPI_Allreduce(&big_m, &total_rows, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* only form when total_rows > 0 */
   if( total_rows > 0 )
   {
      /* now create S */
      /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
      {
         JXF_Int global_start;
         col_starts = jxf_CTAlloc(JXF_Int,2);
         jxf_MPI_Scan( &big_m, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
         col_starts[0] = global_start - m;
         col_starts[1] = global_start;
      }
#else
      col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

      jxf_MPI_Allgather(&big_m,1,JXF_MPI_INT,&col_starts[1],
            1,JXF_MPI_INT,comm);

      for (i=2; i < num_procs+1; i++)
         col_starts[i] += col_starts[i-1];
#endif

      /* We did nothing to A_offd, so all the data kept, just reorder them
       * The create function takes comm, global num rows/cols,
       *    row/col start, num cols offd, nnz diag, nnz offd
       */
      S_offd_nnz = jxf_CSRMatrixNumNonzeros(A_offd);
      S_offd_ncols = jxf_CSRMatrixNumCols(A_offd);

      matS = jxf_ParCSRMatrixCreate( comm,
            total_rows,
            total_rows,
            col_starts,
            col_starts,
            S_offd_ncols,
            S_diag_i[m],
            S_offd_nnz);

      /* S owns different start/end */
      jxf_ParCSRMatrixSetColStartsOwner(matS,1);
      jxf_ParCSRMatrixSetRowStartsOwner(matS,0);/* square matrix, use same row and col start */

      /* first put diagonal data in */
      S_diag = jxf_ParCSRMatrixDiag(matS);

      jxf_CSRMatrixI(S_diag) = S_diag_i;
      jxf_CSRMatrixData(S_diag) = S_diag_data;
      jxf_CSRMatrixJ(S_diag) = S_diag_j;

      /* now start to construct offdiag of S */
      S_offd = jxf_ParCSRMatrixOffd(matS);
      S_offd_i = jxf_TAlloc(JXF_Int, m+1);
      S_offd_j = jxf_TAlloc(JXF_Int, S_offd_nnz);
      S_offd_data = jxf_TAlloc(JXF_Real, S_offd_nnz);
      S_offd_colmap = jxf_CTAlloc(JXF_Int, S_offd_ncols);

      /* simply use a loop to copy data from A_offd */
      S_offd_i[0] = 0;
      k3 = 0;
      for(i = 1 ; i <= e ; i ++)
      {
         S_offd_i[i+1] = k3;
      }
      for(i = 0 ; i < m_e ; i ++)
      {
         col = perm[i + nI];
         k1 = A_offd_i[col];
         k2 = A_offd_i[col+1];
         for(j = k1 ; j < k2 ; j ++)
         {
            S_offd_j[k3] = A_offd_j[j];
            S_offd_data[k3++] = A_offd_data[j];
         }
         S_offd_i[i+e+1] = k3;
      }

      /* give I, J, DATA to S_offd */
      jxf_CSRMatrixI(S_offd) = S_offd_i;
      jxf_CSRMatrixJ(S_offd) = S_offd_j;
      jxf_CSRMatrixData(S_offd) = S_offd_data;

      /* now we need to update S_offd_colmap */

      /* get total num of send */
      JXF_Int num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      JXF_Int begin = jxf_ParCSRCommPkgSendMapStart(comm_pkg,0);
      JXF_Int end = jxf_ParCSRCommPkgSendMapStart(comm_pkg,num_sends);
      send_buf = jxf_TAlloc(JXF_Int, end - begin);
      /* copy new index into send_buf */
#ifdef JXF_NO_GLOBAL_PARTITION
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[0];
      }
#else
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[my_id];
      }
#endif

      /* main communication */
      comm_handle = jxf_ParCSRCommHandleCreate(22, comm_pkg, send_buf, S_offd_colmap);
      jxf_ParCSRCommHandleDestroy(comm_handle);

      /* setup index */
      jxf_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

      jxf_ILUSortOffdColmap(matS);

      /* free */
      jxf_TFree(send_buf);
   }/* end of forming S */

   /* Assemble LDU matrices */
   /* zero out unfactored rows */
   for(k=nLU; k<n; k++)
   {
      D_data[k] = 1.;
   }

   matL = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0 /* num_cols_offd */,
         L_diag_i[n],
         0 /* num_nonzeros_offd */);

   /* Have A own coarse_partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[n]>0)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jxf_TFree(L_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (L_diag_i[n]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         U_diag_i[n],
         0 );

   /* Have A own coarse_partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);

   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[n]>0)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jxf_TFree(U_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (U_diag_i[n]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free */
   jxf_TFree(iw);
   if(!matS)
   {
      /* we allocate some memory for S, need to free if unused */
      jxf_TFree(S_diag_i);
      //      jxf_TFree(col_starts);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;

   return jxf_error_flag;
}

/* ILUT
 * A: input matrix
 * lfil: maximum nnz per row in L and U
 * tol: droptol array in ILUT
 *    tol[0]: matrix B
 *    tol[1]: matrix E and F
 *    tol[2]: matrix S
 * perm: permutation array indicating ordering of factorization. Perm could come from a
 *    CF_marker: array or a reordering routine.
 * qperm: permutation array for column
 * nLU: size of computed LDU factorization. If nLU < n, Schur compelemnt will be formed
 * nI: number of interial unknowns. nLU should obey nLU <= nI.
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sptr: Schur complement
 *
 * Keep the largest lfil entries that is greater than some tol relative
 *    to the input tol and the norm of that row in both L and U
 */
JXF_Int
jxf_ILUSetupILUT(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Real *tol,
      JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI, jxf_ParCSRMatrix **Lptr,
      JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end)
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
   JXF_Real               local_nnz, total_nnz;
   JXF_Int                i, ii, j, k, k1, k2, k3, kl, ku, col, icol, lenl, lenu, lenhu, lenhlr, lenhll, jpos, jrow;
   JXF_Real               inorm, itolb, itolef, itols, dpiv, lxu;
   JXF_Int                *iw,*iL;
   JXF_Real               *w;

   /* memory management */
   JXF_Int                ctrL;
   JXF_Int                ctrU;
   JXF_Int                initial_alloc = 0;
   JXF_Int                capacity_L;
   JXF_Int                capacity_U;
   JXF_Int                ctrS;
   JXF_Int                capacity_S = 0;
   JXF_Int                nnz_A;

   /* communication stuffs for S */
   MPI_Comm                 comm             = jxf_ParCSRMatrixComm(A);
   JXF_Int                S_offd_nnz, S_offd_ncols;
   jxf_ParCSRCommPkg      *comm_pkg;
   jxf_ParCSRCommHandle   *comm_handle;
   JXF_Int                num_procs, my_id;
   JXF_Int                *col_starts;
   JXF_Int                total_rows;
   JXF_Int                num_sends;
   JXF_Int                begin, end;

   /* data objects for A */
   jxf_CSRMatrix          *A_diag          = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix          *A_offd          = jxf_ParCSRMatrixOffd(A);
   JXF_Real               *A_diag_data     = jxf_CSRMatrixData(A_diag);
   JXF_Int                *A_diag_i        = jxf_CSRMatrixI(A_diag);
   JXF_Int                *A_diag_j        = jxf_CSRMatrixJ(A_diag);
   JXF_Int                *A_offd_i        = jxf_CSRMatrixI(A_offd);
   JXF_Int                *A_offd_j        = jxf_CSRMatrixJ(A_offd);
   JXF_Real               *A_offd_data     = jxf_CSRMatrixData(A_offd);

   /* data objects for L, D, U */
   jxf_ParCSRMatrix       *matL;
   jxf_ParCSRMatrix       *matU;
   jxf_CSRMatrix          *L_diag;
   jxf_CSRMatrix          *U_diag;
   JXF_Real               *D_data;
   JXF_Real               *L_diag_data     = NULL;
   JXF_Int                *L_diag_i;
   JXF_Int                *L_diag_j        = NULL;
   JXF_Real               *U_diag_data     = NULL;
   JXF_Int                *U_diag_i;
   JXF_Int                *U_diag_j        = NULL;

   /* data objects for S */
   jxf_ParCSRMatrix       *matS            = NULL;
   jxf_CSRMatrix          *S_diag;
   jxf_CSRMatrix          *S_offd;
   JXF_Real               *S_diag_data     = NULL;
   JXF_Int                *S_diag_i        = NULL;
   JXF_Int                *S_diag_j        = NULL;
   JXF_Int                *S_offd_i        = NULL;
   JXF_Int                *S_offd_j        = NULL;
   JXF_Int                *S_offd_colmap   = NULL;
   JXF_Real               *S_offd_data;
   JXF_Int                *send_buf        = NULL;
   JXF_Int                *u_end_array;
   JXF_Int                u_end_location;

   /* reverse permutation */
   JXF_Int                *rperm;

   /* problem size
    * m is n - nLU, num of rows of local Schur system
    * m_e is the size of interface nodes
    * e is the number of interial rows in local Schur Complement
    */
   JXF_Int                n;
   JXF_Int                m;
   JXF_Int                e;
   JXF_Int                m_e;

   /* start setup
    * check input first
    */
   n = jxf_CSRMatrixNumRows(A_diag);
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }
   m = n - nLU;
   e = nI - nLU;
   m_e = n - nI;
   if(e < 0)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU should not exceed nI.\n");
   }

   u_end_array = jxf_TAlloc(JXF_Int, nLU);

   /* start set up
    * setup communication stuffs first
    */
   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   /* create if not yet built */
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* setup initial memory, in ILUT, just guess with max nnz per row */
   nnz_A = A_diag_i[nLU];
   if(n > 0)
   {
      initial_alloc = jxf_min(nLU + ceil((nnz_A / 2.0) * nLU / n), nLU * lfil);
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   D_data = jxf_CTAlloc(JXF_Real, n);
   L_diag_i = jxf_CTAlloc(JXF_Int, (n+1));
   U_diag_i = jxf_CTAlloc(JXF_Int, (n+1));

   L_diag_j = jxf_CTAlloc(JXF_Int, capacity_L);
   U_diag_j = jxf_CTAlloc(JXF_Int, capacity_U);
   L_diag_data = jxf_CTAlloc(JXF_Real, capacity_L);
   U_diag_data = jxf_CTAlloc(JXF_Real, capacity_U);

   ctrL = ctrU = 0;

   ctrS = 0;
   S_diag_i = jxf_CTAlloc(JXF_Int, (m + 1));
   S_diag_i[0] = 0;
   /* only setup S part when n > nLU */
   if(m > 0)
   {
      capacity_S = jxf_min(m + ceil((nnz_A / 2.0) * m / n), m * lfil);
      S_diag_j = jxf_CTAlloc(JXF_Int, capacity_S);
      S_diag_data = jxf_CTAlloc(JXF_Real, capacity_S);
   }

   /* setting up working array */
   iw = jxf_CTAlloc(JXF_Int,3*n);
   iL = iw + n;
   w = jxf_CTAlloc(JXF_Real,n);
   for(i = 0 ; i < n ; i ++)
   {
      iw[i] = -1;
   }
   L_diag_i[0] = U_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    * rperm[old] -> new
    * perm[new]  -> old
    */
   rperm = iw + 2*n;
   for(i = 0 ; i < n ; i ++)
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
   for(ii = 0 ; ii < nLU ; ii ++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      kl = ii-1;
      /* reset row norm of ith row */
      inorm = .0;
      for(j = k1 ; j < k2 ; j ++)
      {
         inorm += fabs(A_diag_data[j]);
      }
      if(inorm == .0)
      {
         jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: ILUT with zero row.\n");
      }
      inorm /= (JXF_Real)(k2-k1);
      /* set the scaled tol for that row */
      itolb = tol[0] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      w[ii] = 0.0;
      iw[ii] = ii;
      /* copy in data from A */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
         }
         else if(col == ii)
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
      while(lenhll > 0)
      {

         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jxf_ILUMinHeapRemoveIRIi(iL,w,iw,lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jxf_swap2(iL,w,lenhll,kl-lenhlr);
         lenhlr++;
         jxf_ILUMaxrHeapAddRabsI(w+kl,iL+kl,lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow+1];
         for(j = U_diag_i[jrow] ; j < ku ; j ++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = - dpiv*U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if( icol == -1 && ( (col < nLU && fabs(lxu) < itolb) || (col >= nLU && fabs(lxu) < itolef) ) )
            {
               continue;
            }
            if(icol == -1)
            {
               if(col < ii)
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
                  jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
               }
               else if(col == ii)
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
      }/* while loop for the elimination of current row */

      if(fabs(w[ii]) < MAT_TOL)
      {
         w[ii]=1e-06;
      }
      D_data[ii] = 1./w[ii];
      iw[ii] = -1;

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* test if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         ctrL += lenl;
         /* copy large data in */
         for(j = L_diag_i[ii] ; j < ctrL ; j ++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jxf_ILUMaxrHeapRemoveRabsI(w+kl,iL+kl,lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only U part
       */
      ku = lenu+ii;
      for(j = ii + 1 ; j <= ku ; j ++)
      {
         iw[iL[j]] = -1;
      }

      if(lenu < lfil)
      {
         /* we simply keep all of the data, no need to sort */
         lenhu = lenu;
      }
      else
      {
         /* need to sort the first small(hopefully) part of it */
         lenhu = lfil;
         /* quick split, only sort the first small part of the array */
         jxf_ILUMaxQSplitRabsI(w,iL,ii+1,ii+lenhu,ii+lenu);
      }

      U_diag_i[ii+1] = U_diag_i[ii] + lenhu;
      if(lenhu > 0)
      {
         /* test if memory is enough */
         while(ctrU + lenhu > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
         }
         ctrU += lenhu;
         /* copy large data in */
         for(j = U_diag_i[ii] ; j < ctrU ; j ++)
         {
            jpos = ii+1+j-U_diag_i[ii];
            U_diag_j[j] = iL[jpos];
            U_diag_data[j] = w[jpos];
         }
      }
      /* check and build u_end array */
      if(m > 0)
      {
         jxf_qsort1(U_diag_j,U_diag_data,U_diag_i[ii],U_diag_i[ii+1]-1);
         u_end_location = jxf_BinarySearch2(U_diag_j,nLU,U_diag_i[ii],U_diag_i[ii+1]-1,u_end_array + ii);
         if(u_end_location >= 0)
         {
            u_end_array[ii] = u_end_location + 1;
         }
      }
      else
      {
         /* Everything is in U */
         u_end_array[ii] = ctrU;
      }
   }/* end of ii loop from 0 to nLU-1 */


   /* now main loop for Schur comlement part */
   for(ii = nLU ; ii < n ; ii ++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      kl = nLU-1;
      /* reset row norm of ith row */
      inorm = .0;
      for(j = k1 ; j < k2 ; j ++)
      {
         inorm += fabs(A_diag_data[j]);
      }
      if(inorm == .0)
      {
         jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: ILUT with zero row.\n");
      }
      inorm /= (JXF_Real)(k2-k1);
      /* set the scaled tol for that row */
      itols = tol[2] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      /* copy in data from A */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if(col < nLU)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
         }
         else if(col == ii)
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
      while(lenhll > 0)
      {
         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jxf_ILUMinHeapRemoveIRIi(iL,w,iw,lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jxf_swap2(iL,w,lenhll,kl-lenhlr);
         lenhlr++;
         jxf_ILUMaxrHeapAddRabsI(w+kl,iL+kl,lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow+1];
         for(j = U_diag_i[jrow] ; j < ku ; j ++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = - dpiv*U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if(icol == -1 && ( (col < nLU && fabs(lxu) < itolef) || ( col >= nLU && fabs(lxu) < itols ) ) )
            {
               continue;
            }
            if(icol == -1)
            {
               if(col < nLU)
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
                  jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
               }
               else if(col == ii)
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
      }/* while loop for the elimination of current row */

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* test if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         ctrL += lenl;
         /* copy large data in */
         for(j = L_diag_i[ii] ; j < ctrL ; j ++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jxf_ILUMaxrHeapRemoveRabsI(w+kl,iL+kl,lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only S part
       */
      ku = lenu+nLU;
      for(j = nLU ; j <= ku ; j ++)
      {
         iw[iL[j]] = -1;
      }

      /* no dropping at this point of time for S */
      //lenhu = lenu < lfil ? lenu : lfil;
      lenhu = lenu;
      /* quick split, only sort the first small part of the array */
      jxf_ILUMaxQSplitRabsI(w,iL,nLU+1,nLU+lenhu,nLU+lenu);
      /* we have diagonal in S anyway */
      /* test if memory is enough */
      while(ctrS + lenhu + 1 > capacity_S)
      {
         JXF_Int tmp = capacity_S;
         capacity_S = capacity_S * EXPAND_FACT + 1;
         S_diag_j = jxf_TReAlloc_v2(S_diag_j, JXF_Int, tmp, JXF_Int, capacity_S);
         S_diag_data = jxf_TReAlloc_v2(S_diag_data, JXF_Real, tmp, JXF_Real, capacity_S);
      }

      ctrS += (lenhu+1);
      S_diag_i[ii-nLU+1] = ctrS;

      /* copy large data in, diagonal first */
      S_diag_j[S_diag_i[ii-nLU]] = iL[nLU]-nLU;
      S_diag_data[S_diag_i[ii-nLU]] = w[nLU];
      for(j = S_diag_i[ii-nLU] + 1 ; j < ctrS ; j ++)
      {
         jpos = nLU+j-S_diag_i[ii-nLU];
         S_diag_j[j] = iL[jpos]-nLU;
         S_diag_data[j] = w[jpos];
      }
   }/* end of ii loop from nLU to n-1 */

   /*
    * 3: Finishing up and free
    */

   /* First create Schur complement if necessary
    * Check if we need to create Schur complement
    */
   JXF_Int big_m = (JXF_Int)m;
   jxf_MPI_Allreduce(&big_m, &total_rows, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* only form when total_rows > 0 */
   if( total_rows > 0 )
   {
      /* now create S */
      /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
      {
         JXF_Int global_start;
         col_starts = jxf_CTAlloc(JXF_Int,2);
         jxf_MPI_Scan( &big_m, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
         col_starts[0] = global_start - m;
         col_starts[1] = global_start;
      }
#else
      col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

      jxf_MPI_Allgather(&big_m,1,JXF_MPI_INT,&col_starts[1],
            1,JXF_MPI_INT,comm);

      for (i=2; i < num_procs+1; i++)
         col_starts[i] += col_starts[i-1];
#endif
      /* We did nothing to A_offd, so all the data kept, just reorder them
       * The create function takes comm, global num rows/cols,
       *    row/col start, num cols offd, nnz diag, nnz offd
       */
      S_offd_nnz = jxf_CSRMatrixNumNonzeros(A_offd);
      S_offd_ncols = jxf_CSRMatrixNumCols(A_offd);

      matS = jxf_ParCSRMatrixCreate( comm,
            total_rows,
            total_rows,
            col_starts,
            col_starts,
            S_offd_ncols,
            S_diag_i[m],
            S_offd_nnz);

      /* S owns different start/end */
      jxf_ParCSRMatrixSetColStartsOwner(matS,1);
      jxf_ParCSRMatrixSetRowStartsOwner(matS,0);/* square matrix, use same row and col start */

      /* first put diagonal data in */
      S_diag = jxf_ParCSRMatrixDiag(matS);

      jxf_CSRMatrixI(S_diag) = S_diag_i;
      jxf_CSRMatrixData(S_diag) = S_diag_data;
      jxf_CSRMatrixJ(S_diag) = S_diag_j;

      /* now start to construct offdiag of S */
      S_offd = jxf_ParCSRMatrixOffd(matS);
      S_offd_i = jxf_TAlloc(JXF_Int, m+1);
      S_offd_j = jxf_TAlloc(JXF_Int, S_offd_nnz);
      S_offd_data = jxf_TAlloc(JXF_Real, S_offd_nnz);
      S_offd_colmap = jxf_CTAlloc(JXF_Int, S_offd_ncols);

      /* simply use a loop to copy data from A_offd */
      S_offd_i[0] = 0;
      k3 = 0;
      for(i = 1 ; i <= e ; i ++)
      {
         S_offd_i[i] = k3;
      }
      for(i = 0 ; i < m_e ; i ++)
      {
         col = perm[i + nI];
         k1 = A_offd_i[col];
         k2 = A_offd_i[col+1];
         for(j = k1 ; j < k2 ; j ++)
         {
            S_offd_j[k3] = A_offd_j[j];
            S_offd_data[k3++] = A_offd_data[j];
         }
         S_offd_i[i+e+1] = k3;
      }

      /* give I, J, DATA to S_offd */
      jxf_CSRMatrixI(S_offd) = S_offd_i;
      jxf_CSRMatrixJ(S_offd) = S_offd_j;
      jxf_CSRMatrixData(S_offd) = S_offd_data;

      /* now we need to update S_offd_colmap */

      /* get total num of send */
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      begin = jxf_ParCSRCommPkgSendMapStart(comm_pkg,0);
      end = jxf_ParCSRCommPkgSendMapStart(comm_pkg,num_sends);
      send_buf = jxf_TAlloc(JXF_Int, end - begin);
      /* copy new index into send_buf */
#ifdef JXF_NO_GLOBAL_PARTITION
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[0];
      }
#else
      for(i = begin ; i < end ; i ++)
      {
         send_buf[i-begin] = rperm[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,i)] - nLU + col_starts[my_id];
      }
#endif

      /* main communication */
      comm_handle = jxf_ParCSRCommHandleCreate(22, comm_pkg, send_buf, S_offd_colmap);
      /* need this to synchronize, Isend & Irecv used in above functions */
      jxf_ParCSRCommHandleDestroy(comm_handle);

      /* setup index */
      jxf_ParCSRMatrixColMapOffd(matS) = S_offd_colmap;

      jxf_ILUSortOffdColmap(matS);

      /* free */
      jxf_TFree(send_buf);
   }/* end of forming S */

   /* now start to construct L and U */
   for(k=nLU; k<n; k++)
   {
      /* set U after nLU to be 0, and diag to be one */
      U_diag_i[k+1] = U_diag_i[nLU];
      D_data[k] = 1.;
   }

   /* create parcsr matrix */
   matL = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         L_diag_i[n],
         0 );

   /* Have A own coarse_partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[n] > 0)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jxf_TFree(L_diag_j);
      jxf_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (L_diag_i[n]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         U_diag_i[n],
         0 );

   /* Have A own coarse_partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);

   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[n] > 0)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jxf_TFree(U_diag_j);
      jxf_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (U_diag_i[n]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free working array */
   jxf_TFree(iw);
   jxf_TFree(w);

   if(!matS)
   {
      jxf_TFree(S_diag_i);
      //      jxf_TFree(col_starts);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;
   *Sptr = matS;
   *u_end = u_end_array;

   return jxf_error_flag;
}


/* NSH setup */
/* Setup NSH data */
JXF_Int
jxf_NSHSetup( void               *nsh_vdata,
                  jxf_ParCSRMatrix *A,
                  jxf_ParVector    *f,
                  jxf_ParVector    *u )
{
   MPI_Comm             comm              = jxf_ParCSRMatrixComm(A);
   jxf_ParNSHData     *nsh_data         = (jxf_ParNSHData*) nsh_vdata;

   //   JXF_Int            i;
   // JXF_Int            num_threads;
   // JXF_Int            debug_flag = 0;

   /* pointers to NSH data */
   JXF_Int            logging           = jxf_ParNSHDataLogging(nsh_data);
   JXF_Int            print_level       = jxf_ParNSHDataPrintLevel(nsh_data);

   jxf_ParCSRMatrix   *matA             = jxf_ParNSHDataMatA(nsh_data);
   jxf_ParCSRMatrix   *matM             = jxf_ParNSHDataMatM(nsh_data);

   //   JXF_Int            n                 = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A));
   JXF_Int            num_procs,  my_id;

   jxf_ParVector      *Utemp;
   jxf_ParVector      *Ftemp;
   jxf_ParVector      *F_array          = jxf_ParNSHDataF(nsh_data);
   jxf_ParVector      *U_array          = jxf_ParNSHDataU(nsh_data);
   jxf_ParVector      *residual         = jxf_ParNSHDataResidual(nsh_data);
   JXF_Real           *rel_res_norms    = jxf_ParNSHDataRelResNorms(nsh_data);

   /* solver setting */
   JXF_Real           *droptol          = jxf_ParNSHDataDroptol(nsh_data);
   JXF_Real           mr_tol            = jxf_ParNSHDataMRTol(nsh_data);
   JXF_Int            mr_max_row_nnz    = jxf_ParNSHDataMRMaxRowNnz(nsh_data);
   JXF_Int            mr_max_iter       = jxf_ParNSHDataMRMaxIter(nsh_data);
   JXF_Int            mr_col_version    = jxf_ParNSHDataMRColVersion(nsh_data);
   JXF_Real           nsh_tol           = jxf_ParNSHDataNSHTol(nsh_data);
   JXF_Int            nsh_max_row_nnz   = jxf_ParNSHDataNSHMaxRowNnz(nsh_data);
   JXF_Int            nsh_max_iter      = jxf_ParNSHDataNSHMaxIter(nsh_data);

   /* ----- begin -----*/

   //num_threads = jxf_NumThreads();

   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);

   /* Free Previously allocated data, if any not destroyed */
   if(matM)
   {
      jxf_TFree(matM);
      matM = NULL;
   }

   /* clear old l1_norm data, if created */
   if(jxf_ParNSHDataL1Norms(nsh_data))
   {
      jxf_TFree(jxf_ParNSHDataL1Norms(nsh_data));
      jxf_ParNSHDataL1Norms(nsh_data) = NULL;
   }

   /* setup temporary storage
    * first check is they've already here
    */
   if (jxf_ParNSHDataUTemp(nsh_data))
   {
      jxf_ParVectorDestroy(jxf_ParNSHDataUTemp(nsh_data));
      jxf_ParNSHDataUTemp(nsh_data) = NULL;
   }
   if (jxf_ParNSHDataFTemp(nsh_data))
   {
      jxf_ParVectorDestroy(jxf_ParNSHDataFTemp(nsh_data));
      jxf_ParNSHDataFTemp(nsh_data) = NULL;
   }
   if (jxf_ParNSHDataResidual(nsh_data))
   {
      jxf_ParVectorDestroy(jxf_ParNSHDataResidual(nsh_data));
      jxf_ParNSHDataResidual(nsh_data) = NULL;
   }
   if (jxf_ParNSHDataRelResNorms(nsh_data))
   {
      jxf_TFree(jxf_ParNSHDataRelResNorms(nsh_data));
      jxf_ParNSHDataRelResNorms(nsh_data) = NULL;
   }

   /* start to create working vectors */
   Utemp = jxf_ParVectorCreate(jxf_ParCSRMatrixComm(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A));
   jxf_ParVectorInitialize(Utemp);
   jxf_ParVectorSetPartitioningOwner(Utemp,0);
   jxf_ParNSHDataUTemp(nsh_data) = Utemp;

   Ftemp = jxf_ParVectorCreate(jxf_ParCSRMatrixComm(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A));
   jxf_ParVectorInitialize(Ftemp);
   jxf_ParVectorSetPartitioningOwner(Ftemp,0);
   jxf_ParNSHDataFTemp(nsh_data) = Ftemp;
   /* set matrix, solution and rhs pointers */
   matA = A;
   F_array = f;
   U_array = u;

   /* NSH compute approximate inverse, see par_ilu.c */
   jxf_ILUParCSRInverseNSH(matA, &matM, droptol, mr_tol, nsh_tol, DIVIDE_TOL, mr_max_row_nnz,
         nsh_max_row_nnz, mr_max_iter, nsh_max_iter, mr_col_version, print_level);

   /* set pointers to NSH data */
   jxf_ParNSHDataMatA(nsh_data) = matA;
   jxf_ParNSHDataF(nsh_data) = F_array;
   jxf_ParNSHDataU(nsh_data) = U_array;
   jxf_ParNSHDataMatM(nsh_data) = matM;

   /* compute operator complexity */
   jxf_ParCSRMatrixSetDNumNonzeros(matA);
   jxf_ParCSRMatrixSetDNumNonzeros(matM);
   /* compute complexity */
   jxf_ParNSHDataOperatorComplexity(nsh_data) =  jxf_ParCSRMatrixDNumNonzeros(matM)/jxf_ParCSRMatrixDNumNonzeros(matA);
   if (my_id == 0)
   {
      jxf_printf("NSH SETUP: operator complexity = %f  \n", jxf_ParNSHDataOperatorComplexity(nsh_data));
   }

   if ( logging > 1 ) {
      residual =
         jxf_ParVectorCreate(jxf_ParCSRMatrixComm(matA),
               jxf_ParCSRMatrixGlobalNumRows(matA),
               jxf_ParCSRMatrixRowStarts(matA) );
      jxf_ParVectorInitialize(residual);
      jxf_ParVectorSetPartitioningOwner(residual,0);
      jxf_ParNSHDataResidual(nsh_data)= residual;
   }
   else{
      jxf_ParNSHDataResidual(nsh_data) = NULL;
   }
   rel_res_norms = jxf_CTAlloc(JXF_Real, jxf_ParNSHDataMaxIter(nsh_data));
   jxf_ParNSHDataRelResNorms(nsh_data) = rel_res_norms;

   return jxf_error_flag;
}


/* ILU(0) for RAS, has some external rows
 * A = input matrix
 * perm = permutation array indicating ordering of factorization. Perm could come from a
 *    CF_marker array or a reordering routine.
 * nLU = size of computed LDU factorization.
 * Lptr, Dptr, Uptr, Sptr = L, D, U, S factors.
 * will form global Schur Matrix if nLU < n
 */
JXF_Int
jxf_ILUSetupILU0RAS(jxf_ParCSRMatrix *A, JXF_Int *perm, JXF_Int nLU,
      jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr)
{
   JXF_Int                i, ii, j, k, k1, k2, ctrU, ctrL, lenl, lenu, jpiv, col, jpos;
   JXF_Int                *iw, *iL, *iU;
   JXF_Real               dd, t, dpiv, lxu, *wU, *wL;

   /* communication stuffs for S */
   MPI_Comm                 comm          = jxf_ParCSRMatrixComm(A);
   JXF_Int                num_procs;
   //   JXF_Int                S_offd_nnz, S_offd_ncols;
   jxf_ParCSRCommPkg      *comm_pkg;
   //   jxf_ParCSRCommHandle   *comm_handle;
   //   JXF_Int                num_sends, begin, end;
   //   JXF_Int                *send_buf     = NULL;

   /* data objects for A */
   jxf_CSRMatrix          *A_diag       = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix          *A_offd       = jxf_ParCSRMatrixOffd(A);
   JXF_Real               *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int                *A_diag_i     = jxf_CSRMatrixI(A_diag);
   JXF_Int                *A_diag_j     = jxf_CSRMatrixJ(A_diag);
   JXF_Real               *A_offd_data  = jxf_CSRMatrixData(A_offd);
   JXF_Int                *A_offd_i     = jxf_CSRMatrixI(A_offd);
   JXF_Int                *A_offd_j     = jxf_CSRMatrixJ(A_offd);

   /* size of problem and external matrix */
   JXF_Int                n             =  jxf_CSRMatrixNumRows(A_diag);
   //   JXF_Int                m             = n - nLU;
   JXF_Int                ext           = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int                total_rows    = n + ext;
   JXF_Int             *col_starts;
   JXF_Int             global_num_rows;
   JXF_Real               local_nnz, total_nnz;

   /* data objects for L, D, U */
   jxf_ParCSRMatrix       *matL;
   jxf_ParCSRMatrix       *matU;
   jxf_CSRMatrix          *L_diag;
   jxf_CSRMatrix          *U_diag;
   JXF_Real               *D_data;
   JXF_Real               *L_diag_data;
   JXF_Int                *L_diag_i;
   JXF_Int                *L_diag_j;
   JXF_Real               *U_diag_data;
   JXF_Int                *U_diag_i;
   JXF_Int                *U_diag_j;

   /* data objects for E, external matrix */
   JXF_Int                *E_i;
   JXF_Int                *E_j;
   JXF_Real               *E_data;

   /* memory management */
   JXF_Int                initial_alloc = 0;
   JXF_Int                capacity_L;
   JXF_Int                capacity_U;
   JXF_Int                nnz_A = A_diag_i[n];

   /* reverse permutation array */
   JXF_Int                *rperm;
   /* the original permutation array */
   JXF_Int                *perm_old;

   /* start setup
    * get communication stuffs first
    */
   jxf_MPI_Comm_size(comm,&num_procs);
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   /* setup if not yet built */
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* check for correctness */
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }

   /* Allocate memory for L,D,U,S factors */
   if(n > 0)
   {
      initial_alloc = (n + ext) + ceil((nnz_A / 2.0)*total_rows/n);
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   D_data      = jxf_TAlloc(JXF_Real, total_rows);
   L_diag_i    = jxf_TAlloc(JXF_Int, total_rows+1);
   L_diag_j    = jxf_TAlloc(JXF_Int, capacity_L);
   L_diag_data = jxf_TAlloc(JXF_Real, capacity_L);
   U_diag_i    = jxf_TAlloc(JXF_Int, total_rows+1);
   U_diag_j    = jxf_TAlloc(JXF_Int, capacity_U);
   U_diag_data = jxf_TAlloc(JXF_Real, capacity_U);

   /* allocate working arrays */
   iw          = jxf_TAlloc(JXF_Int, 4*total_rows);
   iL          = iw+total_rows;
   rperm       = iw + 2 * total_rows;
   perm_old    = perm;
   perm        = iw + 3 * total_rows;
   wL          = jxf_TAlloc(JXF_Real, total_rows);
   ctrU = ctrL = 0;
   L_diag_i[0] = U_diag_i[0] = 0;
   /* set marker array iw to -1 */
   for( i = 0 ; i < total_rows ; i++ )
   {
      iw[i] = -1;
   }

   /* expand perm to suit extra data, remember to free */
   for( i = 0 ; i < n ; i ++)
   {
      perm[i] = perm_old[i];
   }
   for( i = n ; i < total_rows ; i ++)
   {
      perm[i] = i;
   }

   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    */
   for(i=0 ; i < total_rows ; i++)
   {
      rperm[perm[i]] = i;
   }

   /* get external rows */
   jxf_ILUBuildRASExternalMatrix(A, rperm, &E_i, &E_j, &E_data);

   /*---------  Begin Factorization. Work in permuted space  ----
    * this is the first part, without offd
    */
   for( ii = 0; ii < nLU; ii++ )
   {

      // get row i
      i = perm[ii];
      // get extents of row i
      k1=A_diag_i[i];
      k2=A_diag_i[i+1];

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL+ii;
      wU = wL+ii;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl  = lenu = 0;
      iw[ii] = ii;
      /*-------------------- scan & unwrap column */
      for(j=k1; j < k2; j++)
      {
         col = rperm[A_diag_j[j]];
         t = A_diag_data[j];
         if( col < ii )
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
            dd=t;
         }
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jxf_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jxf_qsort3ir(iL, wL, iw, 0, (lenl-1));
      for(j=0; j<lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for(k=U_diag_i[jpiv]; k<U_diag_i[jpiv+1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if(jpos < 0)
            {
               continue;
            }

            lxu = - U_diag_data[k] * dpiv;
            if(col < ii)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if(col > ii)
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
      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for( j = 0; j < lenu; j++ )
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      while((ctrL+lenl) > capacity_L)
      {
         JXF_Int tmp = capacity_L;
         capacity_L = capacity_L * EXPAND_FACT + 1;
         L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
      }
      jxf_TMemcpy(&(L_diag_j)[ctrL], iL, JXF_Int, lenl);
      jxf_TMemcpy(&(L_diag_data)[ctrL], wL, JXF_Real, lenl);
      L_diag_i[ii+1] = (ctrL+=lenl);

      /* diagonal part (we store the inverse) */
      if(fabs(dd) < MAT_TOL)
      {
         dd = 1.0e-6;
      }
      D_data[ii] = 1./dd;

      /* U part */
      /* Check that memory is sufficient */
      while((ctrU+lenu) > capacity_U)
      {
         JXF_Int tmp = capacity_U;
         capacity_U = capacity_U * EXPAND_FACT + 1;
         U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
         U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
      }
      jxf_TMemcpy(&(U_diag_j)[ctrU], iU, JXF_Int, lenu);
      jxf_TMemcpy(&(U_diag_data)[ctrU], wU, JXF_Real, lenu);
      U_diag_i[ii+1] = (ctrU+=lenu);

   }
   /*---------  Begin Factorization in lower part  ----
    * here we need to get off diagonals in
    */
   for( ii = nLU ; ii < n ; ii++ )
   {

      // get row i
      i = perm[ii];
      // get extents of row i
      k1=A_diag_i[i];
      k2=A_diag_i[i+1];

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL+ii;
      wU = wL+ii;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl  = lenu = 0;
      iw[ii] = ii;
      /*-------------------- scan & unwrap column */
      for(j=k1; j < k2; j++)
      {
         col = rperm[A_diag_j[j]];
         t = A_diag_data[j];
         if( col < ii )
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
            dd=t;
         }
      }

      /*------------------ sjcan offd*/
      k1=A_offd_i[i];
      k2=A_offd_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         /* add offd to U part, all offd are U for this part */
         col = A_offd_j[j] + n;
         t = A_offd_data[j];
         iw[col] = lenu;
         iU[lenu] = col;
         wU[lenu++] = t;
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jxf_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jxf_qsort3ir(iL, wL, iw, 0, (lenl-1));
      for(j=0; j<lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for(k=U_diag_i[jpiv]; k<U_diag_i[jpiv+1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if(jpos < 0)
            {
               continue;
            }

            lxu = - U_diag_data[k] * dpiv;
            if(col < ii)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if(col > ii)
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
      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for( j = 0; j < lenu; j++ )
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      while((ctrL+lenl) > capacity_L)
      {
         JXF_Int tmp = capacity_L;
         capacity_L = capacity_L * EXPAND_FACT + 1;
         L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
      }
      jxf_TMemcpy(&(L_diag_j)[ctrL], iL, JXF_Int, lenl);
      jxf_TMemcpy(&(L_diag_data)[ctrL], wL, JXF_Real, lenl);
      L_diag_i[ii+1] = (ctrL+=lenl);

      /* diagonal part (we store the inverse) */
      if(fabs(dd) < MAT_TOL)
      {
         dd = 1.0e-6;
      }
      D_data[ii] = 1./dd;

      /* U part */
      /* Check that memory is sufficient */
      while((ctrU+lenu) > capacity_U)
      {
         JXF_Int tmp = capacity_U;
         capacity_U = capacity_U * EXPAND_FACT + 1;
         U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
         U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
      }
      jxf_TMemcpy(&(U_diag_j)[ctrU], iU, JXF_Int, lenu);
      jxf_TMemcpy(&(U_diag_data)[ctrU], wU, JXF_Real, lenu);
      U_diag_i[ii+1] = (ctrU+=lenu);

   }

   /*---------  Begin Factorization in external part  ----
    * here we need to get off diagonals in
    */
   for( ii = n ; ii < total_rows ; ii++ )
   {

      // get row i
      i = ii-n;
      // get extents of row i
      k1=E_i[i];
      k2=E_i[i+1];

      /*-------------------- unpack L & U-parts of row of A in arrays w */
      iU = iL+ii;
      wU = wL+ii;
      /*--------------------  diagonal entry */
      dd = 0.0;
      lenl  = lenu = 0;
      iw[ii] = ii;
      /*-------------------- scan & unwrap column */
      for(j=k1; j < k2; j++)
      {
         col = rperm[E_j[j]];
         t = E_data[j];
         if( col < ii )
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
            dd=t;
         }
      }

      /* eliminate row */
      /*-------------------------------------------------------------------------
       *  In order to do the elimination in the correct order we must select the
       *  smallest column index among iL[k], k = j, j+1, ..., lenl-1. For ILU(0),
       *  no new fill-ins are expect, so we can pre-sort iL and wL prior to the
       *  entering the elimination loop.
       *-----------------------------------------------------------------------*/
      //      jxf_quickSortIR(iL, wL, iw, 0, (lenl-1));
      jxf_qsort3ir(iL, wL, iw, 0, (lenl-1));
      for(j=0; j<lenl; j++)
      {
         jpiv = iL[j];
         /* get factor/ pivot element */
         dpiv = wL[j] * D_data[jpiv];
         /* store entry in L */
         wL[j] = dpiv;

         /* zero out element - reset pivot */
         iw[jpiv] = -1;
         /* combine current row and pivot row */
         for(k=U_diag_i[jpiv]; k<U_diag_i[jpiv+1]; k++)
         {
            col = U_diag_j[k];
            jpos = iw[col];

            /* Only fill-in nonzero pattern (jpos != 0) */
            if(jpos < 0)
            {
               continue;
            }

            lxu = - U_diag_data[k] * dpiv;
            if(col < ii)
            {
               /* dealing with L part */
               wL[jpos] += lxu;
            }
            else if(col > ii)
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
      /* restore iw (only need to restore diagonal and U part */
      iw[ii] = -1;
      for( j = 0; j < lenu; j++ )
      {
         iw[iU[j]] = -1;
      }

      /* Update LDU factors */
      /* L part */
      /* Check that memory is sufficient */
      while((ctrL+lenl) > capacity_L)
      {
         JXF_Int tmp = capacity_L;
         capacity_L = capacity_L * EXPAND_FACT + 1;
         L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
      }
      jxf_TMemcpy(&(L_diag_j)[ctrL], iL, JXF_Int, lenl);
      jxf_TMemcpy(&(L_diag_data)[ctrL], wL, JXF_Real, lenl);
      L_diag_i[ii+1] = (ctrL+=lenl);

      /* diagonal part (we store the inverse) */
      if(fabs(dd) < MAT_TOL)
      {
         dd = 1.0e-6;
      }
      D_data[ii] = 1./dd;

      /* U part */
      /* Check that memory is sufficient */
      while((ctrU+lenu) > capacity_U)
      {
         JXF_Int tmp = capacity_U;
         capacity_U = capacity_U * EXPAND_FACT + 1;
         U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
         U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
      }
      jxf_TMemcpy(&(U_diag_j)[ctrU], iU, JXF_Int, lenu);
      jxf_TMemcpy(&(U_diag_data)[ctrU], wU, JXF_Real, lenu);
      U_diag_i[ii+1] = (ctrU+=lenu);

   }
   JXF_Int big_total_rows = (JXF_Int)total_rows;
   jxf_MPI_Allreduce( &big_total_rows, &global_num_rows, 1, JXF_MPI_INT, MPI_SUM, comm);

   /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
   {
      JXF_Int global_start;
      col_starts = jxf_CTAlloc(JXF_Int,2);
      jxf_MPI_Scan( &big_total_rows, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
      col_starts[0] = global_start - total_rows;
      col_starts[1] = global_start;
   }
#else
   col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

   jxf_MPI_Allgather(&big_total_rows,1,JXF_MPI_INT,&col_starts[1],
         1,JXF_MPI_INT,comm);

   for (i=2; i < num_procs+1; i++)
      col_starts[i] += col_starts[i-1];
#endif

   matL = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0,
         ctrL,
         0 );

   /* Have A own row/col partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,1);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (ctrL)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jxf_TFree(L_diag_j);
      jxf_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) ctrL;
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0,
         ctrU,
         0 );

   /* Have A own row/col partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);
   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (ctrU)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we've allocated some memory, so free if not used */
      jxf_TFree(U_diag_j);
      jxf_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) ctrU;
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;
   /* free memory */
   jxf_TFree(wL);
   jxf_TFree(iw);

   /* free external data */
   if(E_i)
   {
      jxf_TFree(E_i);
   }
   if(E_j)
   {
      jxf_TFree(E_j);
      jxf_TFree(E_data);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;

   return jxf_error_flag;
}



/* ILU(k) symbolic factorization for RAS
 * n = total rows of input
 * lfil = level of fill-in, the k in ILU(k)
 * perm = permutation array indicating ordering of factorization. Perm could come from a
 * rperm = reverse permutation array, used here to avoid duplicate memory allocation
 * iw = working array, used here to avoid duplicate memory allocation
 * nLU = size of computed LDU factorization.
 * A/L/U/E_i = the I slot of A, L, U and E
 * A/L/U/E_j = the J slot of A, L, U and E
 * will form global Schur Matrix if nLU < n
 */
JXF_Int
jxf_ILUSetupILUKRASSymbolic(JXF_Int n, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Int *A_offd_i, JXF_Int *A_offd_j,
                              JXF_Int *E_i, JXF_Int *E_j, JXF_Int ext,
                              JXF_Int lfil, JXF_Int *perm,
                              JXF_Int *rperm,   JXF_Int *iw,   JXF_Int nLU,
                              JXF_Int *L_diag_i, JXF_Int *U_diag_i,
                              JXF_Int **L_diag_j, JXF_Int **U_diag_j)
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

   JXF_Int      *temp_L_diag_j, *temp_U_diag_j, *u_levels;
   JXF_Int      *iL, *iLev;
   JXF_Int      ii, i, j, k, ku, lena, lenl, lenu, lenh, ilev, lev, col, icol;
   //   JXF_Int      m = n - nLU;
   JXF_Int      total_rows = ext + n;

   /* memory management */
   JXF_Int      ctrL;
   JXF_Int      ctrU;
   JXF_Int      capacity_L;
   JXF_Int      capacity_U;
   JXF_Int      initial_alloc = 0;
   JXF_Int      nnz_A;

   /* set iL and iLev to right place in iw array */
   iL             = iw + total_rows;
   iLev           = iw + 2*total_rows;

   /* setup initial memory used */
   nnz_A          = A_diag_i[n];
   if(n > 0)
   {
      initial_alloc  = (n + ext) + ceil((nnz_A / 2.0) * total_rows / n);
   }
   capacity_L     = initial_alloc;
   capacity_U     = initial_alloc;

   /* allocate other memory for L and U struct */
   temp_L_diag_j  = jxf_CTAlloc(JXF_Int, capacity_L);
   temp_U_diag_j  = jxf_CTAlloc(JXF_Int, capacity_U);

   u_levels       = jxf_CTAlloc(JXF_Int, capacity_U);
   ctrL = ctrU = 0;

   /* set initial value for working array */
   for(ii = 0 ; ii < total_rows ; ii ++)
   {
      iw[ii] = -1;
   }

   /*
    * 2: Start of main loop
    * those in iL are NEW col index (after permutation)
    */
   for(ii = 0 ; ii < nLU ; ii ++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0;/* this is the current length of heap */
      lenu = ii;
      lena = A_diag_i[i+1];
      /* put those already inside original pattern, and set their level to 0 */
      for(j = A_diag_i[i] ; j < lena ; j ++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
         }
         else if(col > ii)
         {
            /* this is an entry in U */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      }/* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while(lenh > 0)
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
         jxf_ILUMinHeapRemoveIIIi(iL,iLev,iw,lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k]=-1;
         jxf_swap2i(iL,iLev,ii-lenl,lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k+1];
         for(j = U_diag_i[k] ; j < ku ; j ++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if(lev > lfil)
            {
               continue;
            }
            if(icol < 0)
            {
               /* not yet in */
               if(col < ii)
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
                  jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
               }
               else if(col > ii)
               {
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jxf_min(lev, iLev[icol]);
            }
         }/* end of loop j for level update */
      }/* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* check if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jxf_TReAlloc_v2(temp_L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for(j = 0 ; j < lenl ; j ++)
         {
            temp_L_diag_j[ctrL+j] = iL[ii-j-1];
         }
         ctrL += lenl;
      }
      k = lenu - ii;
      U_diag_i[ii+1] = U_diag_i[ii] + k;
      if(k > 0)
      {
         /* check if memory is enough */
         while(ctrU + k > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            temp_U_diag_j = jxf_TReAlloc_v2(temp_U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            u_levels = jxf_TReAlloc_v2(u_levels, JXF_Int, tmp, JXF_Int, capacity_U);
         }
         jxf_TMemcpy(temp_U_diag_j+ctrU,iL+ii,JXF_Int,k);
         jxf_TMemcpy(u_levels+ctrU,iLev+ii,JXF_Int,k);
         ctrU += k;
      }

      /* reset iw */
      for(j = ii ; j < lenu ; j ++)
      {
         iw[iL[j]] = -1;
      }

   }/* end of main loop ii from 0 to nLU-1 */

   /*
    * Offd part
    */
   for(ii = nLU ; ii < n ; ii ++)
   {
      i = perm[ii];
      lenl = 0;
      lenh = 0;/* this is the current length of heap */
      lenu = ii;
      lena = A_diag_i[i+1];
      /* put those already inside original pattern, and set their level to 0 */
      for(j = A_diag_i[i] ; j < lena ; j ++)
      {
         /* get the neworder of that col */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
         }
         else if(col > ii)
         {
            /* this is an entry in U */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      }/* end of j loop for adding pattern in original matrix */

      /* put those already inside offd pattern in, and set their level to 0 */
      lena = A_offd_i[i+1];
      for( j = A_offd_i[i] ; j < lena ; j ++ )
      {
         /* the offd cols are in order */
         col = A_offd_j[j] + n;
         /* col for sure to be greater than ii */
         iL[lenu] = col;
         iLev[lenu] = 0;
         iw[col] = lenu++;
      }

      /*
       * search lower part of current row and update pattern based on level
       */
      while(lenh > 0)
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
         jxf_ILUMinHeapRemoveIIIi(iL,iLev,iw,lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k]=-1;
         jxf_swap2i(iL,iLev,ii-lenl,lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k+1];
         for(j = U_diag_i[k] ; j < ku ; j ++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if(lev > lfil)
            {
               continue;
            }
            if(icol < 0)
            {
               /* not yet in */
               if(col < ii)
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
                  jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
               }
               else if(col > ii)
               {
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jxf_min(lev, iLev[icol]);
            }
         }/* end of loop j for level update */
      }/* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* check if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jxf_TReAlloc_v2(temp_L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for(j = 0 ; j < lenl ; j ++)
         {
            temp_L_diag_j[ctrL+j] = iL[ii-j-1];
         }
         ctrL += lenl;
      }
      k = lenu - ii;
      U_diag_i[ii+1] = U_diag_i[ii] + k;
      if(k > 0)
      {
         /* check if memory is enough */
         while(ctrU + k > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            temp_U_diag_j = jxf_TReAlloc_v2(temp_U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            u_levels = jxf_TReAlloc_v2(u_levels, JXF_Int, tmp, JXF_Int, capacity_U);
         }
         jxf_TMemcpy(temp_U_diag_j+ctrU,iL+ii,JXF_Int,k);
         jxf_TMemcpy(u_levels+ctrU,iLev+ii,JXF_Int,k);
         ctrU += k;
      }

      /* reset iw */
      for(j = ii ; j < lenu ; j ++)
      {
         iw[iL[j]] = -1;
      }

   }/* end of main loop ii from nLU to n */

   /* external part matrix */
   for(ii = n ; ii < total_rows ; ii ++)
   {
      i = ii - n;
      lenl = 0;
      lenh = 0;/* this is the current length of heap */
      lenu = ii;
      lena = E_i[i+1];
      /* put those already inside original pattern, and set their level to 0 */
      for(j = E_i[i] ; j < lena ; j ++)
      {
         /* get the neworder of that col */
         col = E_j[j];
         if(col < ii)
         {
            /*
             * this is an entry in L
             * we maintain a heap structure for L part
             */
            iL[lenh] = col;
            iLev[lenh] = 0;
            iw[col] = lenh++;
            /*now miantian a heap structure*/
            jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
         }
         else if(col > ii)
         {
            /* this is an entry in U */
            iL[lenu] = col;
            iLev[lenu] = 0;
            iw[col] = lenu++;
         }
      }/* end of j loop for adding pattern in original matrix */

      /*
       * search lower part of current row and update pattern based on level
       */
      while(lenh > 0)
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
         jxf_ILUMinHeapRemoveIIIi(iL,iLev,iw,lenh);
         lenh--;
         /* copy to the end of array */
         lenl++;
         /* reset iw for that, not using anymore */
         iw[k]=-1;
         jxf_swap2i(iL,iLev,ii-lenl,lenh);
         /*
          * now the elimination on current row could start.
          * eliminate row k (new index) from current row
          */
         ku = U_diag_i[k+1];
         for(j = U_diag_i[k] ; j < ku ; j ++)
         {
            col = temp_U_diag_j[j];
            lev = u_levels[j] + ilev + 1;
            /* ignore large level */
            icol = iw[col];
            /* skill large level */
            if(lev > lfil)
            {
               continue;
            }
            if(icol < 0)
            {
               /* not yet in */
               if(col < ii)
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
                  jxf_ILUMinHeapAddIIIi(iL,iLev,iw,lenh);
               }
               else if(col > ii)
               {
                  iL[lenu] = col;
                  iLev[lenu] = lev;
                  iw[col] = lenu++;
               }
            }
            else
            {
               iLev[icol] = jxf_min(lev, iLev[icol]);
            }
         }/* end of loop j for level update */
      }/* end of while loop for iith row */

      /* now update everything, indices, levels and so */
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* check if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            temp_L_diag_j = jxf_TReAlloc_v2(temp_L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
         }
         /* now copy L data, reverse order */
         for(j = 0 ; j < lenl ; j ++)
         {
            temp_L_diag_j[ctrL+j] = iL[ii-j-1];
         }
         ctrL += lenl;
      }
      k = lenu - ii;
      U_diag_i[ii+1] = U_diag_i[ii] + k;
      if(k > 0)
      {
         /* check if memory is enough */
         while(ctrU + k > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            temp_U_diag_j = jxf_TReAlloc_v2(temp_U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            u_levels = jxf_TReAlloc_v2(u_levels, JXF_Int, tmp, JXF_Int, capacity_U);
         }
         jxf_TMemcpy(temp_U_diag_j+ctrU,iL+ii,JXF_Int,k);
         jxf_TMemcpy(u_levels+ctrU,iLev+ii,JXF_Int,k);
         ctrU += k;
      }

      /* reset iw */
      for(j = ii ; j < lenu ; j ++)
      {
         iw[iL[j]] = -1;
      }

   }/* end of main loop ii from n to total_rows */

   /*
    * 3: Finishing up and free memory
    */
   jxf_TFree(u_levels);

   *L_diag_j = temp_L_diag_j;
   *U_diag_j = temp_U_diag_j;

   return jxf_error_flag;
}

/* ILU(k) for RAS
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * perm: permutation array indicating ordering of factorization. Perm could come from a
 * CF_marker: array or a reordering routine.
 * nLU: size of computed LDU factorization.
 * Lptr, Dptr, Uptr: L, D, U factors.
 */
JXF_Int
jxf_ILUSetupILUKRAS(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Int *perm, JXF_Int nLU,
      jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr)
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
   if(lfil == 0)
   {
      return jxf_ILUSetupILU0RAS(A,perm,nLU,Lptr,Dptr,Uptr);
   }
   JXF_Int               i, ii, j, k, k1, k2, kl, ku, jpiv, col, icol;
   JXF_Int               *iw;
   MPI_Comm                comm           = jxf_ParCSRMatrixComm(A);
   JXF_Int               num_procs;

   /* data objects for A */
   jxf_CSRMatrix         *A_diag        = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix         *A_offd        = jxf_ParCSRMatrixOffd(A);
   JXF_Real              *A_diag_data   = jxf_CSRMatrixData(A_diag);
   JXF_Int               *A_diag_i      = jxf_CSRMatrixI(A_diag);
   JXF_Int               *A_diag_j      = jxf_CSRMatrixJ(A_diag);
   JXF_Real              *A_offd_data   = jxf_CSRMatrixData(A_offd);
   JXF_Int               *A_offd_i      = jxf_CSRMatrixI(A_offd);
   JXF_Int               *A_offd_j      = jxf_CSRMatrixJ(A_offd);

   /* data objects for L, D, U */
   jxf_ParCSRMatrix      *matL;
   jxf_ParCSRMatrix      *matU;
   jxf_CSRMatrix         *L_diag;
   jxf_CSRMatrix         *U_diag;
   JXF_Real              *D_data;
   JXF_Real              *L_diag_data   = NULL;
   JXF_Int               *L_diag_i;
   JXF_Int               *L_diag_j      = NULL;
   JXF_Real              *U_diag_data   = NULL;
   JXF_Int               *U_diag_i;
   JXF_Int               *U_diag_j      = NULL;

   /* size of problem and external matrix */
   JXF_Int               n              = jxf_CSRMatrixNumRows(A_diag);
   //   JXF_Int               m              = n - nLU;
   JXF_Int               ext            = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int               total_rows     = n + ext;
   JXF_Int            global_num_rows;
   JXF_Int               *col_starts;
   JXF_Real              local_nnz, total_nnz;

   /* data objects for E, external matrix */
   JXF_Int               *E_i;
   JXF_Int               *E_j;
   JXF_Real              *E_data;

   /* communication */
   jxf_ParCSRCommPkg     *comm_pkg;
   jxf_MPI_Comm_size(comm, &num_procs);
   //   jxf_ParCSRCommHandle  *comm_handle;
   //   JXF_Int               *send_buf      = NULL;

   /* reverse permutation array */
   JXF_Int               *rperm;
   /* temp array for old permutation */
   JXF_Int               *perm_old;

   /* start setup */
   /* check input and get problem size */
   n =  jxf_CSRMatrixNumRows(A_diag);
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }

   /* Init I array anyway. S's might be freed later */
   D_data   = jxf_CTAlloc(JXF_Real, total_rows);
   L_diag_i = jxf_CTAlloc(JXF_Int, (total_rows+1));
   U_diag_i = jxf_CTAlloc(JXF_Int, (total_rows+1));

   /* set Comm_Pkg if not yet built */
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /*
    * 2: Symbolic factorization
    * setup iw and rperm first
    */
   /* allocate work arrays */
   iw          = jxf_CTAlloc(JXF_Int, 5*total_rows);
   rperm       = iw + 3*total_rows;
   perm_old    = perm;
   perm        = iw + 4*total_rows;
   L_diag_i[0] = U_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    */
   for(i=0; i<n; i++)
   {
      perm[i] = perm_old[i];
   }
   for(i=n; i<total_rows; i++)
   {
      perm[i] = i;
   }
   for(i=0; i<total_rows; i++)
   {
      rperm[perm[i]] = i;
   }

   /* get external rows */
   jxf_ILUBuildRASExternalMatrix(A,rperm,&E_i,&E_j,&E_data);
   /* do symbolic factorization */
   jxf_ILUSetupILUKRASSymbolic(n, A_diag_i, A_diag_j, A_offd_i, A_offd_j, E_i, E_j, ext, lfil, perm, rperm, iw,
         nLU, L_diag_i, U_diag_i, &L_diag_j, &U_diag_j);

   /*
    * after this, we have our I,J for L, U and S ready, and L sorted
    * iw are still -1 after symbolic factorization
    * now setup helper array here
    */
   if(L_diag_i[total_rows])
   {
      L_diag_data = jxf_CTAlloc(JXF_Real, L_diag_i[total_rows]);
   }
   if(U_diag_i[total_rows])
   {
      U_diag_data = jxf_CTAlloc(JXF_Real, U_diag_i[total_rows]);
   }

   /*
    * 3: Begin real factorization
    * we already have L and U structure ready, so no extra working array needed
    */
   /* first loop for upper part */
   for( ii = 0; ii < nLU; ii++ )
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii+1];
      ku = U_diag_i[ii+1];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      /* set up working arrays */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      D_data[ii] = 0.0;
      iw[ii] = ii;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, D and U */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if(col < ii)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else if(col == ii)
         {
            D_data[ii] = A_diag_data[j];
         }
         else
         {
            U_diag_data[icol] = A_diag_data[j];
         }
      }
      /* elimination */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv+1];

         for(k = U_diag_i[jpiv] ; k < ku ; k ++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if(icol < 0)
            {
               /* not in partern */
               continue;
            }
            if(col < ii)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else if(col == ii)
            {
               /* diag part */
               D_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else
            {
               /* U part */
               U_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
         }
      }
      /* reset working array */
      ku = U_diag_i[ii+1];
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      iw[ii] = -1;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = -1;
      }

      /* diagonal part (we store the inverse) */
      if(fabs(D_data[ii]) < MAT_TOL)
      {
         D_data[ii] = 1e-06;
      }
      D_data[ii] = 1./ D_data[ii];

   }/* end of loop for upper part */

   /* first loop for upper part */
   for( ii = nLU; ii < n; ii++ )
   {
      // get row i
      i = perm[ii];
      kl = L_diag_i[ii+1];
      ku = U_diag_i[ii+1];
      /* set up working arrays */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      D_data[ii] = 0.0;
      iw[ii] = ii;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = j;
      }
      /* copy data from A into L, D and U */
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = rperm[A_diag_j[j]];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if(col < ii)
         {
            L_diag_data[icol] = A_diag_data[j];
         }
         else if(col == ii)
         {
            D_data[ii] = A_diag_data[j];
         }
         else
         {
            U_diag_data[icol] = A_diag_data[j];
         }
      }
      /* copy data from A_offd into L, D and U */
      k1 = A_offd_i[i];
      k2 = A_offd_i[i+1];
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = A_offd_j[j] + n;
         icol = iw[col];
         U_diag_data[icol] = A_offd_data[j];
      }
      /* elimination */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv+1];

         for(k = U_diag_i[jpiv] ; k < ku ; k ++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if(icol < 0)
            {
               /* not in partern */
               continue;
            }
            if(col < ii)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else if(col == ii)
            {
               /* diag part */
               D_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else
            {
               /* U part */
               U_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
         }
      }
      /* reset working array */
      ku = U_diag_i[ii+1];
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      iw[ii] = -1;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = -1;
      }

      /* diagonal part (we store the inverse) */
      if(fabs(D_data[ii]) < MAT_TOL)
      {
         D_data[ii] = 1e-06;
      }
      D_data[ii] = 1./ D_data[ii];

   }/* end of loop for lower part */

   /* last loop through external */
   for( ii = n; ii < total_rows; ii++ )
   {
      // get row i
      i = ii - n;
      kl = L_diag_i[ii+1];
      ku = U_diag_i[ii+1];
      k1 = E_i[i];
      k2 = E_i[i+1];
      /* set up working arrays */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = j;
      }
      D_data[ii] = 0.0;
      iw[ii] = ii;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = j;
      }
      /* copy data from E into L, D and U */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* compute everything in new index */
         col = E_j[j];
         icol = iw[col];
         /* A for sure to be inside the pattern */
         if(col < ii)
         {
            L_diag_data[icol] = E_data[j];
         }
         else if(col == ii)
         {
            D_data[ii] = E_data[j];
         }
         else
         {
            U_diag_data[icol] = E_data[j];
         }
      }
      /* elimination */
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         jpiv = L_diag_j[j];
         L_diag_data[j] *= D_data[jpiv];
         ku = U_diag_i[jpiv+1];

         for(k = U_diag_i[jpiv] ; k < ku ; k ++)
         {
            col = U_diag_j[k];
            icol = iw[col];
            if(icol < 0)
            {
               /* not in partern */
               continue;
            }
            if(col < ii)
            {
               /* L part */
               L_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else if(col == ii)
            {
               /* diag part */
               D_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
            else
            {
               /* U part */
               U_diag_data[icol] -= L_diag_data[j]*U_diag_data[k];
            }
         }
      }
      /* reset working array */
      ku = U_diag_i[ii+1];
      for(j = L_diag_i[ii] ; j < kl ; j ++)
      {
         col = L_diag_j[j];
         iw[col] = -1;
      }
      iw[ii] = -1;
      for(j = U_diag_i[ii] ; j < ku ; j ++)
      {
         col = U_diag_j[j];
         iw[col] = -1;
      }

      /* diagonal part (we store the inverse) */
      if(fabs(D_data[ii]) < MAT_TOL)
      {
         D_data[ii] = 1e-06;
      }
      D_data[ii] = 1./ D_data[ii];

   }/* end of loop for external loop */

   /*
    * 4: Finishing up and free
    */
   JXF_Int big_total_rows = (JXF_Int)total_rows;
   jxf_MPI_Allreduce( &big_total_rows, &global_num_rows, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
   {
      JXF_Int global_start;
      col_starts = jxf_CTAlloc(JXF_Int,2);
      jxf_MPI_Scan( &big_total_rows, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
      col_starts[0] = global_start - total_rows;
      col_starts[1] = global_start;
   }
#else
   col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

   jxf_MPI_Allgather(&big_total_rows,1,JXF_MPI_INT,&col_starts[1],
         1,JXF_MPI_INT,comm);

   for (i=2; i < num_procs+1; i++)
      col_starts[i] += col_starts[i-1];

#endif
   /* Assemble LDU matrices */
   matL = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0 /* num_cols_offd */,
         L_diag_i[total_rows],
         0 /* num_nonzeros_offd */);

   /* Have A own coarse_partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,1);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[total_rows]>0)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jxf_TFree(L_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (L_diag_i[total_rows]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0,
         U_diag_i[total_rows],
         0 );

   /* Have A own coarse_partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);

   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[n]>0)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we allocated some initial length, so free them */
      jxf_TFree(U_diag_j);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (U_diag_i[total_rows]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free */
   jxf_TFree(iw);

   /* free external data */
   if(E_i)
   {
      jxf_TFree(E_i);
   }
   if(E_j)
   {
      jxf_TFree(E_j);
      jxf_TFree(E_data);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;

   return jxf_error_flag;
}

/* ILUT for RAS
 * A: input matrix
 * lfil: level of fill-in, the k in ILU(k)
 * tol: droptol array in ILUT
 *    tol[0]: matrix B
 *    tol[1]: matrix E and F
 *    tol[2]: matrix S
 * perm: permutation array indicating ordering of factorization. Perm could come from a
 * CF_marker: array or a reordering routine.
 * nLU: size of computed LDU factorization. If nLU < n, Schur compelemnt will be formed
 * Lptr, Dptr, Uptr: L, D, U factors.
 * Sptr: Schur complement
 *
 * Keep the largest lfil entries that is greater than some tol relative
 *    to the input tol and the norm of that row in both L and U
 */
JXF_Int
jxf_ILUSetupILUTRAS(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Real *tol,
      JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix **Lptr,
      JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr)
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
   JXF_Real               local_nnz, total_nnz;
   JXF_Int                i, ii, j, k1, k2, k12, k22, kl, ku, col, icol, lenl, lenu, lenhu, lenhlr, lenhll, jpos, jrow;
   JXF_Real               inorm, itolb, itolef, dpiv, lxu;
   JXF_Int                *iw,*iL;
   JXF_Real               *w;

   /* memory management */
   JXF_Int                ctrL;
   JXF_Int                ctrU;
   JXF_Int                initial_alloc = 0;
   JXF_Int                capacity_L;
   JXF_Int                capacity_U;
   JXF_Int                nnz_A;

   /* communication stuffs for S */
   MPI_Comm                 comm          = jxf_ParCSRMatrixComm(A);
   JXF_Int                num_procs;
   jxf_ParCSRCommPkg      *comm_pkg;
   //   jxf_ParCSRCommHandle   *comm_handle;
   JXF_Int                *col_starts;
   //   JXF_Int                num_sends;
   //   JXF_Int                begin, end;

   /* data objects for A */
   jxf_CSRMatrix          *A_diag       = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix          *A_offd       = jxf_ParCSRMatrixOffd(A);
   JXF_Real               *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int                *A_diag_i     = jxf_CSRMatrixI(A_diag);
   JXF_Int                *A_diag_j     = jxf_CSRMatrixJ(A_diag);
   JXF_Int                *A_offd_i     = jxf_CSRMatrixI(A_offd);
   JXF_Int                *A_offd_j     = jxf_CSRMatrixJ(A_offd);
   JXF_Real               *A_offd_data  = jxf_CSRMatrixData(A_offd);

   /* data objects for L, D, U */
   jxf_ParCSRMatrix       *matL;
   jxf_ParCSRMatrix       *matU;
   jxf_CSRMatrix          *L_diag;
   jxf_CSRMatrix          *U_diag;
   JXF_Real               *D_data;
   JXF_Real               *L_diag_data  = NULL;
   JXF_Int                *L_diag_i;
   JXF_Int                *L_diag_j     = NULL;
   JXF_Real               *U_diag_data  = NULL;
   JXF_Int                *U_diag_i;
   JXF_Int                *U_diag_j     = NULL;

   /* size of problem and external matrix */
   JXF_Int                n             = jxf_CSRMatrixNumRows(A_diag);
   //   JXF_Int                m             = n - nLU;
   JXF_Int                ext           = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int                total_rows    = n + ext;
   JXF_Int              global_num_rows;

   /* data objects for E, external matrix */
   JXF_Int                *E_i;
   JXF_Int                *E_j;
   JXF_Real               *E_data;

   /* reverse permutation */
   JXF_Int                *rperm;
   /* old permutation */
   JXF_Int                *perm_old;

   /* start setup
    * check input first
    */
   n = jxf_CSRMatrixNumRows(A_diag);
   if(nLU < 0 || nLU > n)
   {
      jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: nLU out of range.\n");
   }

   /* start set up
    * setup communication stuffs first
    */
   jxf_MPI_Comm_size(comm, &num_procs);
   comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   /* create if not yet built */
   if(!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* setup initial memory */
   nnz_A = A_diag_i[nLU];
   if(n > 0)
   {
      initial_alloc = nLU + ceil(nnz_A / 2.0);
   }
   capacity_L = initial_alloc;
   capacity_U = initial_alloc;

   D_data = jxf_CTAlloc(JXF_Real, total_rows);
   L_diag_i = jxf_CTAlloc(JXF_Int, (total_rows+1));
   U_diag_i = jxf_CTAlloc(JXF_Int, (total_rows+1));

   L_diag_j = jxf_CTAlloc(JXF_Int, capacity_L);
   U_diag_j = jxf_CTAlloc(JXF_Int, capacity_U);
   L_diag_data = jxf_CTAlloc(JXF_Real, capacity_L);
   U_diag_data = jxf_CTAlloc(JXF_Real, capacity_U);

   ctrL = ctrU = 0;

   /* setting up working array */
   iw = jxf_CTAlloc(JXF_Int,4*total_rows);
   iL = iw + total_rows;
   w = jxf_CTAlloc(JXF_Real,total_rows);
   for(i = 0 ; i < total_rows ; i ++)
   {
      iw[i] = -1;
   }
   L_diag_i[0] = U_diag_i[0] = 0;
   /* get reverse permutation (rperm).
    * rperm holds the reordered indexes.
    * rperm[old] -> new
    * perm[new]  -> old
    */
   rperm = iw + 2*total_rows;
   perm_old = perm;
   perm = iw + 3*total_rows;
   for(i = 0 ; i < n ; i ++)
   {
      perm[i] = perm_old[i];
   }
   for(i = n ; i < total_rows ; i ++)
   {
      perm[i] = i;
   }
   for(i = 0 ; i < total_rows ; i ++)
   {
      rperm[perm[i]] = i;
   }
   /* get external matrix */
   jxf_ILUBuildRASExternalMatrix(A,rperm,&E_i,&E_j,&E_data);

   /*
    * 2: Main loop of elimination
    * maintain two heaps
    * |----->*********<-----|-----*********|
    * |col heap***value heap|value in U****|
    */

   /* main outer loop for upper part */
   for(ii = 0 ; ii < nLU ; ii ++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      kl = ii-1;
      /* reset row norm of ith row */
      inorm = .0;
      for(j = k1 ; j < k2 ; j ++)
      {
         inorm += fabs(A_diag_data[j]);
      }
      if(inorm == .0)
      {
         jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: ILUT with zero row.\n");
      }
      inorm /= (JXF_Real)(k2-k1);
      /* set the scaled tol for that row */
      itolb = tol[0] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      w[ii] = 0.0;
      iw[ii] = ii;
      /* copy in data from A */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
         }
         else if(col == ii)
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
      while(lenhll > 0)
      {

         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jxf_ILUMinHeapRemoveIRIi(iL,w,iw,lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jxf_swap2(iL,w,lenhll,kl-lenhlr);
         lenhlr++;
         jxf_ILUMaxrHeapAddRabsI(w+kl,iL+kl,lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow+1];
         for(j = U_diag_i[jrow] ; j < ku ; j ++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = - dpiv*U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if( icol == -1 && ( (col < nLU && fabs(lxu) < itolb) || (col >= nLU && fabs(lxu) < itolef) ) )
            {
               continue;
            }
            if(icol == -1)
            {
               if(col < ii)
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
                  jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
               }
               else if(col == ii)
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
      }/* while loop for the elimination of current row */

      if(fabs(w[ii]) < MAT_TOL)
      {
         w[ii]=1e-06;
      }
      D_data[ii] = 1./w[ii];
      iw[ii] = -1;

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* test if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         ctrL += lenl;
         /* copy large data in */
         for(j = L_diag_i[ii] ; j < ctrL ; j ++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jxf_ILUMaxrHeapRemoveRabsI(w+kl,iL+kl,lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only U part
       */
      ku = lenu+ii;
      for(j = ii + 1 ; j <= ku ; j ++)
      {
         iw[iL[j]] = -1;
      }

      if(lenu < lfil)
      {
         /* we simply keep all of the data, no need to sort */
         lenhu = lenu;
      }
      else
      {
         /* need to sort the first small(hopefully) part of it */
         lenhu = lfil;
         /* quick split, only sort the first small part of the array */
         jxf_ILUMaxQSplitRabsI(w,iL,ii+1,ii+lenhu,ii+lenu);
      }

      U_diag_i[ii+1] = U_diag_i[ii] + lenhu;
      if(lenhu > 0)
      {
         /* test if memory is enough */
         while(ctrU + lenhu > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
         }
         ctrU += lenhu;
         /* copy large data in */
         for(j = U_diag_i[ii] ; j < ctrU ; j ++)
         {
            jpos = ii+1+j-U_diag_i[ii];
            U_diag_j[j] = iL[jpos];
            U_diag_data[j] = w[jpos];
         }
      }
   }/* end of ii loop from 0 to nLU-1 */


   /* second outer loop for lower part */
   for(ii = nLU ; ii < n ; ii ++)
   {
      /* get real row with perm */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      k12 = A_offd_i[i];
      k22 = A_offd_i[i+1];
      kl = ii-1;
      /* reset row norm of ith row */
      inorm = .0;
      for(j = k1 ; j < k2 ; j ++)
      {
         inorm += fabs(A_diag_data[j]);
      }
      for(j = k12 ; j < k22 ; j ++)
      {
         inorm += fabs(A_offd_data[j]);
      }
      if(inorm == .0)
      {
         jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: ILUT with zero row.\n");
      }
      inorm /= (JXF_Real)(k2+k22-k1-k12);
      /* set the scaled tol for that row */
      itolb = tol[0] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      w[ii] = 0.0;
      iw[ii] = ii;
      /* copy in data from A_diag */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* get now col number */
         col = rperm[A_diag_j[j]];
         if(col < ii)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = A_diag_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
         }
         else if(col == ii)
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
      /* copy in data from A_offd */
      for(j = k12 ; j < k22 ; j ++)
      {
         /* get now col number */
         col = A_offd_j[j] + n;
         /* all should greater than ii in lower part */
         lenu++;
         jpos = lenu + ii;
         iL[jpos] = col;
         w[jpos] = A_offd_data[j];
         iw[col] = jpos;
      }

      /*
       * main elimination
       * need to maintain 2 heaps for L, one heap for col and one heaps for value
       * maintian an array for U, and do qsplit with quick sort after that
       * while the heap of col is greater than zero
       */
      while(lenhll > 0)
      {

         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jxf_ILUMinHeapRemoveIRIi(iL,w,iw,lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jxf_swap2(iL,w,lenhll,kl-lenhlr);
         lenhlr++;
         jxf_ILUMaxrHeapAddRabsI(w+kl,iL+kl,lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow+1];
         for(j = U_diag_i[jrow] ; j < ku ; j ++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = - dpiv*U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if( icol == -1 && ( (col < nLU && fabs(lxu) < itolb) || (col >= nLU && fabs(lxu) < itolef) ) )
            {
               continue;
            }
            if(icol == -1)
            {
               if(col < ii)
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
                  jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
               }
               else if(col == ii)
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
      }/* while loop for the elimination of current row */

      if(fabs(w[ii]) < MAT_TOL)
      {
         w[ii]=1e-06;
      }
      D_data[ii] = 1./w[ii];
      iw[ii] = -1;

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* test if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         ctrL += lenl;
         /* copy large data in */
         for(j = L_diag_i[ii] ; j < ctrL ; j ++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jxf_ILUMaxrHeapRemoveRabsI(w+kl,iL+kl,lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only U part
       */
      ku = lenu+ii;
      for(j = ii + 1 ; j <= ku ; j ++)
      {
         iw[iL[j]] = -1;
      }

      if(lenu < lfil)
      {
         /* we simply keep all of the data, no need to sort */
         lenhu = lenu;
      }
      else
      {
         /* need to sort the first small(hopefully) part of it */
         lenhu = lfil;
         /* quick split, only sort the first small part of the array */
         jxf_ILUMaxQSplitRabsI(w,iL,ii+1,ii+lenhu,ii+lenu);
      }

      U_diag_i[ii+1] = U_diag_i[ii] + lenhu;
      if(lenhu > 0)
      {
         /* test if memory is enough */
         while(ctrU + lenhu > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
         }
         ctrU += lenhu;
         /* copy large data in */
         for(j = U_diag_i[ii] ; j < ctrU ; j ++)
         {
            jpos = ii+1+j-U_diag_i[ii];
            U_diag_j[j] = iL[jpos];
            U_diag_data[j] = w[jpos];
         }
      }
   }/* end of ii loop from nLU to n */


   /* main outer loop for upper part */
   for(ii = n ; ii < total_rows ; ii ++)
   {
      /* get real row with perm */
      i = ii-n;
      k1 = E_i[i];
      k2 = E_i[i+1];
      kl = ii-1;
      /* reset row norm of ith row */
      inorm = .0;
      for(j = k1 ; j < k2 ; j ++)
      {
         inorm += fabs(E_data[j]);
      }
      if(inorm == .0)
      {
         jxf_error_w_msg(JXF_ERROR_ARG,"WARNING: ILUT with zero row.\n");
      }
      inorm /= (JXF_Real)(k2-k1);
      /* set the scaled tol for that row */
      itolb = tol[0] * inorm;
      itolef = tol[1] * inorm;

      /* reset displacement */
      lenhll = lenhlr = lenu = 0;
      w[ii] = 0.0;
      iw[ii] = ii;
      /* copy in data from A */
      for(j = k1 ; j < k2 ; j ++)
      {
         /* get now col number */
         col = rperm[E_j[j]];
         if(col < ii)
         {
            /* L part of it */
            iL[lenhll] = col;
            w[lenhll] = E_data[j];
            iw[col] = lenhll++;
            /* add to heap, by col number */
            jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
         }
         else if(col == ii)
         {
            w[ii] = E_data[j];
         }
         else
         {
            lenu++;
            jpos = lenu + ii;
            iL[jpos] = col;
            w[jpos] = E_data[j];
            iw[col] = jpos;
         }
      }

      /*
       * main elimination
       * need to maintain 2 heaps for L, one heap for col and one heaps for value
       * maintian an array for U, and do qsplit with quick sort after that
       * while the heap of col is greater than zero
       */
      while(lenhll > 0)
      {

         /* get the next row from top of the heap */
         jrow = iL[0];
         dpiv = w[0] * D_data[jrow];
         w[0] = dpiv;
         /* now remove it from the top of the heap */
         jxf_ILUMinHeapRemoveIRIi(iL,w,iw,lenhll);
         lenhll--;
         /*
          * reset the drop part to -1
          * we don't need this iw anymore
          */
         iw[jrow] = -1;
         /* need to keep this one, move to the end of the heap */
         /* no longer need to maintain iw */
         jxf_swap2(iL,w,lenhll,kl-lenhlr);
         lenhlr++;
         jxf_ILUMaxrHeapAddRabsI(w+kl,iL+kl,lenhlr);
         /* loop for elimination */
         ku = U_diag_i[jrow+1];
         for(j = U_diag_i[jrow] ; j < ku ; j ++)
         {
            col = U_diag_j[j];
            icol = iw[col];
            lxu = - dpiv*U_diag_data[j];
            /* we don't want to fill small number to empty place */
            if( icol == -1 && ( (col < nLU && fabs(lxu) < itolb) || (col >= nLU && fabs(lxu) < itolef) ) )
            {
               continue;
            }
            if(icol == -1)
            {
               if(col < ii)
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
                  jxf_ILUMinHeapAddIRIi(iL,w,iw,lenhll);
               }
               else if(col == ii)
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
      }/* while loop for the elimination of current row */

      if(fabs(w[ii]) < MAT_TOL)
      {
         w[ii]=1e-06;
      }
      D_data[ii] = 1./w[ii];
      iw[ii] = -1;

      /*
       * now pick up the largest lfil from L
       * L part is guarantee to be larger than itol
       */

      lenl = lenhlr < lfil ? lenhlr : lfil;
      L_diag_i[ii+1] = L_diag_i[ii] + lenl;
      if(lenl > 0)
      {
         /* test if memory is enough */
         while(ctrL + lenl > capacity_L)
         {
            JXF_Int tmp = capacity_L;
            capacity_L = capacity_L * EXPAND_FACT + 1;
            L_diag_j = jxf_TReAlloc_v2(L_diag_j, JXF_Int, tmp, JXF_Int, capacity_L);
            L_diag_data = jxf_TReAlloc_v2(L_diag_data, JXF_Real, tmp, JXF_Real, capacity_L);
         }
         ctrL += lenl;
         /* copy large data in */
         for(j = L_diag_i[ii] ; j < ctrL ; j ++)
         {
            L_diag_j[j] = iL[kl];
            L_diag_data[j] = w[kl];
            jxf_ILUMaxrHeapRemoveRabsI(w+kl,iL+kl,lenhlr);
            lenhlr--;
         }
      }
      /*
       * now reset working array
       * L part already reset when move out of heap, only U part
       */
      ku = lenu+ii;
      for(j = ii + 1 ; j <= ku ; j ++)
      {
         iw[iL[j]] = -1;
      }

      if(lenu < lfil)
      {
         /* we simply keep all of the data, no need to sort */
         lenhu = lenu;
      }
      else
      {
         /* need to sort the first small(hopefully) part of it */
         lenhu = lfil;
         /* quick split, only sort the first small part of the array */
         jxf_ILUMaxQSplitRabsI(w,iL,ii+1,ii+lenhu,ii+lenu);
      }

      U_diag_i[ii+1] = U_diag_i[ii] + lenhu;
      if(lenhu > 0)
      {
         /* test if memory is enough */
         while(ctrU + lenhu > capacity_U)
         {
            JXF_Int tmp = capacity_U;
            capacity_U = capacity_U * EXPAND_FACT + 1;
            U_diag_j = jxf_TReAlloc_v2(U_diag_j, JXF_Int, tmp, JXF_Int, capacity_U);
            U_diag_data = jxf_TReAlloc_v2(U_diag_data, JXF_Real, tmp, JXF_Real, capacity_U);
         }
         ctrU += lenhu;
         /* copy large data in */
         for(j = U_diag_i[ii] ; j < ctrU ; j ++)
         {
            jpos = ii+1+j-U_diag_i[ii];
            U_diag_j[j] = iL[jpos];
            U_diag_data[j] = w[jpos];
         }
      }
   }/* end of ii loop from nLU to total_rows */

   /*
    * 3: Finishing up and free
    */
   JXF_Int big_total_rows = (JXF_Int)total_rows;
   jxf_MPI_Allreduce( &big_total_rows, &global_num_rows, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* need to get new column start */
#ifdef JXF_NO_GLOBAL_PARTITION
   {
      JXF_Int global_start;
      col_starts = jxf_CTAlloc(JXF_Int,2);
      jxf_MPI_Scan( &big_total_rows, &global_start, 1, JXF_MPI_INT, MPI_SUM, comm);
      col_starts[0] = global_start - total_rows;
      col_starts[1] = global_start;
   }
#else
   col_starts = jxf_CTAlloc(JXF_Int,num_procs+1);

   jxf_MPI_Allgather(&big_total_rows,1,JXF_MPI_INT,&col_starts[1],
         1,JXF_MPI_INT,comm);

   for (i=2; i < num_procs+1; i++)
      col_starts[i] += col_starts[i-1];
#endif

   /* create parcsr matrix */
   matL = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0,
         L_diag_i[total_rows],
         0 );

   /* Have A own coarse_partitioning instead of L */
   jxf_ParCSRMatrixSetColStartsOwner(matL,1);
   jxf_ParCSRMatrixSetRowStartsOwner(matL,0);
   L_diag = jxf_ParCSRMatrixDiag(matL);
   jxf_CSRMatrixI(L_diag) = L_diag_i;
   if (L_diag_i[total_rows] > 0)
   {
      jxf_CSRMatrixData(L_diag) = L_diag_data;
      jxf_CSRMatrixJ(L_diag) = L_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jxf_TFree(L_diag_j);
      jxf_TFree(L_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (L_diag_i[total_rows]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matL) = total_nnz;

   matU = jxf_ParCSRMatrixCreate( comm,
         global_num_rows,
         global_num_rows,
         col_starts,
         col_starts,
         0,
         U_diag_i[total_rows],
         0 );

   /* Have A own coarse_partitioning instead of U */
   jxf_ParCSRMatrixSetColStartsOwner(matU,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matU,0);

   U_diag = jxf_ParCSRMatrixDiag(matU);
   jxf_CSRMatrixI(U_diag) = U_diag_i;
   if (U_diag_i[total_rows] > 0)
   {
      jxf_CSRMatrixData(U_diag) = U_diag_data;
      jxf_CSRMatrixJ(U_diag) = U_diag_j;
   }
   else
   {
      /* we initialized some anyway, so remove if unused */
      jxf_TFree(U_diag_j);
      jxf_TFree(U_diag_data);
   }
   /* store (global) total number of nonzeros */
   local_nnz = (JXF_Real) (U_diag_i[total_rows]);
   jxf_MPI_Allreduce(&local_nnz, &total_nnz, 1, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_ParCSRMatrixDNumNonzeros(matU) = total_nnz;

   /* free working array */
   jxf_TFree(iw);
   jxf_TFree(w);

   /* free external data */
   if(E_i)
   {
      jxf_TFree(E_i);
   }
   if(E_j)
   {
      jxf_TFree(E_j);
      jxf_TFree(E_data);
   }

   /* set matrix pointers */
   *Lptr = matL;
   *Dptr = D_data;
   *Uptr = matU;

   return jxf_error_flag;
}
