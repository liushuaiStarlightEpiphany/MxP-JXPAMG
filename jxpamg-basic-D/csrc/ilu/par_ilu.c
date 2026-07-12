//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu.c
 *
 */

#include "jx_ilu.h"
#include "jx_krylov.h"

JX_Int
JX_ILUCreate(JX_Solver *solver)
{
   if (!solver)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   *solver = ((JX_Solver)jx_ILUCreate());
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * JX_ILUDestroy
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUDestroy(JX_Solver solver)
{
   return (jx_ILUDestroy((void *)solver));
}

/*--------------------------------------------------------------------------
 * JX_ILUSetup
 *--------------------------------------------------------------------------*/

// JX_Int
// JX_ILUSetup(JX_Solver solver, JX_ParCSRMatrix A)
// {
//    return (jx_ILUSetup((void *)solver, (jx_ParCSRMatrix *)A));
// }

// JX_Int
// JX_ILUSetup(JX_Solver solver,
//             JX_ParCSRMatrix A,
//             JX_ParVector b,
//             JX_ParVector x)
// {
//    return (jx_ILUSetup((void *)solver,
//                        (jx_ParCSRMatrix *)A,
//                        (jx_ParVector *)b,
//                        (jx_ParVector *)x));
// }

JX_Int
JX_ILUSetup(JX_Solver solver,
            jx_hpCSRMatrix *A,
            JX_ParVector *b,
            JX_ParVector *x)
{
   jx_ParCSRMatrix *par_A = jx_hpCSRMatrixPar(A);

   return (jx_ILUSetup((void *)solver,
                       (jx_ParCSRMatrix *)par_A,
                       (jx_ParVector *)b,
                       (jx_ParVector *)x));
}

/*--------------------------------------------------------------------------
 * JX_ILUSolve
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSolve(JX_Solver solver,
            jx_hpCSRMatrix *A,
            JX_ParVector *b,
            JX_ParVector *x)
{
   jx_ParCSRMatrix *par_A = jx_hpCSRMatrixPar(A);
   return (jx_ILUSolve((void *)solver,
                       (jx_ParCSRMatrix *)par_A,
                       (jx_ParVector *)b,
                       (jx_ParVector *)x));
}

/*--------------------------------------------------------------------------
 * JX_ILUSetPrintLevel
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetPrintLevel(JX_Solver solver, JX_Int print_level)
{
   return jx_ILUSetPrintLevel(solver, print_level);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLogging
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLogging(JX_Solver solver, JX_Int logging)
{
   return jx_ILUSetLogging(solver, logging);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetMaxIter
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetMaxIter(JX_Solver solver, JX_Int max_iter)
{
   return jx_ILUSetMaxIter(solver, max_iter);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetIterativeSetupType
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetIterativeSetupType(JX_Solver solver, JX_Int iter_setup_type)
{
   return jx_ILUSetIterativeSetupType(solver, iter_setup_type);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetIterativeSetupOption
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetIterativeSetupOption(JX_Solver solver, JX_Int iter_setup_option)
{
   return jx_ILUSetIterativeSetupOption(solver, iter_setup_option);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetIterativeSetupMaxIter
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetIterativeSetupMaxIter(JX_Solver solver, JX_Int iter_setup_max_iter)
{
   return jx_ILUSetIterativeSetupMaxIter(solver, iter_setup_max_iter);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetIterativeSetupTolerance
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetIterativeSetupTolerance(JX_Solver solver, JX_Real iter_setup_tolerance)
{
   return jx_ILUSetIterativeSetupTolerance(solver, iter_setup_tolerance);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetTriSolve
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetTriSolve(JX_Solver solver, JX_Int tri_solve)
{
   return jx_ILUSetTriSolve(solver, tri_solve);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLowerJacobiIters
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLowerJacobiIters(JX_Solver solver, JX_Int lower_jacobi_iters)
{
   return jx_ILUSetLowerJacobiIters(solver, lower_jacobi_iters);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetUpperJacobiIters
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetUpperJacobiIters(JX_Solver solver, JX_Int upper_jacobi_iters)
{
   return jx_ILUSetUpperJacobiIters(solver, upper_jacobi_iters);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetTol
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetTol(JX_Solver solver, JX_Real tol)
{
   return jx_ILUSetTol(solver, tol);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetDropThreshold
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetDropThreshold(JX_Solver solver, JX_Real threshold)
{
   return jx_ILUSetDropThreshold(solver, threshold);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetDropThresholdArray
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetDropThresholdArray(JX_Solver solver, JX_Real *threshold)
{
   return jx_ILUSetDropThresholdArray(solver, threshold);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetNSHDropThreshold
 *--------------------------------------------------------------------------*/

// JX_Int
// JX_ILUSetNSHDropThreshold(JX_Solver solver, JX_Real threshold)
// {
//    return jx_ILUSetSchurNSHDropThreshold(solver, threshold);
// }

/*--------------------------------------------------------------------------
 * JX_ILUSetNSHDropThresholdArray
 *--------------------------------------------------------------------------*/

// JX_Int
// JX_ILUSetNSHDropThresholdArray(JX_Solver solver, JX_Real *threshold)
// {
//    return jx_ILUSetSchurNSHDropThresholdArray(solver, threshold);
// }

/*--------------------------------------------------------------------------
 * JX_ILUSetSchurMaxIter
 *--------------------------------------------------------------------------*/

// JX_Int
// JX_ILUSetSchurMaxIter(JX_Solver solver, JX_Int ss_max_iter)
// {
//    return jx_jx_ILUSetSchurSolverMaxIter(solver, ss_max_iter);
// }

/*--------------------------------------------------------------------------
 * JX_ILUSetMaxNnzPerRow
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetMaxNnzPerRow(JX_Solver solver, JX_Int nzmax)
{
   return jx_ILUSetMaxNnzPerRow(solver, nzmax);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLevelOfFill
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLevelOfFill(JX_Solver solver, JX_Int lfil)
{
   return jx_ILUSetLevelOfFill(solver, lfil);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetType
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetType(JX_Solver solver, JX_Int ilu_type)
{
   return jx_ILUSetType(solver, ilu_type);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLocalReordering
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLocalReordering(JX_Solver solver, JX_Int ordering_type)
{
   return jx_ILUSetLocalReordering(solver, ordering_type);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLocalResidualCorrectionType
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLocalResidualCorrectionType(JX_Solver solver, JX_Int lrc_type)
{
   return jx_ILUSetLocalResidualCorrectionType(solver, lrc_type);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLocalResidualCorrectionMaxIter
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLocalResidualCorrectionMaxIter(JX_Solver solver, JX_Int lrc_max_iter)
{
   return jx_ILUSetLocalResidualCorrectionMaxIter(solver, lrc_max_iter);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLocalResidualCorrectionThreshold
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLocalResidualCorrectionThreshold(JX_Solver solver, JX_Real lrc_threshold)
{
   return jx_ILUSetLocalResidualCorrectionThreshold(solver, lrc_threshold);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetLocalResidualCorrectionOmega
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUSetLocalResidualCorrectionOmega(JX_Solver solver, JX_Real lrc_omega)
{
   return jx_ILUSetLocalResidualCorrectionOmega(solver, lrc_omega);
}

JX_Int
JX_ILUSetLocalResidualCorrectionSelectType(JX_Solver solver, JX_Int lrc_select_type)
{
   return jx_ILUSetLocalResidualCorrectionSelectType(solver, lrc_select_type);
}

JX_Int
JX_ILUSetLocalResidualCorrectionEnergyTarget(JX_Solver solver, JX_Real lrc_energy_target)
{
   return jx_ILUSetLocalResidualCorrectionEnergyTarget(solver, lrc_energy_target);
}

/*--------------------------------------------------------------------------
 * JX_ILUGetNumIterations
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUGetNumIterations(JX_Solver solver, JX_Int *num_iterations)
{
   return jx_ILUGetNumIterations(solver, num_iterations);
}

/*--------------------------------------------------------------------------
 * JX_ILUGetFinalRelativeResidualNorm
 *--------------------------------------------------------------------------*/

JX_Int
JX_ILUGetFinalRelativeResidualNorm(JX_Solver solver, JX_Real *res_norm)
{
   return jx_ILUGetFinalRelativeResidualNorm(solver, res_norm);
}

/*--------------------------------------------------------------------------
 * JX_ILUSetIRIters （张娜新增）
 *--------------------------------------------------------------------------*/

// JX_Int
// JX_ILUSetIRIters(JX_Solver solver, JX_Int IR_iters)
// {
//    return jx_ILUSetIRIters(solver, IR_iters);
// }

/* ======================================
 * ====== 具体实现  ======
   ======================================*/
/* Create */
void *
jx_ILUCreate(void)
{
   jx_ParILUData *ilu_data;

   ilu_data = jx_CTAlloc(jx_ParILUData, 1);

   /* general data */
   jx_ParILUDataGlobalSolver(ilu_data) = 0;
   jx_ParILUDataMatA(ilu_data) = NULL;
   jx_ParILUDataMatL(ilu_data) = NULL;
   jx_ParILUDataMatD(ilu_data) = NULL;
   jx_ParILUDataMatU(ilu_data) = NULL;
   jx_ParILUDataMatS(ilu_data) = NULL;
   jx_ParILUDataSchurSolver(ilu_data) = NULL;
   jx_ParILUDataSchurPrecond(ilu_data) = NULL;
   jx_ParILUDataRhs(ilu_data) = NULL;
   jx_ParILUDataX(ilu_data) = NULL;

   /* TODO (VPM): Transform this into a stack array */
   jx_ParILUDataDroptol(ilu_data) = jx_TAlloc(JX_Real, 3);
   jx_ParILUDataDroptol(ilu_data)[0] = 1.0e-02; /* droptol for B */
   jx_ParILUDataDroptol(ilu_data)[1] = 1.0e-02; /* droptol for E and F */
   jx_ParILUDataDroptol(ilu_data)[2] = 1.0e-02; /* droptol for S */
   jx_ParILUDataLfil(ilu_data) = 0;
   jx_ParILUDataMaxRowNnz(ilu_data) = 1000;
   jx_ParILUDataCFMarkerArray(ilu_data) = NULL;
   jx_ParILUDataPerm(ilu_data) = NULL;
   jx_ParILUDataQPerm(ilu_data) = NULL;
   jx_ParILUDataTolDDPQ(ilu_data) = 1.0e-01;
   jx_ParILUDataF(ilu_data) = NULL;
   jx_ParILUDataU(ilu_data) = NULL;
   jx_ParILUDataFTemp(ilu_data) = NULL;
   jx_ParILUDataUTemp(ilu_data) = NULL;
   jx_ParILUDataXTemp(ilu_data) = NULL;
   jx_ParILUDataYTemp(ilu_data) = NULL;
   jx_ParILUDataZTemp(ilu_data) = NULL;
   jx_ParILUDataUExt(ilu_data) = NULL;
   jx_ParILUDataFExt(ilu_data) = NULL;
   jx_ParILUDataResidual(ilu_data) = NULL;
   jx_ParILUDataRelResNorms(ilu_data) = NULL;
   jx_ParILUDataNumIterations(ilu_data) = 0;
   jx_ParILUDataMaxIter(ilu_data) = 20;
   jx_ParILUDataTriSolve(ilu_data) = 1;
   jx_ParILUDataLowerJacobiIters(ilu_data) = 5;
   jx_ParILUDataUpperJacobiIters(ilu_data) = 5;
   jx_ParILUDataLocalResidualCorrectionType(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionMaxIter(ilu_data) = 1;
   jx_ParILUDataLocalResidualCorrectionThreshold(ilu_data) = 0.5;
   jx_ParILUDataLocalResidualCorrectionOmega(ilu_data) = 0.3;
   jx_ParILUDataLocalResidualCorrectionSelectType(ilu_data) = 0;
   jx_ParILUDataLocalResidualCorrectionEnergyTarget(ilu_data) = 0.5;
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
   jx_ParILUDataTol(ilu_data) = 1.0e-7;
   jx_ParILUDataLogging(ilu_data) = 0;
   jx_ParILUDataPrintLevel(ilu_data) = 0;
   jx_ParILUDataL1Norms(ilu_data) = NULL;
   jx_ParILUDataOperatorComplexity(ilu_data) = 0.;
   jx_ParILUDataIluType(ilu_data) = 0;
   jx_ParILUDataNLU(ilu_data) = 0;
   jx_ParILUDataNI(ilu_data) = 0;
   jx_ParILUDataUEnd(ilu_data) = NULL;

   /* Iterative setup variables */
   jx_ParILUDataIterativeSetupType(ilu_data) = 0;
   jx_ParILUDataIterativeSetupOption(ilu_data) = 0;
   jx_ParILUDataIterativeSetupMaxIter(ilu_data) = 100;
   jx_ParILUDataIterativeSetupNumIter(ilu_data) = 0;
   jx_ParILUDataIterativeSetupTolerance(ilu_data) = 1.e-6;
   jx_ParILUDataIterativeSetupHistory(ilu_data) = NULL;

   /* reordering_type default to use local RCM */
   jx_ParILUDataReorderingType(ilu_data) = 1;

   /* see jx_ILUSetType for more default values */
   jx_ParILUDataTestOption(ilu_data) = 0;

   /* -> General slots */
   jx_ParILUDataSchurSolverLogging(ilu_data) = 0;
   jx_ParILUDataSchurSolverPrintLevel(ilu_data) = 0;

   /* -> Schur-GMRES */
   jx_ParILUDataSchurGMRESKDim(ilu_data) = 5;
   jx_ParILUDataSchurGMRESMaxIter(ilu_data) = 5;
   jx_ParILUDataSchurGMRESTol(ilu_data) = 0.0;
   jx_ParILUDataSchurGMRESAbsoluteTol(ilu_data) = 0.0;
   jx_ParILUDataSchurGMRESRelChange(ilu_data) = 0;

   /* -> Schur precond data */
   jx_ParILUDataSchurPrecondIluType(ilu_data) = 0;
   jx_ParILUDataSchurPrecondIluLfil(ilu_data) = 0;
   jx_ParILUDataSchurPrecondIluMaxRowNnz(ilu_data) = 100;
   jx_ParILUDataSchurPrecondIluDroptol(ilu_data) = NULL;
   jx_ParILUDataSchurPrecondPrintLevel(ilu_data) = 0;
   jx_ParILUDataSchurPrecondMaxIter(ilu_data) = 1;
   jx_ParILUDataSchurPrecondTriSolve(ilu_data) = 1;
   jx_ParILUDataSchurPrecondLowerJacobiIters(ilu_data) = 5;
   jx_ParILUDataSchurPrecondUpperJacobiIters(ilu_data) = 5;
   jx_ParILUDataSchurPrecondTol(ilu_data) = 0.0;

   /* -> Schur-NSH */
   jx_ParILUDataSchurNSHSolveMaxIter(ilu_data) = 5;
   jx_ParILUDataSchurNSHSolveTol(ilu_data) = 0.0;
   jx_ParILUDataSchurNSHDroptol(ilu_data) = NULL;
   jx_ParILUDataSchurNSHMaxNumIter(ilu_data) = 2;
   jx_ParILUDataSchurNSHMaxRowNnz(ilu_data) = 1000;
   jx_ParILUDataSchurNSHTol(ilu_data) = 1e-09;
   jx_ParILUDataSchurMRMaxIter(ilu_data) = 2;
   jx_ParILUDataSchurMRColVersion(ilu_data) = 0;
   jx_ParILUDataSchurMRMaxRowNnz(ilu_data) = 200;
   jx_ParILUDataSchurMRTol(ilu_data) = 1e-09;

   return ilu_data;
}

/* Destroy */
JX_Int
jx_ILUDestroy(void *data)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)data;

   if (ilu_data)
   {
      /* final residual vector */
      jx_ParVectorDestroy(jx_ParILUDataResidual(ilu_data));
      jx_TFree(jx_ParILUDataRelResNorms(ilu_data));

      /* temp vectors for solve phase */
      jx_ParVectorDestroy(jx_ParILUDataUTemp(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataFTemp(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataXTemp(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataYTemp(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataZTemp(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataRhs(ilu_data));
      jx_ParVectorDestroy(jx_ParILUDataX(ilu_data));
      jx_TFree(jx_ParILUDataUExt(ilu_data));
      jx_TFree(jx_ParILUDataFExt(ilu_data));
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

      /* l1_norms */
      jx_TFree(jx_ParILUDataL1Norms(ilu_data));

      /* u_end */
      jx_TFree(jx_ParILUDataUEnd(ilu_data));

      /* Factors */
      jx_ParCSRMatrixDestroy(jx_ParILUDataMatS(ilu_data));
      jx_ParCSRMatrixDestroy(jx_ParILUDataMatL(ilu_data));
      jx_ParCSRMatrixDestroy(jx_ParILUDataMatU(ilu_data));
      jx_ParCSRMatrixDestroy(jx_ParILUDataMatLModified(ilu_data));
      jx_ParCSRMatrixDestroy(jx_ParILUDataMatUModified(ilu_data));
      jx_TFree(jx_ParILUDataMatD(ilu_data));
      jx_TFree(jx_ParILUDataMatDModified(ilu_data));

      // /* Schur solver */
      // if (jx_ParILUDataSchurSolver(ilu_data))
      // {
      //    switch (jx_ParILUDataIluType(ilu_data))
      //    {
      //    case 10:
      //    case 11:
      //    case 40:
      //    case 41:
      //    case 50:
      //       /* GMRES for Schur */
      //       JX_ParCSRGMRESDestroy(jx_ParILUDataSchurSolver(ilu_data));
      //       break;

      //    case 20:
      //    case 21:
      //       /* NSH for Schur */
      //       jx_NSHDestroy(jx_ParILUDataSchurSolver(ilu_data));
      //       break;

      //    default:
      //       break;
      //    }
      // }

      // /* ILU as precond for Schur */
      // if (jx_ParILUDataSchurPrecond(ilu_data) &&
      //     (jx_ParILUDataIluType(ilu_data) == 10 ||
      //      jx_ParILUDataIluType(ilu_data) == 11 ||
      //      jx_ParILUDataIluType(ilu_data) == 40 ||
      //      jx_ParILUDataIluType(ilu_data) == 41))
      // {
      //    JX_ILUDestroy(jx_ParILUDataSchurPrecond(ilu_data));
      // }

      /* CF marker array */
      jx_TFree(jx_ParILUDataCFMarkerArray(ilu_data));

      /* permutation array */
      jx_TFree(jx_ParILUDataPerm(ilu_data));
      jx_TFree(jx_ParILUDataQPerm(ilu_data));

      /* Iterative ILU data */
      jx_TFree(jx_ParILUDataIterativeSetupHistory(ilu_data));

      /* droptol array */
      jx_TFree(jx_ParILUDataDroptol(ilu_data));
      jx_TFree(jx_ParILUDataSchurPrecondIluDroptol(ilu_data));
      jx_TFree(jx_ParILUDataSchurNSHDroptol(ilu_data));
   }

   /* ILU data */
   jx_TFree(ilu_data);

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLevelOfFill
 *
 * Set fill level for ILUK
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLevelOfFill(void *ilu_vdata,
                     JX_Int lfil)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLfil(ilu_data) = lfil;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetMaxNnzPerRow
 *
 * Set max non-zeros per row in factors for ILUT
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetMaxNnzPerRow(void *ilu_vdata,
                      JX_Int nzmax)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataMaxRowNnz(ilu_data) = nzmax;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetDropThreshold
 *
 * Set threshold for dropping in LU factors for ILUT
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetDropThreshold(void *ilu_vdata,
                       JX_Real threshold)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   if (!(jx_ParILUDataDroptol(ilu_data)))
   {
      jx_ParILUDataDroptol(ilu_data) = jx_TAlloc(JX_Real, 3);
   }
   jx_ParILUDataDroptol(ilu_data)[0] = threshold;
   jx_ParILUDataDroptol(ilu_data)[1] = threshold;
   jx_ParILUDataDroptol(ilu_data)[2] = threshold;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetDropThresholdArray
 *
 * Set array of threshold for dropping in LU factors for ILUT
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetDropThresholdArray(void *ilu_vdata,
                            JX_Real *threshold)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   if (!(jx_ParILUDataDroptol(ilu_data)))
   {
      jx_ParILUDataDroptol(ilu_data) = jx_TAlloc(JX_Real, 3);
   }

   jx_ParILUDataDroptol(ilu_data)[0] = threshold[0];
   jx_ParILUDataDroptol(ilu_data)[1] = threshold[1];
   jx_ParILUDataDroptol(ilu_data)[2] = threshold[2];

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetType
 *
 * Set ILU factorization type
 *--------------------------------------------------------------------------*/
/* ILU(0) 可以走简化版本 */
JX_Int
jx_ILUSetType(void *ilu_vdata,
              JX_Int ilu_type)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;
   jx_ParILUDataIluType(ilu_data) = ilu_type;

   return jx_error_flag;
}

// JX_Int
// jx_ILUSetType(void *ilu_vdata, JX_Int ilu_type)
// {
//    jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

//    /* Destroy Schur solver if it already exists */
//    if (jx_ParILUDataSchurSolver(ilu_data))
//    {
//       switch (jx_ParILUDataIluType(ilu_data))
//       {
//       case 10:
//       case 11:
//       case 40:
//       case 41:
//       case 50:
//          /* GMRES for Schur */
//          JX_ParCSRGMRESDestroy(jx_ParILUDataSchurSolver(ilu_data));
//          break;

//       case 20:
//       case 21:
//          /* NSH for Schur */
//          jx_NSHDestroy(jx_ParILUDataSchurSolver(ilu_data));
//          break;

//       default:
//          break;
//       }

//       jx_ParILUDataSchurSolver(ilu_data) = NULL;
//    }

//    /* Destroy ILU preconditioner for Schur if it already exists */
//    if (jx_ParILUDataSchurPrecond(ilu_data) &&
//        (jx_ParILUDataIluType(ilu_data) == 10 ||
//         jx_ParILUDataIluType(ilu_data) == 11 ||
//         jx_ParILUDataIluType(ilu_data) == 40 ||
//         jx_ParILUDataIluType(ilu_data) == 41))
//    {
//       JX_ILUDestroy(jx_ParILUDataSchurPrecond(ilu_data));
//       jx_ParILUDataSchurPrecond(ilu_data) = NULL;
//    }

//    jx_ParILUDataIluType(ilu_data) = ilu_type;

//    /* Reset default values for Schur-related types */
//    switch (ilu_type)
//    {
//    case 20:
//    case 21:
//       if (!(jx_ParILUDataSchurNSHDroptol(ilu_data)))
//       {
//          jx_ParILUDataSchurNSHDroptol(ilu_data) =
//              jx_TAlloc(JX_Real, 2);

//          jx_ParILUDataSchurNSHDroptol(ilu_data)[0] = 1.0e-02;
//          jx_ParILUDataSchurNSHDroptol(ilu_data)[1] = 1.0e-02;
//          jx_ParILUDataSchurNSHOwnDroptolData(ilu_data) = 1;
//       }
//       break;

//    case 10:
//    case 11:
//    case 40:
//    case 41:
//    case 50:
//       if (!(jx_ParILUDataSchurPrecondIluDroptol(ilu_data)))
//       {
//          jx_ParILUDataSchurPrecondIluDroptol(ilu_data) =
//              jx_TAlloc(JX_Real, 3);

//          jx_ParILUDataSchurPrecondIluDroptol(ilu_data)[0] = 1.0e-02;
//          jx_ParILUDataSchurPrecondIluDroptol(ilu_data)[1] = 1.0e-02;
//          jx_ParILUDataSchurPrecondIluDroptol(ilu_data)[2] = 1.0e-02;
//          jx_ParILUDataSchurPrecondOwnDroptolData(ilu_data) = 1;
//       }
//       break;

//    default:
//       break;
//    }

//    return jx_error_flag;
// }

/*--------------------------------------------------------------------------
 * jx_ILUSetMaxIter
 *
 * Set max number of iterations for ILU solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetMaxIter(void *ilu_vdata,
                 JX_Int max_iter)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataMaxIter(ilu_data) = max_iter;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetIterativeSetupType
 *
 * Set iterative ILU setup algorithm
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetIterativeSetupType(void *ilu_vdata,
                            JX_Int iter_setup_type)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataIterativeSetupType(ilu_data) = iter_setup_type;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetIterativeSetupOption
 *
 * Set iterative ILU compute option
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetIterativeSetupOption(void *ilu_vdata,
                              JX_Int iter_setup_option)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   /* Compute residuals when using the stopping criteria,
    * if not chosen by the user.
    */
   iter_setup_option |= ((iter_setup_option & 0x02) && !(iter_setup_option & 0x0C)) ? 0x08 : 0;

   /* Compute residuals when asking for convergence history,
    * if not chosen by the user.
    */
   iter_setup_option |= ((iter_setup_option & 0x10) && !(iter_setup_option & 0x08)) ? 0x08 : 0;

   /* Zero out first bit of option. In JX this disables rocSPARSE logging.
    * Keep the behavior for compatibility even if the JX backend does not use it.
    */
   iter_setup_option &= ~0x01;

   jx_ParILUDataIterativeSetupOption(ilu_data) = iter_setup_option;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetIterativeSetupMaxIter
 *
 * Set maximum number of iterations for iterative ILU setup
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetIterativeSetupMaxIter(void *ilu_vdata,
                               JX_Int iter_setup_max_iter)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataIterativeSetupMaxIter(ilu_data) = iter_setup_max_iter;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetIterativeSetupTolerance
 *
 * Set dropping tolerance for iterative ILU setup
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetIterativeSetupTolerance(void *ilu_vdata,
                                 JX_Real iter_setup_tolerance)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataIterativeSetupTolerance(ilu_data) = iter_setup_tolerance;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUGetIterativeSetupHistory
 *
 * Get array of corrections and/or residual norms computed during ILU's
 * iterative setup algorithm.
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetIterativeSetupHistory(void *ilu_vdata,
                               JX_Real **iter_setup_history)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   *iter_setup_history = jx_ParILUDataIterativeSetupHistory(ilu_data);

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetTriSolve
 *
 * Set ILU triangular solver type
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetTriSolve(void *ilu_vdata,
                  JX_Int tri_solve)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataTriSolve(ilu_data) = tri_solve;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLowerJacobiIters
 *
 * Set Lower Jacobi iterations for iterative triangular solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLowerJacobiIters(void *ilu_vdata,
                          JX_Int lower_jacobi_iters)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLowerJacobiIters(ilu_data) = lower_jacobi_iters;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetUpperJacobiIters
 *
 * Set Upper Jacobi iterations for iterative triangular solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetUpperJacobiIters(void *ilu_vdata,
                          JX_Int upper_jacobi_iters)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataUpperJacobiIters(ilu_data) = upper_jacobi_iters;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetTol
 *
 * Set convergence tolerance for ILU solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetTol(void *ilu_vdata,
             JX_Real tol)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataTol(ilu_data) = tol;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetPrintLevel
 *
 * Set print level for ILU solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetPrintLevel(void *ilu_vdata,
                    JX_Int print_level)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataPrintLevel(ilu_data) = print_level;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLogging
 *
 * Set print level for ilu solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLogging(void *ilu_vdata,
                 JX_Int logging)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLogging(ilu_data) = logging;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLocalReordering
 *
 * Set type of reordering for local matrix
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLocalReordering(void *ilu_vdata,
                         JX_Int ordering_type)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataReorderingType(ilu_data) = ordering_type;

   return jx_error_flag;
}

// Schur 补类型的待补充 -----------------------------------------------------
// Schur 补类型的待补充 -----------------------------------------------------
/*--------------------------------------------------------------------------
 * jx_ILUSetLocalResidualCorrectionType
 *
 * 0: off, 1: global linear correction, 2: fixed-mask local linear correction,
 * 3: postprocess fixed Top-ratio local ILU(0) residual correction,
 * 4: interleaved fixed Top-ratio local ILU(0) residual correction.
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLocalResidualCorrectionType(void *ilu_vdata,
                                     JX_Int lrc_type)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionType(ilu_data) = lrc_type;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLocalResidualCorrectionMaxIter
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLocalResidualCorrectionMaxIter(void *ilu_vdata,
                                        JX_Int lrc_max_iter)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionMaxIter(ilu_data) = lrc_max_iter;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLocalResidualCorrectionThreshold
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLocalResidualCorrectionThreshold(void *ilu_vdata,
                                          JX_Real lrc_threshold)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionThreshold(ilu_data) = lrc_threshold;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetLocalResidualCorrectionOmega
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetLocalResidualCorrectionOmega(void *ilu_vdata,
                                      JX_Real lrc_omega)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionOmega(ilu_data) = lrc_omega;

   return jx_error_flag;
}

JX_Int
jx_ILUSetLocalResidualCorrectionSelectType(void *ilu_vdata,
                                           JX_Int lrc_select_type)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionSelectType(ilu_data) = lrc_select_type;

   return jx_error_flag;
}

JX_Int
jx_ILUSetLocalResidualCorrectionEnergyTarget(void *ilu_vdata,
                                             JX_Real lrc_energy_target)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataLocalResidualCorrectionEnergyTarget(ilu_data) = lrc_energy_target;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetSchurSolverKDIM
 *
 * Set KDim (for GMRES) for Solver of Schur System
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetSchurSolverKDIM(void *ilu_vdata,
                         JX_Int ss_kDim)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_ParILUDataSchurGMRESKDim(ilu_data) = ss_kDim;

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUSetSchurSolverMaxIter
 *
 * Set max iteration for Solver of Schur System
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUSetSchurSolverMaxIter(void *ilu_vdata,
                            JX_Int ss_max_iter)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   /* for the GMRES solve, the max iter is same as kdim by default */
   jx_ParILUDataSchurGMRESKDim(ilu_data) = ss_max_iter;
   jx_ParILUDataSchurGMRESMaxIter(ilu_data) = ss_max_iter;

   /* also set this value for NSH solve */
   jx_ParILUDataSchurNSHSolveMaxIter(ilu_data) = ss_max_iter;

   return jx_error_flag;
}

// /*--------------------------------------------------------------------------
//  * jx_ILUSetSchurSolverKDIM
//  *
//  * Set KDim (for GMRES) for Solver of Schur System
//  *--------------------------------------------------------------------------*/

// JX_Int
// jx_ILUSetSchurSolverKDIM(void *ilu_vdata,
//                          JX_Int ss_kDim)
// {
//    jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

//    jx_ParILUDataSchurGMRESKDim(ilu_data) = ss_kDim;

//    return jx_error_flag;
// }

// Schur 补类型的待补充 -----------------------------------------------------
// Schur 补类型的待补充 -----------------------------------------------------

/*--------------------------------------------------------------------------
 * jx_ILUGetNumIterations
 *
 * Get number of iterations for ILU solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetNumIterations(void *ilu_vdata,
                       JX_Int *num_iterations)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   if (!ilu_data)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }

   *num_iterations = jx_ParILUDataNumIterations(ilu_data);

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUGetFinalRelativeResidualNorm
 *
 * Get residual norms for ILU solver
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetFinalRelativeResidualNorm(void *ilu_vdata,
                                   JX_Real *res_norm)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   if (!ilu_data)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }

   *res_norm = jx_ParILUDataFinalRelResidualNorm(ilu_data);

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUWriteSolverParams
 *
 * Print solver params
 *
 * TODO (VPM): check runtime switch to decide whether running on host or device
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUWriteSolverParams(void *ilu_vdata)
{
   jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

   jx_printf("ILU Setup parameters: \n");
   jx_printf("ILU factorization type: %d : ",
             jx_ParILUDataIluType(ilu_data));

   switch (jx_ParILUDataIluType(ilu_data))
   {
   case 0:
      jx_printf("Block Jacobi with ILU(%d) \n",
                jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 1:
      jx_printf("Block Jacobi with ILUT \n");
      jx_printf("drop tolerance for B = %e, E&F = %e, S = %e \n",
                jx_ParILUDataDroptol(ilu_data)[0],
                jx_ParILUDataDroptol(ilu_data)[1],
                jx_ParILUDataDroptol(ilu_data)[2]);
      jx_printf("Max nnz per row = %d \n",
                jx_ParILUDataMaxRowNnz(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 10:
      jx_printf("ILU-GMRES with ILU(%d) \n",
                jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 11:
      jx_printf("ILU-GMRES with ILUT \n");
      jx_printf("drop tolerance for B = %e, E&F = %e, S = %e \n",
                jx_ParILUDataDroptol(ilu_data)[0],
                jx_ParILUDataDroptol(ilu_data)[1],
                jx_ParILUDataDroptol(ilu_data)[2]);
      jx_printf("Max nnz per row = %d \n",
                jx_ParILUDataMaxRowNnz(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 20:
      jx_printf("Newton-Schulz-Hotelling with ILU(%d) \n",
                jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 21:
      jx_printf("Newton-Schulz-Hotelling with ILUT \n");
      jx_printf("drop tolerance for B = %e, E&F = %e, S = %e \n",
                jx_ParILUDataDroptol(ilu_data)[0],
                jx_ParILUDataDroptol(ilu_data)[1],
                jx_ParILUDataDroptol(ilu_data)[2]);
      jx_printf("Max nnz per row = %d \n",
                jx_ParILUDataMaxRowNnz(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 30:
      jx_printf("RAS with ILU(%d) \n",
                jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n", jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 31:
      jx_printf("RAS with ILUT \n");
      jx_printf("drop tolerance for B = %e, E&F = %e, S = %e \n",
                jx_ParILUDataDroptol(ilu_data)[0],
                jx_ParILUDataDroptol(ilu_data)[1],
                jx_ParILUDataDroptol(ilu_data)[2]);
      jx_printf("Max nnz per row = %d \n", jx_ParILUDataMaxRowNnz(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n", jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 40:
      jx_printf("ddPQ-ILU-GMRES with ILU(%d) \n",
                jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 41:
      jx_printf("ddPQ-ILU-GMRES with ILUT \n");
      jx_printf("drop tolerance for B = %e, E&F = %e, S = %e \n",
                jx_ParILUDataDroptol(ilu_data)[0],
                jx_ParILUDataDroptol(ilu_data)[1],
                jx_ParILUDataDroptol(ilu_data)[2]);
      jx_printf("Max nnz per row = %d \n",
                jx_ParILUDataMaxRowNnz(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n",
                jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   case 50:
      jx_printf("RAP-Modified-ILU with ILU(%d) \n", jx_ParILUDataLfil(ilu_data));
      jx_printf("Operator Complexity (Fill factor) = %f \n", jx_ParILUDataOperatorComplexity(ilu_data));
      break;

   default:
      jx_printf("Unknown type \n");
      break;
   }

   jx_printf("\n ILU Solver Parameters: \n");
   jx_printf("Max number of iterations: %d\n", jx_ParILUDataMaxIter(ilu_data));

   if (jx_ParILUDataTriSolve(ilu_data))
   {
      jx_printf("  Triangular solver type: exact (1)\n");
   }
   else
   {
      jx_printf("  Triangular solver type: iterative (0)\n");
      jx_printf(" Lower Jacobi Iterations: %d\n", jx_ParILUDataLowerJacobiIters(ilu_data));
      jx_printf(" Upper Jacobi Iterations: %d\n", jx_ParILUDataUpperJacobiIters(ilu_data));
   }

   jx_printf("      Stopping tolerance: %e\n", jx_ParILUDataTol(ilu_data));

   return jx_error_flag;
}

/******************************************************************************
 *
 * ILU helper functions
 *
 * TODO (VPM): move these to a new "par_ilu_utils.c" file
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * jx_ILUMinHeapAddI
 *
 * Add an element to the heap
 * I means JX_Int
 * R means JX_Real
 * max/min heap
 * r means heap goes from 0 to -1, -2 instead of 0 1 2
 * Ii and Ri means orderd by value of heap, like iw for ILU
 * heap: array of that heap
 * len: the current length of the heap
 * WARNING: You should first put that element to the end of the heap
 *    and add the length of heap by one before call this function.
 * the reason is that we don't want to change something outside the
 *    heap, so left it to the user
 *--------------------------------------------------------------------------*/
JX_Int
jx_ILUMinHeapAddI(JX_Int *heap, JX_Int len)
{
   /* parent, left, right */
   JX_Int p;

   len--; /* now len is the current index */
   while (len > 0)
   {
      /* get the parent index */
      p = (len - 1) / 2;
      if (heap[p] > heap[len])
      {
         /* this is smaller */
         jx_swap(heap, p, len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUMinHeapAddIIIi
 *
 * See jx_ILUMinHeapAddI for detail instructions
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUMinHeapAddIIIi(JX_Int *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p;

   len--; /* now len is the current index */
   while (len > 0)
   {
      /* get the parent index */
      p = (len - 1) / 2;
      if (heap[p] > heap[len])
      {
         /* this is smaller */
         jx_swap(Ii1, heap[p], heap[len]);
         jx_swap2i(heap, I1, p, len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUMinHeapAddIRIi
 *
 * see jx_ILUMinHeapAddI for detail instructions
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUMinHeapAddIRIi(JX_Int *heap, JX_Real *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p;

   len--; /* now len is the current index */
   while (len > 0)
   {
      /* get the parent index */
      p = (len - 1) / 2;
      if (heap[p] > heap[len])
      {
         /* this is smaller */
         jx_swap(Ii1, heap[p], heap[len]);
         jx_swap2(heap, I1, p, len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/* see jx_ILUMinHeapAddI for detail instructions */
JX_Int
jx_ILUMaxHeapAddRabsIIi(JX_Real *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p;
   len--; /* now len is the current index */
   while (len > 0)
   {
      /* get the parent index */
      p = (len - 1) / 2;
      if (jx_abs(heap[p]) < jx_abs(heap[len]))
      {
         /* this is smaller */
         jx_swap(Ii1, heap[p], heap[len]);
         jx_swap2(I1, heap, p, len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/* see jx_ILUMinHeapAddI for detail instructions */
JX_Int
jx_ILUMaxrHeapAddRabsI(JX_Real *heap, JX_Int *I1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p;
   len--; /* now len is the current index */
   while (len > 0)
   {
      /* get the parent index */
      p = (len - 1) / 2;
      if (jx_abs(heap[-p]) < jx_abs(heap[-len]))
      {
         /* this is smaller */
         jx_swap2(I1, heap, -p, -len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUMinHeapRemoveI
 *
 * Swap the first element with the last element of the heap,
 *    reduce size by one, and maintain the heap structure
 * I means JX_Int
 * R means JX_Real
 * max/min heap
 * r means heap goes from 0 to -1, -2 instead of 0 1 2
 * Ii and Ri means orderd by value of heap, like iw for ILU
 * heap: aray of that heap
 * len: current length of the heap
 * WARNING: Remember to change the len yourself
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUMinHeapRemoveI(JX_Int *heap, JX_Int len)
{
   /* parent, left, right */
   JX_Int p, l, r;
   len--; /* now len is the max index */
   /* swap the first element to last */
   jx_swap(heap, 0, len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while (l < len)
   {
      r = 2 * p + 2;

      /* two childs, pick the smaller one */
      l = r >= len || heap[l] < heap[r] ? l : r;
      if (heap[l] < heap[p])
      {
         jx_swap(heap, l, p);
         p = l;
         l = 2 * p + 1;
      }
      else
      {
         break;
      }
   }

   return jx_error_flag;
}

/* see jx_ILUMinHeapRemoveI for detail instructions */
JX_Int
jx_ILUMinHeapRemoveIIIi(JX_Int *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p, l, r;
   len--; /* now len is the max index */
   /* swap the first element to last */
   jx_swap(Ii1, heap[0], heap[len]);
   jx_swap2i(heap, I1, 0, len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while (l < len)
   {
      r = 2 * p + 2;
      /* two childs, pick the smaller one */
      l = r >= len || heap[l] < heap[r] ? l : r;
      if (heap[l] < heap[p])
      {
         jx_swap(Ii1, heap[p], heap[l]);
         jx_swap2i(heap, I1, l, p);
         p = l;
         l = 2 * p + 1;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/* see jx_ILUMinHeapRemoveI for detail instructions */
JX_Int
jx_ILUMinHeapRemoveIRIi(JX_Int *heap, JX_Real *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p, l, r;
   len--; /* now len is the max index */
   /* swap the first element to last */
   jx_swap(Ii1, heap[0], heap[len]);
   jx_swap2(heap, I1, 0, len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while (l < len)
   {
      r = 2 * p + 2;
      /* two childs, pick the smaller one */
      l = r >= len || heap[l] < heap[r] ? l : r;
      if (heap[l] < heap[p])
      {
         jx_swap(Ii1, heap[p], heap[l]);
         jx_swap2(heap, I1, l, p);
         p = l;
         l = 2 * p + 1;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/* see jx_ILUMinHeapRemoveI for detail instructions */
JX_Int
jx_ILUMaxHeapRemoveRabsIIi(JX_Real *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p, l, r;
   len--; /* now len is the max index */
   /* swap the first element to last */
   jx_swap(Ii1, heap[0], heap[len]);
   jx_swap2(I1, heap, 0, len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while (l < len)
   {
      r = 2 * p + 2;
      /* two childs, pick the smaller one */
      l = r >= len || jx_abs(heap[l]) > jx_abs(heap[r]) ? l : r;
      if (jx_abs(heap[l]) > jx_abs(heap[p]))
      {
         jx_swap(Ii1, heap[p], heap[l]);
         jx_swap2(I1, heap, l, p);
         p = l;
         l = 2 * p + 1;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/* see jx_ILUMinHeapRemoveI for detail instructions */
JX_Int
jx_ILUMaxrHeapRemoveRabsI(JX_Real *heap, JX_Int *I1, JX_Int len)
{
   /* parent, left, right */
   JX_Int p, l, r;
   len--; /* now len is the max index */
   /* swap the first element to last */
   jx_swap2(I1, heap, 0, -len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while (l < len)
   {
      r = 2 * p + 2;
      /* two childs, pick the smaller one */
      l = r >= len || jx_abs(heap[-l]) > jx_abs(heap[-r]) ? l : r;
      if (jx_abs(heap[-l]) > jx_abs(heap[-p]))
      {
         jx_swap2(I1, heap, -l, -p);
         p = l;
         l = 2 * p + 1;
      }
      else
      {
         break;
      }
   }
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUMaxQSplitRabsI
 *
 * Split based on quick sort algorithm (avoid sorting the entire array)
 * find the largest k elements out of original array
 *
 * arrayR: input array for compare
 * arrayI: integer array bind with array
 * k: largest k elements
 * len: length of the array
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUMaxQSplitRabsI(JX_Real *array, JX_Int *I, JX_Int left, JX_Int bound, JX_Int right)
{
   JX_Int i, last;
   if (left >= right)
   {
      return jx_error_flag;
   }
   jx_swap2(I, array, left, (left + right) / 2);
   last = left;
   for (i = left + 1; i <= right; i++)
   {
      if (jx_abs(array[i]) > jx_abs(array[left]))
      {
         jx_swap2(I, array, ++last, i);
      }
   }
   jx_swap2(I, array, left, last);
   jx_ILUMaxQSplitRabsI(array, I, left, bound, last - 1);
   if (bound > last)
   {
      jx_ILUMaxQSplitRabsI(array, I, last + 1, bound, right);
   }

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUMaxRabs
 *
 * Helper function to search max value from a row
 * array: the array we work on
 * start: the start of the search range
 * end: the end of the search range
 * nLU: ignore rows (new row index) after nLU
 * rperm: reverse permutation array rperm[old] = new.
 *        if rperm set to NULL, ingore nLU and rperm
 * value: return the value ge get (absolute value)
 * index: return the index of that value, could be NULL which means not return
 * l1_norm: return the l1_norm of the array, could be NULL which means no return
 * nnz: return the number of nonzeros inside this array, could be NULL which means no return
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUMaxRabs(JX_Real *array_data, JX_Int *array_j, JX_Int start, JX_Int end, JX_Int nLU, JX_Int *rperm, JX_Real *value, JX_Int *index, JX_Real *l1_norm, JX_Int *nnz)
{
   JX_Int i, idx, col;
   JX_Real val, max_value, norm, nz;

   nz = 0;
   norm = 0.0;
   max_value = -1.0;
   idx = -1;
   if (rperm)
   {
      /* apply rperm and nLU */
      for (i = start; i < end; i++)
      {
         col = rperm[array_j[i]];
         if (col > nLU)
         {
            /* this old column is in new external part */
            continue;
         }
         nz++;
         val = jx_abs(array_data[i]);
         norm += val;
         if (max_value < val)
         {
            max_value = val;
            idx = i;
         }
      }
   }
   else
   {
      /* basic search */
      for (i = start; i < end; i++)
      {
         val = jx_abs(array_data[i]);
         norm += val;
         if (max_value < val)
         {
            max_value = val;
            idx = i;
         }
      }
      nz = end - start;
   }

   *value = max_value;
   if (index)
   {
      *index = idx;
   }
   if (l1_norm)
   {
      *l1_norm = norm;
   }
   if (nnz)
   {
      *nnz = nz;
   }

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUGetPermddPQPre
 *
 * Pre selection for ddPQ, this is the basic version considering row sparsity
 * n: size of matrix
 * nLU: size we consider ddPQ reorder, only first nLU*nLU block is considered
 * A_diag_i/j/data: information of A
 * tol: tol for ddPQ, normally between 0.1-0.3
 * *perm: current row order
 * *rperm: current column order
 * *pperm_pre: output ddPQ pre row roder
 * *qperm_pre: output ddPQ pre column order
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetPermddPQPre(JX_Int n, JX_Int nLU, JX_Int *A_diag_i, JX_Int *A_diag_j, JX_Real *A_diag_data, JX_Real tol, JX_Int *perm, JX_Int *rperm,
                     JX_Int *pperm_pre, JX_Int *qperm_pre, JX_Int *nB)
{
   JX_Int i, ii, nB_pre, k1, k2;
   JX_Real gtol, max_value, norm;

   JX_Int *jcol, *jnnz;
   JX_Real *weight;

   weight = jx_TAlloc(JX_Real, nLU + 1);
   jcol = jx_TAlloc(JX_Int, nLU + 1);
   jnnz = jx_TAlloc(JX_Int, nLU + 1);

   max_value = -1.0;
   /* first need to build gtol */
   for (ii = 0; ii < nLU; ii++)
   {
      /* find real row */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i + 1];
      /* find max|a| of that row and its index */
      jx_ILUMaxRabs(A_diag_data, A_diag_j, k1, k2, nLU, rperm, weight + ii, jcol + ii, &norm, jnnz + ii);
      weight[ii] /= norm;
      if (weight[ii] > max_value)
      {
         max_value = weight[ii];
      }
   }

   gtol = tol * max_value;

   /* second loop to pre select B */
   nB_pre = 0;
   for (ii = 0; ii < nLU; ii++)
   {
      /* keep this row */
      if (weight[ii] > gtol)
      {
         weight[nB_pre] /= (JX_Real)(jnnz[ii]);
         pperm_pre[nB_pre] = perm[ii];
         qperm_pre[nB_pre++] = A_diag_j[jcol[ii]];
      }
   }

   *nB = nB_pre;

   /* sort from small to large */
   jx_qsort3(weight, pperm_pre, qperm_pre, 0, nB_pre - 1);

   jx_TFree(weight);
   jx_TFree(jcol);
   jx_TFree(jnnz);

   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUGetPermddPQ
 *
 * Get ddPQ version perm array for ParCSR matrices. ddPQ is a two-side
 * permutation for diagonal dominance. Greedy matching selection
 *
 * Parameters:
 *   A: the input matrix
 *   pperm: row permutation (lives at memory_location_A)
 *   qperm: col permutation (lives at memory_location_A)
 *   nB: the size of B block
 *   nI: number of interial nodes
 *   tol: the dropping tolorance for ddPQ
 *   reordering_type: Type of reordering for the interior nodes.
 *
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 *
 * TODO (VPM): Change permutation arrays types to jx_IntArray
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetPermddPQ(jx_ParCSRMatrix *A, JX_Int **io_pperm, JX_Int **io_qperm, JX_Real tol, JX_Int *nB, JX_Int *nI, JX_Int reordering_type)
{
   JX_Int i, nB_pre, irow, jcol, nLU;
   JX_Int *pperm, *qperm;
   JX_Int *rpperm, *rqperm, *pperm_pre, *qperm_pre;

   /* data objects for A */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);

   /* problem size */
   JX_Int n = jx_CSRMatrixNumRows(A_diag);

   /* 1: Setup and create memory
    */

   pperm = NULL;
   qperm = jx_TAlloc(JX_Int, n);
   rpperm = jx_TAlloc(JX_Int, n);
   rqperm = jx_TAlloc(JX_Int, n);

   /* 2: Find interior nodes first
    */
   jx_ILUGetInteriorExteriorPerm(A, &pperm, &nLU, 0);
   *nI = nLU;

   /* 3: Pre selection on interial nodes
    * this pre selection puts external nodes to the last
    * also provide candidate rows for B block
    */

   /* build reverse permutation array
    * rperm[old] = new
    */
   for (i = 0; i < n; i++)
   {
      rpperm[pperm[i]] = i;
   }

   /* build place holder for pre selection pairs */
   pperm_pre = jx_TAlloc(JX_Int, nLU);
   qperm_pre = jx_TAlloc(JX_Int, nLU);

   /* pre selection */
   jx_ILUGetPermddPQPre(n, nLU, A_diag_i, A_diag_j, A_diag_data, tol, pperm, rpperm, pperm_pre, qperm_pre, &nB_pre);

   /* 4: Build B block
    * Greedy selection
    */

   /* rperm[old] = new */
   for (i = 0; i < nLU; i++)
   {
      rpperm[pperm[i]] = -1;
   }

   jx_TMemcpy(rqperm, rpperm, JX_Int, n);
   jx_TMemcpy(qperm, pperm, JX_Int, n);

   /* we sort from small to large, so we need to go from back to start
    * we only need nB_pre to start the loop, after that we could use it for size of B
    */
   for (i = nB_pre - 1, nB_pre = 0; i >= 0; i--)
   {
      irow = pperm_pre[i];
      jcol = qperm_pre[i];

      /* this col is not yet taken */
      if (rqperm[jcol] < 0)
      {
         rpperm[irow] = nB_pre;
         rqperm[jcol] = nB_pre;
         pperm[nB_pre] = irow;
         qperm[nB_pre++] = jcol;
      }
   }

   /* 5: Complete the permutation
    * rperm[old] = new
    * those still mapped to a new index means not yet covered
    */
   nLU = nB_pre;
   for (i = 0; i < n; i++)
   {
      if (rpperm[i] < 0)
      {
         pperm[nB_pre++] = i;
      }
   }
   nB_pre = nLU;
   for (i = 0; i < n; i++)
   {
      if (rqperm[i] < 0)
      {
         qperm[nB_pre++] = i;
      }
   }

   /* Finishing up and free
    */

   switch (reordering_type)
   {
   case 0:
      /* no RCM in this case */
      break;
   case 1:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, nLU, &pperm, &qperm, 0);
      break;
   default:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, nLU, &pperm, &qperm, 0);
      break;
   }

   *nB = nLU;
   *io_pperm = pperm;
   *io_qperm = qperm;

   jx_TFree(rpperm);
   jx_TFree(rqperm);
   jx_TFree(pperm_pre);
   jx_TFree(qperm_pre);
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_ILUGetInteriorExteriorPerm
 *
 * Get perm array from parcsr matrix based on diag and offdiag matrix
 * Just simply loop through the rows of offd of A, check for nonzero rows
 * Put interior nodes at the beginning
 *
 * Parameters:
 *   A: parcsr matrix
 *   perm: permutation array
 *   nLU: number of interial nodes
 *   reordering_type: Type of (additional) reordering for the interior nodes.
 *
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUGetInteriorExteriorPerm(jx_ParCSRMatrix *A, JX_Int **perm, JX_Int *nLU, JX_Int reordering_type)
{
   /* get basic information of A */
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int i, j, first, last, start, end;
   JX_Int num_sends, send_map_start, send_map_end, col;
   jx_CSRMatrix *A_offd;
   JX_Int *A_offd_i;
   A_offd = jx_ParCSRMatrixOffd(A);
   A_offd_i = jx_CSRMatrixI(A_offd);
   first = 0;
   last = n - 1;
   JX_Int *temp_perm = jx_TAlloc(JX_Int, n);
   JX_Int *marker = jx_CTAlloc(JX_Int, n);

   /* first get col nonzero from com_pkg */
   /* get comm_pkg, craete one if we not yet have one */
   jx_ParCSRCommPkg *comm_pkg = jx_ParCSRMatrixCommPkg(A);
   if (!comm_pkg)
   {
      jx_MatvecCommPkgCreate(A);
      comm_pkg = jx_ParCSRMatrixCommPkg(A);
   }

   /* now directly take adavantage of comm_pkg */
   num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
   for (i = 0; i < num_sends; i++)
   {
      send_map_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
      send_map_end = jx_ParCSRCommPkgSendMapStart(comm_pkg, i + 1);
      for (j = send_map_start; j < send_map_end; j++)
      {
         col = jx_ParCSRCommPkgSendMapElmt(comm_pkg, j);
         if (marker[col] == 0)
         {
            temp_perm[last--] = col;
            marker[col] = -1;
         }
      }
   }

   /* now deal with the row */
   for (i = 0; i < n; i++)
   {
      if (marker[i] == 0)
      {
         start = A_offd_i[i];
         end = A_offd_i[i + 1];
         if (start == end)
         {
            temp_perm[first++] = i;
         }
         else
         {
            temp_perm[last--] = i;
         }
      }
   }
   switch (reordering_type)
   {
   case 0:
      /* no RCM in this case */
      break;
   case 1:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, first, &temp_perm, &temp_perm, 1);
      break;
   default:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, first, &temp_perm, &temp_perm, 1);
      break;
   }

   /* set out values */
   *nLU = first;
   if ((*perm) != NULL)
      jx_TFree(*perm);
   *perm = temp_perm;

   jx_TFree(marker);
   return jx_error_flag;
}

/* 对比修改到这里 ---------------------*/
/*
 * Get the (local) ordering of the diag (local) matrix (no permutation). This is the permutation used for the block-jacobi case
 * A: parcsr matrix
 * perm: permutation array
 * nLU: number of interior nodes
 * reordering_type: Type of (additional) reordering for the nodes.
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 */
JX_Int
jx_ILUGetLocalPerm(jx_ParCSRMatrix *A, JX_Int **perm, JX_Int *nLU, JX_Int reordering_type)
{
   /* get basic information of A */
   JX_Int n = jx_ParCSRMatrixNumRows(A);
   JX_Int i;
   JX_Int *temp_perm = jx_TAlloc(JX_Int, n);

   /* set perm array */
   for (i = 0; i < n; i++)
   {
      temp_perm[i] = i;
   }
   switch (reordering_type)
   {
   case 0:
      /* no RCM in this case */
      break;
   case 1:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, n, &temp_perm, &temp_perm, 1);
      break;
   default:
      /* RCM */
      jx_ILULocalRCM(jx_ParCSRMatrixDiag(A), 0, n, &temp_perm, &temp_perm, 1);
      break;
   }
   *nLU = n;
   if ((*perm) != NULL)
      jx_TFree(*perm);
   *perm = temp_perm;

   return jx_error_flag;
}

/* Build the expanded matrix for RAS-1
 * A: input ParCSR matrix
 * E_i, E_j, E_data: information for external matrix
 * rperm: reverse permutation to build real index, rperm[old] = new
 */
JX_Int
jx_ILUBuildRASExternalMatrix(jx_ParCSRMatrix *A, JX_Int *rperm, JX_Int **E_i, JX_Int **E_j, JX_Real **E_data)
{
   JX_Int i, j, idx;
   JX_Int big_col;

   /* data objects for communication */
   MPI_Comm comm = jx_ParCSRMatrixComm(A);
   JX_Int my_id;

   /* data objects for A */
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_col_starts = jx_ParCSRMatrixColStarts(A);
   JX_Int *A_offd_colmap = jx_ParCSRMatrixColMapOffd(A);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);

   /* data objects for external A matrix */
   // Need to check the new version of jx_ParcsrGetExternalRows
   jx_CSRMatrix *A_ext = NULL;
   // # up to local offd cols, no need to be JX_Int
   JX_Int *A_ext_i = NULL;
   // Return global index, JX_Int required
   JX_Int *A_ext_j = NULL;
   JX_Real *A_ext_data = NULL;

   /* data objects for output */
   JX_Int E_nnz;
   JX_Int *E_ext_i = NULL;
   // Local index, no need to use JX_Int
   JX_Int *E_ext_j = NULL;
   JX_Real *E_ext_data = NULL;

   // guess non-zeros for E before start
   JX_Int E_init_alloc;

   /* size */
   JX_Int n = jx_CSRMatrixNumCols(A_diag);
   JX_Int m = jx_CSRMatrixNumCols(A_offd);
   JX_Int A_diag_nnz = A_diag_i[n];
   JX_Int A_offd_nnz = A_offd_i[n];

   /* 1: Set up phase and get external rows
    * Use the JX build-in function
    */

   /* MPI stuff */
   // jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);

   /* Param of jx_ParcsrGetExternalRows:
    * jx_ParCSRMatrix   *A          [in]  -> Input parcsr matrix.
    * JX_Int            indies_len  [in]  -> Input length of indices_len array
    * JX_Int            *indices    [in]  -> Input global indices of rows we want to get
    * jx_CSRMatrix      **A_ext     [out] -> Return the external CSR matrix.
    * jx_ParCSRCommPkg  commpkg_out [out] -> Return commpkg if set to a point. Use NULL here since we don't want it.
    */
   //   jx_ParcsrGetExternalRows( A, m, A_offd_colmap, &A_ext, NULL );
   A_ext = jx_ParCSRMatrixExtractBExt(A, A, 1);

   A_ext_i = jx_CSRMatrixI(A_ext);
   // This should be JX_Int since this is global index, use big_j in csr */
   A_ext_j = jx_CSRMatrixJ(A_ext);
   A_ext_data = jx_CSRMatrixData(A_ext);

   /* guess memory we need to allocate to E_j */
   E_init_alloc = jx_max((JX_Int)(A_diag_nnz / (JX_Real)n / (JX_Real)n * (JX_Real)m * (JX_Real)m + A_offd_nnz), 1);

   /* Initial guess */
   E_ext_i = jx_TAlloc(JX_Int, m + 1);
   E_ext_j = jx_TAlloc(JX_Int, E_init_alloc);
   E_ext_data = jx_TAlloc(JX_Real, E_init_alloc);

   /* 2: Discard unecessary cols
    * Search A_ext_j, discard those cols not belong to current proc
    * First check diag, and search in offd_col_map
    */

   E_nnz = 0;
   E_ext_i[0] = 0;

   for (i = 0; i < m; i++)
   {
      E_ext_i[i] = E_nnz;
      for (j = A_ext_i[i]; j < A_ext_i[i + 1]; j++)
      {
         big_col = A_ext_j[j];
         /* First check if that belongs to the diagonal part */
#ifdef JX_NO_GLOBAL_PARTITION

         if (big_col >= A_col_starts[0] && big_col < A_col_starts[1])
         {
            /* this is a diagonal entry, rperm (map old to new) and shift it */

            /* Note here, the result of big_col - A_col_starts[0] in no longer a JX_Int */
            idx = (JX_Int)(big_col - A_col_starts[0]);
            E_ext_j[E_nnz] = rperm[idx];
            E_ext_data[E_nnz++] = A_ext_data[j];
         }

#else
         if (big_col >= A_col_starts[my_id] && big_col < A_col_starts[my_id + 1])
         {
            /* this is a diagonal entry, rperm (map old to new) and shift it */

            /* Note here, the result of big_col - A_col_starts[0] in no longer a JX_Int */
            idx = (JX_Int)(big_col - A_col_starts[my_id]);
            E_ext_j[E_nnz] = rperm[idx];
            E_ext_data[E_nnz++] = A_ext_data[j];
         }
#endif
         /* If not, apply binary search to check if is offdiagonal */
         else
         {
            /* Search, result is not JX_Int */
            E_ext_j[E_nnz] = jx_BinarySearch(A_offd_colmap, big_col, m);
            if (E_ext_j[E_nnz] >= 0)
            {
               /* this is an offdiagonal entry */
               E_ext_j[E_nnz] = E_ext_j[E_nnz] + n;
               E_ext_data[E_nnz++] = A_ext_data[j];
            }
            else
            {
               /* skip capacity check */
               continue;
            }
         }
         /* capacity check, allocate new memory when full */
         if (E_nnz >= E_init_alloc)
         {
            E_init_alloc = E_init_alloc * EXPAND_FACT + 1;
            E_ext_j = jx_TReAlloc(E_ext_j, JX_Int, E_init_alloc);
            E_ext_data = jx_TReAlloc(E_ext_data, JX_Real, E_init_alloc);
         }
      }
   }
   E_ext_i[m] = E_nnz;

   /* 3: Free and finish up
    * Free memory, set E_i, E_j and E_data
    */

   *E_i = E_ext_i;
   *E_j = E_ext_j;
   *E_data = E_ext_data;

   jx_CSRMatrixDestroy(A_ext);

   return jx_error_flag;
}

/* This function sort offdiagonal map as well as J array for offdiagonal part
 * A: The input CSR matrix
 */
JX_Int
jx_ILUSortOffdColmap(jx_ParCSRMatrix *A)
{
   JX_Int i;
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);
   JX_Int *A_offd_colmap = jx_ParCSRMatrixColMapOffd(A);
   JX_Int len = jx_CSRMatrixNumCols(A_offd);
   JX_Int nnz = jx_CSRMatrixNumNonzeros(A_offd);
   JX_Int *perm = jx_TAlloc(JX_Int, len);
   JX_Int *rperm = jx_TAlloc(JX_Int, len);

   for (i = 0; i < len; i++)
   {
      perm[i] = i;
   }

   jx_qsort2i(A_offd_colmap, perm, 0, len - 1);

   for (i = 0; i < len; i++)
   {
      rperm[perm[i]] = i;
   }

   for (i = 0; i < nnz; i++)
   {
      A_offd_j[i] = rperm[A_offd_j[i]];
   }

   jx_TFree(perm);
   jx_TFree(rperm);

   return jx_error_flag;
}

/* This function computes the RCM ordering of a sub matrix of
 * sparse matrix B = A(perm,perm)
 * For nonsymmetrix problem, is the RCM ordering of B + B'
 * A: The input CSR matrix
 * start:      the start position of the submatrix in B
 * end:        the end position of the submatrix in B ( exclude end, [start,end) )
 * permp:      pointer to the row permutation array such that B = A(perm, perm)
 *             point to NULL if you want to work directly on A
 *             on return, permp will point to the new permutation where
 *             in [start, end) the matrix will reordered
 * qpermp:     pointer to the col permutation array such that B = A(perm, perm)
 *             point to NULL or equal to permp if you want symmetric order
 *             on return, qpermp will point to the new permutation where
 *             in [start, end) the matrix will reordered
 * sym:        set to nonzero to work on A only(symmetric), otherwise A + A'.
 *             WARNING: if you use non-symmetric reordering, that is,
 *             different row and col reordering, the resulting A might be non-symmetric.
 *             Be careful if you are using non-symmetric reordering
 */
JX_Int
jx_ILULocalRCM(jx_CSRMatrix *A, JX_Int start, JX_Int end,
               JX_Int **permp, JX_Int **qpermp, JX_Int sym)
{
   JX_Int i, j, row, col, r1, r2;

   JX_Int num_nodes = end - start;
   JX_Int n = jx_CSRMatrixNumRows(A);
   JX_Int ncol = jx_CSRMatrixNumCols(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   jx_CSRMatrix *GT = NULL;
   jx_CSRMatrix *GGT = NULL;
   //    JX_Int               *AAT_i         = NULL;
   //    JX_Int               *AAT_j         = NULL;
   JX_Int A_nnz = jx_CSRMatrixNumNonzeros(A);
   jx_CSRMatrix *G = NULL;
   JX_Int *G_i = NULL;
   JX_Int *G_j = NULL;
   JX_Real *G_data = NULL;
   JX_Int *G_perm = NULL;
   JX_Int G_nnz;
   JX_Int G_capacity;
   JX_Int *perm_temp = NULL;
   JX_Int *perm = *permp;
   JX_Int *qperm = *qpermp;
   JX_Int *rqperm = NULL;

   /* 1: Preprosessing
    * Check error in input, set some parameters
    */
   if (num_nodes <= 0)
   {
      /* don't do this if we are too small */
      return jx_error_flag;
   }
   if (n != ncol || end > n || start < 0)
   {
      /* don't do this if the input has error */
      jx_printf("Error input, abort RCM\n");
      return jx_error_flag;
   }
   if (!perm)
   {
      /* create permutation array if we don't have one yet */
      perm = jx_TAlloc(JX_Int, n);
      for (i = 0; i < n; i++)
      {
         perm[i] = i;
      }
   }
   if (!qperm)
   {
      /* symmetric reordering, just point it to row reordering */
      qperm = perm;
   }
   rqperm = jx_TAlloc(JX_Int, n);
   for (i = 0; i < n; i++)
   {
      rqperm[qperm[i]] = i;
   }
   /* 2: Build Graph
    * Build Graph for RCM ordering
    */
   G = jx_CSRMatrixCreate(num_nodes, num_nodes, 0);
   jx_CSRMatrixInitialize(G);
   jx_CSRMatrixOwnsData(G) = 1;
   G_i = jx_CSRMatrixI(G);
   if (sym)
   {
      /* Directly use A */
      G_nnz = 0;
      G_capacity = jx_max(A_nnz * n * n / num_nodes / num_nodes - num_nodes, 1);
      G_j = jx_TAlloc(JX_Int, G_capacity);
      for (i = 0; i < num_nodes; i++)
      {
         G_i[i] = G_nnz;
         row = perm[i + start];
         r1 = A_i[row];
         r2 = A_i[row + 1];
         for (j = r1; j < r2; j++)
         {
            col = rqperm[A_j[j]];
            if (col != row && col >= start && col < end)
            {
               /* this is an entry in G */
               G_j[G_nnz++] = col - start;
               if (G_nnz >= G_capacity)
               {
                  JX_Int tmp = G_capacity;
                  G_capacity = G_capacity * EXPAND_FACT + 1;
                  G_j = jx_TReAlloc_v2(G_j, JX_Int, tmp, JX_Int, G_capacity);
               }
            }
         }
      }
      G_i[num_nodes] = G_nnz;
      if (G_nnz == 0)
      {
         // G has only diagonal, no need to do any kind of RCM
         jx_TFree(G_j);
         jx_TFree(rqperm);
         *permp = perm;
         *qpermp = qperm;
         jx_CSRMatrixDestroy(G);
         return jx_error_flag;
      }
      jx_CSRMatrixJ(G) = G_j;
      jx_CSRMatrixNumNonzeros(G) = G_nnz;
   }
   else
   {
      /* Use A + A' */
      G_nnz = 0;
      G_capacity = jx_max(A_nnz * n * n / num_nodes / num_nodes - num_nodes, 1);
      G_j = jx_TAlloc(JX_Int, G_capacity);
      for (i = 0; i < num_nodes; i++)
      {
         G_i[i] = G_nnz;
         row = perm[i + start];
         r1 = A_i[row];
         r2 = A_i[row + 1];
         for (j = r1; j < r2; j++)
         {
            col = rqperm[A_j[j]];
            if (col != row && col >= start && col < end)
            {
               /* this is an entry in G */
               G_j[G_nnz++] = col - start;
               if (G_nnz >= G_capacity)
               {
                  JX_Int tmp = G_capacity;
                  G_capacity = G_capacity * EXPAND_FACT + 1;
                  G_j = jx_TReAlloc_v2(G_j, JX_Int, tmp, JX_Int, G_capacity);
               }
            }
         }
      }
      G_i[num_nodes] = G_nnz;
      if (G_nnz == 0)
      {
         // G has only diagonal, no need to do any kind of RCM
         jx_TFree(G_j);
         jx_TFree(rqperm);
         *permp = perm;
         *qpermp = qperm;
         jx_CSRMatrixDestroy(G);
         return jx_error_flag;
      }
      jx_CSRMatrixJ(G) = G_j;
      G_data = jx_CTAlloc(JX_Real, G_nnz);
      jx_CSRMatrixData(G) = G_data;
      jx_CSRMatrixNumNonzeros(G) = G_nnz;

      /* now sum G with G' */
      jx_CSRMatrixTranspose(G, &GT, 1);
      GGT = jx_CSRMatrixAdd(G, GT);
      jx_CSRMatrixDestroy(G);
      jx_CSRMatrixDestroy(GT);
      G = GGT;
      GGT = NULL;
   }

   /* 3: Build Graph
    * Build RCM
    */
   /* no need to be shared, but perm should be shared */
   G_perm = jx_TAlloc(JX_Int, num_nodes);
   jx_ILULocalRCMOrder(G, G_perm);

   /* 4: Post processing
    * Free, set value, return
    */

   /* update to new index */
   perm_temp = jx_TAlloc(JX_Int, num_nodes);
   for (i = 0; i < num_nodes; i++)
   {
      perm_temp[i] = perm[i + start];
   }
   for (i = 0; i < num_nodes; i++)
   {
      perm[i + start] = perm_temp[G_perm[i]];
   }
   if (perm != qperm)
   {
      for (i = 0; i < num_nodes; i++)
      {
         perm_temp[i] = qperm[i + start];
      }
      for (i = 0; i < num_nodes; i++)
      {
         qperm[i + start] = perm_temp[G_perm[i]];
      }
   }
   *permp = perm;
   *qpermp = qperm;
   jx_CSRMatrixDestroy(G);

   jx_TFree(G_perm);
   jx_TFree(perm_temp);
   jx_TFree(rqperm);

   return jx_error_flag;
}

/* This function finds the unvisited node with the minimum degree
 */
JX_Int
jx_ILULocalRCMMindegree(JX_Int n, JX_Int *degree, JX_Int *marker, JX_Int *rootp)
{
   JX_Int i;
   JX_Int min_degree = n + 1;
   JX_Int root = 0;
   for (i = 0; i < n; i++)
   {
      if (marker[i] < 0)
      {
         if (degree[i] < min_degree)
         {
            root = i;
            min_degree = degree[i];
         }
      }
   }
   *rootp = root;
   return 0;
}

/* This function actually does the RCM ordering of a symmetric csr matrix (entire)
 * A: the csr matrix, A_data is not needed
 * perm: the permutation array, space should be allocated outside
 */
JX_Int
jx_ILULocalRCMOrder(jx_CSRMatrix *A, JX_Int *perm)
{
   JX_Int i, root;
   JX_Int *degree = NULL;
   JX_Int *marker = NULL;
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int n = jx_CSRMatrixNumRows(A);
   JX_Int current_num;
   /* get the degree for each node */
   degree = jx_TAlloc(JX_Int, n);
   marker = jx_TAlloc(JX_Int, n);
   for (i = 0; i < n; i++)
   {
      degree[i] = A_i[i + 1] - A_i[i];
      marker[i] = -1;
   }

   /* start RCM loop */
   current_num = 0;
   while (current_num < n)
   {
      jx_ILULocalRCMMindegree(n, degree, marker, &root);
      /* This is a new connect component */
      jx_ILULocalRCMFindPPNode(A, &root, marker);

      /* Numbering of this component */
      jx_ILULocalRCMNumbering(A, root, marker, perm, &current_num);
   }

   /* free */
   jx_TFree(degree);
   jx_TFree(marker);
   return jx_error_flag;
}

/* This function find a pseudo-peripheral node start from root
 * A: the csr matrix, A_data is not needed
 * rootp: pointer to the root, on return will be a end of the pseudo-peripheral
 * marker: the marker array for unvisited node
 */
JX_Int
jx_ILULocalRCMFindPPNode(jx_CSRMatrix *A, JX_Int *rootp, JX_Int *marker)
{
   JX_Int i, r1, r2, row, min_degree, lev_degree, nlev, newnlev;

   JX_Int root = *rootp;
   JX_Int n = jx_CSRMatrixNumRows(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   /* at most n levels */
   JX_Int *level_i = jx_TAlloc(JX_Int, n + 1);
   JX_Int *level_j = jx_TAlloc(JX_Int, n);

   /* build initial level structure from root */
   jx_ILULocalRCMBuildLevel(A, root, marker, level_i, level_j, &newnlev);

   nlev = newnlev - 1;
   while (nlev < newnlev)
   {
      nlev = newnlev;
      r1 = level_i[nlev - 1];
      r2 = level_i[nlev];
      min_degree = n;
      for (i = r1; i < r2; i++)
      {
         /* select the last level, pick min-degree node */
         row = level_j[i];
         lev_degree = A_i[row + 1] - A_i[row];
         if (min_degree > lev_degree)
         {
            min_degree = lev_degree;
            root = row;
         }
      }
      jx_ILULocalRCMBuildLevel(A, root, marker, level_i, level_j, &newnlev);
   }

   *rootp = root;
   /* free */
   jx_TFree(level_i);
   jx_TFree(level_j);
   return jx_error_flag;
}

/* This function build level structure start from root
 * A: the csr matrix, A_data is not needed
 * root: pointer to the root
 * marker: the marker array for unvisited node
 * level_i: points to the start/end of position on level_j, similar to CSR Matrix
 * level_j: store node number on each level
 * nlevp: return the number of level on this level structure
 */
JX_Int
jx_ILULocalRCMBuildLevel(jx_CSRMatrix *A, JX_Int root, JX_Int *marker,
                         JX_Int *level_i, JX_Int *level_j, JX_Int *nlevp)
{
   JX_Int i, j, l1, l2, l_current, r1, r2, rowi, rowj, nlev;
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);

   /* set first level first */
   level_i[0] = 0;
   level_j[0] = root;
   marker[root] = 0;
   nlev = 1;
   l1 = 0;
   l2 = 1;
   l_current = l2;

   // explore nbhds of all nodes in current level
   while (l2 > l1)
   {
      level_i[nlev++] = l2;
      /* loop through last level */
      for (i = l1; i < l2; i++)
      {
         /* the node to explore */
         rowi = level_j[i];
         r1 = A_i[rowi];
         r2 = A_i[rowi + 1];
         for (j = r1; j < r2; j++)
         {
            rowj = A_j[j];
            if (marker[rowj] < 0)
            {
               /* Aha, an unmarked row */
               marker[rowj] = 0;
               level_j[l_current++] = rowj;
            }
         }
      }
      l1 = l2;
      l2 = l_current;
   }
   /* after this we always have a "ghost" last level */
   nlev--;

   /* reset marker */
   for (i = 0; i < l2; i++)
   {
      marker[level_j[i]] = -1;
   }

   *nlevp = nlev;

   return jx_error_flag;
}

/* This function generate numbering for a connect component
 * A: the csr matrix, A_data is not needed
 * root: pointer to the root
 * marker: the marker array for unvisited node
 * perm: permutation array
 * current_nump: number of nodes already have a perm value
 */

JX_Int
jx_ILULocalRCMNumbering(jx_CSRMatrix *A, JX_Int root, JX_Int *marker, JX_Int *perm, JX_Int *current_nump)
{
   JX_Int i, j, l1, l2, r1, r2, rowi, rowj, row_start, row_end;
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int current_num = *current_nump;

   marker[root] = 0;
   l1 = current_num;
   perm[current_num++] = root;
   l2 = current_num;

   while (l2 > l1)
   {
      /* loop through all nodes is current level */
      for (i = l1; i < l2; i++)
      {
         rowi = perm[i];
         r1 = A_i[rowi];
         r2 = A_i[rowi + 1];
         row_start = current_num;
         for (j = r1; j < r2; j++)
         {
            rowj = A_j[j];
            if (marker[rowj] < 0)
            {
               /* save the degree in marker and add it to perm */
               marker[rowj] = A_i[rowj + 1] - A_i[rowj];
               perm[current_num++] = rowj;
            }
         }
         row_end = current_num;
         jx_ILULocalRCMQsort(perm, row_start, row_end - 1, marker);
      }
      l1 = l2;
      l2 = current_num;
   }

   // reverse
   jx_ILULocalRCMReverse(perm, *current_nump, current_num - 1);
   *current_nump = current_num;
   return jx_error_flag;
}

/* This qsort is very specialized, not worth to put into utilities
 * Sort a part of array perm based on degree value (ascend)
 * That is, if degree[perm[i]] < degree[perm[j]], we should have i < j
 * perm: the perm array
 * start: start in perm
 * end: end in perm
 * degree: degree array
 */

JX_Int
jx_ILULocalRCMQsort(JX_Int *perm, JX_Int start, JX_Int end, JX_Int *degree)
{

   JX_Int i, mid;
   if (start >= end)
   {
      return jx_error_flag;
   }

   jx_swap(perm, start, (start + end) / 2);
   mid = start;
   // loop to split
   for (i = start + 1; i <= end; i++)
   {
      if (degree[perm[i]] < degree[perm[start]])
      {
         jx_swap(perm, ++mid, i);
      }
   }
   jx_swap(perm, start, mid);
   jx_ILULocalRCMQsort(perm, mid + 1, end, degree);
   jx_ILULocalRCMQsort(perm, start, mid - 1, degree);
   return jx_error_flag;
}

/* Last step in RCM, reverse it
 * perm: perm array
 * srart: start position
 * end: end position
 */

JX_Int
jx_ILULocalRCMReverse(JX_Int *perm, JX_Int start, JX_Int end)
{
   JX_Int i, j;
   JX_Int mid = (start + end + 1) / 2;

   for (i = start, j = end; i < mid; i++, j--)
   {
      jx_swap(perm, i, j);
   }
   return jx_error_flag;
}

/* NSH create and solve and help functions */

/* Create */
void *
jx_NSHCreate()
{
   jx_ParNSHData *nsh_data;

   nsh_data = jx_CTAlloc(jx_ParNSHData, 1);

   /* general data */
   jx_ParNSHDataMatA(nsh_data) = NULL;
   jx_ParNSHDataMatM(nsh_data) = NULL;
   jx_ParNSHDataF(nsh_data) = NULL;
   jx_ParNSHDataU(nsh_data) = NULL;
   jx_ParNSHDataResidual(nsh_data) = NULL;
   jx_ParNSHDataRelResNorms(nsh_data) = NULL;
   jx_ParNSHDataNumIterations(nsh_data) = 0;
   jx_ParNSHDataL1Norms(nsh_data) = NULL;
   jx_ParNSHDataFinalRelResidualNorm(nsh_data) = 0.0;
   jx_ParNSHDataTol(nsh_data) = 1e-09;
   jx_ParNSHDataLogging(nsh_data) = 2;
   jx_ParNSHDataPrintLevel(nsh_data) = 2;
   jx_ParNSHDataMaxIter(nsh_data) = 5;

   jx_ParNSHDataOperatorComplexity(nsh_data) = 0.0;
   jx_ParNSHDataDroptol(nsh_data) = jx_TAlloc(JX_Real, 2);
   jx_ParNSHDataOwnDroptolData(nsh_data) = 1;
   jx_ParNSHDataDroptol(nsh_data)[0] = 1.0e-02; /* droptol for MR */
   jx_ParNSHDataDroptol(nsh_data)[1] = 1.0e-02; /* droptol for NSH */
   jx_ParNSHDataUTemp(nsh_data) = NULL;
   jx_ParNSHDataFTemp(nsh_data) = NULL;

   /* MR data */
   jx_ParNSHDataMRMaxIter(nsh_data) = 2;
   jx_ParNSHDataMRTol(nsh_data) = 1e-09;
   jx_ParNSHDataMRMaxRowNnz(nsh_data) = 800;
   jx_ParNSHDataMRColVersion(nsh_data) = 0;

   /* NSH data */
   jx_ParNSHDataNSHMaxIter(nsh_data) = 2;
   jx_ParNSHDataNSHTol(nsh_data) = 1e-09;
   jx_ParNSHDataNSHMaxRowNnz(nsh_data) = 1000;

   return (void *)nsh_data;
}

/* Destroy */
JX_Int
jx_NSHDestroy(void *data)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)data;

   /* residual */
   if (jx_ParNSHDataResidual(nsh_data))
   {
      jx_ParVectorDestroy(jx_ParNSHDataResidual(nsh_data));
      jx_ParNSHDataResidual(nsh_data) = NULL;
   }

   /* residual norms */
   if (jx_ParNSHDataRelResNorms(nsh_data))
   {
      jx_TFree(jx_ParNSHDataRelResNorms(nsh_data));
      jx_ParNSHDataRelResNorms(nsh_data) = NULL;
   }

   /* l1 norms */
   if (jx_ParNSHDataL1Norms(nsh_data))
   {
      jx_TFree(jx_ParNSHDataL1Norms(nsh_data));
      jx_ParNSHDataL1Norms(nsh_data) = NULL;
   }

   /* temp arrays */
   if (jx_ParNSHDataUTemp(nsh_data))
   {
      jx_ParVectorDestroy(jx_ParNSHDataUTemp(nsh_data));
      jx_ParNSHDataUTemp(nsh_data) = NULL;
   }
   if (jx_ParNSHDataFTemp(nsh_data))
   {
      jx_ParVectorDestroy(jx_ParNSHDataFTemp(nsh_data));
      jx_ParNSHDataFTemp(nsh_data) = NULL;
   }

   /* approx inverse matrix */
   if (jx_ParNSHDataMatM(nsh_data))
   {
      jx_ParCSRMatrixDestroy(jx_ParNSHDataMatM(nsh_data));
      jx_ParNSHDataMatM(nsh_data) = NULL;
   }

   /* droptol array */
   if (jx_ParNSHDataOwnDroptolData(nsh_data))
   {
      jx_TFree(jx_ParNSHDataDroptol(nsh_data));
      jx_ParNSHDataOwnDroptolData(nsh_data) = 0;
      jx_ParNSHDataDroptol(nsh_data) = NULL;
   }

   /* nsh data */
   jx_TFree(nsh_data);

   return jx_error_flag;
}

/* Print solver params */
JX_Int
jx_NSHWriteSolverParams(void *nsh_vdata)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_printf("Newton–Schulz–Hotelling Setup parameters: \n");
   jx_printf("NSH max iterations = %d \n", jx_ParNSHDataNSHMaxIter(nsh_data));
   jx_printf("NSH drop tolerance = %e \n", jx_ParNSHDataDroptol(nsh_data)[1]);
   jx_printf("NSH max nnz per row = %d \n", jx_ParNSHDataNSHMaxRowNnz(nsh_data));
   jx_printf("MR max iterations = %d \n", jx_ParNSHDataMRMaxIter(nsh_data));
   jx_printf("MR drop tolerance = %e \n", jx_ParNSHDataDroptol(nsh_data)[0]);
   jx_printf("MR max nnz per row = %d \n", jx_ParNSHDataMRMaxRowNnz(nsh_data));
   jx_printf("Operator Complexity (Fill factor) = %f \n", jx_ParNSHDataOperatorComplexity(nsh_data));
   jx_printf("\n Newton–Schulz–Hotelling Solver Parameters: \n");
   jx_printf("Max number of iterations: %d\n", jx_ParNSHDataMaxIter(nsh_data));
   jx_printf("Stopping tolerance: %e\n", jx_ParNSHDataTol(nsh_data));

   return jx_error_flag;
}

/* set print level */
JX_Int
jx_NSHSetPrintLevel(void *nsh_vdata, JX_Int print_level)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataPrintLevel(nsh_data) = print_level;
   return jx_error_flag;
}
/* set logging level */
JX_Int
jx_NSHSetLogging(void *nsh_vdata, JX_Int logging)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataLogging(nsh_data) = logging;
   return jx_error_flag;
}
/* set max iteration */
JX_Int
jx_NSHSetMaxIter(void *nsh_vdata, JX_Int max_iter)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataMaxIter(nsh_data) = max_iter;
   return jx_error_flag;
}
/* set solver iteration tol */
JX_Int
jx_NSHSetTol(void *nsh_vdata, JX_Real tol)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataTol(nsh_data) = tol;
   return jx_error_flag;
}
/* set global solver */
JX_Int
jx_NSHSetGlobalSolver(void *nsh_vdata, JX_Int global_solver)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataGlobalSolver(nsh_data) = global_solver;
   return jx_error_flag;
}
/* set all droptols */
JX_Int
jx_NSHSetDropThreshold(void *nsh_vdata, JX_Real droptol)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataDroptol(nsh_data)[0] = droptol;
   jx_ParNSHDataDroptol(nsh_data)[1] = droptol;
   return jx_error_flag;
}
/* set array of droptols */
JX_Int
jx_NSHSetDropThresholdArray(void *nsh_vdata, JX_Real *droptol)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   if (jx_ParNSHDataOwnDroptolData(nsh_data))
   {
      jx_TFree(jx_ParNSHDataDroptol(nsh_data));
      jx_ParNSHDataOwnDroptolData(nsh_data) = 0;
   }
   jx_ParNSHDataDroptol(nsh_data) = droptol;
   return jx_error_flag;
}
/* set own data */
JX_Int
jx_NSHSetOwnDroptolData(void *nsh_vdata, JX_Int own_droptol_data)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataOwnDroptolData(nsh_data) = own_droptol_data;
   return jx_error_flag;
}
/* set MR max iter */
JX_Int
jx_NSHSetMRMaxIter(void *nsh_vdata, JX_Int mr_max_iter)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataMRMaxIter(nsh_data) = mr_max_iter;
   return jx_error_flag;
}
/* set MR tol */
JX_Int
jx_NSHSetMRTol(void *nsh_vdata, JX_Real mr_tol)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataMRTol(nsh_data) = mr_tol;
   return jx_error_flag;
}
/* set MR max nonzeros of a row */
JX_Int
jx_NSHSetMRMaxRowNnz(void *nsh_vdata, JX_Int mr_max_row_nnz)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataMRMaxRowNnz(nsh_data) = mr_max_row_nnz;
   return jx_error_flag;
}
/* set MR version, column version or global version */
JX_Int
jx_NSHSetColVersion(void *nsh_vdata, JX_Int mr_col_version)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataMRColVersion(nsh_data) = mr_col_version;
   return jx_error_flag;
}
/* set NSH max iter */
JX_Int
jx_NSHSetNSHMaxIter(void *nsh_vdata, JX_Int nsh_max_iter)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataNSHMaxIter(nsh_data) = nsh_max_iter;
   return jx_error_flag;
}
/* set NSH tol */
JX_Int
jx_NSHSetNSHTol(void *nsh_vdata, JX_Real nsh_tol)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataNSHTol(nsh_data) = nsh_tol;
   return jx_error_flag;
}
/* set NSH max nonzeros of a row */
JX_Int
jx_NSHSetNSHMaxRowNnz(void *nsh_vdata, JX_Int nsh_max_row_nnz)
{
   jx_ParNSHData *nsh_data = (jx_ParNSHData *)nsh_vdata;
   jx_ParNSHDataNSHMaxRowNnz(nsh_data) = nsh_max_row_nnz;
   return jx_error_flag;
}

/* Compute the F norm of CSR matrix
 * A: the target CSR matrix
 * norm_io: output
 */
JX_Int
jx_CSRMatrixNormFro(jx_CSRMatrix *A, JX_Real *norm_io)
{
   JX_Real norm = 0.0;
   JX_Real *data = jx_CSRMatrixData(A);
   JX_Int i, k;
   k = jx_CSRMatrixNumNonzeros(A);
   /* main loop */
   for (i = 0; i < k; i++)
   {
      norm += data[i] * data[i];
   }
   *norm_io = sqrt(norm);
   return jx_error_flag;
}

/* Compute the norm of I-A where I is identity matrix and A is a CSR matrix
 * A: the target CSR matrix
 * norm_io: the output
 */
JX_Int
jx_CSRMatrixResNormFro(jx_CSRMatrix *A, JX_Real *norm_io)
{
   JX_Real norm = 0.0, value;
   JX_Int i, j, k1, k2, n;
   JX_Int *idx = jx_CSRMatrixI(A);
   JX_Int *cols = jx_CSRMatrixJ(A);
   JX_Real *data = jx_CSRMatrixData(A);

   n = jx_CSRMatrixNumRows(A);
   /* main loop to sum up data */
   for (i = 0; i < n; i++)
   {
      k1 = idx[i];
      k2 = idx[i + 1];
      /* check if we have diagonal in A */
      if (k2 > k1)
      {
         if (cols[k1] == i)
         {
            /* reduce 1 on diagonal */
            value = data[k1] - 1.0;
            norm += value * value;
         }
         else
         {
            /* we don't have diagonal in A, so we need to add 1 to norm */
            norm += 1.0;
            norm += data[k1] * data[k1];
         }
      }
      else
      {
         /* we don't have diagonal in A, so we need to add 1 to norm */
         norm += 1.0;
      }
      /* and the rest of the code */
      for (j = k1 + 1; j < k2; j++)
      {
         norm += data[j] * data[j];
      }
   }
   *norm_io = sqrt(norm);
   return jx_error_flag;
}

/* Compute the F norm of ParCSR matrix
 * A: the target CSR matrix
 */
JX_Int
jx_ParCSRMatrixNormFro(jx_ParCSRMatrix *A, JX_Real *norm_io)
{
   JX_Real local_norm = 0.0;
   JX_Real global_norm;
   MPI_Comm comm = jx_ParCSRMatrixComm(A);

   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);

   jx_CSRMatrixNormFro(A_diag, &local_norm);
   /* use global_norm to store offd for now */
   jx_CSRMatrixNormFro(A_offd, &global_norm);

   /* square and sum them */
   local_norm *= local_norm;
   local_norm += global_norm * global_norm;

   /* do communication to get global total sum */
   jx_MPI_Allreduce(&local_norm, &global_norm, 1, JX_MPI_REAL, MPI_SUM, comm);

   *norm_io = sqrt(global_norm);
   return jx_error_flag;
}

/* Compute the F norm of ParCSR matrix
 * Norm of I-A
 * A: the target CSR matrix
 */
JX_Int
jx_ParCSRMatrixResNormFro(jx_ParCSRMatrix *A, JX_Real *norm_io)
{
   JX_Real local_norm = 0.0;
   JX_Real global_norm;
   MPI_Comm comm = jx_ParCSRMatrixComm(A);

   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);

   /* compute I-A for diagonal */
   jx_CSRMatrixResNormFro(A_diag, &local_norm);
   /* use global_norm to store offd for now */
   jx_CSRMatrixNormFro(A_offd, &global_norm);

   /* square and sum them */
   local_norm *= local_norm;
   local_norm += global_norm * global_norm;

   /* do communication to get global total sum */
   jx_MPI_Allreduce(&local_norm, &global_norm, 1, JX_MPI_REAL, MPI_SUM, comm);

   *norm_io = sqrt(global_norm);
   return jx_error_flag;
}

/* Compute the trace of CSR matrix
 * A: the target CSR matrix
 * trace_io: the output trace
 */
JX_Int
jx_CSRMatrixTrace(jx_CSRMatrix *A, JX_Real *trace_io)
{
   JX_Real trace = 0.0;
   JX_Int *idx = jx_CSRMatrixI(A);
   JX_Int *cols = jx_CSRMatrixJ(A);
   JX_Real *data = jx_CSRMatrixData(A);
   JX_Int i, k1, k2, n;
   n = jx_CSRMatrixNumRows(A);
   for (i = 0; i < n; i++)
   {
      k1 = idx[i];
      k2 = idx[i + 1];
      if (cols[k1] == i && k2 > k1)
      {
         /* only add when diagonal is nonzero */
         trace += data[k1];
      }
   }

   *trace_io = trace;
   return jx_error_flag;
}

/* Scale CSR matrix A = scalar * A
 * A: the target CSR matrix
 * scalar: real number
 */
JX_Int
jx_CSRMatrixScaleH(jx_CSRMatrix *A, JX_Real scalar)
{
   JX_Real *data = jx_CSRMatrixData(A);
   JX_Int i, k;
   k = jx_CSRMatrixNumNonzeros(A);
   for (i = 0; i < k; i++)
   {
      data[i] *= scalar;
   }
   return jx_error_flag;
}

/* Scale ParCSR matrix A = scalar * A
 * A: the target CSR matrix
 * scalar: real number
 */
JX_Int
jx_ParCSRMatrixScaleH(jx_ParCSRMatrix *A, JX_Real scalar)
{
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   /* each thread scale local diag and offd */
   jx_CSRMatrixScaleH(A_diag, scalar);
   jx_CSRMatrixScaleH(A_offd, scalar);
   return jx_error_flag;
}

/*--------------------------------------------------------------------------
 * jx_CSRMatrixDropInplace
 *
 * Apply dropping to CSR matrix
 * A: the target CSR matrix
 * droptol: all entries have smaller absolute value than this will be dropped
 * max_row_nnz: max nonzeros allowed for each row, only largest max_row_nnz kept
 * we NEVER drop diagonal entry if exists
 *
 * TODO (VPM): Move this function to seq_mv
 *--------------------------------------------------------------------------*/

JX_Int
jx_CSRMatrixDropInplace(jx_CSRMatrix *A, JX_Real droptol, JX_Int max_row_nnz)
{
   jx_printf("[JX-ILU] WARNING: jx_CSRMatrixDropInplace is called. \n");
}

/*--------------------------------------------------------------------------
 * jx_ILUCSRMatrixInverseSelfPrecondMRGlobal
 *
 * Compute the inverse with MR of original CSR matrix
 * Global(not by each column) and out place version
 * A: the input matrix
 * M: the output matrix
 * droptol: the dropping tolorance
 * tol: when to stop the iteration
 * eps_tol: to avoid divide by 0
 * max_row_nnz: max number of nonzeros per row
 * max_iter: max number of iterations
 * print_level: the print level of this algorithm
 *--------------------------------------------------------------------------*/

JX_Int
jx_ILUCSRMatrixInverseSelfPrecondMRGlobal(jx_CSRMatrix *matA,
                                          jx_CSRMatrix **M,
                                          JX_Real droptol,
                                          JX_Real tol,
                                          JX_Real eps_tol,
                                          JX_Int max_row_nnz,
                                          JX_Int max_iter,
                                          JX_Int print_level)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUCSRMatrixInverseSelfPrecondMRGlobal is called. \n");
}

/*--------------------------------------------------------------------------
 * jx_ILUParCSRInverseNSH
 *
 * Compute inverse with NSH method
 * Use MR to get local initial guess
 * A: input matrix
 * M: output matrix
 * droptol: droptol array. droptol[0] for MR and droptol[1] for NSH.
 * mr_tol: tol for stop iteration for MR
 * nsh_tol: tol for stop iteration for NSH
 * esp_tol: tol for avoid divide by 0
 * mr_max_row_nnz: max number of nonzeros for MR
 * nsh_max_row_nnz: max number of nonzeros for NSH
 * mr_max_iter: max number of iterations for MR
 * nsh_max_iter: max number of iterations for NSH
 * mr_col_version: column version of global version
 *--------------------------------------------------------------------------*/
JX_Int
jx_ILUParCSRInverseNSH(jx_ParCSRMatrix *A,
                       jx_ParCSRMatrix **M,
                       JX_Real *droptol,
                       JX_Real mr_tol,
                       JX_Real nsh_tol,
                       JX_Real eps_tol,
                       JX_Int mr_max_row_nnz,
                       JX_Int nsh_max_row_nnz,
                       JX_Int mr_max_iter,
                       JX_Int nsh_max_iter,
                       JX_Int mr_col_version,
                       JX_Int print_level)
{
   jx_printf("[JX-ILU] WARNING: jx_ILUParCSRInverseNSH is called. \n");
}

// 张娜新增 -----------------------------------------------------------------
/*--------------------------------------------------------------------------
 * jx_ILUSetIRIters
 *
 * Set IR-Iteration times for ILU solver
 *--------------------------------------------------------------------------*/

// JX_Int
// jx_ILUSetIRIters(void *ilu_vdata,
//                  JX_Int IR_iters)
// {
//    jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

//    jx_ParILUDataIRIters(ilu_data) = IR_iters;

//    return jx_error_flag;
// }

/*--------------------------------------------------------------------------
 * jx_ILUSetsweep
 *
 * Set sweep times for ILU solver
 *--------------------------------------------------------------------------*/

// JX_Int
// jx_ILUSetsweep(void *ilu_vdata,
//                JX_Real sweep)
// {
//    jx_ParILUData *ilu_data = (jx_ParILUData *)ilu_vdata;

//    jx_ParILUDatasweep(ilu_data) = sweep;

//    return jx_error_flag;
// }
