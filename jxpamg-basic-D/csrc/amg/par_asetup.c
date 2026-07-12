//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_asetup.c
 *  Date: 2020/03/10
 */ 

#include "jx_asetup.h"

extern JX_Int *Pmarkers_global_rap;    /* Yue Xiaoqiang 2012/10/17 */
extern JX_Int *Pmarkers_global_interp; /* Yue Xiaoqiang 2012/10/17 */

/*!
 * \fn JX_Int JX_ParAdpSetupAMGCreate
 * \date 2020/03/10
 */
JX_Int
JX_ParAdpSetupAMGCreate( MPI_Comm comm, JX_Solver *solver )
{
  *solver = (JX_Solver) jx_ParAdpSetupAMGCreate( comm );
   if (!solver)
   {
      jx_error_in_arg(1);
   }
   return jx_error_flag;
}

/*!
 * \fn JX_Int JX_ParAdpSetupAMGDestroy
 * \date 2020/03/10
 */
JX_Int
JX_ParAdpSetupAMGDestroy( JX_Solver solver )
{
   return( jx_ParAdpSetupAMGDestroy( (void *) solver ) );
}

/*!
 * \fn JX_Int JX_ParAdpSetupAMGSetup
 * \date 2020/03/10
 */
JX_Int
JX_ParAdpSetupAMGSetup( JX_Solver solver, JX_Matrix matA, JX_Vector vecB, JX_Vector vecX )
{
   return( jx_ParAdpSetupAMGSetup( solver, matA, vecB, vecX ) );
}

/*!
 * \fn JX_Int JX_ParAdpSetupAMGSolve
 * \date 2020/03/10
 */
JX_Int
JX_ParAdpSetupAMGSolve( JX_Solver solver, JX_Matrix preA, JX_Matrix matA, JX_Vector vecB, JX_Vector vecX )
{
   return( jx_ParAdpSetupAMGSolve( (void *) solver, (void *) preA, (void *) matA, (void *) vecB, (void *) vecX ) );
}

JX_Int
JX_ParAdpSetupAMGSetIsAdaptive( JX_Solver solver, JX_Bool is_adaptive )
{
   return( jx_ParAdpSetupAMGSetIsAdaptive( (void *) solver, is_adaptive ) );
}

JX_Int
JX_ParAdpSetupAMGSetIsAdapSetup( JX_Solver solver, JX_Bool is_adap_setup )
{
   return( jx_ParAdpSetupAMGSetIsAdapSetup( (void *) solver, is_adap_setup ) );
}

JX_Int
JX_ParAdpSetupAMGSetTol( JX_Solver solver, JX_Real tol )
{
   return( jx_ParAdpSetupAMGSetTol( (void *) solver, tol ) );
}

JX_Int
JX_ParAdpSetupAMGSetMaxRowSum( JX_Solver solver, JX_Real max_row_sum )
{
   return( jx_ParAdpSetupAMGSetMaxRowSum( (void *) solver, max_row_sum ) );
}

JX_Int
JX_ParAdpSetupAMGSetStrongThreshold( JX_Solver solver, JX_Real strong_threshold )
{
   return( jx_ParAdpSetupAMGSetStrongThreshold( (void *) solver, strong_threshold ) );
}

JX_Int
JX_ParAdpSetupAMGSetMaxLevels( JX_Solver solver, JX_Int  max_levels )
{
   return( jx_ParAdpSetupAMGSetMaxLevels( (void *) solver, max_levels ) );
}

JX_Int
JX_ParAdpSetupAMGSetCycleType( JX_Solver solver,JX_Int cycle_type )
{
   return( jx_ParAdpSetupAMGSetCycleType( (void *) solver, cycle_type ) );
}

JX_Int
JX_ParAdpSetupAMGSetRelaxType( JX_Solver solver, JX_Int relax_type )
{
   return( jx_ParAdpSetupAMGSetRelaxType( (void *) solver, relax_type ) );
}

JX_Int
JX_ParAdpSetupAMGSetRelaxOrder( JX_Solver solver, JX_Int relax_order )
{
   return( jx_ParAdpSetupAMGSetRelaxOrder( (void *) solver, relax_order ) );
}

JX_Int
JX_ParAdpSetupAMGSetNsDown( JX_Solver solver, JX_Int ns_down )
{
   return( jx_ParAdpSetupAMGSetNsDown( (void *) solver, ns_down ) );
}

JX_Int
JX_ParAdpSetupAMGSetNsUp( JX_Solver solver, JX_Int ns_up )
{
   return( jx_ParAdpSetupAMGSetNsUp( (void *) solver, ns_up ) );
}

JX_Int
JX_ParAdpSetupAMGSetNsCoarse( JX_Solver solver, JX_Int ns_coarse )
{
   return( jx_ParAdpSetupAMGSetNsCoarse( (void *) solver, ns_coarse ) );
}

JX_Int
JX_ParAdpSetupAMGSetCoarsenType( JX_Solver solver, JX_Int coarsen_type )
{
   return( jx_ParAdpSetupAMGSetCoarsenType( (void *) solver, coarsen_type ) );
}

JX_Int
JX_ParAdpSetupAMGSetInterpType( JX_Solver solver, JX_Int interp_type )
{
   return( jx_ParAdpSetupAMGSetInterpType( (void *) solver, interp_type ) );
}

JX_Int
JX_ParAdpSetupAMGSetPMaxElmts( JX_Solver solver, JX_Int P_max_elmts )
{
   return( jx_ParAdpSetupAMGSetPMaxElmts( (void *) solver, P_max_elmts ) );
}

JX_Int
JX_ParAdpSetupAMGSetAggNumLevels( JX_Solver solver, JX_Int agg_num_levels )
{
   return( jx_ParAdpSetupAMGSetAggNumLevels( (void *) solver, agg_num_levels ) );
}

JX_Int
JX_ParAdpSetupAMGSetCoarseThreshold( JX_Solver solver, JX_Int coarse_threshold )
{
   return (jx_ParAdpSetupAMGSetCoarseThreshold ( (void *) solver, coarse_threshold ) );
}

JX_Int
JX_ParAdpSetupAMGSetAMGPrintLevel( JX_Solver solver, JX_Int amg_print_level )
{
   return( jx_ParAdpSetupAMGSetAMGPrintLevel( (void *) solver, amg_print_level ) );
}

JX_Int
JX_ParAdpSetupAMGSetKDim( JX_Solver solver, JX_Int k_dim )
{
   return( jx_ParAdpSetupAMGSetKDim( (void *) solver, k_dim ) );
}

JX_Int
JX_ParAdpSetupAMGSetTTest( JX_Solver solver, JX_Int TTest )
{
   return( jx_ParAdpSetupAMGSetTTest( (void *) solver, TTest ) );
}

JX_Int
JX_ParAdpSetupAMGSetMaxIter( JX_Solver solver, JX_Int max_iter )
{
   return( jx_ParAdpSetupAMGSetMaxIter( (void *) solver, max_iter ) );
}

JX_Int
JX_ParAdpSetupAMGSetTwoNorm( JX_Solver solver, JX_Int two_norm )
{
   return( jx_ParAdpSetupAMGSetTwoNorm( (void *) solver, two_norm ) );
}

JX_Int
JX_ParAdpSetupAMGSetSolverID( JX_Solver solver, JX_Int solver_id )
{
   return( jx_ParAdpSetupAMGSetSolverID( (void *) solver, solver_id ) );
}

JX_Int
JX_ParAdpSetupAMGSetPrintLevel( JX_Solver solver, JX_Int print_level )
{
   return( jx_ParAdpSetupAMGSetPrintLevel( (void *) solver, print_level ) );
}

JX_Int
JX_ParAdpSetupAMGSetMaxIterSL( JX_Solver solver, JX_Int max_iter_sl )
{
   return( jx_ParAdpSetupAMGSetMaxIterSL( (void *) solver, max_iter_sl ) );
}

JX_Int
JX_ParAdpSetupAMGSetMaxIterReUse( JX_Solver solver, JX_Int max_iter_reuse )
{
   return( jx_ParAdpSetupAMGSetMaxIterReUse( (void *) solver, max_iter_reuse ) );
}

JX_Int
JX_ParAdpSetupAMGSetIsCheckRestarted( JX_Solver solver, JX_Int is_check_restarted )
{
   return( jx_ParAdpSetupAMGSetIsCheckRestarted( (void *) solver, is_check_restarted ) );
}

/*!
 * \fn void *jx_ParAdpSetupAMGCreate
 * \date 2020/03/10
 */
void *
jx_ParAdpSetupAMGCreate( MPI_Comm comm )
{
   jx_ParAdpSetupAMGData *amg_data = jx_CTAlloc(jx_ParAdpSetupAMGData, 1);

   JX_Bool is_adaptive = JX_FALSE;
   JX_Bool is_adap_setup = JX_FALSE;

   JX_Real tol = 1.0e-7;
   JX_Real max_row_sum = 0.9;
   JX_Real strong_threshold = 0.25;

   JX_Int max_levels = 25;
   JX_Int cycle_type = 1;
   JX_Int relax_type = 3;
   JX_Int relax_order = 1;
   JX_Int ns_down = 1;
   JX_Int ns_up = 1;
   JX_Int ns_coarse = 1;
   JX_Int coarsen_type = 6;
   JX_Int interp_type = 0;
   JX_Int P_max_elmts = 0;
   JX_Int agg_num_levels = 0;
   JX_Int coarse_threshold = 9;
   JX_Int amg_print_level = 0;

   JX_Int k_dim = 5;
   JX_Int TTest = 0;
   JX_Int max_iter = 200;
   JX_Int two_norm = 0;
   JX_Int solver_id = 2;
   JX_Int print_level = 0;
   JX_Int max_iter_sl = 1000;
   JX_Int max_iter_reuse = 1000;
   JX_Int is_check_restarted = 0;

   jx_ParAdpSetupAMGDataComm(amg_data) = comm;
   jx_ParAdpSetupAMGDataIterSLLvls(amg_data) = jx_CTAlloc(JX_Int, 25);

   jx_ParAdpSetupAMGSetIsAdaptive(amg_data, is_adaptive);
   jx_ParAdpSetupAMGSetIsAdapSetup(amg_data, is_adap_setup);

   jx_ParAdpSetupAMGSetTol(amg_data, tol);
   jx_ParAdpSetupAMGSetMaxRowSum(amg_data, max_row_sum);
   jx_ParAdpSetupAMGSetStrongThreshold(amg_data, strong_threshold);

   jx_ParAdpSetupAMGSetMaxLevels(amg_data, max_levels);
   jx_ParAdpSetupAMGSetCycleType(amg_data, cycle_type);
   jx_ParAdpSetupAMGSetRelaxType(amg_data, relax_type);
   jx_ParAdpSetupAMGSetRelaxOrder(amg_data, relax_order);
   jx_ParAdpSetupAMGSetNsDown(amg_data, ns_down);
   jx_ParAdpSetupAMGSetNsUp(amg_data, ns_up);
   jx_ParAdpSetupAMGSetNsCoarse(amg_data, ns_coarse);
   jx_ParAdpSetupAMGSetCoarsenType(amg_data, coarsen_type);
   jx_ParAdpSetupAMGSetInterpType(amg_data, interp_type);
   jx_ParAdpSetupAMGSetPMaxElmts(amg_data, P_max_elmts);
   jx_ParAdpSetupAMGSetAggNumLevels(amg_data, agg_num_levels);
   jx_ParAdpSetupAMGSetCoarseThreshold(amg_data, coarse_threshold);
   jx_ParAdpSetupAMGSetAMGPrintLevel(amg_data, amg_print_level);

   jx_ParAdpSetupAMGSetKDim(amg_data, k_dim);
   jx_ParAdpSetupAMGSetTTest(amg_data, TTest);
   jx_ParAdpSetupAMGSetMaxIter(amg_data, max_iter);
   jx_ParAdpSetupAMGSetTwoNorm(amg_data, two_norm);
   jx_ParAdpSetupAMGSetSolverID(amg_data, solver_id);
   jx_ParAdpSetupAMGSetPrintLevel(amg_data, print_level);
   jx_ParAdpSetupAMGSetMaxIterSL(amg_data, max_iter_sl);
   jx_ParAdpSetupAMGSetMaxIterReUse(amg_data, max_iter_reuse);
   jx_ParAdpSetupAMGSetIsCheckRestarted(amg_data, is_check_restarted);

   jx_ParAdpSetupAMGDataNumSL(amg_data) = 0;
   jx_ParAdpSetupAMGDataCanReUse(amg_data) = 0;
   jx_ParAdpSetupAMGDataNumReUse(amg_data) = 0;
   jx_ParAdpSetupAMGDataNumSLAMG(amg_data) = 0;
   jx_ParAdpSetupAMGDataSetupTime(amg_data) = 0.0;
   jx_ParAdpSetupAMGDataSolveTime(amg_data) = 0.0;
   jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
   jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) = 0;
   jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) = 0;
   jx_ParAdpSetupAMGDataTotalConvergeIter(amg_data) = 0;
   jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) = 0;

   return (void *) amg_data;
}

/*!
 * \fn JX_Int jx_ParAdpSetupAMGDestroy
 * \date 2020/03/10
 */
JX_Int
jx_ParAdpSetupAMGDestroy( void *data )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (jx_ParAdpSetupAMGDataAMGSolver(amg_data))
   {
      JX_PAMGDestroy(jx_ParAdpSetupAMGDataAMGSolver(amg_data));
   }
   jx_TFree(jx_ParAdpSetupAMGDataIterSLLvls(amg_data));
   jx_TFree(amg_data);
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetup( void *data, void *matA, void *vecB, void *vecX )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   MPI_Comm comm = jx_ParAdpSetupAMGDataComm(amg_data);
   JX_Real tol = jx_ParAdpSetupAMGDataTol(amg_data);
   JX_Int k_dim = jx_ParAdpSetupAMGDataKDim(amg_data);
   JX_Int TTest = jx_ParAdpSetupAMGDataTTest(amg_data);
   JX_Int max_iter = jx_ParAdpSetupAMGDataMaxIter(amg_data);
   JX_Int two_norm = jx_ParAdpSetupAMGDataTwoNorm(amg_data);
   JX_Int solver_id = jx_ParAdpSetupAMGDataSolverID(amg_data);
   JX_Int print_level = jx_ParAdpSetupAMGDataPrintLevel(amg_data);
   JX_Int is_check_restarted = jx_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JX_Solver solver;
   JX_Solver jac_solver;
   JX_Real starttime, endtime;

   if (TTest) starttime = jx_MPI_Wtime();
   JX_JacobiCreate(&jac_solver);
   JX_JacobiSetMaxIter(jac_solver, 1);
   JX_JacobiSetup(jac_solver, (JX_hpCSRMatrix) matA);
   jx_ParAdpSetupAMGDataJACSolver(amg_data) = jac_solver;
   if (solver_id == 1)
   {
      JX_PCGCreate(comm, &solver);
      JX_PCGSetMaxIter(solver, max_iter);
      JX_PCGSetTol(solver, tol);
      JX_PCGSetTwoNorm(solver, two_norm);
      JX_PCGSetLogging(solver, 1);
      JX_PCGSetPrintLevel(solver, print_level);
      JX_PCGSetPrecond(solver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_PCGSetup(solver, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
   }
   else if (solver_id == 2)
   {
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level);
      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_GMRESSetup(solver, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
   }
   else if (solver_id == 3)
   {
      JX_BiCGSTABCreate(comm, &solver);
      JX_BiCGSTABSetMaxIter(solver, max_iter);
      JX_BiCGSTABSetTol(solver, tol);
      JX_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JX_BiCGSTABSetConvCriteria(solver, 0);
      JX_BiCGSTABSetLogging(solver, 1);
      JX_BiCGSTABSetPrintLevel(solver, print_level);
      JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_BiCGSTABSetup(solver, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
   }
   jx_ParAdpSetupAMGDataSolver(amg_data) = solver;
   if (TTest)
   {
      endtime = jx_MPI_Wtime();
      jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtime - starttime);
   }

   return 0;
}

JX_Int
jx_ParAdpSetupAMGSolve( void *data, void *preA, void *matA, void *vecB, void *vecX )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   JX_Bool is_adaptive = jx_ParAdpSetupAMGDataIsAdaptive(amg_data);
   JX_Bool is_adap_setup = jx_ParAdpSetupAMGDataIsAdapSetup(amg_data);
   JX_Real tol = jx_ParAdpSetupAMGDataTol(amg_data);
   JX_Real max_row_sum = jx_ParAdpSetupAMGDataMaxRowSum(amg_data);
   JX_Real strong_threshold = jx_ParAdpSetupAMGDataStrongThreshold(amg_data);
   JX_Int k_dim = jx_ParAdpSetupAMGDataKDim(amg_data);
   JX_Int max_iter = jx_ParAdpSetupAMGDataMaxIter(amg_data);
   JX_Int two_norm = jx_ParAdpSetupAMGDataTwoNorm(amg_data);
   JX_Int is_check_restarted = jx_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JX_Int TTest = jx_ParAdpSetupAMGDataTTest(amg_data);
   JX_Int solver_id = jx_ParAdpSetupAMGDataSolverID(amg_data);
   JX_Int max_iter_sl = jx_ParAdpSetupAMGDataMaxIterSL(amg_data);
   JX_Int max_iter_reuse = jx_ParAdpSetupAMGDataMaxIterReUse(amg_data);
   JX_Solver solver = jx_ParAdpSetupAMGDataSolver(amg_data);
   JX_Solver amg_solver = jx_ParAdpSetupAMGDataAMGSolver(amg_data);
   JX_Int max_levels = jx_ParAdpSetupAMGDataMaxLevels(amg_data);
   JX_Int cycle_type = jx_ParAdpSetupAMGDataCycleType(amg_data);
   JX_Int relax_type = jx_ParAdpSetupAMGDataRelaxType(amg_data);
   JX_Int relax_order = jx_ParAdpSetupAMGDataRelaxOrder(amg_data);
   JX_Int ns_down = jx_ParAdpSetupAMGDataNsDown(amg_data);
   JX_Int ns_up = jx_ParAdpSetupAMGDataNsUp(amg_data);
   JX_Int ns_coarse = jx_ParAdpSetupAMGDataNsCoarse(amg_data);
   JX_Int coarsen_type = jx_ParAdpSetupAMGDataCoarsenType(amg_data);
   JX_Int interp_type = jx_ParAdpSetupAMGDataInterpType(amg_data);
   JX_Int P_max_elmts = jx_ParAdpSetupAMGDataPMaxElmts(amg_data);
   JX_Int agg_num_levels = jx_ParAdpSetupAMGDataAggNumLevels(amg_data);
   JX_Int coarse_threshold = jx_ParAdpSetupAMGDataCoarseThreshold(amg_data);
   JX_Int amg_print_level = jx_ParAdpSetupAMGDataAMGPrintLevel(amg_data);
   JX_Int *iter_sl_lvls = jx_ParAdpSetupAMGDataIterSLLvls(amg_data);
   JX_Bool ok_ok = JX_TRUE;
   JX_Int iter_sl_vvv = 0;
   JX_Int iter_sl_used = 0;
   JX_Int iter_amg_used = 0;
   JX_Real iter_rest;
   JX_Real res_norm;
   JX_Real res_norm_old;
   JX_Real reduce_factor;
   JX_Real reduce_factor_thrd = 1.0;
   JX_Real starttime, endtime;
   JX_Real starttimeT = 0.0, endtimeT = 0.0;
   JX_Int iter_sl_oused;
   JX_Int num_levels;
   JX_Int num_sl_amg;
   JX_Int num_iterations;

   if (TTest) starttime = jx_MPI_Wtime();
   if (solver_id == 1)
   {
      JX_PCGSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_PCGSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
         JX_PCGGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jx_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JX_JacobiDestroy(jx_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JX_FALSE;
         JX_PCGDestroy(solver);
         jx_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jx_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JX_PCGSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) NULL, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JX_PCGSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
            JX_PCGGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JX_FALSE;
            jx_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JX_PCGDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_PCGSetMaxIter(solver, max_iter);
         JX_PCGSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) NULL, amg_solver);
         JX_PAMGAdapSetup(amg_solver, (JX_hpCSRMatrix) matA, (JX_ParVector) vecB, (JX_ParVector) vecX, 1, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jx_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_PCGSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_PCGGetNumIterations(solver, &num_iterations);
         JX_PCGDestroy(solver);
      }
      else if (ok_ok)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_PCGSetMaxIter(solver, max_iter);
         JX_PCGSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) JX_PAMGSetup, amg_solver);
         JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix) matA);
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_PCGSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_PCGGetNumIterations(solver, &num_iterations);
         JX_PCGDestroy(solver);
      }
   }
   else if (solver_id == 2)
   {
      JX_GMRESSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_GMRESSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
         JX_GMRESGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jx_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JX_JacobiDestroy(jx_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JX_FALSE;
         JX_GMRESDestroy(solver);
         jx_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jx_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) JX_PAMGSetup, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JX_GMRESSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
            JX_GMRESGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JX_FALSE;
            jx_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JX_GMRESDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_GMRESSetMaxIter(solver, max_iter);
         JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) NULL, amg_solver);
         JX_PAMGAdapSetup(amg_solver, (JX_hpCSRMatrix) matA, (JX_ParVector) vecB, (JX_ParVector) vecX, 2, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jx_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_GMRESSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_GMRESGetNumIterations(solver, &num_iterations);
         JX_GMRESDestroy(solver);
      }
      else if (ok_ok)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_GMRESSetMaxIter(solver, max_iter);
         JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) JX_PAMGSetup, amg_solver);
         JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix) matA);
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_GMRESSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_GMRESGetNumIterations(solver, &num_iterations);
         JX_GMRESDestroy(solver);
      }
   }
   else if (solver_id == 3)
   {
      JX_BiCGSTABSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_BiCGSTABSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
         JX_BiCGSTABGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jx_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JX_JacobiDestroy(jx_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JX_FALSE;
         JX_BiCGSTABDestroy(solver);
         jx_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jx_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) JX_PAMGSetup, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JX_BiCGSTABSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector)vecX);
            JX_BiCGSTABGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JX_FALSE;
            jx_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JX_BiCGSTABDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_BiCGSTABSetMaxIter(solver, max_iter);
         JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) NULL, amg_solver);
         JX_PAMGAdapSetup(amg_solver, (JX_hpCSRMatrix) matA, (JX_ParVector) vecB, (JX_ParVector) vecX, 3, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jx_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_BiCGSTABSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_BiCGSTABGetNumIterations(solver, &num_iterations);
         JX_BiCGSTABDestroy(solver);
      }
      else if (ok_ok)
      {
         jx_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jx_MPI_Wtime();
         if (jx_ParAdpSetupAMGDataCanReUse(amg_data)) jx_hpCSRMatrixDestroy((jx_hpCSRMatrix *)jx_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JX_PAMGDestroy(amg_solver);
         JX_PAMGCreate(&amg_solver);
         JX_PAMGSetMaxIter(amg_solver, 1);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JX_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JX_BiCGSTABSetMaxIter(solver, max_iter);
         JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn) JX_PAMGPrecond, (JX_PtrToSolverFcn) JX_PAMGSetup, amg_solver);
         JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix) matA);
         jx_ParAdpSetupAMGDataPreMat(amg_data) = (JX_hpCSRMatrix) matA;
         jx_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jx_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jx_MPI_Wtime();
            jx_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JX_BiCGSTABSolve(solver, (JX_Matrix) preA, (JX_Matrix) matA, (JX_Vector) vecB, (JX_Vector) vecX);
         JX_PAMGGetNumLevels(amg_solver, &num_levels);
         jx_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JX_BiCGSTABGetNumIterations(solver, &num_iterations);
         JX_BiCGSTABDestroy(solver);
      }
   }
   jx_ParAdpSetupAMGDataTotalConvergeIter(amg_data) += num_iterations;
   if (TTest)
   {
      endtime = jx_MPI_Wtime();
      jx_ParAdpSetupAMGDataSolveTime(amg_data) += ((endtime - starttime) - (endtimeT - starttimeT));
   }

   return 0;
}

JX_Int
jx_ParAdpSetupAMGSetIsAdaptive( void *data, JX_Bool is_adaptive )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataIsAdaptive(amg_data) = is_adaptive;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetIsAdapSetup( void *data, JX_Bool is_adap_setup )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataIsAdapSetup(amg_data) = is_adap_setup;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetTol( void *data, JX_Real tol )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (tol < 0 || tol > 1)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataTol(amg_data) = tol;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetMaxRowSum( void *data, JX_Real max_row_sum )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (max_row_sum <= 0 || max_row_sum > 1)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataMaxRowSum(amg_data) = max_row_sum;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetStrongThreshold( void *data, JX_Real strong_threshold )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (strong_threshold < 0 || strong_threshold > 1)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataStrongThreshold(amg_data) = strong_threshold;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetMaxLevels( void *data, JX_Int max_levels )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (max_levels < 1)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataMaxLevels(amg_data) = max_levels;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetCycleType( void *data, JX_Int cycle_type )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (cycle_type < 0 || cycle_type > 2)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataCycleType(amg_data) = cycle_type;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetRelaxType( void *data, JX_Int relax_type )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (relax_type < 0)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataRelaxType(amg_data) = relax_type;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetRelaxOrder( void *data, JX_Int relax_order )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataRelaxOrder(amg_data) = relax_order;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetNsDown( void *data, JX_Int ns_down )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataNsDown(amg_data) = ns_down;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetNsUp( void *data, JX_Int ns_up )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataNsUp(amg_data) = ns_up;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetNsCoarse( void *data, JX_Int ns_coarse )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataNsCoarse(amg_data) = ns_coarse;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetCoarsenType( void *data, JX_Int coarsen_type )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataCoarsenType(amg_data) = coarsen_type;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetInterpType( void *data, JX_Int interp_type )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (interp_type < 0 || interp_type > 100)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataInterpType(amg_data) = interp_type;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetPMaxElmts( void *data, JX_Int P_max_elmts )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (P_max_elmts < 0)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataPMaxElmts(amg_data) = P_max_elmts;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetAggNumLevels( void *data, JX_Int agg_num_levels )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (agg_num_levels < 0)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataAggNumLevels(amg_data) = agg_num_levels;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetCoarseThreshold( void *data, JX_Int coarse_threshold )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataCoarseThreshold(amg_data) = coarse_threshold;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetAMGPrintLevel( void *data, JX_Int amg_print_level )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataAMGPrintLevel(amg_data) = amg_print_level;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetKDim( void *data, JX_Int k_dim )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataKDim(amg_data) = k_dim;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetTTest( void *data, JX_Int TTest )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataTTest(amg_data) = TTest;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetMaxIter( void *data, JX_Int max_iter )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataMaxIter(amg_data) = max_iter;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetTwoNorm( void *data, JX_Int two_norm )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataTwoNorm(amg_data) = two_norm;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetSolverID( void *data, JX_Int solver_id )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataSolverID(amg_data) = solver_id;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetPrintLevel( void *data, JX_Int print_level )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataPrintLevel(amg_data) = print_level;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetMaxIterSL( void *data, JX_Int max_iter_sl )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataMaxIterSL(amg_data) = max_iter_sl;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetMaxIterReUse( void *data, JX_Int max_iter_reuse )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataMaxIterReUse(amg_data) = max_iter_reuse;
   return jx_error_flag;
}

JX_Int
jx_ParAdpSetupAMGSetIsCheckRestarted( void *data, JX_Int is_check_restarted )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   jx_ParAdpSetupAMGDataIsCheckRestarted(amg_data) = is_check_restarted;
   return jx_error_flag;
}

JX_Int
JX_JacobiCreate( JX_Solver *solver )
{
  *solver = (JX_Solver) jx_JacobiCreate();
   if (!solver)
   {
      jx_error_in_arg(1);
   }
   return jx_error_flag;
}

void *
jx_JacobiCreate()
{
   jx_JacobiData *amg_data = jx_CTAlloc(jx_JacobiData, 1);
   jx_JacobiDataMaxIter(amg_data) = 100;
   return (void *) amg_data;
}

JX_Int
JX_JacobiDestroy( JX_Solver solver )
{
   return( jx_JacobiDestroy( (void *) solver ) );
}

JX_Int
jx_JacobiDestroy( void *data )
{
   jx_JacobiData *amg_data = data;
   if (jx_JacobiDataTmpVec(amg_data))
   {
      jx_ParVectorDestroy(jx_JacobiDataTmpVec(amg_data));
   }
   jx_TFree(amg_data);
   return jx_error_flag;
}

JX_Int
JX_JacobiSetMaxIter( JX_Solver solver, JX_Int max_iter )
{
   return( jx_JacobiSetMaxIter( (void *) solver, max_iter ) );
}

JX_Int
jx_JacobiSetMaxIter( void *data, JX_Int max_iter )
{
   jx_JacobiData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! Jacobi object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (max_iter < 0)
   {
      jx_error_in_arg(2);
      return jx_error_flag;
   }
   jx_JacobiDataMaxIter(amg_data) = max_iter;
   return jx_error_flag;
}

JX_Int 
JX_JacobiSetup( JX_Solver solver, JX_hpCSRMatrix hp_matrix )
{
   return( jx_JacobiSetup( (void *) solver, jx_hpCSRMatrixPar((jx_hpCSRMatrix *) hp_matrix )) );
}

JX_Int
jx_JacobiSetup( void *data, jx_ParCSRMatrix *par_matrix )
{
   jx_JacobiData *amg_data = data;
   jx_JacobiDataTmpVec(amg_data) = jx_ParVectorCreate(jx_ParCSRMatrixComm(par_matrix),
                                                      jx_ParCSRMatrixGlobalNumRows(par_matrix),
                                                      jx_ParCSRMatrixRowStarts(par_matrix));
   jx_ParVectorInitialize(jx_JacobiDataTmpVec(amg_data));
   jx_ParVectorSetPartitioningOwner(jx_JacobiDataTmpVec(amg_data), 0);
   return 0;
}

JX_Int
JX_JacobiPrecond( JX_Solver solver, JX_hpCSRMatrix hp_matrix, JX_ParVector par_rhs, JX_ParVector par_app )
{
   return( jx_JacobiPrecond( (void *) solver, jx_hpCSRMatrixPar((jx_hpCSRMatrix *) hp_matrix), (jx_ParVector *) par_rhs, (jx_ParVector *) par_app ) );
}

JX_Int
jx_JacobiPrecond( void *data, jx_ParCSRMatrix *par_matrix, jx_ParVector *par_rhs, jx_ParVector *par_app )
{
   jx_JacobiData *amg_data = data;
   jx_ParVector *tmp_vec = jx_JacobiDataTmpVec(amg_data);
   JX_Int max_iter = jx_JacobiDataMaxIter(amg_data);
   JX_Int cycle_count = 0;
   while (cycle_count < max_iter)
   {
      jx_PAMGRelax0(par_matrix, par_rhs, NULL, 0, 1.0, 0.0, par_app, tmp_vec);
      cycle_count ++;
   }
   return 0;
}

JX_Int
JX_PAMGAdapSetup( JX_Solver solver, JX_hpCSRMatrix matA, JX_ParVector vecB, JX_ParVector vecX,
                  JX_Int sid, JX_Int pl, JX_Real tol, JX_Int k_dim, JX_Int two_norm, JX_Int is_check_restarted,
                  JX_Bool is_adaptive, JX_Int max_iter_sl, JX_Int *iter_sl_used, JX_Int *num_sl, JX_Int *iter_sl_lvls )
{
   return( jx_PAMGAdapSetup( (void *) solver, (jx_hpCSRMatrix *) matA, (jx_ParVector *) vecB, (jx_ParVector *) vecX,
                                      sid, pl, tol, k_dim, two_norm, is_check_restarted,
                                      is_adaptive, max_iter_sl, iter_sl_used, num_sl, iter_sl_lvls ) );
}

JX_Int
jx_PAMGAdapSetup( void *amg_vdata, jx_hpCSRMatrix *hp_matrix, jx_ParVector *par_rhs, jx_ParVector *par_app,
                  JX_Int solver_id, JX_Int yprint_level, JX_Real ytol, JX_Int yk_dim, JX_Int ytwo_norm, JX_Int yis_check_restarted,
                  JX_Bool yis_adaptive, JX_Int ymax_iter_sl, JX_Int *yiter_sl_used, JX_Int *ynum_sl, JX_Int *iter_sl_lvls )
{
   MPI_Comm comm = jx_hpCSRMatrixComm(hp_matrix);
   jx_ParAMGData *amg_data = amg_vdata;

   /* Data Structure variables */
   jx_hpCSRMatrix **A_array;
   jx_ParVector    **F_array;
   jx_ParVector    **U_array;
   jx_ParVector     *Vtemp;
   jx_ParCSRMatrix **P_array;
   jx_ParCSRMatrix **R_array;
   jx_ParVector     *Residual_array;
   JX_Real          **AI_measure_array;
   JX_Int             **CF_marker_array;
   JX_Int             **dof_func_array;
   JX_Int              *dof_func;
   JX_Int              *col_offd_S_to_A;
   JX_Int              *col_offd_Sabs_to_A = NULL;
   JX_Int            *AIR_maxsize_ls;
   JX_Real            strong_threshold;
   JX_Real            AIR_strong_th;
   JX_Real            max_row_sum;
   JX_Real            trunc_factor;
   JX_Real            S_commpkg_switch;

   JX_Int      max_levels;
   JX_Int      amg_logging;
   JX_Int      amg_print_level;
   JX_Int      debug_flag;
   JX_Int      local_num_vars;
   JX_Int      P_max_elmts;
   JX_Int      R_max_size;

   JX_Solver *smoother = NULL;
   JX_Int        smooth_type = jx_ParAMGDataSmoothType(amg_data);
   JX_Int        smooth_num_levels = jx_ParAMGDataSmoothNumLevels(amg_data);
   char      *euclidfile;
   JX_Int	      eu_level;
   JX_Int	      eu_bj;
   JX_Real     eu_sparse_A;

   JX_Solver ysolver;
   JX_Solver jac_solver;
   JX_Int znum_sl = 0;
   JX_Int iter_sl_used = 0;
   JX_Int ziter_sl_used = 0;
   JX_Real iter_rest;
   JX_Real res_norm;
   JX_Real res_norm_old;
   JX_Real reduce_factor;
   JX_Real reduce_factor_thrd = 1.0;

   /* Local variables */
   JX_Real           *AI_measure;
   JX_Real           *AI_measure2;
   JX_Int              *CF_marker;
   JX_Int              *CFN_marker;
   jx_ParCSRMatrix  *par_S;
   jx_ParCSRMatrix  *Sabs = NULL;
   jx_ParCSRMatrix  *par_S2;
   jx_ParCSRMatrix  *par_P;
   jx_ParCSRMatrix  *R;
   jx_hpCSRMatrix  *A_H;

   JX_Int       wall_time_option;         /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_coarsen = 0.0;  /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_rap = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_interp = 0.0;   /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    tmp_wall_time = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */

   JX_Int       old_num_levels, num_levels;
   JX_Int       level;
   JX_Int       local_size, i;
   JX_Int       first_local_row;
   JX_Int       coarse_size_last = -1;
   JX_Int       coarse_size;
   JX_Int       coarsen_type;
   JX_Int       interp_type;
   JX_Int       restri_type;
   JX_Int       measure_type;
   JX_Int       setup_type;
   JX_Int       fine_size;
   JX_Int       rest, tms, indx;
   JX_Int       not_finished_coarsening;
   JX_Int       Setup_err_flag = 0;
   JX_Int       coarse_threshold; /* we don't fix it to be 9.  peghoty 2010/04/14 */
   JX_Int       coarse_threshold_2; /* by xu: 2013/11/25 */
   JX_Int       j, k;
   JX_Int       num_procs, my_id, num_threads;
   JX_Int      *grid_relax_type = jx_ParAMGDataGridRelaxType(amg_data);
   JX_Int       num_functions = jx_ParAMGDataNumFunctions(amg_data);
   JX_Int       spmt_rap_type = jx_ParAMGDataSpMtRapType(amg_data);
   JX_Int       ai_measure_type = jx_ParAMGDataAIMeasureType(amg_data);
   JX_Int       num_paths = jx_ParAMGDataNumPaths(amg_data);
   JX_Int       agg_num_levels = jx_ParAMGDataAggNumLevels(amg_data);
   JX_Int       agg_interp_type = jx_ParAMGDataAggInterpType(amg_data);
   JX_Int       agg_P_max_elmts = jx_ParAMGDataAggPMaxElmts(amg_data);
   //JX_Int       agg_P12_max_elmts = jx_ParAMGDataAggP12MaxElmts(amg_data);
   JX_Real    agg_trunc_factor = jx_ParAMGDataAggTruncFactor(amg_data);
   //JX_Real    agg_P12_trunc_factor = jx_ParAMGDataAggP12TruncFactor(amg_data);
   JX_Int       rap2 = jx_ParAMGDataRAP2(amg_data);
   JX_Int       keepTranspose = jx_ParAMGDataKeepTranspose(amg_data);
   JX_Int       print_coarse_matrix = jx_ParAMGDataPrintCoarseSystem(amg_data);
   JX_Int      *opt_icor;
   JX_Int      *coarse_dof_func;
   JX_Int      *coarse_pnts_global;
   JX_Int	    *coarse_pnts_global1;
   JX_Real    size;
   JX_Real    coarse_ratio;           /* newly added peghoty 2010/04/14 */
   JX_Real    wall_time = 0.0;        /* for debugging instrumentation */
   JX_Int       measure_type_rlx;       // newly added peghoty 2010/05/29
   JX_Int       number_syn_rlx;         // newly added peghoty 2010/05/29
   JX_Real    measure_threshold_rlx;  // newly added peghoty 2010/05/29

   /* coarse matrices */
   char FileNameCoaMat[256];

   /* ai statistic information*/
   JX_Int       num_vars_local = 0, num_vars_global;
   JX_Int       num_ai_local = 0, num_ai_global;
   JX_Int       num_ai_local_valid = 0, num_ai_global_valid;
   JX_Int       num_ai_c_local = 0, num_ai_c_global;
   JX_Real    mai_local = 0.0, mai_global;
   JX_Real    mai_local_valid = 0.0, mai_global_valid;
   JX_Real    mai_c_local = 0.0, mai_c_global;

   JX_Int       num_vars_local_0 = 0;
   JX_Int       num_ai_local_0 = 0;
   JX_Int       num_ai_local_valid_0 = 0;
   JX_Int       num_ai_c_local_0 = 0;
   JX_Real    mai_local_0 = 0.0;
   JX_Real    mai_local_valid_0 = 0.0;
   JX_Real    mai_c_local_0 = 0.0;
   JX_Real    mai_threshold = 0.1;

   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);
   if ((num_procs > 1) && (spmt_rap_type != 1))  /* Yue Xiaoqiang 2012/10/12 */
   {
      spmt_rap_type = 1;
   }
   num_threads = jx_NumThreads();               /* Yue Xiaoqiang 2012/10/12 */
   opt_icor = jx_CTAlloc(JX_Int, 5*num_threads+2); /* Yue Xiaoqiang 2012/10/12 */

   wall_time_option = jx_ParAMGDataWallTimeOption(amg_data);
   old_num_levels   = jx_ParAMGDataNumLevels(amg_data);
   max_levels       = jx_ParAMGDataMaxLevels(amg_data);
   amg_logging      = jx_ParAMGDataLogging(amg_data);
   amg_print_level  = jx_ParAMGDataPrintLevel(amg_data);
   coarsen_type     = jx_ParAMGDataCoarsenType(amg_data);
   measure_type     = jx_ParAMGDataMeasureType(amg_data);
   setup_type       = jx_ParAMGDataSetupType(amg_data);
   debug_flag       = jx_ParAMGDataDebugFlag(amg_data);
   dof_func         = jx_ParAMGDataDofFunc(amg_data);
   interp_type      = jx_ParAMGDataInterpType(amg_data);
   restri_type      = jx_ParAMGDataRestriction(amg_data);
   AIR_strong_th    = jx_ParAMGDataAIRStrongTh(amg_data);
   euclidfile       = jx_ParAMGDataEuclidFile(amg_data);
   eu_level         = jx_ParAMGDataEuLevel(amg_data);
   eu_bj            = jx_ParAMGDataEuBJ(amg_data);
   eu_sparse_A      = jx_ParAMGDataEuSparseA(amg_data);
   coarse_threshold = jx_ParAMGDataCoarseThreshold(amg_data);           /* newly added peghoty 2010/04/14 */
   coarse_ratio     = jx_ParAMGDataCoarseRatio(amg_data);               /* newly added peghoty 2010/04/14 */
   measure_type_rlx = jx_ParAMGDataMeasureTypeRlx(amg_data);            /* newly added peghoty 2010/05/29 */
   number_syn_rlx   = jx_ParAMGDataNumberSynRlx(amg_data);              /* newly added peghoty 2010/05/29 */
   measure_threshold_rlx = jx_ParAMGDataMeasureThresholdRlx(amg_data);  /* newly added peghoty 2010/05/29 */

   coarse_threshold_2 = -1;

   jx_ParCSRMatrixSetNumNonzeros(jx_hpCSRMatrixPar(hp_matrix));
   jx_ParCSRMatrixSetDNumNonzeros(jx_hpCSRMatrixPar(hp_matrix));
   jx_ParAMGDataNumVariables(amg_data) = jx_hpCSRMatrixNumRows(hp_matrix);

   if (setup_type == 0) 
   {
      return(Setup_err_flag);
   }

   par_S = NULL;
   A_H = NULL;

   A_array = jx_hpAMGDataAArray(amg_data);
   P_array = jx_ParAMGDataPArray(amg_data);
   R_array = jx_ParAMGDataRArray(amg_data);
   AIR_maxsize_ls = jx_ParAMGDataAIRMaxSizeLS(amg_data);
   AI_measure_array = jx_ParAMGDataAIMeasureArray(amg_data);
   CF_marker_array = jx_ParAMGDataCFMarkerArray(amg_data);
   dof_func_array  = jx_ParAMGDataDofFuncArray(amg_data);
   local_size = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(hp_matrix));

   grid_relax_type[3] = jx_ParAMGDataUserCoarseRelaxType(amg_data); 

   if (A_array || P_array || R_array || AIR_maxsize_ls || AI_measure_array || CF_marker_array || dof_func_array)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (A_array[j])
         {
            jx_hpCSRMatrixDestroy(A_array[j]);
            A_array[j] = NULL;
         }
       
         if (dof_func_array[j])
         {
            jx_TFree(dof_func_array[j]);
            dof_func_array[j] = NULL;
         }
      }

      for (j = 0; j < old_num_levels - 1; j ++)
      {
         if (P_array[j])
         {
            jx_ParCSRMatrixDestroy(P_array[j]);
            P_array[j] = NULL;
         }
         if (R_array[j])
         {
            jx_ParCSRMatrixDestroy(R_array[j]);
            R_array[j] = NULL;
         }
      }

      jx_TFree(AIR_maxsize_ls);
      AIR_maxsize_ls = NULL;

     /*-------------------------------------------------------------------
      *  Special case use of CF_marker_array when old_num_levels = 1
      *  requires us to attempt this deallocation every time
      *------------------------------------------------------------------*/
      
      if (CF_marker_array[0])
      {
         jx_TFree(CF_marker_array[0]);
         CF_marker_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (CF_marker_array[j])
         {
            jx_TFree(CF_marker_array[j]);
            CF_marker_array[j] = NULL;
         }
      }
      
      /* for AI_measure: added by xwxu */
      if (AI_measure_array[0])
      {
         jx_TFree(AI_measure_array[0]);
         AI_measure_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (AI_measure_array[j])
         {
            jx_TFree(AI_measure_array[j]);
            AI_measure_array[j] = NULL;
         }
      }
   }

   if (A_array == NULL)
   {
      A_array = jx_CTAlloc(jx_hpCSRMatrix*, max_levels);
   }

   if (P_array == NULL && max_levels > 1)
   {
      P_array = jx_CTAlloc(jx_ParCSRMatrix*, max_levels - 1);
   }

   /* If retri_type != 0, R != P^T, allocate R matrices */
   if (restri_type)
   {
      if (R_array == NULL && max_levels > 1)
      {
         R_array = jx_CTAlloc(jx_ParCSRMatrix*, max_levels-1);
         AIR_maxsize_ls = jx_CTAlloc(JX_Int, max_levels-1);
      }
   }

   if (AI_measure_array == NULL)
   {
      AI_measure_array = jx_CTAlloc(JX_Real*, max_levels);
   }
   if (CF_marker_array == NULL)
   {
      CF_marker_array = jx_CTAlloc(JX_Int*, max_levels);
   }
   if (dof_func_array == NULL)
   {
      dof_func_array = jx_CTAlloc(JX_Int*, max_levels);
   }

   if (num_functions > 1 && dof_func == NULL)
   {
      first_local_row = jx_hpCSRMatrixFirstRowIndex(hp_matrix);
      dof_func = jx_CTAlloc(JX_Int, local_size);
      rest = first_local_row-((first_local_row / num_functions)*num_functions);
      indx = num_functions - rest;
      if (rest == 0) 
      {
         indx = 0;
      }
      k = num_functions - 1;
      for (j = indx - 1; j > -1; j --)
      {
         dof_func[j] = k --;
      }
      tms = local_size / num_functions;
      if (tms*num_functions + indx > local_size) 
      {
         tms --;
      }
      for (j = 0; j < tms; j ++)
      {
         for (k = 0; k < num_functions; k ++)
         {
            dof_func[indx++] = k;
         }
      }
      k = 0;
      while (indx < local_size)
      {
         dof_func[indx++] = k ++;
      }
      jx_ParAMGDataDofFunc(amg_data) = dof_func;
   }

   A_array[0] = hp_matrix;

   dof_func_array[0] = dof_func;
   jx_ParAMGDataAIMeasureArray(amg_data) = AI_measure_array;
   jx_ParAMGDataCFMarkerArray(amg_data) = CF_marker_array;
   jx_ParAMGDataDofFuncArray(amg_data) = dof_func_array;
   jx_hpAMGDataAArray(amg_data) = A_array;
   jx_ParAMGDataPArray(amg_data) = P_array;
   /* If R != P^T */
   if (restri_type)
   {
      jx_ParAMGDataRArray(amg_data) = R_array;
      jx_ParAMGDataAIRMaxSizeLS(amg_data) = AIR_maxsize_ls;
   }
   else
   {
      jx_ParAMGDataRArray(amg_data) = P_array;
   }

   Vtemp = jx_ParAMGDataVtemp(amg_data);

   if (Vtemp != NULL)
   {
      jx_ParVectorDestroy(Vtemp);
      Vtemp = NULL;
   }

   Vtemp = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[0]),
                              jx_hpCSRMatrixGlobalNumRows(A_array[0]),
                              jx_hpCSRMatrixRowStarts(A_array[0]));
   jx_ParVectorInitialize(Vtemp);
   jx_ParVectorSetPartitioningOwner(Vtemp,0);
   jx_ParAMGDataVtemp(amg_data) = Vtemp;

   F_array = jx_ParAMGDataFArray(amg_data);
   U_array = jx_ParAMGDataUArray(amg_data);

   if (F_array != NULL || U_array != NULL)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (F_array[j] != NULL)
         {
            jx_ParVectorDestroy(F_array[j]);
            F_array[j] = NULL;
         }
         if (U_array[j] != NULL)
         {
            jx_ParVectorDestroy(U_array[j]);
            U_array[j] = NULL;
         }
      }
   }

   if (F_array == NULL)
   {
      F_array = jx_CTAlloc(jx_ParVector*, max_levels);
   }
   if (U_array == NULL)
   {
      U_array = jx_CTAlloc(jx_ParVector*, max_levels);
   }

   F_array[0] = par_rhs;
   U_array[0] = par_app;

   jx_ParAMGDataFArray(amg_data) = F_array;
   jx_ParAMGDataUArray(amg_data) = U_array;

  /*----------------------------------------------------------
   *   Initialize jx_ParAMGData
   *---------------------------------------------------------*/
   not_finished_coarsening = 1;
   level = 0;
   strong_threshold = jx_ParAMGDataStrongThreshold(amg_data);
   max_row_sum = jx_ParAMGDataMaxRowSum(amg_data);
   trunc_factor = jx_ParAMGDataTruncFactor(amg_data);
   P_max_elmts = jx_ParAMGDataPMaxElmts(amg_data);
   S_commpkg_switch = jx_ParAMGDataSCommPkgSwitch(amg_data);
   if (smooth_num_levels > level)
   {
      smoother = jx_CTAlloc(JX_Solver, smooth_num_levels);
      jx_ParAMGDataSmoother(amg_data) = smoother;
   }

   /* if AI-based coarsening is used, set ai_measure_type = 1. */
   if (coarsen_type == 990 || coarsen_type == 991 || coarsen_type == 993 || 
       coarsen_type == 96 || coarsen_type == 98 || coarsen_type == 910 || 
       coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || 
       coarsen_type == 938 || coarsen_type == 968 || 
       amg_print_level > 0)
   {
      ai_measure_type = 1;
   }

  /*-----------------------------------------------------
   *   Enter Coarsening Loop
   *-----------------------------------------------------*/

   if (amg_print_level > 0 && my_id ==0)
   {
      jx_printf(" \n");
      jx_printf("================================================================== \n");
      jx_printf("+++++++++++++ multi-scale/AI infomation for levels +++++++++++++++ \n");
   }

   while (not_finished_coarsening)
   {
      fine_size = jx_hpCSRMatrixGlobalNumRows(A_array[level]);

      if (level > 0)
      {
         F_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                             jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jx_hpCSRMatrixRowStarts(A_array[level]));
         jx_ParVectorInitialize(F_array[level]);
         jx_ParVectorSetPartitioningOwner(F_array[level], 0);

         U_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                             jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jx_hpCSRMatrixRowStarts(A_array[level]));
         jx_ParVectorInitialize(U_array[level]);
         jx_ParVectorSetPartitioningOwner(U_array[level], 0);
      }

      if (level > 0)
      {
   iter_sl_used = 0;
   reduce_factor_thrd = 1.0;
         jx_ParVectorCopy(F_array[level-1], Vtemp);
         jx_ParCSRMatrixMatvec(-1.0, jx_hpCSRMatrixPar(A_array[level-1]), U_array[level-1], 1.0, Vtemp);
         jx_ParCSRMatrixMatvecT(1.0, P_array[level-1], Vtemp, 0.0, F_array[level]);
   JX_JacobiCreate(&jac_solver);
   JX_JacobiSetMaxIter(jac_solver, 1);
   JX_JacobiSetup(jac_solver, (JX_hpCSRMatrix) (A_array[level]));
         if (solver_id == 1)
         {
      JX_PCGCreate(comm, &ysolver);
      JX_PCGSetMaxIter(ysolver, 1);
      JX_PCGSetTol(ysolver, ytol);
      JX_PCGSetTwoNorm(ysolver, ytwo_norm);
      JX_PCGSetPrintLevel(ysolver, yprint_level);
      JX_PCGSetLogging(ysolver, 1);
      JX_PCGSetPrecond(ysolver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_PCGSetup(ysolver, (JX_Matrix) A_array[level], (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_PCGSolve(ysolver, (JX_Matrix) A_array[level], (JX_Matrix) A_array[level],
                              (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
         JX_PCGGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JX_PCGDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JX_PCGDestroy(ysolver);
         break;
      }
         }
         else if (solver_id == 2)
         {
      JX_GMRESCreate(comm, &ysolver);
      JX_GMRESSetMaxIter(ysolver, 1);
      JX_GMRESSetPrintLevel(ysolver, yprint_level);
      JX_GMRESSetKDim(ysolver, yk_dim);
      JX_GMRESSetIsCheckRestarted(ysolver, yis_check_restarted);
      JX_GMRESSetTol(ysolver, ytol);
      JX_GMRESSetLogging(ysolver, 1);
      JX_GMRESSetPrecond(ysolver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_GMRESSetup(ysolver, (JX_Matrix) A_array[level], (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_GMRESSolve(ysolver, (JX_Matrix) A_array[level], (JX_Matrix) A_array[level],
                                (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
         JX_GMRESGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JX_GMRESDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JX_GMRESDestroy(ysolver);
         break;
      }
         }
         else if (solver_id == 3)
         {
      JX_BiCGSTABCreate(comm, &ysolver);
      JX_BiCGSTABSetMaxIter(ysolver, 1);
      JX_BiCGSTABSetPrintLevel(ysolver, yprint_level);
      JX_BiCGSTABSetTol(ysolver, ytol);
      JX_BiCGSTABSetAbsoluteTol(ysolver, 0.0);
      JX_BiCGSTABSetConvCriteria(ysolver, 0);
      JX_BiCGSTABSetLogging(ysolver, 1);
      JX_BiCGSTABSetPrecond(ysolver, (JX_PtrToSolverFcn) JX_JacobiPrecond, (JX_PtrToSolverFcn) JX_JacobiSetup, jac_solver);
      JX_BiCGSTABSetup(ysolver, (JX_Matrix) A_array[level], (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JX_BiCGSTABSolve(ysolver, (JX_Matrix) A_array[level], (JX_Matrix) A_array[level],
                                   (JX_Vector) F_array[level], (JX_Vector) U_array[level]);
         JX_BiCGSTABGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JX_BiCGSTABDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JX_BiCGSTABDestroy(ysolver);
         break;
      }
         }
         JX_JacobiDestroy(jac_solver);
      }

     /*--------------------------------------------------------------
      *  Select coarse-grid points on 'level' : returns CF_marker
      *  for the level.  Returns strength matrix, par_S 
      *  Returns AI_Measure. 
      *--------------------------------------------------------------*/
     
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();
      if (debug_flag == 3)
      {
         jx_printf("\n ===== Proc = %d     Level = %d  =====\n",my_id, level);
         fflush(NULL);
      }

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      if (max_levels > 1)
      {
         local_num_vars = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level]));

        /*--------------------------------------------------------------------------
         *  Get the Strength Matrix 
         *-------------------------------------------------------------------------*/        
         jx_PAMGCreateS(jx_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }
         /* for AIR, need absolute value SOC */
         if (restri_type)
         {
            jx_PAMGCreateSabs(jx_hpCSRMatrixPar(A_array[level]),AIR_strong_th,1.0,num_functions,dof_func_array[level],&Sabs);
            col_offd_Sabs_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]), Sabs, &col_offd_Sabs_to_A);
            }
         }

         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jx_PAMGMeasureAI(par_S, jx_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

        /*--------------------------------------------------------------------------
         *    Do the appropriate coarsening as follows:
         *
         *      coarsen_type =  0: CLJP 
         *      coarsen_type =  1: Ruge
         *      coarsen_type = 11: Ruge 1st pass only
         *      coarsen_type =  2: Ruge2B
         *      coarsen_type =  3: Ruge3
         *      coarsen_type =  4: Ruge3c
         *      coarsen_type =  5: Ruge relax special points
         *      coarsen_type =  6: Falgout
         *      coarsen_type =  8: PMIS
         *      coarsen_type = 10: HMIS
         *      coarsen_type = 90: RCLJP 
         *      coarsen_type = 91: RRS0
         *      coarsen_type = 990: CLJP_AI 
         *      coarsen_type = 991: Ruge_AI
         *      coarsen_type = 993: Ruge3_AI
         *      coarsen_type = 96: Falgout_AI
         *      coarsen_type = 98: PMIS_AI
         *      coarsen_type = 910: HMIS_AI
         *      coarsen_type = 908, 918, 928, 938, 968: AI-TYPE.
         *                                          peghoty 2010/05/29
         *-------------------------------------------------------------------------*/ 

         //jx_printf("=========== coarsen_type = %d \n", coarsen_type); 

         if (coarsen_type == 6)  
         {
            jx_PAMGCoarsenFalgout(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 96) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenFalgoutAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 8) 
         {
            jx_PAMGCoarsenPMIS(par_S, jx_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 10) 
         {
            jx_PAMGCoarsenHMIS(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }      
         else if ((coarsen_type == 910) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenHMISAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, debug_flag, &CF_marker);
         }      
         else if (coarsen_type == 1 || coarsen_type == 2 ||
                  coarsen_type == 3 || coarsen_type == 4 ||
                  coarsen_type == 5 || coarsen_type == 11  ) 
         {
            jx_PAMGCoarsenRuge(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, coarsen_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 991) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 11, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 992) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 1, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 993) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 3, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 0)    
         {
            jx_PAMGCoarsen(par_S,jx_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 990) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 90)
         {
            jx_PAMGCoarsenRCLJP(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type_rlx, number_syn_rlx,
                                measure_threshold_rlx, 0, debug_flag, &CF_marker);
         }         
         else if (coarsen_type == 91)
         {
            jx_PAMGCoarsenRRS0(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, measure_type_rlx,
                               number_syn_rlx, measure_threshold_rlx, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 98) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenPMISAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || coarsen_type == 938 || 
                  coarsen_type == 968) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenXML(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, coarsen_type, debug_flag, &CF_marker);
         }      

         if (level < agg_num_levels)
         {
            jx_PAMGCoarseParms(comm, local_num_vars, 1,
                               dof_func_array[level], CF_marker,
                               &coarse_dof_func, &coarse_pnts_global1);
            jx_PAMGCreate2ndS(par_S, CF_marker, num_paths, coarse_pnts_global1, &par_S2);

           /*--------------------------------------------------------------------------
            *  Get the AI Measure
            *-------------------------------------------------------------------------*/
            if (ai_measure_type == 1)
            {
               jx_PAMGMeasureAI(par_S2, par_S2, debug_flag, &AI_measure2);
            }

            if (coarsen_type == 10)
            {
               jx_PAMGCoarsenHMIS(par_S2, par_S2, measure_type+3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 910) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenHMISAI(par_S2, par_S2, AI_measure2, 0, measure_type+3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 8)
            {
               jx_PAMGCoarsenPMIS(par_S2, par_S2, 3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 98) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenPMISAI(par_S2, par_S2, AI_measure2, 3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 6)
            {
               jx_PAMGCoarsenFalgout(par_S2, par_S2, measure_type, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 96) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenFalgoutAI(par_S2, par_S2, AI_measure2, measure_type, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 0)
            {
               jx_PAMGCoarsen(par_S2, par_S2, 0, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 990) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenAI(par_S2, par_S2, AI_measure2, 0, debug_flag, &CFN_marker);
            }
            jx_ParCSRMatrixDestroy(par_S2);
            if (ai_measure_type == 1)
            {
               jx_TFree(AI_measure2);
            }
         }

         /* Here for changes of min_coarse_size */

         if (level < agg_num_levels)
         {
            if (agg_interp_type == 4)
            {
               jx_PAMGCorrectCFMarker(CF_marker, local_num_vars, CFN_marker);
               jx_TFree(coarse_pnts_global1);
               jx_TFree(CFN_marker);
            }
         }

        /*--------------------------------------------
         *  store the CF array
         *-------------------------------------------*/ 
         CF_marker_array[level] = CF_marker;

         if (debug_flag == 1)
         {
            wall_time = jx_time_getWallclockSeconds() - wall_time;
            jx_printf("Proc = %d    Level = %d    Coarsen Time = %f\n",my_id,level, wall_time); 
            fflush(NULL);
         }

         if (wall_time_option == 1)
         {
            wall_time_coarsen += (jx_time_getWallclockSeconds() - tmp_wall_time);
         }

        /*-------------------------------------------------------------------------
         *  Get the coarse parameters
         *-------------------------------------------------------------------------*/ 
         jx_PAMGCoarseParms(comm, local_num_vars, num_functions, 
                            dof_func_array[level], CF_marker,
                            &coarse_dof_func, &coarse_pnts_global);

         dof_func_array[level+1] = NULL;
         if (num_functions > 1) dof_func_array[level+1] = coarse_dof_func;

#ifdef JX_NO_GLOBAL_PARTITION
         if (my_id == (num_procs -1)) 
         {
            coarse_size = coarse_pnts_global[1];
         }
         jx_MPI_Bcast(&coarse_size, 1, JX_MPI_INT, num_procs-1, comm);
#else
         coarse_size = coarse_pnts_global[num_procs];
#endif
      
        if ( coarse_size <= coarse_threshold && coarse_size_last >= coarse_size ) {
           break;  
        } 

      }
      else  /* max_levels = 1 */
      {

        /*--------------------------------------------------------------------------
         *  Get the Strength Matrix 
         *-------------------------------------------------------------------------*/        
         jx_PAMGCreateS(jx_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }

         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jx_PAMGMeasureAI(par_S, jx_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         par_S = NULL;
         coarse_pnts_global = NULL;
         CF_marker = jx_CTAlloc(JX_Int, local_size );
         for (i = 0; i < local_size ; i ++) 
         {
            CF_marker[i] = 1;
         }
         CF_marker_array[level] = CF_marker;
         coarse_size = fine_size;

         local_num_vars = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level]));
      }

     /*-------------------------------------------------------------------------
      *  compute ai-information for all-pnts and C-pnts!
      *-------------------------------------------------------------------------*/ 
     
      if ((max_levels == 1 || level != max_levels-1) && (ai_measure_type == 1))
      { 
         for (i = 0; i < local_num_vars; i++) {
             num_vars_local++ ;
             //if (AI_measure[i] > 0.0 && AI_measure[i] < 1.0)
             if (AI_measure[i] > 1.0e-10) 
             {
                mai_local += AI_measure[i];
                num_ai_local++ ;
                if (AI_measure[i] >= mai_threshold) {
                   mai_local_valid += AI_measure[i];
                   num_ai_local_valid++ ;
                   if (CF_marker[i] == 1) {
                      //if (level == 0 ) jx_printf("i_c = %d, ai_measure= %f \n", i, AI_measure[i]);
                      mai_c_local += AI_measure[i];
                      num_ai_c_local++;
                   }
                }
             }
         }
      }

      if (level == 0) {
        num_ai_local_0 = num_ai_local;
        num_ai_local_valid_0 = num_ai_local_valid;
        num_vars_local_0 = num_vars_local;
        num_ai_c_local_0 = num_ai_c_local;
        mai_local_0 = mai_local;
        mai_local_valid_0 = mai_local_valid;
        mai_c_local_0 = mai_c_local;
      }


     /*-------------------------------------------------------------------------
      *  if no coarse-grid, stop coarsening, and set the
      *  coarsest solve to be a single sweep of Jacobi !
      *-------------------------------------------------------------------------*/ 

      if ( (coarse_size == 0) || (coarse_size == fine_size) )
      {
         JX_Int  *num_grid_sweeps   = jx_ParAMGDataNumGridSweeps(amg_data);
         JX_Int **grid_relax_points = jx_ParAMGDataGridRelaxPoints(amg_data);
         if (grid_relax_type[3] == 9)
	 {
	    grid_relax_type[3] = grid_relax_type[0];
	    num_grid_sweeps[3] = 1;
	    if (grid_relax_points)
	    {
	       grid_relax_points[3][0] = 0; 
	    }
	 }
         if (par_S) 
         {
            jx_ParCSRMatrixDestroy(par_S);
         }
         jx_TFree(coarse_pnts_global);
         if (level > 0)
         {
           /*-------------------------------------------------------------
            * Note special case treatment of CF_marker is necessary
            * to do CF relaxation correctly when num_levels = 1
            *------------------------------------------------------------*/ 
            jx_TFree(CF_marker_array[level]);
            jx_ParVectorDestroy(F_array[level]);
            jx_ParVectorDestroy(U_array[level]);
         }
         break; 
      }


      /* Build restriction */
      if (restri_type)
      {
         /* !!! Ensure that CF_marker contains -1 or 1 !!! */
         for (i = 0; i < jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level])); i ++)
         {
            CF_marker[i] = CF_marker[i] > 0 ? 1 : -1;
         }
         if (restri_type == 1) /* distance-1 AIR */
         {
            jx_PAMGBuildRestrAIR(jx_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                 Sabs, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], 
                                 debug_flag, trunc_factor, P_max_elmts, 
                                 col_offd_Sabs_to_A, &R, &R_max_size );
         }
         else /* distance-2 AIR */
         {
            jx_PAMGBuildRestrDist2AIR(jx_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                      Sabs, coarse_pnts_global, num_functions, 
                                      dof_func_array[level], 
                                      debug_flag, trunc_factor, P_max_elmts, 
                                      col_offd_Sabs_to_A, &R, &R_max_size );
         }
         if (Sabs)
         {
            jx_ParCSRMatrixDestroy(Sabs);
            Sabs = NULL;
         }
         jx_TFree(col_offd_Sabs_to_A);
      }


     /*-----------------------------------------------------------------------------------
      *   Build prolongation matrix, P, and place in P_array[level]
      *
      *     interp_type =  0: modified classical interpolation
      *     interp_type =  3: direct interpolation (with separation of weights)
      *     interp_type =  4: multipass interpolation
      *     interp_type =  5: multipass interpolation (with separation of weights)
      *     interp_type =  6: extended classical modified interpolation
      *     interp_type =  7: extended (if no common C neighbor) classical 
      *                       modified interpolation
      *     interp_type =  8: standard interpolation
      *     interp_type =  9: standard interpolation (with separation of weights)
      *-----------------------------------------------------------------------------------*/
       
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      if (level < agg_num_levels)
      {
         if (agg_interp_type == 4)
         {
            jx_PAMGBuildMultipass(jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S, coarse_pnts_global,
                                  num_functions, dof_func_array[level], debug_flag,
                                  agg_trunc_factor, agg_P_max_elmts, 0, col_offd_S_to_A, &par_P);
         }
      }
      else
      {
         if (interp_type == 0)  // classical modified interpolation
         {
            if (spmt_rap_type == 1)
            {
               jx_PAMGBuildInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, trunc_factor,
                                P_max_elmts, col_offd_S_to_A, &par_P ); 
            }
            else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
            {
               jx_PAMGBuildInterp1( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                 coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag, trunc_factor,
                                 P_max_elmts, col_offd_S_to_A, &par_P, opt_icor );
            }
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 3)  // direct interpolation (with separation of weights)
         {
            jx_PAMGBuildDirInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, 
                                trunc_factor, P_max_elmts, 
                                col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         } 
         else if (interp_type == 4)  // multipass interpolation
         {
            jx_PAMGBuildMultipass( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 5) // multipass interpolation (with separation of weights)
         {
            jx_PAMGBuildMultipass( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 6)  // extended classical modified interpolation
         {
             jx_PAMGBuildExtPIInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                   par_S, coarse_pnts_global, num_functions, 
                                   dof_func_array[level], debug_flag,
                                   trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
             jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 7) // extended (if no common C neighbor) classical modified interpolation
         {
            jx_PAMGBuildExtPICCInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                     par_S, coarse_pnts_global, num_functions, 
                                     dof_func_array[level], debug_flag,
                                     trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 8) // standard interpolation
         {
            jx_PAMGBuildStdInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 9) // standard interpolation (with separation of weights)
         {
            jx_PAMGBuildStdInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 100) /* 1pt interpolation */
         {
            jx_PAMGBuildInterpOnePnt( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                    par_S, coarse_pnts_global, num_functions,
                                    dof_func_array[level], debug_flag, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
      }
      
      if (debug_flag == 1)
      {
         wall_time = jx_time_getWallclockSeconds() - wall_time;
         jx_printf("Proc = %d    Level = %d    Build Interp Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_interp += (jx_time_getWallclockSeconds() - tmp_wall_time);
      }

      P_array[level] = par_P; 
      if (restri_type)
      {
         R_array[level] = R;
         AIR_maxsize_ls[level] = R_max_size;
      }

      if (par_S) 
      {
         jx_ParCSRMatrixDestroy(par_S);
      }
      par_S = NULL;


     /*--------------------------------------------------------------------------------
      *   Build coarse-grid operator, A_array[level+1] by R*A*P
      *-------------------------------------------------------------------------------*/
       
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      if (restri_type)
      {
         /* Use two matrix products to generate A_H */
         jx_ParCSRMatrix *Q = jx_ParMatmul(jx_hpCSRMatrixPar(A_array[level]), P_array[level]);
         A_H = jx_hpInithpCSRMatrix();
         jx_hpCSRMatrixPar(A_H) = jx_ParMatmul(R_array[level], Q);
         jx_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jx_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jx_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         jx_ParCSRMatrixOwnsRowStarts(R_array[level]) = 0;
         if (num_procs > 1) jx_MatvecCommPkgCreate(jx_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jx_ParCSRMatrixDestroy(Q);
      }
      else if (rap2)
      {
         /* Use two matrix products to generate A_H */
         jx_ParCSRMatrix *Q = jx_ParMatmul(jx_hpCSRMatrixPar(A_array[level]), P_array[level]);
         A_H = jx_hpInithpCSRMatrix();
         jx_hpCSRMatrixPar(A_H) = jx_ParTMatmul(P_array[level], Q);
         jx_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jx_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jx_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         if (num_procs > 1) jx_MatvecCommPkgCreate(jx_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jx_ParCSRMatrixDestroy(Q);
      }
      else if (spmt_rap_type == 1)
      {
         A_H = jx_hpInithpCSRMatrix();
         jx_PAMGBuildCoarseOperatorKT( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &jx_hpCSRMatrixPar(A_H) );
      }
      else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
      {
         A_H = jx_hpInithpCSRMatrix();
         jx_PAMGBuildCoarseOperatorOMP( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], &jx_hpCSRMatrixPar(A_H), opt_icor );     
      }
      
      if (debug_flag == 1)
      {
         wall_time = jx_time_getWallclockSeconds() - wall_time;
         jx_printf("Proc = %d    Level = %d    Build Coarse Operator Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_rap += (jx_time_getWallclockSeconds() - tmp_wall_time);
      }

      ++ level;

      jx_ParCSRMatrixSetNumNonzeros(jx_hpCSRMatrixPar(A_H));
      jx_ParCSRMatrixSetDNumNonzeros(jx_hpCSRMatrixPar(A_H));
      A_array[level] = A_H;

      if (coarse_size <= coarse_threshold_2) {

         coarse_size_last = coarse_size;

      }

      /* print coarser operator. */
	  if (print_coarse_matrix == 1) {

         jx_sprintf(FileNameCoaMat, "cmat_%d", level);
         jx_hpCSRMatrixPrint(A_array[level], FileNameCoaMat);

      }
      
     /*-------------------------------------------------------------------------------
      *   Switch to CLJP when coarsening slows
      *                                            peghoty  2009/07/09
      *------------------------------------------------------------------------------*/ 
      
      size = ((JX_Real) fine_size )*coarse_ratio;   /* peghoty 2010/04/14 */
      if (coarsen_type > 0 && coarse_size >= (JX_Int) size)
      {
	 coarsen_type = 0;      
      }

      /* How to stop the loop "while" */
      if ( (level == max_levels-1) || (coarse_size <= coarse_threshold) )
      {
         not_finished_coarsening = 0;

#if 0
         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the Strength Matrix 
            *-------------------------------------------------------------------------*/        
            jx_PAMGCreateS(A_array[level],strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
            col_offd_S_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jx_PAMGCreateSCommPkg(A_array[level], par_S, &col_offd_S_to_A);
            }

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jx_PAMGMeasureAI(par_S, A_array[level], debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         /*-------------------------------------------------------------------------
          *  compute ai-information for all-pnts and C-pnts!
          *-------------------------------------------------------------------------*/ 
         local_num_vars = jx_CSRMatrixNumRows(jx_ParCSRMatrixDiag(A_array[level]));

          for (i = 0; i < local_num_vars; i++) {
              num_vars_local++ ;
              //if (AI_measure[i] > 0.0 && AI_measure[i] < 1.0)
              if (AI_measure[i] > 1.0e-10) 
              {
                 mai_local += AI_measure[i];
                 num_ai_local++ ;
                 if (AI_measure[i] > 0.1) {
                    mai_local_valid += AI_measure[i];
                    num_ai_local_valid++ ;
                    mai_c_local += AI_measure[i];
                    num_ai_c_local++;
                 }
              }
          }
#endif
      }

      /*jx_printf("level = %d, not_finished_coarsening = %d \n", level, not_finished_coarsening);
      if (level == 1) {
        num_ai_local_0 = num_ai_local;
        num_ai_local_valid_0 = num_ai_local_valid;
        num_vars_local_0 = num_vars_local;
        num_ai_c_local_0 = num_ai_c_local;
        mai_local_0 = mai_local;
        mai_local_valid_0 = mai_local_valid;
        mai_c_local_0 = mai_c_local;
      }*/

   } // end while loop

  *yiter_sl_used = ziter_sl_used;
  *ynum_sl = znum_sl;

   if (amg_print_level > 0)
   {

      jx_MPI_Allreduce(&num_ai_local_0,&num_ai_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_local_valid_0,&num_ai_global_valid,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_vars_local_0,&num_vars_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_c_local_0,&num_ai_c_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_0,&mai_global,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_valid_0,&mai_global_valid,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_c_local_0,&mai_c_global,1,JX_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JX_Real num_ai = num_ai_global;
         JX_Real num_ai_valid = num_ai_global_valid;
         JX_Real num_var = num_vars_global;
         JX_Real num_ai_c = num_ai_c_global;
         jx_printf(" \n");
         jx_printf("Level 0: \n");
         jx_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jx_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jx_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai) jx_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid) jx_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jx_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jx_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jx_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
      }

      jx_MPI_Allreduce(&num_ai_local,&num_ai_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_local_valid,&num_ai_global_valid,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_vars_local,&num_vars_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_c_local,&num_ai_c_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local,&mai_global,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_valid,&mai_global_valid,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_c_local,&mai_c_global,1,JX_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JX_Real num_ai = num_ai_global;
         JX_Real num_ai_valid = num_ai_global_valid;
         JX_Real num_var = num_vars_global;
         JX_Real num_ai_c = num_ai_c_global;
         jx_printf(" \n");
         jx_printf("All levels (except coarsest level): \n");
         jx_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jx_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jx_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai>0) jx_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid>0) jx_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jx_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jx_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jx_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
         jx_printf("================================================================== \n");
      }
   }

   if ((not_finished_coarsening == 0) && (level > 0))
   {
      F_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jx_hpCSRMatrixRowStarts(A_array[level]));
      jx_ParVectorInitialize(F_array[level]);
      jx_ParVectorSetPartitioningOwner(F_array[level], 0);

      U_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jx_hpCSRMatrixRowStarts(A_array[level]));
      jx_ParVectorInitialize(U_array[level]);
      jx_ParVectorSetPartitioningOwner(U_array[level], 0);
   }

   if (not_finished_coarsening == 0)
   {
      jx_ParVectorCopy(F_array[level-1], Vtemp);
      jx_hpCSRMatrixMatvec(-1.0, A_array[level-1], U_array[level-1], 1.0, Vtemp);
      jx_ParCSRMatrixMatvecT(1.0, P_array[level-1], Vtemp, 0.0, F_array[level]);
      jx_PAMGRelax9(jx_hpCSRMatrixPar(A_array[level]), F_array[level], NULL, 0, 1.0, 0.0, U_array[level], NULL);
   }

   for (i = level; i > 0; i --)
   {
      jx_ParCSRMatrixMatvec(1.0, P_array[i-1], U_array[i], 1.0, U_array[i-1]);
      jx_PAMGRelax0(jx_hpCSRMatrixPar(A_array[i-1]), F_array[i-1], NULL, 0, 1.0, 0.0, U_array[i-1], Vtemp);
   }

  /*-----------------------------------------------------------------------
   * enter all the stuff created, A[level], P[level], CF_marker[level],
   * for levels 1 through coarsest, into amg_data data structure
   *-----------------------------------------------------------------------*/

   num_levels = level + 1;
   jx_ParAMGDataNumLevels(amg_data) = num_levels;
   if (jx_ParAMGDataSmoothNumLevels(amg_data) > level)
   {
      jx_ParAMGDataSmoothNumLevels(amg_data) = level;
   }
   smooth_num_levels = jx_ParAMGDataSmoothNumLevels(amg_data);
   
   for (j = 0; j < smooth_num_levels; j++) // Euclid smoothers
   {
      if (smooth_type == 9 || smooth_type == 19)
      {
         JX_EuclidCreate(comm, &smoother[j]);
         if (euclidfile)
         {
            JX_EuclidSetParamsFromFile(smoother[j], euclidfile);
         }
         JX_EuclidSetLevel(smoother[j], eu_level);
         if (eu_bj)
         {
            JX_EuclidSetBJ(smoother[j], eu_bj);
         }
         if (eu_sparse_A)
         {
            JX_EuclidSetSparseA(smoother[j], eu_sparse_A);
         }
         JX_EuclidSetup(smoother[j], (JX_hpCSRMatrix)A_array[j]);
      }
   }

   if (amg_logging > 1) 
   {
      Residual_array = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[0]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[0]),
                                          jx_hpCSRMatrixRowStarts(A_array[0]));
      jx_ParVectorInitialize(Residual_array);
      jx_ParVectorSetPartitioningOwner(Residual_array,0);
      jx_ParAMGDataResidual(amg_data) = Residual_array;
   }
   else
   {
      jx_ParAMGDataResidual(amg_data) = NULL;
   }

   jx_TFree(opt_icor); /* Yue Xiaoqiang 2012/10/12 */
   if (Pmarkers_global_rap)
   {
      jx_TFree(Pmarkers_global_rap); /* Yue Xiaoqiang 2012/10/17 */
   }
   if (Pmarkers_global_interp)
   {
      jx_TFree(Pmarkers_global_interp); /* Yue Xiaoqiang 2012/10/17 */
   }

  /*--------------------------------------------------------------
   *    Print some stuff
   *-------------------------------------------------------------*/
   
   //if (amg_print_level == 1 || amg_print_level == 3)
   if (amg_print_level > 1)
   {
      /* Write the SETUP parameters */
      jx_PAMGSetupStatus(amg_data, jx_hpCSRMatrixPar(hp_matrix));
   }

   if (wall_time_option == 1)
   {
      jx_printf("\n\nProc = %d, Coarsen Time = %f\n", my_id, wall_time_coarsen);
      jx_printf("Proc = %d, Build Coarse Operator Time = %f\n", my_id, wall_time_rap);
      jx_printf("Proc = %d, Build Interp Time = %f\n\n", my_id, wall_time_interp);
   }

   if (my_id == 0 && amg_print_level > 1)
   {
      /* Write the SOLVE parameters */
      jx_PAMGSolveStatus(amg_data); /* Yue Xiaoqiang 2014/04/12 */
   }

   return(Setup_err_flag);
}

JX_Int
JX_ParAdpSetupAMGGetNumSL( JX_Solver solver, JX_Int *num_sl )
{
   return( jx_ParAdpSetupAMGGetNumSL( (void *) solver, num_sl ) );
}

JX_Int
jx_ParAdpSetupAMGGetNumSL( void *data, JX_Int *num_sl )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *num_sl = jx_ParAdpSetupAMGDataNumSL(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetNumReUse( JX_Solver solver, JX_Int *num_reuse )
{
   return( jx_ParAdpSetupAMGGetNumReUse( (void *) solver, num_reuse ) );
}

JX_Int
jx_ParAdpSetupAMGGetNumReUse( void *data, JX_Int *num_reuse )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *num_reuse = jx_ParAdpSetupAMGDataNumReUse(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetCanReUse( JX_Solver solver, JX_Int *can_reuse )
{
   return( jx_ParAdpSetupAMGGetCanReUse( (void *) solver, can_reuse ) );
}

JX_Int
jx_ParAdpSetupAMGGetCanReUse( void *data, JX_Int *can_reuse )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *can_reuse = jx_ParAdpSetupAMGDataCanReUse(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetNumSLAMG( JX_Solver solver, JX_Int *num_sl_amg )
{
   return( jx_ParAdpSetupAMGGetNumSLAMG( (void *) solver, num_sl_amg ) );
}

JX_Int
jx_ParAdpSetupAMGGetNumSLAMG( void *data, JX_Int *num_sl_amg )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *num_sl_amg = jx_ParAdpSetupAMGDataNumSLAMG(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetCurrentSL( JX_Solver solver, JX_Int *iter_sl, JX_Int *converge, JX_Int *num_lvls, JX_Int **iter_sl_lvls, JX_Int *update )
{
   return( jx_ParAdpSetupAMGGetCurrentSL( (void *) solver, iter_sl, converge, num_lvls, iter_sl_lvls, update ) );
}

JX_Int
jx_ParAdpSetupAMGGetCurrentSL( void *data, JX_Int *iter_sl, JX_Int *converge, JX_Int *num_lvls, JX_Int **iter_sl_lvls, JX_Int *update )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *iter_sl = jx_ParAdpSetupAMGDataCurrentIterSL(amg_data);
  *converge = jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data);
  *num_lvls = jx_ParAdpSetupAMGDataCurrentNumLvls(amg_data);
  *iter_sl_lvls = jx_ParAdpSetupAMGDataIterSLLvls(amg_data);
  *update = jx_ParAdpSetupAMGDataCurrentUpdate(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetTotalNumLevels( JX_Solver solver, JX_Int *total_num_levels )
{
   return( jx_ParAdpSetupAMGGetTotalNumLevels( (void *) solver, total_num_levels ) );
}

JX_Int
jx_ParAdpSetupAMGGetTotalNumLevels( void *data, JX_Int *total_num_levels )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *total_num_levels = jx_ParAdpSetupAMGDataTotalNumLevels(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetTotalRedunIterSL( JX_Solver solver, JX_Int *total_redun_iter_sl )
{
   return( jx_ParAdpSetupAMGGetTotalRedunIterSL( (void *) solver, total_redun_iter_sl ) );
}

JX_Int
jx_ParAdpSetupAMGGetTotalRedunIterSL( void *data, JX_Int *total_redun_iter_sl )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *total_redun_iter_sl = jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetTotalConvergeIter( JX_Solver solver, JX_Int *total_converge_iter )
{
   return( jx_ParAdpSetupAMGGetTotalConvergeIter( (void *) solver, total_converge_iter ) );
}

JX_Int
jx_ParAdpSetupAMGGetTotalConvergeIter( void *data, JX_Int *total_converge_iter )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *total_converge_iter = jx_ParAdpSetupAMGDataTotalConvergeIter(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetTotalRedunIterAMG( JX_Solver solver, JX_Int *total_redun_iter_amg )
{
   return( jx_ParAdpSetupAMGGetTotalRedunIterAMG( (void *) solver, total_redun_iter_amg ) );
}

JX_Int
jx_ParAdpSetupAMGGetTotalRedunIterAMG( void *data, JX_Int *total_redun_iter_amg )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *total_redun_iter_amg = jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetSetupTime( JX_Solver solver, JX_Real *setup_time )
{
   return( jx_ParAdpSetupAMGGetSetupTime( (void *) solver, setup_time ) );
}

JX_Int
jx_ParAdpSetupAMGGetSetupTime( void *data, JX_Real *setup_time )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *setup_time = jx_ParAdpSetupAMGDataSetupTime(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetSolveTime( JX_Solver solver, JX_Real *solve_time )
{
   return( jx_ParAdpSetupAMGGetSolveTime( (void *) solver, solve_time ) );
}

JX_Int
jx_ParAdpSetupAMGGetSolveTime( void *data, JX_Real *solve_time )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jx_printf("Warning! ParAdpSetupAMG object empty!\n");
      jx_error_in_arg(1);
      return jx_error_flag;
   }
  *solve_time = jx_ParAdpSetupAMGDataSetupTime(amg_data);
   return jx_error_flag;
}

JX_Int
JX_ParAdpSetupAMGGetStatus( JX_Solver solver )
{
   return( jx_ParAdpSetupAMGGetStatus( (void *) solver ) );
}

JX_Int
jx_ParAdpSetupAMGGetStatus( void *data )
{
   jx_ParAdpSetupAMGData *amg_data = data;
   MPI_Comm comm = jx_ParAdpSetupAMGDataComm(amg_data);
   JX_Bool is_adaptive = jx_ParAdpSetupAMGDataIsAdaptive(amg_data);
   JX_Bool is_adap_setup = jx_ParAdpSetupAMGDataIsAdapSetup(amg_data);
   JX_Real tol = jx_ParAdpSetupAMGDataTol(amg_data);
   JX_Real max_row_sum = jx_ParAdpSetupAMGDataMaxRowSum(amg_data);
   JX_Real strong_threshold = jx_ParAdpSetupAMGDataStrongThreshold(amg_data);
   JX_Int max_levels = jx_ParAdpSetupAMGDataMaxLevels(amg_data);
   JX_Int cycle_type = jx_ParAdpSetupAMGDataCycleType(amg_data);
   JX_Int relax_type = jx_ParAdpSetupAMGDataRelaxType(amg_data);
   JX_Int relax_order = jx_ParAdpSetupAMGDataRelaxOrder(amg_data);
   JX_Int ns_down = jx_ParAdpSetupAMGDataNsDown(amg_data);
   JX_Int ns_up = jx_ParAdpSetupAMGDataNsUp(amg_data);
   JX_Int ns_coarse = jx_ParAdpSetupAMGDataNsCoarse(amg_data);
   JX_Int coarsen_type = jx_ParAdpSetupAMGDataCoarsenType(amg_data);
   JX_Int interp_type = jx_ParAdpSetupAMGDataInterpType(amg_data);
   JX_Int P_max_elmts = jx_ParAdpSetupAMGDataPMaxElmts(amg_data);
   JX_Int agg_num_levels = jx_ParAdpSetupAMGDataAggNumLevels(amg_data);
   JX_Int coarse_threshold = jx_ParAdpSetupAMGDataCoarseThreshold(amg_data);
   JX_Int amg_print_level = jx_ParAdpSetupAMGDataAMGPrintLevel(amg_data);
   JX_Int k_dim = jx_ParAdpSetupAMGDataKDim(amg_data);
   JX_Int TTest = jx_ParAdpSetupAMGDataTTest(amg_data);
   JX_Int max_iter = jx_ParAdpSetupAMGDataMaxIter(amg_data);
   JX_Int two_norm = jx_ParAdpSetupAMGDataTwoNorm(amg_data);
   JX_Int solver_id = jx_ParAdpSetupAMGDataSolverID(amg_data);
   JX_Int print_level = jx_ParAdpSetupAMGDataPrintLevel(amg_data);
   JX_Int max_iter_sl = jx_ParAdpSetupAMGDataMaxIterSL(amg_data);
   JX_Int max_iter_reuse = jx_ParAdpSetupAMGDataMaxIterReUse(amg_data);
   JX_Int is_check_restarted = jx_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JX_Int my_id;

   jx_MPI_Comm_rank(comm, &my_id);
   if (my_id == 0)
   {
      if (is_adaptive) jx_printf(" is_adaptive = TRUE\n");
      else jx_printf(" is_adaptive = FALSE\n");
      if (is_adap_setup) jx_printf(" is_adap_setup = TRUE\n");
      else jx_printf(" is_adap_setup = FALSE\n");
      jx_printf(" tol = %lf\n", tol);
      jx_printf(" max_row_sum = %lf\n", max_row_sum);
      jx_printf(" strong_threshold = %lf\n", strong_threshold);
      jx_printf(" max_levels = %d\n", max_levels);
      jx_printf(" cycle_type = %d\n", cycle_type);
      jx_printf(" relax_type = %d\n", relax_type);
      jx_printf(" relax_order = %d\n", relax_order);
      jx_printf(" ns_down, up, coarse = %d, %d, %d\n", ns_down, ns_up, ns_coarse);
      jx_printf(" coarsen_type = %d\n", coarsen_type);
      jx_printf(" interp_type = %d\n", interp_type);
      jx_printf(" P_max_elmts = %d\n", P_max_elmts);
      jx_printf(" agg_num_levels = %d\n", agg_num_levels);
      jx_printf(" coarse_threshold = %d\n", coarse_threshold);
      jx_printf(" amg_print_level = %d\n", amg_print_level);
      jx_printf(" k_dim = %d\n", k_dim);
      jx_printf(" TTest = %d\n", TTest);
      jx_printf(" max_iter = %d\n", max_iter);
      jx_printf(" two_norm = %d\n", two_norm);
      jx_printf(" solver_id = %d\n", solver_id);
      jx_printf(" print_level = %d\n", print_level);
      jx_printf(" max_iter_sl = %d\n", max_iter_sl);
      jx_printf(" max_iter_reuse = %d\n", max_iter_reuse);
      jx_printf(" is_check_restarted = %d\n", is_check_restarted);
   }

   return 0;
}
