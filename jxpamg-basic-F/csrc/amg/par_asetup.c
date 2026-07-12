//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_asetup.c
 *  Date: 2020/03/10
 */ 

#include "jxf_asetup.h"

extern JXF_Int *Jxf_Pmarkers_global_rap;    /* Yue Xiaoqiang 2012/10/17 */
extern JXF_Int *Jxf_Pmarkers_global_interp; /* Yue Xiaoqiang 2012/10/17 */

/*!
 * \fn JXF_Int JXF_ParAdpSetupAMGCreate
 * \date 2020/03/10
 */
JXF_Int
JXF_ParAdpSetupAMGCreate( MPI_Comm comm, JXF_Solver *solver )
{
  *solver = (JXF_Solver) jxf_ParAdpSetupAMGCreate( comm );
   if (!solver)
   {
      jxf_error_in_arg(1);
   }
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_ParAdpSetupAMGDestroy
 * \date 2020/03/10
 */
JXF_Int
JXF_ParAdpSetupAMGDestroy( JXF_Solver solver )
{
   return( jxf_ParAdpSetupAMGDestroy( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_ParAdpSetupAMGSetup
 * \date 2020/03/10
 */
JXF_Int
JXF_ParAdpSetupAMGSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_ParAdpSetupAMGSetup( solver, matA, vecB, vecX ) );
}

/*!
 * \fn JXF_Int JXF_ParAdpSetupAMGSolve
 * \date 2020/03/10
 */
JXF_Int
JXF_ParAdpSetupAMGSolve( JXF_Solver solver, JXF_Matrix preA, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_ParAdpSetupAMGSolve( (void *) solver, (void *) preA, (void *) matA, (void *) vecB, (void *) vecX ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetIsAdaptive( JXF_Solver solver, JXF_Bool is_adaptive )
{
   return( jxf_ParAdpSetupAMGSetIsAdaptive( (void *) solver, is_adaptive ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetIsAdapSetup( JXF_Solver solver, JXF_Bool is_adap_setup )
{
   return( jxf_ParAdpSetupAMGSetIsAdapSetup( (void *) solver, is_adap_setup ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_ParAdpSetupAMGSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetMaxRowSum( JXF_Solver solver, JXF_Real max_row_sum )
{
   return( jxf_ParAdpSetupAMGSetMaxRowSum( (void *) solver, max_row_sum ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold )
{
   return( jxf_ParAdpSetupAMGSetStrongThreshold( (void *) solver, strong_threshold ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetMaxLevels( JXF_Solver solver, JXF_Int  max_levels )
{
   return( jxf_ParAdpSetupAMGSetMaxLevels( (void *) solver, max_levels ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetCycleType( JXF_Solver solver,JXF_Int cycle_type )
{
   return( jxf_ParAdpSetupAMGSetCycleType( (void *) solver, cycle_type ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetRelaxType( JXF_Solver solver, JXF_Int relax_type )
{
   return( jxf_ParAdpSetupAMGSetRelaxType( (void *) solver, relax_type ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetRelaxOrder( JXF_Solver solver, JXF_Int relax_order )
{
   return( jxf_ParAdpSetupAMGSetRelaxOrder( (void *) solver, relax_order ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetNsDown( JXF_Solver solver, JXF_Int ns_down )
{
   return( jxf_ParAdpSetupAMGSetNsDown( (void *) solver, ns_down ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetNsUp( JXF_Solver solver, JXF_Int ns_up )
{
   return( jxf_ParAdpSetupAMGSetNsUp( (void *) solver, ns_up ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetNsCoarse( JXF_Solver solver, JXF_Int ns_coarse )
{
   return( jxf_ParAdpSetupAMGSetNsCoarse( (void *) solver, ns_coarse ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type )
{
   return( jxf_ParAdpSetupAMGSetCoarsenType( (void *) solver, coarsen_type ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetInterpType( JXF_Solver solver, JXF_Int interp_type )
{
   return( jxf_ParAdpSetupAMGSetInterpType( (void *) solver, interp_type ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetPMaxElmts( JXF_Solver solver, JXF_Int P_max_elmts )
{
   return( jxf_ParAdpSetupAMGSetPMaxElmts( (void *) solver, P_max_elmts ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetAggNumLevels( JXF_Solver solver, JXF_Int agg_num_levels )
{
   return( jxf_ParAdpSetupAMGSetAggNumLevels( (void *) solver, agg_num_levels ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetCoarseThreshold( JXF_Solver solver, JXF_Int coarse_threshold )
{
   return (jxf_ParAdpSetupAMGSetCoarseThreshold ( (void *) solver, coarse_threshold ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetAMGPrintLevel( JXF_Solver solver, JXF_Int amg_print_level )
{
   return( jxf_ParAdpSetupAMGSetAMGPrintLevel( (void *) solver, amg_print_level ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetKDim( JXF_Solver solver, JXF_Int k_dim )
{
   return( jxf_ParAdpSetupAMGSetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetTTest( JXF_Solver solver, JXF_Int TTest )
{
   return( jxf_ParAdpSetupAMGSetTTest( (void *) solver, TTest ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_ParAdpSetupAMGSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetTwoNorm( JXF_Solver solver, JXF_Int two_norm )
{
   return( jxf_ParAdpSetupAMGSetTwoNorm( (void *) solver, two_norm ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetSolverID( JXF_Solver solver, JXF_Int solver_id )
{
   return( jxf_ParAdpSetupAMGSetSolverID( (void *) solver, solver_id ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetPrintLevel( JXF_Solver solver, JXF_Int print_level )
{
   return( jxf_ParAdpSetupAMGSetPrintLevel( (void *) solver, print_level ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetMaxIterSL( JXF_Solver solver, JXF_Int max_iter_sl )
{
   return( jxf_ParAdpSetupAMGSetMaxIterSL( (void *) solver, max_iter_sl ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetMaxIterReUse( JXF_Solver solver, JXF_Int max_iter_reuse )
{
   return( jxf_ParAdpSetupAMGSetMaxIterReUse( (void *) solver, max_iter_reuse ) );
}

JXF_Int
JXF_ParAdpSetupAMGSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted )
{
   return( jxf_ParAdpSetupAMGSetIsCheckRestarted( (void *) solver, is_check_restarted ) );
}

/*!
 * \fn void *jxf_ParAdpSetupAMGCreate
 * \date 2020/03/10
 */
void *
jxf_ParAdpSetupAMGCreate( MPI_Comm comm )
{
   jxf_ParAdpSetupAMGData *amg_data = jxf_CTAlloc(jxf_ParAdpSetupAMGData, 1);

   JXF_Bool is_adaptive = JXF_FALSE;
   JXF_Bool is_adap_setup = JXF_FALSE;

   JXF_Real tol = 1.0e-7;
   JXF_Real max_row_sum = 0.9;
   JXF_Real strong_threshold = 0.25;

   JXF_Int max_levels = 25;
   JXF_Int cycle_type = 1;
   JXF_Int relax_type = 3;
   JXF_Int relax_order = 1;
   JXF_Int ns_down = 1;
   JXF_Int ns_up = 1;
   JXF_Int ns_coarse = 1;
   JXF_Int coarsen_type = 6;
   JXF_Int interp_type = 0;
   JXF_Int P_max_elmts = 0;
   JXF_Int agg_num_levels = 0;
   JXF_Int coarse_threshold = 9;
   JXF_Int amg_print_level = 0;

   JXF_Int k_dim = 5;
   JXF_Int TTest = 0;
   JXF_Int max_iter = 200;
   JXF_Int two_norm = 0;
   JXF_Int solver_id = 2;
   JXF_Int print_level = 0;
   JXF_Int max_iter_sl = 1000;
   JXF_Int max_iter_reuse = 1000;
   JXF_Int is_check_restarted = 0;

   jxf_ParAdpSetupAMGDataComm(amg_data) = comm;
   jxf_ParAdpSetupAMGDataIterSLLvls(amg_data) = jxf_CTAlloc(JXF_Int, 25);

   jxf_ParAdpSetupAMGSetIsAdaptive(amg_data, is_adaptive);
   jxf_ParAdpSetupAMGSetIsAdapSetup(amg_data, is_adap_setup);

   jxf_ParAdpSetupAMGSetTol(amg_data, tol);
   jxf_ParAdpSetupAMGSetMaxRowSum(amg_data, max_row_sum);
   jxf_ParAdpSetupAMGSetStrongThreshold(amg_data, strong_threshold);

   jxf_ParAdpSetupAMGSetMaxLevels(amg_data, max_levels);
   jxf_ParAdpSetupAMGSetCycleType(amg_data, cycle_type);
   jxf_ParAdpSetupAMGSetRelaxType(amg_data, relax_type);
   jxf_ParAdpSetupAMGSetRelaxOrder(amg_data, relax_order);
   jxf_ParAdpSetupAMGSetNsDown(amg_data, ns_down);
   jxf_ParAdpSetupAMGSetNsUp(amg_data, ns_up);
   jxf_ParAdpSetupAMGSetNsCoarse(amg_data, ns_coarse);
   jxf_ParAdpSetupAMGSetCoarsenType(amg_data, coarsen_type);
   jxf_ParAdpSetupAMGSetInterpType(amg_data, interp_type);
   jxf_ParAdpSetupAMGSetPMaxElmts(amg_data, P_max_elmts);
   jxf_ParAdpSetupAMGSetAggNumLevels(amg_data, agg_num_levels);
   jxf_ParAdpSetupAMGSetCoarseThreshold(amg_data, coarse_threshold);
   jxf_ParAdpSetupAMGSetAMGPrintLevel(amg_data, amg_print_level);

   jxf_ParAdpSetupAMGSetKDim(amg_data, k_dim);
   jxf_ParAdpSetupAMGSetTTest(amg_data, TTest);
   jxf_ParAdpSetupAMGSetMaxIter(amg_data, max_iter);
   jxf_ParAdpSetupAMGSetTwoNorm(amg_data, two_norm);
   jxf_ParAdpSetupAMGSetSolverID(amg_data, solver_id);
   jxf_ParAdpSetupAMGSetPrintLevel(amg_data, print_level);
   jxf_ParAdpSetupAMGSetMaxIterSL(amg_data, max_iter_sl);
   jxf_ParAdpSetupAMGSetMaxIterReUse(amg_data, max_iter_reuse);
   jxf_ParAdpSetupAMGSetIsCheckRestarted(amg_data, is_check_restarted);

   jxf_ParAdpSetupAMGDataNumSL(amg_data) = 0;
   jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 0;
   jxf_ParAdpSetupAMGDataNumReUse(amg_data) = 0;
   jxf_ParAdpSetupAMGDataNumSLAMG(amg_data) = 0;
   jxf_ParAdpSetupAMGDataSetupTime(amg_data) = 0.0;
   jxf_ParAdpSetupAMGDataSolveTime(amg_data) = 0.0;
   jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
   jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) = 0;
   jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) = 0;
   jxf_ParAdpSetupAMGDataTotalConvergeIter(amg_data) = 0;
   jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) = 0;

   return (void *) amg_data;
}

/*!
 * \fn JXF_Int jxf_ParAdpSetupAMGDestroy
 * \date 2020/03/10
 */
JXF_Int
jxf_ParAdpSetupAMGDestroy( void *data )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (jxf_ParAdpSetupAMGDataAMGSolver(amg_data))
   {
      JXF_PAMGDestroy(jxf_ParAdpSetupAMGDataAMGSolver(amg_data));
   }
   jxf_TFree(jxf_ParAdpSetupAMGDataIterSLLvls(amg_data));
   jxf_TFree(amg_data);
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetup( void *data, void *matA, void *vecB, void *vecX )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   MPI_Comm comm = jxf_ParAdpSetupAMGDataComm(amg_data);
   JXF_Real tol = jxf_ParAdpSetupAMGDataTol(amg_data);
   JXF_Int k_dim = jxf_ParAdpSetupAMGDataKDim(amg_data);
   JXF_Int TTest = jxf_ParAdpSetupAMGDataTTest(amg_data);
   JXF_Int max_iter = jxf_ParAdpSetupAMGDataMaxIter(amg_data);
   JXF_Int two_norm = jxf_ParAdpSetupAMGDataTwoNorm(amg_data);
   JXF_Int solver_id = jxf_ParAdpSetupAMGDataSolverID(amg_data);
   JXF_Int print_level = jxf_ParAdpSetupAMGDataPrintLevel(amg_data);
   JXF_Int is_check_restarted = jxf_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JXF_Solver solver;
   JXF_Solver jac_solver;
   JXF_Real starttime, endtime;

   if (TTest) starttime = jxf_MPI_Wtime();
   JXF_JacobiCreate(&jac_solver);
   JXF_JacobiSetMaxIter(jac_solver, 1);
   JXF_JacobiSetup(jac_solver, (JXF_hpCSRMatrix) matA);
   jxf_ParAdpSetupAMGDataJACSolver(amg_data) = jac_solver;
   if (solver_id == 1)
   {
      JXF_PCGCreate(comm, &solver);
      JXF_PCGSetMaxIter(solver, max_iter);
      JXF_PCGSetTol(solver, tol);
      JXF_PCGSetTwoNorm(solver, two_norm);
      JXF_PCGSetLogging(solver, 1);
      JXF_PCGSetPrintLevel(solver, print_level);
      JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_PCGSetup(solver, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
   }
   else if (solver_id == 2)
   {
      JXF_GMRESCreate(comm, &solver);
      JXF_GMRESSetMaxIter(solver, max_iter);
      JXF_GMRESSetKDim(solver, k_dim);
      JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted);
      JXF_GMRESSetTol(solver, tol);
      JXF_GMRESSetLogging(solver, 1);
      JXF_GMRESSetPrintLevel(solver, print_level);
      JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_GMRESSetup(solver, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
   }
   else if (solver_id == 3)
   {
      JXF_BiCGSTABCreate(comm, &solver);
      JXF_BiCGSTABSetMaxIter(solver, max_iter);
      JXF_BiCGSTABSetTol(solver, tol);
      JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JXF_BiCGSTABSetConvCriteria(solver, 0);
      JXF_BiCGSTABSetLogging(solver, 1);
      JXF_BiCGSTABSetPrintLevel(solver, print_level);
      JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_BiCGSTABSetup(solver, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
   }
   jxf_ParAdpSetupAMGDataSolver(amg_data) = solver;
   if (TTest)
   {
      endtime = jxf_MPI_Wtime();
      jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtime - starttime);
   }

   return 0;
}

JXF_Int
jxf_ParAdpSetupAMGSolve( void *data, void *preA, void *matA, void *vecB, void *vecX )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   JXF_Bool is_adaptive = jxf_ParAdpSetupAMGDataIsAdaptive(amg_data);
   JXF_Bool is_adap_setup = jxf_ParAdpSetupAMGDataIsAdapSetup(amg_data);
   JXF_Real tol = jxf_ParAdpSetupAMGDataTol(amg_data);
   JXF_Real max_row_sum = jxf_ParAdpSetupAMGDataMaxRowSum(amg_data);
   JXF_Real strong_threshold = jxf_ParAdpSetupAMGDataStrongThreshold(amg_data);
   JXF_Int k_dim = jxf_ParAdpSetupAMGDataKDim(amg_data);
   JXF_Int max_iter = jxf_ParAdpSetupAMGDataMaxIter(amg_data);
   JXF_Int two_norm = jxf_ParAdpSetupAMGDataTwoNorm(amg_data);
   JXF_Int is_check_restarted = jxf_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JXF_Int TTest = jxf_ParAdpSetupAMGDataTTest(amg_data);
   JXF_Int solver_id = jxf_ParAdpSetupAMGDataSolverID(amg_data);
   JXF_Int max_iter_sl = jxf_ParAdpSetupAMGDataMaxIterSL(amg_data);
   JXF_Int max_iter_reuse = jxf_ParAdpSetupAMGDataMaxIterReUse(amg_data);
   JXF_Solver solver = jxf_ParAdpSetupAMGDataSolver(amg_data);
   JXF_Solver amg_solver = jxf_ParAdpSetupAMGDataAMGSolver(amg_data);
   JXF_Int max_levels = jxf_ParAdpSetupAMGDataMaxLevels(amg_data);
   JXF_Int cycle_type = jxf_ParAdpSetupAMGDataCycleType(amg_data);
   JXF_Int relax_type = jxf_ParAdpSetupAMGDataRelaxType(amg_data);
   JXF_Int relax_order = jxf_ParAdpSetupAMGDataRelaxOrder(amg_data);
   JXF_Int ns_down = jxf_ParAdpSetupAMGDataNsDown(amg_data);
   JXF_Int ns_up = jxf_ParAdpSetupAMGDataNsUp(amg_data);
   JXF_Int ns_coarse = jxf_ParAdpSetupAMGDataNsCoarse(amg_data);
   JXF_Int coarsen_type = jxf_ParAdpSetupAMGDataCoarsenType(amg_data);
   JXF_Int interp_type = jxf_ParAdpSetupAMGDataInterpType(amg_data);
   JXF_Int P_max_elmts = jxf_ParAdpSetupAMGDataPMaxElmts(amg_data);
   JXF_Int agg_num_levels = jxf_ParAdpSetupAMGDataAggNumLevels(amg_data);
   JXF_Int coarse_threshold = jxf_ParAdpSetupAMGDataCoarseThreshold(amg_data);
   JXF_Int amg_print_level = jxf_ParAdpSetupAMGDataAMGPrintLevel(amg_data);
   JXF_Int *iter_sl_lvls = jxf_ParAdpSetupAMGDataIterSLLvls(amg_data);
   JXF_Bool ok_ok = JXF_TRUE;
   JXF_Int iter_sl_vvv = 0;
   JXF_Int iter_sl_used = 0;
   JXF_Int iter_amg_used = 0;
   JXF_Real iter_rest;
   JXF_Real res_norm;
   JXF_Real res_norm_old;
   JXF_Real reduce_factor;
   JXF_Real reduce_factor_thrd = 1.0;
   JXF_Real starttime, endtime;
   JXF_Real starttimeT = 0.0, endtimeT = 0.0;
   JXF_Int iter_sl_oused;
   JXF_Int num_levels;
   JXF_Int num_sl_amg;
   JXF_Int num_iterations;

   if (TTest) starttime = jxf_MPI_Wtime();
   if (solver_id == 1)
   {
      JXF_PCGSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_PCGSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
         JXF_PCGGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jxf_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JXF_JacobiDestroy(jxf_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JXF_FALSE;
         JXF_PCGDestroy(solver);
         jxf_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jxf_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) NULL, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JXF_PCGSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
            JXF_PCGGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JXF_FALSE;
            jxf_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JXF_PCGDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) NULL, amg_solver);
         JXF_PAMGAdapSetup(amg_solver, (JXF_hpCSRMatrix) matA, (JXF_ParVector) vecB, (JXF_ParVector) vecX, 1, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jxf_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_PCGSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGDestroy(solver);
      }
      else if (ok_ok)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) matA);
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_PCGSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGDestroy(solver);
      }
   }
   else if (solver_id == 2)
   {
      JXF_GMRESSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_GMRESSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
         JXF_GMRESGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jxf_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JXF_JacobiDestroy(jxf_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JXF_FALSE;
         JXF_GMRESDestroy(solver);
         jxf_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jxf_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JXF_GMRESSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
            JXF_GMRESGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JXF_FALSE;
            jxf_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JXF_GMRESDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_GMRESSetMaxIter(solver, max_iter);
         JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) NULL, amg_solver);
         JXF_PAMGAdapSetup(amg_solver, (JXF_hpCSRMatrix) matA, (JXF_ParVector) vecB, (JXF_ParVector) vecX, 2, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jxf_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_GMRESSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_GMRESGetNumIterations(solver, &num_iterations);
         JXF_GMRESDestroy(solver);
      }
      else if (ok_ok)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_GMRESSetMaxIter(solver, max_iter);
         JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) matA);
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_GMRESSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_GMRESGetNumIterations(solver, &num_iterations);
         JXF_GMRESDestroy(solver);
      }
   }
   else if (solver_id == 3)
   {
      JXF_BiCGSTABSetMaxIter(solver, 1);
      res_norm = 1.0;
      do {
         if (is_adaptive)
         {
            iter_rest = max_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(tol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_BiCGSTABSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
         JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &res_norm);
         iter_sl_vvv ++;
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

      jxf_ParAdpSetupAMGDataCurrentIterSL(amg_data) = iter_sl_vvv;
      JXF_JacobiDestroy(jxf_ParAdpSetupAMGDataJACSolver(amg_data));
      if (res_norm > tol)
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 0;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_used;
      }
      else
      {
         jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) = 1;
         num_iterations = iter_sl_used;
         ok_ok = JXF_FALSE;
         JXF_BiCGSTABDestroy(solver);
         jxf_ParAdpSetupAMGDataNumSL(amg_data) ++;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) ++;
      }

      if (ok_ok && jxf_ParAdpSetupAMGDataCanReUse(amg_data))
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 0;
         JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         res_norm = 1.0;
         do {
            if (is_adaptive)
            {
               iter_rest = max_iter_reuse - iter_amg_used;
               iter_rest = 1.0 / iter_rest;
               reduce_factor_thrd = pow(tol/res_norm, iter_rest);
            }
            res_norm_old = res_norm;
            JXF_BiCGSTABSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector)vecX);
            JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &res_norm);
            iter_amg_used ++;
            reduce_factor = res_norm / res_norm_old;
         } while ((res_norm > tol) && (reduce_factor < reduce_factor_thrd));

         if (res_norm > tol)
         {
            jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) += iter_amg_used;
         }
         else
         {
            num_iterations = iter_amg_used;
            ok_ok = JXF_FALSE;
            jxf_ParAdpSetupAMGDataNumReUse(amg_data) ++;
            JXF_BiCGSTABDestroy(solver);
         }
      }

      if (ok_ok && is_adap_setup)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_BiCGSTABSetMaxIter(solver, max_iter);
         JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) NULL, amg_solver);
         JXF_PAMGAdapSetup(amg_solver, (JXF_hpCSRMatrix) matA, (JXF_ParVector) vecB, (JXF_ParVector) vecX, 3, amg_print_level,
                          tol, k_dim, two_norm, is_check_restarted, is_adaptive, max_iter_sl, &iter_sl_oused, &num_sl_amg, iter_sl_lvls);
         jxf_ParAdpSetupAMGDataNumSLAMG(amg_data) += num_sl_amg;
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data) += iter_sl_oused;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_BiCGSTABSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataCurrentNumLvls(amg_data) = num_levels;
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
         JXF_BiCGSTABDestroy(solver);
      }
      else if (ok_ok)
      {
         jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data) = 1;
         if (TTest) starttimeT = jxf_MPI_Wtime();
         if (jxf_ParAdpSetupAMGDataCanReUse(amg_data)) jxf_hpCSRMatrixDestroy((jxf_hpCSRMatrix *)jxf_ParAdpSetupAMGDataPreMat(amg_data));
         if (amg_solver) JXF_PAMGDestroy(amg_solver);
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxIter(amg_solver, 1);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetCycleType(amg_solver, cycle_type);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
         JXF_PAMGSetRelaxOrder(amg_solver, relax_order);
         if (ns_down > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);
         if (ns_up > -1) JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);
         JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JXF_PAMGSetInterpType(amg_solver, interp_type);
         JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_BiCGSTABSetMaxIter(solver, max_iter);
         JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond, (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) matA);
         jxf_ParAdpSetupAMGDataPreMat(amg_data) = (JXF_hpCSRMatrix) matA;
         jxf_ParAdpSetupAMGDataCanReUse(amg_data) = 1;
         jxf_ParAdpSetupAMGDataAMGSolver(amg_data) = amg_solver;
         if (TTest)
         {
            endtimeT = jxf_MPI_Wtime();
            jxf_ParAdpSetupAMGDataSetupTime(amg_data) += (endtimeT - starttimeT);
         }
         JXF_BiCGSTABSolve(solver, (JXF_Matrix) preA, (JXF_Matrix) matA, (JXF_Vector) vecB, (JXF_Vector) vecX);
         JXF_PAMGGetNumLevels(amg_solver, &num_levels);
         jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data) += num_levels;
         JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
         JXF_BiCGSTABDestroy(solver);
      }
   }
   jxf_ParAdpSetupAMGDataTotalConvergeIter(amg_data) += num_iterations;
   if (TTest)
   {
      endtime = jxf_MPI_Wtime();
      jxf_ParAdpSetupAMGDataSolveTime(amg_data) += ((endtime - starttime) - (endtimeT - starttimeT));
   }

   return 0;
}

JXF_Int
jxf_ParAdpSetupAMGSetIsAdaptive( void *data, JXF_Bool is_adaptive )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataIsAdaptive(amg_data) = is_adaptive;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetIsAdapSetup( void *data, JXF_Bool is_adap_setup )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataIsAdapSetup(amg_data) = is_adap_setup;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetTol( void *data, JXF_Real tol )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (tol < 0 || tol > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataTol(amg_data) = tol;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetMaxRowSum( void *data, JXF_Real max_row_sum )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (max_row_sum <= 0 || max_row_sum > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataMaxRowSum(amg_data) = max_row_sum;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetStrongThreshold( void *data, JXF_Real strong_threshold )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (strong_threshold < 0 || strong_threshold > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataStrongThreshold(amg_data) = strong_threshold;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetMaxLevels( void *data, JXF_Int max_levels )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (max_levels < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataMaxLevels(amg_data) = max_levels;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetCycleType( void *data, JXF_Int cycle_type )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (cycle_type < 0 || cycle_type > 2)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataCycleType(amg_data) = cycle_type;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetRelaxType( void *data, JXF_Int relax_type )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (relax_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataRelaxType(amg_data) = relax_type;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetRelaxOrder( void *data, JXF_Int relax_order )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataRelaxOrder(amg_data) = relax_order;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetNsDown( void *data, JXF_Int ns_down )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataNsDown(amg_data) = ns_down;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetNsUp( void *data, JXF_Int ns_up )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataNsUp(amg_data) = ns_up;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetNsCoarse( void *data, JXF_Int ns_coarse )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataNsCoarse(amg_data) = ns_coarse;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetCoarsenType( void *data, JXF_Int coarsen_type )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataCoarsenType(amg_data) = coarsen_type;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetInterpType( void *data, JXF_Int interp_type )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (interp_type < 0 || interp_type > 100)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataInterpType(amg_data) = interp_type;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetPMaxElmts( void *data, JXF_Int P_max_elmts )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (P_max_elmts < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataPMaxElmts(amg_data) = P_max_elmts;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetAggNumLevels( void *data, JXF_Int agg_num_levels )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (agg_num_levels < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataAggNumLevels(amg_data) = agg_num_levels;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetCoarseThreshold( void *data, JXF_Int coarse_threshold )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataCoarseThreshold(amg_data) = coarse_threshold;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetAMGPrintLevel( void *data, JXF_Int amg_print_level )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataAMGPrintLevel(amg_data) = amg_print_level;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetKDim( void *data, JXF_Int k_dim )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataKDim(amg_data) = k_dim;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetTTest( void *data, JXF_Int TTest )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataTTest(amg_data) = TTest;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetMaxIter( void *data, JXF_Int max_iter )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataMaxIter(amg_data) = max_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetTwoNorm( void *data, JXF_Int two_norm )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataTwoNorm(amg_data) = two_norm;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetSolverID( void *data, JXF_Int solver_id )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataSolverID(amg_data) = solver_id;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetPrintLevel( void *data, JXF_Int print_level )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataPrintLevel(amg_data) = print_level;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetMaxIterSL( void *data, JXF_Int max_iter_sl )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataMaxIterSL(amg_data) = max_iter_sl;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetMaxIterReUse( void *data, JXF_Int max_iter_reuse )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataMaxIterReUse(amg_data) = max_iter_reuse;
   return jxf_error_flag;
}

JXF_Int
jxf_ParAdpSetupAMGSetIsCheckRestarted( void *data, JXF_Int is_check_restarted )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAdpSetupAMGDataIsCheckRestarted(amg_data) = is_check_restarted;
   return jxf_error_flag;
}

JXF_Int
JXF_JacobiCreate( JXF_Solver *solver )
{
  *solver = (JXF_Solver) jxf_JacobiCreate();
   if (!solver)
   {
      jxf_error_in_arg(1);
   }
   return jxf_error_flag;
}

void *
jxf_JacobiCreate()
{
   jxf_JacobiData *amg_data = jxf_CTAlloc(jxf_JacobiData, 1);
   jxf_JacobiDataMaxIter(amg_data) = 100;
   return (void *) amg_data;
}

JXF_Int
JXF_JacobiDestroy( JXF_Solver solver )
{
   return( jxf_JacobiDestroy( (void *) solver ) );
}

JXF_Int
jxf_JacobiDestroy( void *data )
{
   jxf_JacobiData *amg_data = data;
   if (jxf_JacobiDataTmpVec(amg_data))
   {
      jxf_ParVectorDestroy(jxf_JacobiDataTmpVec(amg_data));
   }
   jxf_TFree(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_JacobiSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_JacobiSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
jxf_JacobiSetMaxIter( void *data, JXF_Int max_iter )
{
   jxf_JacobiData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! Jacobi object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (max_iter < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_JacobiDataMaxIter(amg_data) = max_iter;
   return jxf_error_flag;
}

JXF_Int 
JXF_JacobiSetup( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix )
{
   return( jxf_JacobiSetup( (void *) solver, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) hp_matrix )) );
}

JXF_Int
jxf_JacobiSetup( void *data, jxf_ParCSRMatrix *par_matrix )
{
   jxf_JacobiData *amg_data = data;
   jxf_JacobiDataTmpVec(amg_data) = jxf_ParVectorCreate(jxf_ParCSRMatrixComm(par_matrix),
                                                      jxf_ParCSRMatrixGlobalNumRows(par_matrix),
                                                      jxf_ParCSRMatrixRowStarts(par_matrix));
   jxf_ParVectorInitialize(jxf_JacobiDataTmpVec(amg_data));
   jxf_ParVectorSetPartitioningOwner(jxf_JacobiDataTmpVec(amg_data), 0);
   return 0;
}

JXF_Int
JXF_JacobiPrecond( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix, JXF_ParVector par_rhs, JXF_ParVector par_app )
{
   return( jxf_JacobiPrecond( (void *) solver, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) hp_matrix), (jxf_ParVector *) par_rhs, (jxf_ParVector *) par_app ) );
}

JXF_Int
jxf_JacobiPrecond( void *data, jxf_ParCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app )
{
   jxf_JacobiData *amg_data = data;
   jxf_ParVector *tmp_vec = jxf_JacobiDataTmpVec(amg_data);
   JXF_Int max_iter = jxf_JacobiDataMaxIter(amg_data);
   JXF_Int cycle_count = 0;
   while (cycle_count < max_iter)
   {
      jxf_PAMGRelax0(par_matrix, par_rhs, NULL, 0, 1.0, 0.0, par_app, tmp_vec);
      cycle_count ++;
   }
   return 0;
}

JXF_Int
JXF_PAMGAdapSetup( JXF_Solver solver, JXF_hpCSRMatrix matA, JXF_ParVector vecB, JXF_ParVector vecX,
                  JXF_Int sid, JXF_Int pl, JXF_Real tol, JXF_Int k_dim, JXF_Int two_norm, JXF_Int is_check_restarted,
                  JXF_Bool is_adaptive, JXF_Int max_iter_sl, JXF_Int *iter_sl_used, JXF_Int *num_sl, JXF_Int *iter_sl_lvls )
{
   return( jxf_PAMGAdapSetup( (void *) solver, (jxf_hpCSRMatrix *) matA, (jxf_ParVector *) vecB, (jxf_ParVector *) vecX,
                                      sid, pl, tol, k_dim, two_norm, is_check_restarted,
                                      is_adaptive, max_iter_sl, iter_sl_used, num_sl, iter_sl_lvls ) );
}

JXF_Int
jxf_PAMGAdapSetup( void *amg_vdata, jxf_hpCSRMatrix *hp_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app,
                  JXF_Int solver_id, JXF_Int yprint_level, JXF_Real ytol, JXF_Int yk_dim, JXF_Int ytwo_norm, JXF_Int yis_check_restarted,
                  JXF_Bool yis_adaptive, JXF_Int ymax_iter_sl, JXF_Int *yiter_sl_used, JXF_Int *ynum_sl, JXF_Int *iter_sl_lvls )
{
   MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   jxf_ParAMGData *amg_data = amg_vdata;

   /* Data Structure variables */
   jxf_hpCSRMatrix **A_array;
   jxf_ParVector    **F_array;
   jxf_ParVector    **U_array;
   jxf_ParVector     *Vtemp;
   jxf_ParCSRMatrix **P_array;
   jxf_ParCSRMatrix **R_array;
   jxf_ParVector     *Residual_array;
   JXF_Real          **AI_measure_array;
   JXF_Int             **CF_marker_array;
   JXF_Int             **dof_func_array;
   JXF_Int              *dof_func;
   JXF_Int              *col_offd_S_to_A;
   JXF_Int              *col_offd_Sabs_to_A = NULL;
   JXF_Int            *AIR_maxsize_ls;
   JXF_Real            strong_threshold;
   JXF_Real            AIR_strong_th;
   JXF_Real            max_row_sum;
   JXF_Real            trunc_factor;
   JXF_Real            S_commpkg_switch;

   JXF_Int      max_levels;
   JXF_Int      amg_logging;
   JXF_Int      amg_print_level;
   JXF_Int      debug_flag;
   JXF_Int      local_num_vars;
   JXF_Int      P_max_elmts;
   JXF_Int      R_max_size;

   JXF_Solver *smoother = NULL;
   JXF_Int        smooth_type = jxf_ParAMGDataSmoothType(amg_data);
   JXF_Int        smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data);
   char      *euclidfile;
   JXF_Int	      eu_level;
   JXF_Int	      eu_bj;
   JXF_Real     eu_sparse_A;

   JXF_Solver ysolver;
   JXF_Solver jac_solver;
   JXF_Int znum_sl = 0;
   JXF_Int iter_sl_used = 0;
   JXF_Int ziter_sl_used = 0;
   JXF_Real iter_rest;
   JXF_Real res_norm;
   JXF_Real res_norm_old;
   JXF_Real reduce_factor;
   JXF_Real reduce_factor_thrd = 1.0;

   /* Local variables */
   JXF_Real           *AI_measure;
   JXF_Real           *AI_measure2;
   JXF_Int              *CF_marker;
   JXF_Int              *CFN_marker;
   jxf_ParCSRMatrix  *par_S;
   jxf_ParCSRMatrix  *Sabs = NULL;
   jxf_ParCSRMatrix  *par_S2;
   jxf_ParCSRMatrix  *par_P;
   jxf_ParCSRMatrix  *R;
   jxf_hpCSRMatrix  *A_H;

   JXF_Int       wall_time_option;         /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real    wall_time_coarsen = 0.0;  /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real    wall_time_rap = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real    wall_time_interp = 0.0;   /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real    tmp_wall_time = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */

   JXF_Int       old_num_levels, num_levels;
   JXF_Int       level;
   JXF_Int       local_size, i;
   JXF_Int       first_local_row;
   JXF_Int       coarse_size_last = -1;
   JXF_Int       coarse_size;
   JXF_Int       coarsen_type;
   JXF_Int       interp_type;
   JXF_Int       restri_type;
   JXF_Int       measure_type;
   JXF_Int       setup_type;
   JXF_Int       fine_size;
   JXF_Int       rest, tms, indx;
   JXF_Int       not_finished_coarsening;
   JXF_Int       Setup_err_flag = 0;
   JXF_Int       coarse_threshold; /* we don't fix it to be 9.  peghoty 2010/04/14 */
   JXF_Int       coarse_threshold_2; /* by xu: 2013/11/25 */
   JXF_Int       j, k;
   JXF_Int       num_procs, my_id, num_threads;
   JXF_Int      *grid_relax_type = jxf_ParAMGDataGridRelaxType(amg_data);
   JXF_Int       num_functions = jxf_ParAMGDataNumFunctions(amg_data);
   JXF_Int       spmt_rap_type = jxf_ParAMGDataSpMtRapType(amg_data);
   JXF_Int       ai_measure_type = jxf_ParAMGDataAIMeasureType(amg_data);
   JXF_Int       num_paths = jxf_ParAMGDataNumPaths(amg_data);
   JXF_Int       agg_num_levels = jxf_ParAMGDataAggNumLevels(amg_data);
   JXF_Int       agg_interp_type = jxf_ParAMGDataAggInterpType(amg_data);
   JXF_Int       agg_P_max_elmts = jxf_ParAMGDataAggPMaxElmts(amg_data);
   //JXF_Int       agg_P12_max_elmts = jxf_ParAMGDataAggP12MaxElmts(amg_data);
   JXF_Real    agg_trunc_factor = jxf_ParAMGDataAggTruncFactor(amg_data);
   //JXF_Real    agg_P12_trunc_factor = jxf_ParAMGDataAggP12TruncFactor(amg_data);
   JXF_Int       rap2 = jxf_ParAMGDataRAP2(amg_data);
   JXF_Int       keepTranspose = jxf_ParAMGDataKeepTranspose(amg_data);
   JXF_Int       print_coarse_matrix = jxf_ParAMGDataPrintCoarseSystem(amg_data);
   JXF_Int      *opt_icor;
   JXF_Int      *coarse_dof_func;
   JXF_Int      *coarse_pnts_global;
   JXF_Int	    *coarse_pnts_global1;
   JXF_Real    size;
   JXF_Real    coarse_ratio;           /* newly added peghoty 2010/04/14 */
   JXF_Real    wall_time = 0.0;        /* for debugging instrumentation */
   JXF_Int       measure_type_rlx;       // newly added peghoty 2010/05/29
   JXF_Int       number_syn_rlx;         // newly added peghoty 2010/05/29
   JXF_Real    measure_threshold_rlx;  // newly added peghoty 2010/05/29

   /* coarse matrices */
   char FileNameCoaMat[256];

   /* ai statistic information*/
   JXF_Int       num_vars_local = 0, num_vars_global;
   JXF_Int       num_ai_local = 0, num_ai_global;
   JXF_Int       num_ai_local_valid = 0, num_ai_global_valid;
   JXF_Int       num_ai_c_local = 0, num_ai_c_global;
   JXF_Real    mai_local = 0.0, mai_global;
   JXF_Real    mai_local_valid = 0.0, mai_global_valid;
   JXF_Real    mai_c_local = 0.0, mai_c_global;

   JXF_Int       num_vars_local_0 = 0;
   JXF_Int       num_ai_local_0 = 0;
   JXF_Int       num_ai_local_valid_0 = 0;
   JXF_Int       num_ai_c_local_0 = 0;
   JXF_Real    mai_local_0 = 0.0;
   JXF_Real    mai_local_valid_0 = 0.0;
   JXF_Real    mai_c_local_0 = 0.0;
   JXF_Real    mai_threshold = 0.1;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   if ((num_procs > 1) && (spmt_rap_type != 1))  /* Yue Xiaoqiang 2012/10/12 */
   {
      spmt_rap_type = 1;
   }
   num_threads = jxf_NumThreads();               /* Yue Xiaoqiang 2012/10/12 */
   opt_icor = jxf_CTAlloc(JXF_Int, 5*num_threads+2); /* Yue Xiaoqiang 2012/10/12 */

   wall_time_option = jxf_ParAMGDataWallTimeOption(amg_data);
   old_num_levels   = jxf_ParAMGDataNumLevels(amg_data);
   max_levels       = jxf_ParAMGDataMaxLevels(amg_data);
   amg_logging      = jxf_ParAMGDataLogging(amg_data);
   amg_print_level  = jxf_ParAMGDataPrintLevel(amg_data);
   coarsen_type     = jxf_ParAMGDataCoarsenType(amg_data);
   measure_type     = jxf_ParAMGDataMeasureType(amg_data);
   setup_type       = jxf_ParAMGDataSetupType(amg_data);
   debug_flag       = jxf_ParAMGDataDebugFlag(amg_data);
   dof_func         = jxf_ParAMGDataDofFunc(amg_data);
   interp_type      = jxf_ParAMGDataInterpType(amg_data);
   restri_type      = jxf_ParAMGDataRestriction(amg_data);
   AIR_strong_th    = jxf_ParAMGDataAIRStrongTh(amg_data);
   euclidfile       = jxf_ParAMGDataEuclidFile(amg_data);
   eu_level         = jxf_ParAMGDataEuLevel(amg_data);
   eu_bj            = jxf_ParAMGDataEuBJ(amg_data);
   eu_sparse_A      = jxf_ParAMGDataEuSparseA(amg_data);
   coarse_threshold = jxf_ParAMGDataCoarseThreshold(amg_data);           /* newly added peghoty 2010/04/14 */
   coarse_ratio     = jxf_ParAMGDataCoarseRatio(amg_data);               /* newly added peghoty 2010/04/14 */
   measure_type_rlx = jxf_ParAMGDataMeasureTypeRlx(amg_data);            /* newly added peghoty 2010/05/29 */
   number_syn_rlx   = jxf_ParAMGDataNumberSynRlx(amg_data);              /* newly added peghoty 2010/05/29 */
   measure_threshold_rlx = jxf_ParAMGDataMeasureThresholdRlx(amg_data);  /* newly added peghoty 2010/05/29 */

   coarse_threshold_2 = -1;

   jxf_ParCSRMatrixSetNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix));
   jxf_ParCSRMatrixSetDNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix));
   jxf_ParAMGDataNumVariables(amg_data) = jxf_hpCSRMatrixNumRows(hp_matrix);

   if (setup_type == 0) 
   {
      return(Setup_err_flag);
   }

   par_S = NULL;
   A_H = NULL;

   A_array = jxf_hpAMGDataAArray(amg_data);
   P_array = jxf_ParAMGDataPArray(amg_data);
   R_array = jxf_ParAMGDataRArray(amg_data);
   AIR_maxsize_ls = jxf_ParAMGDataAIRMaxSizeLS(amg_data);
   AI_measure_array = jxf_ParAMGDataAIMeasureArray(amg_data);
   CF_marker_array = jxf_ParAMGDataCFMarkerArray(amg_data);
   dof_func_array  = jxf_ParAMGDataDofFuncArray(amg_data);
   local_size = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix));

   grid_relax_type[3] = jxf_ParAMGDataUserCoarseRelaxType(amg_data); 

   if (A_array || P_array || R_array || AIR_maxsize_ls || AI_measure_array || CF_marker_array || dof_func_array)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (A_array[j])
         {
            jxf_hpCSRMatrixDestroy(A_array[j]);
            A_array[j] = NULL;
         }
       
         if (dof_func_array[j])
         {
            jxf_TFree(dof_func_array[j]);
            dof_func_array[j] = NULL;
         }
      }

      for (j = 0; j < old_num_levels - 1; j ++)
      {
         if (P_array[j])
         {
            jxf_ParCSRMatrixDestroy(P_array[j]);
            P_array[j] = NULL;
         }
         if (R_array[j])
         {
            jxf_ParCSRMatrixDestroy(R_array[j]);
            R_array[j] = NULL;
         }
      }

      jxf_TFree(AIR_maxsize_ls);
      AIR_maxsize_ls = NULL;

     /*-------------------------------------------------------------------
      *  Special case use of CF_marker_array when old_num_levels = 1
      *  requires us to attempt this deallocation every time
      *------------------------------------------------------------------*/
      
      if (CF_marker_array[0])
      {
         jxf_TFree(CF_marker_array[0]);
         CF_marker_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (CF_marker_array[j])
         {
            jxf_TFree(CF_marker_array[j]);
            CF_marker_array[j] = NULL;
         }
      }
      
      /* for AI_measure: added by xwxu */
      if (AI_measure_array[0])
      {
         jxf_TFree(AI_measure_array[0]);
         AI_measure_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (AI_measure_array[j])
         {
            jxf_TFree(AI_measure_array[j]);
            AI_measure_array[j] = NULL;
         }
      }
   }

   if (A_array == NULL)
   {
      A_array = jxf_CTAlloc(jxf_hpCSRMatrix*, max_levels);
   }

   if (P_array == NULL && max_levels > 1)
   {
      P_array = jxf_CTAlloc(jxf_ParCSRMatrix*, max_levels - 1);
   }

   /* If retri_type != 0, R != P^T, allocate R matrices */
   if (restri_type)
   {
      if (R_array == NULL && max_levels > 1)
      {
         R_array = jxf_CTAlloc(jxf_ParCSRMatrix*, max_levels-1);
         AIR_maxsize_ls = jxf_CTAlloc(JXF_Int, max_levels-1);
      }
   }

   if (AI_measure_array == NULL)
   {
      AI_measure_array = jxf_CTAlloc(JXF_Real*, max_levels);
   }
   if (CF_marker_array == NULL)
   {
      CF_marker_array = jxf_CTAlloc(JXF_Int*, max_levels);
   }
   if (dof_func_array == NULL)
   {
      dof_func_array = jxf_CTAlloc(JXF_Int*, max_levels);
   }

   if (num_functions > 1 && dof_func == NULL)
   {
      first_local_row = jxf_hpCSRMatrixFirstRowIndex(hp_matrix);
      dof_func = jxf_CTAlloc(JXF_Int, local_size);
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
      jxf_ParAMGDataDofFunc(amg_data) = dof_func;
   }

   A_array[0] = hp_matrix;

   dof_func_array[0] = dof_func;
   jxf_ParAMGDataAIMeasureArray(amg_data) = AI_measure_array;
   jxf_ParAMGDataCFMarkerArray(amg_data) = CF_marker_array;
   jxf_ParAMGDataDofFuncArray(amg_data) = dof_func_array;
   jxf_hpAMGDataAArray(amg_data) = A_array;
   jxf_ParAMGDataPArray(amg_data) = P_array;
   /* If R != P^T */
   if (restri_type)
   {
      jxf_ParAMGDataRArray(amg_data) = R_array;
      jxf_ParAMGDataAIRMaxSizeLS(amg_data) = AIR_maxsize_ls;
   }
   else
   {
      jxf_ParAMGDataRArray(amg_data) = P_array;
   }

   Vtemp = jxf_ParAMGDataVtemp(amg_data);

   if (Vtemp != NULL)
   {
      jxf_ParVectorDestroy(Vtemp);
      Vtemp = NULL;
   }

   Vtemp = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[0]),
                              jxf_hpCSRMatrixGlobalNumRows(A_array[0]),
                              jxf_hpCSRMatrixRowStarts(A_array[0]));
   jxf_ParVectorInitialize(Vtemp);
   jxf_ParVectorSetPartitioningOwner(Vtemp,0);
   jxf_ParAMGDataVtemp(amg_data) = Vtemp;

   F_array = jxf_ParAMGDataFArray(amg_data);
   U_array = jxf_ParAMGDataUArray(amg_data);

   if (F_array != NULL || U_array != NULL)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (F_array[j] != NULL)
         {
            jxf_ParVectorDestroy(F_array[j]);
            F_array[j] = NULL;
         }
         if (U_array[j] != NULL)
         {
            jxf_ParVectorDestroy(U_array[j]);
            U_array[j] = NULL;
         }
      }
   }

   if (F_array == NULL)
   {
      F_array = jxf_CTAlloc(jxf_ParVector*, max_levels);
   }
   if (U_array == NULL)
   {
      U_array = jxf_CTAlloc(jxf_ParVector*, max_levels);
   }

   F_array[0] = par_rhs;
   U_array[0] = par_app;

   jxf_ParAMGDataFArray(amg_data) = F_array;
   jxf_ParAMGDataUArray(amg_data) = U_array;

  /*----------------------------------------------------------
   *   Initialize jxf_ParAMGData
   *---------------------------------------------------------*/
   not_finished_coarsening = 1;
   level = 0;
   strong_threshold = jxf_ParAMGDataStrongThreshold(amg_data);
   max_row_sum = jxf_ParAMGDataMaxRowSum(amg_data);
   trunc_factor = jxf_ParAMGDataTruncFactor(amg_data);
   P_max_elmts = jxf_ParAMGDataPMaxElmts(amg_data);
   S_commpkg_switch = jxf_ParAMGDataSCommPkgSwitch(amg_data);
   if (smooth_num_levels > level)
   {
      smoother = jxf_CTAlloc(JXF_Solver, smooth_num_levels);
      jxf_ParAMGDataSmoother(amg_data) = smoother;
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
      jxf_printf(" \n");
      jxf_printf("================================================================== \n");
      jxf_printf("+++++++++++++ multi-scale/AI infomation for levels +++++++++++++++ \n");
   }

   while (not_finished_coarsening)
   {
      fine_size = jxf_hpCSRMatrixGlobalNumRows(A_array[level]);

      if (level > 0)
      {
         F_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                             jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jxf_hpCSRMatrixRowStarts(A_array[level]));
         jxf_ParVectorInitialize(F_array[level]);
         jxf_ParVectorSetPartitioningOwner(F_array[level], 0);

         U_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                             jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jxf_hpCSRMatrixRowStarts(A_array[level]));
         jxf_ParVectorInitialize(U_array[level]);
         jxf_ParVectorSetPartitioningOwner(U_array[level], 0);
      }

      if (level > 0)
      {
   iter_sl_used = 0;
   reduce_factor_thrd = 1.0;
         jxf_ParVectorCopy(F_array[level-1], Vtemp);
         jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar(A_array[level-1]), U_array[level-1], 1.0, Vtemp);
         jxf_ParCSRMatrixMatvecT(1.0, P_array[level-1], Vtemp, 0.0, F_array[level]);
   JXF_JacobiCreate(&jac_solver);
   JXF_JacobiSetMaxIter(jac_solver, 1);
   JXF_JacobiSetup(jac_solver, (JXF_hpCSRMatrix) (A_array[level]));
         if (solver_id == 1)
         {
      JXF_PCGCreate(comm, &ysolver);
      JXF_PCGSetMaxIter(ysolver, 1);
      JXF_PCGSetTol(ysolver, ytol);
      JXF_PCGSetTwoNorm(ysolver, ytwo_norm);
      JXF_PCGSetPrintLevel(ysolver, yprint_level);
      JXF_PCGSetLogging(ysolver, 1);
      JXF_PCGSetPrecond(ysolver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_PCGSetup(ysolver, (JXF_Matrix) A_array[level], (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_PCGSolve(ysolver, (JXF_Matrix) A_array[level], (JXF_Matrix) A_array[level],
                              (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
         JXF_PCGGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JXF_PCGDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JXF_PCGDestroy(ysolver);
         break;
      }
         }
         else if (solver_id == 2)
         {
      JXF_GMRESCreate(comm, &ysolver);
      JXF_GMRESSetMaxIter(ysolver, 1);
      JXF_GMRESSetPrintLevel(ysolver, yprint_level);
      JXF_GMRESSetKDim(ysolver, yk_dim);
      JXF_GMRESSetIsCheckRestarted(ysolver, yis_check_restarted);
      JXF_GMRESSetTol(ysolver, ytol);
      JXF_GMRESSetLogging(ysolver, 1);
      JXF_GMRESSetPrecond(ysolver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_GMRESSetup(ysolver, (JXF_Matrix) A_array[level], (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_GMRESSolve(ysolver, (JXF_Matrix) A_array[level], (JXF_Matrix) A_array[level],
                                (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
         JXF_GMRESGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JXF_GMRESDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JXF_GMRESDestroy(ysolver);
         break;
      }
         }
         else if (solver_id == 3)
         {
      JXF_BiCGSTABCreate(comm, &ysolver);
      JXF_BiCGSTABSetMaxIter(ysolver, 1);
      JXF_BiCGSTABSetPrintLevel(ysolver, yprint_level);
      JXF_BiCGSTABSetTol(ysolver, ytol);
      JXF_BiCGSTABSetAbsoluteTol(ysolver, 0.0);
      JXF_BiCGSTABSetConvCriteria(ysolver, 0);
      JXF_BiCGSTABSetLogging(ysolver, 1);
      JXF_BiCGSTABSetPrecond(ysolver, (JXF_PtrToSolverFcn) JXF_JacobiPrecond, (JXF_PtrToSolverFcn) JXF_JacobiSetup, jac_solver);
      JXF_BiCGSTABSetup(ysolver, (JXF_Matrix) A_array[level], (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
      res_norm = 1.0;
      do {
         if (yis_adaptive)
         {
            iter_rest = ymax_iter_sl - iter_sl_used;
            iter_rest = 1.0 / iter_rest;
            reduce_factor_thrd = pow(ytol/res_norm, iter_rest);
         }
         res_norm_old = res_norm;
         JXF_BiCGSTABSolve(ysolver, (JXF_Matrix) A_array[level], (JXF_Matrix) A_array[level],
                                   (JXF_Vector) F_array[level], (JXF_Vector) U_array[level]);
         JXF_BiCGSTABGetFinalRelativeResidualNorm(ysolver, &res_norm);
         iter_sl_used ++;
         reduce_factor = res_norm / res_norm_old;
      } while ((res_norm > ytol) && (reduce_factor < reduce_factor_thrd));

      iter_sl_lvls[level-1] = iter_sl_used;
      if (res_norm > ytol)
      {
         ziter_sl_used += iter_sl_used;
         JXF_BiCGSTABDestroy(ysolver);
      }
      else
      {
         znum_sl ++;
         JXF_BiCGSTABDestroy(ysolver);
         break;
      }
         }
         JXF_JacobiDestroy(jac_solver);
      }

     /*--------------------------------------------------------------
      *  Select coarse-grid points on 'level' : returns CF_marker
      *  for the level.  Returns strength matrix, par_S 
      *  Returns AI_Measure. 
      *--------------------------------------------------------------*/
     
      if (debug_flag == 1) wall_time = jxf_time_getWallclockSeconds();
      if (debug_flag == 3)
      {
         jxf_printf("\n ===== Proc = %d     Level = %d  =====\n",my_id, level);
         fflush(NULL);
      }

      if (wall_time_option == 1) tmp_wall_time = jxf_time_getWallclockSeconds();

      if (max_levels > 1)
      {
         local_num_vars = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level]));

        /*--------------------------------------------------------------------------
         *  Get the Strength Matrix 
         *-------------------------------------------------------------------------*/        
         jxf_PAMGCreateS(jxf_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }
         /* for AIR, need absolute value SOC */
         if (restri_type)
         {
            jxf_PAMGCreateSabs(jxf_hpCSRMatrixPar(A_array[level]),AIR_strong_th,1.0,num_functions,dof_func_array[level],&Sabs);
            col_offd_Sabs_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), Sabs, &col_offd_Sabs_to_A);
            }
         }

         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jxf_PAMGMeasureAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

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

         //jxf_printf("=========== coarsen_type = %d \n", coarsen_type); 

         if (coarsen_type == 6)  
         {
            jxf_PAMGCoarsenFalgout(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 96) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenFalgoutAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 8) 
         {
            jxf_PAMGCoarsenPMIS(par_S, jxf_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 10) 
         {
            jxf_PAMGCoarsenHMIS(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }      
         else if ((coarsen_type == 910) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenHMISAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, debug_flag, &CF_marker);
         }      
         else if (coarsen_type == 1 || coarsen_type == 2 ||
                  coarsen_type == 3 || coarsen_type == 4 ||
                  coarsen_type == 5 || coarsen_type == 11  ) 
         {
            jxf_PAMGCoarsenRuge(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, coarsen_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 991) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 11, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 992) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 1, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 993) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 3, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 0)    
         {
            jxf_PAMGCoarsen(par_S,jxf_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 990) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 90)
         {
            jxf_PAMGCoarsenRCLJP(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type_rlx, number_syn_rlx,
                                measure_threshold_rlx, 0, debug_flag, &CF_marker);
         }         
         else if (coarsen_type == 91)
         {
            jxf_PAMGCoarsenRRS0(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, measure_type_rlx,
                               number_syn_rlx, measure_threshold_rlx, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 98) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenPMISAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || coarsen_type == 938 || 
                  coarsen_type == 968) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenXML(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, coarsen_type, debug_flag, &CF_marker);
         }      

         if (level < agg_num_levels)
         {
            jxf_PAMGCoarseParms(comm, local_num_vars, 1,
                               dof_func_array[level], CF_marker,
                               &coarse_dof_func, &coarse_pnts_global1);
            jxf_PAMGCreate2ndS(par_S, CF_marker, num_paths, coarse_pnts_global1, &par_S2);

           /*--------------------------------------------------------------------------
            *  Get the AI Measure
            *-------------------------------------------------------------------------*/
            if (ai_measure_type == 1)
            {
               jxf_PAMGMeasureAI(par_S2, par_S2, debug_flag, &AI_measure2);
            }

            if (coarsen_type == 10)
            {
               jxf_PAMGCoarsenHMIS(par_S2, par_S2, measure_type+3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 910) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenHMISAI(par_S2, par_S2, AI_measure2, 0, measure_type+3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 8)
            {
               jxf_PAMGCoarsenPMIS(par_S2, par_S2, 3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 98) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenPMISAI(par_S2, par_S2, AI_measure2, 3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 6)
            {
               jxf_PAMGCoarsenFalgout(par_S2, par_S2, measure_type, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 96) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenFalgoutAI(par_S2, par_S2, AI_measure2, measure_type, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 0)
            {
               jxf_PAMGCoarsen(par_S2, par_S2, 0, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 990) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenAI(par_S2, par_S2, AI_measure2, 0, debug_flag, &CFN_marker);
            }
            jxf_ParCSRMatrixDestroy(par_S2);
            if (ai_measure_type == 1)
            {
               jxf_TFree(AI_measure2);
            }
         }

         /* Here for changes of min_coarse_size */

         if (level < agg_num_levels)
         {
            if (agg_interp_type == 4)
            {
               jxf_PAMGCorrectCFMarker(CF_marker, local_num_vars, CFN_marker);
               jxf_TFree(coarse_pnts_global1);
               jxf_TFree(CFN_marker);
            }
         }

        /*--------------------------------------------
         *  store the CF array
         *-------------------------------------------*/ 
         CF_marker_array[level] = CF_marker;

         if (debug_flag == 1)
         {
            wall_time = jxf_time_getWallclockSeconds() - wall_time;
            jxf_printf("Proc = %d    Level = %d    Coarsen Time = %f\n",my_id,level, wall_time); 
            fflush(NULL);
         }

         if (wall_time_option == 1)
         {
            wall_time_coarsen += (jxf_time_getWallclockSeconds() - tmp_wall_time);
         }

        /*-------------------------------------------------------------------------
         *  Get the coarse parameters
         *-------------------------------------------------------------------------*/ 
         jxf_PAMGCoarseParms(comm, local_num_vars, num_functions, 
                            dof_func_array[level], CF_marker,
                            &coarse_dof_func, &coarse_pnts_global);

         dof_func_array[level+1] = NULL;
         if (num_functions > 1) dof_func_array[level+1] = coarse_dof_func;

#ifdef JXF_NO_GLOBAL_PARTITION
         if (my_id == (num_procs -1)) 
         {
            coarse_size = coarse_pnts_global[1];
         }
         jxf_MPI_Bcast(&coarse_size, 1, JXF_MPI_INT, num_procs-1, comm);
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
         jxf_PAMGCreateS(jxf_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }

         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jxf_PAMGMeasureAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         par_S = NULL;
         coarse_pnts_global = NULL;
         CF_marker = jxf_CTAlloc(JXF_Int, local_size );
         for (i = 0; i < local_size ; i ++) 
         {
            CF_marker[i] = 1;
         }
         CF_marker_array[level] = CF_marker;
         coarse_size = fine_size;

         local_num_vars = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level]));
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
                      //if (level == 0 ) jxf_printf("i_c = %d, ai_measure= %f \n", i, AI_measure[i]);
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
         JXF_Int  *num_grid_sweeps   = jxf_ParAMGDataNumGridSweeps(amg_data);
         JXF_Int **grid_relax_points = jxf_ParAMGDataGridRelaxPoints(amg_data);
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
            jxf_ParCSRMatrixDestroy(par_S);
         }
         jxf_TFree(coarse_pnts_global);
         if (level > 0)
         {
           /*-------------------------------------------------------------
            * Note special case treatment of CF_marker is necessary
            * to do CF relaxation correctly when num_levels = 1
            *------------------------------------------------------------*/ 
            jxf_TFree(CF_marker_array[level]);
            jxf_ParVectorDestroy(F_array[level]);
            jxf_ParVectorDestroy(U_array[level]);
         }
         break; 
      }


      /* Build restriction */
      if (restri_type)
      {
         /* !!! Ensure that CF_marker contains -1 or 1 !!! */
         for (i = 0; i < jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level])); i ++)
         {
            CF_marker[i] = CF_marker[i] > 0 ? 1 : -1;
         }
         if (restri_type == 1) /* distance-1 AIR */
         {
            jxf_PAMGBuildRestrAIR(jxf_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                 Sabs, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], 
                                 debug_flag, trunc_factor, P_max_elmts, 
                                 col_offd_Sabs_to_A, &R, &R_max_size );
         }
         else /* distance-2 AIR */
         {
            jxf_PAMGBuildRestrDist2AIR(jxf_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                      Sabs, coarse_pnts_global, num_functions, 
                                      dof_func_array[level], 
                                      debug_flag, trunc_factor, P_max_elmts, 
                                      col_offd_Sabs_to_A, &R, &R_max_size );
         }
         if (Sabs)
         {
            jxf_ParCSRMatrixDestroy(Sabs);
            Sabs = NULL;
         }
         jxf_TFree(col_offd_Sabs_to_A);
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
       
      if (debug_flag == 1) wall_time = jxf_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jxf_time_getWallclockSeconds();

      if (level < agg_num_levels)
      {
         if (agg_interp_type == 4)
         {
            jxf_PAMGBuildMultipass(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S, coarse_pnts_global,
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
               jxf_PAMGBuildInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, trunc_factor,
                                P_max_elmts, col_offd_S_to_A, &par_P ); 
            }
            else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
            {
               jxf_PAMGBuildInterp1( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                 coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag, trunc_factor,
                                 P_max_elmts, col_offd_S_to_A, &par_P, opt_icor );
            }
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 3)  // direct interpolation (with separation of weights)
         {
            jxf_PAMGBuildDirInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, 
                                trunc_factor, P_max_elmts, 
                                col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         } 
         else if (interp_type == 4)  // multipass interpolation
         {
            jxf_PAMGBuildMultipass( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 5) // multipass interpolation (with separation of weights)
         {
            jxf_PAMGBuildMultipass( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 6)  // extended classical modified interpolation
         {
             jxf_PAMGBuildExtPIInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                   par_S, coarse_pnts_global, num_functions, 
                                   dof_func_array[level], debug_flag,
                                   trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
             jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 7) // extended (if no common C neighbor) classical modified interpolation
         {
            jxf_PAMGBuildExtPICCInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                     par_S, coarse_pnts_global, num_functions, 
                                     dof_func_array[level], debug_flag,
                                     trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 8) // standard interpolation
         {
            jxf_PAMGBuildStdInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 9) // standard interpolation (with separation of weights)
         {
            jxf_PAMGBuildStdInterp( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 100) /* 1pt interpolation */
         {
            jxf_PAMGBuildInterpOnePnt( jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                    par_S, coarse_pnts_global, num_functions,
                                    dof_func_array[level], debug_flag, col_offd_S_to_A, &par_P );
            jxf_TFree(col_offd_S_to_A);
         }
      }
      
      if (debug_flag == 1)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         jxf_printf("Proc = %d    Level = %d    Build Interp Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_interp += (jxf_time_getWallclockSeconds() - tmp_wall_time);
      }

      P_array[level] = par_P; 
      if (restri_type)
      {
         R_array[level] = R;
         AIR_maxsize_ls[level] = R_max_size;
      }

      if (par_S) 
      {
         jxf_ParCSRMatrixDestroy(par_S);
      }
      par_S = NULL;


     /*--------------------------------------------------------------------------------
      *   Build coarse-grid operator, A_array[level+1] by R*A*P
      *-------------------------------------------------------------------------------*/
       
      if (debug_flag == 1) wall_time = jxf_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jxf_time_getWallclockSeconds();

      if (restri_type)
      {
         /* Use two matrix products to generate A_H */
         jxf_ParCSRMatrix *Q = jxf_ParMatmul(jxf_hpCSRMatrixPar(A_array[level]), P_array[level]);
         A_H = jxf_hpInithpCSRMatrix();
         jxf_hpCSRMatrixPar(A_H) = jxf_ParMatmul(R_array[level], Q);
         jxf_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jxf_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jxf_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         jxf_ParCSRMatrixOwnsRowStarts(R_array[level]) = 0;
         if (num_procs > 1) jxf_MatvecCommPkgCreate(jxf_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jxf_ParCSRMatrixDestroy(Q);
      }
      else if (rap2)
      {
         /* Use two matrix products to generate A_H */
         jxf_ParCSRMatrix *Q = jxf_ParMatmul(jxf_hpCSRMatrixPar(A_array[level]), P_array[level]);
         A_H = jxf_hpInithpCSRMatrix();
         jxf_hpCSRMatrixPar(A_H) = jxf_ParTMatmul(P_array[level], Q);
         jxf_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jxf_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jxf_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         if (num_procs > 1) jxf_MatvecCommPkgCreate(jxf_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jxf_ParCSRMatrixDestroy(Q);
      }
      else if (spmt_rap_type == 1)
      {
         A_H = jxf_hpInithpCSRMatrix();
         jxf_PAMGBuildCoarseOperatorKT( P_array[level], jxf_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &jxf_hpCSRMatrixPar(A_H) );
      }
      else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
      {
         A_H = jxf_hpInithpCSRMatrix();
         jxf_PAMGBuildCoarseOperatorOMP( P_array[level], jxf_hpCSRMatrixPar(A_array[level]), P_array[level], &jxf_hpCSRMatrixPar(A_H), opt_icor );     
      }
      
      if (debug_flag == 1)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         jxf_printf("Proc = %d    Level = %d    Build Coarse Operator Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_rap += (jxf_time_getWallclockSeconds() - tmp_wall_time);
      }

      ++ level;

      jxf_ParCSRMatrixSetNumNonzeros(jxf_hpCSRMatrixPar(A_H));
      jxf_ParCSRMatrixSetDNumNonzeros(jxf_hpCSRMatrixPar(A_H));
      A_array[level] = A_H;

      if (coarse_size <= coarse_threshold_2) {

         coarse_size_last = coarse_size;

      }

      /* print coarser operator. */
	  if (print_coarse_matrix == 1) {

         jxf_sprintf(FileNameCoaMat, "cmat_%d", level);
         jxf_hpCSRMatrixPrint(A_array[level], FileNameCoaMat);

      }
      
     /*-------------------------------------------------------------------------------
      *   Switch to CLJP when coarsening slows
      *                                            peghoty  2009/07/09
      *------------------------------------------------------------------------------*/ 
      
      size = ((JXF_Real) fine_size )*coarse_ratio;   /* peghoty 2010/04/14 */
      if (coarsen_type > 0 && coarse_size >= (JXF_Int) size)
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
            jxf_PAMGCreateS(A_array[level],strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
            col_offd_S_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jxf_PAMGCreateSCommPkg(A_array[level], par_S, &col_offd_S_to_A);
            }

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jxf_PAMGMeasureAI(par_S, A_array[level], debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         /*-------------------------------------------------------------------------
          *  compute ai-information for all-pnts and C-pnts!
          *-------------------------------------------------------------------------*/ 
         local_num_vars = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A_array[level]));

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

      /*jxf_printf("level = %d, not_finished_coarsening = %d \n", level, not_finished_coarsening);
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

      jxf_MPI_Allreduce(&num_ai_local_0,&num_ai_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_ai_local_valid_0,&num_ai_global_valid,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_vars_local_0,&num_vars_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_ai_c_local_0,&num_ai_c_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_local_0,&mai_global,1,JXF_MPI_REAL,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_local_valid_0,&mai_global_valid,1,JXF_MPI_REAL,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_c_local_0,&mai_c_global,1,JXF_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JXF_Real num_ai = num_ai_global;
         JXF_Real num_ai_valid = num_ai_global_valid;
         JXF_Real num_var = num_vars_global;
         JXF_Real num_ai_c = num_ai_c_global;
         jxf_printf(" \n");
         jxf_printf("Level 0: \n");
         jxf_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jxf_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jxf_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai) jxf_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid) jxf_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jxf_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jxf_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jxf_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
      }

      jxf_MPI_Allreduce(&num_ai_local,&num_ai_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_ai_local_valid,&num_ai_global_valid,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_vars_local,&num_vars_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&num_ai_c_local,&num_ai_c_global,1,JXF_MPI_INT,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_local,&mai_global,1,JXF_MPI_REAL,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_local_valid,&mai_global_valid,1,JXF_MPI_REAL,MPI_SUM,comm);
      jxf_MPI_Allreduce(&mai_c_local,&mai_c_global,1,JXF_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JXF_Real num_ai = num_ai_global;
         JXF_Real num_ai_valid = num_ai_global_valid;
         JXF_Real num_var = num_vars_global;
         JXF_Real num_ai_c = num_ai_c_global;
         jxf_printf(" \n");
         jxf_printf("All levels (except coarsest level): \n");
         jxf_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jxf_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jxf_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai>0) jxf_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid>0) jxf_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jxf_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jxf_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jxf_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
         jxf_printf("================================================================== \n");
      }
   }

   if ((not_finished_coarsening == 0) && (level > 0))
   {
      F_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jxf_hpCSRMatrixRowStarts(A_array[level]));
      jxf_ParVectorInitialize(F_array[level]);
      jxf_ParVectorSetPartitioningOwner(F_array[level], 0);

      U_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jxf_hpCSRMatrixRowStarts(A_array[level]));
      jxf_ParVectorInitialize(U_array[level]);
      jxf_ParVectorSetPartitioningOwner(U_array[level], 0);
   }

   if (not_finished_coarsening == 0)
   {
      jxf_ParVectorCopy(F_array[level-1], Vtemp);
      jxf_hpCSRMatrixMatvec(-1.0, A_array[level-1], U_array[level-1], 1.0, Vtemp);
      jxf_ParCSRMatrixMatvecT(1.0, P_array[level-1], Vtemp, 0.0, F_array[level]);
      jxf_PAMGRelax9(jxf_hpCSRMatrixPar(A_array[level]), F_array[level], NULL, 0, 1.0, 0.0, U_array[level], NULL);
   }

   for (i = level; i > 0; i --)
   {
      jxf_ParCSRMatrixMatvec(1.0, P_array[i-1], U_array[i], 1.0, U_array[i-1]);
      jxf_PAMGRelax0(jxf_hpCSRMatrixPar(A_array[i-1]), F_array[i-1], NULL, 0, 1.0, 0.0, U_array[i-1], Vtemp);
   }

  /*-----------------------------------------------------------------------
   * enter all the stuff created, A[level], P[level], CF_marker[level],
   * for levels 1 through coarsest, into amg_data data structure
   *-----------------------------------------------------------------------*/

   num_levels = level + 1;
   jxf_ParAMGDataNumLevels(amg_data) = num_levels;
   if (jxf_ParAMGDataSmoothNumLevels(amg_data) > level)
   {
      jxf_ParAMGDataSmoothNumLevels(amg_data) = level;
   }
   smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data);
   
   for (j = 0; j < smooth_num_levels; j++) // Euclid smoothers
   {
      if (smooth_type == 9 || smooth_type == 19)
      {
         JXF_EuclidCreate(comm, &smoother[j]);
         if (euclidfile)
         {
            JXF_EuclidSetParamsFromFile(smoother[j], euclidfile);
         }
         JXF_EuclidSetLevel(smoother[j], eu_level);
         if (eu_bj)
         {
            JXF_EuclidSetBJ(smoother[j], eu_bj);
         }
         if (eu_sparse_A)
         {
            JXF_EuclidSetSparseA(smoother[j], eu_sparse_A);
         }
         JXF_EuclidSetup(smoother[j], (JXF_hpCSRMatrix)A_array[j]);
      }
   }

   if (amg_logging > 1) 
   {
      Residual_array = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[0]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[0]),
                                          jxf_hpCSRMatrixRowStarts(A_array[0]));
      jxf_ParVectorInitialize(Residual_array);
      jxf_ParVectorSetPartitioningOwner(Residual_array,0);
      jxf_ParAMGDataResidual(amg_data) = Residual_array;
   }
   else
   {
      jxf_ParAMGDataResidual(amg_data) = NULL;
   }

   jxf_TFree(opt_icor); /* Yue Xiaoqiang 2012/10/12 */
   if (Jxf_Pmarkers_global_rap)
   {
      jxf_TFree(Jxf_Pmarkers_global_rap); /* Yue Xiaoqiang 2012/10/17 */
   }
   if (Jxf_Pmarkers_global_interp)
   {
      jxf_TFree(Jxf_Pmarkers_global_interp); /* Yue Xiaoqiang 2012/10/17 */
   }

  /*--------------------------------------------------------------
   *    Print some stuff
   *-------------------------------------------------------------*/
   
   //if (amg_print_level == 1 || amg_print_level == 3)
   if (amg_print_level > 1)
   {
      /* Write the SETUP parameters */
      jxf_PAMGSetupStatus(amg_data, jxf_hpCSRMatrixPar(hp_matrix));
   }

   if (wall_time_option == 1)
   {
      jxf_printf("\n\nProc = %d, Coarsen Time = %f\n", my_id, wall_time_coarsen);
      jxf_printf("Proc = %d, Build Coarse Operator Time = %f\n", my_id, wall_time_rap);
      jxf_printf("Proc = %d, Build Interp Time = %f\n\n", my_id, wall_time_interp);
   }

   if (my_id == 0 && amg_print_level > 1)
   {
      /* Write the SOLVE parameters */
      jxf_PAMGSolveStatus(amg_data); /* Yue Xiaoqiang 2014/04/12 */
   }

   return(Setup_err_flag);
}

JXF_Int
JXF_ParAdpSetupAMGGetNumSL( JXF_Solver solver, JXF_Int *num_sl )
{
   return( jxf_ParAdpSetupAMGGetNumSL( (void *) solver, num_sl ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetNumSL( void *data, JXF_Int *num_sl )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *num_sl = jxf_ParAdpSetupAMGDataNumSL(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetNumReUse( JXF_Solver solver, JXF_Int *num_reuse )
{
   return( jxf_ParAdpSetupAMGGetNumReUse( (void *) solver, num_reuse ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetNumReUse( void *data, JXF_Int *num_reuse )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *num_reuse = jxf_ParAdpSetupAMGDataNumReUse(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetCanReUse( JXF_Solver solver, JXF_Int *can_reuse )
{
   return( jxf_ParAdpSetupAMGGetCanReUse( (void *) solver, can_reuse ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetCanReUse( void *data, JXF_Int *can_reuse )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *can_reuse = jxf_ParAdpSetupAMGDataCanReUse(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetNumSLAMG( JXF_Solver solver, JXF_Int *num_sl_amg )
{
   return( jxf_ParAdpSetupAMGGetNumSLAMG( (void *) solver, num_sl_amg ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetNumSLAMG( void *data, JXF_Int *num_sl_amg )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *num_sl_amg = jxf_ParAdpSetupAMGDataNumSLAMG(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetCurrentSL( JXF_Solver solver, JXF_Int *iter_sl, JXF_Int *converge, JXF_Int *num_lvls, JXF_Int **iter_sl_lvls, JXF_Int *update )
{
   return( jxf_ParAdpSetupAMGGetCurrentSL( (void *) solver, iter_sl, converge, num_lvls, iter_sl_lvls, update ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetCurrentSL( void *data, JXF_Int *iter_sl, JXF_Int *converge, JXF_Int *num_lvls, JXF_Int **iter_sl_lvls, JXF_Int *update )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *iter_sl = jxf_ParAdpSetupAMGDataCurrentIterSL(amg_data);
  *converge = jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data);
  *num_lvls = jxf_ParAdpSetupAMGDataCurrentNumLvls(amg_data);
  *iter_sl_lvls = jxf_ParAdpSetupAMGDataIterSLLvls(amg_data);
  *update = jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetTotalNumLevels( JXF_Solver solver, JXF_Int *total_num_levels )
{
   return( jxf_ParAdpSetupAMGGetTotalNumLevels( (void *) solver, total_num_levels ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetTotalNumLevels( void *data, JXF_Int *total_num_levels )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *total_num_levels = jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetTotalRedunIterSL( JXF_Solver solver, JXF_Int *total_redun_iter_sl )
{
   return( jxf_ParAdpSetupAMGGetTotalRedunIterSL( (void *) solver, total_redun_iter_sl ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetTotalRedunIterSL( void *data, JXF_Int *total_redun_iter_sl )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *total_redun_iter_sl = jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetTotalConvergeIter( JXF_Solver solver, JXF_Int *total_converge_iter )
{
   return( jxf_ParAdpSetupAMGGetTotalConvergeIter( (void *) solver, total_converge_iter ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetTotalConvergeIter( void *data, JXF_Int *total_converge_iter )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *total_converge_iter = jxf_ParAdpSetupAMGDataTotalConvergeIter(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetTotalRedunIterAMG( JXF_Solver solver, JXF_Int *total_redun_iter_amg )
{
   return( jxf_ParAdpSetupAMGGetTotalRedunIterAMG( (void *) solver, total_redun_iter_amg ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetTotalRedunIterAMG( void *data, JXF_Int *total_redun_iter_amg )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *total_redun_iter_amg = jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetSetupTime( JXF_Solver solver, JXF_Real *setup_time )
{
   return( jxf_ParAdpSetupAMGGetSetupTime( (void *) solver, setup_time ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetSetupTime( void *data, JXF_Real *setup_time )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *setup_time = jxf_ParAdpSetupAMGDataSetupTime(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetSolveTime( JXF_Solver solver, JXF_Real *solve_time )
{
   return( jxf_ParAdpSetupAMGGetSolveTime( (void *) solver, solve_time ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetSolveTime( void *data, JXF_Real *solve_time )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! ParAdpSetupAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *solve_time = jxf_ParAdpSetupAMGDataSetupTime(amg_data);
   return jxf_error_flag;
}

JXF_Int
JXF_ParAdpSetupAMGGetStatus( JXF_Solver solver )
{
   return( jxf_ParAdpSetupAMGGetStatus( (void *) solver ) );
}

JXF_Int
jxf_ParAdpSetupAMGGetStatus( void *data )
{
   jxf_ParAdpSetupAMGData *amg_data = data;
   MPI_Comm comm = jxf_ParAdpSetupAMGDataComm(amg_data);
   JXF_Bool is_adaptive = jxf_ParAdpSetupAMGDataIsAdaptive(amg_data);
   JXF_Bool is_adap_setup = jxf_ParAdpSetupAMGDataIsAdapSetup(amg_data);
   JXF_Real tol = jxf_ParAdpSetupAMGDataTol(amg_data);
   JXF_Real max_row_sum = jxf_ParAdpSetupAMGDataMaxRowSum(amg_data);
   JXF_Real strong_threshold = jxf_ParAdpSetupAMGDataStrongThreshold(amg_data);
   JXF_Int max_levels = jxf_ParAdpSetupAMGDataMaxLevels(amg_data);
   JXF_Int cycle_type = jxf_ParAdpSetupAMGDataCycleType(amg_data);
   JXF_Int relax_type = jxf_ParAdpSetupAMGDataRelaxType(amg_data);
   JXF_Int relax_order = jxf_ParAdpSetupAMGDataRelaxOrder(amg_data);
   JXF_Int ns_down = jxf_ParAdpSetupAMGDataNsDown(amg_data);
   JXF_Int ns_up = jxf_ParAdpSetupAMGDataNsUp(amg_data);
   JXF_Int ns_coarse = jxf_ParAdpSetupAMGDataNsCoarse(amg_data);
   JXF_Int coarsen_type = jxf_ParAdpSetupAMGDataCoarsenType(amg_data);
   JXF_Int interp_type = jxf_ParAdpSetupAMGDataInterpType(amg_data);
   JXF_Int P_max_elmts = jxf_ParAdpSetupAMGDataPMaxElmts(amg_data);
   JXF_Int agg_num_levels = jxf_ParAdpSetupAMGDataAggNumLevels(amg_data);
   JXF_Int coarse_threshold = jxf_ParAdpSetupAMGDataCoarseThreshold(amg_data);
   JXF_Int amg_print_level = jxf_ParAdpSetupAMGDataAMGPrintLevel(amg_data);
   JXF_Int k_dim = jxf_ParAdpSetupAMGDataKDim(amg_data);
   JXF_Int TTest = jxf_ParAdpSetupAMGDataTTest(amg_data);
   JXF_Int max_iter = jxf_ParAdpSetupAMGDataMaxIter(amg_data);
   JXF_Int two_norm = jxf_ParAdpSetupAMGDataTwoNorm(amg_data);
   JXF_Int solver_id = jxf_ParAdpSetupAMGDataSolverID(amg_data);
   JXF_Int print_level = jxf_ParAdpSetupAMGDataPrintLevel(amg_data);
   JXF_Int max_iter_sl = jxf_ParAdpSetupAMGDataMaxIterSL(amg_data);
   JXF_Int max_iter_reuse = jxf_ParAdpSetupAMGDataMaxIterReUse(amg_data);
   JXF_Int is_check_restarted = jxf_ParAdpSetupAMGDataIsCheckRestarted(amg_data);
   JXF_Int my_id;

   jxf_MPI_Comm_rank(comm, &my_id);
   if (my_id == 0)
   {
      if (is_adaptive) jxf_printf(" is_adaptive = TRUE\n");
      else jxf_printf(" is_adaptive = FALSE\n");
      if (is_adap_setup) jxf_printf(" is_adap_setup = TRUE\n");
      else jxf_printf(" is_adap_setup = FALSE\n");
      jxf_printf(" tol = %lf\n", tol);
      jxf_printf(" max_row_sum = %lf\n", max_row_sum);
      jxf_printf(" strong_threshold = %lf\n", strong_threshold);
      jxf_printf(" max_levels = %d\n", max_levels);
      jxf_printf(" cycle_type = %d\n", cycle_type);
      jxf_printf(" relax_type = %d\n", relax_type);
      jxf_printf(" relax_order = %d\n", relax_order);
      jxf_printf(" ns_down, up, coarse = %d, %d, %d\n", ns_down, ns_up, ns_coarse);
      jxf_printf(" coarsen_type = %d\n", coarsen_type);
      jxf_printf(" interp_type = %d\n", interp_type);
      jxf_printf(" P_max_elmts = %d\n", P_max_elmts);
      jxf_printf(" agg_num_levels = %d\n", agg_num_levels);
      jxf_printf(" coarse_threshold = %d\n", coarse_threshold);
      jxf_printf(" amg_print_level = %d\n", amg_print_level);
      jxf_printf(" k_dim = %d\n", k_dim);
      jxf_printf(" TTest = %d\n", TTest);
      jxf_printf(" max_iter = %d\n", max_iter);
      jxf_printf(" two_norm = %d\n", two_norm);
      jxf_printf(" solver_id = %d\n", solver_id);
      jxf_printf(" print_level = %d\n", print_level);
      jxf_printf(" max_iter_sl = %d\n", max_iter_sl);
      jxf_printf(" max_iter_reuse = %d\n", max_iter_reuse);
      jxf_printf(" is_check_restarted = %d\n", is_check_restarted);
   }

   return 0;
}
