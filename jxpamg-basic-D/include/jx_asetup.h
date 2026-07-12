//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jx_asetup.h -- head files for parallel asetup-AMG solver
 *  Date: 2020/03/10
 */ 


#ifndef JX_ASETUP_HEADER
#define JX_ASETUP_HEADER

#ifndef JX_PAMG_HEADER
#include "jx_pamg.h"
#endif
 
#ifndef JX_KRYLOV_HEADER
#include "jx_krylov.h"
#endif

/*----------------------------------------------------------------*
*                   Struct Declaration                           *
*----------------------------------------------------------------*/

/*!
 * \struct jx_JacobiData
 */
typedef struct
{
   JX_Int max_iter;

   jx_ParVector *tmp_vec;
    
} jx_JacobiData;

#define jx_JacobiDataMaxIter(amg_data)  ((amg_data)->max_iter)
#define jx_JacobiDataTmpVec(amg_data)   ((amg_data)->tmp_vec)

/*!
 * \struct jx_ParAdpSetupAMGData
 */
typedef struct
{
   MPI_Comm comm;

   JX_Bool is_adaptive;
   JX_Bool is_adap_setup;

   JX_Real tol;
   JX_Real setup_time;
   JX_Real solve_time;
   JX_Real max_row_sum;
   JX_Real strong_threshold;

   JX_Int max_levels;
   JX_Int cycle_type;
   JX_Int relax_type;
   JX_Int relax_order;
   JX_Int ns_down;
   JX_Int ns_up;
   JX_Int ns_coarse;
   JX_Int coarsen_type;
   JX_Int interp_type;
   JX_Int P_max_elmts;
   JX_Int agg_num_levels;
   JX_Int coarse_threshold;
   JX_Int amg_print_level;

   JX_Int k_dim;
   JX_Int TTest;
   JX_Int num_sl;
   JX_Int max_iter;
   JX_Int two_norm;
   JX_Int can_reuse;
   JX_Int num_reuse;
   JX_Int solver_id;
   JX_Int num_sl_amg;
   JX_Int print_level;
   JX_Int max_iter_sl;
   JX_Int max_iter_reuse;
   JX_Int current_update;
   JX_Int current_iter_sl;
   JX_Int current_num_lvls;
   JX_Int total_num_levels;
   JX_Int is_check_restarted;
   JX_Int current_sl_converge;
   JX_Int total_redun_iter_sl;
   JX_Int total_converge_iter;
   JX_Int total_redun_iter_amg;

   JX_Int *iter_sl_lvls;

   JX_hpCSRMatrix pre_mat;

   JX_Solver solver;
   JX_Solver amg_solver;
   JX_Solver jac_solver;
   
} jx_ParAdpSetupAMGData;

#define jx_ParAdpSetupAMGDataComm(amg_data)              ((amg_data)->comm)
#define jx_ParAdpSetupAMGDataIsAdaptive(amg_data)        ((amg_data)->is_adaptive)
#define jx_ParAdpSetupAMGDataIsAdapSetup(amg_data)       ((amg_data)->is_adap_setup)
#define jx_ParAdpSetupAMGDataTol(amg_data)               ((amg_data)->tol)
#define jx_ParAdpSetupAMGDataSetupTime(amg_data)         ((amg_data)->setup_time)
#define jx_ParAdpSetupAMGDataSolveTime(amg_data)         ((amg_data)->solve_time)
#define jx_ParAdpSetupAMGDataMaxRowSum(amg_data)         ((amg_data)->max_row_sum)
#define jx_ParAdpSetupAMGDataStrongThreshold(amg_data)   ((amg_data)->strong_threshold)
#define jx_ParAdpSetupAMGDataMaxLevels(amg_data)         ((amg_data)->max_levels)
#define jx_ParAdpSetupAMGDataCycleType(amg_data)         ((amg_data)->cycle_type)
#define jx_ParAdpSetupAMGDataRelaxType(amg_data)         ((amg_data)->relax_type)
#define jx_ParAdpSetupAMGDataRelaxOrder(amg_data)        ((amg_data)->relax_order)
#define jx_ParAdpSetupAMGDataNsDown(amg_data)            ((amg_data)->ns_down)
#define jx_ParAdpSetupAMGDataNsUp(amg_data)              ((amg_data)->ns_up)
#define jx_ParAdpSetupAMGDataNsCoarse(amg_data)          ((amg_data)->ns_coarse)
#define jx_ParAdpSetupAMGDataCoarsenType(amg_data)       ((amg_data)->coarsen_type)
#define jx_ParAdpSetupAMGDataInterpType(amg_data)        ((amg_data)->interp_type)
#define jx_ParAdpSetupAMGDataPMaxElmts(amg_data)         ((amg_data)->P_max_elmts)
#define jx_ParAdpSetupAMGDataAggNumLevels(amg_data)      ((amg_data)->agg_num_levels)
#define jx_ParAdpSetupAMGDataCoarseThreshold(amg_data)   ((amg_data)->coarse_threshold)
#define jx_ParAdpSetupAMGDataAMGPrintLevel(amg_data)     ((amg_data)->amg_print_level)
#define jx_ParAdpSetupAMGDataKDim(amg_data)              ((amg_data)->k_dim)
#define jx_ParAdpSetupAMGDataTTest(amg_data)             ((amg_data)->TTest)
#define jx_ParAdpSetupAMGDataNumSL(amg_data)             ((amg_data)->num_sl)
#define jx_ParAdpSetupAMGDataMaxIter(amg_data)           ((amg_data)->max_iter)
#define jx_ParAdpSetupAMGDataTwoNorm(amg_data)           ((amg_data)->two_norm)
#define jx_ParAdpSetupAMGDataCanReUse(amg_data)          ((amg_data)->can_reuse)
#define jx_ParAdpSetupAMGDataNumReUse(amg_data)          ((amg_data)->num_reuse)
#define jx_ParAdpSetupAMGDataSolverID(amg_data)          ((amg_data)->solver_id)
#define jx_ParAdpSetupAMGDataNumSLAMG(amg_data)          ((amg_data)->num_sl_amg)
#define jx_ParAdpSetupAMGDataPrintLevel(amg_data)        ((amg_data)->print_level)
#define jx_ParAdpSetupAMGDataMaxIterSL(amg_data)         ((amg_data)->max_iter_sl)
#define jx_ParAdpSetupAMGDataMaxIterReUse(amg_data)      ((amg_data)->max_iter_reuse)
#define jx_ParAdpSetupAMGDataCurrentUpdate(amg_data)     ((amg_data)->current_update)
#define jx_ParAdpSetupAMGDataCurrentIterSL(amg_data)     ((amg_data)->current_iter_sl)
#define jx_ParAdpSetupAMGDataCurrentNumLvls(amg_data)    ((amg_data)->current_num_lvls)
#define jx_ParAdpSetupAMGDataTotalNumLevels(amg_data)    ((amg_data)->total_num_levels)
#define jx_ParAdpSetupAMGDataIsCheckRestarted(amg_data)  ((amg_data)->is_check_restarted)
#define jx_ParAdpSetupAMGDataCurrentSLConverge(amg_data) ((amg_data)->current_sl_converge)
#define jx_ParAdpSetupAMGDataTotalRedunIterSL(amg_data)  ((amg_data)->total_redun_iter_sl)
#define jx_ParAdpSetupAMGDataTotalConvergeIter(amg_data) ((amg_data)->total_converge_iter)
#define jx_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) ((amg_data)->total_redun_iter_amg)
#define jx_ParAdpSetupAMGDataIterSLLvls(amg_data)        ((amg_data)->iter_sl_lvls)
#define jx_ParAdpSetupAMGDataPreMat(amg_data)            ((amg_data)->pre_mat)
#define jx_ParAdpSetupAMGDataSolver(amg_data)            ((amg_data)->solver)
#define jx_ParAdpSetupAMGDataAMGSolver(amg_data)         ((amg_data)->amg_solver)
#define jx_ParAdpSetupAMGDataJACSolver(amg_data)         ((amg_data)->jac_solver)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/amg/par_asetup.c */
JX_Int JX_ParAdpSetupAMGCreate( MPI_Comm comm, JX_Solver *solver );
JX_Int JX_ParAdpSetupAMGDestroy( JX_Solver solver );
JX_Int JX_ParAdpSetupAMGSetup( JX_Solver solver, JX_Matrix matA, JX_Vector vecB, JX_Vector vecX );
JX_Int JX_ParAdpSetupAMGSolve( JX_Solver solver, JX_Matrix preA, JX_Matrix matA, JX_Vector vecB, JX_Vector vecX );
JX_Int JX_ParAdpSetupAMGSetIsAdaptive( JX_Solver solver, JX_Bool is_adaptive );
JX_Int JX_ParAdpSetupAMGSetIsAdapSetup( JX_Solver solver, JX_Bool is_adap_setup );
JX_Int JX_ParAdpSetupAMGSetTol( JX_Solver solver, JX_Real tol );
JX_Int JX_ParAdpSetupAMGSetMaxRowSum( JX_Solver solver, JX_Real max_row_sum );
JX_Int JX_ParAdpSetupAMGSetStrongThreshold( JX_Solver solver, JX_Real strong_threshold );
JX_Int JX_ParAdpSetupAMGSetMaxLevels( JX_Solver solver, JX_Int  max_levels );
JX_Int JX_ParAdpSetupAMGSetCycleType( JX_Solver solver,JX_Int cycle_type );
JX_Int JX_ParAdpSetupAMGSetRelaxType( JX_Solver solver, JX_Int relax_type );
JX_Int JX_ParAdpSetupAMGSetRelaxOrder( JX_Solver solver, JX_Int relax_order );
JX_Int JX_ParAdpSetupAMGSetNsDown( JX_Solver solver, JX_Int ns_down );
JX_Int JX_ParAdpSetupAMGSetNsUp( JX_Solver solver, JX_Int ns_up );
JX_Int JX_ParAdpSetupAMGSetNsCoarse( JX_Solver solver, JX_Int ns_coarse );
JX_Int JX_ParAdpSetupAMGSetCoarsenType( JX_Solver solver, JX_Int coarsen_type );
JX_Int JX_ParAdpSetupAMGSetInterpType( JX_Solver solver, JX_Int interp_type );
JX_Int JX_ParAdpSetupAMGSetPMaxElmts( JX_Solver solver, JX_Int P_max_elmts );
JX_Int JX_ParAdpSetupAMGSetAggNumLevels( JX_Solver solver, JX_Int agg_num_levels );
JX_Int JX_ParAdpSetupAMGSetCoarseThreshold( JX_Solver solver, JX_Int coarse_threshold );
JX_Int JX_ParAdpSetupAMGSetAMGPrintLevel( JX_Solver solver, JX_Int amg_print_level );
JX_Int JX_ParAdpSetupAMGSetKDim( JX_Solver solver, JX_Int k_dim );
JX_Int JX_ParAdpSetupAMGSetTTest( JX_Solver solver, JX_Int TTest );
JX_Int JX_ParAdpSetupAMGSetMaxIter( JX_Solver solver, JX_Int max_iter );
JX_Int JX_ParAdpSetupAMGSetTwoNorm( JX_Solver solver, JX_Int two_norm );
JX_Int JX_ParAdpSetupAMGSetSolverID( JX_Solver solver, JX_Int solver_id );
JX_Int JX_ParAdpSetupAMGSetPrintLevel( JX_Solver solver, JX_Int print_level );
JX_Int JX_ParAdpSetupAMGSetMaxIterReUse( JX_Solver solver, JX_Int max_iter_reuse );
JX_Int JX_ParAdpSetupAMGSetIsCheckRestarted( JX_Solver solver, JX_Int is_check_restarted );
void *jx_ParAdpSetupAMGCreate( MPI_Comm comm );
JX_Int jx_ParAdpSetupAMGDestroy( void *data );
JX_Int jx_ParAdpSetupAMGSetup( void *data, void *matA, void *vecB, void *vecX );
JX_Int jx_ParAdpSetupAMGSolve( void *data, void *preA, void *matA, void *vecB, void *vecX );
JX_Int jx_ParAdpSetupAMGSetIsAdaptive( void *data, JX_Bool is_adaptive );
JX_Int jx_ParAdpSetupAMGSetIsAdapSetup( void *data, JX_Bool is_adap_setup );
JX_Int jx_ParAdpSetupAMGSetTol( void *data, JX_Real tol );
JX_Int jx_ParAdpSetupAMGSetMaxRowSum( void *data, JX_Real max_row_sum );
JX_Int jx_ParAdpSetupAMGSetStrongThreshold( void *data, JX_Real strong_threshold );
JX_Int jx_ParAdpSetupAMGSetMaxLevels( void *data, JX_Int max_levels );
JX_Int jx_ParAdpSetupAMGSetCycleType( void *data, JX_Int cycle_type );
JX_Int jx_ParAdpSetupAMGSetRelaxType( void *data, JX_Int relax_type );
JX_Int jx_ParAdpSetupAMGSetRelaxOrder( void *data, JX_Int relax_order );
JX_Int jx_ParAdpSetupAMGSetNsDown( void *data, JX_Int ns_down );
JX_Int jx_ParAdpSetupAMGSetNsUp( void *data, JX_Int ns_up );
JX_Int jx_ParAdpSetupAMGSetNsCoarse( void *data, JX_Int ns_coarse );
JX_Int jx_ParAdpSetupAMGSetCoarsenType( void *data, JX_Int coarsen_type );
JX_Int jx_ParAdpSetupAMGSetInterpType( void *data, JX_Int interp_type );
JX_Int jx_ParAdpSetupAMGSetPMaxElmts( void *data, JX_Int P_max_elmts );
JX_Int jx_ParAdpSetupAMGSetAggNumLevels( void *data, JX_Int agg_num_levels );
JX_Int jx_ParAdpSetupAMGSetCoarseThreshold( void *data, JX_Int coarse_threshold );
JX_Int jx_ParAdpSetupAMGSetAMGPrintLevel( void *data, JX_Int amg_print_level );
JX_Int jx_ParAdpSetupAMGSetKDim( void *data, JX_Int k_dim );
JX_Int jx_ParAdpSetupAMGSetTTest( void *data, JX_Int TTest );
JX_Int jx_ParAdpSetupAMGSetMaxIter( void *data, JX_Int max_iter );
JX_Int jx_ParAdpSetupAMGSetTwoNorm( void *data, JX_Int two_norm );
JX_Int jx_ParAdpSetupAMGSetSolverID( void *data, JX_Int solver_id );
JX_Int jx_ParAdpSetupAMGSetPrintLevel( void *data, JX_Int print_level );
JX_Int JX_ParAdpSetupAMGSetMaxIterSL( JX_Solver solver, JX_Int max_iter_sl );
JX_Int jx_ParAdpSetupAMGSetMaxIterSL( void *data, JX_Int max_iter_sl );
JX_Int jx_ParAdpSetupAMGSetMaxIterReUse( void *data, JX_Int max_iter_reuse );
JX_Int jx_ParAdpSetupAMGSetIsCheckRestarted( void *data, JX_Int is_check_restarted );
JX_Int JX_JacobiCreate( JX_Solver *solver );
void *jx_JacobiCreate();
JX_Int JX_JacobiDestroy( JX_Solver solver );
JX_Int jx_JacobiDestroy( void *data );
JX_Int JX_JacobiSetMaxIter( JX_Solver solver, JX_Int max_iter );
JX_Int jx_JacobiSetMaxIter( void *data, JX_Int max_iter );
JX_Int JX_JacobiSetup( JX_Solver solver, JX_hpCSRMatrix par_matrix );
JX_Int jx_JacobiSetup( void *solver, jx_ParCSRMatrix *par_matrix );
JX_Int JX_JacobiPrecond( JX_Solver solver, JX_hpCSRMatrix par_matrix, JX_ParVector par_rhs, JX_ParVector par_app );
JX_Int jx_JacobiPrecond( void *solver, jx_ParCSRMatrix *par_matrix, jx_ParVector *par_rhs, jx_ParVector *par_app );
JX_Int JX_PAMGAdapSetup( JX_Solver solver, JX_hpCSRMatrix matA, JX_ParVector vecB, JX_ParVector vecX,
                         JX_Int sid, JX_Int pl, JX_Real tol, JX_Int k_dim, JX_Int two_norm, JX_Int is_check_restarted,
                         JX_Bool is_adaptive, JX_Int max_iter_sl, JX_Int *iter_sl_used, JX_Int *num_sl, JX_Int *iter_sl_lvls );
JX_Int jx_PAMGAdapSetup( void *amg_vdata, jx_hpCSRMatrix *par_matrix, jx_ParVector *par_rhs, jx_ParVector *par_app,
                         JX_Int solver_id, JX_Int yprint_level, JX_Real ytol, JX_Int yk_dim, JX_Int ytwo_norm, JX_Int yis_check_restarted,
                         JX_Bool yis_adaptive, JX_Int ymax_iter_sl, JX_Int *yiter_sl_used, JX_Int *ynum_sl, JX_Int *iter_sl_lvls );
JX_Int JX_ParAdpSetupAMGGetNumSL( JX_Solver solver, JX_Int *num_sl );
JX_Int jx_ParAdpSetupAMGGetNumSL( void *data, JX_Int *num_sl );
JX_Int JX_ParAdpSetupAMGGetNumReUse( JX_Solver solver, JX_Int *num_reuse );
JX_Int jx_ParAdpSetupAMGGetNumReUse( void *data, JX_Int *num_reuse );
JX_Int JX_ParAdpSetupAMGGetCanReUse( JX_Solver solver, JX_Int *can_reuse );
JX_Int jx_ParAdpSetupAMGGetCanReUse( void *data, JX_Int *can_reuse );
JX_Int JX_ParAdpSetupAMGGetNumSLAMG( JX_Solver solver, JX_Int *num_sl_amg );
JX_Int jx_ParAdpSetupAMGGetNumSLAMG( void *data, JX_Int *num_sl_amg );
JX_Int JX_ParAdpSetupAMGGetCurrentSL( JX_Solver solver, JX_Int *iter_sl, JX_Int *converge, JX_Int *num_lvls, JX_Int **iter_sl_lvls, JX_Int *update );
JX_Int jx_ParAdpSetupAMGGetCurrentSL( void *data, JX_Int *iter_sl, JX_Int *converge, JX_Int *num_lvls, JX_Int **iter_sl_lvls, JX_Int *update );
JX_Int JX_ParAdpSetupAMGGetTotalNumLevels( JX_Solver solver, JX_Int *total_num_levels );
JX_Int jx_ParAdpSetupAMGGetTotalNumLevels( void *data, JX_Int *total_num_levels );
JX_Int JX_ParAdpSetupAMGGetTotalRedunIterSL( JX_Solver solver, JX_Int *total_redun_iter_sl );
JX_Int jx_ParAdpSetupAMGGetTotalRedunIterSL( void *data, JX_Int *total_redun_iter_sl );
JX_Int JX_ParAdpSetupAMGGetTotalConvergeIter( JX_Solver solver, JX_Int *total_converge_iter );
JX_Int jx_ParAdpSetupAMGGetTotalConvergeIter( void *data, JX_Int *total_converge_iter );
JX_Int JX_ParAdpSetupAMGGetTotalRedunIterAMG( JX_Solver solver, JX_Int *total_redun_iter_amg );
JX_Int jx_ParAdpSetupAMGGetTotalRedunIterAMG( void *data, JX_Int *total_redun_iter_amg );
JX_Int JX_ParAdpSetupAMGGetSetupTime( JX_Solver solver, JX_Real *setup_time );
JX_Int jx_ParAdpSetupAMGGetSetupTime( void *data, JX_Real *setup_time );
JX_Int JX_ParAdpSetupAMGGetSolveTime( JX_Solver solver, JX_Real *solve_time );
JX_Int jx_ParAdpSetupAMGGetSolveTime( void *data, JX_Real *solve_time );
JX_Int JX_ParAdpSetupAMGGetStatus( JX_Solver solver );
JX_Int jx_ParAdpSetupAMGGetStatus( void *data );

#endif
