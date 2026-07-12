//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_asetup.h -- head files for parallel asetup-AMG solver
 *  Date: 2020/03/10
 */ 


#ifndef JXF_ASETUP_HEADER
#define JXF_ASETUP_HEADER

#ifndef JXF_PAMG_HEADER
#include "jxf_pamg.h"
#endif
 
#ifndef JXF_KRYLOV_HEADER
#include "jxf_krylov.h"
#endif

/*----------------------------------------------------------------*
*                   Struct Declaration                           *
*----------------------------------------------------------------*/

/*!
 * \struct jxf_JacobiData
 */
typedef struct
{
   JXF_Int max_iter;

   jxf_ParVector *tmp_vec;
    
} jxf_JacobiData;

#define jxf_JacobiDataMaxIter(amg_data)  ((amg_data)->max_iter)
#define jxf_JacobiDataTmpVec(amg_data)   ((amg_data)->tmp_vec)

/*!
 * \struct jxf_ParAdpSetupAMGData
 */
typedef struct
{
   MPI_Comm comm;

   JXF_Bool is_adaptive;
   JXF_Bool is_adap_setup;

   JXF_Real tol;
   JXF_Real setup_time;
   JXF_Real solve_time;
   JXF_Real max_row_sum;
   JXF_Real strong_threshold;

   JXF_Int max_levels;
   JXF_Int cycle_type;
   JXF_Int relax_type;
   JXF_Int relax_order;
   JXF_Int ns_down;
   JXF_Int ns_up;
   JXF_Int ns_coarse;
   JXF_Int coarsen_type;
   JXF_Int interp_type;
   JXF_Int P_max_elmts;
   JXF_Int agg_num_levels;
   JXF_Int coarse_threshold;
   JXF_Int amg_print_level;

   JXF_Int k_dim;
   JXF_Int TTest;
   JXF_Int num_sl;
   JXF_Int max_iter;
   JXF_Int two_norm;
   JXF_Int can_reuse;
   JXF_Int num_reuse;
   JXF_Int solver_id;
   JXF_Int num_sl_amg;
   JXF_Int print_level;
   JXF_Int max_iter_sl;
   JXF_Int max_iter_reuse;
   JXF_Int current_update;
   JXF_Int current_iter_sl;
   JXF_Int current_num_lvls;
   JXF_Int total_num_levels;
   JXF_Int is_check_restarted;
   JXF_Int current_sl_converge;
   JXF_Int total_redun_iter_sl;
   JXF_Int total_converge_iter;
   JXF_Int total_redun_iter_amg;

   JXF_Int *iter_sl_lvls;

   JXF_hpCSRMatrix pre_mat;

   JXF_Solver solver;
   JXF_Solver amg_solver;
   JXF_Solver jac_solver;
   
} jxf_ParAdpSetupAMGData;

#define jxf_ParAdpSetupAMGDataComm(amg_data)              ((amg_data)->comm)
#define jxf_ParAdpSetupAMGDataIsAdaptive(amg_data)        ((amg_data)->is_adaptive)
#define jxf_ParAdpSetupAMGDataIsAdapSetup(amg_data)       ((amg_data)->is_adap_setup)
#define jxf_ParAdpSetupAMGDataTol(amg_data)               ((amg_data)->tol)
#define jxf_ParAdpSetupAMGDataSetupTime(amg_data)         ((amg_data)->setup_time)
#define jxf_ParAdpSetupAMGDataSolveTime(amg_data)         ((amg_data)->solve_time)
#define jxf_ParAdpSetupAMGDataMaxRowSum(amg_data)         ((amg_data)->max_row_sum)
#define jxf_ParAdpSetupAMGDataStrongThreshold(amg_data)   ((amg_data)->strong_threshold)
#define jxf_ParAdpSetupAMGDataMaxLevels(amg_data)         ((amg_data)->max_levels)
#define jxf_ParAdpSetupAMGDataCycleType(amg_data)         ((amg_data)->cycle_type)
#define jxf_ParAdpSetupAMGDataRelaxType(amg_data)         ((amg_data)->relax_type)
#define jxf_ParAdpSetupAMGDataRelaxOrder(amg_data)        ((amg_data)->relax_order)
#define jxf_ParAdpSetupAMGDataNsDown(amg_data)            ((amg_data)->ns_down)
#define jxf_ParAdpSetupAMGDataNsUp(amg_data)              ((amg_data)->ns_up)
#define jxf_ParAdpSetupAMGDataNsCoarse(amg_data)          ((amg_data)->ns_coarse)
#define jxf_ParAdpSetupAMGDataCoarsenType(amg_data)       ((amg_data)->coarsen_type)
#define jxf_ParAdpSetupAMGDataInterpType(amg_data)        ((amg_data)->interp_type)
#define jxf_ParAdpSetupAMGDataPMaxElmts(amg_data)         ((amg_data)->P_max_elmts)
#define jxf_ParAdpSetupAMGDataAggNumLevels(amg_data)      ((amg_data)->agg_num_levels)
#define jxf_ParAdpSetupAMGDataCoarseThreshold(amg_data)   ((amg_data)->coarse_threshold)
#define jxf_ParAdpSetupAMGDataAMGPrintLevel(amg_data)     ((amg_data)->amg_print_level)
#define jxf_ParAdpSetupAMGDataKDim(amg_data)              ((amg_data)->k_dim)
#define jxf_ParAdpSetupAMGDataTTest(amg_data)             ((amg_data)->TTest)
#define jxf_ParAdpSetupAMGDataNumSL(amg_data)             ((amg_data)->num_sl)
#define jxf_ParAdpSetupAMGDataMaxIter(amg_data)           ((amg_data)->max_iter)
#define jxf_ParAdpSetupAMGDataTwoNorm(amg_data)           ((amg_data)->two_norm)
#define jxf_ParAdpSetupAMGDataCanReUse(amg_data)          ((amg_data)->can_reuse)
#define jxf_ParAdpSetupAMGDataNumReUse(amg_data)          ((amg_data)->num_reuse)
#define jxf_ParAdpSetupAMGDataSolverID(amg_data)          ((amg_data)->solver_id)
#define jxf_ParAdpSetupAMGDataNumSLAMG(amg_data)          ((amg_data)->num_sl_amg)
#define jxf_ParAdpSetupAMGDataPrintLevel(amg_data)        ((amg_data)->print_level)
#define jxf_ParAdpSetupAMGDataMaxIterSL(amg_data)         ((amg_data)->max_iter_sl)
#define jxf_ParAdpSetupAMGDataMaxIterReUse(amg_data)      ((amg_data)->max_iter_reuse)
#define jxf_ParAdpSetupAMGDataCurrentUpdate(amg_data)     ((amg_data)->current_update)
#define jxf_ParAdpSetupAMGDataCurrentIterSL(amg_data)     ((amg_data)->current_iter_sl)
#define jxf_ParAdpSetupAMGDataCurrentNumLvls(amg_data)    ((amg_data)->current_num_lvls)
#define jxf_ParAdpSetupAMGDataTotalNumLevels(amg_data)    ((amg_data)->total_num_levels)
#define jxf_ParAdpSetupAMGDataIsCheckRestarted(amg_data)  ((amg_data)->is_check_restarted)
#define jxf_ParAdpSetupAMGDataCurrentSLConverge(amg_data) ((amg_data)->current_sl_converge)
#define jxf_ParAdpSetupAMGDataTotalRedunIterSL(amg_data)  ((amg_data)->total_redun_iter_sl)
#define jxf_ParAdpSetupAMGDataTotalConvergeIter(amg_data) ((amg_data)->total_converge_iter)
#define jxf_ParAdpSetupAMGDataTotalRedunIterAMG(amg_data) ((amg_data)->total_redun_iter_amg)
#define jxf_ParAdpSetupAMGDataIterSLLvls(amg_data)        ((amg_data)->iter_sl_lvls)
#define jxf_ParAdpSetupAMGDataPreMat(amg_data)            ((amg_data)->pre_mat)
#define jxf_ParAdpSetupAMGDataSolver(amg_data)            ((amg_data)->solver)
#define jxf_ParAdpSetupAMGDataAMGSolver(amg_data)         ((amg_data)->amg_solver)
#define jxf_ParAdpSetupAMGDataJACSolver(amg_data)         ((amg_data)->jac_solver)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/amg/par_asetup.c */
JXF_Int JXF_ParAdpSetupAMGCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int JXF_ParAdpSetupAMGDestroy( JXF_Solver solver );
JXF_Int JXF_ParAdpSetupAMGSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_ParAdpSetupAMGSolve( JXF_Solver solver, JXF_Matrix preA, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_ParAdpSetupAMGSetIsAdaptive( JXF_Solver solver, JXF_Bool is_adaptive );
JXF_Int JXF_ParAdpSetupAMGSetIsAdapSetup( JXF_Solver solver, JXF_Bool is_adap_setup );
JXF_Int JXF_ParAdpSetupAMGSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_ParAdpSetupAMGSetMaxRowSum( JXF_Solver solver, JXF_Real max_row_sum );
JXF_Int JXF_ParAdpSetupAMGSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold );
JXF_Int JXF_ParAdpSetupAMGSetMaxLevels( JXF_Solver solver, JXF_Int  max_levels );
JXF_Int JXF_ParAdpSetupAMGSetCycleType( JXF_Solver solver,JXF_Int cycle_type );
JXF_Int JXF_ParAdpSetupAMGSetRelaxType( JXF_Solver solver, JXF_Int relax_type );
JXF_Int JXF_ParAdpSetupAMGSetRelaxOrder( JXF_Solver solver, JXF_Int relax_order );
JXF_Int JXF_ParAdpSetupAMGSetNsDown( JXF_Solver solver, JXF_Int ns_down );
JXF_Int JXF_ParAdpSetupAMGSetNsUp( JXF_Solver solver, JXF_Int ns_up );
JXF_Int JXF_ParAdpSetupAMGSetNsCoarse( JXF_Solver solver, JXF_Int ns_coarse );
JXF_Int JXF_ParAdpSetupAMGSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type );
JXF_Int JXF_ParAdpSetupAMGSetInterpType( JXF_Solver solver, JXF_Int interp_type );
JXF_Int JXF_ParAdpSetupAMGSetPMaxElmts( JXF_Solver solver, JXF_Int P_max_elmts );
JXF_Int JXF_ParAdpSetupAMGSetAggNumLevels( JXF_Solver solver, JXF_Int agg_num_levels );
JXF_Int JXF_ParAdpSetupAMGSetCoarseThreshold( JXF_Solver solver, JXF_Int coarse_threshold );
JXF_Int JXF_ParAdpSetupAMGSetAMGPrintLevel( JXF_Solver solver, JXF_Int amg_print_level );
JXF_Int JXF_ParAdpSetupAMGSetKDim( JXF_Solver solver, JXF_Int k_dim );
JXF_Int JXF_ParAdpSetupAMGSetTTest( JXF_Solver solver, JXF_Int TTest );
JXF_Int JXF_ParAdpSetupAMGSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_ParAdpSetupAMGSetTwoNorm( JXF_Solver solver, JXF_Int two_norm );
JXF_Int JXF_ParAdpSetupAMGSetSolverID( JXF_Solver solver, JXF_Int solver_id );
JXF_Int JXF_ParAdpSetupAMGSetPrintLevel( JXF_Solver solver, JXF_Int print_level );
JXF_Int JXF_ParAdpSetupAMGSetMaxIterReUse( JXF_Solver solver, JXF_Int max_iter_reuse );
JXF_Int JXF_ParAdpSetupAMGSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted );
void *jxf_ParAdpSetupAMGCreate( MPI_Comm comm );
JXF_Int jxf_ParAdpSetupAMGDestroy( void *data );
JXF_Int jxf_ParAdpSetupAMGSetup( void *data, void *matA, void *vecB, void *vecX );
JXF_Int jxf_ParAdpSetupAMGSolve( void *data, void *preA, void *matA, void *vecB, void *vecX );
JXF_Int jxf_ParAdpSetupAMGSetIsAdaptive( void *data, JXF_Bool is_adaptive );
JXF_Int jxf_ParAdpSetupAMGSetIsAdapSetup( void *data, JXF_Bool is_adap_setup );
JXF_Int jxf_ParAdpSetupAMGSetTol( void *data, JXF_Real tol );
JXF_Int jxf_ParAdpSetupAMGSetMaxRowSum( void *data, JXF_Real max_row_sum );
JXF_Int jxf_ParAdpSetupAMGSetStrongThreshold( void *data, JXF_Real strong_threshold );
JXF_Int jxf_ParAdpSetupAMGSetMaxLevels( void *data, JXF_Int max_levels );
JXF_Int jxf_ParAdpSetupAMGSetCycleType( void *data, JXF_Int cycle_type );
JXF_Int jxf_ParAdpSetupAMGSetRelaxType( void *data, JXF_Int relax_type );
JXF_Int jxf_ParAdpSetupAMGSetRelaxOrder( void *data, JXF_Int relax_order );
JXF_Int jxf_ParAdpSetupAMGSetNsDown( void *data, JXF_Int ns_down );
JXF_Int jxf_ParAdpSetupAMGSetNsUp( void *data, JXF_Int ns_up );
JXF_Int jxf_ParAdpSetupAMGSetNsCoarse( void *data, JXF_Int ns_coarse );
JXF_Int jxf_ParAdpSetupAMGSetCoarsenType( void *data, JXF_Int coarsen_type );
JXF_Int jxf_ParAdpSetupAMGSetInterpType( void *data, JXF_Int interp_type );
JXF_Int jxf_ParAdpSetupAMGSetPMaxElmts( void *data, JXF_Int P_max_elmts );
JXF_Int jxf_ParAdpSetupAMGSetAggNumLevels( void *data, JXF_Int agg_num_levels );
JXF_Int jxf_ParAdpSetupAMGSetCoarseThreshold( void *data, JXF_Int coarse_threshold );
JXF_Int jxf_ParAdpSetupAMGSetAMGPrintLevel( void *data, JXF_Int amg_print_level );
JXF_Int jxf_ParAdpSetupAMGSetKDim( void *data, JXF_Int k_dim );
JXF_Int jxf_ParAdpSetupAMGSetTTest( void *data, JXF_Int TTest );
JXF_Int jxf_ParAdpSetupAMGSetMaxIter( void *data, JXF_Int max_iter );
JXF_Int jxf_ParAdpSetupAMGSetTwoNorm( void *data, JXF_Int two_norm );
JXF_Int jxf_ParAdpSetupAMGSetSolverID( void *data, JXF_Int solver_id );
JXF_Int jxf_ParAdpSetupAMGSetPrintLevel( void *data, JXF_Int print_level );
JXF_Int JXF_ParAdpSetupAMGSetMaxIterSL( JXF_Solver solver, JXF_Int max_iter_sl );
JXF_Int jxf_ParAdpSetupAMGSetMaxIterSL( void *data, JXF_Int max_iter_sl );
JXF_Int jxf_ParAdpSetupAMGSetMaxIterReUse( void *data, JXF_Int max_iter_reuse );
JXF_Int jxf_ParAdpSetupAMGSetIsCheckRestarted( void *data, JXF_Int is_check_restarted );
JXF_Int JXF_JacobiCreate( JXF_Solver *solver );
void *jxf_JacobiCreate();
JXF_Int JXF_JacobiDestroy( JXF_Solver solver );
JXF_Int jxf_JacobiDestroy( void *data );
JXF_Int JXF_JacobiSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int jxf_JacobiSetMaxIter( void *data, JXF_Int max_iter );
JXF_Int JXF_JacobiSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int jxf_JacobiSetup( void *solver, jxf_ParCSRMatrix *par_matrix );
JXF_Int JXF_JacobiPrecond( JXF_Solver solver, JXF_hpCSRMatrix par_matrix, JXF_ParVector par_rhs, JXF_ParVector par_app );
JXF_Int jxf_JacobiPrecond( void *solver, jxf_ParCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app );
JXF_Int JXF_PAMGAdapSetup( JXF_Solver solver, JXF_hpCSRMatrix matA, JXF_ParVector vecB, JXF_ParVector vecX,
                         JXF_Int sid, JXF_Int pl, JXF_Real tol, JXF_Int k_dim, JXF_Int two_norm, JXF_Int is_check_restarted,
                         JXF_Bool is_adaptive, JXF_Int max_iter_sl, JXF_Int *iter_sl_used, JXF_Int *num_sl, JXF_Int *iter_sl_lvls );
JXF_Int jxf_PAMGAdapSetup( void *amg_vdata, jxf_hpCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app,
                         JXF_Int solver_id, JXF_Int yprint_level, JXF_Real ytol, JXF_Int yk_dim, JXF_Int ytwo_norm, JXF_Int yis_check_restarted,
                         JXF_Bool yis_adaptive, JXF_Int ymax_iter_sl, JXF_Int *yiter_sl_used, JXF_Int *ynum_sl, JXF_Int *iter_sl_lvls );
JXF_Int JXF_ParAdpSetupAMGGetNumSL( JXF_Solver solver, JXF_Int *num_sl );
JXF_Int jxf_ParAdpSetupAMGGetNumSL( void *data, JXF_Int *num_sl );
JXF_Int JXF_ParAdpSetupAMGGetNumReUse( JXF_Solver solver, JXF_Int *num_reuse );
JXF_Int jxf_ParAdpSetupAMGGetNumReUse( void *data, JXF_Int *num_reuse );
JXF_Int JXF_ParAdpSetupAMGGetCanReUse( JXF_Solver solver, JXF_Int *can_reuse );
JXF_Int jxf_ParAdpSetupAMGGetCanReUse( void *data, JXF_Int *can_reuse );
JXF_Int JXF_ParAdpSetupAMGGetNumSLAMG( JXF_Solver solver, JXF_Int *num_sl_amg );
JXF_Int jxf_ParAdpSetupAMGGetNumSLAMG( void *data, JXF_Int *num_sl_amg );
JXF_Int JXF_ParAdpSetupAMGGetCurrentSL( JXF_Solver solver, JXF_Int *iter_sl, JXF_Int *converge, JXF_Int *num_lvls, JXF_Int **iter_sl_lvls, JXF_Int *update );
JXF_Int jxf_ParAdpSetupAMGGetCurrentSL( void *data, JXF_Int *iter_sl, JXF_Int *converge, JXF_Int *num_lvls, JXF_Int **iter_sl_lvls, JXF_Int *update );
JXF_Int JXF_ParAdpSetupAMGGetTotalNumLevels( JXF_Solver solver, JXF_Int *total_num_levels );
JXF_Int jxf_ParAdpSetupAMGGetTotalNumLevels( void *data, JXF_Int *total_num_levels );
JXF_Int JXF_ParAdpSetupAMGGetTotalRedunIterSL( JXF_Solver solver, JXF_Int *total_redun_iter_sl );
JXF_Int jxf_ParAdpSetupAMGGetTotalRedunIterSL( void *data, JXF_Int *total_redun_iter_sl );
JXF_Int JXF_ParAdpSetupAMGGetTotalConvergeIter( JXF_Solver solver, JXF_Int *total_converge_iter );
JXF_Int jxf_ParAdpSetupAMGGetTotalConvergeIter( void *data, JXF_Int *total_converge_iter );
JXF_Int JXF_ParAdpSetupAMGGetTotalRedunIterAMG( JXF_Solver solver, JXF_Int *total_redun_iter_amg );
JXF_Int jxf_ParAdpSetupAMGGetTotalRedunIterAMG( void *data, JXF_Int *total_redun_iter_amg );
JXF_Int JXF_ParAdpSetupAMGGetSetupTime( JXF_Solver solver, JXF_Real *setup_time );
JXF_Int jxf_ParAdpSetupAMGGetSetupTime( void *data, JXF_Real *setup_time );
JXF_Int JXF_ParAdpSetupAMGGetSolveTime( JXF_Solver solver, JXF_Real *solve_time );
JXF_Int jxf_ParAdpSetupAMGGetSolveTime( void *data, JXF_Real *solve_time );
JXF_Int JXF_ParAdpSetupAMGGetStatus( JXF_Solver solver );
JXF_Int jxf_ParAdpSetupAMGGetStatus( void *data );

#endif
