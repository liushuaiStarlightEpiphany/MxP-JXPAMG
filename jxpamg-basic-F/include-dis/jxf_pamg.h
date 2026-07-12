//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_pamg.h -- head files for parallel amg solver
 *  Date: 2011/09/08
 *
 *  Created by peghoty
 */

#ifndef JXF_PAMG_HEADER
#define JXF_PAMG_HEADER

#ifndef JXF_UTIL_HEADER
#include "jxf_util.h"
#endif

#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                      Macro Definition                          *
 *----------------------------------------------------------------*/

#define LIST_HEAD -1
#define LIST_TAIL -2

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_double_linked_list
 */
struct jxf_double_linked_list
{
   JXF_Int data;
   struct jxf_double_linked_list *next_elt;
   struct jxf_double_linked_list *prev_elt;
   JXF_Int head;
   JXF_Int tail;
};

typedef struct jxf_double_linked_list jxf_ListElement;
typedef jxf_ListElement *jxf_LinkList;

/*!
 * \struct jxf_ParAMGData
 */
#define CUMNUMIT
typedef struct
{
   /* setup params */
   JXF_Int max_levels;
   JXF_Real strong_threshold;
   JXF_Real max_row_sum;
   JXF_Real trunc_factor;
   JXF_Real jacobi_trunc_threshold;
   JXF_Real S_commpkg_switch;
   JXF_Real CR_rate;
   JXF_Real CR_strong_th;
   JXF_Int measure_type;
   JXF_Int setup_type;
   JXF_Int coarsen_type;
   JXF_Int P_max_elmts;
   JXF_Int interp_type;
   JXF_Int restr_par;
   JXF_Real AIR_strong_th;
   JXF_Int agg_num_levels;
   JXF_Int agg_interp_type;
   JXF_Int agg_P_max_elmts;
   JXF_Int agg_P12_max_elmts;
   JXF_Real agg_trunc_factor;
   JXF_Real agg_P12_trunc_factor;
   JXF_Int num_paths;
   JXF_Int post_interp_type;
   JXF_Int num_CR_relax_steps;
   JXF_Int IS_type;
   JXF_Int CR_use_CG;
   JXF_Int cgc_its;
   JXF_Int spmt_rap_type;          /* RAP Type for single processor with multi-threads Yue Xiaoqiang 2012/10/12 */
   JXF_Int wall_time_option;       // run-time of AMG component, Yue Xiaoqiang 2015/09/30
   JXF_Int ai_measure_type;        // Algebraic interface, Yue Xiaoqiang 2014/02/26
   JXF_Int ai_relax_type;          // AI-smoothing, Yue Xiaoqiang 2014/07/06
   JXF_Int measure_type_rlx;       // peghoty 2010/05/29
   JXF_Int number_syn_rlx;         // peghoty 2010/05/29
   JXF_Real measure_threshold_rlx; // peghoty 2010/05/29

   /* solve params */
   JXF_Int max_iter;
   JXF_Int min_iter;
   JXF_Int cycle_type;
   JXF_Int *num_grid_sweeps;
   JXF_Int *grid_relax_type;
   JXF_Int **grid_relax_points;
   JXF_Int relax_order;
   JXF_Int user_relax_type;
   JXF_Int user_coarse_relax_type;
   JXF_Real *relax_weight;
   JXF_Real *omega;
   JXF_Real tol;
   JXF_Real rhsnrm_threshold; // for the V-cycle control, the following strategy is used
                             // when rhsnrm_threshold != 0:
                             // 1. if ||rhs|| > rhsnrm_threshold, use relative norm control
                             // 2. otherwise, use absolute norm control
                             // default of rhsnrm_threshold: 0.0
   /* problem data */
   jxf_ParCSRMatrix *A;
   JXF_Int num_variables;
   JXF_Int num_functions;
   JXF_Int nodal;
   JXF_Int nodal_diag;
   JXF_Int num_points;
   JXF_Int *dof_func;
   JXF_Int *dof_point;
   JXF_Int *point_dof_map;

   /* data generated in the setup phase */
   jxf_hpCSRMatrix **hp_A_array;
   jxf_ParVector **F_array;
   jxf_ParVector **U_array;
   jxf_ParCSRMatrix **P_array;
   jxf_ParCSRMatrix **R_array;
   JXF_Real **AI_measure_array;
   JXF_Int **CF_marker_array;
   JXF_Int **dof_func_array;
   JXF_Int **dof_point_array;
   JXF_Int **point_dof_map_array;
   JXF_Int num_levels;
   JXF_Real **l1_norms;
   JXF_Int *AIR_maxsize_ls;

   /* Block data */
   jxf_ParCSRBlockMatrix **A_block_array;
   jxf_ParCSRBlockMatrix **P_block_array;
   jxf_ParCSRBlockMatrix **R_block_array;

   JXF_Int block_mode;

   /* data for more complex smoothers */
   JXF_Int smooth_num_levels;
   JXF_Int smooth_type;
   JXF_Solver *smoother;
   JXF_Int smooth_num_sweeps;
   JXF_Int schw_variant;
   JXF_Int schw_overlap;
   JXF_Int schw_domain_type;
   JXF_Real schwarz_rlx_weight;
   JXF_Int schwarz_use_nonsymm;
   JXF_Int ps_sym;
   JXF_Int ps_level;
   JXF_Int pi_max_nz_per_row;
   JXF_Int eu_level;
   JXF_Int eu_bj;
   JXF_Real ps_threshold;
   JXF_Real ps_filter;
   JXF_Real pi_drop_tol;
   JXF_Real eu_sparse_A;
   char *euclidfile;

   /* data generated in the solve phase */
   jxf_ParVector *Vtemp;
   jxf_Vector *Vtemp_local;
   JXF_Real *Vtemp_local_data;
   JXF_Real cycle_op_count;
   jxf_ParVector *Rtemp;
   jxf_ParVector *Ptemp;
   jxf_ParVector *Ztemp;

   /* fields used by GSMG and LS interpolation */
   JXF_Int gsmg;        /* nonzero indicates use of GSMG */
   JXF_Int num_samples; /* number of sample vectors */

   /* log info */
   JXF_Int logging;
   JXF_Int num_iterations;
#ifdef CUMNUMIT
   JXF_Int cum_num_iterations;
#endif
   JXF_Real rel_resid_norm;
   jxf_ParVector *residual; /* available if logging > 1 */

   /* output params */
   JXF_Int print_level;
   JXF_Int print_coarse_matrix;
   char log_file_name[256];
   JXF_Int debug_flag;

   /* whether to print the constructed coarse grids BM Oct 22, 2006 */
   JXF_Int plot_grids;
   char plot_filename[251];

   /* coordinate data BM Oct 17, 2006 */
   JXF_Int coorddim;
   float *coordinates;

   /* peghoty 2009/07/27 */
   JXF_Int conv_criteria;
   JXF_Int coarsestsolverid;

   /* peghoty 2010/04/14 */
   JXF_Int coarse_threshold;   // the allowed number of gridpoints in the coarsest level
   JXF_Real coarse_ratio;      // when coarse_size >= fine_size*coarse_ratio,
                              // use CLJP instead of the current coarsening strategy.
   JXF_Real convfac_threshold; // if the current convergence factor > convfac_threshold,
                              // stop the iteration

   /* Use 2 mat-mat-muls instead of triple product*/
   JXF_Int rap2;

   JXF_Int keepTranspose;

} jxf_ParAMGData;

/* setup params */
#define jxf_ParAMGDataRestriction(amg_data) ((amg_data)->restr_par)
#define jxf_ParAMGDataAIRStrongTh(amg_data) ((amg_data)->AIR_strong_th)
#define jxf_ParAMGDataMaxLevels(amg_data) ((amg_data)->max_levels)
#define jxf_ParAMGDataStrongThreshold(amg_data) ((amg_data)->strong_threshold)
#define jxf_ParAMGDataMaxRowSum(amg_data) ((amg_data)->max_row_sum)
#define jxf_ParAMGDataTruncFactor(amg_data) ((amg_data)->trunc_factor)
#define jxf_ParAMGDataJacobiTruncThreshold(amg_data) ((amg_data)->jacobi_trunc_threshold)
#define jxf_ParAMGDataSCommPkgSwitch(amg_data) ((amg_data)->S_commpkg_switch)
#define jxf_ParAMGDataInterpType(amg_data) ((amg_data)->interp_type)
#define jxf_ParAMGDataCoarsenType(amg_data) ((amg_data)->coarsen_type)
#define jxf_ParAMGDataMeasureType(amg_data) ((amg_data)->measure_type)
#define jxf_ParAMGDataSetupType(amg_data) ((amg_data)->setup_type)
#define jxf_ParAMGDataPMaxElmts(amg_data) ((amg_data)->P_max_elmts)
#define jxf_ParAMGDataNumPaths(amg_data) ((amg_data)->num_paths)
#define jxf_ParAMGDataAggNumLevels(amg_data) ((amg_data)->agg_num_levels)
#define jxf_ParAMGDataAggInterpType(amg_data) ((amg_data)->agg_interp_type)
#define jxf_ParAMGDataAggPMaxElmts(amg_data) ((amg_data)->agg_P_max_elmts)
#define jxf_ParAMGDataAggP12MaxElmts(amg_data) ((amg_data)->agg_P12_max_elmts)
#define jxf_ParAMGDataAggTruncFactor(amg_data) ((amg_data)->agg_trunc_factor)
#define jxf_ParAMGDataAggP12TruncFactor(amg_data) ((amg_data)->agg_P12_trunc_factor)
#define jxf_ParAMGDataPostInterpType(amg_data) ((amg_data)->post_interp_type)
#define jxf_ParAMGDataNumCRRelaxSteps(amg_data) ((amg_data)->num_CR_relax_steps)
#define jxf_ParAMGDataCRRate(amg_data) ((amg_data)->CR_rate)
#define jxf_ParAMGDataCRStrongTh(amg_data) ((amg_data)->CR_strong_th)
#define jxf_ParAMGDataISType(amg_data) ((amg_data)->IS_type)
#define jxf_ParAMGDataCRUseCG(amg_data) ((amg_data)->CR_use_CG)
#define jxf_ParAMGDataL1Norms(amg_data) ((amg_data)->l1_norms)
#define jxf_ParAMGDataCGCIts(amg_data) ((amg_data)->cgc_its)
#define jxf_ParAMGDataSpMtRapType(amg_data) ((amg_data)->spmt_rap_type)                 // Yue Xiaoqiang 2012/10/13
#define jxf_ParAMGDataWallTimeOption(amg_data) ((amg_data)->wall_time_option)           // Yue Xiaoqiang 2015/09/30
#define jxf_ParAMGDataAIMeasureType(amg_data) ((amg_data)->ai_measure_type)             // Yue Xiaoqiang 2014/02/26
#define jxf_ParAMGDataAIRelaxType(amg_data) ((amg_data)->ai_relax_type)                 // Yue Xiaoqiang 2014/07/06
#define jxf_ParAMGDataMeasureTypeRlx(amg_data) ((amg_data)->measure_type_rlx)           // peghoty 2010/05/29
#define jxf_ParAMGDataNumberSynRlx(amg_data) ((amg_data)->number_syn_rlx)               // peghoty 2010/05/29
#define jxf_ParAMGDataMeasureThresholdRlx(amg_data) ((amg_data)->measure_threshold_rlx) // peghoty 2010/05/29
#define jxf_ParAMGDataAIRMaxSizeLS(amg_data) ((amg_data)->AIR_maxsize_ls)               // Yue Xiaoqiang 2018/11/25
/* solve params */
#define jxf_ParAMGDataMinIter(amg_data) ((amg_data)->min_iter)
#define jxf_ParAMGDataMaxIter(amg_data) ((amg_data)->max_iter)
#define jxf_ParAMGDataCycleType(amg_data) ((amg_data)->cycle_type)
#define jxf_ParAMGDataNumGridSweeps(amg_data) ((amg_data)->num_grid_sweeps)
#define jxf_ParAMGDataUserCoarseRelaxType(amg_data) ((amg_data)->user_coarse_relax_type)
#define jxf_ParAMGDataGridRelaxType(amg_data) ((amg_data)->grid_relax_type)
#define jxf_ParAMGDataGridRelaxPoints(amg_data) ((amg_data)->grid_relax_points)
#define jxf_ParAMGDataRelaxOrder(amg_data) ((amg_data)->relax_order)
#define jxf_ParAMGDataUserRelaxType(amg_data) ((amg_data)->user_relax_type)
#define jxf_ParAMGDataRelaxWeight(amg_data) ((amg_data)->relax_weight)
#define jxf_ParAMGDataOmega(amg_data) ((amg_data)->omega)
#define jxf_ParAMGDataTol(amg_data) ((amg_data)->tol)
#define jxf_ParAMGDataRhsNrmThreshold(amg_data) ((amg_data)->rhsnrm_threshold)
/* problem data parameters */
#define jxf_ParAMGDataNumVariables(amg_data) ((amg_data)->num_variables)
#define jxf_ParAMGDataNumFunctions(amg_data) ((amg_data)->num_functions)
#define jxf_ParAMGDataNodal(amg_data) ((amg_data)->nodal)
#define jxf_ParAMGDataNodalDiag(amg_data) ((amg_data)->nodal_diag)
#define jxf_ParAMGDataNumPoints(amg_data) ((amg_data)->num_points)
#define jxf_ParAMGDataDofFunc(amg_data) ((amg_data)->dof_func)
#define jxf_ParAMGDataDofPoint(amg_data) ((amg_data)->dof_point)
#define jxf_ParAMGDataPointDofMap(amg_data) ((amg_data)->point_dof_map)
/* data generated by the setup phase */
#define jxf_ParAMGDataAIMeasureArray(amg_data) ((amg_data)->AI_measure_array)
#define jxf_ParAMGDataCFMarkerArray(amg_data) ((amg_data)->CF_marker_array)
#define jxf_hpAMGDataAArray(amg_data) ((amg_data)->hp_A_array)
#define jxf_ParAMGDataFArray(amg_data) ((amg_data)->F_array)
#define jxf_ParAMGDataUArray(amg_data) ((amg_data)->U_array)
#define jxf_ParAMGDataPArray(amg_data) ((amg_data)->P_array)
#define jxf_ParAMGDataRArray(amg_data) ((amg_data)->R_array)
#define jxf_ParAMGDataDofFuncArray(amg_data) ((amg_data)->dof_func_array)
#define jxf_ParAMGDataDofPointArray(amg_data) ((amg_data)->dof_point_array)
#define jxf_ParAMGDataPointDofMapArray(amg_data) ((amg_data)->point_dof_map_array)
#define jxf_ParAMGDataNumLevels(amg_data) ((amg_data)->num_levels)
#define jxf_ParAMGDataSmoothType(amg_data) ((amg_data)->smooth_type)
#define jxf_ParAMGDataSmoothNumLevels(amg_data) ((amg_data)->smooth_num_levels)
#define jxf_ParAMGDataSmoothNumSweeps(amg_data) ((amg_data)->smooth_num_sweeps)
#define jxf_ParAMGDataSmoother(amg_data) ((amg_data)->smoother)
#define jxf_ParAMGDataVariant(amg_data) ((amg_data)->schw_variant)
#define jxf_ParAMGDataOverlap(amg_data) ((amg_data)->schw_overlap)
#define jxf_ParAMGDataDomainType(amg_data) ((amg_data)->schw_domain_type)
#define jxf_ParAMGDataSchwarzRlxWeight(amg_data) ((amg_data)->schwarz_rlx_weight)
#define jxf_ParAMGDataSchwarzUseNonSymm(amg_data) ((amg_data)->schwarz_use_nonsymm)
#define jxf_ParAMGDataSym(amg_data) ((amg_data)->ps_sym)
#define jxf_ParAMGDataLevel(amg_data) ((amg_data)->ps_level)
#define jxf_ParAMGDataMaxNzPerRow(amg_data) ((amg_data)->pi_max_nz_per_row)
#define jxf_ParAMGDataThreshold(amg_data) ((amg_data)->ps_threshold)
#define jxf_ParAMGDataFilter(amg_data) ((amg_data)->ps_filter)
#define jxf_ParAMGDataDropTol(amg_data) ((amg_data)->pi_drop_tol)
#define jxf_ParAMGDataEuclidFile(amg_data) ((amg_data)->euclidfile)
#define jxf_ParAMGDataEuLevel(amg_data) ((amg_data)->eu_level)
#define jxf_ParAMGDataEuSparseA(amg_data) ((amg_data)->eu_sparse_A)
#define jxf_ParAMGDataEuBJ(amg_data) ((amg_data)->eu_bj)
/* block */
#define jxf_ParAMGDataABlockArray(amg_data) ((amg_data)->A_block_array)
#define jxf_ParAMGDataPBlockArray(amg_data) ((amg_data)->P_block_array)
#define jxf_ParAMGDataRBlockArray(amg_data) ((amg_data)->R_block_array)
#define jxf_ParAMGDataBlockMode(amg_data) ((amg_data)->block_mode)
/* data generated in the solve phase */
#define jxf_ParAMGDataVtemp(amg_data) ((amg_data)->Vtemp)
#define jxf_ParAMGDataVtempLocal(amg_data) ((amg_data)->Vtemp_local)
#define jxf_ParAMGDataVtemplocalData(amg_data) ((amg_data)->Vtemp_local_data)
#define jxf_ParAMGDataCycleOpCount(amg_data) ((amg_data)->cycle_op_count)
#define jxf_ParAMGDataRtemp(amg_data) ((amg_data)->Rtemp)
#define jxf_ParAMGDataPtemp(amg_data) ((amg_data)->Ptemp)
#define jxf_ParAMGDataZtemp(amg_data) ((amg_data)->Ztemp)
/* fields used by GSMG */
#define jxf_ParAMGDataGSMG(amg_data) ((amg_data)->gsmg)
#define jxf_ParAMGDataNumSamples(amg_data) ((amg_data)->num_samples)
/* log info data */
#define jxf_ParAMGDataLogging(amg_data) ((amg_data)->logging)
#define jxf_ParAMGDataNumIterations(amg_data) ((amg_data)->num_iterations)
#ifdef CUMNUMIT
#define jxf_ParAMGDataCumNumIterations(amg_data) ((amg_data)->cum_num_iterations)
#endif
#define jxf_ParAMGDataRelativeResidualNorm(amg_data) ((amg_data)->rel_resid_norm)
#define jxf_ParAMGDataResidual(amg_data) ((amg_data)->residual)
/* output parameters */
#define jxf_ParAMGDataPrintLevel(amg_data) ((amg_data)->print_level)
#define jxf_ParAMGDataPrintCoarseSystem(amg_data) ((amg_data)->print_coarse_matrix)
#define jxf_ParAMGDataLogFileName(amg_data) ((amg_data)->log_file_name)
#define jxf_ParAMGDataDebugFlag(amg_data) ((amg_data)->debug_flag)
/* BM Oct 22, 2006 */
#define jxf_ParAMGDataPlotGrids(amg_data) ((amg_data)->plot_grids)
#define jxf_ParAMGDataPlotFileName(amg_data) ((amg_data)->plot_filename)
/* coordinates BM Oct 17, 2006 */
#define jxf_ParAMGDataCoordDim(amg_data) ((amg_data)->coorddim)
#define jxf_ParAMGDataCoordinates(amg_data) ((amg_data)->coordinates)
/* peghoty 2009/07/27 */
#define jxf_ParAMGDataConvCriteria(amg_data) ((amg_data)->conv_criteria)
#define jxf_ParAMGDataCoarsestSolverID(amg_data) ((amg_data)->coarsestsolverid)
/* peghoty 2010/04/14 */
#define jxf_ParAMGDataCoarseThreshold(amg_data) ((amg_data)->coarse_threshold)
#define jxf_ParAMGDataCoarseRatio(amg_data) ((amg_data)->coarse_ratio)

#define jxf_ParAMGDataConvFacThreshold(amg_data) ((amg_data)->convfac_threshold)

#define jxf_ParAMGDataRAP2(amg_data) ((amg_data)->rap2)
#define jxf_ParAMGDataKeepTranspose(amg_data) ((amg_data)->keepTranspose)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/amg/amg_linklist.c */
void jxf_dispose_elt(jxf_LinkList element_ptr);
void jxf_remove_point(jxf_LinkList *LoL_head_ptr,
                     jxf_LinkList *LoL_tail_ptr,
                     JXF_Int measure,
                     JXF_Int index,
                     JXF_Int *lists,
                     JXF_Int *where);
jxf_LinkList jxf_create_elt(JXF_Int Item);
void jxf_enter_on_lists(jxf_LinkList *LoL_head_ptr,
                       jxf_LinkList *LoL_tail_ptr,
                       JXF_Int measure,
                       JXF_Int index,
                       JXF_Int *lists,
                       JXF_Int *where);

/* csrc/amg/par_amg.c */

/*!
 * Create a solver or preconditioner object.
 */
JXF_Int JXF_PAMGCreate(JXF_Solver *solver);

/*!
 * Destroy a solver or preconditioner object.
 */
JXF_Int JXF_PAMGDestroy(JXF_Solver solver);
// JXF_Int JXF_hpPAMGDestroy( JXF_Solver solver);
/*!
 * Set up the PAMG solver or preconditioner.
 * If used as a preconditioner, this function should be passed
 * to the iterative solver SetPrecond function.
 * This page has been automatically generated with DOC++
 * Parameters:
 *    solver — [IN] object to be set up.
 *    par_matrix — [IN] ParCSR matrix used to construct the solver/preconditioner.
 */
// JXF_Int JXF_PAMGSetup( JXF_Solver solver, JXF_ParCSRMatrix  par_matrix );
JXF_Int JXF_PAMGSetup(JXF_Solver solver, JXF_hpCSRMatrix hp_matrix);
// JXF_Int JXF_hpPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix  hp_matrix );
/*!
 * Solve the system or apply PAMG as a preconditioner.
 * If used as a preconditioner, this function should be
 * passed to the iterative solver SetPrecond function.
 * Parameters:
 *  solver — [IN] solver or preconditioner object to be applied.
 *  par_matrix — [IN] ParCSR matrix, matrix of the linear system to be solved
 *  par_rhs — [IN] right hand side of the linear system to be solved
 *  par_app — [OUT] approximated solution of the linear system to be solved
 */
JXF_Int JXF_PAMGSolve(JXF_Solver solver,
                    JXF_hpCSRMatrix par_matrix,
                    JXF_ParVector par_rhs,
                    JXF_ParVector par_app);

/* peghoty 2011/09/04 */
JXF_Int
JXF_PAMGPrecond(JXF_Solver solver,
               JXF_hpCSRMatrix par_matrix,
               JXF_ParVector par_rhs,
               JXF_ParVector par_app);

JXF_Int JXF_PAMGSetRestriction(JXF_Solver solver, JXF_Int restr_par);
JXF_Int JXF_PAMGSetAIRStrongTh(JXF_Solver solver, JXF_Real AIR_strong_th);

/*!
 * Set maximum number of multigrid levels. The default is 25.
 */
JXF_Int JXF_PAMGSetMaxLevels(JXF_Solver solver, JXF_Int max_levels);
JXF_Int JXF_PAMGGetMaxLevels(JXF_Solver solver, JXF_Int *max_levels);

/*!
 * Set AMG strength threshold. The default is 0.25. For 2d Laplace operators,
 * 0.25 is a good value, for 3d Laplace operators, 0.5 or 0.6 is a better value.
 * For elasticity problems, a large strength threshold,
 * such as 0.9, is often better.
 */
JXF_Int JXF_PAMGSetStrongThreshold(JXF_Solver solver, JXF_Real strong_threshold);

JXF_Int JXF_PAMGGetStrongThreshold(JXF_Solver solver, JXF_Real *strong_threshold);

/*!
 * Set a parameter to modify the definition of strength for diagonal dominant
 * portions of the matrix. The default is 0.9. If max row sum is 1,
 * no checking for diagonally dominant rows is performed.
 */
JXF_Int JXF_PAMGSetMaxRowSum(JXF_Solver solver, JXF_Real max_row_sum);
JXF_Int JXF_PAMGGetMaxRowSum(JXF_Solver solver, JXF_Real *max_row_sum);

/*!
 * Define a truncation factor for the interpolation.
 * The default is 0.
 */
JXF_Int JXF_PAMGSetTruncFactor(JXF_Solver solver, JXF_Real trunc_factor);
JXF_Int JXF_PAMGGetTruncFactor(JXF_Solver solver, JXF_Real *trunc_factor);

/*!
 * Define the maximal number of elements per row for the interpolation.
 * The default is 0.
 */
JXF_Int JXF_PAMGSetPMaxElmts(JXF_Solver solver, JXF_Int P_max_elmts);
JXF_Int JXF_PAMGGetPMaxElmts(JXF_Solver solver, JXF_Int *P_max_elmts);

JXF_Int JXF_PAMGSetJacobiTruncThreshold(JXF_Solver solver, JXF_Real jacobi_trunc_threshold);
JXF_Int JXF_PAMGGetJacobiTruncThreshold(JXF_Solver solver, JXF_Real *jacobi_trunc_threshold);
JXF_Int JXF_PAMGSetPostInterpType(JXF_Solver solver, JXF_Int post_interp_type);
JXF_Int JXF_PAMGGetPostInterpType(JXF_Solver solver, JXF_Int *post_interp_type);

/*!
 * Define the largest strength threshold for which the strength
 * matrix S uses the communication package of the operator A.
 * If the strength threshold is larger than this values, a communication
 * package is generated for S. This can save memory and decrease the amount of data
 * that needs to be communicated, if S is substantially sparser than A.
 * The default is 1.0.
 */
JXF_Int JXF_PAMGSetSCommPkgSwitch(JXF_Solver solver, JXF_Real S_commpkg_switch);

/*!
 * Define which parallel interpolation operator is used.
 * There are the following options for interp type:
 * 0    classical modified interpolation
 * 3    direct interpolation (with separation of weights)
 * 4    multipass interpolation
 * 5    multipass interpolation (with separation of weights)
 * 6    extended classical modified interpolation
 * 7    extended (if no common C neighbor) classical modified interpolation
 * 8    standard interpolation
 * 9    standard interpolation (with separation of weights)
 * The default is 0.
 */
JXF_Int JXF_PAMGSetInterpType(JXF_Solver solver, JXF_Int interp_type);

/*!
 * Set minimum number of iterations.
 * The default is 0.
 */
JXF_Int JXF_PAMGSetMinIter(JXF_Solver solver, JXF_Int min_iter);

/*!
 * Set maximum number of iterations, if PAMG is used as a solver. If it is used as a
 * preconditioner, this function has no effect. The default is 20.
 */
JXF_Int JXF_PAMGSetMaxIter(JXF_Solver solver, JXF_Int max_iter);
JXF_Int JXF_PAMGGetMaxIter(JXF_Solver solver, JXF_Int *max_iter);

/*!
 * Define which parallel coarsening algorithm is used. There are the following options for
 * coarsen type:
 * 0  CLJP-coarsening (a parallel coarsening algorithm using independent sets.
 * 1  classical Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
 * 11 one-pass Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
 * 2  Ruge2B
 * 3  classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse
 *    points on the boundaries
 * 4  Ruge3c third pass on boundary,keep C points
 * 5  Ruge relax special points
 * 6  Falgout coarsening (uses 1 first, followed by CLJP using the interior coarse points
 *    generated by 1 as its first independent set)
 * 8  PMIS-coarsening (a parallel coarsening algorithm using independent sets, generating
 *    lower complexities than CLJP, might also lead to slower convergence)
 * 10 HMIS-coarsening (uses one pass Ruge-Stueben on each processor independently, followed
 *    by PMIS using the interior C-points generated as its first independent set)
 * 91 RRS0 relaxed RS0
 * 90 RCLJP relaxed CLJP
 *  The default is 6.
 */
JXF_Int JXF_PAMGSetCoarsenType(JXF_Solver solver, JXF_Int coarsen_type);
JXF_Int JXF_PAMGGetCoarsenType(JXF_Solver solver, JXF_Int *coarsen_type);

/*!
 * Set the measure type
 * 0: measures are determined locally
 * 1: measures are determined globally
 * The default is 0.
 */
JXF_Int JXF_PAMGSetMeasureType(JXF_Solver solver, JXF_Int measure_type);
JXF_Int JXF_PAMGGetMeasureType(JXF_Solver solver, JXF_Int *measure_type);

/*!
 * Set the setup type
 * 0: don't do the AMG setup
 * 1: do the AMG setup
 * The default is 1.
 * peghoty  2009/12/11
 */
JXF_Int JXF_PAMGSetSetupType(JXF_Solver solver, JXF_Int setup_type);

/*!
 * Define the type of cycle. For a V-cycle, set cycle type to 1,
 * for a W-cycle set cycle type to 2.
 * The default is 1.
 */
JXF_Int JXF_PAMGSetCycleType(JXF_Solver solver, JXF_Int cycle_type);
JXF_Int JXF_PAMGGetCycleType(JXF_Solver solver, JXF_Int *cycle_type);

/*!
 * Set the convergence tolerance, if PAMG is used as a solver.
 * If it is used as a preconditioner,
 * this function has no effect.
 * The default is 1.e-7.
 */
JXF_Int JXF_PAMGSetTol(JXF_Solver solver, JXF_Real tol);
JXF_Int JXF_PAMGGetTol(JXF_Solver solver, JXF_Real *tol);

/* peghoty, 2011/09/28 */
JXF_Int JXF_PAMGSetRhsNrmThreshold(JXF_Solver solver, JXF_Real rhsnrm_threshold);

/*!
 * Define the number of sweeps for the fine and coarse grid, the up and down cycle.
 * Note:  This routine will be phased out!!!!
 * Use JXF_PAMGSetNumSweeps or JXF_PAMGSetCycleNumSweeps instead.
 */
JXF_Int JXF_PAMGSetNumGridSweeps(JXF_Solver solver, JXF_Int *num_grid_sweeps);

/*!
 * Set the number of sweeps. On the finest level, the up and
 * the down cycle the number of sweeps are set to
 * num sweeps and on the coarsest level to 1. The default is 1.
 */
JXF_Int JXF_PAMGSetNumSweeps(JXF_Solver solver, JXF_Int num_sweeps);

/*!
 * Set the number of sweeps at a specified cycle.
 * There are the following options for k:
 * the finest level    if k=0
 * the down cycle      if k=1
 * the up cycle        if k=2
 * the coarsest level  if k=3.
 */
JXF_Int JXF_PAMGSetCycleNumSweeps(JXF_Solver solver, JXF_Int num_sweeps, JXF_Int k);
JXF_Int JXF_PAMGGetCycleNumSweeps(JXF_Solver solver, JXF_Int *num_sweeps, JXF_Int k);

/*!
 * Define which smoother is used on the fine and coarse grid, the up and down cycle.
 * Note: This routine will be phased out!!!!
 * Use JXF_PAMGSetRelaxType or JXF_PAMGSetCycleRelaxType instead.
 */
JXF_Int JXF_PAMGSetGridRelaxType(JXF_Solver solver, JXF_Int *grid_relax_type);

/*!
 * Define the smoother to be used. It uses the given smoother on the fine grid,
 * the up and the down cycle and sets the solver on the coarsest level to Gaussian
 * elimination (9). The default is Gauss-  Seidel(3).
 * There are the following options for relax type:
 * 0   Jacobi
 * 1   Gauss-Seidel, sequential (very slow!)
 * 2   Gauss-Seidel, interior points in parallel, boundary sequential (slow!)
 * 3   hybrid Gauss-Seidel or SOR, forward solve
 * 4   hybrid Gauss-Seidel or SOR, backward solve
 * 5   hybrid chaotic Gauss-Seidel (works only with OpenMP)
 * 6   hybrid symmetric Gauss-Seidel or SSOR
 * 7   Jacobi (uses Matvec), only needed in CGNR
 * 9   Gaussian elimination (only on coarsest level)
 * The default is 3.
 */
JXF_Int JXF_PAMGSetRelaxType(JXF_Solver solver, JXF_Int relax_type);

/*!
 * Define the smoother at a given cycle. For options of relax type see description of
 * JXF_PAMGSetRelaxType. Options for k are
 * the finest level    if k=0
 * the down cycle      if k=1
 * the up cycle        if k=2
 * the coarsest level  if k=3.
 */
JXF_Int JXF_PAMGSetCycleRelaxType(JXF_Solver solver, JXF_Int relax_type, JXF_Int k);
JXF_Int JXF_PAMGGetCycleRelaxType(JXF_Solver solver, JXF_Int *relax_type, JXF_Int k);

/*!
 * Define in which order the points are relaxed.
 * There are the following options for relax order:
 * 0  the points are relaxed in natural or lexicographic order on each processor
 * 1  CF-relaxation is used, i.e on the fine grid and the down cycle the
 *    coarse points are relaxed first, followed by the fine points;
 *    on the up cycle the F-points are relaxed first, followed by the C-points.
 *    On the coarsest level, if an iterative scheme is used,
 *    the points are relaxed in lexicographic order.
 * The default is 1 (CF-relaxation).
 */
JXF_Int JXF_PAMGSetRelaxOrder(JXF_Solver solver, JXF_Int relax_order);

/*!
 * Define in which order the points are relaxed.
 * Note: This routine will be phased out!!!!
 * Use JXF_PAMGSetRelaxOrder instead.
 */
JXF_Int JXF_PAMGSetGridRelaxPoints(JXF_Solver solver, JXF_Int **grid_relax_points);

/*!
 * Define the relaxation weight for smoothed Jacobi and hybrid SOR.
 * Note:  This routine will be phased out!!!!
 * Use JXF_PAMGSetRelaxWt or JXF_PAMGSetLevelRelaxWt instead.
 */
JXF_Int JXF_PAMGSetRelaxWeight(JXF_Solver solver, JXF_Real *relax_weight);

/*!
 * Define the relaxation weight for smoothed Jacobi and hybrid SOR on all levels.
 * relax weight > 0  this assigns the given relaxation weight on all levels
 * relax weight = 0  the weight is determined on each level with the estimate
 *                   3 / (4*||D^{-1/2}AD^{-1/2}||), where D is the diagonal matrix
 *                   of A (this should only be used with Jacobi)
 * relax weight = -k the relaxation weight is determined with at most k CG steps on each level
 *                   this should only be used for symmetric positive definite problems
 * The default is 1.
 */
JXF_Int JXF_PAMGSetRelaxWt(JXF_Solver solver, JXF_Real relax_wt);

/*!
 * Define the relaxation weight for smoothed Jacobi and hybrid SOR on the user defined level.
 * Note that the finest level is denoted 0, the next coarser level 1, etc.
 ＊ For nonpositive relax weight, the parameter is determined on
 * the given level as described for JXF_PAMGSetRelaxWt.
 ＊ The default is 1.
 */
JXF_Int JXF_PAMGSetLevelRelaxWt(JXF_Solver solver, JXF_Real relax_wt, JXF_Int level);

/*!
 * Define the outer relaxation weight for hybrid SOR.
 * Note: This routine will be phased out!!!!
 * Use JXF_PAMGSetOuterWt or JXF_PAMGSetLevelOuterWt instead.
 */
JXF_Int JXF_PAMGSetOmega(JXF_Solver solver, JXF_Real *omega);

/*!
 * Define the outer relaxation weight for hybrid SOR and SSOR on all levels.
 * omega > 0   this assigns the same outer relaxation weight omega on each level
 * omega = -k  an outer relaxation weight is determined with at most k CG steps on each level
               (this only makes sense for symmetric positive definite problems and smoothers, e.g. SSOR)
 * The default is 1.
 */
JXF_Int JXF_PAMGSetOuterWt(JXF_Solver solver, JXF_Real outer_wt);

/*!
 * Define the outer relaxation weight for hybrid SOR or SSOR on the user defined level.
 * Note that the finest level is denoted 0, the next coarser level 1, etc.
 * For nonpositive omega, the parameter is determined
 * on the given level as described for JXF_PAMGSetOuterWt.
 * The default is 1.
 */
JXF_Int JXF_PAMGSetLevelOuterWt(JXF_Solver solver, JXF_Real outer_wt, JXF_Int level);

/*!
 * Enables the use of more complex smoothers.
 * The following options exist for smooth type:
 *
 * value   smoother            routines needed to set smoother parameters
 *
 *  6     Schwarz smoothers     JXF_PAMGSetDomainType, JXF_PAMGSetOverlap,
 *                              JXF_PAMGSetVariant,    JXF_PAMGAMGSetSchwarzRlxWeight
 *  7     Pilut                 JXF_PAMGAMGSetDropTol, JXF_PAMGAMGSetMaxNzPerRow
 *  8     ParaSails             JXF_PAMGAMGSetSym,     JXF_PAMGAMGSetLevel,
 *                              JXF_PAMGAMGSetFilter,  JXF_PAMGAMGSetThreshold
 *  9     Euclid                JXF_PAMGAMGSetEuclidFile
 *
 * The default is 6.
 * Also, if no smoother parameters are set via the routines mentioned in the table above,
 * default values are used.
 */
JXF_Int JXF_PAMGSetSmoothType(JXF_Solver solver, JXF_Int smooth_type);
JXF_Int JXF_PAMGGetSmoothType(JXF_Solver solver, JXF_Int *smooth_type);

/*!
 * Set the number of levels for more complex smoothers.
 * The smoothers, as defined
 * by JXF_PAMGSetSmoothType, will be used on level 0 (the finest level) through level
 * smooth num levels-1.
 * The default is 0, i.e. no complex smoothers are used.
 */
JXF_Int JXF_PAMGSetSmoothNumLevels(JXF_Solver solver, JXF_Int smooth_num_levels);
JXF_Int JXF_PAMGGetSmoothNumLevels(JXF_Solver solver, JXF_Int *smooth_num_levels);

/*!
 * Set the number of sweeps for more complex smoothers.
 * The default is 1.
 */
JXF_Int JXF_PAMGSetSmoothNumSweeps(JXF_Solver solver, JXF_Int smooth_num_sweeps);
JXF_Int JXF_PAMGGetSmoothNumSweeps(JXF_Solver solver, JXF_Int *smooth_num_sweeps);

/*!
 * Requests additional computations for diagnostic and similar data to be logged by the user.
 * Default to 0 for do nothing. The latest residual will be available if logging > 1.
 */
JXF_Int JXF_PAMGSetLogging(JXF_Solver solver, JXF_Int logging);
JXF_Int JXF_PAMGGetLogging(JXF_Solver solver, JXF_Int *logging);

/*!
 * Requests automatic printing of setup and solve information.
 * 0   no printout (default)
 * 1   print setup information
 * 2   print solve information
 * 3   print both setup and solve information
 * Note, that if one desires to print information and uses PAMG as a preconditioner,
 * suggested print level is 1 to avoid excessive output,
 * and use print level of solver for solve phase information.
 */
JXF_Int JXF_PAMGSetPrintLevel(JXF_Solver solver, JXF_Int print_level);
JXF_Int JXF_PAMGGetPrintLevel(JXF_Solver solver, JXF_Int *print_level);

JXF_Int JXF_PAMGSetPrintCoarseMatrix(JXF_Solver solver, JXF_Int print_coarse_matrix);

JXF_Int JXF_PAMGSetPrintFileName(JXF_Solver solver, const char *print_file_name);
JXF_Int JXF_PAMGSetDebugFlag(JXF_Solver solver, JXF_Int debug_flag);
JXF_Int JXF_PAMGGetDebugFlag(JXF_Solver solver, JXF_Int *debug_flag);
JXF_Int JXF_PAMGGetNumLevels(JXF_Solver solver, JXF_Int *num_levels);
JXF_Int JXF_PAMGGetNumIterations(JXF_Solver solver, JXF_Int *num_iterations);
JXF_Int JXF_PAMGGetCumNumIterations(JXF_Solver solver, JXF_Int *cum_num_iterations);
JXF_Int JXF_PAMGGetResidual(JXF_Solver solver, JXF_ParVector *residual);
JXF_Int JXF_PAMGGetFinalRelativeResidualNorm(JXF_Solver solver, JXF_Real *rel_resid_norm);
JXF_Int JXF_PAMGSetVariant(JXF_Solver solver, JXF_Int variant);
JXF_Int JXF_PAMGGetVariant(JXF_Solver solver, JXF_Int *variant);
JXF_Int JXF_PAMGSetOverlap(JXF_Solver solver, JXF_Int overlap);
JXF_Int JXF_PAMGGetOverlap(JXF_Solver solver, JXF_Int *overlap);
JXF_Int JXF_PAMGSetDomainType(JXF_Solver solver, JXF_Int domain_type);
JXF_Int JXF_PAMGGetDomainType(JXF_Solver solver, JXF_Int *domain_type);
JXF_Int JXF_PAMGSetSchwarzRlxWeight(JXF_Solver solver, JXF_Real schwarz_rlx_weight);
JXF_Int JXF_PAMGGetSchwarzRlxWeight(JXF_Solver solver, JXF_Real *schwarz_rlx_weight);
JXF_Int JXF_PAMGSetSchwarzUseNonSymm(JXF_Solver solver, JXF_Int use_nonsymm);
JXF_Int JXF_PAMGSetSym(JXF_Solver solver, JXF_Int sym);
JXF_Int JXF_PAMGSetLevel(JXF_Solver solver, JXF_Int level);
JXF_Int JXF_PAMGSetThreshold(JXF_Solver solver, JXF_Real threshold);
JXF_Int JXF_PAMGSetFilter(JXF_Solver solver, JXF_Real filter);
JXF_Int JXF_PAMGSetDropTol(JXF_Solver solver, JXF_Real drop_tol);
JXF_Int JXF_PAMGSetMaxNzPerRow(JXF_Solver solver, JXF_Int max_nz_per_row);
JXF_Int JXF_PAMGSetEuclidFile(JXF_Solver solver, char *euclidfile);
JXF_Int JXF_PAMGSetEuLevel(JXF_Solver solver, JXF_Int eu_level);
JXF_Int JXF_PAMGSetEuSparseA(JXF_Solver solver, JXF_Real eu_sparse_A);
JXF_Int JXF_PAMGSetEuBJ(JXF_Solver solver, JXF_Int eu_bj);
JXF_Int JXF_PAMGSetNumFunctions(JXF_Solver solver, JXF_Int num_functions);
JXF_Int JXF_PAMGGetNumFunctions(JXF_Solver solver, JXF_Int *num_functions);
JXF_Int JXF_PAMGSetNodal(JXF_Solver solver, JXF_Int nodal);
JXF_Int JXF_PAMGSetNodalDiag(JXF_Solver solver, JXF_Int nodal);
JXF_Int JXF_PAMGSetDofFunc(JXF_Solver solver, JXF_Int *dof_func);
JXF_Int JXF_PAMGSetNumPaths(JXF_Solver solver, JXF_Int num_paths);
JXF_Int JXF_PAMGSetAggNumLevels(JXF_Solver solver, JXF_Int agg_num_levels);
JXF_Int JXF_PAMGSetAggInterpType(JXF_Solver solver, JXF_Int agg_interp_type);
JXF_Int JXF_PAMGSetAggPMaxElmts(JXF_Solver solver, JXF_Int agg_P_max_elmts);
JXF_Int JXF_PAMGSetAggP12MaxElmts(JXF_Solver solver, JXF_Int agg_P12_max_elmts);
JXF_Int JXF_PAMGSetAggTruncFactor(JXF_Solver solver, JXF_Real agg_trunc_factor);
JXF_Int JXF_PAMGSetAggP12TruncFactor(JXF_Solver solver, JXF_Real agg_P12_trunc_factor);
JXF_Int JXF_PAMGSetNumCRRelaxSteps(JXF_Solver solver, JXF_Int num_CR_relax_steps);
JXF_Int JXF_PAMGSetCRRate(JXF_Solver solver, JXF_Real CR_rate);
JXF_Int JXF_PAMGSetCRStrongTh(JXF_Solver solver, JXF_Real CR_strong_th);
JXF_Int JXF_PAMGSetISType(JXF_Solver solver, JXF_Int IS_type);
JXF_Int JXF_PAMGSetCRUseCG(JXF_Solver solver, JXF_Int CR_use_CG);
JXF_Int JXF_PAMGSetGSMG(JXF_Solver solver, JXF_Int gsmg);
JXF_Int JXF_PAMGSetNumSamples(JXF_Solver solver, JXF_Int gsmg);
JXF_Int JXF_PAMGSetCGCIts(JXF_Solver solver, JXF_Int its);
// Yue Xiaoqiang 2012/10/22
JXF_Int JXF_PAMGSetSpMtRapType(JXF_Solver solver, JXF_Int spmt_rap_type);
// Yue Xiaoqiang 2015/09/30
JXF_Int JXF_PAMGSetWallTimeOption(JXF_Solver solver, JXF_Int wall_time_option);
// Yue Xiaoqiang 2014/02/26
JXF_Int JXF_PAMGSetAIMeasureType(JXF_Solver solver, JXF_Int ai_measure_type);
// Yue Xiaoqiang 2014/07/06
JXF_Int JXF_PAMGSetAIRelaxType(JXF_Solver solver, JXF_Int ai_relax_type);
// peghoty 2010/05/29
JXF_Int JXF_PAMGSetRelaxedCoarsenMeasureType(JXF_Solver solver, JXF_Int measure_type_rlx);
// peghoty 2010/05/29
JXF_Int JXF_PAMGSetRelaxedCoarsenNumberSyn(JXF_Solver solver, JXF_Int number_syn_rlx);
// peghoty 2010/05/29
JXF_Int JXF_PAMGSetRelaxedCoarsenMeasureThreshold(JXF_Solver solver, JXF_Real measure_threshold_rlx);
JXF_Int JXF_PAMGSetPlotGrids(JXF_Solver solver, JXF_Int plotgrids);
JXF_Int JXF_PAMGSetPlotFileName(JXF_Solver solver, const char *plotfilename);
JXF_Int JXF_PAMGSetCoordDim(JXF_Solver solver, JXF_Int coorddim);
JXF_Int JXF_PAMGSetCoordinates(JXF_Solver solver, float *coordinates);

/*!
 * Set the stopping criteria for convergence
 * There are the following options for conv_criteria:
 * ( u^{k} is the k-th approxiated solution, r^{k}=b-A*u^{k} )
 *     0:  ||r^{k}||_2 / ||r^{0}||_2                 (Relative Residual Norm-2 (default))
 *     1:  ||u^{k} - u^{k-1}||_1 / ||u^{k}||_1       (Relative Error Norm-1)
 *    11:  max { |u_i^{k} - u_i^{k-1}| / |u_i^{k}| } (Relative Point-wise Error Norm-1)
 *     2:  ||u^{k} - u^{k-1}||_2 / ||u^{k}||_2       (Relative Error Norm-2)
 *   where
 *    ||u||_1 = max{|u_i|}， ||u||_2 = sqrt((u_1)^2 + ... + (u_n)^2)
 * The default is 0.
 * peghoty 2009/12/11
 */
JXF_Int JXF_PAMGSetConvCriteria(JXF_Solver solver, JXF_Int conv_criteria);

/*!
 * Set the solver on the coarsest level(including the case when max_levels=1)
 * There are the following options for coarsestsolverid:
 *   0 -> Jacobi or CF-Jacobi
 *   1 -> Gauss-Seidel <--- very slow, sequential
 *   2 -> Gauss_Seidel: interior points in parallel ,
 *        boundary sequential
 *   3 -> hybrid: SOR-J mix off-processor, SOR on-processor
 *                with outer relaxation parameters (forward solve)
 *   4 -> hybrid: SOR-J mix off-processor, SOR on-processor
 *                with outer relaxation parameters (backward solve)
 *   5 -> hybrid: GS-J mix off-processor, chaotic GS on-node
 *   6 -> hybrid: SSOR-J mix off-processor, SSOR on-processor
 *                with outer relaxation parameters
 *   9 -> Gaussian Elimination (default)
 *  10 -> Gaussian Elimination - with pivoting
 *  The default is 9.
 *  peghoty 2010/02/27
 */
JXF_Int JXF_PAMGSetCoarsestSolverID(JXF_Solver solver, JXF_Int coarsestsolverid);

/*!
 * Set the allowed number of gridpoints in the coarsest level.
 * The default is 100.
 * peghoty  2010/04/14
 */
JXF_Int JXF_PAMGSetCoarseThreshold(JXF_Solver solver, JXF_Int coarse_threshold);

/*!
 * Set the parameter which indicates when to modify the
 * current coarsening strategy.
 *   If coarse_size >= fine_size*coarse_ratio, we use CLJP
 * instead of the current coarsening strategy.
 * The default is 0.75.
 * peghoty  2010/04/14
 */
JXF_Int JXF_PAMGSetCoarseRatio(JXF_Solver solver, JXF_Real coarse_ratio);

/*
 *  * Set the parameter convfac_threshold.
 *   *   If the current convegence factor > convfac_threshold, the iteration
 *    * should be terminated at once.
 *     *
 *      * The default is 1.0.
 *       *
 *        * peghoty  2010/06/22
 *         */
JXF_Int JXF_PAMGSetConvFacThreshold(JXF_Solver solver, JXF_Real convfac_threshold);
JXF_Int JXF_PAMGSetRAP2(JXF_Solver solver, JXF_Int rap2);
JXF_Int JXF_PAMGSetKeepTranspose(JXF_Solver solver, JXF_Int keepTranspose);

void *jxf_PAMGCreate();
JXF_Int jxf_PAMGDestroy(void *data);
JXF_Int jxf_hpPAMGDestroy(void *data);
JXF_Int jxf_PAMGSetRestriction(void *data, JXF_Int restr_par);
JXF_Int jxf_PAMGSetAIRStrongTh(void *data, JXF_Real AIR_strong_th);
JXF_Int jxf_PAMGSetMaxLevels(void *data, JXF_Int max_levels);
JXF_Int jxf_PAMGGetMaxLevels(void *data, JXF_Int *max_levels);
JXF_Int jxf_PAMGSetStrongThreshold(void *data, JXF_Real strong_threshold);
JXF_Int jxf_PAMGGetStrongThreshold(void *data, JXF_Real *strong_threshold);
JXF_Int jxf_PAMGSetMaxRowSum(void *data, JXF_Real max_row_sum);
JXF_Int jxf_PAMGGetMaxRowSum(void *data, JXF_Real *max_row_sum);
JXF_Int jxf_PAMGSetTruncFactor(void *data, JXF_Real trunc_factor);
JXF_Int jxf_PAMGGetTruncFactor(void *data, JXF_Real *trunc_factor);
JXF_Int jxf_PAMGSetPMaxElmts(void *data, JXF_Int P_max_elmts);
JXF_Int jxf_PAMGGetPMaxElmts(void *data, JXF_Int *P_max_elmts);
JXF_Int jxf_PAMGSetJacobiTruncThreshold(void *data, JXF_Real jacobi_trunc_threshold);
JXF_Int jxf_PAMGGetJacobiTruncThreshold(void *data, JXF_Real *jacobi_trunc_threshold);
JXF_Int jxf_PAMGSetPostInterpType(void *data, JXF_Int post_interp_type);
JXF_Int jxf_PAMGGetPostInterpType(void *data, JXF_Int *post_interp_type);
JXF_Int jxf_PAMGSetSCommPkgSwitch(void *data, JXF_Real S_commpkg_switch);
JXF_Int jxf_PAMGGetSCommPkgSwitch(void *data, JXF_Real *S_commpkg_switch);
JXF_Int jxf_PAMGSetInterpType(void *data, JXF_Int interp_type);
JXF_Int jxf_PAMGGetInterpType(void *data, JXF_Int *interp_type);
JXF_Int jxf_PAMGSetMinIter(void *data, JXF_Int min_iter);
JXF_Int jxf_PAMGGetMinIter(void *data, JXF_Int *min_iter);
JXF_Int jxf_PAMGSetMaxIter(void *data, JXF_Int max_iter);
JXF_Int jxf_PAMGGetMaxIter(void *data, JXF_Int *max_iter);
JXF_Int jxf_PAMGSetCoarsenType(void *data, JXF_Int coarsen_type);
JXF_Int jxf_PAMGGetCoarsenType(void *data, JXF_Int *coarsen_type);
JXF_Int jxf_PAMGSetMeasureType(void *data, JXF_Int measure_type);
JXF_Int jxf_PAMGGetMeasureType(void *data, JXF_Int *measure_type);
JXF_Int jxf_PAMGSetSetupType(void *data, JXF_Int setup_type);
JXF_Int jxf_PAMGGetSetupType(void *data, JXF_Int *setup_type);
JXF_Int jxf_PAMGSetCycleType(void *data, JXF_Int cycle_type);
JXF_Int jxf_PAMGGetCycleType(void *data, JXF_Int *cycle_type);
JXF_Int jxf_PAMGSetTol(void *data, JXF_Real tol);
JXF_Int jxf_PAMGSetRhsNrmThreshold(void *data, JXF_Real rhsnrm_threshold);
JXF_Int jxf_PAMGGetTol(void *data, JXF_Real *tol);
JXF_Int jxf_PAMGSetNumSweeps(void *data, JXF_Int num_sweeps);
JXF_Int jxf_PAMGSetCycleNumSweeps(void *data, JXF_Int num_sweeps, JXF_Int k);
JXF_Int jxf_PAMGGetCycleNumSweeps(void *data, JXF_Int *num_sweeps, JXF_Int k);
JXF_Int jxf_PAMGSetNumGridSweeps(void *data, JXF_Int *num_grid_sweeps);
JXF_Int jxf_PAMGGetNumGridSweeps(void *data, JXF_Int **num_grid_sweeps);
JXF_Int jxf_PAMGSetRelaxType(void *data, JXF_Int relax_type);
JXF_Int jxf_PAMGSetCycleRelaxType(void *data, JXF_Int relax_type, JXF_Int k);
JXF_Int jxf_PAMGGetCycleRelaxType(void *data, JXF_Int *relax_type, JXF_Int k);
JXF_Int jxf_PAMGSetRelaxOrder(void *data, JXF_Int relax_order);
JXF_Int jxf_PAMGGetRelaxOrder(void *data, JXF_Int *relax_order);
JXF_Int jxf_PAMGSetGridRelaxType(void *data, JXF_Int *grid_relax_type);
JXF_Int jxf_PAMGGetGridRelaxType(void *data, JXF_Int **grid_relax_type);
JXF_Int jxf_PAMGSetGridRelaxPoints(void *data, JXF_Int **grid_relax_points);
JXF_Int jxf_PAMGGetGridRelaxPoints(void *data, JXF_Int ***grid_relax_points);
JXF_Int jxf_PAMGSetRelaxWeight(void *data, JXF_Real *relax_weight);
JXF_Int jxf_PAMGGetRelaxWeight(void *data, JXF_Real **relax_weight);
JXF_Int jxf_PAMGSetRelaxWt(void *data, JXF_Real relax_weight);
JXF_Int jxf_PAMGSetLevelRelaxWt(void *data, JXF_Real relax_weight, JXF_Int level);
JXF_Int jxf_PAMGGetLevelRelaxWt(void *data, JXF_Real *relax_weight, JXF_Int level);
JXF_Int jxf_PAMGSetOmega(void *data, JXF_Real *omega);
JXF_Int jxf_PAMGGetOmega(void *data, JXF_Real **omega);
JXF_Int jxf_PAMGSetOuterWt(void *data, JXF_Real omega);
JXF_Int jxf_PAMGSetLevelOuterWt(void *data, JXF_Real omega, JXF_Int level);
JXF_Int jxf_PAMGGetLevelOuterWt(void *data, JXF_Real *omega, JXF_Int level);
JXF_Int jxf_PAMGSetSmoothType(void *data, JXF_Int smooth_type);
JXF_Int jxf_PAMGGetSmoothType(void *data, JXF_Int *smooth_type);
JXF_Int jxf_PAMGSetSmoothNumLevels(void *data, JXF_Int smooth_num_levels);
JXF_Int jxf_PAMGGetSmoothNumLevels(void *data, JXF_Int *smooth_num_levels);
JXF_Int jxf_PAMGSetSmoothNumSweeps(void *data, JXF_Int smooth_num_sweeps);
JXF_Int jxf_PAMGGetSmoothNumSweeps(void *data, JXF_Int *smooth_num_sweeps);
JXF_Int jxf_PAMGSetLogging(void *data, JXF_Int logging);
JXF_Int jxf_PAMGGetLogging(void *data, JXF_Int *logging);
JXF_Int jxf_PAMGSetPrintLevel(void *data, JXF_Int print_level);
JXF_Int jxf_PAMGGetPrintLevel(void *data, JXF_Int *print_level);
JXF_Int jxf_PAMGSetPrintCoarseMatrix(void *data, JXF_Int print_coarse_matrix);
JXF_Int jxf_PAMGSetPrintFileName(void *data, const char *print_file_name);
JXF_Int jxf_PAMGGetPrintFileName(void *data, char **print_file_name);
JXF_Int jxf_PAMGSetNumIterations(void *data, JXF_Int num_iterations);
JXF_Int jxf_PAMGSetDebugFlag(void *data, JXF_Int debug_flag);
JXF_Int jxf_PAMGGetDebugFlag(void *data, JXF_Int *debug_flag);
JXF_Int jxf_PAMGSetGSMG(void *data, JXF_Int par);
JXF_Int jxf_PAMGSetNumSamples(void *data, JXF_Int par);
JXF_Int jxf_PAMGSetCGCIts(void *data, JXF_Int its);
// Yue Xiaoqiang 2012/10/22
JXF_Int jxf_PAMGSetSpMtRapType(void *data, JXF_Int spmt_rap_type);
// Yue Xiaoqiang 2015/09/30
JXF_Int jxf_PAMGSetWallTimeOption(void *data, JXF_Int wall_time_option);
// Yue Xiaoqiang 2014/02/26
JXF_Int jxf_PAMGSetAIMeasureType(void *amg_data, JXF_Int ai_measure_type);
// Yue Xiaoqiang 2014/07/06
JXF_Int jxf_PAMGSetAIRelaxType(void *data, JXF_Int ai_relax_type);
// peghoty 2010/05/29
JXF_Int jxf_PAMGSetRelaxedCoarsenMeasureType(void *data, JXF_Int measure_type_rlx);
// peghoty 2010/05/29
JXF_Int jxf_PAMGSetRelaxedCoarsenNumberSyn(void *data, JXF_Int number_syn_rlx);
// peghoty 2010/05/29
JXF_Int jxf_PAMGSetRelaxedCoarsenMeasureThreshold(void *data, JXF_Real measure_threshold_rlx);
JXF_Int jxf_PAMGSetPlotGrids(void *data, JXF_Int plotgrids);
JXF_Int jxf_PAMGSetPlotFileName(void *data, const char *plot_file_name);
JXF_Int jxf_PAMGSetCoordDim(void *data, JXF_Int coorddim);
JXF_Int jxf_PAMGSetCoordinates(void *data, float *coordinates);
JXF_Int jxf_PAMGSetNumFunctions(void *data, JXF_Int num_functions);
JXF_Int jxf_PAMGGetNumFunctions(void *data, JXF_Int *num_functions);
JXF_Int jxf_PAMGSetNodal(void *data, JXF_Int nodal);
JXF_Int jxf_PAMGSetNodalDiag(void *data, JXF_Int nodal);
JXF_Int jxf_PAMGSetNumPaths(void *data, JXF_Int num_paths);
JXF_Int jxf_PAMGSetAggNumLevels(void *data, JXF_Int agg_num_levels);
JXF_Int jxf_PAMGSetAggInterpType(void *data, JXF_Int agg_interp_type);
JXF_Int jxf_PAMGSetAggPMaxElmts(void *data, JXF_Int agg_P_max_elmts);
JXF_Int jxf_PAMGSetAggP12MaxElmts(void *data, JXF_Int agg_P12_max_elmts);
JXF_Int jxf_PAMGSetAggTruncFactor(void *data, JXF_Real agg_trunc_factor);
JXF_Int jxf_PAMGSetAggP12TruncFactor(void *data, JXF_Real agg_P12_trunc_factor);
JXF_Int jxf_PAMGSetNumCRRelaxSteps(void *data, JXF_Int num_CR_relax_steps);
JXF_Int jxf_PAMGSetCRRate(void *data, JXF_Real CR_rate);
JXF_Int jxf_PAMGSetCRStrongTh(void *data, JXF_Real CR_strong_th);
JXF_Int jxf_PAMGSetISType(void *data, JXF_Int IS_type);
JXF_Int jxf_PAMGSetCRUseCG(void *data, JXF_Int CR_use_CG);
JXF_Int jxf_PAMGSetNumPoints(void *data, JXF_Int num_points);
JXF_Int jxf_PAMGSetDofFunc(void *data, JXF_Int *dof_func);
JXF_Int jxf_PAMGSetPointDofMap(void *data, JXF_Int *point_dof_map);
JXF_Int jxf_PAMGSetDofPoint(void *data, JXF_Int *dof_point);
JXF_Int jxf_PAMGGetNumLevels(void *data, JXF_Int *num_levels);
JXF_Int jxf_PAMGGetNumIterations(void *data, JXF_Int *num_iterations);
JXF_Int jxf_PAMGGetCumNumIterations(void *data, JXF_Int *cum_num_iterations);
JXF_Int jxf_PAMGGetResidual(void *data, jxf_ParVector **resid);
JXF_Int jxf_PAMGGetRelResidualNorm(void *data, JXF_Real *rel_resid_norm);
JXF_Int jxf_PAMGSetVariant(void *data, JXF_Int variant);
JXF_Int jxf_PAMGGetVariant(void *data, JXF_Int *variant);
JXF_Int jxf_PAMGSetOverlap(void *data, JXF_Int overlap);
JXF_Int jxf_PAMGGetOverlap(void *data, JXF_Int *overlap);
JXF_Int jxf_PAMGSetDomainType(void *data, JXF_Int domain_type);
JXF_Int jxf_PAMGGetDomainType(void *data, JXF_Int *domain_type);
JXF_Int jxf_PAMGSetSchwarzRlxWeight(void *data, JXF_Real schwarz_rlx_weight);
JXF_Int jxf_PAMGGetSchwarzRlxWeight(void *data, JXF_Real *schwarz_rlx_weight);
JXF_Int jxf_PAMGSetSchwarzUseNonSymm(void *data, JXF_Int use_nonsymm);
JXF_Int jxf_PAMGSetSym(void *data, JXF_Int sym);
JXF_Int jxf_PAMGSetLevel(void *data, JXF_Int level);
JXF_Int jxf_PAMGSetThreshold(void *data, JXF_Real thresh);
JXF_Int jxf_PAMGSetFilter(void *data, JXF_Real filter);
JXF_Int jxf_PAMGSetDropTol(void *data, JXF_Real drop_tol);
JXF_Int jxf_PAMGSetMaxNzPerRow(void *data, JXF_Int max_nz_per_row);
JXF_Int jxf_PAMGSetEuclidFile(void *data, char *euclidfile);
JXF_Int jxf_PAMGSetEuLevel(void *data, JXF_Int eu_level);
JXF_Int jxf_PAMGSetEuSparseA(void *data, JXF_Real eu_sparse_A);
JXF_Int jxf_PAMGSetEuBJ(void *data, JXF_Int eu_bj);
/* peghoty 2009/07/27 */
JXF_Int jxf_PAMGSetConvCriteria(void *data, JXF_Int conv_criteria);
JXF_Int jxf_PAMGSetCoarsestSolverID(void *data, JXF_Int coarsestsolverid);
/* peghoty 2010/04/14 */
JXF_Int jxf_PAMGSetCoarseThreshold(void *data, JXF_Int coarse_threshold);
/* peghoty 2010/04/14 */
JXF_Int jxf_PAMGSetCoarseRatio(void *data, JXF_Real coarse_ratio);

JXF_Int jxf_PAMGSetConvFacThreshold(void *data, JXF_Real convfac_threshold);

JXF_Int jxf_PAMGSetRAP2(void *data, JXF_Int rap2);
JXF_Int jxf_PAMGSetKeepTranspose(void *data, JXF_Int keepTranspose);

/* csrc/amg/par_cparam.c */
JXF_Int
jxf_PAMGCoarseParms(MPI_Comm comm,
                   JXF_Int local_num_variables,
                   JXF_Int num_functions,
                   JXF_Int *dof_func,
                   JXF_Int *CF_marker,
                   JXF_Int **coarse_dof_func_ptr,
                   JXF_Int **coarse_pnts_global_ptr);

/* csrc/amg/par_amg_setup.c */
JXF_Int
jxf_PAMGSetup(void *amg_vdata, jxf_ParCSRMatrix *par_matrix);
JXF_Int
jxf_hpPAMGSetup(void *amg_vdata,
               jxf_hpCSRMatrix *hp_matrix);

/* csrc/amg/par_amg_solve.c */
JXF_Int jxf_PAMGSolve(void *amg_vdata, jxf_ParCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app);
JXF_Int jxf_hpPAMGSolve(void *amg_vdata, jxf_hpCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app);
// peghoty 2011/09/04
// JXF_Int jxf_PAMGPrecond( void *amg_vdata, jxf_ParCSRMatrix *par_matrix, jxf_ParVector *par_rhs, jxf_ParVector *par_app );
JXF_Int
jxf_PAMGPrecond(void *amg_vdata,
               jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector *par_rhs,
               jxf_ParVector *par_app);

/* csrc/amg/coarsen_cljp.c */
JXF_Int
jxf_PAMGCoarsen(jxf_ParCSRMatrix *par_S,
               jxf_ParCSRMatrix *par_matrix,
               JXF_Int CF_init,
               JXF_Int debug_flag,
               JXF_Int **CF_marker_ptr);
JXF_Int
jxf_hpPAMGCoarsen(jxf_hpCSRMatrix *par_S,
                 jxf_hpCSRMatrix *par_A,
                 JXF_Int CF_init,
                 JXF_Int debug_flag,
                 JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_cljp_ai.c */
JXF_Int
jxf_PAMGCoarsenAI(jxf_ParCSRMatrix *par_S,
                 jxf_ParCSRMatrix *par_matrix,
                 JXF_Real *AI_measure,
                 JXF_Int CF_init,
                 JXF_Int debug_flag,
                 JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_falgout.c */
JXF_Int
jxf_PAMGCoarsenFalgout(jxf_ParCSRMatrix *par_S,
                      jxf_ParCSRMatrix *par_matrix,
                      JXF_Int measure_type,
                      JXF_Int debug_flag,
                      JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_falgout_ai.c */
JXF_Int
jxf_PAMGCoarsenFalgoutAI(jxf_ParCSRMatrix *par_S,
                        jxf_ParCSRMatrix *par_matrix,
                        JXF_Real *AI_measure,
                        JXF_Int measure_type,
                        JXF_Int debug_flag,
                        JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_hmis.c */
JXF_Int
jxf_PAMGCoarsenHMIS(jxf_ParCSRMatrix *par_S,
                   jxf_ParCSRMatrix *par_matrix,
                   JXF_Int measure_type,
                   JXF_Int debug_flag,
                   JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_hmis_ai.c */
JXF_Int
jxf_PAMGCoarsenHMISAI(jxf_ParCSRMatrix *par_S,
                     jxf_ParCSRMatrix *par_matrix,
                     JXF_Real *AI_measure,
                     JXF_Int CF_init,
                     JXF_Int measure_type,
                     JXF_Int debug_flag,
                     JXF_Int **CF_marker_ptr);
/* csrc/amg/coarsen_pmis.c */
JXF_Int
jxf_PAMGCoarsenPMIS(jxf_ParCSRMatrix *par_S,
                   jxf_ParCSRMatrix *par_matrix,
                   JXF_Int CF_init,
                   JXF_Int debug_flag,
                   JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_pmis_ai.c */
JXF_Int
jxf_PAMGCoarsenPMISAI(jxf_ParCSRMatrix *par_S,
                     jxf_ParCSRMatrix *par_matrix,
                     JXF_Real *AI_measure,
                     JXF_Int CF_init,
                     JXF_Int debug_flag,
                     JXF_Int **CF_marker_ptr);

/* csrc/amg/splitting_ai.c */
JXF_Int
jxf_PAMGMeasureAI(jxf_ParCSRMatrix *par_S,
                 jxf_ParCSRMatrix *par_matrix,
                 JXF_Int debug_flag,
                 JXF_Real **AI_measure_ptr);

JXF_Int
jxf_PAMGMeasureGAI_Base(jxf_ParCSRMatrix *par_S,
                       jxf_ParCSRMatrix *par_A,
                       JXF_Int level,
                       JXF_Int debug_flag,
                       JXF_Real gai_threshold,
                       JXF_Real **measure_gai_ptr,
                       JXF_Int **marker_gai_ptr,
                       JXF_Int **side_gai_ptr,
                       JXF_Real *tau_out);

/* csrc/amg/coarsen_rcljp.c */
JXF_Int
jxf_PAMGCoarsenRCLJP(jxf_ParCSRMatrix *S,
                    jxf_ParCSRMatrix *A,
                    JXF_Int measure_type_rlx,
                    JXF_Int number_syn_rlx,
                    JXF_Real measure_threshold_rlx,
                    JXF_Int CF_init,
                    JXF_Int debug_flag,
                    JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_rrs0.c */
JXF_Int
jxf_PAMGCoarsenRRS0(jxf_ParCSRMatrix *S,
                   jxf_ParCSRMatrix *A,
                   JXF_Int measure_type,
                   JXF_Int measure_type_rlx,
                   JXF_Int number_syn_rlx,
                   JXF_Real measure_threshold_rlx,
                   JXF_Int debug_flag,
                   JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_ruge.c */
JXF_Int
jxf_PAMGCoarsenRuge(jxf_ParCSRMatrix *par_S,
                   jxf_ParCSRMatrix *par_matrix,
                   JXF_Int measure_type,
                   JXF_Int coarsen_type,
                   JXF_Int debug_flag,
                   JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_ruge_ai.c */
JXF_Int
jxf_PAMGCoarsenRugeAI(jxf_ParCSRMatrix *par_S,
                     jxf_ParCSRMatrix *par_matrix,
                     JXF_Real *AI_measure,
                     JXF_Int CF_init,
                     JXF_Int measure_type,
                     JXF_Int coarsen_type,
                     JXF_Int debug_flag,
                     JXF_Int **CF_marker_ptr);

/* csrc/amg/coarsen_xml.c */
JXF_Int
jxf_PAMGCoarsenXML(jxf_ParCSRMatrix *par_S,
                  jxf_ParCSRMatrix *par_matrix,
                  JXF_Real *AI_measure,
                  JXF_Int measure_type,
                  JXF_Int coarsen_type,
                  JXF_Int debug_flag,
                  JXF_Int **CF_marker_ptr);

/* csrc/amg/par_csolver.c */
JXF_Int
jxf_CoarsestSolver(void *amg_vdata,
                  jxf_hpCSRMatrix *par_matrix,
                  jxf_ParVector *par_rhs,
                  jxf_ParVector *par_app);
JXF_Int
jxf_check_convergence(JXF_Int *relax_marker,
                     jxf_ParCSRMatrix *par_matrix,
                     jxf_ParVector *par_rhs,
                     JXF_Real alpha,
                     JXF_Real beta,
                     JXF_Int iter,
                     jxf_ParVector *par_app,
                     jxf_ParVector *Vtemp,
                     jxf_ParVector *Utemp,
                     JXF_Int my_id,
                     JXF_Int print_level,
                     JXF_Int conv_criteria,
                     JXF_Real rhs_norm,
                     JXF_Real resid_nrm_init,
                     JXF_Real *resid_nrm,
                     JXF_Real *relative,
                     JXF_Real *conv_factor);
JXF_Int
jxf_hp_check_convergence(JXF_Int *relax_marker,
                        jxf_hpCSRMatrix *par_matrix,
                        jxf_ParVector *par_rhs,
                        JXF_Real alpha,
                        JXF_Real beta,
                        JXF_Int iter,
                        jxf_ParVector *par_app,
                        jxf_ParVector *Vtemp,
                        jxf_ParVector *Utemp,
                        JXF_Int my_id,
                        JXF_Int print_level,
                        JXF_Int conv_criteria,
                        JXF_Real rhs_norm,
                        JXF_Real resid_nrm_init,
                        JXF_Real *resid_nrm_old,
                        JXF_Real *relative_ret,
                        JXF_Real *conv_factor_old);
JXF_Int jxf_CoarsestSolverInfo(MPI_Comm comm, void *amg_vdata);

/* csrc/amg/par_cycle.c */
JXF_Int jxf_PAMGCycle(void *amg_vdata, jxf_ParVector **F_array, jxf_ParVector **U_array);
JXF_Int jxf_hpPAMGCycle(void *amg_vdata, jxf_ParVector **F_array, jxf_ParVector **U_array);

/* csrc/amg/par_gesolver.c */
JXF_Int
jxf_GaussElimination(JXF_Int solverid,
                    jxf_ParCSRMatrix *par_matrix,
                    jxf_ParVector *par_rhs,
                    jxf_ParVector *par_app);

/* csrc/amg/indepset.c */
JXF_Int
jxf_PAMGIndepSetInit(jxf_ParCSRMatrix *par_S, JXF_Real *measure_array, JXF_Int seq_rand);

JXF_Int
jxf_PAMGIndepSet(jxf_ParCSRMatrix *par_S,
                JXF_Real *measure_array,
                JXF_Int *graph_array,
                JXF_Int graph_array_size,
                JXF_Int *graph_array_offd,
                JXF_Int graph_array_offd_size,
                JXF_Int *IS_marker,
                JXF_Int *IS_marker_offd);

/* csrc/amg/par_indepsetrlx.c */
JXF_Int jxf_PAMGIndepSetRelaxInit(jxf_ParCSRMatrix *S, JXF_Real *measure_array);
JXF_Int
jxf_PAMGIndepSetRelax(jxf_ParCSRMatrix *S,
                     jxf_CSRMatrix *S_ext,
                     JXF_Real *measure_array,
                     JXF_Int *graph_array,
                     JXF_Int graph_array_size,
                     JXF_Int *graph_array_offd,
                     JXF_Int graph_array_offd_size,
                     JXF_Int *IS_marker,
                     JXF_Int *IS_marker_offd,
                     JXF_Int nc_per,
                     JXF_Real measure_th);

/* csrc/amg/par_interp_0.c */
JXF_Int
jxf_PAMGBuildInterp(jxf_ParCSRMatrix *par_A,
                   JXF_Int *CF_marker,
                   jxf_ParCSRMatrix *par_S,
                   JXF_Int *num_cpts_global,
                   JXF_Int num_functions,
                   JXF_Int *dof_func,
                   JXF_Int debug_flag,
                   JXF_Real trunc_factor,
                   JXF_Int max_elmts,
                   JXF_Int *col_offd_S_to_A,
                   jxf_ParCSRMatrix **P_ptr);
JXF_Int jxf_PAMGInterpTruncation(jxf_ParCSRMatrix *par_P, JXF_Real trunc_factor, JXF_Int max_elmts);
JXF_Int jxf_GetCommPkgRTFromCommPkgA(jxf_ParCSRMatrix *RT, jxf_ParCSRMatrix *par_A, JXF_Int *fine_to_coarse_offd);
JXF_Int jxf_PAMGBuildInterp1(jxf_ParCSRMatrix *par_A,
                           JXF_Int *CF_marker,
                           jxf_ParCSRMatrix *par_S,
                           JXF_Int *num_cpts_global,
                           JXF_Int num_functions,
                           JXF_Int *dof_func,
                           JXF_Int debug_flag,
                           JXF_Real trunc_factor,
                           JXF_Int max_elmts,
                           JXF_Int *col_offd_S_to_A,
                           jxf_ParCSRMatrix **P_ptr,
                           JXF_Int *opt_icor);
/* csrc/amg/par_interp_100.c */
JXF_Int
jxf_PAMGBuildInterpOnePnt(jxf_ParCSRMatrix *A,
                         JXF_Int *CF_marker,
                         jxf_ParCSRMatrix *S,
                         JXF_Int *num_cpts_global,
                         JXF_Int num_functions,
                         JXF_Int *dof_func,
                         JXF_Int debug_flag,
                         JXF_Int *col_offd_S_to_A,
                         jxf_ParCSRMatrix **P_ptr);

/* csrc/amg/par_interp_200.c */
JXF_Int
jxf_PAMGBuildInterpDinvRATrans(jxf_ParCSRMatrix *A,
                              jxf_ParCSRMatrix *R,
                              jxf_ParCSRMatrix **P_ptr);

/* csrc/amg/par_interp_3.c */
JXF_Int
jxf_PAMGBuildDirInterp(jxf_ParCSRMatrix *par_A,
                      JXF_Int *CF_marker,
                      jxf_ParCSRMatrix *par_S,
                      JXF_Int *num_cpts_global,
                      JXF_Int num_functions,
                      JXF_Int *dof_func,
                      JXF_Int debug_flag,
                      JXF_Real trunc_factor,
                      JXF_Int max_elmts,
                      JXF_Int *col_offd_S_to_A,
                      jxf_ParCSRMatrix **P_ptr);
/* csrc/amg/par_interp_4.c */
JXF_Int
jxf_PAMGBuildMultipass(jxf_ParCSRMatrix *par_A,
                      JXF_Int *CF_marker,
                      jxf_ParCSRMatrix *par_S,
                      JXF_Int *num_cpts_global,
                      JXF_Int num_functions,
                      JXF_Int *dof_func,
                      JXF_Int debug_flag,
                      JXF_Real trunc_factor,
                      JXF_Int P_max_elmts,
                      JXF_Int weight_option,
                      JXF_Int *col_offd_S_to_A,
                      jxf_ParCSRMatrix **P_ptr);

/* csrc/amg/par_interp_6.c */
JXF_Int
jxf_PAMGBuildExtPIInterp(jxf_ParCSRMatrix *par_A,
                        JXF_Int *CF_marker,
                        jxf_ParCSRMatrix *par_S,
                        JXF_Int *num_cpts_global,
                        JXF_Int num_functions,
                        JXF_Int *dof_func,
                        JXF_Int debug_flag,
                        JXF_Real trunc_factor,
                        JXF_Int max_elmts,
                        JXF_Int *col_offd_S_to_A,
                        jxf_ParCSRMatrix **P_ptr);
JXF_Int
jxf_new_offd_nodes(JXF_Int **found,
                  JXF_Int num_cols_A_offd,
                  JXF_Int *A_ext_i,
                  JXF_Int *A_ext_j,
                  JXF_Int num_cols_S_offd,
                  JXF_Int *col_map_offd,
                  JXF_Int col_1,
                  JXF_Int col_n,
                  JXF_Int *Sop_i,
                  JXF_Int *Sop_j,
                  JXF_Int *CF_marker,
                  jxf_ParCSRCommPkg *comm_pkg);
void jxf_initialize_vecs(JXF_Int diag_n,
                        JXF_Int offd_n,
                        JXF_Int *diag_ftc,
                        JXF_Int *offd_ftc,
                        JXF_Int *diag_pm,
                        JXF_Int *offd_pm,
                        JXF_Int *tmp_CF);
JXF_Int jxf_ssort(JXF_Int *data, JXF_Int size);
JXF_Int jxf_index_of_minimum(JXF_Int *data, JXF_Int size);
void jxf_swap_int(JXF_Int *data, JXF_Int va, JXF_Int vb);
JXF_Int
jxf_alt_insert_new_nodes(jxf_ParCSRCommPkg *comm_pkg,
                        jxf_ParCSRCommPkg *extend_comm_pkg,
                        JXF_Int *IN_marker,
                        JXF_Int full_off_procNodes,
                        JXF_Int *OUT_marker);
JXF_Int
jxf_ParCSRFindExtendCommPkg(jxf_ParCSRMatrix *par_A,
                           JXF_Int newoff,
                           JXF_Int *found,
                           jxf_ParCSRCommPkg **extend_comm_pkg);

/* csrc/amg/par_interp_7.c */
JXF_Int
jxf_PAMGBuildExtPICCInterp(jxf_ParCSRMatrix *par_A,
                          JXF_Int *CF_marker,
                          jxf_ParCSRMatrix *par_S,
                          JXF_Int *num_cpts_global,
                          JXF_Int num_functions,
                          JXF_Int *dof_func,
                          JXF_Int debug_flag,
                          JXF_Real trunc_factor,
                          JXF_Int max_elmts,
                          JXF_Int *col_offd_S_to_A,
                          jxf_ParCSRMatrix **P_ptr);

/* csrc/amg/par_interp_8.c */
JXF_Int
jxf_PAMGBuildStdInterp(jxf_ParCSRMatrix *par_A,
                      JXF_Int *CF_marker,
                      jxf_ParCSRMatrix *par_S,
                      JXF_Int *num_cpts_global,
                      JXF_Int num_functions,
                      JXF_Int *dof_func,
                      JXF_Int debug_flag,
                      JXF_Real trunc_factor,
                      JXF_Int max_elmts,
                      JXF_Int sep_weight,
                      JXF_Int *col_offd_S_to_A,
                      jxf_ParCSRMatrix **P_ptr);

/* csrc/amg/par_itersolver.c */
JXF_Int
jxf_CoarsestIterativeMethod(jxf_hpCSRMatrix *par_matrix,
                           jxf_ParVector *par_rhs,
                           JXF_Int solverid,
                           JXF_Real weight,
                           JXF_Real omega,
                           jxf_ParVector *par_app,
                           jxf_ParVector *Vtemp);

/* csrc/amg/par_itersolver_ai.c */
JXF_Int
jxf_CoarsestIterativeMethodAI(jxf_hpCSRMatrix *par_matrix,
                             jxf_ParVector *par_rhs,
                             JXF_Int *relax_marker,
                             JXF_Int solverid,
                             JXF_Real weight,
                             JXF_Real omega,
                             jxf_ParVector *par_app,
                             jxf_ParVector *Vtemp);

/* csrc/amg/par_relax.c */
JXF_Int
jxf_PAMGRelax(jxf_ParCSRMatrix *par_matrix,
             jxf_ParVector *par_rhs,
             JXF_Int *cf_marker,
             JXF_Int relax_type,
             JXF_Int relax_points,
             JXF_Real relax_weight,
             JXF_Real omega,
             JXF_Real *l1_norms,
             jxf_ParVector *par_app,
             jxf_ParVector *Vtemp,
             jxf_ParVector *Ztemp);
JXF_Int
jxf_hpPAMGRelax(jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_type,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp);
/* csrc/amg/par_relax_if.c */
JXF_Int
jxf_PAMGRelaxIF(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int relax_type,
               JXF_Int relax_order,
               JXF_Int cycle_type,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp);
JXF_Int
jxf_hpPAMGRelaxIF(jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int relax_type,
                 JXF_Int relax_order,
                 JXF_Int cycle_type,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 JXF_Real *l1_norms,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp,
                 jxf_ParVector *Ztemp);
/* csrc/amg/par_relax_ai.c */
JXF_Int
jxf_PAMGRelaxAI(jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector *par_rhs,
               JXF_Int *cf_marker,
               JXF_Int *relax_marker,
               /*JXF_Real             *ai_measure,
               JXF_Real           ai_measure_th,
               JXF_Real           ai_measure_th_h, */
               JXF_Int relax_type,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               jxf_ParVector *par_app,
               jxf_ParVector *Vtemp);
JXF_Int
jxf_hpPAMGRelaxAI(jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 JXF_Int relax_type,
                 JXF_Int relax_points,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);
JXF_Int
jxf_PAMGRelaxIFAI(jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector *par_rhs,
                 JXF_Int *cf_marker,
                 JXF_Int *relax_marker,
                 /*JXF_Real             *ai_measure,
                 JXF_Real           ai_measure_th,
                 JXF_Real           ai_measure_th_h, */
                 JXF_Int relax_type,
                 JXF_Int relax_order,
                 JXF_Int cycle_type,
                 JXF_Real relax_weight,
                 JXF_Real omega,
                 jxf_ParVector *par_app,
                 jxf_ParVector *Vtemp);
JXF_Int
jxf_hpPAMGRelaxIFAI(jxf_hpCSRMatrix *hp_matrix,
                   jxf_ParVector *par_rhs,
                   JXF_Int *cf_marker,
                   JXF_Int *relax_marker,
                   JXF_Int relax_type,
                   JXF_Int relax_order,
                   JXF_Int cycle_type,
                   JXF_Real relax_weight,
                   JXF_Real omega,
                   jxf_ParVector *par_app,
                   jxf_ParVector *Vtemp);

/* csrc/amg/par_status.c */
JXF_Int jxf_PAMGSetupStatus(void *amg_vdata, jxf_ParCSRMatrix *par_matrix);
JXF_Int jxf_hpPAMGSetupStatus(void *amg_vdata, jxf_hpCSRMatrix *par_matrix);
JXF_Int jxf_PAMGSolveStatus(void *data);

/* csrc/amg/par_restr.c */
JXF_Int
jxf_PAMGBuildRestrAIR(jxf_ParCSRMatrix *A,
                     JXF_Int *CF_marker,
                     jxf_ParCSRMatrix *S,
                     JXF_Int *num_cpts_global,
                     JXF_Int num_functions,
                     JXF_Int *dof_func,
                     JXF_Int debug_flag,
                     JXF_Real trunc_factor,
                     JXF_Int max_elmts,
                     JXF_Int *col_offd_S_to_A,
                     jxf_ParCSRMatrix **R_ptr,
                     JXF_Int *R_max_size_ptr);

JXF_Int
jxf_PAMGBuildRestrDist2AIR(jxf_ParCSRMatrix *A,
                          JXF_Int *CF_marker,
                          jxf_ParCSRMatrix *S,
                          JXF_Int *num_cpts_global,
                          JXF_Int num_functions,
                          JXF_Int *dof_func,
                          JXF_Int debug_flag,
                          JXF_Real trunc_factor,
                          JXF_Int max_elmts,
                          JXF_Int *col_offd_S_to_A,
                          jxf_ParCSRMatrix **R_ptr,
                          JXF_Int *R_max_size_ptr);
JXF_Int
jxf_PAMGBuildOPntRestr(jxf_ParCSRMatrix *A, JXF_Int *CF_marker, JXF_Int *num_cpts_global, jxf_ParCSRMatrix **R_ptr, JXF_Int *R_max_size_ptr);

/* csrc/amg/par_sthength.c */
JXF_Int
jxf_PAMGCreateS(jxf_ParCSRMatrix *par_A,
               JXF_Real strength_threshold,
               JXF_Real max_row_sum,
               JXF_Int num_functions,
               JXF_Int *dof_func,
               jxf_ParCSRMatrix **S_ptr);

JXF_Int
jxf_PAMGCreateSabs(jxf_ParCSRMatrix *A,
                  JXF_Real strength_threshold,
                  JXF_Real max_row_sum,
                  JXF_Int num_functions,
                  JXF_Int *dof_func,
                  jxf_ParCSRMatrix **S_ptr);

JXF_Int jxf_PAMGCreateSCommPkg(jxf_ParCSRMatrix *hp_A, jxf_ParCSRMatrix *par_S, JXF_Int **col_offd_S_to_A_ptr);

JXF_Int
jxf_PAMGCreate2ndS(jxf_ParCSRMatrix *S,
                  JXF_Int *CF_marker,
                  JXF_Int num_paths,
                  JXF_Int *coarse_row_starts,
                  jxf_ParCSRMatrix **C_ptr);

JXF_Int jxf_PAMGCorrectCFMarker(JXF_Int *CF_marker, JXF_Int num_var, JXF_Int *new_CF_marker);

JXF_Int jxf_Mumps(jxf_hpCSRMatrix *pA, jxf_ParVector *pb, jxf_ParVector *px);

#if JXF_USING_MUMPS
void jxfvv_mumps_(int *n, int *nnz, int *ia, int *ja, double *aa, double *rhs, double *sol);
#endif

#endif
