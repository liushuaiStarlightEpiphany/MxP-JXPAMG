//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jx_ilu.h -- head files for ILU Preconditioner
 *  Date: 2014/03/24
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JX_ILU_HEADER
#define JX_ILU_HEADER

#ifndef JX_UTIL_HEADER
#include "jx_util.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jx_GridPartitionData
 */
typedef struct
{
    MPI_Comm comm;

    JX_Int part_type;
    JX_Int x_part_len;
    JX_Int y_part_len;
    JX_Int x_lower_idx;
    JX_Int x_upper_idx;
    JX_Int y_lower_idx;
    JX_Int y_upper_idx;
    JX_Int num_sideproc;
    JX_Int num_smallside;
    JX_Int num_largeside;
    JX_Int num_smallcross;
    JX_Int num_largecross;
    JX_Int num_nocrossside;

    JX_Int *xlo_array;
    JX_Int *xup_array;
    JX_Int *ylo_array;
    JX_Int *yup_array;
    JX_Int *sideprocs;
    JX_Int *sideprcpos;
    JX_Int *sideprcxsrt;
    JX_Int *sideprcysrt;
    JX_Int *sideprclength;

} jx_GridPartitionData;

/*!
 * \struct jx_ILUZeroFactorData
 */
typedef struct
{
    MPI_Comm comm;

    JX_Int nx;
    JX_Int ny;
    JX_Int npx;
    JX_Int npy;
    JX_Int num_equns;

    JX_Int ex_len;
    JX_Int ey_len;
    JX_Int dx_len;
    JX_Int kx_len;
    JX_Int dy_len;
    JX_Int lx_len;
    JX_Int ly_len;
    JX_Int length;
    JX_Int max_iter;
    JX_Int pos_int_end;
    JX_Int pos_dwn_end;
    JX_Int pos_lft_end;
    JX_Int num_fill_in_drop;

    JX_Real drop_tol;

    JX_Int *index;
    JX_Int *indexA;
    JX_Int *indexD;
    JX_Int *permute;

    JX_Real *value;
    JX_Real *senddown;

    MPI_Status *status;

    jx_CSRMatrix *matA;

    jx_Vector *aux_vec;
    jx_Vector *res_vec;
    jx_Vector *tmp_vec;

    jx_ParVector *par_aux_vec;
    jx_ParVector *par_res_vec;

    jx_GridPartitionData *par_grid;

} jx_ILUZeroFactorData;

typedef struct jx_ParILUData_struct
{
    /* General data */
    JX_Int global_solver;
    jx_ParCSRMatrix *matA;
    jx_ParCSRMatrix *matL;
    JX_Real *matD;
    jx_ParCSRMatrix *matU;
    jx_ParCSRMatrix *matmL;
    JX_Real *matmD;
    jx_ParCSRMatrix *matmU;
    jx_ParCSRMatrix *matS;
    JX_Real *droptol; /* Array of 3 elements, for B, (E and F), S respectively */
    JX_Int lfil;
    JX_Int maxRowNnz;
    JX_Int *CF_marker_array;
    JX_Int *perm;
    JX_Int *qperm;
    JX_Real tol_ddPQ;
    jx_ParVector *F;
    jx_ParVector *U;
    jx_ParVector *residual;
    JX_Real *rel_res_norms;
    JX_Int num_iterations;
    JX_Real *l1_norms;
    JX_Real final_rel_residual_norm;
    JX_Real tol;
    JX_Real operator_complexity;
    JX_Int logging;
    JX_Int print_level;
    JX_Int max_iter;
    JX_Int tri_solve;
    JX_Int lower_jacobi_iters;
    JX_Int upper_jacobi_iters;
    JX_Int lrc_type;
    JX_Int lrc_max_iter;
    JX_Real lrc_threshold;
    JX_Real lrc_omega;
    JX_Int lrc_select_type;
    JX_Real lrc_energy_target;
    JX_Int *lrc_marker;
    JX_Int lrc_num_selected;
    JX_Int lrc_core_size;
    JX_Int *lrc_core_indices;
    JX_Int *lrc_local_indices;
    JX_Int *lrc_local_map;
    jx_CSRMatrix *lrc_local_mat;
    JX_Int *lrc_local_indexD;
    JX_Int *lrc_local_indexLU;
    JX_Real *lrc_local_valueLU;
    JX_Real *lrc_local_rhs;
    JX_Real *lrc_local_sol;
    JX_Real *lrc_local_work;
    JX_Real lrc_global_core_size;
    JX_Real lrc_global_expanded_size;
    JX_Real lrc_global_local_mat_nnz;
    JX_Real lrc_global_local_ilu_nnz;
    JX_Real lrc_local_mat_complexity;
    JX_Real lrc_local_ilu_complexity;
    JX_Real lrc_extra_complexity;
    JX_Real lrc_storage_complexity;
    JX_Real lrc_setup_work_complexity;
    JX_Real lrc_effective_complexity;
    JX_Int lrc_apply_count;
    JX_Real lrc_apply_time;
    JX_Real lrc_matvec_time;
    JX_Real lrc_local_solve_time;
    JX_Real lrc_probe_core_energy_ratio;
    JX_Real lrc_probe_expanded_energy_ratio;
    JX_Real lrc_offd_coupling_ratio;
    JX_Real lrc_delta_norm_ratio_sum;
    JX_Real lrc_residual_ratio_sum;
    JX_Real lrc_residual_ratio_min;
    JX_Real lrc_residual_ratio_max;
    JX_Int lrc_accept_like_count;
    JX_Int lrc_reject_like_count;
    JX_Int ilu_type;
    JX_Int nLU;
    JX_Int nI;
    JX_Int *u_end; /* used when schur block is formed */

    /* Iterative ILU parameters */
    JX_Int iter_setup_type;
    JX_Int iter_setup_option;
    JX_Int setup_max_iter;
    JX_Int setup_num_iter;
    JX_Real setup_tolerance;
    JX_Complex *setup_history;

    /* temp vectors for solve phase */
    jx_ParVector *Utemp;
    jx_ParVector *Ftemp;
    jx_ParVector *Xtemp;
    jx_ParVector *Ytemp;
    jx_ParVector *Ztemp;
    JX_Real *uext;
    JX_Real *fext;

    /* data structure sor solving Schur System */
    JX_Solver schur_solver;
    JX_Solver schur_precond;
    jx_ParVector *rhs;
    jx_ParVector *x;

    /* Schur solver data */
    JX_Int ss_logging;
    JX_Int ss_print_level;

    /* Schur-GMRES */
    JX_Int ss_kDim;          /* max number of iterations for GMRES */
    JX_Int ss_max_iter;      /* max number of iterations for GMRES solve */
    JX_Real ss_tol;          /* stop iteration tol for GMRES */
    JX_Real ss_absolute_tol; /* absolute tol for GMRES or tol for NSH solve */
    JX_Int ss_rel_change;

    /* Schur-NSH */
    JX_Int ss_nsh_setup_max_iter; /* number of iterations for NSH inverse */
    JX_Int ss_nsh_solve_max_iter; /* max number of iterations for NSH solve */
    JX_Real ss_nsh_setup_tol;     /* stop iteration tol for NSH inverse */
    JX_Real ss_nsh_solve_tol;     /* absolute tol for NSH solve */
    JX_Int ss_nsh_max_row_nnz;    /* max rows of nonzeros for NSH */
    JX_Int ss_nsh_mr_col_version; /* MR column version setting in NSH */
    JX_Int ss_nsh_mr_max_row_nnz; /* max rows for MR  */
    JX_Real *ss_nsh_droptol;      /* droptol array for NSH */
    JX_Int ss_nsh_mr_max_iter;    /* max MR iteration */
    JX_Real ss_nsh_mr_tol;

    /* Schur precond data */
    JX_Int sp_ilu_type;        /* ilu type is use ILU */
    JX_Int sp_ilu_lfil;        /* level of fill in for ILUK */
    JX_Int sp_ilu_max_row_nnz; /* max rows for ILUT  */
    /* droptol for ILUT or MR
     * ILUT: [0], [1], [2] B, E&F, S respectively
     * NSH: [0] for MR, [1] for NSH
     */
    JX_Real *sp_ilu_droptol; /* droptol array for ILUT */
    JX_Int sp_print_level;
    JX_Int sp_max_iter; /* max precond iter or max MR iteration */
    JX_Int sp_tri_solve;
    JX_Int sp_lower_jacobi_iters;
    JX_Int sp_upper_jacobi_iters;
    JX_Real sp_tol;
    JX_Int test_opt; /* TODO (VPM): change this to something more descriptive*/

    /* local reordering */
    JX_Int reordering_type;
} jx_ParILUData;

#define jx_ParILUDataTestOption(ilu_data) ((ilu_data)->test_opt)
#define jx_ParILUDataGlobalSolver(ilu_data) ((ilu_data)->global_solver)
#define jx_ParILUDataMatA(ilu_data) ((ilu_data)->matA)
#define jx_ParILUDataMatL(ilu_data) ((ilu_data)->matL)
#define jx_ParILUDataMatD(ilu_data) ((ilu_data)->matD)
#define jx_ParILUDataMatU(ilu_data) ((ilu_data)->matU)
#define jx_ParILUDataMatLModified(ilu_data) ((ilu_data)->matmL)
#define jx_ParILUDataMatDModified(ilu_data) ((ilu_data)->matmD)
#define jx_ParILUDataMatUModified(ilu_data) ((ilu_data)->matmU)
#define jx_ParILUDataMatS(ilu_data) ((ilu_data)->matS)
#define jx_ParILUDataDroptol(ilu_data) ((ilu_data)->droptol)
#define jx_ParILUDataLfil(ilu_data) ((ilu_data)->lfil)
#define jx_ParILUDataMaxRowNnz(ilu_data) ((ilu_data)->maxRowNnz)
#define jx_ParILUDataCFMarkerArray(ilu_data) ((ilu_data)->CF_marker_array)
#define jx_ParILUDataPerm(ilu_data) ((ilu_data)->perm)
#define jx_ParILUDataPPerm(ilu_data) ((ilu_data)->perm)
#define jx_ParILUDataQPerm(ilu_data) ((ilu_data)->qperm)
#define jx_ParILUDataTolDDPQ(ilu_data) ((ilu_data)->tol_ddPQ)
#define jx_ParILUDataF(ilu_data) ((ilu_data)->F)
#define jx_ParILUDataU(ilu_data) ((ilu_data)->U)
#define jx_ParILUDataResidual(ilu_data) ((ilu_data)->residual)
#define jx_ParILUDataRelResNorms(ilu_data) ((ilu_data)->rel_res_norms)
#define jx_ParILUDataNumIterations(ilu_data) ((ilu_data)->num_iterations)
#define jx_ParILUDataL1Norms(ilu_data) ((ilu_data)->l1_norms)
#define jx_ParILUDataFinalRelResidualNorm(ilu_data) ((ilu_data)->final_rel_residual_norm)
#define jx_ParILUDataTol(ilu_data) ((ilu_data)->tol)
#define jx_ParILUDataOperatorComplexity(ilu_data) ((ilu_data)->operator_complexity)
#define jx_ParILUDataLogging(ilu_data) ((ilu_data)->logging)
#define jx_ParILUDataPrintLevel(ilu_data) ((ilu_data)->print_level)
#define jx_ParILUDataMaxIter(ilu_data) ((ilu_data)->max_iter)
#define jx_ParILUDataTriSolve(ilu_data) ((ilu_data)->tri_solve)
#define jx_ParILUDataLowerJacobiIters(ilu_data) ((ilu_data)->lower_jacobi_iters)
#define jx_ParILUDataUpperJacobiIters(ilu_data) ((ilu_data)->upper_jacobi_iters)
#define jx_ParILUDataLocalResidualCorrectionType(ilu_data) ((ilu_data)->lrc_type)
#define jx_ParILUDataLocalResidualCorrectionMaxIter(ilu_data) ((ilu_data)->lrc_max_iter)
#define jx_ParILUDataLocalResidualCorrectionThreshold(ilu_data) ((ilu_data)->lrc_threshold)
#define jx_ParILUDataLocalResidualCorrectionOmega(ilu_data) ((ilu_data)->lrc_omega)
#define jx_ParILUDataLocalResidualCorrectionSelectType(ilu_data) ((ilu_data)->lrc_select_type)
#define jx_ParILUDataLocalResidualCorrectionEnergyTarget(ilu_data) ((ilu_data)->lrc_energy_target)
#define jx_ParILUDataLocalResidualCorrectionMarker(ilu_data) ((ilu_data)->lrc_marker)
#define jx_ParILUDataLocalResidualCorrectionNumSelected(ilu_data) ((ilu_data)->lrc_num_selected)
#define jx_ParILUDataLocalResidualCorrectionCoreSize(ilu_data) ((ilu_data)->lrc_core_size)
#define jx_ParILUDataLocalResidualCorrectionCoreIndices(ilu_data) ((ilu_data)->lrc_core_indices)
#define jx_ParILUDataLocalResidualCorrectionLocalIndices(ilu_data) ((ilu_data)->lrc_local_indices)
#define jx_ParILUDataLocalResidualCorrectionLocalMap(ilu_data) ((ilu_data)->lrc_local_map)
#define jx_ParILUDataLocalResidualCorrectionLocalMat(ilu_data) ((ilu_data)->lrc_local_mat)
#define jx_ParILUDataLocalResidualCorrectionLocalIndexD(ilu_data) ((ilu_data)->lrc_local_indexD)
#define jx_ParILUDataLocalResidualCorrectionLocalIndexLU(ilu_data) ((ilu_data)->lrc_local_indexLU)
#define jx_ParILUDataLocalResidualCorrectionLocalValueLU(ilu_data) ((ilu_data)->lrc_local_valueLU)
#define jx_ParILUDataLocalResidualCorrectionLocalRhs(ilu_data) ((ilu_data)->lrc_local_rhs)
#define jx_ParILUDataLocalResidualCorrectionLocalSol(ilu_data) ((ilu_data)->lrc_local_sol)
#define jx_ParILUDataLocalResidualCorrectionLocalWork(ilu_data) ((ilu_data)->lrc_local_work)
#define jx_ParILUDataLocalResidualCorrectionGlobalCoreSize(ilu_data) ((ilu_data)->lrc_global_core_size)
#define jx_ParILUDataLocalResidualCorrectionGlobalExpandedSize(ilu_data) ((ilu_data)->lrc_global_expanded_size)
#define jx_ParILUDataLocalResidualCorrectionGlobalLocalMatNnz(ilu_data) ((ilu_data)->lrc_global_local_mat_nnz)
#define jx_ParILUDataLocalResidualCorrectionGlobalLocalIluNnz(ilu_data) ((ilu_data)->lrc_global_local_ilu_nnz)
#define jx_ParILUDataLocalResidualCorrectionLocalMatComplexity(ilu_data) ((ilu_data)->lrc_local_mat_complexity)
#define jx_ParILUDataLocalResidualCorrectionLocalIluComplexity(ilu_data) ((ilu_data)->lrc_local_ilu_complexity)
#define jx_ParILUDataLocalResidualCorrectionExtraComplexity(ilu_data) ((ilu_data)->lrc_extra_complexity)
#define jx_ParILUDataLocalResidualCorrectionStorageComplexity(ilu_data) ((ilu_data)->lrc_storage_complexity)
#define jx_ParILUDataLocalResidualCorrectionSetupWorkComplexity(ilu_data) ((ilu_data)->lrc_setup_work_complexity)
#define jx_ParILUDataLocalResidualCorrectionEffectiveComplexity(ilu_data) ((ilu_data)->lrc_effective_complexity)
#define jx_ParILUDataLocalResidualCorrectionApplyCount(ilu_data) ((ilu_data)->lrc_apply_count)
#define jx_ParILUDataLocalResidualCorrectionApplyTime(ilu_data) ((ilu_data)->lrc_apply_time)
#define jx_ParILUDataLocalResidualCorrectionMatvecTime(ilu_data) ((ilu_data)->lrc_matvec_time)
#define jx_ParILUDataLocalResidualCorrectionLocalSolveTime(ilu_data) ((ilu_data)->lrc_local_solve_time)
#define jx_ParILUDataLocalResidualCorrectionProbeCoreEnergyRatio(ilu_data) ((ilu_data)->lrc_probe_core_energy_ratio)
#define jx_ParILUDataLocalResidualCorrectionProbeExpandedEnergyRatio(ilu_data) ((ilu_data)->lrc_probe_expanded_energy_ratio)
#define jx_ParILUDataLocalResidualCorrectionOffdCouplingRatio(ilu_data) ((ilu_data)->lrc_offd_coupling_ratio)
#define jx_ParILUDataLocalResidualCorrectionDeltaNormRatioSum(ilu_data) ((ilu_data)->lrc_delta_norm_ratio_sum)
#define jx_ParILUDataLocalResidualCorrectionResidualRatioSum(ilu_data) ((ilu_data)->lrc_residual_ratio_sum)
#define jx_ParILUDataLocalResidualCorrectionResidualRatioMin(ilu_data) ((ilu_data)->lrc_residual_ratio_min)
#define jx_ParILUDataLocalResidualCorrectionResidualRatioMax(ilu_data) ((ilu_data)->lrc_residual_ratio_max)
#define jx_ParILUDataLocalResidualCorrectionAcceptLikeCount(ilu_data) ((ilu_data)->lrc_accept_like_count)
#define jx_ParILUDataLocalResidualCorrectionRejectLikeCount(ilu_data) ((ilu_data)->lrc_reject_like_count)
#define jx_ParILUDataIluType(ilu_data) ((ilu_data)->ilu_type)
#define jx_ParILUDataNLU(ilu_data) ((ilu_data)->nLU)
#define jx_ParILUDataNI(ilu_data) ((ilu_data)->nI)
#define jx_ParILUDataUEnd(ilu_data) ((ilu_data)->u_end)
#define jx_ParILUDataUTemp(ilu_data) ((ilu_data)->Utemp)
#define jx_ParILUDataFTemp(ilu_data) ((ilu_data)->Ftemp)
#define jx_ParILUDataXTemp(ilu_data) ((ilu_data)->Xtemp)
#define jx_ParILUDataYTemp(ilu_data) ((ilu_data)->Ytemp)
#define jx_ParILUDataZTemp(ilu_data) ((ilu_data)->Ztemp)
#define jx_ParILUDataUExt(ilu_data) ((ilu_data)->uext)
#define jx_ParILUDataFExt(ilu_data) ((ilu_data)->fext)
#define jx_ParILUDataSchurSolver(ilu_data) ((ilu_data)->schur_solver)
#define jx_ParILUDataSchurPrecond(ilu_data) ((ilu_data)->schur_precond)
#define jx_ParILUDataRhs(ilu_data) ((ilu_data)->rhs)
#define jx_ParILUDataX(ilu_data) ((ilu_data)->x)
#define jx_ParILUDataReorderingType(ilu_data) ((ilu_data)->reordering_type)

/* Iterative ILU setup */
#define jx_ParILUDataIterativeSetupType(ilu_data) ((ilu_data)->iter_setup_type)
#define jx_ParILUDataIterativeSetupOption(ilu_data) ((ilu_data)->iter_setup_option)
#define jx_ParILUDataIterativeSetupMaxIter(ilu_data) ((ilu_data)->setup_max_iter)
#define jx_ParILUDataIterativeSetupNumIter(ilu_data) ((ilu_data)->setup_num_iter)
#define jx_ParILUDataIterativeSetupTolerance(ilu_data) ((ilu_data)->setup_tolerance)
#define jx_ParILUDataIterativeSetupHistory(ilu_data) ((ilu_data)->setup_history)
#define jx_ParILUDataIterSetupCorrectionNorm(ilu_data, i) ((ilu_data)->setup_history[i])
#define jx_ParILUDataIterSetupResidualNorm(ilu_data, i) (((ilu_data)->setup_history + \
                                                          (ilu_data)->setup_num_iter)[i])
/* Schur System */
#define jx_ParILUDataSchurGMRESKDim(ilu_data) ((ilu_data)->ss_kDim)
#define jx_ParILUDataSchurGMRESMaxIter(ilu_data) ((ilu_data)->ss_max_iter)
#define jx_ParILUDataSchurGMRESTol(ilu_data) ((ilu_data)->ss_tol)
#define jx_ParILUDataSchurGMRESAbsoluteTol(ilu_data) ((ilu_data)->ss_absolute_tol)
#define jx_ParILUDataSchurGMRESRelChange(ilu_data) ((ilu_data)->ss_rel_change)
#define jx_ParILUDataSchurPrecondIluType(ilu_data) ((ilu_data)->sp_ilu_type)
#define jx_ParILUDataSchurPrecondIluLfil(ilu_data) ((ilu_data)->sp_ilu_lfil)
#define jx_ParILUDataSchurPrecondIluMaxRowNnz(ilu_data) ((ilu_data)->sp_ilu_max_row_nnz)
#define jx_ParILUDataSchurPrecondIluDroptol(ilu_data) ((ilu_data)->sp_ilu_droptol)
#define jx_ParILUDataSchurPrecondPrintLevel(ilu_data) ((ilu_data)->sp_print_level)
#define jx_ParILUDataSchurPrecondMaxIter(ilu_data) ((ilu_data)->sp_max_iter)
#define jx_ParILUDataSchurPrecondTriSolve(ilu_data) ((ilu_data)->sp_tri_solve)
#define jx_ParILUDataSchurPrecondLowerJacobiIters(ilu_data) ((ilu_data)->sp_lower_jacobi_iters)
#define jx_ParILUDataSchurPrecondUpperJacobiIters(ilu_data) ((ilu_data)->sp_upper_jacobi_iters)
#define jx_ParILUDataSchurPrecondTol(ilu_data) ((ilu_data)->sp_tol)

#define jx_ParILUDataSchurNSHMaxNumIter(ilu_data) ((ilu_data)->ss_nsh_setup_max_iter)
#define jx_ParILUDataSchurNSHSolveMaxIter(ilu_data) ((ilu_data)->ss_nsh_solve_max_iter)
#define jx_ParILUDataSchurNSHTol(ilu_data) ((ilu_data)->ss_nsh_setup_tol)
#define jx_ParILUDataSchurNSHSolveTol(ilu_data) ((ilu_data)->ss_nsh_solve_tol)
#define jx_ParILUDataSchurNSHMaxRowNnz(ilu_data) ((ilu_data)->ss_nsh_max_row_nnz)
#define jx_ParILUDataSchurMRColVersion(ilu_data) ((ilu_data)->ss_nsh_mr_col_version)
#define jx_ParILUDataSchurMRMaxRowNnz(ilu_data) ((ilu_data)->ss_nsh_mr_max_row_nnz)
#define jx_ParILUDataSchurNSHDroptol(ilu_data) ((ilu_data)->ss_nsh_droptol)
#define jx_ParILUDataSchurMRMaxIter(ilu_data) ((ilu_data)->ss_nsh_mr_max_iter)
#define jx_ParILUDataSchurMRTol(ilu_data) ((ilu_data)->ss_nsh_mr_tol)

#define jx_ParILUDataSchurSolverLogging(ilu_data) ((ilu_data)->ss_logging)
#define jx_ParILUDataSchurSolverPrintLevel(ilu_data) ((ilu_data)->ss_print_level)

// ===========================================================

#define FMRK -1
#define CMRK 1
#define UMRK 0
#define S_CMRK 2

#define FPT(i, bsize) (((i) % (bsize)) == FMRK)
#define CPT(i, bsize) (((i) % (bsize)) == CMRK)

#define MAT_TOL 1e-14
#define EXPAND_FACT 1.3

typedef struct jx_ParNSHData_struct
{
    /* solver information */
    JX_Int global_solver;
    jx_ParCSRMatrix *matA;
    jx_ParCSRMatrix *matM;
    jx_ParVector *F;
    jx_ParVector *U;
    jx_ParVector *residual;
    JX_Real *rel_res_norms;
    JX_Int num_iterations;
    JX_Real *l1_norms;
    JX_Real final_rel_residual_norm;
    JX_Real tol;
    JX_Real operator_complexity;

    JX_Int logging;
    JX_Int print_level;
    JX_Int max_iter;

    /* common data slots */
    /* droptol[0]: droptol for MR
     * droptol[1]: droptol for NSH
     */
    JX_Real *droptol;
    JX_Int own_droptol_data;

    /* temp vectors for solve phase */
    jx_ParVector *Utemp;
    jx_ParVector *Ftemp;

    /* data slots for local MR */
    JX_Int mr_max_iter;
    JX_Real mr_tol;
    JX_Int mr_max_row_nnz;
    JX_Int mr_col_version; /* global version or column version MR */

    /* data slots for global NSH */
    JX_Int nsh_max_iter;
    JX_Real nsh_tol;
    JX_Int nsh_max_row_nnz;

} jx_ParNSHData;

#define jx_ParNSHDataGlobalSolver(nsh_data) ((nsh_data)->global_solver)
#define jx_ParNSHDataMatA(nsh_data) ((nsh_data)->matA)
#define jx_ParNSHDataMatM(nsh_data) ((nsh_data)->matM)
#define jx_ParNSHDataF(nsh_data) ((nsh_data)->F)
#define jx_ParNSHDataU(nsh_data) ((nsh_data)->U)
#define jx_ParNSHDataResidual(nsh_data) ((nsh_data)->residual)
#define jx_ParNSHDataRelResNorms(nsh_data) ((nsh_data)->rel_res_norms)
#define jx_ParNSHDataNumIterations(nsh_data) ((nsh_data)->num_iterations)
#define jx_ParNSHDataL1Norms(nsh_data) ((nsh_data)->l1_norms)
#define jx_ParNSHDataFinalRelResidualNorm(nsh_data) ((nsh_data)->final_rel_residual_norm)
#define jx_ParNSHDataTol(nsh_data) ((nsh_data)->tol)
#define jx_ParNSHDataOperatorComplexity(nsh_data) ((nsh_data)->operator_complexity)
#define jx_ParNSHDataLogging(nsh_data) ((nsh_data)->logging)
#define jx_ParNSHDataPrintLevel(nsh_data) ((nsh_data)->print_level)
#define jx_ParNSHDataMaxIter(nsh_data) ((nsh_data)->max_iter)
#define jx_ParNSHDataDroptol(nsh_data) ((nsh_data)->droptol)
#define jx_ParNSHDataOwnDroptolData(nsh_data) ((nsh_data)->own_droptol_data)
#define jx_ParNSHDataUTemp(nsh_data) ((nsh_data)->Utemp)
#define jx_ParNSHDataFTemp(nsh_data) ((nsh_data)->Ftemp)
#define jx_ParNSHDataMRMaxIter(nsh_data) ((nsh_data)->mr_max_iter)
#define jx_ParNSHDataMRTol(nsh_data) ((nsh_data)->mr_tol)
#define jx_ParNSHDataMRMaxRowNnz(nsh_data) ((nsh_data)->mr_max_row_nnz)
#define jx_ParNSHDataMRColVersion(nsh_data) ((nsh_data)->mr_col_version)
#define jx_ParNSHDataNSHMaxIter(nsh_data) ((nsh_data)->nsh_max_iter)
#define jx_ParNSHDataNSHTol(nsh_data) ((nsh_data)->nsh_tol)
#define jx_ParNSHDataNSHMaxRowNnz(nsh_data) ((nsh_data)->nsh_max_row_nnz)

#define DIVIDE_TOL 1e-32

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/ilu/csr_matrix.c */
jx_CSRMatrix *jx_CSRMatrixTDMGReorderByNodes(jx_CSRMatrix *A, JX_Int num_equns);
jx_CSRMatrix *jx_CSRMatrixTDMGReorderByVariables(jx_CSRMatrix *A, JX_Int num_equns);
jx_CSRMatrix *jx_CSRMatrixTDMGReorderAlongX(jx_CSRMatrix *A, JX_Int num_equns, JX_Int nx, JX_Int ny);
jx_CSRMatrix *jx_CSRMatrixTDMGReorderAlongY(jx_CSRMatrix *A, JX_Int num_equns, JX_Int nx, JX_Int ny);

/* csrc/ilu/decomposition.c */
JX_Int jx_ILUZeroDecompositionA(jx_CSRMatrix *A,
                                JX_Int **indexDP_ptr,
                                JX_Int **indexLU_ptr,
                                JX_Real **valueLU_ptr);
JX_Int jx_ILUZeroDecompositionB(jx_CSRMatrix *A,
                                JX_Real drop_tol,
                                JX_Int **indexAP_ptr,
                                JX_Int **indexDP_ptr,
                                JX_Int **indexLU_ptr,
                                JX_Real **valueLU_ptr);
JX_Int *jx_ILUZeroParallelDecompositionA(jx_ParCSRMatrix *par_A,
                                         JX_Int **indexDP_ptr,
                                         JX_Int **indexLU_ptr,
                                         JX_Real **valueLU_ptr,
                                         JX_Int *num_nonzeros,
                                         JX_Int *fill_in_drop,
                                         jx_ILUZeroFactorData *ilu_data);
JX_Int
jx_ILUZeroLocalDecompositionIntURPntsA(JX_Int *IA,
                                       JX_Int *JA,
                                       JX_Real *AA,
                                       JX_Int *indexLU,
                                       JX_Int *indexDP,
                                       JX_Int *placeRC,
                                       JX_Real *valueLU,
                                       JX_Int int_uprgt_pnt,
                                       JX_Int first_row_idx);
void jx_ILUZeroLocalDecompositionDPntsA(JX_Int *IA,
                                        JX_Int *JA,
                                        JX_Real *AA,
                                        JX_Int *indexLU,
                                        JX_Int *indexDP,
                                        JX_Int *placeRC,
                                        JX_Real *valueLU,
                                        JX_Int int_uprgt_pnt,
                                        JX_Int first_row_idx,
                                        JX_Int num_rows,
                                        JX_Int ex_len,
                                        JX_Int *fill_in_drop);
void jx_ILUZeroLocalDecompositionLPntsA(JX_Int *IA,
                                        JX_Int *JA,
                                        JX_Real *AA,
                                        JX_Int *indexLU,
                                        JX_Int *indexDP,
                                        JX_Int *placeRC,
                                        JX_Real *valueLU,
                                        JX_Int int_uprgtdwn_pnt,
                                        JX_Int first_row_idx,
                                        JX_Int num_rows,
                                        JX_Int dwn_num_rows,
                                        JX_Int ex_len,
                                        JX_Int *fill_in_drop);
void jx_ILUZeroLocalDecompositionLDPntsA(JX_Int *IA,
                                         JX_Int *JA,
                                         JX_Real *AA,
                                         JX_Int *indexLU,
                                         JX_Int *indexDP,
                                         JX_Int *placeRC,
                                         JX_Real *valueLU,
                                         JX_Int int_urdl_pnt,
                                         JX_Int first_row_idx,
                                         JX_Int fst_row_idx,
                                         JX_Int num_rows,
                                         JX_Int lft_cnum_rows,
                                         JX_Int dwn_cnum_rows,
                                         JX_Int ex_len,
                                         JX_Int *fill_in_drop);
void jx_ILUZeroFactorDataParallelUPartIntURPntsA(JX_Int *indexLU,
                                                 JX_Int *indexAP,
                                                 JX_Int *indexDP,
                                                 JX_Real *valueLU,
                                                 JX_Int int_uprgt_pnt,
                                                 JX_Int recv_downsrt,
                                                 JX_Int recv_leftsrt,
                                                 JX_Int *permute,
                                                 JX_Real *aux_data,
                                                 JX_Real *res_data,
                                                 JX_Real *app_data,
                                                 JX_Int first_row_idx,
                                                 JX_Int next_row_idx,
                                                 JX_Int nnxt_row_idx);
void jx_ILUZeroFactorDataParallelUPartIntURPntsB(JX_Int *indexLU,
                                                 JX_Int *indexAP,
                                                 JX_Int *indexDP,
                                                 JX_Real *valueLU,
                                                 JX_Int int_uprgt_pnt,
                                                 JX_Int recv_downsrt,
                                                 JX_Int *permute,
                                                 JX_Real *aux_data,
                                                 JX_Real *res_data,
                                                 JX_Real *app_data,
                                                 JX_Int first_row_idx,
                                                 JX_Int next_row_idx);
void jx_ILUZeroFactorDataParallelUPartIntURPntsC(JX_Int *indexLU,
                                                 JX_Int *indexAP,
                                                 JX_Int *indexDP,
                                                 JX_Real *valueLU,
                                                 JX_Int int_uprgt_pnt,
                                                 JX_Int recv_leftsrt,
                                                 JX_Int *permute,
                                                 JX_Real *aux_data,
                                                 JX_Real *res_data,
                                                 JX_Real *app_data,
                                                 JX_Int first_row_idx,
                                                 JX_Int next_row_idx);
void jx_ILUZeroFactorDataParallelUPartIntURPntsD(JX_Int *indexLU,
                                                 JX_Int *indexAP,
                                                 JX_Int *indexDP,
                                                 JX_Real *valueLU,
                                                 JX_Int int_uprgt_pnt,
                                                 JX_Int *permute,
                                                 JX_Real *aux_data,
                                                 JX_Real *res_data,
                                                 JX_Real *app_data,
                                                 JX_Int first_row_idx);

/* csrc/ilu/grid.c */
JX_Int JX_GridPartitionDataCreate(JX_Solver *solver, MPI_Comm comm, JX_Int nx, JX_Int ny, JX_Int npx, JX_Int npy);
JX_Int JX_GridPartitionDataSetEachSides4Comm(JX_Solver solver);
JX_Int JX_GridPartitionDataDestroy(JX_Solver solver);
void *jx_GridPartitionDataInitialize(MPI_Comm comm, JX_Int nx, JX_Int ny, JX_Int npx, JX_Int npy);
JX_Int jx_GridPartitionDataSetEachSides4Comm(void *grid_vdata);
JX_Int jx_GridPartitionDataFinalize(void *grid_vdata);

/* csrc/ilu/ilu.c */
JX_Int JX_ILUZeroFactorDataCreate(JX_Solver *solver, MPI_Comm comm);
JX_Int JX_ILUZeroFactorDataSetMaxIter(JX_Solver solver, JX_Int max_iter);
JX_Int JX_ILUZeroFactorDataSetNxy(JX_Solver solver, JX_Int nx, JX_Int ny);
JX_Int JX_ILUZeroFactorDataSetNpxy(JX_Solver solver, JX_Int npx, JX_Int npy);
JX_Int JX_ILUZeroFactorDataSetNumEquns(JX_Solver solver, JX_Int num_equns);
JX_Int JX_ILUZeroFactorDataSetDropTol(JX_Solver solver, JX_Real drop_tol);
JX_Int JX_ILUZeroFactorDataSetMatA(JX_Solver solver, jx_CSRMatrix *matA);
JX_Int JX_ILUZeroFactorDataGetLULength(JX_Solver solver, JX_Int *lu_length);
JX_Int JX_ILUZeroFactorDataDestroy(JX_Solver solver);
JX_Int JX_ILUZeroFactorDataGenerateParGrid(JX_Solver solver);
JX_Int JX_ILUZeroFactorDataSetup(JX_Solver solver, JX_hpCSRMatrix par_matrix);
JX_Int JX_ILUZeroFactorDataPrecond(JX_Solver solver,
                                   JX_hpCSRMatrix par_matrix,
                                   JX_ParVector par_rhs,
                                   JX_ParVector par_app);
void *jx_ILUZeroFactorDataInitialize(MPI_Comm comm);
JX_Int jx_ILUZeroFactorDataSetMaxIter(void *ilu_vdata, JX_Int max_iter);
JX_Int jx_ILUZeroFactorDataSetNxy(void *ilu_vdata, JX_Int nx, JX_Int ny);
JX_Int jx_ILUZeroFactorDataSetNpxy(void *ilu_vdata, JX_Int npx, JX_Int npy);
JX_Int jx_ILUZeroFactorDataSetNumEquns(void *ilu_vdata, JX_Int num_equns);
JX_Int jx_ILUZeroFactorDataSetDropTol(void *ilu_vdata, JX_Real drop_tol);
JX_Int jx_ILUZeroFactorDataSetMatA(void *ilu_vdata, jx_CSRMatrix *matA);
JX_Int jx_ILUZeroFactorDataGetLULength(void *ilu_vdata, JX_Int *lu_length);
JX_Int jx_ILUZeroFactorDataFinalize(void *ilu_vdata);
JX_Int jx_ILUZeroFactorDataGenerateParGrid(void *ilu_vdata);

/* csrc/ilu/ilucycle.c */
JX_Int jx_ILUZeroFactorDataPrecond(void *ilu_vdata,
                                   jx_ParCSRMatrix *par_A,
                                   jx_ParVector *par_b,
                                   jx_ParVector *par_x);
void jx_ILUZeroFactorDataCycleA(void *ilu_vdata, jx_CSRMatrix *A, jx_Vector *f, jx_Vector *u);
void jx_ILUZeroFactorDataCycleB(void *ilu_vdata, jx_CSRMatrix *A, jx_Vector *f, jx_Vector *u);
void jx_ILUZeroFactorDataParallelCycleA(jx_ILUZeroFactorData *ilu_data,
                                        jx_ParCSRMatrix *par_A,
                                        jx_ParVector *par_b,
                                        jx_ParVector *par_x);
void jx_ILUZeroFactorDataParallelLPartIntURPntsA(JX_Int *indexLU,
                                                 JX_Int *indexAP,
                                                 JX_Int *indexDP,
                                                 JX_Real *valueLU,
                                                 JX_Int *permute,
                                                 JX_Int int_uprgt_pnt,
                                                 JX_Real *aux_data,
                                                 JX_Real *rhs_data,
                                                 JX_Int first_row_idx);
void jx_ILUZeroFactorDataParallelLPartDPntsA(JX_Int *indexLU,
                                             JX_Int *indexAP,
                                             JX_Int *indexDP,
                                             JX_Real *valueLU,
                                             JX_Int int_uprgt_pnt,
                                             JX_Int num_rows,
                                             JX_Int *permute,
                                             JX_Real *aux_data,
                                             JX_Real *res_data,
                                             JX_Real *rhs_data,
                                             JX_Int first_row_idx);
void jx_ILUZeroFactorDataParallelLPartLPntsA(JX_Int *indexLU,
                                             JX_Int *indexAP,
                                             JX_Int *indexDP,
                                             JX_Real *valueLU,
                                             JX_Int int_uprgtdwn_pnt,
                                             JX_Int num_rows,
                                             JX_Int dwn_num_rows,
                                             JX_Int *permute,
                                             JX_Real *aux_data,
                                             JX_Real *res_data,
                                             JX_Real *rhs_data,
                                             JX_Int first_row_idx);
void jx_ILUZeroFactorDataParallelLPartLDPntsA(JX_Int *indexLU,
                                              JX_Int *indexAP,
                                              JX_Int *indexDP,
                                              JX_Real *valueLU,
                                              JX_Int int_urdl_pnt,
                                              JX_Int num_rows,
                                              JX_Int fst_row_idx,
                                              JX_Int *permute,
                                              JX_Int lft_cnum_rows,
                                              JX_Int dwn_dnum_rows,
                                              JX_Real *aux_data,
                                              JX_Real *res_data,
                                              JX_Real *rhs_data,
                                              JX_Int first_row_idx);
void jx_ILUZeroFactorDataParallelUPartLDPntsA(JX_Int *indexLU,
                                              JX_Int *indexAP,
                                              JX_Int *indexDP,
                                              JX_Real *valueLU,
                                              JX_Int int_urdl_pnt,
                                              JX_Int num_rows,
                                              JX_Int *permute,
                                              JX_Real *aux_data,
                                              JX_Real *app_data,
                                              JX_Int first_row_idx);
void jx_ILUZeroFactorDataParallelUPartLPntsA(JX_Int *indexLU,
                                             JX_Int *indexAP,
                                             JX_Int *indexDP,
                                             JX_Real *valueLU,
                                             JX_Int int_urdl_pnt,
                                             JX_Int num_rows,
                                             JX_Int recv_downsrt,
                                             JX_Int *permute,
                                             JX_Real *aux_data,
                                             JX_Real *res_data,
                                             JX_Real *app_data,
                                             JX_Int first_row_idx,
                                             JX_Int next_row_idx);
void jx_ILUZeroFactorDataParallelUPartDPntsA(JX_Int *indexLU,
                                             JX_Int *indexAP,
                                             JX_Int *indexDP,
                                             JX_Real *valueLU,
                                             JX_Int int_uprgt_pnt,
                                             JX_Int num_rows,
                                             JX_Int recv_leftsrt,
                                             JX_Int *permute,
                                             JX_Real *aux_data,
                                             JX_Real *res_data,
                                             JX_Real *app_data,
                                             JX_Int first_row_idx,
                                             JX_Int next_row_idx);

/* csrc/ilu/ilusetup.c */
JX_Int jx_ILUZeroFactorDataSetup(void *ilu_vdata, jx_ParCSRMatrix *par_A);

/* csrc/ilu/par_csr_matrix.c */
jx_ParCSRMatrix *jx_BuildMatParFromOneFile2(char *filename,
                                            JX_Int file_base,
                                            JX_Int *row_part,
                                            JX_Int *col_part,
                                            JX_Int num_equns,
                                            JX_Int nx,
                                            JX_Int ny);
jx_ParCSRMatrix *jx_BuildMatParFromOneFile3(char *filename, JX_Int file_base, JX_Int num_equns, JX_Int nx, JX_Int ny);
jx_ParCSRMatrix *
jx_ParCSRMatrixCreate2(MPI_Comm comm,
                       JX_Int global_num_rows,
                       JX_Int global_num_cols,
                       JX_Int *row_starts,
                       JX_Int *col_starts,
                       JX_Int num_cols_offd,
                       JX_Int num_nonzeros_diag,
                       JX_Int num_nonzeros_offd,
                       JX_Int num_equns,
                       JX_Int nx,
                       JX_Int ny);
jx_ParCSRMatrix *jx_CSRMatrixToParCSRMatrix2(MPI_Comm comm, jx_CSRMatrix *A, JX_Int num_equns, JX_Int nx, JX_Int ny);
jx_CSRMatrix *jx_MergeDiagAndOffdDropSmall(jx_ParCSRMatrix *par_matrix, JX_Real drop_tol);
jx_CSRMatrix *
jx_CSRMatrixMergeReorderIntUpRtDnLtQuasiBdy(MPI_Comm comm,
                                            JX_Int *row_starts,
                                            JX_Int *permute,
                                            jx_CSRMatrix *ser_B,
                                            JX_Int nz_srt,
                                            JX_Int ng_pt,
                                            jx_GridPartitionData *grid_data,
                                            JX_Int *ex_len,
                                            JX_Int *ey_len,
                                            JX_Int *dx_len,
                                            JX_Int *kx_len,
                                            JX_Int *dy_len,
                                            JX_Int *lx_len,
                                            JX_Int *ly_len,
                                            JX_Int *postn_a,
                                            JX_Int *postn_b,
                                            JX_Int *postn_c);

/* csrc/ilu/par_vector.c */
jx_ParVector *jx_BuildRhsParFromOneFile2(char *filename, jx_ParCSRMatrix *A, JX_Int num_equns, JX_Int nx, JX_Int ny);
JX_Int jx_GeneratePartitioning2(JX_Int num_procs, JX_Int **part_ptr, JX_Int num_equns, JX_Int length, JX_Int ny);

/* csrc/ilu/vector.c */
jx_Vector *jx_VectorTDMGReorderByNodes(jx_Vector *x, JX_Int num_equns);
jx_Vector *jx_VectorTDMGReorderByVariables(jx_Vector *x, JX_Int num_equns);
jx_Vector *jx_VectorTDMGReorderAlongX(jx_Vector *x, JX_Int num_equns, JX_Int nx, JX_Int ny);
jx_Vector *jx_VectorTDMGReorderAlongY(jx_Vector *x, JX_Int num_equns, JX_Int nx, JX_Int ny);

/* csrc/ilu/par_ilu.c */
JX_Int JX_ILUCreate(JX_Solver *solver);
JX_Int JX_ILUDestroy(JX_Solver solver);
// JX_Int JX_ILUSetup(JX_Solver solver, JX_ParCSRMatrix A, JX_ParVector b, JX_ParVector x);
JX_Int JX_ILUSetup(JX_Solver solver, jx_hpCSRMatrix *A, JX_ParVector *b, JX_ParVector *x);
// JX_Int JX_ILUSolve(JX_Solver solver, JX_ParCSRMatrix *A, JX_ParVector *b, JX_ParVector *x);
JX_Int JX_ILUSolve(JX_Solver solver, jx_hpCSRMatrix *A, JX_ParVector *b, JX_ParVector *x);
JX_Int JX_ILUSetPrintLevel(JX_Solver solver, JX_Int print_level);
JX_Int JX_ILUSetLogging(JX_Solver solver, JX_Int logging);
JX_Int JX_ILUSetMaxIter(JX_Solver solver, JX_Int max_iter);
JX_Int JX_ILUSetIterativeSetupType(JX_Solver solver, JX_Int iter_setup_type);
JX_Int JX_ILUSetTol(JX_Solver solver, JX_Real tol);
JX_Int JX_ILUSetDropThreshold(JX_Solver solver, JX_Real threshold);
JX_Int JX_ILUSetDropThresholdArray(JX_Solver solver, JX_Real *threshold);
JX_Int JX_ILUSetNSHDropThreshold(JX_Solver solver, JX_Real threshold);
JX_Int JX_ILUSetNSHDropThresholdArray(JX_Solver solver, JX_Real *threshold);
JX_Int JX_ILUSetSchurMaxIter(JX_Solver solver, JX_Int ss_max_iter);
JX_Int JX_ILUSetMaxNnzPerRow(JX_Solver solver, JX_Int nzmax);
JX_Int JX_ILUSetLevelOfFill(JX_Solver solver, JX_Int lfil);
JX_Int JX_ILUSetType(JX_Solver solver, JX_Int ilu_type);
JX_Int JX_ILUSetLocalResidualCorrectionType(JX_Solver solver, JX_Int lrc_type);
JX_Int JX_ILUSetLocalResidualCorrectionMaxIter(JX_Solver solver, JX_Int lrc_max_iter);
JX_Int JX_ILUSetLocalResidualCorrectionThreshold(JX_Solver solver, JX_Real lrc_threshold);
JX_Int JX_ILUSetLocalResidualCorrectionOmega(JX_Solver solver, JX_Real lrc_omega);
JX_Int JX_ILUSetLocalResidualCorrectionSelectType(JX_Solver solver, JX_Int lrc_select_type);
JX_Int JX_ILUSetLocalResidualCorrectionEnergyTarget(JX_Solver solver, JX_Real lrc_energy_target);
JX_Int JX_ILUGetNumIterations(JX_Solver solver, JX_Int *num_iterations);
JX_Int JX_ILUGetFinalRelativeResidualNorm(JX_Solver solver, JX_Real *res_norm);
JX_Int JX_ILUSetLocalReordering(JX_Solver solver, JX_Int ordering_type);
void *jx_ILUCreate();
JX_Int jx_ILUDestroy(void *data);
JX_Int jx_ILUSetLevelOfFill(void *ilu_vdata, JX_Int lfil);
JX_Int jx_ILUSetMaxNnzPerRow(void *ilu_vdata, JX_Int nzmax);
JX_Int jx_ILUSetDropThreshold(void *ilu_vdata, JX_Real threshold);
JX_Int jx_ILUSetDropThresholdArray(void *ilu_vdata, JX_Real *threshold);
// JX_Int jx_ILUSetOwnDropThreshold(void *ilu_vdata, JX_Int own_droptol_data);
JX_Int jx_ILUSetType(void *ilu_vdata, JX_Int ilu_type);
JX_Int jx_ILUSetLocalResidualCorrectionType(void *ilu_vdata, JX_Int lrc_type);
JX_Int jx_ILUSetLocalResidualCorrectionMaxIter(void *ilu_vdata, JX_Int lrc_max_iter);
JX_Int jx_ILUSetLocalResidualCorrectionThreshold(void *ilu_vdata, JX_Real lrc_threshold);
JX_Int jx_ILUSetLocalResidualCorrectionOmega(void *ilu_vdata, JX_Real lrc_omega);
JX_Int jx_ILUSetLocalResidualCorrectionSelectType(void *ilu_vdata, JX_Int lrc_select_type);
JX_Int jx_ILUSetLocalResidualCorrectionEnergyTarget(void *ilu_vdata, JX_Real lrc_energy_target);
JX_Int jx_ILUSetMaxIter(void *ilu_vdata, JX_Int max_iter);
JX_Int jx_ILUSetIterativeSetupType(void *ilu_vdata, JX_Int iter_setup_type);
JX_Int jx_ILUSetTol(void *ilu_vdata, JX_Real tol);
JX_Int jx_ILUSetPrintLevel(void *ilu_vdata, JX_Int print_level);
JX_Int jx_ILUSetLogging(void *ilu_vdata, JX_Int logging);
JX_Int jx_ILUSetLocalReordering(void *ilu_vdata, JX_Int ordering_type);
JX_Int jx_ILUSetSchurSolverKDIM(void *ilu_vdata, JX_Int ss_kDim);
JX_Int jx_ILUSetSchurSolverMaxIter(void *ilu_vdata, JX_Int ss_max_iter);
JX_Int jx_ILUSetSchurSolverTol(void *ilu_vdata, JX_Real ss_tol);
JX_Int jx_ILUSetSchurSolverAbsoluteTol(void *ilu_vdata, JX_Real ss_absolute_tol);
JX_Int jx_ILUSetSchurSolverLogging(void *ilu_vdata, JX_Int ss_logging);
JX_Int jx_ILUSetSchurSolverPrintLevel(void *ilu_vdata, JX_Int ss_print_level);
JX_Int jx_ILUSetSchurSolverRelChange(void *ilu_vdata, JX_Int ss_rel_change);
JX_Int jx_ILUSetSchurPrecondILUType(void *ilu_vdata, JX_Int sp_ilu_type);
JX_Int jx_ILUSetSchurPrecondILULevelOfFill(void *ilu_vdata, JX_Int sp_ilu_lfil);
JX_Int jx_ILUSetSchurPrecondILUMaxNnzPerRow(void *ilu_vdata, JX_Int sp_ilu_max_row_nnz);
JX_Int jx_ILUSetSchurPrecondILUDropThreshold(void *ilu_vdata, JX_Real sp_ilu_droptol);
JX_Int jx_ILUSetSchurPrecondILUDropThresholdArray(void *ilu_vdata, JX_Real *sp_ilu_droptol);
JX_Int jx_ILUSetSchurPrecondILUOwnDropThreshold(void *ilu_vdata, JX_Int sp_own_droptol_data);
JX_Int jx_ILUSetSchurPrecondPrintLevel(void *ilu_vdata, JX_Int sp_print_level);
JX_Int jx_ILUSetSchurPrecondMaxIter(void *ilu_vdata, JX_Int sp_max_iter);
JX_Int jx_ILUSetSchurPrecondTol(void *ilu_vdata, JX_Int sp_tol);
JX_Int jx_ILUSetSchurNSHDropThreshold(void *ilu_vdata, JX_Real threshold);
JX_Int jx_ILUSetSchurNSHDropThresholdArray(void *ilu_vdata, JX_Real *threshold);
JX_Int jx_ILUGetNumIterations(void *ilu_vdata, JX_Int *num_iterations);
JX_Int jx_ILUGetFinalRelativeResidualNorm(void *ilu_vdata, JX_Real *res_norm);
JX_Int jx_ILUWriteSolverParams(void *ilu_vdata);
JX_Int jx_ILUMinHeapAddI(JX_Int *heap, JX_Int len);
JX_Int jx_ILUMinHeapAddIIIi(JX_Int *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMinHeapAddIRIi(JX_Int *heap, JX_Real *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMaxHeapAddRabsIIi(JX_Real *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMaxrHeapAddRabsI(JX_Real *heap, JX_Int *I1, JX_Int len);
JX_Int jx_ILUMinHeapRemoveI(JX_Int *heap, JX_Int len);
JX_Int jx_ILUMinHeapRemoveIIIi(JX_Int *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMinHeapRemoveIRIi(JX_Int *heap, JX_Real *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMaxHeapRemoveRabsIIi(JX_Real *heap, JX_Int *I1, JX_Int *Ii1, JX_Int len);
JX_Int jx_ILUMaxrHeapRemoveRabsI(JX_Real *heap, JX_Int *I1, JX_Int len);
JX_Int jx_ILUMaxQSplitRabsI(JX_Real *array, JX_Int *I1, JX_Int left, JX_Int bound, JX_Int right);
JX_Int jx_ILUMaxRabs(JX_Real *array_data, JX_Int *array_j, JX_Int start, JX_Int end, JX_Int nLU, JX_Int *rperm, JX_Real *value, JX_Int *index, JX_Real *l1_norm, JX_Int *nnz);
JX_Int jx_ILUGetPermddPQPre(JX_Int n, JX_Int nLU, JX_Int *A_diag_i, JX_Int *A_diag_j, JX_Real *A_diag_data, JX_Real tol, JX_Int *perm, JX_Int *rperm, JX_Int *pperm_pre, JX_Int *qperm_pre, JX_Int *nB);
JX_Int jx_ILUGetPermddPQ(jx_ParCSRMatrix *A, JX_Int **pperm, JX_Int **qperm, JX_Real tol, JX_Int *nB, JX_Int *nI, JX_Int reordering_type);
JX_Int jx_ILUGetInteriorExteriorPerm(jx_ParCSRMatrix *A, JX_Int **perm, JX_Int *nLU, JX_Int reordering_type);
JX_Int jx_ILUGetLocalPerm(jx_ParCSRMatrix *A, JX_Int **perm, JX_Int *nLU, JX_Int reordering_type);
JX_Int jx_ILUBuildRASExternalMatrix(jx_ParCSRMatrix *A, JX_Int *rperm, JX_Int **E_i, JX_Int **E_j, JX_Real **E_data);
JX_Int jx_ILUSortOffdColmap(jx_ParCSRMatrix *A);
JX_Int jx_ILULocalRCM(jx_CSRMatrix *A, JX_Int start, JX_Int end, JX_Int **permp, JX_Int **qpermp, JX_Int sym);
JX_Int jx_ILULocalRCMMindegree(JX_Int n, JX_Int *degree, JX_Int *marker, JX_Int *rootp);
JX_Int jx_ILULocalRCMOrder(jx_CSRMatrix *A, JX_Int *perm);
JX_Int jx_ILULocalRCMFindPPNode(jx_CSRMatrix *A, JX_Int *rootp, JX_Int *marker);
JX_Int jx_ILULocalRCMBuildLevel(jx_CSRMatrix *A, JX_Int root, JX_Int *marker, JX_Int *level_i, JX_Int *level_j, JX_Int *nlevp);
JX_Int jx_ILULocalRCMNumbering(jx_CSRMatrix *A, JX_Int root, JX_Int *marker, JX_Int *perm, JX_Int *current_nump);
JX_Int jx_ILULocalRCMQsort(JX_Int *perm, JX_Int start, JX_Int end, JX_Int *degree);
JX_Int jx_ILULocalRCMReverse(JX_Int *perm, JX_Int start, JX_Int end);
// Newton-Schultz-Hotelling (NSH) functions
void *jx_NSHCreate();
JX_Int jx_NSHDestroy(void *data);
JX_Int jx_NSHWriteSolverParams(void *nsh_vdata);
JX_Int jx_NSHSetPrintLevel(void *nsh_vdata, JX_Int print_level);
JX_Int jx_NSHSetLogging(void *nsh_vdata, JX_Int logging);
JX_Int jx_NSHSetMaxIter(void *nsh_vdata, JX_Int max_iter);
JX_Int jx_NSHSetTol(void *nsh_vdata, JX_Real tol);
JX_Int jx_NSHSetGlobalSolver(void *nsh_vdata, JX_Int global_solver);
JX_Int jx_NSHSetDropThreshold(void *nsh_vdata, JX_Real droptol);
JX_Int jx_NSHSetDropThresholdArray(void *nsh_vdata, JX_Real *droptol);
JX_Int jx_NSHSetOwnDroptolData(void *nsh_vdata, JX_Int own_droptol_data);
JX_Int jx_NSHSetMRMaxIter(void *nsh_vdata, JX_Int mr_max_iter);
JX_Int jx_NSHSetMRTol(void *nsh_vdata, JX_Real mr_tol);
JX_Int jx_NSHSetMRMaxRowNnz(void *nsh_vdata, JX_Int mr_max_row_nnz);
JX_Int jx_NSHSetColVersion(void *nsh_vdata, JX_Int mr_col_version);
JX_Int jx_NSHSetNSHMaxIter(void *nsh_vdata, JX_Int nsh_max_iter);
JX_Int jx_NSHSetNSHTol(void *nsh_vdata, JX_Real nsh_tol);
JX_Int jx_NSHSetNSHMaxRowNnz(void *nsh_vdata, JX_Int nsh_max_row_nnz);
JX_Int jx_CSRMatrixNormFro(jx_CSRMatrix *A, JX_Real *norm_io);
JX_Int jx_CSRMatrixResNormFro(jx_CSRMatrix *A, JX_Real *norm_io);
JX_Int jx_ParCSRMatrixNormFro(jx_ParCSRMatrix *A, JX_Real *norm_io);
JX_Int jx_ParCSRMatrixResNormFro(jx_ParCSRMatrix *A, JX_Real *norm_io);
JX_Int jx_CSRMatrixTrace(jx_CSRMatrix *A, JX_Real *trace_io);
JX_Int jx_CSRMatrixScaleH(jx_CSRMatrix *A, JX_Real scalar);
JX_Int jx_ParCSRMatrixScaleH(jx_ParCSRMatrix *A, JX_Real scalar);
JX_Int jx_CSRMatrixDropInplace(jx_CSRMatrix *A, JX_Real droptol, JX_Int max_row_nnz);
JX_Int jx_ILUCSRMatrixInverseSelfPrecondMRGlobal(jx_CSRMatrix *matA, jx_CSRMatrix **M, JX_Real droptol, JX_Real tol, JX_Real eps_tol, JX_Int max_row_nnz, JX_Int max_iter, JX_Int print_level);
JX_Int jx_ILUParCSRInverseNSH(jx_ParCSRMatrix *A, jx_ParCSRMatrix **M, JX_Real *droptol, JX_Real mr_tol, JX_Real nsh_tol, JX_Real eps_tol, JX_Int mr_max_row_nnz, JX_Int nsh_max_row_nnz, JX_Int mr_max_iter, JX_Int nsh_max_iter, JX_Int mr_col_version, JX_Int print_level);

/* csrc/ilu/par_ilu_setup.c */
JX_Int jx_ILUSetup(void *ilu_vdata, jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u);
JX_Int jx_ILUSetupILU0(jx_ParCSRMatrix *A, JX_Int *perm, JX_Int *qperm, JX_Int nLU, JX_Int nI, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr, jx_ParCSRMatrix **Sptr, JX_Int **u_end);
JX_Int jx_ILUSetupMILU0(jx_ParCSRMatrix *A, JX_Int *permp, JX_Int *qpermp, JX_Int nLU, JX_Int nI, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr, jx_ParCSRMatrix **Sptr, JX_Int **u_end, JX_Int modified);
JX_Int jx_ILUSetupILUKSymbolic(JX_Int n, JX_Int *A_diag_i, JX_Int *A_diag_j, JX_Int lfil, JX_Int *perm, JX_Int *rperm, JX_Int *iw, JX_Int nLU, JX_Int *L_diag_i, JX_Int *U_diag_i, JX_Int *S_diag_i, JX_Int **L_diag_j, JX_Int **U_diag_j, JX_Int **S_diag_j, JX_Int **u_end);
JX_Int jx_ILUSetupILUK(jx_ParCSRMatrix *A, JX_Int lfil, JX_Int *perm, JX_Int *qperm, JX_Int nLU, JX_Int nI, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr, jx_ParCSRMatrix **Sptr, JX_Int **u_end);
JX_Int jx_ILUSetupILUT(jx_ParCSRMatrix *A, JX_Int lfil, JX_Real *tol, JX_Int *perm, JX_Int *qperm, JX_Int nLU, JX_Int nI, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr, jx_ParCSRMatrix **Sptr, JX_Int **u_end);
JX_Int jx_NSHSetup(void *nsh_vdata, jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u);
JX_Int jx_ILUSetupILU0RAS(jx_ParCSRMatrix *A, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr);
JX_Int jx_ILUSetupILUKRASSymbolic(JX_Int n, JX_Int *A_diag_i, JX_Int *A_diag_j, JX_Int *A_offd_i, JX_Int *A_offd_j, JX_Int *E_i, JX_Int *E_j, JX_Int ext, JX_Int lfil, JX_Int *perm, JX_Int *rperm, JX_Int *iw, JX_Int nLU, JX_Int *L_diag_i, JX_Int *U_diag_i, JX_Int **L_diag_j, JX_Int **U_diag_j);
JX_Int jx_ILUSetupILUKRAS(jx_ParCSRMatrix *A, JX_Int lfil, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr);
JX_Int jx_ILUSetupILUTRAS(jx_ParCSRMatrix *A, JX_Int lfil, JX_Real *tol, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix **Lptr, JX_Real **Dptr, jx_ParCSRMatrix **Uptr);

/* csrc/ilu/par_ilu_solve.c */
JX_Int jx_ILUSolve(void *ilu_vdata, jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u);
JX_Int jx_ILUSolveSchurGMRES(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, JX_Int *perm, JX_Int *qperm, JX_Int nLU, jx_ParCSRMatrix *L, JX_Real *D, jx_ParCSRMatrix *U, jx_ParCSRMatrix *S, jx_ParVector *ftemp, jx_ParVector *utemp, JX_Solver schur_solver, JX_Solver schur_precond, jx_ParVector *rhs, jx_ParVector *x, JX_Int *u_end);
JX_Int jx_ILUSolveSchurNSH(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix *L, JX_Real *D, jx_ParCSRMatrix *U, jx_ParCSRMatrix *S, jx_ParVector *ftemp, jx_ParVector *utemp, JX_Solver schur_solver, jx_ParVector *rhs, jx_ParVector *x, JX_Int *u_end);
JX_Int jx_ILUSolveLU(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix *L, JX_Real *D, jx_ParCSRMatrix *U, jx_ParVector *utemp, jx_ParVector *ftemp);
JX_Int jx_ILUSolveLUIter(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, JX_Int *perm, JX_Int nLU, jx_ParCSRMatrix *L, JX_Real *D, jx_ParCSRMatrix *U, jx_ParVector *ftemp, jx_ParVector *utemp, JX_Int lower_jacobi_iters, JX_Int upper_jacobi_iters);
JX_Int jx_ILUSolveLURAS(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, JX_Int *perm, jx_ParCSRMatrix *L, JX_Real *D, jx_ParCSRMatrix *U, jx_ParVector *ftemp, jx_ParVector *utemp, JX_Real *fext, JX_Real *uext);
JX_Int jx_NSHSolve(void *nsh_vdata, jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *us);
JX_Int jx_NSHSolveInverse(jx_ParCSRMatrix *A, jx_ParVector *f, jx_ParVector *u, jx_ParCSRMatrix *M, jx_ParVector *ftemp, jx_ParVector *utemp);

/* hypre 版本 */
JX_Int jx_ILUSetIterativeSetupType(void *ilu_vdata, JX_Int iter_setup_type);
JX_Int JX_ILUSetIterativeSetupType(JX_Solver solver, JX_Int iter_setup_type);
JX_Int jx_ILUSetIterativeSetupOption(void *ilu_vdata, JX_Int iter_setup_option);
JX_Int JX_ILUSetIterativeSetupOption(JX_Solver solver, JX_Int iter_setup_option);
JX_Int jx_ILUSetIterativeSetupMaxIter(void *ilu_vdata, JX_Int iter_setup_max_iter);
JX_Int JX_ILUSetIterativeSetupMaxIter(JX_Solver solver, JX_Int iter_setup_max_iter);
JX_Int jx_ILUSetIterativeSetupTolerance(void *ilu_vdata, JX_Real iter_setup_tolerance);
JX_Int JX_ILUSetIterativeSetupTolerance(JX_Solver solver, JX_Real iter_setup_tolerance);
JX_Int jx_ILUSetLowerJacobiIters(void *ilu_vdata, JX_Int lower_jacobi_iters);
JX_Int JX_ILUSetLowerJacobiIters(JX_Solver solver, JX_Int lower_jacobi_iterations);
JX_Int jx_ILUSetTriSolve(void *ilu_vdata, JX_Int tri_solve);
JX_Int JX_ILUSetTriSolve(JX_Solver solver, JX_Int tri_solve);
JX_Int jx_ILUSetUpperJacobiIters(void *ilu_vdata, JX_Int upper_jacobi_iters);
JX_Int JX_ILUSetUpperJacobiIters(JX_Solver solver, JX_Int upper_jacobi_iterations);
// JX_Int jx_ILUSetIRIters(void *ilu_vdata, JX_Int IR_iters);
// JX_Int JX_ILUSetIRIters(JX_Solver solver, JX_Int IR_iters);

#endif
