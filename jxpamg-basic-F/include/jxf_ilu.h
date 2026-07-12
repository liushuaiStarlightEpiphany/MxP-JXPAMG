//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_ilu.h -- head files for ILU Preconditioner
 *  Date: 2014/03/24
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JXF_ILU_HEADER
#define JXF_ILU_HEADER
 
#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

#ifndef JXF_HPCSRMV_HEADER 
#include "jxf_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_GridPartitionData
 */
typedef struct
{
    MPI_Comm comm;
    
    JXF_Int part_type;
    JXF_Int x_part_len;
    JXF_Int y_part_len;
    JXF_Int x_lower_idx;
    JXF_Int x_upper_idx;
    JXF_Int y_lower_idx;
    JXF_Int y_upper_idx;
    JXF_Int num_sideproc;
    JXF_Int num_smallside;
    JXF_Int num_largeside;
    JXF_Int num_smallcross;
    JXF_Int num_largecross;
    JXF_Int num_nocrossside;
    
    JXF_Int *xlo_array;
    JXF_Int *xup_array;
    JXF_Int *ylo_array;
    JXF_Int *yup_array;
    JXF_Int *sideprocs;
    JXF_Int *sideprcpos;
    JXF_Int *sideprcxsrt;
    JXF_Int *sideprcysrt;
    JXF_Int *sideprclength;
    
} jxf_GridPartitionData;

/*!
 * \struct jxf_ILUZeroFactorData
 */
typedef struct
{
    MPI_Comm comm;
    
    JXF_Int nx;
    JXF_Int ny;
    JXF_Int npx;
    JXF_Int npy;
    JXF_Int num_equns;
    
    JXF_Int ex_len;
    JXF_Int ey_len;
    JXF_Int dx_len;
    JXF_Int kx_len;
    JXF_Int dy_len;
    JXF_Int lx_len;
    JXF_Int ly_len;
    JXF_Int length;
    JXF_Int max_iter;
    JXF_Int pos_int_end;
    JXF_Int pos_dwn_end;
    JXF_Int pos_lft_end;
    JXF_Int num_fill_in_drop;
    
    JXF_Real drop_tol;
    
    JXF_Int *index;
    JXF_Int *indexA;
    JXF_Int *indexD;
    JXF_Int *permute;
    
    JXF_Real *value;
    JXF_Real *senddown;
    
    MPI_Status *status;
    
    jxf_CSRMatrix *matA;
    
    jxf_Vector *aux_vec;
    jxf_Vector *res_vec;
    jxf_Vector *tmp_vec;
    
    jxf_ParVector *par_aux_vec;
    jxf_ParVector *par_res_vec;
    
    jxf_GridPartitionData *par_grid;
    
} jxf_ILUZeroFactorData;


typedef struct jxf_ParILUData_struct
{
   //general data
   JXF_Int            global_solver;
   jxf_ParCSRMatrix   *matA;
   jxf_ParCSRMatrix   *matL;
   JXF_Real           *matD;
   jxf_ParCSRMatrix   *matU;
   jxf_ParCSRMatrix   *matS;
   JXF_Real           *droptol;/* should be an array of 3 element, for B, (E and F), S respectively */
   JXF_Int            own_droptol_data;/* should I free droptols */
   JXF_Int            lfil;
   JXF_Int            maxRowNnz;
   JXF_Int            *CF_marker_array;
   JXF_Int            *perm;
   JXF_Int            *qperm;
   JXF_Real           tol_ddPQ;
   jxf_ParVector      *F;
   jxf_ParVector      *U;
   jxf_ParVector      *residual;
   JXF_Real           *rel_res_norms;  
   JXF_Int            num_iterations;
   JXF_Real           *l1_norms;
   JXF_Real           final_rel_residual_norm;
   JXF_Real           tol;
   JXF_Real           operator_complexity;
   
   JXF_Int            logging;
   JXF_Int            print_level;
   JXF_Int            max_iter;
   
   JXF_Int            ilu_type;
   JXF_Int            nLU;
   JXF_Int            nI;
   
   /* used when schur block is formed */
   JXF_Int            *u_end;
 
   /* temp vectors for solve phase */
   jxf_ParVector      *Utemp;
   jxf_ParVector      *Ftemp;
   JXF_Real           *uext;
   JXF_Real           *fext;
   
   /* data structure sor solving Schur System */
   JXF_Solver         schur_solver;
   JXF_Solver         schur_precond;
   jxf_ParVector      *rhs;
   jxf_ParVector      *x;
   
   /* schur solver data */
   JXF_Int            ss_kDim;/* dim and max number of iterations for GMRES or max number of iterations for NSH inverse */
   JXF_Int            ss_max_iter;/* max number of iterations for NSH solve */
   JXF_Real           ss_tol;/* stop iteration tol for GMRES or NSH inverse */
   JXF_Real           ss_absolute_tol;/* absolute tol for GMRES or tol for NSH solve */
   JXF_Int            ss_logging;
   JXF_Int            ss_print_level;
   JXF_Int            ss_rel_change;
   
   /* schur precond data */
   JXF_Int            sp_ilu_type;/* ilu type is use ILU, or max rows of nonzeros for NSH */
   JXF_Int            sp_ilu_lfil;/* level of fill in for ILUK or MR column version setting*/
   JXF_Int            sp_ilu_max_row_nnz;/* max rows for ILUT or MR  */
   /* droptol for ILUT or MR 
    * ILUT: [0], [1], [2] B, E&F, S respectively
    * NSH: [0] for MR, [1] for NSH
    */
   JXF_Real           *sp_ilu_droptol;/* droptol array for ILUT or NSH */
   JXF_Int            sp_own_droptol_data;
   JXF_Int            sp_print_level;
   JXF_Int            sp_max_iter;/* max precond iter or max MR iteration */
   JXF_Real           sp_tol;
   
   /* local reordering */
   JXF_Int 	reordering_type;
   
} jxf_ParILUData;

#define jxf_ParILUDataGlobalSolver(ilu_data)                 ((ilu_data) -> global_solver)
#define jxf_ParILUDataMatA(ilu_data)                         ((ilu_data) -> matA)
#define jxf_ParILUDataMatL(ilu_data)                         ((ilu_data) -> matL)
#define jxf_ParILUDataMatD(ilu_data)                         ((ilu_data) -> matD)
#define jxf_ParILUDataMatU(ilu_data)                         ((ilu_data) -> matU)
#define jxf_ParILUDataMatS(ilu_data)                         ((ilu_data) -> matS)
#define jxf_ParILUDataDroptol(ilu_data)                      ((ilu_data) -> droptol)
#define jxf_ParILUDataOwnDroptolData(ilu_data)               ((ilu_data) -> own_droptol_data)
#define jxf_ParILUDataLfil(ilu_data)                         ((ilu_data) -> lfil)
#define jxf_ParILUDataMaxRowNnz(ilu_data)                    ((ilu_data) -> maxRowNnz)
#define jxf_ParILUDataCFMarkerArray(ilu_data)                ((ilu_data) -> CF_marker_array)
#define jxf_ParILUDataPerm(ilu_data)                         ((ilu_data) -> perm)
#define jxf_ParILUDataPPerm(ilu_data)                        ((ilu_data) -> perm)
#define jxf_ParILUDataQPerm(ilu_data)                        ((ilu_data) -> qperm)
#define jxf_ParILUDataTolDDPQ(ilu_data)                      ((ilu_data) -> tol_ddPQ)
#define jxf_ParILUDataF(ilu_data)                            ((ilu_data) -> F)
#define jxf_ParILUDataU(ilu_data)                            ((ilu_data) -> U)
#define jxf_ParILUDataResidual(ilu_data)                     ((ilu_data) -> residual)
#define jxf_ParILUDataRelResNorms(ilu_data)                  ((ilu_data) -> rel_res_norms)
#define jxf_ParILUDataNumIterations(ilu_data)                ((ilu_data) -> num_iterations)
#define jxf_ParILUDataL1Norms(ilu_data)                      ((ilu_data) -> l1_norms)
#define jxf_ParILUDataFinalRelResidualNorm(ilu_data)         ((ilu_data) -> final_rel_residual_norm)
#define jxf_ParILUDataTol(ilu_data)                          ((ilu_data) -> tol)
#define jxf_ParILUDataOperatorComplexity(ilu_data)           ((ilu_data) -> operator_complexity)
#define jxf_ParILUDataLogging(ilu_data)                      ((ilu_data) -> logging)
#define jxf_ParILUDataPrintLevel(ilu_data)                   ((ilu_data) -> print_level)
#define jxf_ParILUDataMaxIter(ilu_data)                      ((ilu_data) -> max_iter)
#define jxf_ParILUDataIluType(ilu_data)                      ((ilu_data) -> ilu_type)
#define jxf_ParILUDataNLU(ilu_data)                          ((ilu_data) -> nLU)
#define jxf_ParILUDataNI(ilu_data)                           ((ilu_data) -> nI)
#define jxf_ParILUDataUEnd(ilu_data)                         ((ilu_data) -> u_end)
#define jxf_ParILUDataUTemp(ilu_data)                        ((ilu_data) -> Utemp)
#define jxf_ParILUDataFTemp(ilu_data)                        ((ilu_data) -> Ftemp)
#define jxf_ParILUDataUExt(ilu_data)                         ((ilu_data) -> uext)
#define jxf_ParILUDataFExt(ilu_data)                         ((ilu_data) -> fext)
#define jxf_ParILUDataSchurSolver(ilu_data)                  ((ilu_data) -> schur_solver)
#define jxf_ParILUDataSchurPrecond(ilu_data)                 ((ilu_data) -> schur_precond)
#define jxf_ParILUDataRhs(ilu_data)                          ((ilu_data) -> rhs)
#define jxf_ParILUDataX(ilu_data)                            ((ilu_data) -> x)
#define jxf_ParILUDataReorderingType(ilu_data)                            ((ilu_data) -> reordering_type)
/* Schur System */
#define jxf_ParILUDataSchurGMRESKDim(ilu_data)               ((ilu_data) -> ss_kDim)
#define jxf_ParILUDataSchurNSHMaxNumIter(ilu_data)           ((ilu_data) -> ss_kDim)
#define jxf_ParILUDataSchurGMRESMaxIter(ilu_data)            ((ilu_data) -> ss_kDim)
#define jxf_ParILUDataSchurNSHSolveMaxIter(ilu_data)         ((ilu_data) -> ss_max_iter)
#define jxf_ParILUDataSchurGMRESTol(ilu_data)                ((ilu_data) -> ss_tol)
#define jxf_ParILUDataSchurNSHTol(ilu_data)                  ((ilu_data) -> ss_tol)
#define jxf_ParILUDataSchurGMRESAbsoluteTol(ilu_data)        ((ilu_data) -> ss_absolute_tol)
#define jxf_ParILUDataSchurNSHSolveTol(ilu_data)             ((ilu_data) -> ss_absolute_tol)
#define jxf_ParILUDataSchurSolverLogging(ilu_data)           ((ilu_data) -> ss_logging)
#define jxf_ParILUDataSchurSolverPrintLevel(ilu_data)        ((ilu_data) -> ss_print_level)
#define jxf_ParILUDataSchurGMRESRelChange(ilu_data)          ((ilu_data) -> ss_rel_change)
#define jxf_ParILUDataSchurPrecondIluType(ilu_data)          ((ilu_data) -> sp_ilu_type)
#define jxf_ParILUDataSchurNSHMaxRowNnz(ilu_data)            ((ilu_data) -> sp_ilu_type)
#define jxf_ParILUDataSchurPrecondIluLfil(ilu_data)          ((ilu_data) -> sp_ilu_lfil)
#define jxf_ParILUDataSchurMRColVersion(ilu_data)            ((ilu_data) -> sp_ilu_lfil)
#define jxf_ParILUDataSchurPrecondIluMaxRowNnz(ilu_data)     ((ilu_data) -> sp_ilu_max_row_nnz)
#define jxf_ParILUDataSchurMRMaxRowNnz(ilu_data)             ((ilu_data) -> sp_ilu_max_row_nnz)
#define jxf_ParILUDataSchurPrecondIluDroptol(ilu_data)       ((ilu_data) -> sp_ilu_droptol)
#define jxf_ParILUDataSchurNSHDroptol(ilu_data)              ((ilu_data) -> sp_ilu_droptol)
#define jxf_ParILUDataSchurPrecondOwnDroptolData(ilu_data)   ((ilu_data) -> sp_own_droptol_data)
#define jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data)       ((ilu_data) -> sp_own_droptol_data)
#define jxf_ParILUDataSchurPrecondPrintLevel(ilu_data)       ((ilu_data) -> sp_print_level)
#define jxf_ParILUDataSchurPrecondMaxIter(ilu_data)          ((ilu_data) -> sp_max_iter)
#define jxf_ParILUDataSchurMRMaxIter(ilu_data)               ((ilu_data) -> sp_max_iter)
#define jxf_ParILUDataSchurPrecondTol(ilu_data)              ((ilu_data) -> sp_tol)
#define jxf_ParILUDataSchurMRTol(ilu_data)                   ((ilu_data) -> sp_tol)


#define FMRK  -1
#define CMRK  1
#define UMRK  0
#define S_CMRK  2

#define FPT(i, bsize) (((i) % (bsize)) == FMRK)
#define CPT(i, bsize) (((i) % (bsize)) == CMRK)

#define MAT_TOL 1e-14
#define EXPAND_FACT 1.3


typedef struct jxf_ParNSHData_struct
{
   /* solver information */
   JXF_Int             global_solver;
   jxf_ParCSRMatrix    *matA;
   jxf_ParCSRMatrix    *matM;
   jxf_ParVector       *F;
   jxf_ParVector       *U;
   jxf_ParVector       *residual;
   JXF_Real            *rel_res_norms;  
   JXF_Int             num_iterations;
   JXF_Real            *l1_norms;
   JXF_Real            final_rel_residual_norm;
   JXF_Real            tol;
   JXF_Real            operator_complexity;
   
   JXF_Int             logging;
   JXF_Int             print_level;
   JXF_Int             max_iter;
   
   /* common data slots */
   /* droptol[0]: droptol for MR 
    * droptol[1]: droptol for NSH
    */
   JXF_Real            *droptol;
   JXF_Int             own_droptol_data;
   
   /* temp vectors for solve phase */
   jxf_ParVector       *Utemp;
   jxf_ParVector       *Ftemp;
   
   /* data slots for local MR */
   JXF_Int             mr_max_iter;
   JXF_Real            mr_tol;
   JXF_Int             mr_max_row_nnz;
   JXF_Int             mr_col_version;/* global version or column version MR */
   
   /* data slots for global NSH */
   JXF_Int             nsh_max_iter;
   JXF_Real            nsh_tol;
   JXF_Int             nsh_max_row_nnz;

}jxf_ParNSHData;

#define jxf_ParNSHDataGlobalSolver(nsh_data)           ((nsh_data) -> global_solver)
#define jxf_ParNSHDataMatA(nsh_data)                   ((nsh_data) -> matA)
#define jxf_ParNSHDataMatM(nsh_data)                   ((nsh_data) -> matM)
#define jxf_ParNSHDataF(nsh_data)                      ((nsh_data) -> F)
#define jxf_ParNSHDataU(nsh_data)                      ((nsh_data) -> U)
#define jxf_ParNSHDataResidual(nsh_data)               ((nsh_data) -> residual)
#define jxf_ParNSHDataRelResNorms(nsh_data)            ((nsh_data) -> rel_res_norms)
#define jxf_ParNSHDataNumIterations(nsh_data)          ((nsh_data) -> num_iterations)
#define jxf_ParNSHDataL1Norms(nsh_data)                ((nsh_data) -> l1_norms)
#define jxf_ParNSHDataFinalRelResidualNorm(nsh_data)   ((nsh_data) -> final_rel_residual_norm)
#define jxf_ParNSHDataTol(nsh_data)                    ((nsh_data) -> tol)
#define jxf_ParNSHDataOperatorComplexity(nsh_data)     ((nsh_data) -> operator_complexity)
#define jxf_ParNSHDataLogging(nsh_data)                ((nsh_data) -> logging)
#define jxf_ParNSHDataPrintLevel(nsh_data)             ((nsh_data) -> print_level)
#define jxf_ParNSHDataMaxIter(nsh_data)                ((nsh_data) -> max_iter)
#define jxf_ParNSHDataDroptol(nsh_data)                ((nsh_data) -> droptol)
#define jxf_ParNSHDataOwnDroptolData(nsh_data)         ((nsh_data) -> own_droptol_data)
#define jxf_ParNSHDataUTemp(nsh_data)                  ((nsh_data) -> Utemp)
#define jxf_ParNSHDataFTemp(nsh_data)                  ((nsh_data) -> Ftemp)
#define jxf_ParNSHDataMRMaxIter(nsh_data)              ((nsh_data) -> mr_max_iter)
#define jxf_ParNSHDataMRTol(nsh_data)                  ((nsh_data) -> mr_tol)
#define jxf_ParNSHDataMRMaxRowNnz(nsh_data)            ((nsh_data) -> mr_max_row_nnz)
#define jxf_ParNSHDataMRColVersion(nsh_data)           ((nsh_data) -> mr_col_version)
#define jxf_ParNSHDataNSHMaxIter(nsh_data)             ((nsh_data) -> nsh_max_iter)
#define jxf_ParNSHDataNSHTol(nsh_data)                 ((nsh_data) -> nsh_tol)
#define jxf_ParNSHDataNSHMaxRowNnz(nsh_data)           ((nsh_data) -> nsh_max_row_nnz)

#define DIVIDE_TOL 1e-32

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/ilu/csr_matrix.c */
jxf_CSRMatrix *jxf_CSRMatrixTDMGReorderByNodes( jxf_CSRMatrix *A, JXF_Int num_equns );
jxf_CSRMatrix *jxf_CSRMatrixTDMGReorderByVariables( jxf_CSRMatrix *A, JXF_Int num_equns );
jxf_CSRMatrix *jxf_CSRMatrixTDMGReorderAlongX( jxf_CSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );
jxf_CSRMatrix *jxf_CSRMatrixTDMGReorderAlongY( jxf_CSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );

/* csrc/ilu/decomposition.c */
JXF_Int jxf_ILUZeroDecompositionA( jxf_CSRMatrix *A,
                              JXF_Int **indexDP_ptr,
                              JXF_Int **indexLU_ptr,
                              JXF_Real **valueLU_ptr );
JXF_Int jxf_ILUZeroDecompositionB( jxf_CSRMatrix *A,
                              JXF_Real drop_tol,
                              JXF_Int **indexAP_ptr,
                              JXF_Int **indexDP_ptr,
                              JXF_Int **indexLU_ptr,
                              JXF_Real **valueLU_ptr );
JXF_Int *jxf_ILUZeroParallelDecompositionA( jxf_ParCSRMatrix *par_A,
                                       JXF_Int **indexDP_ptr,
                                       JXF_Int **indexLU_ptr,
                                       JXF_Real **valueLU_ptr,
                                       JXF_Int *num_nonzeros,
                                       JXF_Int *fill_in_drop,
                                       jxf_ILUZeroFactorData *ilu_data );
JXF_Int
jxf_ILUZeroLocalDecompositionIntURPntsA( JXF_Int *IA,
                                        JXF_Int *JA,
                                        JXF_Real *AA,
                                        JXF_Int *indexLU,
                                        JXF_Int *indexDP,
                                        JXF_Int *placeRC,
                                        JXF_Real *valueLU,
                                        JXF_Int int_uprgt_pnt,
                                        JXF_Int first_row_idx );
void
jxf_ILUZeroLocalDecompositionDPntsA( JXF_Int *IA,
                                    JXF_Int *JA,
                                    JXF_Real *AA,
                                    JXF_Int *indexLU,
                                    JXF_Int *indexDP,
                                    JXF_Int *placeRC,
                                    JXF_Real *valueLU,
                                    JXF_Int int_uprgt_pnt,
                                    JXF_Int first_row_idx,
                                    JXF_Int num_rows,
                                    JXF_Int ex_len,
                                    JXF_Int *fill_in_drop );
void
jxf_ILUZeroLocalDecompositionLPntsA( JXF_Int *IA,
                                    JXF_Int *JA,
                                    JXF_Real *AA,
                                    JXF_Int *indexLU,
                                    JXF_Int *indexDP,
                                    JXF_Int *placeRC,
                                    JXF_Real *valueLU,
                                    JXF_Int int_uprgtdwn_pnt,
                                    JXF_Int first_row_idx,
                                    JXF_Int num_rows,
                                    JXF_Int dwn_num_rows,
                                    JXF_Int ex_len,
                                    JXF_Int *fill_in_drop );
void
jxf_ILUZeroLocalDecompositionLDPntsA( JXF_Int *IA,
                                     JXF_Int *JA,
                                     JXF_Real *AA,
                                     JXF_Int *indexLU,
                                     JXF_Int *indexDP,
                                     JXF_Int *placeRC,
                                     JXF_Real *valueLU,
                                     JXF_Int int_urdl_pnt,
                                     JXF_Int first_row_idx,
                                     JXF_Int fst_row_idx,
                                     JXF_Int num_rows,
                                     JXF_Int lft_cnum_rows,
                                     JXF_Int dwn_cnum_rows,
                                     JXF_Int ex_len,
                                     JXF_Int *fill_in_drop );
void jxf_ILUZeroFactorDataParallelUPartIntURPntsA( JXF_Int *indexLU,
                                                  JXF_Int *indexAP,
                                                  JXF_Int *indexDP,
                                                  JXF_Real *valueLU,
                                                  JXF_Int int_uprgt_pnt,
                                                  JXF_Int recv_downsrt,
                                                  JXF_Int recv_leftsrt,
                                                  JXF_Int *permute,
                                                  JXF_Real *aux_data,
                                                  JXF_Real *res_data,
                                                  JXF_Real *app_data,
                                                  JXF_Int first_row_idx,
                                                  JXF_Int next_row_idx,
                                                  JXF_Int nnxt_row_idx );
void jxf_ILUZeroFactorDataParallelUPartIntURPntsB( JXF_Int *indexLU,
                                                  JXF_Int *indexAP,
                                                  JXF_Int *indexDP,
                                                  JXF_Real *valueLU,
                                                  JXF_Int int_uprgt_pnt,
                                                  JXF_Int recv_downsrt,
                                                  JXF_Int *permute,
                                                  JXF_Real *aux_data,
                                                  JXF_Real *res_data,
                                                  JXF_Real *app_data,
                                                  JXF_Int first_row_idx,
                                                  JXF_Int next_row_idx );
void jxf_ILUZeroFactorDataParallelUPartIntURPntsC( JXF_Int *indexLU,
                                                  JXF_Int *indexAP,
                                                  JXF_Int *indexDP,
                                                  JXF_Real *valueLU,
                                                  JXF_Int int_uprgt_pnt,
                                                  JXF_Int recv_leftsrt,
                                                  JXF_Int *permute,
                                                  JXF_Real *aux_data,
                                                  JXF_Real *res_data,
                                                  JXF_Real *app_data,
                                                  JXF_Int first_row_idx,
                                                  JXF_Int next_row_idx );
void jxf_ILUZeroFactorDataParallelUPartIntURPntsD( JXF_Int *indexLU,
                                                  JXF_Int *indexAP,
                                                  JXF_Int *indexDP,
                                                  JXF_Real *valueLU,
                                                  JXF_Int int_uprgt_pnt,
                                                  JXF_Int *permute,
                                                  JXF_Real *aux_data,
                                                  JXF_Real *res_data,
                                                  JXF_Real *app_data,
                                                  JXF_Int first_row_idx );

/* csrc/ilu/grid.c */
JXF_Int JXF_GridPartitionDataCreate( JXF_Solver *solver, MPI_Comm comm, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy );
JXF_Int JXF_GridPartitionDataSetEachSides4Comm( JXF_Solver solver );
JXF_Int JXF_GridPartitionDataDestroy( JXF_Solver solver );
void *jxf_GridPartitionDataInitialize( MPI_Comm comm, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy );
JXF_Int jxf_GridPartitionDataSetEachSides4Comm( void *grid_vdata );
JXF_Int jxf_GridPartitionDataFinalize( void *grid_vdata );

/* csrc/ilu/ilu.c */
JXF_Int JXF_ILUZeroFactorDataCreate( JXF_Solver *solver, MPI_Comm comm );
JXF_Int JXF_ILUZeroFactorDataSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_ILUZeroFactorDataSetNxy( JXF_Solver solver, JXF_Int nx, JXF_Int ny );
JXF_Int JXF_ILUZeroFactorDataSetNpxy( JXF_Solver solver, JXF_Int npx, JXF_Int npy );
JXF_Int JXF_ILUZeroFactorDataSetNumEquns( JXF_Solver solver, JXF_Int num_equns );
JXF_Int JXF_ILUZeroFactorDataSetDropTol( JXF_Solver solver, JXF_Real drop_tol );
JXF_Int JXF_ILUZeroFactorDataSetMatA( JXF_Solver solver, jxf_CSRMatrix *matA );
JXF_Int JXF_ILUZeroFactorDataGetLULength( JXF_Solver solver, JXF_Int *lu_length );
JXF_Int JXF_ILUZeroFactorDataDestroy( JXF_Solver solver );
JXF_Int JXF_ILUZeroFactorDataGenerateParGrid( JXF_Solver solver );
JXF_Int JXF_ILUZeroFactorDataSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int JXF_ILUZeroFactorDataPrecond( JXF_Solver       solver,
                                 JXF_hpCSRMatrix par_matrix,
                                 JXF_ParVector    par_rhs,
                                 JXF_ParVector    par_app  );
void *jxf_ILUZeroFactorDataInitialize( MPI_Comm comm );
JXF_Int jxf_ILUZeroFactorDataSetMaxIter( void *ilu_vdata, JXF_Int max_iter );
JXF_Int jxf_ILUZeroFactorDataSetNxy( void *ilu_vdata, JXF_Int nx, JXF_Int ny );
JXF_Int jxf_ILUZeroFactorDataSetNpxy( void *ilu_vdata, JXF_Int npx, JXF_Int npy );
JXF_Int jxf_ILUZeroFactorDataSetNumEquns( void *ilu_vdata, JXF_Int num_equns );
JXF_Int jxf_ILUZeroFactorDataSetDropTol( void *ilu_vdata, JXF_Real drop_tol );
JXF_Int jxf_ILUZeroFactorDataSetMatA( void *ilu_vdata, jxf_CSRMatrix *matA );
JXF_Int jxf_ILUZeroFactorDataGetLULength( void *ilu_vdata, JXF_Int *lu_length );
JXF_Int jxf_ILUZeroFactorDataFinalize( void *ilu_vdata );
JXF_Int jxf_ILUZeroFactorDataGenerateParGrid( void *ilu_vdata );

/* csrc/ilu/ilucycle.c */
JXF_Int jxf_ILUZeroFactorDataPrecond( void            *ilu_vdata,
                                 jxf_ParCSRMatrix *par_A,
                                 jxf_ParVector    *par_b,
                                 jxf_ParVector    *par_x );
void jxf_ILUZeroFactorDataCycleA( void *ilu_vdata, jxf_CSRMatrix *A, jxf_Vector *f, jxf_Vector *u );
void jxf_ILUZeroFactorDataCycleB( void *ilu_vdata, jxf_CSRMatrix *A, jxf_Vector *f, jxf_Vector *u );
void jxf_ILUZeroFactorDataParallelCycleA( jxf_ILUZeroFactorData *ilu_data,
                                         jxf_ParCSRMatrix *par_A,
                                         jxf_ParVector *par_b,
                                         jxf_ParVector *par_x );
void jxf_ILUZeroFactorDataParallelLPartIntURPntsA( JXF_Int *indexLU,
                                                  JXF_Int *indexAP,
                                                  JXF_Int *indexDP,
                                                  JXF_Real *valueLU,
                                                  JXF_Int *permute,
                                                  JXF_Int int_uprgt_pnt,
                                                  JXF_Real *aux_data,
                                                  JXF_Real *rhs_data,
                                                  JXF_Int first_row_idx );
void jxf_ILUZeroFactorDataParallelLPartDPntsA( JXF_Int *indexLU,
                                              JXF_Int *indexAP,
                                              JXF_Int *indexDP,
                                              JXF_Real *valueLU,
                                              JXF_Int int_uprgt_pnt,
                                              JXF_Int num_rows,
                                              JXF_Int *permute,
                                              JXF_Real *aux_data,
                                              JXF_Real *res_data,
                                              JXF_Real *rhs_data,
                                              JXF_Int first_row_idx );
void jxf_ILUZeroFactorDataParallelLPartLPntsA( JXF_Int *indexLU,
                                              JXF_Int *indexAP,
                                              JXF_Int *indexDP,
                                              JXF_Real *valueLU,
                                              JXF_Int int_uprgtdwn_pnt,
                                              JXF_Int num_rows,
                                              JXF_Int dwn_num_rows,
                                              JXF_Int *permute,
                                              JXF_Real *aux_data,
                                              JXF_Real *res_data,
                                              JXF_Real *rhs_data,
                                              JXF_Int first_row_idx );
void jxf_ILUZeroFactorDataParallelLPartLDPntsA( JXF_Int *indexLU,
                                               JXF_Int *indexAP,
                                               JXF_Int *indexDP,
                                               JXF_Real *valueLU,
                                               JXF_Int int_urdl_pnt,
                                               JXF_Int num_rows,
                                               JXF_Int fst_row_idx,
                                               JXF_Int *permute,
                                               JXF_Int lft_cnum_rows,
                                               JXF_Int dwn_dnum_rows,
                                               JXF_Real *aux_data,
                                               JXF_Real *res_data,
                                               JXF_Real *rhs_data,
                                               JXF_Int first_row_idx );
void jxf_ILUZeroFactorDataParallelUPartLDPntsA( JXF_Int *indexLU,
                                               JXF_Int *indexAP,
                                               JXF_Int *indexDP,
                                               JXF_Real *valueLU,
                                               JXF_Int int_urdl_pnt,
                                               JXF_Int num_rows,
                                               JXF_Int *permute,
                                               JXF_Real *aux_data,
                                               JXF_Real *app_data,
                                               JXF_Int first_row_idx );
void jxf_ILUZeroFactorDataParallelUPartLPntsA( JXF_Int *indexLU,
                                              JXF_Int *indexAP,
                                              JXF_Int *indexDP,
                                              JXF_Real *valueLU,
                                              JXF_Int int_urdl_pnt,
                                              JXF_Int num_rows,
                                              JXF_Int recv_downsrt,
                                              JXF_Int *permute,
                                              JXF_Real *aux_data,
                                              JXF_Real *res_data,
                                              JXF_Real *app_data,
                                              JXF_Int first_row_idx,
                                              JXF_Int next_row_idx );
void jxf_ILUZeroFactorDataParallelUPartDPntsA( JXF_Int *indexLU,
                                              JXF_Int *indexAP,
                                              JXF_Int *indexDP,
                                              JXF_Real *valueLU,
                                              JXF_Int int_uprgt_pnt,
                                              JXF_Int num_rows,
                                              JXF_Int recv_leftsrt,
                                              JXF_Int *permute,
                                              JXF_Real *aux_data,
                                              JXF_Real *res_data,
                                              JXF_Real *app_data,
                                              JXF_Int first_row_idx,
                                              JXF_Int next_row_idx );

/* csrc/ilu/ilusetup.c */
JXF_Int jxf_ILUZeroFactorDataSetup( void *ilu_vdata, jxf_ParCSRMatrix *par_A );

/* csrc/ilu/par_csr_matrix.c */
jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile2( char *filename, 
                                             JXF_Int   file_base,
                                             JXF_Int  *row_part, 
                                             JXF_Int  *col_part,
                                             JXF_Int   num_equns,
                                             JXF_Int   nx,
                                             JXF_Int   ny );
jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile3( char *filename, JXF_Int file_base, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );
jxf_ParCSRMatrix *
jxf_ParCSRMatrixCreate2( MPI_Comm   comm,
                        JXF_Int        global_num_rows,
                        JXF_Int        global_num_cols,
                        JXF_Int       *row_starts,
                        JXF_Int       *col_starts,
                        JXF_Int        num_cols_offd,
                        JXF_Int        num_nonzeros_diag,
                        JXF_Int        num_nonzeros_offd,
                        JXF_Int        num_equns,
                        JXF_Int        nx,
                        JXF_Int        ny );
jxf_ParCSRMatrix *jxf_CSRMatrixToParCSRMatrix2( MPI_Comm comm, jxf_CSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );
jxf_CSRMatrix *jxf_MergeDiagAndOffdDropSmall( jxf_ParCSRMatrix *par_matrix, JXF_Real drop_tol );
jxf_CSRMatrix *
jxf_CSRMatrixMergeReorderIntUpRtDnLtQuasiBdy( MPI_Comm              comm,
                                             JXF_Int                  *row_starts,
                                             JXF_Int                  *permute,
                                             jxf_CSRMatrix         *ser_B,
                                             JXF_Int                   nz_srt,
                                             JXF_Int                   ng_pt,
                                             jxf_GridPartitionData *grid_data,
                                             JXF_Int                  *ex_len,
                                             JXF_Int                  *ey_len,
                                             JXF_Int                  *dx_len,
                                             JXF_Int                  *kx_len,
                                             JXF_Int                  *dy_len,
                                             JXF_Int                  *lx_len,
                                             JXF_Int                  *ly_len,
                                             JXF_Int                  *postn_a,
                                             JXF_Int                  *postn_b,
                                             JXF_Int                  *postn_c );

/* csrc/ilu/par_vector.c */
jxf_ParVector *jxf_BuildRhsParFromOneFile2( char *filename, jxf_ParCSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );
JXF_Int jxf_GeneratePartitioning2(JXF_Int num_procs, JXF_Int **part_ptr, JXF_Int num_equns, JXF_Int length, JXF_Int ny);

/* csrc/ilu/vector.c */
jxf_Vector *jxf_VectorTDMGReorderByNodes( jxf_Vector *x, JXF_Int num_equns );
jxf_Vector *jxf_VectorTDMGReorderByVariables( jxf_Vector *x, JXF_Int num_equns );
jxf_Vector *jxf_VectorTDMGReorderAlongX( jxf_Vector *x, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );
jxf_Vector *jxf_VectorTDMGReorderAlongY( jxf_Vector *x, JXF_Int num_equns, JXF_Int nx, JXF_Int ny );

/* csrc/ilu/par_ilu.c */
JXF_Int JXF_ILUCreate( JXF_Solver *solver );
JXF_Int JXF_ILUDestroy( JXF_Solver solver );
JXF_Int JXF_ILUSetup( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int JXF_ILUSolve( JXF_Solver solver, JXF_ParCSRMatrix A, JXF_ParVector b, JXF_ParVector x );
JXF_Int JXF_ILUSetPrintLevel( JXF_Solver solver, JXF_Int print_level );
JXF_Int JXF_ILUSetLogging( JXF_Solver solver, JXF_Int logging );
JXF_Int JXF_ILUSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_ILUSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_ILUSetDropThreshold( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_ILUSetDropThresholdArray( JXF_Solver solver, JXF_Real *threshold );
JXF_Int JXF_ILUSetNSHDropThreshold( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_ILUSetNSHDropThresholdArray( JXF_Solver solver, JXF_Real *threshold );
JXF_Int JXF_ILUSetSchurMaxIter( JXF_Solver solver, JXF_Int ss_max_iter );
JXF_Int JXF_ILUSetMaxNnzPerRow( JXF_Solver solver, JXF_Int nzmax );
JXF_Int JXF_ILUSetLevelOfFill( JXF_Solver solver, JXF_Int lfil );
JXF_Int JXF_ILUSetType( JXF_Solver solver, JXF_Int ilu_type );
JXF_Int JXF_ILUGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_ILUGetFinalRelativeResidualNorm(  JXF_Solver solver, JXF_Real *res_norm );
JXF_Int JXF_ILUSetLocalReordering(  JXF_Solver solver, JXF_Int ordering_type );
void *jxf_ILUCreate();
JXF_Int jxf_ILUDestroy( void *data );
JXF_Int jxf_ILUSetLevelOfFill( void *ilu_vdata, JXF_Int lfil );
JXF_Int jxf_ILUSetMaxNnzPerRow( void *ilu_vdata, JXF_Int nzmax );
JXF_Int jxf_ILUSetDropThreshold( void *ilu_vdata, JXF_Real threshold );
JXF_Int jxf_ILUSetDropThresholdArray( void *ilu_vdata, JXF_Real *threshold );
JXF_Int jxf_ILUSetOwnDropThreshold( void *ilu_vdata, JXF_Int own_droptol_data );
JXF_Int jxf_ILUSetType( void *ilu_vdata, JXF_Int ilu_type );
JXF_Int jxf_ILUSetMaxIter( void *ilu_vdata, JXF_Int max_iter );
JXF_Int jxf_ILUSetTol( void *ilu_vdata, JXF_Real tol );
JXF_Int jxf_ILUSetPrintLevel( void *ilu_vdata, JXF_Int print_level );
JXF_Int jxf_ILUSetLogging( void *ilu_vdata, JXF_Int logging );
JXF_Int jxf_ILUSetLocalReordering( void *ilu_vdata, JXF_Int ordering_type );
JXF_Int jxf_ILUSetSchurSolverKDIM( void *ilu_vdata, JXF_Int ss_kDim );
JXF_Int jxf_ILUSetSchurSolverMaxIter( void *ilu_vdata, JXF_Int ss_max_iter );
JXF_Int jxf_ILUSetSchurSolverTol( void *ilu_vdata, JXF_Real ss_tol );
JXF_Int jxf_ILUSetSchurSolverAbsoluteTol( void *ilu_vdata, JXF_Real ss_absolute_tol );
JXF_Int jxf_ILUSetSchurSolverLogging( void *ilu_vdata, JXF_Int ss_logging );
JXF_Int jxf_ILUSetSchurSolverPrintLevel( void *ilu_vdata, JXF_Int ss_print_level );
JXF_Int jxf_ILUSetSchurSolverRelChange( void *ilu_vdata, JXF_Int ss_rel_change );
JXF_Int jxf_ILUSetSchurPrecondILUType( void *ilu_vdata, JXF_Int sp_ilu_type );
JXF_Int jxf_ILUSetSchurPrecondILULevelOfFill( void *ilu_vdata, JXF_Int sp_ilu_lfil );
JXF_Int jxf_ILUSetSchurPrecondILUMaxNnzPerRow( void *ilu_vdata, JXF_Int sp_ilu_max_row_nnz );
JXF_Int jxf_ILUSetSchurPrecondILUDropThreshold( void *ilu_vdata, JXF_Real sp_ilu_droptol );
JXF_Int jxf_ILUSetSchurPrecondILUDropThresholdArray( void *ilu_vdata, JXF_Real *sp_ilu_droptol );
JXF_Int jxf_ILUSetSchurPrecondILUOwnDropThreshold( void *ilu_vdata, JXF_Int sp_own_droptol_data );
JXF_Int jxf_ILUSetSchurPrecondPrintLevel( void *ilu_vdata, JXF_Int sp_print_level );
JXF_Int jxf_ILUSetSchurPrecondMaxIter( void *ilu_vdata, JXF_Int sp_max_iter );
JXF_Int jxf_ILUSetSchurPrecondTol( void *ilu_vdata, JXF_Int sp_tol );
JXF_Int jxf_ILUSetSchurNSHDropThreshold( void *ilu_vdata, JXF_Real threshold);
JXF_Int jxf_ILUSetSchurNSHDropThresholdArray( void *ilu_vdata, JXF_Real *threshold);
JXF_Int jxf_ILUGetNumIterations( void *ilu_vdata, JXF_Int *num_iterations );
JXF_Int jxf_ILUGetFinalRelativeResidualNorm( void *ilu_vdata, JXF_Real *res_norm );
JXF_Int jxf_ILUWriteSolverParams(void *ilu_vdata);
JXF_Int jxf_ILUMinHeapAddI(JXF_Int *heap, JXF_Int len);
JXF_Int jxf_ILUMinHeapAddIIIi(JXF_Int *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMinHeapAddIRIi(JXF_Int *heap, JXF_Real *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMaxHeapAddRabsIIi(JXF_Real *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMaxrHeapAddRabsI(JXF_Real *heap, JXF_Int *I1, JXF_Int len);
JXF_Int jxf_ILUMinHeapRemoveI(JXF_Int *heap, JXF_Int len);
JXF_Int jxf_ILUMinHeapRemoveIIIi(JXF_Int *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMinHeapRemoveIRIi(JXF_Int *heap, JXF_Real *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMaxHeapRemoveRabsIIi(JXF_Real *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len);
JXF_Int jxf_ILUMaxrHeapRemoveRabsI(JXF_Real *heap, JXF_Int *I1, JXF_Int len);
JXF_Int jxf_ILUMaxQSplitRabsI(JXF_Real *array, JXF_Int *I1, JXF_Int left, JXF_Int bound, JXF_Int right);
JXF_Int jxf_ILUMaxRabs(JXF_Real *array_data, JXF_Int *array_j, JXF_Int start, JXF_Int end, JXF_Int nLU, JXF_Int *rperm, JXF_Real *value, JXF_Int *index, JXF_Real *l1_norm, JXF_Int *nnz);
JXF_Int jxf_ILUGetPermddPQPre(JXF_Int n, JXF_Int nLU, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Real *A_diag_data, JXF_Real tol, JXF_Int *perm, JXF_Int *rperm, JXF_Int *pperm_pre, JXF_Int *qperm_pre, JXF_Int *nB);
JXF_Int jxf_ILUGetPermddPQ(jxf_ParCSRMatrix *A, JXF_Int **pperm, JXF_Int **qperm, JXF_Real tol, JXF_Int *nB, JXF_Int *nI, JXF_Int reordering_type);
JXF_Int jxf_ILUGetInteriorExteriorPerm(jxf_ParCSRMatrix *A, JXF_Int **perm, JXF_Int *nLU, JXF_Int reordering_type);
JXF_Int jxf_ILUGetLocalPerm(jxf_ParCSRMatrix *A, JXF_Int **perm, JXF_Int *nLU, JXF_Int reordering_type);
JXF_Int jxf_ILUBuildRASExternalMatrix(jxf_ParCSRMatrix *A, JXF_Int *rperm, JXF_Int **E_i, JXF_Int **E_j, JXF_Real **E_data);
JXF_Int jxf_ILUSortOffdColmap(jxf_ParCSRMatrix *A);
JXF_Int jxf_ILULocalRCM( jxf_CSRMatrix *A, JXF_Int start, JXF_Int end, JXF_Int **permp, JXF_Int **qpermp, JXF_Int sym);
JXF_Int jxf_ILULocalRCMMindegree(JXF_Int n, JXF_Int *degree, JXF_Int *marker, JXF_Int *rootp);
JXF_Int jxf_ILULocalRCMOrder( jxf_CSRMatrix *A, JXF_Int *perm);
JXF_Int jxf_ILULocalRCMFindPPNode( jxf_CSRMatrix *A, JXF_Int *rootp, JXF_Int *marker);
JXF_Int jxf_ILULocalRCMBuildLevel(jxf_CSRMatrix *A, JXF_Int root, JXF_Int *marker, JXF_Int *level_i, JXF_Int *level_j, JXF_Int *nlevp);
JXF_Int jxf_ILULocalRCMNumbering(jxf_CSRMatrix *A, JXF_Int root, JXF_Int *marker, JXF_Int *perm, JXF_Int *current_nump);
JXF_Int jxf_ILULocalRCMQsort(JXF_Int *perm, JXF_Int start, JXF_Int end, JXF_Int *degree);
JXF_Int jxf_ILULocalRCMReverse(JXF_Int *perm, JXF_Int start, JXF_Int end);
// Newton-Schultz-Hotelling (NSH) functions
void * jxf_NSHCreate();
JXF_Int jxf_NSHDestroy( void *data );
JXF_Int jxf_NSHWriteSolverParams(void *nsh_vdata);
JXF_Int jxf_NSHSetPrintLevel( void *nsh_vdata, JXF_Int print_level );
JXF_Int jxf_NSHSetLogging( void *nsh_vdata, JXF_Int logging );
JXF_Int jxf_NSHSetMaxIter( void *nsh_vdata, JXF_Int max_iter );
JXF_Int jxf_NSHSetTol( void *nsh_vdata, JXF_Real tol );
JXF_Int jxf_NSHSetGlobalSolver( void *nsh_vdata, JXF_Int global_solver );
JXF_Int jxf_NSHSetDropThreshold( void *nsh_vdata, JXF_Real droptol );
JXF_Int jxf_NSHSetDropThresholdArray( void *nsh_vdata, JXF_Real *droptol );
JXF_Int jxf_NSHSetOwnDroptolData( void *nsh_vdata, JXF_Int own_droptol_data );
JXF_Int jxf_NSHSetMRMaxIter( void *nsh_vdata, JXF_Int mr_max_iter );
JXF_Int jxf_NSHSetMRTol( void *nsh_vdata, JXF_Real mr_tol );
JXF_Int jxf_NSHSetMRMaxRowNnz( void *nsh_vdata, JXF_Int mr_max_row_nnz );
JXF_Int jxf_NSHSetColVersion( void *nsh_vdata, JXF_Int mr_col_version );
JXF_Int jxf_NSHSetNSHMaxIter( void *nsh_vdata, JXF_Int nsh_max_iter );
JXF_Int jxf_NSHSetNSHTol( void *nsh_vdata, JXF_Real nsh_tol );
JXF_Int jxf_NSHSetNSHMaxRowNnz( void *nsh_vdata, JXF_Int nsh_max_row_nnz );
JXF_Int jxf_CSRMatrixNormFro(jxf_CSRMatrix *A, JXF_Real *norm_io);
JXF_Int jxf_CSRMatrixResNormFro(jxf_CSRMatrix *A, JXF_Real *norm_io);
JXF_Int jxf_ParCSRMatrixNormFro(jxf_ParCSRMatrix *A, JXF_Real *norm_io);
JXF_Int jxf_ParCSRMatrixResNormFro(jxf_ParCSRMatrix *A, JXF_Real *norm_io);
JXF_Int jxf_CSRMatrixTrace(jxf_CSRMatrix *A, JXF_Real *trace_io);
JXF_Int jxf_CSRMatrixScaleH(jxf_CSRMatrix *A, JXF_Real scalar);
JXF_Int jxf_ParCSRMatrixScaleH(jxf_ParCSRMatrix *A, JXF_Real scalar);
JXF_Int jxf_CSRMatrixDropInplace(jxf_CSRMatrix *A, JXF_Real droptol, JXF_Int max_row_nnz);
JXF_Int jxf_ILUCSRMatrixInverseSelfPrecondMRGlobal(jxf_CSRMatrix *matA, jxf_CSRMatrix **M, JXF_Real droptol, JXF_Real tol, JXF_Real eps_tol, JXF_Int max_row_nnz, JXF_Int max_iter, JXF_Int print_level );
JXF_Int jxf_ILUParCSRInverseNSH(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix **M, JXF_Real *droptol, JXF_Real mr_tol, JXF_Real nsh_tol, JXF_Real eps_tol, JXF_Int mr_max_row_nnz, JXF_Int nsh_max_row_nnz, JXF_Int mr_max_iter, JXF_Int nsh_max_iter, JXF_Int mr_col_version, JXF_Int print_level);

/* csrc/ilu/par_ilu_setup.c */
JXF_Int jxf_ILUSetup( void *ilu_vdata, jxf_ParCSRMatrix *A );
JXF_Int jxf_ILUSetupILU0(jxf_ParCSRMatrix *A, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end);
JXF_Int jxf_ILUSetupILUKSymbolic(JXF_Int n, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Int lfil, JXF_Int *perm, JXF_Int *rperm, JXF_Int *iw, JXF_Int nLU, JXF_Int *L_diag_i, JXF_Int *U_diag_i, JXF_Int *S_diag_i, JXF_Int **L_diag_j, JXF_Int **U_diag_j, JXF_Int **S_diag_j, JXF_Int **u_end);
JXF_Int jxf_ILUSetupILUK(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end);
JXF_Int jxf_ILUSetupILUT(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Real *tol, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, JXF_Int nI, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr, jxf_ParCSRMatrix **Sptr, JXF_Int **u_end);
JXF_Int jxf_NSHSetup( void *nsh_vdata, jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u );
JXF_Int jxf_ILUSetupILU0RAS(jxf_ParCSRMatrix *A, JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr);
JXF_Int jxf_ILUSetupILUKRASSymbolic(JXF_Int n, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Int *A_offd_i, JXF_Int *A_offd_j, JXF_Int *E_i, JXF_Int *E_j, JXF_Int ext, JXF_Int lfil, JXF_Int *perm, JXF_Int *rperm, JXF_Int *iw, JXF_Int nLU, JXF_Int *L_diag_i, JXF_Int *U_diag_i, JXF_Int **L_diag_j, JXF_Int **U_diag_j);
JXF_Int jxf_ILUSetupILUKRAS(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr);
JXF_Int jxf_ILUSetupILUTRAS(jxf_ParCSRMatrix *A, JXF_Int lfil, JXF_Real *tol, JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix **Lptr, JXF_Real** Dptr, jxf_ParCSRMatrix **Uptr);

/* csrc/ilu/par_ilu_solve.c */
JXF_Int jxf_ILUSolve( void *ilu_vdata, jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u );
JXF_Int jxf_ILUSolveSchurGMRES(jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u, JXF_Int *perm, JXF_Int *qperm, JXF_Int nLU, jxf_ParCSRMatrix *L, JXF_Real* D, jxf_ParCSRMatrix *U, jxf_ParCSRMatrix *S, jxf_ParVector *ftemp, jxf_ParVector *utemp, JXF_Solver schur_solver, JXF_Solver schur_precond, jxf_ParVector *rhs, jxf_ParVector *x, JXF_Int *u_end);
JXF_Int jxf_ILUSolveSchurNSH(jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u, JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix *L, JXF_Real* D, jxf_ParCSRMatrix *U, jxf_ParCSRMatrix *S, jxf_ParVector *ftemp, jxf_ParVector *utemp, JXF_Solver schur_solver, jxf_ParVector *rhs, jxf_ParVector *x, JXF_Int *u_end);
JXF_Int jxf_ILUSolveLU(jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u, JXF_Int *perm, JXF_Int nLU, jxf_ParCSRMatrix *L, JXF_Real* D, jxf_ParCSRMatrix *U, jxf_ParVector *utemp, jxf_ParVector *ftemp);
JXF_Int jxf_ILUSolveLURAS(jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u, JXF_Int *perm, jxf_ParCSRMatrix *L, JXF_Real* D, jxf_ParCSRMatrix *U, jxf_ParVector *ftemp, jxf_ParVector *utemp, JXF_Real *fext, JXF_Real *uext);
JXF_Int jxf_NSHSolve( void *nsh_vdata, jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u );
JXF_Int jxf_NSHSolveInverse(jxf_ParCSRMatrix *A, jxf_ParVector *f, jxf_ParVector *u, jxf_ParCSRMatrix *M, jxf_ParVector *ftemp, jxf_ParVector *utemp);

#endif
