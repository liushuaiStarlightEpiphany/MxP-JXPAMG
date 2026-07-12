//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_apctl.h -- head files for 3t solver and preconditioners
 *  Date: 2011/09/08
 *
 *  The 3T matrices are of the following form:
 *    /           \      /           \  
 *   | A11 A12  0  |    | Arr Are  0  |  
 *   | A21 A22 A23 | or | Aer Aee Aei |
 *   |  0  A32 A33 |    |  0  Aie Aii |
 *    \           /      \           /   
 * 
 *  Created by peghoty
 */ 

#ifndef JXF_3T_HEADER
#define JXF_3T_HEADER 

#ifndef JXF_KRYLOV_HEADER
#include "jxf_krylov.h"
#endif

/* relax type for a submatrix */
#define RELAX_AMG      0
#define RELAX_WJACOBI  1

/* solver type for a submatrix */
#define SOLVER_AMG      0
#define SOLVER_AMGGMRES 1

/* block smoothing type for a 3t matrix */
#define BLOCKSMOOTH_GS 0
#define BLOCKSMOOTH_BD 1

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_3tAPCTLData
 */
typedef struct
{
   JXF_Real   tol;                    // tolerance    
   JXF_Int      max_iter;               // maximal number of iteration
   JXF_Int      num_relax_pre;          // number of pre-relaxations 
   JXF_Int      num_relax_post;         // number of post-relaxations    
   JXF_Int      print_level;            // 1: print CPU time 2: print number of iterations 

   JXF_Int      num_iterations;         // number of iterations
   JXF_Real   last_rel_nrm;           // last relative residual norm
   JXF_Real   ave_conv_factor;        // average convergence factor
   
   /* for processors-grouping */
   MPI_Comm comm;
   MPI_Comm comm_x;
   MPI_Comm comm_y;
   
   JXF_Int      groupid_x; 
   JXF_Int      groupid_y;
 
   JXF_Int      np_R;     // number of processors used for ARR
   JXF_Int      np_E;     // number of processors used for AEE
   JXF_Int      np_I;     // number of processors used for AII      

   JXF_Int      num_group;

   JXF_Int      reset_zero;

   /* flag to indicate whether diagonal elements of the three DiagonalBlock
      matrices are firstly stored in each row for the CSR format. peghoty, 2012/03/06 */
   JXF_Int    is_diagelm_first;

   /* the 3t matrix to be solved */
   jxf_ParCSRMatrix *A;
   
   /* Blocks of the 3t matrix A (基于进程组) */
   jxf_ParCSRMatrix *ARR;
   jxf_ParCSRMatrix *AEE;
   jxf_ParCSRMatrix *AII;

   jxf_ParVector    *VRE;
   jxf_ParVector    *VER;
   jxf_ParVector    *VEI;
   jxf_ParVector    *VIE;   

   jxf_ParVector   **VER2;

   /* Blocks of the 3t matrix A (基于所有进程) peghoty, 2012/02/26 */
   jxf_ParCSRMatrix *ARR_all;
   jxf_ParCSRMatrix *AEE_all;
   jxf_ParCSRMatrix *AII_all;

   jxf_ParCSRMatrix **ARR_all2;

   jxf_ParVector    *VRE_all;
   jxf_ParVector    *VER_all;
   jxf_ParVector    *VEI_all;
   jxf_ParVector    *VIE_all; 

   jxf_ParVector   **VRE_all2;
   jxf_ParVector   **VER_all2;

   /* prolongation operator */
   jxf_ParCSRMatrix *P;  
   jxf_ParVector    *PRR;
   jxf_ParVector    *PII;
       
   /* coarse matrix ACC */
   jxf_ParCSRMatrix *ACC; 
   
   /* relaxation type of each subblock */
   JXF_Int     ARR_relax_type;
   JXF_Int     AEE_relax_type;
   JXF_Int     AII_relax_type;           
     
   /* max_iter for each subblock */
   JXF_Int     ARR_interp_maxit;     
   JXF_Int     AII_interp_maxit;
   JXF_Int     ARR_relax_maxit;
   JXF_Int     AEE_relax_maxit;
   JXF_Int     AII_relax_maxit;   
   JXF_Int     ACC_relax_maxit;

   JXF_Int     maxit_default;

   /* tolerace for each subblock */
   JXF_Real  ARR_interp_tol;   
   JXF_Real  AII_interp_tol;
   JXF_Real  ARR_relax_tol;
   JXF_Real  AEE_relax_tol;
   JXF_Real  AII_relax_tol;
   JXF_Real  ACC_relax_tol;

   JXF_Real  tol_default;

   /* solver type for interpolation-building of PRR.  peghoty,2011/10/29 */
   JXF_Int ARR_solver_id; // ARR sometimes is very hard to converge by PAMG, PAMG-GMRES is provided as an alternative

   /* solver type for PCTL iteration.  peghoty,2011/10/29 */
   JXF_Int ACC_solver_id; // ACC sometimes is very hard to converge by PAMG, PAMG-GMRES is provided as an alternative
   
   /* restart parameters for GMRES solver. peghoty,2011/10/29 */
   JXF_Int ARR_kdim;
   JXF_Int ACC_kdim;
    
   /* AMG solver data for each subblock */
   jxf_ParAMGData    *ARR_amg_solver;
   jxf_ParAMGData    *AEE_amg_solver;
   jxf_ParAMGData    *AII_amg_solver;
   jxf_ParAMGData    *ACC_amg_solver;
   
   /* gmres solver data for ARR and ACC. peghoty,2011/10/29 */
   jxf_GMRESData     *ARR_gmres_solver;
   jxf_GMRESData     *ACC_gmres_solver;   

   /* Auxialiary vectors */
   jxf_ParVector     *WRR;
   jxf_ParVector     *WEE;
   jxf_ParVector     *WII;
   jxf_ParVector     *WCC;

   jxf_ParVector     *GCC;

   jxf_ParVector     *RES;
   jxf_ParVector     *RHS;
   
   jxf_ParVector     *JAC;  // used in Jacobi relaxation

   JXF_Real *tTEMP;

   /* parameters to describe the weaking coupling */
   JXF_Real theta_wc_R;
   JXF_Real theta_wc_E;
   JXF_Real theta_wc_I;
   JXF_Real threshold_wc_R;
   JXF_Real threshold_wc_E;
   JXF_Real threshold_wc_I;
   JXF_Int    IS_WC_R;
   JXF_Int    IS_WC_E;    
   JXF_Int    IS_WC_I;
   
   /* parameters to describe the diagonal dominance */
   JXF_Real theta_dd_R;
   JXF_Real theta_dd_E;
   JXF_Real theta_dd_I;
   JXF_Real threshold_dd_R;
   JXF_Real threshold_dd_E;
   JXF_Real threshold_dd_I; 
   JXF_Int    IS_DD_R;
   JXF_Int    IS_DD_E;    
   JXF_Int    IS_DD_I;
   
   /* flag to indicate whether the Coarse-grid Correction is necessary */
   JXF_Int    Need_CC;
   
   /* type of the block smoothing */
   JXF_Int    blocksmooth_type;
   
   /* fixed number of iterations and tolerance for submatrices */ 
   JXF_Int    fixit_pctl_R;
   JXF_Int    fixit_pctl_E;
   JXF_Int    fixit_pctl_I;
   JXF_Int    fixit_brlx_R;
   JXF_Int    fixit_brlx_E;
   JXF_Int    fixit_brlx_I;
   
   /* whether employ the fixed-number-of-iterations mode? peghoty, 2012/02/15 */
   JXF_Int    use_fixedmode_R;
   JXF_Int    use_fixedmode_E;
   JXF_Int    use_fixedmode_I;
   
   /* Whether use the pure PCTL? */
   JXF_Int    use_ppctl;
   
   /* Whether test the number of iterations for each sub-linearsystem solution? */
   JXF_Int    test_subls_iter;

   /* variables to count the number of iteration */
   JXF_Int    num_iter_Ai_pctl_setup;
   JXF_Int    num_iter_Ar_pctl_setup;
   JXF_Int    num_iter_Ae_pctl_precond;
   JXF_Int    num_iter_Ai_pctl_precond;
   JXF_Int    num_iter_Ar_pctl_precond;
   JXF_Int    num_iter_Ac_pctl_precond;    

   JXF_Int    debug_flag;

   JXF_Real strong_threshold;

   JXF_Int interp_type;
   JXF_Int coarsen_type;
   JXF_Int agg_num_levels;
   JXF_Int coarse_threshold;
   JXF_Int print_level_amg;
    
} jxf_3tAPCTLData;

#define jxf_3tAPCTLDatatTEMP(pre_3tapctl_data)               ((pre_3tapctl_data) -> tTEMP)
#define jxf_3tAPCTLDataTol(pre_3tapctl_data)                 ((pre_3tapctl_data) -> tol)
#define jxf_3tAPCTLDataMaxIter(pre_3tapctl_data)             ((pre_3tapctl_data) -> max_iter)
#define jxf_3tAPCTLDataNumRlxPre(pre_3tapctl_data)           ((pre_3tapctl_data) -> num_relax_pre)
#define jxf_3tAPCTLDataNumRlxPost(pre_3tapctl_data)          ((pre_3tapctl_data) -> num_relax_post)
#define jxf_3tAPCTLDataPrintLevel(pre_3tapctl_data)          ((pre_3tapctl_data) -> print_level)
#define jxf_3tAPCTLDataNumIterations(pre_3tapctl_data)       ((pre_3tapctl_data) -> num_iterations)
#define jxf_3tAPCTLDataLastRelNrm(pre_3tapctl_data)          ((pre_3tapctl_data) -> last_rel_nrm)
#define jxf_3tAPCTLDataAveConvFactor(pre_3tapctl_data)       ((pre_3tapctl_data) -> ave_conv_factor)
#define jxf_3tAPCTLDataComm(pre_3tapctl_data)                ((pre_3tapctl_data) -> comm)
#define jxf_3tAPCTLDataCommX(pre_3tapctl_data)               ((pre_3tapctl_data) -> comm_x)
#define jxf_3tAPCTLDataCommY(pre_3tapctl_data)               ((pre_3tapctl_data) -> comm_y)
#define jxf_3tAPCTLDataGroupIdX(pre_3tapctl_data)            ((pre_3tapctl_data) -> groupid_x)
#define jxf_3tAPCTLDataGroupIdY(pre_3tapctl_data)            ((pre_3tapctl_data) -> groupid_y)
#define jxf_3tAPCTLDataNpR(pre_3tapctl_data)                 ((pre_3tapctl_data) -> np_R)
#define jxf_3tAPCTLDataNpE(pre_3tapctl_data)                 ((pre_3tapctl_data) -> np_E)
#define jxf_3tAPCTLDataNpI(pre_3tapctl_data)                 ((pre_3tapctl_data) -> np_I)
#define jxf_3tAPCTLDataNumGroup(pre_3tapctl_data)            ((pre_3tapctl_data) -> num_group)
#define jxf_3tAPCTLDataResetZero(pre_3tapctl_data)           ((pre_3tapctl_data) -> reset_zero)
#define jxf_3tAPCTLDataIsDiagElmFirst(pre_3tapctl_data)      ((pre_3tapctl_data) -> is_diagelm_first)
#define jxf_3tAPCTLDataA(pre_3tapctl_data)                   ((pre_3tapctl_data) -> A)
#define jxf_3tAPCTLDataARR(pre_3tapctl_data)                 ((pre_3tapctl_data) -> ARR)
#define jxf_3tAPCTLDataAEE(pre_3tapctl_data)                 ((pre_3tapctl_data) -> AEE)
#define jxf_3tAPCTLDataAII(pre_3tapctl_data)                 ((pre_3tapctl_data) -> AII)
#define jxf_3tAPCTLDataVRE(pre_3tapctl_data)                 ((pre_3tapctl_data) -> VRE)
#define jxf_3tAPCTLDataVER(pre_3tapctl_data)                 ((pre_3tapctl_data) -> VER)
#define jxf_3tAPCTLDataVER2(pre_3tapctl_data)                ((pre_3tapctl_data) -> VER2)
#define jxf_3tAPCTLDataVEI(pre_3tapctl_data)                 ((pre_3tapctl_data) -> VEI)
#define jxf_3tAPCTLDataVIE(pre_3tapctl_data)                 ((pre_3tapctl_data) -> VIE)
#define jxf_3tAPCTLDataARRAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> ARR_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataARRAll2(pre_3tapctl_data)             ((pre_3tapctl_data) -> ARR_all2)
#define jxf_3tAPCTLDataAEEAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> AEE_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataAIIAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> AII_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataVREAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> VRE_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataVREAll2(pre_3tapctl_data)             ((pre_3tapctl_data) -> VRE_all2)
#define jxf_3tAPCTLDataVERAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> VER_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataVERAll2(pre_3tapctl_data)             ((pre_3tapctl_data) -> VER_all2)
#define jxf_3tAPCTLDataVEIAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> VEI_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataVIEAll(pre_3tapctl_data)              ((pre_3tapctl_data) -> VIE_all) /* peghoty, 2012/02/26 */
#define jxf_3tAPCTLDataP(pre_3tapctl_data)                   ((pre_3tapctl_data) -> P)
#define jxf_3tAPCTLDataPRR(pre_3tapctl_data)                 ((pre_3tapctl_data) -> PRR)
#define jxf_3tAPCTLDataPII(pre_3tapctl_data)                 ((pre_3tapctl_data) -> PII)
#define jxf_3tAPCTLDataACC(pre_3tapctl_data)                 ((pre_3tapctl_data) -> ACC)
#define jxf_3tAPCTLDataARRRelaxType(pre_3tapctl_data)        ((pre_3tapctl_data) -> ARR_relax_type)
#define jxf_3tAPCTLDataAEERelaxType(pre_3tapctl_data)        ((pre_3tapctl_data) -> AEE_relax_type)
#define jxf_3tAPCTLDataAIIRelaxType(pre_3tapctl_data)        ((pre_3tapctl_data) -> AII_relax_type)
#define jxf_3tAPCTLDataARRInterpMaxIt(pre_3tapctl_data)      ((pre_3tapctl_data) -> ARR_interp_maxit)
#define jxf_3tAPCTLDataAIIInterpMaxIt(pre_3tapctl_data)      ((pre_3tapctl_data) -> AII_interp_maxit)
#define jxf_3tAPCTLDataARRRelaxMaxIt(pre_3tapctl_data)       ((pre_3tapctl_data) -> ARR_relax_maxit)
#define jxf_3tAPCTLDataAEERelaxMaxIt(pre_3tapctl_data)       ((pre_3tapctl_data) -> AEE_relax_maxit)
#define jxf_3tAPCTLDataAIIRelaxMaxIt(pre_3tapctl_data)       ((pre_3tapctl_data) -> AII_relax_maxit)
#define jxf_3tAPCTLDataACCRelaxMaxIt(pre_3tapctl_data)       ((pre_3tapctl_data) -> ACC_relax_maxit)
#define jxf_3tAPCTLDataMaxItDefault(pre_3tapctl_data)        ((pre_3tapctl_data) -> maxit_default)
#define jxf_3tAPCTLDataARRInterpTol(pre_3tapctl_data)        ((pre_3tapctl_data) -> ARR_interp_tol)
#define jxf_3tAPCTLDataAIIInterpTol(pre_3tapctl_data)        ((pre_3tapctl_data) -> AII_interp_tol)
#define jxf_3tAPCTLDataARRRelaxTol(pre_3tapctl_data)         ((pre_3tapctl_data) -> ARR_relax_tol)
#define jxf_3tAPCTLDataAEERelaxTol(pre_3tapctl_data)         ((pre_3tapctl_data) -> AEE_relax_tol)
#define jxf_3tAPCTLDataAIIRelaxTol(pre_3tapctl_data)         ((pre_3tapctl_data) -> AII_relax_tol)
#define jxf_3tAPCTLDataACCRelaxTol(pre_3tapctl_data)         ((pre_3tapctl_data) -> ACC_relax_tol)
#define jxf_3tAPCTLDataTolDefault(pre_3tapctl_data)          ((pre_3tapctl_data) -> tol_default)
#define jxf_3tAPCTLDataARRSolverID(pre_3tapctl_data)         ((pre_3tapctl_data) -> ARR_solver_id)
#define jxf_3tAPCTLDataACCSolverID(pre_3tapctl_data)         ((pre_3tapctl_data) -> ACC_solver_id)
#define jxf_3tAPCTLDataARRKDim(pre_3tapctl_data)             ((pre_3tapctl_data) -> ARR_kdim)
#define jxf_3tAPCTLDataACCKDim(pre_3tapctl_data)             ((pre_3tapctl_data) -> ACC_kdim)  
#define jxf_3tAPCTLDataARRAMGSolver(pre_3tapctl_data)        ((pre_3tapctl_data) -> ARR_amg_solver)  
#define jxf_3tAPCTLDataAEEAMGSolver(pre_3tapctl_data)        ((pre_3tapctl_data) -> AEE_amg_solver)
#define jxf_3tAPCTLDataAIIAMGSolver(pre_3tapctl_data)        ((pre_3tapctl_data) -> AII_amg_solver)
#define jxf_3tAPCTLDataACCAMGSolver(pre_3tapctl_data)        ((pre_3tapctl_data) -> ACC_amg_solver) 
#define jxf_3tAPCTLDataARRGMRESSolver(pre_3tapctl_data)      ((pre_3tapctl_data) -> ARR_gmres_solver)
#define jxf_3tAPCTLDataACCGMRESSolver(pre_3tapctl_data)      ((pre_3tapctl_data) -> ACC_gmres_solver)
#define jxf_3tAPCTLDataWRR(pre_3tapctl_data)                 ((pre_3tapctl_data) -> WRR)
#define jxf_3tAPCTLDataWEE(pre_3tapctl_data)                 ((pre_3tapctl_data) -> WEE)
#define jxf_3tAPCTLDataWII(pre_3tapctl_data)                 ((pre_3tapctl_data) -> WII)
#define jxf_3tAPCTLDataWCC(pre_3tapctl_data)                 ((pre_3tapctl_data) -> WCC)
#define jxf_3tAPCTLDataGCC(pre_3tapctl_data)                 ((pre_3tapctl_data) -> GCC)
#define jxf_3tAPCTLDataRES(pre_3tapctl_data)                 ((pre_3tapctl_data) -> RES)
#define jxf_3tAPCTLDataRHS(pre_3tapctl_data)                 ((pre_3tapctl_data) -> RHS)
#define jxf_3tAPCTLDataJAC(pre_3tapctl_data)                 ((pre_3tapctl_data) -> JAC)
#define jxf_3tAPCTLDataThetaWCR(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_wc_R)
#define jxf_3tAPCTLDataThetaWCE(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_wc_E)
#define jxf_3tAPCTLDataThetaWCI(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_wc_I)
#define jxf_3tAPCTLDataThresholdWCR(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_wc_R)
#define jxf_3tAPCTLDataThresholdWCE(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_wc_E)
#define jxf_3tAPCTLDataThresholdWCI(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_wc_I)
#define jxf_3tAPCTLDataISWCR(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_WC_R)
#define jxf_3tAPCTLDataISWCE(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_WC_E)
#define jxf_3tAPCTLDataISWCI(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_WC_I)
#define jxf_3tAPCTLDataThetaDDR(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_dd_R)
#define jxf_3tAPCTLDataThetaDDE(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_dd_E)
#define jxf_3tAPCTLDataThetaDDI(pre_3tapctl_data)            ((pre_3tapctl_data) -> theta_dd_I)
#define jxf_3tAPCTLDataThresholdDDR(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_dd_R)
#define jxf_3tAPCTLDataThresholdDDE(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_dd_E)
#define jxf_3tAPCTLDataThresholdDDI(pre_3tapctl_data)        ((pre_3tapctl_data) -> threshold_dd_I)
#define jxf_3tAPCTLDataISDDR(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_DD_R)
#define jxf_3tAPCTLDataISDDE(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_DD_E)
#define jxf_3tAPCTLDataISDDI(pre_3tapctl_data)               ((pre_3tapctl_data) -> IS_DD_I)
#define jxf_3tAPCTLDataNeedCC(pre_3tapctl_data)              ((pre_3tapctl_data) -> Need_CC)
#define jxf_3tAPCTLDataBlockSmoothType(pre_3tapctl_data)     ((pre_3tapctl_data) -> blocksmooth_type)
#define jxf_3tAPCTLDataFixItPCTLR(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_pctl_R)
#define jxf_3tAPCTLDataFixItPCTLE(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_pctl_E)
#define jxf_3tAPCTLDataFixItPCTLI(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_pctl_I)
#define jxf_3tAPCTLDataFixItBRLXR(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_brlx_R)
#define jxf_3tAPCTLDataFixItBRLXE(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_brlx_E)
#define jxf_3tAPCTLDataFixItBRLXI(pre_3tapctl_data)          ((pre_3tapctl_data) -> fixit_brlx_I)
#define jxf_3tAPCTLDataUseFixedModeR(pre_3tapctl_data)       ((pre_3tapctl_data) -> use_fixedmode_R)
#define jxf_3tAPCTLDataUseFixedModeE(pre_3tapctl_data)       ((pre_3tapctl_data) -> use_fixedmode_E)
#define jxf_3tAPCTLDataUseFixedModeI(pre_3tapctl_data)       ((pre_3tapctl_data) -> use_fixedmode_I)
#define jxf_3tAPCTLDataUsePPCTL(pre_3tapctl_data)            ((pre_3tapctl_data) -> use_ppctl)
// peghoty, 2012/03/23
#define jxf_3tAPCTLDataTestSubLSIter(pre_3tapctl_data)       ((pre_3tapctl_data) -> test_subls_iter)
#define jxf_3tAPCTLDataNumIterAiSetup(pre_3tapctl_data)      ((pre_3tapctl_data) -> num_iter_Ai_pctl_setup)
#define jxf_3tAPCTLDataNumIterArSetup(pre_3tapctl_data)      ((pre_3tapctl_data) -> num_iter_Ar_pctl_setup)
#define jxf_3tAPCTLDataNumIterAePrecond(pre_3tapctl_data)    ((pre_3tapctl_data) -> num_iter_Ae_pctl_precond)
#define jxf_3tAPCTLDataNumIterAiPrecond(pre_3tapctl_data)    ((pre_3tapctl_data) -> num_iter_Ai_pctl_precond)
#define jxf_3tAPCTLDataNumIterArPrecond(pre_3tapctl_data)    ((pre_3tapctl_data) -> num_iter_Ar_pctl_precond)
#define jxf_3tAPCTLDataNumIterAcPrecond(pre_3tapctl_data)    ((pre_3tapctl_data) -> num_iter_Ac_pctl_precond)
#define jxf_3tAPCTLDataDebugFlag(pre_3tapctl_data)           ((pre_3tapctl_data) -> debug_flag)
#define jxf_3tAPCTLDataStrongThreshold(pre_3tapctl_data)     ((pre_3tapctl_data) -> strong_threshold)
#define jxf_3tAPCTLDataInterpType(pre_3tapctl_data)          ((pre_3tapctl_data) -> interp_type)
#define jxf_3tAPCTLDataCoarsenType(pre_3tapctl_data)         ((pre_3tapctl_data) -> coarsen_type)
#define jxf_3tAPCTLDataAggNumLevels(pre_3tapctl_data)        ((pre_3tapctl_data) -> agg_num_levels)
#define jxf_3tAPCTLDataCoarseThreshold(pre_3tapctl_data)     ((pre_3tapctl_data) -> coarse_threshold)
#define jxf_3tAPCTLDataPrintLevelAMG(pre_3tapctl_data)     ((pre_3tapctl_data) -> print_level_amg)


/*!
 * \struct jxf_APCTLKrylovParam
 * \brief parameters for APCTL-Krylov method
 */
typedef struct
{
   JXF_Int       num_group;

   JXF_Int       reset_zero;

   //----------------------------------------------------------------------------
   //  parameters for P-Krylov methods
   //--------------------------------------------------------------------------
   JXF_Int       solver_id;           // 1: PCG; 2: PGMRES; 3: PBicgStab
   JXF_Int       precond_id;           // 1: APCTL; 2: Schur1; 3: Schur2
   JXF_Real    tol;                 // tolerance of the APCTL-Krylov method
   JXF_Int       max_iter;            // maximal number of iteration
   JXF_Int       k_dim;               // number of restart
   JXF_Int       is_check_restarted;  // peghoty, 2011/11/08
   JXF_Int       two_norm;            // for PCG
   JXF_Int       print_level;         // how many info to be output?
   JXF_Int       TTest;               // whether timing the program?
   JXF_Int       keepsol;             // whether save the solution?
   JXF_Int       num_iterations;      // number of iterations of the APCTL-Krylov method
   JXF_Real    last_rel_nrm;        // last relative residual norm     peghoty,2012/03/22

   //----------------------------------------------------------------------------
   //  parameters for PCTL iteration or preconditioner
   //--------------------------------------------------------------------------   
   JXF_Int    print_level_apctl; // how much info to be output in apctl?
                             // 1: CPU information
                             // 2: inner iteration information
                             // 3: both CPU and inner iteration information
   JXF_Int    blocksmooth_type;  // BD or GS type preconditioner when Coarse Correction is not needed?
                             // BLOCKSMOOTH_BD:
                             // BLOCKSMOOTH_GS:
   JXF_Int    num_relax_pre;     // number of sweeps for pre-block-smoothing
   JXF_Int    num_relax_post;    // number of sweeps for post-block-smoothing

   /* solver type and restart number for interpolation-building of PRR */
   JXF_Int    interp_solver_ARR; // 0: SOLVER_AMG; 1: SOLVER_AMGGMRES;   
   JXF_Int    interp_kdim_ARR;   // restart parameters for GMRES solver

   JXF_Int    interp_maxit_ARR;  // maximal number of iteration for ARR to build PRR
   JXF_Int    interp_maxit_AII;  // maximal number of iteration for AII to build PII
   JXF_Real interp_tol_ARR;    // tolerance for ARR to build PRR
   JXF_Real interp_tol_AII;    // tolerance for AII to build PII
   
   JXF_Real  ARR_relax_tol;
   JXF_Real  AEE_relax_tol;
   JXF_Real  AII_relax_tol;
   JXF_Real  ACC_relax_tol;

   JXF_Int    ARR_relax_type;
   JXF_Int    AEE_relax_type;
   JXF_Int    AII_relax_type;

   JXF_Int    fixit_pctl_R;      // fixed number of iterations for ARR in PCTL
   JXF_Int    fixit_pctl_E;      // fixed number of iterations for AEE in PCTL
   JXF_Int    fixit_pctl_I;      // fixed number of iterations for AII in PCTL
   JXF_Int    fixit_brlx_R;      // fixed number of iterations for ARR in Block Relaxation
   JXF_Int    fixit_brlx_E;      // fixed number of iterations for AEE in Block Relaxation
   JXF_Int    fixit_brlx_I;      // fixed number of iterations for AII in Block Relaxation     

   /* whether employ the fixed-number-of-iterations mode? peghoty, 2012/02/15 */
   JXF_Int    use_fixedmode_R;
   JXF_Int    use_fixedmode_E;
   JXF_Int    use_fixedmode_I; 
   
   /* parameters to describe the weaking coupling between AEE and VER, VEI */
   JXF_Real theta_wc_E;
   JXF_Real threshold_wc_E;
   
   /* parameters to describe the diagonal dominance */
   JXF_Real theta_dd_R;
   JXF_Real theta_dd_E;
   JXF_Real theta_dd_I;
   JXF_Real threshold_dd_R;
   JXF_Real threshold_dd_E;
   JXF_Real threshold_dd_I;
   
   /* flag to indicate whether diagonal elements of the three DiagonalBlock
      matrices are firstly stored in each row for the CSR format. peghoty, 2012/03/06 */
   JXF_Int    is_diagelm_first;
   
   /* Whether use the pure PCTL? */
   JXF_Int    use_ppctl;

   /* Whether test the number of iterations for each sub-linearsystem solution? */
   JXF_Int    test_subls_iter;
   
   /* variables to count the number of iteration */
   JXF_Int    num_iter_Ai_pctl_setup;
   JXF_Int    num_iter_Ar_pctl_setup;
   JXF_Int    num_iter_Ae_pctl_precond;
   JXF_Int    num_iter_Ai_pctl_precond;
   JXF_Int    num_iter_Ar_pctl_precond;
   JXF_Int    num_iter_Ac_pctl_precond;
   
   /* whether need coarse-grid correction in apctl preconditioner? */
   JXF_Int    Need_CC;
   
   /* CPU time of each phase */
   JXF_Real cpu_trans;
   JXF_Real cpu_setup;
   JXF_Real cpu_solve;
   JXF_Real cpu_total;                

   JXF_Int debug_flag;

   JXF_Real strong_threshold;

   JXF_Int interp_type;
   JXF_Int coarsen_type;
   JXF_Int agg_num_levels;
   JXF_Int coarse_threshold;
   JXF_Int print_level_amg;

} jxf_APCTLKrylovParam;

#define jxf_APCTLKrylovParamNumGroup(apctlkrylov_param)           ((apctlkrylov_param) -> num_group)
#define jxf_APCTLKrylovParamResetZero(apctlkrylov_param)          ((apctlkrylov_param) -> reset_zero)
#define jxf_APCTLKrylovParamSolverID(apctlkrylov_param)           ((apctlkrylov_param) -> solver_id)
#define jxf_APCTLKrylovParamPrecondID(apctlkrylov_param)          ((apctlkrylov_param) -> precond_id)
#define jxf_APCTLKrylovParamTol(apctlkrylov_param)                ((apctlkrylov_param) -> tol)
#define jxf_APCTLKrylovParamMaxIter(apctlkrylov_param)            ((apctlkrylov_param) -> max_iter)
#define jxf_APCTLKrylovParamKDim(apctlkrylov_param)               ((apctlkrylov_param) -> k_dim)
#define jxf_APCTLKrylovParamISCheckRestarted(apctlkrylov_param)   ((apctlkrylov_param) -> is_check_restarted)
#define jxf_APCTLKrylovParamTwoNorm(apctlkrylov_param)            ((apctlkrylov_param) -> two_norm)
#define jxf_APCTLKrylovParamPrintLevel(apctlkrylov_param)         ((apctlkrylov_param) -> print_level)
#define jxf_APCTLKrylovParamTTest(apctlkrylov_param)              ((apctlkrylov_param) -> TTest)
#define jxf_APCTLKrylovParamKeepSol(apctlkrylov_param)            ((apctlkrylov_param) -> keepsol)
#define jxf_APCTLKrylovParamNumIterations(apctlkrylov_param)      ((apctlkrylov_param) -> num_iterations)
#define jxf_APCTLKrylovParamLastRelNrm(apctlkrylov_param)         ((apctlkrylov_param) -> last_rel_nrm) // 2012/03/22
#define jxf_APCTLKrylovParamPrintLevelAPCTL(apctlkrylov_param)    ((apctlkrylov_param) -> print_level_apctl)
#define jxf_APCTLKrylovParamBlockSmoothType(apctlkrylov_param)    ((apctlkrylov_param) -> blocksmooth_type)
#define jxf_APCTLKrylovParamNumRelaxPre(apctlkrylov_param)        ((apctlkrylov_param) -> num_relax_pre)
#define jxf_APCTLKrylovParamNumRelaxPost(apctlkrylov_param)       ((apctlkrylov_param) -> num_relax_post)
#define jxf_APCTLKrylovParamInterpSolverARR(apctlkrylov_param)    ((apctlkrylov_param) -> interp_solver_ARR)
#define jxf_APCTLKrylovParamInterpKdimARR(apctlkrylov_param)      ((apctlkrylov_param) -> interp_kdim_ARR)
#define jxf_APCTLKrylovParamInterpMaxitARR(apctlkrylov_param)     ((apctlkrylov_param) -> interp_maxit_ARR)
#define jxf_APCTLKrylovParamInterpMaxitAII(apctlkrylov_param)     ((apctlkrylov_param) -> interp_maxit_AII)
#define jxf_APCTLKrylovParamInterpTolARR(apctlkrylov_param)       ((apctlkrylov_param) -> interp_tol_ARR)
#define jxf_APCTLKrylovParamInterpTolAII(apctlkrylov_param)       ((apctlkrylov_param) -> interp_tol_AII)
#define jxf_APCTLKrylovParamARRRelaxTol(apctlkrylov_param)        ((apctlkrylov_param) -> ARR_relax_tol)
#define jxf_APCTLKrylovParamAEERelaxTol(apctlkrylov_param)        ((apctlkrylov_param) -> AEE_relax_tol)
#define jxf_APCTLKrylovParamAIIRelaxTol(apctlkrylov_param)        ((apctlkrylov_param) -> AII_relax_tol)
#define jxf_APCTLKrylovParamACCRelaxTol(apctlkrylov_param)        ((apctlkrylov_param) -> ACC_relax_tol)
#define jxf_APCTLKrylovParamFixitPCTLR(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_pctl_R)
#define jxf_APCTLKrylovParamFixitPCTLE(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_pctl_E)
#define jxf_APCTLKrylovParamFixitPCTLI(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_pctl_I)
#define jxf_APCTLKrylovParamFixitBrlxR(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_brlx_R)
#define jxf_APCTLKrylovParamFixitBrlxE(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_brlx_E)
#define jxf_APCTLKrylovParamFixitBrlxI(apctlkrylov_param)         ((apctlkrylov_param) -> fixit_brlx_I)
#define jxf_APCTLKrylovParamARRRelaxType(apctlkrylov_param)       ((apctlkrylov_param) -> ARR_relax_type)
#define jxf_APCTLKrylovParamAEERelaxType(apctlkrylov_param)       ((apctlkrylov_param) -> AEE_relax_type)
#define jxf_APCTLKrylovParamAIIRelaxType(apctlkrylov_param)       ((apctlkrylov_param) -> AII_relax_type)
#define jxf_APCTLKrylovParamUseFixedModeR(apctlkrylov_param)      ((apctlkrylov_param) -> use_fixedmode_R)
#define jxf_APCTLKrylovParamUseFixedModeE(apctlkrylov_param)      ((apctlkrylov_param) -> use_fixedmode_E)
#define jxf_APCTLKrylovParamUseFixedModeI(apctlkrylov_param)      ((apctlkrylov_param) -> use_fixedmode_I)
#define jxf_APCTLKrylovParamThetaWCE(apctlkrylov_param)           ((apctlkrylov_param) -> theta_wc_E)
#define jxf_APCTLKrylovParamThresholdWCE(apctlkrylov_param)       ((apctlkrylov_param) -> threshold_wc_E)
#define jxf_APCTLKrylovParamThetaDDR(apctlkrylov_param)           ((apctlkrylov_param) -> theta_dd_R)
#define jxf_APCTLKrylovParamThetaDDE(apctlkrylov_param)           ((apctlkrylov_param) -> theta_dd_E)
#define jxf_APCTLKrylovParamThetaDDI(apctlkrylov_param)           ((apctlkrylov_param) -> theta_dd_I)
#define jxf_APCTLKrylovParamThresholdDDR(apctlkrylov_param)       ((apctlkrylov_param) -> threshold_dd_R)
#define jxf_APCTLKrylovParamThresholdDDE(apctlkrylov_param)       ((apctlkrylov_param) -> threshold_dd_E)
#define jxf_APCTLKrylovParamThresholdDDI(apctlkrylov_param)       ((apctlkrylov_param) -> threshold_dd_I)
#define jxf_APCTLKrylovParamISDiagElmFirst(apctlkrylov_param)     ((apctlkrylov_param) -> is_diagelm_first)
#define jxf_APCTLKrylovParamUsePPCTL(apctlkrylov_param)           ((apctlkrylov_param) -> use_ppctl)
// peghoty, 2012/03/23
#define jxf_APCTLKrylovParamTestSubLSIter(apctlkrylov_param)      ((apctlkrylov_param) -> test_subls_iter)
#define jxf_APCTLKrylovParamNumIterAiSetup(apctlkrylov_param)     ((apctlkrylov_param) -> num_iter_Ai_pctl_setup)
#define jxf_APCTLKrylovParamNumIterArSetup(apctlkrylov_param)     ((apctlkrylov_param) -> num_iter_Ar_pctl_setup)
#define jxf_APCTLKrylovParamNumIterAePrecond(apctlkrylov_param)   ((apctlkrylov_param) -> num_iter_Ae_pctl_precond)
#define jxf_APCTLKrylovParamNumIterAiPrecond(apctlkrylov_param)   ((apctlkrylov_param) -> num_iter_Ai_pctl_precond)
#define jxf_APCTLKrylovParamNumIterArPrecond(apctlkrylov_param)   ((apctlkrylov_param) -> num_iter_Ar_pctl_precond)
#define jxf_APCTLKrylovParamNumIterAcPrecond(apctlkrylov_param)   ((apctlkrylov_param) -> num_iter_Ac_pctl_precond)
#define jxf_APCTLKrylovParamNeedCC(apctlkrylov_param)             ((apctlkrylov_param) -> Need_CC)
// peghoty, 2012/03/25
#define jxf_APCTLKrylovParamCPUTrans(apctlkrylov_param)           ((apctlkrylov_param) -> cpu_trans)
#define jxf_APCTLKrylovParamCPUSetup(apctlkrylov_param)           ((apctlkrylov_param) -> cpu_setup)
#define jxf_APCTLKrylovParamCPUSolve(apctlkrylov_param)           ((apctlkrylov_param) -> cpu_solve)
#define jxf_APCTLKrylovParamCPUTotal(apctlkrylov_param)           ((apctlkrylov_param) -> cpu_total)
#define jxf_APCTLKrylovParamDebugFlag(apctlkrylov_param)          ((apctlkrylov_param) -> debug_flag)
#define jxf_APCTLKrylovParamStrongThreshold(apctlkrylov_param)    ((apctlkrylov_param) -> strong_threshold)
#define jxf_APCTLKrylovParamInterpType(apctlkrylov_param)         ((apctlkrylov_param) -> interp_type)
#define jxf_APCTLKrylovParamCoarsenType(apctlkrylov_param)        ((apctlkrylov_param) -> coarsen_type)
#define jxf_APCTLKrylovParamAggNumLevels(apctlkrylov_param)       ((apctlkrylov_param) -> agg_num_levels)
#define jxf_APCTLKrylovParamCoarseThreshold(apctlkrylov_param)    ((apctlkrylov_param) -> coarse_threshold)
#define jxf_APCTLKrylovParamPrintLevelAMG(apctlkrylov_param)    ((apctlkrylov_param) -> print_level_amg)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/
 
/* /csrc/3t/subvec.c */
JXF_Int
jxf_3tGetSubVecs( jxf_Vector  *f, 
                 jxf_Vector **fR_ptr, 
                 jxf_Vector **fE_ptr, 
                 jxf_Vector **fI_ptr );
JXF_Int
jxf_mgGetSubVecs( jxf_Vector  *f,
                 jxf_Vector **fR,
                 jxf_Vector **fE_ptr,
                 jxf_Vector **fI_ptr,
                 JXF_Int      ng );

/* /csrc/3t/submat.c */                   
JXF_Int
jxf_3tGetSubBlocks( jxf_CSRMatrix  *A, 
                   jxf_CSRMatrix **ARR_ptr, 
                   jxf_CSRMatrix **AEE_ptr, 
                   jxf_CSRMatrix **AII_ptr, 
                   jxf_CSRMatrix **ARE_ptr, 
                   jxf_CSRMatrix **AER_ptr, 
                   jxf_CSRMatrix **AEI_ptr, 
                   jxf_CSRMatrix **AIE_ptr );
JXF_Int
jxf_3tGetSubBlocks_REIV( jxf_CSRMatrix  *A, 
                        jxf_CSRMatrix **ARR_ptr, 
                        jxf_CSRMatrix **AEE_ptr, 
                        jxf_CSRMatrix **AII_ptr, 
                        jxf_Vector    **VRE_ptr, 
                        jxf_Vector    **VER_ptr, 
                        jxf_Vector    **VEI_ptr, 
                        jxf_Vector    **VIE_ptr );
JXF_Int
jxf_mgGetSubBlocks_REIV( jxf_CSRMatrix  *A,
                        jxf_CSRMatrix **ARRarray,
                        jxf_CSRMatrix **AEE_ptr,
                        jxf_CSRMatrix **AII_ptr,
                        jxf_Vector    **VREarray,
                        jxf_Vector    **VERarray,
                        jxf_Vector    **VEI_ptr,
                        jxf_Vector    **VIE_ptr,
                        JXF_Int         ng );
JXF_Int
jxf_3tGetSubDiagBlocks( jxf_CSRMatrix  *A, 
                       jxf_CSRMatrix **ARR_ptr, 
                       jxf_CSRMatrix **AEE_ptr, 
                       jxf_CSRMatrix **AII_ptr );                                        
JXF_Int
jxf_3tGetSubBlocks_R( jxf_CSRMatrix  *A, 
                     jxf_CSRMatrix **ARR_ptr, 
                     jxf_CSRMatrix **ARE_ptr );
JXF_Int
jxf_3tGetSubBlocks_RV( jxf_CSRMatrix  *A, 
                      jxf_CSRMatrix **ARR_ptr, 
                      jxf_Vector    **VRE_ptr );                                   
JXF_Int
jxf_3tGetSubBlocks_E( jxf_CSRMatrix  *A, 
                     jxf_CSRMatrix **AEE_ptr, 
                     jxf_CSRMatrix **AER_ptr, 
                     jxf_CSRMatrix **AEI_ptr );
JXF_Int
jxf_3tGetSubBlocks_EV( jxf_CSRMatrix  *A, 
                      jxf_CSRMatrix **AEE_ptr, 
                      jxf_Vector    **VER_ptr, 
                      jxf_Vector    **VEI_ptr );                   
JXF_Int
jxf_3tGetSubBlocks_I( jxf_CSRMatrix  *A, 
                     jxf_CSRMatrix **AII_ptr, 
                     jxf_CSRMatrix **AIE_ptr );
JXF_Int
jxf_3tGetSubBlocks_IV( jxf_CSRMatrix  *A, 
                      jxf_CSRMatrix **AII_ptr, 
                      jxf_Vector    **VIE_ptr );
JXF_Int jxf_3tGetSubBlocks_ARR( jxf_CSRMatrix *A, jxf_CSRMatrix **ARR_ptr );                 
JXF_Int jxf_3tGetSubBlocks_AEE( jxf_CSRMatrix *A, jxf_CSRMatrix **AEE_ptr );
JXF_Int jxf_3tGetSubBlocks_AII( jxf_CSRMatrix *A, jxf_CSRMatrix **AII_ptr );                    
JXF_Int jxf_3tGetSubBlocks_DII( jxf_CSRMatrix *A, jxf_Vector **DII_ptr );  
jxf_ParCSRMatrix *jxf_3tGetAD( jxf_ParCSRMatrix *A, JXF_Int AII_all );
JXF_Int jxf_3tAPCTLDDCheck( JXF_Real theta_dd, JXF_Real threshold_dd, jxf_ParCSRMatrix *A );
JXF_Int 
jxf_3tDDCheckParallel( jxf_CSRMatrix  *A_diag, 
                      jxf_CSRMatrix  *A_offd, 
                      JXF_Real         theta_dd, 
                      JXF_Real         threshold_dd );                                            

JXF_Int 
jxf_3tAPCTLWeakCouplingE( JXF_Real            theta_wc_E, 
                         JXF_Real            threshold_wc_E, 
                         jxf_ParCSRMatrix  *AEE, 
                         jxf_ParVector     *VER, 
                         jxf_ParVector     *VEI );
JXF_Int 
jxf_3tAPCTLmgWeakCouplingE( JXF_Real            theta_wc_E, 
                           JXF_Real            threshold_wc_E, 
                           JXF_Int            ng,
                           jxf_ParCSRMatrix  *AEE, 
                           jxf_ParVector    **VER, 
                           jxf_ParVector     *VEI );
JXF_Int 
jxf_3tWeakCouplingParallel( jxf_ParCSRMatrix  *A, 
                           jxf_ParVector     *V,
                           JXF_Real            theta_wc );                     
JXF_Int
jxf_3tAPCTLWCDD( jxf_3tAPCTLData   *pre_3tapctl_data,
                MPI_Comm          comm,
                jxf_CSRMatrix      *ARR, 
                jxf_CSRMatrix      *AEE,
                jxf_CSRMatrix      *AII,
                jxf_Vector         *VRE,
                jxf_Vector         *VER,
                jxf_Vector         *VEI,
                jxf_Vector         *VIE  );
JXF_Int
jxf_3tWeakCoupling( jxf_CSRMatrix  *A, 
                   jxf_Vector     *V, 
                   JXF_Real         theta_wc, 
                   JXF_Real         threshold_wc );
JXF_Int 
jxf_3tDiagDominant( jxf_CSRMatrix  *A,
                   JXF_Real         theta_dd,
                   JXF_Real         threshold_dd );
JXF_Int
jxf_BSCModifySubMat( jxf_ParCSRMatrix *A, jxf_ParVector *V, JXF_Real *TEMP, JXF_Int hsize );
JXF_Int
jxf_BSCModify2SubMat( jxf_ParCSRMatrix *A, jxf_ParVector *V, JXF_Real *sTEMP, JXF_Real *vTEMP, JXF_Int hsize );

/* /csrc/3t/transf.c */                 
JXF_Int 
jxf_APCTLKrylovSolBack4Jasmin( MPI_Comm      comm,
                              MPI_Comm      comm_x, 
                              JXF_Int           groupid_x, 
                              jxf_ParVector *par_sol, 
                              jxf_ParVector *uR_p, 
                              jxf_ParVector *uE_p, 
                              jxf_ParVector *uI_p );
JXF_Int 
jxf_APCTLKrylovSolBack4mgJasmin( MPI_Comm      comm,
                                MPI_Comm      comm_x,
                                JXF_Int           groupid_x,
                                JXF_Int           ng,
                                jxf_ParVector  *par_sol,
                                jxf_ParVector **uR_p,
                                jxf_ParVector  *uE_p,
                                jxf_ParVector  *uI_p );
JXF_Int 
jxf_APCTLKrylovSolBack4Jasmin2( MPI_Comm      comm,
                               jxf_ParVector *par_sol, 
                               jxf_ParVector *uR_p, 
                               jxf_ParVector *uE_p, 
                               jxf_ParVector *uI_p );
JXF_Int 
jxf_MatVecGroupR( MPI_Comm           comm_x,
                 jxf_ParCSRMatrix   *A, 
                 jxf_ParCSRMatrix  **ARR_ptr, 
                 jxf_ParVector     **VRE_ptr );
JXF_Int 
jxf_MatGroupR( MPI_Comm           comm_x,
              jxf_ParCSRMatrix   *A, 
              jxf_ParCSRMatrix  **ARR_ptr );
JXF_Int 
jxf_MatVecGroupI( MPI_Comm           comm_x,
                 jxf_ParCSRMatrix   *A, 
                 jxf_ParCSRMatrix  **AII_ptr, 
                 jxf_ParVector     **VIE_ptr );
JXF_Int 
jxf_MatGroupI( MPI_Comm           comm_x,
              jxf_ParCSRMatrix   *A, 
              jxf_ParCSRMatrix  **AII_ptr );                                
JXF_Int 
jxf_MatVecGroupE( MPI_Comm           comm_x,
                 jxf_ParCSRMatrix   *A, 
                 jxf_ParCSRMatrix  **AEE_ptr, 
                 jxf_ParVector     **VER_ptr,
                 jxf_ParVector     **VEI_ptr );
JXF_Int 
jxf_DataCombine4ApctlKrylov( // input:
                            jxf_ParCSRMatrix   *ARR_p, 
                            jxf_ParCSRMatrix   *AEE_p, 
                            jxf_ParCSRMatrix   *AII_p, 
                            jxf_ParVector      *VRE_p, 
                            jxf_ParVector      *VER_p, 
                            jxf_ParVector      *VEI_p, 
                            jxf_ParVector      *VIE_p, 
                            jxf_ParVector      *fR_p, 
                            jxf_ParVector      *fE_p, 
                            jxf_ParVector      *fI_p,
                            jxf_ParVector      *uR_p, 
                            jxf_ParVector      *uE_p, 
                            jxf_ParVector      *uI_p,
                            // output:    
                            jxf_ParCSRMatrix  **A_ptr,  
                            jxf_ParVector     **f_ptr, 
                            jxf_ParVector     **u_ptr );
JXF_Int 
jxf_3tGlobalSystem( jxf_CSRMatrix  *ARR, 
                   jxf_CSRMatrix  *AEE, 
                   jxf_CSRMatrix  *AII, 
                   jxf_Vector     *VRE, 
                   jxf_Vector     *VER, 
                   jxf_Vector     *VEI, 
                   jxf_Vector     *VIE, 
                   jxf_Vector     *fR, 
                   jxf_Vector     *fE, 
                   jxf_Vector     *fI,
                   jxf_Vector     *uR, 
                   jxf_Vector     *uE, 
                   jxf_Vector     *uI,                       
                   jxf_CSRMatrix **A_ptr, 
                   jxf_Vector    **f_ptr,
                   jxf_Vector    **u_ptr );
JXF_Int 
jxf_ParaDataTrans4ApctlKrylov( // input:
                              MPI_Comm           comm, 
                              MPI_Comm           comm_x, 
                              JXF_Int                groupid_x,
                              jxf_ParCSRMatrix   *ARR_p, 
                              jxf_ParCSRMatrix   *AEE_p, 
                              jxf_ParCSRMatrix   *AII_p, 
                              jxf_ParVector      *VRE_p, 
                              jxf_ParVector      *VER_p, 
                              jxf_ParVector      *VEI_p, 
                              jxf_ParVector      *VIE_p, 
                              jxf_ParVector      *fR_p, 
                              jxf_ParVector      *fE_p, 
                              jxf_ParVector      *fI_p,
                              jxf_ParVector      *uR_p, 
                              jxf_ParVector      *uE_p, 
                              jxf_ParVector      *uI_p,
                              // output:    
                              jxf_ParCSRMatrix  **ARR_ptr,
                              jxf_ParCSRMatrix  **AEE_ptr,
                              jxf_ParCSRMatrix  **AII_ptr, 
                              jxf_ParVector     **VRE_ptr, 
                              jxf_ParVector     **VER_ptr, 
                              jxf_ParVector     **VEI_ptr, 
                              jxf_ParVector     **VIE_ptr,
                              jxf_ParCSRMatrix  **A_ptr,  
                              jxf_ParVector     **f_ptr, 
                              jxf_ParVector     **u_ptr );
JXF_Int 
jxf_ParaDataTrans4ApctlmgKrylov( // input:
                                MPI_Comm           comm, 
                                MPI_Comm           comm_x, 
                                JXF_Int                groupid_x,
                                jxf_ParCSRMatrix  **ARR_p, 
                                jxf_ParCSRMatrix   *AEE_p, 
                                jxf_ParCSRMatrix   *AII_p, 
                                jxf_ParVector     **VRE_p, 
                                jxf_ParVector     **VER_p, 
                                jxf_ParVector      *VEI_p, 
                                jxf_ParVector      *VIE_p, 
                                jxf_ParVector     **fR_p, 
                                jxf_ParVector      *fE_p, 
                                jxf_ParVector      *fI_p,
                                jxf_ParVector     **uR_p, 
                                jxf_ParVector      *uE_p, 
                                jxf_ParVector      *uI_p,
                                JXF_Int             ng,
                                // output:    
                                jxf_ParCSRMatrix  **ARR_ptr,
                                jxf_ParCSRMatrix  **AEE_ptr,
                                jxf_ParCSRMatrix  **AII_ptr, 
                                jxf_ParVector     **VRE_ptr, 
                                jxf_ParVector    ***VER_ptr, 
                                jxf_ParVector     **VEI_ptr, 
                                jxf_ParVector     **VIE_ptr,
                                jxf_ParCSRMatrix  **A_ptr,  
                                jxf_ParVector     **f_ptr, 
                                jxf_ParVector     **u_ptr );
JXF_Int 
jxf_MatVecGroup2All( MPI_Comm          comm, 
                    JXF_Int               groupid_x,
                    jxf_ParCSRMatrix  *ARR, 
                    jxf_ParCSRMatrix  *AEE, 
                    jxf_ParCSRMatrix  *AII, 
                    jxf_ParVector     *VRE, 
                    jxf_ParVector     *VER, 
                    jxf_ParVector     *VEI, 
                    jxf_ParVector     *VIE,
                    jxf_ParCSRMatrix **ARR_all_ptr, 
                    jxf_ParCSRMatrix **AEE_all_ptr, 
                    jxf_ParCSRMatrix **AII_all_ptr, 
                    jxf_ParVector    **VRE_all_ptr, 
                    jxf_ParVector    **VER_all_ptr, 
                    jxf_ParVector    **VEI_all_ptr, 
                    jxf_ParVector    **VIE_all_ptr );
JXF_Int
jxf_3tDataTransFromSeq2SubPara( JXF_Int               iniguess,
                               MPI_Comm          comm,
                               jxf_CSRMatrix     *A_s, 
                               jxf_Vector        *f_s, 
                               jxf_Vector        *u_s, 
                               jxf_ParCSRMatrix **ARR_p_ptr, 
                               jxf_ParCSRMatrix **AEE_p_ptr, 
                               jxf_ParCSRMatrix **AII_p_ptr, 
                               jxf_ParVector    **VRE_p_ptr, 
                               jxf_ParVector    **VER_p_ptr, 
                               jxf_ParVector    **VEI_p_ptr, 
                               jxf_ParVector    **VIE_p_ptr, 
                               jxf_ParVector    **fR_p_ptr, 
                               jxf_ParVector    **fE_p_ptr, 
                               jxf_ParVector    **fI_p_ptr,
                               jxf_ParVector    **uR_p_ptr, 
                               jxf_ParVector    **uE_p_ptr, 
                               jxf_ParVector    **uI_p_ptr );
JXF_Int
jxf_3tDataTransFromSeq2SubPar0( JXF_Int               iniguess,
                               MPI_Comm          comm,
                               jxf_CSRMatrix     *A_s, 
                               jxf_Vector        *f_s, 
                               jxf_Vector        *u_s, 
                               jxf_ParCSRMatrix **ARR_p_ptr, 
                               jxf_ParCSRMatrix **AEE_p_ptr, 
                               jxf_ParCSRMatrix **AII_p_ptr, 
                               jxf_ParVector    **VRE_p_ptr, 
                               jxf_ParVector    **VER_p_ptr, 
                               jxf_ParVector    **VEI_p_ptr, 
                               jxf_ParVector    **VIE_p_ptr, 
                               jxf_ParVector    **fR_p_ptr, 
                               jxf_ParVector    **fE_p_ptr, 
                               jxf_ParVector    **fI_p_ptr,
                               jxf_ParVector    **uR_p_ptr, 
                               jxf_ParVector    **uE_p_ptr, 
                               jxf_ParVector    **uI_p_ptr );                              
JXF_Int  
jxf_ParaDataTransSEQIF_mp( // input:
                          MPI_Comm           comm, 
                          MPI_Comm           comm_x, 
                          JXF_Int                groupid_x,
                          JXF_Real             theta_wc_E, 
                          JXF_Real             threshold_wc_E,
                          /* 串行离散系统数据 */
                          jxf_CSRMatrix      *A_s, 
                          jxf_Vector         *f_s, 
                          jxf_Vector         *u_s, 
                          // output: 
                          /* 系数矩阵各子块在整个进程组上的并行数据 */
                          jxf_ParCSRMatrix  **ARR_p_ptr, 
                          jxf_ParCSRMatrix  **AEE_p_ptr, 
                          jxf_ParCSRMatrix  **AII_p_ptr, 
                          jxf_ParVector     **VRE_p_ptr, 
                          jxf_ParVector     **VER_p_ptr, 
                          jxf_ParVector     **VEI_p_ptr, 
                          jxf_ParVector     **VIE_p_ptr, 
                          /* 系数矩阵各子块基于进程组的并行数据 */                         
                          jxf_ParCSRMatrix  **ARR_ptr,
                          jxf_ParCSRMatrix  **AEE_ptr,
                          jxf_ParCSRMatrix  **AII_ptr, 
                          jxf_ParVector     **VRE_ptr, 
                          jxf_ParVector     **VER_ptr, 
                          jxf_ParVector     **VEI_ptr, 
                          jxf_ParVector     **VIE_ptr,
                          /* 并行离散系统数据 */
                          jxf_ParCSRMatrix  **A_ptr,  
                          jxf_ParVector     **f_ptr, 
                          jxf_ParVector     **u_ptr,
                          /* 是否需要粗网格校正的标志变量 */
                          JXF_Int               *Need_CC_ptr );
JXF_Int  
jxf_ParaDataTransSEQIF_sp( // input:
                          MPI_Comm           comm, 
                          /* 串行离散系统数据 */
                          jxf_CSRMatrix      *A_s, 
                          jxf_Vector         *f_s, 
                          jxf_Vector         *u_s, 
                          // output:
                          /* 系数矩阵各子块基于进程组的并行数据 */                         
                          jxf_ParCSRMatrix  **ARR_ptr,
                          jxf_ParCSRMatrix  **AEE_ptr,
                          jxf_ParCSRMatrix  **AII_ptr, 
                          jxf_ParVector     **VRE_ptr, 
                          jxf_ParVector     **VER_ptr, 
                          jxf_ParVector     **VEI_ptr, 
                          jxf_ParVector     **VIE_ptr,
                          /* 并行离散系统数据 */
                          jxf_ParCSRMatrix  **A_ptr,  
                          jxf_ParVector     **f_ptr, 
                          jxf_ParVector     **u_ptr );
JXF_Int jxf_APCTLKrylovSolGather( jxf_ParVector *par_sol, jxf_Vector *u_s );
JXF_Int 
jxf_mgGenerateSubBlocks( MPI_Comm           comm,
                        MPI_Comm           comm_x,
                        JXF_Int             groupid_x,
                        JXF_Int             ng,
                        jxf_ParCSRMatrix   *par_mat,
                        jxf_ParCSRMatrix  **ARR_ptr,
                        jxf_ParCSRMatrix  **AEE_ptr,
                        jxf_ParCSRMatrix  **AII_ptr,
                        jxf_ParVector     **VRE_ptr,
                        jxf_ParVector    ***VER_ptr,
                        jxf_ParVector     **VEI_ptr,
                        jxf_ParVector     **VIE_ptr );
                          
                          

/* /csrc/3t/partition.c */
JXF_Int jxf_3tCheckNprocs( MPI_Comm comm, JXF_Int print_level, JXF_Int *np_R_ptr, JXF_Int *np_E_ptr, JXF_Int *np_I_ptr );
JXF_Int jxf_3tGetMyPartition( JXF_Int nprocs, JXF_Int np_R, JXF_Int np_E, JXF_Int np_I, JXF_Int N, JXF_Int **mypartition_ptr );   



/* /csrc/3t/apctl.c */
JXF_Int JXF_3tAPCTLSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetMaxiter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetNumRlxPre( JXF_Solver solver, JXF_Int num_relax );
JXF_Int JXF_3tAPCTLSetNumRlxPost( JXF_Solver solver, JXF_Int num_relax );
JXF_Int JXF_3tAPCTLSetPrintLevel( JXF_Solver solver, JXF_Int print_level );
JXF_Int JXF_3tAPCTLSetIsDiagElmFirst( JXF_Solver solver, JXF_Int is_diagelm_first );
JXF_Int JXF_3tAPCTLSetNumGroup( JXF_Solver solver, JXF_Int num_group );
JXF_Int JXF_3tAPCTLSetNumIterations( JXF_Solver solver, JXF_Int num_iterations );
JXF_Int JXF_3tAPCTLSetLastRelNrm( JXF_Solver solver, JXF_Real last_rel_nrm );
JXF_Int JXF_3tAPCTLSetAveConvFactor( JXF_Solver solver, JXF_Real ave_conv_factor );
JXF_Int JXF_3tAPCTLGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_3tAPCTLGetLastRelNrm( JXF_Solver solver, JXF_Real *rel_resid_norm );
JXF_Int JXF_3tAPCTLSetNpR( JXF_Solver solver, JXF_Int np );
JXF_Int JXF_3tAPCTLSetNpE( JXF_Solver solver, JXF_Int np );
JXF_Int JXF_3tAPCTLSetNpI( JXF_Solver solver, JXF_Int np );
JXF_Int JXF_3tAPCTLSetARRInterpMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetAIIInterpMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetARRRelaxMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetAEERelaxMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetAIIRelaxMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetACCRelaxMaxIt( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_3tAPCTLSetARRInterpTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetAIIInterpTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetARRRelaxTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetAEERelaxTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetAIIRelaxTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetACCRelaxTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_3tAPCTLSetARRRelaxType( JXF_Solver solver, JXF_Int ARR_relax_type );
JXF_Int JXF_3tAPCTLSetAEERelaxType( JXF_Solver solver, JXF_Int AEE_relax_type );
JXF_Int JXF_3tAPCTLSetAIIRelaxType( JXF_Solver solver, JXF_Int AII_relax_type );
JXF_Int JXF_3tAPCTLSetThetaWCR( JXF_Solver solver, JXF_Real theta );
JXF_Int JXF_3tAPCTLSetThetaWCE( JXF_Solver solver, JXF_Real theta ); 
JXF_Int JXF_3tAPCTLSetThetaWCI( JXF_Solver solver, JXF_Real theta ); 
JXF_Int JXF_3tAPCTLSetThresholdWCR( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_3tAPCTLSetThresholdWCE( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_3tAPCTLSetThresholdWCI( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_3tAPCTLSetISWCR( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetISWCE( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetISWCI( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetThetaDDR( JXF_Solver solver, JXF_Real theta ); 
JXF_Int JXF_3tAPCTLSetThetaDDE( JXF_Solver solver, JXF_Real theta );  
JXF_Int JXF_3tAPCTLSetThetaDDI( JXF_Solver solver, JXF_Real theta );  
JXF_Int JXF_3tAPCTLSetThresholdDDR( JXF_Solver solver, JXF_Real threshold ); 
JXF_Int JXF_3tAPCTLSetThresholdDDE( JXF_Solver solver, JXF_Real threshold );
JXF_Int JXF_3tAPCTLSetThresholdDDI( JXF_Solver solver, JXF_Real threshold );   
JXF_Int JXF_3tAPCTLSetISDDR( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetISDDE( JXF_Solver solver, JXF_Int flag ); 
JXF_Int JXF_3tAPCTLSetISDDI( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetNeedCC( JXF_Solver solver, JXF_Int flag );
JXF_Int JXF_3tAPCTLSetBlockSmoothType( JXF_Solver solver, JXF_Int blocksmooth_type );
JXF_Int JXF_3tAPCTLSetFixItPCTLR( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetFixItPCTLE( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetFixItPCTLI( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetFixItBRLXR( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetFixItBRLXE( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetFixItBRLXI( JXF_Solver solver, JXF_Int fixit );
JXF_Int JXF_3tAPCTLSetUseFixedModeR( JXF_Solver solver, JXF_Int usefixedmode ); /* peghoty  2012/02/15 */
JXF_Int JXF_3tAPCTLSetUseFixedModeE( JXF_Solver solver, JXF_Int usefixedmode ); /* peghoty  2012/02/15 */
JXF_Int JXF_3tAPCTLSetUseFixedModeI( JXF_Solver solver, JXF_Int usefixedmode ); /* peghoty  2012/02/15 */
JXF_Int JXF_3tAPCTLSetUsePPCTL( JXF_Solver solver, JXF_Int use_ppctl );   /* peghoty  2012/03/06 */
JXF_Int JXF_3tAPCTLSetARRSolverID( JXF_Solver solver, JXF_Int solverid ); /* peghoty  2011/10/29 */
JXF_Int JXF_3tAPCTLSetACCSolverID( JXF_Solver solver, JXF_Int solverid ); /* peghoty  2011/10/29 */
JXF_Int JXF_3tAPCTLSetARRKDim( JXF_Solver solver, JXF_Int kdim ); /* peghoty  2011/10/29 */
JXF_Int JXF_3tAPCTLSetACCKDim( JXF_Solver solver, JXF_Int kdim ); /* peghoty  2011/10/29 */
JXF_Int JXF_3tAPCTLSetTestSubLSIter( JXF_Solver solver, JXF_Int test_subls_iter ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetDebugFlag( JXF_Solver solver, JXF_Int debug_flag );
JXF_Int JXF_3tAPCTLSetResetZero( JXF_Solver solver, JXF_Int reset_zero );
JXF_Int JXF_3tAPCTLSetInterpType( JXF_Solver solver, JXF_Int interp_type );
JXF_Int JXF_3tAPCTLSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type );
JXF_Int JXF_3tAPCTLSetAggNumLevels( JXF_Solver solver, JXF_Int agg_num_levels );
JXF_Int JXF_3tAPCTLSetCoarseThreshold( JXF_Solver solver, JXF_Int coarse_threshold );
JXF_Int JXF_3tAPCTLSetPrintLevelAMG( JXF_Solver solver, JXF_Int print_level_amg );
JXF_Int JXF_3tAPCTLSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold );
JXF_Int JXF_3tAPCTLSetNumIterAiSetup( JXF_Solver solver, JXF_Int num_iter_Ai_pctl_setup ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetNumIterArSetup( JXF_Solver solver, JXF_Int num_iter_Ar_pctl_setup ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetNumIterAePrecond( JXF_Solver solver, JXF_Int num_iter_Ae_pctl_precond ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetNumIterAiPrecond( JXF_Solver solver, JXF_Int num_iter_Ai_pctl_precond ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetNumIterArPrecond( JXF_Solver solver, JXF_Int num_iter_Ar_pctl_precond ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetNumIterAcPrecond( JXF_Solver solver, JXF_Int num_iter_Ac_pctl_precond ); // peghoty, 2012/02/23
JXF_Int JXF_3tAPCTLSetA( JXF_Solver solver, jxf_ParCSRMatrix *A );
JXF_Int JXF_3tAPCTLSetSubBlocks( JXF_Solver solver );
JXF_Int JXF_3tAPCTLSetARR( JXF_Solver solver, jxf_ParCSRMatrix *ARR );
JXF_Int JXF_3tAPCTLSetAEE( JXF_Solver solver, jxf_ParCSRMatrix *AEE );
JXF_Int JXF_3tAPCTLSetAII( JXF_Solver solver, jxf_ParCSRMatrix *AII );
JXF_Int JXF_3tAPCTLSetVRE( JXF_Solver solver, jxf_ParVector *VRE );
JXF_Int JXF_3tAPCTLSetVER( JXF_Solver solver, jxf_ParVector *VER );
JXF_Int JXF_3tAPCTLSetVER2( JXF_Solver solver, jxf_ParVector **VER );
JXF_Int JXF_3tAPCTLSetVEI( JXF_Solver solver, jxf_ParVector *VEI );
JXF_Int JXF_3tAPCTLSetVIE( JXF_Solver solver, jxf_ParVector *VIE );
JXF_Int JXF_3tAPCTLSetARRAll( JXF_Solver solver, jxf_ParCSRMatrix *ARR );
JXF_Int JXF_3tAPCTLSetARRAll2( JXF_Solver solver, jxf_ParCSRMatrix **ARR );
JXF_Int JXF_3tAPCTLSetAEEAll( JXF_Solver solver, jxf_ParCSRMatrix *AEE );
JXF_Int JXF_3tAPCTLSetAIIAll( JXF_Solver solver, jxf_ParCSRMatrix *AII );
JXF_Int JXF_3tAPCTLSetVREAll( JXF_Solver solver, jxf_ParVector *VRE );
JXF_Int JXF_3tAPCTLSetVREAll2( JXF_Solver solver, jxf_ParVector **VRE );
JXF_Int JXF_3tAPCTLSetVERAll( JXF_Solver solver, jxf_ParVector *VER );
JXF_Int JXF_3tAPCTLSetVERAll2( JXF_Solver solver, jxf_ParVector **VER );
JXF_Int JXF_3tAPCTLSetVEIAll( JXF_Solver solver, jxf_ParVector *VEI );
JXF_Int JXF_3tAPCTLSetVIEAll( JXF_Solver solver, jxf_ParVector *VIE );
JXF_Int JXF_3tAPCTLIterCount( JXF_Solver solver );
JXF_Int JXF_3tAPCTLSetComm( JXF_Solver solver, MPI_Comm comm );
JXF_Int JXF_3tAPCTLSetCommX( JXF_Solver solver, MPI_Comm comm );
JXF_Int JXF_3tAPCTLSetCommY( JXF_Solver solver, MPI_Comm comm );
JXF_Int JXF_3tAPCTLSetGroupIdX( JXF_Solver solver, JXF_Int groupid );
JXF_Int JXF_3tAPCTLSetGroupIdY( JXF_Solver solver, JXF_Int groupid );
JXF_Int JXF_3tAPCTLDataInitialize( JXF_Solver *solver);
JXF_Int JXF_3tAPCTLDestroy( JXF_Solver solver );
JXF_Int jxf_3tAPCTLSetTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetMaxiter( void *solver, JXF_Int max_iter );       
JXF_Int jxf_3tAPCTLSetNumRlxPre( void *solver, JXF_Int num_relax );
JXF_Int jxf_3tAPCTLSetNumRlxPost( void *solver, JXF_Int num_relax );
JXF_Int jxf_3tAPCTLSetPrintLevel( void *solver, JXF_Int print_level );
JXF_Int jxf_3tAPCTLSetIsDiagElmFirst( void *solver, JXF_Int is_diagelm_first );
JXF_Int jxf_3tAPCTLSetNumGroup( void *solver, JXF_Int num_group );
JXF_Int jxf_3tAPCTLSetNumIterations( void *solver, JXF_Int num_iterations );
JXF_Int jxf_3tAPCTLSetLastRelNrm( void *solver, JXF_Real last_rel_nrm );
JXF_Int jxf_3tAPCTLSetAveConvFactor( void *solver, JXF_Real ave_conv_factor );
JXF_Int jxf_3tAPCTLGetNumIterations( void *solver, JXF_Int *num_iterations );
JXF_Int jxf_3tAPCTLGetLastRelNrm( void *solver, JXF_Real *rel_resid_norm );
JXF_Int jxf_3tAPCTLSetNpR( void *solver, JXF_Int np );
JXF_Int jxf_3tAPCTLSetNpE( void *solver, JXF_Int np );
JXF_Int jxf_3tAPCTLSetNpI( void *solver, JXF_Int np );
JXF_Int jxf_3tAPCTLSetARRInterpMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetAIIInterpMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetARRRelaxMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetAEERelaxMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetAIIRelaxMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetACCRelaxMaxIt( void *solver, JXF_Int max_iter );
JXF_Int jxf_3tAPCTLSetMaxItDefault( void *solver, JXF_Int maxit_default );
JXF_Int jxf_3tAPCTLSetARRInterpTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetAIIInterpTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetARRRelaxTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetAEERelaxTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetAIIRelaxTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetACCRelaxTol( void *solver, JXF_Real tol );
JXF_Int jxf_3tAPCTLSetARRRelaxType( void *solver, JXF_Int ARR_relax_type );
JXF_Int jxf_3tAPCTLSetAEERelaxType( void *solver, JXF_Int AEE_relax_type );
JXF_Int jxf_3tAPCTLSetAIIRelaxType( void *solver, JXF_Int AII_relax_type );
JXF_Int jxf_3tAPCTLSetTolDefault( void *solver, JXF_Real tol_default );
JXF_Int jxf_3tAPCTLSetThetaWCR( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThetaWCE( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThetaWCI( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThresholdWCR( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetThresholdWCE( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetThresholdWCI( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetISWCR( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetISWCE( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetISWCI( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetThetaDDR( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThetaDDE( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThetaDDI( void *solver, JXF_Real theta );
JXF_Int jxf_3tAPCTLSetThresholdDDR( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetThresholdDDE( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetThresholdDDI( void *solver, JXF_Real threshold );
JXF_Int jxf_3tAPCTLSetISDDR( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetISDDE( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetISDDI( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetNeedCC( void *solver, JXF_Int flag );
JXF_Int jxf_3tAPCTLSetBlockSmoothType( void *solver, JXF_Int blocksmooth_type );
JXF_Int jxf_3tAPCTLSetFixItPCTLR( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetFixItPCTLE( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetFixItPCTLI( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetFixItBRLXR( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetFixItBRLXE( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetFixItBRLXI( void *solver, JXF_Int fixit );
JXF_Int jxf_3tAPCTLSetUseFixedModeR( void *solver, JXF_Int usefixedmode ); // peghoty, 2012/02/15
JXF_Int jxf_3tAPCTLSetUseFixedModeE( void *solver, JXF_Int usefixedmode ); // peghoty, 2012/02/15
JXF_Int jxf_3tAPCTLSetUseFixedModeI( void *solver, JXF_Int usefixedmode ); // peghoty, 2012/02/15
JXF_Int jxf_3tAPCTLSetUsePPCTL( void *solver, JXF_Int use_ppctl ); // peghoty, 2012/03/06
JXF_Int jxf_APCTLKrylovParamSetDebugFlag( void *param, JXF_Int debug_flag );
JXF_Int jxf_APCTLKrylovParamSetResetZero( void *param, JXF_Int reset_zero );
JXF_Int jxf_APCTLKrylovParamSetStrongThreshold( void *param, JXF_Real strong_threshold );
JXF_Int jxf_APCTLKrylovParamSetInterpType( void *param, JXF_Int interp_type );
JXF_Int jxf_APCTLKrylovParamSetCoarsenType( void *param, JXF_Int coarsen_type );
JXF_Int jxf_APCTLKrylovParamSetAggNumLevels( void *param, JXF_Int agg_num_levels );
JXF_Int jxf_APCTLKrylovParamSetCoarseThreshold( void *param, JXF_Int coarse_threshold );
JXF_Int jxf_APCTLKrylovParamSetPrintLevelAMG( void *param, JXF_Int print_level_amg );
JXF_Int jxf_3tAPCTLSetARRSolverID( void *solver, JXF_Int solverid );
JXF_Int jxf_3tAPCTLSetACCSolverID( void *solver, JXF_Int solverid );
JXF_Int jxf_3tAPCTLSetARRKDim( void *solver, JXF_Int kdim );
JXF_Int jxf_3tAPCTLSetACCKDim( void *solver, JXF_Int kdim );
JXF_Int jxf_3tAPCTLSetTestSubLSIter( void *solver, JXF_Int test_subls_iter );  // peghoty, 2012/02/23           
JXF_Int jxf_3tAPCTLSetDebugFlag( void *solver, JXF_Int debug_flag );
JXF_Int jxf_3tAPCTLSetResetZero( void *solver, JXF_Int reset_zero );
JXF_Int jxf_3tAPCTLSetInterpType( void *solver, JXF_Int interp_type );
JXF_Int jxf_3tAPCTLSetCoarsenType( void *solver, JXF_Int coarsen_type );
JXF_Int jxf_3tAPCTLSetAggNumLevels( void *solver, JXF_Int agg_num_levels );
JXF_Int jxf_3tAPCTLSetCoarseThreshold( void *solver, JXF_Int coarse_threshold );
JXF_Int jxf_3tAPCTLSetPrintLevelAMG( void *solver, JXF_Int print_level_amg );
JXF_Int jxf_3tAPCTLSetStrongThreshold( void *solver, JXF_Real strong_threshold );
JXF_Int jxf_3tAPCTLSetNumIterAiSetup( void *solver, JXF_Int num_iter_Ai_pctl_setup ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetNumIterArSetup( void *solver, JXF_Int num_iter_Ar_pctl_setup ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetNumIterAePrecond( void *solver, JXF_Int num_iter_Ae_pctl_precond ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetNumIterAiPrecond( void *solver, JXF_Int num_iter_Ai_pctl_precond ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetNumIterArPrecond( void *solver, JXF_Int num_iter_Ar_pctl_precond ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetNumIterAcPrecond( void *solver, JXF_Int num_iter_Ac_pctl_precond ); // peghoty, 2012/02/23
JXF_Int jxf_3tAPCTLSetA( void *solver, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetSubBlocks( void *solver );
JXF_Int jxf_3tAPCTLSetARR( void *solver, jxf_ParCSRMatrix *ARR );
JXF_Int jxf_3tAPCTLSetAEE( void *solver, jxf_ParCSRMatrix *AEE );
JXF_Int jxf_3tAPCTLSetAII( void *solver, jxf_ParCSRMatrix *AII );
JXF_Int jxf_3tAPCTLSetVRE( void *solver, jxf_ParVector *VRE );
JXF_Int jxf_3tAPCTLSetVER( void *solver, jxf_ParVector *VER );
JXF_Int jxf_3tAPCTLSetVER2( void *solver, jxf_ParVector **VER );
JXF_Int jxf_3tAPCTLSetVEI( void *solver, jxf_ParVector *VEI );
JXF_Int jxf_3tAPCTLSetVIE( void *solver, jxf_ParVector *VIE );
JXF_Int jxf_3tAPCTLSetARRAll( void *solver, jxf_ParCSRMatrix *ARR );
JXF_Int jxf_3tAPCTLSetARRAll2( void *solver, jxf_ParCSRMatrix **ARR );
JXF_Int jxf_3tAPCTLSetAEEAll( void *solver, jxf_ParCSRMatrix *AEE );
JXF_Int jxf_3tAPCTLSetAIIAll( void *solver, jxf_ParCSRMatrix *AII );
JXF_Int jxf_3tAPCTLSetVREAll( void *solver, jxf_ParVector *VRE );
JXF_Int jxf_3tAPCTLSetVREAll2( void *solver, jxf_ParVector **VRE );
JXF_Int jxf_3tAPCTLSetVERAll( void *solver, jxf_ParVector *VER );
JXF_Int jxf_3tAPCTLSetVERAll2( void *solver, jxf_ParVector **VER );
JXF_Int jxf_3tAPCTLSetVEIAll( void *solver, jxf_ParVector *VEI );
JXF_Int jxf_3tAPCTLSetVIEAll( void *solver, jxf_ParVector *VIE );
JXF_Int jxf_3tAPCTLSetComm( void *solver, MPI_Comm comm );
JXF_Int jxf_3tAPCTLSetCommX( void *solver, MPI_Comm comm );
JXF_Int jxf_3tAPCTLSetCommY( void *solver, MPI_Comm comm );
JXF_Int jxf_3tAPCTLSetGroupIdX( void *solver, JXF_Int groupid );
JXF_Int jxf_3tAPCTLSetGroupIdY( void *solver, JXF_Int groupid );
void *jxf_3tAPCTLDataInitialize();                  
JXF_Int jxf_3tAPCTLDestroy( void *vdata );
JXF_Int jxf_3tAPCTLDestroy_sp( void *vdata );           
JXF_Int jxf_3tAPCTLDestroy_mp( void *vdata ); 
JXF_Int jxf_3tAPCTLIterCount( void *solver );
JXF_Int jxf_3tAPCTLIterCount_sp( void *solver );
JXF_Int jxf_3tAPCTLIterCount_mp( void *solver );
JXF_Int jxf_GetAPCTLNumIterOfSubLS( JXF_Solver precond, jxf_APCTLKrylovParam *apctlkrylov_param ); 
JXF_Int jxf_GetAPCTLmgNumIterOfSubLS( JXF_Solver precond, JXF_Int groupid_x, JXF_Int ng, jxf_APCTLKrylovParam *apctlkrylov_param ); 



/* /csrc/3t/setup.c */
JXF_Int JXF_3tAPCTLSetup( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int jxf_3tAPCTLSetup( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A ); 
JXF_Int jxf_3tAPCTLSetup_sp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );                
JXF_Int jxf_3tAPCTLSetup_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );                    

                              

/* /csrc/3t/solve.c */
JXF_Int 
JXF_3tAPCTLSolve( JXF_Solver        solver,
                 JXF_ParCSRMatrix  A,
                 JXF_ParVector     f,
                 JXF_ParVector     u  );
JXF_Int 
JXF_3tAPCTLPrecond( JXF_Solver        solver,
                   JXF_ParCSRMatrix  A,
                   JXF_ParVector     b,
                   JXF_ParVector     x  );
JXF_Int 
JXF_3tAPCTLmgPrecond( JXF_Solver        solver,
                     JXF_ParCSRMatrix  A,
                     JXF_ParVector     b,
                     JXF_ParVector     x  );
JXF_Int 
JXF_3tABSC1mgPrecond( JXF_Solver        solver,
                     JXF_ParCSRMatrix  A,
                     JXF_ParVector     b,
                     JXF_ParVector     x  );
JXF_Int 
JXF_3tABSC2mgPrecond( JXF_Solver        solver,
                     JXF_ParCSRMatrix  A,
                     JXF_ParVector     b,
                     JXF_ParVector     x  );
JXF_Int
jxf_3tAPCTLOneIteration( jxf_3tAPCTLData    *pre_3tapctl_data,
                        jxf_ParCSRMatrix   *A,
                        jxf_ParVector      *g,
                        jxf_ParVector      *w ); 
JXF_Int
jxf_3tAPCTLmgOneIteration( jxf_3tAPCTLData    *pre_3tapctl_data,
                          jxf_ParCSRMatrix   *A,
                          jxf_ParVector      *g,
                          jxf_ParVector      *w );
JXF_Int
jxf_3tABSC1mgOneIteration( jxf_3tAPCTLData    *pre_3tapctl_data,
                          jxf_ParCSRMatrix   *A,
                          jxf_ParVector      *g,
                          jxf_ParVector      *w );
JXF_Int
jxf_3tABSC2mgOneIteration( jxf_3tAPCTLData    *pre_3tapctl_data,
                          jxf_ParCSRMatrix   *A,
                          jxf_ParVector      *g,
                          jxf_ParVector      *w );
JXF_Int
jxf_3tAPCTLOneIteration_mp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w );  
JXF_Int
jxf_3tAPCTLmgOneIteration_mp( jxf_3tAPCTLData  *pre_3tapctl_data,
                             jxf_ParCSRMatrix *A,
                             jxf_ParVector    *g,
                             jxf_ParVector    *w );
JXF_Int
jxf_3tABSC1mgOneIteration_mp( jxf_3tAPCTLData  *pre_3tapctl_data,
                             jxf_ParCSRMatrix *A,
                             jxf_ParVector    *g,
                             jxf_ParVector    *w );
JXF_Int
jxf_3tABSC2mgOneIteration_mp( jxf_3tAPCTLData  *pre_3tapctl_data,
                             jxf_ParCSRMatrix *A,
                             jxf_ParVector    *g,
                             jxf_ParVector    *w );
JXF_Int
jxf_3tAPCTLOneIteration_sp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w );                                                  
JXF_Int
jxf_3tAPCTLPrecond( jxf_3tAPCTLData    *pre_3tapctl_data,
                   jxf_ParCSRMatrix   *A,
                   jxf_ParVector      *f,
                   jxf_ParVector      *u );                                        
JXF_Int
jxf_3tAPCTLmgPrecond( jxf_3tAPCTLData    *pre_3tapctl_data,
                     jxf_ParCSRMatrix   *A,
                     jxf_ParVector      *f,
                     jxf_ParVector      *u );
JXF_Int
jxf_3tABSC1mgPrecond( jxf_3tAPCTLData    *pre_3tapctl_data,
                     jxf_ParCSRMatrix   *A,
                     jxf_ParVector      *f,
                     jxf_ParVector      *u );
JXF_Int
jxf_3tABSC2mgPrecond( jxf_3tAPCTLData    *pre_3tapctl_data,
                     jxf_ParCSRMatrix   *A,
                     jxf_ParVector      *f,
                     jxf_ParVector      *u );
JXF_Int
jxf_3tAPCTLSolve( jxf_3tAPCTLData    *pre_3tapctl_data,
                 jxf_ParCSRMatrix   *A,
                 jxf_ParVector      *f,
                 jxf_ParVector      *u );



/* /csrc/3t/relax.c */ 
JXF_Int
jxf_3tAPCTLRelax_GSType_mp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w );
JXF_Int
jxf_3tAPCTLmgRelax_GSType_mp( jxf_3tAPCTLData    *pre_3tapctl_data,
                             jxf_ParCSRMatrix   *A,
                             jxf_ParVector      *g,
                             jxf_ParVector      *w );
JXF_Int
jxf_3tAPCTLRelax_BDType_mp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w ); 
JXF_Int
jxf_3tAPCTLmgRelax_BDType_mp( jxf_3tAPCTLData    *pre_3tapctl_data,
                             jxf_ParCSRMatrix   *A,
                             jxf_ParVector      *g,
                             jxf_ParVector      *w );
JXF_Int
jxf_3tAPCTLRelax_GSType_sp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w );
JXF_Int
jxf_3tAPCTLRelax_BDType_sp( jxf_3tAPCTLData    *pre_3tapctl_data,
                           jxf_ParCSRMatrix   *A,
                           jxf_ParVector      *g,
                           jxf_ParVector      *w );         


/* /csrc/3t/param.c */                         
jxf_APCTLKrylovParam *jxf_APCTLKrylovParamInitialize();
JXF_Int jxf_APCTLKrylovParamSetNumGroup( void *param, JXF_Int ng );
JXF_Int jxf_APCTLKrylovParamSetSolverID( void *param, JXF_Int solver_id );
JXF_Int jxf_APCTLKrylovParamSetPrecondID( void *param, JXF_Int precond_id );
JXF_Int jxf_APCTLKrylovParamSetTol( void *param, JXF_Real tol );
JXF_Int jxf_APCTLKrylovParamSetMaxIter( void *param, JXF_Int max_iter );
JXF_Int jxf_APCTLKrylovParamSetKDim( void *param, JXF_Int k_dim );
JXF_Int jxf_APCTLKrylovParamSetISCheckRestarted( void *param, JXF_Int is_check_restarted );
JXF_Int jxf_APCTLKrylovParamSetTwoNorm( void *param, JXF_Int two_norm );
JXF_Int jxf_APCTLKrylovParamSetPrintLevel( void *param, JXF_Int print_level );
JXF_Int jxf_APCTLKrylovParamSetTTest( void *param, JXF_Int TTest );
JXF_Int jxf_APCTLKrylovParamSetKeepSol( void *param, JXF_Int keepsol );
JXF_Int jxf_APCTLKrylovParamSetNumIterations( void *param, JXF_Int num );
JXF_Int jxf_APCTLKrylovParamSetLastRelNrm( void *param, JXF_Real last_rel_nrm );
JXF_Int jxf_APCTLKrylovParamSetPrintLevelAPCTL( void *param, JXF_Int print_level_apctl );
JXF_Int jxf_APCTLKrylovParamSetBlockSmoothType( void *param, JXF_Int blocksmooth_type );
JXF_Int jxf_APCTLKrylovParamSetNumRelaxPre( void *param, JXF_Int num_relax_pre );
JXF_Int jxf_APCTLKrylovParamSetNumRelaxPost( void *param, JXF_Int num_relax_post );
JXF_Int jxf_APCTLKrylovParamSetInterpSolverARR( void *param, JXF_Int interp_solver_ARR );
JXF_Int jxf_APCTLKrylovParamSetInterpKdimARR( void *param, JXF_Int interp_kdim_ARR );
JXF_Int jxf_APCTLKrylovParamSetInterpMaxitARR( void *param, JXF_Int interp_maxit_ARR );
JXF_Int jxf_APCTLKrylovParamSetInterpMaxitAII( void *param, JXF_Int interp_maxit_AII );
JXF_Int jxf_APCTLKrylovParamSetInterpTolARR( void *param, JXF_Real interp_tol_ARR );
JXF_Int jxf_APCTLKrylovParamSetInterpTolAII( void *param, JXF_Real interp_tol_AII );
JXF_Int jxf_APCTLKrylovParamSetARRRelaxTol( void *param, JXF_Real ARR_relax_tol );
JXF_Int jxf_APCTLKrylovParamSetAEERelaxTol( void *param, JXF_Real AEE_relax_tol );
JXF_Int jxf_APCTLKrylovParamSetAIIRelaxTol( void *param, JXF_Real AII_relax_tol );
JXF_Int jxf_APCTLKrylovParamSetACCRelaxTol( void *param, JXF_Real ACC_relax_tol );
JXF_Int jxf_APCTLKrylovParamSetFixitPCTLR( void *param, JXF_Int fixit_pctl_R );
JXF_Int jxf_APCTLKrylovParamSetFixitPCTLE( void *param, JXF_Int fixit_pctl_E );
JXF_Int jxf_APCTLKrylovParamSetFixitPCTLI( void *param, JXF_Int fixit_pctl_I );
JXF_Int jxf_APCTLKrylovParamSetFixitBrlxR( void *param, JXF_Int fixit_brlx_R );
JXF_Int jxf_APCTLKrylovParamSetFixitBrlxE( void *param, JXF_Int fixit_brlx_E );
JXF_Int jxf_APCTLKrylovParamSetFixitBrlxI( void *param, JXF_Int fixit_brlx_I );
JXF_Int jxf_APCTLKrylovParamSetARRRelaxType( void *param, JXF_Int ARR_relax_type );
JXF_Int jxf_APCTLKrylovParamSetAEERelaxType( void *param, JXF_Int AEE_relax_type );
JXF_Int jxf_APCTLKrylovParamSetAIIRelaxType( void *param, JXF_Int AII_relax_type );
JXF_Int jxf_APCTLKrylovParamSetUseFixedModeR( void *param, JXF_Int use_fixedmode_R );
JXF_Int jxf_APCTLKrylovParamSetUseFixedModeE( void *param, JXF_Int use_fixedmode_E );
JXF_Int jxf_APCTLKrylovParamSetUseFixedModeI( void *param, JXF_Int use_fixedmode_I );
JXF_Int jxf_APCTLKrylovParamSetThetaWCE( void *param, JXF_Real theta_wc_E );
JXF_Int jxf_APCTLKrylovParamSetParamThresholdWCE( void *param, JXF_Real threshold_wc_E );
JXF_Int jxf_APCTLKrylovParamSetThetaDDR( void *param, JXF_Real theta_dd_R );
JXF_Int jxf_APCTLKrylovParamSetThetaDDE( void *param, JXF_Real theta_dd_E );
JXF_Int jxf_APCTLKrylovParamSetThetaDDI( void *param, JXF_Real theta_dd_I );
JXF_Int jxf_APCTLKrylovParamSetThresholdDDR( void *param, JXF_Real threshold_dd_R );
JXF_Int jxf_APCTLKrylovParamSetThresholdDDE( void *param, JXF_Real threshold_dd_E );
JXF_Int jxf_APCTLKrylovParamSetThresholdDDI( void *param, JXF_Real threshold_dd_I );
JXF_Int jxf_APCTLKrylovParamSetISDiagElmFirst( void *param, JXF_Int is_diagelm_first );
JXF_Int jxf_APCTLKrylovParamSetUsePPCTL( void *param, JXF_Int use_ppctl );
JXF_Int jxf_APCTLKrylovParamSetTestSubLSIter( void *param, JXF_Int test_subls_iter );
JXF_Int jxf_APCTLKrylovParamSetCPUTrans( void *param, JXF_Real cpu_trans );
JXF_Int jxf_APCTLKrylovParamSetCPUSetup( void *param, JXF_Real cpu_setup );
JXF_Int jxf_APCTLKrylovParamSetCPUSolve( void *param, JXF_Real cpu_solve );
JXF_Int jxf_APCTLKrylovParamSetCPUTotal( void *param, JXF_Real cpu_total );



/* /csrc/3t/jasminIF.c */
JXF_Int 
jxf_ApctlKrylov_JASMIN( jxf_ParCSRMatrix     *ARR_p, 
                       jxf_ParCSRMatrix     *AEE_p, 
                       jxf_ParCSRMatrix     *AII_p, 
                       jxf_ParVector        *VRE_p, 
                       jxf_ParVector        *VER_p, 
                       jxf_ParVector        *VEI_p, 
                       jxf_ParVector        *VIE_p, 
                       jxf_ParVector        *fR_p, 
                       jxf_ParVector        *fE_p, 
                       jxf_ParVector        *fI_p,
                       jxf_ParVector        *uR_p, 
                       jxf_ParVector        *uE_p, 
                       jxf_ParVector        *uI_p,                      
                       jxf_APCTLKrylovParam *apctlkrylov_param );
JXF_Int
jxf_ApctlKrylov_mgJASMIN( jxf_ParCSRMatrix    **ARR_p,
                         jxf_ParCSRMatrix     *AEE_p,
                         jxf_ParCSRMatrix     *AII_p,
                         jxf_ParVector       **VRE_p,
                         jxf_ParVector       **VER_p,
                         jxf_ParVector        *VEI_p,
                         jxf_ParVector        *VIE_p,
                         jxf_ParVector       **fR_p,
                         jxf_ParVector        *fE_p,
                         jxf_ParVector        *fI_p,
                         jxf_ParVector       **uR_p,
                         jxf_ParVector        *uE_p,
                         jxf_ParVector        *uI_p,
                         jxf_APCTLKrylovParam *apctlkrylov_param );
JXF_Int 
jxf_ApctlKrylov_JASMIN_mp( jxf_ParCSRMatrix     *ARR_p, 
                          jxf_ParCSRMatrix     *AEE_p, 
                          jxf_ParCSRMatrix     *AII_p, 
                          jxf_ParVector        *VRE_p, 
                          jxf_ParVector        *VER_p, 
                          jxf_ParVector        *VEI_p, 
                          jxf_ParVector        *VIE_p, 
                          jxf_ParVector        *fR_p, 
                          jxf_ParVector        *fE_p, 
                          jxf_ParVector        *fI_p,
                          jxf_ParVector        *uR_p, 
                          jxf_ParVector        *uE_p, 
                          jxf_ParVector        *uI_p,
                          jxf_APCTLKrylovParam *apctlkrylov_param ); 
JXF_Int
jxf_ApctlKrylov_mgJASMIN_mp( jxf_ParCSRMatrix    **ARR_p, 
                            jxf_ParCSRMatrix     *AEE_p, 
                            jxf_ParCSRMatrix     *AII_p, 
                            jxf_ParVector       **VRE_p, 
                            jxf_ParVector       **VER_p, 
                            jxf_ParVector        *VEI_p, 
                            jxf_ParVector        *VIE_p, 
                            jxf_ParVector       **fR_p, 
                            jxf_ParVector        *fE_p, 
                            jxf_ParVector        *fI_p,
                            jxf_ParVector       **uR_p, 
                            jxf_ParVector        *uE_p, 
                            jxf_ParVector        *uI_p,
                            jxf_APCTLKrylovParam *apctlkrylov_param );
JXF_Int 
jxf_ApctlKrylov_JASMIN_sp( jxf_ParCSRMatrix     *ARR_p, 
                          jxf_ParCSRMatrix     *AEE_p, 
                          jxf_ParCSRMatrix     *AII_p, 
                          jxf_ParVector        *VRE_p, 
                          jxf_ParVector        *VER_p, 
                          jxf_ParVector        *VEI_p, 
                          jxf_ParVector        *VIE_p, 
                          jxf_ParVector        *fR_p, 
                          jxf_ParVector        *fE_p, 
                          jxf_ParVector        *fI_p,
                          jxf_ParVector        *uR_p, 
                          jxf_ParVector        *uE_p, 
                          jxf_ParVector        *uI_p,
                          jxf_APCTLKrylovParam *apctlkrylov_param );                       
JXF_Int JXF_3tAPCTLSetup4Jasmin( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int JXF_3tAPCTLSetup4mgJasmin( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int JXF_3tABSC1Setup4mgJasmin( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int JXF_3tABSC2Setup4mgJasmin( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int jxf_3tAPCTLSetup4Jasmin( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetup4mgJasmin( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tABSC1Setup4mgJasmin( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tABSC2Setup4mgJasmin( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetup4Jasmin_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetup4mgJasmin_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tABSC1Setup4mgJasmin_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tABSC2Setup4mgJasmin_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetup4Jasmin_sp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int JXF_3tAPCTLDestroy4Jasmin( JXF_Solver solver );  
JXF_Int JXF_3tAPCTLDestroy4mgJasmin( JXF_Solver solver );
JXF_Int jxf_3tAPCTLDestroy4Jasmin( void *vdata );
JXF_Int jxf_3tAPCTLDestroy4mgJasmin( void *vdata );
JXF_Int jxf_3tAPCTLDestroy4Jasmin_sp( void *vdata );
JXF_Int jxf_3tAPCTLDestroy4Jasmin_mp( void *vdata );
JXF_Int jxf_3tAPCTLDestroy4mgJasmin_mp( void *vdata );



/* /csrc/3t/interp.c */
jxf_ParCSRMatrix *
jxf_3tApctlInterpolation( MPI_Comm         comm, 
                         MPI_Comm         comm_x,
                         JXF_Int              N,
                         JXF_Int              groupid_x,
                         JXF_Int             *row_starts, 
                         JXF_Int             *col_starts,
                         jxf_ParVector    *WEE );
jxf_ParCSRMatrix *
jxf_3tApctlmgInterpolation( MPI_Comm         comm, 
                           MPI_Comm         comm_x,
                           JXF_Int              N,
                           JXF_Int              groupid_x,
                           JXF_Int              ng,
                           JXF_Int             *row_starts, 
                           JXF_Int             *col_starts,
                           jxf_ParVector    *WEE );
JXF_Int 
jxf_3tAPCTLParaInterpVec( MPI_Comm       comm,
                         MPI_Comm       comm_x, 
                         JXF_Int            groupid_x,
                         JXF_Int           *vec_starts, 
                         jxf_ParVector  *WEE, 
                         jxf_ParVector **PRR_ptr, 
                         jxf_ParVector **PII_ptr );
JXF_Int 
jxf_3tAPCTLmgParaInterpVec( MPI_Comm        comm,
                           MPI_Comm        comm_x, 
                           JXF_Int            groupid_x,
                           JXF_Int           *vec_starts, 
                           JXF_Int          ng,
                           jxf_ParVector   *WEE, 
                           jxf_ParVector ***PRR_ptr, 
                           jxf_ParVector  **PII_ptr );


/* /csrc/3t/coarseop.c */
jxf_ParCSRMatrix *
jxf_3tApctlCoarseOperator_sp( MPI_Comm       comm,
                             jxf_CSRMatrix  *ARR_s, 
                             jxf_CSRMatrix  *AEE_s, 
                             jxf_CSRMatrix  *AII_s, 
                             jxf_Vector     *VRE_s, 
                             jxf_Vector     *VER_s, 
                             jxf_Vector     *VEI_s, 
                             jxf_Vector     *VIE_s,
                             jxf_ParVector  *PRR,
                             jxf_ParVector  *PII  );
jxf_ParCSRMatrix *
jxf_3tApctlCoarseOperator_mp( MPI_Comm          comm,
                             jxf_ParCSRMatrix  *ARR, 
                             jxf_ParCSRMatrix  *AEE, 
                             jxf_ParCSRMatrix  *AII, 
                             jxf_ParVector     *VRE, 
                             jxf_ParVector     *VER, 
                             jxf_ParVector     *VEI, 
                             jxf_ParVector     *VIE,
                             jxf_ParVector     *PRR, 
                             jxf_ParVector     *PII ); 
jxf_ParCSRMatrix *
jxf_3tApctlmgCoarseOperator_mp( MPI_Comm          comm,
                               JXF_Int            ng,
                               jxf_ParCSRMatrix **ARR, 
                               jxf_ParCSRMatrix  *AEE, 
                               jxf_ParCSRMatrix  *AII, 
                               jxf_ParVector    **VRE, 
                               jxf_ParVector    **VER, 
                               jxf_ParVector     *VEI, 
                               jxf_ParVector     *VIE,
                               jxf_ParVector    **PRR, 
                               jxf_ParVector     *PII );


/* /csrc/3t/seqIF.c */
JXF_Int 
jxf_ApctlKrylov_SEQ( MPI_Comm              comm,
                    jxf_CSRMatrix         *A_s,  
                    jxf_Vector            *f_s, 
                    jxf_Vector            *u_s,                      
                    jxf_APCTLKrylovParam  *apctlkrylov_param );
JXF_Int 
jxf_ApctlKrylov_SEQ_mp( MPI_Comm              comm,
                       jxf_CSRMatrix         *A_s,  
                       jxf_Vector            *f_s, 
                       jxf_Vector            *u_s,                      
                       jxf_APCTLKrylovParam  *apctlkrylov_param );
JXF_Int 
jxf_ApctlKrylov_SEQ_sp( MPI_Comm              comm,
                       jxf_CSRMatrix         *A_s,  
                       jxf_Vector            *f_s, 
                       jxf_Vector            *u_s,                      
                       jxf_APCTLKrylovParam  *apctlkrylov_param );                    
JXF_Int JXF_3tAPCTLSetupSEQIF( JXF_Solver solver, JXF_ParCSRMatrix A );
JXF_Int jxf_3tAPCTLSetupSEQIF( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetupSEQIF_mp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int jxf_3tAPCTLSetupSEQIF_sp( jxf_3tAPCTLData *pre_3tapctl_data, jxf_ParCSRMatrix *A );
JXF_Int JXF_3tAPCTLDestroySEQIF( JXF_Solver solver );
JXF_Int jxf_3tAPCTLDestroySEQIF( void *vdata );
JXF_Int jxf_3tAPCTLDestroySEQIF_mp( void *vdata );
JXF_Int jxf_3tAPCTLDestroySEQIF_sp( void *vdata );

                             
#endif                                                     
