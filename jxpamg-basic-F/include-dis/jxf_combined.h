//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_combined.h -- head file for combined preconditioners
 *  Date: 2013/12/11
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JXF_COMBINED_HEADER
#define JXF_COMBINED_HEADER

#ifndef JXF_PAMG_HEADER
#include "jxf_pamg.h"
#endif

#ifndef JXF_EUCLID_HEADER
#include "jxf_euclid.h"
#endif

#ifndef JXF_ILU_HEADER
#include "jxf_ilu.h"
#endif

#define JXF_CP_EUCLID_ERRCHKA \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("", __FUNC__, __FILE__, __LINE__); \
        jxf_printErrorMsg(stderr); \
        jxf_MPI_Abort(jxf_comm_dh, -1); \
      }

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_CombinedPrecondData
 */
typedef struct
{
    JXF_Int nx;
    JXF_Int ny;
    JXF_Int npx;
    JXF_Int npy;
    JXF_Int num_equns;
    
    JXF_Int pre_id;
    JXF_Int theta_psi;
    JXF_Int theta_rho;
    JXF_Int theta_phi;
    JXF_Int actual_psi;
    JXF_Int actual_rho;
    JXF_Int actual_phi;
    JXF_Int euclid_level;
    JXF_Int amg_max_levels;
    JXF_Int amg_relax_type;
    JXF_Int amg_print_level;
    JXF_Int amg_interp_type;
    JXF_Int amg_P_max_elmts;
    JXF_Int amg_measure_type;
    JXF_Int amg_coarsen_type;
    JXF_Int amg_agg_num_levels;
    JXF_Int amg_coarse_threshold;
    
    JXF_Real drop_tol;
    JXF_Real adp_theta;
    JXF_Real theta_dis;
    JXF_Real adp_gamma_3;
    JXF_Real adp_gamma_11;
    JXF_Real amg_strong_threshold;
    
    MPI_Comm comm;
    
    jxf_ParCSRMatrix *pre_mat;
    jxf_ParCSRMatrix *pre_euc_mat;
    
    jxf_ParVector *aux_vec;
    jxf_ParVector *res_vec;
    
    JXF_Solver ilu_data;
    JXF_Solver pamg_data;
    JXF_Solver euclid_data;
    
} jxf_CombinedPrecondData;

#define jxf_CombinedPrecondDataNx(pre_data)                 ((pre_data)->nx)
#define jxf_CombinedPrecondDataNy(pre_data)                 ((pre_data)->ny)
#define jxf_CombinedPrecondDataNpx(pre_data)                ((pre_data)->npx)
#define jxf_CombinedPrecondDataNpy(pre_data)                ((pre_data)->npy)
#define jxf_CombinedPrecondDataNumEquns(pre_data)           ((pre_data)->num_equns)
#define jxf_CombinedPrecondDataPreID(pre_data)              ((pre_data)->pre_id)
#define jxf_CombinedPrecondDataThetaPsi(pre_data)           ((pre_data)->theta_psi)
#define jxf_CombinedPrecondDataThetaRho(pre_data)           ((pre_data)->theta_rho)
#define jxf_CombinedPrecondDataThetaPhi(pre_data)           ((pre_data)->theta_phi)
#define jxf_CombinedPrecondDataActualPsi(pre_data)          ((pre_data)->actual_psi)
#define jxf_CombinedPrecondDataActualRho(pre_data)          ((pre_data)->actual_rho)
#define jxf_CombinedPrecondDataActualPhi(pre_data)          ((pre_data)->actual_phi)
#define jxf_CombinedPrecondDataEuclidLevel(pre_data)        ((pre_data)->euclid_level)
#define jxf_CombinedPrecondDataAMGMaxLevels(pre_data)       ((pre_data)->amg_max_levels)
#define jxf_CombinedPrecondDataAMGRelaxType(pre_data)       ((pre_data)->amg_relax_type)
#define jxf_CombinedPrecondDataAMGPrintLevel(pre_data)      ((pre_data)->amg_print_level)
#define jxf_CombinedPrecondDataAMGInterpType(pre_data)      ((pre_data)->amg_interp_type)
#define jxf_CombinedPrecondDataAMGPMaxElmts(pre_data)       ((pre_data)->amg_P_max_elmts)
#define jxf_CombinedPrecondDataAMGMeasureType(pre_data)     ((pre_data)->amg_measure_type)
#define jxf_CombinedPrecondDataAMGCoarsenType(pre_data)     ((pre_data)->amg_coarsen_type)
#define jxf_CombinedPrecondDataAMGAggNumLevels(pre_data)    ((pre_data)->amg_agg_num_levels)
#define jxf_CombinedPrecondDataAMGCoarseThreshold(pre_data) ((pre_data)->amg_coarse_threshold)
#define jxf_CombinedPrecondDataDropTol(pre_data)            ((pre_data)->drop_tol)
#define jxf_CombinedPrecondDataAdpTheta(pre_data)           ((pre_data)->adp_theta)
#define jxf_CombinedPrecondDataThetaDis(pre_data)           ((pre_data)->theta_dis)
#define jxf_CombinedPrecondDataAdpGammaT(pre_data)          ((pre_data)->adp_gamma_3)
#define jxf_CombinedPrecondDataAdpGammaE(pre_data)          ((pre_data)->adp_gamma_11)
#define jxf_CombinedPrecondDataAMGStrongThreshold(pre_data) ((pre_data)->amg_strong_threshold)
#define jxf_CombinedPrecondDataComm(pre_data)               ((pre_data)->comm)
#define jxf_CombinedPrecondDataPreMat(pre_data)             ((pre_data)->pre_mat)
#define jxf_CombinedPrecondDataPreEucMat(pre_data)          ((pre_data)->pre_euc_mat)
#define jxf_CombinedPrecondDataAuxVector(pre_data)          ((pre_data)->aux_vec)
#define jxf_CombinedPrecondDataResVector(pre_data)          ((pre_data)->res_vec)
#define jxf_CombinedPrecondDataILUZData(pre_data)           ((pre_data)->ilu_data)
#define jxf_CombinedPrecondDataPAMGData(pre_data)           ((pre_data)->pamg_data)
#define jxf_CombinedPrecondDataEuclidData(pre_data)         ((pre_data)->euclid_data)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/combined/amg_euclid.c */
JXF_Int
JXF_PAMGEuclidPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_PAMGEuclidPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_PAMGEuclidPAMGSetupB( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_PAMGEuclidPAMGSetupB( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_EuclidPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_EuclidPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_PAMGEuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_PAMGEuclidSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_EuclidPAMGEuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_EuclidPAMGEuclidSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_PAMGEuclidPAMGSolve( JXF_Solver       solver,
                        JXF_hpCSRMatrix par_matrix,
                        JXF_ParVector    par_b,
                        JXF_ParVector    par_x );
JXF_Int
jxf_PAMGEuclidPAMGSolve( void            *data,
                        jxf_hpCSRMatrix *par_matrix,
                        jxf_ParVector    *par_b,
                        jxf_ParVector    *par_x );
JXF_Int
JXF_EuclidPAMGSolve( JXF_Solver       solver,
                    JXF_hpCSRMatrix par_matrix,
                    JXF_ParVector    par_b,
                    JXF_ParVector    par_x );
JXF_Int
jxf_EuclidPAMGSolve( void            *data,
                    jxf_hpCSRMatrix *par_matrix,
                    jxf_ParVector    *par_b,
                    jxf_ParVector    *par_x );
JXF_Int
JXF_PAMGEuclidSolve( JXF_Solver       solver,
                    JXF_hpCSRMatrix par_matrix,
                    JXF_ParVector    par_b,
                    JXF_ParVector    par_x );
JXF_Int
jxf_PAMGEuclidSolve( void            *data,
                    jxf_hpCSRMatrix *par_matrix,
                    jxf_ParVector    *par_b,
                    jxf_ParVector    *par_x );
JXF_Int
JXF_EuclidPAMGEuclidSolve( JXF_Solver       solver,
                          JXF_hpCSRMatrix par_matrix,
                          JXF_ParVector    par_b,
                          JXF_ParVector    par_x );
JXF_Int
jxf_EuclidPAMGEuclidSolve( void            *data,
                          jxf_hpCSRMatrix *par_matrix,
                          jxf_ParVector    *par_b,
                          jxf_ParVector    *par_x );

/* csrc/combined/amg_ilu.c */
JXF_Int
JXF_ILUZPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
JXF_ILUZPAMGILUZSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix );
JXF_Int
jxf_ILUZPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
jxf_ILUZPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
jxf_ILUZPAMGILUZSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
jxf_PAMGILUZPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix );
JXF_Int
JXF_ILUZPAMGSolve( JXF_Solver       solver,
                  JXF_hpCSRMatrix par_matrix,
                  JXF_ParVector    par_b,
                  JXF_ParVector    par_x );
JXF_Int
jxf_ILUZPAMGSolve( void            *data,
                  jxf_hpCSRMatrix *par_matrix,
                  jxf_ParVector    *par_b,
                  jxf_ParVector    *par_x );
JXF_Int
JXF_ILUZPAMGILUZSolve( JXF_Solver       solver,
                      JXF_hpCSRMatrix par_matrix,
                      JXF_ParVector    par_b,
                      JXF_ParVector    par_x );
JXF_Int
jxf_ILUZPAMGILUZSolve( void            *data,
                      jxf_hpCSRMatrix *par_matrix,
                      jxf_ParVector    *par_b,
                      jxf_ParVector    *par_x );
JXF_Int
JXF_PAMGILUZPAMGSolve( JXF_Solver       solver,
                      JXF_hpCSRMatrix par_matrix,
                      JXF_ParVector    par_b,
                      JXF_ParVector    par_x );
JXF_Int
jxf_PAMGILUZPAMGSolve( void            *data,
                      jxf_hpCSRMatrix *par_matrix,
                      jxf_ParVector    *par_b,
                      jxf_ParVector    *par_x );

/* csrc/combined/combined.c */
JXF_Int
JXF_CombinedPrecondDataCreate( JXF_Solver *solver, MPI_Comm comm );
void *
jxf_CombinedPrecondDataCreate( MPI_Comm comm );
JXF_Int
JXF_CombinedPrecondDataInitialize( JXF_Solver solver, JXF_Int level );
JXF_Int
jxf_CombinedPrecondDataInitialize( void *data, JXF_Int level );
JXF_Int
JXF_CombinedPrecondDataSetEuclidLevel( JXF_Solver solver, JXF_Int level );
JXF_Int
JXF_CombinedPrecondDataInitializeP( JXF_Solver solver,
                                   JXF_Int argc,
                                   char *argv[] );
JXF_Int
jxf_CombinedPrecondDataInitializeP( void *data,
                                   JXF_Int argc,
                                   char *argv[] );
JXF_Int
JXF_CombinedPrecondDataSetEuclidParams( JXF_Solver solver,
                                       JXF_Int argc,
                                       char *argv[] );
JXF_Int
JXF_CombinedPrecondDataSetPreID( JXF_Solver solver, JXF_Int pre_id );
JXF_Int
jxf_CombinedPrecondDataSetPreID( void *data, JXF_Int pre_id );
JXF_Int
JXF_CombinedPrecondDataSetNxyNpxyNequ( JXF_Solver solver, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy, JXF_Int num_equns );
JXF_Int
jxf_CombinedPrecondDataSetNxyNpxyNequ( void *data, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy, JXF_Int num_equns );
JXF_Int
JXF_CombinedPrecondDataSetThetaPsiRhoPhi( JXF_Solver solver,
                                         JXF_Int theta_psi,
                                         JXF_Int theta_rho,
                                         JXF_Int theta_phi,
                                         JXF_Real theta_dis );
JXF_Int
jxf_CombinedPrecondDataSetThetaPsiRhoPhi( void *data,
                                         JXF_Int theta_psi,
                                         JXF_Int theta_rho,
                                         JXF_Int theta_phi,
                                         JXF_Real theta_dis );
JXF_Int
JXF_CombinedPrecondDataSetAMGParameters( JXF_Solver solver,
                                        JXF_Int amg_max_levels,
                                        JXF_Int amg_relax_type,
                                        JXF_Int amg_print_level,
                                        JXF_Int amg_interp_type,
                                        JXF_Int amg_P_max_elmts,
                                        JXF_Int amg_measure_type,
                                        JXF_Int amg_coarsen_type,
                                        JXF_Int amg_agg_num_levels,
                                        JXF_Int amg_coarse_threshold,
                                        JXF_Real amg_strong_threshold );
JXF_Int
jxf_CombinedPrecondDataSetAMGParameters( void *data,
                                        JXF_Int amg_max_levels,
                                        JXF_Int amg_relax_type,
                                        JXF_Int amg_print_level,
                                        JXF_Int amg_interp_type,
                                        JXF_Int amg_P_max_elmts,
                                        JXF_Int amg_measure_type,
                                        JXF_Int amg_coarsen_type,
                                        JXF_Int amg_agg_num_levels,
                                        JXF_Int amg_coarse_threshold,
                                        JXF_Real amg_strong_threshold );
JXF_Int
JXF_CombinedPrecondDataSetPreMat( JXF_Solver solver, jxf_ParCSRMatrix *pre_mat );
JXF_Int
jxf_CombinedPrecondDataSetPreMat( void *data, jxf_ParCSRMatrix *pre_mat );
JXF_Int
JXF_CombinedPrecondDataSetDropTol( JXF_Solver solver, JXF_Real drop_tol );
JXF_Int
jxf_CombinedPrecondDataSetDropTol( void *data, JXF_Real drop_tol );
JXF_Int
JXF_CombinedPrecondDataSetInterpType( JXF_Solver solver, JXF_Int interp_type );
JXF_Int
jxf_CombinedPrecondDataSetInterpType( void *data, JXF_Int interp_type );
JXF_Int
JXF_CombinedPrecondDataSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type );
JXF_Int
jxf_CombinedPrecondDataSetCoarsenType( void *data, JXF_Int coarsen_type );
JXF_Int
JXF_CombinedPrecondDataSetCycleRelaxType( JXF_Solver solver, JXF_Int relax_type );
JXF_Int
jxf_CombinedPrecondDataSetCycleRelaxType( void *data, JXF_Int relax_type );
JXF_Int
JXF_CombinedPrecondDataSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold );
JXF_Int
jxf_CombinedPrecondDataSetStrongThreshold( void *data, JXF_Real strong_threshold );
JXF_Int
JXF_CombinedPrecondDataSetAdpTheta( JXF_Solver solver, JXF_Real adp_theta );
JXF_Int
jxf_CombinedPrecondDataSetAdpTheta( void *data, JXF_Real adp_theta );
JXF_Int
JXF_CombinedPrecondDataSetAdpGammaT( JXF_Solver solver, JXF_Real adp_gamma_3 );
JXF_Int
jxf_CombinedPrecondDataSetAdpGammaT( void *data, JXF_Real adp_gamma_3 );
JXF_Int
JXF_CombinedPrecondDataSetAdpGammaE( JXF_Solver solver, JXF_Real adp_gamma_11 );
JXF_Int
jxf_CombinedPrecondDataSetAdpGammaE( void *data, JXF_Real adp_gamma_11 );
JXF_Int
JXF_CombinedPrecondDataSetILUMatA( JXF_Solver solver, jxf_CSRMatrix *matA );
JXF_Int
jxf_CombinedPrecondDataSetILUMatA( void *data, jxf_CSRMatrix *matA );
JXF_Int
JXF_CombinedPrecondDataSetPreEucMat( JXF_Solver solver, jxf_ParCSRMatrix *matA );
JXF_Int
jxf_CombinedPrecondDataSetPreEucMat( void *data, jxf_ParCSRMatrix *matA );
JXF_Int
JXF_CombinedPrecondDataGetLULength( JXF_Solver solver, JXF_Int *lu_length );
JXF_Int
jxf_CombinedPrecondDataGetLULength( void *data, JXF_Int *lu_length );
// JXF_Int
// JXF_CombinedPrecondDataSetup( JXF_Solver solver, JXF_ParCSRMatrix par_matrix );
JXF_Int
JXF_CombinedPrecondDataSetup( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix );
JXF_Int
jxf_CombinedPrecondDataSetup( void *data, jxf_hpCSRMatrix *par_matrix );
// JXF_Int
// JXF_CombinedPrecondDataAdaptiveSetup2( JXF_Solver solver, JXF_ParCSRMatrix par_matrix );
JXF_Int
JXF_CombinedPrecondDataAdaptiveSetup2( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix );
JXF_Int
jxf_CombinedPrecondDataAdaptiveSetup2( void *data, jxf_hpCSRMatrix *par_matrix );
// JXF_Int
// JXF_CombinedPrecondDataAdaptiveSetup3( JXF_Solver solver, JXF_ParCSRMatrix par_matrix );
JXF_Int
JXF_CombinedPrecondDataAdaptiveSetup3( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix );
JXF_Int
jxf_CombinedPrecondDataAdaptiveSetup3( void *data, jxf_hpCSRMatrix *par_matrix );

JXF_Int
JXF_CombinedPrecondDataSolve( JXF_Solver       solver,
                             JXF_hpCSRMatrix par_matrix,
                             JXF_ParVector    par_b,
                             JXF_ParVector    par_x );
JXF_Int
jxf_CombinedPrecondDataSolve( void            *data,
                             jxf_hpCSRMatrix *par_matrix,
                             jxf_ParVector    *par_b,
                             jxf_ParVector    *par_x );

JXF_Int
JXF_CombinedPrecondDataAdaptiveSolve( JXF_Solver       solver,
                                     JXF_hpCSRMatrix hp_matrix,
                                     JXF_ParVector    par_b,
                                     JXF_ParVector    par_x );
JXF_Int
jxf_CombinedPrecondDataAdaptiveSolve( void            *data,
                                     jxf_hpCSRMatrix *par_matrix,
                                     jxf_ParVector    *par_b,
                                     jxf_ParVector    *par_x );
JXF_Int
JXF_CombinedPrecondDataDestroy( JXF_Solver solver );
JXF_Int
jxf_CombinedPrecondDataDestroy( void *data );

/* csrc/combined/csr_matrix.c */
JXF_Int
jxf_CSRMatrixWeaklyDiagDominant( jxf_CSRMatrix *A,
                                JXF_Int num_equns,
                                JXF_Real theta,
                                JXF_Real gamma_3,
                                JXF_Real gamma_11 );

/* csrc/combined/par_csr_matrix.c */
jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile4( char *filename, JXF_Int file_base, JXF_Int  *row_part, JXF_Int  *col_part );

#endif
