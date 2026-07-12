//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jx_hx.h -- the head file for the auxiliary space Maxwell solver
 *  Date: 2024/01/28
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JX_HX_HEADER
#define JX_HX_HEADER

#ifndef JX_UTIL_HEADER
#include "jx_util.h"
#endif

#ifndef JX_HPCSRMV_HEADER 
#include "jx_hpcsr.h"
#endif

#ifndef JX_KRYLOV_HEADER
#include "jx_krylov.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jx_HXData
 */
typedef struct
{
   JX_Int dim;
   JX_Int beta_is_zero;
   JX_Int owns_Pi;
   JX_Int owns_A_G, owns_A_Pi;
   JX_Int projection_frequency;
   JX_Int solve_counter;
   JX_Int maxit;
   JX_Int cycle_type;
   JX_Int print_level;
   JX_Int A_relax_type;
   JX_Int A_relax_times;
   JX_Int A_cheby_order;
   JX_Int B_G_coarsen_type;
   JX_Int B_G_agg_levels;
   JX_Int B_G_relax_type;
   JX_Int B_G_coarse_relax_type;
   JX_Int B_G_interp_type;
   JX_Int B_G_Pmax;
   JX_Int B_Pi_coarsen_type;
   JX_Int B_Pi_agg_levels;
   JX_Int B_Pi_relax_type;
   JX_Int B_Pi_coarse_relax_type;
   JX_Int B_Pi_interp_type;
   JX_Int B_Pi_Pmax;
   JX_Int num_iterations;

   JX_Real tol;
   JX_Real A_relax_weight;
   JX_Real A_omega;
   JX_Real A_max_eig_est;
   JX_Real A_min_eig_est;
   JX_Real A_cheby_fraction;
   JX_Real B_G_theta;
   JX_Real B_Pi_theta;
   JX_Real rel_resid_norm;

   JX_Real *A_l1_norms;

   jx_ParVector *x, *y, *z;
   jx_ParVector *Gx, *Gy, *Gz;
   jx_ParVector *interior_nodes;
   jx_ParVector *r0, *g0, *r1, *g1, *r2, *g2;

   jx_ParCSRMatrix *A;
   jx_ParCSRMatrix *G;
   jx_hpCSRMatrix *A_G;
   jx_ParCSRMatrix *Pi;
   jx_hpCSRMatrix *A_Pi;
   jx_ParCSRMatrix *Pix, *Piy, *Piz;
   jx_hpCSRMatrix *A_Pix, *A_Piy, *A_Piz;
   jx_ParCSRMatrix *G0;
   jx_hpCSRMatrix *A_G0;

   JX_Solver B_G;
   JX_Solver B_Pi;
   JX_Solver B_Pix, B_Piy, B_Piz;
   JX_Solver B_G0;

} jx_HXData;

#define jx_HXDataDimension(ams_data) ((ams_data)->dim)
#define jx_HXDataA(ams_data) ((ams_data)->A)
#define jx_HXDataDiscreteGradient(ams_data) ((ams_data)->G)
#define jx_HXDataPoissonBeta(ams_data) ((ams_data)->A_G)
#define jx_HXDataPoissonBetaAMG(ams_data) ((ams_data)->B_G)
#define jx_HXDataOwnsPoissonBeta(ams_data) ((ams_data)->owns_A_G)
#define jx_HXDataBetaIsZero(ams_data) ((ams_data)->beta_is_zero)
#define jx_HXDataPiInterpolation(ams_data) ((ams_data)->Pi)
#define jx_HXDataOwnsPiInterpolation(ams_data) ((ams_data)->owns_Pi)
#define jx_HXDataPoissonAlpha(ams_data) ((ams_data)->A_Pi)
#define jx_HXDataPoissonAlphaAMG(ams_data) ((ams_data)->B_Pi)
#define jx_HXDataOwnsPoissonAlpha(ams_data) ((ams_data)->owns_A_Pi)
#define jx_HXDataVertexCoordinateX(ams_data) ((ams_data)->x)
#define jx_HXDataVertexCoordinateY(ams_data) ((ams_data)->y)
#define jx_HXDataVertexCoordinateZ(ams_data) ((ams_data)->z)
#define jx_HXDataEdgeConstantX(ams_data) ((ams_data)->Gx)
#define jx_HXDataEdgeConstantY(ams_data) ((ams_data)->Gy)
#define jx_HXDataEdgeConstantZ(ams_data) ((ams_data)->Gz)
#define jx_HXDataMaxIter(ams_data) ((ams_data)->maxit)
#define jx_HXDataTol(ams_data) ((ams_data)->tol)
#define jx_HXDataCycleType(ams_data) ((ams_data)->cycle_type)
#define jx_HXDataPrintLevel(ams_data) ((ams_data)->print_level)
#define jx_HXDataARelaxType(ams_data) ((ams_data)->A_relax_type)
#define jx_HXDataARelaxTimes(ams_data) ((ams_data)->A_relax_times)
#define jx_HXDataAL1Norms(ams_data) ((ams_data)->A_l1_norms)
#define jx_HXDataARelaxWeight(ams_data) ((ams_data)->A_relax_weight)
#define jx_HXDataAOmega(ams_data) ((ams_data)->A_omega)
#define jx_HXDataAMaxEigEst(ams_data) ((ams_data)->A_max_eig_est)
#define jx_HXDataAMinEigEst(ams_data) ((ams_data)->A_min_eig_est)
#define jx_HXDataAChebyOrder(ams_data) ((ams_data)->A_cheby_order)
#define jx_HXDataAChebyFraction(ams_data) ((ams_data)->A_cheby_fraction)
#define jx_HXDataPoissonAlphaAMGCoarsenType(ams_data) ((ams_data)->B_Pi_coarsen_type)
#define jx_HXDataPoissonAlphaAMGAggLevels(ams_data) ((ams_data)->B_Pi_agg_levels)
#define jx_HXDataPoissonAlphaAMGRelaxType(ams_data) ((ams_data)->B_Pi_relax_type)
#define jx_HXDataPoissonAlphaAMGStrengthThreshold(ams_data) ((ams_data)->B_Pi_theta)
#define jx_HXDataPoissonBetaAMGCoarsenType(ams_data) ((ams_data)->B_G_coarsen_type)
#define jx_HXDataPoissonBetaAMGAggLevels(ams_data) ((ams_data)->B_G_agg_levels)
#define jx_HXDataPoissonBetaAMGRelaxType(ams_data) ((ams_data)->B_G_relax_type)
#define jx_HXDataPoissonBetaAMGStrengthThreshold(ams_data) ((ams_data)->B_G_theta)
#define jx_HXDataTempEdgeVectorR(ams_data) ((ams_data)->r0)
#define jx_HXDataTempEdgeVectorG(ams_data) ((ams_data)->g0)
#define jx_HXDataTempVertexVectorR(ams_data) ((ams_data)->r1)
#define jx_HXDataTempVertexVectorG(ams_data) ((ams_data)->g1)
#define jx_HXDataTempVecVertexVectorR(ams_data) ((ams_data)->r2)
#define jx_HXDataTempVecVertexVectorG(ams_data) ((ams_data)->g2)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/hx/block.c */
JX_Int jx_ParVectorBlockSplit( jx_ParVector *x, jx_ParVector *x_[3], JX_Int dim );
JX_Int jx_ParVectorBlockGather( jx_ParVector *x, jx_ParVector *x_[3], JX_Int dim );
JX_Int jx_PAMGBlockSolve( void *B, jx_hpCSRMatrix *A, jx_ParVector *b, jx_ParVector *x );

/* csrc/hx/comput.c */
JX_Int
jx_HXComputePi( jx_ParCSRMatrix *A,
                jx_ParCSRMatrix *G,
                jx_ParVector *Gx,
                jx_ParVector *Gy,
                jx_ParVector *Gz,
                JX_Int dim,
                jx_ParCSRMatrix **Pi_ptr );
JX_Int
jx_HXComputePixyz( jx_ParCSRMatrix *A,
                   jx_ParCSRMatrix *G,
                   jx_ParVector *Gx,
                   jx_ParVector *Gy,
                   jx_ParVector *Gz,
                   JX_Int dim,
                   jx_ParCSRMatrix **Pix_ptr,
                   jx_ParCSRMatrix **Piy_ptr,
                   jx_ParCSRMatrix **Piz_ptr );
JX_Int
jx_HXComputeGPi( jx_ParCSRMatrix *A,
                 jx_ParCSRMatrix *G,
                 jx_ParVector *Gx,
                 jx_ParVector *Gy,
                 jx_ParVector *Gz,
                 JX_Int dim,
                 jx_ParCSRMatrix **GPi_ptr );

/* csrc/hx/fei.c */
JX_Int
jx_HXFEISetup( void *solver,
               jx_ParCSRMatrix *A,
               jx_ParVector *b,
               jx_ParVector *x,
               JX_Int num_vert,
               JX_Int num_local_vert,
               JX_Int *vert_number,
               JX_Real *vert_coord,
               JX_Int num_edges,
               JX_Int *edge_vertex );
JX_Int jx_HXFEIDestroy( void *solver );
JX_Int
JX_HXFEISetup( JX_Solver solver,
               JX_ParCSRMatrix A,
               JX_ParVector b,
               JX_ParVector x,
               JX_Int *EdgeNodeList_,
               JX_Int *NodeNumbers_,
               JX_Int numEdges_,
               JX_Int numLocalNodes_,
               JX_Int numNodes_,
               JX_Real *NodalCoord_ );
JX_Int JX_HXFEIDestroy( JX_Solver solver );

/* csrc/hx/hx.c */
void *jx_HXCreate();
JX_Int jx_HXDestroy( void *solver );
JX_Int jx_HXSetDimension( void *solver, JX_Int dim );
JX_Int jx_HXSetDiscreteGradient( void *solver, jx_ParCSRMatrix *G );
JX_Int jx_HXSetCoordinateVectors( void *solver, jx_ParVector *x, jx_ParVector *y, jx_ParVector *z );
JX_Int jx_HXSetEdgeConstantVectors( void *solver, jx_ParVector *Gx, jx_ParVector *Gy, jx_ParVector *Gz );
JX_Int jx_HXSetInterpolations( void *solver, jx_ParCSRMatrix *Pi, jx_ParCSRMatrix *Pix, jx_ParCSRMatrix *Piy, jx_ParCSRMatrix *Piz );
JX_Int jx_HXSetAlphaPoissonMatrix( void *solver, jx_hpCSRMatrix *A_Pi );
JX_Int jx_HXSetBetaPoissonMatrix( void *solver, jx_hpCSRMatrix *A_G );
JX_Int jx_HXSetInteriorNodes( void *solver, jx_ParVector *interior_nodes );
JX_Int jx_HXSetProjectionFrequency( void *solver, JX_Int projection_frequency );
JX_Int jx_HXSetMaxIter( void *solver, JX_Int maxit );
JX_Int jx_HXSetTol( void *solver, JX_Real tol );
JX_Int jx_HXSetCycleType( void *solver, JX_Int cycle_type );
JX_Int jx_HXSetPrintLevel( void *solver, JX_Int print_level );
JX_Int jx_HXSetSmoothingOptions( void *solver, JX_Int A_relax_type, JX_Int A_relax_times, JX_Real A_relax_weight, JX_Real A_omega );
JX_Int jx_HXSetChebySmoothingOptions( void *solver, JX_Int A_cheby_order, JX_Int A_cheby_fraction );
JX_Int
jx_HXSetAlphaAMGOptions( void *solver,
                         JX_Int B_Pi_coarsen_type,
                         JX_Int B_Pi_agg_levels,
                         JX_Int B_Pi_relax_type,
                         JX_Real B_Pi_theta,
                         JX_Int B_Pi_interp_type,
                         JX_Int B_Pi_Pmax );
JX_Int jx_HXSetAlphaAMGCoarseRelaxType( void *solver, JX_Int B_Pi_coarse_relax_type );
JX_Int
jx_HXSetBetaAMGOptions( void *solver,
                        JX_Int B_G_coarsen_type,
                        JX_Int B_G_agg_levels,
                        JX_Int B_G_relax_type,
                        JX_Real B_G_theta,
                        JX_Int B_G_interp_type,
                        JX_Int B_G_Pmax );
JX_Int jx_HXSetBetaAMGCoarseRelaxType( void *solver, JX_Int B_G_coarse_relax_type );
JX_Int jx_HXSetup( void *solver, jx_ParCSRMatrix *A );
JX_Int jx_HXSolve( void *solver, jx_ParCSRMatrix *A, jx_ParVector *b, jx_ParVector *x );
JX_Int jx_HXGetNumIterations( void *solver, JX_Int *num_iterations );
JX_Int jx_HXGetFinalRelativeResidualNorm( void *solver, JX_Real *rel_resid_norm );
JX_Int JX_HXCreate( JX_Solver *solver );
JX_Int JX_HXDestroy( JX_Solver solver );
JX_Int JX_HXSetup( JX_Solver solver, JX_hpCSRMatrix A );
JX_Int JX_HXSolve( JX_Solver solver, JX_hpCSRMatrix A, JX_ParVector b, JX_ParVector x );
JX_Int JX_HXSetDimension( JX_Solver solver, JX_Int dim );
JX_Int JX_HXSetDiscreteGradient( JX_Solver solver, JX_hpCSRMatrix G );
JX_Int JX_HXSetCoordinateVectors( JX_Solver solver, JX_ParVector x, JX_ParVector y, JX_ParVector z );
JX_Int JX_HXSetEdgeConstantVectors( JX_Solver solver, JX_ParVector Gx, JX_ParVector Gy, JX_ParVector Gz );
JX_Int JX_HXSetInterpolations( JX_Solver solver, JX_ParCSRMatrix Pi, JX_ParCSRMatrix Pix, JX_ParCSRMatrix Piy, JX_ParCSRMatrix Piz );
JX_Int JX_HXSetAlphaPoissonMatrix( JX_Solver solver, JX_hpCSRMatrix A_alpha );
JX_Int JX_HXSetBetaPoissonMatrix( JX_Solver solver, JX_hpCSRMatrix A_beta );
JX_Int JX_HXSetInteriorNodes( JX_Solver solver, JX_ParVector interior_nodes );
JX_Int JX_HXSetProjectionFrequency( JX_Solver solver, JX_Int projection_frequency );
JX_Int JX_HXSetMaxIter( JX_Solver solver, JX_Int maxit );
JX_Int JX_HXSetTol( JX_Solver solver, JX_Real tol );
JX_Int JX_HXSetCycleType( JX_Solver solver, JX_Int cycle_type );
JX_Int JX_HXSetPrintLevel( JX_Solver solver, JX_Int print_level );
JX_Int JX_HXSetSmoothingOptions( JX_Solver solver, JX_Int relax_type, JX_Int relax_times, JX_Real relax_weight, JX_Real omega );
JX_Int JX_HXSetChebySmoothingOptions( JX_Solver solver, JX_Int cheby_order, JX_Int cheby_fraction );
JX_Int
JX_HXSetAlphaAMGOptions( JX_Solver solver,
                         JX_Int alpha_coarsen_type,
                         JX_Int alpha_agg_levels,
                         JX_Int alpha_relax_type,
                         JX_Real alpha_strength_threshold,
                         JX_Int alpha_interp_type,
                         JX_Int alpha_Pmax );
JX_Int JX_HXSetAlphaAMGCoarseRelaxType( JX_Solver solver, JX_Int alpha_coarse_relax_type );
JX_Int
JX_HXSetBetaAMGOptions( JX_Solver solver,
                        JX_Int beta_coarsen_type,
                        JX_Int beta_agg_levels,
                        JX_Int beta_relax_type,
                        JX_Real beta_strength_threshold,
                        JX_Int beta_interp_type,
                        JX_Int beta_Pmax );
JX_Int JX_HXSetBetaAMGCoarseRelaxType( JX_Solver solver, JX_Int beta_coarse_relax_type );
JX_Int JX_HXGetNumIterations( JX_Solver solver, JX_Int *num_iterations );
JX_Int JX_HXGetFinalRelativeResidualNorm( JX_Solver solver, JX_Real *rel_resid_norm );
JX_Int JX_HXProjectOutGradients( JX_Solver solver, JX_ParVector x );
JX_Int
JX_HXConstructDiscreteGradient( JX_ParCSRMatrix A,
                                JX_ParVector x_coord,
                                JX_Int *edge_vertex,
                                JX_Int edge_orientation,
                                JX_ParCSRMatrix *G );

/* csrc/hx/matrix.c */
jx_CSRMatrix *jx_CSRMatrixDeleteZeros( jx_CSRMatrix *A, JX_Real tol );
JX_Int JX_hpCSRMatrixRead( MPI_Comm comm, const char *file_name, JX_hpCSRMatrix *matrix );
JX_Int JX_hpCSRMatrixDestroy( JX_hpCSRMatrix matrix );

/* csrc/hx/oper.c */
jx_ParVector *jx_ParVectorInRangeOf( jx_ParCSRMatrix *A );
jx_ParVector *jx_ParVectorInDomainOf( jx_ParCSRMatrix *A );
JX_Int jx_ParCSRMatrixFixZeroRows( jx_ParCSRMatrix *A );
JX_Int jx_ParCSRMatrixSetDiagRows( jx_ParCSRMatrix *A, JX_Real d );
JX_Int jx_HXProjectOutGradients( void *solver, jx_ParVector *x );
JX_Int
jx_HXConstructDiscreteGradient( jx_ParCSRMatrix *A,
                                jx_ParVector *x_coord,
                                JX_Int *edge_vertex,
                                JX_Int edge_orientation,
                                jx_ParCSRMatrix **G_ptr );

/* csrc/hx/prec.c */
JX_Int
jx_ParCSRSubspacePrec( jx_ParCSRMatrix *A0,
                       JX_Int A0_relax_type,
                       JX_Int A0_relax_times,
                       JX_Real *A0_l1_norms,
                       JX_Real A0_relax_weight,
                       JX_Real A0_omega,
                       JX_Real A0_max_eig_est,
                       JX_Real A0_min_eig_est,
                       JX_Int A0_cheby_order,
                       JX_Real A0_cheby_fraction,
                       jx_hpCSRMatrix **A,
                       JX_Solver *B,
                       JX_PtrToSolverFcn *HB,
                       jx_ParCSRMatrix **P,
                       jx_ParVector **r,
                       jx_ParVector **g,
                       jx_ParVector *x,
                       jx_ParVector *y,
                       jx_ParVector *r0,
                       jx_ParVector *g0,
                       char *cycle,
                       jx_ParVector *z );

/* csrc/hx/relax.c */
JX_Int
jx_ParCSRRelax( jx_ParCSRMatrix *A,
                jx_ParVector *f,
                JX_Int relax_type,
                JX_Int relax_times,
                JX_Real *l1_norms,
                JX_Real relax_weight,
                JX_Real omega,
                JX_Real max_eig_est,
                JX_Real min_eig_est,
                JX_Int cheby_order,
                JX_Real cheby_fraction,
                jx_ParVector *u,
                jx_ParVector *v,
                jx_ParVector *z );
JX_Int
jx_ParCSRRelaxThreads( jx_ParCSRMatrix *A,
                       jx_ParVector *f,
                       JX_Int relax_type,
                       JX_Int relax_times,
                       JX_Real *l1_norms,
                       JX_Real relax_weight,
                       JX_Real omega,
                       jx_ParVector *u,
                       jx_ParVector *Vtemp,
                       jx_ParVector *z );
JX_Int
jx_ParCSRRelax_Cheby( jx_ParCSRMatrix *A,
                      jx_ParVector *f,
                      JX_Real max_eig,
                      JX_Real min_eig,
                      JX_Real fraction,
                      JX_Int order,
                      JX_Int scale,
                      JX_Int variant,
                      jx_ParVector *u,
                      jx_ParVector *v,
                      jx_ParVector *r );
JX_Int jx_ParCSRMaxEigEstimateCG( jx_ParCSRMatrix *A, JX_Int scale, JX_Int max_iter, JX_Real *max_eig, JX_Real *min_eig );
JX_Int jx_LINPACKcgtql1( JX_Int *n, JX_Real *d, JX_Real *e, JX_Int *ierr );
JX_Real jx_LINPACKcgpthy( JX_Real *a, JX_Real *b );

/* csrc/hx/vector.c */
JX_Int jx_ParVectorSetRandomValues( jx_ParVector *v, JX_Int seed );
JX_Int jx_SeqVectorSetRandomValues( jx_Vector *x, JX_Int seed );
JX_Int JX_ParVectorRead( MPI_Comm comm, const char *file_name, JX_ParVector *vector );
JX_Int JX_ParVectorDestroy( JX_ParVector vector );

#endif
