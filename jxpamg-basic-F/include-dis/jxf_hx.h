//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_hx.h -- the head file for the auxiliary space Maxwell solver
 *  Date: 2024/01/28
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JXF_HX_HEADER
#define JXF_HX_HEADER

#ifndef JXF_UTIL_HEADER
#include "jxf_util.h"
#endif

#ifndef JXF_HPCSRMV_HEADER 
#include "jxf_hpcsr.h"
#endif

#ifndef JXF_KRYLOV_HEADER
#include "jxf_krylov.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_HXData
 */
typedef struct
{
   JXF_Int dim;
   JXF_Int beta_is_zero;
   JXF_Int owns_Pi;
   JXF_Int owns_A_G, owns_A_Pi;
   JXF_Int projection_frequency;
   JXF_Int solve_counter;
   JXF_Int maxit;
   JXF_Int cycle_type;
   JXF_Int print_level;
   JXF_Int A_relax_type;
   JXF_Int A_relax_times;
   JXF_Int A_cheby_order;
   JXF_Int B_G_coarsen_type;
   JXF_Int B_G_agg_levels;
   JXF_Int B_G_relax_type;
   JXF_Int B_G_coarse_relax_type;
   JXF_Int B_G_interp_type;
   JXF_Int B_G_Pmax;
   JXF_Int B_Pi_coarsen_type;
   JXF_Int B_Pi_agg_levels;
   JXF_Int B_Pi_relax_type;
   JXF_Int B_Pi_coarse_relax_type;
   JXF_Int B_Pi_interp_type;
   JXF_Int B_Pi_Pmax;
   JXF_Int num_iterations;

   JXF_Real tol;
   JXF_Real A_relax_weight;
   JXF_Real A_omega;
   JXF_Real A_max_eig_est;
   JXF_Real A_min_eig_est;
   JXF_Real A_cheby_fraction;
   JXF_Real B_G_theta;
   JXF_Real B_Pi_theta;
   JXF_Real rel_resid_norm;

   JXF_Real *A_l1_norms;

   jxf_ParVector *x, *y, *z;
   jxf_ParVector *Gx, *Gy, *Gz;
   jxf_ParVector *interior_nodes;
   jxf_ParVector *r0, *g0, *r1, *g1, *r2, *g2;

   jxf_ParCSRMatrix *A;
   jxf_ParCSRMatrix *G;
   jxf_hpCSRMatrix *A_G;
   jxf_ParCSRMatrix *Pi;
   jxf_hpCSRMatrix *A_Pi;
   jxf_ParCSRMatrix *Pix, *Piy, *Piz;
   jxf_hpCSRMatrix *A_Pix, *A_Piy, *A_Piz;
   jxf_ParCSRMatrix *G0;
   jxf_hpCSRMatrix *A_G0;

   JXF_Solver B_G;
   JXF_Solver B_Pi;
   JXF_Solver B_Pix, B_Piy, B_Piz;
   JXF_Solver B_G0;

} jxf_HXData;

#define jxf_HXDataDimension(ams_data) ((ams_data)->dim)
#define jxf_HXDataA(ams_data) ((ams_data)->A)
#define jxf_HXDataDiscreteGradient(ams_data) ((ams_data)->G)
#define jxf_HXDataPoissonBeta(ams_data) ((ams_data)->A_G)
#define jxf_HXDataPoissonBetaAMG(ams_data) ((ams_data)->B_G)
#define jxf_HXDataOwnsPoissonBeta(ams_data) ((ams_data)->owns_A_G)
#define jxf_HXDataBetaIsZero(ams_data) ((ams_data)->beta_is_zero)
#define jxf_HXDataPiInterpolation(ams_data) ((ams_data)->Pi)
#define jxf_HXDataOwnsPiInterpolation(ams_data) ((ams_data)->owns_Pi)
#define jxf_HXDataPoissonAlpha(ams_data) ((ams_data)->A_Pi)
#define jxf_HXDataPoissonAlphaAMG(ams_data) ((ams_data)->B_Pi)
#define jxf_HXDataOwnsPoissonAlpha(ams_data) ((ams_data)->owns_A_Pi)
#define jxf_HXDataVertexCoordinateX(ams_data) ((ams_data)->x)
#define jxf_HXDataVertexCoordinateY(ams_data) ((ams_data)->y)
#define jxf_HXDataVertexCoordinateZ(ams_data) ((ams_data)->z)
#define jxf_HXDataEdgeConstantX(ams_data) ((ams_data)->Gx)
#define jxf_HXDataEdgeConstantY(ams_data) ((ams_data)->Gy)
#define jxf_HXDataEdgeConstantZ(ams_data) ((ams_data)->Gz)
#define jxf_HXDataMaxIter(ams_data) ((ams_data)->maxit)
#define jxf_HXDataTol(ams_data) ((ams_data)->tol)
#define jxf_HXDataCycleType(ams_data) ((ams_data)->cycle_type)
#define jxf_HXDataPrintLevel(ams_data) ((ams_data)->print_level)
#define jxf_HXDataARelaxType(ams_data) ((ams_data)->A_relax_type)
#define jxf_HXDataARelaxTimes(ams_data) ((ams_data)->A_relax_times)
#define jxf_HXDataAL1Norms(ams_data) ((ams_data)->A_l1_norms)
#define jxf_HXDataARelaxWeight(ams_data) ((ams_data)->A_relax_weight)
#define jxf_HXDataAOmega(ams_data) ((ams_data)->A_omega)
#define jxf_HXDataAMaxEigEst(ams_data) ((ams_data)->A_max_eig_est)
#define jxf_HXDataAMinEigEst(ams_data) ((ams_data)->A_min_eig_est)
#define jxf_HXDataAChebyOrder(ams_data) ((ams_data)->A_cheby_order)
#define jxf_HXDataAChebyFraction(ams_data) ((ams_data)->A_cheby_fraction)
#define jxf_HXDataPoissonAlphaAMGCoarsenType(ams_data) ((ams_data)->B_Pi_coarsen_type)
#define jxf_HXDataPoissonAlphaAMGAggLevels(ams_data) ((ams_data)->B_Pi_agg_levels)
#define jxf_HXDataPoissonAlphaAMGRelaxType(ams_data) ((ams_data)->B_Pi_relax_type)
#define jxf_HXDataPoissonAlphaAMGStrengthThreshold(ams_data) ((ams_data)->B_Pi_theta)
#define jxf_HXDataPoissonBetaAMGCoarsenType(ams_data) ((ams_data)->B_G_coarsen_type)
#define jxf_HXDataPoissonBetaAMGAggLevels(ams_data) ((ams_data)->B_G_agg_levels)
#define jxf_HXDataPoissonBetaAMGRelaxType(ams_data) ((ams_data)->B_G_relax_type)
#define jxf_HXDataPoissonBetaAMGStrengthThreshold(ams_data) ((ams_data)->B_G_theta)
#define jxf_HXDataTempEdgeVectorR(ams_data) ((ams_data)->r0)
#define jxf_HXDataTempEdgeVectorG(ams_data) ((ams_data)->g0)
#define jxf_HXDataTempVertexVectorR(ams_data) ((ams_data)->r1)
#define jxf_HXDataTempVertexVectorG(ams_data) ((ams_data)->g1)
#define jxf_HXDataTempVecVertexVectorR(ams_data) ((ams_data)->r2)
#define jxf_HXDataTempVecVertexVectorG(ams_data) ((ams_data)->g2)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/hx/block.c */
JXF_Int jxf_ParVectorBlockSplit( jxf_ParVector *x, jxf_ParVector *x_[3], JXF_Int dim );
JXF_Int jxf_ParVectorBlockGather( jxf_ParVector *x, jxf_ParVector *x_[3], JXF_Int dim );
JXF_Int jxf_PAMGBlockSolve( void *B, jxf_hpCSRMatrix *A, jxf_ParVector *b, jxf_ParVector *x );

/* csrc/hx/comput.c */
JXF_Int
jxf_HXComputePi( jxf_ParCSRMatrix *A,
                jxf_ParCSRMatrix *G,
                jxf_ParVector *Gx,
                jxf_ParVector *Gy,
                jxf_ParVector *Gz,
                JXF_Int dim,
                jxf_ParCSRMatrix **Pi_ptr );
JXF_Int
jxf_HXComputePixyz( jxf_ParCSRMatrix *A,
                   jxf_ParCSRMatrix *G,
                   jxf_ParVector *Gx,
                   jxf_ParVector *Gy,
                   jxf_ParVector *Gz,
                   JXF_Int dim,
                   jxf_ParCSRMatrix **Pix_ptr,
                   jxf_ParCSRMatrix **Piy_ptr,
                   jxf_ParCSRMatrix **Piz_ptr );
JXF_Int
jxf_HXComputeGPi( jxf_ParCSRMatrix *A,
                 jxf_ParCSRMatrix *G,
                 jxf_ParVector *Gx,
                 jxf_ParVector *Gy,
                 jxf_ParVector *Gz,
                 JXF_Int dim,
                 jxf_ParCSRMatrix **GPi_ptr );

/* csrc/hx/fei.c */
JXF_Int
jxf_HXFEISetup( void *solver,
               jxf_ParCSRMatrix *A,
               jxf_ParVector *b,
               jxf_ParVector *x,
               JXF_Int num_vert,
               JXF_Int num_local_vert,
               JXF_Int *vert_number,
               JXF_Real *vert_coord,
               JXF_Int num_edges,
               JXF_Int *edge_vertex );
JXF_Int jxf_HXFEIDestroy( void *solver );
JXF_Int
JXF_HXFEISetup( JXF_Solver solver,
               JXF_ParCSRMatrix A,
               JXF_ParVector b,
               JXF_ParVector x,
               JXF_Int *EdgeNodeList_,
               JXF_Int *NodeNumbers_,
               JXF_Int numEdges_,
               JXF_Int numLocalNodes_,
               JXF_Int numNodes_,
               JXF_Real *NodalCoord_ );
JXF_Int JXF_HXFEIDestroy( JXF_Solver solver );

/* csrc/hx/hx.c */
void *jxf_HXCreate();
JXF_Int jxf_HXDestroy( void *solver );
JXF_Int jxf_HXSetDimension( void *solver, JXF_Int dim );
JXF_Int jxf_HXSetDiscreteGradient( void *solver, jxf_ParCSRMatrix *G );
JXF_Int jxf_HXSetCoordinateVectors( void *solver, jxf_ParVector *x, jxf_ParVector *y, jxf_ParVector *z );
JXF_Int jxf_HXSetEdgeConstantVectors( void *solver, jxf_ParVector *Gx, jxf_ParVector *Gy, jxf_ParVector *Gz );
JXF_Int jxf_HXSetInterpolations( void *solver, jxf_ParCSRMatrix *Pi, jxf_ParCSRMatrix *Pix, jxf_ParCSRMatrix *Piy, jxf_ParCSRMatrix *Piz );
JXF_Int jxf_HXSetAlphaPoissonMatrix( void *solver, jxf_hpCSRMatrix *A_Pi );
JXF_Int jxf_HXSetBetaPoissonMatrix( void *solver, jxf_hpCSRMatrix *A_G );
JXF_Int jxf_HXSetInteriorNodes( void *solver, jxf_ParVector *interior_nodes );
JXF_Int jxf_HXSetProjectionFrequency( void *solver, JXF_Int projection_frequency );
JXF_Int jxf_HXSetMaxIter( void *solver, JXF_Int maxit );
JXF_Int jxf_HXSetTol( void *solver, JXF_Real tol );
JXF_Int jxf_HXSetCycleType( void *solver, JXF_Int cycle_type );
JXF_Int jxf_HXSetPrintLevel( void *solver, JXF_Int print_level );
JXF_Int jxf_HXSetSmoothingOptions( void *solver, JXF_Int A_relax_type, JXF_Int A_relax_times, JXF_Real A_relax_weight, JXF_Real A_omega );
JXF_Int jxf_HXSetChebySmoothingOptions( void *solver, JXF_Int A_cheby_order, JXF_Int A_cheby_fraction );
JXF_Int
jxf_HXSetAlphaAMGOptions( void *solver,
                         JXF_Int B_Pi_coarsen_type,
                         JXF_Int B_Pi_agg_levels,
                         JXF_Int B_Pi_relax_type,
                         JXF_Real B_Pi_theta,
                         JXF_Int B_Pi_interp_type,
                         JXF_Int B_Pi_Pmax );
JXF_Int jxf_HXSetAlphaAMGCoarseRelaxType( void *solver, JXF_Int B_Pi_coarse_relax_type );
JXF_Int
jxf_HXSetBetaAMGOptions( void *solver,
                        JXF_Int B_G_coarsen_type,
                        JXF_Int B_G_agg_levels,
                        JXF_Int B_G_relax_type,
                        JXF_Real B_G_theta,
                        JXF_Int B_G_interp_type,
                        JXF_Int B_G_Pmax );
JXF_Int jxf_HXSetBetaAMGCoarseRelaxType( void *solver, JXF_Int B_G_coarse_relax_type );
JXF_Int jxf_HXSetup( void *solver, jxf_ParCSRMatrix *A );
JXF_Int jxf_HXSolve( void *solver, jxf_ParCSRMatrix *A, jxf_ParVector *b, jxf_ParVector *x );
JXF_Int jxf_HXGetNumIterations( void *solver, JXF_Int *num_iterations );
JXF_Int jxf_HXGetFinalRelativeResidualNorm( void *solver, JXF_Real *rel_resid_norm );
JXF_Int JXF_HXCreate( JXF_Solver *solver );
JXF_Int JXF_HXDestroy( JXF_Solver solver );
JXF_Int JXF_HXSetup( JXF_Solver solver, JXF_hpCSRMatrix A );
JXF_Int JXF_HXSolve( JXF_Solver solver, JXF_hpCSRMatrix A, JXF_ParVector b, JXF_ParVector x );
JXF_Int JXF_HXSetDimension( JXF_Solver solver, JXF_Int dim );
JXF_Int JXF_HXSetDiscreteGradient( JXF_Solver solver, JXF_hpCSRMatrix G );
JXF_Int JXF_HXSetCoordinateVectors( JXF_Solver solver, JXF_ParVector x, JXF_ParVector y, JXF_ParVector z );
JXF_Int JXF_HXSetEdgeConstantVectors( JXF_Solver solver, JXF_ParVector Gx, JXF_ParVector Gy, JXF_ParVector Gz );
JXF_Int JXF_HXSetInterpolations( JXF_Solver solver, JXF_ParCSRMatrix Pi, JXF_ParCSRMatrix Pix, JXF_ParCSRMatrix Piy, JXF_ParCSRMatrix Piz );
JXF_Int JXF_HXSetAlphaPoissonMatrix( JXF_Solver solver, JXF_hpCSRMatrix A_alpha );
JXF_Int JXF_HXSetBetaPoissonMatrix( JXF_Solver solver, JXF_hpCSRMatrix A_beta );
JXF_Int JXF_HXSetInteriorNodes( JXF_Solver solver, JXF_ParVector interior_nodes );
JXF_Int JXF_HXSetProjectionFrequency( JXF_Solver solver, JXF_Int projection_frequency );
JXF_Int JXF_HXSetMaxIter( JXF_Solver solver, JXF_Int maxit );
JXF_Int JXF_HXSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_HXSetCycleType( JXF_Solver solver, JXF_Int cycle_type );
JXF_Int JXF_HXSetPrintLevel( JXF_Solver solver, JXF_Int print_level );
JXF_Int JXF_HXSetSmoothingOptions( JXF_Solver solver, JXF_Int relax_type, JXF_Int relax_times, JXF_Real relax_weight, JXF_Real omega );
JXF_Int JXF_HXSetChebySmoothingOptions( JXF_Solver solver, JXF_Int cheby_order, JXF_Int cheby_fraction );
JXF_Int
JXF_HXSetAlphaAMGOptions( JXF_Solver solver,
                         JXF_Int alpha_coarsen_type,
                         JXF_Int alpha_agg_levels,
                         JXF_Int alpha_relax_type,
                         JXF_Real alpha_strength_threshold,
                         JXF_Int alpha_interp_type,
                         JXF_Int alpha_Pmax );
JXF_Int JXF_HXSetAlphaAMGCoarseRelaxType( JXF_Solver solver, JXF_Int alpha_coarse_relax_type );
JXF_Int
JXF_HXSetBetaAMGOptions( JXF_Solver solver,
                        JXF_Int beta_coarsen_type,
                        JXF_Int beta_agg_levels,
                        JXF_Int beta_relax_type,
                        JXF_Real beta_strength_threshold,
                        JXF_Int beta_interp_type,
                        JXF_Int beta_Pmax );
JXF_Int JXF_HXSetBetaAMGCoarseRelaxType( JXF_Solver solver, JXF_Int beta_coarse_relax_type );
JXF_Int JXF_HXGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_HXGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *rel_resid_norm );
JXF_Int JXF_HXProjectOutGradients( JXF_Solver solver, JXF_ParVector x );
JXF_Int
JXF_HXConstructDiscreteGradient( JXF_ParCSRMatrix A,
                                JXF_ParVector x_coord,
                                JXF_Int *edge_vertex,
                                JXF_Int edge_orientation,
                                JXF_ParCSRMatrix *G );

/* csrc/hx/matrix.c */
jxf_CSRMatrix *jxf_CSRMatrixDeleteZeros( jxf_CSRMatrix *A, JXF_Real tol );
JXF_Int JXF_hpCSRMatrixRead( MPI_Comm comm, const char *file_name, JXF_hpCSRMatrix *matrix );
JXF_Int JXF_hpCSRMatrixDestroy( JXF_hpCSRMatrix matrix );

/* csrc/hx/oper.c */
jxf_ParVector *jxf_ParVectorInRangeOf( jxf_ParCSRMatrix *A );
jxf_ParVector *jxf_ParVectorInDomainOf( jxf_ParCSRMatrix *A );
JXF_Int jxf_ParCSRMatrixFixZeroRows( jxf_ParCSRMatrix *A );
JXF_Int jxf_ParCSRMatrixSetDiagRows( jxf_ParCSRMatrix *A, JXF_Real d );
JXF_Int jxf_HXProjectOutGradients( void *solver, jxf_ParVector *x );
JXF_Int
jxf_HXConstructDiscreteGradient( jxf_ParCSRMatrix *A,
                                jxf_ParVector *x_coord,
                                JXF_Int *edge_vertex,
                                JXF_Int edge_orientation,
                                jxf_ParCSRMatrix **G_ptr );

/* csrc/hx/prec.c */
JXF_Int
jxf_ParCSRSubspacePrec( jxf_ParCSRMatrix *A0,
                       JXF_Int A0_relax_type,
                       JXF_Int A0_relax_times,
                       JXF_Real *A0_l1_norms,
                       JXF_Real A0_relax_weight,
                       JXF_Real A0_omega,
                       JXF_Real A0_max_eig_est,
                       JXF_Real A0_min_eig_est,
                       JXF_Int A0_cheby_order,
                       JXF_Real A0_cheby_fraction,
                       jxf_hpCSRMatrix **A,
                       JXF_Solver *B,
                       JXF_PtrToSolverFcn *HB,
                       jxf_ParCSRMatrix **P,
                       jxf_ParVector **r,
                       jxf_ParVector **g,
                       jxf_ParVector *x,
                       jxf_ParVector *y,
                       jxf_ParVector *r0,
                       jxf_ParVector *g0,
                       char *cycle,
                       jxf_ParVector *z );

/* csrc/hx/relax.c */
JXF_Int
jxf_ParCSRRelax( jxf_ParCSRMatrix *A,
                jxf_ParVector *f,
                JXF_Int relax_type,
                JXF_Int relax_times,
                JXF_Real *l1_norms,
                JXF_Real relax_weight,
                JXF_Real omega,
                JXF_Real max_eig_est,
                JXF_Real min_eig_est,
                JXF_Int cheby_order,
                JXF_Real cheby_fraction,
                jxf_ParVector *u,
                jxf_ParVector *v,
                jxf_ParVector *z );
JXF_Int
jxf_ParCSRRelaxThreads( jxf_ParCSRMatrix *A,
                       jxf_ParVector *f,
                       JXF_Int relax_type,
                       JXF_Int relax_times,
                       JXF_Real *l1_norms,
                       JXF_Real relax_weight,
                       JXF_Real omega,
                       jxf_ParVector *u,
                       jxf_ParVector *Vtemp,
                       jxf_ParVector *z );
JXF_Int
jxf_ParCSRRelax_Cheby( jxf_ParCSRMatrix *A,
                      jxf_ParVector *f,
                      JXF_Real max_eig,
                      JXF_Real min_eig,
                      JXF_Real fraction,
                      JXF_Int order,
                      JXF_Int scale,
                      JXF_Int variant,
                      jxf_ParVector *u,
                      jxf_ParVector *v,
                      jxf_ParVector *r );
JXF_Int jxf_ParCSRMaxEigEstimateCG( jxf_ParCSRMatrix *A, JXF_Int scale, JXF_Int max_iter, JXF_Real *max_eig, JXF_Real *min_eig );
JXF_Int jxf_LINPACKcgtql1( JXF_Int *n, JXF_Real *d, JXF_Real *e, JXF_Int *ierr );
JXF_Real jxf_LINPACKcgpthy( JXF_Real *a, JXF_Real *b );

/* csrc/hx/vector.c */
JXF_Int jxf_ParVectorSetRandomValues( jxf_ParVector *v, JXF_Int seed );
JXF_Int jxf_SeqVectorSetRandomValues( jxf_Vector *x, JXF_Int seed );
JXF_Int JXF_ParVectorRead( MPI_Comm comm, const char *file_name, JXF_ParVector *vector );
JXF_Int JXF_ParVectorDestroy( JXF_ParVector vector );

#endif
