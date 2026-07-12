//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_krylov.h -- head files for Krylov subspace iterative methods
 *  Date: 2009/10/08
 * 
 *  Created by peghoty
 */ 

#ifndef JXF_KRYLOV_HEADER
#define JXF_KRYLOV_HEADER 

#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

/*----------------------------------------------------------------*
 *                      Macro Definition                          *
 *----------------------------------------------------------------*/ 
                                   
#define jxf_CTAllocF(type, count, funcs) \
( (type *)(*(funcs->CAlloc))((JXF_UInt)(count), (JXF_UInt)sizeof(type)) )
#define jxf_TFreeF( ptr, funcs ) \
( (*(funcs->Free))((char *)ptr), ptr = NULL )

typedef JXF_Int (*JXF_PtrToSolverFcn)(JXF_Solver, JXF_Matrix, JXF_Vector, JXF_Vector);
typedef JXF_Int (*JXF_PtrToModifyPCFcn)(JXF_Solver, JXF_Int, JXF_Real);


/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_PCGFunctions
 * \brief fuctions used in PCG.
 */ 
typedef struct
{
   char * (*CAlloc)        ( size_t count, size_t elt_size );
   JXF_Int    (*Free)          ( char *ptr );
   JXF_Int    (*CommInfo)      ( void *A, JXF_Int *my_id, JXF_Int *num_procs );
   void * (*CreateVector)  ( void *vector );
   JXF_Int    (*DestroyVector) ( void *vector );
   void * (*MatvecCreate)  ( void *A, void *x );
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y );
   JXF_Int    (*MatvecDestroy) ( void *matvec_data );
   JXF_Real (*InnerProd)     ( void *x, void *y );
   JXF_Int    (*CopyVector)    ( void *x, void *y );
   JXF_Int    (*ClearVector)   ( void *x );
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x );
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y );

   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();

} jxf_PCGFunctions;

/*!
 * \struct jxf_PCGData
 */  
typedef struct
{
   JXF_Real   tol;
   JXF_Real   atolf;
   JXF_Real   cf_tol;
   JXF_Real   a_tol;
   JXF_Int      max_iter;
   JXF_Int      two_norm;
   JXF_Int      rel_change;
   JXF_Int      recompute_residual;
   JXF_Int      stop_crit;
   JXF_Int      converged;
   JXF_Int      conv_criteria;      /* peghoty 2010/06/23 */
   JXF_Real   convfac_threshold;  /* peghoty 2010/06/23 */

   void    *A;
   void    *B; /* precond matrix  peghoty 2011/09/04 */  
   void    *p;
   void    *s;
   void    *r; /* ...contains the residual.  This is currently kept permanently.
                  If that is ever changed, it still must be kept if logging > 1 */

   JXF_Int      owns_matvec_data;  /* normally 1; if 0, don't delete it */
   void    *matvec_data;
   void    *precond_data;

   jxf_PCGFunctions *functions;

   /* log info (always logged) */
   JXF_Int      num_iterations;
   JXF_Real   rel_residual_norm;

   JXF_Int      print_level; /* printing when print_level > 0 */
   JXF_Int      logging;     /* extra computations for logging when logging > 0 */
   JXF_Real  *norms;
   JXF_Real  *rel_norms;

} jxf_PCGData;

/*!
 * \struct jxf_GMRESFunctions
 */  
typedef struct
{
   char * (*CAlloc)            ( size_t count, size_t elt_size );
   JXF_Int    (*Free)              ( char *ptr );
   JXF_Int    (*CommInfo)          ( void  *A, JXF_Int   *my_id, JXF_Int   *num_procs );
   void * (*CreateVector)      ( void *vector );
   void * (*CreateVectorArray) ( JXF_Int size, void *vectors );
   JXF_Int    (*DestroyVector)     ( void *vector );
   void * (*MatvecCreate)      ( void *A, void *x );
   JXF_Int    (*Matvec)            ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y );
   JXF_Int    (*MatvecDestroy)     ( void *matvec_data );
   JXF_Real (*InnerProd)         ( void *x, void *y );
   JXF_Int    (*CopyVector)        ( void *x, void *y );
   JXF_Int    (*ClearVector)       ( void *x );
   JXF_Int    (*ScaleVector)       ( JXF_Real alpha, void *x );
   JXF_Int    (*Axpy)              ( JXF_Real alpha, void *x, void *y );

   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();

} jxf_GMRESFunctions;

/*!
 * \struct jxf_GMRESData
 */
typedef struct
{
   JXF_Int      k_dim;
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      rel_change;
   JXF_Int      stop_crit;
   JXF_Int      conv_criteria;       /* peghoty,  2010/06/23 */
   JXF_Real   convfac_threshold;   /* peghoty,  2010/06/23 */
   JXF_Real   resdown_0_threshold; /* Yue Xiaoqiang, 2014/03/24 */
   JXF_Real   convfac_threshold_2; /* Yue Xiaoqiang, 2014/03/24 */
   JXF_Int      is_check_restarted;  /* peghoty,  2011/11/08 */
   JXF_Int      break_adaptive;      /* Yue Xiaoqiang, 2014/10/24 */
   JXF_Int      last_precond_type;   /* Yue Xiaoqiang, 2014/10/24 */
   JXF_Int      comp_amg_iter;       /* Yue Xiaoqiang, 2016/02/20 */
   JXF_Int      comp_ilu_iter;       /* Yue Xiaoqiang, 2016/02/20 */
   JXF_Int      comp_bco_iter;       /* Yue Xiaoqiang, 2016/02/20 */
   JXF_Int      converged;
   JXF_Real   tol;
   JXF_Real   cf_tol;
   JXF_Real   a_tol;
   JXF_Real   rel_residual_norm;

   void    *A;
   void    *B; /* precond matrix  peghoty 2011/09/04 */   
   void    *r;
   void    *w;
   void    *w_2;
   void   **p;

   void    *matvec_data;
   void    *precond_data;

   jxf_GMRESFunctions *functions;

   /* log info (always logged) */
   JXF_Int     num_iterations;
 
   JXF_Int     print_level; /* printing when print_level > 0 */
   JXF_Int     logging;     /* extra computations for logging when logging > 0 */
   JXF_Real *norms;
   char   *log_file_name;

} jxf_GMRESData;

/*!
 * \struct jxf_COGMRESFunctions
 */ 
typedef struct
{
   char *       (*CAlloc)        ( size_t count, size_t elt_size );
   JXF_Int    (*Free)          ( char *ptr );
   JXF_Int    (*CommInfo)      ( void  *A, JXF_Int   *my_id, JXF_Int   *num_procs );
   void *       (*CreateVector)  ( void *vector );
   void *       (*CreateVectorArray)  ( JXF_Int size, void *vectors );
   JXF_Int    (*DestroyVector) ( void *vector );
   void *       (*MatvecCreate)  ( void *A, void *x );
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y );
   JXF_Int    (*MatvecDestroy) ( void *matvec_data );
   JXF_Real   (*InnerProd)     ( void *x, void *y );
   JXF_Int    (*MassInnerProd) ( void *x, void **p, JXF_Int k, JXF_Int unroll, void *result );
   JXF_Int    (*MassDotpTwo)   ( void *x, void *y, void **p, JXF_Int k, JXF_Int unroll, void *result_x, void *result_y );
   JXF_Int    (*CopyVector)    ( void *x, void *y );
   JXF_Int    (*ClearVector)   ( void *x );
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x );
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y );
   JXF_Int    (*MassAxpy)      ( JXF_Real *alpha, void **x, void *y, JXF_Int k, JXF_Int unroll );
   JXF_Int    (*precond)       ();
   JXF_Int    (*precond_setup) ();

   JXF_Int    (*modify_pc)( void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm );

} jxf_COGMRESFunctions;

/*!
 * \struct jxf_COGMRESData
 */
typedef struct
{
   JXF_Int      k_dim;
   JXF_Int      unroll;
   JXF_Int      cgs;
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      rel_change;
   JXF_Int      is_check_restarted;
   JXF_Int      skip_real_r_check;
   JXF_Int      converged;
   JXF_Real   tol;
   JXF_Real   cf_tol;
   JXF_Real   a_tol;
   JXF_Real   rel_residual_norm;

   void  *A;
   void  *r;
   void  *w;
   void  *w_2;
   void  **p;

   void    *matvec_data;
   void    *precond_data;

   jxf_COGMRESFunctions * functions;

   /* log info (always logged) */
   JXF_Int      num_iterations;
 
   JXF_Int     print_level; /* printing when print_level>0 */
   JXF_Int     logging;  /* extra computations for logging when logging>0 */
   JXF_Real  *norms;
   char    *log_file_name;

} jxf_COGMRESData;

/*!
 * \struct jxf_FlexGMRESFunctions
 */
typedef struct
{
   char *       (*CAlloc)        ( size_t count, size_t elt_size );
   JXF_Int    (*Free)          ( char *ptr );
   JXF_Int    (*CommInfo)      ( void  *A, JXF_Int   *my_id, JXF_Int   *num_procs );
   void *       (*CreateVector)  ( void *vector );
   void *       (*CreateVectorArray)  ( JXF_Int size, void *vectors );
   JXF_Int    (*DestroyVector) ( void *vector );
   void *       (*MatvecCreate)  ( void *A, void *x );
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y );
   JXF_Int    (*MatvecDestroy) ( void *matvec_data );
   JXF_Real   (*InnerProd)     ( void *x, void *y );
   JXF_Int    (*CopyVector)    ( void *x, void *y );
   JXF_Int    (*ClearVector)   ( void *x );
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x );
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y );

   JXF_Int    (*precond)(void *vdata, void *A, void *b, void *x );
   JXF_Int    (*precond_setup)(void *vdata, void *A );

   JXF_Int    (*modify_pc)(void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm );

} jxf_FlexGMRESFunctions;

/*!
 * \struct jxf_FlexGMRESData
 */
typedef struct
{
   JXF_Int      k_dim;
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      rel_change;
   JXF_Int      stop_crit;
   JXF_Int      converged;
   JXF_Int      is_check_restarted;
   JXF_Real   tol;
   JXF_Real   cf_tol;
   JXF_Real   a_tol;
   JXF_Real   rel_residual_norm;

   void   **pre_vecs;

   void  *A;
   void  *r;
   void  *w;
   void  *w_2;
   void  **p;

   void    *matvec_data;
   void    *precond_data;

   jxf_FlexGMRESFunctions * functions;

   /* log info (always logged) */
   JXF_Int      num_iterations;
 
   JXF_Int     print_level; /* printing when print_level>0 */
   JXF_Int     logging;  /* extra computations for logging when logging>0 */
   JXF_Real  *norms;
   char    *log_file_name;

} jxf_FlexGMRESData;

/*!
 * \struct jxf_BiCGSTABFunctions
 */ 
typedef struct
{
  void * (*CreateVector)  ( void *vvector );
  JXF_Int    (*DestroyVector) ( void *vvector );
  void * (*MatvecCreate)  ( void *A, void *x );
  JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y );
  JXF_Int    (*MatvecDestroy) ( void *matvec_data );
  JXF_Real (*InnerProd)     ( void *x, void *y );
  JXF_Int    (*CopyVector)    ( void *x, void *y );
  JXF_Int    (*ClearVector)   ( void *x );
  JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x );
  JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x , void *y );
  JXF_Int    (*CommInfo)      ( void *A, JXF_Int *my_id, JXF_Int *num_procs );
  
  JXF_Int    (*precond_setup)();
  JXF_Int    (*precond)();

} jxf_BiCGSTABFunctions;


/*!
 * \struct jxf_BiCGSTABData
 */ 
typedef struct
{
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      conv_criteria;      /* peghoty 2010/03/20 */
   JXF_Real   convfac_threshold;  /* peghoty 2010/06/23 */
   JXF_Int      converged;
   JXF_Real   tol;
   JXF_Real   cf_tol;
   JXF_Real   rel_residual_norm;
   JXF_Real   a_tol;
   
   void    *A;
   void    *B; /* precond matrix  peghoty 2011/09/04 */
   void    *r;
   void    *r0;
   void    *s;
   void    *v;
   void    *p;
   void    *q;
   void    *Xtmp; /* peghoty 2010/03/20 */

   void    *matvec_data;
   void    *precond_data;

   jxf_BiCGSTABFunctions *functions;

   /* log info (always logged) */
   JXF_Int      num_iterations;
 
   /* additional log info (logged when 'logging' > 0) */
   JXF_Int      logging;
   JXF_Int      print_level;
   JXF_Real  *norms;
   char    *log_file_name;

} jxf_BiCGSTABData;


/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/
 
/* /csrc/krylov/krylov.c */
char *jxf_ParKrylovCAlloc( JXF_Int count, JXF_Int elt_size );                              
JXF_Int jxf_ParKrylovFree( char *ptr );                                                  
JXF_Int jxf_ParKrylovCommInfo( void *A, JXF_Int *my_id, JXF_Int *num_procs );                                                 
void *jxf_ParKrylovCreateVector( void *vvector );                           
void *jxf_ParKrylovCreateVectorArray( JXF_Int n, void *vvector );                        
JXF_Int jxf_ParKrylovDestroyVector( void *vvector );                                                                        
void *jxf_ParKrylovMatvecCreate( void *A, void *x );                      
JXF_Int jxf_ParKrylovMatvec( void *matvec_data, JXF_Real alpha, void *A, void *x,JXF_Real beta, void *y );
JXF_Int jxf_ParKrylovMatvecDestroy( void *matvec_data );
JXF_Real jxf_ParKrylovInnerProd( void *x, void *y );
JXF_Int jxf_ParKrylovMassInnerProd( void *x, void **y, JXF_Int k, JXF_Int unroll, void  *result );
JXF_Int jxf_ParKrylovMassDotpTwo( void *x, void *y, void **z, JXF_Int k, JXF_Int unroll, void *result_x, void *result_y );
JXF_Int jxf_ParKrylovCopyVector( void *x, void *y );
JXF_Int jxf_ParKrylovClearVector( void *x );
JXF_Int jxf_ParKrylovScaleVector( JXF_Real alpha, void *x );
JXF_Int jxf_ParKrylovAxpy( JXF_Real alpha, void *x, void *y );
JXF_Int jxf_ParKrylovMassAxpy( JXF_Real *alpha, void **x, void *y, JXF_Int k, JXF_Int unroll );
JXF_Int jxf_ParKrylovIdentitySetup( void *vdata, void *A );            
JXF_Int jxf_ParKrylovIdentity( void *vdata, void *A, void *b, void *x );

/* /csrc/krylov/pcg.c */
//JXF_Int JXF_ParCSRPCGCreate( MPI_Comm comm, JXF_Solver *solver );                
JXF_Int JXF_PCGCreate( MPI_Comm comm, JXF_Solver *solver );
jxf_PCGFunctions *
jxf_PCGFunctionsCreate(
   char * (*CAlloc)        ( size_t count, size_t elt_size ),
   JXF_Int    (*Free)          ( char *ptr ),
   JXF_Int    (*CommInfo)      ( void *A, JXF_Int *my_id, JXF_Int *num_procs ),
   void * (*CreateVector)  ( void *vector ),
   JXF_Int    (*DestroyVector) ( void *vector ),
   void * (*MatvecCreate)  ( void *A, void *x ),
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y ),
   JXF_Int    (*MatvecDestroy) ( void *matvec_data ),
   JXF_Real (*InnerProd)     ( void *x, void *y ),
   JXF_Int    (*CopyVector)    ( void *x, void *y ),
   JXF_Int    (*ClearVector)   ( void *x ),
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x ),
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y ),
   JXF_Int    (*PrecondSetup)  ( void *vdata, void *A ),
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
);                                              
void *jxf_PCGCreate( jxf_PCGFunctions *pcg_functions );                                     
JXF_Int JXF_PCGSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );                     
JXF_Int JXF_PCGSolve( JXF_Solver solver, JXF_Matrix  preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );  

JXF_Int JXF_PCGSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator ); /* peghoty 2011/09/04 */
JXF_Int JXF_PCGSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_PCGGetTol( JXF_Solver solver, JXF_Real *tol );                       
JXF_Int JXF_PCGSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol );
JXF_Int JXF_PCGGetAbsoluteTol( JXF_Solver solver, JXF_Real *a_tol );                        
JXF_Int JXF_PCGSetAbsoluteTolFactor( JXF_Solver solver, JXF_Real abstolf );
JXF_Int JXF_PCGGetAbsoluteTolFactor( JXF_Solver solver, JXF_Real *abstolf );                                        
JXF_Int JXF_PCGSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol );
JXF_Int JXF_PCGGetConvergenceFactorTol( JXF_Solver solver, JXF_Real *cf_tol );                                         
JXF_Int JXF_PCGSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_PCGGetMaxIter( JXF_Solver solver, JXF_Int *max_iter );                                                            
JXF_Int JXF_PCGSetStopCrit( JXF_Solver solver, JXF_Int stop_crit );
JXF_Int JXF_PCGGetStopCrit( JXF_Solver solver, JXF_Int *stop_crit );               
JXF_Int JXF_PCGSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria ); /* peghoty 2010/06/23 */
JXF_Int JXF_PCGSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold );                 
JXF_Int JXF_PCGSetTwoNorm( JXF_Solver solver, JXF_Int two_norm );
JXF_Int JXF_PCGGetTwoNorm( JXF_Solver solver, JXF_Int *two_norm );               
JXF_Int JXF_PCGSetRelChange( JXF_Solver solver, JXF_Int rel_change );
JXF_Int JXF_PCGGetRelChange( JXF_Solver solver, JXF_Int *rel_change );
JXF_Int JXF_PCGSetRecomputeResidual( JXF_Solver solver, JXF_Int recompute_residual );
JXF_Int JXF_PCGGetRecomputeResidual( JXF_Solver solver, JXF_Int *recompute_residual );
JXF_Int 
JXF_PCGSetPrecond( JXF_Solver          solver, 
                  JXF_PtrToSolverFcn  precond,
                  JXF_PtrToSolverFcn  precond_setup, 
                  JXF_Solver          precond_solver );
JXF_Int JXF_PCGGetPrecond( JXF_Solver solver, JXF_Solver *precond_data_ptr );
JXF_Int JXF_PCGSetLogging( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_PCGGetLogging( JXF_Solver solver, JXF_Int *level );
JXF_Int JXF_PCGSetPrintLevel( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_PCGGetPrintLevel( JXF_Solver solver, JXF_Int *level );
JXF_Int JXF_PCGGetNumIterations( JXF_Solver  solver, JXF_Int *num_iterations );                         
JXF_Int JXF_PCGGetConverged( JXF_Solver solver, JXF_Int *converged );                                                    
JXF_Int JXF_PCGGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm );                           
JXF_Int JXF_PCGGetResidual( JXF_Solver solver, void **residual );                                                   
//JXF_Int JXF_ParCSRPCGDestroy( JXF_Solver solver );   
JXF_Int JXF_PCGDestroy( JXF_Solver solver );                                             
JXF_Int jxf_PCGDestroy( void *pcg_vdata );
JXF_Int jxf_PCGGetResidual( void *pcg_vdata, void **residual );
JXF_Int jxf_PCGSetup( void *pcg_vdata, void *matA, void *vecB, void *vecX );
JXF_Int jxf_PCGSolve( void *pcg_vdata, void *preOperator, void *matA, void *vecB, void *vecX );
JXF_Int jxf_PCGSetPreOperator( void *pcg_vdata, void *preOperator ); /* peghoty 2011/09/04 */
JXF_Int jxf_PCGSetTol( void *pcg_vdata, JXF_Real tol );
JXF_Int jxf_PCGGetTol( void *pcg_vdata, JXF_Real *tol );
JXF_Int jxf_PCGSetAbsoluteTol( void *pcg_vdata, JXF_Real a_tol );
JXF_Int jxf_PCGGetAbsoluteTol( void *pcg_vdata, JXF_Real *a_tol );
JXF_Int jxf_PCGSetAbsoluteTolFactor( void *pcg_vdata, JXF_Real atolf );
JXF_Int jxf_PCGGetAbsoluteTolFactor( void *pcg_vdata, JXF_Real *atolf );
JXF_Int jxf_PCGSetConvergenceFactorTol( void *pcg_vdata, JXF_Real cf_tol );
JXF_Int jxf_PCGGetConvergenceFactorTol( void *pcg_vdata, JXF_Real *cf_tol );
JXF_Int jxf_PCGSetMaxIter( void *pcg_vdata, JXF_Int max_iter );
JXF_Int jxf_PCGGetMaxIter( void *pcg_vdata, JXF_Int *max_iter );
JXF_Int jxf_PCGSetTwoNorm( void *pcg_vdata, JXF_Int two_norm );
JXF_Int jxf_PCGGetTwoNorm( void *pcg_vdata, JXF_Int *two_norm );
JXF_Int jxf_PCGSetRelChange( void *pcg_vdata, JXF_Int rel_change );
JXF_Int jxf_PCGGetRelChange( void *pcg_vdata, JXF_Int *rel_change );
JXF_Int jxf_PCGSetRecomputeResidual( void *pcg_vdata, JXF_Int recompute_residual );
JXF_Int jxf_PCGGetRecomputeResidual( void *pcg_vdata, JXF_Int *recompute_residual );
JXF_Int jxf_PCGSetStopCrit( void *pcg_vdata, JXF_Int stop_crit );
JXF_Int jxf_PCGGetStopCrit( void *pcg_vdata, JXF_Int * stop_crit );
JXF_Int jxf_PCGSetConvCriteria( void *pcg_vdata, JXF_Int conv_criteria ); /* peghoty 2010/06/23 */
JXF_Int jxf_PCGSetConvFacThreshold( void *pcg_vdata, JXF_Real convfac_threshold ); /* peghoty 2010/06/23 */
JXF_Int jxf_PCGGetPrecond( void *pcg_vdata, JXF_Solver *precond_data_ptr );
JXF_Int jxf_PCGSetPrecond( void *pcg_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data ); 
JXF_Int jxf_PCGSetPrintLevel( void *pcg_vdata, JXF_Int level );
JXF_Int jxf_PCGGetPrintLevel( void *pcg_vdata, JXF_Int *level );
JXF_Int jxf_PCGSetLogging( void *pcg_vdata, JXF_Int level );
JXF_Int jxf_PCGGetLogging( void *pcg_vdata, JXF_Int *level);
JXF_Int jxf_PCGGetNumIterations( void *pcg_vdata, JXF_Int *num_iterations );
JXF_Int jxf_PCGGetConverged( void *pcg_vdata, JXF_Int *converged );
JXF_Int jxf_PCGPrintLogging( void *pcg_vdata, JXF_Int myid );
JXF_Int jxf_PCGGetFinalRelativeResidualNorm( void *pcg_vdata, JXF_Real *relative_residual_norm );

/* /csrc/krylov/pcogmres.c */
JXF_Int JXF_ParCSRCOGMRESCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int JXF_ParCSRCOGMRESDestroy( JXF_Solver solver );
JXF_Int JXF_COGMRESSetup( JXF_Solver solver, JXF_Matrix A, JXF_Vector b, JXF_Vector x );
JXF_Int JXF_COGMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix A, JXF_Vector b, JXF_Vector x );
JXF_Int JXF_COGMRESSetKDim( JXF_Solver solver, JXF_Int k_dim );
JXF_Int JXF_COGMRESGetKDim( JXF_Solver solver, JXF_Int * k_dim );
JXF_Int JXF_COGMRESSetUnroll( JXF_Solver solver, JXF_Int unroll );
JXF_Int JXF_COGMRESGetUnroll( JXF_Solver solver, JXF_Int * unroll );
JXF_Int JXF_COGMRESSetCGS( JXF_Solver solver, JXF_Int cgs );
JXF_Int JXF_COGMRESGetCGS( JXF_Solver solver, JXF_Int * cgs );
JXF_Int JXF_COGMRESSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_COGMRESGetTol( JXF_Solver solver, JXF_Real * tol );
JXF_Int JXF_COGMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol );
JXF_Int JXF_COGMRESGetAbsoluteTol( JXF_Solver solver, JXF_Real * a_tol );
JXF_Int JXF_COGMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol );
JXF_Int JXF_COGMRESGetConvergenceFactorTol( JXF_Solver solver, JXF_Real * cf_tol );
JXF_Int JXF_COGMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter );
JXF_Int JXF_COGMRESGetMinIter( JXF_Solver solver, JXF_Int * min_iter );
JXF_Int JXF_COGMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_COGMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted );
JXF_Int JXF_COGMRESGetMaxIter( JXF_Solver solver, JXF_Int * max_iter );
JXF_Int
JXF_COGMRESSetPrecond( JXF_Solver          solver,
                      JXF_PtrToSolverFcn  precond,
                      JXF_PtrToSolverFcn  precond_setup,
                      JXF_Solver          precond_solver );
JXF_Int JXF_COGMRESGetPrecond( JXF_Solver  solver, JXF_Solver *precond_data_ptr );
JXF_Int JXF_COGMRESSetPrintLevel( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_COGMRESGetPrintLevel( JXF_Solver solver, JXF_Int * level );
JXF_Int JXF_COGMRESSetLogging( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_COGMRESGetLogging( JXF_Solver solver, JXF_Int * level );
JXF_Int JXF_COGMRESGetNumIterations( JXF_Solver  solver, JXF_Int *num_iterations );
JXF_Int JXF_COGMRESGetConverged( JXF_Solver  solver, JXF_Int *converged );
JXF_Int JXF_COGMRESGetFinalRelativeResidualNorm( JXF_Solver  solver, JXF_Real *norm );
JXF_Int JXF_COGMRESGetResidual( JXF_Solver solver, void *residual );
JXF_Int JXF_COGMRESSetModifyPC( JXF_Solver  solver, JXF_Int (*modify_pc)(JXF_Solver, JXF_Int, JXF_Real) );
jxf_COGMRESFunctions *
jxf_COGMRESFunctionsCreate(
   char *       (*CAlloc)        ( size_t count, size_t elt_size ),
   JXF_Int    (*Free)          ( char *ptr ),
   JXF_Int    (*CommInfo)      ( void  *A, JXF_Int   *my_id, JXF_Int   *num_procs ),
   void *       (*CreateVector)  ( void *vector ),
   void *       (*CreateVectorArray)  ( JXF_Int size, void *vectors ),
   JXF_Int    (*DestroyVector) ( void *vector ),
   void *       (*MatvecCreate)  ( void *A, void *x ),
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y ),
   JXF_Int    (*MatvecDestroy) ( void *matvec_data ),
   JXF_Real   (*InnerProd)     ( void *x, void *y ),
   JXF_Int    (*MassInnerProd) ( void *x, void **y, JXF_Int k, JXF_Int unroll, void *result ),
   JXF_Int    (*MassDotpTwo)   ( void *x, void *y, void **z, JXF_Int k, JXF_Int unroll, void *result_x, void *result_y ),
   JXF_Int    (*CopyVector)    ( void *x, void *y ),
   JXF_Int    (*ClearVector)   ( void *x ),
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x ),
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y ),
   JXF_Int    (*MassAxpy)      ( JXF_Real *alpha, void **x, void *y, JXF_Int k, JXF_Int unroll ),
   JXF_Int    (*PrecondSetup)  ( void *vdata, void *A ),
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
);
void *jxf_COGMRESCreate( jxf_COGMRESFunctions *cogmres_functions );
JXF_Int jxf_COGMRESDestroy( void *cogmres_vdata );
JXF_Int jxf_COGMRESGetResidual( void *cogmres_vdata, void **residual );
JXF_Int jxf_COGMRESSetup( void *cogmres_vdata, void *A, void *b, void *x );
JXF_Int jxf_COGMRESSolve( void  *cogmres_vdata, void *preOperator, void  *A, void  *b, void  *x );
JXF_Int jxf_COGMRESSetKDim( void   *cogmres_vdata, JXF_Int   k_dim );
JXF_Int jxf_COGMRESGetKDim( void   *cogmres_vdata, JXF_Int * k_dim );
JXF_Int jxf_COGMRESSetUnroll( void   *cogmres_vdata, JXF_Int   unroll );
JXF_Int jxf_COGMRESGetUnroll( void   *cogmres_vdata, JXF_Int * unroll );
JXF_Int jxf_COGMRESSetCGS( void   *cogmres_vdata, JXF_Int   cgs );
JXF_Int jxf_COGMRESGetCGS( void   *cogmres_vdata, JXF_Int * cgs );
JXF_Int jxf_COGMRESSetTol( void   *cogmres_vdata, JXF_Real  tol );
JXF_Int jxf_COGMRESGetTol( void   *cogmres_vdata, JXF_Real  * tol );
JXF_Int jxf_COGMRESSetAbsoluteTol( void   *cogmres_vdata, JXF_Real  a_tol );
JXF_Int jxf_COGMRESGetAbsoluteTol( void   *cogmres_vdata, JXF_Real  * a_tol );
JXF_Int jxf_COGMRESSetConvergenceFactorTol( void   *cogmres_vdata, JXF_Real  cf_tol );
JXF_Int jxf_COGMRESGetConvergenceFactorTol( void   *cogmres_vdata, JXF_Real * cf_tol );
JXF_Int jxf_COGMRESSetMinIter( void *cogmres_vdata, JXF_Int   min_iter );
JXF_Int jxf_COGMRESGetMinIter( void *cogmres_vdata, JXF_Int * min_iter );
JXF_Int jxf_COGMRESSetMaxIter( void *cogmres_vdata, JXF_Int   max_iter );
JXF_Int jxf_COGMRESSetIsCheckRestarted( void *cogmres_vdata, JXF_Int is_check_restarted );
JXF_Int jxf_COGMRESGetMaxIter( void *cogmres_vdata, JXF_Int * max_iter );
JXF_Int jxf_COGMRESSetRelChange( void *cogmres_vdata, JXF_Int   rel_change );
JXF_Int jxf_COGMRESGetRelChange( void *cogmres_vdata, JXF_Int * rel_change );
JXF_Int jxf_COGMRESSetSkipRealResidualCheck( void *cogmres_vdata, JXF_Int skip_real_r_check );
JXF_Int jxf_COGMRESGetSkipRealResidualCheck( void *cogmres_vdata, JXF_Int *skip_real_r_check );
JXF_Int
jxf_COGMRESSetPrecond( void  *cogmres_vdata,
        JXF_Int  (*precond)(void*,void*,void*,void*),
        JXF_Int  (*precond_setup)(void*,void*,void*,void*),
        void  *precond_data );
JXF_Int jxf_COGMRESGetPrecond( void *cogmres_vdata, JXF_Solver *precond_data_ptr );
JXF_Int jxf_COGMRESSetPrintLevel( void *cogmres_vdata, JXF_Int   level );
JXF_Int jxf_COGMRESGetPrintLevel( void *cogmres_vdata, JXF_Int * level );
JXF_Int jxf_COGMRESSetLogging( void *cogmres_vdata, JXF_Int   level );
JXF_Int jxf_COGMRESGetLogging( void *cogmres_vdata, JXF_Int * level );
JXF_Int jxf_COGMRESGetNumIterations( void *cogmres_vdata, JXF_Int  *num_iterations );
JXF_Int jxf_COGMRESGetConverged( void *cogmres_vdata, JXF_Int  *converged );
JXF_Int jxf_COGMRESGetFinalRelativeResidualNorm( void   *cogmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int jxf_COGMRESSetModifyPC( void *cogmres_vdata, JXF_Int (*modify_pc)(void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm) );

/* /csrc/krylov/pgmres.c */
//JXF_Int JXF_ParCSRGMRESCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int JXF_GMRESCreate( MPI_Comm comm, JXF_Solver *solver );
//JXF_Int JXF_ParCSRGMRESDestroy( JXF_Solver solver );
JXF_Int JXF_GMRESDestroy( JXF_Solver solver );
JXF_Int JXF_GMRESSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_GMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_GMRESAdaptiveSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_GMRESSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator ); /* peghoty 2011/09/04 */ 
JXF_Int JXF_GMRESSetKDim( JXF_Solver solver, JXF_Int k_dim );
JXF_Int JXF_GMRESGetKDim( JXF_Solver solver, JXF_Int *k_dim );
JXF_Int JXF_GMRESSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_GMRESGetTol( JXF_Solver solver, JXF_Real *tol );
JXF_Int JXF_GMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol );
JXF_Int JXF_GMRESGetAbsoluteTol( JXF_Solver solver, JXF_Real *a_tol );
JXF_Int JXF_GMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol );
JXF_Int JXF_GMRESGetConvergenceFactorTol( JXF_Solver solver, JXF_Real *cf_tol );
JXF_Int JXF_GMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter );
JXF_Int JXF_GMRESGetMinIter( JXF_Solver solver, JXF_Int *min_iter );
JXF_Int JXF_GMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_GMRESGetMaxIter( JXF_Solver solver, JXF_Int *max_iter );
JXF_Int JXF_GMRESSetStopCrit( JXF_Solver solver, JXF_Int stop_crit );
JXF_Int JXF_GMRESGetStopCrit( JXF_Solver solver, JXF_Int *stop_crit );
JXF_Int JXF_GMRESSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria ); /* peghoty 2010/06/23 */
JXF_Int JXF_GMRESSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold ); /* peghoty 2010/06/23 */
JXF_Int JXF_GMRESSetResDownZeroThreshold( JXF_Solver solver, JXF_Real resdown_0_threshold ); /* Yue Xiaoqiang 2014/03/24 */
JXF_Int JXF_GMRESSetConvFacThresholdTwo( JXF_Solver solver, JXF_Real convfac_threshold_2 ); /* Yue Xiaoqiang 2014/03/24 */
JXF_Int JXF_GMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted ); /* peghoty 2011/11/08 */
JXF_Int JXF_GMRESSetRelChange( JXF_Solver solver, JXF_Int rel_change );
JXF_Int JXF_GMRESGetRelChange( JXF_Solver solver, JXF_Int *rel_change );
JXF_Int JXF_GMRESGetLastPrecondType( JXF_Solver solver, JXF_Int *last_precond_type ); /* Yue Xiaoqiang 2014/10/24 */
JXF_Int 
JXF_GMRESSetPrecond( JXF_Solver          solver,
                    JXF_PtrToSolverFcn  precond,
                    JXF_PtrToSolverFcn  precond_setup,
                    JXF_Solver          precond_solver );                         
JXF_Int JXF_GMRESGetPrecond( JXF_Solver solver, JXF_Solver *precond_data_ptr );
JXF_Int JXF_GMRESSetPrintLevel( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_GMRESGetPrintLevel( JXF_Solver solver, JXF_Int *level );
JXF_Int JXF_GMRESGetCompAMGILUBcoIter( JXF_Solver solver, JXF_Int *amg_iter, JXF_Int *ilu_iter, JXF_Int *bco_iter );
JXF_Int JXF_GMRESSetLogging( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_GMRESGetLogging( JXF_Solver solver, JXF_Int *level );
JXF_Int JXF_GMRESGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_GMRESGetConverged( JXF_Solver solver, JXF_Int  *converged );
JXF_Int JXF_GMRESGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm );
JXF_Int JXF_GMRESGetResidual( JXF_Solver solver, void **residual );
void *jxf_ParCSRGMRESCreate( MPI_Comm comm );  /* peghoty 2011/10/26 */
jxf_GMRESFunctions *
jxf_GMRESFunctionsCreate(
   char * (*CAlloc)             ( size_t count, size_t elt_size ),
   JXF_Int    (*Free)               ( char *ptr ),
   JXF_Int    (*CommInfo)           ( void *A, JXF_Int *my_id, JXF_Int *num_procs ),
   void * (*CreateVector)       ( void *vector ),
   void * (*CreateVectorArray)  ( JXF_Int size, void *vectors ),
   JXF_Int    (*DestroyVector)      ( void *vector ),
   void * (*MatvecCreate)       ( void *A, void *x ),
   JXF_Int    (*Matvec)             ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y ),
   JXF_Int    (*MatvecDestroy)      ( void *matvec_data ),
   JXF_Real (*InnerProd)          ( void *x, void *y ),
   JXF_Int    (*CopyVector)         ( void *x, void *y ),
   JXF_Int    (*ClearVector)        ( void *x ),
   JXF_Int    (*ScaleVector)        ( JXF_Real alpha, void *x ),
   JXF_Int    (*Axpy)               ( JXF_Real alpha, void *x, void *y ),
   JXF_Int    (*PrecondSetup)       ( void *vdata, void *A ),
   JXF_Int    (*Precond)            ( void *vdata, void *A, void *b, void *x )
);
void *jxf_GMRESCreate( jxf_GMRESFunctions *gmres_functions );
JXF_Int jxf_GMRESDestroy( void *gmres_vdata );
JXF_Int jxf_GMRESGetResidual( void *gmres_vdata, void **residual );
JXF_Int jxf_GMRESSetup( void *gmres_vdata, void *matA, void *vecB, void *vecX );
JXF_Int jxf_GMRESSetPreOperator( void *gmres_vdata, void *preOperator ); /* peghoty 2011/09/04 */
JXF_Int jxf_GMRESSetKDim( void *gmres_vdata, JXF_Int k_dim );
JXF_Int jxf_GMRESGetKDim( void *gmres_vdata, JXF_Int *k_dim );
JXF_Int jxf_GMRESSetTol( void *gmres_vdata, JXF_Real tol );
JXF_Int jxf_GMRESGetTol( void *gmres_vdata, JXF_Real *tol );
JXF_Int jxf_GMRESSetAbsoluteTol( void *gmres_vdata, JXF_Real a_tol );
JXF_Int jxf_GMRESGetAbsoluteTol( void *gmres_vdata, JXF_Real *a_tol );
JXF_Int jxf_GMRESSetConvergenceFactorTol( void *gmres_vdata, JXF_Real cf_tol );
JXF_Int jxf_GMRESGetConvergenceFactorTol( void *gmres_vdata, JXF_Real *cf_tol );
JXF_Int jxf_GMRESSetMinIter( void *gmres_vdata, JXF_Int min_iter );
JXF_Int jxf_GMRESGetMinIter( void *gmres_vdata, JXF_Int *min_iter );
JXF_Int jxf_GMRESSetMaxIter( void *gmres_vdata, JXF_Int max_iter );
JXF_Int jxf_GMRESGetMaxIter( void *gmres_vdata, JXF_Int *max_iter );
JXF_Int jxf_GMRESSetRelChange( void *gmres_vdata, JXF_Int rel_change );
JXF_Int jxf_GMRESGetRelChange( void *gmres_vdata, JXF_Int *rel_change );
JXF_Int jxf_GMRESGetLastPrecondType( void *gmres_vdata, JXF_Int *last_precond_type ); /* Yue Xiaoqiang 2014/10/24 */
JXF_Int jxf_GMRESSetStopCrit( void *gmres_vdata, JXF_Int stop_crit );
JXF_Int jxf_GMRESGetStopCrit( void *gmres_vdata, JXF_Int *stop_crit );
JXF_Int jxf_GMRESSetConvCriteria( void *gmres_vdata, JXF_Int conv_criteria ); /* peghoty 2010/06/23 */
JXF_Int jxf_GMRESSetConvFacThreshold( void *gmres_vdata, JXF_Real convfac_threshold ); /* peghoty 2010/06/23 */
JXF_Int jxf_GMRESSetResDownZeroThreshold( void *gmres_vdata, JXF_Real resdown_0_threshold ); /* Yue Xiaoqiang 2014/03/24 */
JXF_Int jxf_GMRESSetConvFacThresholdTwo( void *gmres_vdata, JXF_Real convfac_threshold_2 ); /* Yue Xiaoqiang 2014/03/24 */
JXF_Int jxf_GMRESSetIsCheckRestarted( void *gmres_vdata, JXF_Int is_check_restarted ); /* peghoty 2011/11/08 */
JXF_Int jxf_GMRESSetPrecond( void *gmres_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int jxf_GMRESGetPrecond( void *gmres_vdata, JXF_Solver *precond_data_ptr );
JXF_Int jxf_GMRESSetPrintLevel( void *gmres_vdata, JXF_Int level );
JXF_Int jxf_GMRESGetPrintLevel( void *gmres_vdata, JXF_Int *level );
JXF_Int jxf_GMRESGetCompAMGILUBcoIter( void *gmres_vdata, JXF_Int *amg_iter, JXF_Int *ilu_iter, JXF_Int *bco_iter );
JXF_Int jxf_GMRESSetLogging( void *gmres_vdata, JXF_Int level );
JXF_Int jxf_GMRESGetLogging( void *gmres_vdata, JXF_Int *level );
JXF_Int jxf_GMRESGetNumIterations( void *gmres_vdata, JXF_Int *num_iterations );
JXF_Int jxf_GMRESGetConverged( void *gmres_vdata, JXF_Int *converged );
JXF_Int jxf_GMRESGetFinalRelativeResidualNorm( void *gmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int jxf_GMRESSolve( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX );
JXF_Int jxf_GMRESAdaptiveSolve( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX );
JXF_Int jxf_GMRESSolveDefault( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX );

/* /csrc/krylov/pflexgmres.c */
JXF_Int JXF_ParCSRFlexGMRESCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int JXF_ParCSRFlexGMRESDestroy( JXF_Solver solver );
JXF_Int JXF_FlexGMRESSetup( JXF_Solver solver, JXF_Matrix A, JXF_Vector b, JXF_Vector x );
JXF_Int JXF_FlexGMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix A, JXF_Vector b, JXF_Vector x );
JXF_Int JXF_FlexGMRESSetKDim( JXF_Solver solver, JXF_Int k_dim );
JXF_Int JXF_FlexGMRESSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_FlexGMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol );
JXF_Int JXF_FlexGMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol );
JXF_Int JXF_FlexGMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter );
JXF_Int JXF_FlexGMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_FlexGMRESSetPrecond( JXF_Solver solver,
                               JXF_PtrToSolverFcn precond,
                               JXF_PtrToSolverFcn precond_setup,
                               JXF_Solver precond_solver );
JXF_Int JXF_FlexGMRESSetPrintLevel( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_FlexGMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted );
JXF_Int JXF_FlexGMRESSetLogging( JXF_Solver solver, JXF_Int level );
JXF_Int JXF_FlexGMRESGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_FlexGMRESGetConverged( JXF_Solver solver, JXF_Int *converged );
JXF_Int JXF_FlexGMRESGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm );
JXF_Int JXF_FlexGMRESGetResidual( JXF_Solver solver, void *residual );
JXF_Int JXF_FlexGMRESSetModifyPC( JXF_Solver  solver, JXF_Int (*modify_pc)(JXF_Solver, JXF_Int, JXF_Real) );
jxf_FlexGMRESFunctions *jxf_FlexGMRESFunctionsCreate(
   char *       (*CAlloc)        ( size_t count, size_t elt_size ),
   JXF_Int    (*Free)          ( char *ptr ),
   JXF_Int    (*CommInfo)      ( void  *A, JXF_Int   *my_id, JXF_Int   *num_procs ),
   void *       (*CreateVector)  ( void *vector ),
   void *       (*CreateVectorArray)  ( JXF_Int size, void *vectors ),
   JXF_Int    (*DestroyVector) ( void *vector ),
   void *       (*MatvecCreate)  ( void *A, void *x ),
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y ),
   JXF_Int    (*MatvecDestroy) ( void *matvec_data ),
   JXF_Real   (*InnerProd)     ( void *x, void *y ),
   JXF_Int    (*CopyVector)    ( void *x, void *y ),
   JXF_Int    (*ClearVector)   ( void *x ),
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x ),
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y ),
   JXF_Int    (*PrecondSetup)  ( void *vdata, void *A ),
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x ) );
void *jxf_FlexGMRESCreate( jxf_FlexGMRESFunctions *fgmres_functions );
JXF_Int jxf_FlexGMRESDestroy( void *fgmres_vdata );
JXF_Int jxf_FlexGMRESGetResidual( void *fgmres_vdata, void **residual );
JXF_Int jxf_FlexGMRESSetup( void *fgmres_vdata, void *A, void *b, void *x );
JXF_Int jxf_FlexGMRESSolve( void *fgmres_vdata, void *preOperator, void *A, void *b, void *x );
JXF_Int jxf_FlexGMRESSetKDim( void *fgmres_vdata, JXF_Int k_dim );
JXF_Int jxf_FlexGMRESSetTol( void *fgmres_vdata, JXF_Real tol );
JXF_Int jxf_FlexGMRESSetAbsoluteTol( void *fgmres_vdata, JXF_Real a_tol );
JXF_Int jxf_FlexGMRESSetConvergenceFactorTol( void *fgmres_vdata, JXF_Real cf_tol );
JXF_Int jxf_FlexGMRESSetMinIter( void *fgmres_vdata, JXF_Int min_iter );
JXF_Int jxf_FlexGMRESSetMaxIter( void *fgmres_vdata, JXF_Int max_iter );
JXF_Int jxf_FlexGMRESSetStopCrit( void *fgmres_vdata, JXF_Int stop_crit );
JXF_Int jxf_FlexGMRESSetPrecond( void *fgmres_vdata,
                               JXF_Int (*precond)(void*,void*,void*,void*),
                               JXF_Int (*precond_setup)(void*,void*),
                               void *precond_data );
JXF_Int jxf_FlexGMRESSetPrintLevel( void *fgmres_vdata, JXF_Int level );
JXF_Int jxf_FlexGMRESSetIsCheckRestarted( void *fgmres_vdata, JXF_Int is_check_restarted );
JXF_Int jxf_FlexGMRESSetLogging( void *fgmres_vdata, JXF_Int level );
JXF_Int jxf_FlexGMRESGetNumIterations( void *fgmres_vdata, JXF_Int *num_iterations );
JXF_Int jxf_FlexGMRESGetConverged( void *fgmres_vdata, JXF_Int *converged );
JXF_Int jxf_FlexGMRESGetFinalRelativeResidualNorm( void *fgmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int jxf_FlexGMRESSetModifyPC( void *fgmres_vdata, JXF_Int (*modify_pc)(void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm) );
JXF_Int jxf_FlexGMRESModifyPCDefault( void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm );

/* /csrc/krylov/pbicgstab.c */
//JXF_Int JXF_ParCSRBiCGSTABCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int JXF_BiCGSTABCreate( MPI_Comm comm, JXF_Solver *solver );
//JXF_Int JXF_ParCSRBiCGSTABDestroy( JXF_Solver solver );
JXF_Int JXF_BiCGSTABDestroy( JXF_Solver solver );
JXF_Int JXF_BiCGSTABSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_BiCGSTABSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX );
JXF_Int JXF_BiCGSTABSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator ); /* peghoty 2011/09/04 */
JXF_Int JXF_BiCGSTABSetTol( JXF_Solver solver, JXF_Real tol );
JXF_Int JXF_BiCGSTABSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol );
JXF_Int JXF_BiCGSTABSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol );
JXF_Int JXF_BiCGSTABSetMinIter( JXF_Solver solver, JXF_Int min_iter );
JXF_Int JXF_BiCGSTABSetMaxIter( JXF_Solver solver, JXF_Int max_iter );
JXF_Int JXF_BiCGSTABSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria ); /* peghoty 2010/03/20 */
JXF_Int JXF_BiCGSTABSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold ); /* peghoty 2010/03/22 */
JXF_Int
JXF_BiCGSTABSetPrecond( JXF_Solver         solver,
                       JXF_PtrToSolverFcn precond,
                       JXF_PtrToSolverFcn precond_setup,
                       JXF_Solver         precond_solver );
JXF_Int JXF_BiCGSTABGetPrecond( JXF_Solver solver, JXF_Solver *precond_data_ptr );
JXF_Int JXF_BiCGSTABSetLogging( JXF_Solver solver, JXF_Int logging );
JXF_Int JXF_BiCGSTABSetPrintLevel( JXF_Solver solver, JXF_Int print_level );
JXF_Int JXF_BiCGSTABGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations );
JXF_Int JXF_BiCGSTABGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm );
JXF_Int JXF_BiCGSTABGetResidual( JXF_Solver solver, void **residual );
jxf_BiCGSTABFunctions *
jxf_BiCGSTABFunctionsCreate(
   void * (*CreateVector)  ( void *vvector ),
   JXF_Int    (*DestroyVector) ( void *vvector ),
   void * (*MatvecCreate)  ( void *A, void *x ),
   JXF_Int    (*Matvec)        ( void *matvec_data, JXF_Real alpha, void *A, void *x, JXF_Real beta, void *y ),
   JXF_Int    (*MatvecDestroy) ( void *matvec_data ),
   JXF_Real (*InnerProd)     ( void *x, void *y ),
   JXF_Int    (*CopyVector)    ( void *x, void *y ),
   JXF_Int    (*ClearVector)   ( void *x ),
   JXF_Int    (*ScaleVector)   ( JXF_Real alpha, void *x ),
   JXF_Int    (*Axpy)          ( JXF_Real alpha, void *x, void *y ),
   JXF_Int    (*CommInfo)      ( void *A, JXF_Int *my_id, JXF_Int *num_procs ),
   JXF_Int    (*PrecondSetup)  (  void *vdata, void *A ),
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
);
void *jxf_BiCGSTABCreate( jxf_BiCGSTABFunctions * bicgstab_functions );
JXF_Int jxf_BiCGSTABDestroy( void *bicgstab_vdata );
JXF_Int jxf_BiCGSTABSetup( void *bicgstab_vdata, void *matA, void *vecB, void *vecX );
JXF_Int jxf_BiCGSTABSetPreOperator( void *bicgstab_vdata, void *preOperator ); /* peghoty 2011/09/04 */
JXF_Int jxf_BiCGSTABSetTol( void *bicgstab_vdata, JXF_Real tol );
JXF_Int jxf_BiCGSTABSetAbsoluteTol( void *bicgstab_vdata, JXF_Real a_tol );
JXF_Int jxf_BiCGSTABSetConvergenceFactorTol( void *bicgstab_vdata, JXF_Real cf_tol );
JXF_Int jxf_BiCGSTABSetMinIter( void *bicgstab_vdata, JXF_Int min_iter );
JXF_Int jxf_BiCGSTABSetMaxIter( void *bicgstab_vdata, JXF_Int max_iter );
JXF_Int jxf_BiCGSTABSetConvCriteria( void *bicgstab_vdata, JXF_Int conv_criteria ); /* peghoty 2010/03/20 */
JXF_Int jxf_BiCGSTABSetConvFacThreshold( void *bicgstab_vdata, JXF_Real convfac_threshold ); /* peghoty 2010/03/22 */
JXF_Int jxf_BiCGSTABSetPrecond( void *bicgstab_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int jxf_BiCGSTABGetPrecond( void *bicgstab_vdata, JXF_Solver *precond_data_ptr );
JXF_Int jxf_BiCGSTABSetLogging( void *bicgstab_vdata, JXF_Int logging );   
JXF_Int jxf_BiCGSTABSetPrintLevel( void *bicgstab_vdata, JXF_Int print_level );
JXF_Int jxf_BiCGSTABGetConverged( void *bicgstab_vdata, JXF_Int *converged ); 
JXF_Int jxf_BiCGSTABGetNumIterations( void *bicgstab_vdata, JXF_Int *num_iterations );
JXF_Int jxf_BiCGSTABGetFinalRelativeResidualNorm( void *bicgstab_vdata, JXF_Real *relative_residual_norm );
JXF_Int jxf_BiCGSTABGetResidual( void *bicgstab_vdata, void **residual );
JXF_Int jxf_BiCGSTABSolve( void *bicgstab_vdata, void *preOperator, void *matA, void *vecB, void *vecX );

#endif
