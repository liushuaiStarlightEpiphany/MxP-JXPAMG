//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  pcg.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"
#include "jxf_krylov.h"

/*!
 * \fn JXF_Int JXF_ParCSRPCGCreate
 */ 
JXF_Int
JXF_PCGCreate( MPI_Comm comm, JXF_Solver *solver )
{
   /* The function names with a PCG in them are in
      parcsr_ls/pcg_par.c .  These functions do rather little -
      e.g., cast to the correct type - before calling something else.
      These names should be called, e.g., jxf_struct_Free, to reduce the
      chance of name conflicts. */
   jxf_PCGFunctions * pcg_functions =
      jxf_PCGFunctionsCreate(
         jxf_CAlloc, 
         jxf_ParKrylovFree, 
         jxf_ParKrylovCommInfo,
         jxf_ParKrylovCreateVector,
         jxf_ParKrylovDestroyVector, 
         jxf_ParKrylovMatvecCreate,
         jxf_ParKrylovMatvec, 
         jxf_ParKrylovMatvecDestroy,
         jxf_ParKrylovInnerProd, 
         jxf_ParKrylovCopyVector,
         jxf_ParKrylovClearVector,
         jxf_ParKrylovScaleVector, 
         jxf_ParKrylovAxpy,
         jxf_ParKrylovIdentitySetup, 
         jxf_ParKrylovIdentity );

   *solver = ( (JXF_Solver) jxf_PCGCreate( pcg_functions ) );
   
   if (!solver) 
   {
      jxf_error_in_arg(2);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PCGFunctionsCreate
 */
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
)
{
   jxf_PCGFunctions *pcg_functions;
   pcg_functions = (jxf_PCGFunctions *)CAlloc(1, sizeof(jxf_PCGFunctions));

   pcg_functions->CAlloc = CAlloc;
   pcg_functions->Free = Free;
   pcg_functions->CommInfo = CommInfo;
   pcg_functions->CreateVector = CreateVector;
   pcg_functions->DestroyVector = DestroyVector;
   pcg_functions->MatvecCreate = MatvecCreate;
   pcg_functions->Matvec = Matvec;
   pcg_functions->MatvecDestroy = MatvecDestroy;
   pcg_functions->InnerProd = InnerProd;
   pcg_functions->CopyVector = CopyVector;
   pcg_functions->ClearVector = ClearVector;
   pcg_functions->ScaleVector = ScaleVector;
   pcg_functions->Axpy = Axpy;
   /* default preconditioner must be set here but can be changed later... */
   pcg_functions->precond_setup = PrecondSetup;
   pcg_functions->precond       = Precond;

   return pcg_functions;
}

/*!
 * \fn JXF_Int jxf_PCGCreate
 */
void *
jxf_PCGCreate( jxf_PCGFunctions *pcg_functions )
{
   jxf_PCGData *pcg_data;

   pcg_data = jxf_CTAllocF(jxf_PCGData, 1, pcg_functions);

   pcg_data -> functions = pcg_functions;

   /* set defaults */
   (pcg_data -> tol)                = 1.0e-06;
   (pcg_data -> atolf)              = 0.0;
   (pcg_data -> cf_tol)             = 0.0;
   (pcg_data -> a_tol)              = 0.0;
   (pcg_data -> max_iter)           = 1000;
   (pcg_data -> two_norm)           = 0;
   (pcg_data -> rel_change)         = 0;
   (pcg_data -> recompute_residual) = 0;
   (pcg_data -> stop_crit)          = 0;
   (pcg_data -> conv_criteria)      = 0;    // peghoty 2010/06/23
   (pcg_data -> convfac_threshold)  = 1.0e4;  // peghoty 2010/06/23
   (pcg_data -> converged)          = 0;
   (pcg_data -> owns_matvec_data )  = 1;
   (pcg_data -> matvec_data)        = NULL;
   (pcg_data -> precond_data)       = NULL;
   (pcg_data -> print_level)        = 0;
   (pcg_data -> logging)            = 0;
   (pcg_data -> norms)              = NULL;
   (pcg_data -> rel_norms)          = NULL;
   (pcg_data -> B)                  = NULL; // peghoty 2011/10/08
   (pcg_data -> p)                  = NULL;
   (pcg_data -> s)                  = NULL;
   (pcg_data -> r)                  = NULL;

   return (void *) pcg_data;
}

/*!
 * \fn JXF_Int JXF_PCGSetup
 */
JXF_Int 
JXF_PCGSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_PCGSetup( solver, matA, vecB, vecX ) );
}

/*!
 * \fn JXF_Int JXF_PCGSolve
 */
JXF_Int 
JXF_PCGSolve( JXF_Solver  solver, 
             JXF_Matrix  preOperator,
             JXF_Matrix  matA, 
             JXF_Vector  vecB, 
             JXF_Vector  vecX )
{
   return( jxf_PCGSolve( (void *) solver, 
                        (void *) preOperator,
                        (void *) matA, 
                        (void *) vecB, 
                        (void *) vecX ) );
}

/*!
 * \fn JXF_Int JXF_PCGSetPreOperator
 * \brief Set Preconditioner operator.
 * \author peghoty
 * \date  2011/09/04 
 */
JXF_Int
JXF_PCGSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator )
{
   return( jxf_PCGSetPreOperator( (void *) solver, (void *) preOperator ) );
}

/*!
 * \fn JXF_Int JXF_PCGSetXXX
 */

JXF_Int
JXF_PCGSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_PCGSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_PCGGetTol( JXF_Solver solver, JXF_Real *tol )
{
   return( jxf_PCGGetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_PCGSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol )
{
   return( jxf_PCGSetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_PCGGetAbsoluteTol( JXF_Solver solver, JXF_Real *a_tol )
{
   return( jxf_PCGGetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_PCGSetAbsoluteTolFactor( JXF_Solver solver, JXF_Real abstolf )
{
   return( jxf_PCGSetAbsoluteTolFactor( (void *) solver, abstolf ) );
}

JXF_Int
JXF_PCGGetAbsoluteTolFactor( JXF_Solver solver, JXF_Real *abstolf )
{
   return( jxf_PCGGetAbsoluteTolFactor( (void *) solver, abstolf ) );
}

JXF_Int
JXF_PCGSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol )
{
   return jxf_PCGSetConvergenceFactorTol( (void *) solver, cf_tol );
}

JXF_Int
JXF_PCGGetConvergenceFactorTol( JXF_Solver solver, JXF_Real * cf_tol )
{
   return jxf_PCGGetConvergenceFactorTol( (void *) solver, cf_tol );
}

JXF_Int
JXF_PCGSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_PCGSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_PCGGetMaxIter( JXF_Solver solver, JXF_Int *max_iter )
{
   return( jxf_PCGGetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_PCGSetStopCrit( JXF_Solver solver, JXF_Int stop_crit )
{
   return( jxf_PCGSetStopCrit( (void *) solver, stop_crit ) );
}

JXF_Int
JXF_PCGGetStopCrit( JXF_Solver solver, JXF_Int *stop_crit )
{
   return( jxf_PCGGetStopCrit( (void *) solver, stop_crit ) );
}

/*!
 * \fn JXF_Int JXF_PCGSetConvCriteria
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
JXF_PCGSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria )
{
   return( jxf_PCGSetConvCriteria( (void *) solver, conv_criteria ) );
}

/*!
 * \fn JXF_Int JXF_PCGSetConvFacThreshold
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
JXF_PCGSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold )
{
   return( jxf_PCGSetConvFacThreshold( (void *) solver, convfac_threshold ) );
}

JXF_Int
JXF_PCGSetTwoNorm( JXF_Solver solver, JXF_Int two_norm )
{
   return( jxf_PCGSetTwoNorm( (void *) solver, two_norm ) );
}

JXF_Int
JXF_PCGGetTwoNorm( JXF_Solver solver, JXF_Int * two_norm )
{
   return( jxf_PCGGetTwoNorm( (void *) solver, two_norm ) );
}

JXF_Int
JXF_PCGSetRelChange( JXF_Solver solver, JXF_Int rel_change )
{
   return( jxf_PCGSetRelChange( (void *) solver, rel_change ) );
}

JXF_Int
JXF_PCGGetRelChange( JXF_Solver solver, JXF_Int *rel_change )
{
   return( jxf_PCGGetRelChange( (void *) solver, rel_change ) );
}

JXF_Int
JXF_PCGSetRecomputeResidual( JXF_Solver solver, JXF_Int recompute_residual )
{
   return( jxf_PCGSetRecomputeResidual( (void *) solver, recompute_residual ) );
}

JXF_Int
JXF_PCGGetRecomputeResidual( JXF_Solver solver, JXF_Int *recompute_residual )
{
   return( jxf_PCGGetRecomputeResidual( (void *) solver, recompute_residual ) );
}

/*!
 * \fn JXF_Int jxf_PCGSetConvCriteria
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
jxf_PCGSetConvCriteria( void *pcg_vdata, JXF_Int conv_criteria  )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> conv_criteria) = conv_criteria;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PCGSetConvFacThreshold
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
jxf_PCGSetConvFacThreshold( void *pcg_vdata, JXF_Real convfac_threshold )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> convfac_threshold) = convfac_threshold;
   return jxf_error_flag;
}

JXF_Int
JXF_PCGSetPrecond( JXF_Solver         solver,
                  JXF_PtrToSolverFcn precond,
                  JXF_PtrToSolverFcn precond_setup,
                  JXF_Solver         precond_solver )
{
   return( jxf_PCGSetPrecond( (void *) solver, precond, precond_setup, (void *) precond_solver ) );
}

JXF_Int
JXF_PCGGetPrecond( JXF_Solver  solver, JXF_Solver *precond_data_ptr )
{
   return( jxf_PCGGetPrecond( (void *) solver, (JXF_Solver *) precond_data_ptr ) );
}

/*-----------------------------------------------------------------------------
 * JXF_PCGSetLogging, JXF_PCGGetLogging
 * SetLogging sets both the print and log level, for backwards compatibility.
 * Soon the SetPrintLevel call should be deleted.
 *-----------------------------------------------------------------------------*/

JXF_Int
JXF_PCGSetLogging( JXF_Solver solver, JXF_Int level )
{
   return ( jxf_PCGSetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_PCGGetLogging( JXF_Solver solver, JXF_Int *level )
{
   return ( jxf_PCGGetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_PCGSetPrintLevel( JXF_Solver solver, JXF_Int level )
{
   return( jxf_PCGSetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_PCGGetPrintLevel( JXF_Solver solver, JXF_Int *level )
{
   return( jxf_PCGGetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_PCGGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations )
{
   return( jxf_PCGGetNumIterations( (void *) solver, num_iterations ) );
}

JXF_Int
JXF_PCGGetConverged( JXF_Solver solver, JXF_Int *converged )
{
   return( jxf_PCGGetConverged( (void *) solver, converged ) );
}

JXF_Int
JXF_PCGGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm )
{
   return( jxf_PCGGetFinalRelativeResidualNorm( (void *) solver, norm ) );
}

JXF_Int 
JXF_PCGGetResidual( JXF_Solver solver, void **residual )
{
   /* returns a pointer to the residual vector */
   return jxf_PCGGetResidual( (void *) solver, residual );
}

JXF_Int 
JXF_PCGDestroy( JXF_Solver solver )
{
   return( jxf_PCGDestroy( (void *) solver ) );
}

/*!
 * \fn JXF_Int jxf_PCGGetResidual
 */
JXF_Int 
jxf_PCGGetResidual( void *pcg_vdata, void **residual )
{
   /* returns a pointer to the residual vector */
   jxf_PCGData  *pcg_data = pcg_vdata;
   *residual = pcg_data->r;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PCGDestroy
 */
JXF_Int
jxf_PCGDestroy( void *pcg_vdata )
{
   jxf_PCGData      *pcg_data      = pcg_vdata;
   jxf_PCGFunctions *pcg_functions = pcg_data->functions;

   if (pcg_data)
   {
      if ( (pcg_data -> norms) != NULL )
      {
         jxf_TFreeF( pcg_data -> norms, pcg_functions );
         pcg_data -> norms = NULL;
      } 
      if ( (pcg_data -> rel_norms) != NULL )
      {
         jxf_TFreeF( pcg_data -> rel_norms, pcg_functions );
         pcg_data -> rel_norms = NULL;
      }
      if ( pcg_data -> matvec_data != NULL && pcg_data->owns_matvec_data )
      {
         (*(pcg_functions->MatvecDestroy))(pcg_data -> matvec_data);
         pcg_data -> matvec_data = NULL;
      }
      if ( pcg_data -> p != NULL )
      {
         (*(pcg_functions->DestroyVector))(pcg_data -> p);
         pcg_data -> p = NULL;
      }
      if ( pcg_data -> s != NULL )
      {
         (*(pcg_functions->DestroyVector))(pcg_data -> s);
         pcg_data -> s = NULL;
      }
      if ( pcg_data -> r != NULL )
      {
         (*(pcg_functions->DestroyVector))(pcg_data -> r);
         pcg_data -> r = NULL;
      }
      jxf_TFreeF( pcg_data, pcg_functions );
      jxf_TFreeF( pcg_functions, pcg_functions );
   }

   return(jxf_error_flag);
}

/*!
 * \fn JXF_Int jxf_PCGSetup
 */
JXF_Int
jxf_PCGSetup( void *pcg_vdata, void *matA, void *vecB, void *vecX )
{
   jxf_PCGData      *pcg_data          = pcg_vdata;
   jxf_PCGFunctions *pcg_functions     = pcg_data->functions;
   JXF_Int              max_iter          = (pcg_data -> max_iter);

   (pcg_data -> A) = matA;

  /*---------------------------------------------------
   * The arguments for CreateVector are important to
   * maintain consistency between the setup and
   * compute phases of matvec and the preconditioner.
   *--------------------------------------------------*/

   if ( pcg_data -> p != NULL )
   {
      (*(pcg_functions->DestroyVector))(pcg_data -> p);
   }   
   (pcg_data -> p) = (*(pcg_functions->CreateVector))(vecX);

   if ( pcg_data -> s != NULL )
   {
      (*(pcg_functions->DestroyVector))(pcg_data -> s);
   }
   (pcg_data -> s) = (*(pcg_functions->CreateVector))(vecX);

   if ( pcg_data -> r != NULL )
   {
      (*(pcg_functions->DestroyVector))(pcg_data -> r);
   }
   (pcg_data -> r) = (*(pcg_functions->CreateVector))(vecB);

   if ( pcg_data -> matvec_data != NULL && pcg_data->owns_matvec_data )
   {
      (*(pcg_functions->MatvecDestroy))(pcg_data -> matvec_data);
   }
   (pcg_data -> matvec_data) = (*(pcg_functions->MatvecCreate))(matA, vecX);

  /*-------------------------------------------------------------------
   * Allocate space for log info
   *------------------------------------------------------------------*/
   
   if ( (pcg_data->logging) > 0 || (pcg_data->print_level) > 0 ) 
   {
      if ( (pcg_data -> norms) != NULL )
      {
         jxf_TFreeF( pcg_data -> norms, pcg_functions );
      }
      (pcg_data -> norms) = jxf_CTAllocF( JXF_Real, max_iter + 1, pcg_functions);

      if ( (pcg_data -> rel_norms) != NULL )
      {
         jxf_TFreeF( pcg_data -> rel_norms, pcg_functions );
      }
      (pcg_data -> rel_norms) = jxf_CTAllocF( JXF_Real, max_iter + 1, pcg_functions );
   }

   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetTol( void *pcg_vdata, JXF_Real tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> tol) = tol;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetTol( void *pcg_vdata, JXF_Real *tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *tol = (pcg_data -> tol);
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PCGSetPreOperator
 * \author peghoty
 * \date 2011/09/04 
 */
JXF_Int
jxf_PCGSetPreOperator( void *pcg_vdata, void *preOperator )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> B) = preOperator;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetAbsoluteTol( void *pcg_vdata, JXF_Real a_tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> a_tol) = a_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetAbsoluteTol( void *pcg_vdata, JXF_Real *a_tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *a_tol = (pcg_data -> a_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetAbsoluteTolFactor( void *pcg_vdata, JXF_Real atolf )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> atolf) = atolf;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetAbsoluteTolFactor( void *pcg_vdata, JXF_Real *atolf )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *atolf = (pcg_data -> atolf);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetConvergenceFactorTol( void *pcg_vdata, JXF_Real cf_tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> cf_tol) = cf_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetConvergenceFactorTol( void *pcg_vdata, JXF_Real *cf_tol )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *cf_tol = (pcg_data -> cf_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetMaxIter( void *pcg_vdata, JXF_Int max_iter )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetMaxIter( void *pcg_vdata, JXF_Int *max_iter )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *max_iter = (pcg_data -> max_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetTwoNorm( void *pcg_vdata, JXF_Int two_norm )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> two_norm) = two_norm;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetTwoNorm( void *pcg_vdata, JXF_Int *two_norm )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *two_norm = (pcg_data -> two_norm);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetRelChange( void *pcg_vdata, JXF_Int rel_change )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> rel_change) = rel_change;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetRelChange( void *pcg_vdata, JXF_Int *rel_change )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *rel_change = (pcg_data -> rel_change);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetRecomputeResidual( void *pcg_vdata, JXF_Int recompute_residual )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> recompute_residual) = recompute_residual;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetRecomputeResidual( void *pcg_vdata, JXF_Int *recompute_residual )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *recompute_residual = (pcg_data -> recompute_residual);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetStopCrit( void *pcg_vdata, JXF_Int stop_crit )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> stop_crit) = stop_crit;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetStopCrit( void *pcg_vdata, JXF_Int *stop_crit )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *stop_crit = (pcg_data -> stop_crit);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetPrecond( void *pcg_vdata, JXF_Solver *precond_data_ptr )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *precond_data_ptr = (JXF_Solver)(pcg_data -> precond_data);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetPrecond( void  *pcg_vdata,
                  JXF_Int  (*precond)(),
                  JXF_Int  (*precond_setup)(),
                  void  *precond_data )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   jxf_PCGFunctions *pcg_functions = pcg_data->functions;

   (pcg_functions -> precond)       = precond;
   (pcg_functions -> precond_setup) = precond_setup;
   (pcg_data -> precond_data)       = precond_data;
 
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetPrintLevel( void *pcg_vdata, JXF_Int level )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> print_level) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetPrintLevel( void *pcg_vdata, JXF_Int *level )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *level = (pcg_data -> print_level);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGSetLogging( void *pcg_vdata, JXF_Int level )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   (pcg_data -> logging) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetLogging( void *pcg_vdata, JXF_Int *level )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *level = (pcg_data -> logging);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetNumIterations( void *pcg_vdata, JXF_Int *num_iterations )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *num_iterations = (pcg_data -> num_iterations);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetConverged( void *pcg_vdata, JXF_Int *converged )
{
   jxf_PCGData *pcg_data = pcg_vdata;
   *converged = (pcg_data -> converged);
   return jxf_error_flag;
}

JXF_Int
jxf_PCGPrintLogging( void *pcg_vdata, JXF_Int myid )
{
   jxf_PCGData *pcg_data       = pcg_vdata;

   JXF_Int         num_iterations = (pcg_data -> num_iterations);
   JXF_Int         print_level    = (pcg_data -> print_level);
   JXF_Real     *norms          = (pcg_data -> norms);
   JXF_Real     *rel_norms      = (pcg_data -> rel_norms);

   JXF_Int         i;

   if (myid == 0)
   {
      if (print_level > 0)
      {
         for (i = 0; i < num_iterations; i ++)
         {
            jxf_printf("Residual norm[%d] = %e   ", i, norms[i]);
            jxf_printf("Relative residual norm[%d] = %e\n", i, rel_norms[i]);
         }
      }
   }
  
   return jxf_error_flag;
}

JXF_Int
jxf_PCGGetFinalRelativeResidualNorm( void *pcg_vdata, JXF_Real *relative_residual_norm )
{
   jxf_PCGData *pcg_data          = pcg_vdata;
   JXF_Real      rel_residual_norm = (pcg_data -> rel_residual_norm);
  *relative_residual_norm = rel_residual_norm;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PCGSolve
 * \note
 * We use the following convergence test as the default (see Ashby, Holst,
 * Manteuffel, and Saylor):
 *
 *       ||e||_A                           ||r||_C
 *       -------  <=  [kappa_A(C*A)]^(1/2) -------  < tol
 *       ||x||_A                           ||b||_C
 *
 * where we let (for the time being) kappa_A(CA) = 1.
 * We implement the test as:
 *
 *       gamma = <C*r,r>/<C*b,b>  <  (tol^2) = eps 
 */
JXF_Int
jxf_PCGSolve( void *pcg_vdata, void *preOperator, void *matA, void *vecB, void *vecX )
{
   jxf_PCGData      *pcg_data      = pcg_vdata;
   jxf_PCGFunctions *pcg_functions = pcg_data -> functions;

   JXF_Real          r_tol              = (pcg_data -> tol);
   JXF_Real          a_tol              = (pcg_data -> a_tol);
   JXF_Real          atolf              = (pcg_data -> atolf);
   JXF_Real          cf_tol             = (pcg_data -> cf_tol);
   JXF_Int             max_iter           = (pcg_data -> max_iter);
   JXF_Int             two_norm           = (pcg_data -> two_norm);
   JXF_Int             rel_change         = (pcg_data -> rel_change);
   JXF_Int             recompute_residual = (pcg_data -> recompute_residual);
   JXF_Int             stop_crit          = (pcg_data -> stop_crit);
   JXF_Real          convfac_threshold = (pcg_data -> convfac_threshold); // peghoty 2010/06/23
   void           *vecP               = (pcg_data -> p);
   void           *vecS               = (pcg_data -> s);
   void           *vecR               = (pcg_data -> r);
   void           *matvec_data        = (pcg_data -> matvec_data);
   JXF_Int           (*precond)()         = (pcg_functions -> precond);
   void           *precond_data       = (pcg_data -> precond_data);
   JXF_Int             print_level        = (pcg_data -> print_level);
   JXF_Int             logging            = (pcg_data -> logging);
   JXF_Real         *norms              = (pcg_data -> norms);
   JXF_Real         *rel_norms          = (pcg_data -> rel_norms);
                
   JXF_Real          alpha, beta;
   JXF_Real          gamma, gamma_old;
   JXF_Real          i_prod = 0.0;
   JXF_Real          bi_prod, eps;
   JXF_Real          pi_prod, xi_prod;
   JXF_Real          ieee_check = 0.;
                
   JXF_Real          i_prod_0 = 0.0;
   JXF_Real          cf_ave_0 = 0.0;
   JXF_Real          cf_ave_1 = 0.0;
   JXF_Real          weight;
   JXF_Real          ratio;

   JXF_Real          guard_zero_residual, sdotp;
   JXF_Int             tentatively_converged = 0; 

   JXF_Int             i = 0;
   JXF_Int             my_id, num_procs;

   JXF_Real  normup       = 0.0;
   JXF_Real  normdown     = 0.0;
   JXF_Real  conv_factor  = 0.0;

   (pcg_data -> converged) = 0;

   (*(pcg_functions->CommInfo))(matA, &my_id, &num_procs);

  /*-----------------------------------------------------------------------
   * With relative change convergence test on, it is possible to attempt
   * another iteration with a zero residual. This causes the parameter
   * alpha to go NaN. The guard_zero_residual parameter is to circumvent
   * this. Perhaps it should be set to something non-zero (but small).
   *-----------------------------------------------------------------------*/

   guard_zero_residual = 0.0;

  /*-----------------------------------------------------------------------
   * Start pcg solve
   *----------------------------------------------------------------------*/

   /* compute eps */
   if (two_norm)
   {
      /* bi_prod = <b,b> */
      bi_prod = (*(pcg_functions->InnerProd))(vecB, vecB);
      if (print_level > 1 && my_id == 0) 
      {
         jxf_printf(" <b,b>: %e\n",bi_prod);
      }
   }
   else
   {
      /* bi_prod = <C*b,b> */
      (*(pcg_functions->ClearVector))(vecP);
      precond(precond_data, preOperator, vecB, vecP);
      bi_prod = (*(pcg_functions->InnerProd))(vecP, vecB);
      if (print_level > 1 && my_id == 0)
      {
         jxf_printf(" <C*b,b>: %e\n",bi_prod);
      }
   }

   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (bi_prod != 0.) ieee_check = bi_prod / bi_prod; /* INF -> NaN conversion */
   if (ieee_check != ieee_check)
   {
      /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
         for ieee_check self-equality works on all IEEE-compliant compilers/
         machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
         by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
         found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
      if (print_level > 0 || logging > 0)
      {
         jxf_printf("\n\nERROR detected by JXFPAMG ...  BEGIN\n");
         jxf_printf("ERROR -- jxf_PCGSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ...  END\n\n\n");
//         break;
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   eps = r_tol*r_tol; /* note: this may be re-assigned below */
   if ( bi_prod > 0.0 ) 
   {
      if ( stop_crit && !rel_change && atolf <= 0 ) 
      {  
         /* pure absolute tolerance */
         eps = eps / bi_prod;
         /* Note: this section is obsolete.  Aside from backwards comatability
            concerns, we could delete the stop_crit parameter and related code,
            using tol & atolf instead. */
      }
      else if ( atolf > 0 )  /* mixed relative and absolute tolerance */
      {
         bi_prod += atolf;
      }
      else /* DEFAULT (stop_crit and atolf exist for backwards compatibilty and are not in the reference manual) */
      {
        /* convergence criteria:  <C*r,r>  <= max( a_tol^2, r_tol^2 * <C*b,b> )
           note: default for a_tol is 0.0, so relative residual criteria is used unless
           user specifies a_tol, or sets r_tol = 0.0, which means absolute
           tol only is checked  */
         eps = jxf_max(r_tol*r_tol, a_tol*a_tol/bi_prod);
      }
   }
   else /* bi_prod == 0.0: the rhs vector b is zero */
   {
      /* Set x equal to zero and return */
      (*(pcg_functions->CopyVector))(vecB, vecX);
      if (logging > 0 || print_level > 0)
      {
         norms[0]     = 0.0;
         rel_norms[i] = 0.0;
      }

      return jxf_error_flag;
      /* In this case, for the original parcsr pcg, the code would take special
         action to force iterations even though the exact value was known. */
   }

   /* r = b - Ax */
   (*(pcg_functions->CopyVector))(vecB, vecR);
   (*(pcg_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
 
   /* p = C*r */
   (*(pcg_functions->ClearVector))(vecP);
   precond(precond_data, preOperator, vecR, vecP);

   /* gamma = <r,p> */
   gamma = (*(pcg_functions->InnerProd))(vecR,vecP);

   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (gamma != 0.) ieee_check = gamma / gamma; /* INF -> NaN conversion */
   if (ieee_check != ieee_check)
   {
      /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
         for ieee_check self-equality works on all IEEE-compliant compilers/
         machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
         by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
         found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
      if (print_level > 0 || logging > 0)
      {
         jxf_printf("\n\nERROR detected by JXFPAMG ...  BEGIN\n");
         jxf_printf("ERROR -- jxf_PCGSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ...  END\n\n\n");
//         break;
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   if (two_norm)
   {
      normdown = (*(pcg_functions->InnerProd))(vecR,vecR);
   }
   else
   {
      normdown = gamma;
   }


   /* Set initial residual norm */
   if ( logging > 0 || print_level > 0 || cf_tol > 0.0 )
   {
      if (two_norm)
      {
         i_prod_0 = (*(pcg_functions->InnerProd))(vecR, vecR);
      }
      else
      {
         i_prod_0 = gamma;
      }

      if ( logging > 0 || print_level > 0 ) 
      {
         norms[0] = sqrt(i_prod_0);
      }
   }
   
   if ( print_level > 0 && my_id == 0 )
   {
      jxf_printf("\n\n");
      if (two_norm)
      {
         if ( stop_crit && !rel_change && atolf == 0 ) 
         {  
            /* pure absolute tolerance */
            jxf_printf("Iters       ||r||_2     conv.rate\n");
            jxf_printf("-----    ------------   ---------\n");
         }
         else 
         {
            jxf_printf("Iters       ||r||_2     conv.rate  ||r||_2/||b||_2\n");
            jxf_printf("-----    ------------   ---------  ------------ \n");
         }
      }
      else  /* !two_norm */
      {
         jxf_printf("Iters       ||r||_C     conv.rate  ||r||_C/||b||_C\n");
         jxf_printf("-----    ------------    ---------  ------------ \n");
      }
   }

   while ((i+1) <= max_iter)
   {
     /*--------------------------------------------------------------------
      * the core CG calculations...
      *--------------------------------------------------------------------*/
      i ++;

      /* s = A*p */
      (*(pcg_functions->Matvec))(matvec_data, 1.0, matA, vecP, 0.0, vecS);

      /* alpha = gamma / <s,p> */
      sdotp = (*(pcg_functions->InnerProd))(vecS, vecP);
      if (sdotp == 0.0)
      {
         if (i == 1) i_prod = i_prod_0;
         break;
      }
      alpha = gamma / sdotp;

      gamma_old = gamma;

      /* x = x + alpha*p */
      (*(pcg_functions->Axpy))(alpha, vecP, vecX);

      /* r = r - alpha*s */
      (*(pcg_functions->Axpy))(-alpha, vecS, vecR);
         
      /* s = C*r */
      (*(pcg_functions->ClearVector))(vecS);
      precond(precond_data, preOperator, vecR, vecS);

      /* gamma = <r,s> */
      gamma = (*(pcg_functions->InnerProd))(vecR, vecS);

      if (two_norm)
      {
         normup = (*(pcg_functions->InnerProd))(vecR,vecR);
      }
      else
      {
         normup = gamma;
      }
      conv_factor = sqrt(normup / normdown);
      normdown = normup;

      /* set i_prod for convergence test */
      if (two_norm)
      {
         i_prod = (*(pcg_functions->InnerProd))(vecR,vecR);
      }
      else
      {
         i_prod = gamma;
      }
 
      /* print norm info */
      if ( logging > 0 || print_level > 0 )
      {
         norms[i]     = sqrt(i_prod);
         rel_norms[i] = bi_prod ? sqrt(i_prod/bi_prod) : 0;
      }
      if ( print_level > 0 && my_id == 0 )
      {
         if (two_norm)
         {
            if ( stop_crit && !rel_change && atolf==0 ) 
            {  
               /* pure absolute tolerance */
               jxf_printf("% 5d    %e    %f\n", i, norms[i], norms[i]/norms[i-1]);
            }
            else 
            {
               jxf_printf("% 5d    %e    %f    %e\n", i, norms[i], norms[i]/norms[i-1], rel_norms[i]);
            }
         }
         else 
         {
            jxf_printf("% 5d    %e    %f    %e\n", i, norms[i], norms[i]/norms[i-1], rel_norms[i]);
         }
      }

      if (i > 1 && conv_factor > convfac_threshold)
      {
         if (my_id == 0 && print_level > 0)
         {
            jxf_printf("\n  Warning: iteration has terminated because the\n");
            jxf_printf("           current conv_factor > convfac_threshold(%lf) !!\n\n",convfac_threshold);
         }
         break;
      }

      if (i > 1 && conv_factor > convfac_threshold)
      {
         if (my_id == 0 && print_level > 0)
         {
            jxf_printf("\n  Warning: iteration has terminated because the\n");
            jxf_printf("           current conv_factor > convfac_threshold(%lf) !!\n\n",convfac_threshold);
         }
         break;
      }

     /*--------------------------------------------------------------------
      * check for convergence
      *--------------------------------------------------------------------*/
      if (i_prod / bi_prod < eps)  /* the basic convergence test */
      {
         tentatively_converged = 1;
      }

     /* At user request, don't trust the convergence test until we've recomputed
        the residual from scratch.  This is expensive in the usual case where an
        the norm is the energy norm.
        This calculation is coded on the assumption that r's accuracy is only a
        concern for problems where CG takes many iterations. */
      if ( tentatively_converged && recompute_residual )
      {
         /* r = b - Ax */
         (*(pcg_functions->CopyVector))(vecB, vecR);
         (*(pcg_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);

         /* set i_prod for convergence test */
         if (two_norm)
         {
            i_prod = (*(pcg_functions->InnerProd))(vecR, vecR);
         }
         else
         {
            /* s = C*r */
            (*(pcg_functions->ClearVector))(vecS);
            precond(precond_data, preOperator, vecR, vecS);
            /* iprod = gamma = <r,s> */
            i_prod = (*(pcg_functions->InnerProd))(vecR, vecS);
         }
         
         if (i_prod / bi_prod >= eps) 
         {
            tentatively_converged = 0;
         }
      }

     /* At user request, don't treat this as converged unless x didn't change
        much in the last iteration. */ 
      if ( tentatively_converged && rel_change && (i_prod > guard_zero_residual) )
      {
         pi_prod = (*(pcg_functions->InnerProd))(vecP, vecP); 
         xi_prod = (*(pcg_functions->InnerProd))(vecX, vecX);
         ratio = alpha*alpha*pi_prod / xi_prod;
         if (ratio >= eps) 
         {
            tentatively_converged = 0;
         }
      }
      
      /* we've passed all the convergence tests, it's for real */
      if ( tentatively_converged )  
      {
         (pcg_data -> converged) = 1;
         break;
      }

      if ( (gamma<1.0e-292) && ((-gamma)<1.0e-292) ) 
      {
         jxf_error(JXF_ERROR_CONV);  
         break;
      }
      
      /* ... gamma should be >=0.  IEEE subnormal numbers are < 2**(-1022)=2.2e-308
         (and >= 2**(-1074)=4.9e-324).  So a gamma this small means we're getting
         dangerously close to subnormal or zero numbers (usually if gamma is small,
         so will be other variables).  Thus further calculations risk a crash.
         Such small gamma generally means no hope of progress anyway. */

     /*--------------------------------------------------------------------
      * Optional test to see if adequate progress is being made.
      * The average convergence factor is recorded and compared
      * against the tolerance 'cf_tol'. The weighting factor is  
      * intended to pay more attention to the test when an accurate
      * estimate for average convergence factor is available.  
      *--------------------------------------------------------------------*/

      if (cf_tol > 0.0)
      {
         cf_ave_0 = cf_ave_1;

        /* i_prod_0 is zero, or (almost) subnormal, yet i_prod wasn't small
           enough to pass the convergence test.  Therefore initial guess was good,
           and we're just calculating garbage - time to bail out before the
           next step, which will be a divide by zero (or close to it). */         
         if ( i_prod_0 < 1.0e-292 ) 
         {
            jxf_error(JXF_ERROR_CONV);  
            break;
         }
	 cf_ave_1 = pow( i_prod / i_prod_0, 1.0/(2.0*i)); 

         weight   = fabs(cf_ave_1 - cf_ave_0);
         weight   = weight / jxf_max(cf_ave_1, cf_ave_0);
         weight   = 1.0 - weight;
         
         if (weight * cf_ave_1 > cf_tol) break;
      }

     /*--------------------------------------------------------------------
      * back to the core CG calculations
      *--------------------------------------------------------------------*/

      /* beta = gamma / gamma_old */
      beta = gamma / gamma_old;

      /* p = s + beta p */
      (*(pcg_functions->ScaleVector))(beta, vecP);   
      (*(pcg_functions->Axpy))(1.0, vecS, vecP);
   }


   /*--------------------------------------------------------------------
    * Finish up with some outputs.
    *--------------------------------------------------------------------*/

   if ( print_level > 1 && my_id == 0 )
   {
      jxf_printf("\n\n");
   }

   (pcg_data -> num_iterations) = i;
   
   if (bi_prod > 0.0)
   {
      (pcg_data -> rel_residual_norm) = sqrt(i_prod / bi_prod);
   }
   else /* actually, we'll never get here... */
   {
      (pcg_data -> rel_residual_norm) = 0.0;
   }

   return jxf_error_flag;
}
