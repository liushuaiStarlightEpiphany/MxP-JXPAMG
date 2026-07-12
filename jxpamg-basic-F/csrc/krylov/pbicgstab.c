//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  pbicgstab.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"
#include "jxf_krylov.h"

/*!
 * \fn JXF_Int JXF_ParCSRBiCGSTABCreate
 */ 
JXF_Int
JXF_BiCGSTABCreate( MPI_Comm comm, JXF_Solver *solver )
{
   jxf_BiCGSTABFunctions *bicgstab_functions =
      jxf_BiCGSTABFunctionsCreate(
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
         jxf_ParKrylovCommInfo,
         jxf_ParKrylovIdentitySetup, 
         jxf_ParKrylovIdentity );

   *solver = ( (JXF_Solver) jxf_BiCGSTABCreate( bicgstab_functions) );
   
   if (!solver) 
   {
      jxf_error_in_arg(2);
   }
    
   return jxf_error_flag;
}


/*-------------------------------------------------------------------------*
 * JXF_BiCGSTABCreate does not exist.  Call the appropriate function which  *
 * also specifies the vector type, e.g. JXF_ParCSRBiCGSTABCreate            *
 *-------------------------------------------------------------------------*/


/*!
 * \fn JXF_Int JXF_BiCGSTABDestroy
 */ 
JXF_Int 
JXF_BiCGSTABDestroy( JXF_Solver solver )
{
   return( jxf_BiCGSTABDestroy( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetup
 */ 
JXF_Int 
JXF_BiCGSTABSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_BiCGSTABSetup( (void *)solver, (void *)matA, (void *)vecB, (void *)vecX ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSolve
 */ 
JXF_Int 
JXF_BiCGSTABSolve( JXF_Solver solver,
                  JXF_Matrix preOperator,
                  JXF_Matrix matA,
                  JXF_Vector vecB,
                  JXF_Vector vecX  )
{
   return( jxf_BiCGSTABSolve( (void *)solver, (void *)preOperator, (void *)matA, (void *)vecB, (void *)vecX ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetPreOperator
 * \brief Set the preconditioner for BiCGSTAB
 * \author peghoty
 * \date 2011/09/04
 */ 
JXF_Int
JXF_BiCGSTABSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator )
{
   return( jxf_BiCGSTABSetPreOperator( (void *) solver, (void *) preOperator ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetTol
 */
JXF_Int
JXF_BiCGSTABSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_BiCGSTABSetTol( (void *) solver, tol ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetAbsoluteTol
 */
JXF_Int
JXF_BiCGSTABSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol )
{
   return( jxf_BiCGSTABSetAbsoluteTol( (void *) solver, a_tol ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetConvergenceFactorTol
 */
JXF_Int
JXF_BiCGSTABSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol )
{
   return( jxf_BiCGSTABSetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetMinIter
 */
JXF_Int
JXF_BiCGSTABSetMinIter( JXF_Solver solver, JXF_Int min_iter )
{
   return( jxf_BiCGSTABSetMinIter( (void *) solver, min_iter ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetMaxIter
 */
JXF_Int
JXF_BiCGSTABSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_BiCGSTABSetMaxIter( (void *) solver, max_iter ) );
}


/*!
 * \fn JXF_Int JXF_BiCGSTABSetConvCriteria
 * \brief Set the parameter 'conv_criteria'.
 * \author peghoty
 * \date 2010/03/20
 */ 
JXF_Int
JXF_BiCGSTABSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria )
{
   return( jxf_BiCGSTABSetConvCriteria( (void *) solver, conv_criteria ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetConvFacThreshold
 * \brief Set the parameter 'convfac_threshold'.
 * \author peghoty
 * \date 2010/06/22
 */ 
JXF_Int
JXF_BiCGSTABSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold )
{
   return( jxf_BiCGSTABSetConvFacThreshold( (void *) solver, convfac_threshold) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetPrecond
 */ 
JXF_Int
JXF_BiCGSTABSetPrecond( JXF_Solver         solver,
                       JXF_PtrToSolverFcn precond,
                       JXF_PtrToSolverFcn precond_setup,
                       JXF_Solver         precond_solver )
{
   return( jxf_BiCGSTABSetPrecond( (void *) solver,
                                  precond, 
                                  precond_setup,
                                  (void *) precond_solver ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABGetPrecond
 */ 
JXF_Int
JXF_BiCGSTABGetPrecond( JXF_Solver solver, JXF_Solver *precond_data_ptr )
{
   return( jxf_BiCGSTABGetPrecond( (void *) solver, (JXF_Solver *) precond_data_ptr ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetLogging
 */ 
JXF_Int
JXF_BiCGSTABSetLogging( JXF_Solver solver, JXF_Int logging )
{
   return( jxf_BiCGSTABSetLogging( (void *) solver, logging ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABSetPrintLevel
 */ 
JXF_Int
JXF_BiCGSTABSetPrintLevel( JXF_Solver solver, JXF_Int print_level )
{
   return( jxf_BiCGSTABSetPrintLevel( (void *) solver, print_level ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABGetNumIterations
 */
JXF_Int
JXF_BiCGSTABGetNumIterations( JXF_Solver  solver, JXF_Int *num_iterations )
{
   return( jxf_BiCGSTABGetNumIterations( (void *) solver, num_iterations ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABGetFinalRelativeResidualNorm
 */
JXF_Int
JXF_BiCGSTABGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm )
{
   return( jxf_BiCGSTABGetFinalRelativeResidualNorm( (void *) solver, norm ) );
}

/*!
 * \fn JXF_Int JXF_BiCGSTABGetResidual
 */
JXF_Int
JXF_BiCGSTABGetResidual( JXF_Solver solver, void **residual  )
{
   return( jxf_BiCGSTABGetResidual( (void *) solver, residual ) );
}

/*!
 * \fn jxf_BiCGSTABFunctions *jxf_BiCGSTABFunctionsCreate
 */
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
   JXF_Int    (*PrecondSetup)  ( void *vdata, void *A ),
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
)
{
   jxf_BiCGSTABFunctions *bicgstab_functions;
   bicgstab_functions = (jxf_BiCGSTABFunctions *)jxf_CTAlloc(jxf_BiCGSTABFunctions, 1);

   bicgstab_functions->CreateVector = CreateVector;
   bicgstab_functions->DestroyVector = DestroyVector;
   bicgstab_functions->MatvecCreate = MatvecCreate;
   bicgstab_functions->Matvec = Matvec;
   bicgstab_functions->MatvecDestroy = MatvecDestroy;
   bicgstab_functions->InnerProd = InnerProd;
   bicgstab_functions->CopyVector = CopyVector;
   bicgstab_functions->ClearVector = ClearVector;
   bicgstab_functions->ScaleVector = ScaleVector;
   bicgstab_functions->Axpy = Axpy;
   bicgstab_functions->CommInfo = CommInfo;
   bicgstab_functions->precond_setup = PrecondSetup;
   bicgstab_functions->precond = Precond;

   return bicgstab_functions;
}

/*!
 * \fn void *jxf_BiCGSTABCreate
 */
void *
jxf_BiCGSTABCreate( jxf_BiCGSTABFunctions * bicgstab_functions )
{
   jxf_BiCGSTABData *bicgstab_data;
 
   bicgstab_data = jxf_CTAlloc( jxf_BiCGSTABData, 1);
 
   bicgstab_data->functions = bicgstab_functions;

   /* set defaults */
   (bicgstab_data -> tol)                = 1.0e-06;
   (bicgstab_data -> min_iter)           = 0;
   (bicgstab_data -> max_iter)           = 1000;
   (bicgstab_data -> conv_criteria)      = 0;    // peghoty 2010/03/20
   (bicgstab_data -> convfac_threshold)  = 1.0e4;  // peghoty 2010/06/22
   (bicgstab_data -> a_tol)              = 0.0;  
   (bicgstab_data -> cf_tol)             = 0.0;  // peghoty 2010/03/20
   (bicgstab_data -> precond_data)       = NULL;
   (bicgstab_data -> logging)            = 0;
   (bicgstab_data -> print_level)        = 0;
   (bicgstab_data -> B)                  = NULL; // peghoty 2011/10/08
   (bicgstab_data -> p)                  = NULL;
   (bicgstab_data -> q)                  = NULL;
   (bicgstab_data -> r)                  = NULL;
   (bicgstab_data -> r0)                 = NULL;
   (bicgstab_data -> s)                  = NULL;    
   (bicgstab_data -> v)                  = NULL;
   (bicgstab_data -> matvec_data)        = NULL;
   (bicgstab_data -> norms)              = NULL;
   (bicgstab_data -> log_file_name)      = NULL;
 
   return (void *) bicgstab_data;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABDestroy
 */
JXF_Int
jxf_BiCGSTABDestroy( void *bicgstab_vdata )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   jxf_BiCGSTABFunctions *bicgstab_functions = bicgstab_data->functions;

   if (bicgstab_data)
   {
      if ( (bicgstab_data -> norms) != NULL )
      {
         jxf_TFree(bicgstab_data -> norms);
      }
 
      (*(bicgstab_functions->MatvecDestroy))(bicgstab_data -> matvec_data);
 
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> r);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> r0);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> s);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> v);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> p);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> q);
      (*(bicgstab_functions->DestroyVector))(bicgstab_data -> Xtmp); // peghoty 2010/03/20
 
      jxf_TFree(bicgstab_data);
      jxf_TFree(bicgstab_functions);
   }
 
   return(jxf_error_flag);
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetup
 */
JXF_Int
jxf_BiCGSTABSetup( void *bicgstab_vdata, void *matA, void *vecB, void *vecX )
{
   jxf_BiCGSTABData      *bicgstab_data      = bicgstab_vdata;
   jxf_BiCGSTABFunctions *bicgstab_functions = bicgstab_data->functions;

   JXF_Int max_iter = (bicgstab_data -> max_iter);   
   
   (bicgstab_data -> A) = matA;

  /*---------------------------------------------------
   * The arguments for NewVector are important to
   * maintain consistency between the setup and
   * compute phases of matvec and the preconditioner.
   *--------------------------------------------------*/
 
   if ((bicgstab_data -> p) == NULL)
   {
      (bicgstab_data -> p) = (*(bicgstab_functions->CreateVector))(vecB);
   }
   if ((bicgstab_data -> q) == NULL)
   {
      (bicgstab_data -> q) = (*(bicgstab_functions->CreateVector))(vecB);
   }
   if ((bicgstab_data -> r) == NULL)
   {
      (bicgstab_data -> r) = (*(bicgstab_functions->CreateVector))(vecB);
   }
   if ((bicgstab_data -> r0) == NULL)
   {
      (bicgstab_data -> r0) = (*(bicgstab_functions->CreateVector))(vecB);
   }
   if ((bicgstab_data -> s) == NULL)
   {
      (bicgstab_data -> s) = (*(bicgstab_functions->CreateVector))(vecB);
   }
   if ((bicgstab_data -> v) == NULL)
   {
      (bicgstab_data -> v) = (*(bicgstab_functions->CreateVector))(vecB);
   }

   /* peghoty 2010/03/20 */
   if ((bicgstab_data -> Xtmp) == NULL)
   {
      (bicgstab_data -> Xtmp) = (*(bicgstab_functions->CreateVector))(vecB); 
   }     
 
   if ((bicgstab_data -> matvec_data) == NULL)
   { 
      (bicgstab_data -> matvec_data) = (*(bicgstab_functions->MatvecCreate))(matA, vecX);
   }

  /*------------------------------------------------------------------------
   * Allocate space for log info
   *-----------------------------------------------------------------------*/
   if ( (bicgstab_data->logging) > 0 || (bicgstab_data->print_level) > 0 )
   {
      if ((bicgstab_data -> norms) == NULL)
      {
         (bicgstab_data -> norms) = jxf_CTAlloc(JXF_Real, max_iter + 1);
      }
   }
   if ((bicgstab_data -> print_level) > 0)
   {
      if ((bicgstab_data -> log_file_name) == NULL)
      {
         (bicgstab_data -> log_file_name) = "bicgstab.out.log";
      }
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetPreOperator
 * \brief Set the preconditioner operator for BiCGSTAB.
 * \author peghoty
 * \date 2011/09/04
 */ 
JXF_Int
jxf_BiCGSTABSetPreOperator( void *bicgstab_vdata, void *preOperator )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> B) = preOperator;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetTol
 */  
JXF_Int
jxf_BiCGSTABSetTol( void *bicgstab_vdata, JXF_Real tol )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> tol) = tol;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetAbsoluteTol
 */ 
JXF_Int
jxf_BiCGSTABSetAbsoluteTol( void *bicgstab_vdata, JXF_Real a_tol )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> a_tol) = a_tol;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetConvergenceFactorTol
 */ 
JXF_Int
jxf_BiCGSTABSetConvergenceFactorTol( void *bicgstab_vdata, JXF_Real cf_tol )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> cf_tol) = cf_tol;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetMinIter
 */ 
JXF_Int
jxf_BiCGSTABSetMinIter( void *bicgstab_vdata, JXF_Int min_iter )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> min_iter) = min_iter;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetMaxIter
 */ 
JXF_Int
jxf_BiCGSTABSetMaxIter( void *bicgstab_vdata, JXF_Int max_iter )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetConvCriteria
 * \brief Set the parameter 'conv_criteria'.
 * \author peghoty
 * \date 2010/03/20
 */  
JXF_Int
jxf_BiCGSTABSetConvCriteria( void *bicgstab_vdata, JXF_Int conv_criteria )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> conv_criteria) = conv_criteria;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetConvFacThreshold
 * \brief Set the parameter 'convfac_threshold'.
 * \author peghoty
 * \date 2010/06/22
 */
JXF_Int
jxf_BiCGSTABSetConvFacThreshold( void *bicgstab_vdata, JXF_Real convfac_threshold )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> convfac_threshold) = convfac_threshold;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetPrecond
 */ 
JXF_Int
jxf_BiCGSTABSetPrecond( void  *bicgstab_vdata,
                       JXF_Int  (*precond)(),
                       JXF_Int  (*precond_setup)(),
                       void  *precond_data  )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   jxf_BiCGSTABFunctions *bicgstab_functions = bicgstab_data->functions;

   (bicgstab_functions -> precond)       = precond;
   (bicgstab_functions -> precond_setup) = precond_setup;
   (bicgstab_data -> precond_data)       = precond_data;
 
   return jxf_error_flag;
}
 
/*!
 * \fn JXF_Int jxf_BiCGSTABGetPrecond
 */ 
JXF_Int
jxf_BiCGSTABGetPrecond( void *bicgstab_vdata, JXF_Solver *precond_data_ptr )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   *precond_data_ptr = (JXF_Solver)(bicgstab_data -> precond_data);
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetLogging
 */ 
JXF_Int
jxf_BiCGSTABSetLogging( void *bicgstab_vdata, JXF_Int logging )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> logging) = logging;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABSetPrintLevel
 */ 
JXF_Int
jxf_BiCGSTABSetPrintLevel( void *bicgstab_vdata, JXF_Int print_level )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   (bicgstab_data -> print_level) = print_level;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BiCGSTABGetConverged
 */ 
JXF_Int
jxf_BiCGSTABGetConverged( void *bicgstab_vdata, JXF_Int *converged )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   *converged = (bicgstab_data -> converged);
   return jxf_error_flag;
}
 
/*!
 * \fn JXF_Int jxf_BiCGSTABGetNumIterations
 */ 
JXF_Int
jxf_BiCGSTABGetNumIterations( void *bicgstab_vdata, JXF_Int *num_iterations )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   *num_iterations = (bicgstab_data -> num_iterations);
   return jxf_error_flag;
}
 
/*!
 * \fn JXF_Int jxf_BiCGSTABGetFinalRelativeResidualNorm
 */ 
JXF_Int
jxf_BiCGSTABGetFinalRelativeResidualNorm( void *bicgstab_vdata, JXF_Real *relative_residual_norm )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   *relative_residual_norm = (bicgstab_data -> rel_residual_norm);
   return jxf_error_flag;
} 

/*!
 * \fn JXF_Int jxf_BiCGSTABGetResidual
 */ 
JXF_Int
jxf_BiCGSTABGetResidual( void *bicgstab_vdata, void **residual )
{
   jxf_BiCGSTABData *bicgstab_data = bicgstab_vdata;
   *residual = (bicgstab_data -> r);
   return jxf_error_flag;
} 

/*!
 * \fn JXF_Int jxf_BiCGSTABSolve
 */ 
JXF_Int
jxf_BiCGSTABSolve( void *bicgstab_vdata,
                  void *preOperator, 
                  void *matA, 
                  void *vecB, 
                  void *vecX  )
{
   jxf_BiCGSTABData      *bicgstab_data      = bicgstab_vdata;
   jxf_BiCGSTABFunctions *bicgstab_functions = bicgstab_data->functions;

   JXF_Int         min_iter      = (bicgstab_data -> min_iter);
   JXF_Int         max_iter      = (bicgstab_data -> max_iter);
   JXF_Int         conv_criteria = (bicgstab_data -> conv_criteria); // peghoty 2010/03/20
   JXF_Real      convfac_threshold = (bicgstab_data -> convfac_threshold); // peghoty 2010/06/22
   JXF_Real      r_tol         = (bicgstab_data -> tol);
   JXF_Real      cf_tol        = (bicgstab_data -> cf_tol);
   void       *matvec_data   = (bicgstab_data -> matvec_data);
   JXF_Real      a_tol         = (bicgstab_data -> a_tol);
  
   void       *vecR          = (bicgstab_data -> r);
   void       *vecR0         = (bicgstab_data -> r0);
   void       *vecS          = (bicgstab_data -> s);
   void       *vecV          = (bicgstab_data -> v);
   void       *vecP          = (bicgstab_data -> p);
   void       *vecQ          = (bicgstab_data -> q);
   void       *Xtmp          = (bicgstab_data -> Xtmp); // peghoty 2010/03/20

   JXF_Int 	      (*precond)()   = (bicgstab_functions -> precond);
   JXF_Int 	      *precond_data  = (bicgstab_data -> precond_data);

   /* logging variables */
   JXF_Int        logging        = (bicgstab_data -> logging);
   JXF_Int        print_level    = (bicgstab_data -> print_level);
   JXF_Real    *norms          = (bicgstab_data -> norms);

//   JXF_Int        ierr = 0;
   JXF_Int        iter; 
   JXF_Int        my_id, num_procs;
   JXF_Real     alpha, beta, gamma, temp, res, r_norm, b_norm;
   JXF_Real     epsilon = 0.0;
   JXF_Real     epsmac = 1.e-128; 
   JXF_Real     ieee_check = 0.;
   JXF_Real     cf_ave_0 = 0.0;
   JXF_Real     cf_ave_1 = 0.0;
   JXF_Real     weight;
   JXF_Real     r_norm_0;
   JXF_Real     den_norm;
   
   /* peghoty 2010/03/20 */
   JXF_Real     old_relative = 0.0;
   JXF_Real     relative     = 0.0;
   JXF_Real     normup       = 0.0;
   JXF_Real     normdown     = 0.0;
   JXF_Real     conv_factor  = 0.0;
     
   (bicgstab_data -> converged) = 0;

   (*(bicgstab_functions->CommInfo))(matA, &my_id, &num_procs);
   
   if (logging > 0 || print_level > 0)
   {
      norms = (bicgstab_data -> norms);
   }


   /* initialize work arrays */
   (*(bicgstab_functions->CopyVector))(vecB, vecR0);


   /* compute initial residual */
   (*(bicgstab_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR0);
   (*(bicgstab_functions->CopyVector))(vecR0, vecR);
   (*(bicgstab_functions->CopyVector))(vecR0, vecP);

   b_norm = sqrt((*(bicgstab_functions->InnerProd))(vecB, vecB));


   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (b_norm != 0.) ieee_check = b_norm / b_norm; /* INF -> NaN conversion */
   if (ieee_check != ieee_check)
   {
      /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
         for ieee_check self-equality works on all IEEE-compliant compilers/
         machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
         by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
         found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
      if (logging > 0 || print_level > 0)
      {
         jxf_printf("\n\nERROR detected by JXFPAMG ...  BEGIN\n");
         jxf_printf("ERROR -- jxf_BiCGSTABSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ...  END\n\n\n");
//         break;
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }


   res = (*(bicgstab_functions->InnerProd))(vecR0, vecR0);
   r_norm = sqrt(res);
   r_norm_0 = r_norm;
 
 
   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (r_norm != 0.) ieee_check = r_norm / r_norm; /* INF -> NaN conversion */
   if (ieee_check != ieee_check)
   {
      /* ...INFs or NaNs in input can make ieee_check a NaN.  This test
         for ieee_check self-equality works on all IEEE-compliant compilers/
         machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
         by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
         found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF */
      if (logging > 0 || print_level > 0)
      {
         jxf_printf("\n\nERROR detected by JXFPAMG ...  BEGIN\n");
         jxf_printf("ERROR -- jxf_BiCGSTABSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ...  END\n\n\n");
//         break;
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   normdown = r_norm;

   if (logging > 0 || print_level > 0)
   {
      norms[0] = r_norm;
      if (print_level > 0 && my_id == 0)
      {
  	 jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
         {
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         }
         jxf_printf("Initial L2 norm of residual: %e\n\n", r_norm);
      }
   }
   
   
   iter = 0;


   /* peghoty 2010/03/20 */
   if (b_norm > 0.0)
   {
      /* convergence criterion |r_i| <= r_tol*|b| if |b| > 0 */
      den_norm = b_norm;
   }
   else
   {
      /* convergence criterion |r_i| <= r_tol*|r0| if |b| = 0 */
      den_norm = r_norm;
   }
   
   epsilon = jxf_max(a_tol, r_tol*den_norm);
   
   if (conv_criteria) relative = 2*r_tol;  
     
   if (print_level > 0 && my_id == 0)
   {
      if (conv_criteria == 0)
      {
         if (b_norm > 0.0)
         {
             jxf_printf("=============================================\n\n");
             jxf_printf("Iters     resid.norm     conv.rate  rel.res.norm(wrt.b_norm)\n");
             jxf_printf("-----    ------------    ---------- ------------------------\n");
         }
         else
         {
             jxf_printf("=============================================\n\n");
             jxf_printf("Iters     resid.norm     conv.rate  rel.res.norm\n");
             jxf_printf("-----    ------------    ---------- ------------\n");
         }
      }
      else if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
      {
         if (conv_criteria == 1)
         {
            jxf_printf("                                      Relative   \n");
            jxf_printf(" Iters     Error Norm-1    Factor    Error Norm-1\n");
            jxf_printf("-------    ------------   --------  -------------\n");
         }
         else if (conv_criteria == 11)
         {
            jxf_printf("                  Relative Point-wise \n");
            jxf_printf(" Iters   Factor      Error Norm-1     \n");
            jxf_printf("-------  -------    --------------    \n");
         }
         else if (conv_criteria == 2)
         {
            jxf_printf("                                      Relative   \n");
            jxf_printf(" Iters     Error Norm-2    Factor    Error Norm-2\n");
            jxf_printf("-------    ------------   --------  -------------\n");          
         }
      }      
   }

   (bicgstab_data -> num_iterations) = iter;
   
   if (b_norm > 0.0)
   {
      (bicgstab_data -> rel_residual_norm) = r_norm / b_norm;
   }   
      
   
  /*-----------------------------------------------------------------
   * The main 'while' loop 
   *---------------------------------------------------------------*/  
    
   while (iter < max_iter)
   {

        if (r_norm == 0.0)
        {
//	   ierr = 0;
	   return jxf_error_flag;
	}
	
	if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
        {
           (*(bicgstab_functions->CopyVector))(vecX, Xtmp);
        }
	
	/* check for convergence, evaluate actual residual */
	if (conv_criteria == 0)
	{ 
	   if (r_norm <= epsilon && iter >= min_iter) 
           {
	      (*(bicgstab_functions->CopyVector))(vecB,vecR);
              (*(bicgstab_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
	      r_norm = sqrt((*(bicgstab_functions->InnerProd))(vecR, vecR));
	      if (r_norm <= epsilon)
              {
                 if (print_level > 0 && my_id == 0)
                 {
                    jxf_printf("\n\n");
                    jxf_printf(" Final L2 norm of residual: %e\n\n", r_norm);
                 }
                 (bicgstab_data -> converged) = 1;
                 break;
              }
	      else
	      {
	         (*(bicgstab_functions->CopyVector))(vecR, vecP);
	      }
	   }
	}
	else if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
        {           
	   if (relative <= r_tol && iter >= min_iter) 
           {
              /* Is this necessary? peghoty 2010/03/20 */
	      (*(bicgstab_functions->CopyVector))(vecB, vecR);
              (*(bicgstab_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
	      r_norm = sqrt((*(bicgstab_functions->InnerProd))(vecR, vecR));
	      if (r_norm <= epsilon) 
              {
                 if (print_level > 0 && my_id == 0)
                 {
                    jxf_printf("\n\n");
                    jxf_printf(" Final L2 norm of residual: %e\n\n", r_norm);
                 }
                 (bicgstab_data -> converged) = 1;
                 break;
              }
	      else
	      {
	         (*(bicgstab_functions->CopyVector))(vecR, vecP);
	      }
	   }
        }

	
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
           cf_ave_1 = pow( r_norm / r_norm_0, 1.0/(2.0*iter));

           weight   = fabs(cf_ave_1 - cf_ave_0);
           weight   = weight / jxf_max(cf_ave_1, cf_ave_0);
           weight   = 1.0 - weight;
           if (weight * cf_ave_1 > cf_tol) 
           {
              break;
           }
        }


        iter ++;
	
	(*(bicgstab_functions->ClearVector))(vecV);
		
        precond(precond_data, preOperator, vecP, vecV);
	        
        (*(bicgstab_functions->Matvec))(matvec_data, 1.0, matA, vecV, 0.0, vecQ);
	        
      	temp = (*(bicgstab_functions->InnerProd))(vecR0, vecQ);
	      	
      	if (fabs(temp) >= epsmac)
      	{
	   alpha = res / temp;
	}
	else
	{
	   jxf_printf("BiCGSTAB broke down!! divide by near zero\n");
	   return (1);
	}
		
	(*(bicgstab_functions->Axpy))(alpha, vecV, vecX);
	(*(bicgstab_functions->Axpy))(-alpha, vecQ, vecR);
	
	(*(bicgstab_functions->ClearVector))(vecV);
	
        precond(precond_data, preOperator, vecR, vecV);
       
        (*(bicgstab_functions->Matvec))(matvec_data, 1.0, matA, vecV, 0.0, vecS);
        
      	gamma = (*(bicgstab_functions->InnerProd))(vecR, vecS) /
                (*(bicgstab_functions->InnerProd))(vecS, vecS);
          
	(*(bicgstab_functions->Axpy))(gamma, vecV, vecX);

	
	if (conv_criteria == 1)
	{
           old_relative = relative;
           (*(bicgstab_functions->Axpy))(-1.0, vecX, Xtmp);  // Xtmp = Xtmp - vecX
           normup   = jxf_ParVectorNorm1(Xtmp);
           normdown = jxf_ParVectorNorm1(vecX);
           if (normdown < EPS)
           {
              relative = normup;
           }
           else
           {
              relative = normup / normdown;
           }  
           conv_factor = relative / old_relative;
	}
	else if (conv_criteria == 11)
	{
            old_relative = relative;
            (*(bicgstab_functions->Axpy))(-1.0, vecX, Xtmp);  // Xtmp = Xtmp - vecX
            relative = jxf_ParVectorPointWiseRelNorm1(Xtmp, vecX);
            conv_factor = relative / old_relative;
	}
	else if (conv_criteria == 2)
	{
            old_relative = relative;
            (*(bicgstab_functions->Axpy))(-1.0, vecX, Xtmp);  // Xtmp = Xtmp - vecX
            normup = jxf_ParVectorNorm2(Xtmp);
            normdown = jxf_ParVectorNorm2(vecX);
            if (normdown < EPS)
            {
               relative = normup;
            }
            else
            {
               relative = normup / normdown;
            }
            conv_factor = relative / old_relative;
	}
	

	(*(bicgstab_functions->Axpy))(-gamma, vecS, vecR);
      	if (fabs(res) >= epsmac)
        {
           beta = 1.0 / res;
        }
	else
	{
	   jxf_printf("BiCGSTAB broke down!! res = 0 \n");
	   return(2);
	}
        res = (*(bicgstab_functions->InnerProd))(vecR0, vecR);
        beta *= res;    
	(*(bicgstab_functions->Axpy))(-gamma, vecQ, vecP);
      	if (fabs(gamma) >= epsmac)
      	{
           (*(bicgstab_functions->ScaleVector))((beta*alpha / gamma), vecP);
        }
	else
	{
	   jxf_printf("BiCGSTAB broke down!! gamma = 0 \n");
	   return (3);
	}
	(*(bicgstab_functions->Axpy))(1.0, vecR, vecP);

	r_norm = sqrt((*(bicgstab_functions->InnerProd))(vecR, vecR));

        normup = r_norm;
        conv_factor = normup / normdown;
        normdown = normup;
	
	if (logging > 0 || print_level > 0)
	{
	   norms[iter] = r_norm;
	}


        if (print_level > 0 && my_id == 0)
	{
           if (conv_criteria == 0)
           {
              if (b_norm > 0.0)
              {
                 jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
			 norms[iter]/norms[iter-1], norms[iter]/b_norm);
              }
              else
              {
                 jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
		         norms[iter]/norms[iter-1], norms[iter]/r_norm_0);
	      }
	   }
           if (conv_criteria == 1 || conv_criteria == 2)
           {
              if (iter == 1)
              {
                 jxf_printf("   %3d     %e             %e \n", iter, normup, relative);
              }
              else
              {
                 jxf_printf("   %3d     %e   %6f  %e \n", iter, normup, conv_factor, relative);
              }
           }
           if (conv_criteria == 11)
           {
              if (iter == 1)
              {
                 jxf_printf("   %3d               %e \n", iter, relative);
              }
              else
              {
                 jxf_printf("   %3d   %6f    %e \n", iter, conv_factor, relative);
              }
           }
	}

        if (iter > 1 && conv_factor > convfac_threshold)
        {
           if (my_id == 0 && print_level > 0)
           {
              jxf_printf("\n  Warning: iteration has terminated because the\n");
              jxf_printf("           current conv_factor > convfac_threshold(%lf) !!\n\n",convfac_threshold);
           }
           break;
        }

   }

   (bicgstab_data -> num_iterations) = iter;
   
   if (conv_criteria == 0)
   {
      if (b_norm > 0.0)
      {
         (bicgstab_data -> rel_residual_norm) = r_norm / b_norm;
      }
      if (b_norm == 0.0)
      {
         (bicgstab_data -> rel_residual_norm) = r_norm / r_norm_0;
      }
   }
   else
   {
      (bicgstab_data -> rel_residual_norm) = relative;
   }
   
#if 0
   if (iter >= max_iter && r_norm > epsilon) 
   {
      jxf_printf(" Warning: Not reaching the given tolerance in %d iterations!!\n\n", max_iter);
      jxf_error(JXF_ERROR_CONV);
   }
#endif

   return jxf_error_flag;
}
