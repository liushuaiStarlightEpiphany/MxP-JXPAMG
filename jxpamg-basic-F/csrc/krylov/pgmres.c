//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  pgmres.c
 *  Date: 2011/09/03
 */ 

#include "jxf_combined.h"
#include "jxf_krylov.h"

JXF_Int jxf__kdim_memory_error = -1; // Feng Chunsheng & Yue Xiaoqiang 2012/10/26

/*!
 * \fn JXF_Int GMRESCreate
 */ 
JXF_Int
JXF_GMRESCreate( MPI_Comm comm, JXF_Solver *solver )
{
   jxf_GMRESFunctions * gmres_functions =
      jxf_GMRESFunctionsCreate(
         jxf_CAlloc, 
         jxf_ParKrylovFree, 
         jxf_ParKrylovCommInfo,
         jxf_ParKrylovCreateVector,
         jxf_ParKrylovCreateVectorArray,
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

   *solver = ( (JXF_Solver) jxf_GMRESCreate( gmres_functions ) );
   
   if (!solver) 
   {
      jxf_error_in_arg(2);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int GMRESDestroy
 */ 
JXF_Int 
JXF_GMRESDestroy( JXF_Solver solver )
{
   return( jxf_GMRESDestroy( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetup
 */ 
JXF_Int 
JXF_GMRESSetup( JXF_Solver solver, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_GMRESSetup( solver, matA, vecB, vecX ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSolve
 */
JXF_Int 
JXF_GMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_GMRESSolve( solver, preOperator, matA, vecB, vecX ) );
}

/*!
 * \fn JXF_Int JXF_GMRESAdaptiveSolve
 */
JXF_Int
JXF_GMRESAdaptiveSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix matA, JXF_Vector vecB, JXF_Vector vecX )
{
   return( jxf_GMRESAdaptiveSolve( solver, preOperator, matA, vecB, vecX ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetPreOperator
 * \brief Set preconditioner operater.
 * \author peghoty
 * \date 2011/09/04 
 */
JXF_Int
JXF_GMRESSetPreOperator( JXF_Solver solver, JXF_ParCSRMatrix preOperator )
{
   return( jxf_GMRESSetPreOperator( (void *) solver, (void *) preOperator ) );
}

JXF_Int
JXF_GMRESSetKDim( JXF_Solver solver, JXF_Int k_dim )
{
   return( jxf_GMRESSetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_GMRESGetKDim( JXF_Solver solver, JXF_Int *k_dim )
{
   return( jxf_GMRESGetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_GMRESSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_GMRESSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_GMRESGetTol( JXF_Solver solver, JXF_Real *tol )
{
   return( jxf_GMRESGetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_GMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol )
{
   return( jxf_GMRESSetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_GMRESGetAbsoluteTol( JXF_Solver solver, JXF_Real *a_tol )
{
   return( jxf_GMRESGetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_GMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol )
{
   return( jxf_GMRESSetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

JXF_Int
JXF_GMRESGetConvergenceFactorTol( JXF_Solver solver, JXF_Real *cf_tol )
{
   return( jxf_GMRESGetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

JXF_Int
JXF_GMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter )
{
   return( jxf_GMRESSetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_GMRESGetMinIter( JXF_Solver solver, JXF_Int *min_iter )
{
   return( jxf_GMRESGetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_GMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_GMRESSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_GMRESGetMaxIter( JXF_Solver solver, JXF_Int *max_iter )
{
   return( jxf_GMRESGetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_GMRESSetStopCrit( JXF_Solver solver, JXF_Int stop_crit )
{
   return( jxf_GMRESSetStopCrit( (void *) solver, stop_crit ) );
}

JXF_Int
JXF_GMRESGetStopCrit( JXF_Solver solver, JXF_Int *stop_crit )
{
   return( jxf_GMRESGetStopCrit( (void *) solver, stop_crit ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetConvCriteria
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
JXF_GMRESSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria )
{
   return( jxf_GMRESSetConvCriteria( (void *) solver, conv_criteria ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetConvFacThreshold
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
JXF_GMRESSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold )
{
   return( jxf_GMRESSetConvFacThreshold( (void *) solver, convfac_threshold ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetResDownZeroThreshold
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_GMRESSetResDownZeroThreshold( JXF_Solver solver, JXF_Real resdown_0_threshold )
{
   return( jxf_GMRESSetResDownZeroThreshold( (void *) solver, resdown_0_threshold ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetConvFacThresholdTwo
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_GMRESSetConvFacThresholdTwo( JXF_Solver solver, JXF_Real convfac_threshold_2 )
{
   return( jxf_GMRESSetConvFacThresholdTwo( (void *) solver, convfac_threshold_2 ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetIsCheckRestarted
 * \author peghoty
 * \date 2011/11/08 
 */
JXF_Int
JXF_GMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted )
{
   return( jxf_GMRESSetIsCheckRestarted( (void *) solver, is_check_restarted ) );
}

JXF_Int
JXF_GMRESSetRelChange( JXF_Solver solver, JXF_Int rel_change )
{
   return( jxf_GMRESSetRelChange( (void *) solver, rel_change ) );
}

JXF_Int
JXF_GMRESGetRelChange( JXF_Solver solver, JXF_Int *rel_change )
{
   return( jxf_GMRESGetRelChange( (void *) solver, rel_change ) );
}

/*!
 * \fn JXF_Int JXF_GMRESGetLastPrecondType
 * \brief Get last_precond_type
 * \author Yue Xiaoqiang
 * \date 2014/10/24
 */
JXF_Int
JXF_GMRESGetLastPrecondType( JXF_Solver solver, JXF_Int *last_precond_type )
{
   return( jxf_GMRESGetLastPrecondType( (void *) solver, last_precond_type ) );
}

/*!
 * \fn JXF_Int JXF_GMRESSetPrecond
 */
JXF_Int
JXF_GMRESSetPrecond( JXF_Solver          solver,
                    JXF_PtrToSolverFcn  precond,
                    JXF_PtrToSolverFcn  precond_setup,
                    JXF_Solver          precond_solver )
{
   return( jxf_GMRESSetPrecond( (void *) solver,
                               precond, 
                               precond_setup,
                               (void *) precond_solver ) );
}

/*!
 * \fn JXF_Int JXF_GMRESGetPrecond
 */
JXF_Int
JXF_GMRESGetPrecond( JXF_Solver solver, JXF_Solver *precond_data_ptr )
{
   return( jxf_GMRESGetPrecond( (void *) solver, (JXF_Solver *) precond_data_ptr ) );
}

JXF_Int
JXF_GMRESSetPrintLevel( JXF_Solver solver, JXF_Int level )
{
   return( jxf_GMRESSetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_GMRESGetPrintLevel( JXF_Solver solver, JXF_Int *level )
{
   return( jxf_GMRESGetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_GMRESGetCompAMGILUBcoIter( JXF_Solver solver, JXF_Int *amg_iter, JXF_Int *ilu_iter, JXF_Int *bco_iter )
{
   return( jxf_GMRESGetCompAMGILUBcoIter( (void *) solver, amg_iter, ilu_iter, bco_iter ) );
}

JXF_Int
JXF_GMRESSetLogging( JXF_Solver solver, JXF_Int level )
{
   return( jxf_GMRESSetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_GMRESGetLogging( JXF_Solver solver, JXF_Int * level )
{
   return( jxf_GMRESGetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_GMRESGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations )
{
   return( jxf_GMRESGetNumIterations( (void *) solver, num_iterations ) );
}

JXF_Int
JXF_GMRESGetConverged( JXF_Solver solver, JXF_Int *converged )
{
   return( jxf_GMRESGetConverged( (void *) solver, converged ) );
}

JXF_Int
JXF_GMRESGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm )
{
   return( jxf_GMRESGetFinalRelativeResidualNorm( (void *) solver, norm ) );
}

JXF_Int JXF_GMRESGetResidual( JXF_Solver solver, void **residual )
{
   /* returns a pointer to the residual vector */
   return jxf_GMRESGetResidual( (void *) solver, residual );
}

/*!
 * \fn void *jxf_ParCSRGMRESCreate
 * \author peghoty
 * \date 2011/09/03
 */
void *
jxf_ParCSRGMRESCreate( MPI_Comm comm )
{
   jxf_GMRESData *solver = NULL;
   
   jxf_GMRESFunctions * gmres_functions =
      jxf_GMRESFunctionsCreate(
         jxf_CAlloc, 
         jxf_ParKrylovFree, 
         jxf_ParKrylovCommInfo,
         jxf_ParKrylovCreateVector,
         jxf_ParKrylovCreateVectorArray,
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

   solver = ( (jxf_GMRESData *)jxf_GMRESCreate( gmres_functions ) );
   
   if (!solver) 
   {
      jxf_error_in_arg(2);
   }

   return solver;
}

/*!
 * \fn jxf_GMRESFunctions *jxf_GMRESFunctionsCreate
 */
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
)
{
   jxf_GMRESFunctions *gmres_functions;
   gmres_functions = (jxf_GMRESFunctions *)CAlloc( 1, sizeof(jxf_GMRESFunctions) );

   gmres_functions->CAlloc = CAlloc;
   gmres_functions->Free = Free;
   gmres_functions->CommInfo = CommInfo; /* not in PCGFunctionsCreate */
   gmres_functions->CreateVector = CreateVector;
   gmres_functions->CreateVectorArray = CreateVectorArray; /* not in PCGFunctionsCreate */
   gmres_functions->DestroyVector = DestroyVector;
   gmres_functions->MatvecCreate = MatvecCreate;
   gmres_functions->Matvec = Matvec;
   gmres_functions->MatvecDestroy = MatvecDestroy;
   gmres_functions->InnerProd = InnerProd;
   gmres_functions->CopyVector = CopyVector;
   gmres_functions->ClearVector = ClearVector;
   gmres_functions->ScaleVector = ScaleVector;
   gmres_functions->Axpy = Axpy;
   /* default preconditioner must be set here but can be changed later... */
   gmres_functions->precond_setup = PrecondSetup;
   gmres_functions->precond       = Precond;

   return gmres_functions;
}

/*!
 * \fn void *jxf_GMRESCreate
 */
void *
jxf_GMRESCreate( jxf_GMRESFunctions *gmres_functions )
{
   jxf_GMRESData *gmres_data;
 
   gmres_data = jxf_CTAllocF(jxf_GMRESData, 1, gmres_functions);

   gmres_data->functions = gmres_functions;
 
   /* set defaults */
   (gmres_data -> k_dim)               = 5;
   (gmres_data -> tol)                 = 1.0e-06; /* relative residual tol */
   (gmres_data -> cf_tol)              = 0.0;
   (gmres_data -> a_tol)               = 0.0;     /* abs. residual tol */
   (gmres_data -> min_iter)            = 0;
   (gmres_data -> max_iter)            = 1000;
   (gmres_data -> rel_change)          = 0;
   (gmres_data -> stop_crit)           = 0;       /* rel. residual norm  - this is obsolete! */
   (gmres_data -> conv_criteria)       = 0;       /* peghoty 2010/06/23 */
   (gmres_data -> convfac_threshold)   = 1.0e4;   /* peghoty 2010/06/23 */
   (gmres_data -> resdown_0_threshold) = 0.1;     /* Yue Xiaoqiang, 2014/03/24 */
   (gmres_data -> convfac_threshold_2) = 0.1;     /* Yue Xiaoqiang, 2014/03/24 */
   (gmres_data -> comp_amg_iter)       = 0;       /* Yue Xiaoqiang, 2016/02/20 */
   (gmres_data -> comp_ilu_iter)       = 0;       /* Yue Xiaoqiang, 2016/02/20 */
   (gmres_data -> comp_bco_iter)       = 0;       /* Yue Xiaoqiang, 2016/02/20 */
   (gmres_data -> is_check_restarted)  = 0;       /* peghoty 2011/11/08 */
   (gmres_data -> converged)           = 0;
   (gmres_data -> precond_data)        = NULL;
   (gmres_data -> print_level)         = 0;
   (gmres_data -> logging)             = 0;
   (gmres_data -> B)                   = NULL;    /* peghoty 2011/10/08 */
   (gmres_data -> p)                   = NULL;
   (gmres_data -> r)                   = NULL;
   (gmres_data -> w)                   = NULL;
   (gmres_data -> w_2)                 = NULL;
   (gmres_data -> matvec_data)         = NULL;
   (gmres_data -> norms)               = NULL;
   (gmres_data -> log_file_name)       = NULL;
 
   return (void *) gmres_data;
}

/*!
 * \fn JXF_Int jxf_GMRESDestroy
 */
JXF_Int
jxf_GMRESDestroy( void *gmres_vdata )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = gmres_data->functions;
   JXF_Int i;
 
   if (gmres_data)
   {
      if ( (gmres_data->logging>0) || (gmres_data->print_level) > 0 )
      {
         if ( (gmres_data -> norms) != NULL )
         {
            jxf_TFreeF( gmres_data -> norms, gmres_functions );
         }
      }
 
      if ( (gmres_data -> matvec_data) != NULL )
      {
         (*(gmres_functions->MatvecDestroy))(gmres_data -> matvec_data);
      }
 
      if ( (gmres_data -> r) != NULL )
      {
         (*(gmres_functions->DestroyVector))(gmres_data -> r);
      }
      
      if ( (gmres_data -> w) != NULL )
      {
         (*(gmres_functions->DestroyVector))(gmres_data -> w);
      }
      
      if ( (gmres_data -> w_2) != NULL )
      {
         (*(gmres_functions->DestroyVector))(gmres_data -> w_2);
      }

      if ( (gmres_data -> p) != NULL )
      {
         for (i = 0; i < (gmres_data -> k_dim+1); i ++)
         {
            if ( (gmres_data -> p)[i] != NULL )
            {
	       (*(gmres_functions->DestroyVector))( (gmres_data -> p)[i] );
	    }
         }
         jxf_TFreeF( gmres_data->p, gmres_functions );
      }
      jxf_TFreeF( gmres_data, gmres_functions );
      jxf_TFreeF( gmres_functions, gmres_functions );
   }
 
   return jxf_error_flag;
}

JXF_Int 
jxf_GMRESGetResidual( void *gmres_vdata, void **residual )
{
   /* returns a pointer to the residual vector */
   jxf_GMRESData  *gmres_data = gmres_vdata;
   *residual = gmres_data->r;
   return jxf_error_flag;
   
}

/*!
 * \fn JXF_Int jxf_GMRESSetup
 */
JXF_Int
jxf_GMRESSetup( void *gmres_vdata, void *matA, void *vecB, void *vecX )
{
   jxf_GMRESData      *gmres_data      = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = gmres_data->functions;

   JXF_Int            k_dim            = (gmres_data -> k_dim);
   JXF_Int            max_iter         = (gmres_data -> max_iter);
   JXF_Int            rel_change       = (gmres_data -> rel_change);
    
   (gmres_data -> A) = matA;
 
  /*--------------------------------------------------
   * The arguments for NewVector are important to
   * maintain consistency between the setup and
   * compute phases of matvec and the preconditioner.
   *--------------------------------------------------*/
 
 
   if ((gmres_data -> r) == NULL)
   {
      (gmres_data -> r) = (*(gmres_functions->CreateVector))(vecB);
   }
   
   if ((gmres_data -> w) == NULL)
   {
      (gmres_data -> w) = (*(gmres_functions->CreateVector))(vecB);
   }
 
   if (rel_change)
   {  
      if ((gmres_data -> w_2) == NULL)
      {
         (gmres_data -> w_2) = (*(gmres_functions->CreateVector))(vecB);
      }
   }
   
   if ((gmres_data -> matvec_data) == NULL)
   {
      (gmres_data -> matvec_data) = (*(gmres_functions->MatvecCreate))(matA, vecX);
   }

  /*-----------------------------------------------------
   * Allocate space for log info
   *-----------------------------------------------------*/
 
   if ( (gmres_data->logging) > 0 || (gmres_data->print_level) > 0 )
   {
      if ((gmres_data -> norms) == NULL)
      {
         (gmres_data -> norms) = jxf_CTAllocF(JXF_Real, max_iter + 1, gmres_functions);
      }
   }
   if ( (gmres_data->print_level) > 0 ) 
   {
      if ((gmres_data -> log_file_name) == NULL)
      {
         (gmres_data -> log_file_name) = "gmres.out.log";
      }
   }
 
   if ((gmres_data -> p) == NULL)
   {
      (gmres_data -> p) = (*(gmres_functions->CreateVectorArray))(k_dim + 1, vecX);
      if (jxf__kdim_memory_error > 1) // Feng Chunsheng & Yue Xiaoqiang 2012/10/26
      {
         if ( (gmres_data->print_level) > 0 )
         {
            jxf_printf("\n\nWARNING: K_DIM = %d IS TOO BIG, ", k_dim);
            jxf_printf("MODIFIED TO BE THE SUGGESTED VALUE %d\n\n", jxf__kdim_memory_error-1);
         }
         k_dim = jxf__kdim_memory_error - 1;
         gmres_data -> k_dim = k_dim;
         jxf__kdim_memory_error = -1;
         (gmres_data -> p) = (*(gmres_functions->CreateVectorArray))(k_dim + 1, vecX);
      }
   }
 
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetPreOperator
 * \brief Set preconditioner operater.
 * \author peghoty
 * \date 2011/09/04 
 */
JXF_Int
jxf_GMRESSetPreOperator( void *gmres_vdata, void *preOperator )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> B) = preOperator;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetKDim( void *gmres_vdata, JXF_Int k_dim )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> k_dim) = k_dim;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetKDim( void *gmres_vdata, JXF_Int *k_dim )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *k_dim = (gmres_data -> k_dim);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetTol( void *gmres_vdata, JXF_Real tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> tol) = tol;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetTol( void *gmres_vdata, JXF_Real *tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *tol = (gmres_data -> tol);
   return jxf_error_flag;
}
 
JXF_Int
jxf_GMRESSetAbsoluteTol( void *gmres_vdata, JXF_Real a_tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> a_tol) = a_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetAbsoluteTol( void *gmres_vdata, JXF_Real *a_tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *a_tol = (gmres_data -> a_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetConvergenceFactorTol( void *gmres_vdata, JXF_Real cf_tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> cf_tol) = cf_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetConvergenceFactorTol( void *gmres_vdata, JXF_Real *cf_tol )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *cf_tol = (gmres_data -> cf_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetMinIter( void *gmres_vdata, JXF_Int min_iter )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> min_iter) = min_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetMinIter( void *gmres_vdata, JXF_Int *min_iter )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *min_iter = (gmres_data -> min_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetMaxIter( void *gmres_vdata, JXF_Int max_iter )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetMaxIter( void *gmres_vdata, JXF_Int *max_iter )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *max_iter = (gmres_data -> max_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetRelChange( void *gmres_vdata, JXF_Int rel_change )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> rel_change) = rel_change;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetRelChange( void *gmres_vdata, JXF_Int *rel_change )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *rel_change = (gmres_data -> rel_change);
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESGetLastPrecondType
 * \brief Get last_precond_type
 * \author Yue Xiaoqiang
 * \date 2014/10/24
 */
JXF_Int
jxf_GMRESGetLastPrecondType( void *gmres_vdata, JXF_Int *last_precond_type )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   JXF_Int break_adaptive = (gmres_data -> break_adaptive);
   JXF_Int *precond_data = (gmres_data -> precond_data);
   *last_precond_type = break_adaptive * 10 + jxf_CombinedPrecondDataPreID((jxf_CombinedPrecondData *)precond_data);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetStopCrit( void *gmres_vdata, JXF_Int stop_crit )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> stop_crit) = stop_crit;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetStopCrit( void *gmres_vdata, JXF_Int *stop_crit )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *stop_crit = (gmres_data -> stop_crit);
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetConvCriteria
 * \brief Set preconditioner operater.
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
jxf_GMRESSetConvCriteria( void *gmres_vdata, JXF_Int conv_criteria )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> conv_criteria ) = conv_criteria;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetConvFacThreshold
 * \brief Set preconditioner operater.
 * \author peghoty
 * \date 2010/06/23 
 */
JXF_Int
jxf_GMRESSetConvFacThreshold( void *gmres_vdata, JXF_Real convfac_threshold )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> convfac_threshold ) = convfac_threshold;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetResDownZeroThreshold
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
jxf_GMRESSetResDownZeroThreshold( void *gmres_vdata, JXF_Real resdown_0_threshold )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> resdown_0_threshold ) = resdown_0_threshold;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetConvFacThresholdTwo
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
jxf_GMRESSetConvFacThresholdTwo( void *gmres_vdata, JXF_Real convfac_threshold_2 )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> convfac_threshold_2 ) = convfac_threshold_2;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetIsCheckRestarted
 * \brief Set is_check_restarted.
 * \author peghoty
 * \date 2011/11/08 
 */
JXF_Int
jxf_GMRESSetIsCheckRestarted( void *gmres_vdata, JXF_Int is_check_restarted )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> is_check_restarted ) = is_check_restarted;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSetPrecond
 */
JXF_Int
jxf_GMRESSetPrecond( void  *gmres_vdata,
                    JXF_Int  (*precond)(),
                    JXF_Int  (*precond_setup)(),
                    void  *precond_data )
{
   jxf_GMRESData      *gmres_data      = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = gmres_data->functions;
   
   (gmres_functions -> precond)       = precond;
   (gmres_functions -> precond_setup) = precond_setup;
   (gmres_data -> precond_data)       = precond_data;
 
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESGetPrecond
 */
JXF_Int
jxf_GMRESGetPrecond( void *gmres_vdata, JXF_Solver *precond_data_ptr )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *precond_data_ptr = (JXF_Solver)(gmres_data -> precond_data);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetPrintLevel( void *gmres_vdata, JXF_Int level )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> print_level) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetPrintLevel( void *gmres_vdata, JXF_Int *level )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *level = (gmres_data -> print_level);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetCompAMGILUBcoIter( void *gmres_vdata, JXF_Int *amg_iter, JXF_Int *ilu_iter, JXF_Int *bco_iter )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *amg_iter = (gmres_data -> comp_amg_iter);
   *ilu_iter = (gmres_data -> comp_ilu_iter);
   *bco_iter = (gmres_data -> comp_bco_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESSetLogging( void *gmres_vdata, JXF_Int level )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   (gmres_data -> logging) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetLogging( void *gmres_vdata, JXF_Int *level )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *level = (gmres_data -> logging);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetNumIterations( void *gmres_vdata, JXF_Int *num_iterations )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *num_iterations = (gmres_data -> num_iterations);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetConverged( void *gmres_vdata, JXF_Int *converged )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *converged = (gmres_data -> converged);
   return jxf_error_flag;
}

JXF_Int
jxf_GMRESGetFinalRelativeResidualNorm( void *gmres_vdata, JXF_Real *relative_residual_norm )
{
   jxf_GMRESData *gmres_data = gmres_vdata;
   *relative_residual_norm = (gmres_data -> rel_residual_norm);
   return jxf_error_flag;
} 

/*!
 * \fn JXF_Int jxf_GMRESSolve
 */
JXF_Int
jxf_GMRESSolve( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX )
{
   jxf_GMRESData      *gmres_data      = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = (gmres_data -> functions);
   JXF_Int 		     k_dim            = (gmres_data -> k_dim);
   JXF_Int               min_iter         = (gmres_data -> min_iter);
   JXF_Int 		     max_iter         = (gmres_data -> max_iter);
   JXF_Int               rel_change       = (gmres_data -> rel_change);
   JXF_Real 	     r_tol            = (gmres_data -> tol);
   JXF_Real 	     cf_tol           = (gmres_data -> cf_tol);
   JXF_Real            a_tol            = (gmres_data -> a_tol);
   void             *matvec_data      = (gmres_data -> matvec_data);

   JXF_Real            convfac_threshold = (gmres_data -> convfac_threshold); // peghoty 2010/06/23

   void             *vecR             = (gmres_data -> r);
   void             *vecW             = (gmres_data -> w);
   void             *vecW2            = (gmres_data -> w_2); /* vecW2 is only allocated if rel_change = 1 */
   void            **vecP             = (gmres_data -> p);
   JXF_Int 	           (*precond)()       = (gmres_functions -> precond);
   JXF_Int 	            *precond_data     = (gmres_data -> precond_data);
   JXF_Int             print_level        = (gmres_data -> print_level);
   JXF_Int             logging            = (gmres_data -> logging);
   JXF_Real         *norms              = (gmres_data -> norms);
   JXF_Int             is_check_restarted = (gmres_data -> is_check_restarted); /* Xu Xiaowen, 2011/11/08 */
   
   JXF_Int        break_value = 0;
   JXF_Int	      i, j, k;
   JXF_Real     *rs, **hh, *c, *s, *rs_2 = NULL; 
   JXF_Int        iter; 
   JXF_Int        my_id, num_procs;
   JXF_Real     epsilon, gamma, t, r_norm, b_norm, den_norm, x_norm;
   JXF_Real     w_norm;
   
   JXF_Real     epsmac = 1.e-16; 
   JXF_Real     ieee_check = 0.0;

   JXF_Real     guard_zero_residual; 
   JXF_Real     cf_ave_0 = 0.0;
   JXF_Real     cf_ave_1 = 0.0;
   JXF_Real     weight;
   JXF_Real     r_norm_0;
   JXF_Real     relative_error = 1.0;

   JXF_Int        rel_change_passed = 0, num_rel_change_check = 0;

   /* peghoty 2010/06/23 */
   JXF_Real  normup       = 0.0;
   JXF_Real  normdown     = 0.0;
   JXF_Real  conv_factor  = 0.0;

#if 0
   /* 保存离散系统和初始解向量 */
   jxf_CSRMatrix *ser_matrix = NULL;
   ser_matrix = jxf_ParCSRMatrixToCSRMatrixAll((jxf_ParCSRMatrix *)matA);
   jxf_CSRMatrixPrint(ser_matrix,"./A_ser");
   
   jxf_Vector *ser_rhs = NULL;
   ser_rhs = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecB);
   jxf_SeqVectorPrint(ser_rhs, "./f_ser");
   
   jxf_Vector *ser_sol = NULL;
   ser_sol = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecX);
   jxf_SeqVectorPrint(ser_sol, "./u0_ser"); 

   jxf_CSRMatrixDestroy(ser_matrix);
   jxf_SeqVectorDestroy(ser_rhs);
   jxf_SeqVectorDestroy(ser_sol);
   
   jxf_printf("\n\n zzy kdim = %d  data saved!!\n\n", k_dim);
   jxf_printf(" Continue to solve? Type 'enter' to go on!!\n\n"); getchar();
#endif 

   (gmres_data -> converged) = 0;
   
  /*------------------------------------------------------------------------
   * With relative change convergence test on, it is possible to attempt
   * another iteration with a zero residual. This causes the parameter
   * alpha to go NaN. The guard_zero_residual parameter is to circumvent
   * this. Perhaps it should be set to something non-zero (but small).
   *-----------------------------------------------------------------------*/
   
   guard_zero_residual = 0.0;

   (*(gmres_functions->CommInfo))(matA, &my_id, &num_procs);
   
   if ( logging > 0 || print_level > 0 )
   {
      norms = (gmres_data -> norms);
   }

   /* initialize work arrays */
   rs = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions); 
   c  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   s  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   if (rel_change) 
   {
      rs_2 = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions);
   }
   

   hh = jxf_CTAllocF(JXF_Real*, k_dim + 1, gmres_functions); 
   for (i = 0; i < k_dim + 1; i ++)
   {	
      hh[i] = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   }

   (*(gmres_functions->CopyVector))(vecB, vecP[0]);

   /* compute initial residual */
   (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecP[0]);

   b_norm = sqrt((*(gmres_functions->InnerProd))(vecB, vecB));

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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm = sqrt((*(gmres_functions->InnerProd))(vecP[0], vecP[0]));
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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   /* peghoty  2010/06/23 */
   normdown = r_norm;

   if ( logging > 0 || print_level > 0 )
   {
      norms[0] = r_norm;
      if ( print_level > 0 && my_id == 0 )
      {
  	 jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
         {
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         }
         jxf_printf("Initial L2 norm of residual: %e\n", r_norm);
      }
   }
   
   iter = 0;

   if (b_norm > 0.0)
   {
      /* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
      den_norm = b_norm;
   }
   else
   {
      /* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
      den_norm = r_norm;
   }


   /* convergence criteria: |r_i| <= max( a_tol, r_tol * den_norm)
      den_norm = |r_0| or |b|
      note: default for a_tol is 0.0, so relative residual criteria is used unless
            user specifies a_tol, or sets r_tol = 0.0, which means absolute
            tol only is checked  */
      
   epsilon = jxf_max(a_tol, r_tol*den_norm);
   
   /* so now our stop criteria is |r_i| <= epsilon */

   if ( print_level > 0 && my_id == 0 )
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


   /* once the rel. change check has passed, we do not want to check it again */
   rel_change_passed = 0;


   /* outer iteration cycle */
   
   while (iter < max_iter)
   {
        /* initialize first term of hessenberg system */

	rs[0] = r_norm;
        if (r_norm == 0.0)
        {
           jxf_TFreeF(c, gmres_functions); 
           jxf_TFreeF(s, gmres_functions); 
           jxf_TFreeF(rs, gmres_functions);
           if (rel_change)  
           {
              jxf_TFreeF(rs_2,gmres_functions);
           }
           for (i = 0; i < k_dim + 1; i ++) 
           {
              jxf_TFreeF(hh[i], gmres_functions);
           }
           jxf_TFreeF(hh, gmres_functions); 
	   return jxf_error_flag;
	}

        /* see if we are already converged and should print the final norm and exit */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           if (!rel_change) /* shouldn't exit after no iterations if relative change is on */
           {
              (*(gmres_functions->CopyVector))(vecB, vecR);
              (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
              r_norm = sqrt((*(gmres_functions->InnerProd))(vecR, vecR));
              if (r_norm <= epsilon)
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n\n");
                    jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
                 }
                 (gmres_data -> converged) = 1; /* added by peghoty, 2011/11/08 */
                 break;
              }
              else
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf(" false convergence 1:");
                    if (b_norm > 0.0)
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                    }
                    else
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                    }
                 }
                 
                 /* added by Xu Xiaowen, 2011/11/07 */
                 if (is_check_restarted == 0) 
                 {
                    break;
                 }
              }
           }
	}


      	t = 1.0 / r_norm;
	(*(gmres_functions->ScaleVector))(t, vecP[0]);
	i = 0;

        /*** RESTART CYCLE (right-preconditioning) ***/
        
        while (i < k_dim && iter < max_iter)
	{
           i ++;
           iter ++;
           (*(gmres_functions->ClearVector))(vecR);
           //jxf_printf("precond: 1\n");
           precond(precond_data, preOperator, vecP[i-1], vecR);
           (*(gmres_functions->Matvec))(matvec_data, 1.0, matA, vecR, 0.0, vecP[i]);
           
           /* modified Gram_Schmidt */
           for (j = 0; j < i; j ++)
           {
              hh[j][i-1] = (*(gmres_functions->InnerProd))(vecP[j], vecP[i]);
              (*(gmres_functions->Axpy))(-hh[j][i-1],vecP[j], vecP[i]);
           }
           t = sqrt((*(gmres_functions->InnerProd))(vecP[i], vecP[i]));
           hh[i][i-1] = t;
           if (t != 0.0)
           {
              t = 1.0 / t;
              (*(gmres_functions->ScaleVector))(t, vecP[i]);
           }
           
           /* done with modified Gram_schmidt and Arnoldi step. update factorization of hh */
           for (j = 1; j < i; j ++)
           {
              t = hh[j-1][i-1];
              hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
              hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
           }
           t = hh[i][i-1]*hh[i][i-1];
           t += hh[i-1][i-1]*hh[i-1][i-1];
           gamma = sqrt(t);
           if (gamma == 0.0) gamma = epsmac;
           c[i-1] = hh[i-1][i-1]/gamma;
           s[i-1] = hh[i][i-1]/gamma;
           rs[i] = -hh[i][i-1]*rs[i-1];
           rs[i] /=  gamma;
           rs[i-1] = c[i-1]*rs[i-1];
           /* determine residual norm */
           hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
           r_norm = fabs(rs[i]);

           /* peghoty 2010/06/23 */
           normup = r_norm;
           conv_factor = normup / normdown;
           //jxf_printf("\n conv_factor = %lf\n\n",conv_factor);
           normdown = normup;

           /* print ? */
           if ( print_level > 0 )
           {
              norms[iter] = r_norm;
              if ( print_level > 0 && my_id == 0 )
              {
                 if (b_norm > 0.0)
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1],norms[iter]/b_norm);
                 }
                 else
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1], norms[iter]/r_norm_0);
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
 
           /* convergence factor tolerance */
           if (cf_tol > 0.0)
           {
              cf_ave_0 = cf_ave_1;
              cf_ave_1 = pow( r_norm / r_norm_0, 1.0 / (2.0*iter));
              
              weight   = fabs(cf_ave_1 - cf_ave_0);
              weight   = weight / jxf_max(cf_ave_1, cf_ave_0);
              weight   = 1.0 - weight;
              
              if (weight * cf_ave_1 > cf_tol) 
              {
                 break_value = 1;
                 break;
              }
           }
           
           /* should we exit the restart cycle? (conv. check) */
           if (r_norm <= epsilon && iter >= min_iter)
           {
              if (rel_change && !rel_change_passed)
              {
                 
                 /* To decide whether to break here: to actually
                  determine the relative change requires the approx
                  solution (so a triangular solve) and a
                  precond. solve - so if we have to do this many
                  times, it will be expensive...(unlike cg where is
                  is relatively straightforward)

                  previously, the intent (there was a bug), was to
                  exit the restart cycle based on the residual norm
                  and check the relative change outside the cycle.
                  Here we will check the relative here as we don't
                  want to exit the restart cycle prematurely */
                 
                 for (k = 0; k < i; k ++) /* extra copy of rs so we don't need to change the later solve */
                 {
                    rs_2[k] = rs[k];
                 }

                 /* solve tri. system */
                 rs_2[i-1] = rs_2[i-1] / hh[i-1][i-1];
                 for (k = i - 2; k >= 0; k --)
                 {
                    t = 0.0;
                    for (j = k + 1; j < i; j ++)
                    {
                       t -= hh[k][j]*rs_2[j];
                    }
                    t += rs_2[k];
                    rs_2[k] = t / hh[k][k];
                 }
                 
                 (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
                 (*(gmres_functions->ScaleVector))(rs_2[i-1], vecW);
                 for (j = i - 2; j >= 0; j --)
                 {
                    (*(gmres_functions->Axpy))(rs_2[j], vecP[j], vecW);
                 }
                    
                 (*(gmres_functions->ClearVector))(vecR);
                 
                 /* find correction (in r) */
                 //jxf_printf("precond: 2\n");
                 precond(precond_data, preOperator, vecW, vecR);
                 
                 /* copy current solution (x) to w (don't want to over-write x) */
                 (*(gmres_functions->CopyVector))(vecX, vecW);

                 /* add the correction */
                 (*(gmres_functions->Axpy))(1.0, vecR, vecW);

                 /* now w is the approx solution  - get the norm*/
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecW, vecW) );

                 if ( !(x_norm <= guard_zero_residual) ) /* don't divide by zero */
                 {  
                    /* now get x_i - x_i-1 */
                    if (num_rel_change_check)
                    {
                       /* have already checked once so we can avoid another precond. solve */
                       (*(gmres_functions->CopyVector))(vecW, vecR);
                       (*(gmres_functions->Axpy))(-1.0, vecW2, vecR);
                       
                       /* now r contains x_i - x_i-1*/

                       /* save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);
                    }
                    else
                    {
                       /* first time to check rel change*/

                       /* first save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);

                       /* for relative change take x_(i-1) to be 
                          x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                          Now
                          x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                          - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                          = M^{-1} rs_{i-1}{p_{i-1}} */
                       
                       (*(gmres_functions->ClearVector))(vecW);
                       (*(gmres_functions->Axpy))(rs_2[i-1], vecP[i-1], vecW);
                       (*(gmres_functions->ClearVector))(vecR);
                       /* apply the preconditioner */
                       //jxf_printf("precond: 3\n");
                       precond(precond_data, preOperator, vecW, vecR);
                       /* now r contains x_i - x_i-1 */          
                    }
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error = w_norm / x_norm;
                    if (relative_error <= r_tol)
                    {
                       rel_change_passed = 1;
                       break;
                    }
                 }
                 else
                 {
                    rel_change_passed = 1;
                    break;

                 }
                 num_rel_change_check ++;
              }
              else /* no relative change */
              {
                 break;
              }
           }
           
	} /*** end of restart cycle ***/

	/* now compute solution, first solve upper triangular system */

	if (break_value) break;
	
	rs[i-1] = rs[i-1] / hh[i-1][i-1];
	for (k = i - 2; k >= 0; k --)
	{
           t = 0.0;
           for (j = k + 1; j < i; j ++)
           {
              t -= hh[k][j]*rs[j];
           }
           t += rs[k];
           rs[k] = t / hh[k][k];
	}

        (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
        (*(gmres_functions->ScaleVector))(rs[i-1], vecW);
        for (j = i - 2; j >= 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j], vecP[j], vecW);
        }

	(*(gmres_functions->ClearVector))(vecR);
	
	/* find correction (in r) */
        //jxf_printf("precond: 4\n");
        precond(precond_data, preOperator, vecW, vecR);

        /* update current solution x (in x) */
	(*(gmres_functions->Axpy))(1.0, vecR, vecX);
         

        /* check for convergence by evaluating the actual residual */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           /* calculate actual residual norm */
           (*(gmres_functions->CopyVector))(vecB, vecR);
           (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
           r_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );

           if (r_norm <= epsilon)
           {
              if (rel_change && !rel_change_passed) /* calculate the relative change */
              {

                 /* calculate the norm of the solution */
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecX,vecX) );
               
                 if ( !(x_norm <= guard_zero_residual )) /* don't divide by zero */
                 {
                    /* for relative change take x_(i-1) to be 
                       x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                       Now
                       x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                       - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                       = M^{-1} rs_{i-1}{p_{i-1}} */
                    (*(gmres_functions->ClearVector))(vecW);
                    (*(gmres_functions->Axpy))(rs[i-1], vecP[i-1], vecW);
                    (*(gmres_functions->ClearVector))(vecR);
                    
                    /* apply the preconditioner */
                    //jxf_printf("precond: 5\n");
                    precond(precond_data, preOperator, vecW, vecR);
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error= w_norm / x_norm;
                    if ( relative_error < r_tol )
                    {
                       (gmres_data -> converged) = 1;
                       if ( print_level > 0 && my_id == 0 )
                       {
                          jxf_printf("\n");
                          jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                       }
                       break;
                    }
                 }
                 else
                 {
                    (gmres_data -> converged) = 1;
                    if ( print_level > 0 && my_id == 0 )
                    {
                       jxf_printf("\n");
                       jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                    }
                    break;
                 }

              }
              else /* don't need to check rel. change */
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n");
                    jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                 }
                 (gmres_data -> converged) = 1;
                 break;
              }
           }
           else /* conv. has not occurred, according to true residual */ 
           {
              if ( print_level > 0 && my_id == 0)
              {
                 jxf_printf(" false convergence 2:");
                 if (b_norm > 0.0)
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                 }
                 else
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                 }
              }
              
              /* modified by Xu Xiaowen, 2011/11/07 */
              if (is_check_restarted == 1) 
              {
                 (*(gmres_functions->CopyVector))(vecR, vecP[0]);
                 i = 0;
              } 
              else 
              {
                 break;  
              }
           }
	} /* end of convergence check */

        if (iter > 1 && conv_factor > convfac_threshold)
        {
           break;
        }

        /* compute residual vector and continue loop */
	for (j = i; j > 0; j --)
	{
           rs[j-1] = -s[j-1]*rs[j];
           rs[j] = c[j-1]*rs[j];
	}
        
        if (i) 
        {
           (*(gmres_functions->Axpy))(rs[i]-1.0, vecP[i], vecP[i]);
        }
        
        for (j = i - 1 ; j > 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j],vecP[j],vecP[i]);
        }
        
        if (i)
        {
           (*(gmres_functions->Axpy))(rs[0]-1.0, vecP[0], vecP[0]);
           (*(gmres_functions->Axpy))(1.0, vecP[i], vecP[0]);
        }
        
   } /* END of iteration while loop */


   if ( print_level > 1 && my_id == 0 ) jxf_printf("\n\n"); 

   (gmres_data -> num_iterations) = iter;
   
   if (b_norm > 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / b_norm;
   }
   
   if (b_norm == 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / r_norm_0;
   }

#if 0
   if (iter >= max_iter && r_norm > epsilon) 
   {
      jxf_printf(" Warning: Not reaching the given tolerance in %d iterations!!\n\n", max_iter);
      jxf_error(JXF_ERROR_CONV);
   }
#endif 

#if 0
   /* 保存解向量 */   
   jxf_Vector *ser_app = NULL;
   ser_app = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecX);
   jxf_SeqVectorPrint(ser_app, "./u_ser"); 
   jxf_SeqVectorDestroy(ser_app);
   jxf_printf("\n\n zzy solution has been saved!!\n\n");
#endif

#if 0
   /* 修改解向量 */ 
   JXF_Int change_sol;
   jxf_printf("\n");
   jxf_printf(" Do you want the modify the solution with the file './app'?\n");
   jxf_printf(" >>> 0: NO\n");
   jxf_printf(" >>> 1: YES\n");
   jxf_scanf("%d", &change_sol); 
   if (change_sol == 1) 
   {
      FILE *fp = NULL;
      jxf_Vector *tmp_app = NULL;
      JXF_Int ii, num_rows;
      tmp_app = jxf_ParVectorLocalVector((jxf_ParVector *)vecX);
   
      fp = fopen("./app", "r");
      jxf_fscanf(fp, "%d", &num_rows);
      for (ii = 0; ii < num_rows; ii ++)
      {
         jxf_fscanf(fp, "%le", &jxf_VectorData(tmp_app)[ii]);
      }
      jxf_printf("\n\n zzy solution has been changed for test!!\n\n");
   }
#endif
  
   jxf_TFreeF(c,gmres_functions); 
   jxf_TFreeF(s,gmres_functions); 
   jxf_TFreeF(rs,gmres_functions);
   if (rel_change)  
   {
      jxf_TFreeF(rs_2,gmres_functions);
   }

   for (i = 0; i < k_dim + 1; i ++)
   {	
      jxf_TFreeF(hh[i], gmres_functions);
   }
   jxf_TFreeF(hh, gmres_functions); 

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESAdaptiveSolve
 * \brief Adaptive between AMG, ILU(0) and combined preconditioner
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
jxf_GMRESAdaptiveSolve( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX )
{
   jxf_GMRESData      *gmres_data      = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = (gmres_data -> functions);
   JXF_Int 		     k_dim            = (gmres_data -> k_dim);
   JXF_Int               min_iter         = (gmres_data -> min_iter);
   JXF_Int 		     max_iter         = (gmres_data -> max_iter);
   JXF_Int               rel_change       = (gmres_data -> rel_change);
   JXF_Real 	     r_tol            = (gmres_data -> tol);
   JXF_Real 	     cf_tol           = (gmres_data -> cf_tol);
   JXF_Real            a_tol            = (gmres_data -> a_tol);
   void             *matvec_data      = (gmres_data -> matvec_data);

   JXF_Real            convfac_threshold   = (gmres_data -> convfac_threshold); // peghoty 2010/06/23
   JXF_Real            resdown_0_threshold = (gmres_data -> resdown_0_threshold); /* Yue Xiaoqiang, 2014/03/24 */
   JXF_Real            convfac_threshold_2 = (gmres_data -> convfac_threshold_2); /* Yue Xiaoqiang, 2014/03/24 */

   void             *vecR             = (gmres_data -> r);
   void             *vecW             = (gmres_data -> w);
   void             *vecW2            = (gmres_data -> w_2); /* vecW2 is only allocated if rel_change = 1 */
   void            **vecP             = (gmres_data -> p);
   JXF_Int 	           (*precond)()       = (gmres_functions -> precond);
   JXF_Int 	            *precond_data     = (gmres_data -> precond_data);
   JXF_Int             print_level        = (gmres_data -> print_level);
   JXF_Int             logging            = (gmres_data -> logging);
   JXF_Real         *norms              = (gmres_data -> norms);
   JXF_Int             is_check_restarted = (gmres_data -> is_check_restarted); /* Xu Xiaowen, 2011/11/08 */
   
   JXF_Int          *part_aux  = NULL;
   JXF_Int          *part_res  = NULL;
   jxf_ParVector *aux_vec   = NULL;
   jxf_ParVector *res_vec   = NULL;
   JXF_Solver     amg_data  = NULL;
   JXF_Int amg_max_levels, amg_relax_type, amg_print_level;
   JXF_Int amg_P_max_elmts, amg_measure_type, amg_interp_type;
   JXF_Int amg_coarsen_type, amg_agg_num_levels, amg_coarse_threshold;
   JXF_Real amg_strong_threshold;
   
   JXF_Int        break_value = 0;
   JXF_Int        break_adaptive = 0;
   JXF_Real     r_current_norm;
   JXF_Int	      i, j, k;
   JXF_Real     *rs, **hh, *c, *s, *rs_2 = NULL; 
   JXF_Int        iter; 
   JXF_Int        my_id, num_procs;
   JXF_Real     epsilon, gamma, t, r_norm, b_norm, den_norm, x_norm;
   JXF_Real     w_norm;
   
   JXF_Real     epsmac = 1.e-16; 
   JXF_Real     ieee_check = 0.0;

   JXF_Real     guard_zero_residual; 
   JXF_Real     cf_ave_0 = 0.0;
   JXF_Real     cf_ave_1 = 0.0;
   JXF_Real     weight;
   JXF_Real     r_norm_0;
   JXF_Real     relative_error = 1.0;

   JXF_Int        rel_change_passed = 0, num_rel_change_check = 0;

   /* peghoty 2010/06/23 */
   JXF_Real  normup       = 0.0;
   JXF_Real  normdown     = 0.0;
   JXF_Real  conv_factor  = 0.0;

#if 0
   /* 保存离散系统和初始解向量 */
   jxf_CSRMatrix *ser_matrix = NULL;
   ser_matrix = jxf_ParCSRMatrixToCSRMatrixAll((jxf_ParCSRMatrix *)matA);
   jxf_CSRMatrixPrint(ser_matrix,"./A_ser");
   
   jxf_Vector *ser_rhs = NULL;
   ser_rhs = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecB);
   jxf_SeqVectorPrint(ser_rhs, "./f_ser");
   
   jxf_Vector *ser_sol = NULL;
   ser_sol = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecX);
   jxf_SeqVectorPrint(ser_sol, "./u0_ser"); 

   jxf_CSRMatrixDestroy(ser_matrix);
   jxf_SeqVectorDestroy(ser_rhs);
   jxf_SeqVectorDestroy(ser_sol);
   
   jxf_printf("\n\n zzy kdim = %d  data saved!!\n\n", k_dim);
   jxf_printf(" Continue to solve? Type 'enter' to go on!!\n\n"); getchar();
#endif

   if (jxf_CombinedPrecondDataPreID((jxf_CombinedPrecondData *)precond_data) == 3) // regular PGMRES
   {

   (gmres_data -> converged) = 0;
   
  /*------------------------------------------------------------------------
   * With relative change convergence test on, it is possible to attempt
   * another iteration with a zero residual. This causes the parameter
   * alpha to go NaN. The guard_zero_residual parameter is to circumvent
   * this. Perhaps it should be set to something non-zero (but small).
   *-----------------------------------------------------------------------*/
   
   guard_zero_residual = 0.0;

   (*(gmres_functions->CommInfo))(matA, &my_id, &num_procs);
   
   if ( logging > 0 || print_level > 0 )
   {
      norms = (gmres_data -> norms);
   }

   /* initialize work arrays */
   rs = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions); 
   c  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   s  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   if (rel_change) 
   {
      rs_2 = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions);
   }
   

   hh = jxf_CTAllocF(JXF_Real*, k_dim + 1, gmres_functions); 
   for (i = 0; i < k_dim + 1; i ++)
   {	
      hh[i] = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   }

   (*(gmres_functions->CopyVector))(vecB, vecP[0]);

   /* compute initial residual */
   (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecP[0]);

   b_norm = sqrt((*(gmres_functions->InnerProd))(vecB, vecB));

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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm = sqrt((*(gmres_functions->InnerProd))(vecP[0], vecP[0]));
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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   /* peghoty  2010/06/23 */
   normdown = r_norm;

   if ( logging > 0 || print_level > 0 )
   {
      norms[0] = r_norm;
      if ( print_level > 0 && my_id == 0 )
      {
  	 jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
         {
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         }
         jxf_printf("Initial L2 norm of residual: %e\n", r_norm);
      }
   }
   
   iter = 0;

   if (b_norm > 0.0)
   {
      /* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
      den_norm = b_norm;
   }
   else
   {
      /* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
      den_norm = r_norm;
   }


   /* convergence criteria: |r_i| <= max( a_tol, r_tol * den_norm)
      den_norm = |r_0| or |b|
      note: default for a_tol is 0.0, so relative residual criteria is used unless
            user specifies a_tol, or sets r_tol = 0.0, which means absolute
            tol only is checked  */
      
   epsilon = jxf_max(a_tol, r_tol*den_norm);
   
   /* so now our stop criteria is |r_i| <= epsilon */

   if ( print_level > 0 && my_id == 0 )
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


   /* once the rel. change check has passed, we do not want to check it again */
   rel_change_passed = 0;


   /* outer iteration cycle */
   
   while (iter < max_iter)
   {
        /* initialize first term of hessenberg system */

	rs[0] = r_norm;
        if (r_norm == 0.0)
        {
           jxf_TFreeF(c, gmres_functions); 
           jxf_TFreeF(s, gmres_functions); 
           jxf_TFreeF(rs, gmres_functions);
           if (rel_change)  
           {
              jxf_TFreeF(rs_2,gmres_functions);
           }
           for (i = 0; i < k_dim + 1; i ++) 
           {
              jxf_TFreeF(hh[i], gmres_functions);
           }
           jxf_TFreeF(hh, gmres_functions); 
	   return jxf_error_flag;
	}

        /* see if we are already converged and should print the final norm and exit */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           if (!rel_change) /* shouldn't exit after no iterations if relative change is on */
           {
              (*(gmres_functions->CopyVector))(vecB, vecR);
              (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
              r_norm = sqrt((*(gmres_functions->InnerProd))(vecR, vecR));
              if (r_norm <= epsilon)
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n\n");
                    jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
                 }
                 (gmres_data -> converged) = 1; /* added by peghoty, 2011/11/08 */
                 break;
              }
              else
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf(" false convergence 1:");
                    if (b_norm > 0.0)
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                    }
                    else
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                    }
                 }
                 
                 /* added by Xu Xiaowen, 2011/11/07 */
                 if (is_check_restarted == 0) 
                 {
                    break;
                 }
              }
           }
	}


      	t = 1.0 / r_norm;
	(*(gmres_functions->ScaleVector))(t, vecP[0]);
	i = 0;

        /*** RESTART CYCLE (right-preconditioning) ***/
        
        while (i < k_dim && iter < max_iter)
	{
           i ++;
           iter ++;
           (*(gmres_functions->ClearVector))(vecR);
           //jxf_printf("precond: 1\n");
           precond(precond_data, preOperator, vecP[i-1], vecR);
           (*(gmres_functions->Matvec))(matvec_data, 1.0, matA, vecR, 0.0, vecP[i]);
           
           /* modified Gram_Schmidt */
           for (j = 0; j < i; j ++)
           {
              hh[j][i-1] = (*(gmres_functions->InnerProd))(vecP[j], vecP[i]);
              (*(gmres_functions->Axpy))(-hh[j][i-1],vecP[j], vecP[i]);
           }
           t = sqrt((*(gmres_functions->InnerProd))(vecP[i], vecP[i]));
           hh[i][i-1] = t;
           if (t != 0.0)
           {
              t = 1.0 / t;
              (*(gmres_functions->ScaleVector))(t, vecP[i]);
           }
           
           /* done with modified Gram_schmidt and Arnoldi step. update factorization of hh */
           for (j = 1; j < i; j ++)
           {
              t = hh[j-1][i-1];
              hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
              hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
           }
           t = hh[i][i-1]*hh[i][i-1];
           t += hh[i-1][i-1]*hh[i-1][i-1];
           gamma = sqrt(t);
           if (gamma == 0.0) gamma = epsmac;
           c[i-1] = hh[i-1][i-1]/gamma;
           s[i-1] = hh[i][i-1]/gamma;
           rs[i] = -hh[i][i-1]*rs[i-1];
           rs[i] /=  gamma;
           rs[i-1] = c[i-1]*rs[i-1];
           /* determine residual norm */
           hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
           r_norm = fabs(rs[i]);

           /* peghoty 2010/06/23 */
           normup = r_norm;
           conv_factor = normup / normdown;
           //jxf_printf("\n conv_factor = %lf\n\n",conv_factor);
           normdown = normup;

           /* print ? */
           if ( print_level > 0 )
           {
              norms[iter] = r_norm;
              if ( print_level > 0 && my_id == 0 )
              {
                 if (b_norm > 0.0)
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1],norms[iter]/b_norm);
                 }
                 else
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1], norms[iter]/r_norm_0);
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
 
           /* convergence factor tolerance */
           if (cf_tol > 0.0)
           {
              cf_ave_0 = cf_ave_1;
              cf_ave_1 = pow( r_norm / r_norm_0, 1.0 / (2.0*iter));
              
              weight   = fabs(cf_ave_1 - cf_ave_0);
              weight   = weight / jxf_max(cf_ave_1, cf_ave_0);
              weight   = 1.0 - weight;
              
              if (weight * cf_ave_1 > cf_tol) 
              {
                 break_value = 1;
                 break;
              }
           }
           
           /* should we exit the restart cycle? (conv. check) */
           if (r_norm <= epsilon && iter >= min_iter)
           {
              if (rel_change && !rel_change_passed)
              {
                 
                 /* To decide whether to break here: to actually
                  determine the relative change requires the approx
                  solution (so a triangular solve) and a
                  precond. solve - so if we have to do this many
                  times, it will be expensive...(unlike cg where is
                  is relatively straightforward)

                  previously, the intent (there was a bug), was to
                  exit the restart cycle based on the residual norm
                  and check the relative change outside the cycle.
                  Here we will check the relative here as we don't
                  want to exit the restart cycle prematurely */
                 
                 for (k = 0; k < i; k ++) /* extra copy of rs so we don't need to change the later solve */
                 {
                    rs_2[k] = rs[k];
                 }

                 /* solve tri. system */
                 rs_2[i-1] = rs_2[i-1] / hh[i-1][i-1];
                 for (k = i - 2; k >= 0; k --)
                 {
                    t = 0.0;
                    for (j = k + 1; j < i; j ++)
                    {
                       t -= hh[k][j]*rs_2[j];
                    }
                    t += rs_2[k];
                    rs_2[k] = t / hh[k][k];
                 }
                 
                 (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
                 (*(gmres_functions->ScaleVector))(rs_2[i-1], vecW);
                 for (j = i - 2; j >= 0; j --)
                 {
                    (*(gmres_functions->Axpy))(rs_2[j], vecP[j], vecW);
                 }
                    
                 (*(gmres_functions->ClearVector))(vecR);
                 
                 /* find correction (in r) */
                 //jxf_printf("precond: 2\n");
                 precond(precond_data, preOperator, vecW, vecR);
                 
                 /* copy current solution (x) to w (don't want to over-write x) */
                 (*(gmres_functions->CopyVector))(vecX, vecW);

                 /* add the correction */
                 (*(gmres_functions->Axpy))(1.0, vecR, vecW);

                 /* now w is the approx solution  - get the norm*/
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecW, vecW) );

                 if ( !(x_norm <= guard_zero_residual) ) /* don't divide by zero */
                 {  
                    /* now get x_i - x_i-1 */
                    if (num_rel_change_check)
                    {
                       /* have already checked once so we can avoid another precond. solve */
                       (*(gmres_functions->CopyVector))(vecW, vecR);
                       (*(gmres_functions->Axpy))(-1.0, vecW2, vecR);
                       
                       /* now r contains x_i - x_i-1*/

                       /* save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);
                    }
                    else
                    {
                       /* first time to check rel change*/

                       /* first save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);

                       /* for relative change take x_(i-1) to be 
                          x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                          Now
                          x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                          - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                          = M^{-1} rs_{i-1}{p_{i-1}} */
                       
                       (*(gmres_functions->ClearVector))(vecW);
                       (*(gmres_functions->Axpy))(rs_2[i-1], vecP[i-1], vecW);
                       (*(gmres_functions->ClearVector))(vecR);
                       /* apply the preconditioner */
                       //jxf_printf("precond: 3\n");
                       precond(precond_data, preOperator, vecW, vecR);
                       /* now r contains x_i - x_i-1 */          
                    }
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error = w_norm / x_norm;
                    if (relative_error <= r_tol)
                    {
                       rel_change_passed = 1;
                       break;
                    }
                 }
                 else
                 {
                    rel_change_passed = 1;
                    break;

                 }
                 num_rel_change_check ++;
              }
              else /* no relative change */
              {
                 break;
              }
           }
           
	} /*** end of restart cycle ***/

        if (iter > 1 && conv_factor > convfac_threshold)
        {
           break;
        }

	/* now compute solution, first solve upper triangular system */

	if (break_value) break;
	
	rs[i-1] = rs[i-1] / hh[i-1][i-1];
	for (k = i - 2; k >= 0; k --)
	{
           t = 0.0;
           for (j = k + 1; j < i; j ++)
           {
              t -= hh[k][j]*rs[j];
           }
           t += rs[k];
           rs[k] = t / hh[k][k];
	}

        (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
        (*(gmres_functions->ScaleVector))(rs[i-1], vecW);
        for (j = i - 2; j >= 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j], vecP[j], vecW);
        }

	(*(gmres_functions->ClearVector))(vecR);
	
	/* find correction (in r) */
        //jxf_printf("precond: 4\n");
        precond(precond_data, preOperator, vecW, vecR);

        /* update current solution x (in x) */
	(*(gmres_functions->Axpy))(1.0, vecR, vecX);
         

        /* check for convergence by evaluating the actual residual */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           /* calculate actual residual norm */
           (*(gmres_functions->CopyVector))(vecB, vecR);
           (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
           r_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );

           if (r_norm <= epsilon)
           {
              if (rel_change && !rel_change_passed) /* calculate the relative change */
              {

                 /* calculate the norm of the solution */
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecX,vecX) );
               
                 if ( !(x_norm <= guard_zero_residual )) /* don't divide by zero */
                 {
                    /* for relative change take x_(i-1) to be 
                       x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                       Now
                       x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                       - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                       = M^{-1} rs_{i-1}{p_{i-1}} */
                    (*(gmres_functions->ClearVector))(vecW);
                    (*(gmres_functions->Axpy))(rs[i-1], vecP[i-1], vecW);
                    (*(gmres_functions->ClearVector))(vecR);
                    
                    /* apply the preconditioner */
                    //jxf_printf("precond: 5\n");
                    precond(precond_data, preOperator, vecW, vecR);
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error= w_norm / x_norm;
                    if ( relative_error < r_tol )
                    {
                       (gmres_data -> converged) = 1;
                       if ( print_level > 0 && my_id == 0 )
                       {
                          jxf_printf("\n");
                          jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                       }
                       break;
                    }
                 }
                 else
                 {
                    (gmres_data -> converged) = 1;
                    if ( print_level > 0 && my_id == 0 )
                    {
                       jxf_printf("\n");
                       jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                    }
                    break;
                 }

              }
              else /* don't need to check rel. change */
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n");
                    jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                 }
                 (gmres_data -> converged) = 1;
                 break;
              }
           }
           else /* conv. has not occurred, according to true residual */ 
           {
              if ( print_level > 0 && my_id == 0)
              {
                 jxf_printf(" false convergence 2:");
                 if (b_norm > 0.0)
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                 }
                 else
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                 }
              }
              
              /* modified by Xu Xiaowen, 2011/11/07 */
              if (is_check_restarted == 1) 
              {
                 (*(gmres_functions->CopyVector))(vecR, vecP[0]);
                 i = 0;
              } 
              else 
              {
                 break;  
              }
           }
	} /* end of convergence check */

        /* compute residual vector and continue loop */
	for (j = i; j > 0; j --)
	{
           rs[j-1] = -s[j-1]*rs[j];
           rs[j] = c[j-1]*rs[j];
	}
        
        if (i) 
        {
           (*(gmres_functions->Axpy))(rs[i]-1.0, vecP[i], vecP[i]);
        }
        
        for (j = i - 1 ; j > 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j],vecP[j],vecP[i]);
        }
        
        if (i)
        {
           (*(gmres_functions->Axpy))(rs[0]-1.0, vecP[0], vecP[0]);
           (*(gmres_functions->Axpy))(1.0, vecP[i], vecP[0]);
        }
        
   } /* END of iteration while loop */


   if ( print_level > 1 && my_id == 0 ) jxf_printf("\n\n"); 

   (gmres_data -> num_iterations) = iter;
   (gmres_data -> comp_amg_iter)  = iter;
   (gmres_data -> break_adaptive) = break_adaptive;
   
   if (b_norm > 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / b_norm;
   }
   
   if (b_norm == 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / r_norm_0;
   }

   }
   else // adaptive PGMRES
   {

   (gmres_data -> converged) = 0;
   
   MPI_Comm comm = jxf_CombinedPrecondDataComm((jxf_CombinedPrecondData *)precond_data);
   jxf_MPI_Comm_size(comm, &num_procs);

  /*------------------------------------------------------------------------
   * With relative change convergence test on, it is possible to attempt
   * another iteration with a zero residual. This causes the parameter
   * alpha to go NaN. The guard_zero_residual parameter is to circumvent
   * this. Perhaps it should be set to something non-zero (but small).
   *-----------------------------------------------------------------------*/
   
   guard_zero_residual = 0.0;

   (*(gmres_functions->CommInfo))(matA, &my_id, &num_procs);
   
   if ( logging > 0 || print_level > 0 )
   {
      norms = (gmres_data -> norms);
   }

   /* initialize work arrays */
   rs = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions); 
   c  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   s  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   if (rel_change) 
   {
      rs_2 = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions);
   }
   

   hh = jxf_CTAllocF(JXF_Real*, k_dim + 1, gmres_functions); 
   for (i = 0; i < k_dim + 1; i ++)
   {	
      hh[i] = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   }

   (*(gmres_functions->CopyVector))(vecB, vecP[0]);

   /* compute initial residual */
   (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecP[0]);

   b_norm = sqrt((*(gmres_functions->InnerProd))(vecB, vecB));

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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESAdaptiveSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm = sqrt((*(gmres_functions->InnerProd))(vecP[0], vecP[0]));
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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESAdaptiveSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   /* peghoty  2010/06/23 */
   normdown = r_norm;

   if ( logging > 0 || print_level > 0 )
   {
      norms[0] = r_norm;
      if ( print_level > 0 && my_id == 0 )
      {
  	 jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
         {
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         }
         jxf_printf("Initial L2 norm of residual: %e\n", r_norm);
      }
   }
   
   iter = 0;

   if (b_norm > 0.0)
   {
      /* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
      den_norm = b_norm;
   }
   else
   {
      /* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
      den_norm = r_norm;
   }


   /* convergence criteria: |r_i| <= max( a_tol, r_tol * den_norm)
      den_norm = |r_0| or |b|
      note: default for a_tol is 0.0, so relative residual criteria is used unless
            user specifies a_tol, or sets r_tol = 0.0, which means absolute
            tol only is checked  */
      
   epsilon = jxf_max(a_tol, r_tol*den_norm);
   
   /* so now our stop criteria is |r_i| <= epsilon */

   if ( print_level > 0 && my_id == 0 )
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


   /* once the rel. change check has passed, we do not want to check it again */
   rel_change_passed = 0;


   /* outer iteration cycle */
   
   while (iter < max_iter)
   {
        /* initialize first term of hessenberg system */

	rs[0] = r_norm;
	r_current_norm = r_norm;
        if (r_norm == 0.0)
        {
           jxf_TFreeF(c, gmres_functions); 
           jxf_TFreeF(s, gmres_functions); 
           jxf_TFreeF(rs, gmres_functions);
           if (rel_change)  
           {
              jxf_TFreeF(rs_2,gmres_functions);
           }
           for (i = 0; i < k_dim + 1; i ++) 
           {
              jxf_TFreeF(hh[i], gmres_functions);
           }
           jxf_TFreeF(hh, gmres_functions); 
	   return jxf_error_flag;
	}

        /* see if we are already converged and should print the final norm and exit */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           if (!rel_change) /* shouldn't exit after no iterations if relative change is on */
           {
              (*(gmres_functions->CopyVector))(vecB, vecR);
              (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
              r_norm = sqrt((*(gmres_functions->InnerProd))(vecR, vecR));
              if (r_norm <= epsilon)
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n\n");
                    jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
                 }
                 (gmres_data -> converged) = 1; /* added by peghoty, 2011/11/08 */
                 break;
              }
              else
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf(" false convergence 1:");
                    if (b_norm > 0.0)
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                    }
                    else
                    {
                       jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                    }
                 }
                 
                 /* added by Xu Xiaowen, 2011/11/07 */
                 if (is_check_restarted == 0) 
                 {
                    break;
                 }
              }
           }
	}

      	t = 1.0 / r_norm;
	(*(gmres_functions->ScaleVector))(t, vecP[0]);
	i = 0;

        /*** RESTART CYCLE (right-preconditioning) ***/
        
        while (i < k_dim && iter < max_iter)
	{
           // here AMG SETUP phase
           if (break_adaptive == 1)
           {
           amg_max_levels = jxf_CombinedPrecondDataAMGMaxLevels((jxf_CombinedPrecondData *)precond_data);
           amg_relax_type = jxf_CombinedPrecondDataAMGRelaxType((jxf_CombinedPrecondData *)precond_data);
           amg_print_level = jxf_CombinedPrecondDataAMGPrintLevel((jxf_CombinedPrecondData *)precond_data);
           amg_interp_type = jxf_CombinedPrecondDataAMGInterpType((jxf_CombinedPrecondData *)precond_data);
           amg_P_max_elmts = jxf_CombinedPrecondDataAMGPMaxElmts((jxf_CombinedPrecondData *)precond_data);
           amg_measure_type = jxf_CombinedPrecondDataAMGMeasureType((jxf_CombinedPrecondData *)precond_data);
           amg_coarsen_type = jxf_CombinedPrecondDataAMGCoarsenType((jxf_CombinedPrecondData *)precond_data);
           amg_agg_num_levels = jxf_CombinedPrecondDataAMGAggNumLevels((jxf_CombinedPrecondData *)precond_data);
           amg_coarse_threshold = jxf_CombinedPrecondDataAMGCoarseThreshold((jxf_CombinedPrecondData *)precond_data);
           amg_strong_threshold = jxf_CombinedPrecondDataAMGStrongThreshold((jxf_CombinedPrecondData *)precond_data);
              JXF_PAMGCreate(&amg_data);
              JXF_PAMGSetMaxLevels(amg_data, amg_max_levels);
              JXF_PAMGSetMaxIter(amg_data, 1);
              JXF_PAMGSetMeasureType(amg_data, amg_measure_type);
              JXF_PAMGSetCoarsenType(amg_data, amg_coarsen_type);
              JXF_PAMGSetInterpType(amg_data, amg_interp_type);
              JXF_PAMGSetPMaxElmts(amg_data, amg_P_max_elmts);
              JXF_PAMGSetAggNumLevels(amg_data, amg_agg_num_levels);
              JXF_PAMGSetStrongThreshold(amg_data, amg_strong_threshold);
              JXF_PAMGSetPrintLevel(amg_data, amg_print_level);
              JXF_PAMGSetCoarseThreshold(amg_data, amg_coarse_threshold);
              JXF_PAMGSetCycleNumSweeps(amg_data, 1, 1);
              JXF_PAMGSetCycleNumSweeps(amg_data, 1, 2);
              JXF_PAMGSetCycleNumSweeps(amg_data, 1, 3);
              JXF_PAMGSetCycleRelaxType(amg_data, amg_relax_type, 1);
              JXF_PAMGSetCycleRelaxType(amg_data, amg_relax_type, 2);
              JXF_PAMGSetCycleRelaxType(amg_data, 9, 3);
              JXF_PAMGSetup(amg_data, (JXF_hpCSRMatrix)preOperator);
              // Allocate memory for auxiliary and residual vectors
              jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), &part_aux);
              aux_vec = jxf_ParVectorCreate(jxf_hpCSRMatrixComm((jxf_hpCSRMatrix *)preOperator),
                             jxf_hpCSRMatrixGlobalNumRows((jxf_hpCSRMatrix *)preOperator), part_aux);
              jxf_ParVectorInitialize(aux_vec);
              jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), &part_res);
              res_vec = jxf_ParVectorCreate(jxf_hpCSRMatrixComm((jxf_hpCSRMatrix *)preOperator),
                             jxf_hpCSRMatrixGlobalNumRows((jxf_hpCSRMatrix *)preOperator), part_res);
              jxf_ParVectorInitialize(res_vec);
              break_adaptive = 2;
           }
           
           i ++;
           iter ++;
           (*(gmres_functions->ClearVector))(vecR);
           //jxf_printf("precond: 1\n");
           precond(precond_data, preOperator, vecP[i-1], vecR);
           if (break_adaptive == 2)
           {
              jxf_ParVectorCopy(vecP[i-1], res_vec);
              jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), vecR, 1.0, res_vec);
              jxf_ParVectorSetConstantValues(aux_vec, 0.0);
              JXF_PAMGPrecond(amg_data, (JXF_hpCSRMatrix)preOperator, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
              jxf_ParVectorAxpy(1.0, aux_vec, vecR);
           }
           (*(gmres_functions->Matvec))(matvec_data, 1.0, matA, vecR, 0.0, vecP[i]);
           
           /* modified Gram_Schmidt */
           for (j = 0; j < i; j ++)
           {
              hh[j][i-1] = (*(gmres_functions->InnerProd))(vecP[j], vecP[i]);
              (*(gmres_functions->Axpy))(-hh[j][i-1], vecP[j], vecP[i]);
           }
           t = sqrt((*(gmres_functions->InnerProd))(vecP[i], vecP[i]));
           hh[i][i-1] = t;
           if (t != 0.0)
           {
              t = 1.0 / t;
              (*(gmres_functions->ScaleVector))(t, vecP[i]);
           }
           
           /* done with modified Gram_schmidt and Arnoldi step. update factorization of hh */
           for (j = 1; j < i; j ++)
           {
              t = hh[j-1][i-1];
              hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
              hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
           }
           t = hh[i][i-1]*hh[i][i-1];
           t += hh[i-1][i-1]*hh[i-1][i-1];
           gamma = sqrt(t);
           if (gamma == 0.0) gamma = epsmac;
           c[i-1] = hh[i-1][i-1]/gamma;
           s[i-1] = hh[i][i-1]/gamma;
           rs[i] = -hh[i][i-1]*rs[i-1];
           rs[i] /=  gamma;
           rs[i-1] = c[i-1]*rs[i-1];
           /* determine residual norm */
           hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
           r_norm = fabs(rs[i]);

           /* print ? */
           if ( print_level > 0 )
           {
              norms[iter] = r_norm;
              if ( print_level > 0 && my_id == 0 )
              {
                 if (b_norm > 0.0)
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1],norms[iter]/b_norm);
                 }
                 else
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                               norms[iter]/norms[iter-1], norms[iter]/r_norm_0);
                 }
              }
           }

           /* peghoty 2010/06/23 */
           normup = r_norm;
           conv_factor = normup / normdown;
           //jxf_printf("\n conv_factor = %lf\n\n",conv_factor);
           if (b_norm > 0.0)
           {
               if ((break_adaptive == 0) && ((normup/b_norm > resdown_0_threshold)
                                         || ((iter > 1) && (conv_factor > convfac_threshold_2))))
               {
                  (gmres_data -> comp_ilu_iter) = iter;
                  break_adaptive = 1;
                  i = 0;
                  iter = 0;
                  rs[0] = r_current_norm;
                  continue;
               }
           }
           else
           {
               if ((break_adaptive == 0) && ((normup/r_norm_0 > resdown_0_threshold)
                                         || ((iter > 1) && (conv_factor > convfac_threshold_2))))
               {
                  (gmres_data -> comp_ilu_iter) = iter;
                  break_adaptive = 1;
                  i = 0;
                  iter = 0;
                  rs[0] = r_current_norm;
                  continue;
               }
           }
           normdown = normup;
          
           if (iter > 1 && conv_factor > convfac_threshold)
           {
              if (my_id == 0 && print_level > 0)
              {
                 jxf_printf("\n  Warning: iteration has terminated because the\n");
                 jxf_printf("           current conv_factor > convfac_threshold(%lf) !!\n\n",convfac_threshold);
              }
              break;
           }
 
           /* convergence factor tolerance */
           if (cf_tol > 0.0)
           {
              cf_ave_0 = cf_ave_1;
              cf_ave_1 = pow( r_norm / r_norm_0, 1.0 / (2.0*iter));
              
              weight   = fabs(cf_ave_1 - cf_ave_0);
              weight   = weight / jxf_max(cf_ave_1, cf_ave_0);
              weight   = 1.0 - weight;
              
              if (weight * cf_ave_1 > cf_tol) 
              {
                 break_value = 1;
                 break;
              }
           }
           
           /* should we exit the restart cycle? (conv. check) */
           if (r_norm <= epsilon && iter >= min_iter)
           {
              if (rel_change && !rel_change_passed)
              {
                 
                 /* To decide whether to break here: to actually
                  determine the relative change requires the approx
                  solution (so a triangular solve) and a
                  precond. solve - so if we have to do this many
                  times, it will be expensive...(unlike cg where is
                  is relatively straightforward)

                  previously, the intent (there was a bug), was to
                  exit the restart cycle based on the residual norm
                  and check the relative change outside the cycle.
                  Here we will check the relative here as we don't
                  want to exit the restart cycle prematurely */
                 
                 for (k = 0; k < i; k ++) /* extra copy of rs so we don't need to change the later solve */
                 {
                    rs_2[k] = rs[k];
                 }

                 /* solve tri. system */
                 rs_2[i-1] = rs_2[i-1] / hh[i-1][i-1];
                 for (k = i - 2; k >= 0; k --)
                 {
                    t = 0.0;
                    for (j = k + 1; j < i; j ++)
                    {
                       t -= hh[k][j]*rs_2[j];
                    }
                    t += rs_2[k];
                    rs_2[k] = t / hh[k][k];
                 }
                 
                 (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
                 (*(gmres_functions->ScaleVector))(rs_2[i-1], vecW);
                 for (j = i - 2; j >= 0; j --)
                 {
                    (*(gmres_functions->Axpy))(rs_2[j], vecP[j], vecW);
                 }
                    
                 (*(gmres_functions->ClearVector))(vecR);
                 
                 /* find correction (in r) */
                 //jxf_printf("precond: 2\n");
                 precond(precond_data, preOperator, vecW, vecR);
                 if (break_adaptive == 2)
                 {
                    jxf_ParVectorCopy(vecW, res_vec);
                    jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), vecR, 1.0, res_vec);
                    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
                    JXF_PAMGPrecond(amg_data, (JXF_hpCSRMatrix)preOperator,
                                             (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
                    jxf_ParVectorAxpy(1.0, aux_vec, vecR);
                 }
                 
                 /* copy current solution (x) to w (don't want to over-write x) */
                 (*(gmres_functions->CopyVector))(vecX, vecW);

                 /* add the correction */
                 (*(gmres_functions->Axpy))(1.0, vecR, vecW);

                 /* now w is the approx solution  - get the norm*/
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecW, vecW) );

                 if ( !(x_norm <= guard_zero_residual) ) /* don't divide by zero */
                 {  
                    /* now get x_i - x_i-1 */
                    if (num_rel_change_check)
                    {
                       /* have already checked once so we can avoid another precond. solve */
                       (*(gmres_functions->CopyVector))(vecW, vecR);
                       (*(gmres_functions->Axpy))(-1.0, vecW2, vecR);
                       
                       /* now r contains x_i - x_i-1*/

                       /* save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);
                    }
                    else
                    {
                       /* first time to check rel change*/

                       /* first save current soln w in w_2 for next time */
                       (*(gmres_functions->CopyVector))(vecW, vecW2);

                       /* for relative change take x_(i-1) to be 
                          x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                          Now
                          x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                          - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                          = M^{-1} rs_{i-1}{p_{i-1}} */
                       
                       (*(gmres_functions->ClearVector))(vecW);
                       (*(gmres_functions->Axpy))(rs_2[i-1], vecP[i-1], vecW);
                       (*(gmres_functions->ClearVector))(vecR);
                       /* apply the preconditioner */
                       //jxf_printf("precond: 3\n");
                       precond(precond_data, preOperator, vecW, vecR);
                       if (break_adaptive == 2)
                       {
                          jxf_ParVectorCopy(vecW, res_vec);
                          jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), vecR, 1.0, res_vec);
                          jxf_ParVectorSetConstantValues(aux_vec, 0.0);
                          JXF_PAMGPrecond(amg_data, (JXF_hpCSRMatrix)preOperator,
                                                   (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
                          jxf_ParVectorAxpy(1.0, aux_vec, vecR);
                       }
                       /* now r contains x_i - x_i-1 */          
                    }
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error = w_norm / x_norm;
                    if (relative_error <= r_tol)
                    {
                       rel_change_passed = 1;
                       break;
                    }
                 }
                 else
                 {
                    rel_change_passed = 1;
                    break;

                 }
                 num_rel_change_check ++;
              }
              else /* no relative change */
              {
                 break;
              }
           }
           
	} /*** end of restart cycle ***/

        if (iter > 1 && conv_factor > convfac_threshold)
        {
           break;
        }

	/* now compute solution, first solve upper triangular system */

	if (break_value) break;
	
	rs[i-1] = rs[i-1] / hh[i-1][i-1];
	for (k = i - 2; k >= 0; k --)
	{
           t = 0.0;
           for (j = k + 1; j < i; j ++)
           {
              t -= hh[k][j]*rs[j];
           }
           t += rs[k];
           rs[k] = t / hh[k][k];
	}

        (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
        (*(gmres_functions->ScaleVector))(rs[i-1], vecW);
        for (j = i - 2; j >= 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j], vecP[j], vecW);
        }

	(*(gmres_functions->ClearVector))(vecR);
	
	/* find correction (in r) */
        //jxf_printf("precond: 4\n");
        precond(precond_data, preOperator, vecW, vecR);
        if (break_adaptive == 2)
        {
           jxf_ParVectorCopy(vecW, res_vec);
           jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), vecR, 1.0, res_vec);
           jxf_ParVectorSetConstantValues(aux_vec, 0.0);
           JXF_PAMGPrecond(amg_data, (JXF_hpCSRMatrix)preOperator, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
           jxf_ParVectorAxpy(1.0, aux_vec, vecR);
        }

        /* update current solution x (in x) */
	(*(gmres_functions->Axpy))(1.0, vecR, vecX);
         

        /* check for convergence by evaluating the actual residual */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           /* calculate actual residual norm */
           (*(gmres_functions->CopyVector))(vecB, vecR);
           (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
           r_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );

           if (r_norm <= epsilon)
           {
              if (rel_change && !rel_change_passed) /* calculate the relative change */
              {

                 /* calculate the norm of the solution */
                 x_norm = sqrt( (*(gmres_functions->InnerProd))(vecX,vecX) );
               
                 if ( !(x_norm <= guard_zero_residual )) /* don't divide by zero */
                 {
                    /* for relative change take x_(i-1) to be 
                       x + M^{-1}[sum{j=0..i-2} rs_j p_j ]. 
                       Now
                       x_i - x_{i-1}= {x + M^{-1}[sum{j=0..i-1} rs_j p_j ]}
                       - {x + M^{-1}[sum{j=0..i-2} rs_j p_j ]}
                       = M^{-1} rs_{i-1}{p_{i-1}} */
                    (*(gmres_functions->ClearVector))(vecW);
                    (*(gmres_functions->Axpy))(rs[i-1], vecP[i-1], vecW);
                    (*(gmres_functions->ClearVector))(vecR);
                    
                    /* apply the preconditioner */
                    //jxf_printf("precond: 5\n");
                    precond(precond_data, preOperator, vecW, vecR);
                    if (break_adaptive == 2)
                    {
                       jxf_ParVectorCopy(vecW, res_vec);
                       jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)preOperator), vecR, 1.0, res_vec);
                       jxf_ParVectorSetConstantValues(aux_vec, 0.0);
                       JXF_PAMGPrecond(amg_data, (JXF_hpCSRMatrix)preOperator,
                                                (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
                       jxf_ParVectorAxpy(1.0, aux_vec, vecR);
                    }
                    
                    /* find the norm of x_i - x_i-1 */          
                    w_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );
                    relative_error= w_norm / x_norm;
                    if ( relative_error < r_tol )
                    {
                       (gmres_data -> converged) = 1;
                       if ( print_level > 0 && my_id == 0 )
                       {
                          jxf_printf("\n");
                          jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                       }
                       break;
                    }
                 }
                 else
                 {
                    (gmres_data -> converged) = 1;
                    if ( print_level > 0 && my_id == 0 )
                    {
                       jxf_printf("\n");
                       jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                    }
                    break;
                 }

              }
              else /* don't need to check rel. change */
              {
                 if ( print_level > 0 && my_id == 0 )
                 {
                    jxf_printf("\n");
                    jxf_printf("Final L2 norm of residual: %e\n", r_norm);
                 }
                 (gmres_data -> converged) = 1;
                 break;
              }
           }
           else /* conv. has not occurred, according to true residual */ 
           {
              if ( print_level > 0 && my_id == 0)
              {
                 jxf_printf(" false convergence 2:");
                 if (b_norm > 0.0)
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||b|| = %.2le\n", r_norm, r_norm / b_norm);
                 }
                 else
                 {
                    jxf_printf(" ||r|| = %.2le  ||r||/||r_0|| = %.2le\n", r_norm, r_norm / r_norm_0);
                 }
              }
              
              /* modified by Xu Xiaowen, 2011/11/07 */
              if (is_check_restarted == 1) 
              {
                 (*(gmres_functions->CopyVector))(vecR, vecP[0]);
                 i = 0;
              } 
              else 
              {
                 break;  
              }
           }
	} /* end of convergence check */

        /* compute residual vector and continue loop */
	for (j = i; j > 0; j --)
	{
           rs[j-1] = -s[j-1]*rs[j];
           rs[j] = c[j-1]*rs[j];
	}
        
        if (i) 
        {
           (*(gmres_functions->Axpy))(rs[i]-1.0, vecP[i], vecP[i]);
        }
        
        for (j = i - 1 ; j > 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j],vecP[j],vecP[i]);
        }
        
        if (i)
        {
           (*(gmres_functions->Axpy))(rs[0]-1.0, vecP[0], vecP[0]);
           (*(gmres_functions->Axpy))(1.0, vecP[i], vecP[0]);
        }
        
   } /* END of iteration while loop */


   if ( print_level > 1 && my_id == 0 ) jxf_printf("\n\n"); 

   (gmres_data -> num_iterations) = iter;
   if (break_adaptive == 0)
   {
      (gmres_data -> comp_ilu_iter) = iter;
   }
   else
   {
      (gmres_data -> comp_bco_iter) = iter;
   }
   (gmres_data -> break_adaptive) = break_adaptive;
   
   if (b_norm > 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / b_norm;
   }
   
   if (b_norm == 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / r_norm_0;
   }

   if ( amg_data )
   {
      JXF_PAMGDestroy( amg_data );
   }

   if ( aux_vec )
   {
      jxf_ParVectorDestroy( aux_vec );
   }

   if ( res_vec )
   {
      jxf_ParVectorDestroy( res_vec );
   }

   }

#if 0
   if (iter >= max_iter && r_norm > epsilon) 
   {
      jxf_printf(" Warning: Not reaching the given tolerance in %d iterations!!\n\n", max_iter);
      jxf_error(JXF_ERROR_CONV);
   }
#endif

#if 0
   /* 保存解向量 */   
   jxf_Vector *ser_app = NULL;
   ser_app = jxf_ParVectorToVectorAll ((jxf_ParVector *)vecX);
   jxf_SeqVectorPrint(ser_app, "./u_ser"); 
   jxf_SeqVectorDestroy(ser_app);
   jxf_printf("\n\n zzy solution has been saved!!\n\n");
#endif

#if 0
   /* 修改解向量 */ 
   JXF_Int change_sol;
   jxf_printf("\n");
   jxf_printf(" Do you want the modify the solution with the file './app'?\n");
   jxf_printf(" >>> 0: NO\n");
   jxf_printf(" >>> 1: YES\n");
   jxf_scanf("%d", &change_sol); 
   if (change_sol == 1) 
   {
      FILE *fp = NULL;
      jxf_Vector *tmp_app = NULL;
      JXF_Int ii, num_rows;
      tmp_app = jxf_ParVectorLocalVector((jxf_ParVector *)vecX);
   
      fp = fopen("./app", "r");
      jxf_fscanf(fp, "%d", &num_rows);
      for (ii = 0; ii < num_rows; ii ++)
      {
         jxf_fscanf(fp, "%le", &jxf_VectorData(tmp_app)[ii]);
      }
      jxf_printf("\n\n zzy solution has been changed for test!!\n\n");
   }
#endif

   jxf_TFreeF(c,gmres_functions); 
   jxf_TFreeF(s,gmres_functions); 
   jxf_TFreeF(rs,gmres_functions);
   if (rel_change)  
   {
      jxf_TFreeF(rs_2,gmres_functions);
   }
   for (i = 0; i < k_dim + 1; i ++)
   {	
      jxf_TFreeF(hh[i], gmres_functions);
   }
   jxf_TFreeF(hh, gmres_functions);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_GMRESSolveDefault
 * \note This id a simplified version with default parameters. 
 * \author peghoty
 * \date 2011/11/08
 */
JXF_Int
jxf_GMRESSolveDefault( void *gmres_vdata, void *preOperator, void *matA, void *vecB, void *vecX )
{
   jxf_GMRESData      *gmres_data      = gmres_vdata;
   jxf_GMRESFunctions *gmres_functions = gmres_data -> functions;
   JXF_Int 		     k_dim            = (gmres_data -> k_dim);
   JXF_Int               min_iter         = (gmres_data -> min_iter);
   JXF_Int 		     max_iter         = (gmres_data -> max_iter);
   JXF_Real 	     r_tol            = (gmres_data -> tol);
   JXF_Real            a_tol            = (gmres_data -> a_tol);
   void             *matvec_data      = (gmres_data -> matvec_data);
   void             *vecR             = (gmres_data -> r);
   void             *vecW             = (gmres_data -> w);
   void            **vecP             = (gmres_data -> p);
   JXF_Int 	           (*precond)()       = (gmres_functions -> precond);
   JXF_Int 	            *precond_data     = (gmres_data -> precond_data);
   JXF_Int             print_level        = (gmres_data -> print_level);
   JXF_Int             logging            = (gmres_data -> logging);
   JXF_Real         *norms              = (gmres_data -> norms);

   JXF_Int	      i, j, k;
   JXF_Real     *rs, **hh, *c, *s; 
   JXF_Int        iter; 
   JXF_Int        my_id, num_procs;
   JXF_Real     epsilon, gamma, t, r_norm, b_norm, den_norm;
   
   JXF_Real     epsmac = 1.e-16; 
   JXF_Real     ieee_check = 0.0;

   JXF_Real     r_norm_0;

   (gmres_data -> converged) = 0;
   
   (*(gmres_functions->CommInfo))(matA, &my_id, &num_procs);
   
   if ( logging > 0 || print_level > 0 )
   {
      norms = (gmres_data -> norms);
   }

   /* initialize work arrays */
   rs = jxf_CTAllocF(JXF_Real, k_dim + 1, gmres_functions); 
   c  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   s  = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   hh = jxf_CTAllocF(JXF_Real*, k_dim + 1, gmres_functions); 
   for (i = 0; i < k_dim + 1; i ++)
   {	
      hh[i] = jxf_CTAllocF(JXF_Real, k_dim, gmres_functions); 
   }

   (*(gmres_functions->CopyVector))(vecB, vecP[0]);

   /* compute initial residual */
   (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecP[0]);

   b_norm = sqrt((*(gmres_functions->InnerProd))(vecB, vecB));

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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm = sqrt((*(gmres_functions->InnerProd))(vecP[0], vecP[0]));
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
         jxf_printf("\n\nERROR detected by JXFPAMG ... BEGIN\n");
         jxf_printf("ERROR -- jxf_GMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   if ( logging > 0 || print_level > 0 )
   {
      norms[0] = r_norm;
      if ( print_level > 0 && my_id == 0 )
      {
  	 jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
         {
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         }
         jxf_printf("Initial L2 norm of residual: %e\n", r_norm);
      }
   }
   
   iter = 0;

   if (b_norm > 0.0)
   {
      /* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
      den_norm = b_norm;
   }
   else
   {
      /* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
      den_norm = r_norm;
   }


   /* convergence criteria: |r_i| <= max( a_tol, r_tol * den_norm)
      den_norm = |r_0| or |b|
      note: default for a_tol is 0.0, so relative residual criteria is used unless
            user specifies a_tol, or sets r_tol = 0.0, which means absolute
            tol only is checked  */
      
   epsilon = jxf_max(a_tol, r_tol*den_norm);
   
   /* so now our stop criteria is |r_i| <= epsilon */

   if ( print_level > 0 && my_id == 0 )
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


   /* outer iteration cycle */
   
   while (iter < max_iter)
   {
        /* initialize first term of hessenberg system */

	rs[0] = r_norm;
        if (r_norm == 0.0)
        {
           jxf_TFreeF(c, gmres_functions); 
           jxf_TFreeF(s, gmres_functions); 
           jxf_TFreeF(rs, gmres_functions);
           for (i = 0; i < k_dim + 1; i ++) 
           {
              jxf_TFreeF(hh[i], gmres_functions);
           }
           jxf_TFreeF(hh, gmres_functions); 
	   return jxf_error_flag;
	}

        /* see if we are already converged and should print the final norm and exit */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           (*(gmres_functions->CopyVector))(vecB, vecR);
           (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
           r_norm = sqrt((*(gmres_functions->InnerProd))(vecR, vecR));
           if (r_norm <= epsilon)
           {
              if ( print_level > 0 && my_id == 0 )
              {
                 jxf_printf("\n\n");
                 jxf_printf(" Final L2 norm of residual: %e\n\n", r_norm);
              }
              (gmres_data -> converged) = 1; // peghoty, 2011/11/08
              break;
           }
           else
           {
              if ( print_level > 0 && my_id == 0)
              {
                 jxf_printf(" false convergence 1\n");
              }
           }
	}

      	t = 1.0 / r_norm;
	(*(gmres_functions->ScaleVector))(t, vecP[0]);
	i = 0;

        /*** RESTART CYCLE (right-preconditioning) ***/
        
        while (i < k_dim && iter < max_iter)
	{
           i ++;
           iter ++;
           (*(gmres_functions->ClearVector))(vecR);
           precond(precond_data, preOperator, vecP[i-1], vecR);
           (*(gmres_functions->Matvec))(matvec_data, 1.0, matA, vecR, 0.0, vecP[i]);
           
           /* modified Gram_Schmidt */
           for (j = 0; j < i; j ++)
           {
              hh[j][i-1] = (*(gmres_functions->InnerProd))(vecP[j], vecP[i]);
              (*(gmres_functions->Axpy))(-hh[j][i-1],vecP[j], vecP[i]);
           }
           t = sqrt((*(gmres_functions->InnerProd))(vecP[i], vecP[i]));
           hh[i][i-1] = t;	
           if (t != 0.0)
           {
              t = 1.0 / t;
              (*(gmres_functions->ScaleVector))(t, vecP[i]);
           }
           
           /* done with modified Gram_schmidt and Arnoldi step. update factorization of hh */
           for (j = 1; j < i; j ++)
           {
              t = hh[j-1][i-1];
              hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
              hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
           }
           t = hh[i][i-1]*hh[i][i-1];
           t += hh[i-1][i-1]*hh[i-1][i-1];
           gamma = sqrt(t);
           if (gamma == 0.0) gamma = epsmac;
           c[i-1] = hh[i-1][i-1]/gamma;
           s[i-1] = hh[i][i-1]/gamma;
           rs[i] = -hh[i][i-1]*rs[i-1];
           rs[i] /=  gamma;
           rs[i-1] = c[i-1]*rs[i-1];
           /* determine residual norm */
           hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
           r_norm = fabs(rs[i]);

           /* print ? */
           if ( print_level > 0 )
           {
              norms[iter] = r_norm;
              if ( print_level > 0 && my_id == 0 )
              {
                 if (b_norm > 0.0)
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter],
                                              norms[iter]/norms[iter-1],norms[iter]/b_norm);
                 }
                 else
                 {
                    jxf_printf("% 5d    %e    %f   %e\n", iter, norms[iter], norms[iter]/norms[iter-1],norms[iter]/r_norm_0);
                 }
              }
           }
           
           /* should we exit the restart cycle? (conv. check) */
           if (r_norm <= epsilon && iter >= min_iter)
           {
              break;
           }
           
	} /*** end of restart cycle ***/


	/* now compute solution, first solve upper triangular system */

	rs[i-1] = rs[i-1] / hh[i-1][i-1];
	for (k = i - 2; k >= 0; k --)
	{
           t = 0.0;
           for (j = k + 1; j < i; j ++)
           {
              t -= hh[k][j]*rs[j];
           }
           t += rs[k];
           rs[k] = t / hh[k][k];
	}

        (*(gmres_functions->CopyVector))(vecP[i-1], vecW);
        (*(gmres_functions->ScaleVector))(rs[i-1], vecW);
        for (j = i - 2; j >= 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j], vecP[j], vecW);
        }

	(*(gmres_functions->ClearVector))(vecR);
	
	/* find correction (in r) */
        precond(precond_data, preOperator, vecW, vecR);

        /* update current solution x (in x) */
	(*(gmres_functions->Axpy))(1.0, vecR, vecX);
         

        /* check for convergence by evaluating the actual residual */
	if (r_norm <= epsilon && iter >= min_iter) 
        {
           /* calculate actual residual norm */
           (*(gmres_functions->CopyVector))(vecB, vecR);
           (*(gmres_functions->Matvec))(matvec_data, -1.0, matA, vecX, 1.0, vecR);
           r_norm = sqrt( (*(gmres_functions->InnerProd))(vecR,vecR) );

           if (r_norm <= epsilon)
           {
              if ( print_level > 0 && my_id == 0 )
              {
                 jxf_printf("\n");
                 jxf_printf(" Final L2 norm of residual: %e\n", r_norm);
              }
              (gmres_data -> converged) = 1;
              break;
           }
           else /* conv. has not occurred, according to true residual */ 
           {
              if ( print_level > 0 && my_id == 0)
              {
                 jxf_printf(" false convergence 2\n");
              }
              (*(gmres_functions->CopyVector))(vecR, vecP[0]);
              i = 0;
           }
	} /* end of convergence check */

        /* compute residual vector and continue loop */
	for (j = i; j > 0; j --)
	{
           rs[j-1] = -s[j-1]*rs[j];
           rs[j] = c[j-1]*rs[j];
	}
        
        if (i) 
        {
           (*(gmres_functions->Axpy))(rs[i]-1.0, vecP[i], vecP[i]);
        }
        
        for (j = i - 1 ; j > 0; j --)
        {
           (*(gmres_functions->Axpy))(rs[j],vecP[j],vecP[i]);
        }
        
        if (i)
        {
           (*(gmres_functions->Axpy))(rs[0]-1.0, vecP[0], vecP[0]);
           (*(gmres_functions->Axpy))(1.0, vecP[i], vecP[0]);
        }
        
   } /* END of iteration while loop */


   if ( print_level > 1 && my_id == 0 ) jxf_printf("\n\n"); 

   (gmres_data -> num_iterations) = iter;
   
   if (b_norm > 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / b_norm;
   }
   
   if (b_norm == 0.0)
   {
      (gmres_data -> rel_residual_norm) = r_norm / r_norm_0;
   }

#if 0
   if (iter >= max_iter && r_norm > epsilon) 
   {
      jxf_printf(" Warning: Not reaching the given tolerance in %d iterations!!\n\n", max_iter);
      jxf_error(JXF_ERROR_CONV);
   }
#endif
  
   jxf_TFreeF(c,gmres_functions); 
   jxf_TFreeF(s,gmres_functions); 
   jxf_TFreeF(rs,gmres_functions);

   for (i = 0; i < k_dim + 1; i ++)
   {	
      jxf_TFreeF(hh[i], gmres_functions);
   }
   jxf_TFreeF(hh, gmres_functions); 

   return jxf_error_flag;
}
