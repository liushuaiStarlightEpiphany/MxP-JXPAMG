//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  pcogmres.c
 *  Date: 2019/02/12
 */ 

#include "jxf_combined.h"
#include "jxf_krylov.h"

JXF_Int
JXF_ParCSRCOGMRESCreate( MPI_Comm comm, JXF_Solver *solver )
{
   jxf_COGMRESFunctions * cogmres_functions;

   cogmres_functions = jxf_COGMRESFunctionsCreate(
         jxf_CAlloc, jxf_ParKrylovFree, jxf_ParKrylovCommInfo,
         jxf_ParKrylovCreateVector,
         jxf_ParKrylovCreateVectorArray,
         jxf_ParKrylovDestroyVector, jxf_ParKrylovMatvecCreate,
         jxf_ParKrylovMatvec, jxf_ParKrylovMatvecDestroy,
         jxf_ParKrylovInnerProd, jxf_ParKrylovMassInnerProd, 
         jxf_ParKrylovMassDotpTwo, jxf_ParKrylovCopyVector,
         //jxf_ParKrylovCopyVector,
         jxf_ParKrylovClearVector,
         jxf_ParKrylovScaleVector, jxf_ParKrylovAxpy, jxf_ParKrylovMassAxpy,
         jxf_ParKrylovIdentitySetup, jxf_ParKrylovIdentity );

   *solver = ( (JXF_Solver) jxf_COGMRESCreate( cogmres_functions ) );

   if (!solver)
   {
      jxf_error_in_arg(2);
   }

   return jxf_error_flag;
}

JXF_Int
JXF_ParCSRCOGMRESDestroy( JXF_Solver solver )
{
   return( jxf_COGMRESDestroy( (void *) solver ) );
}

JXF_Int
JXF_COGMRESSetup( JXF_Solver solver, JXF_Matrix A, JXF_Vector b, JXF_Vector x )
{
   return( jxf_COGMRESSetup( solver, A, b, x ) );
}

JXF_Int 
JXF_COGMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix A, JXF_Vector b, JXF_Vector x )
{
   return( jxf_COGMRESSolve( solver, preOperator, A, b, x ) );
}

JXF_Int
JXF_COGMRESSetKDim( JXF_Solver solver, JXF_Int k_dim )
{
   return( jxf_COGMRESSetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_COGMRESGetKDim( JXF_Solver solver, JXF_Int * k_dim )
{
   return( jxf_COGMRESGetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_COGMRESSetUnroll( JXF_Solver solver, JXF_Int unroll )
{
   return( jxf_COGMRESSetUnroll( (void *) solver, unroll ) );
}

JXF_Int
JXF_COGMRESGetUnroll( JXF_Solver solver, JXF_Int * unroll )
{
   return( jxf_COGMRESGetUnroll( (void *) solver, unroll ) );
}

JXF_Int
JXF_COGMRESSetCGS( JXF_Solver solver, JXF_Int cgs )
{
   return( jxf_COGMRESSetCGS( (void *) solver, cgs ) );
}

JXF_Int
JXF_COGMRESGetCGS( JXF_Solver solver, JXF_Int * cgs )
{
   return( jxf_COGMRESGetCGS( (void *) solver, cgs ) );
}

JXF_Int
JXF_COGMRESSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_COGMRESSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_COGMRESGetTol( JXF_Solver solver, JXF_Real * tol )
{
   return( jxf_COGMRESGetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_COGMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol )
{
   return( jxf_COGMRESSetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_COGMRESGetAbsoluteTol( JXF_Solver solver, JXF_Real * a_tol )
{
   return( jxf_COGMRESGetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_COGMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol )
{
   return( jxf_COGMRESSetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

JXF_Int
JXF_COGMRESGetConvergenceFactorTol( JXF_Solver solver, JXF_Real * cf_tol )
{
   return( jxf_COGMRESGetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

JXF_Int
JXF_COGMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter )
{
   return( jxf_COGMRESSetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_COGMRESGetMinIter( JXF_Solver solver, JXF_Int * min_iter )
{
   return( jxf_COGMRESGetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_COGMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_COGMRESSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_COGMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted )
{
   return( jxf_COGMRESSetIsCheckRestarted( (void *) solver, is_check_restarted ) );
}

JXF_Int
JXF_COGMRESGetMaxIter( JXF_Solver solver, JXF_Int * max_iter )
{
   return( jxf_COGMRESGetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_COGMRESSetPrecond( JXF_Solver          solver,
                      JXF_PtrToSolverFcn  precond,
                      JXF_PtrToSolverFcn  precond_setup,
                      JXF_Solver          precond_solver )
{
   return( jxf_COGMRESSetPrecond( (void *) solver,
                                 (JXF_Int (*)(void*, void*, void*, void*))precond,
                                 (JXF_Int (*)(void*, void*, void*, void*))precond_setup,
                                 (void *) precond_solver ) );
}

JXF_Int
JXF_COGMRESGetPrecond( JXF_Solver  solver, JXF_Solver *precond_data_ptr )
{
   return( jxf_COGMRESGetPrecond( (void *) solver, (JXF_Solver *) precond_data_ptr ) );
}

JXF_Int
JXF_COGMRESSetPrintLevel( JXF_Solver solver, JXF_Int level )
{
   return( jxf_COGMRESSetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_COGMRESGetPrintLevel( JXF_Solver solver, JXF_Int * level )
{
   return( jxf_COGMRESGetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_COGMRESSetLogging( JXF_Solver solver, JXF_Int level )
{
   return( jxf_COGMRESSetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_COGMRESGetLogging( JXF_Solver solver, JXF_Int * level )
{
   return( jxf_COGMRESGetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_COGMRESGetNumIterations( JXF_Solver  solver, JXF_Int *num_iterations )
{
   return( jxf_COGMRESGetNumIterations( (void *) solver, num_iterations ) );
}

JXF_Int
JXF_COGMRESGetConverged( JXF_Solver  solver, JXF_Int *converged )
{
   return( jxf_COGMRESGetConverged( (void *) solver, converged ) );
}

JXF_Int
JXF_COGMRESGetFinalRelativeResidualNorm( JXF_Solver  solver, JXF_Real *norm )
{
   return( jxf_COGMRESGetFinalRelativeResidualNorm( (void *) solver, norm ) );
}

JXF_Int
JXF_COGMRESGetResidual( JXF_Solver solver, void *residual )
{
   return jxf_COGMRESGetResidual( (void *) solver, (void **) residual );
}

JXF_Int
JXF_COGMRESSetModifyPC( JXF_Solver  solver, JXF_Int (*modify_pc)(JXF_Solver, JXF_Int, JXF_Real) )
{
   return jxf_COGMRESSetModifyPC( (void *) solver, (JXF_Int(*)(void*, JXF_Int, JXF_Real))modify_pc );
}

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
)
{
   jxf_COGMRESFunctions *cogmres_functions;
   cogmres_functions = (jxf_COGMRESFunctions *)CAlloc( 1, sizeof(jxf_COGMRESFunctions) );

   cogmres_functions->CAlloc            = CAlloc;
   cogmres_functions->Free              = Free;
   cogmres_functions->CommInfo          = CommInfo; /* not in PCGFunctionsCreate */
   cogmres_functions->CreateVector      = CreateVector;
   cogmres_functions->CreateVectorArray = CreateVectorArray; /* not in PCGFunctionsCreate */
   cogmres_functions->DestroyVector     = DestroyVector;
   cogmres_functions->MatvecCreate      = MatvecCreate;
   cogmres_functions->Matvec            = Matvec;
   cogmres_functions->MatvecDestroy     = MatvecDestroy;
   cogmres_functions->InnerProd         = InnerProd;
   cogmres_functions->MassInnerProd     = MassInnerProd;
   cogmres_functions->MassDotpTwo       = MassDotpTwo;
   cogmres_functions->CopyVector        = CopyVector;
   cogmres_functions->ClearVector       = ClearVector;
   cogmres_functions->ScaleVector       = ScaleVector;
   cogmres_functions->Axpy              = Axpy;
   cogmres_functions->MassAxpy          = MassAxpy;
   /* default preconditioner must be set here but can be changed later... */
   cogmres_functions->precond_setup     = PrecondSetup;
   cogmres_functions->precond           = Precond;

   return cogmres_functions;
}

void *
jxf_COGMRESCreate( jxf_COGMRESFunctions *cogmres_functions )
{
   jxf_COGMRESData *cogmres_data;

   cogmres_data = jxf_CTAllocF(jxf_COGMRESData, 1, cogmres_functions);

   cogmres_data->functions = cogmres_functions;

   /* set defaults */
   (cogmres_data -> k_dim)              = 5;
   (cogmres_data -> cgs)                = 1; /* if 2 performs reorthogonalization */
   (cogmres_data -> tol)                = 1.0e-06; /* relative residual tol */
   (cogmres_data -> cf_tol)             = 0.0;
   (cogmres_data -> a_tol)              = 0.0; /* abs. residual tol */
   (cogmres_data -> min_iter)           = 0;
   (cogmres_data -> max_iter)           = 1000;
   (cogmres_data -> rel_change)         = 0;
   (cogmres_data -> is_check_restarted) = 0;
   (cogmres_data -> skip_real_r_check)  = 0;
   (cogmres_data -> converged)          = 0;
   (cogmres_data -> precond_data)       = NULL;
   (cogmres_data -> print_level)        = 0;
   (cogmres_data -> logging)            = 0;
   (cogmres_data -> p)                  = NULL;
   (cogmres_data -> r)                  = NULL;
   (cogmres_data -> w)                  = NULL;
   (cogmres_data -> w_2)                = NULL;
   (cogmres_data -> matvec_data)        = NULL;
   (cogmres_data -> norms)              = NULL;
   (cogmres_data -> log_file_name)      = NULL;
   (cogmres_data -> unroll)             = 0;

   return (void *) cogmres_data;
}

JXF_Int
jxf_COGMRESDestroy( void *cogmres_vdata )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   JXF_Int i;

   if (cogmres_data)
   {
      jxf_COGMRESFunctions *cogmres_functions = cogmres_data->functions;
      if ( (cogmres_data->logging>0) || (cogmres_data->print_level) > 0 )
      {
         if ( (cogmres_data -> norms) != NULL )
         jxf_TFreeF( cogmres_data -> norms, cogmres_functions );
      }

      if ( (cogmres_data -> matvec_data) != NULL )
         (*(cogmres_functions->MatvecDestroy))(cogmres_data -> matvec_data);

      if ( (cogmres_data -> r) != NULL )
         (*(cogmres_functions->DestroyVector))(cogmres_data -> r);
      if ( (cogmres_data -> w) != NULL )
         (*(cogmres_functions->DestroyVector))(cogmres_data -> w);
      if ( (cogmres_data -> w_2) != NULL )
         (*(cogmres_functions->DestroyVector))(cogmres_data -> w_2);

      if ( (cogmres_data -> p) != NULL )
      {
         for (i = 0; i < (cogmres_data -> k_dim+1); i++)
         {
            if ( (cogmres_data -> p)[i] != NULL )
            (*(cogmres_functions->DestroyVector))( (cogmres_data -> p) [i]);
         }
         jxf_TFreeF( cogmres_data->p, cogmres_functions );
      }
      jxf_TFreeF( cogmres_data, cogmres_functions );
      jxf_TFreeF( cogmres_functions, cogmres_functions );
   }

   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetResidual( void *cogmres_vdata, void **residual )
{
   jxf_COGMRESData  *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *residual = cogmres_data->r;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetup( void *cogmres_vdata, void *A, void *b, void *x )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   jxf_COGMRESFunctions *cogmres_functions = cogmres_data->functions;

   JXF_Int k_dim            = (cogmres_data -> k_dim);
   JXF_Int max_iter         = (cogmres_data -> max_iter);
   JXF_Int (*precond_setup)(void*,void*,void*,void*) = (cogmres_functions->precond_setup);
   void       *precond_data   = (cogmres_data -> precond_data);
   JXF_Int rel_change       = (cogmres_data -> rel_change);

   (cogmres_data -> A) = A;

   /*--------------------------------------------------
    * The arguments for NewVector are important to
    * maintain consistency between the setup and
    * compute phases of matvec and the preconditioner.
    *--------------------------------------------------*/

   if ((cogmres_data -> p) == NULL)
      (cogmres_data -> p) = (void**)(*(cogmres_functions->CreateVectorArray))(k_dim+1,x);
   if ((cogmres_data -> r) == NULL)
      (cogmres_data -> r) = (*(cogmres_functions->CreateVector))(b);
   if ((cogmres_data -> w) == NULL)
      (cogmres_data -> w) = (*(cogmres_functions->CreateVector))(b);

   if (rel_change)
   {  
      if ((cogmres_data -> w_2) == NULL)
         (cogmres_data -> w_2) = (*(cogmres_functions->CreateVector))(b);
   }

   if ((cogmres_data -> matvec_data) == NULL)
      (cogmres_data -> matvec_data) = (*(cogmres_functions->MatvecCreate))(A, x);

   precond_setup(precond_data, A, b, x);

   /*-----------------------------------------------------
    * Allocate space for log info
    *-----------------------------------------------------*/

   if ( (cogmres_data->logging)>0 || (cogmres_data->print_level) > 0 )
   {
      if ((cogmres_data -> norms) == NULL)
         (cogmres_data -> norms) = jxf_CTAllocF(JXF_Real, max_iter + 1, cogmres_functions);
   }
   if ( (cogmres_data->print_level) > 0 ) 
   {
      if ((cogmres_data -> log_file_name) == NULL)
         (cogmres_data -> log_file_name) = (char*)"cogmres.out.log";
   }

   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSolve( void  *cogmres_vdata, void *preOperator, void  *A, void  *b, void  *x )
{
   jxf_COGMRESData      *cogmres_data      = (jxf_COGMRESData *)cogmres_vdata;
   jxf_COGMRESFunctions *cogmres_functions = cogmres_data->functions;
   JXF_Int     k_dim              = (cogmres_data -> k_dim);
   JXF_Int     unroll             = (cogmres_data -> unroll);
   JXF_Int     cgs                = (cogmres_data -> cgs);
   JXF_Int     min_iter           = (cogmres_data -> min_iter);
   JXF_Int     max_iter           = (cogmres_data -> max_iter);
   JXF_Int     rel_change         = (cogmres_data -> rel_change);
   JXF_Int     is_check_restarted = cogmres_data -> is_check_restarted;
   JXF_Int     skip_real_r_check  = (cogmres_data -> skip_real_r_check);
   JXF_Real    r_tol              = (cogmres_data -> tol);
   JXF_Real    cf_tol             = (cogmres_data -> cf_tol);
   JXF_Real    a_tol              = (cogmres_data -> a_tol);
   void         *matvec_data       = (cogmres_data -> matvec_data);

   void         *r                 = (cogmres_data -> r);
   void         *w                 = (cogmres_data -> w);
   /* note: w_2 is only allocated if rel_change = 1 */
   void         *w_2               = (cogmres_data -> w_2); 

   void        **p                 = (cogmres_data -> p);

   JXF_Int (*precond)(void*,void*,void*,void*) = (cogmres_functions -> precond);
   JXF_Int  *precond_data       = (JXF_Int*)(cogmres_data -> precond_data);

   JXF_Int print_level = (cogmres_data -> print_level);
   JXF_Int logging     = (cogmres_data -> logging);

   JXF_Real     *norms          = (cogmres_data -> norms);
  /* not used yet   char           *log_file_name  = (cogmres_data -> log_file_name);*/
  /*   FILE           *fp; */

   JXF_Int  break_value = 0;
   JXF_Int  i, j, k;
  /*KS: rv is the norm history */
   JXF_Real *rs, *hh, *uu, *c, *s, *rs_2 = NULL, *rv;
  //, *tmp; 
   JXF_Int  iter; 
   JXF_Int  my_id, num_procs;
   JXF_Real epsilon, gamma, t, r_norm, b_norm, den_norm, x_norm;
   JXF_Real w_norm;

   JXF_Real epsmac = 1.e-16; 
   JXF_Real ieee_check = 0.;

   JXF_Real guard_zero_residual; 
   JXF_Real cf_ave_0 = 0.0;
   JXF_Real cf_ave_1 = 0.0;
   JXF_Real weight;
   JXF_Real r_norm_0;
   JXF_Real relative_error = 1.0;

   JXF_Int        rel_change_passed = 0, num_rel_change_check = 0;
   JXF_Int    itmp = 0;

   JXF_Real real_r_norm_old, real_r_norm_new;

   (cogmres_data -> converged) = 0;
   /*-----------------------------------------------------------------------
    * With relative change convergence test on, it is possible to attempt
    * another iteration with a zero residual. This causes the parameter
    * alpha to go NaN. The guard_zero_residual parameter is to circumvent
    * this. Perhaps it should be set to something non-zero (but small).
    *-----------------------------------------------------------------------*/
   guard_zero_residual = 0.0;

   (*(cogmres_functions->CommInfo))(A,&my_id,&num_procs);
   if ( logging>0 || print_level>0 )
   {
      norms          = (cogmres_data -> norms);
   }

   /* initialize work arrays */
   rs = jxf_CTAllocF(JXF_Real,k_dim+1,cogmres_functions);
   c  = jxf_CTAllocF(JXF_Real,k_dim,cogmres_functions);
   s  = jxf_CTAllocF(JXF_Real,k_dim,cogmres_functions);
   if (rel_change) rs_2 = jxf_CTAllocF(JXF_Real, k_dim+1, cogmres_functions); 

   rv = jxf_CTAllocF(JXF_Real, k_dim+1, cogmres_functions);
  
   hh = jxf_CTAllocF(JXF_Real, (k_dim+1)*k_dim, cogmres_functions);
   uu = jxf_CTAllocF(JXF_Real, (k_dim+1)*k_dim, cogmres_functions);

   (*(cogmres_functions->CopyVector))(b,p[0]);

   /* compute initial residual */
   (*(cogmres_functions->Matvec))(matvec_data,-1.0, A, x, 1.0, p[0]);

   b_norm = sqrt((*(cogmres_functions->InnerProd))(b,b));
   real_r_norm_old = b_norm;

   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (b_norm != 0.) ieee_check = b_norm/b_norm; /* INF -> NaN conversion */
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
         jxf_printf("ERROR -- jxf_COGMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied b.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm   = sqrt((*(cogmres_functions->InnerProd))(p[0],p[0]));
   r_norm_0 = r_norm;

   /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
   if (r_norm != 0.) ieee_check = r_norm/r_norm; /* INF -> NaN conversion */
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
         jxf_printf("ERROR -- jxf_COGMRESSolve: INFs and/or NaNs detected in input.\n");
         jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
         jxf_printf("Returning error flag += 101.  Program not terminated.\n");
         jxf_printf("ERROR detected by JXFPAMG ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }
 
   if ( logging>0 || print_level > 0)
   {
      norms[0] = r_norm;
      if ( print_level>1 && my_id == 0 )
      {
         jxf_printf("L2 norm of b: %e\n", b_norm);
         if (b_norm == 0.0)
            jxf_printf("Rel_resid_norm actually contains the residual norm\n");
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
   };

   /* convergence criteria: |r_i| <= max( a_tol, r_tol * den_norm)
      den_norm = |r_0| or |b|
      note: default for a_tol is 0.0, so relative residual criteria is used unless
      user specifies a_tol, or sets r_tol = 0.0, which means absolute
      tol only is checked  */

   epsilon = jxf_max(a_tol,r_tol*den_norm);

   /* so now our stop criteria is |r_i| <= epsilon */

   if ( print_level>1 && my_id == 0 )
   {
      if (b_norm > 0.0)
      {
         jxf_printf("=============================================\n\n");
         jxf_printf("Iters     resid.norm     conv.rate  rel.res.norm\n");
         jxf_printf("-----    ------------    ---------- ------------\n");

      }
      else
      {
         jxf_printf("=============================================\n\n");
         jxf_printf("Iters     resid.norm     conv.rate\n");
         jxf_printf("-----    ------------    ----------\n");
      };
   }


   /* once the rel. change check has passed, we do not want to check it again */
   rel_change_passed = 0;

   while (iter < max_iter)
   {
      /* initialize first term of hessenberg system */
      rs[0] = r_norm;
      if (r_norm == 0.0)
      {
         jxf_TFreeF(c,cogmres_functions); 
         jxf_TFreeF(s,cogmres_functions); 
         jxf_TFreeF(rs,cogmres_functions);
         jxf_TFreeF(rv,cogmres_functions);
         if (rel_change)  jxf_TFreeF(rs_2,cogmres_functions);
         jxf_TFreeF(hh,cogmres_functions); 
         jxf_TFreeF(uu,cogmres_functions); 
         return jxf_error_flag;
      }

      /* see if we are already converged and 
         should print the final norm and exit */

      if (r_norm  <= epsilon && iter >= min_iter) 
      {
         if (!rel_change) /* shouldn't exit after no iterations if
                           * relative change is on*/
         {
            (*(cogmres_functions->CopyVector))(b,r);
            (*(cogmres_functions->Matvec))(matvec_data,-1.0,A,x,1.0,r);
            r_norm = sqrt((*(cogmres_functions->InnerProd))(r,r));
            if (r_norm  <= epsilon)
            {
               if ( print_level>1 && my_id == 0)
               {
                  jxf_printf("\n\n");
                  jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
               }
               break;
            }
            else
            {
               if ( print_level>0 && my_id == 0) jxf_printf("false convergence 1\n");
               if (is_check_restarted == 0) break;
            }
         }
      }

      t = 1.0 / r_norm;
      (*(cogmres_functions->ScaleVector))(t,p[0]);
      i = 0;
      /***RESTART CYCLE (right-preconditioning) ***/
      while (i < k_dim && iter < max_iter)
      {
         i++;
         iter++;
         itmp = (i-1)*(k_dim+1);

         (*(cogmres_functions->ClearVector))(r);
        
         precond(precond_data, preOperator, p[i-1], r);
         (*(cogmres_functions->Matvec))(matvec_data, 1.0, A, r, 0.0, p[i]);
         for (j=0; j<i; j++)
            rv[j]  = 0;

         if (cgs > 1)
         {
            (*(cogmres_functions->MassDotpTwo))((void *) p[i], p[i-1], p, i, unroll, &hh[itmp], &uu[itmp]);
            for (j=0; j<i-1; j++) uu[j*(k_dim+1)+i-1] = uu[itmp+j];
            for (j=0; j<i; j++) rv[j] = hh[itmp+j];
            for (k=0; k < i; k++)
            {
               for (j=0; j < i; j++)
               {
                  hh[itmp+j] -= (uu[k*(k_dim+1)+j]*rv[j]);
               }
            }
            for (j=0; j<i; j++)
               hh[itmp+j]  = -rv[j]-hh[itmp+j];
         }
         else
         {
            (*(cogmres_functions->MassInnerProd))((void *) p[i], p, i, unroll, &hh[itmp]);
            for (j=0; j<i; j++)
               hh[itmp+j]  = -hh[itmp+j];
         }

         (*(cogmres_functions->MassAxpy))(&hh[itmp],p,p[i], i, unroll);
         for (j=0; j<i; j++)
            hh[itmp+j]  = -hh[itmp+j];
         t = sqrt( (*(cogmres_functions->InnerProd))(p[i],p[i]) );
         hh[itmp+i] = t;

         if (hh[itmp+i] != 0.0)
         {
            t = 1.0/t;
            (*(cogmres_functions->ScaleVector))(t,p[i]);
         }
         for (j = 1; j < i; j++)
         {
            t = hh[itmp+j-1];
            hh[itmp+j-1] = s[j-1]*hh[itmp+j] + c[j-1]*t;
            hh[itmp+j] = -s[j-1]*t + c[j-1]*hh[itmp+j];
         }
         t= hh[itmp+i]*hh[itmp+i];
         t+= hh[itmp+i-1]*hh[itmp+i-1];
         gamma = sqrt(t);
         if (gamma == 0.0) gamma = epsmac;
         c[i-1] = hh[itmp+i-1]/gamma;
         s[i-1] = hh[itmp+i]/gamma;
         rs[i] = -hh[itmp+i]*rs[i-1];
         rs[i] /=  gamma;
         rs[i-1] = c[i-1]*rs[i-1];
         // determine residual norm 
         hh[itmp+i-1] = s[i-1]*hh[itmp+i] + c[i-1]*hh[itmp+i-1];
         r_norm = fabs(rs[i]);
         if ( print_level>0 )
         {
            norms[iter] = r_norm;
            if ( print_level>1 && my_id == 0 )
            {
               if (b_norm > 0.0)
                  jxf_printf("% 5d    %e    %f   %e\n", iter, 
                     norms[iter],norms[iter]/norms[iter-1],
                     norms[iter]/b_norm);
               else
                  jxf_printf("% 5d    %e    %f\n", iter, norms[iter],
                     norms[iter]/norms[iter-1]);
            }
         }
         /*convergence factor tolerance */
         if (cf_tol > 0.0)
         {
            cf_ave_0 = cf_ave_1;
            cf_ave_1 = pow( r_norm / r_norm_0, 1.0/(2.0*iter));

            weight = fabs(cf_ave_1 - cf_ave_0);
            weight = weight / jxf_max(cf_ave_1, cf_ave_0);

            weight = 1.0 - weight;
#if 0
           jxf_printf("I = %d: cf_new = %e, cf_old = %e, weight = %e\n",
              i, cf_ave_1, cf_ave_0, weight );
#endif
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
               for (k=0; k<i; k++) /* extra copy of rs so we don't need
                                   to change the later solve */
                  rs_2[k] = rs[k];

               /* solve tri. system*/
               rs_2[i-1] = rs_2[i-1]/hh[itmp+i-1];
               for (k = i-2; k >= 0; k--)
               {
                  t = 0.0;
                  for (j = k+1; j < i; j++)
                  {
                     t -= hh[j*(k_dim+1)+k]*rs_2[j];
                  }
                  t+= rs_2[k];
                  rs_2[k] = t/hh[k*(k_dim+1)+k];
               }
               (*(cogmres_functions->CopyVector))(p[i-1],w);
               (*(cogmres_functions->ScaleVector))(rs_2[i-1],w);
               for (j = i-2; j >=0; j--)
                  (*(cogmres_functions->Axpy))(rs_2[j], p[j], w);

               (*(cogmres_functions->ClearVector))(r);
               /* find correction (in r) */
               precond(precond_data, preOperator, w, r);
               /* copy current solution (x) to w (don't want to over-write x)*/
               (*(cogmres_functions->CopyVector))(x,w);

               /* add the correction */
               (*(cogmres_functions->Axpy))(1.0,r,w);

               /* now w is the approx solution  - get the norm*/
               x_norm = sqrt( (*(cogmres_functions->InnerProd))(w,w) );

               if ( !(x_norm <= guard_zero_residual ))
                  /* don't divide by zero */
               {  /* now get  x_i - x_i-1 */
                  if (num_rel_change_check)
                  {
                     /* have already checked once so we can avoid another precond.
                        solve */
                     (*(cogmres_functions->CopyVector))(w, r);
                     (*(cogmres_functions->Axpy))(-1.0, w_2, r);
                     /* now r contains x_i - x_i-1*/

                     /* save current soln w in w_2 for next time */
                     (*(cogmres_functions->CopyVector))(w, w_2);
                  }
                  else
                  {
                     /* first time to check rel change*/
                     /* first save current soln w in w_2 for next time */
                     (*(cogmres_functions->CopyVector))(w, w_2);

                     (*(cogmres_functions->ClearVector))(w);
                     (*(cogmres_functions->Axpy))(rs_2[i-1], p[i-1], w);
                     (*(cogmres_functions->ClearVector))(r);
                     /* apply the preconditioner */
                     precond(precond_data, preOperator, w, r);
                     /* now r contains x_i - x_i-1 */          
                  }
                  /* find the norm of x_i - x_i-1 */          
                  w_norm = sqrt( (*(cogmres_functions->InnerProd))(r,r) );
                  relative_error = w_norm/x_norm;
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
               num_rel_change_check++;
            }
            else /* no relative change */
            {
               break;
            }
         }
      } /*** end of restart cycle ***/

      /* now compute solution, first solve upper triangular system */
      if (break_value) break;
     
      rs[i-1] = rs[i-1]/hh[itmp+i-1];
      for (k = i-2; k >= 0; k--)
      {
         t = 0.0;
         for (j = k+1; j < i; j++)
         {
            t -= hh[j*(k_dim+1)+k]*rs[j];
         }
         t+= rs[k];
         rs[k] = t/hh[k*(k_dim+1)+k];
      }

      (*(cogmres_functions->CopyVector))(p[i-1],w);
      (*(cogmres_functions->ScaleVector))(rs[i-1],w);
      for (j = i-2; j >=0; j--)
         (*(cogmres_functions->Axpy))(rs[j], p[j], w);

      (*(cogmres_functions->ClearVector))(r);
      /* find correction (in r) */
      precond(precond_data, preOperator, w, r);

      /* update current solution x (in x) */
      (*(cogmres_functions->Axpy))(1.0,r,x);


      /* check for convergence by evaluating the actual residual */
      if (r_norm  <= epsilon && iter >= min_iter)
      {
         if (skip_real_r_check)
         {
            (cogmres_data -> converged) = 1;
            break;
         }

         /* calculate actual residual norm*/
         (*(cogmres_functions->CopyVector))(b,r);
         (*(cogmres_functions->Matvec))(matvec_data,-1.0,A,x,1.0,r);
         real_r_norm_new = r_norm = sqrt( (*(cogmres_functions->InnerProd))(r,r) );

         if (r_norm <= epsilon)
         {
            if (rel_change && !rel_change_passed) /* calculate the relative change */
            {
               /* calculate the norm of the solution */
               x_norm = sqrt( (*(cogmres_functions->InnerProd))(x,x) );

               if ( !(x_norm <= guard_zero_residual ))
               /* don't divide by zero */
               {
                  (*(cogmres_functions->ClearVector))(w);
                  (*(cogmres_functions->Axpy))(rs[i-1], p[i-1], w);
                  (*(cogmres_functions->ClearVector))(r);
                  /* apply the preconditioner */
                  precond(precond_data, preOperator, w, r);
                  /* find the norm of x_i - x_i-1 */          
                  w_norm = sqrt( (*(cogmres_functions->InnerProd))(r,r) );
                  relative_error= w_norm/x_norm;
                  if ( relative_error < r_tol )
                  {
                     (cogmres_data -> converged) = 1;
                     if ( print_level>1 && my_id == 0 )
                     {
                        jxf_printf("\n\n");
                        jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
                     }
                     break;
                  }
               }
               else
               {
                  (cogmres_data -> converged) = 1;
                  if ( print_level>1 && my_id == 0 )
                  {
                     jxf_printf("\n\n");
                     jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
                  }
                  break;
               }
            }
            else /* don't need to check rel. change */
            {
               if ( print_level>1 && my_id == 0 )
               {
                  jxf_printf("\n\n");
                  jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
               }
               (cogmres_data -> converged) = 1;
               break;
            }
         }
         else /* conv. has not occurred, according to true residual */
         {
            /* exit if the real residual norm has not decreased */
            if (real_r_norm_new >= real_r_norm_old)
            {
               if (print_level > 1 && my_id == 0)
               {
                  jxf_printf("\n\n");
                  jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
               }
               (cogmres_data -> converged) = 1;
               break;
            }
            /* report discrepancy between real/COGMRES residuals and restart */
            if ( print_level>0 && my_id == 0)
               jxf_printf("false convergence 2, L2 norm of residual: %e\n", r_norm);
            if (is_check_restarted == 1)
            {
               (*(cogmres_functions->CopyVector))(r,p[0]);
               i = 0;
               real_r_norm_old = real_r_norm_new;
            }
            else
            {
               break;
            }
         }
      } /* end of convergence check */

      /* compute residual vector and continue loop */
      for (j=i ; j > 0; j--)
      {
         rs[j-1] = -s[j-1]*rs[j];
         rs[j] = c[j-1]*rs[j];
      }

      if (i) (*(cogmres_functions->Axpy))(rs[i]-1.0,p[i],p[i]);
      for (j=i-1 ; j > 0; j--)
         (*(cogmres_functions->Axpy))(rs[j],p[j],p[i]);

      if (i)
      {
         (*(cogmres_functions->Axpy))(rs[0]-1.0,p[0],p[0]);
         (*(cogmres_functions->Axpy))(1.0,p[i],p[0]);
      }

   } /* END of iteration while loop */


   (cogmres_data -> num_iterations) = iter;
   if (b_norm > 0.0)
      (cogmres_data -> rel_residual_norm) = r_norm/b_norm;
   if (b_norm == 0.0)
      (cogmres_data -> rel_residual_norm) = r_norm;

   if (iter >= max_iter && r_norm > epsilon && epsilon > 0) jxf_error(JXF_ERROR_CONV);

   jxf_TFreeF(c,cogmres_functions); 
   jxf_TFreeF(s,cogmres_functions); 
   jxf_TFreeF(rs,cogmres_functions);
   jxf_TFreeF(rv,cogmres_functions);
   if (rel_change)  jxf_TFreeF(rs_2,cogmres_functions);

   /*for (i=0; i < k_dim+1; i++)
   {  
      jxf_TFreeF(hh[i],cogmres_functions);
      jxf_TFreeF(uu[i],cogmres_functions);
   }*/
   jxf_TFreeF(hh,cogmres_functions); 
   jxf_TFreeF(uu,cogmres_functions);

   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetKDim( void   *cogmres_vdata, JXF_Int   k_dim )
{
   jxf_COGMRESData *cogmres_data =(jxf_COGMRESData *) cogmres_vdata;
   (cogmres_data -> k_dim) = k_dim;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetKDim( void   *cogmres_vdata, JXF_Int * k_dim )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *k_dim = (cogmres_data -> k_dim);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetUnroll( void   *cogmres_vdata, JXF_Int   unroll )
{
   jxf_COGMRESData *cogmres_data =(jxf_COGMRESData *) cogmres_vdata;
   (cogmres_data -> unroll) = unroll;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetUnroll( void   *cogmres_vdata, JXF_Int * unroll )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *unroll = (cogmres_data -> unroll);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetCGS( void   *cogmres_vdata, JXF_Int   cgs )
{
   jxf_COGMRESData *cogmres_data =(jxf_COGMRESData *) cogmres_vdata;
   (cogmres_data -> cgs) = cgs;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetCGS( void   *cogmres_vdata, JXF_Int * cgs )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *cgs = (cogmres_data -> cgs);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetTol( void   *cogmres_vdata, JXF_Real  tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> tol) = tol;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetTol( void   *cogmres_vdata, JXF_Real  * tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *tol = (cogmres_data -> tol);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetAbsoluteTol( void   *cogmres_vdata, JXF_Real  a_tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> a_tol) = a_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetAbsoluteTol( void   *cogmres_vdata, JXF_Real  * a_tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *a_tol = (cogmres_data -> a_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetConvergenceFactorTol( void   *cogmres_vdata, JXF_Real  cf_tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> cf_tol) = cf_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetConvergenceFactorTol( void   *cogmres_vdata, JXF_Real * cf_tol )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *cf_tol = (cogmres_data -> cf_tol);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetMinIter( void *cogmres_vdata, JXF_Int   min_iter )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> min_iter) = min_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetMinIter( void *cogmres_vdata, JXF_Int * min_iter )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *min_iter = (cogmres_data -> min_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetMaxIter( void *cogmres_vdata, JXF_Int   max_iter )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetIsCheckRestarted( void *cogmres_vdata, JXF_Int is_check_restarted )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> is_check_restarted) = is_check_restarted;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetMaxIter( void *cogmres_vdata, JXF_Int * max_iter )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *max_iter = (cogmres_data -> max_iter);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetRelChange( void *cogmres_vdata, JXF_Int   rel_change )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> rel_change) = rel_change;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetRelChange( void *cogmres_vdata, JXF_Int * rel_change )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *rel_change = (cogmres_data -> rel_change);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetSkipRealResidualCheck( void *cogmres_vdata, JXF_Int skip_real_r_check )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> skip_real_r_check) = skip_real_r_check;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetSkipRealResidualCheck( void *cogmres_vdata, JXF_Int *skip_real_r_check )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *skip_real_r_check = (cogmres_data -> skip_real_r_check);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetPrecond( void  *cogmres_vdata,
        JXF_Int  (*precond)(void*,void*,void*,void*),
        JXF_Int  (*precond_setup)(void*,void*,void*,void*),
        void  *precond_data )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   jxf_COGMRESFunctions *cogmres_functions = cogmres_data->functions;
   (cogmres_functions -> precond)        = precond;
   (cogmres_functions -> precond_setup)  = precond_setup;
   (cogmres_data -> precond_data)   = precond_data;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetPrecond( void *cogmres_vdata, JXF_Solver *precond_data_ptr )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *precond_data_ptr = (JXF_Solver)(cogmres_data -> precond_data);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetPrintLevel( void *cogmres_vdata, JXF_Int   level )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> print_level) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetPrintLevel( void *cogmres_vdata, JXF_Int * level )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *level = (cogmres_data -> print_level);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESSetLogging( void *cogmres_vdata, JXF_Int   level )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   (cogmres_data -> logging) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetLogging( void *cogmres_vdata, JXF_Int * level )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *level = (cogmres_data -> logging);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetNumIterations( void *cogmres_vdata, JXF_Int  *num_iterations )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *num_iterations = (cogmres_data -> num_iterations);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetConverged( void *cogmres_vdata, JXF_Int  *converged )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *converged = (cogmres_data -> converged);
   return jxf_error_flag;
}

JXF_Int
jxf_COGMRESGetFinalRelativeResidualNorm( void   *cogmres_vdata, JXF_Real *relative_residual_norm )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   *relative_residual_norm = (cogmres_data -> rel_residual_norm);
   return jxf_error_flag;
}


JXF_Int
jxf_COGMRESSetModifyPC( void *cogmres_vdata, JXF_Int (*modify_pc)(void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm) )
{
   jxf_COGMRESData *cogmres_data = (jxf_COGMRESData *)cogmres_vdata;
   jxf_COGMRESFunctions *cogmres_functions = cogmres_data->functions;
   (cogmres_functions -> modify_pc)        = modify_pc;
   return jxf_error_flag;
} 
