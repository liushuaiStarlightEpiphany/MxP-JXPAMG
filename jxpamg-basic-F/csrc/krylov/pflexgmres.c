//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  pflexgmres.c
 *  Date: 2019/12/05
 */ 

#include "jxf_pamg.h"
#include "jxf_krylov.h"

JXF_Int
JXF_ParCSRFlexGMRESCreate( MPI_Comm comm, JXF_Solver *solver )
{
   jxf_FlexGMRESFunctions *fgmres_functions;

   if (!solver)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   fgmres_functions =
      jxf_FlexGMRESFunctionsCreate(
         jxf_CAlloc, jxf_ParKrylovFree, jxf_ParKrylovCommInfo,
         jxf_ParKrylovCreateVector,
         jxf_ParKrylovCreateVectorArray,
         jxf_ParKrylovDestroyVector, jxf_ParKrylovMatvecCreate,
         jxf_ParKrylovMatvec, jxf_ParKrylovMatvecDestroy,
         jxf_ParKrylovInnerProd, jxf_ParKrylovCopyVector,
         jxf_ParKrylovClearVector,
         jxf_ParKrylovScaleVector, jxf_ParKrylovAxpy,
         jxf_ParKrylovIdentitySetup, jxf_ParKrylovIdentity );

  *solver = ( (JXF_Solver) jxf_FlexGMRESCreate( fgmres_functions ) );

   return jxf_error_flag;
}

JXF_Int 
JXF_ParCSRFlexGMRESDestroy( JXF_Solver solver )
{
   return( jxf_FlexGMRESDestroy( (void *) solver ) );
}

JXF_Int 
JXF_FlexGMRESSetup( JXF_Solver solver, JXF_Matrix A, JXF_Vector b, JXF_Vector x )
{
   return( jxf_FlexGMRESSetup( solver, A, b, x ) );
}

JXF_Int 
JXF_FlexGMRESSolve( JXF_Solver solver, JXF_Matrix preOperator, JXF_Matrix A, JXF_Vector b, JXF_Vector x )
{
   return( jxf_FlexGMRESSolve( solver, preOperator, A, b, x ) );
}

JXF_Int
JXF_FlexGMRESSetKDim( JXF_Solver solver, JXF_Int k_dim )
{
   return( jxf_FlexGMRESSetKDim( (void *) solver, k_dim ) );
}

JXF_Int
JXF_FlexGMRESSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_FlexGMRESSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_FlexGMRESSetAbsoluteTol( JXF_Solver solver, JXF_Real a_tol )
{
   return( jxf_FlexGMRESSetAbsoluteTol( (void *) solver, a_tol ) );
}

JXF_Int
JXF_FlexGMRESSetConvergenceFactorTol( JXF_Solver solver, JXF_Real cf_tol )
{
   return( jxf_FlexGMRESSetConvergenceFactorTol( (void *) solver, cf_tol ) );
}

JXF_Int
JXF_FlexGMRESSetMinIter( JXF_Solver solver, JXF_Int min_iter )
{
   return( jxf_FlexGMRESSetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_FlexGMRESSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_FlexGMRESSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_FlexGMRESSetPrecond( JXF_Solver solver,
                        JXF_PtrToSolverFcn precond,
                        JXF_PtrToSolverFcn precond_setup,
                        JXF_Solver precond_solver )
{
   return( jxf_FlexGMRESSetPrecond( (void *) solver,
                                   (JXF_Int (*)(void*, void*, void*, void*))precond,
                                   (JXF_Int (*)(void*, void*))precond_setup,
                                   (void *) precond_solver ) );
}

JXF_Int
JXF_FlexGMRESSetPrintLevel( JXF_Solver solver, JXF_Int level )
{
   return( jxf_FlexGMRESSetPrintLevel( (void *) solver, level ) );
}

JXF_Int
JXF_FlexGMRESSetIsCheckRestarted( JXF_Solver solver, JXF_Int is_check_restarted )
{
   return( jxf_FlexGMRESSetIsCheckRestarted( (void *) solver, is_check_restarted ) );
}

JXF_Int
JXF_FlexGMRESSetLogging( JXF_Solver solver, JXF_Int level )
{
   return( jxf_FlexGMRESSetLogging( (void *) solver, level ) );
}

JXF_Int
JXF_FlexGMRESGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations )
{
   return( jxf_FlexGMRESGetNumIterations( (void *) solver, num_iterations ) );
}

JXF_Int
JXF_FlexGMRESGetConverged( JXF_Solver solver, JXF_Int *converged )
{
   return( jxf_FlexGMRESGetConverged( (void *) solver, converged ) );
}

JXF_Int
JXF_FlexGMRESGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *norm )
{
   return( jxf_FlexGMRESGetFinalRelativeResidualNorm( (void *) solver, norm ) );
}

JXF_Int
JXF_FlexGMRESGetResidual( JXF_Solver solver, void *residual )
{
   /* returns a pointer to the residual vector */
   return jxf_FlexGMRESGetResidual( (void *) solver, (void **) residual );
}

JXF_Int
JXF_FlexGMRESSetModifyPC( JXF_Solver  solver, JXF_Int (*modify_pc)(JXF_Solver, JXF_Int, JXF_Real) )
{
   return jxf_FlexGMRESSetModifyPC( (void *) solver, (JXF_Int(*)(void*, JXF_Int, JXF_Real))modify_pc );
}

/*!
 * \fn jxf_FlexGMRESFunctions *jxf_FlexGMRESFunctionsCreate
 */
jxf_FlexGMRESFunctions *
jxf_FlexGMRESFunctionsCreate(
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
   JXF_Int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
)
{
   jxf_FlexGMRESFunctions *fgmres_functions;
   fgmres_functions = (jxf_FlexGMRESFunctions *)CAlloc( 1, sizeof(jxf_FlexGMRESFunctions) );

   fgmres_functions->CAlloc = CAlloc;
   fgmres_functions->Free = Free;
   fgmres_functions->CommInfo = CommInfo; /* not in PCGFunctionsCreate */
   fgmres_functions->CreateVector = CreateVector;
   fgmres_functions->CreateVectorArray = CreateVectorArray; /* not in PCGFunctionsCreate */
   fgmres_functions->DestroyVector = DestroyVector;
   fgmres_functions->MatvecCreate = MatvecCreate;
   fgmres_functions->Matvec = Matvec;
   fgmres_functions->MatvecDestroy = MatvecDestroy;
   fgmres_functions->InnerProd = InnerProd;
   fgmres_functions->CopyVector = CopyVector;
   fgmres_functions->ClearVector = ClearVector;
   fgmres_functions->ScaleVector = ScaleVector;
   fgmres_functions->Axpy = Axpy;
/* default preconditioner must be set here but can be changed later... */
   fgmres_functions->precond_setup = PrecondSetup;
   fgmres_functions->precond       = Precond;

   fgmres_functions->modify_pc     = jxf_FlexGMRESModifyPCDefault;

   return fgmres_functions;
}

/*!
 * \fn void *jxf_FlexGMRESCreate
 */
void *
jxf_FlexGMRESCreate( jxf_FlexGMRESFunctions *fgmres_functions )
{
   jxf_FlexGMRESData *fgmres_data;
 
   fgmres_data = jxf_CTAllocF(jxf_FlexGMRESData, 1, fgmres_functions);

   fgmres_data->functions = fgmres_functions;
 
   /* set defaults */
   (fgmres_data -> k_dim)              = 20;
   (fgmres_data -> tol)                = 1.0e-06;
   (fgmres_data -> cf_tol)             = 0.0;
   (fgmres_data -> a_tol)              = 0.0; /* abs. residual tol */
   (fgmres_data -> min_iter)           = 0;
   (fgmres_data -> max_iter)           = 1000;
   (fgmres_data -> rel_change)         = 0;
   (fgmres_data -> stop_crit)          = 0; /* rel. residual norm */
   (fgmres_data -> converged)          = 0;
   (fgmres_data -> is_check_restarted) = 0;
   (fgmres_data -> precond_data)       = NULL;
   (fgmres_data -> print_level)        = 0;
   (fgmres_data -> logging)            = 0;
   (fgmres_data -> p)                  = NULL;
   (fgmres_data -> r)                  = NULL;
   (fgmres_data -> w)                  = NULL;
   (fgmres_data -> w_2)                = NULL;
   (fgmres_data -> matvec_data)        = NULL;
   (fgmres_data -> norms)              = NULL;
   (fgmres_data -> log_file_name)      = NULL;
 
   return (void *) fgmres_data;
}

/*!
 * \fn JXF_Int jxf_FlexGMRESDestroy
 */
 
JXF_Int
jxf_FlexGMRESDestroy( void *fgmres_vdata )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   JXF_Int i;
 
   if (fgmres_data)
   {
      jxf_FlexGMRESFunctions *fgmres_functions = fgmres_data->functions;
      if ( (fgmres_data->logging>0) || (fgmres_data->print_level) > 0 )
      {
         if ( (fgmres_data -> norms) != NULL ) jxf_TFreeF( fgmres_data -> norms, fgmres_functions );
      }
 
      if ( (fgmres_data -> matvec_data) != NULL ) (*(fgmres_functions->MatvecDestroy))(fgmres_data -> matvec_data);
 
      if ( (fgmres_data -> r) != NULL ) (*(fgmres_functions->DestroyVector))(fgmres_data -> r);
      if ( (fgmres_data -> w) != NULL ) (*(fgmres_functions->DestroyVector))(fgmres_data -> w);
      if ( (fgmres_data -> w_2) != NULL ) (*(fgmres_functions->DestroyVector))(fgmres_data -> w_2);

      if ( (fgmres_data -> p) != NULL )
      {
         for (i = 0; i < (fgmres_data -> k_dim+1); i++)
         {
            if ( (fgmres_data -> p)[i] != NULL ) (*(fgmres_functions->DestroyVector))( (fgmres_data -> p) [i]);
         }
         jxf_TFreeF( fgmres_data->p, fgmres_functions );
      }

      /* fgmres mod  - space for precond. vectors*/
      if ( (fgmres_data -> pre_vecs) != NULL )
      {
         for (i = 0; i < (fgmres_data -> k_dim + 1); i++)
         {
            if ( (fgmres_data -> pre_vecs)[i] != NULL ) (*(fgmres_functions->DestroyVector))( (fgmres_data -> pre_vecs) [i]);
         }
         jxf_TFreeF( fgmres_data->pre_vecs, fgmres_functions );
      }
      /*---*/

      jxf_TFreeF( fgmres_data, fgmres_functions );
      jxf_TFreeF( fgmres_functions, fgmres_functions );
   }
 
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_FlexGMRESGetResidual
 */
JXF_Int
jxf_FlexGMRESGetResidual( void *fgmres_vdata, void **residual )
{
   /* returns a pointer to the residual vector */
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
  *residual = fgmres_data->r;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_FlexGMRESSetup
 */
JXF_Int
jxf_FlexGMRESSetup( void *fgmres_vdata, void *A, void *b, void *x )
{
   jxf_FlexGMRESData *fgmres_data     = (jxf_FlexGMRESData *)fgmres_vdata;
   jxf_FlexGMRESFunctions *fgmres_functions = fgmres_data->functions;

   JXF_Int            k_dim            = (fgmres_data -> k_dim);
   JXF_Int            max_iter         = (fgmres_data -> max_iter);
   //JXF_Int          (*precond_setup)(void*,void*) = (fgmres_functions->precond_setup);
   //void          *precond_data     = (fgmres_data -> precond_data);

   JXF_Int            rel_change       = (fgmres_data -> rel_change);
 
   (fgmres_data -> A) = A;
 
   /*--------------------------------------------------
    * The arguments for NewVector are important to
    * maintain consistency between the setup and
    * compute phases of matvec and the preconditioner.
    *--------------------------------------------------*/
 
   if ((fgmres_data -> p) == NULL) (fgmres_data -> p) = (void**)(*(fgmres_functions->CreateVectorArray))(k_dim+1,x);
   if ((fgmres_data -> r) == NULL) (fgmres_data -> r) = (*(fgmres_functions->CreateVector))(b);
   if ((fgmres_data -> w) == NULL) (fgmres_data -> w) = (*(fgmres_functions->CreateVector))(b);
 
   if (rel_change)
   {  
      if ((fgmres_data -> w_2) == NULL) (fgmres_data -> w_2) = (*(fgmres_functions->CreateVector))(b);
   }
 
   /* fgmres mod */
   (fgmres_data -> pre_vecs) = (void**)(*(fgmres_functions->CreateVectorArray))(k_dim+1,x); 
   /*---*/

   if ((fgmres_data -> matvec_data) == NULL) (fgmres_data -> matvec_data) = (*(fgmres_functions->MatvecCreate))(A, x);
 
   /* precond_setup(precond_data, A); */
 
   /*-----------------------------------------------------
    * Allocate space for log info
    *-----------------------------------------------------*/
 
   if ( (fgmres_data->logging)>0 || (fgmres_data->print_level) > 0 )
   {
      if ((fgmres_data -> norms) == NULL)
         (fgmres_data -> norms) = jxf_CTAllocF(JXF_Real, max_iter + 1,fgmres_functions);
   }

   if ( (fgmres_data->print_level) > 0 )
   {
      if ((fgmres_data -> log_file_name) == NULL) (fgmres_data -> log_file_name) = (char*)"gmres.out.log";
   }
 
   return jxf_error_flag;
}
 
/*!
 * \fn JXF_Int jxf_FlexGMRESSolve
 */
JXF_Int
jxf_FlexGMRESSolve( void *fgmres_vdata, void *preOperator, void *A, void *b, void *x )
{
   jxf_FlexGMRESData  *fgmres_data   = (jxf_FlexGMRESData *)fgmres_vdata;
   jxf_FlexGMRESFunctions *fgmres_functions = fgmres_data->functions;
   JXF_Int 		     k_dim        = (fgmres_data -> k_dim);
   JXF_Int               min_iter     = (fgmres_data -> min_iter);
   JXF_Int 		     max_iter     = (fgmres_data -> max_iter);
   JXF_Real 	     r_tol        = (fgmres_data -> tol);
   JXF_Real 	     cf_tol       = (fgmres_data -> cf_tol);
   JXF_Real        a_tol        = (fgmres_data -> a_tol);
   void             *matvec_data  = (fgmres_data -> matvec_data);

   void             *r            = (fgmres_data -> r);
   void             *w            = (fgmres_data -> w);

   void            **p            = (fgmres_data -> p);

   /* fgmres  mod*/
   void          **pre_vecs       = (fgmres_data ->pre_vecs);
   /*---*/

   JXF_Int 	           (*precond)(void*,void*,void*,void*)   = (fgmres_functions -> precond);
   JXF_Int 	            *precond_data = (JXF_Int*)(fgmres_data -> precond_data);

   JXF_Int             print_level    = (fgmres_data -> print_level);
   JXF_Int             logging        = (fgmres_data -> logging);

   JXF_Int             is_check_restarted = (fgmres_data -> is_check_restarted);

   JXF_Real     *norms          = (fgmres_data -> norms);
   
   JXF_Int        break_value = 0;
   JXF_Int	      i, j, k;
   JXF_Real *rs, **hh, *c, *s; 
   JXF_Int        iter; 
   JXF_Int        my_id, num_procs;
   JXF_Real epsilon, gamma, t, r_norm, b_norm, den_norm;

   JXF_Real epsmac = 1.e-16; 
   JXF_Real ieee_check = 0.;

   JXF_Real cf_ave_0 = 0.0;
   JXF_Real cf_ave_1 = 0.0;
   JXF_Real weight;
   JXF_Real r_norm_0;

   JXF_Int 	      (*modify_pc)(void*,JXF_Int,JXF_Real)   = (fgmres_functions -> modify_pc);

   /* We are not checking rel. change for now... */

   (fgmres_data -> converged) = 0;

   /*-----------------------------------------------------------------------
    * With relative change convergence test on, it is possible to attempt
    * another iteration with a zero residual. This causes the parameter
    * alpha to go NaN. The guard_zero_residual parameter is to circumvent
    * this. Perhaps it should be set to something non-zero (but small).
    *-----------------------------------------------------------------------*/

   (*(fgmres_functions->CommInfo))(A,&my_id,&num_procs);
   if ( logging>0 || print_level>0 )
   {
      norms          = (fgmres_data -> norms);
      /* not used yet      log_file_name  = (fgmres_data -> log_file_name);*/
      /* fp = fopen(log_file_name,"w"); */
   }

   /* initialize work arrays  */
   rs = jxf_CTAllocF(JXF_Real,k_dim+1,fgmres_functions);
   c = jxf_CTAllocF(JXF_Real,k_dim,fgmres_functions);
   s = jxf_CTAllocF(JXF_Real,k_dim,fgmres_functions);

   /* fgmres mod. - need non-modified hessenberg ???? */
   hh = jxf_CTAllocF(JXF_Real*,k_dim+1,fgmres_functions); 
   for (i=0; i < k_dim+1; i++)
   {
      hh[i] = jxf_CTAllocF(JXF_Real,k_dim,fgmres_functions); 
   }

   (*(fgmres_functions->CopyVector))(b,p[0]);

   /* compute initial residual */
   (*(fgmres_functions->Matvec))(matvec_data,-1.0, A, x, 1.0, p[0]);

   b_norm = sqrt((*(fgmres_functions->InnerProd))(b,b));

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
        jxf_printf("\n\nERROR detected by Hypre ... BEGIN\n");
        jxf_printf("ERROR -- jxf_FlexGMRESSolve: INFs and/or NaNs detected in input.\n");
        jxf_printf("User probably placed non-numerics in supplied b.\n");
        jxf_printf("Returning error flag += 101.  Program not terminated.\n");
        jxf_printf("ERROR detected by Hypre ... END\n\n\n");
      }
      jxf_error(JXF_ERROR_GENERIC);
      return jxf_error_flag;
   }

   r_norm = sqrt((*(fgmres_functions->InnerProd))(p[0],p[0]));
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
        jxf_printf("\n\nERROR detected by Hypre ... BEGIN\n");
        jxf_printf("ERROR -- jxf_FlexGMRESSolve: INFs and/or NaNs detected in input.\n");
        jxf_printf("User probably placed non-numerics in supplied A or x_0.\n");
        jxf_printf("Returning error flag += 101.  Program not terminated.\n");
        jxf_printf("ERROR detected by Hypre ... END\n\n\n");
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
         if (b_norm == 0.0) jxf_printf("Rel_resid_norm actually contains the residual norm\n");
         jxf_printf("Initial L2 norm of residual: %e\n", r_norm);
      }
   }
   iter = 0;

   if (b_norm > 0.0)
   {
/* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
     den_norm= b_norm;
   }
   else
   {
/* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
     den_norm= r_norm;
   }

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
      }
   }

   /* outer iteration cycle */
   while (iter < max_iter)
   {
   /* initialize first term of hessenberg system */

      rs[0] = r_norm;
      if (r_norm == 0.0)
      {
         jxf_TFreeF(c,fgmres_functions); 
         jxf_TFreeF(s,fgmres_functions); 
         jxf_TFreeF(rs,fgmres_functions);

         for (i=0; i < k_dim+1; i++)
         {
            jxf_TFreeF(hh[i],fgmres_functions);
         }

         jxf_TFreeF(hh,fgmres_functions); 
	     return jxf_error_flag;
      }

        /* see if we are already converged and should print the final norm and exit */
      if (r_norm  <= epsilon && iter >= min_iter) 
      {
         (*(fgmres_functions->CopyVector))(b,r);
         (*(fgmres_functions->Matvec))(matvec_data,-1.0,A,x,1.0,r);
         r_norm = sqrt((*(fgmres_functions->InnerProd))(r,r));
         if (r_norm <= epsilon)
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

      t = 1.0 / r_norm;

      (*(fgmres_functions->ScaleVector))(t,p[0]);
      i = 0;

      /***RESTART CYCLE (right-preconditioning) ***/
      while (i < k_dim  && iter < max_iter)
      {
         i++;
         iter++;

         (*(fgmres_functions->ClearVector))(pre_vecs[i-1]);

         /* allow some user function here (to change
          * prec. attributes, i.e.tolerances, etc. ? */
         modify_pc(precond_data, iter, r_norm/den_norm );

         /*apply preconditioner and store in pre_vecs */
         precond(precond_data, preOperator, p[i-1], pre_vecs[i-1]);
         /*apply operator and store in p */
         (*(fgmres_functions->Matvec))(matvec_data, 1.0, A, pre_vecs[i-1], 0.0, p[i]);

         /* modified Gram_Schmidt */
         for (j=0; j < i; j++)
         {
            hh[j][i-1] = (*(fgmres_functions->InnerProd))(p[j],p[i]);
            (*(fgmres_functions->Axpy))(-hh[j][i-1],p[j],p[i]);
         }
         t = sqrt((*(fgmres_functions->InnerProd))(p[i],p[i]));
         hh[i][i-1] = t;	
         if (t != 0.0)
         {
            t = 1.0/t;
            (*(fgmres_functions->ScaleVector))(t,p[i]);
         }

         /* done with modified Gram_schmidt and Arnoldi step.
            update factorization of hh */
         for (j = 1; j < i; j++)
         {
            t = hh[j-1][i-1];
            hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
            hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
         }
         t= hh[i][i-1]*hh[i][i-1];
         t+= hh[i-1][i-1]*hh[i-1][i-1];
         gamma = sqrt(t);
         if (gamma == 0.0) gamma = epsmac;
         c[i-1] = hh[i-1][i-1]/gamma;
         s[i-1] = hh[i][i-1]/gamma;
         rs[i] = -hh[i][i-1]*rs[i-1];
         rs[i]/=  gamma;
         rs[i-1] = c[i-1]*rs[i-1];
         /* determine residual norm */
         hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
         r_norm = fabs(rs[i]);

         /* print ? */
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
         if (r_norm  <= epsilon && iter >= min_iter)
         {
            /* no relative change */
            break;
         }
      } /*** end of restart cycle ***/

      /* now compute solution, first solve upper triangular system */

      if (break_value) break;

      rs[i-1] = rs[i-1]/hh[i-1][i-1];
      for (k = i-2; k >= 0; k--)
      {
         t = 0.0;
         for (j = k+1; j < i; j++)
         {
            t -= hh[k][j]*rs[j];
         }
         t+= rs[k];
         rs[k] = t/hh[k][k];
      }

      /* form linear combination of pre_vecs's to get solution */

      (*(fgmres_functions->CopyVector))(pre_vecs[i-1],w);
      (*(fgmres_functions->ScaleVector))(rs[i-1],w);
      for (j = i-2; j >=0; j--) (*(fgmres_functions->Axpy))(rs[j], pre_vecs[j], w);

      /* don't need to un-wind precond... - so now the correction is
       * in w */

      /* update current solution x (in x) */
      (*(fgmres_functions->Axpy))(1.0,w,x);

      /* check for convergence by evaluating the actual residual */
      if (r_norm <= epsilon && iter >= min_iter) 
      {
         /* calculate actual residual norm*/
         (*(fgmres_functions->CopyVector))(b,r);
         (*(fgmres_functions->Matvec))(matvec_data,-1.0,A,x,1.0,r);
         r_norm = sqrt( (*(fgmres_functions->InnerProd))(r,r) );

         if (r_norm <= epsilon)
         {
            if ( print_level>1 && my_id == 0 )
            {
               jxf_printf("\n\n");
               jxf_printf("Final L2 norm of residual: %e\n\n", r_norm);
            }
            (fgmres_data -> converged) = 1;
            break;
         }
         else /* conv. has not occurred, according to true residual */ 
         {
            if ( print_level>0 && my_id == 0) jxf_printf("false convergence 2\n");
            if (is_check_restarted == 1)
            {
               (*(fgmres_functions->CopyVector))(r,p[0]);
               i = 0;
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

      if (i) (*(fgmres_functions->Axpy))(rs[i]-1.0,p[i],p[i]);
      for (j=i-1 ; j > 0; j--) (*(fgmres_functions->Axpy))(rs[j],p[j],p[i]);

      if (i)
      {
         (*(fgmres_functions->Axpy))(rs[0]-1.0,p[0],p[0]);
         (*(fgmres_functions->Axpy))(1.0,p[i],p[0]);
      }
   } /* END of iteration while loop */

   if ( print_level>1 && my_id == 0 ) jxf_printf("\n\n"); 

   (fgmres_data -> num_iterations) = iter;
   if (b_norm > 0.0)
      (fgmres_data -> rel_residual_norm) = r_norm/b_norm;
   if (b_norm == 0.0)
      (fgmres_data -> rel_residual_norm) = r_norm;

   if (iter >= max_iter && r_norm > epsilon && epsilon > 0) jxf_error(JXF_ERROR_CONV);

   jxf_TFreeF(c,fgmres_functions);
   jxf_TFreeF(s,fgmres_functions);
   jxf_TFreeF(rs,fgmres_functions);

   for (i=0; i < k_dim+1; i++)
   {
   	  jxf_TFreeF(hh[i],fgmres_functions);
   }
   jxf_TFreeF(hh,fgmres_functions);

   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetKDim( void *fgmres_vdata, JXF_Int k_dim )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> k_dim) = k_dim;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetTol( void *fgmres_vdata, JXF_Real tol )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> tol) = tol;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetAbsoluteTol( void *fgmres_vdata, JXF_Real a_tol )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> a_tol) = a_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetConvergenceFactorTol( void *fgmres_vdata, JXF_Real cf_tol )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> cf_tol) = cf_tol;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetMinIter( void *fgmres_vdata, JXF_Int min_iter )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> min_iter) = min_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetMaxIter( void *fgmres_vdata, JXF_Int max_iter )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetStopCrit( void *fgmres_vdata, JXF_Int stop_crit )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> stop_crit) = stop_crit;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetPrecond( void *fgmres_vdata,
                        JXF_Int (*precond)(void*,void*,void*,void*),
                        JXF_Int (*precond_setup)(void*,void*),
                        void *precond_data )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   jxf_FlexGMRESFunctions *fgmres_functions = fgmres_data->functions;
   (fgmres_functions -> precond) = precond;
   (fgmres_functions -> precond_setup) = precond_setup;
   (fgmres_data -> precond_data) = precond_data;
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESSetPrintLevel( void *fgmres_vdata, JXF_Int level )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> print_level) = level;
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESSetIsCheckRestarted( void *fgmres_vdata, JXF_Int is_check_restarted )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> is_check_restarted) = is_check_restarted;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESSetLogging( void *fgmres_vdata, JXF_Int level )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   (fgmres_data -> logging) = level;
   return jxf_error_flag;
}

JXF_Int
jxf_FlexGMRESGetNumIterations( void *fgmres_vdata, JXF_Int *num_iterations )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
  *num_iterations = (fgmres_data -> num_iterations);
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESGetConverged( void *fgmres_vdata, JXF_Int *converged )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
  *converged = (fgmres_data -> converged);
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESGetFinalRelativeResidualNorm( void *fgmres_vdata, JXF_Real *relative_residual_norm )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
  *relative_residual_norm = (fgmres_data -> rel_residual_norm);
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESSetModifyPC( void *fgmres_vdata, 
                         JXF_Int (*modify_pc)(void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm) )
{
   jxf_FlexGMRESData *fgmres_data = (jxf_FlexGMRESData *)fgmres_vdata;
   jxf_FlexGMRESFunctions *fgmres_functions = fgmres_data->functions;
   (fgmres_functions -> modify_pc) = modify_pc;
   return jxf_error_flag;
}
 
JXF_Int
jxf_FlexGMRESModifyPCDefault( void *precond_data, JXF_Int iteration, JXF_Real rel_residual_norm )
{
   /* Here would could check the number of its and the current
      residual and make some changes to the preconditioner. */
   return 0;
}
