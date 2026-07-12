//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_ilu.c
 *
 */

#include "jxf_ilu.h"
#include "jxf_krylov.h"

JXF_Int
JXF_ILUCreate( JXF_Solver *solver )
{
   if (!solver)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   *solver = ( (JXF_Solver) jxf_ILUCreate( ) );
   return jxf_error_flag;
}

JXF_Int
JXF_ILUDestroy( JXF_Solver solver )
{
   return( jxf_ILUDestroy( (void *) solver ) );
}

JXF_Int
JXF_ILUSetup( JXF_Solver solver, JXF_ParCSRMatrix A )
{
   return( jxf_ILUSetup( (void *) solver, (jxf_ParCSRMatrix *) A ));
}

JXF_Int
JXF_ILUSolve( JXF_Solver solver, JXF_ParCSRMatrix A, JXF_ParVector b, JXF_ParVector x )
{
   return( jxf_ILUSolve( (void *) solver,
                                 (jxf_ParCSRMatrix *) A,
                                 (jxf_ParVector *) b,
                                 (jxf_ParVector *) x ) );
}

JXF_Int
JXF_ILUSetPrintLevel( JXF_Solver solver, JXF_Int print_level )
{
   return jxf_ILUSetPrintLevel( solver, print_level );
}

JXF_Int
JXF_ILUSetLogging( JXF_Solver solver, JXF_Int logging )
{
   return jxf_ILUSetLogging(solver, logging );
}

JXF_Int
JXF_ILUSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return jxf_ILUSetMaxIter( solver, max_iter );
}

JXF_Int
JXF_ILUSetTol( JXF_Solver solver, JXF_Real tol )
{
   return jxf_ILUSetTol( solver, tol );
}

JXF_Int
JXF_ILUSetDropThreshold( JXF_Solver solver, JXF_Real threshold )
{
   return jxf_ILUSetDropThreshold( solver, threshold );
}

JXF_Int
JXF_ILUSetDropThresholdArray( JXF_Solver solver, JXF_Real *threshold )
{
   return jxf_ILUSetDropThresholdArray( solver, threshold );
}

JXF_Int
JXF_ILUSetNSHDropThreshold( JXF_Solver solver, JXF_Real threshold )
{
   return jxf_ILUSetSchurNSHDropThreshold( solver, threshold );
}

JXF_Int
JXF_ILUSetNSHDropThresholdArray( JXF_Solver solver, JXF_Real *threshold )
{
   return jxf_ILUSetSchurNSHDropThresholdArray( solver, threshold );
}

JXF_Int
JXF_ILUSetSchurMaxIter( JXF_Solver solver, JXF_Int ss_max_iter )
{
   return jxf_ILUSetSchurSolverMaxIter( solver, ss_max_iter );
}

JXF_Int
JXF_ILUSetMaxNnzPerRow( JXF_Solver solver, JXF_Int nzmax )
{
   return jxf_ILUSetMaxNnzPerRow( solver, nzmax );
}

JXF_Int
JXF_ILUSetLevelOfFill( JXF_Solver solver, JXF_Int lfil )
{
   return jxf_ILUSetLevelOfFill( solver, lfil );
}

JXF_Int
JXF_ILUSetType( JXF_Solver solver, JXF_Int ilu_type )
{
   return jxf_ILUSetType( solver, ilu_type );
}

JXF_Int
JXF_ILUGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations )
{
   return jxf_ILUGetNumIterations( solver, num_iterations );
}

JXF_Int
JXF_ILUGetFinalRelativeResidualNorm(  JXF_Solver solver, JXF_Real *res_norm )
{
   return jxf_ILUGetFinalRelativeResidualNorm(solver, res_norm);
}

JXF_Int
JXF_ILUSetLocalReordering(  JXF_Solver solver, JXF_Int ordering_type )
{
   return jxf_ILUSetLocalReordering(solver, ordering_type);
}

/* Create */
void *
jxf_ILUCreate()
{
   jxf_ParILUData *ilu_data = jxf_CTAlloc(jxf_ParILUData, 1);

   /* general data */
   (ilu_data -> global_solver)         = 0;
   (ilu_data -> matA)                  = NULL;
   (ilu_data -> matL)                  = NULL;
   (ilu_data -> matD)                  = NULL;
   (ilu_data -> matU)                  = NULL;
   (ilu_data -> matS)                  = NULL;
   (ilu_data -> schur_solver)          = NULL;
   (ilu_data -> schur_precond)         = NULL;
   (ilu_data -> rhs)                   = NULL;
   (ilu_data -> x)                     = NULL;

   (ilu_data -> droptol)               = jxf_TAlloc(JXF_Real, 3);
   (ilu_data -> own_droptol_data)      = 1;
   (ilu_data -> droptol)[0]            = 1.0e-02;/* droptol for B */
   (ilu_data -> droptol)[1]            = 1.0e-02;/* droptol for E and F */
   (ilu_data -> droptol)[2]            = 1.0e-02;/* droptol for S */
   (ilu_data -> lfil)                  = 0;
   (ilu_data -> maxRowNnz)             = 1000;
   (ilu_data -> CF_marker_array)       = NULL;
   (ilu_data -> perm)                  = NULL;
   (ilu_data -> qperm)                 = NULL;
   (ilu_data -> tol_ddPQ)              = 1.0e-01;

   (ilu_data -> F)                     = NULL;
   (ilu_data -> U)                     = NULL;
   (ilu_data -> Utemp)                 = NULL;
   (ilu_data -> Ftemp)                 = NULL;
   (ilu_data -> uext)                  = NULL;
   (ilu_data -> fext)                  = NULL;
   (ilu_data -> residual)              = NULL;
   (ilu_data -> rel_res_norms)         = NULL;

   (ilu_data -> num_iterations)        = 0;

   (ilu_data -> max_iter)              = 20;
   (ilu_data -> tol)                   = 1.0e-7;

   (ilu_data -> logging)               = 0;
   (ilu_data -> print_level)           = 0;

   (ilu_data -> l1_norms)              = NULL;

   (ilu_data -> operator_complexity)   = 0.;

   (ilu_data -> ilu_type)              = 0;
   (ilu_data -> nLU)                   = 0;
   (ilu_data -> nI)                    = 0;
   (ilu_data -> u_end)                 = NULL;

   /* reordering_type default to use local RCM */
   (ilu_data -> reordering_type) = 1;

   /* see jxf_ILUSetType for more default values */

   return (void *) ilu_data;
}

/* Destroy */
JXF_Int
jxf_ILUDestroy( void *data )
{
   jxf_ParILUData * ilu_data = (jxf_ParILUData*) data;

   /* final residual vector */
   if((ilu_data -> residual))
   {
      jxf_ParVectorDestroy( (ilu_data -> residual) );
      (ilu_data -> residual) = NULL;
   }
   if((ilu_data -> rel_res_norms))
   {
      jxf_TFree( (ilu_data -> rel_res_norms) );
      (ilu_data -> rel_res_norms) = NULL;
   }
   /* temp vectors for solve phase */
   if((ilu_data -> Utemp))
   {
      jxf_ParVectorDestroy( (ilu_data -> Utemp) );
      (ilu_data -> Utemp) = NULL;
   }
   if((ilu_data -> Ftemp))
   {
      jxf_ParVectorDestroy( (ilu_data -> Ftemp) );
      (ilu_data -> Ftemp) = NULL;
   }
   if(jxf_ParILUDataUExt(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataUExt(ilu_data));
      jxf_ParILUDataUExt(ilu_data) = NULL;
   }
   if(jxf_ParILUDataFExt(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataFExt(ilu_data));
      jxf_ParILUDataFExt(ilu_data) = NULL;
   }
   if((ilu_data -> rhs))
   {
      jxf_ParVectorDestroy( (ilu_data -> rhs) );
      (ilu_data -> rhs) = NULL;
   }
   if((ilu_data -> x))
   {
      jxf_ParVectorDestroy( (ilu_data -> x) );
      (ilu_data -> x) = NULL;
   }
   /* l1_norms */
   if((ilu_data -> l1_norms))
   {
      jxf_TFree((ilu_data -> l1_norms));
      (ilu_data -> l1_norms) = NULL;
   }

   /* u_end */
   if(jxf_ParILUDataUEnd(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataUEnd(ilu_data));
      jxf_ParILUDataUEnd(ilu_data) = NULL;
   }

   /* Factors */
   if(ilu_data -> matL)
   {
      jxf_ParCSRMatrixDestroy((ilu_data -> matL));
      (ilu_data -> matL) = NULL;
   }
   if(ilu_data -> matU)
   {
      jxf_ParCSRMatrixDestroy((ilu_data -> matU));
      (ilu_data -> matU) = NULL;
   }
   if(ilu_data -> matD)
   {
      jxf_TFree((ilu_data -> matD));
      (ilu_data -> matD) = NULL;
   }
   if(ilu_data -> matS)
   {
      jxf_ParCSRMatrixDestroy((ilu_data -> matS));
      (ilu_data -> matS) = NULL;
   }
   if(ilu_data -> schur_solver)
   {
      switch(ilu_data -> ilu_type){
         case 10: case 11: case 40: case 41:
            JXF_GMRESDestroy(ilu_data -> schur_solver); //GMRES for Schur
            break;
         case 20: case 21:
            jxf_NSHDestroy(jxf_ParILUDataSchurSolver(ilu_data));//NSH for Schur
            break;
         default:
            break;
      }
      (ilu_data -> schur_solver) = NULL;
   }
   if(ilu_data -> schur_precond)
   {
      switch(ilu_data -> ilu_type){
         case 10: case 11: case 40: case 41:
            JXF_ILUDestroy(ilu_data -> schur_precond); //ILU as precond for Schur
            break;
         default:
            break;
      }
      (ilu_data -> schur_precond) = NULL;
   }
   /* CF marker array */
   if((ilu_data -> CF_marker_array))
   {
      jxf_TFree((ilu_data -> CF_marker_array));
      (ilu_data -> CF_marker_array) = NULL;
   }
   /* permutation array */
   if((ilu_data -> perm))
   {
      jxf_TFree((ilu_data -> perm));
      (ilu_data -> perm) = NULL;
   }
   if((ilu_data -> qperm))
   {
      jxf_TFree((ilu_data -> qperm));
      (ilu_data -> qperm) = NULL;
   }
   /* droptol array */
   if((ilu_data -> own_droptol_data))
   {
      jxf_TFree((ilu_data -> droptol));
      (ilu_data -> own_droptol_data) = 0;
      (ilu_data -> droptol) = NULL;
   }
   if((ilu_data -> sp_own_droptol_data))
   {
      jxf_TFree((ilu_data -> sp_ilu_droptol));
      (ilu_data -> sp_own_droptol_data) = 0;
      (ilu_data -> sp_ilu_droptol) = NULL;
   }
   /* ilu data */
   jxf_TFree(ilu_data);

   return jxf_error_flag;
}

/* set fill level (for ilu(k)) */
JXF_Int
jxf_ILUSetLevelOfFill( void *ilu_vdata, JXF_Int lfil )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> lfil) = lfil;
   return jxf_error_flag;
}

/* set max non-zeros per row in factors (for ilut) */
JXF_Int
jxf_ILUSetMaxNnzPerRow( void *ilu_vdata, JXF_Int nzmax )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> maxRowNnz) = nzmax;
   return jxf_error_flag;
}

/* set threshold for dropping in LU factors (for ilut) */
JXF_Int
jxf_ILUSetDropThreshold( void *ilu_vdata, JXF_Real threshold )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> droptol)[0] = threshold;
   (ilu_data -> droptol)[1] = threshold;
   (ilu_data -> droptol)[2] = threshold;
   return jxf_error_flag;
}

/* set array of threshold for dropping in LU factors (for ilut) */
JXF_Int
jxf_ILUSetDropThresholdArray( void *ilu_vdata, JXF_Real *threshold )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   /* need to free memory if we own droptol array before */
   if((ilu_data -> own_droptol_data))
   {
      jxf_TFree((ilu_data -> droptol));
      (ilu_data -> own_droptol_data) = 0;
   }
   (ilu_data -> droptol) = threshold;
   return jxf_error_flag;
}

/* set if owns threshold data (for ilut) */
JXF_Int
jxf_ILUSetOwnDropThreshold( void *ilu_vdata, JXF_Int own_droptol_data )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> own_droptol_data) = own_droptol_data;
   return jxf_error_flag;
}

/* set ILU factorization type */
JXF_Int
jxf_ILUSetType( void *ilu_vdata, JXF_Int ilu_type )
{
   jxf_ParILUData *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ilu_type) = ilu_type;
   /* reset default value, not a large cost
    * assume we won't change back from
    */
   switch(ilu_type)
   {
      /* NSH type */
      case 20: case 21:
         /* set NSH Solver parameters */
         jxf_ParILUDataSchurNSHSolveMaxIter(ilu_data)        = 5;
         jxf_ParILUDataSchurNSHSolveTol(ilu_data)            = 1.0e-02;
         jxf_ParILUDataSchurSolverLogging(ilu_data)          = 0;
         jxf_ParILUDataSchurSolverPrintLevel(ilu_data)       = 0;
         if(jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data))
         {
            jxf_TFree(jxf_ParILUDataSchurNSHDroptol(ilu_data));
         }
         jxf_ParILUDataSchurNSHDroptol(ilu_data)             = jxf_ParILUDataDroptol(ilu_data);
         jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data)      = 0;

         /* set NHS inverse parameters */
         jxf_ParILUDataSchurNSHMaxNumIter(ilu_data)          = 2;/* kDim */
         jxf_ParILUDataSchurNSHMaxRowNnz(ilu_data)           = 1000;
         jxf_ParILUDataSchurNSHTol(ilu_data)                 = 1.0e-09;

         /* set MR inverse parameters */
         jxf_ParILUDataSchurMRMaxIter(ilu_data)              = 2;
         jxf_ParILUDataSchurMRColVersion(ilu_data)           = 0;/* sp_lfil */
         jxf_ParILUDataSchurMRMaxRowNnz(ilu_data)            = 200;
         jxf_ParILUDataSchurMRTol(ilu_data)                  = 1.0e-09;
         break;
      case 10: case 11: case 40: case 41:
         /* default data for schur solver */
         jxf_ParILUDataSchurGMRESKDim(ilu_data)              = 5;
         jxf_ParILUDataSchurGMRESTol(ilu_data)               = 1.0e-02;
         jxf_ParILUDataSchurGMRESAbsoluteTol(ilu_data)       = 0.0;
         jxf_ParILUDataSchurSolverLogging(ilu_data)          = 0;
         jxf_ParILUDataSchurSolverPrintLevel(ilu_data)       = 0;
         jxf_ParILUDataSchurGMRESRelChange(ilu_data)         = 0;

         /* default data for schur precond
          * default ILU0
          */
         jxf_ParILUDataSchurPrecondIluType(ilu_data)         = 0;
         jxf_ParILUDataSchurPrecondIluLfil(ilu_data)         = 0;
         jxf_ParILUDataSchurPrecondIluMaxRowNnz(ilu_data)    = 1000;
         if(jxf_ParILUDataSchurPrecondOwnDroptolData(ilu_data))
         {
            jxf_TFree(jxf_ParILUDataSchurPrecondIluDroptol(ilu_data));
         }
         jxf_ParILUDataSchurPrecondIluDroptol(ilu_data)      = jxf_ParILUDataDroptol(ilu_data);/* use same droptol */
         jxf_ParILUDataSchurPrecondOwnDroptolData(ilu_data)  = 0;
         jxf_ParILUDataSchurPrecondPrintLevel(ilu_data)      = 0;
         jxf_ParILUDataSchurPrecondMaxIter(ilu_data)         = 1;
         jxf_ParILUDataSchurPrecondTol(ilu_data)             = 1.0e-09;
         break;
      default:
         break;
   }

   return jxf_error_flag;
}

/* Set max number of iterations for ILU solver */
JXF_Int
jxf_ILUSetMaxIter( void *ilu_vdata, JXF_Int max_iter )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> max_iter) = max_iter;
   return jxf_error_flag;
}

/* Set convergence tolerance for ILU solver */
JXF_Int
jxf_ILUSetTol( void *ilu_vdata, JXF_Real tol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> tol) = tol;
   return jxf_error_flag;
}

/* Set print level for ilu solver */
JXF_Int
jxf_ILUSetPrintLevel( void *ilu_vdata, JXF_Int print_level )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> print_level) = print_level;
   return jxf_error_flag;
}

/* Set print level for ilu solver */
JXF_Int
jxf_ILUSetLogging( void *ilu_vdata, JXF_Int logging )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> logging) = logging;
   return jxf_error_flag;
}

/* Set type of reordering for local matrix */
JXF_Int
jxf_ILUSetLocalReordering( void *ilu_vdata, JXF_Int ordering_type )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> reordering_type) = ordering_type;
   return jxf_error_flag;
}

/* Set KDim (for GMRES) for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverKDIM( void *ilu_vdata, JXF_Int ss_kDim )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_kDim) = ss_kDim;
   return jxf_error_flag;
}

/* Set max iteration for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverMaxIter( void *ilu_vdata, JXF_Int ss_max_iter )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   switch(jxf_ParILUDataIluType(ilu_data))
   {
      case 10: case 11: case 40: case 41:
         /* GMRES
          * To avoid restart, GMRES kDim is equal to max num iter
          */
         jxf_ParILUDataSchurGMRESKDim(ilu_data) = ss_max_iter;
         break;
      case 20: case 21:
         /* set max num iter if use NSH solve */
         jxf_ParILUDataSchurNSHSolveMaxIter(ilu_data) = ss_max_iter;
         break;
      default:
         /* warning - not open yet
          *jxf_printf("Current type has no Schur System\n");
          */
         break;
   }
   return jxf_error_flag;
}

/* Set convergence tolerance for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverTol( void *ilu_vdata, JXF_Real ss_tol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_tol) = ss_tol;
   return jxf_error_flag;
}

/* Set absolute tolerance for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverAbsoluteTol( void *ilu_vdata, JXF_Real ss_absolute_tol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_absolute_tol) = ss_absolute_tol;
   return jxf_error_flag;
}

/* Set logging for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverLogging( void *ilu_vdata, JXF_Int ss_logging )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_logging) = ss_logging;
   return jxf_error_flag;
}

/* Set print level for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverPrintLevel( void *ilu_vdata, JXF_Int ss_print_level )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_print_level) = ss_print_level;
   return jxf_error_flag;
}

/* Set rel change (for GMRES) for Solver of Schur System */
JXF_Int
jxf_ILUSetSchurSolverRelChange( void *ilu_vdata, JXF_Int ss_rel_change )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> ss_rel_change) = ss_rel_change;
   return jxf_error_flag;
}

/* Set IUL type for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondILUType( void *ilu_vdata, JXF_Int sp_ilu_type )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_ilu_type) = sp_ilu_type;
   return jxf_error_flag;
}

/* Set IUL level of fill for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondILULevelOfFill( void *ilu_vdata, JXF_Int sp_ilu_lfil )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_ilu_lfil) = sp_ilu_lfil;
   return jxf_error_flag;
}

/* Set IUL max nonzeros per row for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondILUMaxNnzPerRow( void *ilu_vdata, JXF_Int sp_ilu_max_row_nnz )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_ilu_max_row_nnz) = sp_ilu_max_row_nnz;
   return jxf_error_flag;
}

/* Set IUL drop threshold for ILUT for Precond of Schur System
 * We don't want to influence the original ILU, so create new array if not own data
 */
JXF_Int
jxf_ILUSetSchurPrecondILUDropThreshold( void *ilu_vdata, JXF_Real sp_ilu_droptol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   if(jxf_ParILUDataSchurPrecondOwnDroptolData(ilu_data))
   {
      /* if we own data, just change our own data */
      (ilu_data -> sp_ilu_droptol)[0] = sp_ilu_droptol;
      (ilu_data -> sp_ilu_droptol)[1] = sp_ilu_droptol;
      (ilu_data -> sp_ilu_droptol)[2] = sp_ilu_droptol;
   }
   else
   {
      /* if we share data with other, create new one
       * becuase as default we use data from ILU, so we don't want to change it
       */
      jxf_ParILUDataSchurPrecondIluDroptol(ilu_data) = jxf_TAlloc(JXF_Real, 3);
      jxf_ParILUDataSchurPrecondOwnDroptolData(ilu_data)  = 1;
      jxf_ParILUDataSchurPrecondIluDroptol(ilu_data)[0]   = sp_ilu_droptol;
      jxf_ParILUDataSchurPrecondIluDroptol(ilu_data)[1]   = sp_ilu_droptol;
      jxf_ParILUDataSchurPrecondIluDroptol(ilu_data)[2]   = sp_ilu_droptol;
   }
   return jxf_error_flag;
}

/* Set array of IUL drop threshold for ILUT for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondILUDropThresholdArray( void *ilu_vdata, JXF_Real *sp_ilu_droptol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   /* need to free memory if we own droptol array before */
   if((ilu_data -> sp_own_droptol_data))
   {
      jxf_TFree((ilu_data -> sp_ilu_droptol));
      (ilu_data -> sp_own_droptol_data) = 0;
   }
   (ilu_data -> sp_ilu_droptol) = sp_ilu_droptol;
   return jxf_error_flag;
}

/* Set if owns drop threshold array for ILUT for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondILUOwnDropThreshold( void *ilu_vdata, JXF_Int sp_own_droptol_data )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_own_droptol_data) = sp_own_droptol_data;
   return jxf_error_flag;
}

/* Set print level for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondPrintLevel( void *ilu_vdata, JXF_Int sp_print_level )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_print_level) = sp_print_level;
   return jxf_error_flag;
}

/* Set max number of iterations for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondMaxIter( void *ilu_vdata, JXF_Int sp_max_iter )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_max_iter) = sp_max_iter;
   return jxf_error_flag;
}

/* Set onvergence tolerance for Precond of Schur System */
JXF_Int
jxf_ILUSetSchurPrecondTol( void *ilu_vdata, JXF_Int sp_tol )
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   (ilu_data -> sp_tol) = sp_tol;
   return jxf_error_flag;
}

/* Set tolorance for NSH for Schur System
 * We don't want to influence the original ILU, so create new array if not own data
 */
JXF_Int
jxf_ILUSetSchurNSHDropThreshold( void *ilu_vdata, JXF_Real threshold)
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   if(jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data))
   {
      jxf_ParILUDataSchurNSHDroptol(ilu_data)          = jxf_TAlloc(JXF_Real, 2);
      jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data)   = 1;
      jxf_ParILUDataSchurNSHDroptol(ilu_data)[0]       = threshold;
      jxf_ParILUDataSchurNSHDroptol(ilu_data)[1]       = threshold;
   }
   else
   {
      jxf_ParILUDataSchurNSHDroptol(ilu_data)[0]       = threshold;
      jxf_ParILUDataSchurNSHDroptol(ilu_data)[1]       = threshold;
   }
   return jxf_error_flag;
}

/* Set tolorance array for NSH for Schur System */
JXF_Int
jxf_ILUSetSchurNSHDropThresholdArray( void *ilu_vdata, JXF_Real *threshold)
{
   jxf_ParILUData   *ilu_data = (jxf_ParILUData*) ilu_vdata;
   if(jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data))
   {
      jxf_TFree(jxf_ParILUDataSchurNSHDroptol(ilu_data));
      jxf_ParILUDataSchurNSHOwnDroptolData(ilu_data) = 0;
   }
   jxf_ParILUDataSchurNSHDroptol(ilu_data) = threshold;
   return jxf_error_flag;
}

/* Get number of iterations for ILU solver */
JXF_Int
jxf_ILUGetNumIterations( void *ilu_vdata, JXF_Int *num_iterations )
{
   jxf_ParILUData  *ilu_data = (jxf_ParILUData*) ilu_vdata;

   if (!ilu_data)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   *num_iterations = ilu_data->num_iterations;

   return jxf_error_flag;
}

/* Get residual norms for ILU solver */
JXF_Int
jxf_ILUGetFinalRelativeResidualNorm( void *ilu_vdata, JXF_Real *res_norm )
{
   jxf_ParILUData  *ilu_data = (jxf_ParILUData*) ilu_vdata;

   if (!ilu_data)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   *res_norm = ilu_data->final_rel_residual_norm;

   return jxf_error_flag;
}

/* Print solver params */
JXF_Int
jxf_ILUWriteSolverParams(void *ilu_vdata)
{
   jxf_ParILUData  *ilu_data = (jxf_ParILUData*) ilu_vdata;
   jxf_printf("ILU Setup parameters: \n");
   jxf_printf("ILU factorization type: %d : ", (ilu_data -> ilu_type));
   switch(ilu_data -> ilu_type){
      case 0:
         jxf_printf("Block Jacobi with ILU(%d) \n", (ilu_data -> lfil));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 1:
         jxf_printf("Block Jacobi with ILUT \n");
         jxf_printf("drop tolerance for B = %e, E&F = %e, S = %e \n", jxf_ParILUDataDroptol(ilu_data)[0],jxf_ParILUDataDroptol(ilu_data)[1],jxf_ParILUDataDroptol(ilu_data)[2]);
         jxf_printf("Max nnz per row = %d \n", (ilu_data -> maxRowNnz));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 10:
         jxf_printf("ILU-GMRES with ILU(%d) \n", (ilu_data -> lfil));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 11:
         jxf_printf("ILU-GMRES with ILUT \n");
         jxf_printf("drop tolerance for B = %e, E&F = %e, S = %e \n", jxf_ParILUDataDroptol(ilu_data)[0],jxf_ParILUDataDroptol(ilu_data)[1],jxf_ParILUDataDroptol(ilu_data)[2]);
         jxf_printf("Max nnz per row = %d \n", (ilu_data -> maxRowNnz));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 20:
         jxf_printf("Newton-Schulz-Hotelling with ILU(%d) \n", (ilu_data -> lfil));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 21:
         jxf_printf("Newton-Schulz-Hotelling with ILUT \n");
         jxf_printf("drop tolerance for B = %e, E&F = %e, S = %e \n", jxf_ParILUDataDroptol(ilu_data)[0],jxf_ParILUDataDroptol(ilu_data)[1],jxf_ParILUDataDroptol(ilu_data)[2]);
         jxf_printf("Max nnz per row = %d \n", (ilu_data -> maxRowNnz));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 30:
         jxf_printf("RAS with ILU(%d) \n", (ilu_data -> lfil));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 31:
         jxf_printf("RAS with ILUT \n");
         jxf_printf("drop tolerance for B = %e, E&F = %e, S = %e \n", jxf_ParILUDataDroptol(ilu_data)[0],jxf_ParILUDataDroptol(ilu_data)[1],jxf_ParILUDataDroptol(ilu_data)[2]);
         jxf_printf("Max nnz per row = %d \n", (ilu_data -> maxRowNnz));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 40:
         jxf_printf("ddPQ-ILU-GMRES with ILU(%d) \n", (ilu_data -> lfil));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      case 41:
         jxf_printf("ddPQ-ILU-GMRES with ILUT \n");
         jxf_printf("drop tolerance for B = %e, E&F = %e, S = %e \n", jxf_ParILUDataDroptol(ilu_data)[0],jxf_ParILUDataDroptol(ilu_data)[1],jxf_ParILUDataDroptol(ilu_data)[2]);
         jxf_printf("Max nnz per row = %d \n", (ilu_data -> maxRowNnz));
         jxf_printf("Operator Complexity (Fill factor) = %f \n", (ilu_data -> operator_complexity));
         break;
      default: jxf_printf("Unknown type \n");
               break;
   }

   jxf_printf("\n ILU Solver Parameters: \n");
   jxf_printf("Max number of iterations: %d\n", (ilu_data -> max_iter));
   jxf_printf("Stopping tolerance: %e\n", (ilu_data -> tol));

   return jxf_error_flag;
}

/* helper functions */
/*
 * Add an element to the heap
 * I means JXF_Int
 * R means JXF_Real
 * max/min heap
 * r means heap goes from 0 to -1, -2 instead of 0 1 2
 * Ii and Ri means orderd by value of heap, like iw for ILU
 * heap: array of that heap
 * len: the current length of the heap
 * WARNING: You should first put that element to the end of the heap
 *    and add the length of heap by one before call this function.
 * the reason is that we don't want to change something outside the
 *    heap, so left it to the user
 */
JXF_Int
jxf_ILUMinHeapAddI(JXF_Int *heap, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p;
   len--;/* now len is the current index */
   while(len > 0)
   {
      /* get the parent index */
      p = (len-1)/2;
      if(heap[p] > heap[len])
      {
         /* this is smaller */
         jxf_swap(heap,p,len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapAddI for detail instructions */
JXF_Int
jxf_ILUMinHeapAddIIIi(JXF_Int *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p;
   len--;/* now len is the current index */
   while(len > 0)
   {
      /* get the parent index */
      p = (len-1)/2;
      if(heap[p] > heap[len])
      {
         /* this is smaller */
         jxf_swap(Ii1,heap[p],heap[len]);
         jxf_swap2i(heap,I1,p,len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapAddI for detail instructions */
JXF_Int
jxf_ILUMinHeapAddIRIi(JXF_Int *heap, JXF_Real *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p;
   len--;/* now len is the current index */
   while(len > 0)
   {
      /* get the parent index */
      p = (len-1)/2;
      if(heap[p] > heap[len])
      {
         /* this is smaller */
         jxf_swap(Ii1,heap[p],heap[len]);
         jxf_swap2(heap,I1,p,len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapAddI for detail instructions */
JXF_Int
jxf_ILUMaxHeapAddRabsIIi(JXF_Real *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p;
   len--;/* now len is the current index */
   while(len > 0)
   {
      /* get the parent index */
      p = (len-1)/2;
      if(jxf_abs(heap[p]) < jxf_abs(heap[len]))
      {
         /* this is smaller */
         jxf_swap(Ii1,heap[p],heap[len]);
         jxf_swap2(I1,heap,p,len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapAddI for detail instructions */
JXF_Int
jxf_ILUMaxrHeapAddRabsI(JXF_Real *heap, JXF_Int *I1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p;
   len--;/* now len is the current index */
   while(len > 0)
   {
      /* get the parent index */
      p = (len-1)/2;
      if(jxf_abs(heap[-p]) < jxf_abs(heap[-len]))
      {
         /* this is smaller */
         jxf_swap2(I1,heap,-p,-len);
         len = p;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/*
 * Swap the first element with the last element of the heap,
 *    reduce size by one, and maintain the heap structure
 * I means JXF_Int
 * R means JXF_Real
 * max/min heap
 * r means heap goes from 0 to -1, -2 instead of 0 1 2
 * Ii and Ri means orderd by value of heap, like iw for ILU
 * heap: aray of that heap
 * len: current length of the heap
 * WARNING: Remember to change the len youself
 */
JXF_Int
jxf_ILUMinHeapRemoveI(JXF_Int *heap, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p,l,r;
   len--;/* now len is the max index */
   /* swap the first element to last */
   jxf_swap(heap,0,len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while(l < len)
   {
      r = 2*p+2;
      /* two childs, pick the smaller one */
      l = r >= len || heap[l]<heap[r] ? l : r;
      if(heap[l]<heap[p])
      {
         jxf_swap(heap,l,p);
         p = l;
         l = 2*p+1;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapRemoveI for detail instructions */
JXF_Int
jxf_ILUMinHeapRemoveIIIi(JXF_Int *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p,l,r;
   len--;/* now len is the max index */
   /* swap the first element to last */
   jxf_swap(Ii1,heap[0],heap[len]);
   jxf_swap2i(heap,I1,0,len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while(l < len)
   {
      r = 2*p+2;
      /* two childs, pick the smaller one */
      l = r >= len || heap[l]<heap[r] ? l : r;
      if(heap[l]<heap[p])
      {
         jxf_swap(Ii1,heap[p],heap[l]);
         jxf_swap2i(heap,I1,l,p);
         p = l;
         l = 2*p+1;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapRemoveI for detail instructions */
JXF_Int
jxf_ILUMinHeapRemoveIRIi(JXF_Int *heap, JXF_Real *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p,l,r;
   len--;/* now len is the max index */
   /* swap the first element to last */
   jxf_swap(Ii1,heap[0],heap[len]);
   jxf_swap2(heap,I1,0,len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while(l < len)
   {
      r = 2*p+2;
      /* two childs, pick the smaller one */
      l = r >= len || heap[l]<heap[r] ? l : r;
      if(heap[l]<heap[p])
      {
         jxf_swap(Ii1,heap[p],heap[l]);
         jxf_swap2(heap,I1,l,p);
         p = l;
         l = 2*p+1;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapRemoveI for detail instructions */
JXF_Int
jxf_ILUMaxHeapRemoveRabsIIi(JXF_Real *heap, JXF_Int *I1, JXF_Int *Ii1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p,l,r;
   len--;/* now len is the max index */
   /* swap the first element to last */
   jxf_swap(Ii1,heap[0],heap[len]);
   jxf_swap2(I1,heap,0,len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while(l < len)
   {
      r = 2*p+2;
      /* two childs, pick the smaller one */
      l = r >= len || jxf_abs(heap[l])>jxf_abs(heap[r]) ? l : r;
      if(jxf_abs(heap[l])>jxf_abs(heap[p]))
      {
         jxf_swap(Ii1,heap[p],heap[l]);
         jxf_swap2(I1,heap,l,p);
         p = l;
         l = 2*p+1;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* see jxf_ILUMinHeapRemoveI for detail instructions */
JXF_Int
jxf_ILUMaxrHeapRemoveRabsI(JXF_Real *heap, JXF_Int *I1, JXF_Int len)
{
   /* parent, left, right */
   JXF_Int p,l,r;
   len--;/* now len is the max index */
   /* swap the first element to last */
   jxf_swap2(I1,heap,0,-len);
   p = 0;
   l = 1;
   /* while I'm still in the heap */
   while(l < len)
   {
      r = 2*p+2;
      /* two childs, pick the smaller one */
      l = r >= len || jxf_abs(heap[-l])>jxf_abs(heap[-r]) ? l : r;
      if(jxf_abs(heap[-l])>jxf_abs(heap[-p]))
      {
         jxf_swap2(I1,heap,-l,-p);
         p = l;
         l = 2*p+1;
      }
      else
      {
         break;
      }
   }
   return jxf_error_flag;
}

/* Split based on quick sort algorithm (avoid sorting the entire array)
 * find the largest k elements out of original array
 * array: input array for compare
 * I: integer array bind with array
 * k: largest k elements
 * len: length of the array
 */
JXF_Int
jxf_ILUMaxQSplitRabsI(JXF_Real *array, JXF_Int *I, JXF_Int left, JXF_Int bound, JXF_Int right)
{
   JXF_Int i, last;
   if (left >= right)
   {
      return jxf_error_flag;
   }
   jxf_swap2(I,array,left,(left+right)/2);
   last = left;
   for(i = left + 1 ; i <= right ; i ++)
   {
      if(jxf_abs(array[i]) > jxf_abs(array[left]))
      {
         jxf_swap2(I,array,++last,i);
      }
   }
   jxf_swap2(I,array,left,last);
   jxf_ILUMaxQSplitRabsI(array,I,left,bound,last-1);
   if(bound > last)
   {
       jxf_ILUMaxQSplitRabsI(array,I,last+1,bound,right);
   }

   return jxf_error_flag;
}

/* Helper function to search max value from a row
 * array: the array we work on
 * start: the start of the search range
 * end: the end of the search range
 * nLU: ignore rows (new row index) after nLU
 * rperm: reverse permutation array rperm[old] = new.
 *        if rperm set to NULL, ingore nLU and rperm
 * value: return the value ge get (absolute value)
 * index: return the index of that value, could be NULL which means not return
 * l1_norm: return the l1_norm of the array, could be NULL which means no return
 * nnz: return the number of nonzeros inside this array, could be NULL which means no return
 */
JXF_Int
jxf_ILUMaxRabs(JXF_Real *array_data, JXF_Int *array_j, JXF_Int start, JXF_Int end, JXF_Int nLU, JXF_Int *rperm, JXF_Real *value, JXF_Int *index, JXF_Real *l1_norm, JXF_Int *nnz)
{
   JXF_Int i, idx, col;
   JXF_Real val, max_value, norm, nz;

   nz = 0;
   norm = 0.0;
   max_value = -1.0;
   idx = -1;
   if(rperm)
   {
      /* apply rperm and nLU */
      for(i = start ; i < end ; i ++)
      {
         col = rperm[array_j[i]];
         if(col > nLU)
         {
            /* this old column is in new external part */
            continue;
         }
         nz ++;
         val = jxf_abs(array_data[i]);
         norm += val;
         if(max_value < val)
         {
            max_value = val;
            idx = i;
         }
      }
   }
   else
   {
      /* basic search */
      for(i = start ; i < end ; i ++)
      {
         val = jxf_abs(array_data[i]);
         norm += val;
         if(max_value < val)
         {
            max_value = val;
            idx = i;
         }
      }
      nz = end - start;
   }

   *value = max_value;
   if(index)
   {
      *index = idx;
   }
   if(l1_norm)
   {
      *l1_norm = norm;
   }
   if(nnz)
   {
      *nnz = nz;
   }

   return jxf_error_flag;
}

/* Pre selection for ddPQ, this is the basic version considering row sparsity
 * n: size of matrix
 * nLU: size we consider ddPQ reorder, only first nLU*nLU block is considered
 * A_diag_i/j/data: information of A
 * tol: tol for ddPQ, normally between 0.1-0.3
 * *perm: current row order
 * *rperm: current column order
 * *pperm_pre: output ddPQ pre row roder
 * *qperm_pre: output ddPQ pre column order
 */
JXF_Int
jxf_ILUGetPermddPQPre(JXF_Int n, JXF_Int nLU, JXF_Int *A_diag_i, JXF_Int *A_diag_j, JXF_Real *A_diag_data, JXF_Real tol, JXF_Int *perm, JXF_Int *rperm,
                        JXF_Int *pperm_pre, JXF_Int *qperm_pre, JXF_Int *nB)
{
   JXF_Int   i, ii, nB_pre, k1, k2;
   JXF_Real  gtol, max_value, norm;

   JXF_Int   *jcol, *jnnz;
   JXF_Real  *weight;

   weight      = jxf_TAlloc(JXF_Real, nLU + 1);
   jcol        = jxf_TAlloc(JXF_Int, nLU + 1);
   jnnz        = jxf_TAlloc(JXF_Int, nLU + 1);

   max_value   = -1.0;
   /* first need to build gtol */
   for( ii = 0 ; ii < nLU ; ii ++)
   {
      /* find real row */
      i = perm[ii];
      k1 = A_diag_i[i];
      k2 = A_diag_i[i+1];
      /* find max|a| of that row and its index */
      jxf_ILUMaxRabs(A_diag_data, A_diag_j, k1, k2, nLU, rperm, weight + ii, jcol + ii, &norm, jnnz + ii);
      weight[ii] /= norm;
      if(weight[ii] > max_value)
      {
         max_value = weight[ii];
      }
   }

   gtol = tol * max_value;

   /* second loop to pre select B */
   nB_pre = 0;
   for( ii = 0 ; ii < nLU ; ii ++)
   {
      /* keep this row */
      if(weight[ii] > gtol)
      {
         weight[nB_pre] /= (JXF_Real)(jnnz[ii]);
         pperm_pre[nB_pre] = perm[ii];
         qperm_pre[nB_pre++] = A_diag_j[jcol[ii]];
      }
   }

   *nB = nB_pre;

   /* sort from small to large */
   jxf_qsort3(weight, pperm_pre, qperm_pre, 0, nB_pre-1);

   jxf_TFree(weight);
   jxf_TFree(jcol);
   jxf_TFree(jnnz);

   return jxf_error_flag;
}

/* Get ddPQ version perm array for ParCSR
 * Greedy matching selection
 * ddPQ is a two-side permutation for diagonal dominance
 * A: the input matrix
 * pperm: row permutation
 * qperm: col permutation
 * nB: the size of B block
 * nI: number of interial nodes
 * tol: the dropping tolorance for ddPQ
 * reordering_type: Type of reordering for the interior nodes.
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 */

JXF_Int
jxf_ILUGetPermddPQ(jxf_ParCSRMatrix *A, JXF_Int **io_pperm, JXF_Int **io_qperm, JXF_Real tol, JXF_Int *nB, JXF_Int *nI, JXF_Int reordering_type)
{
   JXF_Int         i, nB_pre, irow, jcol, nLU;
   JXF_Int         *pperm, *qperm;
   JXF_Int         *rpperm, *rqperm, *pperm_pre, *qperm_pre;

   /* data objects for A */
   jxf_CSRMatrix   *A_diag = jxf_ParCSRMatrixDiag(A);
   JXF_Int         *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int         *A_diag_j = jxf_CSRMatrixJ(A_diag);
   JXF_Real        *A_diag_data = jxf_CSRMatrixData(A_diag);

   /* problem size */
   JXF_Int         n = jxf_CSRMatrixNumRows(A_diag);

   /* 1: Setup and create memory
    */

   pperm             = NULL;
   qperm             = jxf_TAlloc(JXF_Int, n);
   rpperm            = jxf_TAlloc(JXF_Int, n);
   rqperm            = jxf_TAlloc(JXF_Int, n);

   /* 2: Find interior nodes first
    */
   jxf_ILUGetInteriorExteriorPerm( A, &pperm, &nLU, 0);
   *nI = nLU;

   /* 3: Pre selection on interial nodes
    * this pre selection puts external nodes to the last
    * also provide candidate rows for B block
    */

   /* build reverse permutation array
    * rperm[old] = new
    */
   for(i = 0 ; i < n ; i ++)
   {
      rpperm[pperm[i]] = i;
   }

   /* build place holder for pre selection pairs */
   pperm_pre = jxf_TAlloc(JXF_Int, nLU);
   qperm_pre = jxf_TAlloc(JXF_Int, nLU);

   /* pre selection */
   jxf_ILUGetPermddPQPre(n, nLU, A_diag_i, A_diag_j, A_diag_data, tol, pperm, rpperm, pperm_pre, qperm_pre, &nB_pre);

   /* 4: Build B block
    * Greedy selection
    */

   /* rperm[old] = new */
   for(i = 0 ; i < nLU ; i ++)
   {
      rpperm[pperm[i]] = -1;
   }

   jxf_TMemcpy( rqperm, rpperm, JXF_Int, n);
   jxf_TMemcpy( qperm, pperm, JXF_Int, n);

   /* we sort from small to large, so we need to go from back to start
    * we only need nB_pre to start the loop, after that we could use it for size of B
    */
   for(i = nB_pre-1, nB_pre = 0 ; i >=0 ; i --)
   {
      irow = pperm_pre[i];
      jcol = qperm_pre[i];

      /* this col is not yet taken */
      if(rqperm[jcol] < 0)
      {
         rpperm[irow] = nB_pre;
         rqperm[jcol] = nB_pre;
         pperm[nB_pre] = irow;
         qperm[nB_pre++] = jcol;
      }
   }

   /* 5: Complete the permutation
    * rperm[old] = new
    * those still mapped to a new index means not yet covered
    */
   nLU = nB_pre;
   for(i = 0 ; i < n ; i ++)
   {
      if(rpperm[i] < 0)
      {
         pperm[nB_pre++] = i;
      }
   }
   nB_pre = nLU;
   for(i = 0 ; i < n ; i ++)
   {
      if(rqperm[i] < 0)
      {
         qperm[nB_pre++] = i;
      }
   }

   /* Finishing up and free
    */

   switch(reordering_type)
   {
      case 0:
         /* no RCM in this case */
         break;
      case 1:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, nLU, &pperm, &qperm, 0);
         break;
      default:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, nLU, &pperm, &qperm, 0);
         break;
   }

   *nB = nLU;
   *io_pperm = pperm;
   *io_qperm = qperm;

   jxf_TFree( rpperm);
   jxf_TFree( rqperm);
   jxf_TFree( pperm_pre);
   jxf_TFree( qperm_pre);
   return jxf_error_flag;
}

/*
 * Get perm array from parcsr matrix based on diag and offdiag matrix
 * Just simply loop through the rows of offd of A, check for nonzero rows
 * Put interior nodes at the beginning
 * A: parcsr matrix
 * perm: permutation array
 * nLU: number of interial nodes
 * reordering_type: Type of (additional) reordering for the interior nodes.
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 */
JXF_Int
jxf_ILUGetInteriorExteriorPerm(jxf_ParCSRMatrix *A, JXF_Int **perm, JXF_Int *nLU, JXF_Int reordering_type)
{
   /* get basic information of A */
   JXF_Int            n = jxf_ParCSRMatrixNumRows(A);
   JXF_Int            i, j, first, last, start, end;
   JXF_Int            num_sends, send_map_start, send_map_end, col;
   jxf_CSRMatrix      *A_offd;
   JXF_Int            *A_offd_i;
   A_offd               = jxf_ParCSRMatrixOffd(A);
   A_offd_i             = jxf_CSRMatrixI(A_offd);
   first                = 0;
   last                 = n - 1;
   JXF_Int            *temp_perm = jxf_TAlloc(JXF_Int, n);
   JXF_Int            *marker = jxf_CTAlloc(JXF_Int, n);

   /* first get col nonzero from com_pkg */
   /* get comm_pkg, craete one if we not yet have one */
   jxf_ParCSRCommPkg  *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   /* now directly take adavantage of comm_pkg */
   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   for( i = 0 ; i < num_sends ; i ++ )
   {
      send_map_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg,i);
      send_map_end = jxf_ParCSRCommPkgSendMapStart(comm_pkg,i+1);
      for ( j = send_map_start ; j < send_map_end ; j ++)
      {
         col = jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j);
         if(marker[col] == 0)
         {
            temp_perm[last--] = col;
            marker[col] = -1;
         }
      }
   }

   /* now deal with the row */
   for( i = 0 ; i < n ; i ++)
   {
      if(marker[i] == 0)
      {
         start = A_offd_i[i];
         end = A_offd_i[i+1];
         if(start == end)
         {
            temp_perm[first++] = i;
         }
         else
         {
            temp_perm[last--] = i;
         }
      }
   }
   switch(reordering_type)
   {
      case 0:
         /* no RCM in this case */
         break;
      case 1:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, first, &temp_perm, &temp_perm, 1);
         break;
      default:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, first, &temp_perm, &temp_perm, 1);
         break;
   }

   /* set out values */
   *nLU = first;
   if((*perm) != NULL) jxf_TFree(*perm);
   *perm = temp_perm;

   jxf_TFree(marker);
   return jxf_error_flag;
}

/*
 * Get the (local) ordering of the diag (local) matrix (no permutation). This is the permutation used for the block-jacobi case
 * A: parcsr matrix
 * perm: permutation array
 * nLU: number of interior nodes
 * reordering_type: Type of (additional) reordering for the nodes.
 * Currently only supports RCM reordering. Set to 0 for no reordering.
 */
JXF_Int
jxf_ILUGetLocalPerm(jxf_ParCSRMatrix *A, JXF_Int **perm, JXF_Int *nLU, JXF_Int reordering_type)
{
   /* get basic information of A */
   JXF_Int            n = jxf_ParCSRMatrixNumRows(A);
   JXF_Int            i;
   JXF_Int            *temp_perm = jxf_TAlloc(JXF_Int, n);

   /* set perm array */
   for( i = 0 ; i < n ; i ++ )
   {
      temp_perm[i] = i;
   }
   switch(reordering_type)
   {
      case 0:
         /* no RCM in this case */
         break;
      case 1:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, n, &temp_perm, &temp_perm, 1);
         break;
      default:
         /* RCM */
         jxf_ILULocalRCM( jxf_ParCSRMatrixDiag(A), 0, n, &temp_perm, &temp_perm, 1);
         break;
   }
   *nLU = n;
   if((*perm) != NULL) jxf_TFree(*perm);
   *perm = temp_perm;

   return jxf_error_flag;
}

/* Build the expanded matrix for RAS-1
 * A: input ParCSR matrix
 * E_i, E_j, E_data: information for external matrix
 * rperm: reverse permutation to build real index, rperm[old] = new
 */
JXF_Int
jxf_ILUBuildRASExternalMatrix(jxf_ParCSRMatrix *A, JXF_Int *rperm, JXF_Int **E_i, JXF_Int **E_j, JXF_Real **E_data)
{
   JXF_Int                i, j, idx;
   JXF_Int   big_col;

   /* data objects for communication */
   MPI_Comm                 comm = jxf_ParCSRMatrixComm(A);
   JXF_Int                my_id;

   /* data objects for A */
   jxf_CSRMatrix          *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix          *A_offd = jxf_ParCSRMatrixOffd(A);
   JXF_Int   *A_col_starts = jxf_ParCSRMatrixColStarts(A);
   JXF_Int   *A_offd_colmap = jxf_ParCSRMatrixColMapOffd(A);
   JXF_Int                *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int                *A_offd_i = jxf_CSRMatrixI(A_offd);

   /* data objects for external A matrix */
   // Need to check the new version of jxf_ParcsrGetExternalRows
   jxf_CSRMatrix          *A_ext = NULL;
   // # up to local offd cols, no need to be JXF_Int
   JXF_Int                *A_ext_i = NULL;
   // Return global index, JXF_Int required
   JXF_Int   *A_ext_j = NULL;
   JXF_Real               *A_ext_data = NULL;

   /* data objects for output */
   JXF_Int                E_nnz;
   JXF_Int                *E_ext_i = NULL;
   // Local index, no need to use JXF_Int
   JXF_Int                *E_ext_j = NULL;
   JXF_Real               *E_ext_data = NULL;

   //guess non-zeros for E before start
   JXF_Int                E_init_alloc;

   /* size */
   JXF_Int                n = jxf_CSRMatrixNumCols(A_diag);
   JXF_Int                m = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int                A_diag_nnz = A_diag_i[n];
   JXF_Int                A_offd_nnz = A_offd_i[n];

   /* 1: Set up phase and get external rows
    * Use the HYPRE build-in function
    */

   /* MPI stuff */
   //jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   /* Param of jxf_ParcsrGetExternalRows:
    * jxf_ParCSRMatrix   *A          [in]  -> Input parcsr matrix.
    * JXF_Int            indies_len  [in]  -> Input length of indices_len array
    * JXF_Int            *indices    [in]  -> Input global indices of rows we want to get
    * jxf_CSRMatrix      **A_ext     [out] -> Return the external CSR matrix.
    * jxf_ParCSRCommPkg  commpkg_out [out] -> Return commpkg if set to a point. Use NULL here since we don't want it.
    */
   //   jxf_ParcsrGetExternalRows( A, m, A_offd_colmap, &A_ext, NULL );
   A_ext = jxf_ParCSRMatrixExtractBExt(A, A, 1);

   A_ext_i              = jxf_CSRMatrixI(A_ext);
   //This should be JXF_Int since this is global index, use big_j in csr */
   A_ext_j = jxf_CSRMatrixJ(A_ext);
   A_ext_data           = jxf_CSRMatrixData(A_ext);

   /* guess memory we need to allocate to E_j */
   E_init_alloc =  jxf_max( (JXF_Int) ( A_diag_nnz / (JXF_Real) n / (JXF_Real) n * (JXF_Real) m * (JXF_Real) m + A_offd_nnz), 1);

   /* Initial guess */
   E_ext_i     = jxf_TAlloc(JXF_Int, m + 1 );
   E_ext_j     = jxf_TAlloc(JXF_Int, E_init_alloc );
   E_ext_data  = jxf_TAlloc(JXF_Real, E_init_alloc );

   /* 2: Discard unecessary cols
    * Search A_ext_j, discard those cols not belong to current proc
    * First check diag, and search in offd_col_map
    */

   E_nnz       = 0;
   E_ext_i[0]  = 0;

   for( i = 0 ;  i < m ; i ++)
   {
      E_ext_i[i] = E_nnz;
      for( j = A_ext_i[i] ; j < A_ext_i[i+1] ; j ++)
      {
         big_col = A_ext_j[j];
         /* First check if that belongs to the diagonal part */
#ifdef JXF_NO_GLOBAL_PARTITION

         if( big_col >= A_col_starts[0] && big_col < A_col_starts[1] )
         {
            /* this is a diagonal entry, rperm (map old to new) and shift it */

            /* Note here, the result of big_col - A_col_starts[0] in no longer a JXF_Int */
            idx = (JXF_Int)(big_col - A_col_starts[0]);
            E_ext_j[E_nnz]       = rperm[idx];
            E_ext_data[E_nnz++]  = A_ext_data[j];
         }

#else
         if( big_col >= A_col_starts[my_id] && big_col < A_col_starts[my_id+1] )
         {
            /* this is a diagonal entry, rperm (map old to new) and shift it */

            /* Note here, the result of big_col - A_col_starts[0] in no longer a JXF_Int */
            idx = (JXF_Int)(big_col - A_col_starts[my_id]);
            E_ext_j[E_nnz]       = rperm[idx];
            E_ext_data[E_nnz++]  = A_ext_data[j];
         }
#endif
         /* If not, apply binary search to check if is offdiagonal */
         else
         {
            /* Search, result is not JXF_Int */
            E_ext_j[E_nnz] = jxf_BinarySearch( A_offd_colmap, big_col, m);
            if( E_ext_j[E_nnz] >= 0)
            {
               /* this is an offdiagonal entry */
               E_ext_j[E_nnz]      = E_ext_j[E_nnz] + n;
               E_ext_data[E_nnz++] = A_ext_data[j];
            }
            else
            {
               /* skip capacity check */
               continue;
            }
         }
         /* capacity check, allocate new memory when full */
         if(E_nnz >= E_init_alloc)
         {
            E_init_alloc   = E_init_alloc * EXPAND_FACT + 1;
            E_ext_j        = jxf_TReAlloc(E_ext_j, JXF_Int, E_init_alloc);
            E_ext_data     = jxf_TReAlloc(E_ext_data, JXF_Real, E_init_alloc);
         }
      }
   }
   E_ext_i[m] = E_nnz;

   /* 3: Free and finish up
    * Free memory, set E_i, E_j and E_data
    */

   *E_i     = E_ext_i;
   *E_j     = E_ext_j;
   *E_data  = E_ext_data;

   jxf_CSRMatrixDestroy(A_ext);

   return jxf_error_flag;

}

/* This function sort offdiagonal map as well as J array for offdiagonal part
 * A: The input CSR matrix
 */
JXF_Int
jxf_ILUSortOffdColmap(jxf_ParCSRMatrix *A)
{
   JXF_Int i;
   jxf_CSRMatrix *A_offd    = jxf_ParCSRMatrixOffd(A);
   JXF_Int *A_offd_j        = jxf_CSRMatrixJ(A_offd);
   JXF_Int *A_offd_colmap   = jxf_ParCSRMatrixColMapOffd(A);
   JXF_Int len              = jxf_CSRMatrixNumCols(A_offd);
   JXF_Int nnz              = jxf_CSRMatrixNumNonzeros(A_offd);
   JXF_Int *perm            = jxf_TAlloc(JXF_Int,len);
   JXF_Int *rperm           = jxf_TAlloc(JXF_Int,len);

   for(i = 0 ; i < len ; i ++)
   {
      perm[i] = i;
   }

   jxf_qsort2i(A_offd_colmap,perm,0,len-1);

   for(i = 0 ; i < len ; i ++)
   {
      rperm[perm[i]] = i;
   }

   for(i = 0 ; i < nnz ; i ++)
   {
      A_offd_j[i] = rperm[A_offd_j[i]];
   }

   jxf_TFree(perm);
   jxf_TFree(rperm);

   return jxf_error_flag;
}

/* This function computes the RCM ordering of a sub matrix of
 * sparse matrix B = A(perm,perm)
 * For nonsymmetrix problem, is the RCM ordering of B + B'
 * A: The input CSR matrix
 * start:      the start position of the submatrix in B
 * end:        the end position of the submatrix in B ( exclude end, [start,end) )
 * permp:      pointer to the row permutation array such that B = A(perm, perm)
 *             point to NULL if you want to work directly on A
 *             on return, permp will point to the new permutation where
 *             in [start, end) the matrix will reordered
 * qpermp:     pointer to the col permutation array such that B = A(perm, perm)
 *             point to NULL or equal to permp if you want symmetric order
 *             on return, qpermp will point to the new permutation where
 *             in [start, end) the matrix will reordered
 * sym:        set to nonzero to work on A only(symmetric), otherwise A + A'.
 *             WARNING: if you use non-symmetric reordering, that is,
 *             different row and col reordering, the resulting A might be non-symmetric.
 *             Be careful if you are using non-symmetric reordering
 */
JXF_Int
jxf_ILULocalRCM( jxf_CSRMatrix *A, JXF_Int start, JXF_Int end,
                     JXF_Int **permp, JXF_Int **qpermp, JXF_Int sym)
{
   JXF_Int               i, j, row, col, r1, r2;

   JXF_Int               num_nodes      = end - start;
   JXF_Int               n              = jxf_CSRMatrixNumRows(A);
   JXF_Int               ncol           = jxf_CSRMatrixNumCols(A);
   JXF_Int               *A_i           = jxf_CSRMatrixI(A);
   JXF_Int               *A_j           = jxf_CSRMatrixJ(A);
   jxf_CSRMatrix         *GT            = NULL;
   jxf_CSRMatrix         *GGT           = NULL;
   //    JXF_Int               *AAT_i         = NULL;
   //    JXF_Int               *AAT_j         = NULL;
   JXF_Int               A_nnz          = jxf_CSRMatrixNumNonzeros(A);
   jxf_CSRMatrix         *G             = NULL;
   JXF_Int               *G_i           = NULL;
   JXF_Int               *G_j           = NULL;
   JXF_Real              *G_data           = NULL;
   JXF_Int               *G_perm        = NULL;
   JXF_Int               G_nnz;
   JXF_Int               G_capacity;
   JXF_Int               *perm_temp     = NULL;
   JXF_Int               *perm          = *permp;
   JXF_Int               *qperm         = *qpermp;
   JXF_Int               *rqperm        = NULL;

   /* 1: Preprosessing
    * Check error in input, set some parameters
    */
   if(num_nodes <= 0)
   {
      /* don't do this if we are too small */
      return jxf_error_flag;
   }
   if(n!=ncol || end > n || start < 0)
   {
      /* don't do this if the input has error */
      jxf_printf("Error input, abort RCM\n");
      return jxf_error_flag;
   }
   if(!perm)
   {
      /* create permutation array if we don't have one yet */
      perm = jxf_TAlloc( JXF_Int, n);
      for(i = 0 ; i < n ; i ++)
      {
         perm[i] = i;
      }
   }
   if(!qperm)
   {
      /* symmetric reordering, just point it to row reordering */
      qperm = perm;
   }
   rqperm = jxf_TAlloc(JXF_Int, n);
   for(i = 0 ; i < n ; i ++)
   {
      rqperm[qperm[i]] = i;
   }
   /* 2: Build Graph
    * Build Graph for RCM ordering
    */
   G = jxf_CSRMatrixCreate(num_nodes, num_nodes, 0);
   jxf_CSRMatrixInitialize(G);
   jxf_CSRMatrixOwnsData(G) = 1;
   G_i = jxf_CSRMatrixI(G);
   if(sym)
   {
      /* Directly use A */
      G_nnz = 0;
      G_capacity = jxf_max(A_nnz * n * n / num_nodes / num_nodes - num_nodes, 1);
      G_j = jxf_TAlloc(JXF_Int, G_capacity);
      for(i = 0 ; i < num_nodes ; i ++)
      {
         G_i[i] = G_nnz;
         row = perm[i + start];
         r1 = A_i[row];
         r2 = A_i[row+1];
         for(j = r1 ; j < r2 ; j ++)
         {
            col = rqperm[A_j[j]];
            if(col != row && col >= start && col < end)
            {
               /* this is an entry in G */
               G_j[G_nnz++] = col - start;
               if(G_nnz >= G_capacity)
               {
                  JXF_Int tmp = G_capacity;
                  G_capacity = G_capacity * EXPAND_FACT + 1;
                  G_j = jxf_TReAlloc_v2(G_j, JXF_Int, tmp, JXF_Int, G_capacity);
               }
            }
         }
      }
      G_i[num_nodes] = G_nnz;
      if(G_nnz == 0)
      {
         //G has only diagonal, no need to do any kind of RCM
         jxf_TFree(G_j);
         jxf_TFree(rqperm);
         *permp   = perm;
         *qpermp  = qperm;
         jxf_CSRMatrixDestroy(G);
         return jxf_error_flag;
      }
      jxf_CSRMatrixJ(G) = G_j;
      jxf_CSRMatrixNumNonzeros(G) = G_nnz;
   }
   else
   {
      /* Use A + A' */
      G_nnz = 0;
      G_capacity = jxf_max(A_nnz * n * n / num_nodes / num_nodes - num_nodes, 1);
      G_j = jxf_TAlloc(JXF_Int, G_capacity);
      for(i = 0 ; i < num_nodes ; i ++)
      {
         G_i[i] = G_nnz;
         row = perm[i + start];
         r1 = A_i[row];
         r2 = A_i[row+1];
         for(j = r1 ; j < r2 ; j ++)
         {
            col = rqperm[A_j[j]];
            if(col != row && col >= start && col < end)
            {
               /* this is an entry in G */
               G_j[G_nnz++] = col - start;
               if(G_nnz >= G_capacity)
               {
                  JXF_Int tmp = G_capacity;
                  G_capacity = G_capacity * EXPAND_FACT + 1;
                  G_j = jxf_TReAlloc_v2(G_j, JXF_Int, tmp, JXF_Int, G_capacity);
               }
            }
         }
      }
      G_i[num_nodes] = G_nnz;
      if(G_nnz == 0)
      {
         //G has only diagonal, no need to do any kind of RCM
         jxf_TFree(G_j);
         jxf_TFree(rqperm);
         *permp   = perm;
         *qpermp  = qperm;
         jxf_CSRMatrixDestroy(G);
         return jxf_error_flag;
      }
      jxf_CSRMatrixJ(G) = G_j;
      G_data = jxf_CTAlloc(JXF_Real, G_nnz);
      jxf_CSRMatrixData(G) = G_data;
      jxf_CSRMatrixNumNonzeros(G) = G_nnz;

      /* now sum G with G' */
      jxf_CSRMatrixTranspose(G, &GT, 1);
      GGT = jxf_CSRMatrixAdd(G, GT);
      jxf_CSRMatrixDestroy(G);
      jxf_CSRMatrixDestroy(GT);
      G = GGT;
      GGT = NULL;
   }

   /* 3: Build Graph
    * Build RCM
    */
   /* no need to be shared, but perm should be shared */
   G_perm = jxf_TAlloc(JXF_Int, num_nodes);
   jxf_ILULocalRCMOrder( G, G_perm);

   /* 4: Post processing
    * Free, set value, return
    */

   /* update to new index */
   perm_temp = jxf_TAlloc(JXF_Int, num_nodes);
   for( i = 0 ; i < num_nodes ; i ++)
   {
      perm_temp[i] = perm[i + start];
   }
   for( i = 0 ; i < num_nodes ; i ++)
   {
      perm[i+start] = perm_temp[G_perm[i]];
   }
   if(perm != qperm)
   {
      for( i = 0 ; i < num_nodes ; i ++)
      {
         perm_temp[i] = qperm[i + start];
      }
      for( i = 0 ; i < num_nodes ; i ++)
      {
         qperm[i+start] = perm_temp[G_perm[i]];
      }
   }
   *permp   = perm;
   *qpermp  = qperm;
   jxf_CSRMatrixDestroy(G);

   jxf_TFree(G_perm);
   jxf_TFree(perm_temp);
   jxf_TFree(rqperm);

   return jxf_error_flag;
}

/* This function finds the unvisited node with the minimum degree
 */
JXF_Int
jxf_ILULocalRCMMindegree(JXF_Int n, JXF_Int *degree, JXF_Int *marker, JXF_Int *rootp)
{
    JXF_Int i;
    JXF_Int min_degree = n+1;
    JXF_Int root = 0;
    for(i = 0 ; i < n ; i ++)
    {
        if(marker[i] < 0)
        {
            if(degree[i] < min_degree)
            {
                root = i;
                min_degree = degree[i];
            }
        }
    }
    *rootp = root;
    return 0;
}

/* This function actually does the RCM ordering of a symmetric csr matrix (entire)
 * A: the csr matrix, A_data is not needed
 * perm: the permutation array, space should be allocated outside
 */
JXF_Int
jxf_ILULocalRCMOrder( jxf_CSRMatrix *A, JXF_Int *perm)
{
   JXF_Int      i, root;
   JXF_Int      *degree     = NULL;
   JXF_Int      *marker     = NULL;
   JXF_Int      *A_i        = jxf_CSRMatrixI(A);
   JXF_Int      n           = jxf_CSRMatrixNumRows(A);
   JXF_Int      current_num;
   /* get the degree for each node */
   degree = jxf_TAlloc(JXF_Int, n);
   marker = jxf_TAlloc(JXF_Int, n);
   for(i = 0 ; i < n ; i ++)
   {
      degree[i] = A_i[i+1] - A_i[i];
      marker[i] = -1;
   }

   /* start RCM loop */
   current_num = 0;
   while(current_num < n)
   {
      jxf_ILULocalRCMMindegree( n, degree, marker, &root);
      /* This is a new connect component */
      jxf_ILULocalRCMFindPPNode(A, &root, marker);

      /* Numbering of this component */
      jxf_ILULocalRCMNumbering(A, root, marker, perm, &current_num);
   }

   /* free */
   jxf_TFree(degree);
   jxf_TFree(marker);
   return jxf_error_flag;
}

/* This function find a pseudo-peripheral node start from root
 * A: the csr matrix, A_data is not needed
 * rootp: pointer to the root, on return will be a end of the pseudo-peripheral
 * marker: the marker array for unvisited node
 */
JXF_Int
jxf_ILULocalRCMFindPPNode( jxf_CSRMatrix *A, JXF_Int *rootp, JXF_Int *marker)
{
   JXF_Int      i, r1, r2, row, min_degree, lev_degree, nlev, newnlev;

   JXF_Int      root           = *rootp;
   JXF_Int      n              = jxf_CSRMatrixNumRows(A);
   JXF_Int      *A_i           = jxf_CSRMatrixI(A);
   /* at most n levels */
   JXF_Int      *level_i       = jxf_TAlloc(JXF_Int, n+1);
   JXF_Int      *level_j       = jxf_TAlloc(JXF_Int, n);

   /* build initial level structure from root */
   jxf_ILULocalRCMBuildLevel( A, root, marker, level_i, level_j, &newnlev);

   nlev  = newnlev - 1;
   while(nlev < newnlev)
   {
      nlev = newnlev;
      r1 =  level_i[nlev-1];
      r2 =  level_i[nlev];
      min_degree = n;
      for(i = r1 ; i < r2 ; i ++)
      {
         /* select the last level, pick min-degree node */
         row = level_j[i];
         lev_degree = A_i[row+1] - A_i[row];
         if(min_degree > lev_degree)
         {
            min_degree = lev_degree;
            root = row;
         }
      }
      jxf_ILULocalRCMBuildLevel( A, root, marker, level_i, level_j, &newnlev);
   }

   *rootp = root;
   /* free */
   jxf_TFree(level_i);
   jxf_TFree(level_j);
   return jxf_error_flag;
}

/* This function build level structure start from root
 * A: the csr matrix, A_data is not needed
 * root: pointer to the root
 * marker: the marker array for unvisited node
 * level_i: points to the start/end of position on level_j, similar to CSR Matrix
 * level_j: store node number on each level
 * nlevp: return the number of level on this level structure
 */
JXF_Int
jxf_ILULocalRCMBuildLevel(jxf_CSRMatrix *A, JXF_Int root, JXF_Int *marker,
                              JXF_Int *level_i, JXF_Int *level_j, JXF_Int *nlevp)
{
   JXF_Int      i, j, l1, l2, l_current, r1, r2, rowi, rowj, nlev;
   JXF_Int      *A_i = jxf_CSRMatrixI(A);
   JXF_Int      *A_j = jxf_CSRMatrixJ(A);

   /* set first level first */
   level_i[0] = 0;
   level_j[0] = root;
   marker[root] = 0;
   nlev = 1;
   l1 = 0;
   l2 = 1;
   l_current = l2;

   //explore nbhds of all nodes in current level
   while(l2 > l1)
   {
      level_i[nlev++] = l2;
      /* loop through last level */
      for(i = l1 ; i < l2 ; i ++)
      {
         /* the node to explore */
         rowi = level_j[i];
         r1 = A_i[rowi];
         r2 = A_i[rowi + 1];
         for(j = r1 ; j < r2 ; j ++)
         {
            rowj = A_j[j];
            if( marker[rowj] < 0 )
            {
               /* Aha, an unmarked row */
               marker[rowj] = 0;
               level_j[l_current++] = rowj;
            }
         }
      }
      l1 = l2;
      l2 = l_current;
   }
   /* after this we always have a "ghost" last level */
   nlev --;

   /* reset marker */
   for(i = 0 ; i < l2 ; i ++)
   {
      marker[level_j[i]] = -1;
   }

   *nlevp = nlev;

   return jxf_error_flag;
}

/* This function generate numbering for a connect component
 * A: the csr matrix, A_data is not needed
 * root: pointer to the root
 * marker: the marker array for unvisited node
 * perm: permutation array
 * current_nump: number of nodes already have a perm value
 */

JXF_Int
jxf_ILULocalRCMNumbering(jxf_CSRMatrix *A, JXF_Int root, JXF_Int *marker, JXF_Int *perm, JXF_Int *current_nump)
{
    JXF_Int        i, j, l1, l2, r1, r2, rowi, rowj, row_start, row_end;
    JXF_Int        *A_i        = jxf_CSRMatrixI(A);
    JXF_Int        *A_j        = jxf_CSRMatrixJ(A);
    JXF_Int        current_num = *current_nump;


    marker[root]        = 0;
    l1                  = current_num;
    perm[current_num++] = root;
    l2                  = current_num;

    while(l2 > l1)
    {
       /* loop through all nodes is current level */
       for(i = l1 ; i < l2 ; i ++)
       {
          rowi = perm[i];
          r1 = A_i[rowi];
          r2 = A_i[rowi+1];
          row_start = current_num;
          for(j = r1 ; j < r2 ; j ++)
          {
             rowj = A_j[j];
             if(marker[rowj] < 0)
             {
                /* save the degree in marker and add it to perm */
                marker[rowj] = A_i[rowj+1] - A_i[rowj];
                perm[current_num++] = rowj;
             }
          }
          row_end = current_num;
          jxf_ILULocalRCMQsort(perm, row_start, row_end-1, marker);
       }
       l1 = l2;
       l2 = current_num;
    }

    //reverse
    jxf_ILULocalRCMReverse(perm, *current_nump, current_num-1);
    *current_nump = current_num;
    return jxf_error_flag;
}

/* This qsort is very specialized, not worth to put into utilities
 * Sort a part of array perm based on degree value (ascend)
 * That is, if degree[perm[i]] < degree[perm[j]], we should have i < j
 * perm: the perm array
 * start: start in perm
 * end: end in perm
 * degree: degree array
 */

JXF_Int
jxf_ILULocalRCMQsort(JXF_Int *perm, JXF_Int start, JXF_Int end, JXF_Int *degree)
{

    JXF_Int i, mid;
    if(start >= end)
    {
        return jxf_error_flag;
    }

    jxf_swap(perm, start, (start + end) / 2);
    mid = start;
    //loop to split
    for(i = start + 1 ; i <= end ; i ++)
    {
        if(degree[perm[i]] < degree[perm[start]])
        {
            jxf_swap(perm, ++mid, i);
        }
    }
    jxf_swap(perm, start, mid);
    jxf_ILULocalRCMQsort(perm, mid+1, end, degree);
    jxf_ILULocalRCMQsort(perm, start, mid-1, degree);
    return jxf_error_flag;
}

/* Last step in RCM, reverse it
 * perm: perm array
 * srart: start position
 * end: end position
 */

JXF_Int
jxf_ILULocalRCMReverse(JXF_Int *perm, JXF_Int start, JXF_Int end)
{
    JXF_Int     i, j;
    JXF_Int     mid = (start + end + 1) / 2;

    for(i = start, j = end ; i < mid ; i ++, j--)
    {
        jxf_swap(perm, i, j);
    }
   return jxf_error_flag;
}

/* NSH create and solve and help functions */

/* Create */
void *
jxf_NSHCreate()
{
   jxf_ParNSHData  *nsh_data;

   nsh_data = jxf_CTAlloc(jxf_ParNSHData,  1);

   /* general data */
   jxf_ParNSHDataMatA(nsh_data)                  = NULL;
   jxf_ParNSHDataMatM(nsh_data)                  = NULL;
   jxf_ParNSHDataF(nsh_data)                     = NULL;
   jxf_ParNSHDataU(nsh_data)                     = NULL;
   jxf_ParNSHDataResidual(nsh_data)              = NULL;
   jxf_ParNSHDataRelResNorms(nsh_data)           = NULL;
   jxf_ParNSHDataNumIterations(nsh_data)         = 0;
   jxf_ParNSHDataL1Norms(nsh_data)               = NULL;
   jxf_ParNSHDataFinalRelResidualNorm(nsh_data)  = 0.0;
   jxf_ParNSHDataTol(nsh_data)                   = 1e-09;
   jxf_ParNSHDataLogging(nsh_data)               = 2;
   jxf_ParNSHDataPrintLevel(nsh_data)            = 2;
   jxf_ParNSHDataMaxIter(nsh_data)               = 5;

   jxf_ParNSHDataOperatorComplexity(nsh_data)    = 0.0;
   jxf_ParNSHDataDroptol(nsh_data)               = jxf_TAlloc(JXF_Real,2);
   jxf_ParNSHDataOwnDroptolData(nsh_data)        = 1;
   jxf_ParNSHDataDroptol(nsh_data)[0]            = 1.0e-02;/* droptol for MR */
   jxf_ParNSHDataDroptol(nsh_data)[1]            = 1.0e-02;/* droptol for NSH */
   jxf_ParNSHDataUTemp(nsh_data)                 = NULL;
   jxf_ParNSHDataFTemp(nsh_data)                 = NULL;

   /* MR data */
   jxf_ParNSHDataMRMaxIter(nsh_data)             = 2;
   jxf_ParNSHDataMRTol(nsh_data)                 = 1e-09;
   jxf_ParNSHDataMRMaxRowNnz(nsh_data)           = 800;
   jxf_ParNSHDataMRColVersion(nsh_data)          = 0;

   /* NSH data */
   jxf_ParNSHDataNSHMaxIter(nsh_data)            = 2;
   jxf_ParNSHDataNSHTol(nsh_data)                = 1e-09;
   jxf_ParNSHDataNSHMaxRowNnz(nsh_data)          = 1000;

   return (void *) nsh_data;
}

/* Destroy */
JXF_Int
jxf_NSHDestroy( void *data )
{
   jxf_ParNSHData * nsh_data = (jxf_ParNSHData*) data;

   /* residual */
   if(jxf_ParNSHDataResidual(nsh_data))
   {
      jxf_ParVectorDestroy( jxf_ParNSHDataResidual(nsh_data) );
      jxf_ParNSHDataResidual(nsh_data) = NULL;
   }

   /* residual norms */
   if(jxf_ParNSHDataRelResNorms(nsh_data))
   {
      jxf_TFree( jxf_ParNSHDataRelResNorms(nsh_data) );
      jxf_ParNSHDataRelResNorms(nsh_data) = NULL;
   }

   /* l1 norms */
   if(jxf_ParNSHDataL1Norms(nsh_data))
   {
      jxf_TFree( jxf_ParNSHDataL1Norms(nsh_data) );
      jxf_ParNSHDataL1Norms(nsh_data) = NULL;
   }

   /* temp arrays */
   if(jxf_ParNSHDataUTemp(nsh_data))
   {
      jxf_ParVectorDestroy( jxf_ParNSHDataUTemp(nsh_data) );
      jxf_ParNSHDataUTemp(nsh_data) = NULL;
   }
   if(jxf_ParNSHDataFTemp(nsh_data))
   {
      jxf_ParVectorDestroy( jxf_ParNSHDataFTemp(nsh_data) );
      jxf_ParNSHDataFTemp(nsh_data) = NULL;
   }

   /* approx inverse matrix */
   if(jxf_ParNSHDataMatM(nsh_data))
   {
      jxf_ParCSRMatrixDestroy( jxf_ParNSHDataMatM(nsh_data) );
      jxf_ParNSHDataMatM(nsh_data) = NULL;
   }

   /* droptol array */
  if(jxf_ParNSHDataOwnDroptolData(nsh_data))
  {
     jxf_TFree(jxf_ParNSHDataDroptol(nsh_data));
     jxf_ParNSHDataOwnDroptolData(nsh_data) = 0;
     jxf_ParNSHDataDroptol(nsh_data) = NULL;
  }

   /* nsh data */
   jxf_TFree(nsh_data);

   return jxf_error_flag;
}

/* Print solver params */
JXF_Int
jxf_NSHWriteSolverParams(void *nsh_vdata)
{
   jxf_ParNSHData  *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_printf("Newton–Schulz–Hotelling Setup parameters: \n");
   jxf_printf("NSH max iterations = %d \n", jxf_ParNSHDataNSHMaxIter(nsh_data));
   jxf_printf("NSH drop tolerance = %e \n", jxf_ParNSHDataDroptol(nsh_data)[1]);
   jxf_printf("NSH max nnz per row = %d \n", jxf_ParNSHDataNSHMaxRowNnz(nsh_data));
   jxf_printf("MR max iterations = %d \n", jxf_ParNSHDataMRMaxIter(nsh_data));
   jxf_printf("MR drop tolerance = %e \n", jxf_ParNSHDataDroptol(nsh_data)[0]);
   jxf_printf("MR max nnz per row = %d \n", jxf_ParNSHDataMRMaxRowNnz(nsh_data));
   jxf_printf("Operator Complexity (Fill factor) = %f \n", jxf_ParNSHDataOperatorComplexity(nsh_data));
   jxf_printf("\n Newton–Schulz–Hotelling Solver Parameters: \n");
   jxf_printf("Max number of iterations: %d\n", jxf_ParNSHDataMaxIter(nsh_data));
   jxf_printf("Stopping tolerance: %e\n", jxf_ParNSHDataTol(nsh_data));

   return jxf_error_flag;
}

/* set print level */
JXF_Int
jxf_NSHSetPrintLevel( void *nsh_vdata, JXF_Int print_level )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataPrintLevel(nsh_data) = print_level;
   return jxf_error_flag;
}
/* set logging level */
JXF_Int
jxf_NSHSetLogging( void *nsh_vdata, JXF_Int logging )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataLogging(nsh_data) = logging;
   return jxf_error_flag;
}
/* set max iteration */
JXF_Int
jxf_NSHSetMaxIter( void *nsh_vdata, JXF_Int max_iter )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataMaxIter(nsh_data) = max_iter;
   return jxf_error_flag;
}
/* set solver iteration tol */
JXF_Int
jxf_NSHSetTol( void *nsh_vdata, JXF_Real tol )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataTol(nsh_data) = tol;
   return jxf_error_flag;
}
/* set global solver */
JXF_Int
jxf_NSHSetGlobalSolver( void *nsh_vdata, JXF_Int global_solver )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataGlobalSolver(nsh_data) = global_solver;
   return jxf_error_flag;
}
/* set all droptols */
JXF_Int
jxf_NSHSetDropThreshold( void *nsh_vdata, JXF_Real droptol )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataDroptol(nsh_data)[0] = droptol;
   jxf_ParNSHDataDroptol(nsh_data)[1] = droptol;
   return jxf_error_flag;
}
/* set array of droptols */
JXF_Int
jxf_NSHSetDropThresholdArray( void *nsh_vdata, JXF_Real *droptol )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   if(jxf_ParNSHDataOwnDroptolData(nsh_data))
   {
      jxf_TFree(jxf_ParNSHDataDroptol(nsh_data));
      jxf_ParNSHDataOwnDroptolData(nsh_data) = 0;
   }
   jxf_ParNSHDataDroptol(nsh_data) = droptol;
   return jxf_error_flag;
}
/* set own data */
JXF_Int
jxf_NSHSetOwnDroptolData( void *nsh_vdata, JXF_Int own_droptol_data )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataOwnDroptolData(nsh_data) = own_droptol_data;
   return jxf_error_flag;
}
/* set MR max iter */
JXF_Int
jxf_NSHSetMRMaxIter( void *nsh_vdata, JXF_Int mr_max_iter )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataMRMaxIter(nsh_data) = mr_max_iter;
   return jxf_error_flag;
}
/* set MR tol */
JXF_Int
jxf_NSHSetMRTol( void *nsh_vdata, JXF_Real mr_tol )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataMRTol(nsh_data) = mr_tol;
   return jxf_error_flag;
}
/* set MR max nonzeros of a row */
JXF_Int
jxf_NSHSetMRMaxRowNnz( void *nsh_vdata, JXF_Int mr_max_row_nnz )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataMRMaxRowNnz(nsh_data) = mr_max_row_nnz;
   return jxf_error_flag;
}
/* set MR version, column version or global version */
JXF_Int
jxf_NSHSetColVersion( void *nsh_vdata, JXF_Int mr_col_version )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataMRColVersion(nsh_data) = mr_col_version;
   return jxf_error_flag;
}
/* set NSH max iter */
JXF_Int
jxf_NSHSetNSHMaxIter( void *nsh_vdata, JXF_Int nsh_max_iter )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataNSHMaxIter(nsh_data) = nsh_max_iter;
   return jxf_error_flag;
}
/* set NSH tol */
JXF_Int
jxf_NSHSetNSHTol( void *nsh_vdata, JXF_Real nsh_tol )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataNSHTol(nsh_data) = nsh_tol;
   return jxf_error_flag;
}
/* set NSH max nonzeros of a row */
JXF_Int
jxf_NSHSetNSHMaxRowNnz( void *nsh_vdata, JXF_Int nsh_max_row_nnz )
{
   jxf_ParNSHData   *nsh_data = (jxf_ParNSHData*) nsh_vdata;
   jxf_ParNSHDataNSHMaxRowNnz(nsh_data) = nsh_max_row_nnz;
   return jxf_error_flag;
}

/* Compute the F norm of CSR matrix
 * A: the target CSR matrix
 * norm_io: output
 */
JXF_Int
jxf_CSRMatrixNormFro(jxf_CSRMatrix *A, JXF_Real *norm_io)
{
   JXF_Real norm = 0.0;
   JXF_Real *data = jxf_CSRMatrixData(A);
   JXF_Int i,k;
   k = jxf_CSRMatrixNumNonzeros(A);
   /* main loop */
   for(i = 0 ; i < k ; i ++)
   {
      norm += data[i] * data[i];
   }
   *norm_io = sqrt(norm);
   return jxf_error_flag;

}

/* Compute the norm of I-A where I is identity matrix and A is a CSR matrix
 * A: the target CSR matrix
 * norm_io: the output
 */
JXF_Int
jxf_CSRMatrixResNormFro(jxf_CSRMatrix *A, JXF_Real *norm_io)
{
   JXF_Real        norm = 0.0, value;
   JXF_Int         i, j, k1, k2, n;
   JXF_Int         *idx  = jxf_CSRMatrixI(A);
   JXF_Int         *cols = jxf_CSRMatrixJ(A);
   JXF_Real        *data = jxf_CSRMatrixData(A);

   n = jxf_CSRMatrixNumRows(A);
   /* main loop to sum up data */
   for(i = 0 ; i < n ; i ++)
   {
      k1 = idx[i];
      k2 = idx[i+1];
      /* check if we have diagonal in A */
      if(k2 > k1)
      {
         if(cols[k1] == i)
         {
            /* reduce 1 on diagonal */
            value = data[k1] - 1.0;
            norm += value * value;
         }
         else
         {
            /* we don't have diagonal in A, so we need to add 1 to norm */
            norm += 1.0;
            norm += data[k1] * data[k1];
         }
      }
      else
      {
         /* we don't have diagonal in A, so we need to add 1 to norm */
         norm += 1.0;
      }
      /* and the rest of the code */
      for(j = k1 + 1 ; j < k2 ; j ++)
      {
         norm += data[j] * data[j];
      }
   }
   *norm_io = sqrt(norm);
   return jxf_error_flag;
}

/* Compute the F norm of ParCSR matrix
 * A: the target CSR matrix
 */
JXF_Int
jxf_ParCSRMatrixNormFro(jxf_ParCSRMatrix *A, JXF_Real *norm_io)
{
   JXF_Real        local_norm = 0.0;
   JXF_Real        global_norm;
   MPI_Comm          comm = jxf_ParCSRMatrixComm(A);

   jxf_CSRMatrix   *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix   *A_offd = jxf_ParCSRMatrixOffd(A);

   jxf_CSRMatrixNormFro(A_diag, &local_norm);
   /* use global_norm to store offd for now */
   jxf_CSRMatrixNormFro(A_offd, &global_norm);

   /* square and sum them */
   local_norm *= local_norm;
   local_norm += global_norm*global_norm;

   /* do communication to get global total sum */
   jxf_MPI_Allreduce(&local_norm, &global_norm, 1, JXF_MPI_REAL, MPI_SUM, comm);

   *norm_io = sqrt(global_norm);
   return jxf_error_flag;

}

/* Compute the F norm of ParCSR matrix
 * Norm of I-A
 * A: the target CSR matrix
 */
JXF_Int
jxf_ParCSRMatrixResNormFro(jxf_ParCSRMatrix *A, JXF_Real *norm_io)
{
   JXF_Real        local_norm = 0.0;
   JXF_Real        global_norm;
   MPI_Comm          comm = jxf_ParCSRMatrixComm(A);

   jxf_CSRMatrix   *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix   *A_offd = jxf_ParCSRMatrixOffd(A);

   /* compute I-A for diagonal */
   jxf_CSRMatrixResNormFro(A_diag, &local_norm);
   /* use global_norm to store offd for now */
   jxf_CSRMatrixNormFro(A_offd, &global_norm);

   /* square and sum them */
   local_norm *= local_norm;
   local_norm += global_norm*global_norm;

   /* do communication to get global total sum */
   jxf_MPI_Allreduce(&local_norm, &global_norm, 1, JXF_MPI_REAL, MPI_SUM, comm);

   *norm_io = sqrt(global_norm);
   return jxf_error_flag;

}

/* Compute the trace of CSR matrix
 * A: the target CSR matrix
 * trace_io: the output trace
 */
JXF_Int
jxf_CSRMatrixTrace(jxf_CSRMatrix *A, JXF_Real *trace_io)
{
   JXF_Real  trace = 0.0;
   JXF_Int   *idx = jxf_CSRMatrixI(A);
   JXF_Int   *cols = jxf_CSRMatrixJ(A);
   JXF_Real  *data = jxf_CSRMatrixData(A);
   JXF_Int i,k1,k2,n;
   n = jxf_CSRMatrixNumRows(A);
   for(i = 0 ; i < n ; i ++)
   {
      k1 = idx[i];
      k2 = idx[i+1];
      if(cols[k1] == i && k2 > k1)
      {
         /* only add when diagonal is nonzero */
         trace += data[k1];
      }
   }

   *trace_io = trace;
   return jxf_error_flag;

}

/* Scale CSR matrix A = scalar * A
 * A: the target CSR matrix
 * scalar: real number
 */
JXF_Int
jxf_CSRMatrixScaleH(jxf_CSRMatrix *A, JXF_Real scalar)
{
   JXF_Real  *data = jxf_CSRMatrixData(A);
   JXF_Int   i,k;
   k = jxf_CSRMatrixNumNonzeros(A);
   for(i = 0 ; i < k ; i ++)
   {
      data[i] *= scalar;
   }
   return jxf_error_flag;
}

/* Scale ParCSR matrix A = scalar * A
 * A: the target CSR matrix
 * scalar: real number
 */
JXF_Int
jxf_ParCSRMatrixScaleH(jxf_ParCSRMatrix *A, JXF_Real scalar)
{
   jxf_CSRMatrix   *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix   *A_offd = jxf_ParCSRMatrixOffd(A);
   /* each thread scale local diag and offd */
   jxf_CSRMatrixScaleH(A_diag, scalar);
   jxf_CSRMatrixScaleH(A_offd, scalar);
   return jxf_error_flag;
}

/* Apply dropping to CSR matrix
 * A: the target CSR matrix
 * droptol: all entries have smaller absolute value than this will be dropped
 * max_row_nnz: max nonzeros allowed for each row, only largest max_row_nnz kept
 * we NEVER drop diagonal entry if exists
 */
JXF_Int
jxf_CSRMatrixDropInplace(jxf_CSRMatrix *A, JXF_Real droptol, JXF_Int max_row_nnz)
{
   JXF_Int      i, j, k1, k2;
   JXF_Int      *idx, len, drop_len;
   JXF_Real     *data, value, itol, norm;

   /* info of matrix A */
   JXF_Int      n = jxf_CSRMatrixNumRows(A);
   JXF_Int      m = jxf_CSRMatrixNumCols(A);
   JXF_Int      *A_i = jxf_CSRMatrixI(A);
   JXF_Int      *A_j = jxf_CSRMatrixJ(A);
   JXF_Real     *A_data = jxf_CSRMatrixData(A);
   JXF_Real     nnzA = jxf_CSRMatrixNumNonzeros(A);

   /* new data */
   JXF_Int      *new_i;
   JXF_Int      *new_j;
   JXF_Real     *new_data;

   /* memory */
   JXF_Int      capacity;
   JXF_Int      ctrA;

   /* setup */
   capacity = nnzA*0.3+1;
   ctrA = 0;
   new_i = jxf_TAlloc(JXF_Int, n+1);
   new_j = jxf_TAlloc(JXF_Int, capacity);
   new_data = jxf_TAlloc(JXF_Real, capacity);

   idx = jxf_TAlloc(JXF_Int, m);
   data = jxf_TAlloc(JXF_Real, m);

   /* start of main loop */
   new_i[0] = 0;
   for(i = 0 ; i < n ; i ++)
   {
      len = 0;
      k1 = A_i[i];
      k2 = A_i[i+1];
      /* compute droptol for current row */
      norm = 0.0;
      for(j = k1 ; j < k2 ; j ++)
      {
         norm += jxf_abs(A_data[j]);
      }
      if(k2 > k1)
      {
         norm /= (JXF_Real)(k2 - k1);
      }
      itol = droptol * norm;
      /* we don't want to drop the diagonal entry, so use an if statement here */
      if(A_j[k1] == i)
      {
         /* we have diagonal entry, skip it */
         idx[len] = A_j[k1];
         data[len++] = A_data[k1];
         for(j = k1 + 1 ; j < k2 ; j ++)
         {
            value = A_data[j];
            if(jxf_abs(value) < itol)
            {
               /* skip small element */
               continue;
            }
            idx[len] = A_j[j];
            data[len++] = A_data[j];
         }

         /* now apply drop on length */
         if(len > max_row_nnz)
         {
            drop_len = max_row_nnz;
            jxf_ILUMaxQSplitRabsI( data + 1, idx + 1, 0, drop_len - 1 , len - 2);
         }
         else
         {
            /* don't need to sort, we keep all of them */
            drop_len = len;
         }
         /* copy data */
         while(ctrA + drop_len > capacity)
         {
            JXF_Int tmp = capacity;
            capacity = capacity * EXPAND_FACT + 1;
            new_j = jxf_TReAlloc_v2(new_j, JXF_Int, tmp, JXF_Int, capacity);
            new_data = jxf_TReAlloc_v2(new_data, JXF_Real, tmp, JXF_Real, capacity);
         }
         jxf_TMemcpy( new_j + ctrA, idx,JXF_Int, drop_len);
         jxf_TMemcpy( new_data + ctrA, data,JXF_Real, drop_len);
         ctrA += drop_len;
         new_i[i+1] = ctrA;
      }
      else
      {
         /* we don't have diagonal entry */
         for(j = k1 ; j < k2 ; j ++)
         {
            value = A_data[j];
            if(jxf_abs(value)<itol)
            {
               /* skip small element */
               continue;
            }
            idx[len] = A_j[j];
            data[len++] = A_data[j];
         }

         /* now apply drop on length */
         if(len > max_row_nnz)
         {
            drop_len = max_row_nnz;
            jxf_ILUMaxQSplitRabsI( data, idx, 0, drop_len, len - 1);
         }
         else
         {
            /* don't need to sort, we keep all of them */
            drop_len = len;
         }

         /* copy data */
         while(ctrA + drop_len > capacity)
         {
            JXF_Int tmp = capacity;
            capacity = capacity * EXPAND_FACT + 1;
            new_j = jxf_TReAlloc_v2(new_j, JXF_Int, tmp, JXF_Int, capacity);
            new_data = jxf_TReAlloc_v2(new_data, JXF_Real, tmp, JXF_Real, capacity);
         }
         jxf_TMemcpy( new_j + ctrA, idx,JXF_Int, drop_len);
         jxf_TMemcpy( new_data + ctrA, data,JXF_Real, drop_len);
         ctrA += drop_len;
         new_i[i+1] = ctrA;
      }
   }/* end of main loop */
   /* destory data if A own them */
   if(jxf_CSRMatrixOwnsData(A))
   {
      jxf_TFree(A_i);
      jxf_TFree(A_j);
      jxf_TFree(A_data);
   }

   jxf_CSRMatrixI(A) = new_i;
   jxf_CSRMatrixJ(A) = new_j;
   jxf_CSRMatrixData(A) = new_data;
   jxf_CSRMatrixNumNonzeros(A) = ctrA;
   jxf_CSRMatrixOwnsData(A) = 1;

   jxf_TFree(idx);
   jxf_TFree(data);

   return jxf_error_flag;
}

/* Compute the inverse with MR of original CSR matrix
 * Global(not by each column) and out place version
 * A: the input matrix
 * M: the output matrix
 * droptol: the dropping tolorance
 * tol: when to stop the iteration
 * eps_tol: to avoid divide by 0
 * max_row_nnz: max number of nonzeros per row
 * max_iter: max number of iterations
 * print_level: the print level of this algorithm
 */
JXF_Int
jxf_ILUCSRMatrixInverseSelfPrecondMRGlobal(jxf_CSRMatrix *matA, jxf_CSRMatrix **M, JXF_Real droptol,
                                               JXF_Real tol, JXF_Real eps_tol, JXF_Int max_row_nnz, JXF_Int max_iter,
                                               JXF_Int print_level )
{
   JXF_Int         i, k1, k2;
   JXF_Real        value, trace1, trace2, alpha, r_norm;

   /* martix A */
   JXF_Int         *A_i = jxf_CSRMatrixI(matA);
   JXF_Int         *A_j = jxf_CSRMatrixJ(matA);
   JXF_Real        *A_data = jxf_CSRMatrixData(matA);

   /* complexity */
   JXF_Real        nnzA = jxf_CSRMatrixNumNonzeros(matA);
   JXF_Real        nnzM = 0;

   /* inverse matrix */
   jxf_CSRMatrix   *inM = *M;
   jxf_CSRMatrix   *matM;
   JXF_Int         *M_i;
   JXF_Int         *M_j;
   JXF_Real        *M_data;

   /* idendity matrix */
   jxf_CSRMatrix   *matI;
   JXF_Int         *I_i;
   JXF_Int         *I_j;
   JXF_Real        *I_data;

   /* helper matrices */
   jxf_CSRMatrix   *matR;
   jxf_CSRMatrix   *matR_temp;
   jxf_CSRMatrix   *matZ;
   jxf_CSRMatrix   *matC;
   jxf_CSRMatrix   *matW;

   JXF_Real        time_s, time_e;

   JXF_Int         n = jxf_CSRMatrixNumRows(matA);

   /* create initial guess and matrix I */
   matM = jxf_CSRMatrixCreate(n,n,n);
   M_i = jxf_TAlloc(JXF_Int, n+1);
   M_j = jxf_TAlloc(JXF_Int, n);
   M_data = jxf_TAlloc(JXF_Real, n);

   matI = jxf_CSRMatrixCreate(n,n,n);
   I_i = jxf_TAlloc(JXF_Int, n+1);
   I_j = jxf_TAlloc(JXF_Int, n);
   I_data = jxf_TAlloc(JXF_Real, n);

   /* now loop to create initial guess */
   M_i[0] = 0;
   I_i[0] = 0;
   for(i = 0 ; i < n ; i ++)
   {
      M_i[i+1] = i+1;
      M_j[i] = i;
      k1 = A_i[i];
      k2 = A_i[i+1];
      if(k2 > k1)
      {
         if(A_j[k1] == i)
         {
            value = A_data[k1];
            if(jxf_abs(value) < MAT_TOL)
            {
               value = 1.0;
            }
            M_data[i] = 1.0/value;
         }
         else
         {
            M_data[i] = 1.0;
         }
      }
      else
      {
         M_data[i] = 1.0;
      }
      I_i[i+1] = i+1;
      I_j[i] = i;
      I_data[i] = 1.0;
   }

   jxf_CSRMatrixI(matM) = M_i;
   jxf_CSRMatrixJ(matM) = M_j;
   jxf_CSRMatrixData(matM) = M_data;
   jxf_CSRMatrixOwnsData(matM) = 1;

   jxf_CSRMatrixI(matI) = I_i;
   jxf_CSRMatrixJ(matI) = I_j;
   jxf_CSRMatrixData(matI) = I_data;
   jxf_CSRMatrixOwnsData(matI) = 1;

   /* now start the main loop */
   if(print_level > 1)
   {
      /* time the iteration */
      time_s = jxf_MPI_Wtime();
   }

   /* main loop */
   for(i = 0 ; i < max_iter ; i ++)
   {
      nnzM = jxf_CSRMatrixNumNonzeros(matM);
      /* R = I - AM */
      matR_temp = jxf_CSRMatrixMultiply(matA,matM);

      jxf_CSRMatrixScaleH(matR_temp, -1.0);

      matR = jxf_CSRMatrixAdd(matI,matR_temp);
      jxf_CSRMatrixDestroy(matR_temp);

      /* r_norm */
      jxf_CSRMatrixNormFro(matR, &r_norm);
      if(r_norm < tol)
      {
         break;
      }

      /* Z = MR and dropping */
      matZ = jxf_CSRMatrixMultiply(matM, matR);
      //jxf_CSRMatrixNormFro(matZ, &z_norm);
      jxf_CSRMatrixDropInplace(matZ, droptol, max_row_nnz);

      /* C = A*Z */
      matC = jxf_CSRMatrixMultiply(matA, matZ);

      /* W = R' * C */
      jxf_CSRMatrixTranspose(matR,&matR_temp,1);
      matW = jxf_CSRMatrixMultiply(matR_temp,matC);

      /* trace and alpha */
      jxf_CSRMatrixTrace(matW, &trace1);
      jxf_CSRMatrixNormFro(matC, &trace2);
      trace2 *= trace2;

      if(jxf_abs(trace2) < eps_tol)
      {
         break;
      }

      alpha = trace1 / trace2;

      /* M - M + alpha * Z */
      jxf_CSRMatrixScaleH(matZ, alpha);

      jxf_CSRMatrixDestroy(matR);
      matR = jxf_CSRMatrixAdd(matM, matZ);
      jxf_CSRMatrixDestroy(matM);
      matM = matR;

      jxf_CSRMatrixDestroy(matZ);
      jxf_CSRMatrixDestroy(matW);
      jxf_CSRMatrixDestroy(matC);
      jxf_CSRMatrixDestroy(matR_temp);

   }/* end of main loop i for compute inverse matrix */

   /* time if we need to print */
   if(print_level > 1)
   {
      time_e = jxf_MPI_Wtime();
      if(i == 0)
      {
         i = 1;
      }
      jxf_printf("matrix size %5d\nfinal norm at loop %5d is %16.12f, time per iteration is %16.12f, complexity is %16.12f out of maximum %16.12f\n",n,i,r_norm, (time_e-time_s)/i, nnzM/nnzA, n/nnzA*n);
   }

   jxf_CSRMatrixDestroy(matI);
   if(inM)
   {
      jxf_CSRMatrixDestroy(inM);
   }
   *M = matM;

   return jxf_error_flag;

}

/* Compute inverse with NSH method
 * Use MR to get local initial guess
 * A: input matrix
 * M: output matrix
 * droptol: droptol array. droptol[0] for MR and droptol[1] for NSH.
 * mr_tol: tol for stop iteration for MR
 * nsh_tol: tol for stop iteration for NSH
 * esp_tol: tol for avoid divide by 0
 * mr_max_row_nnz: max number of nonzeros for MR
 * nsh_max_row_nnz: max number of nonzeros for NSH
 * mr_max_iter: max number of iterations for MR
 * nsh_max_iter: max number of iterations for NSH
 * mr_col_version: column version of global version
 */
JXF_Int
jxf_ILUParCSRInverseNSH(jxf_ParCSRMatrix *A, jxf_ParCSRMatrix **M, JXF_Real *droptol, JXF_Real mr_tol,
                            JXF_Real nsh_tol, JXF_Real eps_tol, JXF_Int mr_max_row_nnz, JXF_Int nsh_max_row_nnz,
                            JXF_Int mr_max_iter, JXF_Int nsh_max_iter, JXF_Int mr_col_version,
                            JXF_Int print_level)
{
   JXF_Int               i;

   /* data slots for matrices */
   jxf_ParCSRMatrix      *matM = NULL;
   jxf_ParCSRMatrix      *inM = *M;
   jxf_ParCSRMatrix      *AM,*MAM;
   JXF_Real              norm, s_norm;
   MPI_Comm                comm = jxf_ParCSRMatrixComm(A);
   JXF_Int               myid;


   jxf_CSRMatrix         *A_diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix         *M_diag = NULL;
   jxf_CSRMatrix         *M_offd;
   JXF_Int               *M_offd_i;

   JXF_Real              time_s, time_e;

   JXF_Int               n = jxf_CSRMatrixNumRows(A_diag);

   /* setup */
   jxf_MPI_Comm_rank(comm, &myid);

   M_offd_i = jxf_TAlloc(JXF_Int, n+1);

   if(mr_col_version)
   {
      jxf_printf("Column version is not yet support, switch to global version\n");
   }

   /* call MR to build loacl initial matrix
    * droptol here should be larger
    * we want same number for MR and NSH to let user set them eaiser
    * but we don't want a too dense MR initial guess
    */
   jxf_ILUCSRMatrixInverseSelfPrecondMRGlobal(A_diag, &M_diag, droptol[0] * 10.0, mr_tol, eps_tol, mr_max_row_nnz, mr_max_iter, print_level );

   /* create empty offdiagonal */
   for(i = 0 ; i <= n ; i ++)
   {
      M_offd_i[i] = 0;
   }

   /* create parCSR matM */
   matM = jxf_ParCSRMatrixCreate( comm,
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixGlobalNumRows(A),
         jxf_ParCSRMatrixRowStarts(A),
         jxf_ParCSRMatrixColStarts(A),
         0,
         jxf_CSRMatrixNumNonzeros(M_diag),
         0 );

   jxf_CSRMatrixDestroy(jxf_ParCSRMatrixDiag(matM));
   jxf_ParCSRMatrixDiag(matM) = M_diag;

   M_offd = jxf_ParCSRMatrixOffd(matM);
   jxf_CSRMatrixI(M_offd) = M_offd_i;
   jxf_CSRMatrixOwnsData(M_offd) = 1;

   jxf_ParCSRMatrixSetColStartsOwner(matM,0);
   jxf_ParCSRMatrixSetRowStartsOwner(matM,0);

   /* now start NSH
    * Mj+1 = 2Mj - MjAMj
    */

   AM = jxf_ParMatmul(A, matM);
   jxf_ParCSRMatrixResNormFro(AM, &norm);
   s_norm = norm;
   jxf_ParCSRMatrixDestroy(AM);
   if(print_level > 1)
   {
      if(myid == 0)
      {
         jxf_printf("before NSH the norm is %16.12f\n", norm);
      }
      time_s = jxf_MPI_Wtime();
   }

   for(i = 0 ; i < nsh_max_iter ; i ++)
   {
      /* compute XjAXj */
      AM = jxf_ParMatmul(A, matM);
      jxf_ParCSRMatrixResNormFro(AM, &norm);
      if(norm < nsh_tol)
      {
         break;
      }
      MAM = jxf_ParMatmul(matM, AM);
      jxf_ParCSRMatrixDestroy(AM);

      /* apply dropping */
      //jxf_ParCSRMatrixNormFro(MAM, &norm);
      /* drop small entries based on 2-norm */
      jxf_ParCSRMatrixDropSmallEntries(MAM, droptol[1], 2);

      /* update Mj+1 = 2Mj - MjAMj
       * the result holds it own start/end data!
       */
      jxf_ParcsrAdd(2.0, matM,-1.0, MAM, &AM);
      jxf_ParCSRMatrixDestroy(matM);
      matM = AM;

      /* destroy */
      jxf_ParCSRMatrixDestroy(MAM);
   }

   if(print_level > 1)
   {
      time_e = jxf_MPI_Wtime();
      /* at this point of time, norm has to be already computed */
      if(i == 0)
      {
         i = 1;
      }
      if(myid == 0)
      {
         jxf_printf("after %5d NSH iterations the norm is %16.12f, time per iteration is %16.12f\n", i, norm, (time_e-time_s)/i);
      }
   }

   if(s_norm < norm)
   {
      /* the residual norm increase after NSH iteration, need to let user know */
      if(myid == 0)
      {
         jxf_printf("Warning: NSH divergence, probably bad approximate invese matrix.\n");
      }
   }

   if(inM)
   {
      jxf_ParCSRMatrixDestroy(inM);
   }
   *M = matM;

   return jxf_error_flag;
}
