//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csolver.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_CoarsestSolver
 * \brief Solver on coarsest grid for PAMG, which can 
 *        also be taken as an iterative solver.
 * \author peghoty
 * \date 2009/07/25
 */
JXF_Int
jxf_CoarsestSolver( void             *amg_vdata,
                   jxf_hpCSRMatrix  *hp_matrix, 
                   jxf_ParVector     *par_rhs,
                   jxf_ParVector     *par_app )
{
   //MPI_Comm        comm     = MPI_COMM_WORLD;  // ERROR in Clone cases.
   MPI_Comm         comm     = jxf_hpCSRMatrixComm(hp_matrix); 
   jxf_ParAMGData  *amg_data = amg_vdata;

   JXF_Real tol           = jxf_ParAMGDataTol(amg_data);
   JXF_Real weight        = (jxf_ParAMGDataRelaxWeight(amg_data))[0];
   JXF_Real omega         = (jxf_ParAMGDataOmega(amg_data))[0];
   JXF_Int    solverid      = jxf_ParAMGDataCoarsestSolverID(amg_data);
   JXF_Int    max_iter      = jxf_ParAMGDataMaxIter(amg_data); 
   JXF_Int    min_iter      = jxf_ParAMGDataMinIter(amg_data);  
   JXF_Int    conv_criteria = jxf_ParAMGDataConvCriteria(amg_data);
   JXF_Int    print_level   = jxf_ParAMGDataPrintLevel(amg_data);

   JXF_Real convfac_threshold = jxf_ParAMGDataConvFacThreshold(amg_data);

   JXF_Real **AI_measure_array  = jxf_ParAMGDataAIMeasureArray(amg_data);
   JXF_Int      *relax_marker_ai;
   JXF_Int      *relax_marker_ess;
   JXF_Int       n, num_ai_th_step, num_ai_th, num_ess_th;
   JXF_Real    ai_measure_th_min = 0.2;
   JXF_Real    ai_measure_th_max = 2.01;
//   JXF_Real    ai_factor = 0.5;
//   JXF_Real    ai_measure_th, ai_measure_th_l, ai_measure_th_h;
   JXF_Real    ai_measure_th_l, ai_measure_th_h;
   JXF_Real    ai_measure;
//   JXF_Int       i,j,nstep;
   JXF_Int       i;

   JXF_Int    iter, size;
   JXF_Int    num_procs, my_id;
   JXF_Int    iter_ai, iter_ess;

   JXF_Real   alpha          = 1.0;
   JXF_Real   beta           = -1.0;
//   JXF_Real   normup         = 0.0;
//   JXF_Real   normdown       = 0.0;
   JXF_Real   conv_factor    = 0.0;
   JXF_Real   resid_nrm      = 0.0;
   JXF_Real   resid_nrm_init = 0.0;
   JXF_Real   resid_nrm_0 = 0.0;
   JXF_Real   relative       = 0.0;
//   JXF_Real   old_relative   = 0.0;
   JXF_Real   rhs_norm       = 0.0;
//   JXF_Real   old_resid      = 0.0;
  

   jxf_Vector               *local_vector_ai = NULL;
   JXF_Real                  *local_data_ai = NULL;
   jxf_Vector               *local_vector_ess = NULL;
   JXF_Real                  *local_data_ess = NULL;
   JXF_Int                      local_size;
   JXF_Real   relative_ai       = 0.0;
   JXF_Real   rhs_norm_ai       = 0.0;
   JXF_Real   resid_nrm_ai      = 0.0;
   JXF_Real   resid_nrm_init_ai = 0.0;
   JXF_Real   relative_ess       = 0.0;
   JXF_Real   rhs_norm_ess       = 0.0;
   JXF_Real   resid_nrm_ess      = 0.0;
   JXF_Real   resid_nrm_init_ess = 0.0;
 
   JXF_Int     *partitioning   = NULL;

   jxf_ParVector  *Vtemp    = NULL;
   jxf_ParVector  *Utemp    = NULL;

   JXF_Int      print_level_ai = print_level;
   JXF_Int      print_level_ess = print_level;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == 0 && print_level)
   {
      jxf_CoarsestSolverInfo( comm, amg_data );    
   }

   size = jxf_ParVectorGlobalSize(par_rhs);
   partitioning = jxf_ParVectorPartitioning(par_rhs);


   Vtemp = jxf_ParVectorCreate(comm, size, partitioning);
   jxf_ParVectorSetPartitioningOwner(Vtemp, 0);
   jxf_ParVectorInitialize(Vtemp);

   if ( solverid == 9 || solverid == 10 )   /* Direct Solve: Gaussian Elimination */
   {
      jxf_GaussElimination(solverid, jxf_hpCSRMatrixPar(hp_matrix), par_rhs, par_app);

      jxf_ParVectorCopy(par_rhs, Vtemp);
      jxf_hpCSRMatrixMatvec(alpha, hp_matrix, par_app, beta, Vtemp);
      resid_nrm = jxf_ParVectorNorm2(Vtemp);

      if (my_id == 0 && print_level) 
      {
         jxf_printf(" Residual Norm-2 = %e\n\n\n", resid_nrm);
      }
   }
   else
   {
      Utemp = jxf_ParVectorCreate(comm, size, partitioning);
      jxf_ParVectorSetPartitioningOwner(Utemp, 0);
      jxf_ParVectorInitialize(Utemp);

      iter = 0;

      if (conv_criteria == 0)
      {
         jxf_ParVectorCopy(par_rhs, Vtemp);
         jxf_hpCSRMatrixMatvec(alpha, hp_matrix, par_app, beta, Vtemp);
         resid_nrm = jxf_ParVectorNorm2(Vtemp);

         resid_nrm_init = resid_nrm;
         rhs_norm = jxf_ParVectorNorm2(par_rhs);
         if (rhs_norm)
         {
            relative = resid_nrm_init / rhs_norm;
         }
         else
         {
            relative = resid_nrm_init;
         }
      }
      else if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
      {
         relative = 2*tol;  /* to make sure entering the "while loop" once at leaset peghoty 2009/07/24 */
      }

      if (my_id == 0 && print_level)
      {     
         if (conv_criteria == 0)
         { 
            jxf_printf("                                             Relative  \n");
            jxf_printf("                Res Norm-2      Factor      Res Norm-2 \n");
            jxf_printf("                ----------     --------    ------------\n");
            jxf_printf("    Initial    %e                %e\n", resid_nrm_init, relative);
         }
         if (conv_criteria == 1)
         {
            jxf_printf("                                           Relative   \n");
            jxf_printf("                Error Norm-1    Factor    Error Norm-1\n");
            jxf_printf("                ------------   --------  -------------\n");
         }
         if (conv_criteria == 11)
         {
            jxf_printf("                         Relative Point-wise \n");
            jxf_printf("                Factor      Error Norm-1     \n");
            jxf_printf("                -------    --------------    \n");
         }
         if (conv_criteria == 2)
         {
            jxf_printf("                                           Relative   \n");
            jxf_printf("                Error Norm-2    Factor    Error Norm-2\n");
            jxf_printf("                ------------   --------  -------------\n");           
         }
      }

      /* solverid >=90: AI-prior smoothing strategy. */
      if (solverid >= 90) {

         solverid = solverid-90;

         n = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix));
         relax_marker_ai = jxf_CTAlloc(JXF_Int, n);
         relax_marker_ess = jxf_CTAlloc(JXF_Int, n);
         for (i=0; i<n; i++) {
            relax_marker_ai[i] = 0;
            relax_marker_ess[i] = 1;
         }

         /* following code need to be optimized. */
         //ai_measure_th = ai_measure_th_max;
         ai_measure_th_h = ai_measure_th_max;
         ai_measure_th_l = ai_measure_th_min;
         num_ai_th = 0;
         num_ess_th = n;
         //while (ai_measure_th >= ai_measure_th_min) {

         num_ai_th_step = 0;
         for (i=0; i<n; i++) {
            ai_measure = AI_measure_array[0][i];
            //jxf_printf("ai_measure=%f, ai_measure_th_l=%f, ai_measure_th_h=%f \n", ai_measure, ai_measure_th_l, ai_measure_th_h);
            if  ( ai_measure > ai_measure_th_l && ai_measure <= ai_measure_th_h 
                  //&& CF_marker_array[level][i] == 1
                ) {
                relax_marker_ai[i] = 1;
                relax_marker_ess[i] = 0;
                num_ai_th++;
                num_ai_th_step++;
                num_ess_th--;
            }
         }

         //ai_measure_th = ai_measure_th_l;
         //ai_measure_th_h = ai_measure_th_l;
         //}

         /* compute initial rhs_norm. */
         if (conv_criteria == 0)
         {
            /* local rhs_norm */
            jxf_ParVectorCopy(par_rhs, Vtemp);
            jxf_ParVectorCopy(par_rhs, Utemp);

            local_vector_ai = jxf_ParVectorLocalVector(Vtemp);
            local_data_ai   = jxf_VectorData(local_vector_ai);
            local_size   = jxf_VectorSize(local_vector_ai);

            local_vector_ess = jxf_ParVectorLocalVector(Utemp);
            local_data_ess   = jxf_VectorData(local_vector_ess);

            for (i=0; i<local_size; i++) {
                if (relax_marker_ai[i] == 0) {
                   local_data_ai[i] = 0.0;
                } else if (relax_marker_ess[i] == 0){
                   local_data_ess[i] = 0.0;
                } else {
                   jxf_printf("ERROR: relax_marker_ai does not match relax_marker_ess.");
                   exit(0);
                }
            }
            rhs_norm_ai = jxf_ParVectorNorm2(Vtemp);
            rhs_norm_ess = jxf_ParVectorNorm2(Utemp);

         } else {
            jxf_printf("NOT SUPPORT for conv_criteria != 0.");
            exit(0);
         }

         jxf_printf("++++++++++++++++++++++++++++++++++++++\n");
         jxf_printf("number of ai-points = %d\n", num_ai_th);
         jxf_printf("number of non-ai-points = %d\n", num_ess_th);
         jxf_printf("++++++++++++++++++++++++++++++++++++++\n");
         iter_ai = 100; 
         iter_ess = 100;
         //while ((relative >= tol || iter_ai < min_iter) && iter < max_iter && (iter_ai > 1 || iter_ess > 1))
         while ( (relative >= tol || iter < min_iter) && iter < max_iter && conv_factor <= convfac_threshold)
         {

            ++iter;
            resid_nrm_0 = resid_nrm;

            /*******************************************************
             * AI-prior smoothing process: divide into two steps. *
             ******************************************************/
            /* first step: AI-parts. */
            /* compute initial norm. */
            if (conv_criteria == 0)
            {
               // local resid_nrm.
               jxf_ParVectorCopy(par_rhs, Vtemp);
               jxf_hpCSRMatrixMatvec(alpha, hp_matrix, par_app, beta, Vtemp);
               jxf_ParVectorCopy(Vtemp, Utemp);

               local_vector_ai = jxf_ParVectorLocalVector(Vtemp);
               local_data_ai   = jxf_VectorData(local_vector_ai);
               local_size   = jxf_VectorSize(local_vector_ai);

               for (i=0; i<local_size; i++) {
                   if (relax_marker_ai[i] == 0) {
                      local_data_ai[i] = 0.0;
                   }
               }
               resid_nrm_ai = jxf_ParVectorNorm2(Vtemp);
               resid_nrm_init_ai = resid_nrm_ai;

               if (rhs_norm_ai)
               {
                  relative_ai = resid_nrm_init_ai / rhs_norm_ai;
               }
               else
               {
                  relative_ai = resid_nrm_init_ai;
               }

            } else {
               jxf_printf("NOT SUPPORT for conv_criteria != 0.");
               exit(0);
            }

            /* AI-parts: relaxation. */
            conv_factor = 0.0;
            iter_ai = 0;
            print_level_ai = 0;
            if (print_level_ai) jxf_printf("\n     Initial_ai    %e                %e\n", resid_nrm_init_ai, relative_ai);
            //while ( (relative >= tol || iter_ai < min_iter) && iter_ai < max_iter && num_ai_th>0 && conv_factor <= 0.9)
            while ( iter_ai < 10 && (relative_ai >= tol || iter_ai < min_iter) && num_ai_th>0 && conv_factor <= 0.6*convfac_threshold)
            //while ( num_ai_th>0 && conv_factor <= 0.95)
            {

               if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
               {
                  jxf_ParVectorCopy(par_app, Utemp);
               }

               jxf_CoarsestIterativeMethodAI( hp_matrix, par_rhs, relax_marker_ai, 
                                             solverid, weight, omega, par_app, Vtemp );

               ++ iter_ai;
               
               jxf_hp_check_convergence( relax_marker_ai, hp_matrix, par_rhs, alpha, beta, iter_ai, par_app, Vtemp, Utemp, 
                                     my_id, print_level_ai, conv_criteria, rhs_norm_ai, resid_nrm_init_ai,
                                     &resid_nrm_ai, &relative_ai, &conv_factor );
               //if (conv_factor > 1.0) exit(0);

            } // end for while.

            /* second step: ESS-parts. */
            /* compute initial norm. */
            if (conv_criteria == 0)
            {
               // local resid_nrm.
               jxf_ParVectorCopy(par_rhs, Utemp);
               jxf_hpCSRMatrixMatvec(alpha, hp_matrix, par_app, beta, Utemp);

               local_vector_ess = jxf_ParVectorLocalVector(Utemp);
               local_size   = jxf_VectorSize(local_vector_ess);
               local_data_ess   = jxf_VectorData(local_vector_ess);

               for (i=0; i<local_size; i++) {
                   if (relax_marker_ess[i] == 0){
                      local_data_ess[i] = 0.0;
                   }
               }
               resid_nrm_ess = jxf_ParVectorNorm2(Utemp);
               resid_nrm_init_ess = resid_nrm_ess;

               if (rhs_norm_ess)
               {
                  relative_ess = resid_nrm_init_ess / rhs_norm_ess;
               }
               else
               {
                  relative_ess = resid_nrm_init_ess;
               }

            } else {
               jxf_printf("NOT SUPPORT for conv_criteria != 0.");
               exit(0);
            }
 
            /* ESS-parts: relaxation. */
            iter_ess = 0;
            conv_factor = 0.0;  
            print_level_ess = 0;
            if (print_level_ess) jxf_printf("      Initial_ess    %e                %e\n", resid_nrm_init_ess, relative_ess);
            //while ( (relative >= tol || iter_ess < min_iter) && iter_ess < max_iter && num_ess_th>0 && conv_factor <= 0.9)
            while ( iter_ess < 10 && (relative_ess >= tol || iter_ess < min_iter) && num_ess_th>0 && conv_factor <= 0.6*convfac_threshold)
            //while ( num_ess_th>0 && conv_factor <= 0.95)
            {

               if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
               {
                  jxf_ParVectorCopy(par_app, Utemp);
               }

               ++ iter_ess;
               jxf_CoarsestIterativeMethodAI( hp_matrix, par_rhs, relax_marker_ess, 
                                             solverid, weight, omega, par_app, Vtemp );
               jxf_hp_check_convergence( relax_marker_ess, hp_matrix, par_rhs, alpha, beta, iter_ess, par_app, Vtemp, Utemp, 
                                     my_id, print_level_ess, conv_criteria, rhs_norm_ess, resid_nrm_init_ess, 
                                     &resid_nrm_ess, &relative_ess, &conv_factor );
               //if (conv_factor > 1.0) exit(0);

            } // end for while ess-parts.

            //jxf_printf("iter = %d, iter_ai = %d, iter_ess = %d\n", iter, iter_ai, iter_ess);
            jxf_hp_check_convergence( NULL, hp_matrix, par_rhs, alpha, beta, iter, par_app, Vtemp, Utemp, 
                                  my_id, print_level, conv_criteria, rhs_norm, resid_nrm_init,
                                  &resid_nrm_0, &relative, &conv_factor );
            resid_nrm = resid_nrm_0;
         }

      } else { 

         /******************************
          * general smoothing methods. *
          ******************************/
         while ( (relative >= tol || iter < min_iter) && iter < max_iter && conv_factor <= convfac_threshold)
         {

            if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
            {
               jxf_ParVectorCopy(par_app, Utemp);
            }

            ++ iter;
            jxf_CoarsestIterativeMethod( hp_matrix, par_rhs, solverid, weight, omega, par_app, Vtemp );

            //jxf_printf("resid_nrm = %f \n ", resid_nrm);
            jxf_hp_check_convergence( NULL, hp_matrix, par_rhs, alpha, beta, iter, par_app, Vtemp, Utemp, 
                                  my_id, print_level, conv_criteria, rhs_norm, resid_nrm_init,
                                  &resid_nrm, &relative, &conv_factor );
         } // end for while.

      }

      jxf_ParAMGDataRelativeResidualNorm(amg_data) = relative;
      jxf_ParAMGDataNumIterations(amg_data) = iter;  

     /*-----------------------------------------------------------------------
      *    Compute closing statistics
      *----------------------------------------------------------------------*/
      if (conv_criteria == 0)
      {
         //jxf_printf("iter = %d, tol = %e, resid_nrm_init= %e, resid_nrm= %e \n", iter, tol, resid_nrm_init, resid_nrm);
         if (iter > 0 && tol >= 0. && resid_nrm_init) 
         {
            conv_factor = pow( (resid_nrm / resid_nrm_init), (1.0 / (JXF_Real)iter) );
         }
         else
         {
            conv_factor = 1.;
         }
      }
      else
      {
         if (iter > 0 && tol >= 0.)
         {
            conv_factor = pow( relative, (1.0 / (JXF_Real)iter) );
         }
         else
         {
            conv_factor = 1.; 
         }       
      }
   
      if (my_id == 0 && iter == max_iter && print_level)
      {
         jxf_printf("\n=======================================================");
         jxf_printf("\n Warning: Convergence tolerance was not achieved     \n");
         jxf_printf("          within the allowed %d iterations!  \n", max_iter);
         jxf_printf("=======================================================\n");
      }

      if (my_id == 0 && tol >= 0. && print_level)
      {
        jxf_printf("\n Average Convergence Factor = %f\n\n", conv_factor); 
      }

      jxf_ParVectorDestroy(Utemp);
   }

   jxf_ParVectorDestroy(Vtemp);

   return 0;
}


/*!
 * \fn JXF_Int jxf_check_convergence
 * \brief check convergence and print iteration inforation.
 * \date 2013/10/18
 */
JXF_Int  
jxf_check_convergence( JXF_Int             *relax_marker,
                      jxf_ParCSRMatrix *par_matrix,
                      jxf_ParVector    *par_rhs,
                      JXF_Real           alpha,
                      JXF_Real           beta,
                      JXF_Int              iter,
                      jxf_ParVector    *par_app,
                      jxf_ParVector    *Vtemp,
                      jxf_ParVector    *Utemp,
                      JXF_Int              my_id,
                      JXF_Int              print_level,
                      JXF_Int              conv_criteria,
                      JXF_Real           rhs_norm,
                      JXF_Real           resid_nrm_init,
                      JXF_Real           *resid_nrm_old,
                      JXF_Real           *relative_ret,
                      JXF_Real           *conv_factor_old )
{

   jxf_Vector               *local_vector = NULL;
   JXF_Real                  *local_data = NULL;
   JXF_Int                      i, local_size;
   JXF_Real   normup         = 0.0;
   JXF_Real   normdown       = 0.0;
   JXF_Real   old_relative   = 0.0;
   JXF_Real   old_resid      = 0.0;
   JXF_Real   resid_nrm      = *resid_nrm_old;
   JXF_Real   conv_factor    = *conv_factor_old;
   JXF_Real   relative       = *relative_ret;

   if (relax_marker == NULL) {
       //jxf_printf("relax_marker == NULL.\n");
   } else {
       //jxf_printf("relax_marker != NULL.\n");
   }

   if (conv_criteria == 0)
   {
      old_resid = resid_nrm;

      jxf_ParVectorCopy(par_rhs, Vtemp);
      jxf_ParCSRMatrixMatvec(alpha, par_matrix, par_app, beta, Vtemp);

      //resid_nrm = jxf_ParVectorNorm2(Vtemp);
      //jxf_printf("resid_nrm_0 = %f \n", resid_nrm);
      if (relax_marker) {
         local_vector = jxf_ParVectorLocalVector(Vtemp);
         local_data   = jxf_VectorData(local_vector);
         local_size   = jxf_VectorSize(local_vector);
         for (i=0; i<local_size; i++) {
             if (relax_marker[i] == 0) local_data[i] = 0;
         }
      } 
      resid_nrm = jxf_ParVectorNorm2(Vtemp);
      //jxf_printf("resid_nrm = %f \n", resid_nrm);
      //jxf_printf("resid_nrm_old = %f \n", old_resid);
            
      if (old_resid) 
      {
         conv_factor = resid_nrm / old_resid;
      }
      else 
      {
         conv_factor = resid_nrm;
      }
            
      if (rhs_norm)
      {
         relative = resid_nrm / rhs_norm;
      }
      else
      {
         relative = resid_nrm / resid_nrm_init;
      }
   }

   if (conv_criteria == 1)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      normup   = jxf_ParVectorNorm1(Utemp);
      normdown = jxf_ParVectorNorm1(par_app);
      if (normdown < EPS)
      {
         relative = normup;
      }
      else
      {
         relative = normup / normdown;
      }
            
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   } 

   if (conv_criteria == 11)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      relative = jxf_ParVectorPointWiseRelNorm1(Utemp,par_app);
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   }

   if (conv_criteria == 2)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      normup   = jxf_ParVectorNorm2(Utemp);
      normdown = jxf_ParVectorNorm2(par_app);
      if (normdown < EPS)
      {
         relative = normup;
      }
      else
      {
         relative = normup / normdown;
      }
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   }

   *resid_nrm_old      = resid_nrm;
   *conv_factor_old    = conv_factor;
   *relative_ret       = relative;

   if (my_id == 0 && print_level)
   { 
      if (conv_criteria == 0)
      {
         jxf_printf("    ITER %3d   %e    %f    %e \n", iter, resid_nrm, conv_factor, relative);
      }
      if (conv_criteria == 1)
      {
         if (iter == 1)
         {
            jxf_printf("    ITER %3d   %e             %e \n", iter, normup, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %e   %6f  %e \n", iter, normup, conv_factor, relative);
         }
      }
      if (conv_criteria == 11)
      { 
         if (iter == 1)
         {
            jxf_printf("    ITER %3d               %e \n", iter, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %6f    %e \n", iter, conv_factor, relative);
         }
      }
      if (conv_criteria == 2)
      {
         if (iter == 1)
         {
            jxf_printf("    ITER %3d   %e             %e \n", iter, normup, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %e   %6f  %e \n", iter, normup, conv_factor, relative);
         }
      }
   }

   return 0;

}

JXF_Int  
jxf_hp_check_convergence( JXF_Int             *relax_marker,
                      jxf_hpCSRMatrix *par_matrix,
                      jxf_ParVector    *par_rhs,
                      JXF_Real           alpha,
                      JXF_Real           beta,
                      JXF_Int              iter,
                      jxf_ParVector    *par_app,
                      jxf_ParVector    *Vtemp,
                      jxf_ParVector    *Utemp,
                      JXF_Int              my_id,
                      JXF_Int              print_level,
                      JXF_Int              conv_criteria,
                      JXF_Real           rhs_norm,
                      JXF_Real           resid_nrm_init,
                      JXF_Real           *resid_nrm_old,
                      JXF_Real           *relative_ret,
                      JXF_Real           *conv_factor_old)
{

   jxf_Vector               *local_vector = NULL;
   JXF_Real                  *local_data = NULL;
   JXF_Int                      i, local_size;
   JXF_Real   normup         = 0.0;
   JXF_Real   normdown       = 0.0;
   JXF_Real   old_relative   = 0.0;
   JXF_Real   old_resid      = 0.0;
   JXF_Real   resid_nrm      = *resid_nrm_old;
   JXF_Real   conv_factor    = *conv_factor_old;
   JXF_Real   relative       = *relative_ret;

   if (relax_marker == NULL) {
       //jxf_printf("relax_marker == NULL.\n");
   } else {
       //jxf_printf("relax_marker != NULL.\n");
   }

   if (conv_criteria == 0)
   {
      old_resid = resid_nrm;

      jxf_ParVectorCopy(par_rhs, Vtemp);
      jxf_hpCSRMatrixMatvec(alpha, par_matrix, par_app, beta, Vtemp);

      //resid_nrm = jxf_ParVectorNorm2(Vtemp);
      //jxf_printf("resid_nrm_0 = %f \n", resid_nrm);
      if (relax_marker) {
         local_vector = jxf_ParVectorLocalVector(Vtemp);
         local_data   = jxf_VectorData(local_vector);
         local_size   = jxf_VectorSize(local_vector);
         for (i=0; i<local_size; i++) {
             if (relax_marker[i] == 0) local_data[i] = 0;
         }
      } 
      resid_nrm = jxf_ParVectorNorm2(Vtemp);
      //jxf_printf("resid_nrm = %f \n", resid_nrm);
      //jxf_printf("resid_nrm_old = %f \n", old_resid);
            
      if (old_resid) 
      {
         conv_factor = resid_nrm / old_resid;
      }
      else 
      {
         conv_factor = resid_nrm;
      }
            
      if (rhs_norm)
      {
         relative = resid_nrm / rhs_norm;
      }
      else
      {
         relative = resid_nrm / resid_nrm_init;
      }
   }

   if (conv_criteria == 1)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      normup   = jxf_ParVectorNorm1(Utemp);
      normdown = jxf_ParVectorNorm1(par_app);
      if (normdown < EPS)
      {
         relative = normup;
      }
      else
      {
         relative = normup / normdown;
      }
            
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   } 

   if (conv_criteria == 11)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      relative = jxf_ParVectorPointWiseRelNorm1(Utemp,par_app);
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   }

   if (conv_criteria == 2)
   {
      if (iter-1 > 0) 
      {
         old_relative = relative;
      }
      jxf_ParVectorAxpy(-1.0, par_app, Utemp);  // Utemp = Utemp - par_app
      normup   = jxf_ParVectorNorm2(Utemp);
      normdown = jxf_ParVectorNorm2(par_app);
      if (normdown < EPS)
      {
         relative = normup;
      }
      else
      {
         relative = normup / normdown;
      }
      if (iter-1 > 0) 
      {
         conv_factor = relative / old_relative;
      }
   }

   *resid_nrm_old      = resid_nrm;
   *conv_factor_old    = conv_factor;
   *relative_ret       = relative;

   if (my_id == 0 && print_level)
   { 
      if (conv_criteria == 0)
      {
         jxf_printf("    ITER %3d   %e    %f    %e \n", iter, resid_nrm, conv_factor, relative);
      }
      if (conv_criteria == 1)
      {
         if (iter == 1)
         {
            jxf_printf("    ITER %3d   %e             %e \n", iter, normup, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %e   %6f  %e \n", iter, normup, conv_factor, relative);
         }
      }
      if (conv_criteria == 11)
      { 
         if (iter == 1)
         {
            jxf_printf("    ITER %3d               %e \n", iter, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %6f    %e \n", iter, conv_factor, relative);
         }
      }
      if (conv_criteria == 2)
      {
         if (iter == 1)
         {
            jxf_printf("    ITER %3d   %e             %e \n", iter, normup, relative);
         }
         else
         {
            jxf_printf("    ITER %3d   %e   %6f  %e \n", iter, normup, conv_factor, relative);
         }
      }
   }

   return 0;

}

/*!
 * \fn JXF_Int jxf_CoarsestSolverInfo
 * \brief Print information for the coarsest solvers.
 * \author peghoty
 * \date 2009/07/25
 */
JXF_Int
jxf_CoarsestSolverInfo( MPI_Comm comm, void *amg_vdata )
{
   //MPI_Comm        comm     = MPI_COMM_WORLD;
   jxf_ParAMGData  *amg_data = amg_vdata;

   JXF_Real tol           = jxf_ParAMGDataTol(amg_data);
   JXF_Real weight        = (jxf_ParAMGDataRelaxWeight(amg_data))[0];
   JXF_Real omega         = (jxf_ParAMGDataOmega(amg_data))[0];
   JXF_Int    solverid      = jxf_ParAMGDataCoarsestSolverID(amg_data);
   JXF_Int    max_iter      = jxf_ParAMGDataMaxIter(amg_data); 
   JXF_Int    min_iter      = jxf_ParAMGDataMinIter(amg_data);  
   JXF_Int    conv_criteria = jxf_ParAMGDataConvCriteria(amg_data);
   JXF_Int    my_id;

   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == 0)
   {
     jxf_printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"); 
     jxf_printf("Using iterative method since max_levels = 1. \n"); 
     jxf_printf("---------------------------------------------\n"); 
     jxf_printf("Coarsest Solver INFO:\n");  
     jxf_printf("---------------------\n"); 
     jxf_printf("Solverid : %d\n", solverid);    
     jxf_printf("Method   : ");   
     switch (solverid)
     {
        case 0:
        jxf_printf("Weighted Jacobi\n");
        break;
        
        case 1:
        jxf_printf("Gauss-Seidel\n");
        jxf_printf("Remark   : very slow, sequential\n");
        break;
        
        case 2:
        jxf_printf("Gauss_Seidel\n");
        jxf_printf("Remark   : interior points in parallel,boundary sequential\n");
        break;
        
        case 3:
        jxf_printf("hybrid\n");
        jxf_printf("Remark   : SOR-J mix off-processor\n");
        jxf_printf("           SOR on-processor with outer relaxation parameters\n");
        jxf_printf("           (forward solve)\n");
        break;
        
        case 4:
        jxf_printf("hybrid\n");
        jxf_printf("Remark   : SOR-J mix off-processor\n");
        jxf_printf("           SOR on-processor with outer relaxation parameters\n");
        jxf_printf("           (backward solve)\n");
        break;
        
        case 5:
        jxf_printf("hybrid\n");
        jxf_printf("Remark   : GS-J mix off-processor, chaotic GS on-node\n");
        break;
        
        case 6:
        jxf_printf("hybrid\n");
        jxf_printf("Remark   : SSOR-J mix off-processor\n");
        jxf_printf("           SSOR on-processor with outer relaxation parameters\n");
        break;
        
        case 9:
        jxf_printf("Gaussian Elimination (default)\n");
        break;
        
        case 10:
        jxf_printf("Gaussian Elimination - with pivoting\n");
        break;
     }
      
     if (solverid != 9 && solverid != 10)
     {
        jxf_printf("tol      : %le\n", tol); 
        jxf_printf("weight   : %.1f\n", weight); 
        jxf_printf("omega    : %.1f\n", omega);
        jxf_printf("max_iter : %d\n", max_iter); 
        jxf_printf("min_iter : %d\n", min_iter); 
        jxf_printf("Convergence Criteria :"); 
        switch (conv_criteria)
        {
           case 0:
           jxf_printf(" Relative Residual Norm-2 (default)\n");
           break;
           
           case 1:
           jxf_printf(" Relative Error Norm-1\n");
           break;
           
           case 11:
           jxf_printf(" Relative Point-wise Error Norm-1\n");
           break;
           
           case 2:
           jxf_printf(" Relative Error Norm-2\n");
           break;
        }
     }
     jxf_printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");      
   }

   return 0;
}
