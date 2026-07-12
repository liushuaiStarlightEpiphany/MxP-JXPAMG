//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_amg_solver.c -- SOLUTION phase of PAMG.
 *  Date: 2011/09/27
 */ 
 
#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGSolve
 * \brief Solution phase of PAMG.
 * \param amg_vdata pointer to the amgdata object
 * \param par_matrix pointer to the coefficient matrix
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector 
 * \remark (peghoty  2009/07/24)
 *   Compared with the hypre code, there are two main alterations in 
 *   this function as follows:
 *     1. dealing with the cases "max_levels = 1" and "max_levels > 1" separately;
 *     2. add a parameter "conv_criteria", the user can choose a proper convergence 
 *        control by providing the value of the parameter
 *          conv_criteria =  0:  ||r^{k}||_2 / ||b||_2
 *          conv_criteria =  1:  ||u^{k} - u^{k-1}||_1 / ||u^{k}||_1
 *          conv_criteria = 11:  max { |u_i^{k} - u_i^{k-1}| / |u_i^{k}| }
 *          conv_criteria =  2:  ||u^{k} - u^{k-1}||_2 / ||u^{k}||_2
 *   where 
 *         ||u||_1 = max{|u_i|}， ||u||_2 = sqrt((u_1)^2 + ... + (u_n)^2)
 * \date 2011/09/28
 */

JXF_Int
jxf_hpPAMGSolve( void            *amg_vdata,
              jxf_hpCSRMatrix *hp_matrix,
              jxf_ParVector    *par_rhs,
              jxf_ParVector    *par_app)
{
   //MPI_Comm        comm     = MPI_COMM_WORLD;
   MPI_Comm         comm     = jxf_hpCSRMatrixComm(hp_matrix); 
   jxf_ParAMGData  *amg_data = amg_vdata;
   /* Data Structure variables */
   JXF_Int    amg_print_level;
   JXF_Int    amg_logging;
   JXF_Int    cycle_count;
   JXF_Int    num_levels;
   JXF_Int    max_levels;
   JXF_Real tol;
   JXF_Real rhsnrm_threshold;

   jxf_hpCSRMatrix  **A_array;
   jxf_ParVector    **F_array;
   jxf_ParVector    **U_array;

   /*  Local variables  */
   JXF_Int      j;
   JXF_Int      Solve_err_flag;
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      num_procs, my_id;

   JXF_Real   alpha = 1.0;
   JXF_Real   beta = -1.0;
   JXF_Real   cycle_op_count;
   JXF_Real   total_coeffs;
   JXF_Real   total_variables;
   JXF_Real  *num_coeffs;
   JXF_Real  *num_variables;
   JXF_Real   cycle_cmplxty  = 0.0;
   JXF_Real   operat_cmplxty = 0.0;
   JXF_Real   grid_cmplxty   = 0.0;
   JXF_Real   conv_factor    = 0.0;
   JXF_Real   resid_nrm      = 0.0;
   JXF_Real   resid_nrm_init = 0.0;
   JXF_Real   relative       = 0.0;
   JXF_Real   old_relative   = 0.0;
   JXF_Real   rhs_norm       = 0.0;
   JXF_Real   old_resid      = 0.0;
   JXF_Real   ieee_check     = 0.0;

   jxf_ParVector  *Vtemp    = NULL;
   jxf_ParVector  *Residual = NULL;
   jxf_ParVector  *Utemp    = NULL;

   JXF_Real  normup   = 0.0;
   JXF_Real  normdown = 0.0;
   JXF_Int    *partitioning;
   JXF_Int     size;
   
   /* peghoty 2009/07/25 */
   JXF_Int conv_criteria;
   JXF_Real convfac_threshold;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   /* peghoty 2009/07/27 */
   max_levels = jxf_ParAMGDataMaxLevels(amg_data);
   num_levels = jxf_ParAMGDataNumLevels(amg_data);
   if (max_levels == 1)
   {
      jxf_CoarsestSolver( amg_data, hp_matrix, par_rhs, par_app );
      return(0);
   }
   
   conv_criteria    = jxf_ParAMGDataConvCriteria(amg_data);
   convfac_threshold = jxf_ParAMGDataConvFacThreshold(amg_data);
   amg_print_level  = jxf_ParAMGDataPrintLevel(amg_data);
   amg_logging      = jxf_ParAMGDataLogging(amg_data);
   tol              = jxf_ParAMGDataTol(amg_data);
   rhsnrm_threshold = jxf_ParAMGDataRhsNrmThreshold(amg_data);
   min_iter         = jxf_ParAMGDataMinIter(amg_data);
   max_iter         = jxf_ParAMGDataMaxIter(amg_data);
   Vtemp            = jxf_ParAMGDataVtemp(amg_data);

#if 0 // not robust. peghoty, 2012/02/23
   size = jxf_ParVectorGlobalSize(par_rhs);
   partitioning = jxf_ParVectorPartitioning(par_rhs);
   Utemp = jxf_ParVectorCreate(comm, size, partitioning);
   jxf_ParVectorSetPartitioningOwner(Utemp, 0);
   jxf_ParVectorInitialize(Utemp);
#else
   size = jxf_hpCSRMatrixGlobalNumRows(hp_matrix);
   partitioning = jxf_hpCSRMatrixRowStarts(hp_matrix);
   Utemp = jxf_ParVectorCreate(comm, size, partitioning);
   jxf_ParVectorInitialize(Utemp);
   jxf_ParVectorSetPartitioningOwner(Utemp, 0); 
#endif

   /* Get the amg_data data */
   A_array = jxf_hpAMGDataAArray(amg_data);
   F_array = jxf_ParAMGDataFArray(amg_data);
   U_array = jxf_ParAMGDataUArray(amg_data);

   /* they are moved from "amg_setup".  peghoty, 2009/07/27 */
   F_array[0] = par_rhs;
   U_array[0] = par_app;

   if (amg_logging > 1)
   {
      Residual = jxf_ParAMGDataResidual(amg_data);
   }
   
   num_coeffs       = jxf_CTAlloc(JXF_Real, num_levels);
   num_variables    = jxf_CTAlloc(JXF_Real, num_levels);
   num_coeffs[0]    = (JXF_Real) jxf_hpCSRMatrixNumNonzeros(A_array[0]);
   num_variables[0] = (JXF_Real) jxf_hpCSRMatrixGlobalNumRows(A_array[0]);

   for (j = 1; j < num_levels; j ++)
   {
      num_coeffs[j]    = (JXF_Real) jxf_hpCSRMatrixNumNonzeros(A_array[j]);
      num_variables[j] = (JXF_Real) jxf_hpCSRMatrixGlobalNumRows(A_array[j]);
   }

   /* Write the solver parameters */ // set the solve status into SETUP phase, Yue Xiaoqiang, 2014/04/12
   // if (my_id == 0 && amg_print_level > 1)
   // { 
   //    jxf_PAMGSolveStatus(amg_data);
   // }

   /* Initialize the solver error flag and assorted bookkeeping variables */
   Solve_err_flag  = 0;
   total_coeffs    = 0;
   total_variables = 0;
   cycle_count     = 0;
   operat_cmplxty  = 0;
   grid_cmplxty    = 0;

   /* write some initial info */
   if (my_id == 0 && amg_print_level > 1 && tol > 0.)
   {
      jxf_printf("\n\nJXFPAMG SOLUTION INFO:\n");
   }


   //-----------------------------------------------------------------------//
   //       Compute initial fine-grid norm and print the table header       //
   //-----------------------------------------------------------------------//
   
   if (conv_criteria == 0)
   {
      if (tol >= 0.)
      {
          if (amg_logging > 1) 
          {
             jxf_ParVectorCopy(F_array[0], Residual);
             jxf_hpCSRMatrixMatvec(alpha, A_array[0], U_array[0], beta, Residual);
             resid_nrm = jxf_ParVectorNorm2(Residual);
          }
          else 
          {
             jxf_ParVectorCopy(F_array[0], Vtemp);
             jxf_hpCSRMatrixMatvec(alpha, A_array[0], U_array[0], beta, Vtemp);
             resid_nrm = jxf_ParVectorNorm2(Vtemp);
          }
         /*--------------------------------------------------------------------------------------
          * Since it does not diminish performance, attempt to return an error flag
          * and notify users when they supply bad input. 
          *-------------------------------------------------------------------------------------*/
          if (resid_nrm != 0.) ieee_check = resid_nrm / resid_nrm; /* INF -> NaN conversion */
          
          if (ieee_check != ieee_check)
          {
            /*----------------------------------------------------------------------------
             * ...INFs or NaNs in input can make ieee_check a NaN.  This test
             *   for ieee_check self-equality works on all IEEE-compliant compilers/
             *   machines, c.f. page 8 of "Lecture Notes on the Status of IEEE 754"
             *   by W. Kahan, May 31, 1996.  Currently (July 2002) this paper may be
             *   found at http://HTTP.CS.Berkeley.EDU/~wkahan/ieee754status/IEEE754.PDF 
             *---------------------------------------------------------------------------*/
             if (amg_print_level > 0)
             {
                jxf_printf("\n\nERROR detected by JXFPAMG ...  BEGIN\n");
                jxf_printf("ERROR -- jxf_PAMGSolve: INFs and/or NaNs detected in input.\n");
                jxf_printf("User probably placed non-numerics in supplied A, x_0, or b.\n");
                jxf_printf("ERROR detected by JXFPAMG ...  END\n\n\n");
//                break;
             }
             jxf_error(JXF_ERROR_GENERIC);
             return jxf_error_flag;
          }

          resid_nrm_init = resid_nrm;
          rhs_norm = jxf_ParVectorNorm2(F_array[0]);
          if (rhsnrm_threshold) /* peghoty, 2011/09/28 */
          {
             if (rhs_norm > rhsnrm_threshold)
             {
                relative = resid_nrm_init / rhs_norm;
             }
             else
             {
                relative = resid_nrm_init;
             }
          }
          else
          {
             if (rhs_norm)
             {
                relative = resid_nrm_init / rhs_norm;
             }
             else
             {
                relative = resid_nrm_init;
             }
          }
      }
      else
      {
         relative = 1.0;
      }
   }
   
   if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
   {
     /*---------------------------------------------------------------------------
      * To guarantee entering the "while loop" once at least peghoty, 2009/07/24 
      *--------------------------------------------------------------------------*/
      relative = 2*tol;  
   }

   if (my_id == 0 && amg_print_level > 1 && tol >= 0.)
   {    
      if (conv_criteria == 0)
      { 
         jxf_printf("                                             Relative  \n");
         jxf_printf("                Res Norm-2       Factor      Res Norm-2\n");
         jxf_printf("                ----------      -------     -----------\n");
         jxf_printf("    Initial     %e                \n",resid_nrm_init);
         //jxf_printf("    Initial     %e                %e\n",resid_nrm_init,relative);
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


   //-----------------------------------------------------------------------//
   //                    M A I N   V  C Y C L E   L O O P                   //
   //-----------------------------------------------------------------------//

//   jxf_ParAMGDataCycleOpCount(amg_data) = 0; 
   while ( (relative >= tol || cycle_count < min_iter) && cycle_count < max_iter )
   {
      jxf_ParAMGDataCycleOpCount(amg_data) = 0;
      /* Op count only needed for one cycle */

      if (conv_criteria == 1 || conv_criteria == 11 || conv_criteria == 2)
      {
         jxf_ParVectorCopy(U_array[0], Utemp);
      }


     /*----------------------------------------------
      *  Excute One Cycle
      *---------------------------------------------*/
      jxf_hpPAMGCycle(amg_data, F_array, U_array);  

     /*------------------------------------------------------------------
      *  Compute the residual norm on fine-grid  peghoty, 2009/07/23
      *-----------------------------------------------------------------*/
     
      if (conv_criteria == 0)
      {
         if (tol >= 0.)
         {
            old_resid = resid_nrm;

            if (amg_logging > 1) 
            {
               jxf_ParVectorCopy(F_array[0], Residual);
               jxf_hpCSRMatrixMatvec(alpha, A_array[0], U_array[0], beta, Residual);
               resid_nrm = jxf_ParVectorNorm2(Residual);
            }
            else 
            {
               jxf_ParVectorCopy(F_array[0], Vtemp);
               jxf_hpCSRMatrixMatvec(alpha, A_array[0], U_array[0], beta, Vtemp);
               resid_nrm = jxf_ParVectorNorm2(Vtemp);
            }

            if (old_resid)
            { 
               conv_factor = resid_nrm / old_resid;
            }   
            else
            { 
               conv_factor = resid_nrm;
            }
            
            if (rhsnrm_threshold) /* peghoty, 2011/09/28 */
            {
               if (rhs_norm > rhsnrm_threshold)
               {
                  relative = resid_nrm / rhs_norm;
               }
               else
               {
                  relative = resid_nrm / resid_nrm_init;
               }
            }
            else
            {
               if (rhs_norm)
               {
                  relative = resid_nrm / rhs_norm;
               }
               else
               {
                  relative = resid_nrm / resid_nrm_init;
               }
            }
         }
      } // end if (conv_criteria == 0)
      
      if (conv_criteria == 1)
      {
         if (cycle_count > 0) old_relative = relative;
         jxf_ParVectorAxpy(-1.0, U_array[0], Utemp);  // Utemp = Utemp - U_array[0]
         normup   = jxf_ParVectorNorm1(Utemp);
         normdown = jxf_ParVectorNorm1(U_array[0]);
         if (normdown < EPS)
         {
            relative = normup;
         }
         else
         {
            relative = normup / normdown;
         }
         if (cycle_count > 0) 
         {
            conv_factor = relative / old_relative;
         }
      } // end if (conv_criteria == 1)

      if (conv_criteria == 11)
      {
         if (cycle_count > 0) 
         {
            old_relative = relative;
         }
         jxf_ParVectorAxpy(-1.0, U_array[0], Utemp);  // Utemp = Utemp - U_array[0]
         relative = jxf_ParVectorPointWiseRelNorm1(Utemp,U_array[0]);
         if (cycle_count > 0) 
         {
            conv_factor = relative / old_relative;
         }
      } // end if (conv_criteria == 11)

      if(conv_criteria == 2)
      {
         if (cycle_count > 0) old_relative = relative;
         jxf_ParVectorAxpy(-1.0, U_array[0], Utemp);  // Utemp = Utemp - U_array[0]
         normup   = jxf_ParVectorNorm2(Utemp);
         normdown = jxf_ParVectorNorm2(U_array[0]);
         if (normdown < EPS)
         {
            relative = normup;
         }
         else
         {
            relative = normup / normdown;
         }
         if (cycle_count > 0) 
         {
            conv_factor = relative / old_relative;
         }
      } // end if (conv_criteria == 2)

      ++ cycle_count;

#ifdef CUMNUMIT
      ++ jxf_ParAMGDataCumNumIterations(amg_data);
#endif

      if (my_id == 0 && amg_print_level > 1 && tol >= 0.)
      { 
         if (conv_criteria == 0)
         {
            jxf_printf("    Cycle %3d   %e    %6f    %e \n", cycle_count, resid_nrm, conv_factor, relative);
         }
         else if (conv_criteria == 1)
         {
            if (cycle_count == 1)
            {
               jxf_printf("    Cycle %3d   %e             %e \n", cycle_count, normup, relative);
            }
            else
            {
               jxf_printf("    Cycle %3d   %e   %6f  %e \n", cycle_count, normup, conv_factor, relative);
            }
         }
         else if (conv_criteria == 11)
         {
            if (cycle_count == 1)
            {
               jxf_printf("    Cycle %3d               %e \n", cycle_count, relative);
            }
            else
            {
               jxf_printf("    Cycle %3d   %6f    %e \n", cycle_count, conv_factor, relative);
            }
         }
         else if (conv_criteria == 2)
         {
            if (cycle_count == 1)
            {
               jxf_printf("    Cycle %3d   %e             %e \n", cycle_count, normup, relative);
            }
            else
            {
               jxf_printf("    Cycle %3d   %e   %6f  %e \n", cycle_count, normup, conv_factor, relative);
            }
         }
      }

      if (conv_factor > convfac_threshold)
      {
         if (my_id == 0 && amg_print_level > 1)
         {
            jxf_printf("\n  Warning: iteration has terminated because the\n");
            jxf_printf("           current conv_factor > convfac_threshold(%lf) !!\n\n",convfac_threshold);
         }
         break;
      }

   }  
   jxf_ParAMGDataRelativeResidualNorm(amg_data) = relative;
   jxf_ParAMGDataNumIterations(amg_data) = cycle_count;

   if (cycle_count == max_iter && tol > 0.) 
   {
      Solve_err_flag = 1;
   }


   //-----------------------------------------------------------------------//
   //                    Compute closing statistics                         //
   //-----------------------------------------------------------------------//    
    
   if (conv_criteria == 0)
   {
      if (cycle_count > 0 && tol >= 0. && resid_nrm_init) 
      {
         conv_factor = pow((resid_nrm/resid_nrm_init),(1.0/(JXF_Real) cycle_count));
      }
      else
      {
         conv_factor = 1.0;
      }
   }
   else
   {
      if (cycle_count > 0 && tol >= 0)
      {
         conv_factor = pow(relative,(1.0/(JXF_Real) cycle_count));
      }
      else
      {
         conv_factor = 1.;
      }        
   }

   for (j = 0; j < num_levels; j ++)
   {
      total_coeffs    += num_coeffs[j];
      total_variables += num_variables[j];
   }

   cycle_op_count = jxf_ParAMGDataCycleOpCount(amg_data);

   if (num_variables[0])
   {
      grid_cmplxty = total_variables / num_variables[0];
   }
   
   if (num_coeffs[0])
   {
      operat_cmplxty = total_coeffs / num_coeffs[0];
      cycle_cmplxty  = cycle_op_count / num_coeffs[0];
   }

   if (my_id == 0 && amg_print_level > 1)
   {
      if (Solve_err_flag == 1)
      {
         jxf_printf("\n==================================================");
         jxf_printf("\n Warning: Convergence tolerance was not achieved\n");
         jxf_printf("          within the allowed %d V-cycles!\n",max_iter);
         jxf_printf("==================================================\n");
      }
      if (tol >= 0.)
      {
         jxf_printf("\n\n Average Convergence Factor = %f",conv_factor);
      }
      jxf_printf("\n\n     Complexity:    grid = %f\n",grid_cmplxty);
      jxf_printf("                operator = %f\n",  operat_cmplxty);
      jxf_printf("                   cycle = %f\n\n", cycle_cmplxty);
   }

   jxf_TFree(num_coeffs);
   jxf_TFree(num_variables);
   jxf_ParVectorDestroy(Utemp);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGPrecond
 * \brief PAMG as a preconditioner.
 * \param amg_vdata pointer to the amgdata object
 * \param par_matrix pointer to the coefficient matrix 
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \author peghoty 
 * \date 2011/09/28
 */
JXF_Int
jxf_PAMGPrecond( void            *amg_vdata,
                jxf_hpCSRMatrix *hp_matrix,
                jxf_ParVector    *par_rhs,
                jxf_ParVector    *par_app)
{
   jxf_ParAMGData  *amg_data = amg_vdata;
   jxf_ParVector  **F_array = jxf_ParAMGDataFArray(amg_data);
   jxf_ParVector  **U_array = jxf_ParAMGDataUArray(amg_data);   

   JXF_Int max_iter    = jxf_ParAMGDataMaxIter(amg_data); 
   JXF_Int max_levels  = jxf_ParAMGDataMaxLevels(amg_data);  
   JXF_Int cycle_count = 0;

   if (max_levels == 1)
   {
      jxf_CoarsestSolver( amg_data, hp_matrix, par_rhs, par_app );
      return(0);
   }    
   F_array[0] = par_rhs;
   U_array[0] = par_app;


  /*----------------------------------------------------------------
   * Generally, max_iter should be set to 1 when AMG is used as  
   * a preconditioner, here, we keep the parameter max_iter for  
   * convenience if max_iter should be larger than 1.
   *---------------------------------------------------------------*/
   while (cycle_count < max_iter)
   {
      jxf_hpPAMGCycle(amg_data, F_array, U_array); 
      cycle_count ++;
   }
   return 0; 
}        
