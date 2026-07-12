//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_status.c -- Information outputing of the Setup and Solution 
 *                  phase of PAMG.
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGSetupStatus
 * \brief Information outputing of the Setup phase of PAMG.
 * \date 2011/09/03
 */ 
JXF_Int
jxf_hpPAMGSetupStatus( void            *amg_vdata,
                    jxf_hpCSRMatrix *hp_matrix )
{
   MPI_Comm       comm     = jxf_hpCSRMatrixComm(hp_matrix);   
   jxf_ParAMGData *amg_data = amg_vdata;

   /* Data Structure variables */

   jxf_hpCSRMatrix  **A_array;
   jxf_ParCSRMatrix  **P_array;
   jxf_ParCSRMatrix  **R_array;

   jxf_CSRMatrix *A_diag;
   JXF_Real       *A_diag_data;
   JXF_Int          *A_diag_i;

   jxf_CSRMatrix *A_offd;   
   JXF_Real       *A_offd_data;
   JXF_Int          *A_offd_i;

   jxf_CSRMatrix *P_diag;
   JXF_Real       *P_diag_data;
   JXF_Int          *P_diag_i;

   jxf_CSRMatrix *P_offd;   
   JXF_Real       *P_offd_data;
   JXF_Int          *P_offd_i;

   JXF_Int	    numrows;
   JXF_Int	   *row_starts;
   JXF_Int      num_levels; 
   JXF_Int      coarsen_type;
   JXF_Int      interp_type;
   JXF_Int      restri_type;
   JXF_Int      measure_type;
   JXF_Int      agg_interp_type;
   JXF_Int      agg_num_levels;
   JXF_Real   global_nonzeros;

   JXF_Real  *send_buff;
   JXF_Real  *gather_buff;
 
   /* Local variables */
   JXF_Int      level;
   JXF_Int      i,j;
   JXF_Int      fine_size;
 
   JXF_Int      min_entries;
   JXF_Int      max_entries;

   JXF_Int      num_procs,my_id;

   JXF_Real   min_rowsum;
   JXF_Real   max_rowsum;
   JXF_Real   sparse;

   JXF_Int      coarse_size;
   JXF_Int      entries;

   JXF_Real   avg_entries;
   JXF_Real   rowsum;

   JXF_Real   min_weight;
   JXF_Real   max_weight;

   JXF_Int      global_min_e;
   JXF_Int      global_max_e;
   JXF_Real   global_min_rsum;
   JXF_Real   global_max_rsum;
   JXF_Real   global_min_wt;
   JXF_Real   global_max_wt;

   JXF_Real  *num_coeffs;
   JXF_Real  *num_variables;
   JXF_Real   total_variables; 
   JXF_Real   operat_cmplxty = 0.0;
   JXF_Real   grid_cmplxty   = 0.0;

   /* amg solve params */
   JXF_Int      max_iter;
   JXF_Int      cycle_type;    
   JXF_Int     *num_grid_sweeps;  
   JXF_Int     *grid_relax_type;   
   JXF_Int      relax_order;
   JXF_Int    **grid_relax_points; 
   JXF_Real  *relax_weight;
   JXF_Real  *omega;
   JXF_Real   tol;
   JXF_Int      conv_criteria; /* peghoty 2009/07/27 */
   //JXF_Int      amg_print_level;
   JXF_Int     *AIR_maxsize_ls;

   JXF_Int one       =  1;
   JXF_Int minus_one = -1;
   JXF_Int zero      =  0;
 
   jxf_MPI_Comm_size(comm, &num_procs);   
   jxf_MPI_Comm_rank(comm, &my_id);

   A_array = jxf_hpAMGDataAArray(amg_data);
   P_array = jxf_ParAMGDataPArray(amg_data);
   R_array = jxf_ParAMGDataRArray(amg_data);
   num_levels = jxf_ParAMGDataNumLevels(amg_data);
   coarsen_type = jxf_ParAMGDataCoarsenType(amg_data);
   interp_type = jxf_ParAMGDataInterpType(amg_data);
   restri_type = jxf_ParAMGDataRestriction(amg_data);
   measure_type = jxf_ParAMGDataMeasureType(amg_data);
   agg_num_levels = jxf_ParAMGDataAggNumLevels(amg_data);
   agg_interp_type = jxf_ParAMGDataAggInterpType(amg_data);
   AIR_maxsize_ls = jxf_ParAMGDataAIRMaxSizeLS(amg_data);

  /*----------------------------------------------------------
   * Get the amg_data data
   *--------------------------------------------------------*/

   num_levels        = jxf_ParAMGDataNumLevels(amg_data);
   max_iter          = jxf_ParAMGDataMaxIter(amg_data);
   cycle_type        = jxf_ParAMGDataCycleType(amg_data);    
   num_grid_sweeps   = jxf_ParAMGDataNumGridSweeps(amg_data);  
   grid_relax_type   = jxf_ParAMGDataGridRelaxType(amg_data);
   grid_relax_points = jxf_ParAMGDataGridRelaxPoints(amg_data);
   relax_weight      = jxf_ParAMGDataRelaxWeight(amg_data); 
   relax_order       = jxf_ParAMGDataRelaxOrder(amg_data); 
   omega             = jxf_ParAMGDataOmega(amg_data); 
   tol               = jxf_ParAMGDataTol(amg_data);
   conv_criteria     = jxf_ParAMGDataConvCriteria(amg_data); /* peghoty 2009/07/27 */
   //amg_print_level   = jxf_ParAMGDataPrintLevel(amg_data);

   send_buff   = jxf_CTAlloc(JXF_Real, 6);
#ifdef JXF_NO_GLOBAL_PARTITION
   gather_buff = jxf_CTAlloc(JXF_Real, 6);    
#else
   gather_buff = jxf_CTAlloc(JXF_Real, 6*num_procs);    
#endif

   if (my_id == 0)
   {
      jxf_printf("\njxf_PAMG SETUP PARAMETERS:\n\n");
      jxf_printf(" Max levels = %d\n",jxf_ParAMGDataMaxLevels(amg_data));
      jxf_printf(" Num levels = %d\n\n",num_levels);
      jxf_printf(" Strength Threshold = %f\n", jxf_ParAMGDataStrongThreshold(amg_data));
      if (restri_type > 0)
      {
     jxf_printf(" AIR Strength Threshold = %f for restriction\n", jxf_ParAMGDataAIRStrongTh(amg_data));
      }
      jxf_printf(" Interpolation Truncation Factor = %f\n", jxf_ParAMGDataTruncFactor(amg_data));
      jxf_printf(" Maximum Row Sum Threshold for Dependency Weakening = %f\n\n", jxf_ParAMGDataMaxRowSum(amg_data));

      if (coarsen_type == 0)
      {
	 jxf_printf(" Coarsening Type = Cleary-Luby-Jones-Plassman\n");
      }
      else if (abs(coarsen_type) == 1) 
      {
	 jxf_printf(" Coarsening Type = Ruge\n");
      }
      else if (abs(coarsen_type) == 2) 
      {
	 jxf_printf(" Coarsening Type = Ruge2B\n");
      }
      else if (abs(coarsen_type) == 3) 
      {
	 jxf_printf(" Coarsening Type = Ruge3\n");
      }
      else if (abs(coarsen_type) == 4) 
      {
	 jxf_printf(" Coarsening Type = Ruge 3c \n");
      }
      else if (abs(coarsen_type) == 5) 
      {
	 jxf_printf(" Coarsening Type = Ruge relax special points \n");
      }
      else if (abs(coarsen_type) == 6) 
      {
	 jxf_printf(" Coarsening Type = Falgout-CLJP \n");
      }
      else if (abs(coarsen_type) == 8) 
      {
	 jxf_printf(" Coarsening Type = PMIS \n");
      }
      else if (abs(coarsen_type) == 10) 
      {
	 jxf_printf(" Coarsening Type = HMIS \n");
      }
      else if (abs(coarsen_type) == 11) 
      {
	 jxf_printf(" Coarsening Type = Ruge 1st pass only \n");
      }
      else if (coarsen_type == 990)
      {
	 jxf_printf(" Coarsening Type = Cleary-Luby-Jones-Plassman - AI\n");
      }
      else if (coarsen_type == 991)
      {
	 jxf_printf(" Coarsening Type = Ruge 1st pass only - AI \n");
      }
      else if (coarsen_type == 992)
      {
	 jxf_printf(" Coarsening Type = Ruge - AI\n");
      }
      else if (coarsen_type == 993)
      {
	 jxf_printf(" Coarsening Type = Ruge3 - AI\n");
      }
      else if (abs(coarsen_type) == 96)
      {
	 jxf_printf(" Coarsening Type = Falgout-CLJP - AI \n");
      }
      else if (abs(coarsen_type) == 98)
      {
	 jxf_printf(" Coarsening Type = PMIS - AI \n");
      }
      else if (abs(coarsen_type) == 910)
      {
	 jxf_printf(" Coarsening Type = HMIS - AI \n");
      }
      else if (coarsen_type == 908)
      {
	 jxf_printf(" Coarsening Type = CLJP -AI + PMIS - AI \n");
      }
      else if (coarsen_type == 918)
      {
         jxf_printf(" Coarsening Type = Ruge 1st pass only - AI + PMIS - AI \n");
      }
      else if (coarsen_type == 928)
      {
         jxf_printf(" Coarsening Type = Ruge - AI + PMIS - AI \n");
      }
      else if (coarsen_type == 938)
      {
         jxf_printf(" Coarsening Type = Ruge3 - AI + PMIS - AI \n");
      }
      else if (coarsen_type == 968)
      {
         jxf_printf(" Coarsening Type = Falgout-CLJP - AI + PMIS - AI \n");
      }

      if (agg_num_levels > 0)
      {
         jxf_printf("\n No. of levels of aggressive coarsening: %d\n\n", agg_num_levels);
         if (agg_interp_type == 4)
            jxf_printf(" Interpolation on agg. levels = multipass interpolation\n");
      }

      if (coarsen_type)
      {
         jxf_printf(" measures are determined %s\n\n", (measure_type ? "globally" : "locally"));
      }

#ifdef JXF_NO_GLOBAL_PARTITION
      jxf_printf("\n No global partition option chosen.\n\n");
#endif

      if (interp_type == 0)
      {
	 jxf_printf(" Interpolation = modified classical interpolation\n");
      }
      else if (interp_type == 3) 
      {
	 jxf_printf(" Interpolation = direct interpolation (with separation of weights)\n");
      }
      else if (interp_type == 4) 
      {
	 jxf_printf(" Interpolation = multipass interpolation\n");
      }
      else if (interp_type == 5) 
      {
	 jxf_printf(" Interpolation = multipass interpolation (with separation of weights)\n");
      }
      else if (interp_type == 6) 
      {
	 jxf_printf(" Interpolation = extended classical modified interpolation\n");
      }
      else if (interp_type == 7) 
      {
	 jxf_printf(" Interpolation = extended (if no common C neighbor) classical modified interpolation\n");
      }
      else if (interp_type == 8) 
      {
	 jxf_printf(" Interpolation = standard interpolation\n");
      }
      else if (interp_type == 9) 
      {
	 jxf_printf(" Interpolation = standard interpolation (with separation of weights)\n");
      }
      else if (interp_type == 100) 
      {
	 jxf_printf(" Interpolation = one-point interpolation (strongest C neighbor)\n");
      }
      else if (interp_type == 20) 
      {
	 jxf_printf(" Interpolation = M^{-1}*A^T*R^T\n");
      }

      if (restri_type == 1)
      {
     jxf_printf(" Restriction = local approximate ideal restriction (AIR-1)\n");
      }
      else if (restri_type == 2)
      {
     jxf_printf(" Restriction = local approximate ideal restriction (AIR-2)\n");
      }
      else if (restri_type)
      {
     jxf_printf(" Restriction = one-point restriction\n");
      }

      jxf_printf("\nOperator Matrix Information:\n\n");

      jxf_printf("            nonzero         entries p");
      jxf_printf("er row        row sums\n");
      jxf_printf("lev   rows  entries  sparse  min  max   ");
      jxf_printf("avg       min         max\n");
      jxf_printf("===================================================================\n");
   }
  
  
  /*-----------------------------------------------------
   *  Enter Statistics Loop
   *-----------------------------------------------------*/

   num_coeffs = jxf_CTAlloc(JXF_Real, num_levels);

   num_variables = jxf_CTAlloc(JXF_Real, num_levels);

   for (level = 0; level < num_levels; level ++)
   { 
       A_diag      = jxf_hpCSRMatrixDiag(A_array[level]);
       A_diag_data = jxf_CSRMatrixData(A_diag);
       A_diag_i    = jxf_CSRMatrixI(A_diag);
         
       A_offd      = jxf_hpCSRMatrixOffd(A_array[level]);   
       A_offd_data = jxf_CSRMatrixData(A_offd);
       A_offd_i    = jxf_CSRMatrixI(A_offd);
         
       row_starts = jxf_hpCSRMatrixRowStarts(A_array[level]);
         
       fine_size         = jxf_hpCSRMatrixGlobalNumRows(A_array[level]);
       global_nonzeros   = jxf_hpCSRMatrixDNumNonzeros(A_array[level]);
       num_coeffs[level] = global_nonzeros;
       num_variables[level] = (JXF_Real) fine_size;
         
       sparse = global_nonzeros / ((JXF_Real) fine_size * (JXF_Real) fine_size);

       min_entries = 0;
       max_entries = 0;
       min_rowsum  = 0.0;
       max_rowsum  = 0.0;
         
       if (jxf_CSRMatrixNumRows(A_diag))
       {
          min_entries = (A_diag_i[1] - A_diag_i[0]) + (A_offd_i[1] - A_offd_i[0]);
          for (j = A_diag_i[0]; j < A_diag_i[1]; j ++)
          {
             min_rowsum += A_diag_data[j];
          }
          for (j = A_offd_i[0]; j < A_offd_i[1]; j ++)
          {
             min_rowsum += A_offd_data[j];
          }
            
          //jxf_printf("rowsum_init = %e\n", min_rowsum);

          max_rowsum = min_rowsum;
            
          for (j = 0; j < jxf_CSRMatrixNumRows(A_diag); j ++)
          {
             entries = (A_diag_i[j+1] - A_diag_i[j]) + (A_offd_i[j+1] - A_offd_i[j]);
             min_entries = jxf_min(entries, min_entries);
             max_entries = jxf_max(entries, max_entries);
               
             rowsum = 0.0;
             for (i = A_diag_i[j]; i < A_diag_i[j+1]; i ++)
             {
                rowsum += A_diag_data[i];
             }
               
             for (i = A_offd_i[j]; i < A_offd_i[j+1]; i ++)
             {
                rowsum += A_offd_data[i];
             }

#if 0            
             if (rowsum < -1.0e-1) {
                jxf_printf("i= %d, rowsum = %e\n", j, rowsum);
                for (i = A_diag_i[j]; i < A_diag_i[j+1]; i ++)
                {
                   jxf_printf("ia= %d, A_diag = %e\n", i, A_diag_data[i]);
                }
               
                for (i = A_offd_i[j]; i < A_offd_i[j+1]; i ++)
                {
                   jxf_printf("ia= %d, A_offd = %e\n", i, A_offd_data[i]);
                }

             }
#endif
               
             min_rowsum = jxf_min(rowsum, min_rowsum);
             max_rowsum = jxf_max(rowsum, max_rowsum);
          }
       }
       avg_entries = global_nonzeros / ((JXF_Real) fine_size);

      
#ifdef JXF_NO_GLOBAL_PARTITION       

       numrows = row_starts[1] - row_starts[0];
       
      /*------------------------------------------------------------- 
       * if we don't have any rows, then don't have this 
       * count toward min row sum or min num entries. 
       *------------------------------------------------------------*/       
       if (!numrows) 
       {
          min_entries = 1000000;
          min_rowsum  = 1.0e7;
       }
       
       send_buff[0] = - (JXF_Real) min_entries;
       send_buff[1] = (JXF_Real) max_entries;
       send_buff[2] = - min_rowsum;
       send_buff[3] = max_rowsum;

       jxf_MPI_Reduce(send_buff, gather_buff, 4, JXF_MPI_REAL, MPI_MAX, 0, comm);
       
       if (my_id == 0)
       {
          global_min_e = - gather_buff[0];
          global_max_e =   gather_buff[1];
          global_min_rsum = - gather_buff[2];
          global_max_rsum =   gather_buff[3];
          
          jxf_printf("%2d %7d %8.0f  %0.3f  %4d %4d",
                 level, fine_size, global_nonzeros, sparse, global_min_e, global_max_e);
          jxf_printf("  %4.1f  %10.3e  %10.3e\n", avg_entries, global_min_rsum, global_max_rsum);
       }
       
#else

       send_buff[0] = (JXF_Real) min_entries;
       send_buff[1] = (JXF_Real) max_entries;
       send_buff[2] = min_rowsum;
       send_buff[3] = max_rowsum;
       
       jxf_MPI_Gather(send_buff, 4, JXF_MPI_REAL, gather_buff, 4, JXF_MPI_REAL, 0, comm);

       if (my_id == 0)
       {
          global_min_e = 1000000;
          global_max_e = 0;
          global_min_rsum = 1.0e7;
          global_max_rsum = 0.0;
          for (j = 0; j < num_procs; j ++)
          {
             numrows = row_starts[j+1] - row_starts[j];
             if (numrows)
             {
                global_min_e    = jxf_min(global_min_e, (JXF_Int) gather_buff[j*4]);
                global_min_rsum = jxf_min(global_min_rsum, gather_buff[j*4 + 2]);
             }
             global_max_e    = jxf_max(global_max_e, (JXF_Int) gather_buff[j*4 + 1]);
             global_max_rsum = jxf_max(global_max_rsum, gather_buff[j*4 + 3]);
          }

          jxf_printf("%2d %7d %8.0f  %0.3f  %4d %4d",
                  level, fine_size, global_nonzeros, sparse, global_min_e, global_max_e);
          jxf_printf("  %4.1f  %10.3e  %10.3e\n", avg_entries, global_min_rsum, global_max_rsum);
       }

#endif
    
   }

   if (my_id == 0)
   {
      jxf_printf("\n\nInterpolation Matrix Information:\n");
      jxf_printf("                 entries/row    min     max");
      jxf_printf("         row sums\n");
      jxf_printf("lev  rows cols    min max  ");
      jxf_printf("   weight   weight     min       max \n");
      jxf_printf("=================================================================\n");
   }
  

  /*-----------------------------------------------------
   *  Enter Statistics Loop
   *-----------------------------------------------------*/

   for (level = 0; level < num_levels-1; level ++)
   {
      P_diag      = jxf_ParCSRMatrixDiag(P_array[level]);
      P_diag_data = jxf_CSRMatrixData(P_diag);
      P_diag_i    = jxf_CSRMatrixI(P_diag);
         
      P_offd      = jxf_ParCSRMatrixOffd(P_array[level]);   
      P_offd_data = jxf_CSRMatrixData(P_offd);
      P_offd_i    = jxf_CSRMatrixI(P_offd);
         
      row_starts = jxf_ParCSRMatrixRowStarts(P_array[level]);
         
      fine_size       = jxf_ParCSRMatrixGlobalNumRows(P_array[level]);
      coarse_size     = jxf_ParCSRMatrixGlobalNumCols(P_array[level]);
      global_nonzeros = jxf_ParCSRMatrixNumNonzeros(P_array[level]);
         
      min_weight  = 1.0;
      max_weight  = 0.0;
      max_rowsum  = 0.0;
      min_rowsum  = 0.0;
      min_entries = 0;
      max_entries = 0;
         
      if (jxf_CSRMatrixNumRows(P_diag))
      {
         if (jxf_CSRMatrixNumCols(P_diag)) 
         {
            min_weight = P_diag_data[0];
         }
         
         for (j = P_diag_i[0]; j < P_diag_i[1]; j ++)
         {
            min_weight = jxf_min(min_weight, P_diag_data[j]);
            if (P_diag_data[j] != 1.0)
            {
               max_weight = jxf_max(max_weight, P_diag_data[j]);
            }
            min_rowsum += P_diag_data[j];
         }
         for (j = P_offd_i[0]; j < P_offd_i[1]; j ++)
         {        
            min_weight = jxf_min(min_weight, P_offd_data[j]); 
            if (P_offd_data[j] != 1.0)
            {
               max_weight = jxf_max(max_weight, P_offd_data[j]);  
            }   
            min_rowsum += P_offd_data[j];
         }
            
         max_rowsum = min_rowsum;
            
         min_entries = (P_diag_i[1] - P_diag_i[0]) + (P_offd_i[1] - P_offd_i[0]); 
         max_entries = 0;
            
         for (j = 0; j < jxf_CSRMatrixNumRows(P_diag); j ++)
         {
            entries = (P_diag_i[j+1] - P_diag_i[j]) + (P_offd_i[j+1] - P_offd_i[j]);
            min_entries = jxf_min(entries, min_entries);
            max_entries = jxf_max(entries, max_entries);
               
            rowsum = 0.0;
            for (i = P_diag_i[j]; i < P_diag_i[j+1]; i ++)
            {
               min_weight = jxf_min(min_weight, P_diag_data[i]);
               if (P_diag_data[i] != 1.0)
               {
                  max_weight = jxf_max(max_weight, P_diag_data[i]);
               }
               rowsum += P_diag_data[i];
            }
               
            for (i = P_offd_i[j]; i < P_offd_i[j+1]; i ++)
            {
               min_weight = jxf_min(min_weight, P_offd_data[i]);
               if (P_offd_data[i] != 1.0)
               { 
                  max_weight = jxf_max(max_weight, P_offd_data[i]);
               }
               rowsum += P_offd_data[i];
            }
               
            min_rowsum = jxf_min(rowsum, min_rowsum);
            max_rowsum = jxf_max(rowsum, max_rowsum);
         }
         
      }
      avg_entries = ((JXF_Real) global_nonzeros) / ((JXF_Real) fine_size);


#ifdef JXF_NO_GLOBAL_PARTITION

      numrows = row_starts[1] - row_starts[0];

      /*------------------------------------------------------------- 
       * if we don't have any rows, then don't have this 
       * count toward min row sum or min num entries. 
       *------------------------------------------------------------*/        
      if (!numrows)
      {
         min_entries = 1000000;
         min_rowsum  = 1.0e7;
         min_weight  = 1.0e7;
       }
       
      send_buff[0] = - (JXF_Real) min_entries;
      send_buff[1] = (JXF_Real) max_entries;
      send_buff[2] = - min_rowsum;
      send_buff[3] = max_rowsum;
      send_buff[4] = - min_weight;
      send_buff[5] = max_weight;

      jxf_MPI_Reduce(send_buff, gather_buff, 6, JXF_MPI_REAL, MPI_MAX, 0, comm);

      if (my_id == 0)
      {
         global_min_e    = -gather_buff[0];
         global_max_e    =  gather_buff[1];
         global_min_rsum = -gather_buff[2];
         global_max_rsum =  gather_buff[3];
         global_min_wt   = -gather_buff[4];
         global_max_wt   =  gather_buff[5];
         jxf_printf("%2d %5d x %-5d %3d %3d", level, fine_size, coarse_size,  global_min_e, global_max_e);
         jxf_printf("  %10.3e %9.3e %9.3e %9.3e\n", global_min_wt, global_max_wt, global_min_rsum, global_max_rsum);
      }

#else
      
      send_buff[0] = (JXF_Real) min_entries;
      send_buff[1] = (JXF_Real) max_entries;
      send_buff[2] = min_rowsum;
      send_buff[3] = max_rowsum;
      send_buff[4] = min_weight;
      send_buff[5] = max_weight;
      
      jxf_MPI_Gather(send_buff, 6, JXF_MPI_REAL, gather_buff, 6, JXF_MPI_REAL, 0, comm);
      
      if (my_id == 0)
      {
         global_min_e    = 1000000;
         global_max_e    = 0;
         global_min_rsum = 1.0e7;
         global_max_rsum = 0.0;
         global_min_wt   = 1.0e7;
         global_max_wt   = 0.0;
         
         for (j = 0; j < num_procs; j ++)
         {
            numrows = row_starts[j+1] - row_starts[j];
            if (numrows)
            {
               global_min_e = jxf_min(global_min_e, (JXF_Int) gather_buff[j*6]);
               global_min_rsum = jxf_min(global_min_rsum, gather_buff[j*6+2]);
               global_min_wt = jxf_min(global_min_wt, gather_buff[j*6+4]);
            }
            global_max_e = jxf_max(global_max_e, (JXF_Int) gather_buff[j*6+1]);
            global_max_rsum = jxf_max(global_max_rsum, gather_buff[j*6+3]);
            global_max_wt = jxf_max(global_max_wt, gather_buff[j*6+5]);
         }
         jxf_printf("%2d %5d x %-5d %3d %3d", level, fine_size, coarse_size,  global_min_e, global_max_e);
         jxf_printf("  %10.3e %9.3e %9.3e %9.3e\n", global_min_wt, global_max_wt, global_min_rsum, global_max_rsum);
      }

#endif

   }

   if (restri_type)
   {

   if (my_id == 0)
   {
      jxf_printf("\n\nRestriction Matrix Information:\n");
      jxf_printf("                 entries/row    min     max");
      jxf_printf("         row sums       size(ls)\n");
      jxf_printf("lev  rows cols    min max  ");
      jxf_printf("   weight   weight     min       max      max\n");
      jxf_printf("=========================================================================\n");
   }


  /*-----------------------------------------------------
   *  Enter Statistics Loop
   *-----------------------------------------------------*/

   for (level = 0; level < num_levels-1; level ++)
   {
      P_diag      = jxf_ParCSRMatrixDiag(R_array[level]);
      P_diag_data = jxf_CSRMatrixData(P_diag);
      P_diag_i    = jxf_CSRMatrixI(P_diag);
         
      P_offd      = jxf_ParCSRMatrixOffd(R_array[level]);   
      P_offd_data = jxf_CSRMatrixData(P_offd);
      P_offd_i    = jxf_CSRMatrixI(P_offd);
         
      row_starts = jxf_ParCSRMatrixRowStarts(R_array[level]);
         
      fine_size       = jxf_ParCSRMatrixGlobalNumRows(R_array[level]);
      coarse_size     = jxf_ParCSRMatrixGlobalNumCols(R_array[level]);
      global_nonzeros = jxf_ParCSRMatrixNumNonzeros(R_array[level]);
         
      min_weight  = 1.0;
      max_weight  = 0.0;
      max_rowsum  = 0.0;
      min_rowsum  = 0.0;
      min_entries = 0;
      max_entries = 0;
         
      if (jxf_CSRMatrixNumRows(P_diag))
      {
         if (jxf_CSRMatrixNumCols(P_diag)) 
         {
            min_weight = P_diag_data[0];
         }
         
         for (j = P_diag_i[0]; j < P_diag_i[1]; j ++)
         {
            min_weight = jxf_min(min_weight, P_diag_data[j]);
            if (P_diag_data[j] != 1.0)
            {
               max_weight = jxf_max(max_weight, P_diag_data[j]);
            }
            min_rowsum += P_diag_data[j];
         }
         for (j = P_offd_i[0]; j < P_offd_i[1]; j ++)
         {        
            min_weight = jxf_min(min_weight, P_offd_data[j]); 
            if (P_offd_data[j] != 1.0)
            {
               max_weight = jxf_max(max_weight, P_offd_data[j]);  
            }   
            min_rowsum += P_offd_data[j];
         }
            
         max_rowsum = min_rowsum;
            
         min_entries = (P_diag_i[1] - P_diag_i[0]) + (P_offd_i[1] - P_offd_i[0]); 
         max_entries = 0;
            
         for (j = 0; j < jxf_CSRMatrixNumRows(P_diag); j ++)
         {
            entries = (P_diag_i[j+1] - P_diag_i[j]) + (P_offd_i[j+1] - P_offd_i[j]);
            min_entries = jxf_min(entries, min_entries);
            max_entries = jxf_max(entries, max_entries);
               
            rowsum = 0.0;
            for (i = P_diag_i[j]; i < P_diag_i[j+1]; i ++)
            {
               min_weight = jxf_min(min_weight, P_diag_data[i]);
               if (P_diag_data[i] != 1.0)
               {
                  max_weight = jxf_max(max_weight, P_diag_data[i]);
               }
               rowsum += P_diag_data[i];
            }
               
            for (i = P_offd_i[j]; i < P_offd_i[j+1]; i ++)
            {
               min_weight = jxf_min(min_weight, P_offd_data[i]);
               if (P_offd_data[i] != 1.0)
               { 
                  max_weight = jxf_max(max_weight, P_offd_data[i]);
               }
               rowsum += P_offd_data[i];
            }
               
            min_rowsum = jxf_min(rowsum, min_rowsum);
            max_rowsum = jxf_max(rowsum, max_rowsum);
         }
         
      }
      avg_entries = ((JXF_Real) global_nonzeros) / ((JXF_Real) fine_size);


#ifdef JXF_NO_GLOBAL_PARTITION

      numrows = row_starts[1] - row_starts[0];

      /*------------------------------------------------------------- 
       * if we don't have any rows, then don't have this 
       * count toward min row sum or min num entries. 
       *------------------------------------------------------------*/        
      if (!numrows)
      {
         min_entries = 1000000;
         min_rowsum  = 1.0e7;
         min_weight  = 1.0e7;
       }
       
      send_buff[0] = - (JXF_Real) min_entries;
      send_buff[1] = (JXF_Real) max_entries;
      send_buff[2] = - min_rowsum;
      send_buff[3] = max_rowsum;
      send_buff[4] = - min_weight;
      send_buff[5] = max_weight;

      jxf_MPI_Reduce(send_buff, gather_buff, 6, JXF_MPI_REAL, MPI_MAX, 0, comm);

      if (my_id == 0)
      {
         global_min_e    = -gather_buff[0];
         global_max_e    =  gather_buff[1];
         global_min_rsum = -gather_buff[2];
         global_max_rsum =  gather_buff[3];
         global_min_wt   = -gather_buff[4];
         global_max_wt   =  gather_buff[5];
         jxf_printf("%2d %5d x %-5d %3d %3d", level, fine_size, coarse_size,  global_min_e, global_max_e);
         jxf_printf("  %10.3e %9.3e %9.3e %9.3e", global_min_wt, global_max_wt, global_min_rsum, global_max_rsum);
         jxf_printf("  %3d\n", AIR_maxsize_ls[level]);
      }

#else
      
      send_buff[0] = (JXF_Real) min_entries;
      send_buff[1] = (JXF_Real) max_entries;
      send_buff[2] = min_rowsum;
      send_buff[3] = max_rowsum;
      send_buff[4] = min_weight;
      send_buff[5] = max_weight;
      
      jxf_MPI_Gather(send_buff, 6, JXF_MPI_REAL, gather_buff, 6, JXF_MPI_REAL, 0, comm);
      
      if (my_id == 0)
      {
         global_min_e    = 1000000;
         global_max_e    = 0;
         global_min_rsum = 1.0e7;
         global_max_rsum = 0.0;
         global_min_wt   = 1.0e7;
         global_max_wt   = 0.0;
         
         for (j = 0; j < num_procs; j ++)
         {
            numrows = row_starts[j+1] - row_starts[j];
            if (numrows)
            {
               global_min_e = jxf_min(global_min_e, (JXF_Int) gather_buff[j*6]);
               global_min_rsum = jxf_min(global_min_rsum, gather_buff[j*6+2]);
               global_min_wt = jxf_min(global_min_wt, gather_buff[j*6+4]);
            }
            global_max_e = jxf_max(global_max_e, (JXF_Int) gather_buff[j*6+1]);
            global_max_rsum = jxf_max(global_max_rsum, gather_buff[j*6+3]);
            global_max_wt = jxf_max(global_max_wt, gather_buff[j*6+5]);
         }
         jxf_printf("%2d %5d x %-5d %3d %3d", level, fine_size, coarse_size,  global_min_e, global_max_e);
         jxf_printf("  %10.3e %9.3e %9.3e %9.3e", global_min_wt, global_max_wt, global_min_rsum, global_max_rsum);
         jxf_printf("  %3d\n", AIR_maxsize_ls[level]);
      }

#endif

   }

   }

   total_variables = 0;
   operat_cmplxty  = 0;
   
   for (j = 0; j < jxf_ParAMGDataNumLevels(amg_data); j ++)
   {
      operat_cmplxty  +=  num_coeffs[j] / num_coeffs[0];
      total_variables += num_variables[j];
   }
   
   if (num_variables[0] != 0)
   {
      grid_cmplxty = total_variables / num_variables[0];
   }
 
   if (my_id == 0)
   {
      jxf_printf("\n\n     Complexity:    grid = %f\n", grid_cmplxty);
      jxf_printf("                operator = %f\n", operat_cmplxty);
   }

   if (my_id == 0) jxf_printf("\n\n");

   //if (my_id == 0 && amg_print_level == 1)
   if (my_id == 0)
   { 
      jxf_printf("\n\njxf_PAMG SOLVER PARAMETERS:\n\n");
      jxf_printf("  Maximum number of cycles:         %d \n",max_iter);
      jxf_printf("  Stopping Tolerance:               %e \n",tol); 
      /* peghoty 2009/07/27 */
      jxf_printf("  Convergence Criteria:");
      switch(conv_criteria)
      {
         case 0:
         jxf_printf("  Relative Residual Norm-2 (default)\n");
         break;
         case 1:
         jxf_printf("  Relative Error Norm-1\n");
         break;
         case 11:
         jxf_printf("  Relative Point-wise Error Norm-1\n");
         break;
         case 2:
         jxf_printf("  Relative Error Norm-2\n");
         break;
      }
      jxf_printf("  Cycle type (1 = V, 2 = W, etc.):  %d\n\n", cycle_type);
      jxf_printf("  Relaxation Parameters:\n");
      jxf_printf("   Visiting Grid:                     down   up  coarse\n");
      jxf_printf("            Number of sweeps:         %4d   %2d  %4d \n",
              num_grid_sweeps[1], num_grid_sweeps[2], num_grid_sweeps[3]);
      jxf_printf("   Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:   %4d   %2d  %4d \n",
              grid_relax_type[1], grid_relax_type[2], grid_relax_type[3]);
      jxf_printf("   Point types, partial sweeps (1=C, -1=F):\n");

      if (grid_relax_points && grid_relax_type[1] != 8)
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[1][j]);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[2][j]);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[3][j]);
         }
         jxf_printf("\n\n");
      }
      else if (relax_order == 1 && grid_relax_type[1] != 8)
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d  %2d", one, minus_one);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d  %2d", minus_one, one);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n\n");
      }
      else 
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n\n");
      }

      for (j = 0; j < num_levels; j ++)
      {
         if (relax_weight[j] != 1)
         {
	    jxf_printf(" Relaxation Weight %f level %d\n", relax_weight[j], j);
	 }
      }
      for (j = 0; j < num_levels; j ++)
      {
         if (omega[j] != 1)
         {
            jxf_printf(" Outer relaxation weight %f level %d\n", omega[j], j);
         }
      }
   }

   jxf_TFree(num_coeffs);
   jxf_TFree(num_variables);
   jxf_TFree(send_buff);
   jxf_TFree(gather_buff);
   
   return(0);
}


/*!
 * \fn JXF_Int jxf_PAMGSolveStatus
 * \brief Information outputing of the Solution phase of PAMG.
 * \date 2011/09/03
 */ 
JXF_Int    
jxf_PAMGSolveStatus( void *data )
{ 
   jxf_ParAMGData  *amg_data = data;
   JXF_Int j;
   JXF_Int one = 1;
   JXF_Int minus_one = -1;
   JXF_Int zero = 0;

  /*-------------------------------------------------------------------------
   * Get the amg_data data
   *------------------------------------------------------------------------*/
   JXF_Int     amg_print_level   = jxf_ParAMGDataPrintLevel(amg_data);
   JXF_Int     max_iter          = jxf_ParAMGDataMaxIter(amg_data);
   JXF_Int     cycle_type        = jxf_ParAMGDataCycleType(amg_data);
   JXF_Int    *num_grid_sweeps   = jxf_ParAMGDataNumGridSweeps(amg_data);  
   JXF_Int    *grid_relax_type   = jxf_ParAMGDataGridRelaxType(amg_data);
   JXF_Int   **grid_relax_points = jxf_ParAMGDataGridRelaxPoints(amg_data);
   JXF_Int     relax_order       = jxf_ParAMGDataRelaxOrder(amg_data);
   JXF_Real  tol               = jxf_ParAMGDataTol(amg_data);
   JXF_Int     conv_criteria     = jxf_ParAMGDataConvCriteria(amg_data); /* peghoty 2009/07/27 */

   if (amg_print_level == 2 || amg_print_level == 3)
   { 
      jxf_printf("\n\nJXF_PAMG SOLVER PARAMETERS:\n\n");
      jxf_printf("  Maximum number of cycles:         %d \n", max_iter);
      jxf_printf("  Stopping Tolerance:               %e \n", tol);
      jxf_printf("  Convergence Criteria:"); /* peghoty 2009/07/27 */
      switch (conv_criteria)
      {
         case 0:
         jxf_printf("  Relative Residual Norm-2 (default)\n");
         break;
         case 1:
         jxf_printf("  Relative Error Norm-1\n");
         break;
         case 11:
         jxf_printf("  Relative Point-wise Error Norm-1\n");
         break;
         case 2:
         jxf_printf("  Relative Error Norm-2\n");
         break;
      }
      jxf_printf("  Cycle type (1 = V, 2 = W, etc.):  %d\n\n", cycle_type);
      jxf_printf("  Relaxation Parameters:\n");
      jxf_printf("   Visiting Grid:                     down   up  coarse\n");
      jxf_printf("            Number of sweeps:         %4d   %2d  %4d \n",
              num_grid_sweeps[1], num_grid_sweeps[2], num_grid_sweeps[3]);
      jxf_printf("   Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:   %4d   %2d  %4d \n",
              grid_relax_type[1], grid_relax_type[2], grid_relax_type[3]);
      jxf_printf("   Point types, partial sweeps (1=C, -1=F):\n");
      if (grid_relax_points)
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[1][j]);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[2][j]);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", grid_relax_points[3][j]);
         }
         jxf_printf("\n\n");
      }
      else if (relax_order == 1)
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d  %2d", one, minus_one);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d  %2d", minus_one, one);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n\n");
      }
      else 
      {
         jxf_printf("                  Pre-CG relaxation (down):");
         for (j = 0; j < num_grid_sweeps[1]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n");
         jxf_printf("                   Post-CG relaxation (up):");
         for (j = 0; j < num_grid_sweeps[2]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n");
         jxf_printf("                             Coarsest grid:");
         for (j = 0; j < num_grid_sweeps[3]; j ++)
         {
            jxf_printf("  %2d", zero);
         }
         jxf_printf("\n\n");
      }
      jxf_printf(" Output flag (print_level): %d \n", amg_print_level);
   }
 
   return 0;
}
