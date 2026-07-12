//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_coarsen_rrs0.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"


/*=======================================================================================
 * 
 * Multi-phase Parallel R-S coarsenning(MPRS), a new parallel coarsening algorithm 
 * for AMG based on R-S algorithm, proposed by Zeyao Mo & Xiaowen Xu, Oct,2003.                  
 * The code was devoloped by X.Xu, Apirl,2004. 
 *
 * 1st modified by X.Xu, Feb,2005.
 * 2nd modified by X.Xu, June,2005. 
 *
 *  (Due to the idea proposed by Z.Mo, the algorithm called RRS(Relaxed RS))
 * 
 *======================================================================================*/
 
#define C_PT  1
#define F_PT -1
#define COMMON_C_PT  2
#define Z_PT -2
#define UNDECIDED 0 

/*!
 * \fn JXF_Int jxf_PAMGCoarsenRRS0
 * \brief MPRS Coarsening routine.
 * \author Xu Xiaowen
 * \date Apirl,2004
 */  
JXF_Int
jxf_PAMGCoarsenRRS0( jxf_ParCSRMatrix  *S,
                    jxf_ParCSRMatrix  *A,
                    JXF_Int               measure_type,
                    JXF_Int               measure_type_rlx,
                    JXF_Int               number_syn_rlx,
                    JXF_Real            measure_threshold_rlx,
                    JXF_Int               debug_flag,
                    JXF_Int             **CF_marker_ptr )
{
   MPI_Comm             comm      = jxf_ParCSRMatrixComm(S);
   jxf_ParCSRCommPkg    *comm_pkg  = jxf_ParCSRMatrixCommPkg(A);
   jxf_ParCSRCommHandle *comm_handle;
   jxf_CSRMatrix *S_diag = jxf_ParCSRMatrixDiag(S);
   jxf_CSRMatrix *S_offd = jxf_ParCSRMatrixOffd(S);
   JXF_Int          *S_i           = jxf_CSRMatrixI(S_diag);
   JXF_Int          *S_j           = jxf_CSRMatrixJ(S_diag);
   JXF_Int          *S_offd_i      = jxf_CSRMatrixI(S_offd);
   JXF_Int          *S_offd_j      = jxf_CSRMatrixJ(S_offd);
   JXF_Int           num_variables = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int           num_cols_offd = jxf_CSRMatrixNumCols(S_offd);
                  
   jxf_CSRMatrix *S_ext   = NULL;
   JXF_Int          *S_ext_i = NULL;
   JXF_Int          *S_ext_j = NULL;
                 
   jxf_CSRMatrix *ST;
   JXF_Int          *ST_i;
   JXF_Int          *ST_j;
                 
   JXF_Int          *CF_marker;
   JXF_Int          *CF_marker_offd;

   JXF_Int          *measure_array;
   JXF_Int 	        *int_buf_data;

   JXF_Int           i, j, k, ji, jS;
   JXF_Int           index;
   JXF_Int           num_per;
   JXF_Int		 num_nonzeros;
   JXF_Int		 num_procs, my_id;
   JXF_Int		 num_sends = 0;
   JXF_Int		 first_col = 0, start;
   JXF_Int		 col_0 = 0, col_n = 0;

   jxf_LinkList   LoL_head;
   jxf_LinkList   LoL_tail;

   JXF_Int          *lists, *where;
   JXF_Int           measure, new_meas, measure_global;
   JXF_Real        measure_th;
   JXF_Int           num_left,global_num_left;
   JXF_Int           nabor, nabor_two;

   JXF_Int           ierr = 0;
   JXF_Int           f_pnt = F_PT;
   JXF_Real	 wall_time = 0.0;
   JXF_Int           *non_finish;
   JXF_Int           iter = 0;
   JXF_Int           iter_coar_dep = 0;
   JXF_Int           iter_coar_inf = 0;
   JXF_Int           iter_fine =0;
   JXF_Int           kk;

  /*----------------------------------------------------------
   * Initialize the C/F marker, LoL_head, LoL_tail  arrays
   *--------------------------------------------------------*/
   LoL_head = NULL;
   LoL_tail = NULL;
   lists = jxf_CTAlloc(JXF_Int, num_variables);
   where = jxf_CTAlloc(JXF_Int, num_variables);

   CF_marker = jxf_CTAlloc(JXF_Int, num_variables); // by xu
   for (j = 0; j < num_variables; j ++)
   {
      CF_marker[j] = UNDECIDED;
   } 

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
  
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A); 
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);

   if (num_cols_offd) S_offd_j = jxf_CSRMatrixJ(S_offd);

   jS = S_i[num_variables];

   ST = jxf_CSRMatrixCreate(num_variables, num_variables, jS);
   ST_i = jxf_CTAlloc(JXF_Int, num_variables + 1);
   ST_j = jxf_CTAlloc(JXF_Int, jS);
   jxf_CSRMatrixI(ST) = ST_i;
   jxf_CSRMatrixJ(ST) = ST_j;

  /*----------------------------------------------------------
   * generate transpose of S, ST
   *----------------------------------------------------------*/
   for (i = 0; i <= num_variables; i ++) ST_i[i] = 0;
 
   for (i = 0; i < jS; i ++)
   {
	 ST_i[S_j[i]+1] ++;
   }
   for (i = 0; i < num_variables; i ++)
   {
      ST_i[i+1] += ST_i[i];
   }
   for (i = 0; i < num_variables; i ++)
   {
      for (j = S_i[i]; j < S_i[i+1]; j ++)
      {
	 index = S_j[j];
       	 ST_j[ST_i[index]] = i;
       	 ST_i[index]++;
      }
   }
   for (i = num_variables; i > 0; i --)
   {
      ST_i[i] = ST_i[i-1];
   }
   ST_i[0] = 0;


  /*----------------------------------------------------------
   * Compute the measures
   *
   * The measures are given by the row sums of ST.
   * Hence, measure_array[i] is the number of influences
   * of variable i.
   * correct actual measures through adding influences from
   * neighbor processors
   *----------------------------------------------------------*/

   measure_array = jxf_CTAlloc(JXF_Int, num_variables);

   for (i = 0; i < num_variables; i ++)
   {
      measure_array[i] = ST_i[i+1] - ST_i[i];
   }

   if (num_procs > 1)
   {
      S_ext = jxf_ParCSRMatrixExtractBExt(S, A, 0);
      S_ext_i = jxf_CSRMatrixI(S_ext);
      S_ext_j = jxf_CSRMatrixJ(S_ext);
      num_nonzeros = S_ext_i[num_cols_offd];
      first_col = jxf_ParCSRMatrixFirstColDiag(S);
      col_0 = first_col - 1;
      col_n = col_0 + num_variables;
      if (measure_type)
      {
         for (i = 0; i < num_nonzeros; i ++)
         {
           index = S_ext_j[i] - first_col;
	   if (index > -1 && index < num_variables)
           measure_array[index] ++;
         }
      }
   }

  /*------------------------------------------------------
   * Loop until all points are either fine or coarse.
   *----------------------------------------------------*/

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();

  /* first coarsening phase */

  /*------------------------------------------------------
   *   Initialize the lists
   *-----------------------------------------------------*/
   num_left = num_variables;
 
   for (j = 0; j < num_variables; j ++) 
   {    
      measure = measure_array[j];
      if (measure > 0) 
      {
         jxf_enter_on_lists(&LoL_head, &LoL_tail, measure, j, lists, where);
      }
      else
      {
         if (measure < 0) jxf_printf("negative measure!\n");
         CF_marker[j] = f_pnt;
         for (k = S_i[j]; k < S_i[j+1]; k ++)
         {
            nabor = S_j[k];
            if (nabor < j)
            {
               new_meas = measure_array[nabor];
	       if (new_meas > 0)
	       {
                  jxf_remove_point(&LoL_head, &LoL_tail, new_meas, nabor, lists, where);
               }

               new_meas = ++(measure_array[nabor]);
                 
               jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas, nabor, lists, where);
            }
	    else
            {
               new_meas = ++(measure_array[nabor]);
            }
         }
         --num_left;
      }
   }
   
   non_finish = jxf_CTAlloc(JXF_Int, num_cols_offd);
   for (i = 0; i < num_cols_offd; i ++)
   {
      non_finish[i] = 1;
   }

  /****************************************************************
   *
   *  Main loop of MPRS coarsening first coloring pass.
   *
   *  WHILE there are still points to classify Parallel DO:
   *     num_per = 0;
   *     measure_th = th_m * global_max_measure;
   *     While ( nc+nf < num_variables && num_per < number_syn_rlx ) Do:
   *        1) find first point, i,  on list with max_measure.
   *           if (measure(i) < measure_th) break;
   *           make i a C-point, remove it from the lists.
   *           num_per++;
   *           nc++;
   *        2) For each point, j,  in S_i^T,
   *           a) Set j to be an F-point.
   *              nf++;
   *           b) For each point, k, in S_j
   *                  move k to the list in LoL with measure one
   *                  greater than it occupies (creating new LoL
   *                  entry if necessary).
   *        3) For each point, j,  in S_i,
   *                  move j to the list in LoL with measure one
   *                  smaller than it occupies (creating new LoL
   *                  entry if necessary).
   *     EndDo
   *     if ( num_per < number_syn_rlx ) measure_th = measure_th/2;
   *     Communication.
   *     Update the measure if necessary, and at the same time,
   *     Treat the 'F-F dependence' which obviate the 'C1'.
   *  ENDDO
   *
   ****************************************************************/

   CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
   int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

   for (j = 0; j < num_cols_offd; j ++)
   {
      CF_marker_offd[j] = UNDECIDED;
   } 

   while (1)
   {
      num_per = 0;
      if (num_left > 0)
      {
         index = LoL_head -> head;       // Feb,24,2005
         measure = measure_array[index]; // Feb,24,2005
      }
      
      if (measure_type_rlx)
      {
         jxf_MPI_Allreduce(&measure,&measure_global,1,JXF_MPI_INT,MPI_MAX,comm); // Feb,24,2005
      }
      else
      {
         measure_global = measure;
      }

      measure_th = measure_threshold_rlx*((JXF_Real) measure_global); // Feb,24,2005

      while (num_left > 0 && num_per < number_syn_rlx)
      {
         index = LoL_head -> head;
         measure = measure_array[index];
         
         if (measure < measure_th) break;  // Feb,24,2005
             
         CF_marker[index] = C_PT;
         measure = measure_array[index];
         measure_array[index] = 0;
         --num_left;
         ++num_per;
	 jxf_remove_point(&LoL_head, &LoL_tail, measure, index, lists, where);
	 
	/*----------------------------------------------
	 * treat the influence point
	 *---------------------------------------------*/
         for (j = ST_i[index]; j < ST_i[index+1]; j ++)
	 {
            nabor = ST_j[j];
            if (CF_marker[nabor] == UNDECIDED)
            {
               CF_marker[nabor] = F_PT;
               measure = measure_array[nabor];
               jxf_remove_point(&LoL_head, &LoL_tail, measure, nabor, lists, where);
               -- num_left;
               for (k = S_i[nabor]; k < S_i[nabor+1]; k ++)
               {
                  nabor_two = S_j[k];
                  if (CF_marker[nabor_two] == UNDECIDED)
                  {
                     measure = measure_array[nabor_two];
                     jxf_remove_point(&LoL_head, &LoL_tail, measure,nabor_two, lists, where);
                     new_meas = ++(measure_array[nabor_two]);      
                     jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas,nabor_two, lists, where);
                  }
               }
            }
         }
	 
	/*-------------------------------------------------------
	 * treat the depedence point
	 *-----------------------------------------------------*/
         for (j = S_i[index]; j < S_i[index+1]; j++)
         {
            nabor = S_j[j];
            if (CF_marker[nabor] == UNDECIDED)
            {
               measure = measure_array[nabor];
               jxf_remove_point(&LoL_head, &LoL_tail, measure, nabor, lists, where);
               measure_array[nabor] = --measure;
               if (measure > 0)
               {
                  jxf_enter_on_lists(&LoL_head, &LoL_tail, measure, nabor,lists, where);
               }
	       else
	       {
                  CF_marker[nabor] = F_PT;
                  -- num_left;
                  for (k = S_i[nabor]; k < S_i[nabor+1]; k ++)
                  {
                     nabor_two = S_j[k];
                     if (CF_marker[nabor_two] == UNDECIDED)
                     {
                        new_meas = measure_array[nabor_two];
                        jxf_remove_point(&LoL_head, &LoL_tail, new_meas,nabor_two, lists, where);
                        new_meas = ++(measure_array[nabor_two]);
                        jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas,nabor_two, lists, where);
                     }
                  }
	       }
            }
         }
      }  //end 'while'
      
      iter ++;
     
     /*****************************************
      * Test for convergence.
      ****************************************/
      
     /********************************************************************
      * Communication-1: exchange dependent boundary data for CF_marker    
      *******************************************************************/
      
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
        start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
        for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
        {
           int_buf_data[index++] = CF_marker[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
        }
      }

      if (num_procs > 1)
      {
         /* obtain my dependent points. */
         comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd); 
         jxf_ParCSRCommHandleDestroy(comm_handle);
      }
   
     /******************************************************************
      * Communication-2: exchange influential boundary data for CF_marker
      ******************************************************************/

     /********************************************************
      * Update the measure if necessary.
      ********************************************************/
      // consider the coarse points and its influence points.
      for (k = 0; k < num_variables; k ++) 
      {	      
	 if (CF_marker[k] == UNDECIDED)
	 {
	    for (j = S_offd_i[k]; j < S_offd_i[k+1]; j ++)
       	    {	     
	      i = S_offd_j[j];
	      if (CF_marker_offd[i] == 1) // consider the coarse points.
	      {
                    iter_coar_inf ++;
                    CF_marker[k] = F_PT;
                    measure = measure_array[k];
                    jxf_remove_point(&LoL_head, &LoL_tail, measure, k, lists, where);
                    -- num_left;
                    for (kk = S_i[k]; kk < S_i[k+1]; kk ++)
                    {
                       nabor = S_j[kk];
                       if (CF_marker[nabor] == UNDECIDED)
                       {
                          measure = measure_array[nabor];
                          jxf_remove_point(&LoL_head, &LoL_tail, measure,nabor, lists, where);
                          new_meas = ++(measure_array[nabor]);      
                          jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas,nabor, lists, where);
                       }
                    }
                    break;
              }
            }
         }
      }
           
      for (i = 0; i < num_cols_offd; i ++)
      {
         if (CF_marker_offd[i] == 1 && non_finish[i]) // consider the coarse points and its dependence points.
	 {
            non_finish[i] = 0;
            for (ji = S_ext_i[i]; ji < S_ext_i[i+1]; ji ++) // treate the depedence points
	    {
              j = S_ext_j[ji];
              if (j > col_0 && j < col_n)
	      {
	         j = j-first_col;
                 if (CF_marker[j] == UNDECIDED)
		 {
                    iter_coar_dep ++;
                    measure = measure_array[j];
                    jxf_remove_point(&LoL_head, &LoL_tail, measure, j, lists, where);
                    measure_array[j] = --measure;
                    if (measure > 0)
                    {
                       jxf_enter_on_lists(&LoL_head, &LoL_tail, measure, j,lists, where);
                    }
	            else
	            {
                       CF_marker[j] = F_PT;
                       -- num_left;
                       for (k = S_i[j]; k < S_i[j+1]; k ++)
                       {
                          nabor = S_j[k];
                          if (CF_marker[nabor] == UNDECIDED)
			  {
                             new_meas = measure_array[nabor];
                             jxf_remove_point(&LoL_head, &LoL_tail, new_meas,nabor, lists, where);
                             new_meas = ++(measure_array[nabor]);
                             jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas,nabor, lists, where);
                          }
                       }
	            }
                 }
	      }
	    }	    
	 }
	 if (CF_marker_offd[i] == -1 && non_finish[i]) // consider the fine points and its dependent points.
	 {
            non_finish[i] = 0;
            for (ji = S_ext_i[i]; ji < S_ext_i[i+1]; ji ++) // treat the dependent points
	    { 
               j = S_ext_j[ji];
	       if (j > col_0 && j < col_n)
	       {
		  j=j-first_col;
                  if (CF_marker[j] == UNDECIDED)
                  {   
                     iter_fine++;
                     measure = measure_array[j];
                     jxf_remove_point(&LoL_head, &LoL_tail, measure,j, lists, where);
                     new_meas = ++(measure_array[j]);      
                     jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas,j, lists, where);
                  }
               }
	    }
	 } //end fine points
      } //end 'update'
    

     /*****************************************
      * Test for convergence.
      ****************************************/

      jxf_MPI_Allreduce(&num_left,&global_num_left,1,JXF_MPI_INT,MPI_SUM,comm);

      if (global_num_left == 0) break;
      if (num_left <= 0) measure = 0;
      
   } //end 'while'
   
   if (debug_flag == 3)
   {
       jxf_printf("my_id= %d,iter= %d,iter_coar_dep= %d,iter_coar_inf= %d,iter_fine= %d\n",
               my_id,iter,iter_coar_dep,iter_coar_inf,iter_fine);
   }
       
   jxf_TFree(measure_array);
   jxf_CSRMatrixDestroy(ST);

   if (debug_flag == 3)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d    Coarsen 1st pass = %f\n", my_id, wall_time);
   }

   jxf_TFree(lists);
   jxf_TFree(where);
   jxf_TFree(LoL_head);
   jxf_TFree(LoL_tail);

  /*---------------------------------------------------
   * Clean up and return
   *---------------------------------------------------*/

   jxf_TFree(CF_marker_offd);
   jxf_TFree(int_buf_data);
   if ( num_procs > 1) 
   {
      jxf_CSRMatrixDestroy(S_ext); 
   }
   
   *CF_marker_ptr = CF_marker;

   return (ierr);
}
