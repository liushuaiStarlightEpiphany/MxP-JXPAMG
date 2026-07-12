//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_coarsen_ruge.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

#define C_PT 1
#define F_PT -1
#define Z_PT -2
#define SF_PT -3  /* special fine points */
#define UNDECIDED 0 

/*!
 * \fn JXF_Int jxf_PAMGCoarsenRuge
 * \brief Ruge's coarsening algorithm.
 * CF_init == 9: for XML coarsening,
 *               reset zero for CF_marker on ESS-parts while preversing CF_marker on AI-parts) 
 * \date 2011/09/03
 */ 
JXF_Int
jxf_PAMGCoarsenRugeAI( jxf_ParCSRMatrix    *par_S,
                    jxf_ParCSRMatrix    *par_A,
                    JXF_Real             *AI_measure,
                    JXF_Int                 CF_init,
                    JXF_Int                 measure_type,
                    JXF_Int                 coarsen_type,
                    JXF_Int                 debug_flag,
                    JXF_Int               **CF_marker_ptr )
{
   MPI_Comm             comm     = jxf_ParCSRMatrixComm(par_S);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_S);
   jxf_ParCSRCommHandle *comm_handle;
   
   jxf_CSRMatrix *S_diag        = jxf_ParCSRMatrixDiag(par_S);
   jxf_CSRMatrix *S_offd        = jxf_ParCSRMatrixOffd(par_S);
   JXF_Int          *S_i           = jxf_CSRMatrixI(S_diag);
   JXF_Int          *S_j           = jxf_CSRMatrixJ(S_diag);
   JXF_Int          *S_offd_i      = jxf_CSRMatrixI(S_offd);
   JXF_Int          *S_offd_j      = NULL;
   JXF_Int           num_variables = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int           num_cols_offd = jxf_CSRMatrixNumCols(S_offd);
   JXF_Int          *col_map_offd  = jxf_ParCSRMatrixColMapOffd(par_S);
                  
   jxf_CSRMatrix *S_ext   = NULL;
   JXF_Int          *S_ext_i = NULL;
   JXF_Int          *S_ext_j = NULL;
                 
   jxf_CSRMatrix *ST;
   JXF_Int          *ST_i;
   JXF_Int          *ST_j;
                 
   JXF_Int          *CF_marker;
   JXF_Int          *CF_marker_offd = NULL;
   JXF_Int           ci_tilde = -1;
   JXF_Int           ci_tilde_mark = -1;
   JXF_Int           ci_tilde_offd = -1;
   JXF_Int           ci_tilde_offd_mark = -1;

   JXF_Int          *measure_array;
   JXF_Int          *graph_array;
   JXF_Int 	        *int_buf_data = NULL;
   JXF_Int 	        *ci_array = NULL;

   JXF_Int           i, j = 0, k, jS;
   JXF_Int           ji, jj, jk, jm, index;
   JXF_Int		 set_empty = 1;
   JXF_Int		 C_i_nonempty = 0;
   JXF_Int		 num_nonzeros;
   JXF_Int		 num_procs, my_id;
   JXF_Int		 num_sends = 0;
   JXF_Int		 first_col = 0;
   JXF_Int		 start = 0;
   JXF_Int		 col_0 = 0;
   JXF_Int		 col_n = 0;

   jxf_LinkList   LoL_head;
   jxf_LinkList   LoL_tail;

   JXF_Int          *lists, *where;
   JXF_Int           measure, new_meas;
   JXF_Int           num_left, elmt;
   JXF_Int           nabor, nabor_two;

   JXF_Int           ierr = 0;
   JXF_Int           use_commpkg_A = 0;
   JXF_Int           break_var = 0;
   JXF_Int           f_pnt = F_PT;
   JXF_Real        wall_time = 0.0;

   if (coarsen_type < 0) coarsen_type = -coarsen_type;

  /*-------------------------------------------------------
   * Initialize the C/F marker, LoL_head, LoL_tail  arrays
   *-------------------------------------------------------*/

   LoL_head = NULL;
   LoL_tail = NULL;
   lists = jxf_CTAlloc(JXF_Int, num_variables);
   where = jxf_CTAlloc(JXF_Int, num_variables);

  /*--------------------------------------------------------------
   * Use a CSR strength matrix, par_S.
   *
   * For now, the "strength" of dependence/influence is defined in
   * the following way: i depends on j if
   *     aij > jxf_max (k != i) aik,    aii < 0
   * or
   *     aij < jxf_min (k != i) aik,    aii >= 0
   * Then S_ij = 1, else S_ij = 0.
   *
   * NOTE: the entries are negative initially, corresponding
   * to "unaccounted-for" dependence.
   *----------------------------------------------------------------*/

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (!comm_pkg)
   {
      use_commpkg_A = 1;
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
   }

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
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

   for (i = 0; i <= num_variables; i ++)
   {
      ST_i[i] = 0;
   }
 
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
       	 ST_i[index] ++;
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

   /* special case for Falgout coarsening */
   if (coarsen_type == 6) 
   {
      f_pnt = Z_PT;
      coarsen_type = 1;
   }
   if (coarsen_type == 10)
   {
      f_pnt = Z_PT;
      coarsen_type = 11;
   }

   if ( (measure_type || (coarsen_type != 1 && coarsen_type != 11)) && num_procs > 1 ) 
   {
      if (use_commpkg_A)
      {
         S_ext = jxf_ParCSRMatrixExtractBExt(par_S, par_A, 0);
      }
      else
      {
         S_ext = jxf_ParCSRMatrixExtractBExt(par_S, par_S, 0);
      }
      S_ext_i = jxf_CSRMatrixI(S_ext);
      S_ext_j = jxf_CSRMatrixJ(S_ext);
      num_nonzeros = S_ext_i[num_cols_offd];
      first_col = jxf_ParCSRMatrixFirstColDiag(par_S);
      col_0 = first_col - 1;
      col_n = col_0 + num_variables;
      if (measure_type)
      {
	 for (i = 0; i < num_nonzeros; i ++)
         {
            index = S_ext_j[i] - first_col;
            if (index > -1 && index < num_variables)
            {
               measure_array[index] ++;
            }
         } 
      } 
   }

  /*---------------------------------------------------
   * Add AI_measure to measure_array.
   *---------------------------------------------------*/
   if (CF_init == 0) {
      for (i = 0; i < num_variables; i ++)
      { 
         if (AI_measure[i] > 0.1 ) {
            measure_array[i] += 1000*AI_measure[i];
            //jxf_printf("i= %d, measure= %d, ai_measure= %f \n", i, measure_array[i], AI_measure[i]);
         }
      }
   } else if (CF_init == 9){
      for (i = 0; i < num_variables; i ++)
      { 
         if (AI_measure[i] > 0.1 ) {
            measure_array[i] = 0;
         }
      }
   } else {
      jxf_printf(" error in jxf_PAMGCoarsenRugeAI for CF_init = %d \n", CF_init);
      exit(0);
   }

  /*---------------------------------------------------
   * Loop until all points are either fine or coarse.
   *---------------------------------------------------*/

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();


   //================================//
   //    first coarsening phase      //
   //================================//

  /*----------------------------------------------------------------
   *   Initialize the lists
   *---------------------------------------------------------------*/

   if (CF_init == 9){
      CF_marker = *CF_marker_ptr;
      num_left = 0;
      for (i = 0; i < num_variables; i ++)
      {
         //if ( (S_i[i+1] - S_i[i]) == 0 && (S_offd_i[i+1] - S_offd_i[i]) == 0 )
	 if ( (S_i[i+1] - S_i[i]) == 0 && (S_offd_i[i+1]-S_offd_i[i]) == 0 && AI_measure[i] <= 0.1)
         {
            CF_marker[i] = SF_PT;
            measure_array[j] = 0;
         } else if (AI_measure[i] <= 0.1) {
            CF_marker[i] = UNDECIDED;
            num_left ++;
         }
      }
   } else {
      CF_marker = jxf_CTAlloc(JXF_Int, num_variables);
   
      num_left = 0;
      for (j = 0; j < num_variables; j ++)
      {
         //if ( (S_i[j+1] - S_i[j]) == 0 && (S_offd_i[j+1] - S_offd_i[j]) == 0 )
	 if ( (S_i[i+1] - S_i[i]) == 0 && (S_offd_i[i+1]-S_offd_i[i]) == 0 && AI_measure[i] <= 0.1)
         {
            CF_marker[j] = SF_PT;
            measure_array[j] = 0;
         }
         else
         {
            CF_marker[j] = UNDECIDED;
            num_left ++;
         }
      } 
   }

   for (j = 0; j < num_variables; j ++) 
   {    
      measure = measure_array[j];
      if (CF_marker[j] != SF_PT)
      {
         if (measure > 0)
         {
            jxf_enter_on_lists(&LoL_head, &LoL_tail, measure, j, lists, where);
         }
         else if (CF_init == 0)
         {
            if (measure < 0) jxf_printf("negative measure!\n");
            CF_marker[j] = f_pnt;
            for (k = S_i[j]; k < S_i[j+1]; k ++)
            {
               nabor = S_j[k];
               if (CF_marker[nabor] != SF_PT)
               {
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
            }
            --num_left;
         }
         else if (CF_init == 9)
         {
            if (measure < 0) jxf_printf("negative measure!\n");
            if (CF_marker[j] == f_pnt ) {
               for (k = S_i[j]; k < S_i[j+1]; k ++)
               {
                  nabor = S_j[k];
                  if (CF_marker[nabor] != SF_PT)
                  {
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
               }
            }
            --num_left;
         }
      }
   }

  
  /***********************************************************************************
   *
   *  Main loop of Ruge-Stueben first coloring pass.
   *
   *  WHILE there are still points to classify DO:
   *        1) find first point, i,  on list with max_measure
   *           make i a C-point, remove it from the lists
   *        2) For each point, j,  in S_i^T,
   *           a) Set j to be an F-point
   *           b) For each point, k, in S_j
   *                  move k to the list in LoL with measure one
   *                  greater than it occupies (creating new LoL
   *                  entry if necessary)
   *        3) For each point, j,  in S_i,
   *                  move j to the list in LoL with measure one
   *                  smaller than it occupies (creating new LoL
   *                  entry if necessary)
   *
   *************************************************************************************/

   while (num_left > 0)
   {
      index = LoL_head -> head;

      CF_marker[index] = C_PT;
      measure = measure_array[index];
      measure_array[index] = 0;
      --num_left;
      
      jxf_remove_point(&LoL_head, &LoL_tail, measure, index, lists, where);
  
      for (j = ST_i[index]; j < ST_i[index+1]; j ++)
      {
         nabor = ST_j[j];
         if (CF_marker[nabor] == UNDECIDED)
         {
            CF_marker[nabor] = F_PT;
            measure = measure_array[nabor];

            jxf_remove_point(&LoL_head, &LoL_tail, measure, nabor, lists, where);
            --num_left;

            for (k = S_i[nabor]; k < S_i[nabor+1]; k ++)
            {
               nabor_two = S_j[k];
               if (CF_marker[nabor_two] == UNDECIDED)
               {
                  measure = measure_array[nabor_two];
                  jxf_remove_point(&LoL_head, &LoL_tail, measure, nabor_two, lists, where);

                  new_meas = ++(measure_array[nabor_two]);
                 
                  jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas, nabor_two, lists, where);
               }
            }
         }
      }
      for (j = S_i[index]; j < S_i[index+1]; j ++)
      {
         nabor = S_j[j];
         if (CF_marker[nabor] == UNDECIDED)
         {
            measure = measure_array[nabor];

            jxf_remove_point(&LoL_head, &LoL_tail, measure, nabor, lists, where);

            measure_array[nabor] = --measure;
	
	    if (measure > 0)
	    {
               jxf_enter_on_lists(&LoL_head, &LoL_tail, measure, nabor, lists, where);
            }
	    else
	    {
               CF_marker[nabor] = F_PT;
               --num_left;

               for (k = S_i[nabor]; k < S_i[nabor+1]; k ++)
               {
                  nabor_two = S_j[k];
                  if (CF_marker[nabor_two] == UNDECIDED)
                  {
                     new_meas = measure_array[nabor_two];
                     jxf_remove_point(&LoL_head, &LoL_tail, new_meas, nabor_two, lists, where);

                     new_meas = ++(measure_array[nabor_two]);
                 
                     jxf_enter_on_lists(&LoL_head, &LoL_tail, new_meas, nabor_two, lists, where);
                  }
               }
	    }
         }
      }
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

   if (coarsen_type == 11)
   {
      *CF_marker_ptr = CF_marker;
      return 0;
   }


   //==============================================================//
   //   second pass, check fine points for coarse neighbors        //
   //   for coarsen_type = 2, the second pass includes             //
   //   off-processore boundary points                             // 
   //==============================================================//
   
  /*---------------------------------------------------
   * Initialize the graph array
   *---------------------------------------------------*/

   graph_array = jxf_CTAlloc(JXF_Int, num_variables);

   for (i = 0; i < num_variables; i ++)
   {
      graph_array[i] = -1;
   }

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();

   if (coarsen_type == 2)
   {
     /*------------------------------------------------
      * Exchange boundary data for CF_marker
      *------------------------------------------------*/
    
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg,num_sends));
    
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
         comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);
         jxf_ParCSRCommHandleDestroy(comm_handle);
      }
    
      ci_array = jxf_CTAlloc(JXF_Int, num_cols_offd);
      for (i = 0; i < num_cols_offd; i ++)
      {
         ci_array[i] = -1;
      }
	
      for (i = 0; i < num_variables; i ++)
      {
	 if (ci_tilde_mark != i) ci_tilde = -1;
	 if (ci_tilde_offd_mark != i) ci_tilde_offd = -1;
         if (CF_marker[i] == -1)
         {
            break_var = 1;
            for (ji = S_i[i]; ji < S_i[i+1]; ji ++)
            {
               j = S_j[ji];
               if (CF_marker[j] > 0)
               {
                  graph_array[j] = i;
               }
            }
            for (ji = S_offd_i[i]; ji < S_offd_i[i+1]; ji ++)
            {
               j = S_offd_j[ji];
               if (CF_marker_offd[j] > 0)
               {
                  ci_array[j] = i;
               }
            }
            for (ji = S_i[i]; ji < S_i[i+1]; ji ++)
            {
               j = S_j[ji];
               if (CF_marker[j] == -1)
               {
                  set_empty = 1;
                  for (jj = S_i[j]; jj < S_i[j+1]; jj ++)
                  {
                     index = S_j[jj];
                     if (graph_array[index] == i)
                     {
                        set_empty = 0;
                        break;
                     }
                  }
		  if (set_empty)
                  {
                     for (jj = S_offd_i[j]; jj < S_offd_i[j+1]; jj ++)
                     {
                        index = S_offd_j[jj];
                        if (ci_array[index] == i)
                        {
                           set_empty = 0;
                           break;
                        }
                     }
                  } 
                  if (set_empty)
                  {
                     if (C_i_nonempty)
                     {
                        CF_marker[i] = 1;
                        if (ci_tilde > -1)
                        {
                           CF_marker[ci_tilde] = -1;
                           ci_tilde = -1;
                        }
                        if (ci_tilde_offd > -1)
                        {
                           CF_marker_offd[ci_tilde_offd] = -1;
                           ci_tilde_offd = -1;
                        }
                        C_i_nonempty = 0;
                        break_var = 0;
                        break;
                     }
                     else
                     {
                        ci_tilde = j;
                        ci_tilde_mark = i;
                        CF_marker[j] = 1;
                        C_i_nonempty = 1;
                        i --;
                        break_var = 0;
                        break;
                     }
                  }
               }
            }
            if (break_var)
            {
               for (ji = S_offd_i[i]; ji < S_offd_i[i+1]; ji ++)
               {
                  j = S_offd_j[ji];
                  if (CF_marker_offd[j] == -1)
                  {
                     set_empty = 1;
                     for (jj = S_ext_i[j]; jj < S_ext_i[j+1]; jj ++)
                     {
                        index = S_ext_j[jj];
                        if (index > col_0 && index < col_n) /* index interior */
                        {
                           if (graph_array[index-first_col] == i)
                           {
                              set_empty = 0;
                              break;
                           }
                        }
                        else
                        {
   		           jk = jxf_BinarySearch(col_map_offd,index,num_cols_offd);
                           if (jk != -1)
                           {
                              if (ci_array[jk] == i)
                              {
                                 set_empty = 0;
                                 break;
                              }
                           }
                        }
                     }
                     if (set_empty)
                     {
                        if (C_i_nonempty)
                        {
                           CF_marker[i] = 1;
                           if (ci_tilde > -1)
                           {
                              CF_marker[ci_tilde] = -1;
                              ci_tilde = -1;
                           }
                           if (ci_tilde_offd > -1)
                           {
                              CF_marker_offd[ci_tilde_offd] = -1;
                              ci_tilde_offd = -1;
                           }
                           C_i_nonempty = 0;
                           break;
                        }
                        else
                        {
                           ci_tilde_offd = j;
                           ci_tilde_offd_mark = i;
                           CF_marker_offd[j] = 1;
                           C_i_nonempty = 1;
                           i --;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else
   {
      for (i = 0; i < num_variables; i ++)
      {
	 if (ci_tilde_mark != i) ci_tilde = -1;
         if (CF_marker[i] == -1)
         {
   	    for (ji = S_i[i]; ji < S_i[i+1]; ji ++)
   	    {
   	       j = S_j[ji];
   	       if (CF_marker[j] > 0)
   	       {
   	          graph_array[j] = i;
   	       }
    	    }
   	    for (ji = S_i[i]; ji < S_i[i+1]; ji ++)
   	    {
   	       j = S_j[ji];
   	       if (CF_marker[j] == -1)
   	       {
   	          set_empty = 1;
   	          for (jj = S_i[j]; jj < S_i[j+1]; jj ++)
   	          {
   		     index = S_j[jj];
   		     if (graph_array[index] == i)
   		     {
   		        set_empty = 0;
   		        break;
   		     }
   	          }
   	          if (set_empty)
   	          {
   		     if (C_i_nonempty)
   		     {
   		        CF_marker[i] = 1;
   		        if (ci_tilde > -1)
   		        {
   			   CF_marker[ci_tilde] = -1;
   		           ci_tilde = -1;
   		        }
   	    		C_i_nonempty = 0;
   		        break;
   		     }
   		     else
   		     {
   		        ci_tilde = j;
   		        ci_tilde_mark = i;
   		        CF_marker[j] = 1;
   		        C_i_nonempty = 1;
		        i --;
		        break;
		     }
	          }
	       }
	    }
	 }
      }
   }

   if (debug_flag == 3 && coarsen_type != 2)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d    Coarsen 2nd pass = %f\n", my_id, wall_time); 
   }



   //===============================================================//
   //  third pass, check boundary fine points for coarse neighbors  //
   //===============================================================//

   if (coarsen_type == 3 || coarsen_type == 4)
   {
      if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();

      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

     /*------------------------------------------------
      * Exchange boundary data for CF_marker
      *------------------------------------------------*/

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
         comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);
         jxf_ParCSRCommHandleDestroy(comm_handle);   
      }

      ci_array = jxf_CTAlloc(JXF_Int, num_cols_offd);
      for (i = 0; i < num_cols_offd; i ++)
      {
         ci_array[i] = -1;
      }
   }

   if (coarsen_type > 1 && coarsen_type < 5)
   { 
      for (i = 0; i < num_variables; i ++)
      {
         graph_array[i] = -1;
      }
      for (i = 0; i < num_cols_offd; i ++)
      {
	 if (ci_tilde_mark != i) ci_tilde = -1;
	 if (ci_tilde_offd_mark != i) ci_tilde_offd = -1;
         if (CF_marker_offd[i] == -1)
         {
   	    for (ji = S_ext_i[i]; ji < S_ext_i[i+1]; ji ++)
   	    {
   	       j = S_ext_j[ji];
   	       if (j > col_0 && j < col_n)
   	       {
   	          j = j - first_col;
   	          if (CF_marker[j] > 0)
   	          {
   	             graph_array[j] = i;
   	          }
   	       }
   	       else
   	       {
   		  jj = jxf_BinarySearch(col_map_offd, j, num_cols_offd);
   		  if (jj != -1 && CF_marker_offd[jj] > 0)
   		  {
   	             ci_array[jj] = i;
   	          }
    	       }	
    	    }
   	    for (ji = S_ext_i[i]; ji < S_ext_i[i+1]; ji ++)
   	    {
   	       j = S_ext_j[ji];
   	       if (j > col_0 && j < col_n)
   	       {
   	          j = j - first_col;
   	          if ( CF_marker[j] == -1)
   	          {
   	             set_empty = 1;
   	             for (jj = S_i[j]; jj < S_i[j+1]; jj ++)
   	             {
   		        index = S_j[jj];
   		        if (graph_array[index] == i)
   		        {
   		           set_empty = 0;
   		           break;
   		        }
   	             }
   	             for (jj = S_offd_i[j]; jj < S_offd_i[j+1]; jj ++)
   	             {
   		        index = S_offd_j[jj];
   		        if (ci_array[index] == i)
   		        {
   		           set_empty = 0;
   		           break;
   		        }
   	             }
   	             if (set_empty)
   	             {
   		        if (C_i_nonempty)
   		        {
   		           CF_marker_offd[i] = 1;
   		           if (ci_tilde > -1)
   		           {
   			      CF_marker[ci_tilde] = -1;
			      ci_tilde = -1;
   		           }
   		           if (ci_tilde_offd > -1)
   		           {
   			      CF_marker_offd[ci_tilde_offd] = -1;
			      ci_tilde_offd = -1;
   		           }
                           C_i_nonempty = 0;
   		           break;
   		        }
   		        else
   		        {
   		           ci_tilde = j;
   		           ci_tilde_mark = i;
   		           CF_marker[j] = 1;
   		           C_i_nonempty = 1;
   		           i --;
   		           break;
   		        }
   	             }
   	          }
   	       }
   	       else
   	       {
   		  jm = jxf_BinarySearch(col_map_offd, j, num_cols_offd);
   		  if (jm != -1 && CF_marker_offd[jm] == -1)
   	          {
   	             set_empty = 1;
   	             for (jj = S_ext_i[jm]; jj < S_ext_i[jm+1]; jj ++)
   	             {
   		        index = S_ext_j[jj];
   		        if (index > col_0 && index < col_n) 
   		  	{
   		           if (graph_array[index-first_col] == i)
   		           {
   		              set_empty = 0;
   		              break;
   		           }
   	                }
   			else
   			{
   		           jk = jxf_BinarySearch(col_map_offd, index, num_cols_offd);
   			   if (jk != -1)
   			   {
   		              if (ci_array[jk] == i)
   		              {
   		                 set_empty = 0;
   		                 break;
   		              }
   		           }
   	                }
   	             }
   	             if (set_empty)
   	             {
   		        if (C_i_nonempty)
   		        {
   		           CF_marker_offd[i] = 1;
   		           if (ci_tilde > -1)
   		           {
   			      CF_marker[ci_tilde] = -1;
   			      ci_tilde = -1;
   		           }
   		           if (ci_tilde_offd > -1)
   		           {
   			      CF_marker_offd[ci_tilde_offd] = -1;
   			      ci_tilde_offd = -1;
   		           }
                           C_i_nonempty = 0;
   		           break;
   		        }
   		        else
   		        {
   		           ci_tilde_offd = jm;
   		           ci_tilde_offd_mark = i;
   		           CF_marker_offd[jm] = 1;
   		           C_i_nonempty = 1;
   		           i --;
   		           break;
   		        }
   		     }
   	          }
   	       }
   	    }
         }
      }
      
     /*---------------------------------------------------------------------------------
      * Send boundary data for CF_marker back
      *--------------------------------------------------------------------------------*/
      if (num_procs > 1)
      {
         comm_handle = jxf_ParCSRCommHandleCreate(12, comm_pkg, CF_marker_offd, int_buf_data);
         jxf_ParCSRCommHandleDestroy(comm_handle);   
      }
   
      
     /* only CF_marker entries from larger procs are accepted  
	if coarsen_type = 4 coarse points are not overwritten  */
 
      index = 0;
      if (coarsen_type != 4)
      {
         for (i = 0; i < num_sends; i ++)
         {
	    start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
            if (jxf_ParCSRCommPkgSendProc(comm_pkg,i) > my_id)
	    {
               for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
               {
                  CF_marker[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)] = int_buf_data[index++]; 
               }
            }
	    else
	    {
	       index += jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - start;
	    }
         }
      }
      else
      {
         for (i = 0; i < num_sends; i ++)
         {
	    start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
            if (jxf_ParCSRCommPkgSendProc(comm_pkg,i) > my_id)
	    {
              for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
              {
                 elmt = jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j);
                 if (CF_marker[elmt] != 1)
                 {
                    CF_marker[elmt] = int_buf_data[index];
                 }
		 index ++; 
              }
            }
	    else
	    {
	       index += jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - start;
	    }
         }
      }
      
      if (debug_flag == 3)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         if (coarsen_type == 4)
         {
            jxf_printf("Proc = %d    Coarsen 3rd pass = %f\n", my_id, wall_time); 
         }
         if (coarsen_type == 3)
         {
            jxf_printf("Proc = %d    Coarsen 3rd pass = %f\n", my_id, wall_time); 
         }
         if (coarsen_type == 2)
         {
            jxf_printf("Proc = %d    Coarsen 2nd pass = %f\n", my_id, wall_time);
         } 
      }
   }

   if (coarsen_type == 5)
   {
     /*------------------------------------------------
      * Exchange boundary data for CF_marker
      *------------------------------------------------*/

      if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();
    
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
    
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
         comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);   
         jxf_ParCSRCommHandleDestroy(comm_handle);
      }
    
      ci_array = jxf_CTAlloc(JXF_Int, num_cols_offd);
      for (i = 0; i < num_cols_offd; i ++)
      {
         ci_array[i] = -1;
      }
      for (i = 0; i < num_variables; i ++)
      {
         graph_array[i] = -1;
      }

      for (i = 0; i < num_variables; i ++)
      {
         if (CF_marker[i] == -1 && (S_offd_i[i+1]-S_offd_i[i]) > 0)
         {
            break_var = 1;
            for (ji = S_i[i]; ji < S_i[i+1]; ji ++)
            {
               j = S_j[ji];
               if (CF_marker[j] > 0)
               {
                  graph_array[j] = i;
               }
            }
            for (ji = S_offd_i[i]; ji < S_offd_i[i+1]; ji ++)
            {
               j = S_offd_j[ji];
               if (CF_marker_offd[j] > 0)
               {
                  ci_array[j] = i;
               }
            }
            for (ji = S_offd_i[i]; ji < S_offd_i[i+1]; ji ++)
            {
               j = S_offd_j[ji];
               if (CF_marker_offd[j] == -1)
               {
                  set_empty = 1;
                  for (jj = S_ext_i[j]; jj < S_ext_i[j+1]; jj ++)
                  {
                     index = S_ext_j[jj];
                     if (index > col_0 && index < col_n) /* index interior */
                     {
                        if (graph_array[index-first_col] == i)
                        {
                           set_empty = 0;
                           break;
                        }
                     }
                     else
                     {
   		        jk = jxf_BinarySearch(col_map_offd, index, num_cols_offd);
                        if (jk != -1)
                        {
                           if (ci_array[jk] == i)
                           {
                              set_empty = 0;
                              break;
                           }
                        }
                     }
                  }
                  if (set_empty)
                  {
                     if (C_i_nonempty)
                     {
                        CF_marker[i] = -2;
                        C_i_nonempty = 0;
                        break;
                     }
                     else
                     {
                        C_i_nonempty = 1;
                        i --;
                        break;
                     }
                  }
               }
            }
         }
      }
      
      if (debug_flag == 3)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         jxf_printf("Proc = %d    Coarsen special points = %f\n", my_id, wall_time); 
      }

   }
   
   
  /*---------------------------------------------------
   * Clean up and return
   *---------------------------------------------------*/

   if (coarsen_type != 1)
   {   
      jxf_TFree(CF_marker_offd);
      jxf_TFree(int_buf_data);
      jxf_TFree(ci_array);
   }   
   jxf_TFree(graph_array);
   if ( (measure_type || (coarsen_type != 1 && coarsen_type != 11)) && num_procs > 1 )
   {
      jxf_CSRMatrixDestroy(S_ext);
   } 
   
   *CF_marker_ptr = CF_marker;
   
   return (ierr);
}
