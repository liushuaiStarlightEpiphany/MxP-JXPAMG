//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_coarsen_pmis.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

#define C_PT  1
#define F_PT -1
#define SF_PT -3
#define COMMON_C_PT  2
#define Z_PT -2

/*!
 * \fn JXF_Int jxf_PAMGCoarsenPMIS
 * \brief PMIS (Modified Independent Set ) coarsening routine.
 * \note Don't worry about strong F-F connections without a common C point.
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGCoarsenPMIS( jxf_ParCSRMatrix  *par_S,
                    jxf_ParCSRMatrix  *par_A,
                    JXF_Int               CF_init,
                    JXF_Int               debug_flag,
                    JXF_Int             **CF_marker_ptr )
{
   MPI_Comm               comm        = jxf_ParCSRMatrixComm(par_S);
   jxf_ParCSRCommPkg      *comm_pkg    = jxf_ParCSRMatrixCommPkg(par_S);
   jxf_ParCSRCommHandle   *comm_handle = NULL;

   jxf_CSRMatrix    *S_diag    = jxf_ParCSRMatrixDiag(par_S);
   JXF_Int             *S_diag_i  = jxf_CSRMatrixI(S_diag);
   JXF_Int             *S_diag_j  = jxf_CSRMatrixJ(S_diag);

   jxf_CSRMatrix    *S_offd    = jxf_ParCSRMatrixOffd(par_S);
   JXF_Int             *S_offd_i  = jxf_CSRMatrixI(S_offd);
   JXF_Int             *S_offd_j  = NULL;

   JXF_Int              num_variables = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int 		    num_cols_offd = 0;

   JXF_Int		    num_sends = 0;
   JXF_Int  	   *int_buf_data;
   JXF_Real	   *buf_data;

   JXF_Int             *CF_marker;
   JXF_Int             *CF_marker_offd;
                      
   JXF_Real          *measure_array;
   JXF_Int             *graph_array;
   JXF_Int             *graph_array_offd;
   JXF_Int              graph_size;
   JXF_Int              graph_offd_size;
   JXF_Int              global_graph_size;
                      
   JXF_Int              i, j, jS, ig;
   JXF_Int              index, start, my_id, num_procs, jrow, cnt, elmt;
                      
   JXF_Int              ierr = 0;

   JXF_Real           wall_time = 0.0;
   JXF_Int              iter = 0;


  /***********************************************************************
   *  BEFORE THE INDEPENDENT SET COARSENING LOOP:
   *   measure_array: calculate the measures, and communicate them
   *     (this array contains measures for both local and external nodes)
   *   CF_marker, CF_marker_offd: initialize CF_marker
   *     (separate arrays for local and external; 
   *      0=unassigned, negative=F point, positive=C point)
   ***********************************************************************/      

  /*-------------------------------------------------------------------------
   * Use the ParCSR strength matrix, par_S.
   *
   * For now, the "strength" of dependence/influence is defined in
   * the following way: i depends on j if
   *     aij > jxf_max (k != i) aik,    aii < 0
   * or
   *     aij < jxf_min (k != i) aik,    aii >= 0
   * Then S_ij = 1, else S_ij = 0.
   *
   * NOTE: S_data is not used; in stead, only strong columns are retained
   *       in S_j, which can then be used like S_data
   *-----------------------------------------------------------------------*/

   if (debug_flag == 3) wall_time = jxf_time_getWallclockSeconds();
   
   jxf_MPI_Comm_size(comm,&num_procs);
   jxf_MPI_Comm_rank(comm,&my_id);

   if (!comm_pkg)
   {
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
   }

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);

   int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   buf_data     = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
 
   num_cols_offd = jxf_CSRMatrixNumCols(S_offd);

   S_diag_j = jxf_CSRMatrixJ(S_diag);

   if (num_cols_offd)
   {
      S_offd_j = jxf_CSRMatrixJ(S_offd);
   }

  /*---------------------------------------------------------------
   * Compute the measures
   *
   * The measures are currently given by the column sums of par_S.
   * Hence, measure_array[i] is the number of influences
   * of variable i.
   *
   * The measures are augmented by a random number
   * between 0 and 1.
   *-------------------------------------------------------------*/

   measure_array = jxf_CTAlloc(JXF_Real, num_variables + num_cols_offd);

   /* first calculate the local part of the sums for the external nodes */
   for (i = 0; i < S_offd_i[num_variables]; i ++)
   { 
      measure_array[num_variables + S_offd_j[i]] += 1.0;
   }

   /* now send those locally calculated values for the external nodes to the neighboring processors */
   if (num_procs > 1)
   {
      comm_handle = jxf_ParCSRCommHandleCreate(2, comm_pkg, &measure_array[num_variables], buf_data);
   }

   /* calculate the local part for the local nodes */
   for (i = 0; i < S_diag_i[num_variables]; i ++)
   { 
      measure_array[S_diag_j[i]] += 1.0;
   }

   /* finish the communication */
   if (num_procs > 1)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle);
   }
       
   /* now add the externally calculated part of the local nodes to the local nodes */
   index = 0;
   for (i = 0; i < num_sends; i ++)
   {
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
      {
         measure_array[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)] += buf_data[index++];
      }
   }

   /* set the measures of the external nodes to zero */
   for (i = num_variables; i < num_variables + num_cols_offd; i ++)
   { 
      measure_array[i] = 0;
   }

   /* this augments the measures with a random number between 0 and 1 */
   /* (only for the local part) */
   /* this augments the measures */
   if (CF_init == 2)
   {
      jxf_PAMGIndepSetInit(par_S, measure_array, 1);
   }
   else
   {
      jxf_PAMGIndepSetInit(par_S, measure_array, 0);
   }

   /*for (i = 0; i < num_variables; i ++)
   { 
      jxf_printf("i= %d,     meausre= %f\n", i, measure_array[i]);
   }*/

  /*---------------------------------------------------
   * Initialize the graph arrays, and CF_marker arrays
   *---------------------------------------------------*/

   /* first the off-diagonal part of the graph array */
   if (num_cols_offd) 
   {
      graph_array_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
   }
   else
   {
      graph_array_offd = NULL;
   }

   for (ig = 0; ig < num_cols_offd; ig ++)
   {
      graph_array_offd[ig] = ig;
   }

   graph_offd_size = num_cols_offd;

   /* now the local part of the graph array, and the local CF_marker array */
   graph_array = jxf_CTAlloc(JXF_Int, num_variables);

   if (CF_init == 1)
   { 
      CF_marker = *CF_marker_ptr;
      cnt = 0;
      for (i = 0; i < num_variables; i ++)
      {
         if ( (S_offd_i[i+1] - S_offd_i[i]) > 0 || CF_marker[i] == -1)
	 {
            CF_marker[i] = 0;
	 }
         if (CF_marker[i] == Z_PT)
         {
            if (measure_array[i] >= 1.0 || (S_diag_i[i+1]-S_diag_i[i]) > 0)
            {
               CF_marker[i] = 0;
               graph_array[cnt++] = i;
            }
            else
            {
               CF_marker[i] = F_PT;
            }
         }
         else if (CF_marker[i] == SF_PT)
         {
            measure_array[i] = 0;
         }
         else
         {
            graph_array[cnt++] = i;
         }
      }
   }
   else
   {
      CF_marker = jxf_CTAlloc(JXF_Int, num_variables);
      cnt = 0;
      for (i = 0; i < num_variables; i ++)
      {
         CF_marker[i] = 0;
         if ( (S_diag_i[i+1] - S_diag_i[i]) == 0 && (S_offd_i[i+1] - S_offd_i[i]) == 0)
         {
            CF_marker[i] = SF_PT;
            measure_array[i] = 0;
         }
         else
         {
            graph_array[cnt++] = i;
         }
      }
   }
   graph_size = cnt;

   /* now the off-diagonal part of CF_marker */
   if (num_cols_offd)
   {
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
   }
   else
   {
      CF_marker_offd = NULL;
   }
   
   for (i = 0; i < num_cols_offd; i ++)
   {
      CF_marker_offd[i] = 0;
   }
  
  
  /*--------------------------------------------------------
   * Communicate the local measures, which are complete,
   * to the external nodes
   *------------------------------------------------------*/
   
   index = 0;
   for (i = 0; i < num_sends; i ++)
   {
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
      {
         jrow = jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j);
         buf_data[index++] = measure_array[jrow];
      }
   }
   
   if (num_procs > 1)
   { 
      comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, buf_data, &measure_array[num_variables]);
      jxf_ParCSRCommHandleDestroy(comm_handle);   
   } 
 
   if (debug_flag == 3)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d    Initialize CLJP phase = %f\n", my_id, wall_time); 
   }

  /*****************************************************************************
   * THE INDEPENDENT SET COARSENING LOOP:
   *****************************************************************************/      

  /*----------------------------------------------------
   * Loop until all points are either fine or coarse.
   *---------------------------------------------------*/
   //choose allreduce
   #if JXF_hoAllreduce
      #define choosen 2
   #else
      #define choosen 3
   #endif

   #if (choosen == 2 || choosen == 3 || choosen == 4)
      global_graph_size = 1;
      MPI_Request *req = jxf_CTAlloc(MPI_Request, 1);
      MPI_Status  *status = jxf_CTAlloc(MPI_Status, 1);
   #endif

   #if (choosen == 2)
      MPI_Comm  Node_Comm  = jxf_hpNodeComm(jxf_hp_hardware);
      JXF_Int tmp_graph_size;
   #endif

   while (1)
   {
      #if (choosen == 2 || choosen ==  3 || choosen == 4)
         if (global_graph_size == 0)
         {
            break;
         }
      #endif
      /* stop the coarsening if nothing left to be coarsened */
      #if (choosen == 0)
         jxf_MPI_Allreduce(&graph_size, &global_graph_size, 1, JXF_MPI_INT, MPI_SUM, comm);
      #elif (choosen == 1) 
         jxf_My_MPI_Allreduce(&graph_size, &global_graph_size, 1, JXF_MPI_INT, MPI_SUM, comm);
      #elif (choosen == 2)
         jxf_My_MPI_Iallreduce_First(&graph_size, &global_graph_size, 1, JXF_MPI_INT, MPI_SUM, comm);
         if(Node_Comm != MPI_COMM_NULL)  MPI_Iallreduce(&global_graph_size, &tmp_graph_size, 1, JXF_MPI_INT, MPI_SUM, Node_Comm, req);
      #elif (choosen == 3)
         MPI_Iallreduce(&graph_size, &global_graph_size, 1, JXF_MPI_INT, MPI_SUM, comm, req);
      #elif (choosen == 4)
         MPI_Ireduce(&graph_size, &global_graph_size, 1, JXF_MPI_INT, MPI_SUM, num_procs - 1, comm, req);
      #endif

      #if JXF_PRINT_MODULE_TIME && JXF_USE_HPCTOOLKIT == 0
         jxf_sprintf(timetag,"%s_l%d", "Allreduce",TimerLevel);
         TimerManagerEnd(timetag);
      #endif

      #if (choosen == 0 || choosen == 1)
         if (global_graph_size == 0)
         {
            break;
         }
      #endif

     /*----------------------------------------------------------------
      * Pick an independent set of points with
      * maximal measure.
      * At the end, CF_marker is complete, but still needs to be
      * communicated to CF_marker_offd
      *---------------------------------------------------------------*/
      if (!CF_init || iter)
      {
         jxf_PAMGIndepSet( par_S, measure_array, graph_array, 
                          graph_size, 
                          graph_array_offd, graph_offd_size, 
                          CF_marker, CF_marker_offd );

      /*-------------------------------------------------------------
       * Exchange boundary data for CF_marker: send internal
       *  points to external points
       *------------------------------------------------------------*/

         if (num_procs > 1)
         {
            comm_handle = jxf_ParCSRCommHandleCreate(12, comm_pkg, CF_marker_offd, int_buf_data);
            jxf_ParCSRCommHandleDestroy(comm_handle);   
         }

         index = 0;
         for (i = 0; i < num_sends; i ++)
         {
            start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
            for (j=start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
            {
               elmt = jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j);
               if (!int_buf_data[index] && CF_marker[elmt] > 0)
               {
                  CF_marker[elmt] = 0; 
                  index ++;
               }
               else
               {
                  int_buf_data[index++] = CF_marker[elmt];
               }
            }
         }
         
 
         if (num_procs > 1)
         {
            comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);
            jxf_ParCSRCommHandleDestroy(comm_handle);   
         }
         
      } // end (!CF_init || iter)

      iter ++;


     /*------------------------------------------------
      * Set C-pts and F-pts.
      *------------------------------------------------*/

      for (ig = 0; ig < graph_size; ig ++) 
      {
      
         i = graph_array[ig];

        /*------------------------------------------------------
	 * If the measure of i is smaller than 1, then
         * make i and F point (because it does not influence
         * any other point), and remove all edges of
	 * equation i.
	 *-----------------------------------------------------*/

	 if (measure_array[i] < 1.0)
	 {
	    /* make point i an F point*/
	    CF_marker[i]= F_PT;

	    /* remove the edges in equation i */
	    /* first the local part */
	    for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS++) 
	    {
	       j = S_diag_j[jS];
	       if (j > -1)
	       { 
	          /* column number is still positive; not accounted for yet */
	          S_diag_j[jS]  = -S_diag_j[jS]-1;
	       }
	    }
	    
	    /* now the external part */
	    for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS++) 
	    {
	       j = S_offd_j[jS];
	       if (j > -1)
	       { 
	          /* column number is still positive; not accounted for yet */
	          S_offd_j[jS]  = -S_offd_j[jS]-1;
	       }
	    }
	    
	 } // end if (measure_array[i] < 1.0)


        /*--------------------------------------------------
         * First treat the case where point i is in the
         * independent set: make i a C point, 
         * take out all the graph edges for
         * equation i.
         *------------------------------------------------*/
       
         if (CF_marker[i] > 0) 
         {
            /* set to be a C-pt */
	    CF_marker[i] = C_PT;

            /* remove the edges in equation i */
	    /* first the local part */
	    for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS ++) 
	    {
               j = S_diag_j[jS];
               if (j > -1)
               { 
                  /* column number is still positive; not accounted for yet */
	          S_diag_j[jS]  = - S_diag_j[jS] - 1;
               }
            }
         
	    /* now the external part */
	    for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS ++) 
	    {
               j = S_offd_j[jS];
               if (j > -1)
               { 
                  /* column number is still positive; not accounted for yet */
                  S_offd_j[jS]  = - S_offd_j[jS] - 1;
               }
            }
         
         }


        /*-------------------------------------------------------
         * Now treat the case where point i is not in the
         * independent set: loop over
         * all the points j that influence equation i; if
         * j is a C point, then make i an F point.
         * If i is a new F point, then remove all the edges
         * from the graph for equation i.
         *-------------------------------------------------------*/

         else 
         {

            /* first the local part */
            for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS ++) 
            {
               /* j is the column number, or the local number of the point influencing i */
               j = S_diag_j[jS];
               
               if (j < 0) 
               {
                  j = - j - 1;
               }

               if (CF_marker[j] > 0)
               { 
                  /* j is a C-point */
                  CF_marker[i] = F_PT;
               }
            }
            
            /* now the external part */
            for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS ++) 
            {
               j = S_offd_j[jS];
               
               if (j < 0)
               {
                  j = - j - 1;
               }
               
               if (CF_marker_offd[j] > 0)
               { 
                  /* j is a C-point */
                  CF_marker[i] = F_PT;
               }
            }

            /* remove all the edges for equation i if i is a new F point */
            if (CF_marker[i] == F_PT)
            {
               /* first the local part */
               for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS ++) 
               {
                  j = S_diag_j[jS];
                  if (j > -1)
                  { 
                     /* column number is still positive; not accounted for yet */
                     S_diag_j[jS]  = - S_diag_j[jS] - 1;
                  }
               }
               
               /* now the external part */
               for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS ++) 
               {
                  j = S_offd_j[jS];
                  if (j > -1)
                  { 
                     /* column number is still positive; not accounted for yet */
                     S_offd_j[jS]  = -S_offd_j[jS] - 1;
                  }
               }
            }   
	    
         } // end if (CF_marker[i] > 0)
       
      } /* end first loop over graph */


     /*****************************************************************
      * now communicate CF_marker to CF_marker_offd, to make sure
      * that new external F points are known on this processor 
      *****************************************************************/

     /*------------------------------------------------
      * Exchange boundary data for CF_marker: send internal
      *   points to external points
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


     /*-------------------------------------------------------
      * Now loop over the points i in the unassigned
      * graph again. For all points i that are no new C or
      * F points, remove the edges in equation i that
      * connect to C or F points.
      * (We have removed the rows for the new C and F
      * points above; now remove the columns.)
      *------------------------------------------------------*/

      for (ig = 0; ig < graph_size; ig ++) 
      {
         i = graph_array[ig];

         if (CF_marker[i] == 0) 
         {
            /* first the local part */
            for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS ++) 
            {
               j = S_diag_j[jS];
               
               if (j < 0) 
               {
                  j = - j - 1;
               }

               if (!CF_marker[j] == 0 && S_diag_j[jS] > -1)
               { 
                  /* connection to C or F point, and column number is still positive; not accounted for yet */
                  S_diag_j[jS]  = -S_diag_j[jS]-1;
               }
            }
            
            /* now the external part */
            for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS ++) 
            {
               j = S_offd_j[jS];
               
               if (j < 0) 
               {
                  j = - j - 1;
               }

               if (!CF_marker_offd[j]==0 && S_offd_j[jS] > -1)
               { 
                  /* connection to C or F point, and column number is still positive; not accounted for yet */
                  S_offd_j[jS]  = -S_offd_j[jS]-1;
               }
            }
            
         } // end if (CF_marker[i] == 0) 
         
      } /* end second loop over graph */


     /*------------------------------------------------
      * Update subgraph
      *------------------------------------------------*/

      for (ig = 0; ig < graph_size; ig ++) 
      {
         i = graph_array[ig];
       
         if (!CF_marker[i]==0) /* C or F point */
	 {
            /* the independent set subroutine needs measure 0 for removed nodes */
            measure_array[i] = 0;
            /* take point out of the subgraph */
            graph_size --;
            graph_array[ig] = graph_array[graph_size];
            graph_array[graph_size] = i;
            ig --;
         }
      }
     
      for (ig = 0; ig < graph_offd_size; ig ++) 
      {
         i = graph_array_offd[ig];
       
         if (!CF_marker_offd[i] == 0) /* C or F point */
         {
            /* the independent set subroutine needs measure 0 for removed nodes */
            measure_array[i+num_variables] = 0;
            /* take point out of the subgraph */
            graph_offd_size --;
            graph_array_offd[ig] = graph_array_offd[graph_offd_size];
            graph_array_offd[graph_offd_size] = i;
            ig --;
         }
      }
   #if (choosen == 2)
      if(Node_Comm != MPI_COMM_NULL)  MPI_Wait(req, status);
      jxf_My_MPI_Iallreduce_Sencond(&tmp_graph_size, 1, MPI_INT, MPI_SUM, comm, req);
      global_graph_size = tmp_graph_size;
   #elif (choosen == 3)
      MPI_Wait(req, status);
   #elif (choosen == 4)
      MPI_Wait(req, status);
      jxf_my_bcast(&global_graph_size, 1, MPI_INT, num_procs - 1, comm);
   #endif
   } /* end while */

#if (choosen == 2 || choosen == 3 || choosen == 4)
   jxf_TFree(req);
   jxf_TFree(status);
#endif

  /*---------------------------------------------------
   * Clean up and return
   *---------------------------------------------------*/

   /* Reset S_matrix */
   for (i = 0; i < S_diag_i[num_variables]; i ++)
   {
      if (S_diag_j[i] < 0)
      {
         S_diag_j[i] = - S_diag_j[i] - 1;
      }
   }
   
   for (i = 0; i < S_offd_i[num_variables]; i ++)
   {
      if (S_offd_j[i] < 0)
      {
         S_offd_j[i] = - S_offd_j[i] - 1;
      }
   }
   
   /*
   for (i = 0; i < num_variables; i ++)
   {
      if (CF_marker[i] == SF_PT) 
      {
         CF_marker[i] = F_PT;
      }
   }
   */

   jxf_TFree(measure_array);
   jxf_TFree(graph_array);
   if (num_cols_offd) 
   {
      jxf_TFree(graph_array_offd);
   }
   jxf_TFree(buf_data);
   jxf_TFree(int_buf_data);
   jxf_TFree(CF_marker_offd);
   
   /*
   if (num_procs > 1) 
   {
      jxf_CSRMatrixDestroy(S_ext);
   }
   */

   /*for (i = 0; i < num_variables; i ++)
   { 
      jxf_printf("i= %d,     CF_marker= %d\n", i, CF_marker[i]);
   }*/

   *CF_marker_ptr = CF_marker;

   return (ierr);
}
