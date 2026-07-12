//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_interp_3.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGBuildDirInterp
 * \brief Direct Interpolation.
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGBuildDirInterp( jxf_ParCSRMatrix   *par_A,
                       JXF_Int               *CF_marker,
                       jxf_ParCSRMatrix   *par_S,
                       JXF_Int               *num_cpts_global,
                       JXF_Int                num_functions,
                       JXF_Int               *dof_func,
                       JXF_Int                debug_flag,
                       JXF_Real             trunc_factor,
                       JXF_Int                max_elmts,
                       JXF_Int               *col_offd_S_to_A,
                       jxf_ParCSRMatrix  **P_ptr )
{
   MPI_Comm 	         comm     = jxf_ParCSRMatrixComm(par_A);   
   jxf_ParCSRCommPkg     *comm_pkg = jxf_ParCSRMatrixCommPkg(par_A);
   jxf_ParCSRCommHandle  *comm_handle;

   jxf_CSRMatrix *A_diag      = jxf_ParCSRMatrixDiag(par_A);
   JXF_Real       *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int          *A_diag_i    = jxf_CSRMatrixI(A_diag);
   JXF_Int          *A_diag_j    = jxf_CSRMatrixJ(A_diag);

   jxf_CSRMatrix *A_offd          = jxf_ParCSRMatrixOffd(par_A);   
   JXF_Real       *A_offd_data     = jxf_CSRMatrixData(A_offd);
   JXF_Int          *A_offd_i        = jxf_CSRMatrixI(A_offd);
   JXF_Int          *A_offd_j        = jxf_CSRMatrixJ(A_offd);
   JXF_Int           num_cols_A_offd = jxf_CSRMatrixNumCols(A_offd);

   jxf_CSRMatrix *S_diag   = jxf_ParCSRMatrixDiag(par_S);
   JXF_Int          *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int          *S_diag_j = jxf_CSRMatrixJ(S_diag);

   jxf_CSRMatrix *S_offd   = jxf_ParCSRMatrixOffd(par_S);   
   JXF_Int          *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int          *S_offd_j = jxf_CSRMatrixJ(S_offd);

   jxf_ParCSRMatrix *par_P;
   JXF_Int		   *col_map_offd_P = NULL;

   JXF_Int             *CF_marker_offd = NULL;
   JXF_Int             *dof_func_offd = NULL;

   jxf_CSRMatrix    *P_diag;
   jxf_CSRMatrix    *P_offd;   

   JXF_Real          *P_diag_data;
   JXF_Int             *P_diag_i;
   JXF_Int             *P_diag_j;
   JXF_Real          *P_offd_data;
   JXF_Int             *P_offd_i;
   JXF_Int             *P_offd_j;

   JXF_Int              P_diag_size, P_offd_size;
   
   JXF_Int             *P_marker, *P_marker_offd;

   JXF_Int              jj_counter, jj_counter_offd;
   JXF_Int             *jj_count, *jj_count_offd;
   JXF_Int              jj_begin_row, jj_begin_row_offd;
   JXF_Int              jj_end_row, jj_end_row_offd;
   
   JXF_Int              start_indexing = 0; /* start indexing for P_data at 0 */

   JXF_Int              n_fine = jxf_CSRMatrixNumRows(A_diag);

   JXF_Int             *fine_to_coarse;
   JXF_Int             *fine_to_coarse_offd;
   JXF_Int             *coarse_counter;
   JXF_Int              coarse_shift;
   JXF_Int              total_global_cpts;
   JXF_Int              num_cols_P_offd, my_first_cpt;

   JXF_Int              i, i1;
   JXF_Int              j, jl, jj;
   JXF_Int              start;
   
   JXF_Real           diagonal;
   JXF_Real           sum_N_pos, sum_P_pos;
   JXF_Real           sum_N_neg, sum_P_neg;
   JXF_Real           alfa = 1.0;
   JXF_Real           beta = 1.0;
   
   JXF_Real           zero = 0.0;
   JXF_Real           one  = 1.0;
   
   JXF_Int              my_id;
   JXF_Int              num_procs;
   JXF_Int              num_threads;
   JXF_Int              num_sends;
   JXF_Int              index;
   JXF_Int              ns, ne, size, rest;
   JXF_Int             *int_buf_data;

   JXF_Real           wall_time = 0.0;  /* for debugging instrumentation  */

   jxf_MPI_Comm_size(comm, &num_procs);   
   jxf_MPI_Comm_rank(comm, &my_id);
   num_threads = jxf_NumThreads();

#ifdef JXF_NO_GLOBAL_PARTITION
   my_first_cpt = num_cpts_global[0];
   if (my_id == (num_procs - 1)) 
   {
      total_global_cpts = num_cpts_global[1];
   }
   jxf_MPI_Bcast(&total_global_cpts, 1, JXF_MPI_INT, num_procs - 1, comm);
#else
   my_first_cpt = num_cpts_global[my_id];
   total_global_cpts = num_cpts_global[num_procs];
#endif

  /*-------------------------------------------------------------------
   * Get the CF_marker data for the off-processor columns
   *-------------------------------------------------------------------*/

   if (debug_flag == 4) wall_time = jxf_time_getWallclockSeconds();

   if (num_cols_A_offd) 
   {
      CF_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);
   }
   
   if (num_functions > 1 && num_cols_A_offd)
   {
      dof_func_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);
   }

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_A); 
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
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
	
   comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, CF_marker_offd);
   jxf_ParCSRCommHandleDestroy(comm_handle);
      
   if (num_functions > 1)
   {
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
	 start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	 for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
	 {
            int_buf_data[index++] = dof_func[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
         }
      }
	
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, dof_func_offd);
      jxf_ParCSRCommHandleDestroy(comm_handle);   
   }

   if (debug_flag == 4)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d     Interp: Comm 1 CF_marker =    %f\n", my_id, wall_time);
      fflush(NULL);
   }

  /*-----------------------------------------------------------------------
   *  First Pass: Determine size of P and fill in fine_to_coarse mapping.
   *-----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------
   *  Intialize counters and allocate mapping vector.
   *-----------------------------------------------------------------------*/

   coarse_counter = jxf_CTAlloc(JXF_Int, num_threads);
   jj_count = jxf_CTAlloc(JXF_Int, num_threads);
   jj_count_offd = jxf_CTAlloc(JXF_Int, num_threads);
   fine_to_coarse = jxf_CTAlloc(JXF_Int, n_fine);
   
#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
   for (i = 0; i < n_fine; i ++) 
   {
      fine_to_coarse[i] = -1;
   }

   jj_counter = start_indexing;
   jj_counter_offd = start_indexing;
      
  /*-----------------------------------------------------------------------
   *  Loop over fine grid.
   *-----------------------------------------------------------------------*/

  /* RDF: this looks a little tricky, but doable */
  
#define JXF_SMP_PRIVATE i,j,i1,jj,ns,ne,size,rest
#include "../../include/jxf_smp_forloop.h"
  
   for (j = 0; j < num_threads; j ++)
   {
      size = n_fine / num_threads;
      rest = n_fine - size*num_threads;
      if (j < rest)
      {
         ns = j*size + j;
         ne = (j+1)*size + j + 1;
      }
      else
      {
         ns = j*size + rest;
         ne = (j+1)*size + rest;
      }
      
      for (i = ns; i < ne; i ++)
      {
      
        /*--------------------------------------------------------------------
         *  If i is a C-point, interpolation is the identity. Also set up
         *  mapping vector.
         *--------------------------------------------------------------------*/

         if (CF_marker[i] >= 0)
         {
            jj_count[j] ++;
            fine_to_coarse[i] = coarse_counter[j];
            coarse_counter[j] ++;
         }
      
        /*--------------------------------------------------------------------
         *  If i is an F-point, interpolation is from the C-points that
         *  strongly influence i.
         *--------------------------------------------------------------------*/

         else
         {
            for (jj = S_diag_i[i]; jj < S_diag_i[i+1]; jj ++)
            {
               i1 = S_diag_j[jj];           
               if (CF_marker[i1] > 0)
               {
                  jj_count[j] ++;
               }
            }

            if (num_procs > 1)
            {
               if (col_offd_S_to_A)
               {
                  for (jj = S_offd_i[i]; jj < S_offd_i[i+1]; jj++)
                  {
                     i1 = col_offd_S_to_A[S_offd_j[jj]];           
                     if (CF_marker_offd[i1] > 0)
                     {
                        jj_count_offd[j] ++;
                     }
                  }
               }
               else
               {
                  for (jj = S_offd_i[i]; jj < S_offd_i[i+1]; jj ++)
                  {
                     i1 = S_offd_j[jj];           
                     if (CF_marker_offd[i1] > 0)
                     {
                        jj_count_offd[j] ++;
                     }
                  }
               }
            } // end if (num_procs > 1)
         } // end if
      } // end for i
   } // end for j



  /*-----------------------------------------------------------------------
   *  Allocate  arrays.
   *-----------------------------------------------------------------------*/

   for (i = 0; i < num_threads - 1; i ++)
   {
      coarse_counter[i+1] += coarse_counter[i];
      jj_count[i+1] += jj_count[i];
      jj_count_offd[i+1] += jj_count_offd[i];
   }
   i = num_threads - 1;
   jj_counter = jj_count[i];
   jj_counter_offd = jj_count_offd[i];

   P_diag_size = jj_counter;

   P_diag_i    = jxf_CTAlloc(JXF_Int, n_fine+1);
   P_diag_j    = jxf_CTAlloc(JXF_Int, P_diag_size);
   P_diag_data = jxf_CTAlloc(JXF_Real, P_diag_size);

   P_diag_i[n_fine] = jj_counter; 

   P_offd_size = jj_counter_offd;

   P_offd_i    = jxf_CTAlloc(JXF_Int, n_fine + 1);
   P_offd_j    = jxf_CTAlloc(JXF_Int, P_offd_size);
   P_offd_data = jxf_CTAlloc(JXF_Real, P_offd_size);

  /*-----------------------------------------------------------------------
   *  Intialize some stuff.
   *-----------------------------------------------------------------------*/

   jj_counter = start_indexing;
   jj_counter_offd = start_indexing;

   if (debug_flag == 4)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d     Interp: Internal work 1 =     %f\n", my_id, wall_time);
      fflush(NULL);
   }

  /*-----------------------------------------------------------------------
   *  Send and receive fine_to_coarse info.
   *-----------------------------------------------------------------------*/ 

   if (debug_flag == 4) wall_time = jxf_time_getWallclockSeconds();

   fine_to_coarse_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd); 

#define JXF_SMP_PRIVATE i,j,ns,ne,size,rest,coarse_shift
#include "../../include/jxf_smp_forloop.h"

   for (j = 0; j < num_threads; j ++)
   {
      coarse_shift = 0;
      
      if (j > 0) 
      {
         coarse_shift = coarse_counter[j-1];
      }
      
      size = n_fine / num_threads;
      rest = n_fine - size*num_threads;
      
      if (j < rest)
      {
         ns = j*size + j;
         ne = (j+1)*size + j + 1;
      }
      else
      {
         ns = j*size + rest;
         ne = (j+1)*size + rest;
      }
      
      for (i = ns; i < ne; i ++)
      {
         fine_to_coarse[i] += my_first_cpt+coarse_shift;
      }
   }
   
   index = 0;
   for (i = 0; i < num_sends; i ++)
   {
      start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
      for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
      {
         int_buf_data[index++] = fine_to_coarse[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
      }
   }
	
   comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, fine_to_coarse_offd);  
   jxf_ParCSRCommHandleDestroy(comm_handle);   

   if (debug_flag == 4)
   {
      wall_time = jxf_time_getWallclockSeconds() - wall_time;
      jxf_printf("Proc = %d     Interp: Comm 4 FineToCoarse = %f\n", my_id, wall_time);
      fflush(NULL);
   }

   if (debug_flag == 4) wall_time = jxf_time_getWallclockSeconds();

#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
   for (i = 0; i < n_fine; i ++) 
   {
      fine_to_coarse[i] -= my_first_cpt;
   }

   
  /*-----------------------------------------------------------------------
   *  Loop over fine grid points.
   *-----------------------------------------------------------------------*/
    
#define JXF_SMP_PRIVATE i,j,jl,i1,jj,ns,ne,size,rest,diagonal,P_marker,P_marker_offd,jj_counter,jj_counter_offd,jj_begin_row,jj_end_row,jj_begin_row_offd,jj_end_row_offd
#include "../../include/jxf_smp_forloop.h"

   for (jl = 0; jl < num_threads; jl ++)
   {
      size = n_fine / num_threads;
      rest = n_fine - size*num_threads;
      if (jl < rest)
      {
         ns = jl*size + jl;
         ne = (jl+1)*size + jl + 1;
      }
      else
      {
         ns = jl*size + rest;
         ne = (jl+1)*size + rest;
      }
      
      jj_counter = 0;
      if (jl > 0) 
      {
         jj_counter = jj_count[jl-1];
      }
      
      jj_counter_offd = 0;
      if (jl > 0) 
      {
         jj_counter_offd = jj_count_offd[jl-1];
      }

      P_marker = jxf_CTAlloc(JXF_Int, n_fine);
      if (num_cols_A_offd)
      {
         P_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_A_offd);
      }
      else
      {
         P_marker_offd = NULL;
      }

      for (i = 0; i < n_fine; i ++)
      {      
         P_marker[i] = -1;
      }
      for (i = 0; i < num_cols_A_offd; i ++)
      {      
         P_marker_offd[i] = -1;
      }
 
      for (i = ns; i < ne; i ++)
      {
             
        /*--------------------------------------------------------------------
         *  If i is a c-point, interpolation is the identity.
         *--------------------------------------------------------------------*/
      
         if (CF_marker[i] >= 0)
         {
            P_diag_i[i] = jj_counter;
            P_diag_j[jj_counter]    = fine_to_coarse[i];
            P_diag_data[jj_counter] = one;
            jj_counter ++;
         }
      
        /*--------------------------------------------------------------------
         *  If i is an F-point, build interpolation.
         *--------------------------------------------------------------------*/

         else
         {         
            /* Diagonal part of P */
            P_diag_i[i] = jj_counter;
            jj_begin_row = jj_counter;

            for (jj = S_diag_i[i]; jj < S_diag_i[i+1]; jj ++)
            {
               i1 = S_diag_j[jj];   

              /*--------------------------------------------------------------
               * If neighbor i1 is a C-point, set column number in P_diag_j
               * and initialize interpolation weight to zero.
               *--------------------------------------------------------------*/

               if (CF_marker[i1] >= 0)
               {
                  P_marker[i1] = jj_counter;
                  P_diag_j[jj_counter]    = fine_to_coarse[i1];
                  P_diag_data[jj_counter] = zero;
                  jj_counter ++;
               }

            } // end for jj
            
            jj_end_row = jj_counter;

            /* Off-Diagonal part of P */
            P_offd_i[i] = jj_counter_offd;
            jj_begin_row_offd = jj_counter_offd;


            if (num_procs > 1)
            {
               if (col_offd_S_to_A)
               {
                  for (jj = S_offd_i[i]; jj < S_offd_i[i+1]; jj ++)
                  {
                     i1 = col_offd_S_to_A[S_offd_j[jj]];   

                    /*-----------------------------------------------------------
                     * If neighbor i1 is a C-point, set column number in P_offd_j
                     * and initialize interpolation weight to zero.
                     *-----------------------------------------------------------*/

                     if (CF_marker_offd[i1] >= 0)
                     {
                        P_marker_offd[i1] = jj_counter_offd;
                        P_offd_j[jj_counter_offd]  = i1;
                        P_offd_data[jj_counter_offd] = zero;
                        jj_counter_offd ++;
                     }
                  } // end for jj
               }
               else
               {
                  for (jj = S_offd_i[i]; jj < S_offd_i[i+1]; jj ++)
                  {
                     i1 = S_offd_j[jj];   

                    /*-------------------------------------------------------------
                     * If neighbor i1 is a C-point, set column number in P_offd_j
                     * and initialize interpolation weight to zero.
                     *-----------------------------------------------------------*/

                     if (CF_marker_offd[i1] >= 0)
                     {
                        P_marker_offd[i1] = jj_counter_offd;
                        P_offd_j[jj_counter_offd]  = i1;
                        P_offd_data[jj_counter_offd] = zero;
                        jj_counter_offd ++;
                     }
                  } // end for jj
               } // end if (col_offd_S_to_A)
            } // end if (num_procs > 1)
      
            jj_end_row_offd = jj_counter_offd;
         
            diagonal = A_diag_data[A_diag_i[i]];

     
            /* Loop over ith row of A.  First, the diagonal part of A */
            sum_N_pos = 0;
            sum_N_neg = 0;
            sum_P_pos = 0;
            sum_P_neg = 0;

            for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
            {
            
               i1 = A_diag_j[jj];
               
               if (num_functions == 1 || dof_func[i1] == dof_func[i])
               { 
                  if (A_diag_data[jj] > 0)
                  {
	             sum_N_pos += A_diag_data[jj];
	          }
                  else
                  {
	             sum_N_neg += A_diag_data[jj];
	          }
               } 
               
              /*--------------------------------------------------------------
               * Case 1: neighbor i1 is a C-point and strongly influences i,
               * accumulate a_{i,i1} into the interpolation weight.
               *--------------------------------------------------------------*/

               if (P_marker[i1] >= jj_begin_row)
               {
                  P_diag_data[P_marker[i1]] += A_diag_data[jj];
                  if (A_diag_data[jj] > 0)
                  {
                     sum_P_pos += A_diag_data[jj];
                  }
                  else
                  {
                     sum_P_neg += A_diag_data[jj];
                  }
               }

            }  // end for jj  
       
           
           /*----------------------------------------------------------------
            * Still looping over ith row of A. Next, loop over the 
            * off-diagonal part of A 
            *---------------------------------------------------------------*/

            if (num_procs > 1)
            {
               for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
               {
                  i1 = A_offd_j[jj];
                  if (num_functions == 1 || dof_func_offd[i1] == dof_func[i])
                  { 
                     if (A_offd_data[jj] > 0)
                     {
                        sum_N_pos += A_offd_data[jj];
                     }
                     else
                     {
                        sum_N_neg += A_offd_data[jj];
                     }
                  } 

                 /*--------------------------------------------------------------
                  * Case 1: neighbor i1 is a C-point and strongly influences i,
                  * accumulate a_{i,i1} into the interpolation weight.
                  *--------------------------------------------------------------*/

                  if (P_marker_offd[i1] >= jj_begin_row_offd)
                  {
                     P_offd_data[P_marker_offd[i1]] += A_offd_data[jj];
                     if (A_offd_data[jj] > 0)
                     {
                        sum_P_pos += A_offd_data[jj];
                     }
                     else
                     {
                        sum_P_neg += A_offd_data[jj];
                     }
                  }
               } // end for jj
            }  // end if (num_procs > 1)
                    
            if (sum_P_neg) 
            {
               alfa = sum_N_neg / sum_P_neg / diagonal;
            }
            if (sum_P_pos) 
            {
               beta = sum_N_pos / sum_P_pos / diagonal;
            }

           /*-----------------------------------------------------------------
            * Set interpolation weight by dividing by the diagonal.
            *-----------------------------------------------------------------*/

            for (jj = jj_begin_row; jj < jj_end_row; jj ++)
            {
               if (P_diag_data[jj]> 0)
               {
                  P_diag_data[jj] *= -beta;
               }
               else
               {
                  P_diag_data[jj] *= -alfa;
               }
            }

            for (jj = jj_begin_row_offd; jj < jj_end_row_offd; jj ++)
            {
               if (P_offd_data[jj] > 0)
               {
                  P_offd_data[jj] *= -beta;
               }
               else
               {
                  P_offd_data[jj] *= -alfa;
               }
            }
           
         } // end if (CF_marker[i] >= 0)

         P_offd_i[i+1] = jj_counter_offd;
         
      } // end for i
     
      jxf_TFree(P_marker);
      jxf_TFree(P_marker_offd);
   
   } // end for jl
                                                                                
   par_P = jxf_ParCSRMatrixCreate(comm,
                                 jxf_ParCSRMatrixGlobalNumRows(par_A),
                                 total_global_cpts,
                                 jxf_ParCSRMatrixColStarts(par_A),
                                 num_cpts_global,
                                 0,
                                 P_diag_i[n_fine],
                                 P_offd_i[n_fine]);
                                                                                                                                                        
   P_diag = jxf_ParCSRMatrixDiag(par_P);
   jxf_CSRMatrixData(P_diag) = P_diag_data;
   jxf_CSRMatrixI(P_diag) = P_diag_i;
   jxf_CSRMatrixJ(P_diag) = P_diag_j;
   P_offd = jxf_ParCSRMatrixOffd(par_P);
   jxf_CSRMatrixData(P_offd) = P_offd_data;
   jxf_CSRMatrixI(P_offd) = P_offd_i;
   jxf_CSRMatrixJ(P_offd) = P_offd_j;
   jxf_ParCSRMatrixOwnsRowStarts(par_P) = 0;

   /* Compress P, removing coefficients smaller than trunc_factor*Max */

   if (trunc_factor != 0.0 || max_elmts > 0)
   {
      jxf_PAMGInterpTruncation(par_P, trunc_factor, max_elmts);
      P_diag_data = jxf_CSRMatrixData(P_diag);
      P_diag_i = jxf_CSRMatrixI(P_diag);
      P_diag_j = jxf_CSRMatrixJ(P_diag);
      P_offd_data = jxf_CSRMatrixData(P_offd);
      P_offd_i = jxf_CSRMatrixI(P_offd);
      P_offd_j = jxf_CSRMatrixJ(P_offd);
      P_diag_size = P_diag_i[n_fine];
      P_offd_size = P_offd_i[n_fine];
   }

   num_cols_P_offd = 0;
   if (P_offd_size)
   {
      P_marker = jxf_CTAlloc(JXF_Int, num_cols_A_offd);

#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols_A_offd; i ++)
      {
         P_marker[i] = 0;
      }
                                                                                
      num_cols_P_offd = 0;
      for (i = 0; i < P_offd_size; i ++)
      {
         index = P_offd_j[i];
         if (!P_marker[index])
         {
            num_cols_P_offd ++;
            P_marker[index] = 1;
         }
      }
                                                                                
      col_map_offd_P = jxf_CTAlloc(JXF_Int, num_cols_P_offd);
                                                                                
      index = 0;
      for (i = 0; i < num_cols_P_offd; i ++)
      {
         while (P_marker[index] == 0) 
         {
            index ++;
         }
         col_map_offd_P[i] = index ++;
      }

#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
      for (i = 0; i < P_offd_size; i ++)
      {
         P_offd_j[i] = jxf_BinarySearch(col_map_offd_P, P_offd_j[i], num_cols_P_offd);
      }
      jxf_TFree(P_marker); 
   }

   for (i = 0; i < n_fine; i ++)
   {
      if (CF_marker[i] == -3) 
      {
         CF_marker[i] = -1;
      }
   }
   
   if (num_cols_P_offd)
   { 
      jxf_ParCSRMatrixColMapOffd(par_P) = col_map_offd_P;
      jxf_CSRMatrixNumCols(P_offd) = num_cols_P_offd;
   } 

   jxf_GetCommPkgRTFromCommPkgA(par_P, par_A, fine_to_coarse_offd); 

   *P_ptr = par_P;

   jxf_TFree(CF_marker_offd);
   jxf_TFree(dof_func_offd);
   jxf_TFree(int_buf_data);
   jxf_TFree(fine_to_coarse);
   jxf_TFree(fine_to_coarse_offd);
   jxf_TFree(coarse_counter);
   jxf_TFree(jj_count);
   jxf_TFree(jj_count_offd);

   return(0);  

}
