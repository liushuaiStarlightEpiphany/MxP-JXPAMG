//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_com.c -- communication operations for parallel matrices and vectors.
 *  Date: 2011/09/07
 */ 

#include "jxf_mv.h"

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

/*!
 * \fn JXF_Int jxf_MatvecCommPkgCreate
 * \brief Generates the comm_pkg for A. 
 *        If no row and / or column partitioning is given, 
 *        the routine determines them with MPE_Decomp1d.
 * \date 2011/09/07
 */
JXF_Int
jxf_MatvecCommPkgCreate (jxf_ParCSRMatrix *par_A)
{
   JXF_Int			 num_sends;
   JXF_Int			*send_procs;
   JXF_Int			*send_map_starts;
   JXF_Int			*send_map_elmts;
 
   JXF_Int			 num_recvs;
   JXF_Int			*recv_procs;
   JXF_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jxf_ParCSRMatrixComm(par_A);
   
   jxf_ParCSRCommPkg	*comm_pkg;

   JXF_Int   first_col_diag = jxf_ParCSRMatrixFirstColDiag(par_A);
   JXF_Int  *col_map_offd   = jxf_ParCSRMatrixColMapOffd(par_A);
   JXF_Int	 num_cols_offd  = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(par_A));


#if JXF_NO_GLOBAL_PARTITION

   JXF_Int  row_start = 0, row_end = 0;
   JXF_Int  col_start = 0, col_end = 0;
   JXF_Int  global_num_cols;
   
   jxf_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jxf_ParCSRMatrixGetLocalRange( par_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jxf_ParCSRMatrixGlobalNumCols(par_A); 

   /* Create the assumed partition */
   if (jxf_ParCSRMatrixAssumedPartition(par_A) == NULL)
   {
      jxf_ParCSRMatrixCreateAssumedPartition(par_A);
   }
   apart = jxf_ParCSRMatrixAssumedPartition(par_A);
   
  /*-----------------------------------------------------------
   * get commpkg info information 
   *----------------------------------------------------------*/

   jxf_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
                             col_start, col_end, 
                             num_cols_offd, global_num_cols,
                             &num_recvs, &recv_procs, &recv_vec_starts,
                             &num_sends, &send_procs, &send_map_starts, 
                             &send_map_elmts, apart );
   
#else
   
   JXF_Int *col_starts    = jxf_ParCSRMatrixColStarts(par_A);
   JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_ParCSRMatrixDiag(par_A));
  

   jxf_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
                                num_cols_diag, num_cols_offd,
                                first_col_diag, col_map_offd, 1,
                                &num_recvs, &recv_procs, &recv_vec_starts,
                                &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
#endif


  /*-----------------------------------------------------------
   *  setup commpkg
   *----------------------------------------------------------*/

   comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);

   jxf_ParCSRCommPkgComm(comm_pkg) = comm;

   jxf_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
   jxf_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
   jxf_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
   jxf_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
   jxf_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
   jxf_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
   jxf_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;

#if JXF_REODER_SEND_RECV
jxf_hpIfReorderMatvecCommPkg(comm_pkg)   = 0;
  /*-----------------------------------------------------------
   *  reoderMatvecCommPkg
   *----------------------------------------------------------*/
         // JXF_Int myid,i;
         // jxf_MPI_Comm_rank(comm, &myid);
         
         // if(myid == 1)
         // {
         //    jxf_printf("raw send_procs is :\n");
         //    for(i = 0; i < num_sends;i++)
         //    {
         //       jxf_printf("%d\t",send_procs[i]);
         //    }
         //    jxf_printf("\n");
         //    jxf_printf("raw SendMapStarts is :\n");
         //    for(i = 0; i < num_sends + 1;i++)
         //    {
         //       jxf_printf("%d\t",send_map_starts[i]);
         //    }
         //    jxf_printf("\n\n");
         // }
   jxf_hpCSRCommPkgReorder(comm_pkg);
         // if(myid == 1)
         // {
         //    jxf_printf("reorder send_procs is :\n");
         //    for(i = 0; i < num_sends;i++)
         //    {
         //       jxf_printf("%d\t",jxf_hpSendReorderProc(comm_pkg, i));
         //    }
         //    jxf_printf("\n");
         //    jxf_printf("send_reorder_map is :\n");
         //    for(i = 0; i < num_sends;i++)
         //    {
         //       jxf_printf("%d\t",jxf_hpSendReorderMapIndex(comm_pkg, i));
         //    }
         //    jxf_printf("\n");
         //    jxf_printf("send_map_starts is :\n");
         //    for(i = 0; i < num_sends + 1;i++)
         //    {
         //       jxf_printf("%d\t",jxf_ParCSRCommPkgSendMapStart(comm_pkg, i));
         //    }
         //    jxf_printf("\n\n");
         // }
#endif

   jxf_ParCSRMatrixCommPkg(par_A) = comm_pkg;

   return jxf_error_flag;
}


/*!
 * \fn jxf_ParCSRCommHandle *jxf_ParCSRCommHandleCreate
 * \brief Sets up a communication handle, posts receives and initiates sends. 
 *        It always requires num_sends, num_recvs, recv_procs and send_procs 
 *        to be set in comm_pkg.
 * \note There are different options for job:
 * job = 1 : is used to initialize communication exchange for the parts
 *           of vector needed to perform a Matvec,  it requires send_data 
 *           and recv_data to be doubles, recv_vec_starts and 
 *           send_map_starts need to be set in comm_pkg.
 * job = 2 : is used to initialize communication exchange for the parts
 *           of vector needed to perform a MatvecT,  it requires send_data 
 *           and recv_data to be doubles, recv_vec_starts and 
 *           send_map_starts need to be set in comm_pkg.
 * job = 11: similar to job = 1, but exchanges data of type JXF_Int (not JXF_Real),
 *           requires send_data and recv_data to be ints
 *           recv_vec_starts and send_map_starts need to be set in comm_pkg.
 * job = 12: similar to job = 1, but exchanges data of type JXF_Int (not JXF_Real),
 *           requires send_data and recv_data to be ints
 *           recv_vec_starts and send_map_starts need to be set in comm_pkg.
 * [default:  ignores send_data and recv_data, requires send_mpi_types
 *           and recv_mpi_types to be set in comm_pkg.
 *           datatypes need to point to absolute
 *           addresses, e.g. generated using jxf_MPI_Address.]  
 * \date 2011/09/07
 */
jxf_ParCSRCommHandle *
jxf_ParCSRCommHandleCreate( JXF_Int 	             job,
			   jxf_ParCSRCommPkg *comm_pkg,
                           void             *send_data, 
                           void             *recv_data )
{
   JXF_Int                  num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   JXF_Int                  num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
   MPI_Comm             comm      = jxf_ParCSRCommPkgComm(comm_pkg);

   jxf_ParCSRCommHandle *comm_handle;
   JXF_Int                  num_requests;
   MPI_Request         *requests;

   JXF_Int                  i, j;
   JXF_Int			my_id, num_procs;
   JXF_Int			ip, vec_start, vec_len;

   num_requests = num_sends + num_recvs;
   requests = jxf_CTAlloc(MPI_Request, num_requests);
 
   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   j = 0;
   switch (job)
   {
      case  1:
      {
         JXF_Real *d_send_data = (JXF_Real *) send_data;
         JXF_Real *d_recv_data = (JXF_Real *) recv_data;
         for (i = 0; i < num_recvs; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               //    vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpRecvReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpRecvVecStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpRecvVecStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
               vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
            #endif
            jxf_MPI_Irecv(&d_recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
         }
         for (i = 0; i < num_sends; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
               //    vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
               //    vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpSendReorderProc(comm_pkg, i);
                  vec_start = jxf_hpSendMapStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpSendMapStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
	            vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	            vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
            #endif
   	    jxf_MPI_Isend(&d_send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
         }
      }
      break;
      
      
      case  2:
      {
         JXF_Real *d_send_data = (JXF_Real *) send_data;
         JXF_Real *d_recv_data = (JXF_Real *) recv_data;
         for (i = 0; i < num_sends; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
               //    vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpSendReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpSendMapStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpSendMapStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
	            vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	            vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
            #endif
            jxf_MPI_Irecv(&d_recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
         }
         for (i = 0; i < num_recvs; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i);
               //    vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               //    vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpRecvReorderProc(comm_pkg, i);
                  vec_start = jxf_hpRecvVecStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpRecvVecStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i);
               vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
            #endif
            jxf_MPI_Isend(&d_send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
         }
      }
      break;
      
      
      case  11:
      {
         JXF_Int *i_send_data = (JXF_Int *) send_data;
         JXF_Int *i_recv_data = (JXF_Int *) recv_data;
         for (i = 0; i < num_recvs; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               //    vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpRecvReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpRecvVecStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpRecvVecStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
                  ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
                  vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
                  vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
            #endif
            jxf_MPI_Irecv(&i_recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
         }
         for (i = 0; i < num_sends; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
               //    vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpSendReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpSendMapStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpSendMapStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
	            vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	            vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
            #endif
            jxf_MPI_Isend(&i_send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
         }
      }
      break;
      
      
      case  12:
      {
         JXF_Int *i_send_data = (JXF_Int *) send_data;
         JXF_Int *i_recv_data = (JXF_Int *) recv_data;
         for (i = 0; i < num_sends; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
               //    vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpSendReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpSendMapStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpSendMapStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
	            vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
	            vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
            #endif
            jxf_MPI_Irecv(&i_recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
         }
         for (i = 0; i < num_recvs; i ++)
         {
            #if JXF_REODER_SEND_RECV 
               // if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
               // {
               //    ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
               //    vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               //    vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
               // }
               // else
               // {
                  ip        = jxf_hpRecvReorderProc(comm_pkg, i); 
                  vec_start = jxf_hpRecvVecStartReorder(comm_pkg,i,0);
                  vec_len   = jxf_hpRecvVecStartReorder(comm_pkg,i,1) - vec_start;
               // }
            #else
               ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
               vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
               vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
            #endif
            jxf_MPI_Isend(&i_send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
         }
         break;
      }
   }
 
 
  /*--------------------------------------------------------------------
   * set up comm_handle and return
   *--------------------------------------------------------------------*/

   comm_handle = jxf_CTAlloc(jxf_ParCSRCommHandle, 1);

   jxf_ParCSRCommHandleCommPkg(comm_handle)     = comm_pkg;
   jxf_ParCSRCommHandleSendData(comm_handle)    = send_data;
   jxf_ParCSRCommHandleRecvData(comm_handle)    = recv_data;
   jxf_ParCSRCommHandleNumRequests(comm_handle) = num_requests;
   jxf_ParCSRCommHandleRequests(comm_handle)    = requests;

   return ( comm_handle );
}


/*!
 * \fn JXF_Int jxf_ParCSRCommHandleDestroy 
 * \date 2011/09/05
 */
JXF_Int
jxf_ParCSRCommHandleDestroy( jxf_ParCSRCommHandle *comm_handle )
{
   MPI_Status   *status   = NULL;
   MPI_Request  *requests = jxf_ParCSRCommHandleRequests(comm_handle);
   JXF_Int num_requests = jxf_ParCSRCommHandleNumRequests(comm_handle);

   if (comm_handle == NULL) 
   {
      return jxf_error_flag;
   }
   
   if (num_requests)
   {
      status = jxf_CTAlloc(MPI_Status, num_requests);
      jxf_MPI_Waitall(num_requests, requests, status);
      jxf_TFree(status);
   }

   jxf_TFree(requests);
   jxf_TFree(comm_handle);

   return jxf_error_flag;
}


/*!
 * \fn void jxf_MatvecCommPkgCreate_core
 * \brief This routine does all the communications and computations for
 *            jxf_MatCommPkgCreate(jxf_ParCSRMatrix *par_A)
 *        and  
 *            jxf_BoolMatCommPkgCreate(jxf_ParCSRBooleanMatrix *par_A)
 *        To support both data types, it has hardly any data structures other than JXF_Int *.
 * \date 2011/09/05
 */
void
jxf_MatvecCommPkgCreate_core(  /* input args: */
                              MPI_Comm comm,      JXF_Int *col_map_offd, 
                              JXF_Int first_col_diag, JXF_Int *col_starts,
                              JXF_Int num_cols_diag,  JXF_Int  num_cols_offd,
                              JXF_Int firstColDiag,   JXF_Int *colMapOffd,
                              /* = 1 for a matrix with floating-point data, = 0 for Boolean matrix */
                              JXF_Int data,  
                              /* pointers to output args: */
                              JXF_Int  *p_num_recvs,  JXF_Int **p_recv_procs,  JXF_Int **p_recv_vec_starts,
                              JXF_Int  *p_num_sends,  JXF_Int **p_send_procs,  JXF_Int **p_send_map_starts,
                              JXF_Int **p_send_map_elmts )
{
   JXF_Int	i, j;
   JXF_Int	num_procs, my_id, proc_num, num_elmts;
   JXF_Int	local_info, offd_col;
   JXF_Int	*proc_mark, *proc_add, *tmp, *recv_buf, *displs, *info;
   
   /* outputs: */
   JXF_Int  num_recvs, * recv_procs, * recv_vec_starts;
   JXF_Int  num_sends, * send_procs, * send_map_starts, * send_map_elmts;
   JXF_Int  ip, vec_start, vec_len, num_requests;

   MPI_Request *requests = NULL;
   MPI_Status  *status   = NULL; 

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);

   proc_mark = jxf_CTAlloc(JXF_Int, num_procs);
   proc_add  = jxf_CTAlloc(JXF_Int, num_procs);
   info      = jxf_CTAlloc(JXF_Int, num_procs);

  /* ----------------------------------------------------------------------
   * determine which processors to receive from (set proc_mark) and num_recvs,
   * at the end of the loop proc_mark[i] contains the number of elements to be
   * received from Proc. i
   * ---------------------------------------------------------------------*/

   for (i = 0; i < num_procs; i ++)
   {
      proc_add[i] = 0;
   }

   proc_num = 0;
   
   if (num_cols_offd) offd_col = col_map_offd[0];
   
   num_recvs = 0;
   j = 0;
   for (i = 0; i < num_cols_offd; i ++)
   {
      if (num_cols_diag) 
      {
         proc_num = jxf_min(num_procs-1, offd_col / num_cols_diag);
      }
	
      while (col_starts[proc_num] > offd_col)
      {
         proc_num = proc_num - 1;
      }
      
      while (col_starts[proc_num+1] - 1 < offd_col)
      {
         proc_num = proc_num + 1;
      }
      
      proc_mark[num_recvs] = proc_num;
      j = i;
      while (col_starts[proc_num+1] > offd_col)
      {
         proc_add[num_recvs] ++;
         if (j < num_cols_offd - 1) 
         {
            j ++;
            offd_col = col_map_offd[j];
         }
         else
         {
            j ++;
            offd_col = col_starts[num_procs];
         }
      }
      num_recvs ++;
	
      if (j < num_cols_offd) 
      {
         i = j - 1;
      }
      else 
      {
         i = j;
      }
   }

   local_info = 2*num_recvs;
			
   jxf_MPI_Allgather(&local_info, 1, JXF_MPI_INT, info, 1, JXF_MPI_INT, comm); 

  /* --------------------------------------------------------------------------
   * generate information to be sent: tmp contains for each recv_proc:
   * id of recv_procs, number of elements to be received for this processor,
   * indices of elements (in this order)
   *---------------------------------------------------------------------------*/

   displs = jxf_CTAlloc(JXF_Int, num_procs + 1);
   displs[0] = 0;
   for (i = 1; i < num_procs + 1; i ++)
   {
      displs[i] = displs[i-1] + info[i-1]; 
   }
   recv_buf = jxf_CTAlloc(JXF_Int, displs[num_procs]); 

   recv_procs = NULL;
   tmp = NULL;
   if (num_recvs)
   {
      recv_procs = jxf_CTAlloc(JXF_Int, num_recvs);
      tmp        = jxf_CTAlloc(JXF_Int, local_info);
   }
   recv_vec_starts = jxf_CTAlloc(JXF_Int, num_recvs + 1);


   j = 0;
   if (num_recvs) recv_vec_starts[0] = 0;
   for (i = 0; i < num_recvs; i ++)
   {
      num_elmts = proc_add[i];
      recv_procs[i] = proc_mark[i];
      recv_vec_starts[i+1] = recv_vec_starts[i] + num_elmts;
      tmp[j++] = proc_mark[i];
      tmp[j++] = num_elmts;
   }

   jxf_MPI_Allgatherv(tmp, local_info, JXF_MPI_INT, recv_buf, info, displs, JXF_MPI_INT, comm);
	

  /*----------------------------------------------------------------------
   * determine num_sends and number of elements to be sent
   *---------------------------------------------------------------------*/

   num_sends = 0;
   num_elmts = 0;
   proc_add[0] = 0;
   for (i = 0; i < num_procs; i ++)
   {
      j = displs[i];
      while (j < displs[i+1])
      {
	 if (recv_buf[j++] == my_id)
	 {
	    proc_mark[num_sends] = i;
	    num_sends ++;
	    proc_add[num_sends] = proc_add[num_sends-1] + recv_buf[j];
	    break;
	 }
	 j ++;
      }	
   }
	
		
  /* ----------------------------------------------------------------------
   * determine send_procs and actual elements to be send (in send_map_elmts)
   * and send_map_starts whose i-th entry points to the beginning of the 
   * elements to be send to proc. i
   * ---------------------------------------------------------------------*/

   send_procs = NULL;
   send_map_elmts = NULL;

   if (num_sends)
   {
      send_procs = jxf_CTAlloc(JXF_Int, num_sends);
      send_map_elmts = jxf_CTAlloc(JXF_Int, proc_add[num_sends]);
   }
   send_map_starts = jxf_CTAlloc(JXF_Int, num_sends+1);
   num_requests = num_recvs+num_sends;
   if (num_requests)
   {
      requests = jxf_CTAlloc(MPI_Request, num_requests);
      status   = jxf_CTAlloc(MPI_Status, num_requests);
   }

   if (num_sends) send_map_starts[0] = 0;
   for (i = 0; i < num_sends; i ++)
   {
      send_map_starts[i+1] = proc_add[i+1];
      send_procs[i] = proc_mark[i];
   }

   j = 0;
   for (i = 0; i < num_sends; i ++)
   {
      vec_start = send_map_starts[i];
      vec_len = send_map_starts[i+1] - vec_start;
      ip = send_procs[i];
      jxf_MPI_Irecv(&send_map_elmts[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
   }
   for (i = 0; i < num_recvs; i ++)
   {
      vec_start = recv_vec_starts[i];
      vec_len = recv_vec_starts[i+1] - vec_start;
      ip = recv_procs[i];
      jxf_MPI_Isend(&col_map_offd[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
   }

   if (num_requests)
   {
      jxf_MPI_Waitall(num_requests, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }

   if (num_sends)
   {
      for (i = 0; i < send_map_starts[num_sends]; i ++)
      {
         send_map_elmts[i] -= first_col_diag;
      }
   }

   // if(my_id == 7)
   // {
   //    jxf_printf("num_sends is: \n");
   //    for(i = 0; i < num_sends; i++)
   //    {
   //       jxf_printf("%d\t",send_procs[i]);
   //    }
   //    jxf_printf("\n");
   //    jxf_printf("send_map_starts:\n");
   //    for(i = 0; i < num_sends+1; i++)
   //    {
   //       jxf_printf("%d\t",send_map_starts[i]);
   //    }
   //    jxf_printf("\n");
   //    // jxf_printf("send_map_elmts:\n");
   //    // for(i = 0; i < proc_add[num_sends]; i++)
   //    // {
   //    //    jxf_printf("%d\t",send_map_elmts[i]);
   //    // }
   //    // jxf_printf("\n");
   // }

   jxf_TFree(proc_add);
   jxf_TFree(proc_mark); 
   jxf_TFree(tmp);
   jxf_TFree(recv_buf);
   jxf_TFree(displs);
   jxf_TFree(info);
 
   /* finish up with the hand-coded call-by-reference... */
   *p_num_recvs = num_recvs;
   *p_recv_procs = recv_procs;
   *p_recv_vec_starts = recv_vec_starts;
   *p_num_sends = num_sends;
   *p_send_procs = send_procs;
   *p_send_map_starts = send_map_starts;
   *p_send_map_elmts = send_map_elmts;
}


/*!
 * \fn JXF_Int jxf_MatvecCommPkgDestroy
 * \date 2011/09/05
 */
JXF_Int
jxf_MatvecCommPkgDestroy( jxf_ParCSRCommPkg *comm_pkg )
{
   if (jxf_ParCSRCommPkgNumSends(comm_pkg))
   {
      jxf_TFree(jxf_ParCSRCommPkgSendProcs(comm_pkg));
      jxf_TFree(jxf_ParCSRCommPkgSendMapElmts(comm_pkg));
   }
   jxf_TFree(jxf_ParCSRCommPkgSendMapStarts(comm_pkg));

   if (jxf_ParCSRCommPkgNumRecvs(comm_pkg))
   {
      jxf_TFree(jxf_ParCSRCommPkgRecvProcs(comm_pkg));
   }
   jxf_TFree(jxf_ParCSRCommPkgRecvVecStarts(comm_pkg));

#if JXF_REODER_SEND_RECV
   if(jxf_hpSendReorderMap(comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderMap(comm_pkg));
   }
   if(jxf_hpRecvReorderMap(comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderMap(comm_pkg));
   }
   if(jxf_hpSendReorderProcs(comm_pkg))
   {
      jxf_TFree(jxf_hpSendReorderProcs(comm_pkg));
   }
   if(jxf_hpRecvReorderProcs(comm_pkg))
   {
      jxf_TFree(jxf_hpRecvReorderProcs(comm_pkg));
   }
#endif

   jxf_TFree(comm_pkg);

   return jxf_error_flag;
}
