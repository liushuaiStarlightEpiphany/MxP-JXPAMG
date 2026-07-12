//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hp_csr_com.c -- communication operations for hierarchy parallel matrices and vectors.
 *  Date: 2011/09/07
 */ 

#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

// JX_Int
// jx_hpMatvecCommPkgCreate (jx_hpCSRMatrix *hp_A, jx_HardwareInfo *hardware)
// {
//    JX_Int			 num_sends;
//    JX_Int			*send_procs;
//    JX_Int			*send_map_starts;
//    JX_Int			*send_map_elmts;
 
//    JX_Int			 num_recvs;
//    JX_Int			*recv_procs;
//    JX_Int			*recv_vec_starts;
   
//    MPI_Comm              comm = jx_hpCSRMatrixComm(hp_A);
//    jx_ParCSRCommPkg	*comm_pkg;

//    JX_Int   first_col_diag = jx_hpCSRMatrixFirstColDiag(hp_A);
//    JX_Int  *col_map_offd   = jx_hpCSRMatrixColMapOffd(hp_A);
//    JX_Int	 num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffd(hp_A));


// #if JX_NO_GLOBAL_PARTITION

//    JX_Int  row_start = 0, row_end = 0;
//    JX_Int  col_start = 0, col_end = 0;
//    JX_Int  global_num_cols;
   
//    jx_IJAssumedPart  *apart;
   
//    /*  get parcsr_A information */
//    jx_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
//    global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A); 

//    /* Create the assumed partition */
//    if (jx_hpCSRMatrixAssumedPartition(hp_A) == NULL)
//    {
//       jx_hpCSRMatrixCreateAssumedPartition(hp_A);
//    }
//    apart = jx_hpCSRMatrixAssumedPartition(hp_A);
   
//   /*-----------------------------------------------------------
//    * get commpkg info information 
//    *----------------------------------------------------------*/

//    jx_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
//                              col_start, col_end, 
//                              num_cols_offd, global_num_cols,
//                              &num_recvs, &recv_procs, &recv_vec_starts,
//                              &num_sends, &send_procs, &send_map_starts, 
//                              &send_map_elmts, apart );
   
// #else
   
//    JX_Int *col_starts    = jx_hpCSRMatrixColStarts(hp_A);
//    JX_Int	num_cols_diag = jx_CSRMatrixNumCols(jx_hpCSRMatrixDiag(hp_A));
  

//    jx_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
//                                 num_cols_diag, num_cols_offd,
//                                 first_col_diag, col_map_offd, 1,
//                                 &num_recvs, &recv_procs, &recv_vec_starts,
//                                 &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
// #endif


//   /*-----------------------------------------------------------
//    *  setup commpkg
//    *----------------------------------------------------------*/

//    comm_pkg = jx_CTAlloc(jx_ParCSRCommPkg, 1);

//    jx_ParCSRCommPkgComm(comm_pkg) = comm;

//   /*-----------------------------------------------------------
//    *  reoderMatvecCommPkg
//    *----------------------------------------------------------*/

//    if(!jx_hpIfReorderMatvecCommPkg(comm_pkg))
//    {
//       jx_hpCSRCommCommPkgSendReorder(num_sends, hardware, &send_procs);
//       jx_hpCSRCommCommPkgRecvReorder(num_recvs, hardware, &recv_procs);
//       jx_hpIfReorderMatvecCommPkg(comm_pkg) = 1;
//    }

//    jx_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
//    jx_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
//    jx_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
//    jx_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
//    jx_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
//    jx_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
//    jx_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;

//    jx_hpCSRMatrixCommPkg(hp_A) = comm_pkg;

//    return jx_error_flag;
// }

// JX_Int
// jx_hpMatvecCommPkgCreateLevel (jx_hpCSRMatrix *hp_A, JX_Int Level, jx_HardwareInfo *hardware)
// {
//    JX_Int			 num_sends;
//    JX_Int			*send_procs;
//    JX_Int			*send_map_starts;
//    JX_Int			*send_map_elmts;
 
//    JX_Int			 num_recvs;
//    JX_Int			*recv_procs;
//    JX_Int			*recv_vec_starts;

//    JX_Int	      num_cols_offd = 0;
//    JX_Int         *col_map_offd = NULL;
   
//    MPI_Comm              comm = jx_hpCSRMatrixComm(hp_A);
//    jx_ParCSRCommPkg	*comm_pkg;

//    if(Level == 1)
//    {
//       col_map_offd   = jx_hpCSRMatrixColMapOffdNode(hp_A);
//       //num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdNode(hp_A));
//    }
//    else if(Level == 2)
//    {
//       col_map_offd   = jx_hpCSRMatrixColMapOffdCpu(hp_A);
//       //num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdCpu(hp_A));
//    }
//    else
//    {
//       jx_printf("level must be 1 or 2\n");
//       return jx_error_flag;
//    }
//    num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffd(hp_A));

//    // JX_Int myid ;
//    // jx_MPI_Comm_rank(comm, &myid);
//    // JX_Int i;
//    // if(myid == 0)
//    // {
//    //    jx_printf("\n myid is %d, col_map_offd is ", myid);
//    //    for(i = 0; i < num_cols_offd; i++)
//    //    {
//    //       jx_printf("%d\t",col_map_offd[i]);
//    //    }
//    //    jx_printf("\n");
//    // }
   
//    JX_Int   first_col_diag = jx_hpCSRMatrixFirstColDiag(hp_A);

// #if JX_NO_GLOBAL_PARTITION

//    JX_Int  row_start = 0, row_end = 0;
//    JX_Int  col_start = 0, col_end = 0;
//    JX_Int  global_num_cols;
   
//    jx_IJAssumedPart  *apart;
   
//    /*  get parcsr_A information */
//    jx_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
//    global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A); 

//    /* Create the assumed partition */
//    if (jx_hpCSRMatrixAssumedPartition(hp_A) == NULL)
//    {
//       jx_hpCSRMatrixCreateAssumedPartition(hp_A);
//    }
//    apart = jx_hpCSRMatrixAssumedPartition(hp_A);
   
//   /*-----------------------------------------------------------
//    * get commpkg info information 
//    *----------------------------------------------------------*/

//    jx_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
//                              col_start, col_end, 
//                              num_cols_offd, global_num_cols,
//                              &num_recvs, &recv_procs, &recv_vec_starts,
//                              &num_sends, &send_procs, &send_map_starts, 
//                              &send_map_elmts, apart );
   
// #else
   
//    // JX_Int *col_starts    = NULL;
//    // if(Level == 1)
//    // {
//    //    col_starts = jx_hpCSRMatrixColNodeStart(hp_A);
//    // }
//    // else if(Level == 2)
//    // {
//    //    col_starts = jx_hpCSRMatrixColCpuStart(hp_A);
//    // }
//    JX_Int *col_starts    = jx_hpCSRMatrixColStarts(hp_A);
//    JX_Int	num_cols_diag = jx_CSRMatrixNumCols(jx_hpCSRMatrixDiag(hp_A));
  

// // jx_printf("before run success\n");

   
//    //if(myid == 0) jx_printf("num_cols_diag is %d, num_cols_offd is %d\n",num_cols_diag, num_cols_offd);

//    jx_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
//                                 num_cols_diag, num_cols_offd,
//                                 first_col_diag, col_map_offd, 
//                                 1,
//                                 &num_recvs, &recv_procs, &recv_vec_starts,
//                                 &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
// #endif

// jx_printf("jx_MatvecCommPkgCreate_core run success\n");
//   /*-----------------------------------------------------------
//    *  setup commpkg
//    *----------------------------------------------------------*/

//    comm_pkg = jx_CTAlloc(jx_ParCSRCommPkg, 1);

//    jx_ParCSRCommPkgComm(comm_pkg) = comm;
   
//    jx_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
//    jx_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
//    jx_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
//    jx_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
//    jx_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
//    jx_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
//    jx_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;

//    //jx_hpCSRMatrixCommPkg(hp_A) = comm_pkg;
//    if(Level == 1)
//    {
//       jx_hpCSRMatrixNodeCommPkg(hp_A) = jx_CTAlloc(jx_ParCSRCommPkg, 1);
//       jx_hpCSRMatrixNodeCommPkg(hp_A) = comm_pkg;
//    }
//    else if(Level == 2) 
//    {
//       jx_hpCSRMatrixCpuCommPkg(hp_A) = jx_CTAlloc(jx_ParCSRCommPkg, 1);
//       jx_hpCSRMatrixCpuCommPkg(hp_A) = comm_pkg;
//    }

//    return jx_error_flag;
// }

// jx_ParCSRCommHandle *
// jx_hpCSRCommHandleCreate(  JX_Int 	             job,
// 			                  jx_ParCSRCommPkg *comm_pkg,
//                            void             *send_data, 
//                            void             *recv_data, 
//                            jx_HardwareInfo *hardware)
// {
//    JX_Int                  num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
//    JX_Int                  num_recvs = jx_ParCSRCommPkgNumRecvs(comm_pkg);
//    MPI_Comm             comm      = jx_ParCSRCommPkgComm(comm_pkg);

//    jx_ParCSRCommHandle *comm_handle;
//    JX_Int                  num_requests;
//    MPI_Request         *requests;

//    JX_Int                  i, j;
//    JX_Int			my_id, num_procs;
//    JX_Int			ip, vec_start, vec_len;

//    num_requests = num_sends + num_recvs;
//    requests = jx_CTAlloc(MPI_Request, num_requests);
 
//    jx_MPI_Comm_size(comm, &num_procs);
//    jx_MPI_Comm_rank(comm, &my_id);

//    j = 0;
//    switch (job)
//    {
//       case  1:
//       {
//          JX_Real *d_send_data = (JX_Real *) send_data;
//          JX_Real *d_recv_data = (JX_Real *) recv_data;
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jx_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jx_MPI_Irecv(&d_recv_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jx_ParCSRCommPkgSendProc(comm_pkg, i);
// 	         vec_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
// 	         vec_len   = jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//    	      jx_MPI_Isend(&d_send_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//       }
//       break;
      
      
//       case  2:
//       {
//          JX_Real *d_send_data = (JX_Real *) send_data;
//          JX_Real *d_recv_data = (JX_Real *) recv_data;
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jx_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jx_MPI_Irecv(&d_recv_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jx_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jx_MPI_Isend(&d_send_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//       }
//       break;
      
//       // case  3:
//       // {
//       //    /*
//       //    1.进程先要重新分为结点内和结点外
//       //    2.结点外先发
//       //    3.结点内，CPU外再发
//       //    4.结点内，CPU内最后发
//       //    */
         
         
//       //    //数组如何按照硬件排序，时间复杂度尽可能低
//       //    //遍历两遍：
//       //    //第一遍得到发送的有多少CPU内，结点内
//       //    //第二遍再重排


//       //    for (i = 0; i < num_recvs; i ++)
//       //    {
//       //       ip        = ip_recv_reorder[i]; 
//       //       vec_start = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//       //       vec_len   = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//       //       jx_MPI_Irecv(&d_recv_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//       //    }

//       //    for (i = 0; i < num_sends; i ++)
//       //    {
//       //       ip        = ip_send_reorder[i];
//       //       vec_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
//       //       vec_len   = jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//       //       jx_MPI_Isend(&d_send_data[vec_start], vec_len, JX_MPI_REAL, ip, 0, comm, &requests[j++]);
//       //    }
//       // }
//       // break;
      
//       case  11:
//       {
//          JX_Int *i_send_data = (JX_Int *) send_data;
//          JX_Int *i_recv_data = (JX_Int *) recv_data;
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jx_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jx_MPI_Irecv(&i_recv_data[vec_start], vec_len, JX_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jx_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jx_MPI_Isend(&i_send_data[vec_start], vec_len, JX_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//       }
//       break;
      
      
//       case  12:
//       {
//          JX_Int *i_send_data = (JX_Int *) send_data;
//          JX_Int *i_recv_data = (JX_Int *) recv_data;
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jx_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jx_MPI_Irecv(&i_recv_data[vec_start], vec_len, JX_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jx_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jx_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jx_MPI_Isend(&i_send_data[vec_start], vec_len, JX_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          break;
//       }
//    }
 
 
//   /*--------------------------------------------------------------------
//    * set up comm_handle and return
//    *--------------------------------------------------------------------*/

//    comm_handle = jx_CTAlloc(jx_ParCSRCommHandle, 1);

//    jx_ParCSRCommHandleCommPkg(comm_handle)     = comm_pkg;
//    jx_ParCSRCommHandleSendData(comm_handle)    = send_data;
//    jx_ParCSRCommHandleRecvData(comm_handle)    = recv_data;
//    jx_ParCSRCommHandleNumRequests(comm_handle) = num_requests;
//    jx_ParCSRCommHandleRequests(comm_handle)    = requests;

//    return ( comm_handle );
// }

//1.进程先要重新分为结点内和结点外
//2.结点外先发
//3.结点内，CPU外再发
//4.结点内，CPU内最后发

//数组如何按照硬件排序，时间复杂度尽可能低
//遍历两遍：
//第一遍得到发送的有多少CPU内，结点内
//第二遍再重排
//只有一个结点或CPU的时候不需要重新排序
#if JX_REODER_SEND_RECV

JX_Int jx_hpCSRCommPkgRecvVecStartReorder(jx_ParCSRCommPkg *comm_pkg ,JX_Int i)
{
   JX_Int VecStart;
   if(jx_hpIfReorderMatvecCommPkg(comm_pkg) )
   {
      JX_Int k = jx_hpRecvReorderMapIndex(comm_pkg, i);
      VecStart = comm_pkg -> recv_vec_starts[k];
   }
   else
   {
      VecStart = comm_pkg -> recv_vec_starts[i];
   }
   return VecStart;
}

JX_Int jx_hpCSRCommPkgRSendMapStartReorder(jx_ParCSRCommPkg *comm_pkg ,JX_Int i)
{
   JX_Int VecStart;
   if(jx_hpIfReorderMatvecCommPkg(comm_pkg) )
   {
      JX_Int k = jx_hpSendReorderMapIndex(comm_pkg, i);
      VecStart = comm_pkg -> send_map_starts[k];
   }
   else
   {
      VecStart = comm_pkg -> send_map_starts[i];
   }
   return VecStart;
}

JX_Int jx_hpCSRCommCommPkgSendReorder(JX_Int num_sends,jx_HardwareInfo  *hardware,
                                       JX_Int *SendProcs,  JX_Int **SendProcsReorder, JX_Int **send_reorder_map)
{
   JX_Int i,ip;
   JX_Int num_ip_cpu = 0;
   JX_Int num_ip_node = 0;
   JX_Int j_cpu = 0;
   JX_Int j_node = 0;
   JX_Int j_out_node = 0;

   JX_Int  proc_cpu_start   = jx_hpProcCpuStart(hardware);
   JX_Int  proc_cpu_end     = jx_hpProcCpuEnd(hardware);
   JX_Int  proc_node_start  = jx_hpProcNodeStart(hardware);
   JX_Int  proc_node_end    = jx_hpProcNodeEnd(hardware);
   JX_Int  bind_type        = jx_hpBindType(hardware);

   JX_Int  *send_procs = SendProcs;
   JX_Int  *ip_send_reorder = jx_CTAlloc(JX_Int,num_sends); 
   JX_Int  *ip_send_reorder_map = jx_CTAlloc(JX_Int,num_sends); 

   for (i = 0; i < num_sends; i ++)
   {
      ip        = send_procs[i]; 
      if(bind_type ==2 && ip >= proc_cpu_start && ip <= proc_cpu_end) num_ip_cpu++;
      else if(bind_type >=1 && ip >= proc_node_start && ip <= proc_node_end) num_ip_node++;
   }

   for (i = 0; i < num_sends; i ++)
   {
      ip        = send_procs[i];
      #if JX_REODER_SEND_CLOSER 
      if(num_ip_cpu !=0 && ip >= proc_cpu_start && ip <= proc_cpu_end)
      {
         // ip_send_reorder[num_sends - num_ip_cpu + j_cpu] = ip;
         // ip_send_reorder_map[num_sends - num_ip_cpu + j_cpu] = i;
         ip_send_reorder[j_cpu] = ip;
         ip_send_reorder_map[j_cpu] = i;
         j_cpu++;
      }
      else if(num_ip_node !=0 && ip >= proc_node_start && ip <= proc_node_end)
      {
         // ip_send_reorder[num_sends - num_ip_cpu - num_ip_node + j_node] = ip;
         // ip_send_reorder_map[num_sends - num_ip_cpu - num_ip_node + j_node] = i;
         ip_send_reorder[num_ip_cpu + j_node] = ip;
         ip_send_reorder_map[num_ip_cpu + j_node] = i;
         j_node++;
      }
      else
      {
         // ip_send_reorder[j_out_node] = ip;
         // ip_send_reorder_map[j_out_node] = i;
         ip_send_reorder[num_ip_cpu + num_ip_node + j_out_node] = ip;
         ip_send_reorder_map[num_ip_cpu + num_ip_node + j_out_node] = i;
         j_out_node++;
      }
      #else
      if(num_ip_cpu !=0 && ip >= proc_cpu_start && ip <= proc_cpu_end)
      {
         ip_send_reorder[num_sends - num_ip_cpu + j_cpu] = ip;
         ip_send_reorder_map[num_sends - num_ip_cpu + j_cpu] = i;
         j_cpu++;
      }
      else if(num_ip_node !=0 && ip >= proc_node_start && ip <= proc_node_end)
      {
         ip_send_reorder[num_sends - num_ip_cpu - num_ip_node + j_node] = ip;
         ip_send_reorder_map[num_sends - num_ip_cpu - num_ip_node + j_node] = i;
         j_node++;
      }
      else
      {
         ip_send_reorder[j_out_node] = ip;
         ip_send_reorder_map[j_out_node] = i;
         j_out_node++;
      }
      #endif
   }
   // jx_ParCSRCommPkgSendProcs(comm_pkg) = ip_send_reorder;
   *SendProcsReorder = ip_send_reorder;
   *send_reorder_map = ip_send_reorder_map;
   return 0;
}

JX_Int jx_hpCSRCommCommPkgRecvReorder(JX_Int num_recvs,jx_HardwareInfo  *hardware, 
                                       JX_Int *RecvProcs, JX_Int **RecvProcsReorder, JX_Int **recv_reorder_map)
{
   JX_Int i,ip;
   JX_Int num_ip_cpu = 0;
   JX_Int num_ip_node = 0;
   JX_Int j_cpu = 0;
   JX_Int j_node = 0;
   JX_Int j_out_node = 0;

   JX_Int  proc_cpu_start   = jx_hpProcCpuStart(hardware);
   JX_Int  proc_cpu_end     = jx_hpProcCpuEnd(hardware);
   JX_Int  proc_node_start  = jx_hpProcNodeStart(hardware);
   JX_Int  proc_node_end    = jx_hpProcNodeEnd(hardware);
   JX_Int  bind_type        = jx_hpBindType(hardware);

   JX_Int  *ip_recv_reorder = jx_CTAlloc(JX_Int,num_recvs);
   JX_Int  *ip_recv_reorder_map = jx_CTAlloc(JX_Int,num_recvs);
   JX_Int  *recv_procs = RecvProcs;

   for (i = 0; i < num_recvs; i ++)
   {
      ip        = recv_procs[i]; 
      if(bind_type ==2 && ip >= proc_cpu_start && ip <= proc_cpu_end) num_ip_cpu++;
      else if(bind_type >=1 && ip >= proc_node_start && ip <= proc_node_end) num_ip_node++;
   }
   for (i = 0; i < num_recvs; i ++)
   {
      ip        = recv_procs[i]; 
      #if JX_REODER_RECV_CLOSER
      if(num_ip_cpu !=0 && ip >= proc_cpu_start && ip <= proc_cpu_end)
      {
         // ip_recv_reorder[num_recvs - num_ip_cpu + j_cpu] = ip;
         // ip_recv_reorder_map[num_recvs - num_ip_cpu + j_cpu] = i;
         ip_recv_reorder[j_cpu] = ip;
         ip_recv_reorder_map[j_cpu] = i;
         j_cpu++;
      }
      else if(num_ip_node !=0 && ip >= proc_node_start && ip <= proc_node_end)
      {
         // ip_recv_reorder[num_recvs - num_ip_cpu - num_ip_node + j_node] = ip;
         // ip_recv_reorder_map[num_recvs - num_ip_cpu - num_ip_node + j_node] = i;
         ip_recv_reorder[num_ip_cpu + j_node] = ip;
         ip_recv_reorder_map[num_ip_cpu + j_node] = i;         
         j_node++;
      }
      else
      {
         ip_recv_reorder[num_ip_cpu + num_ip_node + j_out_node] = ip;
         ip_recv_reorder_map[num_ip_cpu + num_ip_node + j_out_node] = i;
         j_out_node++;
      }      
      #else
      if(num_ip_cpu !=0 && ip >= proc_cpu_start && ip <= proc_cpu_end)
      {
         ip_recv_reorder[num_recvs - num_ip_cpu + j_cpu] = ip;
         ip_recv_reorder_map[num_recvs - num_ip_cpu + j_cpu] = i;
         j_cpu++;
      }
      else if(num_ip_node !=0 && ip >= proc_node_start && ip <= proc_node_end)
      {
         ip_recv_reorder[num_recvs - num_ip_cpu - num_ip_node + j_node] = ip;
         ip_recv_reorder_map[num_recvs - num_ip_cpu - num_ip_node + j_node] = i;
         j_node++;
      }
      else
      {
         ip_recv_reorder[j_out_node] = ip;
         ip_recv_reorder_map[j_out_node] = i;
         j_out_node++;
      }
      #endif
   }
   // jx_ParCSRCommPkgRecvProcs(comm_pkg) = ip_recv_reorder;
   *RecvProcsReorder = ip_recv_reorder;
   *recv_reorder_map = ip_recv_reorder_map;
   return 0;
}


JX_Int jx_hpCSRCommPkgReorder(jx_ParCSRCommPkg *comm_pkg)
{
   if(!jx_hpIfReorderMatvecCommPkg(comm_pkg))
   {
      if(hp_hardware)
      {
         JX_Int *send_reorder;
         JX_Int *recv_reorder;
         JX_Int *send_reorder_map;
         JX_Int *recv_reorder_map;
         jx_hpCSRCommCommPkgSendReorder(jx_ParCSRCommPkgNumSends(comm_pkg), hp_hardware, jx_ParCSRCommPkgSendProcs(comm_pkg), &send_reorder, &send_reorder_map);
         jx_hpCSRCommCommPkgRecvReorder(jx_ParCSRCommPkgNumRecvs(comm_pkg), hp_hardware, jx_ParCSRCommPkgRecvProcs(comm_pkg), &recv_reorder, &recv_reorder_map);
         jx_hpSendReorderMap(comm_pkg) = send_reorder_map;
         jx_hpRecvReorderMap(comm_pkg) = recv_reorder_map;
         jx_hpSendReorderProcs(comm_pkg) = send_reorder;
         jx_hpRecvReorderProcs(comm_pkg) = recv_reorder;
         jx_hpIfReorderMatvecCommPkg(comm_pkg) = 1;
      }
      else
      {
         JX_Int myid;
         jx_MPI_Comm_rank(jx_ParCSRCommPkgComm(comm_pkg), &myid);
         if(myid == 0)
         {
            jx_printf("\nreorder matrix_comm_pkg need hardware info\n");
            jx_printf("This is not erro, but we recommand you to add hardwareinfo into matrix!\n");
         }
         jx_hpSendReorderMap(comm_pkg) = NULL;
         jx_hpRecvReorderMap(comm_pkg) = NULL;
      }
   }
   return 0;
}
#endif

#if JX_COMP_COMM_OVERLAP
/*!
 * \fn JX_Int jx_MatvecCPUCommPkgCreate
 * \brief Generates the comm_pkg for A. 
 *        If no row and / or column partitioning is given, 
 *        the routine determines them with MPE_Decomp1d.
 * \date 2011/09/07
 */
JX_Int
jx_hpMatvecCPUOutsideCommPkgCreate (jx_hpCSRMatrix *hp_A)
{
   JX_Int			 num_sends;
   JX_Int			*send_procs;
   JX_Int			*send_map_starts;
   JX_Int			*send_map_elmts;
 
   JX_Int			 num_recvs;
   JX_Int			*recv_procs;
   JX_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jx_hpCSRMatrixComm(hp_A);
   
   jx_ParCSRCommPkg	*comm_pkg;

   JX_Int   first_col_diag = jx_hpCSRMatrixFirstColDiag(hp_A);
   JX_Int  *col_map_offd   = jx_hpCSRMatrixColMapOffdlevel(jx_hpCSRMatrixCpuOutside(hp_A));
   JX_Int	num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdlevel(jx_hpCSRMatrixCpuOutside(hp_A)));


#if JX_NO_GLOBAL_PARTITION

   JX_Int  row_start = 0, row_end = 0;
   JX_Int  col_start = 0, col_end = 0;
   JX_Int  global_num_cols;
   
   jx_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jx_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jx_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jx_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jx_hpCSRMatrixAssumedPartition(hp_A);
   
  /*-----------------------------------------------------------
   * get commpkg info information 
   *----------------------------------------------------------*/

   jx_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
                             col_start, col_end, 
                             num_cols_offd, global_num_cols,
                             &num_recvs, &recv_procs, &recv_vec_starts,
                             &num_sends, &send_procs, &send_map_starts, 
                             &send_map_elmts, apart );
   
#else
   
   JX_Int *col_starts    = jx_hpCSRMatrixColStarts(hp_A);
   JX_Int	num_cols_diag = jx_CSRMatrixNumCols(jx_hpCSRMatrixDiag(hp_A));
  

   jx_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
                                num_cols_diag, num_cols_offd,
                                first_col_diag, col_map_offd, 1,
                                &num_recvs, &recv_procs, &recv_vec_starts,
                                &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
#endif


  /*-----------------------------------------------------------
   *  setup commpkg
   *----------------------------------------------------------*/

   comm_pkg = jx_CTAlloc(jx_ParCSRCommPkg, 1);

   jx_ParCSRCommPkgComm(comm_pkg) = comm;

   jx_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
   jx_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
   jx_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
   jx_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
   jx_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
   jx_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
   jx_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;
#if JX_REODER_SEND_RECV
   jx_hpIfReorderMatvecCommPkg(comm_pkg)   = 0;
#endif
   jx_hpCSRMatrixlevelCommPkg(jx_hpCSRMatrixCpuOutside(hp_A)) = comm_pkg;

   return jx_error_flag;
}

JX_Int
jx_hpMatvecNodeOutsideCommPkgCreate (jx_hpCSRMatrix *hp_A)
{
   JX_Int			 num_sends;
   JX_Int			*send_procs;
   JX_Int			*send_map_starts;
   JX_Int			*send_map_elmts;
 
   JX_Int			 num_recvs;
   JX_Int			*recv_procs;
   JX_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jx_hpCSRMatrixComm(hp_A);
   
   jx_ParCSRCommPkg	*comm_pkg;

   JX_Int   first_col_diag = jx_hpCSRMatrixFirstColDiag(hp_A);
   JX_Int  *col_map_offd   = jx_hpCSRMatrixColMapOffdlevel(jx_hpCSRMatrixNodeOutside(hp_A));
   //JX_Int	num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffd(hp_A));
   JX_Int	num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdlevel(jx_hpCSRMatrixNodeOutside(hp_A)));


#if JX_NO_GLOBAL_PARTITION

   JX_Int  row_start = 0, row_end = 0;
   JX_Int  col_start = 0, col_end = 0;
   JX_Int  global_num_cols;
   
   jx_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jx_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jx_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jx_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jx_hpCSRMatrixAssumedPartition(hp_A);
   
  /*-----------------------------------------------------------
   * get commpkg info information 
   *----------------------------------------------------------*/

   jx_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
                             col_start, col_end, 
                             num_cols_offd, global_num_cols,
                             &num_recvs, &recv_procs, &recv_vec_starts,
                             &num_sends, &send_procs, &send_map_starts, 
                             &send_map_elmts, apart );
   
#else
   
   JX_Int *col_starts    = jx_hpCSRMatrixColStarts(hp_A);
   JX_Int	num_cols_diag = jx_CSRMatrixNumCols(jx_hpCSRMatrixDiag(hp_A));
  

   jx_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
                                num_cols_diag, num_cols_offd,
                                first_col_diag, col_map_offd, 1,
                                &num_recvs, &recv_procs, &recv_vec_starts,
                                &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
#endif


  /*-----------------------------------------------------------
   *  setup commpkg
   *----------------------------------------------------------*/

   comm_pkg = jx_CTAlloc(jx_ParCSRCommPkg, 1);

   jx_ParCSRCommPkgComm(comm_pkg) = comm;

   jx_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
   jx_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
   jx_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
   jx_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
   jx_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
   jx_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
   jx_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;
#if JX_REODER_SEND_RECV
   jx_hpIfReorderMatvecCommPkg(comm_pkg)   = 0;
#endif
   jx_hpCSRMatrixlevelCommPkg(jx_hpCSRMatrixNodeOutside(hp_A)) = comm_pkg;
   //jx_printf("jx_ParCSRCommPkgNumSends is %d\n",num_sends);
   return jx_error_flag;
}

JX_Int
jx_hpMatvecCoreOutsideCommPkgCreate (jx_hpCSRMatrix *hp_A)
{
   JX_Int			 num_sends;
   JX_Int			*send_procs;
   JX_Int			*send_map_starts;
   JX_Int			*send_map_elmts;
 
   JX_Int			 num_recvs;
   JX_Int			*recv_procs;
   JX_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jx_hpCSRMatrixComm(hp_A);
   
   jx_ParCSRCommPkg	*comm_pkg;

   JX_Int   first_col_diag = jx_hpCSRMatrixFirstColDiag(hp_A);
   JX_Int  *col_map_offd   = jx_hpCSRMatrixColMapOffdlevel(jx_hpCSRMatrixCpu(hp_A));
   JX_Int	num_cols_offd  = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdlevel(jx_hpCSRMatrixCpu(hp_A)));


#if JX_NO_GLOBAL_PARTITION

   JX_Int  row_start = 0, row_end = 0;
   JX_Int  col_start = 0, col_end = 0;
   JX_Int  global_num_cols;
   
   jx_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jx_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jx_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jx_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jx_hpCSRMatrixAssumedPartition(hp_A);
   
  /*-----------------------------------------------------------
   * get commpkg info information 
   *----------------------------------------------------------*/

   jx_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
                             col_start, col_end, 
                             num_cols_offd, global_num_cols,
                             &num_recvs, &recv_procs, &recv_vec_starts,
                             &num_sends, &send_procs, &send_map_starts, 
                             &send_map_elmts, apart );
   
#else
   
   JX_Int *col_starts    = jx_hpCSRMatrixColStarts(hp_A);
   JX_Int	num_cols_diag = jx_CSRMatrixNumCols(jx_hpCSRMatrixDiag(hp_A));
  

   jx_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
                                num_cols_diag, num_cols_offd,
                                first_col_diag, col_map_offd, 1,
                                &num_recvs, &recv_procs, &recv_vec_starts,
                                &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
#endif


  /*-----------------------------------------------------------
   *  setup commpkg
   *----------------------------------------------------------*/

   comm_pkg = jx_CTAlloc(jx_ParCSRCommPkg, 1);

   jx_ParCSRCommPkgComm(comm_pkg) = comm;

   jx_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
   jx_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
   jx_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
   jx_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
   jx_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
   jx_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
   jx_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;
#if JX_REODER_SEND_RECV
   jx_hpIfReorderMatvecCommPkg(comm_pkg)   = 0;
#endif
   jx_hpCSRMatrixlevelCommPkg(jx_hpCSRMatrixCpu(hp_A)) = comm_pkg;
   return jx_error_flag;
}
#endif