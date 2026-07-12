//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hp_csr_com.c -- communication operations for hierarchy parallel matrices and vectors.
 *  Date: 2011/09/07
 */ 

#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

// JXF_Int
// jxf_hpMatvecCommPkgCreate (jxf_hpCSRMatrix *hp_A, jxf_HardwareInfo *hardware)
// {
//    JXF_Int			 num_sends;
//    JXF_Int			*send_procs;
//    JXF_Int			*send_map_starts;
//    JXF_Int			*send_map_elmts;
 
//    JXF_Int			 num_recvs;
//    JXF_Int			*recv_procs;
//    JXF_Int			*recv_vec_starts;
   
//    MPI_Comm              comm = jxf_hpCSRMatrixComm(hp_A);
//    jxf_ParCSRCommPkg	*comm_pkg;

//    JXF_Int   first_col_diag = jxf_hpCSRMatrixFirstColDiag(hp_A);
//    JXF_Int  *col_map_offd   = jxf_hpCSRMatrixColMapOffd(hp_A);
//    JXF_Int	 num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffd(hp_A));


// #if JXF_NO_GLOBAL_PARTITION

//    JXF_Int  row_start = 0, row_end = 0;
//    JXF_Int  col_start = 0, col_end = 0;
//    JXF_Int  global_num_cols;
   
//    jxf_IJAssumedPart  *apart;
   
//    /*  get parcsr_A information */
//    jxf_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
//    global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A); 

//    /* Create the assumed partition */
//    if (jxf_hpCSRMatrixAssumedPartition(hp_A) == NULL)
//    {
//       jxf_hpCSRMatrixCreateAssumedPartition(hp_A);
//    }
//    apart = jxf_hpCSRMatrixAssumedPartition(hp_A);
   
//   /*-----------------------------------------------------------
//    * get commpkg info information 
//    *----------------------------------------------------------*/

//    jxf_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
//                              col_start, col_end, 
//                              num_cols_offd, global_num_cols,
//                              &num_recvs, &recv_procs, &recv_vec_starts,
//                              &num_sends, &send_procs, &send_map_starts, 
//                              &send_map_elmts, apart );
   
// #else
   
//    JXF_Int *col_starts    = jxf_hpCSRMatrixColStarts(hp_A);
//    JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_A));
  

//    jxf_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
//                                 num_cols_diag, num_cols_offd,
//                                 first_col_diag, col_map_offd, 1,
//                                 &num_recvs, &recv_procs, &recv_vec_starts,
//                                 &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
// #endif


//   /*-----------------------------------------------------------
//    *  setup commpkg
//    *----------------------------------------------------------*/

//    comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);

//    jxf_ParCSRCommPkgComm(comm_pkg) = comm;

//   /*-----------------------------------------------------------
//    *  reoderMatvecCommPkg
//    *----------------------------------------------------------*/

//    if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
//    {
//       jxf_hpCSRCommCommPkgSendReorder(num_sends, hardware, &send_procs);
//       jxf_hpCSRCommCommPkgRecvReorder(num_recvs, hardware, &recv_procs);
//       jxf_hpIfReorderMatvecCommPkg(comm_pkg) = 1;
//    }

//    jxf_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
//    jxf_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
//    jxf_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
//    jxf_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
//    jxf_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
//    jxf_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
//    jxf_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;

//    jxf_hpCSRMatrixCommPkg(hp_A) = comm_pkg;

//    return jxf_error_flag;
// }

// JXF_Int
// jxf_hpMatvecCommPkgCreateLevel (jxf_hpCSRMatrix *hp_A, JXF_Int Level, jxf_HardwareInfo *hardware)
// {
//    JXF_Int			 num_sends;
//    JXF_Int			*send_procs;
//    JXF_Int			*send_map_starts;
//    JXF_Int			*send_map_elmts;
 
//    JXF_Int			 num_recvs;
//    JXF_Int			*recv_procs;
//    JXF_Int			*recv_vec_starts;

//    JXF_Int	      num_cols_offd = 0;
//    JXF_Int         *col_map_offd = NULL;
   
//    MPI_Comm              comm = jxf_hpCSRMatrixComm(hp_A);
//    jxf_ParCSRCommPkg	*comm_pkg;

//    if(Level == 1)
//    {
//       col_map_offd   = jxf_hpCSRMatrixColMapOffdNode(hp_A);
//       //num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdNode(hp_A));
//    }
//    else if(Level == 2)
//    {
//       col_map_offd   = jxf_hpCSRMatrixColMapOffdCpu(hp_A);
//       //num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdCpu(hp_A));
//    }
//    else
//    {
//       jxf_printf("level must be 1 or 2\n");
//       return jxf_error_flag;
//    }
//    num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffd(hp_A));

//    // JXF_Int myid ;
//    // jxf_MPI_Comm_rank(comm, &myid);
//    // JXF_Int i;
//    // if(myid == 0)
//    // {
//    //    jxf_printf("\n myid is %d, col_map_offd is ", myid);
//    //    for(i = 0; i < num_cols_offd; i++)
//    //    {
//    //       jxf_printf("%d\t",col_map_offd[i]);
//    //    }
//    //    jxf_printf("\n");
//    // }
   
//    JXF_Int   first_col_diag = jxf_hpCSRMatrixFirstColDiag(hp_A);

// #if JXF_NO_GLOBAL_PARTITION

//    JXF_Int  row_start = 0, row_end = 0;
//    JXF_Int  col_start = 0, col_end = 0;
//    JXF_Int  global_num_cols;
   
//    jxf_IJAssumedPart  *apart;
   
//    /*  get parcsr_A information */
//    jxf_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
//    global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A); 

//    /* Create the assumed partition */
//    if (jxf_hpCSRMatrixAssumedPartition(hp_A) == NULL)
//    {
//       jxf_hpCSRMatrixCreateAssumedPartition(hp_A);
//    }
//    apart = jxf_hpCSRMatrixAssumedPartition(hp_A);
   
//   /*-----------------------------------------------------------
//    * get commpkg info information 
//    *----------------------------------------------------------*/

//    jxf_NewCommPkgCreate_core( comm, col_map_offd, first_col_diag, 
//                              col_start, col_end, 
//                              num_cols_offd, global_num_cols,
//                              &num_recvs, &recv_procs, &recv_vec_starts,
//                              &num_sends, &send_procs, &send_map_starts, 
//                              &send_map_elmts, apart );
   
// #else
   
//    // JXF_Int *col_starts    = NULL;
//    // if(Level == 1)
//    // {
//    //    col_starts = jxf_hpCSRMatrixColNodeStart(hp_A);
//    // }
//    // else if(Level == 2)
//    // {
//    //    col_starts = jxf_hpCSRMatrixColCpuStart(hp_A);
//    // }
//    JXF_Int *col_starts    = jxf_hpCSRMatrixColStarts(hp_A);
//    JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_A));
  

// // jxf_printf("before run success\n");

   
//    //if(myid == 0) jxf_printf("num_cols_diag is %d, num_cols_offd is %d\n",num_cols_diag, num_cols_offd);

//    jxf_MatvecCommPkgCreate_core( comm, col_map_offd, first_col_diag, col_starts,
//                                 num_cols_diag, num_cols_offd,
//                                 first_col_diag, col_map_offd, 
//                                 1,
//                                 &num_recvs, &recv_procs, &recv_vec_starts,
//                                 &num_sends, &send_procs, &send_map_starts, &send_map_elmts );
// #endif

// jxf_printf("jxf_MatvecCommPkgCreate_core run success\n");
//   /*-----------------------------------------------------------
//    *  setup commpkg
//    *----------------------------------------------------------*/

//    comm_pkg = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);

//    jxf_ParCSRCommPkgComm(comm_pkg) = comm;
   
//    jxf_ParCSRCommPkgNumRecvs(comm_pkg)      = num_recvs;
//    jxf_ParCSRCommPkgRecvProcs(comm_pkg)     = recv_procs;
//    jxf_ParCSRCommPkgRecvVecStarts(comm_pkg) = recv_vec_starts;
//    jxf_ParCSRCommPkgNumSends(comm_pkg)      = num_sends;
//    jxf_ParCSRCommPkgSendProcs(comm_pkg)     = send_procs;
//    jxf_ParCSRCommPkgSendMapStarts(comm_pkg) = send_map_starts;
//    jxf_ParCSRCommPkgSendMapElmts(comm_pkg)  = send_map_elmts;

//    //jxf_hpCSRMatrixCommPkg(hp_A) = comm_pkg;
//    if(Level == 1)
//    {
//       jxf_hpCSRMatrixNodeCommPkg(hp_A) = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
//       jxf_hpCSRMatrixNodeCommPkg(hp_A) = comm_pkg;
//    }
//    else if(Level == 2) 
//    {
//       jxf_hpCSRMatrixCpuCommPkg(hp_A) = jxf_CTAlloc(jxf_ParCSRCommPkg, 1);
//       jxf_hpCSRMatrixCpuCommPkg(hp_A) = comm_pkg;
//    }

//    return jxf_error_flag;
// }

// jxf_ParCSRCommHandle *
// jxf_hpCSRCommHandleCreate(  JXF_Int 	             job,
// 			                  jxf_ParCSRCommPkg *comm_pkg,
//                            void             *send_data, 
//                            void             *recv_data, 
//                            jxf_HardwareInfo *hardware)
// {
//    JXF_Int                  num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
//    JXF_Int                  num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
//    MPI_Comm             comm      = jxf_ParCSRCommPkgComm(comm_pkg);

//    jxf_ParCSRCommHandle *comm_handle;
//    JXF_Int                  num_requests;
//    MPI_Request         *requests;

//    JXF_Int                  i, j;
//    JXF_Int			my_id, num_procs;
//    JXF_Int			ip, vec_start, vec_len;

//    num_requests = num_sends + num_recvs;
//    requests = jxf_CTAlloc(MPI_Request, num_requests);
 
//    jxf_MPI_Comm_size(comm, &num_procs);
//    jxf_MPI_Comm_rank(comm, &my_id);

//    j = 0;
//    switch (job)
//    {
//       case  1:
//       {
//          JXF_Real *d_send_data = (JXF_Real *) send_data;
//          JXF_Real *d_recv_data = (JXF_Real *) recv_data;
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jxf_MPI_Irecv(&d_recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
// 	         vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
// 	         vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//    	      jxf_MPI_Isend(&d_send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//       }
//       break;
      
      
//       case  2:
//       {
//          JXF_Real *d_send_data = (JXF_Real *) send_data;
//          JXF_Real *d_recv_data = (JXF_Real *) recv_data;
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jxf_MPI_Irecv(&d_recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jxf_MPI_Isend(&d_send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
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
//       //       vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//       //       vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//       //       jxf_MPI_Irecv(&d_recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
//       //    }

//       //    for (i = 0; i < num_sends; i ++)
//       //    {
//       //       ip        = ip_send_reorder[i];
//       //       vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
//       //       vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//       //       jxf_MPI_Isend(&d_send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
//       //    }
//       // }
//       // break;
      
//       case  11:
//       {
//          JXF_Int *i_send_data = (JXF_Int *) send_data;
//          JXF_Int *i_recv_data = (JXF_Int *) recv_data;
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jxf_MPI_Irecv(&i_recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jxf_MPI_Isend(&i_send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//       }
//       break;
      
      
//       case  12:
//       {
//          JXF_Int *i_send_data = (JXF_Int *) send_data;
//          JXF_Int *i_recv_data = (JXF_Int *) recv_data;
//          for (i = 0; i < num_sends; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgSendProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
//             vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
//             jxf_MPI_Irecv(&i_recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          for (i = 0; i < num_recvs; i ++)
//          {
//             ip        = jxf_ParCSRCommPkgRecvProc(comm_pkg, i); 
//             vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
//             vec_len   = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
//             jxf_MPI_Isend(&i_send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
//          }
//          break;
//       }
//    }
 
 
//   /*--------------------------------------------------------------------
//    * set up comm_handle and return
//    *--------------------------------------------------------------------*/

//    comm_handle = jxf_CTAlloc(jxf_ParCSRCommHandle, 1);

//    jxf_ParCSRCommHandleCommPkg(comm_handle)     = comm_pkg;
//    jxf_ParCSRCommHandleSendData(comm_handle)    = send_data;
//    jxf_ParCSRCommHandleRecvData(comm_handle)    = recv_data;
//    jxf_ParCSRCommHandleNumRequests(comm_handle) = num_requests;
//    jxf_ParCSRCommHandleRequests(comm_handle)    = requests;

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
#if JXF_REODER_SEND_RECV

JXF_Int jxf_hpCSRCommPkgRecvVecStartReorder(jxf_ParCSRCommPkg *comm_pkg ,JXF_Int i)
{
   JXF_Int VecStart;
   if(jxf_hpIfReorderMatvecCommPkg(comm_pkg) )
   {
      JXF_Int k = jxf_hpRecvReorderMapIndex(comm_pkg, i);
      VecStart = comm_pkg -> recv_vec_starts[k];
   }
   else
   {
      VecStart = comm_pkg -> recv_vec_starts[i];
   }
   return VecStart;
}

JXF_Int jxf_hpCSRCommPkgRSendMapStartReorder(jxf_ParCSRCommPkg *comm_pkg ,JXF_Int i)
{
   JXF_Int VecStart;
   if(jxf_hpIfReorderMatvecCommPkg(comm_pkg) )
   {
      JXF_Int k = jxf_hpSendReorderMapIndex(comm_pkg, i);
      VecStart = comm_pkg -> send_map_starts[k];
   }
   else
   {
      VecStart = comm_pkg -> send_map_starts[i];
   }
   return VecStart;
}

JXF_Int jxf_hpCSRCommCommPkgSendReorder(JXF_Int num_sends,jxf_HardwareInfo  *hardware,
                                       JXF_Int *SendProcs,  JXF_Int **SendProcsReorder, JXF_Int **send_reorder_map)
{
   JXF_Int i,ip;
   JXF_Int num_ip_cpu = 0;
   JXF_Int num_ip_node = 0;
   JXF_Int j_cpu = 0;
   JXF_Int j_node = 0;
   JXF_Int j_out_node = 0;

   JXF_Int  proc_cpu_start   = jxf_hpProcCpuStart(hardware);
   JXF_Int  proc_cpu_end     = jxf_hpProcCpuEnd(hardware);
   JXF_Int  proc_node_start  = jxf_hpProcNodeStart(hardware);
   JXF_Int  proc_node_end    = jxf_hpProcNodeEnd(hardware);
   JXF_Int  bind_type        = jxf_hpBindType(hardware);

   JXF_Int  *send_procs = SendProcs;
   JXF_Int  *ip_send_reorder = jxf_CTAlloc(JXF_Int,num_sends); 
   JXF_Int  *ip_send_reorder_map = jxf_CTAlloc(JXF_Int,num_sends); 

   for (i = 0; i < num_sends; i ++)
   {
      ip        = send_procs[i]; 
      if(bind_type ==2 && ip >= proc_cpu_start && ip <= proc_cpu_end) num_ip_cpu++;
      else if(bind_type >=1 && ip >= proc_node_start && ip <= proc_node_end) num_ip_node++;
   }

   for (i = 0; i < num_sends; i ++)
   {
      ip        = send_procs[i];
      #if JXF_REODER_SEND_CLOSER 
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
   // jxf_ParCSRCommPkgSendProcs(comm_pkg) = ip_send_reorder;
   *SendProcsReorder = ip_send_reorder;
   *send_reorder_map = ip_send_reorder_map;
   return 0;
}

JXF_Int jxf_hpCSRCommCommPkgRecvReorder(JXF_Int num_recvs,jxf_HardwareInfo  *hardware, 
                                       JXF_Int *RecvProcs, JXF_Int **RecvProcsReorder, JXF_Int **recv_reorder_map)
{
   JXF_Int i,ip;
   JXF_Int num_ip_cpu = 0;
   JXF_Int num_ip_node = 0;
   JXF_Int j_cpu = 0;
   JXF_Int j_node = 0;
   JXF_Int j_out_node = 0;

   JXF_Int  proc_cpu_start   = jxf_hpProcCpuStart(hardware);
   JXF_Int  proc_cpu_end     = jxf_hpProcCpuEnd(hardware);
   JXF_Int  proc_node_start  = jxf_hpProcNodeStart(hardware);
   JXF_Int  proc_node_end    = jxf_hpProcNodeEnd(hardware);
   JXF_Int  bind_type        = jxf_hpBindType(hardware);

   JXF_Int  *ip_recv_reorder = jxf_CTAlloc(JXF_Int,num_recvs);
   JXF_Int  *ip_recv_reorder_map = jxf_CTAlloc(JXF_Int,num_recvs);
   JXF_Int  *recv_procs = RecvProcs;

   for (i = 0; i < num_recvs; i ++)
   {
      ip        = recv_procs[i]; 
      if(bind_type ==2 && ip >= proc_cpu_start && ip <= proc_cpu_end) num_ip_cpu++;
      else if(bind_type >=1 && ip >= proc_node_start && ip <= proc_node_end) num_ip_node++;
   }
   for (i = 0; i < num_recvs; i ++)
   {
      ip        = recv_procs[i]; 
      #if JXF_REODER_RECV_CLOSER
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
   // jxf_ParCSRCommPkgRecvProcs(comm_pkg) = ip_recv_reorder;
   *RecvProcsReorder = ip_recv_reorder;
   *recv_reorder_map = ip_recv_reorder_map;
   return 0;
}


JXF_Int jxf_hpCSRCommPkgReorder(jxf_ParCSRCommPkg *comm_pkg)
{
   if(!jxf_hpIfReorderMatvecCommPkg(comm_pkg))
   {
      if(jxf_hp_hardware)
      {
         JXF_Int *send_reorder;
         JXF_Int *recv_reorder;
         JXF_Int *send_reorder_map;
         JXF_Int *recv_reorder_map;
         jxf_hpCSRCommCommPkgSendReorder(jxf_ParCSRCommPkgNumSends(comm_pkg), jxf_hp_hardware, jxf_ParCSRCommPkgSendProcs(comm_pkg), &send_reorder, &send_reorder_map);
         jxf_hpCSRCommCommPkgRecvReorder(jxf_ParCSRCommPkgNumRecvs(comm_pkg), jxf_hp_hardware, jxf_ParCSRCommPkgRecvProcs(comm_pkg), &recv_reorder, &recv_reorder_map);
         jxf_hpSendReorderMap(comm_pkg) = send_reorder_map;
         jxf_hpRecvReorderMap(comm_pkg) = recv_reorder_map;
         jxf_hpSendReorderProcs(comm_pkg) = send_reorder;
         jxf_hpRecvReorderProcs(comm_pkg) = recv_reorder;
         jxf_hpIfReorderMatvecCommPkg(comm_pkg) = 1;
      }
      else
      {
         JXF_Int myid;
         jxf_MPI_Comm_rank(jxf_ParCSRCommPkgComm(comm_pkg), &myid);
         if(myid == 0)
         {
            jxf_printf("\nreorder matrix_comm_pkg need hardware info\n");
            jxf_printf("This is not erro, but we recommand you to add hardwareinfo into matrix!\n");
         }
         jxf_hpSendReorderMap(comm_pkg) = NULL;
         jxf_hpRecvReorderMap(comm_pkg) = NULL;
      }
   }
   return 0;
}
#endif

#if JXF_COMP_COMM_OVERLAP
/*!
 * \fn JXF_Int jxf_MatvecCPUCommPkgCreate
 * \brief Generates the comm_pkg for A. 
 *        If no row and / or column partitioning is given, 
 *        the routine determines them with MPE_Decomp1d.
 * \date 2011/09/07
 */
JXF_Int
jxf_hpMatvecCPUOutsideCommPkgCreate (jxf_hpCSRMatrix *hp_A)
{
   JXF_Int			 num_sends;
   JXF_Int			*send_procs;
   JXF_Int			*send_map_starts;
   JXF_Int			*send_map_elmts;
 
   JXF_Int			 num_recvs;
   JXF_Int			*recv_procs;
   JXF_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jxf_hpCSRMatrixComm(hp_A);
   
   jxf_ParCSRCommPkg	*comm_pkg;

   JXF_Int   first_col_diag = jxf_hpCSRMatrixFirstColDiag(hp_A);
   JXF_Int  *col_map_offd   = jxf_hpCSRMatrixColMapOffdlevel(jxf_hpCSRMatrixCpuOutside(hp_A));
   JXF_Int	num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdlevel(jxf_hpCSRMatrixCpuOutside(hp_A)));


#if JXF_NO_GLOBAL_PARTITION

   JXF_Int  row_start = 0, row_end = 0;
   JXF_Int  col_start = 0, col_end = 0;
   JXF_Int  global_num_cols;
   
   jxf_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jxf_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jxf_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jxf_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jxf_hpCSRMatrixAssumedPartition(hp_A);
   
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
   
   JXF_Int *col_starts    = jxf_hpCSRMatrixColStarts(hp_A);
   JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_A));
  

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
#endif
   jxf_hpCSRMatrixlevelCommPkg(jxf_hpCSRMatrixCpuOutside(hp_A)) = comm_pkg;

   return jxf_error_flag;
}

JXF_Int
jxf_hpMatvecNodeOutsideCommPkgCreate (jxf_hpCSRMatrix *hp_A)
{
   JXF_Int			 num_sends;
   JXF_Int			*send_procs;
   JXF_Int			*send_map_starts;
   JXF_Int			*send_map_elmts;
 
   JXF_Int			 num_recvs;
   JXF_Int			*recv_procs;
   JXF_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jxf_hpCSRMatrixComm(hp_A);
   
   jxf_ParCSRCommPkg	*comm_pkg;

   JXF_Int   first_col_diag = jxf_hpCSRMatrixFirstColDiag(hp_A);
   JXF_Int  *col_map_offd   = jxf_hpCSRMatrixColMapOffdlevel(jxf_hpCSRMatrixNodeOutside(hp_A));
   //JXF_Int	num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffd(hp_A));
   JXF_Int	num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdlevel(jxf_hpCSRMatrixNodeOutside(hp_A)));


#if JXF_NO_GLOBAL_PARTITION

   JXF_Int  row_start = 0, row_end = 0;
   JXF_Int  col_start = 0, col_end = 0;
   JXF_Int  global_num_cols;
   
   jxf_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jxf_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jxf_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jxf_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jxf_hpCSRMatrixAssumedPartition(hp_A);
   
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
   
   JXF_Int *col_starts    = jxf_hpCSRMatrixColStarts(hp_A);
   JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_A));
  

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
#endif
   jxf_hpCSRMatrixlevelCommPkg(jxf_hpCSRMatrixNodeOutside(hp_A)) = comm_pkg;
   //jxf_printf("jxf_ParCSRCommPkgNumSends is %d\n",num_sends);
   return jxf_error_flag;
}

JXF_Int
jxf_hpMatvecCoreOutsideCommPkgCreate (jxf_hpCSRMatrix *hp_A)
{
   JXF_Int			 num_sends;
   JXF_Int			*send_procs;
   JXF_Int			*send_map_starts;
   JXF_Int			*send_map_elmts;
 
   JXF_Int			 num_recvs;
   JXF_Int			*recv_procs;
   JXF_Int			*recv_vec_starts;
   
   MPI_Comm              comm = jxf_hpCSRMatrixComm(hp_A);
   
   jxf_ParCSRCommPkg	*comm_pkg;

   JXF_Int   first_col_diag = jxf_hpCSRMatrixFirstColDiag(hp_A);
   JXF_Int  *col_map_offd   = jxf_hpCSRMatrixColMapOffdlevel(jxf_hpCSRMatrixCpu(hp_A));
   JXF_Int	num_cols_offd  = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdlevel(jxf_hpCSRMatrixCpu(hp_A)));


#if JXF_NO_GLOBAL_PARTITION

   JXF_Int  row_start = 0, row_end = 0;
   JXF_Int  col_start = 0, col_end = 0;
   JXF_Int  global_num_cols;
   
   jxf_IJAssumedPart  *apart;
   
   /*  get parcsr_A information */
   jxf_hpCSRMatrixGetLocalRange( hp_A, &row_start, &row_end, &col_start, &col_end );
   global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A); 

   /* Create the assumed partition */
   if (jxf_hpCSRMatrixAssumedPartition(hp_A) == NULL)
   {
      jxf_hpCSRMatrixCreateAssumedPartition(hp_A);
   }
   apart = jxf_hpCSRMatrixAssumedPartition(hp_A);
   
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
   
   JXF_Int *col_starts    = jxf_hpCSRMatrixColStarts(hp_A);
   JXF_Int	num_cols_diag = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_A));
  

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
#endif
   jxf_hpCSRMatrixlevelCommPkg(jxf_hpCSRMatrixCpu(hp_A)) = comm_pkg;
   return jxf_error_flag;
}
#endif