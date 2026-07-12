#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

#ifndef JX_ENV_VARIABLE_H
#include "jx_env_variable.h"
#endif

// add global variable hp_hardware
jx_HardwareInfo *hp_hardware = NULL;

#ifdef USING_HWLOC 
#include "hwloc.h"

/*!
 * \fn jx_HardwareInfo *jx_hpCreateHardwareInfo
 * \brief Create Hardware Info. 
 * \date 2021/10/15 create
 * \date 2022/11/24 change hardware to global variable
 * \author mrz
 */
JX_Int jx_hpCreateHardwareInfo(MPI_Comm comm)
{
   if(hp_hardware == NULL)
   {
      hp_hardware = jx_CTAlloc(jx_HardwareInfo, 1);
      //add process binding operation
      /*
      0:bind to node
      1:bind to cpu
      2:bind to core
      */
      jx_hpBindType(hp_hardware) = jx_getenv(Bind_Type);

      jx_hpComm(hp_hardware)           = comm;
      jx_hpCoreComm(hp_hardware)       = MPI_COMM_NULL;
      jx_hpIntraCpuComm(hp_hardware)	= MPI_COMM_NULL;
      jx_hpNodeComm(hp_hardware)	      = MPI_COMM_NULL;
      //These comm is not used
      // jx_hpCpuComm(hp_hardware)        = MPI_COMM_NULL;
      // jx_hpIntraNodeComm(hp_hardware)	= MPI_COMM_NULL;

      jx_hpBindToType_Hwloc(hp_hardware);

   }

   return 0;
}

/*!
 * \fn JX_Int jx_hpGetLevelBindId
 * \brief use api of mpi according myid bind to type. 
 * \date 2022/11/30
 * \author mrz
 */
JX_Int jx_hpGetLevelBindId(JX_Int myid, JX_Int numprocess, JX_Int num_level, 
                           JX_Int *id_in_next_level, JX_Int *process_id_in_level, JX_Int *numprocess_in_level, JX_Int* start_process_id_in_level, JX_Int *end_process_id_in_level)
{
   JX_Int process_per_in_level;
   JX_Int r_process_per_in_level;
   JX_Int my_id_in_level, my_process_id_in_level, my_numprocess_in_level;
   JX_Int my_start_process_id_in_level, my_end_process_id_in_level;

   process_per_in_level   = numprocess / num_level;
   r_process_per_in_level = numprocess - process_per_in_level * num_level;

   if(myid < r_process_per_in_level *(process_per_in_level + 1))
   {
      my_numprocess_in_level       = process_per_in_level + 1;
      my_id_in_level               = myid / (my_numprocess_in_level);
      my_process_id_in_level       = myid % (my_numprocess_in_level);
   }
   else
   {
      my_numprocess_in_level = process_per_in_level;
      my_id_in_level = r_process_per_in_level + (myid - r_process_per_in_level * (process_per_in_level + 1))/process_per_in_level ;
      my_process_id_in_level = (myid - r_process_per_in_level * (process_per_in_level + 1)) % process_per_in_level;
   }

   if(numprocess < num_level)
   {
      JX_Int remainder = num_level - numprocess;
      JX_Int remainder_in_process = remainder / numprocess  + 1;
      if(remainder % numprocess) remainder_in_process += 1;
      if (myid < remainder)
      {
         my_id_in_level = myid * remainder_in_process;
      }
      else
      {
         my_id_in_level = remainder * remainder_in_process + (myid - remainder) * (num_level/numprocess);
      }
   }

   my_start_process_id_in_level = myid - my_process_id_in_level;
   my_end_process_id_in_level   = my_start_process_id_in_level + my_numprocess_in_level - 1;

   *id_in_next_level    = my_id_in_level;
   if(process_id_in_level) *process_id_in_level = my_process_id_in_level;
   if(numprocess_in_level) *numprocess_in_level = my_numprocess_in_level;
   if(start_process_id_in_level) *start_process_id_in_level = my_start_process_id_in_level;
   if(end_process_id_in_level)   *end_process_id_in_level   = my_end_process_id_in_level;
   return 0;
}

/*!
 * \fn JX_Int jx_hpBindToType_Hwloc
 * \brief use api of hwloc according myid bind to type. 
 * \date 2021/08/26
 * \author mrz
 */

//TODO:Comments need to be removed

JX_Int jx_hpBindToType_Hwloc(jx_HardwareInfo *hardware)
{
   JX_Int myid, numprocess, global_pid;
   MPI_Comm comm = jx_hpComm(hardware);

   JX_Int type = jx_hpBindType(hardware);
   //JX_Int type_id_in_node, num_type_in_node; 
   JX_Int num_numanode_in_machine, num_core_in_numanode, num_machine, num_core_in_total_machine, num_numanode_in_total_machine;
   JX_Int process_per_numanode, process_per_core, my_machine_id, my_numanode_id, pid_in_numanode, pid_in_core, my_core_id;
   JX_Int process_statr_in_machine, process_statr_in_numanode;
   JX_Int process_end_in_machine, process_end_in_numanode;
   hwloc_topology_t topology;
   hwloc_cpuset_t cpuset;
   hwloc_obj_t obj_in_numanode, obj_in_core, obj_in_machine;
   JX_Int err, i, j, k, tmp_index;
   JX_Int ***hardware_ID;

   MPI_Comm shmcomm = MPI_COMM_NULL;
   // get communicator in the same node
   MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
   jx_MPI_Comm_rank(shmcomm, &myid);
   jx_MPI_Comm_size(shmcomm, &numprocess);
   jx_MPI_Comm_rank(comm, &global_pid);

   // char processor_name[MPI_MAX_PROCESSOR_NAME];
   // MPI_Get_processor_name(processor_name,&namelen);

   /*
      1. load topology
   */

   /* create a topology */
   err = hwloc_topology_init(&topology);
   if (err < 0) 
   {
      fprintf(stderr, "failed to initialize the topology\n");
      return EXIT_FAILURE;
   }

#if 0
   // hwloc-1.13.1 is not supporsed
   hwloc_topology_set_all_types_filter(topology, HWLOC_TYPE_FILTER_KEEP_NONE);
   hwloc_topology_set_type_filter(topology, HWLOC_OBJ_NUMANODE, HWLOC_TYPE_FILTER_KEEP_ALL);
   hwloc_topology_set_type_filter(topology, HWLOC_OBJ_CORE, HWLOC_TYPE_FILTER_KEEP_ALL);
   // hwloc_topology_set_type_filter(topology, HWLOC_OBJ_PACKAGE, HWLOC_TYPE_FILTER_KEEP_ALL);
   if(type == 0) hwloc_topology_set_type_filter(topology, HWLOC_OBJ_MACHINE, HWLOC_TYPE_FILTER_KEEP_ALL);
#endif

   err = hwloc_topology_load(topology);
   if (err < 0) 
   {
      fprintf(stderr, "failed to load the topology\n");
      hwloc_topology_destroy(topology);
      return EXIT_FAILURE;
   }

   /*
      2.  get haedware info
         include three level: Machine NumaNode Core
         note: there are all in the same node
   */

   num_machine          = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_MACHINE);
   hardware_ID          = jx_CTAlloc(int**, num_machine + 1);
   hardware_ID[0]       = jx_CTAlloc(int*, 1);
   hardware_ID[0][0]    = jx_CTAlloc(int, 1);
   hardware_ID[0][0][0] = num_machine;
   
   tmp_index                       = 0;
   num_core_in_numanode            = 0;
   num_numanode_in_machine         = 0;
   num_numanode_in_total_machine   = 0;
   num_core_in_total_machine       = 0;

   for(k = 0; k < num_machine; k++)
   {
      obj_in_machine              = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, k);
      num_numanode_in_machine     = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj_in_machine->cpuset, HWLOC_OBJ_NUMANODE);
      num_numanode_in_total_machine = num_numanode_in_total_machine + num_numanode_in_machine;
      hardware_ID[k + 1]          = jx_CTAlloc(int*, num_numanode_in_machine + 1);
      hardware_ID[k + 1][0]       = jx_CTAlloc(int, 1);
      hardware_ID[k + 1][0][0]    = num_numanode_in_machine;
      for(i = 0; i < num_numanode_in_machine; i++)
      {
         obj_in_numanode       = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, i);  
         num_core_in_numanode  = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj_in_numanode->cpuset, HWLOC_OBJ_CORE);
         num_core_in_total_machine = num_core_in_total_machine + num_core_in_numanode;
         hardware_ID[k + 1][i + 1]    = jx_CTAlloc(int, num_core_in_numanode + 1);
         hardware_ID[k + 1][i + 1][0] = num_core_in_numanode;

         for(j = 0; j < num_core_in_numanode; j++) 
         {
            // wrtite CPU index id
            hardware_ID[k + 1][i + 1][j + 1] = tmp_index + j;
         }
         tmp_index  = tmp_index + hardware_ID[k + 1][i + 1][0];
      }
   }


#if 0
   //print hardware_ID
   if(myid == 0)
   {
      for(k = 0; k < hardware_ID[0][0][0] + 1; k++)
      {
         for(i = 0; i < hardware_ID[k][0][0] + 1; i++)
         {
            for(j = 0; j < hardware_ID[k][i][0] + 1; j++) 
            {
               // wrtite CPU index id
               jx_printf("hardware_ID[%d][%d][%d]  is %d\n", k, i, j, hardware_ID[k][i][j] );
               if(k == 0 || i == 0) break;
            }
            if(k == 0) break;
         }
      }
   }
#endif

   /*
      3. process division and bind
   */
   if(type == 0)
   {
      jx_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);
      
      obj_in_machine  = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, my_machine_id);
      
      cpuset = hwloc_bitmap_dup(obj_in_machine->cpuset);
      if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_PROCESS)) 
      {
         char *str;
         int error = errno;
         hwloc_bitmap_asprintf(&str, obj_in_machine->cpuset);
         printf("Couldn't bind to cpuset %s: %s\n", str, strerror(error));
         free(str);
      }
      hwloc_bitmap_free(cpuset);

      jx_hpNumCoreInCpu(hardware)  = hardware_ID[my_machine_id + 1][1][0];
      jx_hpProcCpuStart(hardware)  = 0;
      jx_hpProcCpuEnd(hardware)	  = 0;   
      jx_hpProcNodeStart(hardware) = 0;	    
      jx_hpProcNodeEnd(hardware)	  = 0;  
   }
   else if(type == 1)
   {
      jx_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);
      jx_hpGetLevelBindId(pid_in_numanode, process_per_numanode, hardware_ID[my_machine_id+1][0][0], &my_numanode_id, &pid_in_core, &process_per_core, &process_statr_in_numanode, &process_end_in_numanode);
   #if 0
      jx_printf("myid is %d, start_node_p_id is %d, end_node_p_id is %d\n", myid, process_statr_in_machine, process_end_in_machine);
      jx_printf("myid is %d, start_cpu_p_id is %d, end_cpu_p_id is %d\n", myid, process_statr_in_numanode, process_end_in_numanode);
   #endif
      for(i = 0; i < my_machine_id; i++)
      {
         my_numanode_id = my_numanode_id + hardware_ID[my_machine_id+1][0][0];
      }

      #if 0
         jx_printf("myid is %d, hardware_ID is %d\n", myid, my_numanode_id);
      #endif

      obj_in_numanode  = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, my_numanode_id);
      
      cpuset = hwloc_bitmap_dup(obj_in_numanode->cpuset);
      if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_PROCESS)) 
      {
         char *str;
         int error = errno;
         hwloc_bitmap_asprintf(&str, obj_in_numanode->cpuset);
         printf("Couldn't bind to cpuset %s: %s\n", str, strerror(error));
         free(str);
      }
      hwloc_bitmap_free(cpuset);

      jx_hpNumCoreInCpu(hardware)  = process_per_numanode;
      jx_hpProcCpuStart(hardware)  = 0;
      jx_hpProcCpuEnd(hardware)	  = 0;   
      jx_hpProcNodeStart(hardware) = global_pid - myid;	    
      jx_hpProcNodeEnd(hardware)	  = global_pid - myid + numprocess - 1; 
   }
   else if(type == 2)
   {  //if(myid == 0)jx_printf("myid is %d, numproces is %d, num_machine is %d\n ", myid, numprocess, num_machine);
      jx_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);
      //if(myid == 0)jx_printf("my_machine_id is %d, pid_in_numanode is %d, process_per_numanode is %d\n ", my_machine_id, pid_in_numanode, process_per_numanode);
      jx_hpGetLevelBindId(pid_in_numanode, process_per_numanode, hardware_ID[my_machine_id+1][0][0], &my_numanode_id, &pid_in_core, &process_per_core, &process_statr_in_numanode, &process_end_in_numanode);
      //if(myid == 0)jx_printf("my_numanode_id is %d, pid_in_core is %d, process_per_core is %d\n ", my_numanode_id, pid_in_core, process_per_core);
      //if(myid == 0)jx_printf("pid_in_core is %d, process_per_core is %d, hardware_ID is %d\n ", pid_in_core, process_per_core, hardware_ID[my_machine_id+1][my_numanode_id+1][0]);
      jx_hpGetLevelBindId(pid_in_core, process_per_core, hardware_ID[my_machine_id+1][my_numanode_id+1][0], &my_core_id, NULL, NULL, NULL, NULL);
      //jx_printf("myid is %d, my_machine_id is %d, my_numanode_id is %d, my_core_id %d\n ",myid, my_machine_id, my_numanode_id, my_core_id);

      /*
      // use another bind way
      obj_in_machine  = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, my_machine_id);
      obj_in_numanode = hwloc_get_obj_inside_cpuset_by_type(topology, obj_in_machine->cpuset, HWLOC_OBJ_NUMANODE, my_numanode_id);
      obj_in_core = hwloc_get_obj_inside_cpuset_by_type(topology, obj_in_numanode->cpuset, HWLOC_OBJ_CORE, my_core_id);
      */
      obj_in_core = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, hardware_ID[my_machine_id+1][my_numanode_id+1][my_core_id+1]);
      
      cpuset = hwloc_bitmap_dup(obj_in_core->cpuset);

      if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_PROCESS)) 
      {
         char *str;
         int error = errno;
         hwloc_bitmap_asprintf(&str, obj_in_core->cpuset);
         printf("Couldn't bind to cpuset %s: %s\n", str, strerror(error));
         free(str);
      }
      hwloc_bitmap_free(cpuset);

      jx_hpNumCoreInCpu(hardware)  = process_per_numanode;
      jx_hpProcNodeStart(hardware) = global_pid - myid;	    
      jx_hpProcNodeEnd(hardware)	  = global_pid - myid + numprocess - 1; 
      jx_hpProcCpuStart(hardware)  = jx_hpProcNodeStart(hardware) + process_statr_in_machine + process_statr_in_numanode;
      //if(global_pid == 0) jx_printf("process_statr_in_machine is %d, process_end_in_numanode is %d\n", process_statr_in_machine, process_end_in_numanode);
      jx_hpProcCpuEnd(hardware)	  = jx_hpProcNodeStart(hardware) + process_statr_in_machine + process_end_in_numanode; 
   }
   
   jx_hpNumCpuInNode(hardware)  = num_numanode_in_total_machine;
   jx_hpNumCoreInNode(hardware) = numprocess;
   #if 0
      if(myid == 0)
      {
         jx_printf("myid is %d, hardware_ID is %d, process_per_numanode is %d\n", myid, hardware_ID[my_machine_id+1][my_numanode_id+1][my_core_id+1], process_per_numanode);
         jx_printf("myid is %d, start_node_p_id is %d, end_node_p_id is %d\n", myid, jx_hpProcNodeStart(hardware), jx_hpProcNodeEnd(hardware));
         jx_printf("myid is %d, start_cpu_p_id is %d, end_cpu_p_id is %d\n", myid, jx_hpProcCpuStart(hardware), jx_hpProcCpuEnd(hardware));
      }
      
   #endif

//createhardware Allreduce Test
#if 1
   /**
    * [0, 1, [2],| 3, 4, [[5]], ||6, 7, [8],| 9, 10, [[11]]]
    * || node level, | cpu level]
    * [] cpu represent
    * [[]] node represent
    */

   // definition
   MPI_Comm CORE_COMM, Intra_CPU_COMM, NODE_COMM; //CPU_COMM;

   //init
   CORE_COMM =  MPI_COMM_NULL;
   Intra_CPU_COMM  =  MPI_COMM_NULL;
   NODE_COMM =  MPI_COMM_NULL;
   // CPU_COMM  =  MPI_COMM_NULL;

   // 同一Core内的Comm
   jx_MPI_Comm_split(shmcomm, jx_hpProcCpuStart(hardware), global_pid, &CORE_COMM);
   
   // 结点内CPU代表Comm
   // JX_Int cpu_color = global_pid % process_per_core;
   JX_Int Intra_cpu_color = 1;
   if(pid_in_core != process_per_core - 1)
   {
      Intra_cpu_color = MPI_UNDEFINED;
   }
   jx_MPI_Comm_split(shmcomm, Intra_cpu_color, global_pid, &Intra_CPU_COMM);

   // 所有CPU代表Comm(暂未使用到)
   // JX_Int cpu_color = global_pid % process_per_core;
   // JX_Int cpu_color = 2;
   // if(pid_in_core != process_per_core - 1)
   // {
   //    cpu_color = MPI_UNDEFINED;
   // }
   // jx_MPI_Comm_split(comm, cpu_color, global_pid, &CPU_COMM);


   // Node代表Comm
   // JX_Int node_color = global_pid % process_per_numanode;
   JX_Int node_color = 3;
   if(pid_in_numanode != process_per_numanode - 1)
   {
      node_color = MPI_UNDEFINED;
   }
   jx_MPI_Comm_split(comm, node_color, global_pid, &NODE_COMM);

   jx_hpCoreComm(hardware)      = CORE_COMM;
   jx_hpIntraCpuComm(hardware)  = Intra_CPU_COMM;
   jx_hpNodeComm(hardware)      = NODE_COMM;
   // jx_hpCpuComm(hardware)	= CPU_COMM;
   // jx_hpIntraNodeComm(hardware) = shmcomm;

   #if 0
      JX_Int  mycoreid=-1, myintracpuid=-1, mynodeid=-1, mycpuid=-1;
      if(CORE_COMM != MPI_COMM_NULL) jx_MPI_Comm_rank(CORE_COMM, &mycoreid);
      if(CPU_COMM != MPI_COMM_NULL) jx_MPI_Comm_rank(CPU_COMM, &mycpuid);
      if(Intra_CPU_COMM != MPI_COMM_NULL) jx_MPI_Comm_rank(Intra_CPU_COMM, &myintracpuid);
      if(NODE_COMM != MPI_COMM_NULL) jx_MPI_Comm_rank(NODE_COMM, &mynodeid);
      jx_printf("myid is %d, mycoreid is %d, mycpuid is %d, myintracpuid is %d, mynodeid is %d\n",
         global_pid, mycoreid, mycpuid, myintracpuid, mynodeid);
      jx_MPI_Barrier(comm);
      // exit(0);
   #endif

#endif

   hwloc_topology_destroy(topology); 
   jx_MPI_Comm_free(&shmcomm);

   for(k = 0; k < num_machine + 1; k++)
   {
      num_numanode_in_machine = hardware_ID[k][0][0];
      for(i = 0; i < num_numanode_in_machine + 1; i++)
      {
         jx_TFree(hardware_ID[k][i]);
         if(k == 0) break;
      }
      jx_TFree(hardware_ID[k]);
   }
   jx_TFree(hardware_ID);

   return EXIT_SUCCESS;
}
#endif
/*!
 * \fn jx_hpCSRMatrixSetCpuLevel
 * \brief Get col_start and col_end of offd_cpu 
 * \author mrz 
 * \date 2021/10/15
 */
JX_Int jx_hpCSRMatrixSetCpuLevel(jx_HardwareInfo *hardware)
{
   JX_Int my_id, nprocs;
   MPI_Comm comm = jx_hpComm(hardware);
   JX_Int num_core_in_cpu = jx_hpNumCoreInCpu(hardware);
   JX_Int num_core_in_node = jx_hpNumCoreInNode(hardware);
   JX_Int bind_type = jx_hpBindType(hardware);

   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &nprocs);

   if(bind_type == 0 || bind_type == 1) //绑定到CPU层或node层
   {
      jx_hpProcCpuStart(hardware) = 0;
      jx_hpProcCpuEnd(hardware) = 0;
   }
   else if(num_core_in_cpu == 1)
   {
      jx_hpProcCpuStart(hardware) = my_id;
      jx_hpProcCpuEnd(hardware) = my_id;
   }
   else
   {
      MPI_Comm Cpu_Comm = MPI_COMM_NULL;
      JX_Int core_id_start_in_node = my_id / num_core_in_node * num_core_in_node;
      JX_Int core_id_in_node = my_id % num_core_in_node;
      JX_Int core_id_start_in_cpu = core_id_in_node / num_core_in_cpu * num_core_in_cpu;
      JX_Int process_id_start_in_cpu = core_id_start_in_node + core_id_start_in_cpu; 

      JX_Int end_core_id_start_in_cpu = nprocs / num_core_in_cpu * num_core_in_cpu;
      JX_Int r_cpu = nprocs%num_core_in_cpu;
      if(r_cpu != 0 && my_id >= end_core_id_start_in_cpu)  
      {
         jx_hpProcCpuStart(hardware) = process_id_start_in_cpu;
         jx_hpProcCpuEnd(hardware) = process_id_start_in_cpu + r_cpu - 1;
      }
      else
      {
         jx_hpProcCpuStart(hardware) = process_id_start_in_cpu;
         jx_hpProcCpuEnd(hardware) = process_id_start_in_cpu + num_core_in_cpu - 1;
      }
#if 0
      jx_MPI_Comm_split(comm, process_id_start_in_cpu, my_id, &Cpu_Comm);
#endif
      jx_hpIntraCpuComm(hardware) = Cpu_Comm;  
   }

   return jx_error_flag;
}


/*!
 * \fn jx_hpCSRMatrixGetNodeCol
 * \brief Get col_start and col_end of offd_node 
 * \author mrz 
 * \date 2021/08/30
 */
JX_Int jx_hpCSRMatrixSetNodeLevel(jx_HardwareInfo *hardware)
{
   JX_Int my_id,nprocs;
   MPI_Comm comm = jx_hpComm(hardware);   
   JX_Int num_cpu_in_node = jx_hpNumCpuInNode(hardware);
   JX_Int num_core_in_node = jx_hpNumCoreInNode(hardware);
   JX_Int bind_type = jx_hpBindType(hardware);

   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &nprocs);

   if( bind_type == 0 ) //绑定到node层
   {
      jx_hpProcNodeStart(hardware) = 0;
      jx_hpProcNodeEnd(hardware)   = 0;
   }
   else if(num_cpu_in_node == 1)
   {
      jx_hpProcNodeStart(hardware) = jx_hpProcCpuStart(hardware);
      jx_hpProcNodeEnd(hardware)   = jx_hpProcCpuEnd(hardware);
      jx_hpNodeComm(hardware)      = jx_hpIntraCpuComm(hardware);
   }
   else
   {
      if(bind_type == 2)
      {
         MPI_Comm Node_Comm = MPI_COMM_NULL;
         JX_Int process_id_start_in_node = my_id / num_core_in_node * num_core_in_node;

         JX_Int end_core_id_start_in_node = nprocs / num_core_in_node * num_core_in_node;
         JX_Int r_node = nprocs%num_core_in_node;
         if(r_node != 0 && my_id >= end_core_id_start_in_node)
         {
            jx_hpProcNodeStart(hardware) = process_id_start_in_node;
            jx_hpProcNodeEnd(hardware) = process_id_start_in_node + r_node - 1;
         }
         else
         {
            jx_hpProcNodeStart(hardware) = process_id_start_in_node;
            jx_hpProcNodeEnd(hardware) = process_id_start_in_node + num_core_in_node - 1;
         }
#if 0
         jx_MPI_Comm_split(comm, process_id_start_in_node, my_id, &Node_Comm);
#endif
         jx_hpNodeComm(hardware) = Node_Comm;
      }
      else
      {
         MPI_Comm Node_Comm = MPI_COMM_NULL;
         JX_Int process_id_start_in_node = my_id / num_cpu_in_node * num_cpu_in_node;

         JX_Int end_cpu_id_start_in_node = nprocs / num_cpu_in_node * num_cpu_in_node;
         JX_Int r_node = nprocs%num_cpu_in_node;
         if(r_node != 0 && my_id >= end_cpu_id_start_in_node)
         {
            jx_hpProcNodeStart(hardware) = process_id_start_in_node;
            jx_hpProcNodeEnd(hardware) = process_id_start_in_node + r_node - 1;
         }
         else
         {
            jx_hpProcNodeStart(hardware) = process_id_start_in_node;
            jx_hpProcNodeEnd(hardware) = process_id_start_in_node + num_cpu_in_node - 1;
         }
#if 0
         jx_MPI_Comm_split(comm, process_id_start_in_node, my_id, &Node_Comm);
#endif
         jx_hpNodeComm(hardware) = Node_Comm;
      }
   }
   return jx_error_flag;
}

/*!
 * \fn jx_hpBuildMatParFromOneFile
 * \brief read matrix from files and create hierarchy diag block. 
 * \date 2021/10/27
 * \edit 2022/3/7
 * \author mrz
 */
jx_hpCSRMatrix *
jx_hpBuildMatParFromOneFile( char *filename, JX_Int num_functions, JX_Int file_base)
{
   jx_hpCSRMatrix *A = jx_hpInithpCSRMatrix();
   jx_hpCSRMatrixPar(A)  = jx_BuildMatParFromOneFile(filename, num_functions, file_base);
#ifdef USING_HWLOC
   A = jx_hpCreateMatrixLevelBlock(A);
#endif
   return (A);
}

/*!
 * \fn jx_hpCSRMatrixCreate
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrix * jx_hpCSRMatrixCreate(MPI_Comm   comm,
                       JX_Int        global_num_rows,
                       JX_Int        global_num_cols,
                       JX_Int       *row_starts,
                       JX_Int       *col_starts,
                       JX_Int        num_cols_offd,
                       JX_Int        num_nonzeros_diag,
                       JX_Int        num_nonzeros_offd)
{
   jx_hpCSRMatrix *matrix = jx_hpInithpCSRMatrix();
   jx_hpCSRMatrixPar(matrix) = jx_ParCSRMatrixCreate(comm, global_num_rows, global_num_cols,
      row_starts, col_starts, num_cols_offd, num_nonzeros_diag, num_nonzeros_offd);
#ifdef USING_HWLOC
   matrix = jx_hpCreateMatrixLevelBlock(matrix);
#endif
   return matrix;
}

/*!
 * \fn jx_hpCSRMatrixInitialize
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
JX_Int jx_hpCSRMatrixInitialize(jx_hpCSRMatrix  *A)
{
   if(!jx_hpCSRMatrixPar(A))
   {
      jx_error_in_arg(1);
      return jx_error_flag;      
   }
   jx_ParCSRMatrixInitialize(jx_hpCSRMatrixPar(A));
   return jx_error_flag;
}

/*!
 * \fn jx_hpInithpCSRMatrix
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrix  *jx_hpInithpCSRMatrix()
{
   jx_hpCSRMatrix  *A = jx_CTAlloc(jx_hpCSRMatrix, 1);
   jx_hpCSRMatrixPar(A) = NULL;
   jx_hpCSRMatrixCpu(A) = jx_hpInitMatrixLevel();
   jx_hpCSRMatrixNode(A) = jx_hpInitMatrixLevel();
   jx_hpCSRMatrixCore(A) = jx_hpInitMatrixLevel();
   return A;
}

/*!
 * \fn jx_hpInitMatrixLevel
 * \brief Init Matrix Level.
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrixLevel *jx_hpInitMatrixLevel()
{
   jx_hpCSRMatrixLevel *level_matrix = jx_CTAlloc(jx_hpCSRMatrixLevel, 1);
   jx_hpCSRMatrixOffdlevel(level_matrix) = NULL;
   jx_hpCSRMatrixColMapOffdlevel(level_matrix) = NULL;
   jx_hpCSRMatrixlevelCommPkg(level_matrix) = NULL;
   return level_matrix;
}

/*!
 * \fn jx_createMatrixLevel
 * \brief create block of each Matrix Level.
 * \date 2022/3/7
 * \author mrz
 */
jx_hpCSRMatrix *jx_hpCreateMatrixLevelBlock(jx_hpCSRMatrix *A)
{
   if(!jx_hpCSRMatrixCpu(A))
   jx_hpCSRMatrixCpu(A)  = jx_hpCreateMatrixCpu(A);
   if(!jx_hpCSRMatrixNode(A))
   jx_hpCSRMatrixNode(A) = jx_hpCreateMatrixNode(A);
   if(!jx_hpCSRMatrixCore(A))
   jx_hpCSRMatrixCore(A) = jx_hpCreateMatrixCore(A);
   return A;
}

/*!
 * \fn jx_hpCreateMatrixCpu
 * \brief create cpu level of matrix 
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrixLevel *jx_hpCreateMatrixCore(jx_hpCSRMatrix *A )
{
   JX_Int i;
   JX_Int  num_rows = jx_hpCSRMatrixNumRows(A);
   jx_hpCSRMatrixLevel  *core_matrix  = jx_CTAlloc(jx_hpCSRMatrixLevel, 1);
   jx_CSRMatrix         *offd_core    = jx_CTAlloc(jx_CSRMatrix, 1);
   JX_Int               *offd_core_i  = jx_CTAlloc(JX_Int, num_rows + 1);

   for (i = 0; i < num_rows + 1; i ++)
   {
      offd_core_i[i] = 0;
   }
   jx_CSRMatrixNumCols(offd_core) = 0;
   jx_CSRMatrixI(offd_core) = offd_core_i; 
   jx_hpCSRMatrixOffdlevel(core_matrix) = offd_core;
   return core_matrix;
}

/*!
 * \fn jx_hpCreateMatrixCpu
 * \brief create cpu level of matrix 
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrixLevel *jx_hpCreateMatrixCpu(jx_hpCSRMatrix *A )
{
   if(hp_hardware)
   {
      //jx_hpCSRMatrixLevel *cpu_matrix = NULL;
      //JX_Int bind_type = jx_hpBindType(hardware);
      //JX_Int num_core_in_cpu = jx_hpNumCoreInCpu(hardware);
      

      jx_hpCSRMatrixLevel *cpu_matrix = jx_hpCSRMatrixGenerateOffdCpu(A, hp_hardware);
      return cpu_matrix;
   }
   else
   {
      jx_printf("erro: CreateMatrixCpu, hardwareinfo is empty\n");
      return NULL;
   }
}

/*!
 * \fn JX_Int jx_hpCreateMatrixNode
 * \brief create node level of matrix  
 * \date 2021/11/13
 * \author mrz
 */
jx_hpCSRMatrixLevel *jx_hpCreateMatrixNode(jx_hpCSRMatrix *A )
{
   if(hp_hardware)
   {
      //JX_Int num_cpu_in_node = jx_hpNumCpuInNode(hardware);
      //JX_Int num_core_in_cpu = jx_hpNumCoreInCpu(hardware);
   
      //JX_Int bind_type = jx_hpBindType(hardware);
      jx_hpCSRMatrixLevel *node_matrix = jx_hpCSRMatrixGenerateOffdNode(A, hp_hardware);
      return node_matrix;
   }
   else
   {
      jx_printf("erro: CreateMatrixNode, hardwareinfo is empty\n");
      return NULL;
   }
}

/*!
 * \fn jx_hpCSRMatrixGenerateOffdCpu
 * \brief Generate cpu level matrix 
 * \author mrz 
 * \date 2021/08/31
 */
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdCpu(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo * hardware)
{
   //test need
   // MPI_Comm comm = jx_hpCSRMatrixComm(hp_matrix);
   // JX_Int myid;
   // jx_MPI_Comm_rank(comm, &myid);

   jx_hpCSRMatrixLevel *cpu_matrix = jx_hpInitMatrixLevel();

   JX_Int  i, j, j_cpu, col_index, col_map_index;
   jx_CSRMatrix *offd = jx_hpCSRMatrixOffd(hp_matrix);
   JX_Int  num_rows = jx_CSRMatrixNumRows(offd);
   //JX_Int  num_cols = jx_CSRMatrixNumCols(offd);
   JX_Int  num_cols = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   JX_Real *offd_data = jx_CSRMatrixData(offd);
   JX_Int *offd_i = jx_CSRMatrixI(offd);
   JX_Int *offd_j = jx_CSRMatrixJ(offd);
   JX_Int *offd_col_map = jx_hpCSRMatrixColMapOffd(hp_matrix);

   jx_CSRMatrix *offd_cpu;
   JX_Int *col_map_offd_cpu  = NULL;
   JX_Real *offd_cpu_data;
   JX_Int  *offd_cpu_i ;
   JX_Int  *offd_cpu_j;

   JX_Int  *marker_cpu;
   JX_Int num_cols_offd_cpu = 0;
   JX_Int first_elmt = offd_i[0];
   JX_Int num_nonzeros_cpu;
   JX_Int counter_cpu;

   JX_Int *row_start = jx_hpCSRMatrixRowStarts(hp_matrix);
   JX_Int cpu_colstartnum = row_start[jx_hpProcCpuStart(hardware)];
   JX_Int cpu_colendnum = row_start[jx_hpProcCpuEnd(hardware) + 1] - 1;


   //generate offd_cpu and col_cpu_map
   offd_cpu_i = jx_CTAlloc(JX_Int, num_rows + 1);
   offd_cpu = jx_CTAlloc(jx_CSRMatrix, 1);
   marker_cpu = jx_CTAlloc(JX_Int, num_cols);

   for (i = 0; i < num_cols; i ++)
   {
      marker_cpu[i] = 0;
   }
   
   num_nonzeros_cpu = 0;

   for (i = 0; i < num_rows; i ++)
   {
      offd_cpu_i[i] = num_nonzeros_cpu;

      for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
      {
         // col_index = offd_j[j]; offd_col_map[col_index]
         col_map_index = offd_col_map[offd_j[j]];
         if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
         {
            if (!marker_cpu[col_map_index])
            {
               marker_cpu[col_map_index] = 1;
               num_cols_offd_cpu ++;
            }
               num_nonzeros_cpu ++;
         }
      }
   }

   if(num_nonzeros_cpu == 0)
   {
      for (i = 0; i < num_rows + 1; i ++)
      {
         offd_cpu_i[i] = 0;
      }
      jx_CSRMatrixNumCols(offd_cpu) = 0;
      jx_CSRMatrixI(offd_cpu) = offd_cpu_i;
      jx_hpCSRMatrixOffdlevel(cpu_matrix) = offd_cpu;
   }
   else
   {
      offd_cpu_i[num_rows] = num_nonzeros_cpu;
      col_map_offd_cpu = jx_CTAlloc(JX_Int, num_cols_offd_cpu);

      counter_cpu = 0;
      for (i = 0; i < num_cols; i ++)
      {    
         if (marker_cpu[i])
         {
            col_map_offd_cpu[counter_cpu] = i;
            marker_cpu[i] = counter_cpu;
            counter_cpu ++;
         }              
      }
      jx_CSRMatrixNumRows(offd_cpu) = num_rows;
      jx_CSRMatrixNumNonzeros(offd_cpu) = num_nonzeros_cpu;
      jx_CSRMatrixNumCols(offd_cpu) = num_cols_offd_cpu;
      jx_CSRMatrixJ(offd_cpu) = NULL;
      jx_CSRMatrixData(offd_cpu) = NULL;
      jx_CSRMatrixOwnsData(offd_cpu) = 1;
      jx_CSRMatrixNumRownnz(offd_cpu) = num_rows;
      jx_CSRMatrixI(offd_cpu) = offd_cpu_i;

      jx_CSRMatrixInitialize(offd_cpu);

      offd_cpu_data = jx_CSRMatrixData(offd_cpu);
      offd_cpu_j = jx_CSRMatrixJ(offd_cpu);
      
      j_cpu = 0;
      for (i = 0; i < num_rows; i ++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
         {
            col_index = offd_j[j];
            col_map_index = offd_col_map[col_index];
            if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
            {
               offd_cpu_data[j_cpu] = offd_data[j];
               offd_cpu_j[j_cpu++] = marker_cpu[col_map_index];
            }
         }    
      }
      jx_CSRMatrixOwnsData(offd_cpu) = 1;
      jx_hpCSRMatrixOffdlevel(cpu_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(cpu_matrix) = offd_cpu;
      jx_hpCSRMatrixColMapOffdlevel(cpu_matrix) = jx_CTAlloc(JX_Int, num_cols_offd_cpu);
      jx_hpCSRMatrixColMapOffdlevel(cpu_matrix) = col_map_offd_cpu;   
   }
   
   jx_TFree(marker_cpu);
   return cpu_matrix;
}

/*!
 * \fn jx_hpCSRMatrixGenerateOffdNode
 * \brief Generate node level matrix 
 * \author mrz 
 * \date 2021/08/31
 */
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdNode(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware)
{
      //test need
   // MPI_Comm comm = jx_hpCSRMatrixComm(hp_matrix);
   // JX_Int myid;
   // jx_MPI_Comm_rank(comm, &myid);

   jx_hpCSRMatrixLevel *node_matrix = jx_hpInitMatrixLevel();
   JX_Int  i, j, j_node, col_index, col_map_index;
   jx_CSRMatrix *offd = jx_hpCSRMatrixOffd(hp_matrix);
   JX_Int  num_rows = jx_CSRMatrixNumRows(offd);
   //JX_Int  num_cols = jx_CSRMatrixNumCols(offd);
   JX_Int  num_cols = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   JX_Real *offd_data = jx_CSRMatrixData(offd);
   JX_Int *offd_i = jx_CSRMatrixI(offd);
   JX_Int *offd_j = jx_CSRMatrixJ(offd);
   JX_Int *offd_col_map = jx_hpCSRMatrixColMapOffd(hp_matrix);

   jx_CSRMatrix *offd_node;
   JX_Int  *col_map_offd_node = NULL;
   JX_Real *offd_node_data;
   JX_Int  *offd_node_i;
   JX_Int  *offd_node_j;   

   JX_Int *marker_node;
   JX_Int num_cols_offd_node = 0;
   JX_Int first_elmt = offd_i[0];
   JX_Int num_nonzeros_node;
   JX_Int counter_node;

   JX_Int *row_start = jx_hpCSRMatrixRowStarts(hp_matrix);
   JX_Int  node_colstartnum = row_start[jx_hpProcNodeStart(hardware)];
   JX_Int  node_colendnum   = row_start[jx_hpProcNodeEnd(hardware)+ 1] - 1;
   
   // generate offd
   offd_node_i = jx_CTAlloc(JX_Int, num_rows + 1);
   offd_node = jx_CTAlloc(jx_CSRMatrix, 1);
   marker_node = jx_CTAlloc(JX_Int, num_cols);

   for (i = 0; i < num_cols; i ++)
   {
      marker_node[i] = 0;
   }     

   num_nonzeros_node = 0;

   for (i = 0; i < num_rows; i ++)
   {
      offd_node_i[i] = num_nonzeros_node;
   
      for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
      {
         //col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
         // {
         //    continue;
         // }
         // else 
         if(col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
         {
               // if (!marker_node[col_index])
               // {
               //    marker_node[col_index] = 1;
               //    num_cols_offd_node ++;
               // }
               if (!marker_node[col_map_index])
               {
                  marker_node[col_map_index] = 1;
                  num_cols_offd_node ++;
               }
                  num_nonzeros_node ++;
         }
      }
   }

   if(num_nonzeros_node == 0)
   {
      for (i = 0; i < num_rows + 1; i ++)
      {
         offd_node_i[i] = 0;
      }
      jx_CSRMatrixNumCols(offd_node) = 0;
      jx_CSRMatrixI(offd_node) = offd_node_i;
      jx_hpCSRMatrixOffdlevel(node_matrix) = offd_node;
   }
   else
   {
      offd_node_i[num_rows] = num_nonzeros_node;
      col_map_offd_node = jx_CTAlloc(JX_Int, num_cols_offd_node);

      counter_node = 0;
      for (i = 0; i < num_cols; i ++)
      {
         if (marker_node[i])
         {
            col_map_offd_node[counter_node] = i;
            marker_node[i] = counter_node;
            counter_node ++;
         }                
      }

      jx_CSRMatrixNumRows(offd_node) = num_rows;
      jx_CSRMatrixNumNonzeros(offd_node) = num_nonzeros_node;
      jx_CSRMatrixNumCols(offd_node) = num_cols_offd_node;
      jx_CSRMatrixJ(offd_node) = NULL;
      jx_CSRMatrixData(offd_node) = NULL;
      jx_CSRMatrixOwnsData(offd_node) = 1;
      jx_CSRMatrixNumRownnz(offd_node) = num_rows;
      jx_CSRMatrixI(offd_node) = offd_node_i;

      jx_CSRMatrixInitialize(offd_node);

      offd_node_data = jx_CSRMatrixData(offd_node);
      offd_node_j = jx_CSRMatrixJ(offd_node);

      j_node = 0;
      for (i = 0; i < num_rows; i ++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
         {
            col_index = offd_j[j];
            col_map_index = offd_col_map[col_index];
            // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
            // {
            //    continue;
            // }            
            // else 
            if (col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
            {
               offd_node_data[j_node] = offd_data[j];
               offd_node_j[j_node++] = marker_node[col_map_index];
            }
         }    
      }
      jx_CSRMatrixOwnsData(offd_node) = 1;
      jx_hpCSRMatrixOffdlevel(node_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(node_matrix) = offd_node;
      jx_hpCSRMatrixColMapOffdlevel(node_matrix) = jx_CTAlloc(JX_Int, num_cols_offd_node);
      jx_hpCSRMatrixColMapOffdlevel(node_matrix) = col_map_offd_node;
   }
   jx_TFree(marker_node); 
   return node_matrix;
}

/*!
 * \fn JX_Int jx_hpMatrixLevelToPar
 * \date 2022/03/31
 */
jx_ParCSRMatrix * 
jx_hpMatrixLevelToPar( jx_hpCSRMatrix *hp_matrix, JX_Int Level)
{
   jx_ParCSRMatrix *par_matrix = jx_CTAlloc(jx_ParCSRMatrix, 1);
   jx_ParCSRMatrixComm(par_matrix) = jx_hpCSRMatrixComm(hp_matrix);
   jx_ParCSRMatrixGlobalNumRows(par_matrix) = jx_hpCSRMatrixGlobalNumRows(hp_matrix);
   jx_ParCSRMatrixGlobalNumCols(par_matrix)  = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   jx_ParCSRMatrixFirstRowIndex(par_matrix) = jx_hpCSRMatrixFirstRowIndex(hp_matrix);
   jx_ParCSRMatrixFirstColDiag(par_matrix) = jx_hpCSRMatrixFirstColDiag(hp_matrix);
   jx_ParCSRMatrixLastRowIndex(par_matrix) = jx_hpCSRMatrixLastRowIndex(hp_matrix);
   jx_ParCSRMatrixLastColDiag(par_matrix) = jx_hpCSRMatrixLastRowIndex(hp_matrix);
   jx_ParCSRMatrixDiag(par_matrix) = jx_hpCSRMatrixDiag(hp_matrix);
   if(Level == 0)//core
   {
      jx_ParCSRMatrixOffd(par_matrix) = jx_hpCSRMatrixOffdCore(hp_matrix);
      jx_ParCSRMatrixColMapOffd(par_matrix) = jx_hpCSRMatrixColMapOffdCore(hp_matrix);
      jx_ParCSRMatrixCommPkg(par_matrix) = jx_hpCSRMatrixCoreCommPkg(hp_matrix);
      jx_ParCSRMatrixNumNonzeros(par_matrix) = jx_hpCSRMatrixDNumNonzeros(hp_matrix);    
   }
   else if(Level == 1)//cpu
   {
      jx_ParCSRMatrixOffd(par_matrix) = jx_hpCSRMatrixOffdCpu(hp_matrix);
      jx_ParCSRMatrixColMapOffd(par_matrix) = jx_hpCSRMatrixColMapOffdCpu(hp_matrix);
      jx_ParCSRMatrixCommPkg(par_matrix) = jx_hpCSRMatrixCpuCommPkg(hp_matrix);
      jx_ParCSRMatrixNumNonzeros(par_matrix) = jx_hpCSRMatrixDNumNonzeros(hp_matrix) + jx_CSRMatrixNumNonzeros(jx_hpCSRMatrixOffdCpu(hp_matrix));
   }
   else if(Level == 2)//node
   {
      jx_ParCSRMatrixOffd(par_matrix) = jx_hpCSRMatrixOffdNode(hp_matrix);
      jx_ParCSRMatrixColMapOffd(par_matrix) = jx_hpCSRMatrixColMapOffdNode(hp_matrix);
      jx_ParCSRMatrixCommPkg(par_matrix) = jx_hpCSRMatrixNodeCommPkg(hp_matrix);
      jx_ParCSRMatrixNumNonzeros(par_matrix) = jx_hpCSRMatrixDNumNonzeros(hp_matrix) + jx_CSRMatrixNumNonzeros(jx_hpCSRMatrixOffdNode(hp_matrix));   
   }
   else
   {
      jx_printf("Level must be 1 or 2\n");
      return NULL;
   }
   jx_ParCSRMatrixDiagT(par_matrix) = jx_hpCSRMatrixDiagT(hp_matrix);
   jx_ParCSRMatrixOffdT(par_matrix) = NULL;//暂时先赋值为空
   jx_ParCSRMatrixRowStarts(par_matrix) = jx_hpCSRMatrixRowStarts(hp_matrix);//SPMV中不需要RowStarts
   jx_ParCSRMatrixColStarts(par_matrix) = jx_hpCSRMatrixColStarts(hp_matrix);
   jx_ParCSRMatrixCommPkgT(par_matrix) = NULL;//暂时先赋值为空
   jx_ParCSRMatrixOwnsData(par_matrix) = jx_hpCSRMatrixOwnsData(hp_matrix);
   jx_ParCSRMatrixOwnsRowStarts(par_matrix) = jx_hpCSRMatrixOwnsRowStarts(hp_matrix);
   jx_ParCSRMatrixOwnsColStarts(par_matrix) = jx_hpCSRMatrixOwnsColStarts(hp_matrix);
   jx_ParCSRMatrixDNumNonzeros(par_matrix) = jx_hpCSRMatrixDNumNonzeros(hp_matrix);
   jx_ParCSRMatrixRowindices(par_matrix)  = jx_hpCSRMatrixRowindices(hp_matrix);
   jx_ParCSRMatrixRowvalues(par_matrix) = jx_hpCSRMatrixRowvalues(hp_matrix);
   jx_ParCSRMatrixGetrowactive(par_matrix)  = jx_hpCSRMatrixGetrowactive(hp_matrix);
   jx_ParCSRMatrixAssumedPartition(par_matrix) = jx_hpCSRMatrixAssumedPartition(hp_matrix);
   return par_matrix;
}


/*!
 * \fn JX_Int jx_ParCSRMatrixSetRowStartsOwner
 * \date 2009/07/10
 */
JX_Int 
jx_hpCSRMatrixSetRowStartsOwner( jx_hpCSRMatrix *matrix, JX_Int owns_row_starts )
{
   if (!matrix)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }

   jx_hpCSRMatrixOwnsRowStarts(matrix) = owns_row_starts;

   return jx_error_flag;
}


/*!
 * \fn JX_Int jx_hpCSRMatrixDestroyAssumedPartition
 * \date 2021/09/08
 */
JX_Int 
jx_hpCSRMatrixDestroyAssumedPartition( jx_hpCSRMatrix *matrix )
{
   return jx_ParCSRMatrixDestroyAssumedPartition(jx_hpCSRMatrixPar(matrix));
}

/*!
 * \fn JX_Int jx_hpCSRMatrixDestroy
 * \brief Destroy a hpCSR matrix.
 * \date 2021/09/08
 */
JX_Int 
jx_hpCSRMatrixDestroy(jx_hpCSRMatrix *matrix)
{
   if (matrix)
   {
      if(jx_hpCSRMatrixPar(matrix))
      {
         jx_ParCSRMatrixDestroy(jx_hpCSRMatrixPar(matrix));
      }
      if(jx_hpCSRMatrixNode(matrix))
      {
         jx_hpCSRMatrixLevelDestroy(jx_hpCSRMatrixNode(matrix));
      }      
      if(jx_hpCSRMatrixCpu(matrix))
      {
         jx_hpCSRMatrixLevelDestroy(jx_hpCSRMatrixCpu(matrix));
      }
      if(jx_hpCSRMatrixCore(matrix))
      {
         jx_hpCSRMatrixLevelDestroy(jx_hpCSRMatrixCore(matrix));
      }
      jx_TFree(matrix);
   }

   return jx_error_flag;
}

JX_Int 
jx_hpCSRhardwareDestroy()
{
   if(hp_hardware)
   {

      if(jx_hpCoreComm(hp_hardware)!=MPI_COMM_NULL)
      {
         jx_MPI_Comm_free(&jx_hpCoreComm(hp_hardware));
         jx_hpCoreComm(hp_hardware) = MPI_COMM_NULL;
      }

      if(jx_hpIntraCpuComm(hp_hardware)!=MPI_COMM_NULL)
      {
         jx_MPI_Comm_free(&jx_hpIntraCpuComm(hp_hardware));
         jx_hpIntraCpuComm(hp_hardware) = MPI_COMM_NULL;
      }

      if(jx_hpNodeComm(hp_hardware)!=MPI_COMM_NULL)
      {
         jx_MPI_Comm_free(&jx_hpNodeComm(hp_hardware));
         jx_hpNodeComm(hp_hardware) = MPI_COMM_NULL;
      }
      jx_TFree(hp_hardware);
   }

   return jx_error_flag;
}

JX_Int 
jx_hpCSRMatrixLevelDestroy(jx_hpCSRMatrixLevel *matrix)
{
   if (matrix)
   {
      if(jx_hpCSRMatrixOffdlevel(matrix))
      {
         jx_CSRMatrixDestroy(jx_hpCSRMatrixOffdlevel(matrix));
      }
      if(jx_hpCSRMatrixColMapOffdlevel(matrix))
      {
         jx_TFree(jx_hpCSRMatrixColMapOffdlevel(matrix));
      }
      if(jx_hpCSRMatrixlevelCommPkg(matrix))
      {
         jx_MatvecCommPkgDestroy(jx_hpCSRMatrixlevelCommPkg(matrix));
      }
      jx_TFree(matrix);
   }

   return jx_error_flag;
}

JX_Int
jx_hpCSRMatrixPrint( jx_hpCSRMatrix *hp_matrix, const char *file_name )
{
   MPI_Comm comm;
   JX_Int  global_num_rows;
   JX_Int  global_num_cols;
   JX_Int  num_cols_cpu = 0;
   JX_Int  num_cols_node = 0;
   JX_Int *col_map_offd;
#ifndef JX_NO_GLOBAL_PARTITION
   JX_Int *row_starts;
   JX_Int *col_starts;
#endif
   JX_Int   my_id, i, num_procs;
   char  new_file_offd_cpu[80], new_file_offd_node[80],new_file_o[80], new_file_info[80];
   FILE *fp;
   JX_Int num_cols_offd = 0;
#ifdef JX_NO_GLOBAL_PARTITION
   JX_Int row_s, row_e, col_s, col_e;
#endif
   if (!hp_matrix)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }

   comm = jx_hpCSRMatrixComm(hp_matrix);
   global_num_rows = jx_hpCSRMatrixGlobalNumRows(hp_matrix);
   global_num_cols = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   col_map_offd = jx_hpCSRMatrixColMapOffd(hp_matrix);
   

#ifndef JX_NO_GLOBAL_PARTITION
   row_starts = jx_hpCSRMatrixRowStarts(hp_matrix);
   col_starts = jx_hpCSRMatrixColStarts(hp_matrix);
#endif

   if (jx_hpCSRMatrixOffd(hp_matrix))
   {
      num_cols_offd = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffd(hp_matrix));
   }
   if(jx_hpCSRMatrixOffdCpu(hp_matrix))
   {
      num_cols_cpu = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdCpu(hp_matrix));
   }
   if(jx_hpCSRMatrixOffdNode(hp_matrix))
   {
      num_cols_node = jx_CSRMatrixNumCols(jx_hpCSRMatrixOffdCpu(hp_matrix));
   }


   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &num_procs);


   jx_sprintf(new_file_offd_cpu, "%s.offd_cpu.%d", file_name, my_id);
   jx_sprintf(new_file_offd_node, "%s.offd_node.%d", file_name, my_id);
   jx_sprintf(new_file_o, "%s.O.%d", file_name, my_id);
   jx_sprintf(new_file_info, "%s.INFO.%d", file_name, my_id);
   if(num_cols_cpu != 0) jx_CSRMatrixPrint(jx_hpCSRMatrixOffdCpu(hp_matrix), new_file_offd_cpu);
   if(num_cols_node != 0) jx_CSRMatrixPrint(jx_hpCSRMatrixOffdNode(hp_matrix), new_file_offd_node);   
   if(num_cols_offd != 0)
   {
      jx_CSRMatrixPrint(jx_hpCSRMatrixOffd(hp_matrix),new_file_o);
   }
   
   fp = fopen(new_file_info, "w");
   jx_fprintf(fp, "%d\n", global_num_rows);
   jx_fprintf(fp, "%d\n", global_num_cols);
   jx_fprintf(fp, "%d\n", num_cols_offd);
   
#ifdef JX_NO_GLOBAL_PARTITION
   row_s = jx_hpCSRMatrixFirstRowIndex(hp_matrix);
   row_e = jx_hpCSRMatrixLastRowIndex(hp_matrix);
   col_s =  jx_hpCSRMatrixFirstColDiag(hp_matrix);
   col_e =  jx_hpCSRMatrixLastColDiag(hp_matrix);
   /* add 1 to the ends because this is a starts partition */
   jx_fprintf(fp, "%d %d %d %d\n", row_s, row_e + 1, col_s, col_e + 1);
#else
   jx_fprintf(fp, "row_start col_start: \n");
   for (i = 0; i < num_procs; i ++)
   {
      jx_fprintf(fp, "%d %d\n", row_starts[i], col_starts[i]);
   }
#endif

   jx_fprintf(fp, "col_map_offd: \n");
   for (i = 0; i < num_cols_offd; i ++)
   {
      jx_fprintf(fp, "%d\n", col_map_offd[i]);
   }

   fclose(fp);

   return jx_error_flag;
}

/*!
 * \fn jx_hpCSRMatrixGenerateOffdNodeOutside
 * \brief Generate outside node level matrix 
 * \author mrz 
 * \date 2021/08/31
 */
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdNodeOutside(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware)
{
      //test need
   MPI_Comm comm = jx_hpCSRMatrixComm(hp_matrix);
   JX_Int myid;
   jx_MPI_Comm_rank(comm, &myid);
   
   jx_hpCSRMatrixLevel *node_outside_matrix = jx_hpInitMatrixLevel();
   JX_Int  i, j, j_node_outside, col_index, col_map_index;
   jx_CSRMatrix *offd = jx_hpCSRMatrixOffd(hp_matrix);
   JX_Int  num_rows = jx_CSRMatrixNumRows(offd);
   //JX_Int  num_cols = jx_CSRMatrixNumCols(offd);
   JX_Int  num_cols = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   JX_Real *offd_data = jx_CSRMatrixData(offd);
   JX_Int *offd_i = jx_CSRMatrixI(offd);
   JX_Int *offd_j = jx_CSRMatrixJ(offd);
   JX_Int *offd_col_map = jx_hpCSRMatrixColMapOffd(hp_matrix);

   jx_CSRMatrix *offd_node_outside;
   JX_Int  *col_map_offd_node_outside = NULL;
   JX_Real *offd_node_outside_data;
   JX_Int  *offd_node_outside_i;
   JX_Int  *offd_node_outside_j;   

   JX_Int *marker_node_outside;
   JX_Int num_cols_offd_node_outside = 0;
   JX_Int first_elmt = offd_i[0];
   JX_Int num_nonzeros_node_outside;
   JX_Int counter_node_outside;

   JX_Int *row_start = jx_hpCSRMatrixRowStarts(hp_matrix);
   JX_Int  node_colstartnum = row_start[jx_hpProcNodeStart(hardware)];
   JX_Int  node_colendnum   = row_start[jx_hpProcNodeEnd(hardware)+ 1] - 1;
   // //TODO printf待注释
   // if(myid == 0)
   // {
   //    for(i = 0; i < 3; i++)
   //    {
   //       jx_printf("row_start[%d] is %d\n",i, row_start[i]);
   //    }
   //    jx_printf("node_colstartnum is %d, node_colendnum is %d\n",
   //    node_colstartnum, node_colendnum);
   // }

   // generate offd
   offd_node_outside_i = jx_CTAlloc(JX_Int, num_rows + 1);
   offd_node_outside   = jx_CTAlloc(jx_CSRMatrix, 1);
   marker_node_outside = jx_CTAlloc(JX_Int, num_cols);

   for (i = 0; i < num_cols; i ++)
   {
      marker_node_outside[i] = 0;
   }     

   num_nonzeros_node_outside = 0;

   for (i = 0; i < num_rows; i ++)
   {
      offd_node_outside_i[i] = num_nonzeros_node_outside;
   
      for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
      {
         //col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
         // {
         //    continue;
         // }
         // else 
         if(col_map_index < node_colstartnum || col_map_index > node_colendnum)
         {
               // if (!marker_node[col_index])
               // {
               //    marker_node[col_index] = 1;
               //    num_cols_offd_node ++;
               // }
               if (!marker_node_outside[col_map_index])
               {
                  marker_node_outside[col_map_index] = 1;
                  num_cols_offd_node_outside ++;
               }
                  num_nonzeros_node_outside ++;
         }
      }
   }

   if(num_nonzeros_node_outside == 0)
   {
      for (i = 0; i < num_rows + 1; i ++)
      {
         offd_node_outside_i[i] = 0;
      }
      jx_CSRMatrixNumCols(offd_node_outside) = 0;
      jx_CSRMatrixI(offd_node_outside) = offd_node_outside_i;
      jx_hpCSRMatrixOffdlevel(node_outside_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(node_outside_matrix) = offd_node_outside;
   }
   else
   {
      offd_node_outside_i[num_rows] = num_nonzeros_node_outside;
      col_map_offd_node_outside = jx_CTAlloc(JX_Int, num_cols_offd_node_outside);

      counter_node_outside = 0;
      for (i = 0; i < num_cols; i ++)
      {
         if (marker_node_outside[i])
         {
            col_map_offd_node_outside[counter_node_outside] = i;
            marker_node_outside[i] = counter_node_outside;
            counter_node_outside ++;
         }                
      }

      jx_CSRMatrixNumRows(offd_node_outside) = num_rows;
      jx_CSRMatrixNumNonzeros(offd_node_outside) = num_nonzeros_node_outside;
      jx_CSRMatrixNumCols(offd_node_outside) = num_cols_offd_node_outside;
      jx_CSRMatrixJ(offd_node_outside) = NULL;
      jx_CSRMatrixData(offd_node_outside) = NULL;
      jx_CSRMatrixOwnsData(offd_node_outside) = 1;
      jx_CSRMatrixNumRownnz(offd_node_outside) = num_rows;
      jx_CSRMatrixI(offd_node_outside) = offd_node_outside_i;

      jx_CSRMatrixInitialize(offd_node_outside);

      offd_node_outside_data = jx_CSRMatrixData(offd_node_outside);
      offd_node_outside_j = jx_CSRMatrixJ(offd_node_outside);

      j_node_outside = 0;
      for (i = 0; i < num_rows; i ++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
         {
            col_index = offd_j[j];
            col_map_index = offd_col_map[col_index];
            // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
            // {
            //    continue;
            // }            
            // else 
            if (col_map_index < node_colstartnum || col_map_index > node_colendnum)
            {
               offd_node_outside_data[j_node_outside] = offd_data[j];
               offd_node_outside_j[j_node_outside++] = marker_node_outside[col_map_index];
            }
         }    
      }
      jx_CSRMatrixOwnsData(offd_node_outside) = 1;
      jx_hpCSRMatrixOffdlevel(node_outside_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(node_outside_matrix) = offd_node_outside;
      jx_hpCSRMatrixColMapOffdlevel(node_outside_matrix) = jx_CTAlloc(JX_Int, num_cols_offd_node_outside);
      jx_hpCSRMatrixColMapOffdlevel(node_outside_matrix) = col_map_offd_node_outside;      
   }
   
   jx_TFree(marker_node_outside);
   return node_outside_matrix;
}

/*!
 * \fn jx_hpCSRMatrixGenerateOffdNode
 * \brief Generate node level matrix 
 * \author mrz 
 * \date 2021/08/31
 */
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdCpuOutside(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware)
{
      //test need
   // MPI_Comm comm = jx_hpCSRMatrixComm(hp_matrix);
   // JX_Int myid;
   // jx_MPI_Comm_rank(comm, &myid);

   jx_hpCSRMatrixLevel *cpu_outside_matrix = jx_hpInitMatrixLevel();
   JX_Int  i, j, j_cpu_outside, col_index, col_map_index;
   jx_CSRMatrix *offd = jx_hpCSRMatrixOffd(hp_matrix);
   JX_Int  num_rows = jx_CSRMatrixNumRows(offd);
   //JX_Int  num_cols = jx_CSRMatrixNumCols(offd);
   JX_Int  num_cols = jx_hpCSRMatrixGlobalNumCols(hp_matrix);
   JX_Real *offd_data = jx_CSRMatrixData(offd);
   JX_Int *offd_i = jx_CSRMatrixI(offd);
   JX_Int *offd_j = jx_CSRMatrixJ(offd);
   JX_Int *offd_col_map = jx_hpCSRMatrixColMapOffd(hp_matrix);

   jx_CSRMatrix *offd_cpu_outside;
   JX_Int  *col_map_offd_cpu_outside = NULL;
   JX_Real *offd_cpu_outside_data;
   JX_Int  *offd_cpu_outside_i;
   JX_Int  *offd_cpu_outside_j;   

   JX_Int *marker_cpu_outside;
   JX_Int num_cols_offd_cpu_outside = 0;
   JX_Int first_elmt = offd_i[0];
   JX_Int num_nonzeros_cpu_outside;
   JX_Int counter_cpu_outside;

   JX_Int *row_start = jx_hpCSRMatrixRowStarts(hp_matrix);
   JX_Int  node_colstartnum = row_start[jx_hpProcNodeStart(hardware)];
   JX_Int  node_colendnum   = row_start[jx_hpProcNodeEnd(hardware)+ 1] - 1;
   JX_Int cpu_colstartnum = row_start[jx_hpProcCpuStart(hardware)];
   JX_Int cpu_colendnum = row_start[jx_hpProcCpuEnd(hardware) + 1] - 1;
   offd_cpu_outside_i = jx_CTAlloc(JX_Int, num_rows + 1);
   offd_cpu_outside = jx_CTAlloc(jx_CSRMatrix, 1);
   marker_cpu_outside = jx_CTAlloc(JX_Int, num_cols);

   for (i = 0; i < num_cols; i ++)
   {
      marker_cpu_outside[i] = 0;
   }     

   num_nonzeros_cpu_outside = 0;

   for (i = 0; i < num_rows; i ++)
   {
      offd_cpu_outside_i[i] = num_nonzeros_cpu_outside;
   
      for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
      {
         //col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
         {
            continue;
         }
         else if(col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
         {
               // if (!marker_node[col_index])
               // {
               //    marker_node[col_index] = 1;
               //    num_cols_offd_node ++;
               // }
               if (!marker_cpu_outside[col_map_index])
               {
                  marker_cpu_outside[col_map_index] = 1;
                  num_cols_offd_cpu_outside ++;
               }
                  num_nonzeros_cpu_outside ++;
         }
      }
   }

   if(num_nonzeros_cpu_outside == 0)
   {
      for (i = 0; i < num_rows + 1; i ++)
      {
         offd_cpu_outside_i[i] = 0;
      }
      jx_CSRMatrixNumCols(offd_cpu_outside) = 0;
      jx_CSRMatrixI(offd_cpu_outside) = offd_cpu_outside_i;
      jx_hpCSRMatrixOffdlevel(cpu_outside_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(cpu_outside_matrix) = offd_cpu_outside;
   }
   else
   {
      offd_cpu_outside_i[num_rows] = num_nonzeros_cpu_outside;
      col_map_offd_cpu_outside = jx_CTAlloc(JX_Int, num_cols_offd_cpu_outside);

      counter_cpu_outside = 0;
      for (i = 0; i < num_cols; i ++)
      {
         if (marker_cpu_outside[i])
         {
            col_map_offd_cpu_outside[counter_cpu_outside] = i;
            marker_cpu_outside[i] = counter_cpu_outside;
            counter_cpu_outside ++;
         }                
      }

      jx_CSRMatrixNumRows(offd_cpu_outside) = num_rows;
      jx_CSRMatrixNumNonzeros(offd_cpu_outside) = num_nonzeros_cpu_outside;
      jx_CSRMatrixNumCols(offd_cpu_outside) = num_cols_offd_cpu_outside;
      jx_CSRMatrixJ(offd_cpu_outside) = NULL;
      jx_CSRMatrixData(offd_cpu_outside) = NULL;
      jx_CSRMatrixOwnsData(offd_cpu_outside) = 1;
      jx_CSRMatrixNumRownnz(offd_cpu_outside) = num_rows;
      jx_CSRMatrixI(offd_cpu_outside) = offd_cpu_outside_i;

      jx_CSRMatrixInitialize(offd_cpu_outside);

      offd_cpu_outside_data = jx_CSRMatrixData(offd_cpu_outside);
      offd_cpu_outside_j = jx_CSRMatrixJ(offd_cpu_outside);

      j_cpu_outside = 0;
      for (i = 0; i < num_rows; i ++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i+1] - first_elmt; j ++)
         {
            col_index = offd_j[j];
            col_map_index = offd_col_map[col_index];
            if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
            {
               continue;
            }            
            else if (col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
            {
               offd_cpu_outside_data[j_cpu_outside] = offd_data[j];
               offd_cpu_outside_j[j_cpu_outside++] = marker_cpu_outside[col_map_index];
            }
         }    
      }
      jx_CSRMatrixOwnsData(offd_cpu_outside) = 1;
      jx_hpCSRMatrixOffdlevel(cpu_outside_matrix) = jx_CTAlloc(jx_CSRMatrix,  1); 
      jx_hpCSRMatrixOffdlevel(cpu_outside_matrix) = offd_cpu_outside;
      jx_hpCSRMatrixColMapOffdlevel(cpu_outside_matrix) = jx_CTAlloc(JX_Int, num_cols_offd_cpu_outside);
      jx_hpCSRMatrixColMapOffdlevel(cpu_outside_matrix) = col_map_offd_cpu_outside;
   }

   jx_TFree(marker_cpu_outside);
   return cpu_outside_matrix;
}
