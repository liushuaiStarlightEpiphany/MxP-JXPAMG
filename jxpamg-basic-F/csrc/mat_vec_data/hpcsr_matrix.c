#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

#ifndef JXF_ENV_VARIABLE_H
#include "jxf_env_variable.h"
#endif

// add global variable jxf_hp_hardware
jxf_HardwareInfo *jxf_hp_hardware = NULL;

#ifdef USING_HWLOC
#include "hwloc.h"

/*!
 * \fn jxf_HardwareInfo *jxf_hpCreateHardwareInfo
 * \brief Create Hardware Info.
 * \date 2021/10/15 create
 * \date 2022/11/24 change hardware to global variable
 * \author mrz
 */
JXF_Int jxf_hpCreateHardwareInfo(MPI_Comm comm)
{
   if (jxf_hp_hardware == NULL)
   {
      jxf_hp_hardware = jxf_CTAlloc(jxf_HardwareInfo, 1);
      // add process binding operation
      /*
      0:bind to node
      1:bind to cpu
      2:bind to core
      */
      jxf_hpBindType(jxf_hp_hardware) = jxf_getenv(JXF_Bind_Type);

      jxf_hpComm(jxf_hp_hardware) = comm;
      jxf_hpCoreComm(jxf_hp_hardware) = MPI_COMM_NULL;
      jxf_hpIntraCpuComm(jxf_hp_hardware) = MPI_COMM_NULL;
      jxf_hpNodeComm(jxf_hp_hardware) = MPI_COMM_NULL;
      // These comm is not used
      //  jxf_hpCpuComm(jxf_hp_hardware)        = MPI_COMM_NULL;
      //  jxf_hpIntraNodeComm(jxf_hp_hardware)	= MPI_COMM_NULL;

      jxf_hpBindToType_Hwloc(jxf_hp_hardware);
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_hpGetLevelBindId
 * \brief use api of mpi according myid bind to type.
 * \date 2022/11/30
 * \author mrz
 */
JXF_Int jxf_hpGetLevelBindId(JXF_Int myid, JXF_Int numprocess, JXF_Int num_level,
                             JXF_Int *id_in_next_level, JXF_Int *process_id_in_level, JXF_Int *numprocess_in_level, JXF_Int *start_process_id_in_level, JXF_Int *end_process_id_in_level)
{
   JXF_Int process_per_in_level;
   JXF_Int r_process_per_in_level;
   JXF_Int my_id_in_level, my_process_id_in_level, my_numprocess_in_level;
   JXF_Int my_start_process_id_in_level, my_end_process_id_in_level;

   process_per_in_level = numprocess / num_level;
   r_process_per_in_level = numprocess - process_per_in_level * num_level;

   if (myid < r_process_per_in_level * (process_per_in_level + 1))
   {
      my_numprocess_in_level = process_per_in_level + 1;
      my_id_in_level = myid / (my_numprocess_in_level);
      my_process_id_in_level = myid % (my_numprocess_in_level);
   }
   else
   {
      my_numprocess_in_level = process_per_in_level;
      my_id_in_level = r_process_per_in_level + (myid - r_process_per_in_level * (process_per_in_level + 1)) / process_per_in_level;
      my_process_id_in_level = (myid - r_process_per_in_level * (process_per_in_level + 1)) % process_per_in_level;
   }

   if (numprocess < num_level)
   {
      JXF_Int remainder = num_level - numprocess;
      JXF_Int remainder_in_process = remainder / numprocess + 1;
      if (remainder % numprocess)
         remainder_in_process += 1;
      if (myid < remainder)
      {
         my_id_in_level = myid * remainder_in_process;
      }
      else
      {
         my_id_in_level = remainder * remainder_in_process + (myid - remainder) * (num_level / numprocess);
      }
   }

   my_start_process_id_in_level = myid - my_process_id_in_level;
   my_end_process_id_in_level = my_start_process_id_in_level + my_numprocess_in_level - 1;

   *id_in_next_level = my_id_in_level;
   if (process_id_in_level)
      *process_id_in_level = my_process_id_in_level;
   if (numprocess_in_level)
      *numprocess_in_level = my_numprocess_in_level;
   if (start_process_id_in_level)
      *start_process_id_in_level = my_start_process_id_in_level;
   if (end_process_id_in_level)
      *end_process_id_in_level = my_end_process_id_in_level;
   return 0;
}

/*!
 * \fn JXF_Int jxf_hpBindToType_Hwloc
 * \brief use api of hwloc according myid bind to type.
 * \date 2021/08/26
 * \author mrz
 */

// TODO:Comments need to be removed

JXF_Int jxf_hpBindToType_Hwloc(jxf_HardwareInfo *hardware)
{
   JXF_Int myid, numprocess, global_pid;
   MPI_Comm comm = jxf_hpComm(hardware);

   JXF_Int type = jxf_hpBindType(hardware);
   // JXF_Int type_id_in_node, num_type_in_node;
   JXF_Int num_numanode_in_machine, num_core_in_numanode, num_machine, num_core_in_total_machine, num_numanode_in_total_machine;
   JXF_Int process_per_numanode, process_per_core, my_machine_id, my_numanode_id, pid_in_numanode, pid_in_core, my_core_id;
   JXF_Int process_statr_in_machine, process_statr_in_numanode;
   JXF_Int process_end_in_machine, process_end_in_numanode;
   hwloc_topology_t topology;
   hwloc_cpuset_t cpuset;
   hwloc_obj_t obj_in_numanode, obj_in_core, obj_in_machine;
   JXF_Int err, i, j, k, tmp_index;
   JXF_Int ***hardware_ID;

   MPI_Comm shmcomm = MPI_COMM_NULL;
   // get communicator in the same node
   MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
   jxf_MPI_Comm_rank(shmcomm, &myid);
   jxf_MPI_Comm_size(shmcomm, &numprocess);
   jxf_MPI_Comm_rank(comm, &global_pid);

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

   num_machine = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_MACHINE);
   hardware_ID = jxf_CTAlloc(int **, num_machine + 1);
   hardware_ID[0] = jxf_CTAlloc(int *, 1);
   hardware_ID[0][0] = jxf_CTAlloc(int, 1);
   hardware_ID[0][0][0] = num_machine;

   tmp_index = 0;
   num_core_in_numanode = 0;
   num_numanode_in_machine = 0;
   num_numanode_in_total_machine = 0;
   num_core_in_total_machine = 0;

   for (k = 0; k < num_machine; k++)
   {
      obj_in_machine = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, k);
      num_numanode_in_machine = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj_in_machine->cpuset, HWLOC_OBJ_NUMANODE);
      num_numanode_in_total_machine = num_numanode_in_total_machine + num_numanode_in_machine;
      hardware_ID[k + 1] = jxf_CTAlloc(int *, num_numanode_in_machine + 1);
      hardware_ID[k + 1][0] = jxf_CTAlloc(int, 1);
      hardware_ID[k + 1][0][0] = num_numanode_in_machine;
      for (i = 0; i < num_numanode_in_machine; i++)
      {
         obj_in_numanode = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, i);
         num_core_in_numanode = hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj_in_numanode->cpuset, HWLOC_OBJ_CORE);
         num_core_in_total_machine = num_core_in_total_machine + num_core_in_numanode;
         hardware_ID[k + 1][i + 1] = jxf_CTAlloc(int, num_core_in_numanode + 1);
         hardware_ID[k + 1][i + 1][0] = num_core_in_numanode;

         for (j = 0; j < num_core_in_numanode; j++)
         {
            // wrtite CPU index id
            hardware_ID[k + 1][i + 1][j + 1] = tmp_index + j;
         }
         tmp_index = tmp_index + hardware_ID[k + 1][i + 1][0];
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
               jxf_printf("hardware_ID[%d][%d][%d]  is %d\n", k, i, j, hardware_ID[k][i][j] );
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
   if (type == 0)
   {
      jxf_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);

      obj_in_machine = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, my_machine_id);

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

      jxf_hpNumCoreInCpu(hardware) = hardware_ID[my_machine_id + 1][1][0];
      jxf_hpProcCpuStart(hardware) = 0;
      jxf_hpProcCpuEnd(hardware) = 0;
      jxf_hpProcNodeStart(hardware) = 0;
      jxf_hpProcNodeEnd(hardware) = 0;
   }
   else if (type == 1)
   {
      jxf_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);
      jxf_hpGetLevelBindId(pid_in_numanode, process_per_numanode, hardware_ID[my_machine_id + 1][0][0], &my_numanode_id, &pid_in_core, &process_per_core, &process_statr_in_numanode, &process_end_in_numanode);
#if 0
      jxf_printf("myid is %d, start_node_p_id is %d, end_node_p_id is %d\n", myid, process_statr_in_machine, process_end_in_machine);
      jxf_printf("myid is %d, start_cpu_p_id is %d, end_cpu_p_id is %d\n", myid, process_statr_in_numanode, process_end_in_numanode);
#endif
      for (i = 0; i < my_machine_id; i++)
      {
         my_numanode_id = my_numanode_id + hardware_ID[my_machine_id + 1][0][0];
      }

#if 0
         jxf_printf("myid is %d, hardware_ID is %d\n", myid, my_numanode_id);
#endif

      obj_in_numanode = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, my_numanode_id);

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

      jxf_hpNumCoreInCpu(hardware) = process_per_numanode;
      jxf_hpProcCpuStart(hardware) = 0;
      jxf_hpProcCpuEnd(hardware) = 0;
      jxf_hpProcNodeStart(hardware) = global_pid - myid;
      jxf_hpProcNodeEnd(hardware) = global_pid - myid + numprocess - 1;
   }
   else if (type == 2)
   { // if(myid == 0)jxf_printf("myid is %d, numproces is %d, num_machine is %d\n ", myid, numprocess, num_machine);
      jxf_hpGetLevelBindId(myid, numprocess, num_machine, &my_machine_id, &pid_in_numanode, &process_per_numanode, &process_statr_in_machine, &process_end_in_machine);
      // if(myid == 0)jxf_printf("my_machine_id is %d, pid_in_numanode is %d, process_per_numanode is %d\n ", my_machine_id, pid_in_numanode, process_per_numanode);
      jxf_hpGetLevelBindId(pid_in_numanode, process_per_numanode, hardware_ID[my_machine_id + 1][0][0], &my_numanode_id, &pid_in_core, &process_per_core, &process_statr_in_numanode, &process_end_in_numanode);
      // if(myid == 0)jxf_printf("my_numanode_id is %d, pid_in_core is %d, process_per_core is %d\n ", my_numanode_id, pid_in_core, process_per_core);
      // if(myid == 0)jxf_printf("pid_in_core is %d, process_per_core is %d, hardware_ID is %d\n ", pid_in_core, process_per_core, hardware_ID[my_machine_id+1][my_numanode_id+1][0]);
      jxf_hpGetLevelBindId(pid_in_core, process_per_core, hardware_ID[my_machine_id + 1][my_numanode_id + 1][0], &my_core_id, NULL, NULL, NULL, NULL);
      // jxf_printf("myid is %d, my_machine_id is %d, my_numanode_id is %d, my_core_id %d\n ",myid, my_machine_id, my_numanode_id, my_core_id);

      /*
      // use another bind way
      obj_in_machine  = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, my_machine_id);
      obj_in_numanode = hwloc_get_obj_inside_cpuset_by_type(topology, obj_in_machine->cpuset, HWLOC_OBJ_NUMANODE, my_numanode_id);
      obj_in_core = hwloc_get_obj_inside_cpuset_by_type(topology, obj_in_numanode->cpuset, HWLOC_OBJ_CORE, my_core_id);
      */
      obj_in_core = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, hardware_ID[my_machine_id + 1][my_numanode_id + 1][my_core_id + 1]);

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

      jxf_hpNumCoreInCpu(hardware) = process_per_numanode;
      jxf_hpProcNodeStart(hardware) = global_pid - myid;
      jxf_hpProcNodeEnd(hardware) = global_pid - myid + numprocess - 1;
      jxf_hpProcCpuStart(hardware) = jxf_hpProcNodeStart(hardware) + process_statr_in_machine + process_statr_in_numanode;
      // if(global_pid == 0) jxf_printf("process_statr_in_machine is %d, process_end_in_numanode is %d\n", process_statr_in_machine, process_end_in_numanode);
      jxf_hpProcCpuEnd(hardware) = jxf_hpProcNodeStart(hardware) + process_statr_in_machine + process_end_in_numanode;
   }

   jxf_hpNumCpuInNode(hardware) = num_numanode_in_total_machine;
   jxf_hpNumCoreInNode(hardware) = numprocess;
#if 0
      if(myid == 0)
      {
         jxf_printf("myid is %d, hardware_ID is %d, process_per_numanode is %d\n", myid, hardware_ID[my_machine_id+1][my_numanode_id+1][my_core_id+1], process_per_numanode);
         jxf_printf("myid is %d, start_node_p_id is %d, end_node_p_id is %d\n", myid, jxf_hpProcNodeStart(hardware), jxf_hpProcNodeEnd(hardware));
         jxf_printf("myid is %d, start_cpu_p_id is %d, end_cpu_p_id is %d\n", myid, jxf_hpProcCpuStart(hardware), jxf_hpProcCpuEnd(hardware));
      }

#endif

// createhardware Allreduce Test
#if 1
   /**
    * [0, 1, [2],| 3, 4, [[5]], ||6, 7, [8],| 9, 10, [[11]]]
    * || node level, | cpu level]
    * [] cpu represent
    * [[]] node represent
    */

   // definition
   MPI_Comm CORE_COMM, Intra_CPU_COMM, NODE_COMM; // CPU_COMM;

   // init
   CORE_COMM = MPI_COMM_NULL;
   Intra_CPU_COMM = MPI_COMM_NULL;
   NODE_COMM = MPI_COMM_NULL;
   // CPU_COMM  =  MPI_COMM_NULL;

   // 同一Core内的Comm
   jxf_MPI_Comm_split(shmcomm, jxf_hpProcCpuStart(hardware), global_pid, &CORE_COMM);

   // 结点内CPU代表Comm
   // JXF_Int cpu_color = global_pid % process_per_core;
   JXF_Int Intra_cpu_color = 1;
   if (pid_in_core != process_per_core - 1)
   {
      Intra_cpu_color = MPI_UNDEFINED;
   }
   jxf_MPI_Comm_split(shmcomm, Intra_cpu_color, global_pid, &Intra_CPU_COMM);

   // 所有CPU代表Comm(暂未使用到)
   // JXF_Int cpu_color = global_pid % process_per_core;
   // JXF_Int cpu_color = 2;
   // if(pid_in_core != process_per_core - 1)
   // {
   //    cpu_color = MPI_UNDEFINED;
   // }
   // jxf_MPI_Comm_split(comm, cpu_color, global_pid, &CPU_COMM);

   // Node代表Comm
   // JXF_Int node_color = global_pid % process_per_numanode;
   JXF_Int node_color = 3;
   if (pid_in_numanode != process_per_numanode - 1)
   {
      node_color = MPI_UNDEFINED;
   }
   jxf_MPI_Comm_split(comm, node_color, global_pid, &NODE_COMM);

   jxf_hpCoreComm(hardware) = CORE_COMM;
   jxf_hpIntraCpuComm(hardware) = Intra_CPU_COMM;
   jxf_hpNodeComm(hardware) = NODE_COMM;
   // jxf_hpCpuComm(hardware)	= CPU_COMM;
   // jxf_hpIntraNodeComm(hardware) = shmcomm;

#if 0
      JXF_Int  mycoreid=-1, myintracpuid=-1, mynodeid=-1, mycpuid=-1;
      if(CORE_COMM != MPI_COMM_NULL) jxf_MPI_Comm_rank(CORE_COMM, &mycoreid);
      if(CPU_COMM != MPI_COMM_NULL) jxf_MPI_Comm_rank(CPU_COMM, &mycpuid);
      if(Intra_CPU_COMM != MPI_COMM_NULL) jxf_MPI_Comm_rank(Intra_CPU_COMM, &myintracpuid);
      if(NODE_COMM != MPI_COMM_NULL) jxf_MPI_Comm_rank(NODE_COMM, &mynodeid);
      jxf_printf("myid is %d, mycoreid is %d, mycpuid is %d, myintracpuid is %d, mynodeid is %d\n",
         global_pid, mycoreid, mycpuid, myintracpuid, mynodeid);
      jxf_MPI_Barrier(comm);
      // exit(0);
#endif

#endif

   hwloc_topology_destroy(topology);
   jxf_MPI_Comm_free(&shmcomm);

   for (k = 0; k < num_machine + 1; k++)
   {
      num_numanode_in_machine = hardware_ID[k][0][0];
      for (i = 0; i < num_numanode_in_machine + 1; i++)
      {
         jxf_TFree(hardware_ID[k][i]);
         if (k == 0)
            break;
      }
      jxf_TFree(hardware_ID[k]);
   }
   jxf_TFree(hardware_ID);

   return EXIT_SUCCESS;
}
#endif
/*!
 * \fn jxf_hpCSRMatrixSetCpuLevel
 * \brief Get col_start and col_end of offd_cpu
 * \author mrz
 * \date 2021/10/15
 */
JXF_Int jxf_hpCSRMatrixSetCpuLevel(jxf_HardwareInfo *hardware)
{
   JXF_Int my_id, nprocs;
   MPI_Comm comm = jxf_hpComm(hardware);
   JXF_Int num_core_in_cpu = jxf_hpNumCoreInCpu(hardware);
   JXF_Int num_core_in_node = jxf_hpNumCoreInNode(hardware);
   JXF_Int bind_type = jxf_hpBindType(hardware);

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &nprocs);

   if (bind_type == 0 || bind_type == 1) // 绑定到CPU层或node层
   {
      jxf_hpProcCpuStart(hardware) = 0;
      jxf_hpProcCpuEnd(hardware) = 0;
   }
   else if (num_core_in_cpu == 1)
   {
      jxf_hpProcCpuStart(hardware) = my_id;
      jxf_hpProcCpuEnd(hardware) = my_id;
   }
   else
   {
      MPI_Comm Cpu_Comm = MPI_COMM_NULL;
      JXF_Int core_id_start_in_node = my_id / num_core_in_node * num_core_in_node;
      JXF_Int core_id_in_node = my_id % num_core_in_node;
      JXF_Int core_id_start_in_cpu = core_id_in_node / num_core_in_cpu * num_core_in_cpu;
      JXF_Int process_id_start_in_cpu = core_id_start_in_node + core_id_start_in_cpu;

      JXF_Int end_core_id_start_in_cpu = nprocs / num_core_in_cpu * num_core_in_cpu;
      JXF_Int r_cpu = nprocs % num_core_in_cpu;
      if (r_cpu != 0 && my_id >= end_core_id_start_in_cpu)
      {
         jxf_hpProcCpuStart(hardware) = process_id_start_in_cpu;
         jxf_hpProcCpuEnd(hardware) = process_id_start_in_cpu + r_cpu - 1;
      }
      else
      {
         jxf_hpProcCpuStart(hardware) = process_id_start_in_cpu;
         jxf_hpProcCpuEnd(hardware) = process_id_start_in_cpu + num_core_in_cpu - 1;
      }
#if 0
      jxf_MPI_Comm_split(comm, process_id_start_in_cpu, my_id, &Cpu_Comm);
#endif
      jxf_hpIntraCpuComm(hardware) = Cpu_Comm;
   }

   return jxf_error_flag;
}

/*!
 * \fn jxf_hpCSRMatrixGetNodeCol
 * \brief Get col_start and col_end of offd_node
 * \author mrz
 * \date 2021/08/30
 */
JXF_Int jxf_hpCSRMatrixSetNodeLevel(jxf_HardwareInfo *hardware)
{
   JXF_Int my_id, nprocs;
   MPI_Comm comm = jxf_hpComm(hardware);
   JXF_Int num_cpu_in_node = jxf_hpNumCpuInNode(hardware);
   JXF_Int num_core_in_node = jxf_hpNumCoreInNode(hardware);
   JXF_Int bind_type = jxf_hpBindType(hardware);

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &nprocs);

   if (bind_type == 0) // 绑定到node层
   {
      jxf_hpProcNodeStart(hardware) = 0;
      jxf_hpProcNodeEnd(hardware) = 0;
   }
   else if (num_cpu_in_node == 1)
   {
      jxf_hpProcNodeStart(hardware) = jxf_hpProcCpuStart(hardware);
      jxf_hpProcNodeEnd(hardware) = jxf_hpProcCpuEnd(hardware);
      jxf_hpNodeComm(hardware) = jxf_hpIntraCpuComm(hardware);
   }
   else
   {
      if (bind_type == 2)
      {
         MPI_Comm Node_Comm = MPI_COMM_NULL;
         JXF_Int process_id_start_in_node = my_id / num_core_in_node * num_core_in_node;

         JXF_Int end_core_id_start_in_node = nprocs / num_core_in_node * num_core_in_node;
         JXF_Int r_node = nprocs % num_core_in_node;
         if (r_node != 0 && my_id >= end_core_id_start_in_node)
         {
            jxf_hpProcNodeStart(hardware) = process_id_start_in_node;
            jxf_hpProcNodeEnd(hardware) = process_id_start_in_node + r_node - 1;
         }
         else
         {
            jxf_hpProcNodeStart(hardware) = process_id_start_in_node;
            jxf_hpProcNodeEnd(hardware) = process_id_start_in_node + num_core_in_node - 1;
         }
#if 0
         jxf_MPI_Comm_split(comm, process_id_start_in_node, my_id, &Node_Comm);
#endif
         jxf_hpNodeComm(hardware) = Node_Comm;
      }
      else
      {
         MPI_Comm Node_Comm = MPI_COMM_NULL;
         JXF_Int process_id_start_in_node = my_id / num_cpu_in_node * num_cpu_in_node;

         JXF_Int end_cpu_id_start_in_node = nprocs / num_cpu_in_node * num_cpu_in_node;
         JXF_Int r_node = nprocs % num_cpu_in_node;
         if (r_node != 0 && my_id >= end_cpu_id_start_in_node)
         {
            jxf_hpProcNodeStart(hardware) = process_id_start_in_node;
            jxf_hpProcNodeEnd(hardware) = process_id_start_in_node + r_node - 1;
         }
         else
         {
            jxf_hpProcNodeStart(hardware) = process_id_start_in_node;
            jxf_hpProcNodeEnd(hardware) = process_id_start_in_node + num_cpu_in_node - 1;
         }
#if 0
         jxf_MPI_Comm_split(comm, process_id_start_in_node, my_id, &Node_Comm);
#endif
         jxf_hpNodeComm(hardware) = Node_Comm;
      }
   }
   return jxf_error_flag;
}

/*!
 * \fn jxf_hpBuildMatParFromOneFile
 * \brief read matrix from files and create hierarchy diag block.
 * \date 2021/10/27
 * \edit 2022/3/7
 * \author mrz
 */
jxf_hpCSRMatrix *
jxf_hpBuildMatParFromOneFile(char *filename, JXF_Int num_functions, JXF_Int file_base)
{
   jxf_hpCSRMatrix *A = jxf_hpInithpCSRMatrix();
   jxf_hpCSRMatrixPar(A) = jxf_BuildMatParFromOneFile(filename, num_functions, file_base);
#ifdef USING_HWLOC
   A = jxf_hpCreateMatrixLevelBlock(A);
#endif
   return (A);
}

/*!
 * \fn jxf_hpCSRMatrixCreate
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrix *jxf_hpCSRMatrixCreate(MPI_Comm comm,
                                       JXF_Int global_num_rows,
                                       JXF_Int global_num_cols,
                                       JXF_Int *row_starts,
                                       JXF_Int *col_starts,
                                       JXF_Int num_cols_offd,
                                       JXF_Int num_nonzeros_diag,
                                       JXF_Int num_nonzeros_offd)
{
   jxf_hpCSRMatrix *matrix = jxf_hpInithpCSRMatrix();
   jxf_hpCSRMatrixPar(matrix) = jxf_ParCSRMatrixCreate(comm, global_num_rows, global_num_cols,
                                                       row_starts, col_starts, num_cols_offd, num_nonzeros_diag, num_nonzeros_offd);
#ifdef USING_HWLOC
   matrix = jxf_hpCreateMatrixLevelBlock(matrix);
#endif
   return matrix;
}

/*!
 * \fn jxf_hpCSRMatrixInitialize
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
JXF_Int jxf_hpCSRMatrixInitialize(jxf_hpCSRMatrix *A)
{
   if (!jxf_hpCSRMatrixPar(A))
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParCSRMatrixInitialize(jxf_hpCSRMatrixPar(A));
   return jxf_error_flag;
}

/*!
 * \fn jxf_hpInithpCSRMatrix
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrix *jxf_hpInithpCSRMatrix()
{
   jxf_hpCSRMatrix *A = jxf_CTAlloc(jxf_hpCSRMatrix, 1);
   jxf_hpCSRMatrixPar(A) = NULL;
   jxf_hpCSRMatrixCpu(A) = jxf_hpInitMatrixLevel();
   jxf_hpCSRMatrixNode(A) = jxf_hpInitMatrixLevel();
   jxf_hpCSRMatrixCore(A) = jxf_hpInitMatrixLevel();
   return A;
}

/*!
 * \fn jxh_hpInithpCSRMatrix
 * \brief Init hpcsr matrix.
 * \date 2021/11/13
 * \author mrz
 */
jxh_hpCSRMatrix *jxh_hpInithpCSRMatrix()
{
   jxh_hpCSRMatrix *A = jxf_CTAlloc(jxh_hpCSRMatrix, 1);

   jxh_hpCSRMatrixPar(A) = NULL;
   jxh_hpCSRMatrixCpu(A) = NULL;
   jxh_hpCSRMatrixNode(A) = NULL;
   jxh_hpCSRMatrixCore(A) = NULL;

   return A;
}

/*!
 * \fn jxf_hpInitMatrixLevel
 * \brief Init Matrix Level.
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrixLevel *jxf_hpInitMatrixLevel()
{
   jxf_hpCSRMatrixLevel *level_matrix = jxf_CTAlloc(jxf_hpCSRMatrixLevel, 1);
   jxf_hpCSRMatrixOffdlevel(level_matrix) = NULL;
   jxf_hpCSRMatrixColMapOffdlevel(level_matrix) = NULL;
   jxf_hpCSRMatrixlevelCommPkg(level_matrix) = NULL;
   return level_matrix;
}

/*!
 * \fn jxf_createMatrixLevel
 * \brief create block of each Matrix Level.
 * \date 2022/3/7
 * \author mrz
 */
jxf_hpCSRMatrix *jxf_hpCreateMatrixLevelBlock(jxf_hpCSRMatrix *A)
{
   if (!jxf_hpCSRMatrixCpu(A))
      jxf_hpCSRMatrixCpu(A) = jxf_hpCreateMatrixCpu(A);
   if (!jxf_hpCSRMatrixNode(A))
      jxf_hpCSRMatrixNode(A) = jxf_hpCreateMatrixNode(A);
   if (!jxf_hpCSRMatrixCore(A))
      jxf_hpCSRMatrixCore(A) = jxf_hpCreateMatrixCore(A);
   return A;
}

/*!
 * \fn jxf_hpCreateMatrixCpu
 * \brief create cpu level of matrix
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixCore(jxf_hpCSRMatrix *A)
{
   JXF_Int i;
   JXF_Int num_rows = jxf_hpCSRMatrixNumRows(A);
   jxf_hpCSRMatrixLevel *core_matrix = jxf_CTAlloc(jxf_hpCSRMatrixLevel, 1);
   jxf_CSRMatrix *offd_core = jxf_CTAlloc(jxf_CSRMatrix, 1);
   JXF_Int *offd_core_i = jxf_CTAlloc(JXF_Int, num_rows + 1);

   for (i = 0; i < num_rows + 1; i++)
   {
      offd_core_i[i] = 0;
   }
   jxf_CSRMatrixNumCols(offd_core) = 0;
   jxf_CSRMatrixI(offd_core) = offd_core_i;
   jxf_hpCSRMatrixOffdlevel(core_matrix) = offd_core;
   return core_matrix;
}

/*!
 * \fn jxf_hpCreateMatrixCpu
 * \brief create cpu level of matrix
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixCpu(jxf_hpCSRMatrix *A)
{
   if (jxf_hp_hardware)
   {
      // jxf_hpCSRMatrixLevel *cpu_matrix = NULL;
      // JXF_Int bind_type = jxf_hpBindType(hardware);
      // JXF_Int num_core_in_cpu = jxf_hpNumCoreInCpu(hardware);

      jxf_hpCSRMatrixLevel *cpu_matrix = jxf_hpCSRMatrixGenerateOffdCpu(A, jxf_hp_hardware);
      return cpu_matrix;
   }
   else
   {
      jxf_printf("erro: CreateMatrixCpu, hardwareinfo is empty\n");
      return NULL;
   }
}

/*!
 * \fn JXF_Int jxf_hpCreateMatrixNode
 * \brief create node level of matrix
 * \date 2021/11/13
 * \author mrz
 */
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixNode(jxf_hpCSRMatrix *A)
{
   if (jxf_hp_hardware)
   {
      // JXF_Int num_cpu_in_node = jxf_hpNumCpuInNode(hardware);
      // JXF_Int num_core_in_cpu = jxf_hpNumCoreInCpu(hardware);

      // JXF_Int bind_type = jxf_hpBindType(hardware);
      jxf_hpCSRMatrixLevel *node_matrix = jxf_hpCSRMatrixGenerateOffdNode(A, jxf_hp_hardware);
      return node_matrix;
   }
   else
   {
      jxf_printf("erro: CreateMatrixNode, hardwareinfo is empty\n");
      return NULL;
   }
}

/*!
 * \fn jxf_hpCSRMatrixGenerateOffdCpu
 * \brief Generate cpu level matrix
 * \author mrz
 * \date 2021/08/31
 */
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdCpu(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware)
{
   // test need
   //  MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   //  JXF_Int myid;
   //  jxf_MPI_Comm_rank(comm, &myid);

   jxf_hpCSRMatrixLevel *cpu_matrix = jxf_hpInitMatrixLevel();

   JXF_Int i, j, j_cpu, col_index, col_map_index;
   jxf_CSRMatrix *offd = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(offd);
   // JXF_Int  num_cols = jxf_CSRMatrixNumCols(offd);
   JXF_Int num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Int *offd_col_map = jxf_hpCSRMatrixColMapOffd(hp_matrix);

   jxf_CSRMatrix *offd_cpu;
   JXF_Int *col_map_offd_cpu = NULL;
   JXF_Real *offd_cpu_data;
   JXF_Int *offd_cpu_i;
   JXF_Int *offd_cpu_j;

   JXF_Int *marker_cpu;
   JXF_Int num_cols_offd_cpu = 0;
   JXF_Int first_elmt = offd_i[0];
   JXF_Int num_nonzeros_cpu;
   JXF_Int counter_cpu;

   JXF_Int *row_start = jxf_hpCSRMatrixRowStarts(hp_matrix);
   JXF_Int cpu_colstartnum = row_start[jxf_hpProcCpuStart(hardware)];
   JXF_Int cpu_colendnum = row_start[jxf_hpProcCpuEnd(hardware) + 1] - 1;

   // generate offd_cpu and col_cpu_map
   offd_cpu_i = jxf_CTAlloc(JXF_Int, num_rows + 1);
   offd_cpu = jxf_CTAlloc(jxf_CSRMatrix, 1);
   marker_cpu = jxf_CTAlloc(JXF_Int, num_cols);

   for (i = 0; i < num_cols; i++)
   {
      marker_cpu[i] = 0;
   }

   num_nonzeros_cpu = 0;

   for (i = 0; i < num_rows; i++)
   {
      offd_cpu_i[i] = num_nonzeros_cpu;

      for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
      {
         // col_index = offd_j[j]; offd_col_map[col_index]
         col_map_index = offd_col_map[offd_j[j]];
         if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
         {
            if (!marker_cpu[col_map_index])
            {
               marker_cpu[col_map_index] = 1;
               num_cols_offd_cpu++;
            }
            num_nonzeros_cpu++;
         }
      }
   }

   if (num_nonzeros_cpu == 0)
   {
      for (i = 0; i < num_rows + 1; i++)
      {
         offd_cpu_i[i] = 0;
      }
      jxf_CSRMatrixNumCols(offd_cpu) = 0;
      jxf_CSRMatrixI(offd_cpu) = offd_cpu_i;
      jxf_hpCSRMatrixOffdlevel(cpu_matrix) = offd_cpu;
   }
   else
   {
      offd_cpu_i[num_rows] = num_nonzeros_cpu;
      col_map_offd_cpu = jxf_CTAlloc(JXF_Int, num_cols_offd_cpu);

      counter_cpu = 0;
      for (i = 0; i < num_cols; i++)
      {
         if (marker_cpu[i])
         {
            col_map_offd_cpu[counter_cpu] = i;
            marker_cpu[i] = counter_cpu;
            counter_cpu++;
         }
      }
      jxf_CSRMatrixNumRows(offd_cpu) = num_rows;
      jxf_CSRMatrixNumNonzeros(offd_cpu) = num_nonzeros_cpu;
      jxf_CSRMatrixNumCols(offd_cpu) = num_cols_offd_cpu;
      jxf_CSRMatrixJ(offd_cpu) = NULL;
      jxf_CSRMatrixData(offd_cpu) = NULL;
      jxf_CSRMatrixOwnsData(offd_cpu) = 1;
      jxf_CSRMatrixNumRownnz(offd_cpu) = num_rows;
      jxf_CSRMatrixI(offd_cpu) = offd_cpu_i;

      jxf_CSRMatrixInitialize(offd_cpu);

      offd_cpu_data = jxf_CSRMatrixData(offd_cpu);
      offd_cpu_j = jxf_CSRMatrixJ(offd_cpu);

      j_cpu = 0;
      for (i = 0; i < num_rows; i++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
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
      jxf_CSRMatrixOwnsData(offd_cpu) = 1;
      jxf_hpCSRMatrixOffdlevel(cpu_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(cpu_matrix) = offd_cpu;
      jxf_hpCSRMatrixColMapOffdlevel(cpu_matrix) = jxf_CTAlloc(JXF_Int, num_cols_offd_cpu);
      jxf_hpCSRMatrixColMapOffdlevel(cpu_matrix) = col_map_offd_cpu;
   }

   jxf_TFree(marker_cpu);
   return cpu_matrix;
}

/*!
 * \fn jxf_hpCSRMatrixGenerateOffdNode
 * \brief Generate node level matrix
 * \author mrz
 * \date 2021/08/31
 */
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdNode(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware)
{
   // test need
   // MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   // JXF_Int myid;
   // jxf_MPI_Comm_rank(comm, &myid);

   jxf_hpCSRMatrixLevel *node_matrix = jxf_hpInitMatrixLevel();
   JXF_Int i, j, j_node, col_index, col_map_index;
   jxf_CSRMatrix *offd = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(offd);
   // JXF_Int  num_cols = jxf_CSRMatrixNumCols(offd);
   JXF_Int num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Int *offd_col_map = jxf_hpCSRMatrixColMapOffd(hp_matrix);

   jxf_CSRMatrix *offd_node;
   JXF_Int *col_map_offd_node = NULL;
   JXF_Real *offd_node_data;
   JXF_Int *offd_node_i;
   JXF_Int *offd_node_j;

   JXF_Int *marker_node;
   JXF_Int num_cols_offd_node = 0;
   JXF_Int first_elmt = offd_i[0];
   JXF_Int num_nonzeros_node;
   JXF_Int counter_node;

   JXF_Int *row_start = jxf_hpCSRMatrixRowStarts(hp_matrix);
   JXF_Int node_colstartnum = row_start[jxf_hpProcNodeStart(hardware)];
   JXF_Int node_colendnum = row_start[jxf_hpProcNodeEnd(hardware) + 1] - 1;

   // generate offd
   offd_node_i = jxf_CTAlloc(JXF_Int, num_rows + 1);
   offd_node = jxf_CTAlloc(jxf_CSRMatrix, 1);
   marker_node = jxf_CTAlloc(JXF_Int, num_cols);

   for (i = 0; i < num_cols; i++)
   {
      marker_node[i] = 0;
   }

   num_nonzeros_node = 0;

   for (i = 0; i < num_rows; i++)
   {
      offd_node_i[i] = num_nonzeros_node;

      for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
      {
         // col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
         // {
         //    continue;
         // }
         // else
         if (col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
         {
            // if (!marker_node[col_index])
            // {
            //    marker_node[col_index] = 1;
            //    num_cols_offd_node ++;
            // }
            if (!marker_node[col_map_index])
            {
               marker_node[col_map_index] = 1;
               num_cols_offd_node++;
            }
            num_nonzeros_node++;
         }
      }
   }

   if (num_nonzeros_node == 0)
   {
      for (i = 0; i < num_rows + 1; i++)
      {
         offd_node_i[i] = 0;
      }
      jxf_CSRMatrixNumCols(offd_node) = 0;
      jxf_CSRMatrixI(offd_node) = offd_node_i;
      jxf_hpCSRMatrixOffdlevel(node_matrix) = offd_node;
   }
   else
   {
      offd_node_i[num_rows] = num_nonzeros_node;
      col_map_offd_node = jxf_CTAlloc(JXF_Int, num_cols_offd_node);

      counter_node = 0;
      for (i = 0; i < num_cols; i++)
      {
         if (marker_node[i])
         {
            col_map_offd_node[counter_node] = i;
            marker_node[i] = counter_node;
            counter_node++;
         }
      }

      jxf_CSRMatrixNumRows(offd_node) = num_rows;
      jxf_CSRMatrixNumNonzeros(offd_node) = num_nonzeros_node;
      jxf_CSRMatrixNumCols(offd_node) = num_cols_offd_node;
      jxf_CSRMatrixJ(offd_node) = NULL;
      jxf_CSRMatrixData(offd_node) = NULL;
      jxf_CSRMatrixOwnsData(offd_node) = 1;
      jxf_CSRMatrixNumRownnz(offd_node) = num_rows;
      jxf_CSRMatrixI(offd_node) = offd_node_i;

      jxf_CSRMatrixInitialize(offd_node);

      offd_node_data = jxf_CSRMatrixData(offd_node);
      offd_node_j = jxf_CSRMatrixJ(offd_node);

      j_node = 0;
      for (i = 0; i < num_rows; i++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
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
      jxf_CSRMatrixOwnsData(offd_node) = 1;
      jxf_hpCSRMatrixOffdlevel(node_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(node_matrix) = offd_node;
      jxf_hpCSRMatrixColMapOffdlevel(node_matrix) = jxf_CTAlloc(JXF_Int, num_cols_offd_node);
      jxf_hpCSRMatrixColMapOffdlevel(node_matrix) = col_map_offd_node;
   }
   jxf_TFree(marker_node);
   return node_matrix;
}

/*!
 * \fn JXF_Int jxf_hpMatrixLevelToPar
 * \date 2022/03/31
 */
jxf_ParCSRMatrix *
jxf_hpMatrixLevelToPar(jxf_hpCSRMatrix *hp_matrix, JXF_Int Level)
{
   jxf_ParCSRMatrix *par_matrix = jxf_CTAlloc(jxf_ParCSRMatrix, 1);
   jxf_ParCSRMatrixComm(par_matrix) = jxf_hpCSRMatrixComm(hp_matrix);
   jxf_ParCSRMatrixGlobalNumRows(par_matrix) = jxf_hpCSRMatrixGlobalNumRows(hp_matrix);
   jxf_ParCSRMatrixGlobalNumCols(par_matrix) = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   jxf_ParCSRMatrixFirstRowIndex(par_matrix) = jxf_hpCSRMatrixFirstRowIndex(hp_matrix);
   jxf_ParCSRMatrixFirstColDiag(par_matrix) = jxf_hpCSRMatrixFirstColDiag(hp_matrix);
   jxf_ParCSRMatrixLastRowIndex(par_matrix) = jxf_hpCSRMatrixLastRowIndex(hp_matrix);
   jxf_ParCSRMatrixLastColDiag(par_matrix) = jxf_hpCSRMatrixLastRowIndex(hp_matrix);
   jxf_ParCSRMatrixDiag(par_matrix) = jxf_hpCSRMatrixDiag(hp_matrix);
   if (Level == 0) // core
   {
      jxf_ParCSRMatrixOffd(par_matrix) = jxf_hpCSRMatrixOffdCore(hp_matrix);
      jxf_ParCSRMatrixColMapOffd(par_matrix) = jxf_hpCSRMatrixColMapOffdCore(hp_matrix);
      jxf_ParCSRMatrixCommPkg(par_matrix) = jxf_hpCSRMatrixCoreCommPkg(hp_matrix);
      jxf_ParCSRMatrixNumNonzeros(par_matrix) = jxf_hpCSRMatrixDNumNonzeros(hp_matrix);
   }
   else if (Level == 1) // cpu
   {
      jxf_ParCSRMatrixOffd(par_matrix) = jxf_hpCSRMatrixOffdCpu(hp_matrix);
      jxf_ParCSRMatrixColMapOffd(par_matrix) = jxf_hpCSRMatrixColMapOffdCpu(hp_matrix);
      jxf_ParCSRMatrixCommPkg(par_matrix) = jxf_hpCSRMatrixCpuCommPkg(hp_matrix);
      jxf_ParCSRMatrixNumNonzeros(par_matrix) = jxf_hpCSRMatrixDNumNonzeros(hp_matrix) + jxf_CSRMatrixNumNonzeros(jxf_hpCSRMatrixOffdCpu(hp_matrix));
   }
   else if (Level == 2) // node
   {
      jxf_ParCSRMatrixOffd(par_matrix) = jxf_hpCSRMatrixOffdNode(hp_matrix);
      jxf_ParCSRMatrixColMapOffd(par_matrix) = jxf_hpCSRMatrixColMapOffdNode(hp_matrix);
      jxf_ParCSRMatrixCommPkg(par_matrix) = jxf_hpCSRMatrixNodeCommPkg(hp_matrix);
      jxf_ParCSRMatrixNumNonzeros(par_matrix) = jxf_hpCSRMatrixDNumNonzeros(hp_matrix) + jxf_CSRMatrixNumNonzeros(jxf_hpCSRMatrixOffdNode(hp_matrix));
   }
   else
   {
      jxf_printf("Level must be 1 or 2\n");
      return NULL;
   }
   jxf_ParCSRMatrixDiagT(par_matrix) = jxf_hpCSRMatrixDiagT(hp_matrix);
   jxf_ParCSRMatrixOffdT(par_matrix) = NULL;                                    // 暂时先赋值为空
   jxf_ParCSRMatrixRowStarts(par_matrix) = jxf_hpCSRMatrixRowStarts(hp_matrix); // SPMV中不需要RowStarts
   jxf_ParCSRMatrixColStarts(par_matrix) = jxf_hpCSRMatrixColStarts(hp_matrix);
   jxf_ParCSRMatrixCommPkgT(par_matrix) = NULL; // 暂时先赋值为空
   jxf_ParCSRMatrixOwnsData(par_matrix) = jxf_hpCSRMatrixOwnsData(hp_matrix);
   jxf_ParCSRMatrixOwnsRowStarts(par_matrix) = jxf_hpCSRMatrixOwnsRowStarts(hp_matrix);
   jxf_ParCSRMatrixOwnsColStarts(par_matrix) = jxf_hpCSRMatrixOwnsColStarts(hp_matrix);
   jxf_ParCSRMatrixDNumNonzeros(par_matrix) = jxf_hpCSRMatrixDNumNonzeros(hp_matrix);
   jxf_ParCSRMatrixRowindices(par_matrix) = jxf_hpCSRMatrixRowindices(hp_matrix);
   jxf_ParCSRMatrixRowvalues(par_matrix) = jxf_hpCSRMatrixRowvalues(hp_matrix);
   jxf_ParCSRMatrixGetrowactive(par_matrix) = jxf_hpCSRMatrixGetrowactive(hp_matrix);
   jxf_ParCSRMatrixAssumedPartition(par_matrix) = jxf_hpCSRMatrixAssumedPartition(hp_matrix);
   return par_matrix;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixSetRowStartsOwner
 * \date 2009/07/10
 */
JXF_Int
jxf_hpCSRMatrixSetRowStartsOwner(jxf_hpCSRMatrix *matrix, JXF_Int owns_row_starts)
{
   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_hpCSRMatrixOwnsRowStarts(matrix) = owns_row_starts;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_hpCSRMatrixDestroyAssumedPartition
 * \date 2021/09/08
 */
JXF_Int
jxf_hpCSRMatrixDestroyAssumedPartition(jxf_hpCSRMatrix *matrix)
{
   return jxf_ParCSRMatrixDestroyAssumedPartition(jxf_hpCSRMatrixPar(matrix));
}

/*!
 * \fn JXF_Int jxf_hpCSRMatrixDestroy
 * \brief Destroy a hpCSR matrix.
 * \date 2021/09/08
 */
JXF_Int
jxf_hpCSRMatrixDestroy(jxf_hpCSRMatrix *matrix)
{
   if (matrix)
   {
      if (jxf_hpCSRMatrixPar(matrix))
      {
         jxf_ParCSRMatrixDestroy(jxf_hpCSRMatrixPar(matrix));
      }
      if (jxf_hpCSRMatrixNode(matrix))
      {
         jxf_hpCSRMatrixLevelDestroy(jxf_hpCSRMatrixNode(matrix));
      }
      if (jxf_hpCSRMatrixCpu(matrix))
      {
         jxf_hpCSRMatrixLevelDestroy(jxf_hpCSRMatrixCpu(matrix));
      }
      if (jxf_hpCSRMatrixCore(matrix))
      {
         jxf_hpCSRMatrixLevelDestroy(jxf_hpCSRMatrixCore(matrix));
      }
      jxf_TFree(matrix);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_hpCSRhardwareDestroy()
{
   if (jxf_hp_hardware)
   {

      if (jxf_hpCoreComm(jxf_hp_hardware) != MPI_COMM_NULL)
      {
         jxf_MPI_Comm_free(&jxf_hpCoreComm(jxf_hp_hardware));
         jxf_hpCoreComm(jxf_hp_hardware) = MPI_COMM_NULL;
      }

      if (jxf_hpIntraCpuComm(jxf_hp_hardware) != MPI_COMM_NULL)
      {
         jxf_MPI_Comm_free(&jxf_hpIntraCpuComm(jxf_hp_hardware));
         jxf_hpIntraCpuComm(jxf_hp_hardware) = MPI_COMM_NULL;
      }

      if (jxf_hpNodeComm(jxf_hp_hardware) != MPI_COMM_NULL)
      {
         jxf_MPI_Comm_free(&jxf_hpNodeComm(jxf_hp_hardware));
         jxf_hpNodeComm(jxf_hp_hardware) = MPI_COMM_NULL;
      }
      jxf_TFree(jxf_hp_hardware);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_hpCSRMatrixLevelDestroy(jxf_hpCSRMatrixLevel *matrix)
{
   if (matrix)
   {
      if (jxf_hpCSRMatrixOffdlevel(matrix))
      {
         jxf_CSRMatrixDestroy(jxf_hpCSRMatrixOffdlevel(matrix));
      }
      if (jxf_hpCSRMatrixColMapOffdlevel(matrix))
      {
         jxf_TFree(jxf_hpCSRMatrixColMapOffdlevel(matrix));
      }
      if (jxf_hpCSRMatrixlevelCommPkg(matrix))
      {
         jxf_MatvecCommPkgDestroy(jxf_hpCSRMatrixlevelCommPkg(matrix));
      }
      jxf_TFree(matrix);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_hpCSRMatrixPrint(jxf_hpCSRMatrix *hp_matrix, const char *file_name)
{
   MPI_Comm comm;
   JXF_Int global_num_rows;
   JXF_Int global_num_cols;
   JXF_Int num_cols_cpu = 0;
   JXF_Int num_cols_node = 0;
   JXF_Int *col_map_offd;
#ifndef JXF_NO_GLOBAL_PARTITION
   JXF_Int *row_starts;
   JXF_Int *col_starts;
#endif
   JXF_Int my_id, i, num_procs;
   char new_file_offd_cpu[80], new_file_offd_node[80], new_file_o[80], new_file_info[80];
   FILE *fp;
   JXF_Int num_cols_offd = 0;
#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int row_s, row_e, col_s, col_e;
#endif
   if (!hp_matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   comm = jxf_hpCSRMatrixComm(hp_matrix);
   global_num_rows = jxf_hpCSRMatrixGlobalNumRows(hp_matrix);
   global_num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   col_map_offd = jxf_hpCSRMatrixColMapOffd(hp_matrix);

#ifndef JXF_NO_GLOBAL_PARTITION
   row_starts = jxf_hpCSRMatrixRowStarts(hp_matrix);
   col_starts = jxf_hpCSRMatrixColStarts(hp_matrix);
#endif

   if (jxf_hpCSRMatrixOffd(hp_matrix))
   {
      num_cols_offd = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffd(hp_matrix));
   }
   if (jxf_hpCSRMatrixOffdCpu(hp_matrix))
   {
      num_cols_cpu = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdCpu(hp_matrix));
   }
   if (jxf_hpCSRMatrixOffdNode(hp_matrix))
   {
      num_cols_node = jxf_CSRMatrixNumCols(jxf_hpCSRMatrixOffdCpu(hp_matrix));
   }

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   jxf_sprintf(new_file_offd_cpu, "%s.offd_cpu.%d", file_name, my_id);
   jxf_sprintf(new_file_offd_node, "%s.offd_node.%d", file_name, my_id);
   jxf_sprintf(new_file_o, "%s.O.%d", file_name, my_id);
   jxf_sprintf(new_file_info, "%s.INFO.%d", file_name, my_id);
   if (num_cols_cpu != 0)
      jxf_CSRMatrixPrint(jxf_hpCSRMatrixOffdCpu(hp_matrix), new_file_offd_cpu);
   if (num_cols_node != 0)
      jxf_CSRMatrixPrint(jxf_hpCSRMatrixOffdNode(hp_matrix), new_file_offd_node);
   if (num_cols_offd != 0)
   {
      jxf_CSRMatrixPrint(jxf_hpCSRMatrixOffd(hp_matrix), new_file_o);
   }

   fp = fopen(new_file_info, "w");
   jxf_fprintf(fp, "%d\n", global_num_rows);
   jxf_fprintf(fp, "%d\n", global_num_cols);
   jxf_fprintf(fp, "%d\n", num_cols_offd);

#ifdef JXF_NO_GLOBAL_PARTITION
   row_s = jxf_hpCSRMatrixFirstRowIndex(hp_matrix);
   row_e = jxf_hpCSRMatrixLastRowIndex(hp_matrix);
   col_s = jxf_hpCSRMatrixFirstColDiag(hp_matrix);
   col_e = jxf_hpCSRMatrixLastColDiag(hp_matrix);
   /* add 1 to the ends because this is a starts partition */
   jxf_fprintf(fp, "%d %d %d %d\n", row_s, row_e + 1, col_s, col_e + 1);
#else
   jxf_fprintf(fp, "row_start col_start: \n");
   for (i = 0; i < num_procs; i++)
   {
      jxf_fprintf(fp, "%d %d\n", row_starts[i], col_starts[i]);
   }
#endif

   jxf_fprintf(fp, "col_map_offd: \n");
   for (i = 0; i < num_cols_offd; i++)
   {
      jxf_fprintf(fp, "%d\n", col_map_offd[i]);
   }

   fclose(fp);

   return jxf_error_flag;
}

/*!
 * \fn jxf_hpCSRMatrixGenerateOffdNodeOutside
 * \brief Generate outside node level matrix
 * \author mrz
 * \date 2021/08/31
 */
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdNodeOutside(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware)
{
   // test need
   MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   JXF_Int myid;
   jxf_MPI_Comm_rank(comm, &myid);

   jxf_hpCSRMatrixLevel *node_outside_matrix = jxf_hpInitMatrixLevel();
   JXF_Int i, j, j_node_outside, col_index, col_map_index;
   jxf_CSRMatrix *offd = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(offd);
   // JXF_Int  num_cols = jxf_CSRMatrixNumCols(offd);
   JXF_Int num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Int *offd_col_map = jxf_hpCSRMatrixColMapOffd(hp_matrix);

   jxf_CSRMatrix *offd_node_outside;
   JXF_Int *col_map_offd_node_outside = NULL;
   JXF_Real *offd_node_outside_data;
   JXF_Int *offd_node_outside_i;
   JXF_Int *offd_node_outside_j;

   JXF_Int *marker_node_outside;
   JXF_Int num_cols_offd_node_outside = 0;
   JXF_Int first_elmt = offd_i[0];
   JXF_Int num_nonzeros_node_outside;
   JXF_Int counter_node_outside;

   JXF_Int *row_start = jxf_hpCSRMatrixRowStarts(hp_matrix);
   JXF_Int node_colstartnum = row_start[jxf_hpProcNodeStart(hardware)];
   JXF_Int node_colendnum = row_start[jxf_hpProcNodeEnd(hardware) + 1] - 1;
   // //TODO printf待注释
   // if(myid == 0)
   // {
   //    for(i = 0; i < 3; i++)
   //    {
   //       jxf_printf("row_start[%d] is %d\n",i, row_start[i]);
   //    }
   //    jxf_printf("node_colstartnum is %d, node_colendnum is %d\n",
   //    node_colstartnum, node_colendnum);
   // }

   // generate offd
   offd_node_outside_i = jxf_CTAlloc(JXF_Int, num_rows + 1);
   offd_node_outside = jxf_CTAlloc(jxf_CSRMatrix, 1);
   marker_node_outside = jxf_CTAlloc(JXF_Int, num_cols);

   for (i = 0; i < num_cols; i++)
   {
      marker_node_outside[i] = 0;
   }

   num_nonzeros_node_outside = 0;

   for (i = 0; i < num_rows; i++)
   {
      offd_node_outside_i[i] = num_nonzeros_node_outside;

      for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
      {
         // col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         // if (offd_col_map[col_index] >= cpu_colstartnum && offd_col_map[col_index] <= cpu_colendnum)
         // {
         //    continue;
         // }
         // else
         if (col_map_index < node_colstartnum || col_map_index > node_colendnum)
         {
            // if (!marker_node[col_index])
            // {
            //    marker_node[col_index] = 1;
            //    num_cols_offd_node ++;
            // }
            if (!marker_node_outside[col_map_index])
            {
               marker_node_outside[col_map_index] = 1;
               num_cols_offd_node_outside++;
            }
            num_nonzeros_node_outside++;
         }
      }
   }

   if (num_nonzeros_node_outside == 0)
   {
      for (i = 0; i < num_rows + 1; i++)
      {
         offd_node_outside_i[i] = 0;
      }
      jxf_CSRMatrixNumCols(offd_node_outside) = 0;
      jxf_CSRMatrixI(offd_node_outside) = offd_node_outside_i;
      jxf_hpCSRMatrixOffdlevel(node_outside_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(node_outside_matrix) = offd_node_outside;
   }
   else
   {
      offd_node_outside_i[num_rows] = num_nonzeros_node_outside;
      col_map_offd_node_outside = jxf_CTAlloc(JXF_Int, num_cols_offd_node_outside);

      counter_node_outside = 0;
      for (i = 0; i < num_cols; i++)
      {
         if (marker_node_outside[i])
         {
            col_map_offd_node_outside[counter_node_outside] = i;
            marker_node_outside[i] = counter_node_outside;
            counter_node_outside++;
         }
      }

      jxf_CSRMatrixNumRows(offd_node_outside) = num_rows;
      jxf_CSRMatrixNumNonzeros(offd_node_outside) = num_nonzeros_node_outside;
      jxf_CSRMatrixNumCols(offd_node_outside) = num_cols_offd_node_outside;
      jxf_CSRMatrixJ(offd_node_outside) = NULL;
      jxf_CSRMatrixData(offd_node_outside) = NULL;
      jxf_CSRMatrixOwnsData(offd_node_outside) = 1;
      jxf_CSRMatrixNumRownnz(offd_node_outside) = num_rows;
      jxf_CSRMatrixI(offd_node_outside) = offd_node_outside_i;

      jxf_CSRMatrixInitialize(offd_node_outside);

      offd_node_outside_data = jxf_CSRMatrixData(offd_node_outside);
      offd_node_outside_j = jxf_CSRMatrixJ(offd_node_outside);

      j_node_outside = 0;
      for (i = 0; i < num_rows; i++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
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
      jxf_CSRMatrixOwnsData(offd_node_outside) = 1;
      jxf_hpCSRMatrixOffdlevel(node_outside_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(node_outside_matrix) = offd_node_outside;
      jxf_hpCSRMatrixColMapOffdlevel(node_outside_matrix) = jxf_CTAlloc(JXF_Int, num_cols_offd_node_outside);
      jxf_hpCSRMatrixColMapOffdlevel(node_outside_matrix) = col_map_offd_node_outside;
   }

   jxf_TFree(marker_node_outside);
   return node_outside_matrix;
}

/*!
 * \fn jxf_hpCSRMatrixGenerateOffdNode
 * \brief Generate node level matrix
 * \author mrz
 * \date 2021/08/31
 */
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdCpuOutside(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware)
{
   // test need
   // MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   // JXF_Int myid;
   // jxf_MPI_Comm_rank(comm, &myid);

   jxf_hpCSRMatrixLevel *cpu_outside_matrix = jxf_hpInitMatrixLevel();
   JXF_Int i, j, j_cpu_outside, col_index, col_map_index;
   jxf_CSRMatrix *offd = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(offd);
   // JXF_Int  num_cols = jxf_CSRMatrixNumCols(offd);
   JXF_Int num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_matrix);
   JXF_Real *offd_data = jxf_CSRMatrixData(offd);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *offd_j = jxf_CSRMatrixJ(offd);
   JXF_Int *offd_col_map = jxf_hpCSRMatrixColMapOffd(hp_matrix);

   jxf_CSRMatrix *offd_cpu_outside;
   JXF_Int *col_map_offd_cpu_outside = NULL;
   JXF_Real *offd_cpu_outside_data;
   JXF_Int *offd_cpu_outside_i;
   JXF_Int *offd_cpu_outside_j;

   JXF_Int *marker_cpu_outside;
   JXF_Int num_cols_offd_cpu_outside = 0;
   JXF_Int first_elmt = offd_i[0];
   JXF_Int num_nonzeros_cpu_outside;
   JXF_Int counter_cpu_outside;

   JXF_Int *row_start = jxf_hpCSRMatrixRowStarts(hp_matrix);
   JXF_Int node_colstartnum = row_start[jxf_hpProcNodeStart(hardware)];
   JXF_Int node_colendnum = row_start[jxf_hpProcNodeEnd(hardware) + 1] - 1;
   JXF_Int cpu_colstartnum = row_start[jxf_hpProcCpuStart(hardware)];
   JXF_Int cpu_colendnum = row_start[jxf_hpProcCpuEnd(hardware) + 1] - 1;
   offd_cpu_outside_i = jxf_CTAlloc(JXF_Int, num_rows + 1);
   offd_cpu_outside = jxf_CTAlloc(jxf_CSRMatrix, 1);
   marker_cpu_outside = jxf_CTAlloc(JXF_Int, num_cols);

   for (i = 0; i < num_cols; i++)
   {
      marker_cpu_outside[i] = 0;
   }

   num_nonzeros_cpu_outside = 0;

   for (i = 0; i < num_rows; i++)
   {
      offd_cpu_outside_i[i] = num_nonzeros_cpu_outside;

      for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
      {
         // col_index = offd_j[j];
         col_map_index = offd_col_map[offd_j[j]];
         if (col_map_index >= cpu_colstartnum && col_map_index <= cpu_colendnum)
         {
            continue;
         }
         else if (col_map_index >= node_colstartnum && col_map_index <= node_colendnum)
         {
            // if (!marker_node[col_index])
            // {
            //    marker_node[col_index] = 1;
            //    num_cols_offd_node ++;
            // }
            if (!marker_cpu_outside[col_map_index])
            {
               marker_cpu_outside[col_map_index] = 1;
               num_cols_offd_cpu_outside++;
            }
            num_nonzeros_cpu_outside++;
         }
      }
   }

   if (num_nonzeros_cpu_outside == 0)
   {
      for (i = 0; i < num_rows + 1; i++)
      {
         offd_cpu_outside_i[i] = 0;
      }
      jxf_CSRMatrixNumCols(offd_cpu_outside) = 0;
      jxf_CSRMatrixI(offd_cpu_outside) = offd_cpu_outside_i;
      jxf_hpCSRMatrixOffdlevel(cpu_outside_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(cpu_outside_matrix) = offd_cpu_outside;
   }
   else
   {
      offd_cpu_outside_i[num_rows] = num_nonzeros_cpu_outside;
      col_map_offd_cpu_outside = jxf_CTAlloc(JXF_Int, num_cols_offd_cpu_outside);

      counter_cpu_outside = 0;
      for (i = 0; i < num_cols; i++)
      {
         if (marker_cpu_outside[i])
         {
            col_map_offd_cpu_outside[counter_cpu_outside] = i;
            marker_cpu_outside[i] = counter_cpu_outside;
            counter_cpu_outside++;
         }
      }

      jxf_CSRMatrixNumRows(offd_cpu_outside) = num_rows;
      jxf_CSRMatrixNumNonzeros(offd_cpu_outside) = num_nonzeros_cpu_outside;
      jxf_CSRMatrixNumCols(offd_cpu_outside) = num_cols_offd_cpu_outside;
      jxf_CSRMatrixJ(offd_cpu_outside) = NULL;
      jxf_CSRMatrixData(offd_cpu_outside) = NULL;
      jxf_CSRMatrixOwnsData(offd_cpu_outside) = 1;
      jxf_CSRMatrixNumRownnz(offd_cpu_outside) = num_rows;
      jxf_CSRMatrixI(offd_cpu_outside) = offd_cpu_outside_i;

      jxf_CSRMatrixInitialize(offd_cpu_outside);

      offd_cpu_outside_data = jxf_CSRMatrixData(offd_cpu_outside);
      offd_cpu_outside_j = jxf_CSRMatrixJ(offd_cpu_outside);

      j_cpu_outside = 0;
      for (i = 0; i < num_rows; i++)
      {
         for (j = offd_i[i] - first_elmt; j < offd_i[i + 1] - first_elmt; j++)
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
      jxf_CSRMatrixOwnsData(offd_cpu_outside) = 1;
      jxf_hpCSRMatrixOffdlevel(cpu_outside_matrix) = jxf_CTAlloc(jxf_CSRMatrix, 1);
      jxf_hpCSRMatrixOffdlevel(cpu_outside_matrix) = offd_cpu_outside;
      jxf_hpCSRMatrixColMapOffdlevel(cpu_outside_matrix) = jxf_CTAlloc(JXF_Int, num_cols_offd_cpu_outside);
      jxf_hpCSRMatrixColMapOffdlevel(cpu_outside_matrix) = col_map_offd_cpu_outside;
   }

   jxf_TFree(marker_cpu_outside);
   return cpu_outside_matrix;
}
