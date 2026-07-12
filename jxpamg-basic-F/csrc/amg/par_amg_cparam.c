//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_amg_cparam.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGCoarseParms
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGCoarseParms( MPI_Comm   comm,
                    JXF_Int        local_num_variables,
                    JXF_Int        num_functions,
                    JXF_Int       *dof_func,
                    JXF_Int       *CF_marker, 
                    JXF_Int      **coarse_dof_func_ptr, 
                    JXF_Int      **coarse_pnts_global_ptr) 
{
   JXF_Int  i;
   JXF_Int  ierr = 0;
   JXF_Int	num_procs;
   JXF_Int  local_coarse_size = 0;

   JXF_Int *coarse_dof_func;
   JXF_Int *coarse_pnts_global;

   jxf_MPI_Comm_size(comm, &num_procs);

   for (i = 0; i < local_num_variables; i ++)
   {
      if (CF_marker[i] == 1) 
      {
         local_coarse_size++;
      }
   }
   
   if (num_functions > 1)
   {
      coarse_dof_func = jxf_CTAlloc(JXF_Int, local_coarse_size);
      local_coarse_size = 0;
      for (i = 0; i < local_num_variables; i ++)
      {
         if (CF_marker[i] == 1)
         {
            coarse_dof_func[local_coarse_size++] = dof_func[i];
         }
      }
      *coarse_dof_func_ptr = coarse_dof_func;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int scan_recv;
   coarse_pnts_global = jxf_CTAlloc(JXF_Int, 2);
   jxf_MPI_Scan(&local_coarse_size, &scan_recv, 1, JXF_MPI_INT, MPI_SUM, comm);
   /* first point in my range */ 
   coarse_pnts_global[0] = scan_recv - local_coarse_size;
   /* first point in next proc's range */
   coarse_pnts_global[1] = scan_recv;  
#else
   coarse_pnts_global = jxf_CTAlloc(JXF_Int, num_procs + 1);
   jxf_MPI_Allgather(&local_coarse_size, 1, JXF_MPI_INT, &coarse_pnts_global[1], 1, JXF_MPI_INT, comm);
   for (i = 2; i < num_procs + 1; i ++)
   {
      coarse_pnts_global[i] += coarse_pnts_global[i-1];
   }
#endif

   *coarse_pnts_global_ptr = coarse_pnts_global;

   return (ierr);
}
