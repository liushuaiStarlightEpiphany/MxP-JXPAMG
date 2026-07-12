//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  grid.c
 *
 * Created by Yue Xiaoqiang 2014/08/03
 * Xiangtan University
 *
 */

#include "jxf_ilu.h"

/*!
 * \fn JXF_Int JXF_GridPartitionDataCreate
 * \author Yue Xiaoqiang
 * \date 2014/10/20
 */
JXF_Int
JXF_GridPartitionDataCreate( JXF_Solver *solver, MPI_Comm comm, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy )
{
  *solver = (JXF_Solver)jxf_GridPartitionDataInitialize(comm, nx, ny, npx, npy);
   if (!solver)
   {
      jxf_error_in_arg(1);
   }
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_GridPartitionDataSetEachSides4Comm
 * \author Yue Xiaoqiang
 * \date 2014/08/03
 */
JXF_Int
JXF_GridPartitionDataSetEachSides4Comm( JXF_Solver solver )
{
   return( jxf_GridPartitionDataSetEachSides4Comm( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_GridPartitionDataDestroy
 * \author Yue Xiaoqiang
 * \date 2014/08/03
 */
JXF_Int
JXF_GridPartitionDataDestroy( JXF_Solver solver )
{
   return( jxf_GridPartitionDataFinalize( (void *) solver ) );
}

/*!
 * \fn void *jxf_GridPartitionDataInitialize
 * \brief Initialize
 * \author Yue Xiaoqiang
 * \date 2014/10/20
 */
void *
jxf_GridPartitionDataInitialize( MPI_Comm comm, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy )
{
    jxf_GridPartitionData *grid_data = jxf_CTAlloc(jxf_GridPartitionData, 1);
    JXF_Int *nx_part = NULL;
    JXF_Int *ny_part = NULL;
    JXF_Int num_procs, my_id, iy, ix, cnt = 0;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    grid_data->comm = comm;
    grid_data->xlo_array = jxf_CTAlloc(JXF_Int, num_procs);
    grid_data->xup_array = jxf_CTAlloc(JXF_Int, num_procs);
    grid_data->ylo_array = jxf_CTAlloc(JXF_Int, num_procs);
    grid_data->yup_array = jxf_CTAlloc(JXF_Int, num_procs);
    jxf_GeneratePartitioning(nx, npx, &nx_part);
    jxf_GeneratePartitioning(ny, npy, &ny_part);
    if (0) // Processors in x-direction first
    {
        for (iy = 0; iy < npy; iy ++)
        {
            for (ix = 0; ix < npx; ix ++)
            {
                grid_data->xlo_array[cnt] = nx_part[ix];
                grid_data->xup_array[cnt] = nx_part[ix+1];
                grid_data->ylo_array[cnt] = ny_part[iy];
                grid_data->yup_array[cnt] = ny_part[iy+1];
                cnt ++;
            }
        }
    }
    else // Processors in y-direction first
    {
        for (ix = 0; ix < npx; ix ++)
        {
            for (iy = 0; iy < npy; iy ++)
            {
                grid_data->xlo_array[cnt] = nx_part[ix];
                grid_data->xup_array[cnt] = nx_part[ix+1];
                grid_data->ylo_array[cnt] = ny_part[iy];
                grid_data->yup_array[cnt] = ny_part[iy+1];
                cnt ++;
            }
        }
    }
    jxf_TFree(nx_part);
    jxf_TFree(ny_part);
    grid_data->x_lower_idx = grid_data->xlo_array[my_id];
    grid_data->x_upper_idx = grid_data->xup_array[my_id];
    grid_data->y_lower_idx = grid_data->ylo_array[my_id];
    grid_data->y_upper_idx = grid_data->yup_array[my_id];
    grid_data->x_part_len = grid_data->xup_array[my_id] - grid_data->xlo_array[my_id];
    grid_data->y_part_len = grid_data->yup_array[my_id] - grid_data->ylo_array[my_id];
    
    return (void *)grid_data;
}

/*!
 * \fn JXF_Int jxf_GridPartitionDataSetEachSides4Comm
 * \brief Set sideprocs and sideprcfgstln
 * \note not for the case of a street intersection
 * \author Yue Xiaoqiang
 * \date 2014/08/03
 */
JXF_Int
jxf_GridPartitionDataSetEachSides4Comm( void *grid_vdata )
{
    jxf_GridPartitionData *grid_data = grid_vdata;
    JXF_Int x_lower_idx = grid_data->x_lower_idx;
    JXF_Int x_upper_idx = grid_data->x_upper_idx;
    JXF_Int y_lower_idx = grid_data->y_lower_idx;
    JXF_Int y_upper_idx = grid_data->y_upper_idx;
    JXF_Int *tmpsidepro = NULL;
    JXF_Int myid, nprocs, i, x_cnt, y_cnt, pid;
    JXF_Int num_sideproc, num_smallside, num_largeside;
    JXF_Int num_small_dlr, num_large_ulr, num_upside, num_rightside;
    
    jxf_MPI_Comm_size(grid_data->comm, &nprocs);
    jxf_MPI_Comm_rank(grid_data->comm, &myid);
    tmpsidepro = jxf_CTAlloc(JXF_Int, nprocs);
    // store all the side processors, Yue Xiaoqiang, 2014/08/03
    num_sideproc = 0;
    for (i = 0; i < nprocs; i ++)
    {
        if (i == myid)
        {
            continue;
        }
        x_cnt = 0;
        if ((grid_data->xlo_array[i] >= x_lower_idx) && (grid_data->xlo_array[i] <= x_upper_idx))
        {
            x_cnt ++;
        }
        if ((grid_data->xup_array[i] >= x_lower_idx) && (grid_data->xup_array[i] <= x_upper_idx))
        {
            x_cnt ++;
        }
        if (x_cnt < 1) // x_cnt = 0
        {
            // symmetric case, Yue Xiaoqiang, 2014/08/03
            if ((x_lower_idx >= grid_data->xlo_array[i]) && (x_lower_idx <= grid_data->xup_array[i]))
            {
                x_cnt ++;
            }
            if ((x_upper_idx >= grid_data->xlo_array[i]) && (x_upper_idx <= grid_data->xup_array[i]))
            {
                x_cnt ++;
            }
            if (x_cnt < 1)
            {
                continue;
            }
        }
        y_cnt = 0;
        if ((grid_data->ylo_array[i] >= y_lower_idx) && (grid_data->ylo_array[i] <= y_upper_idx))
        {
            y_cnt ++;
        }
        if ((grid_data->yup_array[i] >= y_lower_idx) && (grid_data->yup_array[i] <= y_upper_idx))
        {
            y_cnt ++;
        }
        if (y_cnt > 0)
        {
            tmpsidepro[num_sideproc++] = i;
            continue;
        }
        // symmetric case, Yue Xiaoqiang, 2014/08/03
        // y_cnt = 0
        if ((y_lower_idx >= grid_data->ylo_array[i]) && (y_lower_idx <= grid_data->yup_array[i]))
        {
            y_cnt ++;
        }
        if ((y_upper_idx >= grid_data->ylo_array[i]) && (y_upper_idx <= grid_data->yup_array[i]))
        {
            y_cnt ++;
        }
        if (y_cnt > 0)
        {
            tmpsidepro[num_sideproc++] = i;
        }
    }
    grid_data->num_sideproc = num_sideproc;
    num_smallside = 0;
    num_largeside = 0;
    num_small_dlr = 0;
    num_large_ulr = 0;
    num_upside = 0;
    num_rightside = 0;
    for (i = 0; i < num_sideproc; i ++)
    {
        pid = tmpsidepro[i];
        // determinations to deal with possible cases, Yue Xiaoqiang, 2014/08/03
        if ((grid_data->xlo_array[pid] == x_lower_idx) && (grid_data->xup_array[pid] == x_upper_idx))
        {
            if (y_lower_idx < grid_data->ylo_array[pid])
            {
                num_largeside ++; // UP
                num_upside ++;
            }
            else // won't be equal
            {
                num_smallside ++; // DOWN
            }
        }
        else if ((grid_data->ylo_array[pid] == y_lower_idx) && (grid_data->yup_array[pid] == y_upper_idx))
        {
            if (x_lower_idx < grid_data->xlo_array[pid])
            {
                num_largeside ++; // RIGHT
                num_rightside ++;
            }
            else // won't be equal
            {
                num_smallside ++; // LEFT
            }
        }
        else if (grid_data->xlo_array[pid] == x_upper_idx) // case 3
        {
            if (grid_data->ylo_array[pid] == y_upper_idx)
            {
                num_large_ulr ++; // UP - RIGHT
            }
            else if (grid_data->yup_array[pid] == y_lower_idx)
            {
                num_small_dlr ++; // DOWN - RIGHT
            }
            else
            {
                num_largeside ++; // RIGHT
                num_rightside ++;
            }
        }
        else if (grid_data->xup_array[pid] == x_lower_idx) // case 4
        {
            if (grid_data->ylo_array[pid] == y_upper_idx)
            {
                num_large_ulr ++; // UP - LEFT
            }
            else if (grid_data->yup_array[pid] == y_lower_idx)
            {
                num_small_dlr ++; // DOWN - LEFT
            }
            else
            {
                num_smallside ++; // LEFT
            }
        }
        else if (grid_data->ylo_array[pid] == y_upper_idx) // case 5
        {
            if (grid_data->xlo_array[pid] == x_upper_idx)
            {
                num_large_ulr ++; // UP - RIGHT
            }
            else if (grid_data->xup_array[pid] == x_lower_idx)
            {
                num_large_ulr ++; // UP - LEFT
            }
            else
            {
                num_largeside ++; // UP
                num_upside ++;
            }
        }
        else if (grid_data->yup_array[pid] == y_lower_idx) // case 6
        {
            if (grid_data->xlo_array[pid] == x_upper_idx)
            {
                num_small_dlr ++; // DOWN - RIGHT
            }
            else if (grid_data->xup_array[pid] == x_lower_idx)
            {
                num_small_dlr ++; // DOWN - LEFT
            }
            else
            {
                num_smallside ++; // DOWN
            }
        }
    }
    if (num_upside > 0)
    {
        if (num_rightside > 0)
        {
            grid_data->part_type = 1; // UP and RIGHT both exist
        }
        else
        {
            grid_data->part_type = 2; // Only UP exists
        }
    }
    else
    {
        if (num_rightside > 0)
        {
            grid_data->part_type = 3; // Only RIGHT exists
        }
        else
        {
            grid_data->part_type = 4; // None of UP or RIGHT
        }
    }
    grid_data->num_smallside = num_smallside;
    grid_data->num_largeside = num_largeside;
    grid_data->num_smallcross = num_small_dlr;
    grid_data->num_largecross = num_large_ulr;
    grid_data->num_nocrossside = num_smallside + num_largeside;
    //jxf_assert(num_smallside+num_largeside+num_small_dlr+num_large_ulr == num_sideproc);
    num_large_ulr = num_smallside + num_largeside + num_small_dlr;
    num_small_dlr = num_smallside + num_largeside;
    num_largeside = num_smallside;
    num_smallside = 0;
    grid_data->sideprocs = jxf_CTAlloc(JXF_Int, num_sideproc);
    grid_data->sideprcpos = jxf_CTAlloc(JXF_Int, num_sideproc);
    grid_data->sideprcxsrt = jxf_CTAlloc(JXF_Int, num_sideproc);
    grid_data->sideprcysrt = jxf_CTAlloc(JXF_Int, num_sideproc);
    grid_data->sideprclength = jxf_CTAlloc(JXF_Int, num_sideproc);
    for (i = 0; i < num_sideproc; i ++)
    {
        pid = tmpsidepro[i];
        // determinations to deal with possible cases, Yue Xiaoqiang, 2014/08/03
        if ((grid_data->xlo_array[pid] == x_lower_idx) && (grid_data->xup_array[pid] == x_upper_idx))
        {
            if (y_lower_idx < grid_data->ylo_array[pid])
            {
                grid_data->sideprocs[num_largeside] = pid;
                grid_data->sideprcpos[num_largeside] = 1;
                grid_data->sideprcxsrt[num_largeside] = x_lower_idx;
                grid_data->sideprcysrt[num_largeside] = y_upper_idx - 1;
                grid_data->sideprclength[num_largeside] = x_upper_idx - x_lower_idx;
                num_largeside ++;
            }
            else // won't be equal
            {
                grid_data->sideprocs[num_smallside] = pid;
                grid_data->sideprcpos[num_smallside] = 2;
                grid_data->sideprcxsrt[num_smallside] = x_lower_idx;
                grid_data->sideprcysrt[num_smallside] = y_lower_idx;
                grid_data->sideprclength[num_smallside] = x_upper_idx - x_lower_idx;
                num_smallside ++;
            }
        }
        else if ((grid_data->ylo_array[pid] == y_lower_idx) && (grid_data->yup_array[pid] == y_upper_idx))
        {
            if (x_lower_idx < grid_data->xlo_array[pid])
            {
                grid_data->sideprocs[num_largeside] = pid;
                grid_data->sideprcpos[num_largeside] = 3;
                grid_data->sideprcxsrt[num_largeside] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_largeside] = y_lower_idx;
                grid_data->sideprclength[num_largeside] = y_upper_idx - y_lower_idx;
                num_largeside ++;
            }
            else // won't be equal
            {
                grid_data->sideprocs[num_smallside] = pid;
                grid_data->sideprcpos[num_smallside] = 4;
                grid_data->sideprcxsrt[num_smallside] = x_lower_idx;
                grid_data->sideprcysrt[num_smallside] = y_lower_idx;
                grid_data->sideprclength[num_smallside] = y_upper_idx - y_lower_idx;
                num_smallside ++;
            }
        }
        else if (grid_data->xlo_array[pid] == x_upper_idx) // case 3
        {
            if (grid_data->ylo_array[pid] == y_upper_idx)
            {
                grid_data->sideprocs[num_large_ulr] = pid;
                grid_data->sideprcpos[num_large_ulr] = 6;
                grid_data->sideprcxsrt[num_large_ulr] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_large_ulr] = y_upper_idx - 1;
                grid_data->sideprclength[num_large_ulr] = 1;
                num_large_ulr ++;
            }
            else if (grid_data->yup_array[pid] == y_lower_idx)
            {
                grid_data->sideprocs[num_small_dlr] = pid;
                grid_data->sideprcpos[num_small_dlr] = 8;
                grid_data->sideprcxsrt[num_small_dlr] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_small_dlr] = y_lower_idx;
                grid_data->sideprclength[num_small_dlr] = 1;
                num_small_dlr ++;
            }
            else
            {
                grid_data->sideprocs[num_largeside] = pid;
                grid_data->sideprcpos[num_largeside] = 3;
                grid_data->sideprcxsrt[num_largeside] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_largeside] = jxf_max(y_lower_idx, grid_data->ylo_array[pid]);
                grid_data->sideprclength[num_largeside] = jxf_min(y_upper_idx,
                                   grid_data->yup_array[pid]) - jxf_max(y_lower_idx, grid_data->ylo_array[pid]);
                num_largeside ++;
            }
        }
        else if (grid_data->xup_array[pid] == x_lower_idx) // case 4
        {
            if (grid_data->ylo_array[pid] == y_upper_idx)
            {
                grid_data->sideprocs[num_large_ulr] = pid;
                grid_data->sideprcpos[num_large_ulr] = 5;
                grid_data->sideprcxsrt[num_large_ulr] = x_lower_idx;
                grid_data->sideprcysrt[num_large_ulr] = y_upper_idx - 1;
                grid_data->sideprclength[num_large_ulr] = 1;
                num_large_ulr ++;
            }
            else if (grid_data->yup_array[pid] == y_lower_idx)
            {
                grid_data->sideprocs[num_small_dlr] = pid;
                grid_data->sideprcpos[num_small_dlr] = 7;
                grid_data->sideprcxsrt[num_small_dlr] = x_lower_idx;
                grid_data->sideprcysrt[num_small_dlr] = y_lower_idx;
                grid_data->sideprclength[num_small_dlr] = 1;
                num_small_dlr ++;
            }
            else
            {
                grid_data->sideprocs[num_smallside] = pid;
                grid_data->sideprcpos[num_smallside] = 4;
                grid_data->sideprcxsrt[num_smallside] = x_lower_idx;
                grid_data->sideprcysrt[num_smallside] = jxf_max(y_lower_idx, grid_data->ylo_array[pid]);
                grid_data->sideprclength[num_smallside] = jxf_min(y_upper_idx,
                                   grid_data->yup_array[pid]) - jxf_max(y_lower_idx, grid_data->ylo_array[pid]);
                num_smallside ++;
            }
        }
        else if (grid_data->ylo_array[pid] == y_upper_idx) // case 5
        {
            if (grid_data->xlo_array[pid] == x_upper_idx)
            {
                grid_data->sideprocs[num_large_ulr] = pid;
                grid_data->sideprcpos[num_large_ulr] = 6;
                grid_data->sideprcxsrt[num_large_ulr] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_large_ulr] = y_upper_idx - 1;
                grid_data->sideprclength[num_large_ulr] = 1;
                num_large_ulr ++;
            }
            else if (grid_data->xup_array[pid] == x_lower_idx)
            {
                grid_data->sideprocs[num_large_ulr] = pid;
                grid_data->sideprcpos[num_large_ulr] = 5;
                grid_data->sideprcxsrt[num_large_ulr] = x_lower_idx;
                grid_data->sideprcysrt[num_large_ulr] = y_upper_idx - 1;
                grid_data->sideprclength[num_large_ulr] = 1;
                num_large_ulr ++;
            }
            else
            {
                grid_data->sideprocs[num_largeside] = pid;
                grid_data->sideprcpos[num_largeside] = 1;
                grid_data->sideprcxsrt[num_largeside] = jxf_max(x_lower_idx, grid_data->xlo_array[pid]);
                grid_data->sideprcysrt[num_largeside] = y_upper_idx - 1;
                grid_data->sideprclength[num_largeside] = jxf_min(x_upper_idx,
                                   grid_data->xup_array[pid]) - jxf_max(x_lower_idx, grid_data->xlo_array[pid]);
                num_largeside ++;
            }
        }
        else if (grid_data->yup_array[pid] == y_lower_idx) // case 6
        {
            if (grid_data->xlo_array[pid] == x_upper_idx)
            {
                grid_data->sideprocs[num_small_dlr] = pid;
                grid_data->sideprcpos[num_small_dlr] = 8;
                grid_data->sideprcxsrt[num_small_dlr] = x_upper_idx - 1;
                grid_data->sideprcysrt[num_small_dlr] = y_lower_idx;
                grid_data->sideprclength[num_small_dlr] = 1;
                num_small_dlr ++;
            }
            else if (grid_data->xup_array[pid] == x_lower_idx)
            {
                grid_data->sideprocs[num_small_dlr] = pid;
                grid_data->sideprcpos[num_small_dlr] = 7;
                grid_data->sideprcxsrt[num_small_dlr] = x_lower_idx;
                grid_data->sideprcysrt[num_small_dlr] = y_lower_idx;
                grid_data->sideprclength[num_small_dlr] = 1;
                num_small_dlr ++;
            }
            else
            {
                grid_data->sideprocs[num_smallside] = pid;
                grid_data->sideprcpos[num_smallside] = 2;
                grid_data->sideprcxsrt[num_smallside] = jxf_max(x_lower_idx, grid_data->xlo_array[pid]);
                grid_data->sideprcysrt[num_smallside] = y_lower_idx;
                grid_data->sideprclength[num_smallside] = jxf_min(x_upper_idx,
                                   grid_data->xup_array[pid]) - jxf_max(x_lower_idx, grid_data->xlo_array[pid]);
                num_smallside ++;
            }
        }
    }
    //jxf_assert(num_large_ulr == num_sideproc);
    jxf_TFree(tmpsidepro);
    
    return 0;
}

/*!
 * \fn JXF_Int jxf_GridPartitionDataFinalize
 * \brief Finalize
 * \author Yue Xiaoqiang
 * \date 2014/08/03
 */
JXF_Int
jxf_GridPartitionDataFinalize( void *grid_vdata )
{
    jxf_GridPartitionData *grid_data = grid_vdata;
    
    if (grid_data)
    {
        if (grid_data->xlo_array)
        {
           jxf_TFree(grid_data->xlo_array);
        }
        if (grid_data->xup_array)
        {
           jxf_TFree(grid_data->xup_array);
        }
        if (grid_data->ylo_array)
        {
           jxf_TFree(grid_data->ylo_array);
        }
        if (grid_data->yup_array)
        {
           jxf_TFree(grid_data->yup_array);
        }
        if (grid_data->sideprocs)
        {
            jxf_TFree(grid_data->sideprocs);
        }
        if (grid_data->sideprcpos)
        {
           jxf_TFree(grid_data->sideprcpos);
        }
        if (grid_data->sideprcxsrt)
        {
           jxf_TFree(grid_data->sideprcxsrt);
        }
        if (grid_data->sideprcysrt)
        {
           jxf_TFree(grid_data->sideprcysrt);
        }
        if (grid_data->sideprclength)
        {
           jxf_TFree(grid_data->sideprclength);
        }
        jxf_TFree(grid_data);
    }
    
    return 0;
}
