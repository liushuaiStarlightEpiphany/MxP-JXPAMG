//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ilucycle.c
 *
 * Created by Yue Xiaoqiang 2014/03/24
 * Xiangtan University
 *
 */

#include "jxf_ilu.h"

/*!
 * \fn jxf_ILUZeroFactorDataPrecond
 * \brief Precond
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
jxf_ILUZeroFactorDataPrecond( void            *ilu_vdata,
                             jxf_ParCSRMatrix *par_A,
                             jxf_ParVector    *par_b,
                             jxf_ParVector    *par_x )
{
    MPI_Comm comm = jxf_ParCSRMatrixComm(par_A);
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    jxf_CSRMatrix *A = NULL;
    jxf_Vector *f = NULL;
    jxf_Vector *u = NULL;
    JXF_Int num_procs, my_id;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    if (num_procs == 1)
    {
        A = jxf_ParCSRMatrixDiag(par_A);
        f = jxf_ParVectorLocalVector(par_b);
        u = jxf_ParVectorLocalVector(par_x);
        //jxf_printf(" >>> Precond: Drop-Tol = %le on ID = %d of %d\n", ilu_data->drop_tol, my_id, num_procs);
        if (ilu_data->drop_tol == 0.0)
        {
            jxf_ILUZeroFactorDataCycleA(ilu_data, A, f, u);
        }
        else
        {
            jxf_ILUZeroFactorDataCycleB(ilu_data, A, f, u);
        }
    }
    else if (ilu_data->matA != NULL) /* Still Serial ILU(0), Yue Xiaoqiang, 2014/04/18 */
    {
        f = jxf_ParVectorToVectorAll(par_b);
        if (my_id == 0)
        {
            A = ilu_data->matA;
            //jxf_printf(" >>> Precond: Drop-Tol = %le on ID = %d of %d\n", ilu_data->drop_tol, my_id, num_procs);
            if (ilu_data->drop_tol == 0.0)
            {
                jxf_ILUZeroFactorDataCycleA(ilu_data, A, f, ilu_data->tmp_vec);
            }
            else
            {
                jxf_ILUZeroFactorDataCycleB(ilu_data, A, f, ilu_data->tmp_vec);
            }
        }
        jxf_VectorToParVector_Allocated2(comm, ilu_data->tmp_vec, jxf_ParVectorPartitioning(par_x), par_x);
        jxf_SeqVectorDestroy(f);
    }
    else if (ilu_data->drop_tol == 0.0)
    {
        jxf_ILUZeroFactorDataParallelCycleA(ilu_data, par_A, par_b, par_x);
    }
    else
    {
        jxf_printf("\n >>> Wrong parameters in ILU Cycle\n");
        exit(0);
    }
    
    return 0;
}

/*!
 * \fn void jxf_ILUZeroFactorDataCycleA
 * \brief Cycle
 * \author Yue Xiaoqiang
 * \date 2014/03/16
 */
void
jxf_ILUZeroFactorDataCycleA( void *ilu_vdata, jxf_CSRMatrix *A, jxf_Vector *f, jxf_Vector *u )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    JXF_Int num_rows = jxf_VectorSize(u);
    JXF_Int *IA = jxf_CSRMatrixI(A);
    JXF_Real *f_data = jxf_VectorData(f);
    JXF_Real *u_data = jxf_VectorData(u);
    JXF_Int *indexDP = ilu_data->indexD;
    JXF_Int *indexLU = ilu_data->index;
    JXF_Real *work = jxf_VectorData(ilu_data->aux_vec);
    JXF_Real *valueLU = ilu_data->value;
    JXF_Int i, j, row_end = 0;
    JXF_Real tmp;
    
    work[0] = f_data[0];
    for (i = 1; i < num_rows; i ++)
    {
        tmp = f_data[i];
        row_end = indexDP[i];
        for (j = IA[i]; j < row_end; j ++)
        {
            tmp -= valueLU[j] * work[indexLU[j]];
        }
        work[i] = tmp;
    }
    u_data[num_rows-1] = work[num_rows-1] * valueLU[row_end];
    for (i = num_rows-2; i > -1; i --)
    {
        tmp = work[i];
        row_end = IA[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            tmp -= valueLU[j] * u_data[indexLU[j]];
        }
        u_data[i] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataCycleB
 * \brief Cycle
 * \author Yue Xiaoqiang
 * \date 2014/03/16
 */
void
jxf_ILUZeroFactorDataCycleB( void *ilu_vdata, jxf_CSRMatrix *A, jxf_Vector *f, jxf_Vector *u )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    JXF_Int num_rows = jxf_VectorSize(u);
    JXF_Int *indexAP = ilu_data->indexA;
    JXF_Real *f_data = jxf_VectorData(f);
    JXF_Real *u_data = jxf_VectorData(u);
    JXF_Int *indexDP = ilu_data->indexD;
    JXF_Int *indexLU = ilu_data->index;
    JXF_Real *work = jxf_VectorData(ilu_data->aux_vec);
    JXF_Real *valueLU = ilu_data->value;
    JXF_Int i, j, row_end = 0;
    JXF_Real tmp;
    
    work[0] = f_data[0];
    for (i = 1; i < num_rows; i ++)
    {
        tmp = f_data[i];
        row_end = indexDP[i];
        for (j = indexAP[i]; j < row_end; j ++)
        {
            tmp -= valueLU[j] * work[indexLU[j]];
        }
        work[i] = tmp;
    }
    u_data[num_rows-1] = work[num_rows-1] * valueLU[row_end];
    for (i = num_rows-2; i > -1; i --)
    {
        tmp = work[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            tmp -= valueLU[j] * u_data[indexLU[j]];
        }
        u_data[i] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelCycleA
 * \brief Parallel Cycle
 * \author Yue Xiaoqiang
 * \date 2014/10/28
 */
void
jxf_ILUZeroFactorDataParallelCycleA( jxf_ILUZeroFactorData *ilu_data,
                                    jxf_ParCSRMatrix *par_A,
                                    jxf_ParVector *par_b,
                                    jxf_ParVector *par_x )
{
    MPI_Comm comm = jxf_ParCSRMatrixComm(par_A);
    JXF_Int num_rows = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(par_A));
    JXF_Int *row_starts = jxf_ParCSRMatrixRowStarts(par_A);
    jxf_Vector *rhs_vec = jxf_ParVectorLocalVector(par_b);
    jxf_Vector *app_vec = jxf_ParVectorLocalVector(par_x);
    JXF_Int ex_len = ilu_data->ex_len;
    JXF_Int dx_len = ilu_data->dx_len;
    JXF_Int kx_len = ilu_data->kx_len;
    JXF_Int dy_len = ilu_data->dy_len;
    JXF_Int lx_len = ilu_data->lx_len;
    JXF_Int ly_len = ilu_data->ly_len;
    JXF_Int postn_b = ilu_data->pos_dwn_end;
    JXF_Int postn_c = ilu_data->pos_lft_end;
    JXF_Int ng_pt = ilu_data->num_equns;
    JXF_Int *indexLU = ilu_data->index;
    JXF_Int *indexAP = ilu_data->indexA + ex_len;
    JXF_Int *indexDP = ilu_data->indexD + ex_len;
    JXF_Int *permute = ilu_data->permute;
    JXF_Real *valueLU = ilu_data->value;
    JXF_Real *senddown = ilu_data->senddown;
    MPI_Status *status = ilu_data->status;
    jxf_Vector *aux_vec = jxf_ParVectorLocalVector(ilu_data->par_aux_vec);
    jxf_Vector *res_vec = jxf_ParVectorLocalVector(ilu_data->par_res_vec);
    jxf_GridPartitionData *grid_data = ilu_data->par_grid;
    JXF_Int num_smallside = grid_data->num_smallside;
    JXF_Int num_largeside = grid_data->num_largeside; // YUE: Needn't to consider the cross ones in 5-point stencil
    JXF_Int *smallprocs = grid_data->sideprocs;
    JXF_Int *largeprocs = smallprocs + num_smallside;
    JXF_Int *smallprcpos = grid_data->sideprcpos;
    JXF_Int *largeprcpos = smallprcpos + num_smallside;
    JXF_Int *smallprclength = grid_data->sideprclength;
    JXF_Int *largeprclength = smallprclength + num_smallside;
    JXF_Int ng_pt_mx = ng_pt * grid_data->x_part_len;
    //JXF_Int ng_pt_my = ng_pt * grid_data->y_part_len;
    JXF_Int num_procs, my_id, j, k, pid, jend, jjend, cnt, dnt;
    JXF_Int int_uprgt_pnt, recv_downsrt, recv_leftsrt, recv_dnftsrt;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    int_uprgt_pnt = num_rows - ilu_data->ey_len;
    /* TYPE 1. Back substitution of L parts for interior, up and right quasi-boundary points */
    jxf_ILUZeroFactorDataParallelLPartIntURPntsA(indexLU, indexAP, indexDP, valueLU,
                                                permute, int_uprgt_pnt, jxf_VectorData(aux_vec),
                                                jxf_VectorData(rhs_vec), row_starts[my_id]);
    /* Send back substitution of L parts for UP and RIGHT quasi-boundary points. Xiaoqiang Yue, 2014/04/28 */
    recv_downsrt = 0;
    /* jxf_MPI_Irecv and jxf_MPI_Isend shouldn't be used here as summations of each offset */
    /* Receive data from down processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 2) /* DOWN */
        {
            jend = (smallprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_downsrt,
                jend, JXF_MPI_REAL, smallprocs[pid], my_id*233, comm, &status[pid]);
            recv_downsrt += jend; // YUE: aims to loop for pid
        }
    }
    recv_leftsrt = recv_downsrt;
    /* Receive data from left processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 4) /* LEFT */
        {
            jend = (smallprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_leftsrt,
                jend, JXF_MPI_REAL, smallprocs[pid], my_id*235, comm, &status[pid]);
            recv_leftsrt += jend; // YUE: aims to loop for pid
        }
    }
    /* Fill data, Send data to up processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 1) /* UP */
        {
            jend = (largeprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Send(jxf_VectorData(aux_vec)+dy_len, jend, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*233, comm);
        }
    }
    /* Fill data, Send data to right processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 3) /* RIGHT */
        {
            jjend = largeprclength[pid] - 1; // YUE: exclude the first block
            jend = jjend * ng_pt;
            cnt = 0;
            dnt = dx_len;
            for (j = 0; j < jjend; j ++)
            {
                for (k = 0; k < ng_pt; k ++)
                {
                    senddown[cnt++] = jxf_VectorData(aux_vec)[dnt+k];
                }
                dnt += kx_len;
            }
            jxf_MPI_Send(senddown, jend, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*235, comm);
        }
    }
    //jxf_MPI_Barrier(comm);
    //exit(0);
    //jxf_assert(postn_a == int_uprgt_pnt);
    /* TYPE 2. Back substitution of L parts for down quasi-boundary points */
    jxf_ILUZeroFactorDataParallelLPartDPntsA(indexLU, indexAP, indexDP, valueLU,
                                            int_uprgt_pnt, postn_b, permute,
                                            jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                            jxf_VectorData(rhs_vec), row_starts[my_id]);
    /* TYPE 3. Back substitution of L parts for left quasi-boundary points */
    jxf_ILUZeroFactorDataParallelLPartLPntsA(indexLU, indexAP, indexDP, valueLU,
                                            postn_b, postn_c, recv_downsrt, permute,
                                            jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                            jxf_VectorData(rhs_vec), row_starts[my_id]);
    recv_dnftsrt = recv_downsrt = recv_leftsrt;
    /* Receive data from down processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 2) /* DOWN */
        {
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_downsrt,
                ng_pt, JXF_MPI_REAL, smallprocs[pid], my_id*237, comm, &status[pid]);
            recv_downsrt += ng_pt; // YUE: aims to loop for pid
        }
    }
    recv_leftsrt = recv_downsrt;
    /* Receive data from left processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 4) /* LEFT */
        {
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_leftsrt,
                ng_pt, JXF_MPI_REAL, smallprocs[pid], my_id*231, comm, &status[pid]);
            recv_leftsrt += ng_pt; // YUE: aims to loop for pid
        }
    }
    /* Fill data, Send data to up processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 1) /* UP */
        {
            jxf_MPI_Send(jxf_VectorData(aux_vec)+ly_len, ng_pt, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*237, comm);
        }
    }
    /* Fill data, Send data to right processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 3) /* RIGHT */
        {
            jxf_MPI_Send(jxf_VectorData(aux_vec)+lx_len, ng_pt, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*231, comm);
        }
    }
    //jxf_MPI_Barrier(comm);
    //exit(0);
    /* TYPE 4. Back substitution of L parts for LD-Corner quasi-boundary points */
    if (my_id != 0)
    {
        jxf_ILUZeroFactorDataParallelLPartLDPntsA(indexLU, indexAP, indexDP, valueLU,
                                                 postn_c, num_rows, row_starts[my_id-1],
                                                 permute, recv_dnftsrt, recv_downsrt,
                                                 jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                                 jxf_VectorData(rhs_vec), row_starts[my_id]);
    }
    /* TYPE 4. Back substitution of U parts for LD-Corner quasi-boundary points */
    if (my_id != 0)
    {
        jxf_ILUZeroFactorDataParallelUPartLDPntsA(indexLU, indexAP, indexDP,
                                                 valueLU, postn_c, num_rows,
                                                 permute, jxf_VectorData(aux_vec),
                                                 jxf_VectorData(app_vec), row_starts[my_id]);
    }
    recv_downsrt = 0;
    /* Receive data from UP processor, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 1) /* UP */
        {
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_downsrt,
                ng_pt, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*234, comm, &status[pid]);
            recv_downsrt += ng_pt; // YUE: aims to loop for pid
        }
    }
    recv_leftsrt = ng_pt_mx; // YUE: pay attention to the offset
    /* Receive data from RIGHT processor, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 3) /* RIGHT */
        {
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_leftsrt,
                ng_pt, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*236, comm, &status[pid]);
            recv_leftsrt += ng_pt; // YUE: aims to loop for pid
        }
    }
    /* Fill data, Send data to down processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 2) /* DOWN */
        {
            jxf_MPI_Send(jxf_VectorData(app_vec), ng_pt, JXF_MPI_REAL, smallprocs[pid], my_id*234, comm);
        }
    }
    /* Fill data, Send data to left processors, 1 block, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 4) /* LEFT */
        {
            jxf_MPI_Send(jxf_VectorData(app_vec), ng_pt, JXF_MPI_REAL, smallprocs[pid], my_id*236, comm);
        }
    }
    //jxf_MPI_Barrier(comm);
    //exit(0);
    /* TYPE 3. Back substitution of U parts for left quasi-boundary points */
    jxf_ILUZeroFactorDataParallelUPartLPntsA(indexLU, indexAP, indexDP, valueLU,
                                            postn_b, postn_c, recv_downsrt, permute,
                                            jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                            jxf_VectorData(app_vec), row_starts[my_id], row_starts[my_id+1]);
    /* TYPE 2. Back substitution of U parts for down quasi-boundary points */
    jxf_ILUZeroFactorDataParallelUPartDPntsA(indexLU, indexAP, indexDP, valueLU,
                                            int_uprgt_pnt, postn_b, recv_leftsrt,
                                            permute, jxf_VectorData(aux_vec),
                                            jxf_VectorData(res_vec), jxf_VectorData(app_vec),
                                            row_starts[my_id], row_starts[my_id+1]);
    //jxf_MPI_Barrier(comm);
    //exit(0);
    /* Receive data from UP processor, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 1) /* UP */
        {
            jend = (largeprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_downsrt,
                jend, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*232, comm, &status[pid]);
            recv_downsrt += jend; // YUE: aims to loop for pid
        }
    }
    /* Receive data from RIGHT processor, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_largeside; pid ++)
    {
        if (largeprcpos[pid] == 3) /* RIGHT */
        {
            jend = (largeprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Recv(jxf_VectorData(res_vec)+recv_leftsrt,
                jend, JXF_MPI_REAL, largeprocs[pid], largeprocs[pid]*238, comm, &status[pid]);
            recv_leftsrt += jend; // YUE: aims to loop for pid
        }
    }
    /* Fill data, Send data to down processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 2) /* DOWN */
        {
            jend = (smallprclength[pid] - 1) * ng_pt; // YUE: exclude the first block
            jxf_MPI_Send(jxf_VectorData(app_vec)+ng_pt, jend, JXF_MPI_REAL, smallprocs[pid], my_id*232, comm);
        }
    }
    /* Fill data, Send data to left processors, Yue Xiaoqiang, 2014/08/06 */
    for (pid = 0; pid < num_smallside; pid ++)
    {
        if (smallprcpos[pid] == 4) /* LEFT */
        {
            jjend = smallprclength[pid] - 1; // YUE: exclude the first block
            jend = jjend * ng_pt;
            dnt = cnt = 0;
            for (j = 0; j < jjend; j ++)
            {
                dnt += ng_pt_mx;
                for (k = 0; k < ng_pt; k ++)
                {
                    senddown[cnt++] = jxf_VectorData(app_vec)[dnt+k];
                }
            }
            jxf_MPI_Send(senddown, jend, JXF_MPI_REAL, smallprocs[pid], my_id*238, comm);
        }
    }
    //jxf_MPI_Barrier(comm);
    //exit(0);
    /* TYPE 1: Back substitution of U parts for interior, up and right quasi-boundary points */
    if (grid_data->part_type == 1) /* UP and RIGHT both exist */
    {
        //jxf_printf("\n >>> 1 Proc %d of %d\n\n", my_id, num_procs);
        jxf_ILUZeroFactorDataParallelUPartIntURPntsA(indexLU, indexAP, indexDP, valueLU,
                                                    int_uprgt_pnt, recv_downsrt, recv_leftsrt,
                                                    permute, jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                                    jxf_VectorData(app_vec), row_starts[my_id],
                                                    row_starts[my_id+1], row_starts[my_id+2]);
    }
    else if (grid_data->part_type == 2) /* Only UP exists */
    {
        //jxf_printf("\n >>> 2 Proc %d of %d\n\n", my_id, num_procs);
        jxf_ILUZeroFactorDataParallelUPartIntURPntsB(indexLU, indexAP, indexDP, valueLU,
                                                    int_uprgt_pnt, recv_downsrt, permute,
                                                    jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                                    jxf_VectorData(app_vec), row_starts[my_id], row_starts[my_id+1]);
    }
    else if (grid_data->part_type == 3) /* Only RIGHT exists */
    {
        //jxf_printf("\n >>> 3 Proc %d of %d\n\n", my_id, num_procs);
        jxf_ILUZeroFactorDataParallelUPartIntURPntsC(indexLU, indexAP, indexDP, valueLU,
                                                    int_uprgt_pnt, recv_leftsrt, permute,
                                                    jxf_VectorData(aux_vec), jxf_VectorData(res_vec),
                                                    jxf_VectorData(app_vec), row_starts[my_id], row_starts[my_id+1]);
    }
    else if (grid_data->part_type == 4) /* None of UP or RIGHT */
    {
        //jxf_printf("\n >>> 4 Proc %d of %d\n\n", my_id, num_procs);
        jxf_ILUZeroFactorDataParallelUPartIntURPntsD(indexLU, indexAP, indexDP, valueLU,
                                                    int_uprgt_pnt, permute, jxf_VectorData(aux_vec),
                                                    jxf_VectorData(res_vec), jxf_VectorData(app_vec),
                                                    row_starts[my_id]);
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelLPartIntURPntsA
 * \brief Parallel Cycle of L parts for Int-U-R points
 * \author Yue Xiaoqiang
 * \date 2014/08/26
 */
void
jxf_ILUZeroFactorDataParallelLPartIntURPntsA( JXF_Int *indexLU,
                                             JXF_Int *indexAP,
                                             JXF_Int *indexDP,
                                             JXF_Real *valueLU,
                                             JXF_Int *permute,
                                             JXF_Int int_uprgt_pnt,
                                             JXF_Real *aux_data,
                                             JXF_Real *rhs_data,
                                             JXF_Int first_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    aux_data[0] = rhs_data[permute[0]];
    for (i = 1; i < int_uprgt_pnt; i ++)
    {
        tmp = rhs_data[permute[i]];
        row_end = indexDP[i];
        for (j = indexAP[i]; j < row_end; j ++)
        {
            tmp -= valueLU[j] * aux_data[indexLU[j]-first_row_idx];
        }
        aux_data[i] = tmp;
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelLPartDPntsA
 * \brief Parallel Cycle of L parts for D points
 * \author Yue Xiaoqiang
 * \date 2014/08/28
 */
void
jxf_ILUZeroFactorDataParallelLPartDPntsA( JXF_Int *indexLU,
                                         JXF_Int *indexAP,
                                         JXF_Int *indexDP,
                                         JXF_Real *valueLU,
                                         JXF_Int int_uprgt_pnt,
                                         JXF_Int num_rows,
                                         JXF_Int *permute,
                                         JXF_Real *aux_data,
                                         JXF_Real *res_data,
                                         JXF_Real *rhs_data,
                                         JXF_Int first_row_idx )
{
    JXF_Int i, j, l, row_end, sstmp = 0;
    JXF_Real tmp;
    
    for (i = int_uprgt_pnt; i < num_rows; i ++)
    {
        tmp = rhs_data[permute[i]];
        row_end = indexDP[i];
        for (j = indexAP[i]; j < row_end; j ++)
        {
            l = indexLU[j] - first_row_idx;
            if (l < 0)
            {
                tmp -= valueLU[j] * res_data[sstmp];
                sstmp ++;
            }
            else
            {
                tmp -= valueLU[j] * aux_data[l];
            }
        }
        aux_data[i] = tmp;
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelLPartLPntsA
 * \brief Parallel Cycle of L parts for L points
 * \author Yue Xiaoqiang
 * \date 2014/08/28
 */
void
jxf_ILUZeroFactorDataParallelLPartLPntsA( JXF_Int *indexLU,
                                         JXF_Int *indexAP,
                                         JXF_Int *indexDP,
                                         JXF_Real *valueLU,
                                         JXF_Int int_uprgtdwn_pnt,
                                         JXF_Int num_rows,
                                         JXF_Int dwn_num_rows,
                                         JXF_Int *permute,
                                         JXF_Real *aux_data,
                                         JXF_Real *res_data,
                                         JXF_Real *rhs_data,
                                         JXF_Int first_row_idx )
{
    JXF_Int i, j, l, row_end, sstmp = dwn_num_rows;
    JXF_Real tmp;
    
    for (i = int_uprgtdwn_pnt; i < num_rows; i ++)
    {
        tmp = rhs_data[permute[i]];
        row_end = indexDP[i];
        for (j = indexAP[i]; j < row_end; j ++)
        {
            l = indexLU[j] - first_row_idx;
            if (l < 0)
            {
                tmp -= valueLU[j] * res_data[sstmp];
                sstmp ++;
            }
            else
            {
                tmp -= valueLU[j] * aux_data[l];
            }
        }
        aux_data[i] = tmp;
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelLPartLDPntsA
 * \brief Parallel Cycle of L parts for LD-Corner points
 * \author Yue Xiaoqiang
 * \date 2014/08/28
 */
void
jxf_ILUZeroFactorDataParallelLPartLDPntsA( JXF_Int *indexLU,
                                          JXF_Int *indexAP,
                                          JXF_Int *indexDP,
                                          JXF_Real *valueLU,
                                          JXF_Int int_urdl_pnt,
                                          JXF_Int num_rows,
                                          JXF_Int fst_row_idx,
                                          JXF_Int *permute,
                                          JXF_Int lft_cnum_rows,
                                          JXF_Int dwn_dnum_rows,
                                          JXF_Real *aux_data,
                                          JXF_Real *res_data,
                                          JXF_Real *rhs_data,
                                          JXF_Int first_row_idx )
{
    JXF_Int sstmp = lft_cnum_rows, kktmp = dwn_dnum_rows;
    JXF_Int i, j, l, row_end;
    JXF_Real tmp;
    
    for (i = int_urdl_pnt; i < num_rows; i ++)
    {
        tmp = rhs_data[permute[i]];
        row_end = indexDP[i];
        for (j = indexAP[i]; j < row_end; j ++)
        {
            l = indexLU[j] - first_row_idx;
            if (l < 0)
            {
                if (l+first_row_idx < fst_row_idx) // DOWN
                {
                    tmp -= valueLU[j] * res_data[sstmp];
                    sstmp ++;
                }
                else // LEFT
                {
                    tmp -= valueLU[j] * res_data[kktmp];
                    kktmp ++;
                }
            }
            else
            {
                tmp -= valueLU[j] * aux_data[l];
            }
        }
        aux_data[i] = tmp;
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartLDPntsA
 * \brief Parallel Cycle of U parts for LD-Corner points
 * \author Yue Xiaoqiang
 * \date 2014/08/30
 */
void
jxf_ILUZeroFactorDataParallelUPartLDPntsA( JXF_Int *indexLU,
                                          JXF_Int *indexAP,
                                          JXF_Int *indexDP,
                                          JXF_Real *valueLU,
                                          JXF_Int int_urdl_pnt,
                                          JXF_Int num_rows,
                                          JXF_Int *permute,
                                          JXF_Real *aux_data,
                                          JXF_Real *app_data,
                                          JXF_Int first_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    app_data[permute[num_rows-1]] = aux_data[num_rows-1] * valueLU[indexDP[num_rows-1]];
    for (i = num_rows-2; i >= int_urdl_pnt; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartLPntsA
 * \brief Parallel Cycle of U parts for L points
 * \author Yue Xiaoqiang
 * \date 2014/08/30
 */
void
jxf_ILUZeroFactorDataParallelUPartLPntsA( JXF_Int *indexLU,
                                         JXF_Int *indexAP,
                                         JXF_Int *indexDP,
                                         JXF_Real *valueLU,
                                         JXF_Int int_urdl_pnt,
                                         JXF_Int num_rows,
                                         JXF_Int recv_downsrt,
                                         JXF_Int *permute,
                                         JXF_Real *aux_data,
                                         JXF_Real *res_data,
                                         JXF_Real *app_data,
                                         JXF_Int first_row_idx,
                                         JXF_Int next_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = num_rows-1; i >= int_urdl_pnt; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            if (indexLU[j] < next_row_idx)
            {
                tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
            }
            else
            {
                recv_downsrt --;
                tmp -= valueLU[j] * res_data[recv_downsrt];
            }
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartDPntsA
 * \brief Parallel Cycle of U parts for D points
 * \author Yue Xiaoqiang
 * \date 2014/08/30
 */
void
jxf_ILUZeroFactorDataParallelUPartDPntsA( JXF_Int *indexLU,
                                         JXF_Int *indexAP,
                                         JXF_Int *indexDP,
                                         JXF_Real *valueLU,
                                         JXF_Int int_uprgt_pnt,
                                         JXF_Int num_rows,
                                         JXF_Int recv_leftsrt,
                                         JXF_Int *permute,
                                         JXF_Real *aux_data,
                                         JXF_Real *res_data,
                                         JXF_Real *app_data,
                                         JXF_Int first_row_idx,
                                         JXF_Int next_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = num_rows-1; i >= int_uprgt_pnt; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            if (indexLU[j] < next_row_idx)
            {
                tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
            }
            else
            {
                recv_leftsrt --;
                tmp -= valueLU[j] * res_data[recv_leftsrt];
            }
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartIntURPntsA
 * \brief Parallel Cycle of U parts for Int-U-R points
 * \author Yue Xiaoqiang
 * \date 2014/09/03
 */
void
jxf_ILUZeroFactorDataParallelUPartIntURPntsA( JXF_Int *indexLU,
                                             JXF_Int *indexAP,
                                             JXF_Int *indexDP,
                                             JXF_Real *valueLU,
                                             JXF_Int int_uprgt_pnt,
                                             JXF_Int recv_downsrt,
                                             JXF_Int recv_leftsrt,
                                             JXF_Int *permute,
                                             JXF_Real *aux_data,
                                             JXF_Real *res_data,
                                             JXF_Real *app_data,
                                             JXF_Int first_row_idx,
                                             JXF_Int next_row_idx,
                                             JXF_Int nnxt_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = int_uprgt_pnt-1; i > -1; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            if (indexLU[j] < next_row_idx)
            {
                tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
            }
            else if (indexLU[j] < nnxt_row_idx)
            {
                recv_leftsrt --;
                tmp -= valueLU[j] * res_data[recv_leftsrt];
            }
            else
            {
                recv_downsrt --;
                tmp -= valueLU[j] * res_data[recv_downsrt];
            }
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartIntURPntsB
 * \brief Parallel Cycle of U parts for Int-U-R points
 * \author Yue Xiaoqiang
 * \date 2014/09/11
 */
void
jxf_ILUZeroFactorDataParallelUPartIntURPntsB( JXF_Int *indexLU,
                                             JXF_Int *indexAP,
                                             JXF_Int *indexDP,
                                             JXF_Real *valueLU,
                                             JXF_Int int_uprgt_pnt,
                                             JXF_Int recv_downsrt,
                                             JXF_Int *permute,
                                             JXF_Real *aux_data,
                                             JXF_Real *res_data,
                                             JXF_Real *app_data,
                                             JXF_Int first_row_idx,
                                             JXF_Int next_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = int_uprgt_pnt-1; i > -1; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            if (indexLU[j] < next_row_idx)
            {
                tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
            }
            else
            {
                recv_downsrt --;
                tmp -= valueLU[j] * res_data[recv_downsrt];
            }
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartIntURPntsC
 * \brief Parallel Cycle of U parts for Int-U-R points
 * \author Yue Xiaoqiang
 * \date 2014/09/17
 */
void
jxf_ILUZeroFactorDataParallelUPartIntURPntsC( JXF_Int *indexLU,
                                             JXF_Int *indexAP,
                                             JXF_Int *indexDP,
                                             JXF_Real *valueLU,
                                             JXF_Int int_uprgt_pnt,
                                             JXF_Int recv_leftsrt,
                                             JXF_Int *permute,
                                             JXF_Real *aux_data,
                                             JXF_Real *res_data,
                                             JXF_Real *app_data,
                                             JXF_Int first_row_idx,
                                             JXF_Int next_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = int_uprgt_pnt-1; i > -1; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            if (indexLU[j] < next_row_idx)
            {
                tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
            }
            else
            {
                recv_leftsrt --;
                tmp -= valueLU[j] * res_data[recv_leftsrt];
            }
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*!
 * \fn void jxf_ILUZeroFactorDataParallelUPartIntURPntsD
 * \brief Parallel Cycle of U parts for Int-U-R points
 * \author Yue Xiaoqiang
 * \date 2014/09/23
 */
void
jxf_ILUZeroFactorDataParallelUPartIntURPntsD( JXF_Int *indexLU,
                                             JXF_Int *indexAP,
                                             JXF_Int *indexDP,
                                             JXF_Real *valueLU,
                                             JXF_Int int_uprgt_pnt,
                                             JXF_Int *permute,
                                             JXF_Real *aux_data,
                                             JXF_Real *res_data,
                                             JXF_Real *app_data,
                                             JXF_Int first_row_idx )
{
    JXF_Int i, j, row_end;
    JXF_Real tmp;
    
    for (i = int_uprgt_pnt-1; i > -1; i --)
    {
        tmp = aux_data[i];
        row_end = indexAP[i+1];
        for (j = indexDP[i]+1; j < row_end; j ++)
        {
            tmp -= valueLU[j] * app_data[permute[indexLU[j]-first_row_idx]];
        }
        app_data[permute[i]] = tmp * valueLU[indexDP[i]];
    }
}

/*------------------------------------------------------------*/
/*------------------------ End of File -----------------------*/
/*------------------------------------------------------------*/
