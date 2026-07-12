//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  mat_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

static void jxf_setup_matvec_sends_private( jxf_Mat_dh mat, JXF_Int *inlist );
static void jxf_setup_matvec_receives_private( jxf_Mat_dh mat, JXF_Int *beg_rows,
                       JXF_Int *end_rows, JXF_Int reqlen, JXF_Int *reqind, JXF_Int *outlist );

static jxf_bool jxf_commsOnly = jxf_false;

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhCreate"
void jxf_Mat_dhCreate( jxf_Mat_dh *mat )
{
    JXF_START_FUNC_DH
    struct _jxf_mat_dh *tmp = (struct _jxf_mat_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_mat_dh)); JXF_CHECK_V_ERROR;
    
   *mat = tmp;
    jxf_commsOnly = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-jxf_commsOnly");
    if (jxf_myid_dh == 0 && jxf_commsOnly == jxf_true)
    {
        fflush(stdout);
    }
    tmp->m = 0;
    tmp->n = 0;
    tmp->beg_row = 0;
    tmp->bs = 1;
    tmp->rp = NULL;
    tmp->len = NULL;
    tmp->cval = NULL;
    tmp->aval = NULL;
    tmp->diag = NULL;
    tmp->fill = NULL;
    tmp->owner = jxf_true;
    tmp->len_private = 0;
    tmp->rowCheckedOut = -1;
    tmp->cval_private = NULL;
    tmp->aval_private = NULL;
    tmp->row_perm = NULL;
    tmp->num_recv = 0;
    tmp->num_send = 0;
    tmp->recv_req = NULL;
    tmp->send_req = NULL;
    tmp->status = NULL;
    tmp->recvbuf = NULL;
    tmp->sendbuf = NULL;
    tmp->sendind = NULL;
    tmp->sendlen = 0;
    tmp->recvlen = 0;
    tmp->numb = NULL;
    tmp->matvecIsSetup = jxf_false;
    jxf_Mat_dhZeroTiming(tmp); JXF_CHECK_V_ERROR;
    tmp->matvec_timing = jxf_true;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_Mat");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhDestroy"
void jxf_Mat_dhDestroy( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
    if (mat->owner)
    {
        if (mat->rp != NULL)
        {
            JXF_FREE_DH(mat->rp); JXF_CHECK_V_ERROR;
        }
        if (mat->len != NULL)
        {
            JXF_FREE_DH(mat->len); JXF_CHECK_V_ERROR;
        }
        if (mat->cval != NULL)
        {
            JXF_FREE_DH(mat->cval); JXF_CHECK_V_ERROR;
        }
        if (mat->aval != NULL)
        {
            JXF_FREE_DH(mat->aval); JXF_CHECK_V_ERROR;
        }
        if (mat->diag != NULL)
        {
            JXF_FREE_DH(mat->diag); JXF_CHECK_V_ERROR;
        }
        if (mat->fill != NULL)
        {
            JXF_FREE_DH(mat->fill); JXF_CHECK_V_ERROR;
        }
        if (mat->cval_private != NULL)
        {
            JXF_FREE_DH(mat->cval_private); JXF_CHECK_V_ERROR;
        }
        if (mat->aval_private != NULL)
        {
            JXF_FREE_DH(mat->aval_private); JXF_CHECK_V_ERROR;
        }
        if (mat->row_perm != NULL)
        {
            JXF_FREE_DH(mat->row_perm); JXF_CHECK_V_ERROR;
        }
    }
    for (i = 0; i < mat->num_recv; i ++) jxf_MPI_Request_free(&mat->recv_req[i]);
    for (i = 0; i < mat->num_send; i ++) jxf_MPI_Request_free(&mat->send_req[i]);
    if (mat->recv_req != NULL)
    {
        JXF_FREE_DH(mat->recv_req); JXF_CHECK_V_ERROR;
    }
    if (mat->send_req != NULL)
    {
        JXF_FREE_DH(mat->send_req); JXF_CHECK_V_ERROR;
    }
    if (mat->status != NULL)
    {
        JXF_FREE_DH(mat->status); JXF_CHECK_V_ERROR;
    }
    if (mat->recvbuf != NULL)
    {
        JXF_FREE_DH(mat->recvbuf); JXF_CHECK_V_ERROR;
    }
    if (mat->sendbuf != NULL)
    {
        JXF_FREE_DH(mat->sendbuf); JXF_CHECK_V_ERROR;
    }
    if (mat->sendind != NULL)
    {
        JXF_FREE_DH(mat->sendind); JXF_CHECK_V_ERROR;
    }
    if (mat->matvecIsSetup)
    {
        jxf_Mat_dhMatVecSetdown(mat); JXF_CHECK_V_ERROR;
    }
    if (mat->numb != NULL)
    {
        jxf_Numbering_dhDestroy(mat->numb); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(mat); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVecSetDown"
void jxf_Mat_dhMatVecSetdown( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    if (jxf_ignoreMe) JXF_SET_V_ERROR("not implemented");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVecSetup"
void jxf_Mat_dhMatVecSetup( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh == 1)
    {
        goto DO_NOTHING;
    }
    else
    {
        JXF_Int *outlist, *inlist;
        JXF_Int ierr, i, row, *rp = mat->rp, *cval = mat->cval;
        jxf_Numbering_dh numb;
        JXF_Int m = mat->m;
        JXF_Int firstLocal = mat->beg_row;
        JXF_Int lastLocal = firstLocal+m;
        JXF_Int *beg_rows, *end_rows;
        
        mat->recv_req = (MPI_Request *)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
        mat->send_req = (MPI_Request *)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
        mat->status = (MPI_Status *)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
        beg_rows = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        end_rows = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        if (jxf_np_dh == 1)
        {
            beg_rows[0] = 0;
            end_rows[0] = m;
        }
        else
        {
            ierr = jxf_MPI_Allgather(&firstLocal, 1, JXF_MPI_INT, beg_rows, 1, JXF_MPI_INT, jxf_comm_dh); JXF_CHECK_MPI_V_ERROR(ierr);
            ierr = jxf_MPI_Allgather(&lastLocal, 1, JXF_MPI_INT, end_rows, 1, JXF_MPI_INT, jxf_comm_dh); JXF_CHECK_MPI_V_ERROR(ierr);
        }
        outlist = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        inlist = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        for (i = 0; i < jxf_np_dh; ++ i)
        {
            outlist[i] = 0;
            inlist[i] = 0;
        }
        jxf_Numbering_dhCreate(&(mat->numb)); JXF_CHECK_V_ERROR;
        numb = mat->numb;
        jxf_Numbering_dhSetup(numb, mat); JXF_CHECK_V_ERROR;
        jxf_setup_matvec_receives_private(mat, beg_rows, end_rows,
                          numb->num_ext, numb->idx_ext, outlist); JXF_CHECK_V_ERROR;
        if (jxf_np_dh == 1)
        {
            inlist[0] = outlist[0];
        }
        else
        {
            ierr = jxf_MPI_Alltoall(outlist, 1, JXF_MPI_INT, inlist, 1, JXF_MPI_INT, jxf_comm_dh); JXF_CHECK_MPI_V_ERROR(ierr);
        }
        jxf_setup_matvec_sends_private(mat, inlist); JXF_CHECK_V_ERROR;
        for (row = 0; row < m; row ++)
        {
            JXF_Int len = rp[row+1] - rp[row];
            JXF_Int *ind = cval+rp[row];
            
            jxf_Numbering_dhGlobalToLocal(numb, len, ind, ind); JXF_CHECK_V_ERROR;
        }
        JXF_FREE_DH(outlist); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(inlist); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(beg_rows); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(end_rows); JXF_CHECK_V_ERROR;
    }
    
DO_NOTHING: ;
    
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_setup_matvec_receives_private"
void jxf_setup_matvec_receives_private( jxf_Mat_dh mat,
                                    JXF_Int *beg_rows,
                                    JXF_Int *end_rows,
                                    JXF_Int reqlen,
                                    JXF_Int *reqind,
                                    JXF_Int *outlist )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, i, j, this_pe;
    MPI_Request request;
    JXF_Int m = mat->m;
    
    mat->num_recv = 0;
    mat->recvbuf = (JXF_Real *)JXF_MALLOC_DH((reqlen+m)*sizeof(JXF_Real));
    for (i = 0; i < reqlen; i = j)
    {
        this_pe = jxf_mat_jxf_find_owner(beg_rows, end_rows, reqind[i]); JXF_CHECK_V_ERROR;
        for (j = i+1; j < reqlen; j ++)
        {
            if (reqind[j] < beg_rows[this_pe] || reqind[j] > end_rows[this_pe]) break;
        }
        ierr = jxf_MPI_Isend(&reqind[i], j-i, JXF_MPI_INT, this_pe, 444, jxf_comm_dh, &request); JXF_CHECK_MPI_V_ERROR(ierr);
        ierr = jxf_MPI_Request_free(&request); JXF_CHECK_MPI_V_ERROR(ierr);
        outlist[this_pe] = j - i;
        ierr = jxf_MPI_Recv_init(&mat->recvbuf[i+m], j-i,
                 JXF_MPI_REAL, this_pe, 555, jxf_comm_dh, &mat->recv_req[mat->num_recv]); JXF_CHECK_MPI_V_ERROR(ierr);
        mat->num_recv ++;
        mat->recvlen += j - i;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_setup_matvec_sends_private"
void jxf_setup_matvec_sends_private( jxf_Mat_dh mat, JXF_Int *inlist )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, i, j, sendlen, first = mat->beg_row;
    MPI_Request *requests;
    MPI_Status  *statuses;
    
    requests = (MPI_Request *)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    statuses = (MPI_Status *)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
    sendlen = 0;
    for (i = 0; i < jxf_np_dh; i ++) sendlen += inlist[i];
    mat->sendlen = sendlen;
    mat->sendbuf = (JXF_Real *)JXF_MALLOC_DH(sendlen*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    mat->sendind = (JXF_Int *)JXF_MALLOC_DH(sendlen*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    j = 0;
    mat->num_send = 0;
    for (i = 0; i < jxf_np_dh; i ++)
    {
        if (inlist[i] != 0)
        {
            ierr = jxf_MPI_Irecv(&mat->sendind[j], inlist[i],
                      JXF_MPI_INT, i, 444, jxf_comm_dh, &requests[mat->num_send]); JXF_CHECK_MPI_V_ERROR(ierr);
            ierr = jxf_MPI_Send_init(&mat->sendbuf[j], inlist[i],
                    JXF_MPI_REAL, i, 555, jxf_comm_dh, &mat->send_req[mat->num_send]); JXF_CHECK_MPI_V_ERROR(ierr);
            mat->num_send ++;
            j += inlist[i];
        }
    }
    mat->time[JXF_MATVEC_WORDS] = j;
    ierr = jxf_MPI_Waitall(mat->num_send, requests, statuses); JXF_CHECK_MPI_V_ERROR(ierr);
    for (i = 0; i < mat->sendlen; i ++) mat->sendind[i] -= first;
    JXF_FREE_DH(requests);
    JXF_FREE_DH(statuses);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVec"
void jxf_Mat_dhMatVec( jxf_Mat_dh mat, JXF_Real *x, JXF_Real *b )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh == 1)
    {
        jxf_Mat_dhMatVec_uni(mat, x, b); JXF_CHECK_V_ERROR;
    }
    else
    {
        JXF_Int ierr, i, row, m = mat->m;
        JXF_Int *rp = mat->rp, *cval = mat->cval;
        JXF_Real *aval = mat->aval;
        JXF_Int *sendind = mat->sendind;
        JXF_Int sendlen = mat->sendlen;
        JXF_Real *sendbuf = mat->sendbuf;
        JXF_Real *recvbuf = mat->recvbuf;
        JXF_Real t1 = 0, t2 = 0, t3 = 0, t4 = 0;
        jxf_bool timeFlag = mat->matvec_timing;
        
        if (timeFlag) t1 = jxf_MPI_Wtime();
        if (!jxf_commsOnly)
        {
            for (i = 0; i < sendlen; i ++) sendbuf[i] = x[sendind[i]];
        }
        if (timeFlag)
        {
            t2 = jxf_MPI_Wtime();
            mat->time[JXF_MATVEC_TIME] += (t2 - t1);
        }
        ierr = jxf_MPI_Startall(mat->num_recv, mat->recv_req); JXF_CHECK_MPI_V_ERROR(ierr);
        ierr = jxf_MPI_Startall(mat->num_send, mat->send_req); JXF_CHECK_MPI_V_ERROR(ierr);
        ierr = jxf_MPI_Waitall(mat->num_recv, mat->recv_req, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
        ierr = jxf_MPI_Waitall(mat->num_send, mat->send_req, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
        if (timeFlag)
        {
            t3 = jxf_MPI_Wtime();
            mat->time[JXF_MATVEC_MPI_TIME] += (t3 - t2);
        }
        if (!jxf_commsOnly)
        {
            for (i = 0; i < m; i ++) recvbuf[i] = x[i];
            for (row = 0; row < m; row ++)
            {
                JXF_Int len = rp[row+1] - rp[row];
                JXF_Int *ind = cval + rp[row];
                JXF_Real *val = aval + rp[row];
                JXF_Real temp = 0.0;
                
                for (i = 0; i < len; i ++)
                {
                    temp += (val[i] * recvbuf[ind[i]]);
                }
                b[row] = temp;
            }
        }
        if (timeFlag)
        {
            t4 = jxf_MPI_Wtime();
            mat->time[JXF_MATVEC_TOTAL_TIME] += (t4 - t1);
            mat->time[JXF_MATVEC_TIME] += (t4 - t3);
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVec_omp"
void jxf_Mat_dhMatVec_omp( jxf_Mat_dh mat, JXF_Real *x, JXF_Real *b )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, i, row, m = mat->m;
    JXF_Int *rp = mat->rp, *cval = mat->cval;
    JXF_Real *aval = mat->aval;
    JXF_Int *sendind = mat->sendind;
    JXF_Int sendlen = mat->sendlen;
    JXF_Real *sendbuf = mat->sendbuf;
    JXF_Real *recvbuf = mat->recvbuf;
    JXF_Real t1 = 0, t2 = 0, t3 = 0, t4 = 0, tx = 0;
    JXF_Real *val, temp;
    JXF_Int len, *ind;
    jxf_bool timeFlag = mat->matvec_timing;
    
    if (timeFlag) t1 = jxf_MPI_Wtime();
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(runtime) private(i)
#endif
    for (i = 0; i < sendlen; i ++) sendbuf[i] = x[sendind[i]];
    if (timeFlag)
    {
        t2 = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_TIME] += (t2 - t1);
    }
    ierr = jxf_MPI_Startall(mat->num_recv, mat->recv_req); JXF_CHECK_MPI_V_ERROR(ierr);
    ierr = jxf_MPI_Startall(mat->num_send, mat->send_req); JXF_CHECK_MPI_V_ERROR(ierr);
    ierr = jxf_MPI_Waitall(mat->num_recv, mat->recv_req, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
    ierr = jxf_MPI_Waitall(mat->num_send, mat->send_req, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
    if (timeFlag)
    {
        t3 = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_MPI_TIME] += (t3 - t2);
    }
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(runtime) private(i)
#endif
    for (i = 0; i < m; i ++) recvbuf[i] = x[i];
    if (timeFlag)
    {
        tx = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_MPI_TIME2] += (tx - t1);
    }
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(runtime) private(row,i,len,ind,val,temp)
#endif
    for (row = 0; row < m; row ++)
    {
        len = rp[row+1] - rp[row];
        ind = cval + rp[row];
        val = aval + rp[row];
        temp = 0.0;
        for (i = 0; i < len; i ++)
        {
            temp += (val[i] * recvbuf[ind[i]]);
        }
        b[row] = temp;
    }
    if (timeFlag)
    {
        t4 = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_TOTAL_TIME] += (t4 - t1);
        mat->time[JXF_MATVEC_TIME] += (t4 - t3);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVec_uni_omp"
void jxf_Mat_dhMatVec_uni_omp( jxf_Mat_dh mat, JXF_Real *x, JXF_Real *b )
{
    JXF_START_FUNC_DH
    JXF_Int i, row, m = mat->m;
    JXF_Int *rp = mat->rp, *cval = mat->cval;
    JXF_Real *aval = mat->aval;
    JXF_Real t1 = 0, t2 = 0;
    jxf_bool timeFlag = mat->matvec_timing;
    
    if (timeFlag)
    {
        t1 = jxf_MPI_Wtime();
    }
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(runtime) private(row,i)
#endif
    for (row = 0; row < m; row ++)
    {
        JXF_Int len = rp[row+1] - rp[row];
        JXF_Int *ind = cval+rp[row];
        JXF_Real *val = aval+rp[row];
        JXF_Real temp = 0.0;
        
        for (i = 0; i < len; i ++)
        {
            temp += (val[i] * x[ind[i]]);
        }
        b[row] = temp;
    }
    if (timeFlag)
    {
        t2 = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_TIME] += (t2 - t1);
        mat->time[JXF_MATVEC_TOTAL_TIME] += (t2 - t1);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMatVec_uni"
void jxf_Mat_dhMatVec_uni( jxf_Mat_dh mat, JXF_Real *x, JXF_Real *b )
{
    JXF_START_FUNC_DH
    JXF_Int i, row, m = mat->m;
    JXF_Int *rp = mat->rp, *cval = mat->cval;
    JXF_Real *aval = mat->aval;
    JXF_Real t1 = 0, t2 = 0;
    jxf_bool timeFlag = mat->matvec_timing;
    
    if (timeFlag) t1 = jxf_MPI_Wtime();
    for (row = 0; row < m; row ++)
    {
        JXF_Int len = rp[row+1] - rp[row];
        JXF_Int *ind = cval + rp[row];
        JXF_Real *val = aval + rp[row];
        JXF_Real temp = 0.0;
        
        for (i = 0; i < len; i ++)
        {
            temp += (val[i] * x[ind[i]]);
        }
        b[row] = temp;
    }
    if (timeFlag)
    {
        t2 = jxf_MPI_Wtime();
        mat->time[JXF_MATVEC_TIME] += (t2 - t1);
        mat->time[JXF_MATVEC_TOTAL_TIME] += (t2 - t1);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhReadNz"
JXF_Int jxf_Mat_dhReadNz( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, retval = mat->rp[mat->m];
    JXF_Int nz = retval;
    
    ierr = jxf_MPI_Allreduce(&nz, &retval, 1, JXF_MPI_INT, MPI_SUM, jxf_comm_dh); JXF_CHECK_MPI_ERROR(ierr);
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhZeroTiming"
void jxf_Mat_dhZeroTiming( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
    for (i = 0; i < JXF_MAT_DH_BINS; ++ i)
    {
        mat->time[i] = 0;
        mat->time_max[i] = 0;
        mat->time_min[i] = 0;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhReduceTiming"
void jxf_Mat_dhReduceTiming( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    if (mat->time[JXF_MATVEC_MPI_TIME])
    {
        mat->time[JXF_MATVEC_RATIO] = mat->time[JXF_MATVEC_TIME] / mat->time[JXF_MATVEC_MPI_TIME];
    }
    jxf_MPI_Allreduce(mat->time, mat->time_min, JXF_MAT_DH_BINS, JXF_MPI_REAL, MPI_MIN, jxf_comm_dh);
    jxf_MPI_Allreduce(mat->time, mat->time_max, JXF_MAT_DH_BINS, JXF_MPI_REAL, MPI_MAX, jxf_comm_dh);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPermute"
void jxf_Mat_dhPermute( jxf_Mat_dh A, JXF_Int *n2o, jxf_Mat_dh *Bout )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh B;
    JXF_Int  i, j, *RP = A->rp, *CVAL = A->cval;
    JXF_Int  *o2n, *rp, *cval, m = A->m, nz = RP[m];
    JXF_Real *aval, *AVAL = A->aval;
    
    jxf_Mat_dhCreate(&B); JXF_CHECK_V_ERROR;
    B->m = B->n = m;
   *Bout = B;
    o2n = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) o2n[n2o[i]] = i;
    rp = B->rp = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = B->cval = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    aval = B->aval = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    rp[0] = 0;
    for (i = 0; i < m; ++ i)
    {
        JXF_Int oldRow = n2o[i];
        rp[i+1] = RP[oldRow+1] - RP[oldRow];
    }
    for (i = 1; i <= m; ++ i) rp[i] = rp[i] + rp[i-1];
    for (i = 0; i < m; ++ i)
    {
        JXF_Int oldRow = n2o[i];
        JXF_Int idx = rp[i];
        
        for (j = RP[oldRow]; j < RP[oldRow+1]; ++ j)
        {
            cval[idx] = o2n[CVAL[j]];
            aval[idx] = AVAL[j];
            ++ idx;
        }
    }
    JXF_FREE_DH(o2n); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintGraph"
void jxf_Mat_dhPrintGraph( jxf_Mat_dh A, jxf_SubdomainGraph_dh sg, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int pe, id = jxf_myid_dh;
    JXF_Int ierr;
    
    if (sg != NULL)
    {
        id = sg->o2n_sub[id];
    }
    for (pe = 0; pe < jxf_np_dh; ++ pe)
    {
        ierr = jxf_MPI_Barrier(jxf_comm_dh); JXF_CHECK_MPI_V_ERROR(ierr);
        if (id == pe)
        {
            if (sg == NULL)
            {
                jxf_mat_dh_print_graph_private(A->m, A->beg_row,
                          A->rp, A->cval, A->aval, NULL, NULL, NULL, fp); JXF_CHECK_V_ERROR;
            }
            else
            {
                JXF_Int beg_row = sg->beg_rowP[jxf_myid_dh];
                
                jxf_mat_dh_print_graph_private(A->m, beg_row,
                      A->rp, A->cval, A->aval, sg->n2o_row, sg->o2n_col, sg->o2n_ext, fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintRows"
void jxf_Mat_dhPrintRows( jxf_Mat_dh A, jxf_SubdomainGraph_dh sg, FILE *fp )
{
    JXF_START_FUNC_DH
    jxf_bool noValues;
    JXF_Int m = A->m, *rp = A->rp, *cval = A->cval;
    JXF_Real *aval = A->aval;
    
    noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
    if (noValues) aval = NULL;
    if (sg == NULL)
    {
        JXF_Int i, j;
        JXF_Int beg_row = A->beg_row;
        
        jxf_fprintf(fp, "\n----- A, unpermuted ------------------------------------\n");
        for (i = 0; i < m; ++ i)
        {
            jxf_fprintf(fp, "%i :: ", 1+i+beg_row);
            for (j = rp[i]; j < rp[i+1]; ++ j)
            {
                if (noValues)
                {
                    jxf_fprintf(fp, "%i ", 1+cval[j]);
                }
                else
                {
                    jxf_fprintf(fp, "%i,%g ; ", 1+cval[j], aval[j]);
                }
            }
            jxf_fprintf(fp, "\n");
        }
    }
    else if (jxf_np_dh == 1)
    {
        JXF_Int i, k, idx = 1;
        JXF_Int oldRow;
        
        for (i = 0; i < sg->blocks; ++ i)
        {
            JXF_Int oldBlock = sg->n2o_sub[i];
            JXF_Int beg_row = sg->beg_row[oldBlock];
            JXF_Int end_row = beg_row + sg->row_count[oldBlock];
            
            jxf_fprintf(fp, "\n");
            jxf_fprintf(fp, "\n----- A, permuted, single mpi task  ------------------\n");
            jxf_fprintf(fp, "---- new subdomain: %i;  old subdomain: %i\n", i, oldBlock);
            jxf_fprintf(fp, "     old beg_row:   %i;  new beg_row:   %i\n",
                                       sg->beg_row[oldBlock], sg->beg_rowP[oldBlock]);
            jxf_fprintf(fp, "     local rows in this block: %i\n", sg->row_count[oldBlock]);
            jxf_fprintf(fp, "     bdry rows in this block:  %i\n", sg->bdry_count[oldBlock]);
            jxf_fprintf(fp, "     1st bdry row= %i \n", 1+end_row-sg->bdry_count[oldBlock]);
            for (oldRow = beg_row; oldRow < end_row; ++ oldRow)
            {
                JXF_Int len = 0, *cval;
                JXF_Real *aval;
                
                jxf_fprintf(fp, "%3i (old= %3i) :: ", idx, 1+oldRow);
                ++ idx;
                jxf_Mat_dhGetRow(A, oldRow, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                for (k = 0; k < len; ++ k)
                {
                    if (noValues)
                    {
                        jxf_fprintf(fp, "%i ", 1+sg->o2n_col[cval[k]]);
                    }
                    else
                    {
                        jxf_fprintf(fp, "%i,%g ; ", 1+sg->o2n_col[cval[k]], aval[k]);
                    }
                }
                jxf_fprintf(fp, "\n");
                jxf_Mat_dhRestoreRow(A, oldRow, &len, &cval, &aval); JXF_CHECK_V_ERROR;
            }
        }
    }
    else
    {
        jxf_Hash_i_dh hash = sg->o2n_ext;
        JXF_Int *o2n_col = sg->o2n_col, *n2o_row = sg->n2o_row;
        JXF_Int beg_row = sg->beg_row[jxf_myid_dh];
        JXF_Int beg_rowP = sg->beg_rowP[jxf_myid_dh];
        JXF_Int i, j;
        
        for (i = 0; i < m; ++ i)
        {
            JXF_Int row = n2o_row[i];
            
            jxf_fprintf(fp, "%3i (old= %3i) :: ", 1+i+beg_rowP, 1+row+beg_row);
            for (j = rp[row]; j < rp[row+1]; ++ j)
            {
                JXF_Int col = cval[j];
                
                if (col >= beg_row && col < beg_row+m)
                {
                    col = o2n_col[col-beg_row] + beg_rowP;
                }
                else
                {
                    JXF_Int tmp = col;
                    
                    tmp = jxf_Hash_i_dhLookup(hash, col); JXF_CHECK_V_ERROR;
                    if (tmp == -1)
                    {
                        jxf_sprintf(jxf_msgBuf_dh, "nonlocal column= %i not in hash table", 1+col);
                        JXF_SET_V_ERROR(jxf_msgBuf_dh);
                    }
                    else
                    {
                        col = tmp;
                    }
                }
                if (noValues)
                {
                    jxf_fprintf(fp, "%i ", 1+col);
                }
                else
                {
                    jxf_fprintf(fp, "%i,%g ; ", 1+col, aval[j]);
                }
            }
            jxf_fprintf(fp, "\n");
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintTriples"
void jxf_Mat_dhPrintTriples( jxf_Mat_dh A, jxf_SubdomainGraph_dh sg, char *filename )
{
    JXF_START_FUNC_DH
    JXF_Int m = A->m, *rp = A->rp, *cval = A->cval;
    JXF_Real *aval = A->aval;
    jxf_bool noValues;
    jxf_bool matlab;
    FILE *fp;
    
    noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
    if (noValues) aval = NULL;
    matlab = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-matlab"));
    if (sg == NULL)
    {
        JXF_Int i, j, pe;
        JXF_Int beg_row = A->beg_row;
        JXF_Real val;
        
        for (pe = 0; pe < jxf_np_dh; ++ pe)
        {
            jxf_MPI_Barrier(jxf_comm_dh);
            if (pe == jxf_myid_dh)
            {
                if (pe == 0)
                {
                    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
                }
                else
                {
                    fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
                }
                for (i = 0; i < m; ++ i)
                {
                    for (j = rp[i]; j < rp[i+1]; ++ j)
                    {
                        if (noValues)
                        {
                            jxf_fprintf(fp, "%i %i\n", 1+i+beg_row, 1+cval[j]);
                        }
                        else
                        {
                            val = aval[j];
                            if (val == 0.0 && matlab) val = _JXF_MATLAB_ZERO_;
                            jxf_fprintf(fp, JXF_TRIPLES_FORMAT, 1+i+beg_row, 1+cval[j], val);
                        }
                    }
                }
                jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    else if (jxf_np_dh == 1)
    {
        JXF_Int i, j, k, idx = 1;
        
        fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
        for (i = 0; i < sg->blocks; ++ i)
        {
            JXF_Int oldBlock = sg->n2o_sub[i];
            JXF_Int beg_row = sg->beg_rowP[oldBlock];
            JXF_Int end_row = beg_row + sg->row_count[oldBlock];
            
            for (j = beg_row; j < end_row; ++ j)
            {
                JXF_Int len = 0, *cval;
                JXF_Real *aval;
                JXF_Int oldRow = sg->n2o_row[j];
                
                jxf_Mat_dhGetRow(A, oldRow, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                if (noValues)
                {
                    for (k = 0; k < len; ++ k)
                    {
                        jxf_fprintf(fp, "%i %i\n", idx, 1+sg->o2n_col[cval[k]]);
                    }
                    ++ idx;
                }
                else
                {
                    for (k = 0; k < len; ++ k)
                    {
                        JXF_Real val = aval[k];
                        
                        if (val == 0.0 && matlab) val = _JXF_MATLAB_ZERO_;
                        jxf_fprintf(fp, JXF_TRIPLES_FORMAT, idx, 1+sg->o2n_col[cval[k]], val);
                    }
                    ++ idx;
                }
                jxf_Mat_dhRestoreRow(A, oldRow, &len, &cval, &aval); JXF_CHECK_V_ERROR;
            }
        }
    }
    else
    {
        jxf_Hash_i_dh hash = sg->o2n_ext;
        JXF_Int *o2n_col = sg->o2n_col, *n2o_row = sg->n2o_row;
        JXF_Int beg_row = sg->beg_row[jxf_myid_dh];
        JXF_Int beg_rowP = sg->beg_rowP[jxf_myid_dh];
        JXF_Int i, j, pe;
        JXF_Int id = sg->o2n_sub[jxf_myid_dh];
        
        for (pe = 0; pe < jxf_np_dh; ++ pe)
        {
            jxf_MPI_Barrier(jxf_comm_dh);
            if (id == pe)
            {
                if (pe == 0)
                {
                    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
                }
                else
                {
                    fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
                }
                for (i = 0; i < m; ++ i)
                {
                    JXF_Int row = n2o_row[i];
                    
                    for (j = rp[row]; j < rp[row+1]; ++ j)
                    {
                        JXF_Int col = cval[j];
                        JXF_Real val = 0.0;
                        
                        if (aval != NULL) val = aval[j];
                        if (val == 0.0 && matlab) val = _JXF_MATLAB_ZERO_;
                        if (col >= beg_row && col < beg_row+m)
                        {
                            col = o2n_col[col-beg_row] + beg_rowP;
                        }
                        else
                        {
                            JXF_Int tmp = col;
                            tmp = jxf_Hash_i_dhLookup(hash, col); JXF_CHECK_V_ERROR;
                            if (tmp == -1)
                            {
                                jxf_sprintf(jxf_msgBuf_dh, "nonlocal column= %i not in hash table", 1+col);
                                JXF_SET_V_ERROR(jxf_msgBuf_dh);
                            }
                            else
                            {
                                col = tmp;
                            }
                        }
                        if (noValues)
                        {
                            jxf_fprintf(fp, "%i %i\n", 1+i+beg_rowP, 1+col);
                        }
                        else
                        {
                            jxf_fprintf(fp, JXF_TRIPLES_FORMAT, 1+i+beg_rowP, 1+col, val);
                        }
                    }
                }
                jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintCSR"
void jxf_Mat_dhPrintCSR( jxf_Mat_dh A, jxf_SubdomainGraph_dh sg, char *filename )
{
    JXF_START_FUNC_DH
    FILE *fp;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single mpi task");
    }
    if (sg != NULL)
    {
        JXF_SET_V_ERROR("not implemented for reordered matrix (jxf_SubdomainGraph_dh should be NULL)");
    }
    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
    if (sg == NULL)
    {
        jxf_mat_dh_print_csr_private(A->m, A->rp, A->cval, A->aval, fp); JXF_CHECK_V_ERROR;
    }
    else
    {
        jxf_mat_dh_print_csr_private(A->m, A->rp, A->cval, A->aval, fp); JXF_CHECK_V_ERROR;
    }
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintBIN"
void jxf_Mat_dhPrintBIN( jxf_Mat_dh A, jxf_SubdomainGraph_dh sg, char *filename )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    if (sg != NULL)
    {
        JXF_SET_V_ERROR("not implemented for reordering; ensure sg=NULL");
    }
    jxf_io_dh_print_ebin_mat_private(A->m, A->beg_row,
             A->rp, A->cval, A->aval, NULL, NULL, NULL, filename); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhReadCSR"
void jxf_Mat_dhReadCSR( jxf_Mat_dh *mat, char *filename )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh A;
    FILE *fp;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    fp = jxf_openFile_dh(filename, "r"); JXF_CHECK_V_ERROR;
    jxf_Mat_dhCreate(&A); JXF_CHECK_V_ERROR;
    jxf_mat_dh_read_csr_private(&A->m, &A->rp, &A->cval, &A->aval, fp); JXF_CHECK_V_ERROR;
    A->n = A->m;
   *mat = A;
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhReadTriples"
void jxf_Mat_dhReadTriples( jxf_Mat_dh *mat, JXF_Int ignore, char *filename )
{
    JXF_START_FUNC_DH
    FILE *fp = NULL;
    jxf_Mat_dh A = NULL;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    fp = jxf_openFile_dh(filename, "r"); JXF_CHECK_V_ERROR;
    jxf_Mat_dhCreate(&A); JXF_CHECK_V_ERROR;
    jxf_mat_dh_read_triples_private(ignore, &A->m, &A->rp, &A->cval, &A->aval, fp); JXF_CHECK_V_ERROR;
    A->n = A->m;
   *mat = A;
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhReadBIN"
void jxf_Mat_dhReadBIN( jxf_Mat_dh *mat, char *filename )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh A;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    jxf_Mat_dhCreate(&A); JXF_CHECK_V_ERROR;
    jxf_io_dh_read_ebin_mat_private(&A->m, &A->rp, &A->cval, &A->aval, filename); JXF_CHECK_V_ERROR;
    A->n = A->m;
   *mat = A;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhTranspose"
void jxf_Mat_dhTranspose( jxf_Mat_dh A, jxf_Mat_dh *Bout )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh B;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only for sequential");
    }
    jxf_Mat_dhCreate(&B); JXF_CHECK_V_ERROR;
   *Bout = B;
    B->m = B->n = A->m;
    jxf_mat_dh_transpose_private(A->m, A->rp, &B->rp, A->cval, &B->cval, A->aval, &B->aval); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhMakeStructurallySymmetric"
void jxf_Mat_dhMakeStructurallySymmetric( jxf_Mat_dh A )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only for sequential");
    }
    jxf_make_symmetric_private(A->m, &A->rp, &A->cval, &A->aval); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

void jxf_insert_diags_private( jxf_Mat_dh A, JXF_Int ct );

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhFixDiags"
void jxf_Mat_dhFixDiags( jxf_Mat_dh A )
{
    JXF_START_FUNC_DH
    JXF_Int i, j;
    JXF_Int *rp = A->rp, *cval = A->cval, m = A->m;
    jxf_bool ct = 0;
    JXF_Real *aval = A->aval;
    
    for (i = 0; i < m; ++ i)
    {
        jxf_bool flag = jxf_true;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int col = cval[j];
            
            if (col == i)
            {
                flag = jxf_false;
                break;
            }
        }
        if (flag) ++ ct;
    }
    if (ct)
    {
        jxf_printf("\njxf_Mat_dhFixDiags:: %i diags not explicitly present; inserting!\n", ct);
        jxf_insert_diags_private(A, ct); JXF_CHECK_V_ERROR;
        rp = A->rp;
        cval = A->cval;
        aval = A->aval;
    }
    for (i = 0; i < m; ++ i)
    {
        JXF_Real sum = 0.0;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            sum += fabs(aval[j]);
        }
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            if (cval[j] == i)
            {
                aval[j] = sum;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_insert_diags_private"
void jxf_insert_diags_private( jxf_Mat_dh A, JXF_Int ct )
{
    JXF_START_FUNC_DH
    JXF_Int *RP = A->rp, *CVAL = A->cval;
    JXF_Int *rp, *cval, m = A->m;
    JXF_Real *aval, *AVAL = A->aval;
    JXF_Int nz = RP[m] + ct;
    JXF_Int i, j, idx = 0;
    
    rp = A->rp = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = A->cval = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    aval = A->aval = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    rp[0] = 0;
    for (i = 0; i < m; ++ i)
    {
        jxf_bool flag = jxf_true;
        
        for (j = RP[i]; j < RP[i+1]; ++ j)
        {
            cval[idx] = CVAL[j];
            aval[idx] = AVAL[j];
            ++ idx;
            if (CVAL[j] == i) flag = jxf_false;
        }
        if (flag)
        {
            cval[idx] = i;
            aval[idx] = 0.0;
            ++ idx;
        }
        rp[i+1] = idx;
    }
    JXF_FREE_DH(RP); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(CVAL); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(AVAL); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPrintDiags"
void jxf_Mat_dhPrintDiags( jxf_Mat_dh A, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, m = A->m;
    JXF_Int *rp = A->rp, *cval = A->cval;
    JXF_Real *aval = A->aval;
    
    jxf_fprintf(fp, "=================== diagonal elements ====================\n");
    for (i = 0; i < m; ++ i)
    {
        jxf_bool flag = jxf_true;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            if (cval[j] == i)
            {
                jxf_fprintf(fp, "%i  %g\n", i+1, aval[j]);
                flag = jxf_false;
                break;
            }
        }
        if (flag)
        {
            jxf_fprintf(fp, "%i  ---------- missing\n", i+1);
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhGetRow"
void jxf_Mat_dhGetRow( jxf_Mat_dh B, JXF_Int globalRow, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_Int row = globalRow - B->beg_row;
    
    if (row > B->m)
    {
        jxf_sprintf(jxf_msgBuf_dh, "requested globalRow= %i, which is local row= %i, but only have %i rows!",
                                                         globalRow, row, B->m);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
   *len = B->rp[row+1] - B->rp[row];
    if (ind != NULL) *ind = B->cval + B->rp[row];
    if (val != NULL) *val = B->aval + B->rp[row];
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhRestoreRow"
void jxf_Mat_dhRestoreRow( jxf_Mat_dh B, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhRowPermute"
void jxf_Mat_dhRowPermute( jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    if (jxf_ignoreMe) JXF_SET_V_ERROR("turned off; compilation problem on blue");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPartition"
void jxf_build_adj_lists_private( jxf_Mat_dh mat, JXF_Int **rpOUT, JXF_Int **cvalOUT )
{
    JXF_START_FUNC_DH
    JXF_Int m = mat->m;
    JXF_Int *RP = mat->rp, *CVAL = mat->cval;
    JXF_Int nz = RP[m];
    JXF_Int i, j, *rp, *cval, idx = 0;
    
    rp = *rpOUT = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = *cvalOUT = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rp[0] = 0;
    for (i = 0; i < m; ++ i)
    {
        for (j = RP[i]; j < RP[i+1]; ++ j)
        {
            JXF_Int col = CVAL[j];
            
            if (col != i)
            {
                cval[idx++] = col;
            }
        }
        rp[i+1] = idx;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mat_dhPartition"
void jxf_Mat_dhPartition( jxf_Mat_dh mat, JXF_Int blocks, JXF_Int **beg_rowOUT, JXF_Int **row_countOUT,  JXF_Int **n2oOUT, JXF_Int **o2nOUT )
{
    JXF_START_FUNC_DH
#ifndef HAVE_METIS_DH
    if (jxf_ignoreMe) JXF_SET_V_ERROR("not compiled for metis!");
#else
    JXF_Int *beg_row, *row_count, *n2o, *o2n, bk, new, *part;
    JXF_Int m = mat->m;
    JXF_Int i, cutEdgeCount;
    JXF_Real zero = 0.0;
    JXF_Int metisOpts[5] = {0, 0, 0, 0, 0};
    JXF_Int *rp, *cval;
    
    beg_row = *beg_rowOUT = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    row_count = *row_countOUT = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
   *n2oOUT = n2o = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
   *o2nOUT = o2n = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    jxf_build_adj_lists_private(mat, &rp, &cval); JXF_CHECK_V_ERROR;
    part = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    METIS_PartGraphKway(&m, rp, cval, NULL, NULL, &zero, &zero, &blocks, metisOpts, &cutEdgeCount, part);
    JXF_FREE_DH(rp); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(cval); JXF_CHECK_V_ERROR;
    if (mat->debug)
    {
        jxf_printf_dh("\nmetis partitioning vector; blocks= %i\n", blocks);
        for (i = 0; i < m; ++ i) jxf_printf_dh("  %i %i\n", i+1, part[i]);
    }
    for (i = 0; i < blocks; ++ i) row_count[i] = 0;
    for (i = 0; i < m; ++ i)
    {
        bk = part[i];
        row_count[bk] += 1;
    }
    beg_row[0] = 0;
    for (i = 1; i < blocks; ++ i) beg_row[i] = beg_row[i-1] + row_count[i-1];
    if (mat->debug)
    {
        jxf_printf_dh("\nrow_counts: ");
        for (i = 0; i < blocks; ++ i) jxf_printf_dh(" %i", row_count[i]);
        jxf_printf_dh("\nbeg_row: ");
        for (i = 0; i < blocks; ++ i) jxf_printf_dh(" %i", beg_row[i]+1);
        jxf_printf_dh("\n");
    }
    JXF_Int *tmp = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    memcpy(tmp, beg_row, blocks*sizeof(JXF_Int));
    for (i = 0; i < m; ++ i)
    {
        bk = part[i];
        new = tmp[bk];
        tmp[bk] += 1;
        o2n[i] = new;
        n2o[new] = i;
    }
    JXF_FREE_DH(tmp);
    JXF_FREE_DH(part); JXF_CHECK_V_ERROR;
#endif
    JXF_END_FUNC_DH
}
