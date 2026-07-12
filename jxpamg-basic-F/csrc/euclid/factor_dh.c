//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  factor_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

void jxf_Factor_dh_junk()
{
}

static void jxf_adjust_bj_private( jxf_Factor_dh mat );
static void jxf_unadjust_bj_private( jxf_Factor_dh mat );

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhCreate"
void jxf_Factor_dhCreate( jxf_Factor_dh *mat )
{
    JXF_START_FUNC_DH
    struct _jxf_factor_dh *tmp;
    if (jxf_np_dh > JXF_MAX_MPI_TASKS)
    {
        JXF_SET_V_ERROR("you must change JXF_MAX_MPI_TASKS and recompile!");
    }
    tmp = (struct _jxf_factor_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_factor_dh)); JXF_CHECK_V_ERROR;
   *mat = tmp;
    tmp->m = 0;
    tmp->n = 0;
    tmp->id = jxf_myid_dh;
    tmp->beg_row = 0;
    tmp->first_bdry = 0;
    tmp->bdry_count = 0;
    tmp->blockJacobi = jxf_false;
    tmp->rp = NULL;
    tmp->cval = NULL;
    tmp->aval = NULL;
    tmp->fill = NULL;
    tmp->diag = NULL;
    tmp->alloc = 0;
    tmp->work_y_lo = tmp->work_x_hi = NULL;
    tmp->sendbufLo = tmp->sendbufHi = NULL;
    tmp->sendindLo = tmp->sendindHi = NULL;
    tmp->num_recvLo = tmp->num_recvHi = 0;
    tmp->num_sendLo = tmp->num_sendHi = 0;
    tmp->sendlenLo = tmp->sendlenHi = 0;
    tmp->solveIsSetup = jxf_false;
    tmp->numbSolve = NULL;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_Factor");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhDestroy"
void jxf_Factor_dhDestroy( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    if (mat->rp != NULL)
    {
        JXF_FREE_DH(mat->rp); JXF_CHECK_V_ERROR;
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
    if (mat->work_y_lo != NULL)
    {
        JXF_FREE_DH(mat->work_y_lo); JXF_CHECK_V_ERROR;
    }
    if (mat->work_x_hi != NULL)
    {
        JXF_FREE_DH(mat->work_x_hi); JXF_CHECK_V_ERROR;
    }
    if (mat->sendbufLo != NULL)
    {
        JXF_FREE_DH(mat->sendbufLo); JXF_CHECK_V_ERROR;
    }
    if (mat->sendbufHi != NULL)
    {
        JXF_FREE_DH(mat->sendbufHi); JXF_CHECK_V_ERROR;
    }
    if (mat->sendindLo != NULL)
    {
        JXF_FREE_DH(mat->sendindLo); JXF_CHECK_V_ERROR;
    }
    if (mat->sendindHi != NULL)
    {
        JXF_FREE_DH(mat->sendindHi); JXF_CHECK_V_ERROR;
    }
    if (mat->numbSolve != NULL)
    {
        jxf_Numbering_dhDestroy(mat->numbSolve); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(mat); JXF_CHECK_V_ERROR; 
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_create_fake_mat_private"
static void jxf_create_fake_mat_private( jxf_Factor_dh mat, jxf_Mat_dh *matFakeIN )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh matFake;
    
    jxf_Mat_dhCreate(matFakeIN); JXF_CHECK_V_ERROR;
    matFake = *matFakeIN;
    matFake->m = mat->m;
    matFake->n = mat->n;
    matFake->rp = mat->rp;
    matFake->cval = mat->cval;
    matFake->aval = mat->aval;
    matFake->beg_row = mat->beg_row;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_destroy_fake_mat_private"
static void jxf_destroy_fake_mat_private( jxf_Mat_dh matFake )
{
    JXF_START_FUNC_DH
    matFake->rp = NULL;
    matFake->cval = NULL;
    matFake->aval = NULL;
    jxf_Mat_dhDestroy(matFake); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhReadNz"
JXF_Int jxf_Factor_dhReadNz( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, retval = mat->rp[mat->m];
    JXF_Int nz = retval;
    
    ierr = jxf_MPI_Allreduce(&nz, &retval, 1, JXF_MPI_INT, MPI_SUM, jxf_comm_dh); JXF_CHECK_MPI_ERROR(ierr);
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhPrintRows"
void jxf_Factor_dhPrintRows( jxf_Factor_dh mat, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int beg_row = mat->beg_row;
    JXF_Int m = mat->m, i, j;
    jxf_bool noValues;
    
    noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
    if (mat->aval == NULL) noValues = jxf_true;
    if (mat->blockJacobi)
    {
        jxf_adjust_bj_private(mat); JXF_CHECK_V_ERROR;
    }
    jxf_fprintf(fp, "\n----------------------- jxf_Factor_dhPrintRows ------------------\n");
    if (mat->blockJacobi)
    {
        jxf_fprintf(fp, "@@@ Block Jacobi ILU; adjusted values from zero-based @@@\n");
    }
    for (i = 0; i < m; ++ i)
    {
        jxf_fprintf(fp, "%i :: ", 1+i+beg_row);
        for (j = mat->rp[i]; j < mat->rp[i+1]; ++ j)
        {
            if (noValues)
            {
                jxf_fprintf(fp, "%i ", 1+mat->cval[j]);
            }
            else
            {
                jxf_fprintf(fp, "%i,%g ; ", 1+mat->cval[j], mat->aval[j]);
            }
        }
        jxf_fprintf(fp, "\n");
    }
    if (mat->blockJacobi)
    {
        jxf_unadjust_bj_private(mat); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhPrintDiags"
void jxf_Factor_dhPrintDiags( jxf_Factor_dh mat, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int beg_row = mat->beg_row;
    JXF_Int m = mat->m, i, pe, *diag = mat->diag;
    JXF_REAL_DH *aval = mat->aval;
    
    jxf_fprintf_dh(fp, "\n----------------------- jxf_Factor_dhPrintDiags ------------------\n");
    jxf_fprintf_dh(fp, "(grep for 'ZERO')\n");
    for (pe = 0; pe < jxf_np_dh; ++ pe)
    {
        jxf_MPI_Barrier(jxf_comm_dh);
        if (mat->id == pe)
        {
            jxf_fprintf(fp, "----- subdomain: %i  processor: %i\n", pe, jxf_myid_dh);
            for (i = 0; i < m; ++ i)
            {
                JXF_REAL_DH val = aval[diag[i]];
                if (val)
                {
                    jxf_fprintf(fp, "%i %g\n", i+1+beg_row, aval[diag[i]]);
                }
                else
                {
                    jxf_fprintf(fp, "%i %g ZERO\n", i+1+beg_row, aval[diag[i]]);
                }
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhPrintGraph"
void jxf_Factor_dhPrintGraph( jxf_Factor_dh mat, char *filename )
{
    JXF_START_FUNC_DH
    FILE *fp;
    JXF_Int i, j, m = mat->m, *work, *rp = mat->rp, *cval = mat->cval;
    
    if (jxf_np_dh > 1) JXF_SET_V_ERROR("only implemented for single mpi task");
    work = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        for (j = 0; j < m; ++ j) work[j] = 0;
        for (j = rp[i]; j < rp[i]; ++ j) work[cval[j]] = 1;
        for (j = 0; j < m; ++ j)
        {
            if (work[j])
            {
                jxf_fprintf(fp, " x ");
            }
            else
            {
                jxf_fprintf(fp, "   ");
            }
        }
        jxf_fprintf(fp, "\n");
    }
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(work);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhPrintTriples"
void jxf_Factor_dhPrintTriples( jxf_Factor_dh mat, char *filename )
{
    JXF_START_FUNC_DH
    JXF_Int pe, i, j;
    JXF_Int m = mat->m, *rp = mat->rp;
    JXF_Int beg_row = mat->beg_row;
    JXF_REAL_DH *aval = mat->aval;
    jxf_bool noValues;
    FILE *fp;
    
    if (mat->blockJacobi)
    {
        jxf_adjust_bj_private(mat); JXF_CHECK_V_ERROR;
    }
    noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
    if (noValues) aval = NULL;
    for (pe = 0; pe < jxf_np_dh; ++ pe)
    {
        jxf_MPI_Barrier(jxf_comm_dh);
        if (mat->id == pe)
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
                        jxf_fprintf(fp, "%i %i\n", 1+i+beg_row, 1+mat->cval[j]);
                    }
                    else
                    {
                        jxf_fprintf(fp, JXF_TRIPLES_FORMAT, 1+i+beg_row, 1+mat->cval[j], aval[j]);
                    }
                }
            }
            jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
        }
    }
    if (mat->blockJacobi)
    {
        jxf_unadjust_bj_private(mat); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_setup_receives_private"
static JXF_Int jxf_setup_receives_private( jxf_Factor_dh mat,
                                   JXF_Int *beg_rows,
                                   JXF_Int *end_rows,
                                   JXF_Real *recvBuf,
                                   MPI_Request *req,
                                   JXF_Int *reqind,
                                   JXF_Int reqlen,
                                   JXF_Int *outlist,
                                   jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, this_pe, num_recv = 0;
    MPI_Request request;
    
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nFACT ========================================================\n");
        jxf_fprintf(jxf_logFile, "FACT STARTING: jxf_setup_receives_private\n");
    }
    for (i = 0; i < reqlen; i = j)
    {
        this_pe = jxf_mat_jxf_find_owner(beg_rows, end_rows, reqind[i]); JXF_CHECK_ERROR(-1);
        for (j = i+1; j < reqlen; j ++)
        {
            JXF_Int idx = reqind[j];
            
            if (idx < beg_rows[this_pe] || idx >= end_rows[this_pe])
            {
                break;
            }
        }
        if (debug)
        {
            JXF_Int k;
            
            jxf_fprintf(jxf_logFile, "FACT need nodes from P_%i: ", this_pe);
            for (k = i; k < j; ++ k) jxf_fprintf(jxf_logFile, "%i ", 1+reqind[k]);
            jxf_fprintf(jxf_logFile,"\n");
        }
        outlist[this_pe] = j - i;
        jxf_MPI_Isend(reqind+i, j-i, JXF_MPI_INT, this_pe, 444, jxf_comm_dh, &request);
        jxf_MPI_Request_free(&request);
        jxf_MPI_Recv_init(recvBuf+i, j-i, JXF_MPI_REAL, this_pe, 555, jxf_comm_dh, req+num_recv);
        ++ num_recv;
    }
    JXF_END_FUNC_VAL(num_recv);
}

#undef __FUNC__
#define __FUNC__ "jxf_setup_sends_private"
static void jxf_setup_sends_private( jxf_Factor_dh mat, JXF_Int *inlist, JXF_Int *o2n_subdomain, jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int i, jLo, jHi, sendlenLo, sendlenHi, first = mat->beg_row;
    MPI_Request *requests = mat->requests, *sendReq;
    MPI_Status *statuses = mat->status;
    jxf_bool isHigher;
    JXF_Int *rcvBuf;
    JXF_Real *sendBuf;
    JXF_Int myidNEW = o2n_subdomain[jxf_myid_dh];
    JXF_Int count;
    
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "FACT \nSTARTING: jxf_setup_sends_private\n");
    }
    sendlenLo = sendlenHi = 0;
    for (i = 0; i < jxf_np_dh; i ++)
    {
        if (inlist[i])
        {
            if (o2n_subdomain[i] < myidNEW)
            {
                sendlenLo += inlist[i];
            }
            else
            {
                sendlenHi += inlist[i];
            }
        }
    }
    mat->sendlenLo = sendlenLo;
    mat->sendlenHi = sendlenHi;
    mat->sendbufLo = (JXF_Real *)JXF_MALLOC_DH(sendlenLo * sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    mat->sendbufHi = (JXF_Real *)JXF_MALLOC_DH(sendlenHi * sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    mat->sendindLo = (JXF_Int *)JXF_MALLOC_DH(sendlenLo * sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    mat->sendindHi = (JXF_Int *)JXF_MALLOC_DH(sendlenHi * sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    count = 0;
    jLo = jHi = 0;
    mat->num_sendLo = 0;
    mat->num_sendHi = 0;
    for (i = 0; i < jxf_np_dh; i ++)
    {
        if (inlist[i])
        {
            isHigher = (o2n_subdomain[i] < myidNEW) ? jxf_false : jxf_true;
            if (isHigher)
            {
                rcvBuf = &mat->sendindHi[jHi];
                sendBuf = &mat->sendbufHi[jHi];
                sendReq = &mat->send_reqHi[mat->num_sendHi];
                mat->num_sendHi ++;
                jHi += inlist[i];
            }
            else
            {
                rcvBuf = &mat->sendindLo[jLo];
                sendBuf = &mat->sendbufLo[jLo];
                sendReq = &mat->send_reqLo[mat->num_sendLo];
                mat->num_sendLo ++;
                jLo += inlist[i];
            }
            jxf_MPI_Irecv(rcvBuf, inlist[i], JXF_MPI_INT, i, 444, jxf_comm_dh, requests+count);
            ++ count;
            jxf_MPI_Send_init(sendBuf, inlist[i], JXF_MPI_REAL, i, 555, jxf_comm_dh, sendReq);
        }
    }
    jxf_MPI_Waitall(count, requests, statuses);
    if (debug)
    {
        JXF_Int j;
        
        jLo = jHi = 0;
        jxf_fprintf(jxf_logFile, "\nFACT columns that I must send to other subdomains:\n");
        for (i = 0; i < jxf_np_dh; i ++)
        {
            if (inlist[i])
            {
                isHigher = (o2n_subdomain[i] < myidNEW) ? jxf_false : jxf_true;
                if (isHigher)
                {
                    rcvBuf = &mat->sendindHi[jHi];
                    jHi += inlist[i];
                }
                else
                {
                    rcvBuf = &mat->sendindLo[jLo];
                    jLo += inlist[i];
                }
                jxf_fprintf(jxf_logFile, "FACT  send to P_%i: ", i);
                for (j = 0; j < inlist[i]; ++ j) jxf_fprintf(jxf_logFile, "%i ", rcvBuf[j]+1);
                jxf_fprintf(jxf_logFile, "\n");
            }
        }
    }
    for (i = 0; i < mat->sendlenLo; i ++) mat->sendindLo[i] -= first;
    for (i = 0; i < mat->sendlenHi; i ++) mat->sendindHi[i] -= first;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhSolveSetup"
void jxf_Factor_dhSolveSetup( jxf_Factor_dh mat, jxf_SubdomainGraph_dh sg )
{
    JXF_START_FUNC_DH
    JXF_Int *outlist, *inlist;
    JXF_Int i, row, *rp = mat->rp, *cval = mat->cval;
    jxf_Numbering_dh numb;
    JXF_Int m = mat->m;
    JXF_Int *beg_rows = sg->beg_rowP, *row_count = sg->row_count, *end_rows;
    jxf_Mat_dh matFake;
    jxf_bool debug = jxf_false;
    JXF_Real *recvBuf;
    
    if (mat->debug && jxf_logFile != NULL) debug = jxf_true;
    end_rows = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    outlist = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    inlist  = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < jxf_np_dh; ++ i)
    {
        inlist[i] = 0;
        outlist[i] = 0;
        end_rows[i] = beg_rows[i] + row_count[i];
    }
    jxf_create_fake_mat_private(mat, &matFake); JXF_CHECK_V_ERROR;
    jxf_Numbering_dhCreate(&(mat->numbSolve)); JXF_CHECK_V_ERROR;
    numb = mat->numbSolve;
    jxf_Numbering_dhSetup(numb, matFake); JXF_CHECK_V_ERROR;
    jxf_destroy_fake_mat_private(matFake); JXF_CHECK_V_ERROR;
    if (debug)
    {
        jxf_fprintf(stderr, "jxf_Numbering_dhSetup completed\n");
    }
    i = m + numb->num_ext;
    mat->work_y_lo = (JXF_Real *)JXF_MALLOC_DH(i*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    mat->work_x_hi = (JXF_Real *)JXF_MALLOC_DH(i*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "FACT num_extLo= %i  num_extHi= %i\n", numb->num_extLo, numb->num_extHi);
    }
    mat->num_recvLo = 0;
    mat->num_recvHi = 0;
    if (numb->num_extLo)
    {
        recvBuf = mat->work_y_lo + m;
        mat->num_recvLo = jxf_setup_receives_private(mat, beg_rows, end_rows,
            recvBuf, mat->recv_reqLo, numb->idx_extLo, numb->num_extLo, outlist, debug); JXF_CHECK_V_ERROR;
    }
    if (numb->num_extHi)
    {
        recvBuf = mat->work_x_hi + m + numb->num_extLo;
        mat->num_recvHi = jxf_setup_receives_private(mat, beg_rows, end_rows,
                    recvBuf, mat->recv_reqHi, numb->idx_extHi, numb->num_extHi, outlist, debug); JXF_CHECK_V_ERROR;
    }
    jxf_MPI_Alltoall(outlist, 1, JXF_MPI_INT, inlist, 1, JXF_MPI_INT, jxf_comm_dh);
    jxf_setup_sends_private(mat, inlist, sg->o2n_sub, debug); JXF_CHECK_V_ERROR;
    for (row = 0; row < m; row ++)
    {
        JXF_Int len = rp[row+1] - rp[row];
        JXF_Int *ind = cval + rp[row];
        
        jxf_Numbering_dhGlobalToLocal(numb, len, ind, ind); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(outlist); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(inlist); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(end_rows); JXF_CHECK_V_ERROR;
    if (debug)
    {
        JXF_Int ii, jj;
        
        jxf_fprintf(jxf_logFile, "\n--------- row/col structure, after global to local renumbering\n");
        for (ii = 0; ii < mat->m; ++ ii)
        {
            jxf_fprintf(jxf_logFile, "local row %i :: ", ii+1);
            for (jj = mat->rp[ii]; jj < mat->rp[ii+1]; ++ jj)
            {
                jxf_fprintf(jxf_logFile, "%i ", 1+mat->cval[jj]);
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
        jxf_fprintf(jxf_logFile, "\n");
        fflush(jxf_logFile);
    }
    JXF_END_FUNC_DH
}

static void jxf_forward_solve_private( JXF_Int m, JXF_Int from, JXF_Int to,
          JXF_Int *rp, JXF_Int *cval, JXF_Int *diag, JXF_Real *aval, JXF_Real *rhs, JXF_Real *work_y, jxf_bool debug );
static void jxf_backward_solve_private( JXF_Int m, JXF_Int from, JXF_Int to,
          JXF_Int *rp, JXF_Int *cval, JXF_Int *diag, JXF_Real *aval, JXF_Real *work_y, JXF_Real *work_x, jxf_bool debug );
static JXF_Int jxf_beg_rowG;

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhSolve"
void jxf_Factor_dhSolve( JXF_Real *rhs, JXF_Real *lhs, jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    jxf_Factor_dh mat = ctx->F;
    JXF_Int from, to;
    JXF_Int ierr, i, m = mat->m, first_bdry = mat->first_bdry;
    JXF_Int offsetLo = mat->numbSolve->num_extLo;
    JXF_Int offsetHi = mat->numbSolve->num_extHi;
    JXF_Int *rp = mat->rp, *cval = mat->cval, *diag = mat->diag;
    JXF_Real *aval = mat->aval;
    JXF_Int *sendindLo = mat->sendindLo, *sendindHi = mat->sendindHi;
    JXF_Int sendlenLo = mat->sendlenLo, sendlenHi = mat->sendlenHi;
    JXF_Real *sendbufLo = mat->sendbufLo, *sendbufHi = mat->sendbufHi; 
    JXF_Real *work_y = mat->work_y_lo;
    JXF_Real *work_x = mat->work_x_hi;
    jxf_bool debug = jxf_false;
    
    if (mat->debug && jxf_logFile != NULL) debug = jxf_true;
    if (debug) jxf_beg_rowG = ctx->F->beg_row;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\n=====================================================\n");
        jxf_fprintf(jxf_logFile, "FACT jxf_Factor_dhSolve: num_recvLo= %i num_recvHi = %i\n", mat->num_recvLo, mat->num_recvHi);
    }
    if (mat->num_recvLo)
    {
        jxf_MPI_Startall(mat->num_recvLo, mat->recv_reqLo);
    }
    if (mat->num_recvHi)
    {
        jxf_MPI_Startall(mat->num_recvHi, mat->recv_reqHi);
    }
    from = 0;
    to = first_bdry;
    if (from != to)
    {
        jxf_forward_solve_private(m, from, to, rp, cval, diag, aval, rhs, work_y, debug); JXF_CHECK_V_ERROR;
    }
    if (mat->num_recvLo)
    {
        jxf_MPI_Waitall(mat->num_recvLo, mat->recv_reqLo, mat->status);
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "FACT got 'y' values from lower neighbors; work buffer:\n  ");
            for (i = 0; i < offsetLo; ++ i)
            {
                jxf_fprintf(jxf_logFile, "%g ", work_y[m+i]);
            }
        }
    }
    from = first_bdry;
    to = m;
    if (from != to)
    {
        jxf_forward_solve_private(m, from, to, rp, cval, diag, aval, rhs, work_y, debug); JXF_CHECK_V_ERROR;
    }
    if (mat->num_sendHi)
    {
        for (i = 0; i < sendlenHi; i ++)
        {
            sendbufHi[i] = work_y[sendindHi[i]];
        }
        jxf_MPI_Startall(mat->num_sendHi, mat->send_reqHi);
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "\nFACT sending 'y' values to higher neighbor:\nFACT   ");
            for (i = 0; i < sendlenHi; i ++)
            {
                jxf_fprintf(jxf_logFile, "%g ", sendbufHi[i]);
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    if (mat->num_recvHi)
    {
        ierr = jxf_MPI_Waitall(mat->num_recvHi, mat->recv_reqHi, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "FACT got 'x' values from higher neighbors:\n  ");
            for (i = m+offsetLo; i < m+offsetLo+offsetHi; ++ i)
            {
                jxf_fprintf(jxf_logFile, "%g ", work_x[i]);
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    from = m;
    to = first_bdry;
    if (from != to)
    {
        jxf_backward_solve_private(m, from, to, rp, cval, diag, aval, work_y, work_x, debug); JXF_CHECK_V_ERROR;
    }
    if (mat->num_sendLo)
    {
        for (i = 0; i < sendlenLo; i ++)
        {
            sendbufLo[i] = work_x[sendindLo[i]];
        }
        ierr = jxf_MPI_Startall(mat->num_sendLo, mat->send_reqLo); JXF_CHECK_MPI_V_ERROR(ierr);
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "\nFACT sending 'x' values to lower neighbor:\nFACT   ");
            for (i = 0; i < sendlenLo; i ++)
            {
                jxf_fprintf(jxf_logFile, "%g ", sendbufLo[i]);
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    from = first_bdry;
    to = 0;
    if (from != to)
    {
        jxf_backward_solve_private(m, from, to, rp, cval, diag, aval, work_y, work_x, debug); JXF_CHECK_V_ERROR;
    }
    memcpy(lhs, work_x, m*sizeof(JXF_Real));
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nFACT solution: ");
        for (i = 0; i < m; ++ i)
        {
            jxf_fprintf(jxf_logFile, "%g ", lhs[i]);
        }
        jxf_fprintf(jxf_logFile, "\n");
    }
    if (mat->num_sendLo)
    {
        ierr = jxf_MPI_Waitall(mat->num_sendLo, mat->send_reqLo, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
    }
    if (mat->num_sendHi)
    {
        ierr = jxf_MPI_Waitall(mat->num_sendHi, mat->send_reqHi, mat->status); JXF_CHECK_MPI_V_ERROR(ierr);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_forward_solve_private"
void jxf_forward_solve_private( JXF_Int m,
                            JXF_Int from,
                            JXF_Int to,
                            JXF_Int *rp,
                            JXF_Int *cval,
                            JXF_Int *diag,
                            JXF_Real *aval,
                            JXF_Real *rhs,
                            JXF_Real *work_y,
                            jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, idx;
    
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nFACT starting jxf_forward_solve_private; from= %i; to= %i, m= %i\n", 1+from, 1+to, m);
    }
    if (debug)
    {
        for (i = from; i < to; ++ i)
        {
            JXF_Int len = diag[i] - rp[i];
            JXF_Int *col = cval + rp[i];
            JXF_Real *val  = aval + rp[i];
            JXF_Real sum = rhs[i];
            
            jxf_fprintf(jxf_logFile, "FACT   solving for work_y[%i] (global)\n", i+1+jxf_beg_rowG);
            jxf_fprintf(jxf_logFile, "FACT        sum = %g\n", sum);
            for (j = 0; j < len; ++ j)
            {
                idx = col[j];
                sum -= (val[j] * work_y[idx]);
                jxf_fprintf(jxf_logFile, "FACT        sum(%g) -= val[j] (%g) * work_y[%i] (%g)\n",
                                                                sum, val[j], 1+idx, work_y[idx]);
            }
            work_y[i] = sum;
            jxf_fprintf(jxf_logFile, "FACT  work_y[%i] = %g\n", 1+i+jxf_beg_rowG, work_y[i]);
            jxf_fprintf(jxf_logFile, "-----------\n");
        }
        jxf_fprintf(jxf_logFile, "\nFACT   work vector at end of forward solve:\n");
        for (i = 0; i < to; i ++) jxf_fprintf(jxf_logFile, "    %i %g\n", i+1+jxf_beg_rowG, work_y[i]);
    }
    else
    {
        for (i = from; i < to; ++ i)
        {
            JXF_Int len = diag[i] - rp[i];
            JXF_Int *col = cval + rp[i];
            JXF_Real *val = aval + rp[i];
            JXF_Real sum = rhs[i];
            
            for (j = 0; j < len; ++ j)
            {
                idx = col[j];
                sum -= (val[j] * work_y[idx]);
            }
            work_y[i] = sum;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_backward_solve_private"
void jxf_backward_solve_private( JXF_Int m,
                             JXF_Int from,
                             JXF_Int to,
                             JXF_Int *rp,
                             JXF_Int *cval,
                             JXF_Int *diag,
                             JXF_Real *aval,
                             JXF_Real *work_y,
                             JXF_Real *work_x,
                             jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, idx;
    
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nFACT starting jxf_backward_solve_private; from= %i; to= %i, m= %i\n", 1+from, 1+to, m);
        for (i = from-1; i >= to; -- i)
        {
            JXF_Int len = rp[i+1] - diag[i] - 1;
            JXF_Int *col = cval + diag[i] + 1;
            JXF_Real *val = aval + diag[i] + 1;
            JXF_Real sum = work_y[i];
            
            jxf_fprintf(jxf_logFile, "FACT   solving for work_x[%i]\n", i+1+jxf_beg_rowG);
            for (j = 0; j < len; ++ j)
            {
                idx = col[j];
                sum -= (val[j] * work_x[idx]);
                jxf_fprintf(jxf_logFile, "FACT        sum(%g) -= val[j] (%g) * work_x[idx] (%g)\n",
                                                                        sum, val[j], work_x[idx]);
            }
            work_x[i] = sum * aval[diag[i]];
            jxf_fprintf(jxf_logFile, "FACT   work_x[%i] = %g\n", 1+i, work_x[i]);
            jxf_fprintf(jxf_logFile, "----------\n");
        }
    }
    else
    {
        for (i = from-1; i >= to; -- i)
        {
            JXF_Int len = rp[i+1] - diag[i] - 1;
            JXF_Int *col = cval + diag[i] + 1;
            JXF_Real *val  = aval + diag[i] + 1;
            JXF_Real sum = work_y[i];
            
            for (j = 0; j < len; ++ j)
            {
                idx = col[j];
                sum -= (val[j] * work_x[idx]);
            }
            work_x[i] = sum * aval[diag[i]];
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhInit"
void jxf_Factor_dhInit( void *A, jxf_bool fillFlag, jxf_bool avalFlag, JXF_Real rho, JXF_Int id, JXF_Int beg_rowP, jxf_Factor_dh *Fout )
{
    JXF_START_FUNC_DH
    JXF_Int m, n, beg_row, alloc;
    jxf_Factor_dh F;
    
    jxf_EuclidGetDimensions(A, &beg_row, &m, &n); JXF_CHECK_V_ERROR;
    alloc = rho * m;
    jxf_Factor_dhCreate(&F); JXF_CHECK_V_ERROR;
   *Fout = F;
    F->m = m;
    F->n = n;
    F->beg_row = beg_rowP;
    F->id = id;
    F->alloc = alloc;
    F->rp = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    F->rp[0] = 0;
    F->cval = (JXF_Int *)JXF_MALLOC_DH(alloc*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    F->diag = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (fillFlag)
    {
        F->fill = (JXF_Int *)JXF_MALLOC_DH(alloc*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    }
    if (avalFlag)
    {
        F->aval = (JXF_REAL_DH *)JXF_MALLOC_DH(alloc*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhReallocate"
void jxf_Factor_dhReallocate( jxf_Factor_dh F, JXF_Int used, JXF_Int additional )
{
    JXF_START_FUNC_DH
    JXF_Int alloc = F->alloc;
    
    if (used+additional > F->alloc)
    {
        JXF_Int *tmpI;
        
        while (alloc < used+additional) alloc *= 2.0;
        F->alloc = alloc;
        tmpI = F->cval;
        F->cval = (JXF_Int *)JXF_MALLOC_DH(alloc*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        memcpy(F->cval, tmpI, used*sizeof(JXF_Int));
        JXF_FREE_DH(tmpI); JXF_CHECK_V_ERROR;
        if (F->fill != NULL)
        {
            tmpI = F->fill;
            F->fill = (JXF_Int *)JXF_MALLOC_DH(alloc*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
            memcpy(F->fill, tmpI, used*sizeof(JXF_Int));
            JXF_FREE_DH(tmpI); JXF_CHECK_V_ERROR;
        }
        if (F->aval != NULL)
        {
            JXF_REAL_DH *tmpF = F->aval;
            F->aval = (JXF_REAL_DH *)JXF_MALLOC_DH(alloc*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
            memcpy(F->aval, tmpF, used*sizeof(JXF_REAL_DH));
            JXF_FREE_DH(tmpF); JXF_CHECK_V_ERROR;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhTranspose"
void jxf_Factor_dhTranspose( jxf_Factor_dh A, jxf_Factor_dh *Bout )
{
    JXF_START_FUNC_DH
    jxf_Factor_dh B;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only for sequential");
    }
    jxf_Factor_dhCreate(&B); JXF_CHECK_V_ERROR;
   *Bout = B;
    B->m = B->n = A->m;
    if (B->aval == NULL)
    {
        jxf_mat_dh_transpose_private(A->m, A->rp, &B->rp, A->cval, &B->cval, A->aval, NULL); JXF_CHECK_V_ERROR;
    }
    else
    {
        jxf_mat_dh_transpose_private(A->m, A->rp, &B->rp, A->cval, &B->cval, A->aval, &B->aval); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhSolveSeq"
void jxf_Factor_dhSolveSeq( JXF_Real *rhs, JXF_Real *lhs, jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    jxf_Factor_dh F = ctx->F;
    JXF_Int *rp, *cval, *diag;
    JXF_Int i, j, *vi, nz, m = F->m;
    JXF_REAL_DH *aval, *work;
    JXF_REAL_DH *v, sum;
    jxf_bool debug = jxf_false;
    
    if (ctx->F->debug && jxf_logFile != NULL) debug = jxf_true;
    rp = F->rp;
    cval = F->cval;
    aval = F->aval;
    diag = F->diag;
    work = ctx->work;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nFACT ============================================================\n");
        jxf_fprintf(jxf_logFile, "FACT starting jxf_Factor_dhSolveSeq\n");
        jxf_fprintf(jxf_logFile, "\nFACT   STARTING FORWARD SOLVE\n------------\n");
        work[0] = rhs[0];
        jxf_fprintf(jxf_logFile, "FACT   work[0] = %g\n------------\n", work[0]);
        for (i = 1; i < m; i ++)
        {
            v = aval + rp[i];
            vi = cval + rp[i];
            nz = diag[i] - rp[i];
            jxf_fprintf(jxf_logFile, "FACT   solving for work[%i]\n", i+1);
            sum = rhs[i];
            for (j = 0; j < nz; ++ j)
            {
                sum -= (v[j] * work[vi[j]]);
                jxf_fprintf(jxf_logFile, "FACT         sum (%g) -= v[j] (%g) * work[vi[j]] (%g)\n", sum, v[j], work[vi[j]]);
            }
            work[i] = sum;
            jxf_fprintf(jxf_logFile, "FACT   work[%i] = %g\n------------\n", 1+i, work[i]);
        }
        jxf_fprintf(jxf_logFile, "\nFACT   work vector at end of forward solve:\n");
        for (i = 0; i < m; i ++) jxf_fprintf(jxf_logFile, "    %i %g\n", i+1, work[i]);
        jxf_fprintf(jxf_logFile, "\nFACT   STARTING JXF_BACKWARD SOLVE\n--------------\n");
        for (i = m-1; i >= 0; i --)
        {
            v = aval + diag[i] + 1;
            vi = cval + diag[i] + 1;
            nz = rp[i+1] - diag[i] - 1;
            jxf_fprintf(jxf_logFile, "FACT   solving for lhs[%i]\n", i+1);
            sum = work[i];
            for (j = 0; j < nz; ++ j)
            {
                sum -= (v[j] * work[vi[j]]);
                jxf_fprintf(jxf_logFile, "FACT         sum (%g) -= v[j] (%g) * work[vi[j]] (%g)\n", sum, v[j], work[vi[j]]);
            }
            lhs[i] = work[i] = sum * aval[diag[i]];
            jxf_fprintf(jxf_logFile, "FACT   lhs[%i] = %g\n------------\n", 1+i, lhs[i]);
            jxf_fprintf(jxf_logFile, "FACT   solving for lhs[%i]\n", i+1);
        }
        jxf_fprintf(jxf_logFile, "\nFACT solution: ");
        for (i = 0; i < m; ++ i) jxf_fprintf(jxf_logFile, "%g ", lhs[i]);
        jxf_fprintf(jxf_logFile, "\n");
    }
    else
    {
        work[0] = rhs[0];
        for (i = 1; i < m; i ++)
        {
            v = aval + rp[i];
            vi = cval + rp[i];
            nz = diag[i] - rp[i];
            sum = rhs[i];
            while (nz --) sum -= (*v++ * work[*vi++]);
            work[i] = sum;
        }
        for (i = m-1; i >= 0; i --)
        {
            v = aval + diag[i] + 1;
            vi = cval + diag[i] + 1;
            nz = rp[i+1] - diag[i] - 1;
            sum = work[i];
            while (nz --) sum -= (*v++ * work[*vi++]);
            lhs[i] = work[i] = sum*aval[diag[i]];
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_adjust_bj_private"
void jxf_adjust_bj_private( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int nz = mat->rp[mat->m];
    JXF_Int beg_row = mat->beg_row;
    
    for (i = 0; i < nz; ++ i) mat->cval[i] += beg_row;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_unadjust_bj_private"
void jxf_unadjust_bj_private( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int nz = mat->rp[mat->m];
    JXF_Int beg_row = mat->beg_row;
    
    for (i = 0; i < nz; ++ i) mat->cval[i] -= beg_row;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhMaxPivotInverse"
JXF_Real jxf_Factor_dhMaxPivotInverse( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = mat->m, *diags = mat->diag;
    JXF_REAL_DH *aval = mat->aval;
    JXF_Real minGlobal = 0.0, min = aval[diags[0]];
    JXF_Real retval;
    
    for (i = 0; i < m; ++ i) min = JXF_MIN(min, fabs(aval[diags[i]]));
    if (jxf_np_dh == 1)
    {
        minGlobal = min;
    }
    else
    {
        jxf_MPI_Reduce(&min, &minGlobal, 1, JXF_MPI_REAL, MPI_MIN, 0, jxf_comm_dh);
    }
    if (minGlobal == 0)
    {
        retval = 0;
    }
    else
    {
        retval = 1.0 / minGlobal;
    }
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhMaxValue"
JXF_Real jxf_Factor_dhMaxValue( jxf_Factor_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Real maxGlobal = 0.0, max = 0.0;
    JXF_Int i, nz = mat->rp[mat->m];
    JXF_REAL_DH *aval = mat->aval;
    
    for (i = 0; i < nz; ++ i)
    {
        max = JXF_MAX(max, fabs(aval[i]));
    }
    if (jxf_np_dh == 1)
    {
        maxGlobal = max;
    }
    else
    {
        jxf_MPI_Reduce(&max, &maxGlobal, 1, JXF_MPI_REAL, MPI_MAX, 0, jxf_comm_dh);
    }
    JXF_END_FUNC_VAL(maxGlobal)
}

#undef __FUNC__
#define __FUNC__ "jxf_Factor_dhCondEst"
JXF_Real jxf_Factor_dhCondEst( jxf_Factor_dh mat, jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Real max = 0.0, maxGlobal = 0.0;
    JXF_Real *x;
    JXF_Int i, m = mat->m;
    jxf_Vec_dh lhs, rhs;
    
    jxf_Vec_dhCreate(&lhs); JXF_CHECK_ERROR(-1);
    jxf_Vec_dhInit(lhs, m); JXF_CHECK_ERROR(-1);
    jxf_Vec_dhDuplicate(lhs,&rhs); JXF_CHECK_ERROR(-1);
    jxf_Vec_dhSet(rhs, 1.0); JXF_CHECK_ERROR(-1);
    jxf_Euclid_dhApply(ctx, rhs->vals, lhs->vals); JXF_CHECK_ERROR(-1);
    x = lhs->vals;
    for (i = 0; i < m; ++ i)
    {
        max = JXF_MAX(max, fabs(x[i]));
    }
    if (jxf_np_dh == 1)
    {
        maxGlobal = max;
    }
    else
    {
        jxf_MPI_Reduce(&max, &maxGlobal, 1, JXF_MPI_REAL, MPI_MAX, 0, jxf_comm_dh);
    }
    JXF_END_FUNC_VAL(maxGlobal)
}
