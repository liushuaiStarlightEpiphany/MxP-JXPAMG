//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  external_rows_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

enum{ JXF_ROW_CT_TAG, JXF_NZ_CT_TAG, JXF_ROW_LENGTH_TAG, JXF_ROW_NUMBER_TAG, JXF_CVAL_TAG, JXF_FILL_TAG, JXF_AVAL_TAG };

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhCreate"
void jxf_ExternalRows_dhCreate( jxf_ExternalRows_dh *er )
{
    JXF_START_FUNC_DH
    struct _jxf_extrows_dh *tmp = (struct _jxf_extrows_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_extrows_dh)); JXF_CHECK_V_ERROR;
    
   *er = tmp;
    if (JXF_MAX_MPI_TASKS < jxf_np_dh)
    {
        JXF_SET_V_ERROR("JXF_MAX_MPI_TASKS is too small; change, then recompile!");
    }
    JXF_Int i;
    for (i = 0; i < JXF_MAX_MPI_TASKS; ++ i)
    {
        tmp->rcv_row_lengths[i] = NULL;
        tmp->rcv_row_numbers[i] = NULL;
    }
    tmp->cvalExt = NULL;
    tmp->fillExt = NULL;
    tmp->avalExt = NULL;
    tmp->my_row_counts = NULL;
    tmp->my_row_numbers = NULL;
    tmp->cvalSend = NULL;
    tmp->fillSend = NULL;
    tmp->avalSend = NULL;
    tmp->rowLookup = NULL;
    tmp->sg = NULL;
    tmp->F = NULL;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_ExtRows");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhDestroy"
void jxf_ExternalRows_dhDestroy( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
    for (i = 0; i < JXF_MAX_MPI_TASKS; ++ i)
    {
        if (er->rcv_row_lengths[i] != NULL)
        {
            JXF_FREE_DH(er->rcv_row_lengths[i]); JXF_CHECK_V_ERROR;
        }
        if (er->rcv_row_numbers[i] != NULL)
        {
            JXF_FREE_DH(er->rcv_row_numbers[i]); JXF_CHECK_V_ERROR;
        }
    }
    if (er->cvalExt != NULL)
    {
        JXF_FREE_DH(er->cvalExt); JXF_CHECK_V_ERROR;
    }
    if (er->fillExt != NULL)
    {
        JXF_FREE_DH(er->fillExt); JXF_CHECK_V_ERROR;
    }
    if (er->avalExt != NULL)
    {
        JXF_FREE_DH(er->avalExt); JXF_CHECK_V_ERROR;
    }
    if (er->my_row_counts != NULL)
    {
        JXF_FREE_DH(er->my_row_counts); JXF_CHECK_V_ERROR;
    }
    if (er->my_row_numbers != NULL)
    {
        JXF_FREE_DH(er->my_row_numbers); JXF_CHECK_V_ERROR;
    }
    if (er->cvalSend != NULL)
    {
        JXF_FREE_DH(er->cvalSend); JXF_CHECK_V_ERROR;
    }
    if (er->fillSend != NULL)
    {
        JXF_FREE_DH(er->fillSend); JXF_CHECK_V_ERROR;
    }
    if (er->avalSend != NULL)
    {
        JXF_FREE_DH(er->avalSend); JXF_CHECK_V_ERROR;
    }
    if (er->rowLookup != NULL)
    {
        jxf_Hash_dhDestroy(er->rowLookup); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(er); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhInit"
void jxf_ExternalRows_dhInit( jxf_ExternalRows_dh er, jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    er->sg = ctx->sg;
    er->F = ctx->F;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhGetRow"
void jxf_ExternalRows_dhGetRow( jxf_ExternalRows_dh er, JXF_Int globalRow, JXF_Int *len, JXF_Int **cval, JXF_Int **fill, JXF_REAL_DH **aval )
{
    JXF_START_FUNC_DH
    if (er->rowLookup == NULL)
    {
       *len = 0;
    }
    else
    {
        jxf_HashData *r = NULL;
        r = jxf_Hash_dhLookup(er->rowLookup, globalRow); JXF_CHECK_V_ERROR;
        if (r != NULL)
        {
           *len = r->iData;
            if (cval != NULL) *cval = r->iDataPtr;
            if (fill != NULL) *fill = r->iDataPtr2;
            if (aval != NULL) *aval = r->fDataPtr;
        }
        else
        {
           *len = 0;
        }
    }
    JXF_END_FUNC_DH
}

static void jxf_rcv_ext_storage_private( jxf_ExternalRows_dh extRows );
static void jxf_build_hash_table_private( jxf_ExternalRows_dh er );
static void jxf_rcv_external_rows_private( jxf_ExternalRows_dh er );
static void jxf_allocate_ext_row_storage_private( jxf_ExternalRows_dh er );
static void jxf_print_received_rows_private( jxf_ExternalRows_dh er );

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhRecvRows"
void jxf_ExternalRows_dhRecvRows( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && er->debug) debug = jxf_true;
    if (er->sg->loCount > 0)
    {
        jxf_rcv_ext_storage_private(er); JXF_CHECK_V_ERROR;
        jxf_allocate_ext_row_storage_private(er); JXF_CHECK_V_ERROR;
        jxf_build_hash_table_private(er); JXF_CHECK_V_ERROR;
        jxf_rcv_external_rows_private(er); JXF_CHECK_V_ERROR;
        if (debug)
        {
            jxf_print_received_rows_private(er); JXF_CHECK_V_ERROR;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_rcv_ext_storage_private"
void jxf_rcv_ext_storage_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
    JXF_Int *rcv_row_counts = er->rcv_row_counts;
    JXF_Int *rcv_nz_counts = er->rcv_nz_counts;
    JXF_Int **lengths = er->rcv_row_lengths, **numbers = er->rcv_row_numbers;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && er->debug) debug = jxf_true;
    for (i = 0; i < loCount; ++ i)
    {
        JXF_Int nabor = loNabors[i];
        
        jxf_MPI_Irecv(rcv_row_counts+i, 1, JXF_MPI_INT, nabor, JXF_ROW_CT_TAG, jxf_comm_dh, er->req1+i);
        jxf_MPI_Irecv(rcv_nz_counts+i, 1, JXF_MPI_INT, nabor, JXF_NZ_CT_TAG, jxf_comm_dh, er->req2+i);
    }
    jxf_MPI_Waitall(loCount, er->req1, er->status);
    jxf_MPI_Waitall(loCount, er->req2, er->status);
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nEXR jxf_rcv_ext_storage_private:: <nabor,rowCount,nzCount>\nEXR ");
        for (i = 0; i < loCount; ++ i)
        {
            jxf_fprintf(jxf_logFile, "<%i,%i,%i> ", loNabors[i], rcv_row_counts[i], rcv_nz_counts[i]);
        }
    }
    for (i = 0; i < loCount; ++ i)
    {
        JXF_Int nz = rcv_nz_counts[i];
        JXF_Int nabor = loNabors[i];
        
        lengths[i] = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        numbers[i] = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        jxf_MPI_Irecv(lengths[i], nz, JXF_MPI_INT, nabor, JXF_ROW_LENGTH_TAG, jxf_comm_dh, er->req1+i);
        jxf_MPI_Irecv(numbers[i], nz, JXF_MPI_INT, nabor, JXF_ROW_NUMBER_TAG, jxf_comm_dh, er->req2+i);
    }
    jxf_MPI_Waitall(loCount, er->req1, er->status);
    jxf_MPI_Waitall(loCount, er->req2, er->status);
    if (debug)
    {
        JXF_Int j, nz;
        
        for (i = 0; i < loCount; ++ i)
        {
            jxf_fprintf(jxf_logFile, "\nEXR rows <number,length> to be received from P_%i\nEXR ", loNabors[i]);
            nz = rcv_row_counts[i];
            for (j = 0; j < nz; ++ j) jxf_fprintf(jxf_logFile, "<%i,%i> ", numbers[i][j], lengths[i][j]);
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_allocate_ext_row_storage_private"
void jxf_allocate_ext_row_storage_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int i, nz = 0;
    JXF_Int loCount = er->sg->loCount;
    JXF_Int *rcv_nz_counts = er->rcv_nz_counts;
    
    for (i = 0; i < loCount; ++ i) nz += rcv_nz_counts[i];
    er->cvalExt = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    er->fillExt = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    er->avalExt = (JXF_REAL_DH *)JXF_MALLOC_DH(nz*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_build_hash_table_private"
void jxf_build_hash_table_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int loCount = er->sg->loCount;
    JXF_Int i, j, offset, rowCt = 0;
    jxf_Hash_dh table;
    jxf_HashData record;
    JXF_Int *extRowCval = er->cvalExt, *extRowFill = er->fillExt;
    JXF_REAL_DH *extRowAval = er->avalExt;
    JXF_Int *rcv_row_counts = er->rcv_row_counts;
    JXF_Int **rcv_row_numbers = er->rcv_row_numbers;
    JXF_Int **rcv_row_lengths = er->rcv_row_lengths;
    
    for (i = 0; i < loCount; ++ i) rowCt += rcv_row_counts[i];
    jxf_Hash_dhCreate(&table, rowCt); JXF_CHECK_V_ERROR;
    er->rowLookup = table;
    offset = 0;
    for (i = 0; i < loCount; ++ i)
    {
        JXF_Int rowCount = rcv_row_counts[i];
        
        for (j = 0; j < rowCount; ++ j)
        {
            JXF_Int row = rcv_row_numbers[i][j];
            JXF_Int rowLength = rcv_row_lengths[i][j];
            
            record.iData = rowLength;
            record.iDataPtr = extRowCval + offset;
            record.iDataPtr2 = extRowFill + offset;
            record.fDataPtr = extRowAval + offset;
            jxf_Hash_dhInsert(table, row, &record); JXF_CHECK_V_ERROR;
            offset += rowLength;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_rcv_external_rows_private"
void jxf_rcv_external_rows_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int *rcv_nz_counts = er->rcv_nz_counts;
    JXF_Int i, loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
    JXF_Int nabor, nz = 0, offset = 0;
    JXF_Int *extRowCval = er->cvalExt, *extRowFill = er->fillExt;
    JXF_Real *extRowAval = er->avalExt;
    
    nz = 0;
    for (i = 0; i < loCount; ++ i)
    {
        nabor = loNabors[i];
        nz = rcv_nz_counts[i];
        jxf_MPI_Irecv(extRowCval+offset, nz, JXF_MPI_INT, nabor, JXF_CVAL_TAG, jxf_comm_dh, er->req1+i);
        jxf_MPI_Irecv(extRowFill+offset, nz, JXF_MPI_INT, nabor, JXF_FILL_TAG, jxf_comm_dh, er->req2+i);
        jxf_MPI_Irecv(extRowAval+offset, nz, JXF_MPI_REAL, nabor, JXF_AVAL_TAG, jxf_comm_dh, er->req3+i);
        offset += nz;
    }
    jxf_MPI_Waitall(loCount, er->req1, er->status);
    jxf_MPI_Waitall(loCount, er->req2, er->status);
    jxf_MPI_Waitall(loCount, er->req3, er->status);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_print_received_rows_private"
void jxf_print_received_rows_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    jxf_bool noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
    JXF_Int i, j, k, rwCt, idx = 0, nabor;
    JXF_Int loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
    JXF_Int n = er->F->n;
    
    jxf_fprintf(jxf_logFile, "\nEXR ================= received rows, printed from buffers ==============\n");
    for (i = 0; i < loCount; ++ i)
    {
        rwCt = er->rcv_row_counts[i];
        nabor = loNabors[i];
        jxf_fprintf(jxf_logFile, "\nEXR Rows received from P_%i:\n", nabor);
        for (j = 0; j < rwCt; ++ j)
        {
            JXF_Int rowNum = er->rcv_row_numbers[i][j];
            JXF_Int rowLen  = er->rcv_row_lengths[i][j];
            
            jxf_fprintf(jxf_logFile, "EXR %i :: ", 1+rowNum);
            for (k = 0; k < rowLen; ++ k)
            {
                if (noValues)
                {
                    jxf_fprintf(jxf_logFile, "%i,%i ; ", er->cvalExt[idx], er->fillExt[idx]);
                }
                else
                {
                    jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", er->cvalExt[idx], er->fillExt[idx], er->avalExt[idx]);
                }
                ++ idx;
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    jxf_fprintf(jxf_logFile, "\nEXR =============== received rows, printed from hash table =============\n");
    for (i = 0; i < n; ++ i)
    {
        JXF_Int len, *cval, *fill;
        JXF_REAL_DH *aval;
        
        jxf_ExternalRows_dhGetRow(er, i, &len, &cval, &fill, &aval); JXF_CHECK_V_ERROR;
        if (len > 0)
        {
            jxf_fprintf(jxf_logFile, "EXR %i :: ", i+1);
            for (j = 0; j < len; ++ j)
            {
                if (noValues)
                {
                    jxf_fprintf(jxf_logFile, "%i,%i ; ", cval[j], fill[j]);
                }
                else
                {
                    jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", cval[j], fill[j], aval[j]);
                }
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    JXF_END_FUNC_DH
}

static void jxf_send_ext_storage_private( jxf_ExternalRows_dh er );
static void jxf_send_external_rows_private( jxf_ExternalRows_dh er );
static void jxf_waitfor_sends_private( jxf_ExternalRows_dh er );

#undef __FUNC__
#define __FUNC__ "jxf_ExternalRows_dhSendRows"
void jxf_ExternalRows_dhSendRows( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    if (er->sg->hiCount > 0)
    {
        jxf_send_ext_storage_private(er); JXF_CHECK_V_ERROR;
        jxf_send_external_rows_private(er); JXF_CHECK_V_ERROR;
        jxf_waitfor_sends_private(er); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_send_ext_storage_private"
void jxf_send_ext_storage_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int nz, i, j;
    JXF_Int *nzCounts, *nzNumbers;
    JXF_Int hiCount = er->sg->hiCount, *hiNabors = er->sg->hiNabors;
    JXF_Int *rp = er->F->rp, *diag = er->F->diag;
    JXF_Int m = er->F->m;
    JXF_Int beg_row = er->F->beg_row;
    JXF_Int rowCount = er->F->bdry_count;
    JXF_Int first_bdry = er->F->first_bdry;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && er->debug) debug = jxf_true;
    nzCounts = er->my_row_counts = (JXF_Int *)JXF_MALLOC_DH(rowCount*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    nzNumbers = er->my_row_numbers = (JXF_Int *)JXF_MALLOC_DH(rowCount*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    nz = 0;
    for (i = first_bdry, j = 0; i < m; ++ i, ++ j)
    {
        JXF_Int tmp = (rp[i+1] - diag[i]);
        
        nz += tmp;
        nzCounts[j] = tmp;
    }
    er->nzSend = nz;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "EXR jxf_send_ext_storage_private:: rowCount = %i\n", rowCount);
        jxf_fprintf(jxf_logFile, "EXR jxf_send_ext_storage_private:: nz Count = %i\n", nz);
    }
    for (i = 0; i < hiCount; ++ i)
    {
        JXF_Int nabor = hiNabors[i];
        
        jxf_MPI_Isend(&rowCount, 1, JXF_MPI_INT, nabor, JXF_ROW_CT_TAG, jxf_comm_dh, er->req1+i);
        jxf_MPI_Isend(&nz, 1, JXF_MPI_INT, nabor, JXF_NZ_CT_TAG, jxf_comm_dh, er->req2+i);
    }
    for (i = 0, j = first_bdry; j < m; ++ i, ++ j)
    {
        nzNumbers[i] = j + beg_row;
    }
    for (i = 0; i < hiCount; ++ i)
    {
        JXF_Int nabor = hiNabors[i];
        
        jxf_MPI_Isend(nzNumbers, rowCount, JXF_MPI_INT, nabor, JXF_ROW_NUMBER_TAG, jxf_comm_dh, er->req3+i);
        jxf_MPI_Isend(nzCounts, rowCount, JXF_MPI_INT, nabor, JXF_ROW_LENGTH_TAG, jxf_comm_dh, er->req4+i);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_send_external_rows_private"
void jxf_send_external_rows_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, hiCount = er->sg->hiCount, *hiNabors = er->sg->hiNabors;
    JXF_Int offset, nz = er->nzSend;
    JXF_Int *cvalSend, *fillSend;
    JXF_REAL_DH *avalSend;
    JXF_Int *cval = er->F->cval, *fill = er->F->fill;
    JXF_Int m = er->F->m;
    JXF_Int *rp = er->F->rp, *diag = er->F->diag;
    JXF_Int first_bdry = er->F->first_bdry;
    JXF_REAL_DH *aval = er->F->aval;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && er->debug) debug = jxf_true;
    cvalSend = er->cvalSend = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    fillSend = er->fillSend = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    avalSend = er->avalSend = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    offset = 0;
    for (i = first_bdry, j = 0; i < m; ++ i, ++ j)
    {
        JXF_Int tmp = (rp[i+1] - diag[i]);
        
        memcpy(cvalSend+offset, cval+diag[i], tmp*sizeof(JXF_Int));
        memcpy(fillSend+offset, fill+diag[i], tmp*sizeof(JXF_Int));
        memcpy(avalSend+offset, aval+diag[i], tmp*sizeof(JXF_Real));
        offset += tmp;
    }
    if (debug)
    {
        JXF_Int beg_row = er->F->beg_row;
        JXF_Int idx = 0;
        jxf_bool noValues = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues"));
        
        jxf_fprintf(jxf_logFile, "\nEXR ======================= send buffers ======================\n");
        for (i = first_bdry, j = 0; i < m; ++ i, ++ j)
        {
            JXF_Int tmp = (rp[i+1] - diag[i]);
            
            jxf_fprintf(jxf_logFile, "EXR %i :: ", i+beg_row);
            for (j = 0; j < tmp; ++ j)
            {
                if (noValues)
                {
                    jxf_fprintf(jxf_logFile, "%i,%i ; ", cvalSend[idx], fillSend[idx]);
                }
                else
                {
                    jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", cvalSend[idx], fillSend[idx], avalSend[idx]);
                }
                ++ idx;
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
    }
    for (i = 0; i < hiCount; ++ i)
    {
        JXF_Int nabor = hiNabors[i];
        
        jxf_MPI_Isend(cvalSend, nz, JXF_MPI_INT, nabor, JXF_CVAL_TAG, jxf_comm_dh, er->cval_req+i);
        jxf_MPI_Isend(fillSend, nz, JXF_MPI_INT, nabor, JXF_FILL_TAG, jxf_comm_dh, er->fill_req+i);
        jxf_MPI_Isend(avalSend, nz, JXF_MPI_REAL, nabor, JXF_AVAL_TAG, jxf_comm_dh, er->aval_req+i);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_waitfor_sends_private"
void jxf_waitfor_sends_private( jxf_ExternalRows_dh er )
{
    JXF_START_FUNC_DH
    MPI_Status *status = er->status;
    JXF_Int hiCount = er->sg->hiCount;
    
    if (hiCount)
    {
        jxf_MPI_Waitall(hiCount, er->req1, status);
        jxf_MPI_Waitall(hiCount, er->req2, status);
        jxf_MPI_Waitall(hiCount, er->req3, status);
        jxf_MPI_Waitall(hiCount, er->req4, status);
        jxf_MPI_Waitall(hiCount, er->cval_req, status);
        jxf_MPI_Waitall(hiCount, er->fill_req, status);
        jxf_MPI_Waitall(hiCount, er->aval_req, status);
    }
    JXF_END_FUNC_DH
}
