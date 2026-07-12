//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ilu_mpi_bj_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

JXF_Int jxf_symbolic_row_private( JXF_Int localRow, JXF_Int beg_row, JXF_Int end_row,
      JXF_Int *list, JXF_Int *marker, JXF_Int *tmpFill, JXF_Int len, JXF_Int *CVAL, JXF_Real *AVAL, JXF_Int *o2n_col, jxf_Euclid_dh ctx );
static JXF_Int jxf_numeric_row_private( JXF_Int localRow, JXF_Int beg_row, JXF_Int end_row,
                JXF_Int len, JXF_Int *CVAL, JXF_Real *AVAL, JXF_REAL_DH *work, JXF_Int *o2n_col, jxf_Euclid_dh ctx );

#undef __FUNC__
#define __FUNC__ "jxf_iluk_mpi_bj"
void jxf_iluk_mpi_bj( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int *rp, *cval, *diag;
    JXF_Int *CVAL;
    JXF_Int i, j, len, count, col, idx = 0;
    JXF_Int *list, *marker, *fill, *tmpFill;
    JXF_Int temp, m, from = ctx->from, to = ctx->to;
    JXF_Int *n2o_row, *o2n_col;
    JXF_Int first_row, last_row;
    JXF_Real *AVAL;
    JXF_REAL_DH *work, *aval;
    jxf_Factor_dh F = ctx->F;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    
    if (ctx->F == NULL)
    {
        JXF_SET_V_ERROR("ctx->F is NULL");
    }
    if (ctx->F->rp == NULL)
    {
        JXF_SET_V_ERROR("ctx->F->rp is NULL");
    }
    m = F->m;
    rp = F->rp;
    cval = F->cval;
    fill = F->fill;
    diag = F->diag;
    aval = F->aval;
    work = ctx->work;
    n2o_row = sg->n2o_row;
    o2n_col = sg->o2n_col;
    list = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    tmpFill = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        marker[i] = -1;
        work[i] = 0.0;
    }
    first_row = sg->beg_row[jxf_myid_dh];
    last_row  = first_row + sg->row_count[jxf_myid_dh];
    for (i = from; i < to; ++ i)
    {
        JXF_Int row = n2o_row[i];
        JXF_Int globalRow = row + first_row;
        
        jxf_EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        if (ctx->isScaled)
        {
            jxf_compute_scaling_private(i, len, AVAL, ctx); JXF_CHECK_V_ERROR;
        }
        count = jxf_symbolic_row_private(i, first_row, last_row,
             list, marker, tmpFill, len, CVAL, AVAL, o2n_col, ctx); JXF_CHECK_V_ERROR;
        if (idx + count > F->alloc)
        {
            jxf_Factor_dhReallocate(F, idx, count); JXF_CHECK_V_ERROR;
            JXF_SET_INFO("REALLOCATED from lu_mpi_bj");
            cval = F->cval;
            fill = F->fill;
            aval = F->aval;
        }
        col = list[m];
        while (count --)
        {
            cval[idx] = col;
            fill[idx] = tmpFill[col];
            ++ idx;
            col = list[col];
        }
        rp[i+1] = idx;
        temp = rp[i];
        while (cval[temp] != i) ++ temp;
        diag[i] = temp;
        jxf_numeric_row_private(i, first_row, last_row, len, CVAL, AVAL, work, o2n_col, ctx); JXF_CHECK_V_ERROR
        jxf_EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            col = cval[j];
            aval[j] = work[col];
            work[col] = 0.0;
        }
        if (!aval[diag[i]])
        {
            jxf_sprintf(jxf_msgBuf_dh, "zero diagonal in local row %i", i+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    JXF_FREE_DH(list); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(tmpFill); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_symbolic_row_private"
JXF_Int jxf_symbolic_row_private( JXF_Int localRow,
                          JXF_Int beg_row,
                          JXF_Int end_row,
                          JXF_Int *list,
                          JXF_Int *marker,
                          JXF_Int *tmpFill,
                          JXF_Int len,
                          JXF_Int *CVAL,
                          JXF_Real *AVAL,
                          JXF_Int *o2n_col,
                          jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int level = ctx->level, m = ctx->F->m;
    JXF_Int *cval = ctx->F->cval, *diag = ctx->F->diag, *rp = ctx->F->rp;
    JXF_Int *fill = ctx->F->fill;
    JXF_Int count = 0;
    JXF_Int j, node, tmp, col, head;
    JXF_Int fill1, fill2;
    float val;
    JXF_Real thresh = ctx->sparseTolA;
    JXF_REAL_DH scale;
    
    scale = ctx->scale[localRow];
    ctx->stats[JXF_NZA_STATS] += (JXF_Real)len;
    list[m] = m;
    for (j = 0; j < len; ++ j)
    {
        tmp = m;
        col = *CVAL ++;
        val = *AVAL ++;
        if (col >= beg_row && col < end_row)
        {
            col -= beg_row;
            col = o2n_col[col];
            if (fabs(scale*val) > thresh || col == localRow)
            {
                ++ count;
                while (col > list[tmp]) tmp = list[tmp];
                list[col] = list[tmp];
                list[tmp] = col;
                tmpFill[col] = 0;
                marker[col] = localRow;
            }
        }
    }
    if (marker[localRow] != localRow)
    {
        tmp = m;
        while (localRow > list[tmp]) tmp = list[tmp];
        list[localRow] = list[tmp];
        list[tmp] = localRow;
        tmpFill[localRow] = 0;
        marker[localRow] = localRow;
        ++ count;
    }
    ctx->stats[JXF_NZA_USED_STATS] += (JXF_Real)count;
    head = m;
    if (level > 0)
    {
        while (list[head] < localRow)
        {
            node = list[head];
            fill1 = tmpFill[node];
            if (fill1 < level)
            {
                for (j = diag[node]+1; j < rp[node+1]; ++ j)
                {
                    col = cval[j];
                    fill2 = fill1 + fill[j] + 1;
                    if (fill2 <= level)
                    {
                        if (marker[col] < localRow)
                        {
                            tmp = head;
                            marker[col] = localRow;
                            tmpFill[col] = fill2;
                            while (col > list[tmp]) tmp = list[tmp];
                            list[col] = list[tmp];
                            list[tmp] = col;
                            ++ count;
                        }
                        else
                        {
                            tmpFill[col] = (fill2 < tmpFill[col]) ? fill2 : tmpFill[col];
                        }
                    }
                }
            }
            head = list[head];
        }
    }
    JXF_END_FUNC_VAL(count)
}

#undef __FUNC__
#define __FUNC__ "jxf_numeric_row_private"
JXF_Int jxf_numeric_row_private( JXF_Int localRow,
                         JXF_Int beg_row,
                         JXF_Int end_row,
                         JXF_Int len,
                         JXF_Int *CVAL,
                         JXF_Real *AVAL,
                         JXF_REAL_DH *work,
                         JXF_Int *o2n_col,
                         jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Real pc, pv, multiplier;
    JXF_Int j, k, col, row;
    JXF_Int *rp = ctx->F->rp, *cval = ctx->F->cval;
    JXF_Int *diag = ctx->F->diag;
    JXF_Real val;
    JXF_REAL_DH *aval = ctx->F->aval, scale;
    
    scale = ctx->scale[localRow];
    for (j = rp[localRow]; j < rp[localRow+1]; ++ j)
    {
        col = cval[j];
        work[col] = 0.0;
    }
    for (j = 0; j < len; ++ j)
    {
        col = *CVAL ++;
        val = *AVAL ++;
        if (col >= beg_row && col < end_row)
        {
            col -= beg_row;
            col = o2n_col[col];
            work[col] = val * scale;
        }
    }
    for (j = rp[localRow]; j < diag[localRow]; ++ j)
    {
        row = cval[j];
        pc = work[row];
        if (pc != 0.0)
        {
            pv = aval[diag[row]];
            multiplier = pc / pv;
            work[row] = multiplier;
            for (k = diag[row]+1; k < rp[row+1]; ++ k)
            {
                col = cval[k];
                work[col] -= (multiplier * aval[k]);
            }
        }
    }
    JXF_END_FUNC_VAL(0)
}
