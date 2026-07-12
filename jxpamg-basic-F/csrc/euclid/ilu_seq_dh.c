//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ilu_seq_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

static jxf_bool jxf_check_constraint_private( jxf_Euclid_dh ctx, JXF_Int b, JXF_Int j );
static JXF_Int jxf_symbolic_row_private( JXF_Int localRow, JXF_Int *list, JXF_Int *marker,
             JXF_Int *tmpFill, JXF_Int len, JXF_Int *CVAL, JXF_Real *AVAL, JXF_Int *o2n_col, jxf_Euclid_dh ctx, jxf_bool debug );
static JXF_Int jxf_numeric_row_private( JXF_Int localRow, JXF_Int len,
                JXF_Int *CVAL, JXF_Real *AVAL, JXF_REAL_DH *work, JXF_Int *o2n_col, jxf_Euclid_dh ctx, jxf_bool debug );

#undef __FUNC__
#define __FUNC__ "jxf_compute_scaling_private"
void jxf_compute_scaling_private( JXF_Int row, JXF_Int len, JXF_Real *AVAL, jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Real tmp = 0.0;
    JXF_Int j;
    
    for (j = 0; j < len; ++ j) tmp = JXF_MAX(tmp, fabs(AVAL[j]));
    if (tmp)
    {
        ctx->scale[row] = 1.0 / tmp;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_iluk_seq"
void jxf_iluk_seq( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int *rp, *cval, *diag;
    JXF_Int *CVAL;
    JXF_Int i, j, len, count, col, idx = 0;
    JXF_Int *list, *marker, *fill, *tmpFill;
    JXF_Int temp, m, from = ctx->from, to = ctx->to;
    JXF_Int *n2o_row, *o2n_col, beg_row, beg_rowP;
    JXF_Real *AVAL;
    JXF_REAL_DH *work, *aval;
    jxf_Factor_dh F = ctx->F;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_ilu")) debug = jxf_true;
    m = F->m;
    rp = F->rp;
    cval = F->cval;
    fill = F->fill;
    diag = F->diag;
    aval = F->aval;
    work = ctx->work;
    count = rp[from];
    if (sg == NULL)
    {
        JXF_SET_V_ERROR("subdomain graph is NULL");
    }
    n2o_row = ctx->sg->n2o_row;
    o2n_col = ctx->sg->o2n_col;
    beg_row = ctx->sg->beg_row[jxf_myid_dh];
    beg_rowP = ctx->sg->beg_rowP[jxf_myid_dh];
    list = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    tmpFill = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) marker[i] = -1;
    for (i = 0; i < m; ++ i) work[i] = 0.0;
    for (i = from; i < to; ++ i)
    {
        JXF_Int row = n2o_row[i];
        JXF_Int globalRow = row + beg_row;
        
        if (debug)
        {
            jxf_fprintf(jxf_logFile,
              "ILU_seq ================================= starting local row: %i, (global= %i) level= %i\n",
                                              i+1, i+1+sg->beg_rowP[jxf_myid_dh], ctx->level);
        }
        jxf_EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        if (ctx->isScaled)
        {
            jxf_compute_scaling_private(i, len, AVAL, ctx); JXF_CHECK_V_ERROR;
        }
        count = jxf_symbolic_row_private(i, list, marker, tmpFill, len, CVAL, AVAL, o2n_col, ctx, debug); JXF_CHECK_V_ERROR;
        if (idx + count > F->alloc)
        {
            jxf_Factor_dhReallocate(F, idx, count); JXF_CHECK_V_ERROR;
            JXF_SET_INFO("REALLOCATED from ilu_seq");
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
        jxf_numeric_row_private(i, len, CVAL, AVAL, work, o2n_col, ctx, debug); JXF_CHECK_V_ERROR
        jxf_EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "ILU_seq:  ");
            for (j = rp[i]; j < rp[i+1]; ++ j)
            {
                col = cval[j];
                aval[j] = work[col];
                work[col] = 0.0;
                jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", 1+cval[j], fill[j], aval[j]);
                fflush(jxf_logFile);
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
        else
        {
            for (j = rp[i]; j < rp[i+1]; ++ j)
            {
                col = cval[j];
                aval[j] = work[col];
                work[col] = 0.0;
            }
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
    if (beg_rowP)
    {
        JXF_Int start = rp[from];
        JXF_Int stop = rp[to];
        
        for (i = start; i < stop; ++ i) cval[i] += beg_rowP;
    }
    for (i = to+1; i < m; ++ i) rp[i] = 0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_iluk_seq_block"
void jxf_iluk_seq_block( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int *rp, *cval, *diag;
    JXF_Int *CVAL;
    JXF_Int h, i, j, len, count, col, idx = 0;
    JXF_Int *list, *marker, *fill, *tmpFill;
    JXF_Int temp, m;
    JXF_Int *n2o_row, *o2n_col, *beg_rowP, *n2o_sub, blocks;
    JXF_Int *row_count, *dummy = NULL, dummy2[1];
    JXF_Real *AVAL;
    JXF_REAL_DH *work, *aval;
    jxf_Factor_dh F = ctx->F;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    jxf_bool bj = jxf_false, constrained = jxf_false;
    JXF_Int discard = 0;
    JXF_Int gr = -1;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_ilu")) debug = jxf_true;
    if (!strcmp(ctx->algo_par, "bj")) bj = jxf_true;
    constrained = ! jxf_Parser_dhHasSwitch(jxf_parser_dh, "-unconstrained");
    m = F->m;
    rp = F->rp;
    cval = F->cval;
    fill = F->fill;
    diag = F->diag;
    aval = F->aval;
    work = ctx->work;
    if (sg != NULL)
    {
        n2o_row = sg->n2o_row;
        o2n_col = sg->o2n_col;
        row_count = sg->row_count;
        beg_rowP = sg->beg_rowP;
        n2o_sub = sg->n2o_sub;
        blocks = sg->blocks;
    }
    else
    {
        dummy = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        for (i = 0; i < m; ++ i) dummy[i] = i;
        n2o_row = dummy;
        o2n_col = dummy;
        dummy2[0] = m; row_count = dummy2;
        beg_rowP = dummy;
        n2o_sub = dummy;
        blocks = 1;
    }
    list = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    tmpFill = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) marker[i] = -1;
    for (i = 0; i < m; ++ i) work[i] = 0.0;
    for (h = 0; h < blocks; ++ h)
    {
        JXF_Int curBlock = n2o_sub[h];
        JXF_Int first_row = beg_rowP[curBlock];
        JXF_Int end_row = first_row + row_count[curBlock];
        
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "\n\nILU_seq BLOCK: %i @@@@@@@@@@@@@@@ \n", curBlock);
        }
        for (i = first_row; i < end_row; ++ i)
        {
            JXF_Int row = n2o_row[i];
            
            ++ gr;
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_seq  global: %i  local: %i =================================\n",
                                            1+gr, 1+i-first_row);
            }
            jxf_EuclidGetRow(ctx->A, row, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
            if (ctx->isScaled)
            {
                jxf_compute_scaling_private(i, len, AVAL, ctx); JXF_CHECK_V_ERROR;
            }
            count = jxf_symbolic_row_private(i, list, marker,
                          tmpFill, len, CVAL, AVAL, o2n_col, ctx, debug); JXF_CHECK_V_ERROR;
            if (idx + count > F->alloc)
            {
                jxf_Factor_dhReallocate(F, idx, count); JXF_CHECK_V_ERROR;
                JXF_SET_INFO("REALLOCATED from ilu_seq");
                cval = F->cval;
                fill = F->fill;
                aval = F->aval;
            }
            col = list[m];
            while (count --)
            {
                if (constrained && !bj)
                {
                    if (col >= first_row && col < end_row)
                    {
                        cval[idx] = col;
                        fill[idx] = tmpFill[col];
                        ++ idx;
                    }
                    else
                    {
                        if (jxf_check_constraint_private(ctx, curBlock, col))
                        {
                            cval[idx] = col;
                            fill[idx] = tmpFill[col];
                            ++ idx;
                        }
                        else
                        {
                            ++ discard;
                        }
                    }
                    col = list[col];
                }
                else if (bj)
                {
                    if (col >= first_row && col < end_row)
                    {
                        cval[idx] = col;
                        fill[idx] = tmpFill[col];
                        ++ idx;
                    }
                    else
                    {
                        ++ discard;
                    }
                    col = list[col];
                }
                else
                {
                    cval[idx] = col;
                    fill[idx] = tmpFill[col];
                    ++ idx;
                    col = list[col];
                }
            }
            rp[i+1] = idx;
            temp = rp[i];
            while (cval[temp] != i) ++ temp;
            diag[i] = temp;
            jxf_numeric_row_private(i, len, CVAL, AVAL, work, o2n_col, ctx, debug); JXF_CHECK_V_ERROR
            jxf_EuclidRestoreRow(ctx->A, row, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_seq: ");
                for (j = rp[i]; j < rp[i+1]; ++ j)
                {
                    col = cval[j];
                    aval[j] = work[col];
                    work[col] = 0.0;
                    jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", 1+cval[j], fill[j], aval[j]);
                }
                jxf_fprintf(jxf_logFile, "\n");
            }
            else
            {
                for (j = rp[i]; j < rp[i+1]; ++ j)
                {
                    col = cval[j];
                    aval[j] = work[col];
                    work[col] = 0.0;
                }
            }
            if (!aval[diag[i]])
            {
                jxf_sprintf(jxf_msgBuf_dh, "zero diagonal in local row %i", i+1);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
        }
    }
    if (dummy != NULL)
    {
        JXF_FREE_DH(dummy); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(list); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(tmpFill); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_symbolic_row_private"
JXF_Int jxf_symbolic_row_private( JXF_Int localRow,
                          JXF_Int *list,
                          JXF_Int *marker,
                          JXF_Int *tmpFill,
                          JXF_Int len,
                          JXF_Int *CVAL,
                          JXF_Real *AVAL,
                          JXF_Int *o2n_col,
                          jxf_Euclid_dh ctx,
                          jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int level = ctx->level, m = ctx->F->m;
    JXF_Int *cval = ctx->F->cval, *diag = ctx->F->diag, *rp = ctx->F->rp;
    JXF_Int *fill = ctx->F->fill;
    JXF_Int count = 0;
    JXF_Int j, node, tmp, col, head;
    JXF_Int fill1, fill2, beg_row;
    JXF_Real val;
    JXF_Real thresh = ctx->sparseTolA;
    JXF_REAL_DH scale;
    
    scale = ctx->scale[localRow];
    ctx->stats[JXF_NZA_STATS] += (JXF_Real)len;
    beg_row = ctx->sg->beg_row[jxf_myid_dh];
    list[m] = m;
    for (j = 0; j < len; ++ j)
    {
        tmp = m;
        col = *CVAL ++;
        col -= beg_row;
        col = o2n_col[col];
        val = *AVAL ++;
        val *= scale;
        if (fabs(val) > thresh || col == localRow)
        {
            ++ count;
            while (col > list[tmp]) tmp = list[tmp];
            list[col] = list[tmp];
            list[tmp] = col;
            tmpFill[col] = 0;
            marker[col] = localRow;
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
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_seq   sf updating from row: %i\n", 1+node);
            }
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
                         JXF_Int len,
                         JXF_Int *CVAL,
                         JXF_Real *AVAL,
                         JXF_REAL_DH *work,
                         JXF_Int *o2n_col,
                         jxf_Euclid_dh ctx,
                         jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Real  pc, pv, multiplier;
    JXF_Int j, k, col, row;
    JXF_Int *rp = ctx->F->rp, *cval = ctx->F->cval;
    JXF_Int *diag = ctx->F->diag;
    JXF_Int beg_row;
    JXF_Real val;
    JXF_REAL_DH *aval = ctx->F->aval, scale;
    
    scale = ctx->scale[localRow];
    beg_row  = ctx->sg->beg_row[jxf_myid_dh];
    for (j = rp[localRow]; j < rp[localRow+1]; ++ j)
    {
        col = cval[j];
        work[col] = 0.0;
    }
    for (j = 0; j < len; ++ j)
    {
        col = *CVAL ++;
        col -= beg_row;
        val = *AVAL ++;
        col = o2n_col[col];
        work[col] = val * scale;
    }
    for (j = rp[localRow]; j < diag[localRow]; ++ j)
    {
        row = cval[j];
        pc = work[row];
        pv = aval[diag[row]];
        if (pc != 0.0 && pv != 0.0)
        {
            multiplier = pc / pv;
            work[row] = multiplier;
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_seq   nf updating from row: %i; multiplier= %g\n", 1+row, multiplier);
            }
            for (k = diag[row]+1; k < rp[row+1]; ++ k)
            {
                col = cval[k];
                work[col] -= (multiplier * aval[k]);
            }
        }
        else
        {
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_seq   nf NO UPDATE from row %i; pc = %g; pv = %g\n", 1+row, pc, pv);
            }
        }
    }
    JXF_END_FUNC_VAL(0)
}

JXF_Int jxf_ilut_row_private( JXF_Int localRow, JXF_Int *list, JXF_Int *o2n_col, JXF_Int *marker,
                JXF_Int len, JXF_Int *CVAL, JXF_Real *AVAL, JXF_REAL_DH *work, jxf_Euclid_dh ctx, jxf_bool debug );

#undef __FUNC__
#define __FUNC__ "jxf_ilut_seq"
void jxf_ilut_seq( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int *rp, *cval, *diag, *CVAL;
    JXF_Int i, len, count, col, idx = 0;
    JXF_Int *list, *marker;
    JXF_Int temp, m, from, to;
    JXF_Int *n2o_row, *o2n_col, beg_row, beg_rowP;
    JXF_Real *AVAL, droptol; 
    JXF_REAL_DH *work, *aval, val;
    jxf_Factor_dh F = ctx->F;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    jxf_bool debug = jxf_false;
    
    if (jxf_logFile != NULL && jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_ilu")) debug = jxf_true;
    m = F->m;
    rp = F->rp;
    cval = F->cval;
    diag = F->diag;
    aval = F->aval;
    work = ctx->work;
    from = ctx->from;
    to = ctx->to;
    count = rp[from];
    droptol = ctx->droptol;
    if (sg == NULL)
    {
        JXF_SET_V_ERROR("subdomain graph is NULL");
    }
    n2o_row = ctx->sg->n2o_row;
    o2n_col = ctx->sg->o2n_col;
    beg_row = ctx->sg->beg_row[jxf_myid_dh];
    beg_rowP = ctx->sg->beg_rowP[jxf_myid_dh];
    list = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) marker[i] = -1;
    rp[0] = 0;
    for (i = 0; i < m; ++ i) work[i] = 0.0;
    for (i = from; i < to; ++ i)
    {
        JXF_Int row = n2o_row[i];
        JXF_Int globalRow = row + beg_row;
        
        jxf_EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        jxf_compute_scaling_private(i, len, AVAL, ctx); JXF_CHECK_V_ERROR; 
        count = jxf_ilut_row_private(i, list, o2n_col, marker, len, CVAL, AVAL, work, ctx, debug); JXF_CHECK_V_ERROR;
        jxf_EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        if (idx + count > F->alloc)
        {
            jxf_Factor_dhReallocate(F, idx, count); JXF_CHECK_V_ERROR;
            JXF_SET_INFO("REALLOCATED from ilu_seq");
            cval = F->cval;
            aval = F->aval;
        }
        col = list[m];
        while (count --)
        {
            val = work[col];
            if (col == i || fabs(val) > droptol)
            {
                cval[idx] = col;
                aval[idx++] = val;
                work[col] = 0.0;
            }
            col = list[col];
        }
        rp[i+1] = idx;
        temp = rp[i];
        while (cval[temp] != i) ++ temp;
        diag[i] = temp;
        if (!aval[diag[i]])
        {
            jxf_sprintf(jxf_msgBuf_dh, "zero diagonal in local row %i", i+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    if (beg_rowP)
    {
        JXF_Int start = rp[from];
        JXF_Int stop = rp[to];
        
        for (i = start; i < stop; ++ i) cval[i] += beg_rowP;
    }
    JXF_FREE_DH(list);
    JXF_FREE_DH(marker);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_ilut_row_private"
JXF_Int jxf_ilut_row_private( JXF_Int localRow,
                      JXF_Int *list,
                      JXF_Int *o2n_col,
                      JXF_Int *marker,
                      JXF_Int len,
                      JXF_Int *CVAL,
                      JXF_Real *AVAL,
                      JXF_REAL_DH *work,
                      jxf_Euclid_dh ctx,
                      jxf_bool debug )
{
    JXF_START_FUNC_DH
    jxf_Factor_dh F = ctx->F;
    JXF_Int j, col, m = ctx->m, *rp = F->rp, *cval = F->cval;
    JXF_Int tmp, *diag = F->diag;
    JXF_Int head;
    JXF_Int count = 0, beg_row;
    JXF_Real val;
    JXF_Real mult, *aval = F->aval;
    JXF_Real scale, pv, pc;
    JXF_Real droptol = ctx->droptol;
    JXF_Real thresh = ctx->sparseTolA;
    
    scale = ctx->scale[localRow];
    ctx->stats[JXF_NZA_STATS] += (JXF_Real)len;
    beg_row = ctx->sg->beg_row[jxf_myid_dh];
    list[m] = m;
    for (j = 0; j < len; ++ j)
    {
        tmp = m;
        col = *CVAL ++;
        col -= beg_row;
        col = o2n_col[col];
        val = *AVAL ++;
        val *= scale;
        if (fabs(val) > thresh || col == localRow)
        {
            ++ count;
            while (col > list[tmp]) tmp = list[tmp];
            list[col] = list[tmp];
            list[tmp] = col;
            work[col] = val;
            marker[col] = localRow;
        }
    }
    if (marker[localRow] != localRow)
    {
        tmp = m;
        while (localRow > list[tmp]) tmp = list[tmp];
        list[localRow] = list[tmp];
        list[tmp] = localRow;
        marker[localRow] = localRow;
        ++ count;
    }
    head = m;
    while (list[head] < localRow)
    {
        JXF_Int row = list[head];
        
        pc = work[row];
        if (pc != 0.0)
        {
            pv = aval[diag[row]];
            mult = pc / pv;
            if (fabs(mult) > droptol)
            {
                work[row] = mult;
                for (j = diag[row]+1; j < rp[row+1]; ++ j)
                {
                    col = cval[j];
                    work[col] -= (mult * aval[j]);
                    if (marker[col] < localRow)
                    {
                        marker[col] = localRow;
                        tmp = head;
                        while (col > list[tmp]) tmp = list[tmp];
                        list[col] = list[tmp];
                        list[tmp] = col;
                        ++ count;
                    }
                }
            }
        }
        head = list[head];
    }
    JXF_END_FUNC_VAL(count)
}

#undef __FUNC__
#define __FUNC__ "jxf_check_constraint_private"
jxf_bool jxf_check_constraint_private( jxf_Euclid_dh ctx, JXF_Int p1, JXF_Int j )
{
    JXF_START_FUNC_DH
    jxf_bool retval = jxf_false;
    JXF_Int i, p2;
    JXF_Int *nabors, count;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    
    if (sg == NULL)
    {
        JXF_SET_ERROR(-1, "ctx->sg == NULL");
    }
    p2 = jxf_SubdomainGraph_dhFindOwner(ctx->sg, j, jxf_true);
    nabors = sg->adj + sg->ptrs[p1];
    count = sg->ptrs[p1+1] - sg->ptrs[p1];
    for (i = 0; i < count; ++ i)
    {
        if (nabors[i] == p2)
        {
            retval = jxf_true;
            break;
        }
    }
    JXF_END_FUNC_VAL(retval)
}
