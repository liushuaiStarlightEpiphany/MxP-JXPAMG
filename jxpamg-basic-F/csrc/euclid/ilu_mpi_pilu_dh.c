//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_ilu_mpi_pilu_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

static void jxf_iluk_symbolic_row_private( JXF_Int localRow, JXF_Int len, JXF_Int *CVAL,
             JXF_Real *AVAL, jxf_ExternalRows_dh extRows, jxf_SortedList_dh sList, jxf_Euclid_dh ctx, jxf_bool debug );
static void jxf_iluk_numeric_row_private( JXF_Int new_row,
                  jxf_ExternalRows_dh extRows, jxf_SortedList_dh slist, jxf_Euclid_dh ctx, jxf_bool debug );

#undef __FUNC__
#define __FUNC__ "jxf_iluk_mpi_pilu"
void jxf_iluk_mpi_pilu( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int from = ctx->from, to = ctx->to;
    JXF_Int i, m;
    JXF_Int *n2o_row;
    JXF_Int *rp, *cval, *diag, *fill;
    JXF_Int beg_row, beg_rowP, end_rowP;
    jxf_SubdomainGraph_dh sg = ctx->sg;
    JXF_Int *CVAL, len, idx = 0, count;
    JXF_Real *AVAL;
    JXF_REAL_DH *aval;
    jxf_Factor_dh F = ctx->F;
    jxf_SortedList_dh slist = ctx->slist;
    jxf_ExternalRows_dh extRows = ctx->extRows;
    jxf_bool bj, noValues, debug = jxf_false;
    
    if (jxf_logFile != NULL && jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_ilu")) debug = jxf_true;
    noValues = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noValues");
    bj = ctx->F->blockJacobi;
    m = F->m;
    rp = F->rp;
    cval = F->cval;
    fill = F->fill;
    diag = F->diag;
    aval = F->aval;
    n2o_row = sg->n2o_row;
    if (from != 0) idx = rp[from];
    beg_row = sg->beg_row[jxf_myid_dh];
    beg_rowP = sg->beg_rowP[jxf_myid_dh];
    end_rowP = beg_rowP + sg->row_count[jxf_myid_dh];
    for (i = from; i < to; ++ i)
    {
        JXF_Int row = n2o_row[i];
        JXF_Int globalRow = row + beg_row;
        
        if (debug)
        {
            jxf_fprintf(jxf_logFile,
              "\nILU_pilu global: %i  old_Local: %i =========================================================\n",
                              i+1+beg_rowP, row+1);
        }
        jxf_EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        if (debug)
        {
            JXF_Int h;
            
            jxf_fprintf(jxf_logFile, "ILU_pilu  jxf_EuclidGetRow:\n");
            for (h = 0; h < len; ++ h) jxf_fprintf(jxf_logFile, "    %i   %g\n", 1+CVAL[h], AVAL[h]);
        }
        if (ctx->isScaled)
        {
            jxf_compute_scaling_private(i, len, AVAL, ctx); JXF_CHECK_V_ERROR;
        }
        jxf_SortedList_dhReset(slist, i); JXF_CHECK_V_ERROR;
        jxf_iluk_symbolic_row_private(i, len, CVAL, AVAL, extRows, slist, ctx, debug); JXF_CHECK_V_ERROR;
        jxf_SortedList_dhEnforceConstraint(slist, sg);
        if (!noValues)
        {
            jxf_iluk_numeric_row_private(i, extRows, slist, ctx, debug); JXF_CHECK_V_ERROR;
        }
        jxf_EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); JXF_CHECK_V_ERROR;
        count = jxf_SortedList_dhReadCount(slist); JXF_CHECK_V_ERROR;
        if (idx + count > F->alloc)
        {
            jxf_Factor_dhReallocate(F, idx, count); JXF_CHECK_V_ERROR;
            JXF_SET_INFO("REALLOCATED from jxf_ilu_mpi_pilu");
            cval = F->cval;
            fill = F->fill;
            aval = F->aval;
        }
        if (bj)
        {
            JXF_Int col;
            
            while (count --)
            {
                jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(slist); JXF_CHECK_V_ERROR;
                col = sr->col;
                if (col >= beg_rowP && col < end_rowP)
                {
                    cval[idx] = col;
                    if (noValues)
                    {
                        aval[idx] = 0.0;
                    }
                    else
                    {
                        aval[idx] = sr->val;
                    }
                    fill[idx] = sr->level;
                    ++ idx;
                }
            }
        }
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "ILU_pilu  ");
            while (count --)
            {
                jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(slist); JXF_CHECK_V_ERROR;
                cval[idx] = sr->col;
                aval[idx] = sr->val;
                fill[idx] = sr->level;
                jxf_fprintf(jxf_logFile, "%i,%i,%g ; ", 1+cval[idx], fill[idx], aval[idx]);
                ++ idx;
            }
            jxf_fprintf(jxf_logFile, "\n");
        }
        else
        {
            while (count --)
            {
                jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(slist); JXF_CHECK_V_ERROR;
                cval[idx] = sr->col;
                aval[idx] = sr->val;
                fill[idx] = sr->level;
                ++ idx;
            }
        }
        rp[i+1] = idx;
        JXF_Int temp = rp[i];
        jxf_bool flag = jxf_true;
        while (temp < idx)
        {
            if (cval[temp] == i+beg_rowP)
            {
                diag[i] = temp;
                flag = jxf_false;
                break;
            }
            ++ temp;
        }
        if (flag)
        {
            if (jxf_logFile != NULL)
            {
                JXF_Int k;
                jxf_fprintf(jxf_logFile, "Failed to jxf_find diag in localRow %i (globalRow %i; ct= %i)\n   ",
                                                                   1+i, i+1+beg_rowP, rp[i+1] - rp[i]);
                for (k = rp[i]; k < rp[i+1]; ++ k)
                {
                    jxf_fprintf(jxf_logFile, "%i ", cval[i]+1);
                }
                jxf_fprintf(jxf_logFile, "\n\n");
            }
            jxf_sprintf(jxf_msgBuf_dh, "failed to jxf_find diagonal for localRow: %i", 1+i);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
        if (!aval[diag[i]])
        {
            jxf_sprintf(jxf_msgBuf_dh, "zero diagonal in local row %i", i+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    if (bj)
    {
        JXF_Int nz = rp[m];
        for (i = 0; i < nz; ++ i) cval[i] -= beg_rowP;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_iluk_symbolic_row_private"
void jxf_iluk_symbolic_row_private( JXF_Int localRow,
                                JXF_Int len,
                                JXF_Int *CVAL,
                                JXF_Real *AVAL,
                                jxf_ExternalRows_dh extRows,
                                jxf_SortedList_dh slist,
                                jxf_Euclid_dh ctx,
                                jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int level = ctx->level, m = ctx->m;
    JXF_Int beg_row = ctx->sg->beg_row[jxf_myid_dh];
    JXF_Int beg_rowP = ctx->sg->beg_rowP[jxf_myid_dh];
    JXF_Int *cval = ctx->F->cval, *diag = ctx->F->diag; 
    JXF_Int *rp = ctx->F->rp, *fill = ctx->F->fill;
    JXF_Int j, node, col;
    JXF_Int end_rowP = beg_rowP + m;
    JXF_Int level_1, level_2;
    JXF_Int *cvalPtr, *fillPtr;
    jxf_SRecord sr, *srPtr;
    JXF_REAL_DH scale, *avalPtr;
    JXF_Real thresh = ctx->sparseTolA;
    jxf_bool wasInserted;
    JXF_Int count = 0;
    
    scale = ctx->scale[localRow];
    ctx->stats[JXF_NZA_STATS] += (JXF_Real)len;
    sr.level = 0;
    for (j = 0; j < len; ++ j)
    {
        sr.col = CVAL[j];
        sr.val = scale * AVAL[j];
        wasInserted = jxf_SortedList_dhPermuteAndInsert(slist, &sr, thresh); JXF_CHECK_V_ERROR;
        if (wasInserted) ++ count;
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "ILU_pilu   inserted from A: col= %i  val= %g\n", 1+CVAL[j], sr.val);
        }
    }
    sr.val = 0.0;
    sr.col = localRow + beg_rowP;
    srPtr = jxf_SortedList_dhFind(slist, &sr); JXF_CHECK_V_ERROR;
    if (srPtr == NULL)
    {
        jxf_SortedList_dhInsert(slist, &sr); JXF_CHECK_V_ERROR;
        ++ count;
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "ILU_pilu   inserted missing diagonal: %i\n", 1+localRow+beg_row);
        }
    }
    ctx->stats[JXF_NZA_USED_STATS] += (JXF_Real)count;
    sr.val = 0.0;
    if (level > 0)
    {
        while(1)
        {
            srPtr = jxf_SortedList_dhGetSmallestLowerTri(slist); JXF_CHECK_V_ERROR;
            if (srPtr == NULL) break;
            node = srPtr->col;
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "ILU_pilu   sf updating from row: %i\n", 1+srPtr->col);
            }
            level_1 = srPtr->level;
            if (level_1 < level)
            {
                if (node >= beg_rowP && node < end_rowP)
                {
                    node -= beg_rowP;
                    len = rp[node+1] - diag[node] - 1;
                    cvalPtr = cval + diag[node] + 1;
                    fillPtr = fill + diag[node] + 1;
                }
                else
                {
                    len = 0;
                    jxf_ExternalRows_dhGetRow(extRows, node, &len, &cvalPtr, &fillPtr, &avalPtr); JXF_CHECK_V_ERROR;
                    if (debug && len == 0)
                    {
                        jxf_fprintf(stderr, "ILU_pilu  sf failed to get extern row: %i\n", 1+node);
                    }
                }
                for (j = 0; j < len; ++ j)
                {
                    col = *cvalPtr ++;
                    level_2 = 1+ level_1 + *fillPtr ++;
                    if (level_2 <= level)
                    {
                        sr.col = col;
                        sr.level = level_2;
                        sr.val = 0.0;
                        jxf_SortedList_dhInsertOrUpdate(slist, &sr); JXF_CHECK_V_ERROR;
                    }
                }
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_iluk_numeric_row_private"
void jxf_iluk_numeric_row_private( JXF_Int new_row, jxf_ExternalRows_dh extRows, jxf_SortedList_dh slist, jxf_Euclid_dh ctx, jxf_bool debug )
{
    JXF_START_FUNC_DH
    JXF_Int m = ctx->m;
    JXF_Int beg_rowP = ctx->sg->beg_rowP[jxf_myid_dh];
    JXF_Int end_rowP = beg_rowP + m;
    JXF_Int len, row;
    JXF_Int *rp = ctx->F->rp, *cval = ctx->F->cval, *diag = ctx->F->diag;
    JXF_REAL_DH *avalPtr, *aval = ctx->F->aval;
    JXF_Int *cvalPtr;
    JXF_Real multiplier, pc, pv;
    jxf_SRecord sr, *srPtr;
    
    jxf_SortedList_dhResetGetSmallest(slist); JXF_CHECK_V_ERROR;
    while (1)
    {
        srPtr = jxf_SortedList_dhGetSmallestLowerTri(slist); JXF_CHECK_V_ERROR;
        if (srPtr == NULL) break;
        row = srPtr->col;
        if (row >= beg_rowP && row < end_rowP)
        {
            JXF_Int local_row = row - beg_rowP;
            
            len = rp[local_row+1] - diag[local_row];
            cvalPtr = cval + diag[local_row];
            avalPtr = aval + diag[local_row];
        }
        else
        {
            len = 0;
            jxf_ExternalRows_dhGetRow(extRows, row, &len, &cvalPtr, NULL, &avalPtr); JXF_CHECK_V_ERROR;
            if (debug && len == 0)
            {
                jxf_fprintf(stderr, "ILU_pilu  failed to get extern row: %i\n", 1+row);
            }
        }
        if (len)
        {
            sr.col = row;
            srPtr = jxf_SortedList_dhFind(slist, &sr); JXF_CHECK_V_ERROR;
            if (srPtr == NULL)
            {
                jxf_sprintf(jxf_msgBuf_dh, "jxf_find failed for sr.col = %i while factoring local row= %i \n",
                                                                                    1+sr.col, new_row+1);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
            pc = srPtr->val;
            if (pc != 0.0)
            {
                pv = *avalPtr ++;
                -- len;
                ++ cvalPtr;
                multiplier = pc / pv;
                srPtr->val = multiplier;
                if (debug)
                {
                    jxf_fprintf(jxf_logFile, "ILU_pilu   nf updating from row: %i; multiplier = %g\n",
                                                                           1+srPtr->col, multiplier);
                }
                while (len --)
                {
                    sr.col = *cvalPtr ++;
                    sr.val = *avalPtr ++;
                    srPtr = jxf_SortedList_dhFind(slist, &sr); JXF_CHECK_V_ERROR;
                    if (srPtr != NULL)
                    {
                        srPtr->val -= (multiplier * sr.val);
                    }
                }
            }
            else
            {
                if (debug)
                {
                    jxf_fprintf(jxf_logFile, "ILU_pilu   NO UPDATE from row: %i; srPtr->val = 0.0\n", 1+srPtr->col);
                }
            }
        }
    }
    JXF_END_FUNC_DH
}
