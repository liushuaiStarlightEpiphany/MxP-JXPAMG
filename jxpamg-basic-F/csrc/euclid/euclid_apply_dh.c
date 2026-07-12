//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  euclid_apply_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

static void jxf_scale_rhs_private( jxf_Euclid_dh ctx, JXF_Real *rhs );
static void jxf_permute_vec_n2o_private( jxf_Euclid_dh ctx, JXF_Real *xIN, JXF_Real *xOUT );
static void jxf_permute_vec_o2n_private( jxf_Euclid_dh ctx, JXF_Real *xIN, JXF_Real *xOUT );

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhApply"
void jxf_Euclid_dhApply( jxf_Euclid_dh ctx, JXF_Real *rhs, JXF_Real *lhs )
{
    JXF_START_FUNC_DH
    JXF_Real *rhs_, *lhs_;
    JXF_Real t1, t2;
    
    t1 = jxf_MPI_Wtime();
    ctx->from = 0;
    ctx->to = ctx->m;
    if (!strcmp(ctx->algo_ilu, "none") || !strcmp(ctx->algo_par, "none"))
    {
        JXF_Int i, m = ctx->m;
        for (i = 0; i < m; ++ i) lhs[i] = rhs[i];
        goto END_OF_FUNCTION;
    }
    if (ctx->sg != NULL)
    {
        jxf_permute_vec_n2o_private(ctx, rhs, lhs); JXF_CHECK_V_ERROR;
        rhs_ = lhs;
        lhs_ = ctx->work2;
    }
    else
    {
        rhs_ = rhs;
        lhs_ = lhs;
    }
    if (ctx->isScaled)
    {
        jxf_scale_rhs_private(ctx, rhs_); JXF_CHECK_V_ERROR;
    }
    if (jxf_np_dh == 1 || !strcmp(ctx->algo_par, "bj"))
    {
        jxf_Factor_dhSolveSeq(rhs_, lhs_, ctx); JXF_CHECK_V_ERROR;
    }
    else
    {
        jxf_Factor_dhSolve(rhs_, lhs_, ctx); JXF_CHECK_V_ERROR;
    }
    if (ctx->sg != NULL)
    {
        jxf_permute_vec_o2n_private(ctx, lhs_, lhs); JXF_CHECK_V_ERROR;
    }
    
END_OF_FUNCTION: ;
    
    t2 = jxf_MPI_Wtime();
    ctx->timing[JXF_TRI_SOLVE_T] += (t2 - t1);
    ctx->timing[JXF_TOTAL_SOLVE_TEMP_T] = t2 - ctx->timing[JXF_SOLVE_START_T];
    ctx->its += 1;
    ctx->itsTotal += 1;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_scale_rhs_private"
void jxf_scale_rhs_private( jxf_Euclid_dh ctx, JXF_Real *rhs )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = ctx->m;
    JXF_REAL_DH *scale = ctx->scale;
    
    if (scale != NULL)
    {
#ifdef USING_OPENMP_DH
#pragma omp for schedule(static)
#endif
        for (i = 0; i < m; ++ i)
        {
            rhs[i] *= scale[i];
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_permute_vec_o2n_private"
void jxf_permute_vec_o2n_private( jxf_Euclid_dh ctx, JXF_Real *xIN, JXF_Real *xOUT )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = ctx->m;
    JXF_Int *o2n = ctx->sg->o2n_col;
    
    for (i = 0; i < m; ++ i) xOUT[i] = xIN[o2n[i]];
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_permute_vec_n2o_private"
void jxf_permute_vec_n2o_private( jxf_Euclid_dh ctx, JXF_Real *xIN, JXF_Real *xOUT )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = ctx->m;
    JXF_Int *n2o = ctx->sg->n2o_row;
    
    for (i = 0; i < m; ++ i) xOUT[i] = xIN[n2o[i]];
    JXF_END_FUNC_DH
}
