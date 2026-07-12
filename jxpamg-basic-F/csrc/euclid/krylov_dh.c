//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  krylov_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_bicgstab_euclid"
void jxf_bicgstab_euclid( jxf_Mat_dh A, jxf_Euclid_dh ctx, JXF_Real *x, JXF_Real *b, JXF_Int *itsOUT )
{
    JXF_START_FUNC_DH
    JXF_Int its, m = ctx->m;
    jxf_bool monitor;
    JXF_Int maxIts = ctx->maxIts;
    JXF_Real atol = ctx->atol, rtol = ctx->rtol;
    JXF_Real alpha, alpha_1 = 0., beta_1, widget, widget_1 = 0.;
    JXF_Real rho_1, rho_2 = 0., s_norm, eps, exit_a, b_iprod, r_iprod;
    JXF_Real *t, *s, *s_hat, *v, *p, *p_hat, *r, *r_hat;
    
    monitor = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-monitor");
    t = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    s = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    s_hat = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    v = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    p = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    p_hat = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    r = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    r_hat = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    jxf_Mat_dhMatVec(A, x, s);
    jxf_CopyVec(m, b, r);
    jxf_Axpy(m, -1.0, s, r);
    jxf_CopyVec(m, r, r_hat);
    b_iprod = jxf_InnerProd(m, b, b); JXF_CHECK_V_ERROR;
    exit_a = atol*atol*b_iprod;  JXF_CHECK_V_ERROR;
    eps = rtol * rtol * b_iprod;
    its = 0;
    while (1)
    {
        ++ its;
        rho_1 = jxf_InnerProd(m, r_hat, r);
        if (rho_1 == 0)
        {
            JXF_SET_V_ERROR("(r_hat . r) = 0; method fails");
        }
        if (its == 1)
        {
            jxf_CopyVec(m, r, p); JXF_CHECK_V_ERROR;
        }
        else
        {
            beta_1 = (rho_1 / rho_2) * (alpha_1 / widget_1);
            jxf_Axpy(m, -widget_1, v, p); JXF_CHECK_V_ERROR;
            jxf_ScaleVec(m, beta_1, p); JXF_CHECK_V_ERROR;
            jxf_Axpy(m, 1.0, r, p); JXF_CHECK_V_ERROR;
        }
        jxf_Euclid_dhApply(ctx, p, p_hat); JXF_CHECK_V_ERROR;
        jxf_Mat_dhMatVec(A, p_hat, v); JXF_CHECK_V_ERROR;
        JXF_Real tmp = jxf_InnerProd(m, r_hat, v); JXF_CHECK_V_ERROR;
        alpha = rho_1 / tmp;
        jxf_CopyVec(m, r, s); JXF_CHECK_V_ERROR;
        jxf_Axpy(m, -alpha, v, s); JXF_CHECK_V_ERROR;
        s_norm = jxf_InnerProd(m, s, s);
        if (s_norm < exit_a)
        {
            JXF_SET_INFO("reached absolute stopping criteria");
            break;
        }
        jxf_Euclid_dhApply(ctx, s, s_hat); JXF_CHECK_V_ERROR;
        jxf_Mat_dhMatVec(A, s_hat, t); JXF_CHECK_V_ERROR;
        JXF_Real tmp1, tmp2;
        tmp1 = jxf_InnerProd(m, t, s); JXF_CHECK_V_ERROR;
        tmp2 = jxf_InnerProd(m, t, t); JXF_CHECK_V_ERROR;
        widget = tmp1 / tmp2;
        jxf_Axpy(m, alpha, p_hat, x); JXF_CHECK_V_ERROR;
        jxf_Axpy(m, widget, s_hat, x); JXF_CHECK_V_ERROR;
        jxf_CopyVec(m, s, r); JXF_CHECK_V_ERROR;
        jxf_Axpy(m, -widget, t, r); JXF_CHECK_V_ERROR;
        r_iprod = jxf_InnerProd(m, r, r); JXF_CHECK_V_ERROR;
        if (r_iprod < eps)
        {
            JXF_SET_INFO("stipulated residual reduction achieved");
            break;
        }
        if (monitor && jxf_myid_dh == 0)
        {
            jxf_fprintf(stderr, "[it = %i] %e\n", its, sqrt(r_iprod/b_iprod));
        }
        rho_2 = rho_1;
        widget_1 = widget;
        alpha_1 = alpha;
        if (its >= maxIts)
        {
            its = -its;
            break;
        }
    }
   *itsOUT = its;
    JXF_FREE_DH(t);
    JXF_FREE_DH(s);
    JXF_FREE_DH(s_hat);
    JXF_FREE_DH(v);
    JXF_FREE_DH(p);
    JXF_FREE_DH(p_hat);
    JXF_FREE_DH(r);
    JXF_FREE_DH(r_hat);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_cg_euclid"
void jxf_cg_euclid( jxf_Mat_dh A, jxf_Euclid_dh ctx, JXF_Real *x, JXF_Real *b, JXF_Int *itsOUT )
{
    JXF_START_FUNC_DH
    JXF_Int its, m = A->m;
    JXF_Real *p, *r, *s;
    JXF_Real alpha, beta, gamma, gamma_old, eps, bi_prod, i_prod;
    jxf_bool monitor;
    JXF_Int maxIts = ctx->maxIts;
    JXF_Real rtol = ctx->rtol;
    
    monitor = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-monitor");
    bi_prod = jxf_InnerProd(m, b, b); JXF_CHECK_V_ERROR;
    eps = (rtol * rtol) * bi_prod;
    p = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    s = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    r = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real));
    jxf_Mat_dhMatVec(A, x, r); JXF_CHECK_V_ERROR;
    jxf_ScaleVec(m, -1.0, r); JXF_CHECK_V_ERROR;
    jxf_Axpy(m, 1.0, b, r); JXF_CHECK_V_ERROR;
    jxf_Euclid_dhApply(ctx, r, p); JXF_CHECK_V_ERROR;
    gamma = jxf_InnerProd(m, r, p); JXF_CHECK_V_ERROR;
    its = 0;
    while (1)
    {
        ++ its;
        jxf_Mat_dhMatVec(A, p, s);  JXF_CHECK_V_ERROR;
        JXF_Real tmp = jxf_InnerProd(m, s, p); JXF_CHECK_V_ERROR;
        alpha = gamma / tmp;
        gamma_old = gamma;
        jxf_Axpy(m, alpha, p, x); JXF_CHECK_V_ERROR;
        jxf_Axpy(m, -alpha, s, r); JXF_CHECK_V_ERROR;
        jxf_Euclid_dhApply(ctx, r, s); JXF_CHECK_V_ERROR;
        gamma = jxf_InnerProd(m, r, s); JXF_CHECK_V_ERROR;
        i_prod = jxf_InnerProd(m, r, r); JXF_CHECK_V_ERROR;
        if (monitor && jxf_myid_dh == 0)
        {
            jxf_fprintf(stderr, "iter = %i  rel. resid. norm: %e\n", its, sqrt(i_prod/bi_prod));
        }
        if (i_prod < eps) break;
        beta = gamma / gamma_old;
        jxf_ScaleVec(m, beta, p); JXF_CHECK_V_ERROR;
        jxf_Axpy(m, 1.0, s, p); JXF_CHECK_V_ERROR;
        if (its >= maxIts)
        {
            its = -its;
            break;
        }
    }
   *itsOUT = its;
    JXF_FREE_DH(p);
    JXF_FREE_DH(s);
    JXF_FREE_DH(r);
    JXF_END_FUNC_DH
}
