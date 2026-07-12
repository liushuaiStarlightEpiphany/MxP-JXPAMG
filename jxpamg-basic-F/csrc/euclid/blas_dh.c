//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  blas_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_matvec_euclid_seq"
void jxf_matvec_euclid_seq( JXF_Int n, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, JXF_Real *x, JXF_Real *y )
{
    JXF_START_FUNC_DH
    JXF_Int i, j;
    JXF_Int from, to, col;
    JXF_Real sum;
    
    if (jxf_np_dh > 1) JXF_SET_V_ERROR("only for sequential case!\n");
#ifdef USING_OPENMP_DH
#pragma omp parallel private(j, col, sum, from, to) default(shared) firstprivate(n, rp, cval, aval, x, y)
#endif
    {
#ifdef USING_OPENMP_DH
#pragma omp for schedule(static)
#endif
        for (i = 0; i < n; ++ i)
        {
            sum = 0.0;
            from = rp[i];
            to = rp[i+1];
            for (j = from; j < to; ++ j)
            {
                col = cval[j];
                sum += (aval[j] * x[col]);
            }
            y[i] = sum;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Axpy"
void jxf_Axpy( JXF_Int n, JXF_Real alpha, JXF_Real *x, JXF_Real *y )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(static) firstprivate(alpha, x, y) private(i)
#endif
    for (i = 0; i < n; ++ i)
    {
        y[i] = alpha * x[i] + y[i];
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_CopyVec"
void jxf_CopyVec( JXF_Int n, JXF_Real *xIN, JXF_Real *yOUT )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(static) firstprivate(yOUT, xIN) private(i)
#endif
    for (i = 0; i < n; ++ i)
    {
        yOUT[i] = xIN[i];
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_ScaleVec"
void jxf_ScaleVec( JXF_Int n, JXF_Real alpha, JXF_Real *x )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(static) firstprivate(alpha, x) private(i)
#endif
    for (i = 0; i < n; ++ i)
    {
        x[i] *= alpha;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_InnerProd"
JXF_Real jxf_InnerProd( JXF_Int n, JXF_Real *x, JXF_Real *y )
{
    JXF_START_FUNC_DH
    JXF_Real result, local_result = 0.0;
    JXF_Int i;
    
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(static) firstprivate(x, y) private(i) reduction(+:local_result)
#endif
    for (i = 0; i < n; ++ i)
    {
        local_result += x[i] * y[i];
    }
    if (jxf_np_dh > 1)
    {
        jxf_MPI_Allreduce(&local_result, &result, 1, JXF_MPI_REAL, MPI_SUM, jxf_comm_dh);
    }
    else
    {
        result = local_result;
    }
    JXF_END_FUNC_VAL(result)
}

#undef __FUNC__
#define __FUNC__ "jxf_Norm2"
JXF_Real jxf_Norm2( JXF_Int n, JXF_Real *x )
{
    JXF_START_FUNC_DH
    JXF_Real result, local_result = 0.0;
    JXF_Int i;
    
#ifdef USING_OPENMP_DH
#pragma omp parallel for schedule(static) firstprivate(x) private(i) reduction(+:local_result)
#endif
    for (i = 0; i < n; ++ i)
    {
        local_result += (x[i] * x[i]);
    }
    if (jxf_np_dh > 1)
    {
        jxf_MPI_Allreduce(&local_result, &result, 1, JXF_MPI_REAL, MPI_SUM, jxf_comm_dh);
    }
    else
    {
        result = local_result;
    }
    result = sqrt(result);
    JXF_END_FUNC_VAL(result)
}
