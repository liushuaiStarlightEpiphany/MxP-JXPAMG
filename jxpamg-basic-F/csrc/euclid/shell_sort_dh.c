//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  shell_sort_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_shellSort_int"
void jxf_shellSort_int( const JXF_Int n, JXF_Int *x )
{
    JXF_START_FUNC_DH
    JXF_Int m, max, j, k, itemp;
    
    m = n/2;
    while (m > 0)
    {
        max = n - m;
        for (j = 0; j < max; j ++)
        {
            for (k = j; k >= 0; k -= m)
            {
                if (x[k+m] >= x[k]) break;
                itemp = x[k+m];
                x[k+m] = x[k];
                x[k] = itemp;
            }
        }
        m = m / 2;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_shellSort_float"
void jxf_shellSort_float( const JXF_Int n, JXF_Real *x )
{
    JXF_START_FUNC_DH
    JXF_Int m, max, j, k;
    JXF_Real itemp;
    
    m = n / 2;
    while (m > 0)
    {
        max = n - m;
        for (j = 0; j < max; j ++)
        {
            for (k = j; k >= 0; k -= m)
            {
                if (x[k+m] >= x[k]) break;
                itemp = x[k+m];
                x[k+m] = x[k];
                x[k] = itemp;
            }
        }
        m = m / 2;
    }
    JXF_END_FUNC_DH
}
