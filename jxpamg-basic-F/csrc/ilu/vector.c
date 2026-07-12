//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  vector.c -- operations for sequential vectors.
 *  Date: 2014/07/05
 */

#include "jxf_mv.h"

/*!
 * \fn jxf_Vector *jxf_VectorTDMGReorderByNodes
 * \brief Reorder a vector ordered by variables by nodes
 * \author Yue Xiaoqiang
 * \date 2014/04/09
 */
jxf_Vector *
jxf_VectorTDMGReorderByNodes( jxf_Vector *x, JXF_Int num_equns )
{
    JXF_Int size = jxf_VectorSize(x);
    JXF_Real *x_data = jxf_VectorData(x);
    jxf_Vector *y = jxf_SeqVectorCreate(size);
    JXF_Real *y_data = NULL;
    JXF_Int sub_size = size / num_equns;
    JXF_Int i, k, t_i = 0;
    
    jxf_SeqVectorInitialize(y);
    y_data = jxf_VectorData(y);
    for (i = 0; i < sub_size; i ++)
    {
        for (k = 0; k < num_equns-2; k ++)
        {
            y_data[t_i++] = x_data[k*sub_size+i];
        }
        y_data[t_i++] = x_data[k*sub_size+i];
        k ++;
        y_data[t_i++] = x_data[k*sub_size+i];
    }
    
    return y;
}

/*!
 * \fn jxf_Vector *jxf_VectorTDMGReorderByVariables
 * \brief Reorder a vector ordered by nodes by variables
 * \author Yue Xiaoqiang
 * \date 2014/04/09
 */
jxf_Vector *
jxf_VectorTDMGReorderByVariables( jxf_Vector *x, JXF_Int num_equns )
{
    JXF_Int size = jxf_VectorSize(x);
    JXF_Real *x_data = jxf_VectorData(x);
    jxf_Vector *y = jxf_SeqVectorCreate(size);
    JXF_Real *y_data = NULL;
    JXF_Int sub_size = size / num_equns;
    JXF_Int Rowx = 0, Rowy, mdo;
    
    jxf_SeqVectorInitialize(y);
    y_data = jxf_VectorData(y);
    for (Rowy = 0; Rowy < size; Rowy ++)
    {
        mdo = Rowy / sub_size;
        Rowx = num_equns * (Rowy - mdo * sub_size) + mdo;
        y_data[Rowy] = x_data[Rowx];
    }
    
    return y;
}

/*!
 * \fn jxf_Vector *jxf_VectorTDMGReorderAlongX
 * \brief Reorder a vector ordered by x-axis direction
 * \author Yue Xiaoqiang
 * \date 04/09/2014
 */
jxf_Vector *
jxf_VectorTDMGReorderAlongX( jxf_Vector *x, JXF_Int num_equns, JXF_Int nx, JXF_Int ny )
{
    JXF_Int size = jxf_VectorSize(x);
    JXF_Real *x_data = jxf_VectorData(x);
    jxf_Vector *y = jxf_SeqVectorCreate(size);
    JXF_Real *y_data = NULL;
    JXF_Int i, j, k, fxp, t_i = 0;
    
    jxf_SeqVectorInitialize(y);
    y_data = jxf_VectorData(y);
    for (i = 0; i < ny; i ++)
    {
        for (j = 0; j < nx; j ++)
        {
            fxp = num_equns * (j * ny + i);
            for (k = 0; k < num_equns; k ++)
            {
                y_data[t_i++] = x_data[fxp+k];
            }
        }
    }
    
    return y;
}

/*!
 * \fn jxf_Vector *jxf_VectorTDMGReorderAlongY
 * \brief Reorder a vector ordered by y-axis direction
 * \author Yue Xiaoqiang
 * \date 04/09/2014
 */
jxf_Vector *
jxf_VectorTDMGReorderAlongY( jxf_Vector *x, JXF_Int num_equns, JXF_Int nx, JXF_Int ny )
{
    JXF_Int size = jxf_VectorSize(x);
    JXF_Real *x_data = jxf_VectorData(x);
    jxf_Vector *y = jxf_SeqVectorCreate(size);
    JXF_Real *y_data = NULL;
    JXF_Int i, j, k, fxp, t_i = 0;
    
    jxf_SeqVectorInitialize(y);
    y_data = jxf_VectorData(y);
    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j ++)
        {
            fxp = num_equns * (j * nx + i);
            for (k = 0; k < num_equns; k ++)
            {
                y_data[t_i++] = x_data[fxp++];
            }
        }
    }
    
    return y;
}
