//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ilu.c
 *
 * Created by Yue Xiaoqiang 2014/03/24
 * Xiangtan University
 *
 */

#include "jxf_ilu.h"

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataCreate
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_ILUZeroFactorDataCreate( JXF_Solver *solver, MPI_Comm comm )
{
  *solver = (JXF_Solver)jxf_ILUZeroFactorDataInitialize(comm);
   if (!solver)
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetMaxIter
 * \brief Set max_iter
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
JXF_ILUZeroFactorDataSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_ILUZeroFactorDataSetMaxIter( (void *) solver, max_iter ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetNxy
 * \brief Set nx and ny
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
JXF_ILUZeroFactorDataSetNxy( JXF_Solver solver, JXF_Int nx, JXF_Int ny )
{
   return( jxf_ILUZeroFactorDataSetNxy( (void *) solver, nx, ny ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetNpxy
 * \brief Set npx and npy
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
JXF_ILUZeroFactorDataSetNpxy( JXF_Solver solver, JXF_Int npx, JXF_Int npy )
{
   return( jxf_ILUZeroFactorDataSetNpxy( (void *) solver, npx, npy ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetNumEquns
 * \brief Set num_equns
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
JXF_ILUZeroFactorDataSetNumEquns( JXF_Solver solver, JXF_Int num_equns )
{
   return( jxf_ILUZeroFactorDataSetNumEquns( (void *) solver, num_equns ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetDropTol
 * \brief Set drop_tol
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
JXF_ILUZeroFactorDataSetDropTol( JXF_Solver solver, JXF_Real drop_tol )
{
   return( jxf_ILUZeroFactorDataSetDropTol( (void *) solver, drop_tol ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetMatA
 * \brief Set matA
 * \author Yue Xiaoqiang
 * \date 2014/04/03
 */
JXF_Int
JXF_ILUZeroFactorDataSetMatA( JXF_Solver solver, jxf_CSRMatrix *matA )
{
   return( jxf_ILUZeroFactorDataSetMatA( (void *) solver, matA ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataGetLULength
 * \brief Get lu_length
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
JXF_ILUZeroFactorDataGetLULength( JXF_Solver solver, JXF_Int *lu_length )
{
   return( jxf_ILUZeroFactorDataGetLULength( (void *) solver, lu_length ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataDestroy
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_ILUZeroFactorDataDestroy( JXF_Solver solver )
{
   return( jxf_ILUZeroFactorDataFinalize( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataGenerateParGrid
 * \brief Generate par_grid
 * \author Yue Xiaoqiang
 * \date 2014/10/20
 */
JXF_Int
JXF_ILUZeroFactorDataGenerateParGrid( JXF_Solver solver )
{
   return( jxf_ILUZeroFactorDataGenerateParGrid( (void *) solver ) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataSetup
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_ILUZeroFactorDataSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
   return( jxf_ILUZeroFactorDataSetup( (void *) solver, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) par_matrix )) );
}

/*!
 * \fn JXF_Int JXF_ILUZeroFactorDataPrecond
 * \author Yue Xiaoqiang
 * \date 2014/03/24
 */
JXF_Int
JXF_ILUZeroFactorDataPrecond( JXF_Solver       solver,
                             JXF_hpCSRMatrix par_matrix,
                             JXF_ParVector    par_rhs,
                             JXF_ParVector    par_app  )
{
   return( jxf_ILUZeroFactorDataPrecond( (void *) solver,
                                        jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) par_matrix),
                                        (jxf_ParVector *) par_rhs,
                                        (jxf_ParVector *) par_app ) );
}

/*!
 * \fn void *jxf_ILUZeroFactorDataInitialize
 * \brief Initialize
 * \author Zhiyang Zhou, Yue Xiaoqiang
 * \date 2010/12/11, 2014/03/24
 */
void *
jxf_ILUZeroFactorDataInitialize( MPI_Comm comm )
{
    jxf_ILUZeroFactorData *ilu_data = jxf_CTAlloc(jxf_ILUZeroFactorData, 1);
    JXF_Int max_iter = 200;
    JXF_Real drop_tol = 0.0;
    
    ilu_data->comm = comm;
    ilu_data->max_iter = max_iter;
    ilu_data->drop_tol = drop_tol;
    
    return (void *)ilu_data;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetMaxIter
 * \brief Set max_iter
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
jxf_ILUZeroFactorDataSetMaxIter( void *ilu_vdata, JXF_Int max_iter )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->max_iter = max_iter;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetNxy
 * \brief Set nx and ny
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
jxf_ILUZeroFactorDataSetNxy( void *ilu_vdata, JXF_Int nx, JXF_Int ny )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->nx = nx;
    ilu_data->ny = ny;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetNpxy
 * \brief Set npx and npy
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
jxf_ILUZeroFactorDataSetNpxy( void *ilu_vdata, JXF_Int npx, JXF_Int npy )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->npx = npx;
    ilu_data->npy = npy;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetNumEquns
 * \brief Set num_equns
 * \author Yue Xiaoqiang
 * \date 2014/04/21
 */
JXF_Int
jxf_ILUZeroFactorDataSetNumEquns( void *ilu_vdata, JXF_Int num_equns )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->num_equns = num_equns;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetDropTol
 * \brief Set drop_tol
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
jxf_ILUZeroFactorDataSetDropTol( void *ilu_vdata, JXF_Real drop_tol )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->drop_tol = drop_tol;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetMatA
 * \brief Set matA
 * \author Yue Xiaoqiang
 * \date 2014/04/03
 */
JXF_Int
jxf_ILUZeroFactorDataSetMatA( void *ilu_vdata, jxf_CSRMatrix *matA )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->matA = matA;
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataGetLULength
 * \brief Get lu_length
 * \author Yue Xiaoqiang
 * \date 2014/03/17
 */
JXF_Int
jxf_ILUZeroFactorDataGetLULength( void *ilu_vdata, JXF_Int *lu_length )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    if (!ilu_data)
    {
        jxf_printf(" Warning: ILUZeroFactor object empty!\n");
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
   *lu_length = ilu_data->length;
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataFinalize
 * \brief Finalize
 * \author Zhiyang Zhou, Yue Xiaoqiang
 * \date 2010/12/11, 2014/03/24
 */
JXF_Int
jxf_ILUZeroFactorDataFinalize( void *ilu_vdata )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    
    if (ilu_data)
    {
        if (ilu_data->index)
        {
            jxf_TFree(ilu_data->index);
        }
        if (ilu_data->indexA)
        {
            jxf_TFree(ilu_data->indexA);
        }
        if (ilu_data->indexD)
        {
            jxf_TFree(ilu_data->indexD);
        }
        if (ilu_data->permute)
        {
            jxf_TFree(ilu_data->permute);
        }
        if (ilu_data->value)
        {
            jxf_TFree(ilu_data->value);
        }
        if (ilu_data->senddown)
        {
            jxf_TFree(ilu_data->senddown);
        }
        if (ilu_data->status)
        {
            jxf_TFree(ilu_data->status);
        }
        if (ilu_data->aux_vec)
        {
            jxf_SeqVectorDestroy(ilu_data->aux_vec);
        }
        if (ilu_data->res_vec)
        {
            jxf_SeqVectorDestroy(ilu_data->res_vec);
        }
        if (ilu_data->tmp_vec)
        {
            jxf_SeqVectorDestroy(ilu_data->tmp_vec);
        }
        if (ilu_data->par_aux_vec)
        {
            jxf_ParVectorDestroy(ilu_data->par_aux_vec);
        }
        if (ilu_data->par_res_vec)
        {
            jxf_ParVectorDestroy(ilu_data->par_res_vec);
        }
        if (ilu_data->par_grid)
        {
            jxf_GridPartitionDataFinalize(ilu_data->par_grid);
        }
        jxf_TFree(ilu_data);
    }
    
    return 0;
}

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataGenerateParGrid
 * \brief Generate par_grid
 * \author Yue Xiaoqiang
 * \date 2014/10/20
 */
JXF_Int
jxf_ILUZeroFactorDataGenerateParGrid( void *ilu_vdata )
{
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    ilu_data->par_grid = jxf_GridPartitionDataInitialize(ilu_data->comm,
                     ilu_data->nx, ilu_data->ny, ilu_data->npx, ilu_data->npy);
    jxf_GridPartitionDataSetEachSides4Comm(ilu_data->par_grid);
    return 0;
}
