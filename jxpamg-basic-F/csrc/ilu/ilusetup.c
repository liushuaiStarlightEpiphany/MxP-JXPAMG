//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * ilusetup.c
 *
 * Created by Yue Xiaoqiang 2014/03/24
 * Xiangtan University
 *
 */

#include "jxf_ilu.h"

/*!
 * \fn JXF_Int jxf_ILUZeroFactorDataSetup
 * \brief SETUP phase
 * \author peghoty, Yue Xiaoqiang
 * \date 2010/12/11, 2014/03/24
 */
JXF_Int
jxf_ILUZeroFactorDataSetup( void *ilu_vdata, jxf_ParCSRMatrix *par_A )
{
    MPI_Comm comm = jxf_ParCSRMatrixComm(par_A);
    jxf_ILUZeroFactorData *ilu_data = ilu_vdata;
    JXF_Real drop_tol = ilu_data->drop_tol;
    JXF_Int *index = NULL;
    JXF_Int *indexA = NULL;
    JXF_Int *indexD = NULL;
    JXF_Real *value = NULL;
    jxf_CSRMatrix *A = NULL;
    JXF_Int *part_aux = NULL;
    JXF_Int *part_res = NULL;
    JXF_Int num_procs, my_id, num_rows;
    JXF_Int length, dfill_in, maxlength, maxdfill_in;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    if (num_procs == 1)
    {
        A = jxf_ParCSRMatrixDiag(par_A);
        num_rows = jxf_CSRMatrixNumRows(A);
        //jxf_printf(" >>> SETUP: Drop-Tol = %le on ID = %d of %d\n", drop_tol, my_id, num_procs);
        if (drop_tol == 0.0)
        {
            ilu_data->num_fill_in_drop = jxf_ILUZeroDecompositionA(A, &indexD, &index, &value);
            ilu_data->length = jxf_CSRMatrixNumNonzeros(A);
        }
        else
        {
            ilu_data->num_fill_in_drop = jxf_ILUZeroDecompositionB(A, drop_tol, &indexA, &indexD, &index, &value);
            ilu_data->length = indexA[num_rows];
        }
        ilu_data->index = index;
        ilu_data->indexA = indexA;
        ilu_data->indexD = indexD;
        ilu_data->value = value;
        
        ilu_data->aux_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->aux_vec);
        ilu_data->res_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->res_vec);
        ilu_data->tmp_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->tmp_vec);
    }
    else if ((ilu_data->matA != NULL) && (my_id == 0))
    {
        A = ilu_data->matA;
        num_rows = jxf_CSRMatrixNumRows(A);
        //jxf_printf(" >>> SETUP: Drop-Tol = %le on ID = %d of %d\n", drop_tol, my_id, num_procs);
        if (drop_tol == 0.0)
        {
            ilu_data->num_fill_in_drop = jxf_ILUZeroDecompositionA(A, &indexD, &index, &value);
            ilu_data->length = jxf_CSRMatrixNumNonzeros(A);
        }
        else
        {
            ilu_data->num_fill_in_drop = jxf_ILUZeroDecompositionB(A, drop_tol, &indexA, &indexD, &index, &value);
            ilu_data->length = indexA[num_rows];
        }
        ilu_data->index = index;
        ilu_data->indexA = indexA;
        ilu_data->indexD = indexD;
        ilu_data->value = value;
        
        ilu_data->aux_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->aux_vec);
        ilu_data->res_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->res_vec);
        ilu_data->tmp_vec = jxf_SeqVectorCreate(num_rows);
        jxf_SeqVectorInitialize(ilu_data->tmp_vec);
    }
    else if (drop_tol == 0.0)
    {
        //jxf_printf(" >>> SETUP: Drop-Tol = %le on ID = %d of %d\n", drop_tol, my_id, num_procs);
        indexA = jxf_ILUZeroParallelDecompositionA(par_A, &indexD, &index, &value, &length, &dfill_in, ilu_data);
        jxf_MPI_Reduce(&length, &maxlength, 1, JXF_MPI_INT, MPI_MAX, 0, comm);
        jxf_MPI_Reduce(&dfill_in, &maxdfill_in, 1, JXF_MPI_INT, MPI_MAX, 0, comm);
        ilu_data->index = index;
        ilu_data->indexA = indexA;
        ilu_data->indexD = indexD;
        ilu_data->value = value;
        ilu_data->length = maxlength;
        ilu_data->num_fill_in_drop = maxdfill_in;
        
        // Allocate memory for auxiliary vector
        jxf_ParCSRMatrixGetRowPartitioning(par_A, &part_aux);
        ilu_data->par_aux_vec = jxf_ParVectorCreate(comm, jxf_ParCSRMatrixGlobalNumRows(par_A), part_aux);
        jxf_ParVectorInitialize(ilu_data->par_aux_vec);
        // Allocate memory for residual vector
        jxf_ParCSRMatrixGetRowPartitioning(par_A, &part_res);
        ilu_data->par_res_vec = jxf_ParVectorCreate(comm, jxf_ParCSRMatrixGlobalNumRows(par_A), part_res);
        jxf_ParVectorInitialize(ilu_data->par_res_vec);
    }
    else
    {
        jxf_printf("\n >>> Wrong parameters in ILU Setup\n");
        exit(0);
    }
    
    return 0;
}

/*------------------------------------------------------------*/
/*------------------------ End of File -----------------------*/
/*------------------------------------------------------------*/
