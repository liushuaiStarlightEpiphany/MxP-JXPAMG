//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  get_row_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

#if defined(JXFPAMG_GET_ROW)

#undef __FUNC__
#define __FUNC__ "jxf_EuclidGetRow (JXFPAMG_GET_ROW)"
void jxf_EuclidGetRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_Int ierr;
    jxf_ParCSRMatrix *mat = (jxf_ParCSRMatrix *)A;
    
    ierr = jxf_ParCSRMatrixGetRow(mat, row, len, ind, val);
    if (ierr)
    {
        jxf_sprintf(jxf_msgBuf_dh, "jxf_ParCSRMatrixGetRow(row= %i) returned %i", row+1, ierr);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidRestoreRow (JXFPAMG_GET_ROW)"
void jxf_EuclidRestoreRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_Int ierr;
    jxf_ParCSRMatrix *mat = (jxf_ParCSRMatrix *)A;
    
    ierr = jxf_ParCSRMatrixRestoreRow(mat, row, len, ind, val);
    if (ierr)
    {
        jxf_sprintf(jxf_msgBuf_dh, "jxf_ParCSRMatrixRestoreRow(row= %i) returned %i", row+1, ierr);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidGetDimensions (JXFPAMG_GET_ROW)"
void jxf_EuclidGetDimensions( void *A, JXF_Int *beg_row, JXF_Int *rowsLocal, JXF_Int *rowsGlobal )
{
    JXF_START_FUNC_DH
    JXF_Int ierr, m, n;
    JXF_Int row_start, row_end, col_start, col_end;
    jxf_ParCSRMatrix *mat = (jxf_ParCSRMatrix *)A;
    
    ierr = jxf_ParCSRMatrixGetDims(mat, &m, &n);
    if (ierr)
    {
        jxf_sprintf(jxf_msgBuf_dh, "jxf_ParCSRMatrixGetDims() returned %i", ierr);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    ierr = jxf_ParCSRMatrixGetLocalRange(mat, &row_start, &row_end, &col_start, &col_end);
    if (ierr)
    {
        jxf_sprintf(jxf_msgBuf_dh, "jxf_ParCSRMatrixGetLocalRange() returned %i", ierr);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
   *beg_row = row_start;
   *rowsLocal = (row_end - row_start + 1);
   *rowsGlobal = n;
    JXF_END_FUNC_DH
}

#else

#undef __FUNC__
#define __FUNC__ "jxf_EuclidGetRow (ERROR)"
void jxf_EuclidGetRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_SET_ERROR(EUCLID_ERROR, "Oops; missing XXX_GET_ROW definition!");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidRestoreRow (ERROR)"
void jxf_EuclidRestoreRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val )
{
    JXF_START_FUNC_DH
    JXF_SET_ERROR(EUCLID_ERROR, "Oops; missing XXX_GET_ROW definition!");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidGetDimensions (ERROR)"
void jxf_EuclidGetDimensions( void *A, JXF_Int *beg_row, JXF_Int *rowsLocal, JXF_Int *rowsGlobal )
{
    JXF_START_FUNC_DH
    JXF_SET_ERROR(EUCLID_ERROR, "Oops; missing XXX_GET_ROW definition!");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidReadLocalNz (ERROR)"
JXF_Int jxf_EuclidReadLocalNz( void *A )
{
    JXF_START_FUNC_DH
    JXF_SET_ERROR(EUCLID_ERROR, "Oops; missing XXX_GET_ROW definition!");
    JXF_END_FUNC_DH
}

#endif

#undef __FUNC__
#define __FUNC__ "jxf_PrintMatUsingGetRow"
void jxf_PrintMatUsingGetRow( void* A, JXF_Int beg_row, JXF_Int m, JXF_Int *n2o_row, JXF_Int *n2o_col, char *filename )
{
    JXF_START_FUNC_DH
    FILE *fp;
    JXF_Int *o2n_col = NULL, pe, i, j, *cval, len;
    JXF_Int newCol, newRow;
    JXF_Real *aval;
    
    if (n2o_col != NULL)
    {
        o2n_col = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        for (i = 0; i < m; ++ i) o2n_col[n2o_col[i]] = i;
    }
    for (pe = 0; pe < jxf_np_dh; ++ pe)
    {
        jxf_MPI_Barrier(jxf_comm_dh);
        if (jxf_myid_dh == pe)
        {
            if (pe == 0)
            {
                fp = fopen(filename, "w");
            }
            else
            {
                fp = fopen(filename, "a");
            }
            if (fp == NULL)
            {
                jxf_sprintf(jxf_msgBuf_dh, "can't open %s for writing\n", filename);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
            for (i = 0; i < m; ++ i)
            {
                if (n2o_row == NULL)
                {
                    jxf_EuclidGetRow(A, i+beg_row, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                    for (j = 0; j < len; ++ j)
                    {
                        jxf_fprintf(fp, "%i %i %g\n", i+1, cval[j], aval[j]);
                    }
                    jxf_EuclidRestoreRow(A, i, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                }
                else
                {
                    newRow = n2o_row[i] + beg_row;
                    jxf_EuclidGetRow(A, newRow, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                    for (j = 0; j < len; ++ j)
                    {
                        newCol = o2n_col[cval[j]-beg_row] + beg_row;
                        jxf_fprintf(fp, "%i %i %g\n", i+1, newCol, aval[j]);
                    }
                    jxf_EuclidRestoreRow(A, i, &len, &cval, &aval); JXF_CHECK_V_ERROR;
                }
            }
            fclose(fp);
        }
    }
    if (n2o_col != NULL)
    {
        JXF_FREE_DH(o2n_col); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhInputJXFpamgMat"
void jxf_Euclid_dhInputJXFpamgMat( jxf_Euclid_dh ctx, JXF_ParCSRMatrix A )
{
    JXF_START_FUNC_DH
    JXF_Int M, N;
    JXF_Int beg_row, end_row, junk;
    
    jxf_ParCSRMatrixGetDims((jxf_ParCSRMatrix *)A, &M, &N);
    if (M != N)
    {
        jxf_sprintf(jxf_msgBuf_dh, "Global matrix is not square: M= %i, N= %i", M, N);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    jxf_ParCSRMatrixGetLocalRange((jxf_ParCSRMatrix *)A, &beg_row, &end_row, &junk, &junk);
    ctx->m = end_row - beg_row + 1;
    ctx->n = M;
    ctx->A = (void *)A;
    JXF_END_FUNC_DH
}
