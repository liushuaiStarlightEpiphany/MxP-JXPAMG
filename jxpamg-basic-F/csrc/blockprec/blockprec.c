//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  blockprec.c
 *  Date: 2021/12/10
 *
 *  Created by dyt
 */

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

#ifndef JXF_BLOCKPREC_HEADER
#include "jxf_blockprec.h"
#endif

/*!
 * \fn JXF_Int JXF_BlockPrecCreate
 */
JXF_Int JXF_BlockPrecCreate(MPI_Comm comm, JXF_Solver *solver)
{
    *solver = ((JXF_Solver)jxf_BlockPrecCreate(comm));

    if (!solver)
    {
        jxf_error_in_arg(2);
    }

    return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_BlockPrecCreate
 */
void *jxf_BlockPrecCreate(MPI_Comm comm)
{
    jxf_BlockPrecData *bp_data = jxf_CTAlloc(jxf_BlockPrecData, 1);

    (bp_data->max_iter) = 2;
    (bp_data->inv_iter) = 2;
    (bp_data->level) = 0;

    jxf_BlockPrecDataDf(bp_data) = NULL;
    jxf_BlockPrecDataDBf(bp_data) = NULL;
    jxf_BlockPrecDataVtemp(bp_data) = NULL;
    jxf_BlockPrecDataUtemp(bp_data) = NULL;

    return (void *)bp_data;
}

JXF_Int JXF_BlockPrecSetMaxIter(JXF_Solver solver, JXF_Int max_iter)
{
    return (jxf_BlockPrecSetMaxIter((void *)solver, max_iter));
}

JXF_Int JXF_BlockPrecGetMaxIter(JXF_Solver solver, JXF_Int *max_iter)
{
    return (jxf_BlockPrecGetMaxIter((void *)solver, max_iter));
}

JXF_Int JXF_BlockPrecSetInvIter(JXF_Solver solver, JXF_Int inv_iter)
{
    return (jxf_BlockPrecSetInvIter((void *)solver, inv_iter));
}

JXF_Int JXF_BlockPrecGetInvIter(JXF_Solver solver, JXF_Int *inv_iter)
{
    return (jxf_BlockPrecGetInvIter((void *)solver, inv_iter));
}

JXF_Int JXF_BlockPrecSetLevel(JXF_Solver solver, JXF_Int level)
{
    return (jxf_BlockPrecSetLevel((void *)solver, level));
}

JXF_Int JXF_BlockPrecGetLevel(JXF_Solver solver, JXF_Int *level)
{
    return (jxf_BlockPrecGetLevel((void *)solver, level));
}

JXF_Int
jxf_BlockPrecSetMaxIter(void *solver, JXF_Int max_iter)
{
    jxf_BlockPrecData *bp_data = solver;
    jxf_BlockPrecDataMaxIter(bp_data) = max_iter;

    return jxf_error_flag;
}

JXF_Int
jxf_BlockPrecGetMaxIter(void *solver, JXF_Int *max_iter)
{
    jxf_BlockPrecData *bp_data = solver;
    *max_iter = jxf_BlockPrecDataMaxIter(bp_data);
    return jxf_error_flag;
}

JXF_Int
jxf_BlockPrecSetInvIter(void *solver, JXF_Int inv_iter)
{
    jxf_BlockPrecData *bp_data = solver;
    jxf_BlockPrecDataInvIter(bp_data) = inv_iter;

    return jxf_error_flag;
}

JXF_Int
jxf_BlockPrecGetInvIter(void *solver, JXF_Int *inv_iter)
{
    jxf_BlockPrecData *bp_data = solver;
    *inv_iter = jxf_BlockPrecDataInvIter(bp_data);
    return jxf_error_flag;
}

JXF_Int
jxf_BlockPrecSetLevel(void *solver, JXF_Int level)
{
    jxf_BlockPrecData *bp_data = solver;
    jxf_BlockPrecDataLevel(bp_data) = level;

    return jxf_error_flag;
}

JXF_Int
jxf_BlockPrecGetLevel(void *solver, JXF_Int *level)
{
    jxf_BlockPrecData *bp_data = solver;
    *level = jxf_BlockPrecDataLevel(bp_data);
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_BlockPrecDestroy
 */
JXF_Int JXF_BlockPrecDestroy(JXF_Solver solver)
{
    return (jxf_BlockPrecDestroy((void *)solver));
}

/*!
 * \fn JXF_Int jxf_BlockPrecDestroy
 */
JXF_Int jxf_BlockPrecDestroy(void *solver)
{
    jxf_BlockPrecData *bp_data = solver;

    jxf_ParVectorDestroy(jxf_BlockPrecDataDf(bp_data));
    jxf_ParVectorDestroy(jxf_BlockPrecDataDBf(bp_data));
    jxf_ParVectorDestroy(jxf_BlockPrecDataVtemp(bp_data));
    jxf_ParVectorDestroy(jxf_BlockPrecDataUtemp(bp_data));

    jxf_TFree(bp_data);
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_BlockPrecSetup
 * \date 2021/12/10
 */
JXF_Int
JXF_BlockPrecSetup(JXF_Solver solver,
                   JXF_hpCSRMatrix hp_matrix)
{
    return (jxf_BlockPrecSetup((void *)solver, (jxf_hpCSRMatrix *)hp_matrix));
}

/*!
 * \fn JXF_Int jxf_BlockPrecSetup
 * \brief Setup phase of BLock preconditioner.
 * \date 2021/12/10
 */
JXF_Int
jxf_BlockPrecSetup(void *solver,
                   jxf_hpCSRMatrix *hp_matrix)
{
    MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
    jxf_BlockPrecData *bp_data = solver;
    jxf_ParVector *Df;
    jxf_ParVector *DBf;
    jxf_ParVector *Vtemp;
    jxf_ParVector *Utemp;
    JXF_Int myid;

    jxf_MPI_Comm_rank(comm, &myid);

    if (myid == 0)
    {
        jxf_printf("bds_maxiter = %d ", jxf_BlockPrecDataMaxIter(bp_data));
        jxf_printf("bds_inviter = %d ", jxf_BlockPrecDataInvIter(bp_data));
        jxf_printf("bds_level = %d ", jxf_BlockPrecDataLevel(bp_data));
    }

    Df = jxf_BlockPrecDataDf(bp_data);
    if (Df != NULL)
    {
        jxf_ParVectorDestroy(Df);
        Df = NULL;
    }
    Df = jxf_ParVectorCreate(comm, jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                             jxf_hpCSRMatrixRowStarts(hp_matrix));
    jxf_ParVectorInitialize(Df);
    jxf_ParVectorSetPartitioningOwner(Df, 0);
    jxf_BlockPrecDataDf(bp_data) = Df;

    DBf = jxf_BlockPrecDataDBf(bp_data);
    if (DBf != NULL)
    {
        jxf_ParVectorDestroy(DBf);
        DBf = NULL;
    }
    DBf = jxf_ParVectorCreate(comm, jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                              jxf_hpCSRMatrixRowStarts(hp_matrix));
    jxf_ParVectorInitialize(DBf);
    jxf_ParVectorSetPartitioningOwner(DBf, 0);
    jxf_BlockPrecDataDBf(bp_data) = DBf;

    Vtemp = jxf_BlockPrecDataVtemp(bp_data);
    if (Vtemp != NULL)
    {
        jxf_ParVectorDestroy(Vtemp);
        Vtemp = NULL;
    }
    Vtemp = jxf_ParVectorCreate(comm, jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                                jxf_hpCSRMatrixRowStarts(hp_matrix));
    jxf_ParVectorInitialize(Vtemp);
    jxf_ParVectorSetPartitioningOwner(Vtemp, 0);
    jxf_BlockPrecDataVtemp(bp_data) = Vtemp;

    Utemp = jxf_BlockPrecDataUtemp(bp_data);
    if (Utemp != NULL)
    {
        jxf_ParVectorDestroy(Utemp);
        Utemp = NULL;
    }
    Utemp = jxf_ParVectorCreate(comm, jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                                jxf_hpCSRMatrixRowStarts(hp_matrix));
    jxf_ParVectorInitialize(Utemp);
    jxf_ParVectorSetPartitioningOwner(Utemp, 0);
    jxf_BlockPrecDataUtemp(bp_data) = Utemp;

    return 0;
}

/*!
 * \fn JXF_Int JXF_BlockPrec_JAC
 * \brief Solver phase of Block preconditioner. (using Jacobi as inner iteration)
 * \date 2021/12/10
 */
JXF_Int
JXF_BlockPrec_JAC(JXF_Solver solver,
                  JXF_hpCSRMatrix hp_matrix,
                  JXF_ParVector par_rhs,
                  JXF_ParVector par_app)
{
    return (jxf_BlockPrec_JAC((void *)solver,
                              (jxf_hpCSRMatrix *)hp_matrix,
                              (jxf_ParVector *)par_rhs,
                              (jxf_ParVector *)par_app));
}

/*!
 * \fn JXF_Int JXF_BlockPrec_GS
 * \brief Solver phase of Block preconditioner. (using GS as inner iteration -- core level)
 * \date 2022/5/18
 */
JXF_Int
JXF_BlockPrec_GS(JXF_Solver solver,
                 JXF_hpCSRMatrix hp_matrix,
                 JXF_ParVector par_rhs,
                 JXF_ParVector par_app)
{
    return (jxf_BlockPrec_GS((void *)solver,
                             (jxf_hpCSRMatrix *)hp_matrix,
                             (jxf_ParVector *)par_rhs,
                             (jxf_ParVector *)par_app));
}

/*!
 * \fn JXF_Int JXF_BlockPrec_SGS
 * \brief Solver phase of Block preconditioner. (using SGS as inner iteration -- core level)
 * \date 2022/5/28
 */
JXF_Int
JXF_BlockPrec_SGS(JXF_Solver solver,
                  JXF_hpCSRMatrix hp_matrix,
                  JXF_ParVector par_rhs,
                  JXF_ParVector par_app)
{
    return (jxf_BlockPrec_SGS((void *)solver,
                              (jxf_hpCSRMatrix *)hp_matrix,
                              (jxf_ParVector *)par_rhs,
                              (jxf_ParVector *)par_app));
}

/*!
 * \fn JXF_Int jxf_BlockPrec_JAC
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2021/12/10
 */

JXF_Int
jxf_BlockPrec_JAC(void *solver,
                  jxf_hpCSRMatrix *hp_matrix,
                  jxf_ParVector *par_rhs,
                  jxf_ParVector *par_app)
{
    jxf_BlockPrecData *bp_data = solver;
    JXF_Int bp_maxiter = jxf_BlockPrecDataMaxIter(bp_data);
    JXF_Int bp_inviter = jxf_BlockPrecDataInvIter(bp_data);
    JXF_Int bp_level = jxf_BlockPrecDataLevel(bp_data);

    jxf_CSRMatrix *A_diag = jxf_hpCSRMatrixDiag(hp_matrix);
    JXF_Int *A_diag_i = jxf_CSRMatrixI(A_diag);
    JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);

    jxf_ParCSRMatrix *A_B = jxf_hpMatrixLevelToPar(hp_matrix, bp_level);

    jxf_ParCSRMatrix *A = jxf_hpCSRMatrixPar(hp_matrix);

    jxf_Vector *f_local = jxf_ParVectorLocalVector(par_rhs);
    JXF_Real *f_data = jxf_VectorData(f_local);

    jxf_ParVector *Df = NULL;
    jxf_ParVector *DBf = NULL;
    jxf_ParVector *Vtemp = NULL;
    jxf_ParVector *Utemp = NULL;

    Df = jxf_BlockPrecDataDf(bp_data);
    DBf = jxf_BlockPrecDataDBf(bp_data);
    Vtemp = jxf_BlockPrecDataVtemp(bp_data);
    Utemp = jxf_BlockPrecDataUtemp(bp_data);

    jxf_Vector *Df_local = jxf_ParVectorLocalVector(Df);
    JXF_Real *Df_data = jxf_VectorData(Df_local);
    jxf_Vector *DBf_local = jxf_ParVectorLocalVector(DBf);

    jxf_Vector *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
    JXF_Real *Vtemp_data = jxf_VectorData(Vtemp_local);
    jxf_Vector *Utemp_local = jxf_ParVectorLocalVector(Utemp);
    JXF_Real *Utemp_data = jxf_VectorData(Utemp_local);

    JXF_Int i, it, invit, ierr = 0;
    JXF_Int local_size = jxf_VectorSize(f_local);

#if 1
    /* DBf = D_B^{-1}*f */
    for (i = 0; i < local_size; i++)
    {
        /* Df = D^{-1}*f */
        Df_data[i] = f_data[i] / A_diag_data[A_diag_i[i]];
    }
    jxf_ParVectorCopy(Df, DBf);
    for (invit = 1; invit < bp_inviter; invit++)
    {
        jxf_ParCSRMatrixMatvec(1.0, A_B, DBf, 0.0, Vtemp);
        for (i = 0; i < local_size; i++)
        {
            Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
        }

        jxf_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
        jxf_SeqVectorAxpy(1.0, Df_local, DBf_local);
    }

    jxf_ParVectorCopy(DBf, par_app);
    for (it = 1; it < bp_maxiter; it++)
    {
        /* temp_1 = A * u */
        jxf_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

        /* temp_1 = D_B^{-1} * A*u */
        for (i = 0; i < local_size; i++)
        {
            Df_data[i] = Utemp_data[i] / A_diag_data[A_diag_i[i]]; /* Df = D^{-1} * A*u  */
        }
        jxf_SeqVectorCopy(Df_local, Utemp_local);
        for (invit = 1; invit < bp_inviter; invit++)
        {
            /* Vtemp = A_B * Utemp */
            jxf_ParCSRMatrixMatvec(1.0, A_B, Utemp, 0.0, Vtemp);

            /* Vtemp = D^{-1} * Vtemp */
            for (i = 0; i < local_size; i++)
            {
                Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
            }

            /* Utemp = -Vtemp + Utemp + Df */
            jxf_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
            jxf_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }
        /* u = -Utemp + u */
        jxf_ParVectorAxpy(-1.0, Utemp, par_app);

        /* u = DBf + u */
        jxf_ParVectorAxpy(1.0, DBf, par_app);
    }
#endif

    jxf_TFree(A_B);

// BlockPrec//
#if 0
         /* Df = D^{-1}*f */
         for(i = 0; i < local_size; i++)
         {
             Df_data[i] = f_data[i] / A_diag_data[A_diag_i[i]];
         }
 
         /* temp = Df */
         ierr = jxf_ParVectorCopy(Df, par_app);
 
         for(it = 1; it < bp_maxiter; it++)
         {   
             /* temp = u */
             jxf_ParVectorCopy(par_app, Utemp);
 
             /* temp = A*u */
             ierr = jxf_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);
 
             /* temp = D^{-1}*temp */
             for(i = 0; i < local_size; i++)
             {
                 Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
             }
 
             /* u = -temp + u */
             ierr = jxf_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
 
             /* u = Df + u */
             ierr = jxf_SeqVectorAxpy(1.0, Df_local, Utemp_local);
         }
#endif

    return ierr;
}

/*!
 * \fn JXF_Int jxf_BlockPrec_GS
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2022/5/18
 */
JXF_Int
jxf_BlockPrec_GS(void *solver,
                 jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector *par_rhs,
                 jxf_ParVector *par_app)
{
    jxf_BlockPrecData *bp_data = solver;
    JXF_Int bp_maxiter = jxf_BlockPrecDataMaxIter(bp_data);
    JXF_Int bp_inviter = jxf_BlockPrecDataInvIter(bp_data);

    jxf_CSRMatrix *A_diag = jxf_hpCSRMatrixDiag(hp_matrix);
    JXF_Int *A_diag_i = jxf_CSRMatrixI(A_diag);
    JXF_Int *A_diag_j = jxf_CSRMatrixJ(A_diag);
    JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);

    JXF_Int n = jxf_CSRMatrixNumRows(A_diag);

    jxf_ParCSRMatrix *A = jxf_hpCSRMatrixPar(hp_matrix);

    jxf_Vector *f_local = jxf_ParVectorLocalVector(par_rhs);

    jxf_ParVector *Df = NULL;
    jxf_ParVector *DBf = NULL;
    jxf_ParVector *Vtemp = NULL;
    jxf_ParVector *Utemp = NULL;

    Df = jxf_BlockPrecDataDf(bp_data);
    DBf = jxf_BlockPrecDataDBf(bp_data);
    Vtemp = jxf_BlockPrecDataVtemp(bp_data);
    Utemp = jxf_BlockPrecDataUtemp(bp_data);

    jxf_Vector *Df_local = jxf_ParVectorLocalVector(Df);
    JXF_Real *Df_data = jxf_VectorData(Df_local);
    jxf_Vector *DBf_local = jxf_ParVectorLocalVector(DBf);

    jxf_Vector *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
    JXF_Real *Vtemp_data = jxf_VectorData(Vtemp_local);
    jxf_Vector *Utemp_local = jxf_ParVectorLocalVector(Utemp);
    JXF_Real *Utemp_data = jxf_VectorData(Utemp_local);

    JXF_Int i, it, invit, ierr = 0;
    JXF_Int ii, jj;

    JXF_Real zero = 0.0;
    JXF_Real res;

    /* Df = (D+L)^{-1}*f */
    jxf_SeqVectorCopy(f_local, Df_local);
    for (i = 0; i < n; i++)
    {
        if (A_diag_data[A_diag_i[i]] != zero)
        {
            res = Df_data[i];
            for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
            {
                ii = A_diag_j[jj];
                if (ii < i)
                {
                    res -= A_diag_data[jj] * Df_data[ii];
                }
            }
            Df_data[i] = res / A_diag_data[A_diag_i[i]];
        }
    }
    jxf_SeqVectorCopy(Df_local, DBf_local);
    for (invit = 1; invit < bp_inviter; invit++)
    {
        /*  Vtemp = A*DBf */
        jxf_CSRMatrixMatvec_origin(1.0, A_diag, DBf_local, 0.0, Vtemp_local);
        /*  Vtemp = (D+L)^{-1} * Vtemp */
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Vtemp_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii < i)
                    {
                        res -= A_diag_data[jj] * Vtemp_data[ii];
                    }
                }
                Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }
        /*  DBf = DBf - Vtemp + Df */
        jxf_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
        jxf_SeqVectorAxpy(1.0, Df_local, DBf_local);
    }

    jxf_ParVectorCopy(DBf, par_app);
    for (it = 1; it < bp_maxiter; it++)
    {
        /* Utemp = A * u */
        jxf_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

        /* Utemp = D_B^{-1} * A * u */
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Utemp_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii < i)
                    {
                        res -= A_diag_data[jj] * Df_data[ii];
                    }
                }
                Df_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }
        jxf_SeqVectorCopy(Df_local, Utemp_local);
        for (invit = 1; invit < bp_inviter; invit++)
        {
            /*  Vtemp = A * Utemp */
            jxf_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);
            /*  Vtemp = (D+L)^{-1} * Vtemp */
            for (i = 0; i < n; i++)
            {
                if (A_diag_data[A_diag_i[i]] != zero)
                {
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                    {
                        ii = A_diag_j[jj];
                        if (ii < i)
                        {
                            res -= A_diag_data[jj] * Vtemp_data[ii];
                        }
                    }
                    Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
                }
            }
            /*  Utemp = Utemp - Vtemp + Df */
            jxf_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
            jxf_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }

        /* u = -Utemp + u */
        jxf_ParVectorAxpy(-1.0, Utemp, par_app);

        /* u = DBf + u */
        jxf_ParVectorAxpy(1.0, DBf, par_app);
    }

    return ierr;
}

/*!
 * \fn JXF_Int jxf_BlockPrec_SGS
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2022/6/8
 */
JXF_Int
jxf_BlockPrec_SGS(void *solver,
                  jxf_hpCSRMatrix *hp_matrix,
                  jxf_ParVector *par_rhs,
                  jxf_ParVector *par_app)
{
    jxf_BlockPrecData *bp_data = solver;
    JXF_Int bp_maxiter = jxf_BlockPrecDataMaxIter(bp_data);
    JXF_Int bp_inviter = jxf_BlockPrecDataInvIter(bp_data);

    jxf_CSRMatrix *A_diag = jxf_hpCSRMatrixDiag(hp_matrix);
    JXF_Int *A_diag_i = jxf_CSRMatrixI(A_diag);
    JXF_Int *A_diag_j = jxf_CSRMatrixJ(A_diag);
    JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);

    JXF_Int n = jxf_CSRMatrixNumRows(A_diag);

    jxf_ParCSRMatrix *A = jxf_hpCSRMatrixPar(hp_matrix);

    jxf_Vector *f_local = jxf_ParVectorLocalVector(par_rhs);

    jxf_ParVector *Df = NULL;
    jxf_ParVector *DBf = NULL;
    jxf_ParVector *Vtemp = NULL;
    jxf_ParVector *Utemp = NULL;

    Df = jxf_BlockPrecDataDf(bp_data);
    DBf = jxf_BlockPrecDataDBf(bp_data);
    Vtemp = jxf_BlockPrecDataVtemp(bp_data);
    Utemp = jxf_BlockPrecDataUtemp(bp_data);

    jxf_Vector *Df_local = jxf_ParVectorLocalVector(Df);
    JXF_Real *Df_data = jxf_VectorData(Df_local);
    jxf_Vector *DBf_local = jxf_ParVectorLocalVector(DBf);

    jxf_Vector *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
    JXF_Real *Vtemp_data = jxf_VectorData(Vtemp_local);
    jxf_Vector *Utemp_local = jxf_ParVectorLocalVector(Utemp);
    JXF_Real *Utemp_data = jxf_VectorData(Utemp_local);

    JXF_Int i, it, invit, ierr = 0;
    JXF_Int ii, jj;

    JXF_Real zero = 0.0;
    JXF_Real res;

    /* Df = (D+U)^{-1}*D*(D+L)^{-1}*f */
    /* Df = (D+L)^{-1}*f  */
    jxf_SeqVectorCopy(f_local, Df_local);
    for (i = 0; i < n; i++)
    {
        if (A_diag_data[A_diag_i[i]] != zero)
        {
            res = Df_data[i];
            for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
            {
                ii = A_diag_j[jj];
                if (ii < i)
                {
                    res -= A_diag_data[jj] * Df_data[ii];
                }
            }
            Df_data[i] = res / A_diag_data[A_diag_i[i]];
        }
    }
    /* Df = D * Df */
    for (i = 0; i < n; i++)
    {
        if (A_diag_data[A_diag_i[i]] != zero)
        {
            Df_data[i] *= A_diag_data[A_diag_i[i]];
        }
    }
    /*  Df = (D+U)^{-1} * Df  */
    for (i = n - 1; i >= 0; i--)
    {
        if (A_diag_data[A_diag_i[i]] != zero)
        {
            res = Df_data[i];
            for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
            {
                ii = A_diag_j[jj];
                if (ii > i)
                {
                    res -= A_diag_data[jj] * Df_data[ii];
                }
            }
            Df_data[i] = res / A_diag_data[A_diag_i[i]];
        }
    }
    jxf_SeqVectorCopy(Df_local, DBf_local);

    for (invit = 1; invit < bp_inviter; invit++)
    {
        /*  Vtemp = A*DBf */
        jxf_CSRMatrixMatvec_origin(1.0, A_diag, DBf_local, 0.0, Vtemp_local);
        /*  Vtemp = (D+L)^{-1} * Vtemp */
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Vtemp_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii < i)
                    {
                        res -= A_diag_data[jj] * Vtemp_data[ii];
                    }
                }
                Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }
        /*  Vtemp = D * Vtemp */
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                Vtemp_data[i] *= A_diag_data[A_diag_i[i]];
            }
        }
        /*  Vtemp = (D+U)^{-1} * Vtemp */
        for (i = n - 1; i >= 0; i--)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Vtemp_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii > i)
                    {
                        res -= A_diag_data[jj] * Vtemp_data[ii];
                    }
                }
                Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }

        /*  DBf = DBf - Vtemp + Df */
        jxf_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
        jxf_SeqVectorAxpy(1.0, Df_local, DBf_local);
    }

    jxf_ParVectorCopy(DBf, par_app);
    for (it = 1; it < bp_maxiter; it++)
    {
        /* Utemp = A * u */
        jxf_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

        /* Utemp = D_B^{-1} * A * u */
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Utemp_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii < i)
                    {
                        res -= A_diag_data[jj] * Df_data[ii];
                    }
                }
                Df_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }
        for (i = 0; i < n; i++)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                Df_data[i] *= A_diag_data[A_diag_i[i]];
            }
        }
        for (i = n - 1; i >= 0; i--)
        {
            if (A_diag_data[A_diag_i[i]] != zero)
            {
                res = Df_data[i];
                for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                {
                    ii = A_diag_j[jj];
                    if (ii > i)
                    {
                        res -= A_diag_data[jj] * Df_data[ii];
                    }
                }
                Df_data[i] = res / A_diag_data[A_diag_i[i]];
            }
        }
        jxf_SeqVectorCopy(Df_local, Utemp_local);
        for (invit = 1; invit < bp_inviter; invit++)
        {
            /*  Vtemp = A*Utemp + Vtemp */
            jxf_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);
            /*  Vtemp = (D_B)^{-1} * Vtemp */
            for (i = 0; i < n; i++)
            {
                if (A_diag_data[A_diag_i[i]] != zero)
                {
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                    {
                        ii = A_diag_j[jj];
                        if (ii < i)
                        {
                            res -= A_diag_data[jj] * Vtemp_data[ii];
                        }
                    }
                    Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
                }
            }
            for (i = 0; i < n; i++)
            {
                if (A_diag_data[A_diag_i[i]] != zero)
                {
                    Vtemp_data[i] *= A_diag_data[A_diag_i[i]];
                }
            }
            for (i = n - 1; i >= 0; i--)
            {
                if (A_diag_data[A_diag_i[i]] != zero)
                {
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i] + 1; jj < A_diag_i[i + 1]; jj++)
                    {
                        ii = A_diag_j[jj];
                        if (ii > i)
                        {
                            res -= A_diag_data[jj] * Vtemp_data[ii];
                        }
                    }
                    Vtemp_data[i] = res / A_diag_data[A_diag_i[i]];
                }
            }
            /*  Utemp = Utemp - Vtemp + Df */
            jxf_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
            jxf_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }

        /* u = -Utemp + u */
        jxf_ParVectorAxpy(-1.0, Utemp, par_app);

        /* u = DBf + u */
        jxf_ParVectorAxpy(1.0, DBf, par_app);
    }

    return ierr;
}
