//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  blockprec.c
 *  Date: 2021/12/10
 *  
 *  Created by dyt
 */ 

#ifndef JX_HPCSRMV_HEADER 
#include "jx_hpcsr.h"
#endif

#ifndef JX_BLOCKPREC_HEADER 
#include "jx_blockprec.h"
#endif

/*!
 * \fn JX_Int JX_BlockPrecCreate
 */ 
JX_Int JX_BlockPrecCreate( MPI_Comm comm, JX_Solver *solver )
{
    *solver = ( (JX_Solver)jx_BlockPrecCreate(comm) );

    if (!solver) 
    {
      jx_error_in_arg(2);
    }

    return jx_error_flag;
}

/*!
 * \fn JX_Int jx_BlockPrecCreate
 */
void * jx_BlockPrecCreate(MPI_Comm comm)
{
    jx_BlockPrecData *bp_data = jx_CTAlloc(jx_BlockPrecData, 1);

    (bp_data -> max_iter)  = 2; 
    (bp_data -> inv_iter)  = 2;
    (bp_data -> level)     = 0;

    jx_BlockPrecDataDf(bp_data)    = NULL;
    jx_BlockPrecDataDBf(bp_data)   = NULL;
    jx_BlockPrecDataVtemp(bp_data) = NULL;
    jx_BlockPrecDataUtemp(bp_data) = NULL;

    return (void *) bp_data;
}   

JX_Int JX_BlockPrecSetMaxIter( JX_Solver solver, JX_Int max_iter)
{
    return( jx_BlockPrecSetMaxIter( (void *) solver, max_iter) );
}

JX_Int JX_BlockPrecGetMaxIter( JX_Solver solver, JX_Int *max_iter)
{
    return( jx_BlockPrecGetMaxIter( (void *) solver, max_iter) );
}

JX_Int JX_BlockPrecSetInvIter( JX_Solver solver, JX_Int inv_iter)
{
    return( jx_BlockPrecSetInvIter( (void *) solver, inv_iter) );
}

JX_Int JX_BlockPrecGetInvIter( JX_Solver solver, JX_Int *inv_iter)
{
    return( jx_BlockPrecGetInvIter( (void *) solver, inv_iter) );
}

JX_Int JX_BlockPrecSetLevel( JX_Solver solver, JX_Int level)
{
    return( jx_BlockPrecSetLevel( (void *) solver, level) );
}

JX_Int JX_BlockPrecGetLevel( JX_Solver solver, JX_Int *level)
{
    return( jx_BlockPrecGetLevel( (void *) solver, level) );
}

JX_Int
jx_BlockPrecSetMaxIter( void *solver, JX_Int max_iter )
{
    jx_BlockPrecData       *bp_data  = solver;
    jx_BlockPrecDataMaxIter(bp_data) = max_iter;

    return jx_error_flag;
}

JX_Int
jx_BlockPrecGetMaxIter( void *solver, JX_Int *max_iter )
{
   jx_BlockPrecData *bp_data = solver;
   *max_iter = jx_BlockPrecDataMaxIter(bp_data);
   return jx_error_flag;
}

JX_Int
jx_BlockPrecSetInvIter( void *solver, JX_Int inv_iter )
{
    jx_BlockPrecData       *bp_data  = solver;
    jx_BlockPrecDataInvIter(bp_data) = inv_iter;

    return jx_error_flag;
}

JX_Int
jx_BlockPrecGetInvIter( void *solver, JX_Int *inv_iter )
{
   jx_BlockPrecData *bp_data = solver;
   *inv_iter = jx_BlockPrecDataInvIter(bp_data);
   return jx_error_flag;
}

JX_Int
jx_BlockPrecSetLevel( void *solver, JX_Int level )
{
    jx_BlockPrecData     *bp_data  = solver;
    jx_BlockPrecDataLevel(bp_data) = level;

    return jx_error_flag;
}

JX_Int
jx_BlockPrecGetLevel( void *solver, JX_Int *level )
{
   jx_BlockPrecData *bp_data = solver;
   *level = jx_BlockPrecDataLevel(bp_data);
   return jx_error_flag;
}


/*!
 * \fn JX_Int JX_BlockPrecDestroy
 */
JX_Int JX_BlockPrecDestroy( JX_Solver solver)
{
    return( jx_BlockPrecDestroy( (void *) solver));
}

/*!
 * \fn JX_Int jx_BlockPrecDestroy
 */
JX_Int jx_BlockPrecDestroy( void *solver )
{
    jx_BlockPrecData *bp_data  = solver;

    jx_ParVectorDestroy(jx_BlockPrecDataDf(bp_data));
    jx_ParVectorDestroy(jx_BlockPrecDataDBf(bp_data));
    jx_ParVectorDestroy(jx_BlockPrecDataVtemp(bp_data));
    jx_ParVectorDestroy(jx_BlockPrecDataUtemp(bp_data));

    jx_TFree(bp_data);
    return jx_error_flag;
}

/*!
 * \fn JX_Int JX_BlockPrecSetup
 * \date 2021/12/10
 */
JX_Int 
JX_BlockPrecSetup( JX_Solver       solver, 
                   JX_hpCSRMatrix  hp_matrix )
{
    return ( jx_BlockPrecSetup( (void *) solver, (jx_hpCSRMatrix *) hp_matrix));
}   

/*!
 * \fn JX_Int jx_BlockPrecSetup
 * \brief Setup phase of BLock preconditioner.
 * \date 2021/12/10
 */
JX_Int
jx_BlockPrecSetup( void            *solver, 
                   jx_hpCSRMatrix  *hp_matrix )
{
    MPI_Comm          comm    = jx_hpCSRMatrixComm(hp_matrix);
    jx_BlockPrecData *bp_data = solver;
    jx_ParVector   *Df;
    jx_ParVector   *DBf;
    jx_ParVector   *Vtemp;
    jx_ParVector   *Utemp;
    JX_Int          myid;

    jx_MPI_Comm_rank(comm,&myid);

    if(myid == 0)
    {
        jx_printf("bds_maxiter = %d ", jx_BlockPrecDataMaxIter(bp_data));
        jx_printf("bds_inviter = %d ", jx_BlockPrecDataInvIter(bp_data));
        jx_printf("bds_level = %d ", jx_BlockPrecDataLevel(bp_data));
    }

    Df = jx_BlockPrecDataDf(bp_data);
    if (Df != NULL)
    {
        jx_ParVectorDestroy(Df);
        Df = NULL;
    }
    Df = jx_ParVectorCreate(comm, jx_hpCSRMatrixGlobalNumRows(hp_matrix), 
                                  jx_hpCSRMatrixRowStarts(hp_matrix));
    jx_ParVectorInitialize(Df);
    jx_ParVectorSetPartitioningOwner(Df, 0);
    jx_BlockPrecDataDf(bp_data) = Df;

    DBf = jx_BlockPrecDataDBf(bp_data);
    if (DBf != NULL)
    {
        jx_ParVectorDestroy(DBf);
        DBf = NULL;
    }
    DBf = jx_ParVectorCreate(comm, jx_hpCSRMatrixGlobalNumRows(hp_matrix), 
                                  jx_hpCSRMatrixRowStarts(hp_matrix));
    jx_ParVectorInitialize(DBf);
    jx_ParVectorSetPartitioningOwner(DBf, 0);
    jx_BlockPrecDataDBf(bp_data) = DBf;

    Vtemp = jx_BlockPrecDataVtemp(bp_data);
    if (Vtemp != NULL)
    {
        jx_ParVectorDestroy(Vtemp);
        Vtemp = NULL;
    }
    Vtemp = jx_ParVectorCreate(comm, jx_hpCSRMatrixGlobalNumRows(hp_matrix), 
                                     jx_hpCSRMatrixRowStarts(hp_matrix));
    jx_ParVectorInitialize(Vtemp);
    jx_ParVectorSetPartitioningOwner(Vtemp, 0);
    jx_BlockPrecDataVtemp(bp_data) = Vtemp;

    Utemp = jx_BlockPrecDataUtemp(bp_data);
    if (Utemp != NULL)
    {
        jx_ParVectorDestroy(Utemp);
        Utemp = NULL;
    }
    Utemp = jx_ParVectorCreate(comm, jx_hpCSRMatrixGlobalNumRows(hp_matrix), 
                                     jx_hpCSRMatrixRowStarts(hp_matrix));
    jx_ParVectorInitialize(Utemp);
    jx_ParVectorSetPartitioningOwner(Utemp, 0);
    jx_BlockPrecDataUtemp(bp_data) = Utemp;

    return 0;
}                   

/*!
 * \fn JX_Int JX_BlockPrec_JAC
 * \brief Solver phase of Block preconditioner. (using Jacobi as inner iteration)
 * \date 2021/12/10
 */
JX_Int 
JX_BlockPrec_JAC(JX_Solver       solver,
                 JX_hpCSRMatrix  hp_matrix,
                 JX_ParVector    par_rhs,
                 JX_ParVector    par_app  )
{
    return (jx_BlockPrec_JAC((void *) solver, 
                            (jx_hpCSRMatrix *) hp_matrix,
                            (jx_ParVector *) par_rhs,
                            (jx_ParVector *) par_app));
}

/*!
 * \fn JX_Int JX_BlockPrec_GS
 * \brief Solver phase of Block preconditioner. (using GS as inner iteration -- core level)
 * \date 2022/5/18
 */
JX_Int 
JX_BlockPrec_GS(JX_Solver       solver,
                JX_hpCSRMatrix  hp_matrix,
                JX_ParVector    par_rhs,
                JX_ParVector    par_app  )
{
    return (jx_BlockPrec_GS((void *) solver, 
                            (jx_hpCSRMatrix *) hp_matrix,
                            (jx_ParVector *) par_rhs,
                            (jx_ParVector *) par_app));
}

/*!
 * \fn JX_Int JX_BlockPrec_SGS
 * \brief Solver phase of Block preconditioner. (using SGS as inner iteration -- core level)
 * \date 2022/5/28
 */
JX_Int 
JX_BlockPrec_SGS( JX_Solver       solver,
                     JX_hpCSRMatrix  hp_matrix,
                     JX_ParVector    par_rhs,
                     JX_ParVector    par_app  )
{
    return (jx_BlockPrec_SGS((void *) solver, 
                                (jx_hpCSRMatrix *) hp_matrix,
                                (jx_ParVector *) par_rhs,
                                (jx_ParVector *) par_app));
}

/*!
 * \fn JX_Int jx_BlockPrec_JAC
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix 
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2021/12/10
 */

JX_Int
jx_BlockPrec_JAC(void            *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  )
{
    jx_BlockPrecData *bp_data    = solver;
    JX_Int            bp_maxiter = jx_BlockPrecDataMaxIter(bp_data);
    JX_Int            bp_inviter = jx_BlockPrecDataInvIter(bp_data);
    JX_Int            bp_level   = jx_BlockPrecDataLevel(bp_data);

    jx_CSRMatrix    *A_diag      = jx_hpCSRMatrixDiag(hp_matrix); 
    JX_Int          *A_diag_i    = jx_CSRMatrixI(A_diag);
    JX_Real         *A_diag_data = jx_CSRMatrixData(A_diag);

    jx_ParCSRMatrix *A_B         = jx_hpMatrixLevelToPar(hp_matrix, bp_level);

    jx_ParCSRMatrix *A           = jx_hpCSRMatrixPar(hp_matrix);

    jx_Vector       *f_local = jx_ParVectorLocalVector(par_rhs);
    JX_Real         *f_data  = jx_VectorData(f_local);

    jx_ParVector    *Df      = NULL;
    jx_ParVector    *DBf     = NULL;
    jx_ParVector    *Vtemp   = NULL;
    jx_ParVector    *Utemp   = NULL;

    Df    = jx_BlockPrecDataDf(bp_data);
    DBf   = jx_BlockPrecDataDBf(bp_data);
    Vtemp = jx_BlockPrecDataVtemp(bp_data);
    Utemp = jx_BlockPrecDataUtemp(bp_data);

    jx_Vector      *Df_local = jx_ParVectorLocalVector(Df);
    JX_Real        *Df_data  = jx_VectorData(Df_local);
    jx_Vector      *DBf_local = jx_ParVectorLocalVector(DBf);

    jx_Vector      *Vtemp_local = jx_ParVectorLocalVector(Vtemp);
    JX_Real        *Vtemp_data  = jx_VectorData(Vtemp_local);
    jx_Vector      *Utemp_local = jx_ParVectorLocalVector(Utemp);
    JX_Real        *Utemp_data  = jx_VectorData(Utemp_local);

    JX_Int          i, it, invit, ierr = 0;
    JX_Int          local_size  = jx_VectorSize(f_local);

    #if 1
        /* DBf = D_B^{-1}*f */
        for(i = 0; i < local_size; i++)
        {   
            /* Df = D^{-1}*f */
            Df_data[i] = f_data[i] / A_diag_data[A_diag_i[i]];  
        }
        jx_ParVectorCopy(Df, DBf);
        for(invit = 1; invit < bp_inviter; invit++)
        {   
            jx_ParCSRMatrixMatvec(1.0, A_B, DBf, 0.0, Vtemp);
            for(i = 0; i < local_size; i++)
            {
                Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
            }

            jx_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
            jx_SeqVectorAxpy(1.0, Df_local, DBf_local);
        }
        
        jx_ParVectorCopy(DBf, par_app);
        for(it = 1; it < bp_maxiter; it++) {
            /* temp_1 = A * u */
            jx_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

            /* temp_1 = D_B^{-1} * A*u */
            for(i = 0; i < local_size; i++)
            {
                Df_data[i] = Utemp_data[i] / A_diag_data[A_diag_i[i]];  /* Df = D^{-1} * A*u  */
            }
            jx_SeqVectorCopy(Df_local, Utemp_local);
            for(invit = 1; invit < bp_inviter; invit++)
            {   
                /* Vtemp = A_B * Utemp */
                jx_ParCSRMatrixMatvec(1.0, A_B, Utemp, 0.0, Vtemp);

                /* Vtemp = D^{-1} * Vtemp */
                for(i = 0; i < local_size; i++)
                {
                    Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
                }

                /* Utemp = -Vtemp + Utemp + Df */
                jx_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
                jx_SeqVectorAxpy(1.0, Df_local, Utemp_local);
            }
            /* u = -Utemp + u */
            jx_ParVectorAxpy(-1.0, Utemp, par_app);

            /* u = DBf + u */
            jx_ParVectorAxpy(1.0, DBf, par_app);
        }
    #endif

    jx_TFree(A_B);
    
    //BlockPrec//
    #if 0
        /* Df = D^{-1}*f */
        for(i = 0; i < local_size; i++)
        {
            Df_data[i] = f_data[i] / A_diag_data[A_diag_i[i]];
        }

        /* temp = Df */
        ierr = jx_ParVectorCopy(Df, par_app);

        for(it = 1; it < bp_maxiter; it++)
        {   
            /* temp = u */
            jx_ParVectorCopy(par_app, Utemp);

            /* temp = A*u */
            ierr = jx_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);

            /* temp = D^{-1}*temp */
            for(i = 0; i < local_size; i++)
            {
                Vtemp_data[i] = Vtemp_data[i] / A_diag_data[A_diag_i[i]];
            }

            /* u = -temp + u */
            ierr = jx_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);

            /* u = Df + u */
            ierr = jx_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }
    #endif 

    return ierr;
}

/*!
 * \fn JX_Int jx_BlockPrec_GS
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix 
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2022/5/18
 */
JX_Int
jx_BlockPrec_GS( void            *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  )
{
    jx_BlockPrecData *bp_data    = solver;
    JX_Int            bp_maxiter = jx_BlockPrecDataMaxIter(bp_data);
    JX_Int            bp_inviter = jx_BlockPrecDataInvIter(bp_data);

    jx_CSRMatrix    *A_diag      = jx_hpCSRMatrixDiag(hp_matrix); 
    JX_Int          *A_diag_i    = jx_CSRMatrixI(A_diag);
    JX_Int          *A_diag_j    = jx_CSRMatrixJ(A_diag);
    JX_Real         *A_diag_data = jx_CSRMatrixData(A_diag);

    JX_Int           n           = jx_CSRMatrixNumRows(A_diag);

    jx_ParCSRMatrix *A           = jx_hpCSRMatrixPar(hp_matrix);

    jx_Vector       *f_local = jx_ParVectorLocalVector(par_rhs);

    jx_ParVector    *Df      = NULL;
    jx_ParVector    *DBf     = NULL;
    jx_ParVector    *Vtemp   = NULL;
    jx_ParVector    *Utemp   = NULL;

    Df    = jx_BlockPrecDataDf(bp_data);
    DBf   = jx_BlockPrecDataDBf(bp_data);
    Vtemp = jx_BlockPrecDataVtemp(bp_data);
    Utemp = jx_BlockPrecDataUtemp(bp_data);

    jx_Vector      *Df_local  = jx_ParVectorLocalVector(Df);
    JX_Real        *Df_data   = jx_VectorData(Df_local);
    jx_Vector      *DBf_local = jx_ParVectorLocalVector(DBf);

    jx_Vector      *Vtemp_local = jx_ParVectorLocalVector(Vtemp);
    JX_Real        *Vtemp_data  = jx_VectorData(Vtemp_local);
    jx_Vector      *Utemp_local = jx_ParVectorLocalVector(Utemp);
    JX_Real        *Utemp_data  = jx_VectorData(Utemp_local);

    JX_Int          i, it, invit, ierr = 0;
    JX_Int          ii, jj;

    JX_Real         zero = 0.0;
    JX_Real         res;

    /* Df = (D+L)^{-1}*f */
    jx_SeqVectorCopy(f_local, Df_local);
    for (i = 0; i < n; i++)	
    {
        if ( A_diag_data[A_diag_i[i]] != zero )
        {
            res = Df_data[i];
            for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
    jx_SeqVectorCopy(Df_local, DBf_local);
    for(invit = 1; invit < bp_inviter; invit++)
    {   
        /*  Vtemp = A*DBf */
        jx_CSRMatrixMatvec_origin(1.0, A_diag, DBf_local, 0.0, Vtemp_local);
        /*  Vtemp = (D+L)^{-1} * Vtemp */
        for (i = 0; i < n; i++)	
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {
                res = Vtemp_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
        jx_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
        jx_SeqVectorAxpy(1.0, Df_local, DBf_local);
    }
    
    jx_ParVectorCopy(DBf, par_app);
    for(it = 1; it < bp_maxiter; it++) 
    {
        /* Utemp = A * u */
        jx_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

        /* Utemp = D_B^{-1} * A * u */
        for (i = 0; i < n; i++)	
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {
                res = Utemp_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
        jx_SeqVectorCopy(Df_local, Utemp_local);
        for(invit = 1; invit < bp_inviter; invit++)
        {   
            /*  Vtemp = A * Utemp */
            jx_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);
            /*  Vtemp = (D+L)^{-1} * Vtemp */
            for (i = 0; i < n; i++)	
            {
                if ( A_diag_data[A_diag_i[i]] != zero )
                {
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
            jx_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
            jx_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }
        
        /* u = -Utemp + u */
        jx_ParVectorAxpy(-1.0, Utemp, par_app);

        /* u = DBf + u */
        jx_ParVectorAxpy(1.0, DBf, par_app);
    }

    return ierr;
}


/*!
 * \fn JX_Int jx_BlockPrec_SGS
 * \brief BLock preconditioner.
 * \param solver pointer to NULL
 * \param hp_matrix pointer to the coefficient matrix 
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2022/6/8
 */
JX_Int
jx_BlockPrec_SGS(void            *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  )
{
    jx_BlockPrecData *bp_data    = solver;
    JX_Int            bp_maxiter = jx_BlockPrecDataMaxIter(bp_data);
    JX_Int            bp_inviter = jx_BlockPrecDataInvIter(bp_data);

    jx_CSRMatrix    *A_diag      = jx_hpCSRMatrixDiag(hp_matrix); 
    JX_Int          *A_diag_i    = jx_CSRMatrixI(A_diag);
    JX_Int          *A_diag_j    = jx_CSRMatrixJ(A_diag);
    JX_Real         *A_diag_data = jx_CSRMatrixData(A_diag);

    JX_Int           n           = jx_CSRMatrixNumRows(A_diag);

    jx_ParCSRMatrix *A           = jx_hpCSRMatrixPar(hp_matrix);

    jx_Vector       *f_local = jx_ParVectorLocalVector(par_rhs);

    jx_ParVector    *Df      = NULL;
    jx_ParVector    *DBf     = NULL;
    jx_ParVector    *Vtemp   = NULL;
    jx_ParVector    *Utemp   = NULL;

    Df    = jx_BlockPrecDataDf(bp_data);
    DBf   = jx_BlockPrecDataDBf(bp_data);
    Vtemp = jx_BlockPrecDataVtemp(bp_data);
    Utemp = jx_BlockPrecDataUtemp(bp_data);

    jx_Vector      *Df_local = jx_ParVectorLocalVector(Df);
    JX_Real        *Df_data  = jx_VectorData(Df_local);
    jx_Vector      *DBf_local = jx_ParVectorLocalVector(DBf);

    jx_Vector      *Vtemp_local = jx_ParVectorLocalVector(Vtemp);
    JX_Real        *Vtemp_data  = jx_VectorData(Vtemp_local);
    jx_Vector      *Utemp_local = jx_ParVectorLocalVector(Utemp);
    JX_Real        *Utemp_data  = jx_VectorData(Utemp_local);

    JX_Int          i, it, invit, ierr = 0;
    JX_Int          ii, jj;

    JX_Real         zero = 0.0;
    JX_Real         res;

    /* Df = (D+U)^{-1}*D*(D+L)^{-1}*f */
    /* Df = (D+L)^{-1}*f  */
    jx_SeqVectorCopy(f_local, Df_local);
    for (i = 0; i < n; i++)	
    {
        if ( A_diag_data[A_diag_i[i]] != zero )
        {
            res = Df_data[i];
            for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
        if ( A_diag_data[A_diag_i[i]] != zero )
        {   
            Df_data[i] *= A_diag_data[A_diag_i[i]];
        }
    }
    /*  Df = (D+U)^{-1} * Df  */
    for (i = n-1; i >= 0; i--)
    {
        if ( A_diag_data[A_diag_i[i]] != zero )
        {   
            res = Df_data[i];
            for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
    jx_SeqVectorCopy(Df_local, DBf_local);

    for (invit = 1; invit < bp_inviter; invit++)
    {   
        /*  Vtemp = A*DBf */
        jx_CSRMatrixMatvec_origin(1.0, A_diag, DBf_local, 0.0, Vtemp_local);
        /*  Vtemp = (D+L)^{-1} * Vtemp */
        for (i = 0; i < n; i++)	
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {
                res = Vtemp_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
            if ( A_diag_data[A_diag_i[i]] != zero )
            {   
                Vtemp_data[i] *= A_diag_data[A_diag_i[i]];
            }
        }
        /*  Vtemp = (D+U)^{-1} * Vtemp */
        for (i = n-1; i >= 0; i--)
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {   
                res = Vtemp_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
        jx_SeqVectorAxpy(-1.0, Vtemp_local, DBf_local);
        jx_SeqVectorAxpy(1.0, Df_local, DBf_local);
    }
    
    jx_ParVectorCopy(DBf, par_app);
    for(it = 1; it < bp_maxiter; it++) 
    {
        /* Utemp = A * u */
        jx_ParCSRMatrixMatvec(1.0, A, par_app, 0.0, Utemp);

        /* Utemp = D_B^{-1} * A * u */
        for (i = 0; i < n; i++)	
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {
                res = Utemp_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
            if ( A_diag_data[A_diag_i[i]] != zero )
            {   
                Df_data[i] *= A_diag_data[A_diag_i[i]];
            }
        }
        for (i = n-1; i >= 0; i--)
        {
            if ( A_diag_data[A_diag_i[i]] != zero )
            {   
                res = Df_data[i];
                for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
        jx_SeqVectorCopy(Df_local, Utemp_local);
        for(invit = 1; invit < bp_inviter; invit++)
        {   
            /*  Vtemp = A*Utemp + Vtemp */
            jx_CSRMatrixMatvec_origin(1.0, A_diag, Utemp_local, 0.0, Vtemp_local);
            /*  Vtemp = (D_B)^{-1} * Vtemp */
            for (i = 0; i < n; i++)	
            {
                if ( A_diag_data[A_diag_i[i]] != zero )
                {
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
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
                if ( A_diag_data[A_diag_i[i]] != zero )
                {   
                    Vtemp_data[i] *= A_diag_data[A_diag_i[i]];
                }
            }
            for (i = n-1; i >= 0; i--)
            {
                if ( A_diag_data[A_diag_i[i]] != zero )
                {   
                    res = Vtemp_data[i];
                    for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj++)
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
            jx_SeqVectorAxpy(-1.0, Vtemp_local, Utemp_local);
            jx_SeqVectorAxpy(1.0, Df_local, Utemp_local);
        }
        
        /* u = -Utemp + u */
        jx_ParVectorAxpy(-1.0, Utemp, par_app);

        /* u = DBf + u */
        jx_ParVectorAxpy(1.0, DBf, par_app);
    }

    return ierr;
}
