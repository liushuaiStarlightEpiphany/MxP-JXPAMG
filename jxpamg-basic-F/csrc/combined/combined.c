//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  combined.c
 *  Date: 2013/12/11
 *
 *  Created by Yue Xiaoqiang
 *
 */

#include "jxf_combined.h"

#define JXF_EUCLID_ERRCHKB \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("", __FUNC__, __FILE__, __LINE__); \
        jxf_printErrorMsg(stderr); \
        jxf_MPI_Abort(jxf_comm_dh, -1); \
      }

#undef ENABLE_EUCLID_LOGGING

#if !defined(ENABLE_EUCLID_LOGGING)
#undef JXF_START_FUNC_DH
#undef JXF_END_FUNC_VAL
#undef JXF_END_FUNC_DH
#define JXF_START_FUNC_DH
#define JXF_END_FUNC_DH
#define JXF_END_FUNC_VAL(a) return(a);
#endif

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataCreate
 */
JXF_Int
JXF_CombinedPrecondDataCreate( JXF_Solver *solver, MPI_Comm comm )
{
   *solver = (JXF_Solver) jxf_CombinedPrecondDataCreate(comm);
    if (!solver)
    {
        jxf_error_in_arg(1);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn void *jxf_CombinedPrecondDataCreate
 */
void *
jxf_CombinedPrecondDataCreate( MPI_Comm comm )
{
    jxf_CombinedPrecondData *pre_data;
    JXF_Int pre_id = 1;
    
    pre_data = jxf_CTAlloc(jxf_CombinedPrecondData, 1);
    jxf_CombinedPrecondDataComm(pre_data) = comm;
    jxf_CombinedPrecondDataSetPreID(pre_data, pre_id);
    
    return (void *)pre_data;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataInitialize
 */
JXF_Int
JXF_CombinedPrecondDataInitialize( JXF_Solver solver, JXF_Int level )
{
    return( jxf_CombinedPrecondDataInitialize( (void *) solver, level ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataInitialize
 */
JXF_Int
jxf_CombinedPrecondDataInitialize( void *data, JXF_Int level )
{
    if (!data)
    {
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    jxf_CSRMatrix *ser_mat = NULL;
    JXF_Int nprocs, myid, stg_dd;
    
    jxf_MPI_Comm_size(comm, &nprocs);
    jxf_MPI_Comm_rank(comm, &myid);
    if (pre_id == 11) // PAMG-Euclid-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
    }
    else if (pre_id == 12) // Euclid-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
    }
    else if (pre_id == 13) // PAMG-Euclid
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
    }
    else if (pre_id == 14) // Euclid-PAMG-Euclid
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
    }
    else if (pre_id == 15) // Adaptive
    {
        if (nprocs == 1)
        {
            if (myid == 0)
            {
                stg_dd = jxf_CSRMatrixWeaklyDiagDominant(
                                    jxf_ParCSRMatrixDiag(jxf_CombinedPrecondDataPreMat(pre_data)),
                                        jxf_CombinedPrecondDataNumEquns(pre_data),
                                            jxf_CombinedPrecondDataAdpTheta(pre_data),
                                                jxf_CombinedPrecondDataAdpGammaT(pre_data),
                                                    jxf_CombinedPrecondDataAdpGammaE(pre_data));
                if (stg_dd == 1)
                {
                    JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
                    JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
                    JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                    JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
                    jxf_CombinedPrecondDataPreID(pre_data) = 13; // switch to PAMG-Euclid
                }
                else
                {
                    JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                    JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
                }
            }
        }
        else
        {
            ser_mat = jxf_ParCSRMatrixToCSRMatrixAll(jxf_CombinedPrecondDataPreMat(pre_data));
            if (myid == 0)
            {
                stg_dd = jxf_CSRMatrixWeaklyDiagDominant(ser_mat,
                                        jxf_CombinedPrecondDataNumEquns(pre_data),
                                            jxf_CombinedPrecondDataAdpTheta(pre_data),
                                                jxf_CombinedPrecondDataAdpGammaT(pre_data),
                                                    jxf_CombinedPrecondDataAdpGammaE(pre_data));
            }
            jxf_CSRMatrixDestroy(ser_mat);
            jxf_MPI_Bcast(&stg_dd, 1, JXF_MPI_INT, 0, comm);
            if (stg_dd == 1)
            {
                JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
                JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
                JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
                jxf_CombinedPrecondDataPreID(pre_data) = 13; // switch to PAMG-Euclid
            }
            else
            {
                JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
            }
        }
    }
    else if (pre_id == 16) // ILU-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
    }
    else if (pre_id == 17) // ILU-PAMG-ILU
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
    }
    else if (pre_id == 18) // PAMG-Euclid-PAMG, where Euclid for (A+A^T)/2
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidLevel(jxf_CombinedPrecondDataEuclidData(pre_data), level);
    }
    else if (pre_id == 19) // PAMG-ILU-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
    }
    else
    {
        jxf_printf("\n Warning: Wrong preconditioner ID\n\n");
        exit(-3);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetEuclidLevel
 */
#undef __FUNC__
#define __FUNC__ "JXF_CombinedPrecondDataSetEuclidLevel"
JXF_Int
JXF_CombinedPrecondDataSetEuclidLevel( JXF_Solver solver, JXF_Int level )
{
    char str_level[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_level,"%d",level);
    jxf_Parser_dhInsert(jxf_parser_dh, "-level", str_level); JXF_EUCLID_ERRCHKB;
    JXF_END_FUNC_VAL(0)
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataInitializeP
 */
JXF_Int
JXF_CombinedPrecondDataInitializeP( JXF_Solver solver,
                                   JXF_Int argc,
                                   char *argv[] )
{
    return( jxf_CombinedPrecondDataInitializeP( (void *) solver, argc, argv ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataInitializeP
 */
JXF_Int
jxf_CombinedPrecondDataInitializeP( void *data,
                                   JXF_Int argc,
                                   char *argv[] )
{
    if (!data)
    {
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    jxf_CSRMatrix *ser_mat = NULL;
    JXF_Int nprocs, myid, stg_dd;
    
    jxf_MPI_Comm_size(comm, &nprocs);
    jxf_MPI_Comm_rank(comm, &myid);
    if (pre_id == 11) // PAMG-Euclid-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
    }
    else if (pre_id == 12) // Euclid-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
    }
    else if (pre_id == 13) // PAMG-Euclid
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
    }
    else if (pre_id == 14) // Euclid-PAMG-Euclid
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
    }
    else if (pre_id == 15) // Adaptive
    {
        if (nprocs == 1)
        {
            if (myid == 0)
            {
                stg_dd = jxf_CSRMatrixWeaklyDiagDominant(
                                    jxf_ParCSRMatrixDiag(jxf_CombinedPrecondDataPreMat(pre_data)),
                                        jxf_CombinedPrecondDataNumEquns(pre_data),
                                            jxf_CombinedPrecondDataAdpTheta(pre_data),
                                                jxf_CombinedPrecondDataAdpGammaT(pre_data),
                                                    jxf_CombinedPrecondDataAdpGammaE(pre_data));
                if (stg_dd == 1)
                {
                    JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
                    JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
                    JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                    JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
                    jxf_CombinedPrecondDataPreID(pre_data) = 13; // switch to PAMG-Euclid
                }
                else
                {
                    JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                    JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
                }
            }
        }
        else
        {
            ser_mat = jxf_ParCSRMatrixToCSRMatrixAll(jxf_CombinedPrecondDataPreMat(pre_data));
            if (myid == 0)
            {
                stg_dd = jxf_CSRMatrixWeaklyDiagDominant(ser_mat,
                                        jxf_CombinedPrecondDataNumEquns(pre_data),
                                            jxf_CombinedPrecondDataAdpTheta(pre_data),
                                                jxf_CombinedPrecondDataAdpGammaT(pre_data),
                                                    jxf_CombinedPrecondDataAdpGammaE(pre_data));
            }
            jxf_CSRMatrixDestroy(ser_mat);
            jxf_MPI_Bcast(&stg_dd, 1, JXF_MPI_INT, 0, comm);
            if (stg_dd == 1)
            {
                JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
                JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
                JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
                jxf_CombinedPrecondDataPreID(pre_data) = 13; // switch to PAMG-Euclid
            }
            else
            {
                JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
            }
        }
    }
    else if (pre_id == 16) // ILU-PAMG
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
    }
    else if (pre_id == 17) // ILU-PAMG-ILU
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
    }
    else if (pre_id == 18) // PAMG-Euclid-PAMG, where Euclid for (A+A^T)/2
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetAIMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), 0);
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_CombinedPrecondDataSetEuclidParams(jxf_CombinedPrecondDataEuclidData(pre_data), argc, argv);
    }
    else
    {
        jxf_printf("\n Warning: Wrong preconditioner ID\n\n");
        exit(-3);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetEuclidParams
 */
#undef __FUNC__
#define __FUNC__ "JXF_CombinedPrecondDataSetEuclidParams"
JXF_Int
JXF_CombinedPrecondDataSetEuclidParams( JXF_Solver solver,
                                       JXF_Int argc,
                                       char *argv[] )
{
    jxf_Parser_dhInit(jxf_parser_dh, argc, argv);
    JXF_CP_EUCLID_ERRCHKA;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetPreID
 */
JXF_Int
JXF_CombinedPrecondDataSetPreID( JXF_Solver solver, JXF_Int pre_id )
{
    return( jxf_CombinedPrecondDataSetPreID( (void *) solver, pre_id ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetPreID
 */
JXF_Int
jxf_CombinedPrecondDataSetPreID( void *data, JXF_Int pre_id )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataPreID(pre_data) = pre_id;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetNxyNpxyNequ
 */
JXF_Int
JXF_CombinedPrecondDataSetNxyNpxyNequ( JXF_Solver solver, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy, JXF_Int num_equns )
{
    return( jxf_CombinedPrecondDataSetNxyNpxyNequ( (void *) solver, nx, ny, npx, npy, num_equns ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetNxyNpxyNequ
 */
JXF_Int
jxf_CombinedPrecondDataSetNxyNpxyNequ( void *data, JXF_Int nx, JXF_Int ny, JXF_Int npx, JXF_Int npy, JXF_Int num_equns )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataNx(pre_data) = nx;
    jxf_CombinedPrecondDataNy(pre_data) = ny;
    jxf_CombinedPrecondDataNpx(pre_data) = npx;
    jxf_CombinedPrecondDataNpy(pre_data) = npy;
    jxf_CombinedPrecondDataNumEquns(pre_data) = num_equns;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetThetaPsiRhoPhi
 */
JXF_Int
JXF_CombinedPrecondDataSetThetaPsiRhoPhi( JXF_Solver solver,
                                         JXF_Int theta_psi,
                                         JXF_Int theta_rho,
                                         JXF_Int theta_phi,
                                         JXF_Real theta_dis )
{
    return( jxf_CombinedPrecondDataSetThetaPsiRhoPhi( (void *) solver, theta_psi, theta_rho, theta_phi, theta_dis ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetThetaPsiRhoPhi
 */
JXF_Int
jxf_CombinedPrecondDataSetThetaPsiRhoPhi( void *data,
                                         JXF_Int theta_psi,
                                         JXF_Int theta_rho,
                                         JXF_Int theta_phi,
                                         JXF_Real theta_dis )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataThetaPsi(pre_data) = theta_psi;
    jxf_CombinedPrecondDataThetaRho(pre_data) = theta_rho;
    jxf_CombinedPrecondDataThetaPhi(pre_data) = theta_phi;
    jxf_CombinedPrecondDataThetaDis(pre_data) = theta_dis;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetAMGParameters
 */
JXF_Int
JXF_CombinedPrecondDataSetAMGParameters( JXF_Solver solver,
                                        JXF_Int amg_max_levels,
                                        JXF_Int amg_relax_type,
                                        JXF_Int amg_print_level,
                                        JXF_Int amg_interp_type,
                                        JXF_Int amg_P_max_elmts,
                                        JXF_Int amg_measure_type,
                                        JXF_Int amg_coarsen_type,
                                        JXF_Int amg_agg_num_levels,
                                        JXF_Int amg_coarse_threshold,
                                        JXF_Real amg_strong_threshold )
{
    return( jxf_CombinedPrecondDataSetAMGParameters( (void *) solver, amg_max_levels, amg_relax_type,
                                                    amg_print_level, amg_interp_type, amg_P_max_elmts,
                                                    amg_measure_type, amg_coarsen_type, amg_agg_num_levels,
                                                    amg_coarse_threshold, amg_strong_threshold ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetAMGParameters
 */
JXF_Int
jxf_CombinedPrecondDataSetAMGParameters( void *data,
                                        JXF_Int amg_max_levels,
                                        JXF_Int amg_relax_type,
                                        JXF_Int amg_print_level,
                                        JXF_Int amg_interp_type,
                                        JXF_Int amg_P_max_elmts,
                                        JXF_Int amg_measure_type,
                                        JXF_Int amg_coarsen_type,
                                        JXF_Int amg_agg_num_levels,
                                        JXF_Int amg_coarse_threshold,
                                        JXF_Real amg_strong_threshold )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAMGMaxLevels(pre_data) = amg_max_levels;
    jxf_CombinedPrecondDataAMGRelaxType(pre_data) = amg_relax_type;
    jxf_CombinedPrecondDataAMGPrintLevel(pre_data) = amg_print_level;
    jxf_CombinedPrecondDataAMGInterpType(pre_data) = amg_interp_type;
    jxf_CombinedPrecondDataAMGPMaxElmts(pre_data) = amg_P_max_elmts;
    jxf_CombinedPrecondDataAMGMeasureType(pre_data) = amg_measure_type;
    jxf_CombinedPrecondDataAMGCoarsenType(pre_data) = amg_coarsen_type;
    jxf_CombinedPrecondDataAMGAggNumLevels(pre_data) = amg_agg_num_levels;
    jxf_CombinedPrecondDataAMGCoarseThreshold(pre_data) = amg_coarse_threshold;
    jxf_CombinedPrecondDataAMGStrongThreshold(pre_data) = amg_strong_threshold;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetPreMat
 */
JXF_Int
JXF_CombinedPrecondDataSetPreMat( JXF_Solver solver, jxf_ParCSRMatrix *pre_mat )
{
    return( jxf_CombinedPrecondDataSetPreMat( (void *) solver, pre_mat ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetPreMat
 */
JXF_Int
jxf_CombinedPrecondDataSetPreMat( void *data, jxf_ParCSRMatrix *pre_mat )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataPreMat(pre_data) = pre_mat;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetDropTol
 */
JXF_Int
JXF_CombinedPrecondDataSetDropTol( JXF_Solver solver, JXF_Real drop_tol )
{
    return( jxf_CombinedPrecondDataSetDropTol( (void *) solver, drop_tol ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetDropTol
 */
JXF_Int
jxf_CombinedPrecondDataSetDropTol( void *data, JXF_Real drop_tol )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataDropTol(pre_data) = drop_tol;
    
    return jxf_error_flag;
}

JXF_Int
JXF_CombinedPrecondDataSetInterpType( JXF_Solver solver, JXF_Int interp_type )
{
    return( jxf_CombinedPrecondDataSetInterpType( (void *) solver, interp_type ) );
}

JXF_Int
jxf_CombinedPrecondDataSetInterpType( void *data, JXF_Int interp_type )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAMGInterpType(pre_data) = interp_type;
    
    return jxf_error_flag;
}

JXF_Int
JXF_CombinedPrecondDataSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type )
{
    return( jxf_CombinedPrecondDataSetCoarsenType( (void *) solver, coarsen_type ) );
}

JXF_Int
jxf_CombinedPrecondDataSetCoarsenType( void *data, JXF_Int coarsen_type )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAMGCoarsenType(pre_data) = coarsen_type;
    
    return jxf_error_flag;
}

JXF_Int
JXF_CombinedPrecondDataSetCycleRelaxType( JXF_Solver solver, JXF_Int relax_type )
{
    return( jxf_CombinedPrecondDataSetCycleRelaxType( (void *) solver, relax_type ) );
}

JXF_Int
jxf_CombinedPrecondDataSetCycleRelaxType( void *data, JXF_Int relax_type )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAMGRelaxType(pre_data) = relax_type;
    
    return jxf_error_flag;
}

JXF_Int
JXF_CombinedPrecondDataSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold )
{
    return( jxf_CombinedPrecondDataSetStrongThreshold( (void *) solver, strong_threshold ) );
}

JXF_Int
jxf_CombinedPrecondDataSetStrongThreshold( void *data, JXF_Real strong_threshold )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAMGStrongThreshold(pre_data) = strong_threshold;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetAdpTheta
 */
JXF_Int
JXF_CombinedPrecondDataSetAdpTheta( JXF_Solver solver, JXF_Real adp_theta )
{
    return( jxf_CombinedPrecondDataSetAdpTheta( (void *) solver, adp_theta ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetAdpTheta
 */
JXF_Int
jxf_CombinedPrecondDataSetAdpTheta( void *data, JXF_Real adp_theta )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAdpTheta(pre_data) = adp_theta;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetAdpGammaT
 */
JXF_Int
JXF_CombinedPrecondDataSetAdpGammaT( JXF_Solver solver, JXF_Real adp_gamma_3 )
{
    return( jxf_CombinedPrecondDataSetAdpGammaT( (void *) solver, adp_gamma_3 ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetAdpGammaT
 */
JXF_Int
jxf_CombinedPrecondDataSetAdpGammaT( void *data, JXF_Real adp_gamma_3 )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAdpGammaT(pre_data) = adp_gamma_3;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetAdpGammaE
 */
JXF_Int
JXF_CombinedPrecondDataSetAdpGammaE( JXF_Solver solver, JXF_Real adp_gamma_11 )
{
    return( jxf_CombinedPrecondDataSetAdpGammaE( (void *) solver, adp_gamma_11 ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetAdpGammaE
 */
JXF_Int
jxf_CombinedPrecondDataSetAdpGammaE( void *data, JXF_Real adp_gamma_11 )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    jxf_CombinedPrecondDataAdpGammaE(pre_data) = adp_gamma_11;
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetILUMatA
 */
JXF_Int
JXF_CombinedPrecondDataSetILUMatA( JXF_Solver solver, jxf_CSRMatrix *matA )
{
    return( jxf_CombinedPrecondDataSetILUMatA( (void *) solver, matA ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetILUMatA
 */
JXF_Int
jxf_CombinedPrecondDataSetILUMatA( void *data, jxf_CSRMatrix *matA )
{
    jxf_CombinedPrecondData *pre_data = data;
    if (!pre_data)
    {
        jxf_printf(" Warning: CombinedPrecond object empty!\n");
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
    JXF_ILUZeroFactorDataSetMatA(jxf_CombinedPrecondDataILUZData(pre_data), matA);
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetPreEucMat
 */
JXF_Int
JXF_CombinedPrecondDataSetPreEucMat( JXF_Solver solver, jxf_ParCSRMatrix *matA )
{
    return( jxf_CombinedPrecondDataSetPreEucMat( (void *) solver, matA ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetPreEucMat
 */
JXF_Int
jxf_CombinedPrecondDataSetPreEucMat( void *data, jxf_ParCSRMatrix *matA )
{
    jxf_CombinedPrecondData *pre_data = data;
    if (!pre_data)
    {
        jxf_printf(" Warning: CombinedPrecond object empty!\n");
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
    jxf_CombinedPrecondDataPreEucMat(pre_data) = matA;
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataGetLULength
 */
JXF_Int
JXF_CombinedPrecondDataGetLULength( JXF_Solver solver, JXF_Int *lu_length )
{
    return( jxf_CombinedPrecondDataGetLULength( (void *) solver, lu_length ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataGetLULength
 */
JXF_Int
jxf_CombinedPrecondDataGetLULength( void *data, JXF_Int *lu_length )
{
    jxf_CombinedPrecondData *pre_data = data;
    if (!pre_data)
    {
        jxf_printf(" Warning: CombinedPrecond object empty!\n");
        jxf_error_in_arg(1);
        return jxf_error_flag;
    }
    JXF_ILUZeroFactorDataGetLULength(jxf_CombinedPrecondDataILUZData(pre_data), lu_length);
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSetup
 */
JXF_Int
JXF_CombinedPrecondDataSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_CombinedPrecondDataSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSetup
 */
JXF_Int
jxf_CombinedPrecondDataSetup( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    
    if (pre_id == 11) // PAMG-Euclid-PAMG
    {
        jxf_PAMGEuclidPAMGSetup(data, par_matrix);
    }
    else if (pre_id == 12) // Euclid-PAMG
    {
        jxf_EuclidPAMGSetup(data, par_matrix);
    }
    else if (pre_id == 13) // PAMG-Euclid
    {
        jxf_PAMGEuclidSetup(data, par_matrix);
    }
    else if (pre_id == 14) // Euclid-PAMG-Euclid
    {
        jxf_EuclidPAMGEuclidSetup(data, par_matrix);
    }
    else if (pre_id == 15) // Euclid
    {
        JXF_EuclidSetup(jxf_CombinedPrecondDataEuclidData(pre_data), (JXF_hpCSRMatrix)par_matrix);
    }
    else if (pre_id == 16) // ILU-PAMG
    {
        jxf_ILUZPAMGSetup(data, par_matrix);
    }
    else if (pre_id == 17) // ILU-PAMG-ILU
    {
        jxf_ILUZPAMGILUZSetup(data, par_matrix);
    }
    else if (pre_id == 18) // PAMG-Euclid-PAMG, where Euclid for (A+A^T)/2
    {
        jxf_PAMGEuclidPAMGSetupB(data, par_matrix);
    }
    else if (pre_id == 19) // PAMG-ILU-PAMG
    {
        jxf_PAMGILUZPAMGSetup(data, par_matrix);
    }
    else
    {
        jxf_printf("\n Warning: Wrong preconditioner ID\n\n");
        exit(-3);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataAdaptiveSetup2
 */
JXF_Int
JXF_CombinedPrecondDataAdaptiveSetup2( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix )
{
    return( jxf_CombinedPrecondDataAdaptiveSetup2( (void *) solver, (jxf_hpCSRMatrix *) hp_matrix ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataAdaptiveSetup2
 */
JXF_Int
jxf_CombinedPrecondDataAdaptiveSetup2( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Int nx = jxf_CombinedPrecondDataNx(pre_data);
    JXF_Int ny = jxf_CombinedPrecondDataNy(pre_data);
    JXF_Int npx = jxf_CombinedPrecondDataNpx(pre_data);
    JXF_Int npy = jxf_CombinedPrecondDataNpy(pre_data);
    JXF_Int num_equns = jxf_CombinedPrecondDataNumEquns(pre_data);
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    JXF_Int euclid_level = jxf_CombinedPrecondDataEuclidLevel(pre_data);
    JXF_Real drop_tol = jxf_CombinedPrecondDataDropTol(pre_data);
    
    if (pre_id == 1) // Solve by Euclid firstly, Yue Xiaoqiang 2014/10/24
    {
        JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
        JXF_EuclidSetLevel(jxf_CombinedPrecondDataEuclidData(pre_data), euclid_level);
        JXF_EuclidSetSparseA(jxf_CombinedPrecondDataEuclidData(pre_data), drop_tol);
        JXF_EuclidSetup(jxf_CombinedPrecondDataEuclidData(pre_data), (JXF_hpCSRMatrix)par_matrix);
    }
    else if (pre_id == 4) // Solve by ILU(0) - crossed firstly, Yue Xiaoqiang 2014/12/03
    {
        JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
        JXF_ILUZeroFactorDataSetNxy(jxf_CombinedPrecondDataILUZData(pre_data), nx, ny);
        JXF_ILUZeroFactorDataSetNpxy(jxf_CombinedPrecondDataILUZData(pre_data), npx, npy);
        JXF_ILUZeroFactorDataSetNumEquns(jxf_CombinedPrecondDataILUZData(pre_data), num_equns);
        JXF_ILUZeroFactorDataSetDropTol(jxf_CombinedPrecondDataILUZData(pre_data), drop_tol);
        JXF_ILUZeroFactorDataGenerateParGrid(jxf_CombinedPrecondDataILUZData(pre_data));
        JXF_ILUZeroFactorDataSetup(jxf_CombinedPrecondDataILUZData(pre_data),(JXF_hpCSRMatrix)par_matrix);
    }  
      
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataAdaptiveSetup3
 */
JXF_Int
JXF_CombinedPrecondDataAdaptiveSetup3( JXF_Solver solver, JXF_hpCSRMatrix hp_matrix )
{
    return( jxf_CombinedPrecondDataAdaptiveSetup3( (void *) solver, (jxf_hpCSRMatrix *) hp_matrix) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataAdaptiveSetup3
 */
JXF_Int
jxf_CombinedPrecondDataAdaptiveSetup3( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Int nx = jxf_CombinedPrecondDataNx(pre_data);
    JXF_Int ny = jxf_CombinedPrecondDataNy(pre_data);
    JXF_Int npx = jxf_CombinedPrecondDataNpx(pre_data);
    JXF_Int npy = jxf_CombinedPrecondDataNpy(pre_data);
    JXF_Int num_equns = jxf_CombinedPrecondDataNumEquns(pre_data);
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    JXF_Int num_rows = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(par_matrix));
    JXF_Int theta_psi = jxf_CombinedPrecondDataThetaPsi(pre_data);
    JXF_Int theta_rho = jxf_CombinedPrecondDataThetaRho(pre_data);
    JXF_Int theta_phi = jxf_CombinedPrecondDataThetaPhi(pre_data);
    JXF_Int euclid_level = jxf_CombinedPrecondDataEuclidLevel(pre_data);
    JXF_Int amg_max_levels = jxf_CombinedPrecondDataAMGMaxLevels(pre_data);
    JXF_Int amg_relax_type = jxf_CombinedPrecondDataAMGRelaxType(pre_data);
    JXF_Int amg_print_level = jxf_CombinedPrecondDataAMGPrintLevel(pre_data);
    JXF_Int amg_interp_type = jxf_CombinedPrecondDataAMGInterpType(pre_data);
    JXF_Int amg_P_max_elmts = jxf_CombinedPrecondDataAMGPMaxElmts(pre_data);
    JXF_Int amg_measure_type = jxf_CombinedPrecondDataAMGMeasureType(pre_data);
    JXF_Int amg_coarsen_type = jxf_CombinedPrecondDataAMGCoarsenType(pre_data);
    JXF_Int amg_agg_num_levels = jxf_CombinedPrecondDataAMGAggNumLevels(pre_data);
    JXF_Int amg_coarse_threshold = jxf_CombinedPrecondDataAMGCoarseThreshold(pre_data);
    JXF_Int global_num_rows = jxf_hpCSRMatrixGlobalNumRows(par_matrix);
    JXF_Real drop_tol = jxf_CombinedPrecondDataDropTol(pre_data);
    JXF_Real theta_dis = jxf_CombinedPrecondDataThetaDis(pre_data);
    JXF_Real amg_strong_threshold = jxf_CombinedPrecondDataAMGStrongThreshold(pre_data);
    jxf_CSRMatrix *ser_matrix = NULL;
    JXF_Real *max_dvd_min = NULL;
    JXF_Int *mat_dis_psi = NULL, *max_mat_dis = NULL;
    JXF_Int row, pos_srt, pos_num, mat_flg = 0, tmp_row = 0;
    JXF_Int mat_psi, max_mat_psi, mat_rho, mat_phi;
    JXF_Int num_procs;
    JXF_Real fabs_max, fabs_min;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_CombinedPrecondDataActualPsi(pre_data) = 0; // undecided
    jxf_CombinedPrecondDataActualRho(pre_data) = 0; // undecided
    jxf_CombinedPrecondDataActualPhi(pre_data) = 0; // undecided
    ser_matrix = jxf_MergeDiagAndOffdDropSmall(jxf_hpCSRMatrixPar(par_matrix), 0.0); // Delete zero-elements
    max_dvd_min = jxf_CTAlloc(JXF_Real, num_rows);
    for (row = 0; row < num_rows; row ++)
    {
        pos_srt = jxf_CSRMatrixI(ser_matrix)[row] + 1;
        pos_num = jxf_CSRMatrixI(ser_matrix)[row+1] - pos_srt;
        fabs_max = jxf_DoubleArrayAbsMaxElement(&jxf_CSRMatrixData(ser_matrix)[pos_srt], pos_num);
        fabs_min = jxf_DoubleArrayAbsMinElement(&jxf_CSRMatrixData(ser_matrix)[pos_srt], pos_num);
        max_dvd_min[row] = fabs_max / fabs_min;
    }
    mat_psi = (JXF_Int)log10(jxf_DoubleArrayMaxElement(max_dvd_min, num_rows));
    jxf_MPI_Allreduce(&mat_psi, &max_mat_psi, 1, JXF_MPI_INT, MPI_MAX, comm);
    jxf_CombinedPrecondDataActualPsi(pre_data) = max_mat_psi;
    if (max_mat_psi < theta_psi) // Solve by AMG, Yue Xiaoqiang 2014/10/24
    {
        JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
        JXF_PAMGSetMaxLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_max_levels);
        JXF_PAMGSetMaxIter(jxf_CombinedPrecondDataPAMGData(pre_data), 1);
        JXF_PAMGSetMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_measure_type);
        JXF_PAMGSetCoarsenType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarsen_type);
        JXF_PAMGSetInterpType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_interp_type);
        JXF_PAMGSetPMaxElmts(jxf_CombinedPrecondDataPAMGData(pre_data), amg_P_max_elmts);
        JXF_PAMGSetAggNumLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_agg_num_levels);
        JXF_PAMGSetStrongThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_strong_threshold);
        JXF_PAMGSetPrintLevel(jxf_CombinedPrecondDataPAMGData(pre_data), amg_print_level);
        JXF_PAMGSetCoarseThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarse_threshold);
        JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 1);
        JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 2);
        JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 3);
        JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 1);
        JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 2);
        JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), 9, 3);
        JXF_PAMGSetup(jxf_CombinedPrecondDataPAMGData(pre_data), (JXF_hpCSRMatrix)par_matrix);
        jxf_CombinedPrecondDataPreID(pre_data) = 3;
    }
    else
    {
        mat_psi = max_mat_psi + 1;
        mat_dis_psi = jxf_CTAlloc(JXF_Int, mat_psi);
        memset(mat_dis_psi, 0, mat_psi*sizeof(JXF_Int));
        for (row = 0; row < num_rows; row ++)
        {
            mat_dis_psi[(JXF_Int)log10(max_dvd_min[row])] ++;
        }
        max_mat_dis = jxf_CTAlloc(JXF_Int, mat_psi);
        memset(max_mat_dis, 0, mat_psi*sizeof(JXF_Int));
        jxf_MPI_Allreduce(mat_dis_psi, max_mat_dis, mat_psi, JXF_MPI_INT, MPI_SUM, comm);
        mat_rho = 0;
        for (row = 0; row < mat_psi; row ++)
        {
            if (max_mat_dis[row] >= theta_dis*global_num_rows)
            {
                mat_rho ++;
            }
        }
        jxf_CombinedPrecondDataActualRho(pre_data) = mat_rho;
        if (mat_rho < theta_rho || mat_rho == 1) // Solve by AMG, Yue Xiaoqiang 2014/10/24
        {
            JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
            JXF_PAMGSetMaxLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_max_levels);
            JXF_PAMGSetMaxIter(jxf_CombinedPrecondDataPAMGData(pre_data), 1);
            JXF_PAMGSetMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_measure_type);
            JXF_PAMGSetCoarsenType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarsen_type);
            JXF_PAMGSetInterpType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_interp_type);
            JXF_PAMGSetPMaxElmts(jxf_CombinedPrecondDataPAMGData(pre_data), amg_P_max_elmts);
            JXF_PAMGSetAggNumLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_agg_num_levels);
            JXF_PAMGSetStrongThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_strong_threshold);
            JXF_PAMGSetPrintLevel(jxf_CombinedPrecondDataPAMGData(pre_data), amg_print_level);
            JXF_PAMGSetCoarseThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarse_threshold);
            JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 1);
            JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 2);
            JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 3);
            JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 1);
            JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 2);
            JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), 9, 3);
            JXF_PAMGSetup(jxf_CombinedPrecondDataPAMGData(pre_data), (JXF_hpCSRMatrix)par_matrix);
            jxf_CombinedPrecondDataPreID(pre_data) = 3;
        }
        else
        {
            mat_phi = 0;
            for (row = 0; row < mat_psi; row ++)
            {
                if (max_mat_dis[row] >= theta_dis*global_num_rows)
                {
                    if (mat_flg == 0)
                    {
                        mat_flg = 1;
                        tmp_row = row;
                    }
                    else
                    {
                        mat_phi += (row - tmp_row - 1);
                        tmp_row = row;
                    }
                }
            }
            jxf_CombinedPrecondDataActualPhi(pre_data) = mat_phi;
            if (mat_phi < theta_phi) // Solve by AMG, Yue Xiaoqiang 2014/10/24
            {
                JXF_PAMGCreate(&jxf_CombinedPrecondDataPAMGData(pre_data));
                JXF_PAMGSetMaxLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_max_levels);
                JXF_PAMGSetMaxIter(jxf_CombinedPrecondDataPAMGData(pre_data), 1);
                JXF_PAMGSetMeasureType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_measure_type);
                JXF_PAMGSetCoarsenType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarsen_type);
                JXF_PAMGSetInterpType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_interp_type);
                JXF_PAMGSetPMaxElmts(jxf_CombinedPrecondDataPAMGData(pre_data), amg_P_max_elmts);
                JXF_PAMGSetAggNumLevels(jxf_CombinedPrecondDataPAMGData(pre_data), amg_agg_num_levels);
                JXF_PAMGSetStrongThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_strong_threshold);
                JXF_PAMGSetPrintLevel(jxf_CombinedPrecondDataPAMGData(pre_data), amg_print_level);
                JXF_PAMGSetCoarseThreshold(jxf_CombinedPrecondDataPAMGData(pre_data), amg_coarse_threshold);
                JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 1);
                JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 2);
                JXF_PAMGSetCycleNumSweeps(jxf_CombinedPrecondDataPAMGData(pre_data), 1, 3);
                JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 1);
                JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), amg_relax_type, 2);
                JXF_PAMGSetCycleRelaxType(jxf_CombinedPrecondDataPAMGData(pre_data), 9, 3);
                JXF_PAMGSetup(jxf_CombinedPrecondDataPAMGData(pre_data), (JXF_hpCSRMatrix)par_matrix);
                jxf_CombinedPrecondDataPreID(pre_data) = 3;
            }
            else if (pre_id == 1) // Solve by Euclid firstly, Yue Xiaoqiang 2014/10/24
            {
                JXF_EuclidCreate(comm, &jxf_CombinedPrecondDataEuclidData(pre_data));
                JXF_EuclidSetLevel(jxf_CombinedPrecondDataEuclidData(pre_data), euclid_level);
                JXF_EuclidSetSparseA(jxf_CombinedPrecondDataEuclidData(pre_data), drop_tol);
                JXF_EuclidSetup(jxf_CombinedPrecondDataEuclidData(pre_data), (JXF_hpCSRMatrix)par_matrix);
            }
            else if (pre_id == 4) // Solve by ILU(0) - crossed firstly, Yue Xiaoqiang 2014/12/03
            {
                JXF_ILUZeroFactorDataCreate(&jxf_CombinedPrecondDataILUZData(pre_data), comm);
                JXF_ILUZeroFactorDataSetNxy(jxf_CombinedPrecondDataILUZData(pre_data), nx, ny);
                JXF_ILUZeroFactorDataSetNpxy(jxf_CombinedPrecondDataILUZData(pre_data), npx, npy);
                JXF_ILUZeroFactorDataSetNumEquns(jxf_CombinedPrecondDataILUZData(pre_data), num_equns);
                JXF_ILUZeroFactorDataSetDropTol(jxf_CombinedPrecondDataILUZData(pre_data), drop_tol);
                JXF_ILUZeroFactorDataGenerateParGrid(jxf_CombinedPrecondDataILUZData(pre_data));
                JXF_ILUZeroFactorDataSetup(jxf_CombinedPrecondDataILUZData(pre_data),(JXF_hpCSRMatrix)par_matrix);
            }
        }
    }
    jxf_TFree(max_dvd_min);
    if (mat_dis_psi)
    {
        jxf_TFree(mat_dis_psi);
    }
    if (max_mat_dis)
    {
        jxf_TFree(max_mat_dis);
    }
    jxf_CSRMatrixDestroy(ser_matrix);
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataSolve
 */
JXF_Int
JXF_CombinedPrecondDataSolve( JXF_Solver       solver,
                             JXF_hpCSRMatrix par_matrix,
                             JXF_ParVector    par_b,
                             JXF_ParVector    par_x )
{
    return( jxf_CombinedPrecondDataSolve( (void *) solver,
                                         (jxf_hpCSRMatrix *) par_matrix,
                                         (jxf_ParVector *)par_b,
                                         (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataSolve
 */
JXF_Int
jxf_CombinedPrecondDataSolve( void            *data,
                             jxf_hpCSRMatrix *par_matrix,
                             jxf_ParVector    *par_b,
                             jxf_ParVector    *par_x )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    
    if (pre_id == 11) // PAMG-Euclid-PAMG
    {
        jxf_PAMGEuclidPAMGSolve(data,par_matrix, par_b, par_x);
    }
    else if (pre_id == 12) // Euclid-PAMG
    {
        jxf_EuclidPAMGSolve(data, par_matrix, par_b, par_x);
    }
    else if (pre_id == 13) // PAMG-Euclid
    {
        jxf_PAMGEuclidSolve(data,par_matrix, par_b, par_x);
    }
    else if (pre_id == 14) // Euclid-PAMG-Euclid
    {
        jxf_EuclidPAMGEuclidSolve(data, par_matrix, par_b, par_x);
    }
    else if (pre_id == 15) // Euclid
    {
        JXF_EuclidSolve(jxf_CombinedPrecondDataEuclidData(pre_data),
                (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    }
    else if (pre_id == 16) // ILU-PAMG
    {
        jxf_ILUZPAMGSolve(data, par_matrix, par_b, par_x);
    }
    else if (pre_id == 17) // ILU-PAMG-ILU
    {
        jxf_ILUZPAMGILUZSolve(data, par_matrix, par_b, par_x);
    }
    else if (pre_id == 18) // PAMG-Euclid-PAMG, where Euclid for (A+A^T)/2
    {
        jxf_PAMGEuclidPAMGSolve(data, par_matrix, par_b, par_x);
    }
    else if (pre_id == 19) // PAMG-ILU-PAMG
    {
        jxf_PAMGILUZPAMGSolve(data, par_matrix, par_b, par_x);
    }
    else
    {
        jxf_printf("\n Warning: Wrong preconditioner ID\n\n");
        exit(-3);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataAdaptiveSolve
 */
JXF_Int
JXF_CombinedPrecondDataAdaptiveSolve( JXF_Solver       solver,
                                     JXF_hpCSRMatrix hp_matrix,
                                     JXF_ParVector    par_b,
                                     JXF_ParVector    par_x )
{
    return( jxf_CombinedPrecondDataAdaptiveSolve( (void *) solver,
                                                 (jxf_hpCSRMatrix *) hp_matrix,
                                                 (jxf_ParVector *)par_b,
                                                 (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataAdaptiveSolve
 */
JXF_Int
jxf_CombinedPrecondDataAdaptiveSolve( void            *data,
                                     jxf_hpCSRMatrix *par_matrix,
                                     jxf_ParVector    *par_b,
                                     jxf_ParVector    *par_x )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Int pre_id = jxf_CombinedPrecondDataPreID(pre_data);
    
    if (pre_id == 1) // Euclid
    {
        JXF_EuclidSolve(jxf_CombinedPrecondDataEuclidData(pre_data),
                   (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    }
    else if (pre_id == 3) // AMG
    {
        JXF_PAMGPrecond(jxf_CombinedPrecondDataPAMGData(pre_data),
                   (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    }
    else if (pre_id == 4) // ILU(0) - crossed
    {
        JXF_ILUZeroFactorDataPrecond(jxf_CombinedPrecondDataILUZData(pre_data),
                   (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    }
    else
    {
        jxf_printf("\n Warning: Wrong adaptive preconditioner ID\n\n");
        exit(-3);
    }
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_CombinedPrecondDataDestroy
 */
JXF_Int
JXF_CombinedPrecondDataDestroy( JXF_Solver solver )
{
    return( jxf_CombinedPrecondDataDestroy( (void *) solver ) );
}

/*!
 * \fn JXF_Int jxf_CombinedPrecondDataDestroy
 */
JXF_Int
jxf_CombinedPrecondDataDestroy( void *data )
{
    jxf_CombinedPrecondData *pre_data = data;
    
    if ( jxf_CombinedPrecondDataPAMGData(pre_data) )
    {
        JXF_PAMGDestroy( jxf_CombinedPrecondDataPAMGData(pre_data) );
    }
    if ( jxf_CombinedPrecondDataEuclidData(pre_data) )
    {
        JXF_EuclidDestroy( jxf_CombinedPrecondDataEuclidData(pre_data) );
    }
    if ( jxf_CombinedPrecondDataILUZData(pre_data) )
    {
        JXF_ILUZeroFactorDataDestroy( jxf_CombinedPrecondDataILUZData(pre_data) );
    }
    if ( jxf_CombinedPrecondDataAuxVector(pre_data) )
    {
        jxf_ParVectorDestroy( jxf_CombinedPrecondDataAuxVector(pre_data) );
    }
    if ( jxf_CombinedPrecondDataResVector(pre_data) )
    {
        jxf_ParVectorDestroy( jxf_CombinedPrecondDataResVector(pre_data) );
    }
    
    jxf_TFree(pre_data);
    
    return jxf_error_flag;
}
