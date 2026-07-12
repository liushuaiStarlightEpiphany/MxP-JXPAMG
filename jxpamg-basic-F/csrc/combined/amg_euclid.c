//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  amg_euclid.c
 *  Date: 2013/12/11
 *
 *  Created by Yue Xiaoqiang
 *
 */

#include "jxf_combined.h"

/*!
 * \fn JXF_Int JXF_PAMGEuclidPAMGSetup
 */
JXF_Int
JXF_PAMGEuclidPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_PAMGEuclidPAMGSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_PAMGEuclidPAMGSetup
 */
JXF_Int
jxf_PAMGEuclidPAMGSetup( void *data, jxf_hpCSRMatrix *hp_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, 6);
    JXF_PAMGSetInterpType(amg_solver, 0);
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, 0.25);
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 0);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) hp_matrix);
    
    // SETUP phase for Euclid
    JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar(hp_matrix), &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(hp_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(hp_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_PAMGEuclidPAMGSetupB
 */
JXF_Int
JXF_PAMGEuclidPAMGSetupB( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_PAMGEuclidPAMGSetupB( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_PAMGEuclidPAMGSetupB
 */
JXF_Int
jxf_PAMGEuclidPAMGSetupB( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    //jxf_ParCSRMatrix *pre_euc_mat = jxf_CombinedPrecondDataPreEucMat(pre_data);
    jxf_hpCSRMatrix *hp_pre_euc_mat = jxf_hpInithpCSRMatrix();
    jxf_hpCSRMatrixPar(hp_pre_euc_mat) = jxf_CombinedPrecondDataPreEucMat(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, 6);
    JXF_PAMGSetInterpType(amg_solver, 0);
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, 0.25);
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 0);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)par_matrix);
    
    // SETUP phase for Euclid
    JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_pre_euc_mat);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar(par_matrix), &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_EuclidPAMGSetup
 */
JXF_Int
JXF_EuclidPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_EuclidPAMGSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_EuclidPAMGSetup
 */
JXF_Int
jxf_EuclidPAMGSetup( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, 6);
    JXF_PAMGSetInterpType(amg_solver, 0);
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, 0.25);
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 0);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) par_matrix);
    
    // SETUP phase for Euclid
    JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)par_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar(par_matrix), &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_PAMGEuclidSetup
 */
JXF_Int
JXF_PAMGEuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_PAMGEuclidSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_PAMGEuclidSetup
 */
JXF_Int
jxf_PAMGEuclidSetup( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, 6);
    JXF_PAMGSetInterpType(amg_solver, 0);
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, 0.25);
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 0);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) par_matrix);
    
    // SETUP phase for Euclid
    JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)par_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar(par_matrix), &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_EuclidPAMGEuclidSetup
 */
JXF_Int
JXF_EuclidPAMGEuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_EuclidPAMGEuclidSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_EuclidPAMGEuclidSetup
 */
JXF_Int
jxf_EuclidPAMGEuclidSetup( void *data, jxf_hpCSRMatrix *par_matrix )
{
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, 6);
    JXF_PAMGSetInterpType(amg_solver, 0);
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, 0.25);
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 0);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, 3, 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) par_matrix);
    
    // SETUP phase for Euclid
    JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)par_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(jxf_hpCSRMatrixPar(par_matrix), &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_hpCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_PAMGEuclidPAMGSolve
 */
JXF_Int
JXF_PAMGEuclidPAMGSolve( JXF_Solver       solver,
                        JXF_hpCSRMatrix par_matrix,
                        JXF_ParVector    par_b,
                        JXF_ParVector    par_x )
{
    return( jxf_PAMGEuclidPAMGSolve( (void *) solver,
                                    (jxf_hpCSRMatrix *) par_matrix,
                                    (jxf_ParVector *)par_b,
                                    (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_PAMGEuclidPAMGSolve
 */
JXF_Int
jxf_PAMGEuclidPAMGSolve( void            *data,
                        jxf_hpCSRMatrix *par_matrix,
                        jxf_ParVector    *par_b,
                        jxf_ParVector    *par_x )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One V-cycle for A w = g, where A = par_matrix, g = par_b, w = par_x
    jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar(par_matrix), par_x, 1.0, res_vec);
    // 3. One ILU(k) for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_EuclidSolve(euclid_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    // 5. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar(par_matrix), par_x, 1.0, res_vec);
    // 6. One V-cycle for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 7. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_EuclidPAMGSolve
 */
JXF_Int
JXF_EuclidPAMGSolve( JXF_Solver       solver,
                    JXF_hpCSRMatrix par_matrix,
                    JXF_ParVector    par_b,
                    JXF_ParVector    par_x )
{
    return( jxf_EuclidPAMGSolve( (void *) solver,
                                (jxf_hpCSRMatrix *) par_matrix,
                                (jxf_ParVector *)par_b,
                                (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_EuclidPAMGSolve
 */
JXF_Int
jxf_EuclidPAMGSolve( void            *data,
                    jxf_hpCSRMatrix *par_matrix,
                    jxf_ParVector    *par_b,
                    jxf_ParVector    *par_x )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One ILU(k) for A w = g, where A = par_matrix, g = par_b, w = par_x
    jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_EuclidSolve(euclid_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar(par_matrix), par_x, 1.0, res_vec);
    // 3. One V-cycle for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_PAMGEuclidSolve
 */
JXF_Int
JXF_PAMGEuclidSolve( JXF_Solver       solver,
                    JXF_hpCSRMatrix par_matrix,
                    JXF_ParVector    par_b,
                    JXF_ParVector    par_x )
{
    return( jxf_PAMGEuclidSolve( (void *) solver,
                                (jxf_hpCSRMatrix *) par_matrix,
                                (jxf_ParVector *)par_b,
                                (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_PAMGEuclidSolve
 */
JXF_Int
jxf_PAMGEuclidSolve( void            *data,
                    jxf_hpCSRMatrix *par_matrix,
                    jxf_ParVector    *par_b,
                    jxf_ParVector    *par_x )
{
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One V-cycle for A w = g, where A = par_matrix, g = par_b, w = par_x
    //jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, jxf_hpCSRMatrixPar(par_matrix), par_x, 1.0, res_vec);
    // 3. One ILU(k) for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    //jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_EuclidSolve(euclid_solver, (JXF_hpCSRMatrix)par_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_EuclidPAMGEuclidSolve
 */
JXF_Int
JXF_EuclidPAMGEuclidSolve( JXF_Solver       solver,
                          JXF_hpCSRMatrix par_matrix,
                          JXF_ParVector    par_b,
                          JXF_ParVector    par_x )
{
    return( jxf_EuclidPAMGEuclidSolve( (void *) solver,
                                      (jxf_hpCSRMatrix *) par_matrix,
                                      (jxf_ParVector *)par_b,
                                      (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_EuclidPAMGEuclidSolve
 */
JXF_Int
jxf_EuclidPAMGEuclidSolve( void            *data,
                          jxf_hpCSRMatrix *hp_matrix,
                          jxf_ParVector    *par_b,
                          jxf_ParVector    *par_x )
{
    jxf_ParCSRMatrix* par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Solver euclid_solver = jxf_CombinedPrecondDataEuclidData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One ILU(k) for A w = g, where A = par_matrix, g = par_b, w = par_x
    jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_EuclidSolve(euclid_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_x, 1.0, res_vec);
    // 3. One V-cycle for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    // 5. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_x, 1.0, res_vec);
    // 6. One ILU(k) for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_EuclidSolve(euclid_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 7. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}
