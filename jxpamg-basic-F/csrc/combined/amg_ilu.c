//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  amg_ilu.c
 *  Date: 2014/03/24
 *
 *  Created by Yue Xiaoqiang
 *
 */

#include "jxf_combined.h"

/*!
 * \fn JXF_Int JXF_ILUZPAMGSetup
 */
JXF_Int
JXF_ILUZPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_ILUZPAMGSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_ILUZPAMGSetup
 */
JXF_Int
jxf_ILUZPAMGSetup( void *data, jxf_hpCSRMatrix *hp_matrix )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Real drop_tol = jxf_CombinedPrecondDataDropTol(pre_data);
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
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
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // SETUP phase for ILU
    JXF_ILUZeroFactorDataSetMaxIter(ilu_solver, 1);
    JXF_ILUZeroFactorDataSetDropTol(ilu_solver, drop_tol);
    JXF_ILUZeroFactorDataSetup(ilu_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(par_matrix, &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_ILUZPAMGILUZSetup
 */
JXF_Int
JXF_ILUZPAMGILUZSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_ILUZPAMGILUZSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

/*!
 * \fn JXF_Int jxf_ILUZPAMGILUZSetup
 */
JXF_Int
jxf_ILUZPAMGILUZSetup( void *data, jxf_hpCSRMatrix *hp_matrix )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Real drop_tol = jxf_CombinedPrecondDataDropTol(pre_data);
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, jxf_CombinedPrecondDataAMGCoarsenType(pre_data));
    JXF_PAMGSetInterpType(amg_solver, jxf_CombinedPrecondDataAMGInterpType(pre_data));
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, jxf_CombinedPrecondDataAMGStrongThreshold(pre_data));
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 3);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, jxf_CombinedPrecondDataAMGRelaxType(pre_data), 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, jxf_CombinedPrecondDataAMGRelaxType(pre_data), 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // SETUP phase for ILU
    JXF_ILUZeroFactorDataSetMaxIter(ilu_solver, 1);
    JXF_ILUZeroFactorDataSetDropTol(ilu_solver, drop_tol);
    JXF_ILUZeroFactorDataSetup(ilu_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(par_matrix, &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

JXF_Int
JXF_PAMGILUZPAMGSetup( JXF_Solver solver, JXF_hpCSRMatrix par_matrix )
{
    return( jxf_PAMGILUZPAMGSetup( (void *) solver, (jxf_hpCSRMatrix *) par_matrix ) );
}

JXF_Int
jxf_PAMGILUZPAMGSetup( void *data, jxf_hpCSRMatrix *hp_matrix )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    MPI_Comm comm = jxf_CombinedPrecondDataComm(pre_data);
    JXF_Real drop_tol = jxf_CombinedPrecondDataDropTol(pre_data);
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    JXF_Int *partition = NULL;
    
    // SETUP phase for PAMG
    JXF_PAMGSetMaxLevels(amg_solver, 25);
    JXF_PAMGSetMaxIter(amg_solver, 1);
    JXF_PAMGSetCycleType(amg_solver, 1);
    JXF_PAMGSetMeasureType(amg_solver, 0);
    JXF_PAMGSetCoarsenType(amg_solver, jxf_CombinedPrecondDataAMGCoarsenType(pre_data));
    JXF_PAMGSetInterpType(amg_solver, jxf_CombinedPrecondDataAMGInterpType(pre_data));
    JXF_PAMGSetSpMtRapType(amg_solver, 1);
    JXF_PAMGSetStrongThreshold(amg_solver, jxf_CombinedPrecondDataAMGStrongThreshold(pre_data));
    JXF_PAMGSetMaxRowSum(amg_solver, 0.9);
    JXF_PAMGSetPrintLevel(amg_solver, 3);
    JXF_PAMGSetRelaxWt(amg_solver, 1.0);
    JXF_PAMGSetOuterWt(amg_solver, 1.0);
    JXF_PAMGSetCoarseThreshold(amg_solver, 100);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);
    JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);
    JXF_PAMGSetCycleRelaxType(amg_solver, jxf_CombinedPrecondDataAMGRelaxType(pre_data), 1);
    JXF_PAMGSetCycleRelaxType(amg_solver, jxf_CombinedPrecondDataAMGRelaxType(pre_data), 2);
    JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);
    JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // SETUP phase for ILU
    JXF_ILUZeroFactorDataSetMaxIter(ilu_solver, 1);
    JXF_ILUZeroFactorDataSetDropTol(ilu_solver, drop_tol);
    JXF_ILUZeroFactorDataSetup(ilu_solver, (JXF_hpCSRMatrix)hp_matrix);
    
    // Allocate memory for auxiliary vector
    jxf_ParCSRMatrixGetRowPartitioning(par_matrix, &partition);
    jxf_CombinedPrecondDataAuxVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataAuxVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataAuxVector(pre_data));
    
    // Allocate memory for residual vector
    jxf_CombinedPrecondDataResVector(pre_data) = jxf_ParVectorCreate(comm,
                               jxf_ParCSRMatrixGlobalNumRows(par_matrix), partition);
    jxf_ParVectorSetPartitioningOwner(jxf_CombinedPrecondDataResVector(pre_data), 0);
    jxf_ParVectorInitialize(jxf_CombinedPrecondDataResVector(pre_data));
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_ILUZPAMGSolve
 */
JXF_Int
JXF_ILUZPAMGSolve( JXF_Solver       solver,
                  JXF_hpCSRMatrix par_matrix,
                  JXF_ParVector    par_b,
                  JXF_ParVector    par_x )
{
    return( jxf_ILUZPAMGSolve( (void *) solver,
                              (jxf_hpCSRMatrix *) par_matrix,
                              (jxf_ParVector *)par_b,
                              (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_ILUZPAMGSolve
 */
JXF_Int
jxf_ILUZPAMGSolve( void            *data,
                  jxf_hpCSRMatrix *hp_matrix,
                  jxf_ParVector    *par_b,
                  jxf_ParVector    *par_x )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One ILU(0) for A w = g, where A = par_matrix, g = par_b, w = par_x
    //jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_ILUZeroFactorDataPrecond(ilu_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_x, 1.0, res_vec);
    // 3. One V-cycle for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_ILUZPAMGILUZSolve
 */
JXF_Int
JXF_ILUZPAMGILUZSolve( JXF_Solver       solver,
                      JXF_hpCSRMatrix par_matrix,
                      JXF_ParVector    par_b,
                      JXF_ParVector    par_x )
{
    return( jxf_ILUZPAMGILUZSolve( (void *) solver,
                                  (jxf_hpCSRMatrix *) par_matrix,
                                  (jxf_ParVector *)par_b,
                                  (jxf_ParVector *)par_x ) );
}

/*!
 * \fn JXF_Int jxf_ILUZPAMGILUZSolve
 */
JXF_Int
jxf_ILUZPAMGILUZSolve( void            *data,
                      jxf_hpCSRMatrix *hp_matrix,
                      jxf_ParVector    *par_b,
                      jxf_ParVector    *par_x )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One ILU(0) for A w = g, where A = par_matrix, g = par_b, w = par_x
    //jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_ILUZeroFactorDataPrecond(ilu_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
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
    // 6. One ILU(0) for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    //jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_ILUZeroFactorDataPrecond(ilu_solver, (JXF_hpCSRMatrix)hp_matrix,(JXF_ParVector)res_vec,(JXF_ParVector)aux_vec);
    // 7. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    
    return jxf_error_flag;
}

JXF_Int
JXF_PAMGILUZPAMGSolve( JXF_Solver       solver,
                      JXF_hpCSRMatrix par_matrix,
                      JXF_ParVector    par_b,
                      JXF_ParVector    par_x )
{
    return( jxf_PAMGILUZPAMGSolve( (void *) solver,
                                  (jxf_hpCSRMatrix *) par_matrix,
                                  (jxf_ParVector *)par_b,
                                  (jxf_ParVector *)par_x ) );
}

JXF_Int
jxf_PAMGILUZPAMGSolve( void            *data,
                      jxf_hpCSRMatrix *hp_matrix,
                      jxf_ParVector    *par_b,
                      jxf_ParVector    *par_x )
{
    jxf_ParCSRMatrix *par_matrix = jxf_hpCSRMatrixPar(hp_matrix);
    jxf_CombinedPrecondData *pre_data = data;
    JXF_Solver ilu_solver = jxf_CombinedPrecondDataILUZData(pre_data);
    JXF_Solver amg_solver = jxf_CombinedPrecondDataPAMGData(pre_data);
    jxf_ParVector *aux_vec = jxf_CombinedPrecondDataAuxVector(pre_data);
    jxf_ParVector *res_vec = jxf_CombinedPrecondDataResVector(pre_data);
    
    // 1. One V-cycle for A w = g, where A = par_matrix, g = par_b, w = par_x
    jxf_ParVectorSetConstantValues(par_x, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)par_b, (JXF_ParVector)par_x);
    // 2. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_x, 1.0, res_vec);
    // 3. One ILU(0) for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    //jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_ILUZeroFactorDataPrecond(ilu_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 4. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);
    // 5. r = g - Aw, where r = res_vec, g = par_b, w = par_x
    jxf_ParVectorCopy(par_b, res_vec);
    jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_x, 1.0, res_vec);
    // 6. One V-cycle for A v = r, where A = par_matrix, r = res_vec, v = aux_vec
    jxf_ParVectorSetConstantValues(aux_vec, 0.0);
    JXF_PAMGPrecond(amg_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)res_vec, (JXF_ParVector)aux_vec);
    // 7. w = w + v, where w = par_x, v = aux_vec
    jxf_ParVectorAxpy(1.0, aux_vec, par_x);

    return jxf_error_flag;
}
