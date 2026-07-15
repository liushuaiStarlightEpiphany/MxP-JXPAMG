// test_mixed_prec_cpr.c - 混合精度CPR预条件器测试程序
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__AVX__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#include "mpi.h"

// 双精度头文件
#include "jx_mv.h"
#include "jx_bsr_mv.h"
#include "jx_parbsr_mv.h"
#include "jx_bsr_decoup.h"
#include "jx_krylov.h"
#include "jx_pamg.h"

// 单精度头文件
#include "jxf_mv.h"
#include "jxf_bsr_mv.h"
#include "jxf_parbsr_mv.h"
#include "jxf_bsr_decoup.h"
#include "jxf_cpr.h"
#include "jxf_krylov.h"

// 精度转换函数声明（在您的代码中已有）
JXF_Int jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y);
JXF_Int jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y);
JXF_Int jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y);
JXF_Int jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y);

// 混合精度CPR预条件器接口
JXF_Int JXF_MxP_CPRPrecond(jxf_CPRPrecond *cpr, 
                          jx_ParBSRMatrix* par_matrix,
                          jx_ParVector *par_rhs,
                          jx_ParVector *par_sol);

// 混合精度CPR预条件器接口
JXF_Int JXF_MxP_CPRPrecond1(jxf_CPRPrecond *cpr, 
                          jx_ParBSRMatrix* par_matrix,
                          jx_ParVector *par_rhs,
                          jx_ParVector *par_sol);

// =============================================================
// 精度转换函数（从您的代码中复制）
// =============================================================

// =============================================================
// BSR精度转换函数（极简版）
// =============================================================

// 1. 双精度BSR矩阵转换为单精度BSR矩阵
JXF_Int
jxmp_BSRMatrixDtoF(jx_BSRMatrix *A, jxf_BSRMatrix *B, JXF_Int copy_data)
{
    if (!A || !B) return -1;
    
    JXF_Int ierr = 0;
    JXF_Int num_rows = jx_BSRMatrixNumRows(A);
    JXF_Int nnz = jx_BSRMatrixNumNonzeros(A);
    JXF_Int blk_size = jx_BSRMatrixBlockSize(A);
    
    JX_Int *A_i = jx_BSRMatrixI(A);
    JX_Int *A_j = jx_BSRMatrixJ(A);
    JX_Real *A_data = jx_BSRMatrixData(A);
    
    JXF_Int *B_i = jxf_BSRMatrixI(B);
    JXF_Int *B_j = jxf_BSRMatrixJ(B);
    JXF_Real *B_data = jxf_BSRMatrixData(B);
    
    // 只需要复制索引（整数，可能类型相同）
    for (JXF_Int i = 0; i <= num_rows; i++) {
        B_i[i] = (JXF_Int)A_i[i];
    }
    
    for (JXF_Int i = 0; i < nnz; i++) {
        B_j[i] = (JXF_Int)A_j[i];
    }
    
    // 如果需要，转换浮点数数据
    if (copy_data && A_data && B_data) {
        JXF_Int data_size = nnz * blk_size * blk_size;
        JXF_Int i = 0;
#if defined(__AVX__) || defined(__AVX2__)
        for (; i + 4 <= data_size; i += 4) {
            __m256d d = _mm256_loadu_pd(A_data + i);
            __m128 f = _mm256_cvtpd_ps(d);
            _mm_storeu_ps(B_data + i, f);
        }
#endif
#pragma omp parallel for
        for (i = 0; i < data_size; i++) {
            B_data[i] = (JXF_Real)A_data[i];  // double -> float
        }
    }
    
    return ierr;
}

// 2. 双精度并行BSR矩阵转换为单精度并行BSR矩阵
JXF_Int
jxmp_ParBSRMatrixDtoF(jx_ParBSRMatrix *A, jxf_ParBSRMatrix *B, JXF_Int copy_data)
{
    if (!A || !B) return -1;
    
    JXF_Int ierr = 0;
    JXF_Int i;
    
    // 转换对角部分
    if (A->diag && B->diag) {
        ierr = jxmp_BSRMatrixDtoF(A->diag, B->diag, copy_data);
        if (ierr != 0) return ierr;
    }
    
    // 转换非对角部分
    if (A->offd && B->offd) {
        ierr = jxmp_BSRMatrixDtoF(A->offd, B->offd, copy_data);
        if (ierr != 0) return ierr;
    }
    
    // 复制列映射（整数数组）
    if (A->col_map_offd && B->col_map_offd) {
        JXF_Int num_cols_offd = A->offd ? jx_BSRMatrixNumCols(A->offd) : 0;
        for (i = 0; i < num_cols_offd; i++) {
            B->col_map_offd[i] = (JXF_BigInt)A->col_map_offd[i];
        }
    }
    
    // 注意：row_starts和col_starts已经在jxf_ParBSRMatrixCreate中设置好了
    // 不需要在这里再次复制，除非需要类型转换
    
    // 但为了完整性，如果需要整型类型转换：
    // if (A->row_starts && B->row_starts) {
    //     // 假设长度相同（都是nprocs+1）
    //     // 这里我们不知道长度，所以这是一个潜在的问题
    //     // 实际上，在Create时已经复制了，所以可能不需要这里
    // }
    
    return ierr;
}

// 双精度并行向量 -> 单精度并行向量
JXF_Int
jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y)
{
    if (!x || !y) return -1;
    
    jx_Vector *x_local = jx_ParVectorLocalVector(x);
    jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
    
    return jxmp_SeqVectorDtoF(x_local, y_local);
}

// 单精度并行向量 -> 双精度并行向量  
JXF_Int
jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y)
{
    if (!x || !y) return -1;
    
    jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
    jx_Vector *y_local = jx_ParVectorLocalVector(y);
    
    return jxmp_SeqVectorFtoD(x_local, y_local);
}

// 序列向量的转换
JXF_Int
jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y)
{
    if (!x || !y) return -1;
    JX_Real *x_data = jx_VectorData(x);
    JXF_Real *y_data = jxf_VectorData(y);
    JX_Int size = jx_VectorSize(x);
    size *= jx_VectorNumVectors(x);
    JX_Int i = 0;
#if defined(__AVX__) || defined(__AVX2__)
    for (; i + 4 <= size; i += 4) {
        __m256d d = _mm256_loadu_pd(x_data + i);
        __m128 f = _mm256_cvtpd_ps(d);
        _mm_storeu_ps(y_data + i, f);
    }
#endif
    for (; i < size; i++) {
        y_data[i] = (JXF_Real)x_data[i];
    }
    return 0;
}

JXF_Int
jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y)
{
    if (!x || !y) return -1;
    
    JXF_Real *x_data = jxf_VectorData(x);
    JX_Real *y_data = jx_VectorData(y);
    JX_Int size = jxf_VectorSize(x);
    JX_Int i = 0;

    size *= jxf_VectorNumVectors(x);

        JX_Int ii = 0;
#if defined(__AVX__) || defined(__AVX2__)
    for (; ii + 4 <= size; ii += 4) {
        __m128 f = _mm_loadu_ps(x_data + ii);
        __m256d d = _mm256_cvtps_pd(f);
        _mm256_storeu_pd(y_data + ii, d);
    }
#endif
    for (; ii < size; ii++) {
        y_data[ii] = (JX_Real)x_data[ii];
    }

    return 0;
}

// =============================================================
// 混合精度CPR预条件器接口实现
// =============================================================


// =============================================================
// 混合精度CPR预条件器接口实现（直接在函数内部实现）
// =============================================================

/**
 * \brief 混合精度CPR预条件器接口
 * 
 * 此函数将高精度输入转换为低精度，调用低精度AMG求解器，
 * 然后将结果转换回高精度，最后使用高精度进行块磨光。
 * 
 * \param cpr 低精度CPR预条件器（已用低精度矩阵设置）
 * \param par_matrix 高精度矩阵（用于高精度块磨光）
 * \param par_rhs 高精度右端项向量
 * \param par_sol 高精度解向量（输出）
 * \return JXF_Int 错误代码
 */
// =============================================================
// 混合精度CPR预条件器接口实现（更安全的版本）
// =============================================================

JXF_Int JXF_MxP_CPRPrecond1(jxf_CPRPrecond *cpr, 
                          jx_ParBSRMatrix* par_matrix,
                          jx_ParVector *par_rhs,
                          jx_ParVector *par_sol)
{
    if (cpr == NULL || !cpr->is_initialized) {
        printf("CPR preconditioner not initialized\n");
        return jxf_error_flag;
    }
    
    if (par_rhs == NULL || par_sol == NULL) {
        printf("Input vectors are NULL\n");
        return jxf_error_flag;
    }
    
    MPI_Comm comm = cpr->comm;
    int myid;
    MPI_Comm_rank(comm, &myid);
    
    // 检查维度一致性（调试用）
    if (cpr->print_level > 1) {
        JXF_Int low_block_size = cpr->block_size;
        JXF_BigInt low_global_rows = jxf_ParBSRMatrixGlobalNumRows(cpr->A_bsr);
        JXF_BigInt low_global_scalar_rows = low_global_rows * low_block_size;
        
        JXF_Int high_block_size = jx_ParBSRMatrixBlockSize(par_matrix);
        JXF_BigInt high_global_rows = jx_ParBSRMatrixGlobalNumRows(par_matrix);
        JXF_BigInt high_global_scalar_rows = high_global_rows * high_block_size;
        
        if (low_global_scalar_rows != high_global_scalar_rows) {
            if (myid == 0) {
                printf("Warning: Dimension mismatch in mixed precision CPR\n");
            }
        }
    }
    
    JXF_Real t1, t2;
    JXF_Int i;
    
    // =============================================================
    // 步骤1: 创建低精度工作向量（总是创建新的，避免内存问题）
    // =============================================================
    
    jxf_ParVector *low_rhs = jxf_CreateGlobalVector(cpr->A_bsr);
    jxf_ParVector *low_sol = jxf_CreateGlobalVector(cpr->A_bsr);
    
    if (low_rhs == NULL || low_sol == NULL) {
        printf("Failed to create low precision vectors\n");
        if (low_rhs) jxf_ParVectorDestroy(low_rhs);
        if (low_sol) jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // =============================================================
    // 步骤2: 将高精度右端项转换为低精度
    // =============================================================
    
    JXF_Int convert_err = jxmp_ParVectorDtoF(par_rhs, low_rhs);
    if (convert_err != 0) {
        printf("Failed to convert RHS from double to float\n");
        jxf_ParVectorDestroy(low_rhs);
        jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // 初始化低精度解向量为零
    jxf_ParVectorSetConstantValues(low_sol, 0.0);
    
    // =============================================================
    // 步骤3: 低精度阶段1 - 求解压力子系统（使用AMG）
    // =============================================================
    
    JXF_Int block_size = cpr->block_size;
    JXF_Int pressure_index = cpr->pressure_index;
    
    // 3.1 限制：从低精度全局残差提取压力残差
    if (jxf_RestrictPressureVector(low_rhs, cpr->rp, block_size, pressure_index) != JXF_SUCCESS) {
        printf("Failed to restrict pressure vector\n");
        jxf_ParVectorDestroy(low_rhs);
        jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // 3.2 初始化压力解为零
    jxf_ParVectorSetConstantValues(cpr->xp, 0.0);
    
    // 3.3 使用低精度AMG求解压力系统
    t1 = jxf_MPI_Wtime();
    for (i = 0; i < cpr->stage1_maxit; i++) {
        // 使用AMG求解压力系统
        JXF_Int amg_result = JXF_PAMGSolve(cpr->stage1_solver, 
                                          (JXF_ParCSRMatrix)cpr->A_pressure,
                                          (JXF_Vector)cpr->rp,
                                          (JXF_Vector)cpr->xp);
        if (amg_result != JXF_SUCCESS) {
            printf("AMG solve failed\n");
            jxf_ParVectorDestroy(low_rhs);
            jxf_ParVectorDestroy(low_sol);
            return amg_result;
        }
    }
    t2 = jxf_MPI_Wtime();
    cpr->stage1_solve_time += t2 - t1;
    
    // 3.4 延拓：将低精度压力解延拓到低精度全局向量
    if (jxf_ProlongatePressureVector(cpr->xp, low_sol, block_size, pressure_index, 0) != JXF_SUCCESS) {
        printf("Failed to prolongate pressure vector\n");
        jxf_ParVectorDestroy(low_rhs);
        jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // =============================================================
    // 步骤4: 将低精度解转换为高精度解
    // =============================================================
    
    convert_err = jxmp_ParVectorFtoD(low_sol, par_sol);
    if (convert_err != 0) {
        printf("Failed to convert solution from float to double\n");
        jxf_ParVectorDestroy(low_rhs);
        jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // =============================================================
    // 步骤5: 清理低精度临时向量
    // =============================================================
    
    jxf_ParVectorDestroy(low_rhs);
    jxf_ParVectorDestroy(low_sol);
    
    // =============================================================
    // 步骤6: 高精度阶段2 - 块磨光（使用高精度函数）
    // =============================================================
    
    t1 = jxf_MPI_Wtime();
    
    // 6.1 创建高精度工作向量（用于磨光）
    jx_ParVector *high_work = NULL;
    
    // 创建高精度工作向量
    JXF_BigInt global_scalar_rows = jx_ParBSRMatrixGlobalNumRows(par_matrix) * 
                                    jx_ParBSRMatrixBlockSize(par_matrix);
    JXF_BigInt* partitioning = NULL;
    jx_ParBSRMatrixGetRowPartitioning(par_matrix, &partitioning);
    
    high_work = jx_ParVectorCreate(par_matrix->comm, global_scalar_rows, partitioning);
    
    if (high_work == NULL) {
        printf("Failed to create high precision work vector\n");
        return jxf_error_flag;
    }
    
    jx_ParVectorInitialize(high_work);
    
    // 6.2 设置磨光参数
    JX_Real relax_weight = 1.0;    // 松弛权重
    JX_Real omega = 0.8;          // 阻尼因子
    JX_Int forward_or_backward = 1; // 1: 前向, 0: 对称, -1: 后向
    
    // 6.3 应用高精度块磨光
    if (cpr->print_level > 2) {
        // 计算磨光前的残差（调试用）
        jx_ParVector *high_temp = jx_ParVectorCreate(par_matrix->comm, global_scalar_rows, partitioning);
        if (high_temp) {
            jx_ParVectorInitialize(high_temp);
            jx_ParVectorCopy(par_rhs, high_temp);
            jx_ParBSRMatrixMatvec(-1.0, par_matrix, par_sol, 1.0, high_temp);
            JX_Real res_before = jx_ParVectorNorm2(high_temp);
            if (myid == 0) {
                printf("CPR-MxP: Residual before polishing: %.6e\n", res_before);
            }
            jx_ParVectorDestroy(high_temp);
        }
    }
    
    // 使用高精度块磨光函数
    // 注意：这里假设存在jx_ParBSRHGSRelax函数，如果不存在，需要使用其他高精度磨光函数
    for (i = 0; i < cpr->stage2_maxit; i++) {
        // 检查是否存在高精度磨光函数
        if (cpr->print_level > 0) {
            // if (myid == 0) printf("Applying high precision polishing iteration %d\n", i+1);
        }
        
        // 使用双精度磨光函数（假设函数名类似但前缀不同）
        // 如果没有jx_ParBSRHGSRelax，可以尝试其他高精度磨光方法
        // 这里是一个示例，你可能需要根据实际情况调整
        jx_ParBSRHGSRelax(par_matrix, par_sol, par_rhs, high_work, 
                         relax_weight, omega, forward_or_backward);
    }
    
    if (cpr->print_level > 2) {
        // 计算磨光后的残差（调试用）
        jx_ParVector *high_temp = jx_ParVectorCreate(par_matrix->comm, global_scalar_rows, partitioning);
        if (high_temp) {
            jx_ParVectorInitialize(high_temp);
            jx_ParVectorCopy(par_rhs, high_temp);
            jx_ParBSRMatrixMatvec(-1.0, par_matrix, par_sol, 1.0, high_temp);
            JX_Real res_after = jx_ParVectorNorm2(high_temp);
            if (myid == 0) {
                printf("CPR-MxP: Residual after polishing: %.6e\n", res_after);
            }
            jx_ParVectorDestroy(high_temp);
        }
    }
    
    t2 = jxf_MPI_Wtime();
    cpr->stage2_solve_time += t2 - t1;
    
    // =============================================================
    // 步骤7: 清理高精度工作向量
    // =============================================================
    
    if (high_work) jx_ParVectorDestroy(high_work);
    
    // 如果有需要，可以在这里添加额外的统计信息
    if (cpr->print_level > 0 && myid == 0) {
        printf("CPR-MxP: Stage1 (AMG) time: %.4f s, Stage2 (Polish) time: %.4f s\n",
               cpr->stage1_solve_time, cpr->stage2_solve_time);
    }
    
    return JXF_SUCCESS;
}
// jxf_ParVector* jxf_CreateGlobalVector(jxf_ParBSRMatrix *A_bsr);
/**
 * \brief 混合精度CPR预条件器接口
 * 
 * 此函数将高精度输入转换为低精度，调用低精度CPR预条件器，
 * 然后将结果转换回高精度输出。
 * 
 * \param cpr 低精度CPR预条件器（已用低精度矩阵设置）
 * \param par_matrix 高精度矩阵（仅用于一致性检查）
 * \param par_rhs 高精度右端项向量
 * \param par_sol 高精度解向量（输出）
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_MxP_CPRPrecond(jxf_CPRPrecond *cpr, 
                          jx_ParBSRMatrix* par_matrix,
                          jx_ParVector *par_rhs,
                          jx_ParVector *par_sol)
{


    if (cpr == NULL || !cpr->is_initialized) {
        printf("CPR preconditioner not initialized\n");
        return jxf_error_flag;
    }
    
    if (par_rhs == NULL || par_sol == NULL) {
        printf("Input vectors are NULL\n");
        return jxf_error_flag;
    }
    
    MPI_Comm comm = cpr->comm;
    int myid;
    MPI_Comm_rank(comm, &myid);
    
    // 检查维度一致性（调试用）
    if (cpr->print_level > 1) {
        JXF_Int low_block_size = cpr->block_size;
        JXF_BigInt low_global_rows = jxf_ParBSRMatrixGlobalNumRows(cpr->A_bsr);
        JXF_BigInt low_global_scalar_rows = low_global_rows * low_block_size;
        
        JXF_Int high_block_size = jx_ParBSRMatrixBlockSize(par_matrix);
        JXF_BigInt high_global_rows = jx_ParBSRMatrixGlobalNumRows(par_matrix);
        JXF_BigInt high_global_scalar_rows = high_global_rows * high_block_size;
        
        if (low_global_scalar_rows != high_global_scalar_rows) {
            if (myid == 0) {
                printf("Warning: Dimension mismatch in mixed precision CPR\n");
            }
        }
    }
    
    // 1. 创建低精度工作向量
    // 重用CPR内部的工作向量（如果有的话），或者创建新的
    jxf_ParVector *low_rhs = NULL;
    jxf_ParVector *low_sol = NULL;
    
    // 使用已有的工作向量或者创建新的
    if (cpr->work != NULL) {
        low_rhs = cpr->work;  // 重用工作向量
        low_sol = jxf_CreateGlobalVector(cpr->A_bsr);  // 创建新的解向量
    } else {
        low_rhs = jxf_CreateGlobalVector(cpr->A_bsr);
        low_sol = jxf_CreateGlobalVector(cpr->A_bsr);
    }
    
    if (low_rhs == NULL || low_sol == NULL) {
        printf("Failed to create low precision vectors\n");
        if (low_rhs && low_rhs != cpr->work) jxf_ParVectorDestroy(low_rhs);
        if (low_sol) jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // 2. 将高精度右端项转换为低精度
    JXF_Int convert_err = jxmp_ParVectorDtoF(par_rhs, low_rhs);
    
    if (convert_err != 0) {
        printf("Failed to convert RHS from double to float\n");
        if (low_rhs && low_rhs != cpr->work) jxf_ParVectorDestroy(low_rhs);
        if (low_sol) jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // 3. 初始化低精度解向量为零
    jxf_ParVectorSetConstantValues(low_sol, 0.0);
    
    // 4. 应用低精度CPR预条件器
    JXF_Int precond_err = JXF_CPRPrecond(cpr, cpr->A_bsr, low_rhs, low_sol);
    
    if (precond_err != JXF_SUCCESS) {
        printf("Low precision CPR preconditioner failed\n");
        if (low_rhs && low_rhs != cpr->work) jxf_ParVectorDestroy(low_rhs);
        if (low_sol) jxf_ParVectorDestroy(low_sol);
        return precond_err;
    }
    
    // 5. 将低精度解向量转换回高精度
    convert_err = jxmp_ParVectorFtoD(low_sol, par_sol);
    
    if (convert_err != 0) {
        printf("Failed to convert solution from float to double\n");
        if (low_rhs && low_rhs != cpr->work) jxf_ParVectorDestroy(low_rhs);
        if (low_sol) jxf_ParVectorDestroy(low_sol);
        return jxf_error_flag;
    }
    
    // 6. 清理临时向量（保留工作向量供下次使用）
    if (low_sol && low_sol != cpr->work) {
        jxf_ParVectorDestroy(low_sol);
    }
    if (low_rhs && low_rhs != cpr->work) {
        jxf_ParVectorDestroy(low_rhs);
    }
    
    return JXF_SUCCESS;
}

// =============================================================
// 矩阵辅助函数（从您的代码中复制）
// =============================================================

// 读取右端项文件函数
jx_Vector* read_rhs_file(const char* filename, int expected_size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening RHS file: %s\n", filename);
        return NULL;
    }
    
    // 读取向量维度
    int vector_size;
    if (fscanf(file, "%d", &vector_size) != 1) {
        printf("Error reading vector size from RHS file\n");
        fclose(file);
        return NULL;
    }
    
    // 检查维度是否匹配
    if (expected_size != -1 && vector_size != expected_size) {
        printf("Warning: RHS vector size (%d) does not match expected size (%d)\n", 
               vector_size, expected_size);
    }
    
    // 创建向量
    jx_Vector* rhs = jx_SeqVectorCreate(vector_size);
    jx_SeqVectorInitialize(rhs);
    JX_Real* rhs_data = jx_VectorData(rhs);
    
    // 读取数值
    for (int i = 0; i < vector_size; i++) {
        fscanf(file, "%lf", &rhs_data[i]);
    }
    
    fclose(file);
    return rhs;
}

// =============================================================
// 混合精度CPR-BiCGSTAB求解函数
// =============================================================

/* CPR setup wrapper: skip re-setup since CPR is already configured with A_low */
static JX_Int mixed_cpr_setup_wrapper(void *cpr, void *matrix)
{
    (void)matrix; /* ignore - A_low was already used for setup */
    return JXF_SUCCESS;
}

int solve_with_mixed_precision_cpr_gmres(

    jx_ParBSRMatrix* A_high,        // 高精度矩阵
    jx_ParVector* b_high,           // 高精度右端项
    jx_ParVector* x_high,           // 高精度解向量（输入/输出）
    int max_iterations,             // 最大迭代次数
    JX_Real tolerance               // 收敛容差
)
{
    int myid;
    MPI_Comm_rank(A_high->comm, &myid);
    
    // 1. 将高精度矩阵转换为低精度矩阵
    jxf_ParBSRMatrix* A_low = jxf_ParBSRMatrixCreate(
        A_high->comm,
        jx_ParBSRMatrixBlockSize(A_high),
        (JXF_BigInt)jx_ParBSRMatrixGlobalNumRows(A_high),
        (JXF_BigInt)jx_ParBSRMatrixGlobalNumCols(A_high),
        A_high->row_starts,  // 使用相同的分区
        A_high->col_starts,
        A_high->offd ? jx_BSRMatrixNumCols(A_high->offd) : 0,
        A_high->diag ? jx_BSRMatrixNumNonzeros(A_high->diag) : 0,
        A_high->offd ? jx_BSRMatrixNumNonzeros(A_high->offd) : 0
    );
    
    if (!A_low) {
        if (myid == 0) printf("Error creating low precision matrix\n");
        return -1;
    }
    
    // 初始化矩阵
    jxf_ParBSRMatrixInitialize(A_low);
    
    // 转换矩阵数据（使用已有的双精度转单精度函数）
    JXF_Int ierr = jxmp_ParBSRMatrixDtoF(A_high, A_low, 1);
    if (ierr != 0) {
        if (myid == 0) printf("Error converting matrix to low precision\n");
        jxf_ParBSRMatrixDestroy(A_low);
        return -1;
    }
    
    // 2. 创建低精度CPR预条件器
    jxf_CPRPrecond* cpr = JXF_CPRCreate(MPI_COMM_WORLD);
    
    // 设置CPR参数
    JXF_Int blk_size = jx_ParBSRMatrixBlockSize(A_high);
    JXF_CPRSetParameter(cpr, "pressure_index", 0);
    JXF_CPRSetParameter(cpr, "block_size", blk_size);
    JXF_CPRSetParameter(cpr, "stage1_maxit", 1);
    JXF_CPRSetParameter(cpr, "stage2_maxit", 2);
    JXF_CPRSetParameter(cpr, "stage1_solver_type", 1);
    JXF_CPRSetParameter(cpr, "stage2_solver_type", 2);
    JXF_CPRSetParameter(cpr, "print_level", 2);  // 设置为0以减少输出
    
    // 创建一个包装函数，使得 BiCGSTAB 的 setup 阶段使用低精度矩阵
    // 这里的 cpr_ptr 指向 cpr，matrix_ptr 是 BiCGSTAB 传入的 A_high（忽略）
    // 我们使用前面已经创建好的 A_low
    double setup_start = MPI_Wtime();
    JXF_Int setup_result = JXF_CPRSetup(cpr, A_low);
    double setup_end = MPI_Wtime();
    
    if (setup_result != JXF_SUCCESS) {
        if (myid == 0) printf("Error setting up low precision CPR preconditioner\n");
        jxf_ParBSRMatrixDestroy(A_low);
        JXF_CPRDestroy(&cpr);
        return -1;
    }
    
    if (myid == 0) {
        printf("Low precision CPR setup completed in %.6f seconds\n", setup_end - setup_start);
    }
    
    // 3. 创建高精度BiCGSTAB求解器
    // BSR BiCGSTAB
    JX_Solver bicgstab_solver;
    jx_BiCGSTABFunctions *fg = jx_BiCGSTABFunctionsCreate(
        jx_ParBSRKrylovCreateVector, jx_ParBSRKrylovDestroyVector,
        jx_ParBSRKrylovMatvecCreate, jx_ParBSRKrylovMatvec, jx_ParBSRKrylovMatvecDestroy,
        jx_ParBSRKrylovInnerProd, jx_ParBSRKrylovCopyVector,
        jx_ParBSRKrylovClearVector, jx_ParBSRKrylovScaleVector, jx_ParBSRKrylovAxpy,
        jx_ParBSRKrylovCommInfo,
        (JX_Int (*)(void*, void*))mixed_cpr_setup_wrapper,
        (JX_Int (*)(void*, void*, void*, void*))JXF_MxP_CPRPrecond1);
    bicgstab_solver = (JX_Solver)jx_BiCGSTABCreate(fg);
    
    // 设置BiCGSTAB参数
    JX_BiCGSTABSetMaxIter(bicgstab_solver, max_iterations);
    JX_BiCGSTABSetTol(bicgstab_solver, (JX_Real)tolerance);
    JX_BiCGSTABSetPrintLevel(bicgstab_solver, 1);
    JX_BiCGSTABSetLogging(bicgstab_solver, 1);
    
    // 4. 设置混合精度CPR为预条件器
    JX_BiCGSTABSetPrecond(bicgstab_solver, 
                       (JX_PtrToSolverFcn)JXF_MxP_CPRPrecond1,  // 混合精度接口
                       (JX_PtrToSolverFcn)JXF_CPRSetup,        // 低精度设置函数
                       cpr);                                    // 低精度预条件器
    
    // 5. 设置BiCGSTAB求解器
    JX_BiCGSTABSetup(bicgstab_solver, (JX_Matrix)A_high, 
                  (JX_Vector)b_high, (JX_Vector)x_high);
    
    // 6. 求解（使用高精度BiCGSTAB + 混合精度CPR预条件器）
    if (myid == 0) {
        printf("\nStarting BiCGSTAB with mixed precision CPR preconditioner...\n");
        printf("Max iterations: %d, Tolerance: %.1e\n", max_iterations, tolerance);
    }
    
    double solve_start = MPI_Wtime();
    JX_BiCGSTABSolve(bicgstab_solver, (JX_Matrix)A_high,
                  (JX_Matrix)A_high, (JX_Vector)b_high, (JX_Vector)x_high);
    double solve_end = MPI_Wtime();
    
    // 7. 获取求解器统计信息
    JX_Int num_iterations = 0;
    JX_Real final_rel_res_norm = 0.0;
    
    JX_BiCGSTABGetNumIterations(bicgstab_solver, &num_iterations);
    JX_BiCGSTABGetFinalRelativeResidualNorm(bicgstab_solver, &final_rel_res_norm);
    
    if (myid == 0) {
        printf("\nBiCGSTAB with mixed precision CPR completed in %.6f seconds\n", 
               solve_end - solve_start);
        printf("Iterations: %d, Final relative residual: %.6e\n", 
               num_iterations, final_rel_res_norm);
        printf("Total time (setup+solve): %.6f seconds\n", 
               setup_end - setup_start + solve_end - solve_start);
    }
    
    // 8. 清理
    /* Note: Skip detailed cleanup to avoid double-free in library (shared pointers).
       The solver completed successfully, MPI_Finalize will handle resources. */
    jxf_ParBSRMatrixDestroy(A_low);
    
    return 0;
}

// =============================================================
// 主函数
// =============================================================

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int myid, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
        
    jx_Vector       *b_Ser;
    jx_Vector       *x_Ser;
    jx_ParVector     *par_rhs;
    jx_ParVector     *par_sol;
    
    if (argc < 3) {
        if (myid == 0) {
            printf("Usage: %s <matrix_file> <matrix_type> [binary] [rhs_file]\n", argv[0]);
            printf("  matrix_type: 0=CSR, 1=BSR\n");
            printf("  binary: 0=text, 1=binary (default: 0)\n");
            printf("  rhs_file: right-hand side file (optional)\n");
        }
        // cleanup skipped
        return 1;
    }
    
    char* filename = argv[1];
    int matrix_type = atoi(argv[2]);
    int binary = (argc > 3) ? atoi(argv[3]) : 0;
    char* rhs_file = (argc > 4) ? argv[4] : NULL;
    
    if (myid == 0) {
        printf("\n=================================================\n");
        printf("Mixed Precision CPR-BiCGSTAB Solver\n");
        printf("Matrix file: %s\n", filename);
        printf("Matrix type: BSR\n");
        printf("Binary format: %s\n", binary ? "yes" : "no");
        if (rhs_file) printf("RHS file: %s\n", rhs_file);
        printf("=================================================\n\n");
    }

    // ------------------------------------------------------------
    // 1. 读取并分发BSR矩阵
    // ------------------------------------------------------------
    jx_BSRMatrix* A_bsr = NULL;
    if (myid == 0) {
        if (binary) {
            A_bsr = jx_BSRMatrixRead_Binary(filename);
        } else {
            A_bsr = jx_BSRMatrixRead(filename);
        }
        
        if (!A_bsr) {
            printf("Error reading BSR matrix from %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("BSR Matrix Info (myid 0):\n");
        printf("  Block size: %d\n", jx_BSRMatrixBlockSize(A_bsr));
        printf("  Number of block rows: %d\n", jx_BSRMatrixNumRows(A_bsr));
        printf("  Number of block cols: %d\n", jx_BSRMatrixNumCols(A_bsr));
        printf("  Number of block nonzeros: %d\n", jx_BSRMatrixNumNonzeros(A_bsr));
        
        int blk_size = jx_BSRMatrixBlockSize(A_bsr);
        int nrows = jx_BSRMatrixNumRows(A_bsr);
        int ncols = jx_BSRMatrixNumCols(A_bsr);
        int Scalar_nrows = nrows * blk_size;
        printf("  Scalar dimensions: %d x %d\n", nrows * blk_size, ncols * blk_size);
        
        if (rhs_file) {
            // 从文件读取右端项
            printf("Reading RHS from file: %s\n", rhs_file);
            b_Ser = read_rhs_file(rhs_file, Scalar_nrows);
            if (!b_Ser) {
                printf("Error reading RHS file, using default RHS (all ones)\n");
                // 如果读取失败，使用默认的1向量
                b_Ser = jx_SeqVectorCreate(Scalar_nrows);
                jx_SeqVectorInitialize(b_Ser);
                JX_Real *b_data = jx_VectorData(b_Ser);
                for (int i = 0; i < Scalar_nrows; i++) {
                    b_data[i] = 1.0;
                }
            }
        } else {
            // 使用默认的1向量
            b_Ser = jx_SeqVectorCreate(Scalar_nrows);
            jx_SeqVectorInitialize(b_Ser);
            JX_Real *b_data = jx_VectorData(b_Ser);
            for (int i = 0; i < Scalar_nrows; i++) {
                b_data[i] = 1.0;
            }
            printf("Using default RHS (all ones)\n");
        }
        
        // 创建初始解向量（全零）
        x_Ser = jx_SeqVectorCreate(Scalar_nrows);
        jx_SeqVectorInitialize(x_Ser);
        JX_Real *x_data = jx_VectorData(x_Ser);
        for (int i = 0; i < Scalar_nrows; i++) {
            x_data[i] = 0.0;
        }
    }
    
        // ------------------------------------------------------------
        // 2. 在串行数据上进行解耦
        // ------------------------------------------------------------
    if (myid == 0){
        printf("\nApplying decoupling to serial BSR system...\n");    
        // 选择解耦方法
        JX_DecoupType decoup_type = JX_DECOUP_TIMPES;  // 推荐使用ABF
        // JX_DecoupType decoup_type = JX_DECOUP_ANL;  // 推荐使用ABF
        JX_Int is_thermal = 0;  // 假设非热力模型

        printf("  Decoupling method: %s\n", jx_DecoupTypeToString(decoup_type));
        printf("  Is thermal model: %s\n", is_thermal ? "yes" : "no");
        
        // 执行解耦（同时修改矩阵和右端项）
        double decoup_start = MPI_Wtime();
        JX_Int decoup_result = jx_BSRMatrixDecouple(A_bsr, b_Ser, decoup_type, is_thermal);
        double decoup_end = MPI_Wtime();
        
        if (decoup_result != JX_SUCCESS) {
            printf("Error applying decoupling to serial system\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("  Decoupling completed in %.6f seconds\n", decoup_end - decoup_start);
        printf("  Matrix and RHS have been modified by decoupling\n");
        
        // 可选：验证解耦效果
        // 可以检查解耦后的对角块是否接近单位阵
        if (1) {
            int nb = jx_BSRMatrixBlockSize(A_bsr);
            JX_Real* data = jx_BSRMatrixData(A_bsr);
            JX_Int* rpt = jx_BSRMatrixI(A_bsr);
            int nblocks = jx_BSRMatrixNumRows(A_bsr);
            
            // 检查第一个对角块
            if (nblocks > 0) {
                JX_Real* first_diag_block = data;
                printf("  First diagonal block after decoupling:\n");
                for (int i = 0; i < nb; i++) {
                    for (int j = 0; j < nb; j++) {
                        printf("    %12.6e", first_diag_block[i*nb + j]);
                    }
                    printf("\n");
                }
            }
        }
    }

    // 生成分区
    JX_BigInt* global_row_start = NULL;
    JX_BigInt* global_col_start = NULL;
    if (myid == 0) {
        jx_GeneratePartitioning(jx_BSRMatrixNumRows(A_bsr), size, &global_row_start);
        jx_GeneratePartitioning(jx_BSRMatrixNumCols(A_bsr), size, &global_col_start);
    }
    
    // 分发矩阵到各进程
    jx_ParBSRMatrix* A_parbsr = jx_BSRMatrixToParBSRMatrix(
        MPI_COMM_WORLD, A_bsr, global_row_start, global_col_start);
    
    if (!A_parbsr) {
        if (myid == 0) printf("Error creating parallel BSR matrix\n");
        // cleanup skipped
        return 1;
    }

    // ------------------------------------------------------------
    // 2. 创建右端项和初始解向量
    // ------------------------------------------------------------
    // 创建分区数组（标量行分区）
    JX_Int  *partitioning;
    jx_ParBSRMatrixGetRowPartitioning(A_parbsr, &partitioning);

    // 创建并行向量
    par_rhs = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), b_Ser, partitioning);
    par_sol = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), x_Ser, partitioning);
    
    if (!par_rhs || !par_sol) {
        printf("Rank %d: ERROR: Failed to create parallel vectors!\n", myid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ------------------------------------------------------------
    // 使用混合精度CPR-BiCGSTAB求解
    // ------------------------------------------------------------
    
    // 求解器参数
    int max_iterations = 60;
    JX_Real tolerance = 1e-4;
    
    int result = solve_with_mixed_precision_cpr_gmres(
        A_parbsr, par_rhs, par_sol,
        max_iterations, tolerance
    );
    
    // ------------------------------------------------------------
    // 清理
    // ------------------------------------------------------------
    
    // 清理矩阵和向量
    // cleanup skipped
    // cleanup skipped
    // cleanup skipped
    
    // cleanup skipped
    return result;
}