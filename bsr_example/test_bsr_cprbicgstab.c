// test_cpr_gmres.c - CPR-BiCGSTAB方法测试程序
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#include "jx_mv.h"
#include "jx_bsr_mv.h"
#include "jx_parbsr_mv.h"
#include "jx_cpr.h"
#include "jx_bsr_decoup.h"
#include "jx_pamg.h"
#include "jx_krylov.h"

typedef double Real;
typedef int Int;

jx_Vector* read_rhs_file(const char* filename, int expected_size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening RHS file: %s\n", filename);
        return NULL;
    }
    int vector_size;
    if (fscanf(file, "%d", &vector_size) != 1) {
        printf("Error reading vector size from RHS file\n");
        fclose(file);
        return NULL;
    }
    if (expected_size != -1 && vector_size != expected_size) {
        printf("Warning: RHS vector size (%d) does not match expected size (%d)\n",
               vector_size, expected_size);
    }
    jx_Vector* rhs = jx_SeqVectorCreate(vector_size);
    jx_SeqVectorInitialize(rhs);
    JX_Real* rhs_data = jx_VectorData(rhs);
    for (int i = 0; i < vector_size; i++) {
        fscanf(file, "%lf", &rhs_data[i]);
    }
    fclose(file);
    return rhs;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int myid, size;
    double total_start;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (myid == 0) {
        printf("  running cpr-bicgstab with %d MPI-processes\n", size);
    }

    jx_Vector *b_Ser = NULL, *x_Ser = NULL;
    jx_ParVector *par_rhs = NULL, *par_sol = NULL;

    if (argc < 3) {
        if (myid == 0) {
            printf("Usage: %s <matrix_file> <matrix_type> [binary] [rhs_file] [decoup_type]\n", argv[0]);
            printf("  matrix_type: 0=CSR, 1=BSR\n");
            printf("  binary: 0=text, 1=binary (default: 0)\n");
            printf("  rhs_file: right-hand side file (optional)\n");
            printf("  decoup_type: 0=None 1=ABF 2=ANL 3=SEM 4=QI 5=TIMPES 6=TIMPES2 (default: 1)\n");
        }
        MPI_Finalize();
        return 1;
    }

    char* filename = argv[1];
    int matrix_type = atoi(argv[2]);
    int binary = (argc > 3) ? atoi(argv[3]) : 0;
    char* rhs_file = (argc > 4 && argv[4][0] != '\0') ? argv[4] : NULL;
    JX_DecoupType decoup_type = (argc > 5) ? (JX_DecoupType)atoi(argv[5]) : JX_DECOUP_TIMPES;

    if (matrix_type != 1) {
        if (myid == 0) printf("CSR matrix testing (not supported in this test)\n");
        MPI_Finalize();
        return 1;
    }

    if (myid == 0) {
        printf("\n=================================================\n");
        printf("Testing CPR-BiCGSTAB for BSR Matrix\n");
        printf("Matrix file: %s\n", filename);
        printf("Decoupling: %s (%d)\n", jx_DecoupTypeToString(decoup_type), decoup_type);
        printf("=================================================\n\n");
    }

    // ------------------------------------------------------------
    // 1. 读取BSR矩阵
    // ------------------------------------------------------------
    jx_BSRMatrix* A_bsr = NULL;
    if (myid == 0) {
        if (binary) A_bsr = jx_BSRMatrixRead_Binary(filename);
        else        A_bsr = jx_BSRMatrixRead(filename);

        if (!A_bsr) {
            printf("Error reading BSR matrix from %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("BSR Matrix Info:\n");
        printf("  Block size: %d\n", jx_BSRMatrixBlockSize(A_bsr));
        printf("  Block rows: %d\n", jx_BSRMatrixNumRows(A_bsr));
        printf("  Block cols: %d\n", jx_BSRMatrixNumCols(A_bsr));
        printf("  Block nonzeros: %d\n", jx_BSRMatrixNumNonzeros(A_bsr));
        int blk_size = jx_BSRMatrixBlockSize(A_bsr);
        int nrows = jx_BSRMatrixNumRows(A_bsr);
        printf("  Scalar dimensions: %d x %d\n", nrows * blk_size, jx_BSRMatrixNumCols(A_bsr) * blk_size);

        int Scalar_nrows = nrows * blk_size;
        if (rhs_file) {
            b_Ser = read_rhs_file(rhs_file, Scalar_nrows);
            if (!b_Ser) {
                printf("Error reading RHS file, using default\n");
                b_Ser = jx_SeqVectorCreate(Scalar_nrows);
                jx_SeqVectorInitialize(b_Ser);
                JX_Real *d = jx_VectorData(b_Ser);
                for (int i = 0; i < Scalar_nrows; i++) d[i] = 1.0;
            }
        } else {
            b_Ser = jx_SeqVectorCreate(Scalar_nrows);
            jx_SeqVectorInitialize(b_Ser);
            JX_Real *d = jx_VectorData(b_Ser);
            for (int i = 0; i < Scalar_nrows; i++) d[i] = 1.0;
        }

        x_Ser = jx_SeqVectorCreate(Scalar_nrows);
        jx_SeqVectorInitialize(x_Ser);
        JX_Real *xd = jx_VectorData(x_Ser);
        for (int i = 0; i < Scalar_nrows; i++) xd[i] = 0.0;
    }

    // ------------------------------------------------------------
    // 2. 执行解耦
    // ------------------------------------------------------------
    JX_Int is_thermal = 0;
    if (myid == 0) {
        printf("\nApplying decoupling...\n");
        printf("  Method: %s\n", jx_DecoupTypeToString(decoup_type));

        total_start = MPI_Wtime();
        double decoup_start = MPI_Wtime();
        JX_Int decoup_result = jx_BSRMatrixDecouple(A_bsr, b_Ser, decoup_type, is_thermal);
        double decoup_end = MPI_Wtime();

        if (decoup_result != JX_SUCCESS) {
            printf("  Error applying decoupling\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("  Completed in %.6f seconds\n", decoup_end - decoup_start);

        if (decoup_type != JX_DECOUP_NONE) {
            int nb = jx_BSRMatrixBlockSize(A_bsr);
            JX_Real* data = jx_BSRMatrixData(A_bsr);
            int nblocks = jx_BSRMatrixNumRows(A_bsr);
            if (nblocks > 0) {
                printf("  First diagonal block after decoupling:\n");
                for (int i = 0; i < nb; i++) {
                    printf("   ");
                    for (int j = 0; j < nb; j++)
                        printf(" %12.6e", data[i*nb + j]);
                    printf("\n");
                }
            }
        }
        fflush(stdout);
    }

    // ------------------------------------------------------------
    // 3. 生成分区并分发矩阵
    // ------------------------------------------------------------
    JX_BigInt* global_row_start = NULL;
    JX_BigInt* global_col_start = NULL;
    if (myid == 0) {
        jx_GeneratePartitioning(jx_BSRMatrixNumRows(A_bsr), size, &global_row_start);
        jx_GeneratePartitioning(jx_BSRMatrixNumCols(A_bsr), size, &global_col_start);
    }

    jx_ParBSRMatrix* A_parbsr = jx_BSRMatrixToParBSRMatrix(
        MPI_COMM_WORLD, A_bsr, global_row_start, global_col_start);

    if (!A_parbsr) {
        if (myid == 0) printf("Error creating parallel BSR matrix\n");
        MPI_Finalize();
        return 1;
    }

    int blk_size = jx_ParBSRMatrixBlockSize(A_parbsr);
    int global_scalar_size = jx_ParBSRMatrixGlobalNumRows(A_parbsr) * blk_size;

    JX_Int *partitioning;
    jx_ParBSRMatrixGetRowPartitioning(A_parbsr, &partitioning);

    par_rhs = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), b_Ser, partitioning);
    par_sol = jx_VectorToParVector(jx_ParBSRMatrixComm(A_parbsr), x_Ser, partitioning);

    if (!par_rhs || !par_sol) {
        printf("Rank %d: ERROR: Failed to create parallel vectors!\n", myid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // ------------------------------------------------------------
    // 4. 创建CPR预条件器
    // ------------------------------------------------------------
    if (myid == 0) {
        printf("\nCreating CPR preconditioner...\n");
        fflush(stdout);
    }

    jx_CPRPrecond* cpr = JX_CPRCreate(MPI_COMM_WORLD);
    JX_CPRSetParameter(cpr, "pressure_index", 0);
    JX_CPRSetParameter(cpr, "block_size", blk_size);
    JX_CPRSetParameter(cpr, "stage1_maxit", 1);
    JX_CPRSetParameter(cpr, "stage1_solver_type", 2);
    JX_CPRSetParameter(cpr, "stage2_maxit", 2);
    JX_CPRSetParameter(cpr, "stage2_solver_type", 2);
    JX_CPRSetParameter(cpr, "print_level", 1);
    JX_CPRSetRealParameter(cpr, "threshold", 1e-15);

    double setup_start = MPI_Wtime();
    JX_Int setup_result = JX_CPRSetup(cpr, A_parbsr);
    double setup_end = MPI_Wtime();

    if (setup_result != JX_SUCCESS) {
        if (myid == 0) printf("Error setting up CPR preconditioner\n");
        MPI_Finalize();
        return 1;
    }
    if (myid == 0) printf("CPR setup completed in %.6f seconds\n", setup_end - setup_start);

    // ------------------------------------------------------------
    // 5. 创建GMRES求解器
    // ------------------------------------------------------------
    if (myid == 0) printf("\nCreating BiCGSTAB solver with CPR preconditioner...\n");
    fflush(stdout);

    int max_iter = 200;
    Real tol = 1e-4;

    // 创建BSR版BiCGSTAB求解器
    jx_BiCGSTABFunctions *bicgstab_functions = jx_BiCGSTABFunctionsCreate(
        jx_ParBSRKrylovCreateVector, jx_ParBSRKrylovDestroyVector,
        jx_ParBSRKrylovMatvecCreate, jx_ParBSRKrylovMatvec, jx_ParBSRKrylovMatvecDestroy,
        jx_ParBSRKrylovInnerProd, jx_ParBSRKrylovCopyVector,
        jx_ParBSRKrylovClearVector, jx_ParBSRKrylovScaleVector, jx_ParBSRKrylovAxpy,
        jx_ParBSRKrylovCommInfo,
        (JX_Int (*)(void*, void*))JX_CPRSetup,
        (JX_Int (*)(void*, void*, void*, void*))JX_CPRPrecond);
    JX_Solver bicgstab_solver = (JX_Solver)jx_BiCGSTABCreate(bicgstab_functions);
    JX_BiCGSTABSetMaxIter(bicgstab_solver, max_iter);
    JX_BiCGSTABSetTol(bicgstab_solver, tol);
    JX_BiCGSTABSetLogging(bicgstab_solver, 1);
    JX_BiCGSTABSetPrintLevel(bicgstab_solver, 1);

    JX_BiCGSTABSetPrecond(bicgstab_solver,
                         (JX_PtrToSolverFcn)JX_CPRPrecond,
                         (JX_PtrToSolverFcn)JX_CPRSetup,
                         (JX_Solver)cpr);

    double bicgstab_setup_start = MPI_Wtime();
    JX_BiCGSTABSetup(bicgstab_solver, (JX_Matrix)A_parbsr,
                    (JX_Vector)par_rhs, (JX_Vector)par_sol);
    double bicgstab_setup_end = MPI_Wtime();

    if (myid == 0) {
        printf("BiCGSTAB setup completed in %.6f seconds\n", bicgstab_setup_end - bicgstab_setup_start);
        fflush(stdout);
    }

    // ------------------------------------------------------------
    // 6. 求解
    // ------------------------------------------------------------
    if (myid == 0) {
        printf("\nSolving with CPR-BiCGSTAB (max_iter=%d, tol=%.1e)...\n", max_iter, tol);
        fflush(stdout);
    }

    jx_ParVector* residual = jx_ParVectorCreate(MPI_COMM_WORLD, global_scalar_size, partitioning);
    jx_ParVectorInitialize(residual);

    jx_ParVectorCopy(par_rhs, residual);
    jx_ParBSRMatrixMatvec(-1.0, A_parbsr, par_sol, 1.0, residual);
    Real init_res_norm = jx_ParVectorNorm2(residual);
    if (myid == 0) printf("Initial residual norm: %.6e\n", init_res_norm);

    double solve_start = MPI_Wtime();
    JX_BiCGSTABSolve(bicgstab_solver, (JX_Matrix)A_parbsr,
                    (JX_Matrix)A_parbsr, (JX_Vector)par_rhs, (JX_Vector)par_sol);
    double solve_end = MPI_Wtime();

    jx_ParVectorCopy(par_rhs, residual);
    jx_ParBSRMatrixMatvec(-1.0, A_parbsr, par_sol, 1.0, residual);
    Real final_res_norm = jx_ParVectorNorm2(residual);

    JX_Int num_iterations = 0;
    JX_Real final_rel_res_norm = 0.0;
    JX_BiCGSTABGetNumIterations(bicgstab_solver, &num_iterations);
    JX_BiCGSTABGetFinalRelativeResidualNorm(bicgstab_solver, &final_rel_res_norm);

    Real rhs_norm = jx_ParVectorNorm2(par_rhs);

    if (myid == 0) {
        Real rel_res = final_res_norm / (rhs_norm + 1e-16);
        printf("\n========= Results [%s] =========\n", jx_DecoupTypeToString(decoup_type));
        printf("  Iterations: %d\n", num_iterations);
        printf("  Final residual: %.6e (rel: %.6e)\n", final_res_norm, final_rel_res_norm);
        printf("  Relative residual (wrt rhs): %.6e\n", rel_res);
        printf("  Converged: %s\n", rel_res < tol ? "YES" : "NO");
        printf("  Solve time: %.6f s\n", solve_end - solve_start);
        printf("  Total time (setup+solve): %.6f s\n", solve_end - total_start);
        printf("================================\n\n");
    }

    MPI_Finalize();
    return 0;
}
