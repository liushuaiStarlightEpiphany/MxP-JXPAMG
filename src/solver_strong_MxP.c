//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * solver_strong.c -- This is an example demonstrating how to solve a given
 * linear system in the form of files by calling JXPAMG solver and
 * the JXPAMG-based Krylov iterative methods.
 *
 * One of the following solvers can be used by
 * assigning parameter 'solver_id'
 *  solver_id = 0:  PAMG
 *  solver_id = 11: CG
 *  solver_id = 12: PAMG-CG
 *  solver_id = 13: DS-CG
 *  solver_id = 14: Euclid-CG
 *  solver_id = 21: GMRES
 *  solver_id = 22: PAMG-GMRES
 *  solver_id = 23: DS-GMRES
 *  solver_id = 24: Euclid-GMRES
 *  solver_id = 31: BiCGSTAB
 *  solver_id = 32: PAMG-BiCGSTAB
 *  solver_id = 33: DS-BiCGSTAB
 *  solver_id = 34: Euclid-BiCGSTAB
 *
 * The input data for the linear system to be solved should be
 * provided in the form as follows (the indices are of Fortran
 * style, i.e, starting from 1 not 0):
 *  1. matrix file
 *----------------------------------------
 *     n
 *     ia(i), i = 1(1)n+1
 *     ja(i), i = 1(1)nz
 *      a(i), i = 1(1)nz
 *----------------------------------------
 *  where nz = ia(n+1)-1.
 *  2. rhs file
 *----------------------------------------
 *     n
 *     f(i),i = 1(1)n
 *----------------------------------------
 *
 *  Created by peghoty 2011/05/13
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 *  Modified by Yue Xiaoqiang 2012/10/22
 *
 *  Xiangtan University
 *  yuexq@xtu.edu.cn
 *
 */

#include "jx_pamg.h"
#include "jx_krylov.h"
#include "jx_diagscale.h"
#include "jx_blockprec.h"
#include "jx_euclid.h"
#include "jx_combined.h"

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

#include "jxf_pamg.h"
#include "jxf_krylov.h"
#include "jxf_diagscale.h"
#include "jxf_blockprec.h"
#include "jxf_euclid.h"
#include "jxf_combined.h"

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

JX_Int myid;
JX_Int coreNums = 20;
JX_Int dot_type = 0;
JX_Int jx_spmv_type = 0;
JXF_Int jxf_spmv_type = 0;
// JXF_Int jxh_spmv_type = 0;
// JXF_Int mix_spmv_type = 0;

//++++++++++++++++++++++++++

JX_Int
jx_ParVectorSetRandomValues(jx_ParVector *v, JX_Int seed);
JX_Int
jx_SeqVectorSetRandomValues(jx_Vector *x, JX_Int seed);

/* 低精度预条件 */
JXF_Int
JX_PAMGPrecond_MxP(JXF_Solver solver,
                   JXF_hpCSRMatrix par_matrix,
                   JX_ParVector par_rhs,
                   JX_ParVector par_app);
JX_Int
jx_PAMGPrecond_MxP(void *amg_vdata,
                   jxf_hpCSRMatrix *hp_matrix,
                   jx_ParVector *par_rhs,
                   jx_ParVector *par_app);

/* 高精度矩阵转低精度 */
JXF_Int
jxmp_ParCSRMatrixDtoF(jx_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JX_Int copy_data);
JXF_Int
jxmp_CSRMatrixDtoF(jx_CSRMatrix *A, jxf_CSRMatrix *B, JX_Int copy_data);

/* 高精度向量转低精度 */
JXF_Int
jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y);
JXF_Int
jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y);

/* 低精度向量转高精度 */
JXF_Int
jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y);
JXF_Int
jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y);

JXF_Int
jx_hpPAMGCycle_MxP(void *amg_vdata,
                   jxf_ParVector **F_array,
                   jxf_ParVector **U_array);

int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   // JX_Int myid, nprocs;
   JX_Int nprocs;
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
   JX_Int nthreads;
#endif

   JX_Int arg_index = 0;
   JX_Int print_usage = 0;

   JX_Real starttime, endtime;
   JX_Real starttimeT, endtimeT;

   char *MatFile = NULL;
   char *RhsFile = NULL;
   char *GusFile = NULL;
   char AppFile[120];

   jx_hpCSRMatrix *hp_matrix = NULL;
   jx_ParCSRMatrix *par_matrix = NULL;
   jx_ParVector *par_rhs = NULL;
   jx_ParVector *par_sol = NULL;

   jxf_hpCSRMatrix *hp_matrix_F = NULL;
   jxf_ParCSRMatrix *par_matrix_F = NULL;
   // jxf_ParVector *par_rhs_F = NULL;
   // jxf_ParVector *par_sol_F = NULL;

   JX_Int *partitioning = NULL;

   JX_Int **grid_relax_points = NULL;

   JX_Int build_matrix_arg_index;
   JX_Int build_rhs_arg_index;
   /* DiagScale Precond */
   JX_Solver ds_solver;

   /* Euclid solver */
   JX_Solver euclid_solver;
   JX_Int euclid_level;
   JX_Int euclid_bj;

   /* JXFPAMG solver */
   JXF_Solver amg_solver;
   JXF_Int max_levels;
   JXF_Int cycle_type;
   JXF_Int relax_type;
   JXF_Int measure_type;
   JXF_Int rap2;
   JXF_Int num_functions;
   JXF_Int ns_down;
   JXF_Int ns_up;
   JXF_Int ns_coarse;
   JXF_Int restri_type;
   JXF_Int keepTranspose;
   JXF_Int coarsen_type;
   JXF_Int interp_type;
   JXF_Int P_max_elmts;
   JXF_Int agg_num_levels;
   JXF_Int ai_measure_type;
   JXF_Int ai_relax_type;
   JXF_Real strong_threshold;
   JXF_Real max_row_sum;

   JXF_Real norm1, norm2;

   JXF_Real relax_wt;
   JXF_Real outer_wt;

   JXF_Real S_commpkg_switch;
   JXF_Real AIR_strong_th;

   JXF_Int coarse_threshold;
   JXF_Real coarse_ratio;
   JXF_Int coarsestsolverid;
   JXF_Int conv_criteria;
   JXF_Int amg_print_level;

   /* iterative method */
   JXF_Solver solver;
   JXF_Real tol;
   JXF_Int max_iter;
   JXF_Int k_dim;
   JXF_Int is_check_restarted; /* peghoty, 2011/11/08 */
   JXF_Int twonorm;
   JXF_Int solver_id;
   JXF_Int problem_id;
   JXF_Int file_base;
   JXF_Int print_level;
   JXF_Int keepsol;
   JXF_Int TTest;

   /* other variables */
   JXF_Int glosize, i;
   JXF_Int initguess = 0;
   JXF_Int num_iterations;
   JXF_Real final_res_norm;
   JXF_Real norm;

   //--------------------------
   //  启动 MPI
   //--------------------------
   jx_MPI_Init(&argc, &argv);
   jx_MPI_Comm_rank(comm, &myid);
   jx_MPI_Comm_size(comm, &nprocs);
#ifdef USING_HWLOC
   jx_hpCreateHardwareInfo(comm);
#endif

   //-----------------------
   //  参数设置
   //-----------------------
   build_matrix_arg_index = argc;
   build_rhs_arg_index = argc;
   max_levels = 25;  /* 最大网格层数 */
   cycle_type = 1;   /* Cycle 类型  1: V_Cycle; 2：W_Cycle */
   relax_type = 3;   /* Relax 类型  3: hGS; 6：hSGS */
   measure_type = 0; /* 影响值的计算方式 0：局部；1：全局 */
   rap2 = 0;         /* RAP计算方式  0：RAP；1：先算Q=AP，再算RQ */
   num_functions = 1;
   ns_down = 1;
   ns_up = 1;
   ns_coarse = 1;
   restri_type = 0;         /* 限制算子类型 0                   : P^T, 1 : AIR, 2: AIR-2 */
   keepTranspose = 0;       /* 存放限制算子  0：no；1：yes */
   coarsen_type = 8;        /* 粗化策略 */
   interp_type = 0;         /* 插值策略 */
   P_max_elmts = 0;         /* 插值算子每行最大非零元个数 */
   agg_num_levels = 0;      /* Aggressive粗化的层数 */
   ai_measure_type = 0;     /* AI-策略  0                 : no; 1  : yes */
   ai_relax_type = 0;       /* AI-磨光  0               : no; 1  : yes */
   amg_print_level = 3;     /* work only when AMG as preconditioner */
   strong_threshold = 0.25; /* 强弱连通参数, 0.25 for 2D, 0.5 for 3D is recommended */
   max_row_sum = 0.9;       /* 行和参数 */
   relax_wt = 1.0;
   outer_wt = 1.0;
   S_commpkg_switch = 1.0;
   AIR_strong_th = 0.25;

   coarse_threshold = 1000; /* 最粗网格层上网格节点个数的最大值 */
   coarse_ratio = 0.75;     /* 相邻两个网格层的粗点个数超过细点个数的 coarse_ratio, 则换成 CLJP 粗化 */
   coarsestsolverid = 9;    /* 最粗网格层解法器 */
   conv_criteria = 0;       /* 收敛准则类型 */

   euclid_level = 0; /* level of fill-in */
   euclid_bj = 0;    /* Select PILU (0) or Block Jacobi ILU (1) */

   tol = 1.0e-6;           /* 控制精度 */
   max_iter = 500;         /* 迭代法最大迭代次数 */
   k_dim = 50;             /* 回头数 */
   is_check_restarted = 1; /* peghoty, 2011/11/08 */
   twonorm = 1;            /* PCG 法中的范数控制类型，0   : B 范数; 1: l2 范数 */
   print_level = 3;        /* 0                     : 关闭；1：Setup参数；2：Solve参数；3：Setup+Solve参数 */
   keepsol = 0;            /* 是否保存解向量 */
   TTest = 1;              /* 是否测试时间 */
   solver_id = 22;
   problem_id = -1; /* 默认为-1，需要通过命令行读取矩阵请勿修改此参数*/
   file_base = 0;
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
   nthreads = 1; /* 线程数 */
#endif

   //-----------------------
   //  命令行修改参数
   //-----------------------
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-fromonecsrfile") == 0)
      {
         arg_index++;
         build_matrix_arg_index = arg_index;
      }
      else if (strcmp(argv[arg_index], "-rhsfromfile") == 0)
      {
         arg_index++;
         build_rhs_arg_index = arg_index;
      }
      else if (strcmp(argv[arg_index], "-sid") == 0)
      {
         arg_index++;
         solver_id = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-pid") == 0)
      {
         arg_index++;
         problem_id = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nf") == 0)
      {
         arg_index++;
         num_functions = atoi(argv[arg_index++]);
      }
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
      else if (strcmp(argv[arg_index], "-nts") == 0)
      {
         arg_index++;
         nthreads = atoi(argv[arg_index++]);
      }
#endif
      else if (strcmp(argv[arg_index], "-rap2") == 0)
      {
         arg_index++;
         rap2 = 1;
      }
      else if (strcmp(argv[arg_index], "-air") == 0)
      {
         arg_index++;
         restri_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-kt") == 0)
      {
         arg_index++;
         keepTranspose = 1;
      }
      else if (strcmp(argv[arg_index], "-mxit") == 0)
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-Pmx") == 0)
      {
         arg_index++;
         P_max_elmts = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-agg_nl") == 0)
      {
         arg_index++;
         agg_num_levels = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-tnrm") == 0)
      {
         arg_index++;
         twonorm = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-kdim") == 0)
      {
         arg_index++;
         k_dim = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-euc_lvl") == 0)
      {
         arg_index++;
         euclid_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-euc_bj") == 0)
      {
         arg_index++;
         euclid_bj = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-rlx") == 0)
      {
         arg_index++;
         relax_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ai_rlx") == 0)
      {
         arg_index++;
         ai_relax_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ai_mt") == 0)
      {
         arg_index++;
         ai_measure_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ct") == 0)
      {
         arg_index++;
         coarsen_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ipt") == 0)
      {
         arg_index++;
         interp_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mxct") == 0)
      {
         arg_index++;
         coarse_threshold = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-str") == 0)
      {
         arg_index++;
         strong_threshold = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mxrs") == 0)
      {
         arg_index++;
         max_row_sum = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-scs") == 0)
      {
         arg_index++;
         S_commpkg_switch = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-airst") == 0)
      {
         arg_index++;
         AIR_strong_th = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-amg_ptlv") == 0)
      {
         arg_index++;
         amg_print_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ptlv") == 0)
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-help") == 0)
      {
         print_usage = 1;
         break;
      }
      else
      {
         arg_index++;
      }
   }
   if (print_usage)
   {
      jx_printf("\n");
      jx_printf("  Usage: %s [<options>]\n", argv[0]);
      jx_printf("\n");
      jx_printf("    -bt <val>   : bind type: 0:bind to node, 1:bind to cpu, 2:bind to core\n");
      jx_printf("    -sid <val>   : solver id\n");
      jx_printf("    -pid <val>   : problem id\n");
      jx_printf("  -fromonecsrfile <filename> : ");
      jx_printf("matrix read from a single file (CSR format)\n");
      jx_printf("  -rhsfromonefile        : ");
      jx_printf("rhs read from a single file (CSR format)\n");
      jx_printf("    -nf  <val>   : number of functions\n");
      jx_printf("    -nts <val>   : threads number\n");
      jx_printf("    -rap2        : 2nd implementation of RAP\n");
      jx_printf("    -kt          : keep transpose\n");
      jx_printf("    -Pmx <val>   : maximal number of elements per row for interpolation\n");
      jx_printf(" -agg_nl <val>   : num_levels for aggressive coarsening\n");
      jx_printf("   -tnrm <val>   : two_norm in PCG\n");
      jx_printf("   -kdim <val>   : krylov dimension\n");
      jx_printf("   -dtol <val>   : Drop-tolerance for ILU(0) factorization\n");
      jx_printf("-euc_lvl <val>   : level of fill-in for Euclid\n");
      jx_printf(" -euc_bj <val>   : PILU or Block Jacobi ILU\n");
      jx_printf("  -rlx <val>     : relaxation type\n");
      jx_printf("  -ai_rlx <val>  : AI relaxation type\n");
      jx_printf("  -ai_mt  <val>  : AI measure type\n");
      jx_printf("     -ct <val>   : coarsening type\n");
      jx_printf("    -ipt <val>   : interpolation type\n");
      jx_printf("   -mxct <val>   : max. size on coarsest grid\n");
      jx_printf("    -str <val>   : AMG strength threshold\n");
      jx_printf("   -mxrs <val>   : maximum row sum threshold for dependency weakening\n");
      jx_printf("-amg_ptlv <val>  : print_level of AMG when AMG as preconditioner\n");
      jx_printf("-theta_ms <iiiv> : threshold for multi-scale judgement\n");
      jx_printf("   -ptlv <val>   : print_level\n");
      jx_printf("   -help         : using help message\n\n");
      exit(1);
   }

   //-----------------------
   //  提供线性系统文件
   //-----------------------
   switch (problem_id)
   {
   case 1:
      MatFile = "/vol8/home/xtu_pcy/matrix/Jump/mat_csr_64X64X64_1e3.bin";
      RhsFile = "/vol8/home/xtu_pcy/4/Jump/rhs_64X64X64.bin";
      file_base = 2;
      break;
   case 2:
      MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Anisotropy/mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Anisotropy/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 3:
      MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e3.bin";
      RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 4:
      MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e6.bin";
      RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 5: // 论文测试矩阵
      MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Rand/mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Rand/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 6: // 论文测试矩阵
      MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/test/Checkerboard_mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/test/Checkerboard_rhs_128X128X128.bin";
      file_base = 2;
      break;
   }

   if (build_matrix_arg_index < argc)
   {
      MatFile = argv[build_matrix_arg_index];
      file_base = 0;
   }
   else if (problem_id < 0)
   {
      jx_printf("Error: No filename specified \n");
      exit(1);
   }

   if (build_rhs_arg_index < argc)
   {
      RhsFile = argv[build_rhs_arg_index];
   }
   else if (problem_id < 0)
   {
      jx_printf("Error: No filename specified \n");
      exit(1);
   }

   if (myid == 0)
   {
      jx_printf("\n\n+++++++++++++++++++++\n matrix:\n %s\n+++++++++++++++++++++\n", MatFile);
      jx_printf("\n+++++++++++++++++++++\n rhs:\n %s\n+++++++++++++++++++++\n", RhsFile);
      jx_printf("\n+++++++++++++++++++++");
#if JX_USING_BIG_INT
      jx_printf(" Using BIG_INT,");
#endif
#if JX_USING_BIG_DOUBLE
      jx_printf(" BIG_DOUBLE,");
#endif
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
      jx_printf(" With OpenMP using %d threads,", nthreads);
#endif
      jx_printf(" MPI using %d processors +++++++++++++++++++++\n\n", nprocs);
   }

   //----------------------------------------------------------------
   // 设定线程数
   //----------------------------------------------------------------
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
   omp_set_num_threads(nthreads);
#endif

   if (restri_type) /* Set Restriction to be AIR */
   {
      interp_type = 100; /* 1-pt Interp */
      relax_type = 0;
      ns_down = 0;
      ns_up = 3;
      grid_relax_points = jx_CTAlloc(JX_Int *, 4);
      grid_relax_points[0] = NULL;
      grid_relax_points[1] = jx_CTAlloc(JX_Int, ns_down);
      grid_relax_points[2] = jx_CTAlloc(JX_Int, ns_up);
      grid_relax_points[3] = jx_CTAlloc(JX_Int, ns_coarse);
      for (i = 0; i < ns_down; i++)
         grid_relax_points[1][i] = 0; /* down cycle */
      if (ns_up == 3)                 /* up cycle */
      {
         grid_relax_points[2][0] = -1; // F
         grid_relax_points[2][1] = -1; // F
         grid_relax_points[2][2] = 1;  // C
      }
      else if (ns_up == 2)
      {
         grid_relax_points[2][0] = -1;
         grid_relax_points[2][1] = -1;
      }
      for (i = 0; i < ns_coarse; i++)
         grid_relax_points[3][i] = 0; /* coarse: all */
      coarse_threshold = 100;
      agg_num_levels = 0; /* does not support aggressive coarsening */
   }

   //--------------------------------------------------------------------
   //  利用文件创建并行矩阵和并行右端, 以及并行初始迭代向量(零向量或随机向量)
   //--------------------------------------------------------------------
   if (TTest)
      starttime = jx_MPI_Wtime();

   hp_matrix = jx_hpBuildMatParFromOneFile(MatFile, 1, file_base);
   par_rhs = jx_hpBuildRhsParFromOneFile(RhsFile, hp_matrix, file_base);

   if (!GusFile)
   {
      glosize = jx_ParVectorGlobalSize(par_rhs);
      partitioning = jx_ParVectorPartitioning(par_rhs);
      par_sol = jx_ParVectorCreate(comm, glosize, partitioning);
      jx_ParVectorSetPartitioningOwner(par_sol, 0);
      jx_ParVectorInitialize(par_sol);
      if (initguess == 0)
      {
         jx_ParVectorSetConstantValues(par_sol, 0.0);
      }
      else
      {
         jx_ParVectorSetRandomValues(par_sol, 22775);
         norm = jx_ParVectorInnerProd(par_sol, par_sol);
         norm = 1. / sqrt(norm);
         jx_ParVectorScale(norm, par_sol);
      }
   }
   else
   {
      par_sol = jx_hpBuildRhsParFromOneFile(GusFile, hp_matrix, 0);
   }

   //--------------------------------------------------------------------
   //  将并行矩阵和并行右端以及并行初始迭代向量转成低精度版本
   //--------------------------------------------------------------------
   if (myid == 0)
      jx_printf(" 将并行矩阵和并行右端以及并行初始迭代向量转成低精度版本+++++ \n\n");
   par_matrix = jx_hpCSRMatrixPar(hp_matrix);
   par_matrix_F = jxf_ParCSRMatrixCreate(comm,
                                         jx_ParCSRMatrixGlobalNumRows(par_matrix), // 矩阵 A 的全局行数（所有进程上的总行数）
                                         jx_ParCSRMatrixGlobalNumRows(par_matrix), // 矩阵 A 的全局行数（所有进程上的总行数）
                                         jx_ParCSRMatrixRowStarts(par_matrix),
                                         jx_ParCSRMatrixColStarts(par_matrix),
                                         jx_CSRMatrixNumCols(jx_ParCSRMatrixOffd(par_matrix)),
                                         jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixDiag(par_matrix)),
                                         jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixOffd(par_matrix)));

   jxf_ParCSRMatrixInitialize(par_matrix_F);
   jxmp_ParCSRMatrixDtoF(par_matrix, par_matrix_F, 1);
   hp_matrix_F = jxf_hpInithpCSRMatrix();
   jxf_hpCSRMatrixPar(hp_matrix_F) = par_matrix_F;

   endtime = jx_MPI_Wtime();
   jx_GetWallTime(comm, "BuildParLinearSystem", starttime, endtime, 0, 3);

   // // 待补充
   // par_rhs_F = jxf_ParVectorCreate(comm, glosize, partitioning);

   //--------------------------------------------------------------------
   //  求解线性代数系统
   //--------------------------------------------------------------------
   starttimeT = jx_MPI_Wtime();

   switch (solver_id)
   {
   case 22: /* PAMG-GMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: PAMG-GMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver, max_levels);
      JXF_PAMGSetMaxIter(amg_solver, 1);
      JXF_PAMGSetCycleType(amg_solver, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver, measure_type);
      JXF_PAMGSetRAP2(amg_solver, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);                          /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);                         /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);                        /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);                        /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3); /* relax_type for "coarsest" */ /* relax_type for "coarsest" */

      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                         (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix_F);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Setup", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                    (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Solve", starttime, endtime, 0, 3);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver);
      JX_GMRESDestroy(solver);
   }
   break;

   case 32: /* PAMG-BiCGSTAB */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: PAMG-BiCGSTAB \n\n");

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver, max_levels);
      JXF_PAMGSetMaxIter(amg_solver, 1);
      JXF_PAMGSetCycleType(amg_solver, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver, measure_type);
      JXF_PAMGSetRAP2(amg_solver, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);

      JX_BiCGSTABCreate(comm, &solver);
      JX_BiCGSTABSetMaxIter(solver, max_iter);
      JX_BiCGSTABSetTol(solver, tol);
      JX_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JX_BiCGSTABSetConvCriteria(solver, 0);
      JX_BiCGSTABSetLogging(solver, 1);
      JX_BiCGSTABSetPrintLevel(solver, print_level);

      JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                            (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix_F);

      JX_BiCGSTABSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-BiCGSTAB Setup", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_BiCGSTABSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                       (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-BiCGSTAB Solve", starttime, endtime, 0, 3);
      }

      JX_BiCGSTABGetNumIterations(solver, &num_iterations);
      JX_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver);
      JX_BiCGSTABDestroy(solver);
   }
   break;

   case 42: /* PAMG-FlexGMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: MxP PAMG-FlexGMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver, max_levels);
      JXF_PAMGSetMaxIter(amg_solver, 1);
      JXF_PAMGSetCycleType(amg_solver, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver, measure_type);
      JXF_PAMGSetRAP2(amg_solver, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      JX_ParCSRFlexGMRESCreate(comm, &solver);
      JX_FlexGMRESSetKDim(solver, k_dim);
      JX_FlexGMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_FlexGMRESSetMaxIter(solver, max_iter);
      JX_FlexGMRESSetTol(solver, tol);
      JX_FlexGMRESSetLogging(solver, 1);
      JX_FlexGMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_FlexGMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                             (JX_PtrToSolverFcn)JXF_PAMGSetup, amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix_F);

      /* this is optional - could be a user defined one instead */
      JX_FlexGMRESSetModifyPC(solver, (JX_PtrToModifyPCFcn)jx_FlexGMRESModifyPCDefault);

      JX_FlexGMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Setup", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_FlexGMRESSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                        (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Solve", starttime, endtime, 0, 3);
      }

      JX_FlexGMRESGetNumIterations(solver, &num_iterations);
      JX_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver);
      JX_ParCSRFlexGMRESDestroy(solver);
   }
   break;
   }

   if (TTest)
   {
      endtimeT = jx_MPI_Wtime();
      jx_GetWallTime(comm, "Total Solve Time", starttimeT, endtimeT, 0, 3);
   }

   // rel-residual
   norm2 = jx_ParVectorNorm2(par_rhs);
   jx_hpCSRMatrixMatvec(1.0, hp_matrix, par_sol, -1.0, par_rhs);
   norm1 = jx_ParVectorNorm2(par_rhs);
   norm1 = norm1 / norm2;

   if (myid == 0)
   {
      jx_printf("\n true rel residual is %.4le\n", norm1);
   }

   //-----------------------------------------------------------
   //  将近似解向量保存到文件中
   //-----------------------------------------------------------
   if (keepsol)
   {
      jx_Vector *ser_sol = NULL;
      ser_sol = jx_ParVectorToVectorAll(par_sol);
      if (myid == 0)
      {
         jx_sprintf(AppFile, "%s%06d%s%06d%s%03d%s%03d",
                    "app_", nprocs, "_", problem_id, "_", solver_id, "_", relax_type);
         jx_SeqVectorPrint(ser_sol, AppFile);
      }
      jx_SeqVectorDestroy(ser_sol);
   }

//-----------------------------------------------------------
//  销毁并行矩阵和并行向量
//-----------------------------------------------------------
#ifdef USING_HWLOC
   jx_hpCSRhardwareDestroy();
#endif
   jx_hpCSRMatrixDestroy(hp_matrix);
   // jx_ParVectorDestroy(par_rhs);
   // jx_ParVectorDestroy(par_sol);

   //-----------------------------------------------------------
   //  终止 MPI
   //-----------------------------------------------------------
   jx_MPI_Finalize();

   return 0;
}

JXF_Int
JX_PAMGPrecond_MxP(JXF_Solver solver,
                   JXF_hpCSRMatrix par_matrix,
                   JX_ParVector par_rhs,
                   JX_ParVector par_app)
{

   return (jx_PAMGPrecond_MxP((void *)solver,
                              (jxf_hpCSRMatrix *)par_matrix,
                              (jx_ParVector *)par_rhs,
                              (jx_ParVector *)par_app));
}

JXF_Int
jx_PAMGPrecond_MxP(void *amg_vdata,
                   jxf_hpCSRMatrix *hp_matrix,
                   jx_ParVector *par_rhs,
                   jx_ParVector *par_app)
{
   // jx_printf("\n >>> jx_PAMGPrecond_MxP \n");
   jxf_ParAMGData *amg_data = amg_vdata;
   jxf_ParVector **F_array = jx_ParAMGDataFArray(amg_data);
   jxf_ParVector **U_array = jx_ParAMGDataUArray(amg_data);

   JXF_Int max_iter = jx_ParAMGDataMaxIter(amg_data);
   JXF_Int max_levels = jx_ParAMGDataMaxLevels(amg_data);
   JXF_Int cycle_count = 0;

   jxf_ParVector *par_rhs_F = NULL;
   jxf_ParVector *par_app_F = NULL;

   par_rhs_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
                                   jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                                   jxf_hpCSRMatrixRowStarts(hp_matrix));
   jxf_ParVectorInitialize(par_rhs_F);
   jxmp_ParVectorDtoF(par_rhs, par_rhs_F);
   par_app_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
                                   jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
                                   jxf_hpCSRMatrixRowStarts(hp_matrix));
   jxf_ParVectorInitialize(par_app_F);
   jxmp_ParVectorDtoF(par_app, par_app_F);

   if (max_levels == 1)
   {
      jxf_CoarsestSolver(amg_data, hp_matrix, par_rhs_F, par_app_F);
      return (0);
   }

   F_array[0] = par_rhs_F;
   U_array[0] = par_app_F;

   /*----------------------------------------------------------------
    * Generally, max_iter should be set to 1 when AMG is used as
    * a preconditioner, here, we keep the parameter max_iter for
    * convenience if max_iter should be larger than 1.
    *---------------------------------------------------------------*/
   while (cycle_count < max_iter)
   {
      jxf_hpPAMGCycle(amg_data, F_array, U_array);
      // jx_hpPAMGCycle_MxP(amg_data, F_array, U_array);
      cycle_count++;
   }

   jxmp_ParVectorFtoD(U_array[0], par_app);

   return 0;
}

/*!
 * \fn JX_Int jx_ParVectorSetRandomValues
 */
JX_Int
jx_ParVectorSetRandomValues(jx_ParVector *v, JX_Int seed)
{
   JX_Int my_id;
   jx_Vector *v_local = jx_ParVectorLocalVector(v);

   MPI_Comm comm = jx_ParVectorComm(v);
   jx_MPI_Comm_rank(comm, &my_id);

   seed *= (my_id + 1);

   return jx_SeqVectorSetRandomValues(v_local, seed);
}

/*!
 * \fn JX_Int jx_SeqVectorSetRandomValues
 */
JX_Int
jx_SeqVectorSetRandomValues(jx_Vector *x, JX_Int seed)
{
   JX_Real *vector_data = jx_VectorData(x);
   JX_Int size = jx_VectorSize(x);
   JX_Int i, ierr = 0;

   jx_SeedRand(seed);

   size *= jx_VectorNumVectors(x);

   /* RDF: threading this loop may cause problems because of jx_Rand() */
   for (i = 0; i < size; i++)
   {
      vector_data[i] = 2.0 * jx_Rand() - 1.0;
   }

   return ierr;
}

/* 高精度矩阵转低精度 */
JXF_Int
jxmp_CSRMatrixDtoF(jx_CSRMatrix *A, jxf_CSRMatrix *B, JX_Int copy_data)
{
   jx_printf("\n >>> jxmp_CSRMatrixDtoF \n");
   JX_Int ierr = 0;
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Real *A_data;

   JX_Int *B_i = jxf_CSRMatrixI(B);
   JX_Int *B_j = jxf_CSRMatrixJ(B);
   JXF_Real *B_data;
   JX_Int i, j;

   for (i = 0; i < num_rows; i++)
   {
      B_i[i] = A_i[i];
      for (j = A_i[i]; j < A_i[i + 1]; j++)
      {
         B_j[j] = A_j[j];
      }
   }
   B_i[num_rows] = A_i[num_rows];

   if (copy_data)
   {
      A_data = jx_CSRMatrixData(A);
      B_data = jxf_CSRMatrixData(B);

      for (i = 0; i < num_rows; i++)
      {
         for (j = A_i[i]; j < A_i[i + 1]; j++)
         {

            B_data[j] = (JXF_Real)A_data[j];
         }
      }
   }

   return ierr;
}

JXF_Int
jxmp_ParCSRMatrixDtoF(jx_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JX_Int copy_data)
{
   jx_printf("\n >>> jxmp_ParCSRMatrixDtoF \n");
   jx_CSRMatrix *A_diag;
   jx_CSRMatrix *A_offd;
   JX_Int *col_map_offd_A; // A 非对角块的列映射表
   JX_Int num_cols_offd;
   jxf_CSRMatrix *B_diag;
   jxf_CSRMatrix *B_offd;
   JXF_Int *col_map_offd_B;
   JXF_Int i;

   A_diag = jx_ParCSRMatrixDiag(A);
   A_offd = jx_ParCSRMatrixOffd(A);
   col_map_offd_A = jx_ParCSRMatrixColMapOffd(A);
   num_cols_offd = jx_CSRMatrixNumCols(A_offd);

   B_diag = jxf_ParCSRMatrixDiag(B);
   B_offd = jxf_ParCSRMatrixOffd(B);
   col_map_offd_B = jxf_ParCSRMatrixColMapOffd(B);

   jxmp_CSRMatrixDtoF(A_diag, B_diag, copy_data);

   jxmp_CSRMatrixDtoF(A_offd, B_offd, copy_data);

   if (num_cols_offd && col_map_offd_B == NULL)
   {
      col_map_offd_B = jx_CTAlloc(JX_Int, num_cols_offd);
      jxf_ParCSRMatrixColMapOffd(B) = col_map_offd_B;
   }
   for (i = 0; i < num_cols_offd; i++)
   {
      col_map_offd_B[i] = col_map_offd_A[i];
   }

   return jx_error_flag;
}

/* 高精度向量转低精度 */
JXF_Int
jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y)
{
   jx_Vector *x_local = jx_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);

   return jxmp_SeqVectorDtoF(x_local, y_local);
}

JXF_Int
jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y)
{
   JX_Real *x_data = jx_VectorData(x);
   JX_Int size = jx_VectorSize(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JX_Int i;

   size *= jx_VectorNumVectors(x);

#define JX_SMP_PRIVATE i
#include "../jxpamg-basic-D/include/jx_smp_forloop.h"
   for (i = 0; i < size; i++)
   {
      y_data[i] = (JXF_Real)x_data[i]; // 双精度到单精度转换
   }

   return 0;
}

/* 低精度向量转高精度 */
JXF_Int
jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y)
{
   // jx_printf("\n >>> jxmp_ParVectorFtoD \n");
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jx_Vector *y_local = jx_ParVectorLocalVector(y);

   return jxmp_SeqVectorFtoD(x_local, y_local);
}

JXF_Int
jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y)
{
   // jx_printf("\n >>> jxmp_SeqVectorFtoD \n");
   JXF_Real *x_data = jxf_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int size = jxf_VectorSize(x);
   JX_Int i = 0;

   size *= jxf_VectorNumVectors(x);

#define JX_SMP_PRIVATE i
#include "../jxpamg-basic-D/include/jx_smp_forloop.h"
   for (i = 0; i < size; i++)
   {
      y_data[i] = (JX_Real)x_data[i]; // 双精度到单精度转换
   }

   return 0;
}