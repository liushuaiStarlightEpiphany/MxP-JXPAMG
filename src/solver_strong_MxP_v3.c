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

#include "adp_cycle.h"

// JX_Int
// jx_ParVectorSetRandomValues(jx_ParVector *v, JX_Int seed);
// JX_Int
// jx_SeqVectorSetRandomValues(jx_Vector *x, JX_Int seed);

// /* 低精度预条件 */
// JXF_Int
// JX_PAMGPrecond_MxP(JXF_Solver solver, jxf_hpCSRMatrix par_matrix, JX_ParVector par_rhs, JX_ParVector par_app);
// JX_Int
// jx_PAMGPrecond_MxP(void *amg_vdata, jxf_hpCSRMatrix *hp_matrix, jx_ParVector *par_rhs, jx_ParVector *par_app);

// /* 高精度矩阵转低精度 */
// JXF_Int
// jxmp_CSRMatrixDtoF(jx_CSRMatrix *A, jxf_CSRMatrix *B, JX_Int copy_data);
// JXF_Int
// jxmp_ParCSRMatrixDtoF(jx_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JX_Int copy_data);

// /* 高精度向量转低精度 */
// JXF_Int
// jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y);
// JXF_Int
// jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y);

// /* 低精度向量转高精度 */
// JXF_Int
// jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y);
// JXF_Int
// jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y);

// JXF_Int
// JXF_PAMGSetup_DtoF(void *amg_solver, void *amg_solver_F, jxf_hpCSRMatrix *hp_matrix_F);

// ===================================================
int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   JX_Int myid, nprocs;
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
   jxf_ParVector *par_rhs_F = NULL;
   jxf_ParVector *par_sol_F = NULL;

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
   /* -----------------------------------------------------------------
    * 双精度版本 (Double Precision) - 使用 JX_ 前缀
    * 对应 JX_Real (double)
    * ----------------------------------------------------------------- */
   JX_Solver amg_solver;
   JX_Real strong_threshold;
   JX_Real max_row_sum;
   JX_Real relax_wt;
   JX_Real outer_wt;
   JX_Real AIR_strong_th;
   JX_Real norm1, norm2;     /* 用于计算残差或模量 */
   JX_Real coarse_ratio;     /* 粗网格比例 */
   JX_Real S_commpkg_switch; /* 通信包切换阈值 */

   /* -----------------------------------------------------------------
    * 单精度版本 (Single Precision) - 使用 JXF_ 前缀
    * 对应 JXF_Real (float)
    * ----------------------------------------------------------------- */
   JXF_Solver amg_solver_F;
   JXF_Real strong_threshold_F;
   JXF_Real max_row_sum_F;
   JXF_Real relax_wt_F;
   JXF_Real outer_wt_F;
   JXF_Real AIR_strong_th_F;
   JXF_Real norm1_F, norm2_F;
   JXF_Real coarse_ratio_F;
   JXF_Real S_commpkg_switch_F;

   /* -----------------------------------------------------------------
    * 公用整型变量 (Int is Int)
    * JX_Int 和 JXF_Int 均为 int，内存长度一致，可以共用以节省空间
    * ----------------------------------------------------------------- */
   JX_Int max_levels;
   JX_Int cycle_type;
   JX_Int relax_type;
   JX_Int measure_type;
   JX_Int rap2;
   JX_Int num_functions;
   JX_Int ns_down;
   JX_Int ns_up;
   JX_Int ns_coarse;
   JX_Int restri_type;
   JX_Int keepTranspose;
   JX_Int coarsen_type;
   JX_Int interp_type;
   JX_Int P_max_elmts;
   JX_Int agg_num_levels;
   JX_Int ai_measure_type;
   JX_Int ai_relax_type;
   JX_Int coarse_threshold;
   JX_Int coarsestsolverid;
   JX_Int conv_criteria;
   JX_Int amg_print_level;

   /* iterative method */
   JX_Solver solver;
   JX_Real tol;
   JX_Int max_iter;
   JX_Int k_dim;
   JX_Int is_check_restarted; /* peghoty, 2011/11/08 */
   JX_Int twonorm;
   JX_Int solver_id;
   JX_Int problem_id;
   JX_Int file_base;
   JX_Int print_level;
   JX_Int keepsol;
   JX_Int TTest;

   /* other variables */
   JX_Int glosize, i;
   JX_Int initguess = 0;
   JX_Int num_iterations;
   JX_Real final_res_norm;
   JX_Real norm;
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
   is_check_restarted = 0; /* peghoty, 2011/11/08 */
   twonorm = 1;            /* PCG 法中的范数控制类型，0   : B 范数; 1: l2 范数 */
   print_level = 3;        /* 0 : 关闭；1：Setup参数；2：Solve参数；3：Setup+Solve参数 */
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
      MatFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_A.bin";
      RhsFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_b.bin";
      file_base = 2;
      break;
   case 2:
      MatFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_02/solverchallenge23_02_A.bin";
      RhsFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_02/solverchallenge23_02_b.bin";
      file_base = 2;
      break;
   case 3:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/mat_csr_128X128X128.bin";
      RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 4:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e3.bin";
      RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 5:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e6.bin";
      RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 6:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Checkerboard_mat_csr_128X128X128.bin";
      RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Checkerboard_rhs_128X128X128.bin";
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

   // par_rhs_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
   //                                 jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
   //                                 jxf_hpCSRMatrixRowStarts(hp_matrix));
   // jxf_ParVectorInitialize(par_rhs_F);
   // jxmp_ParVectorDtoF(par_rhs, par_rhs_F);

   // par_sol_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
   //                                 jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
   //                                 jxf_hpCSRMatrixRowStarts(hp_matrix));
   // jxf_ParVectorInitialize(par_sol_F);
   // jxmp_ParVectorDtoF(par_sol, par_sol_F);

   endtime = jx_MPI_Wtime();
   jx_GetWallTime(comm, "BuildParLinearSystem", starttime, endtime, 0, 3);

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

      JX_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JX_PAMGSetMaxLevels(amg_solver, max_levels);
      JX_PAMGSetMaxIter(amg_solver, 1);
      JX_PAMGSetCycleType(amg_solver, cycle_type);
      JX_PAMGSetMeasureType(amg_solver, measure_type);
      JX_PAMGSetRAP2(amg_solver, rap2);
      JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JX_PAMGSetInterpType(amg_solver, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver, relax_wt);
      JX_PAMGSetOuterWt(amg_solver, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond,
                         (JX_PtrToSolverFcn)JX_PAMGSetup, amg_solver);

      JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix)hp_matrix);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                    (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Solve", starttime, endtime, 0, 2);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JX_PAMGDestroy(amg_solver);
      JX_GMRESDestroy(solver);
   }
   break;

   case 2201: /* PAMG-GMRES (Mixed Precision) */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: PAMG-GMRES(%d) - Mixed Precision (Float Precond) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      /* 3. 创建双精度 GMRES 迭代器 (外部始终是双精度) */
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                         (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver_F);

      JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Setup (Mixed)", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Solve (Mixed)", starttime, endtime, 0, 3);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver_F);
      JX_GMRESDestroy(solver);
   }
   break;

   case 2202: /* PAMG-GMRES (Mixed Precision) */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: Adp-PAMG-GMRES(%d) - Adaptive Level-wise Mixed Precision \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      // ==================================================
      JX_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JX_PAMGSetMaxLevels(amg_solver, max_levels);
      JX_PAMGSetMaxIter(amg_solver, 1);
      JX_PAMGSetCycleType(amg_solver, cycle_type);
      JX_PAMGSetMeasureType(amg_solver, measure_type);
      JX_PAMGSetRAP2(amg_solver, rap2);
      JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JX_PAMGSetInterpType(amg_solver, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver, relax_wt);
      JX_PAMGSetOuterWt(amg_solver, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      // ==================================================
      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      /* 3. 创建双精度 GMRES 迭代器 (外部始终是双精度) ===================================== */
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                         (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver_F);

      JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix)hp_matrix);
      // JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);
      JX_PAMGSetup_DtoF(amg_solver, amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Setup (Mixed)", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Solve (Mixed)", starttime, endtime, 0, 3);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver_F);
      JX_PAMGDestroy(amg_solver);
      JX_GMRESDestroy(solver);
   }
   break;

   case 2203: /* PAMG-GMRES (Mixed Precision) */
   {
      JX_Int level;
      JX_Int adp_num_levels;
      JX_Int *adp_level_prec = NULL;
      jx_AdpSolver adp_solver = NULL;

      if (myid == 0)
         jx_printf("\n >>> Solver: 2203 Adp-PAMG-GMRES(%d) - Adaptive Level-wise Mixed Precision \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      /* ==================================================
       * 1. 创建双精度 AMG
       * ================================================== */
      JX_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JX_PAMGSetMaxLevels(amg_solver, max_levels);
      JX_PAMGSetMaxIter(amg_solver, 1);
      JX_PAMGSetCycleType(amg_solver, cycle_type);
      JX_PAMGSetMeasureType(amg_solver, measure_type);
      JX_PAMGSetRAP2(amg_solver, rap2);
      JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JX_PAMGSetInterpType(amg_solver, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver, relax_wt);
      JX_PAMGSetOuterWt(amg_solver, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      /* ==================================================
       * 2. 创建单精度 AMG
       * ================================================== */
      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      /* ==================================================
       * 3. 完成 D/F 两套 AMG setup
       * ================================================== */
      JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix)hp_matrix);
      JX_PAMGSetup_DtoF(amg_solver, amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

      /* ==================================================
       * 4. 创建并绑定 adaptive solver
       * ================================================== */
      if (JX_AdpAMGCreate(&adp_solver))
      {
         jx_printf("ERROR: JX_AdpAMGCreate failed\n");
         exit(1);
      }

      if (JX_AdpAMGSetData(adp_solver,
                           (jx_ParAMGData *)amg_solver,
                           (jxf_ParAMGData *)amg_solver_F))
      {
         jx_printf("ERROR: JX_AdpAMGSetData failed\n");
         jx_printf("D num_levels = %d\n", jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver));
         jx_printf("F num_levels = %d\n", jxf_ParAMGDataNumLevels((jxf_ParAMGData *)amg_solver_F));
         exit(1);
      }

      /* ==================================================
       * 5. 设置每层精度：先做最简单验证
       *    所有层 FP64，只有最粗层 FP32
       * ================================================== */
      adp_num_levels = adp_solver->num_levels;

      if (myid == 0)
      {
         jx_printf("[main] D num_levels = %d\n",
                   jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver));
         jx_printf("[main] F num_levels = %d\n",
                   jxf_ParAMGDataNumLevels((jxf_ParAMGData *)amg_solver_F));
         jx_printf("[main] adp_solver->num_levels = %d\n",
                   adp_solver->num_levels);
         jx_printf("[main] adp_num_levels = %d\n", adp_num_levels);
      }

      adp_level_prec = (JX_Int *)jx_CTAlloc(JX_Int, adp_num_levels);
      // if (adp_level_prec == NULL)
      // {
      //    jx_printf("ERROR: failed to allocate adp_level_prec\n");
      //    exit(1);
      // }

      // for (level = 0; level < adp_num_levels; level++)
      // {
      //    adp_level_prec[level] = JX_PRECISION_FP64;
      // }

      // adp_level_prec[adp_num_levels - 1] = JX_PRECISION_FP32;

      // if (JX_AdpAMGSetLevelPrec(adp_solver, adp_num_levels, adp_level_prec))
      // {
      //    jx_printf("ERROR: JX_AdpAMGSetLevelPrec failed\n");
      //    jx_TFree(adp_level_prec);
      //    adp_level_prec = NULL;
      //    exit(1);
      // }

      // JX_Int num_fp32_levels = adp_num_levels; /* 先测最后 1 层，可改成 2、3、4... */
      JX_Int num_fp32_levels = 1;
      JX_Int num_fp64_levels = 0;
      JX_Int num_fp32_levels_real = 0;

      /* -------------------------------
       * 先默认所有层都是 FP64
       * ------------------------------- */
      for (level = 0; level < adp_num_levels; level++)
      {
         adp_level_prec[level] = JX_PRECISION_FP64;
      }

      /* -------------------------------
       * 最后 num_fp32_levels 层切到 FP32
       * ------------------------------- */
      // for (level = adp_num_levels - num_fp32_levels; level < adp_num_levels; level++)
      // {
      //    if (level >= 0)
      //    {
      //       adp_level_prec[level] = JX_PRECISION_FP32;
      //    }
      // }
      if (adp_num_levels > 0)
      {
         adp_level_prec[0] = JX_PRECISION_FP32;
      }

      /* -------------------------------
       * 打印当前 level-wise 精度配置
       * ------------------------------- */
      if (myid == 0)
      {
         jx_printf("[main] adp_num_levels = %d\n", adp_num_levels);
         jx_printf("[main] requested num_fp32_levels = %d\n", num_fp32_levels);
         jx_printf("[main] level precision schedule:\n");

         for (level = 0; level < adp_num_levels; level++)
         {
            if (adp_level_prec[level] == JX_PRECISION_FP32)
            {
               num_fp32_levels_real++;
               jx_printf("  level %d -> FP32\n", level);
            }
            else
            {
               num_fp64_levels++;
               jx_printf("  level %d -> FP64\n", level);
            }
         }

         jx_printf("[main] summary: FP64 levels = %d, FP32 levels = %d\n",
                   num_fp64_levels, num_fp32_levels_real);
      }

      if (JX_AdpAMGSetLevelPrec(adp_solver, adp_num_levels, adp_level_prec))
      {
         jx_printf("ERROR: JX_AdpAMGSetLevelPrec failed\n");
         jx_TFree(adp_level_prec);
         adp_level_prec = NULL;
         exit(1);
      }

      /* SetLevelPrec 内部已经拷贝了一份，这里释放临时数组 */
      jx_TFree(adp_level_prec);
      adp_level_prec = NULL;

      /* 3. 创建双精度 GMRES 迭代器 (外部始终是双精度) ===================================== */
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      // JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
      //                    (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver_F);

      JX_GMRESSetPrecond(solver,
                         (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP_V2,
                         (JX_PtrToSolverFcn)JX_AdpAMGSetup,
                         (JX_Solver)adp_solver);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Setup (Mixed)", starttime, endtime, 0, 3);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-GMRES Solve (Mixed)", starttime, endtime, 0, 3);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      // jx_TFree(level_prec);
      // level_prec = NULL;
      JX_AdpAMGDestroy(adp_solver);
      JXF_PAMGDestroy(amg_solver_F);
      JX_PAMGDestroy(amg_solver);
      JX_GMRESDestroy(solver);
   }
   break;

   case 42: /* PAMG-FlexGMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: PAMG-FlexGMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JX_PAMGSetMaxLevels(amg_solver, max_levels);
      JX_PAMGSetMaxIter(amg_solver, 1);
      JX_PAMGSetCycleType(amg_solver, cycle_type);
      JX_PAMGSetMeasureType(amg_solver, measure_type);
      JX_PAMGSetRAP2(amg_solver, rap2);
      JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JX_PAMGSetInterpType(amg_solver, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver, relax_wt);
      JX_PAMGSetOuterWt(amg_solver, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      JX_ParCSRFlexGMRESCreate(comm, &solver);
      JX_FlexGMRESSetKDim(solver, k_dim);
      JX_FlexGMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_FlexGMRESSetMaxIter(solver, max_iter);
      JX_FlexGMRESSetTol(solver, tol);
      JX_FlexGMRESSetLogging(solver, 1);
      JX_FlexGMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_FlexGMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond,
                             (JX_PtrToSolverFcn)JX_PAMGSetup, amg_solver);

      JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix)hp_matrix);

      /* this is optional - could be a user defined one instead */
      JX_FlexGMRESSetModifyPC(solver, (JX_PtrToModifyPCFcn)jx_FlexGMRESModifyPCDefault);

      JX_FlexGMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_FlexGMRESSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                        (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Solve", starttime, endtime, 0, 2);
      }

      JX_FlexGMRESGetNumIterations(solver, &num_iterations);
      JX_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JX_PAMGDestroy(amg_solver);
      JX_ParCSRFlexGMRESDestroy(solver);
   }
   break;

   case 4201: /* PAMG-FlexGMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: MxP PAMG-FlexGMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      JX_ParCSRFlexGMRESCreate(comm, &solver);
      JX_FlexGMRESSetKDim(solver, k_dim);
      JX_FlexGMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_FlexGMRESSetMaxIter(solver, max_iter);
      JX_FlexGMRESSetTol(solver, tol);
      JX_FlexGMRESSetLogging(solver, 1);
      JX_FlexGMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_FlexGMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                             (JX_PtrToSolverFcn)JXF_PAMGSetup, amg_solver_F);

      JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

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

      JXF_PAMGDestroy(amg_solver_F);
      JX_ParCSRFlexGMRESDestroy(solver);
   }
   break;

   case 4202: /* PAMG-FlexGMRES(half 版本) */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: MxP half PAMG-FlexGMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      JX_ParCSRFlexGMRESCreate(comm, &solver);
      JX_FlexGMRESSetKDim(solver, k_dim);
      JX_FlexGMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_FlexGMRESSetMaxIter(solver, max_iter);
      JX_FlexGMRESSetTol(solver, tol);
      JX_FlexGMRESSetLogging(solver, 1);
      JX_FlexGMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_FlexGMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                             (JX_PtrToSolverFcn)JXF_PAMGSetup, amg_solver_F);

      JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

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

      JXF_PAMGDestroy(amg_solver_F);
      JX_ParCSRFlexGMRESDestroy(solver);
   }
   break;

   case 4203: /* PAMG-GMRES (Mixed Precision) */
   {
      JX_Int level;
      JX_Int adp_num_levels;
      JX_Int *adp_level_prec = NULL;
      jx_AdpSolver adp_solver = NULL;

      if (myid == 0)
         jx_printf("\n >>> Solver: 2203 Adp-PAMG-GMRES(%d) - Adaptive Level-wise Mixed Precision \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      /* ==================================================
       * 1. 创建双精度 AMG
       * ================================================== */
      JX_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JX_PAMGSetMaxLevels(amg_solver, max_levels);
      JX_PAMGSetMaxIter(amg_solver, 1);
      JX_PAMGSetCycleType(amg_solver, cycle_type);
      JX_PAMGSetMeasureType(amg_solver, measure_type);
      JX_PAMGSetRAP2(amg_solver, rap2);
      JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JX_PAMGSetInterpType(amg_solver, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver, relax_wt);
      JX_PAMGSetOuterWt(amg_solver, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      /* ==================================================
       * 2. 创建单精度 AMG
       * ================================================== */
      JXF_PAMGCreate(&amg_solver_F);
      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver_F, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver_F, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver_F, max_levels);
      JXF_PAMGSetMaxIter(amg_solver_F, 1);
      JXF_PAMGSetCycleType(amg_solver_F, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver_F, measure_type);
      JXF_PAMGSetRAP2(amg_solver_F, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver_F, keepTranspose);
      JXF_PAMGSetCoarsenType(amg_solver_F, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver_F, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver_F, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver_F, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver_F, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver_F, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver_F, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver_F, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver_F, amg_print_level);
      JXF_PAMGSetCoarseThreshold(amg_solver_F, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver_F, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver_F, outer_wt);
      JXF_PAMGSetSCommPkgSwitch(amg_solver_F, S_commpkg_switch);
      JXF_PAMGSetAIRStrongTh(amg_solver_F, AIR_strong_th);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);          /* relax_type for "coarsest" */

      /* ==================================================
       * 3. 完成 D/F 两套 AMG setup
       * ================================================== */
      JX_PAMGSetup(amg_solver, (JX_hpCSRMatrix)hp_matrix);
      JX_PAMGSetup_DtoF(amg_solver, amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

      /* ==================================================
       * 4. 创建并绑定 adaptive solver
       * ================================================== */
      if (JX_AdpAMGCreate(&adp_solver))
      {
         jx_printf("ERROR: JX_AdpAMGCreate failed\n");
         exit(1);
      }

      if (JX_AdpAMGSetData(adp_solver,
                           (jx_ParAMGData *)amg_solver,
                           (jxf_ParAMGData *)amg_solver_F))
      {
         jx_printf("ERROR: JX_AdpAMGSetData failed\n");
         jx_printf("D num_levels = %d\n", jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver));
         jx_printf("F num_levels = %d\n", jxf_ParAMGDataNumLevels((jxf_ParAMGData *)amg_solver_F));
         exit(1);
      }

      /* ==================================================
       * 5. 设置每层精度
       * ================================================== */
      adp_num_levels = adp_solver->num_levels;

      if (myid == 0)
      {
         jx_printf("[main] D num_levels = %d\n",
                   jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver));
         jx_printf("[main] F num_levels = %d\n",
                   jxf_ParAMGDataNumLevels((jxf_ParAMGData *)amg_solver_F));
         jx_printf("[main] adp_solver->num_levels = %d\n",
                   adp_solver->num_levels);
         jx_printf("[main] adp_num_levels = %d\n", adp_num_levels);
      }

      adp_level_prec = (JX_Int *)jx_CTAlloc(JX_Int, adp_num_levels);

      JX_Int num_fp32_levels = adp_num_levels; /* 先测最后 1 层，可改成 2、3、4... */
      JX_Int num_fp64_levels = 0;
      JX_Int num_fp32_levels_real = 0;

      /* -------------------------------
       * 先默认所有层都是 FP64
       * ------------------------------- */
      for (level = 0; level < adp_num_levels; level++)
      {
         adp_level_prec[level] = JX_PRECISION_FP64;
      }

      /* -------------------------------
       * 最后 num_fp32_levels 层切到 FP32
       * ------------------------------- */
      for (level = adp_num_levels - num_fp32_levels; level < adp_num_levels; level++)
      {
         if (level >= 0)
         {
            adp_level_prec[level] = JX_PRECISION_FP32;
         }
      }

      /* -------------------------------
       * 打印当前 level-wise 精度配置
       * ------------------------------- */
      if (myid == 0)
      {
         jx_printf("[main] adp_num_levels = %d\n", adp_num_levels);
         jx_printf("[main] requested num_fp32_levels = %d\n", num_fp32_levels);
         jx_printf("[main] level precision schedule:\n");

         for (level = 0; level < adp_num_levels; level++)
         {
            if (adp_level_prec[level] == JX_PRECISION_FP32)
            {
               num_fp32_levels_real++;
               jx_printf("  level %d -> FP32\n", level);
            }
            else
            {
               num_fp64_levels++;
               jx_printf("  level %d -> FP64\n", level);
            }
         }

         jx_printf("[main] summary: FP64 levels = %d, FP32 levels = %d\n",
                   num_fp64_levels, num_fp32_levels_real);
      }

      if (JX_AdpAMGSetLevelPrec(adp_solver, adp_num_levels, adp_level_prec))
      {
         jx_printf("ERROR: JX_AdpAMGSetLevelPrec failed\n");
         jx_TFree(adp_level_prec);
         adp_level_prec = NULL;
         exit(1);
      }

      /* SetLevelPrec 内部已经拷贝了一份，这里释放临时数组 */
      jx_TFree(adp_level_prec);
      adp_level_prec = NULL;

      /* 3. 创建双精度 GMRES 迭代器 (外部始终是双精度) ===================================== */
      JX_ParCSRFlexGMRESCreate(comm, &solver);
      JX_FlexGMRESSetKDim(solver, k_dim);
      JX_FlexGMRESSetMaxIter(solver, max_iter);
      JX_FlexGMRESSetTol(solver, tol);
      JX_FlexGMRESSetLogging(solver, 1);
      JX_FlexGMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      // JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
      //                    (JX_PtrToSolverFcn)JXF_PAMGSetup, (JX_Solver)amg_solver_F);

      JX_FlexGMRESSetPrecond(solver,
                             (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP_V2,
                             (JX_PtrToSolverFcn)JX_AdpAMGSetup,
                             (JX_Solver)adp_solver);

      JX_FlexGMRESSetModifyPC(solver, (JX_PtrToModifyPCFcn)jx_FlexGMRESModifyPCDefault);
      JX_FlexGMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_FlexGMRESSolve(solver, (JX_Matrix)hp_matrix, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-FlexGMRES Solve", starttime, endtime, 0, 2);
      }

      JX_FlexGMRESGetNumIterations(solver, &num_iterations);
      JX_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JX_AdpAMGDestroy(adp_solver);
      JXF_PAMGDestroy(amg_solver_F);
      JX_PAMGDestroy(amg_solver);
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

// /*!
//  * \fn JX_Int jx_ParVectorSetRandomValues
//  */
// JX_Int
// jx_ParVectorSetRandomValues(jx_ParVector *v, JX_Int seed)
// {
//    JX_Int my_id;
//    jx_Vector *v_local = jx_ParVectorLocalVector(v);

//    MPI_Comm comm = jx_ParVectorComm(v);
//    jx_MPI_Comm_rank(comm, &my_id);

//    seed *= (my_id + 1);

//    return jx_SeqVectorSetRandomValues(v_local, seed);
// }

// /*!
//  * \fn JX_Int jx_SeqVectorSetRandomValues
//  */
// JX_Int
// jx_SeqVectorSetRandomValues(jx_Vector *x, JX_Int seed)
// {
//    JX_Real *vector_data = jx_VectorData(x);
//    JX_Int size = jx_VectorSize(x);
//    JX_Int i, ierr = 0;

//    jx_SeedRand(seed);

//    size *= jx_VectorNumVectors(x);

//    /* RDF: threading this loop may cause problems because of jx_Rand() */
//    for (i = 0; i < size; i++)
//    {
//       vector_data[i] = 2.0 * jx_Rand() - 1.0;
//    }

//    return ierr;
// }

// /* 高精度矩阵转低精度 */
// JXF_Int
// jxmp_CSRMatrixDtoF(jx_CSRMatrix *A, jxf_CSRMatrix *B, JX_Int copy_data)
// {
//    // jx_printf("\n >>> jxmp_CSRMatrixDtoF \n");

//    JX_Int ierr = 0;
//    JX_Int num_rows = jx_CSRMatrixNumRows(A);
//    JX_Int *A_i = jx_CSRMatrixI(A);
//    JX_Int *A_j = jx_CSRMatrixJ(A);
//    JX_Real *A_data;

//    JX_Int *B_i = jxf_CSRMatrixI(B);
//    JX_Int *B_j = jxf_CSRMatrixJ(B);
//    JXF_Real *B_data;
//    JX_Int i, j;

//    for (i = 0; i < num_rows; i++)
//    {
//       B_i[i] = A_i[i];
//       for (j = A_i[i]; j < A_i[i + 1]; j++)
//       {
//          B_j[j] = A_j[j];
//       }
//    }
//    B_i[num_rows] = A_i[num_rows];

//    if (copy_data)
//    {
//       A_data = jx_CSRMatrixData(A);
//       B_data = jxf_CSRMatrixData(B);

//       for (i = 0; i < num_rows; i++)
//       {
//          for (j = A_i[i]; j < A_i[i + 1]; j++)
//          {

//             B_data[j] = (JXF_Real)A_data[j];
//          }
//       }
//    }

//    return ierr;
// }

// JXF_Int
// jxmp_ParCSRMatrixDtoF(jx_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JX_Int copy_data)
// {
//    // jx_printf("\n >>> jxmp_ParCSRMatrixDtoF \n");

//    jx_CSRMatrix *A_diag;
//    jx_CSRMatrix *A_offd;
//    JX_Int *col_map_offd_A; // A 非对角块的列映射表
//    JX_Int num_cols_offd;
//    jxf_CSRMatrix *B_diag;
//    jxf_CSRMatrix *B_offd;
//    JXF_Int *col_map_offd_B;
//    JXF_Int i;

//    A_diag = jx_ParCSRMatrixDiag(A);
//    A_offd = jx_ParCSRMatrixOffd(A);
//    col_map_offd_A = jx_ParCSRMatrixColMapOffd(A);
//    num_cols_offd = jx_CSRMatrixNumCols(A_offd);

//    B_diag = jxf_ParCSRMatrixDiag(B);
//    B_offd = jxf_ParCSRMatrixOffd(B);
//    col_map_offd_B = jxf_ParCSRMatrixColMapOffd(B);

//    jxmp_CSRMatrixDtoF(A_diag, B_diag, copy_data);

//    jxmp_CSRMatrixDtoF(A_offd, B_offd, copy_data);

//    if (num_cols_offd && col_map_offd_B == NULL)
//    {
//       col_map_offd_B = jx_CTAlloc(JX_Int, num_cols_offd);
//       jxf_ParCSRMatrixColMapOffd(B) = col_map_offd_B;
//    }
//    for (i = 0; i < num_cols_offd; i++)
//    {
//       col_map_offd_B[i] = col_map_offd_A[i];
//    }

//    return jx_error_flag;
// }

// /* 高精度向量转低精度 */
// JXF_Int
// jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y)
// {
//    jx_Vector *x_local = jx_ParVectorLocalVector(x);
//    jxf_Vector *y_local = jxf_ParVectorLocalVector(y);

//    return jxmp_SeqVectorDtoF(x_local, y_local);
// }

// JXF_Int
// jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y)
// {
//    JX_Real *x_data = jx_VectorData(x);
//    JX_Int size = jx_VectorSize(x);
//    JXF_Real *y_data = jxf_VectorData(y);
//    JX_Int i;

//    size *= jx_VectorNumVectors(x);

// #define JX_SMP_PRIVATE i
// #include "../jxpamg-basic-DGAI/include/jx_smp_forloop.h"
//    for (i = 0; i < size; i++)
//    {
//       y_data[i] = (JXF_Real)x_data[i]; // 双精度到单精度转换
//    }

//    return 0;
// }

// /* 低精度向量转高精度 */
// JXF_Int
// jxmp_ParVectorFtoD(jxf_ParVector *x, jx_ParVector *y)
// {
//    // jx_printf("\n >>> jxmp_ParVectorFtoD \n");

//    jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
//    jx_Vector *y_local = jx_ParVectorLocalVector(y);

//    return jxmp_SeqVectorFtoD(x_local, y_local);
// }

// JXF_Int
// jxmp_SeqVectorFtoD(jxf_Vector *x, jx_Vector *y)
// {
//    // jx_printf("\n >>> jxmp_SeqVectorFtoD \n");

//    JXF_Real *x_data = jxf_VectorData(x);
//    JX_Real *y_data = jx_VectorData(y);
//    JX_Int size = jxf_VectorSize(x);
//    JX_Int i = 0;

//    size *= jxf_VectorNumVectors(x);

// #define JX_SMP_PRIVATE i
// #include "../jxpamg-basic-DGAI/include/jx_smp_forloop.h"
//    for (i = 0; i < size; i++)
//    {
//       y_data[i] = (JX_Real)x_data[i]; // 双精度到单精度转换
//    }

//    return 0;
// }

// //===================================================
// static JXF_Real *
// jxmp_CopyRealArrayDtoF(const JX_Real *src, JXF_Int n)
// {
//    JXF_Int i;
//    JXF_Real *dst;

//    if (src == NULL || n <= 0)
//    {
//       return NULL;
//    }

//    dst = jxf_CTAlloc(JXF_Real, n);
//    for (i = 0; i < n; i++)
//    {
//       dst[i] = (JXF_Real)src[i];
//    }

//    return dst;
// }

// static JXF_Int *
// jxmp_CopyIntArrayDtoF(const JX_Int *src, JXF_Int n)
// {
//    JXF_Int i;
//    JXF_Int *dst;

//    if (src == NULL || n <= 0)
//    {
//       return NULL;
//    }

//    dst = jxf_CTAlloc(JXF_Int, n);
//    for (i = 0; i < n; i++)
//    {
//       dst[i] = (JXF_Int)src[i];
//    }

//    return dst;
// }

// JXF_Int
// JXF_PAMGSetup_DtoF(void *amg_solver, void *amg_solver_F, jxf_hpCSRMatrix *hp_matrix_F)
// {
//    /* D 端：已经 setup 完成的双精度 AMG */
//    jx_ParAMGData *amg_data_D = amg_solver;
//    /* F 端：准备由 D->F 填充的单精度 AMG */
//    jxf_ParAMGData *amg_data_F = amg_solver_F;

//    MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix_F);

//    /* Data Structure variables */
//    jxf_hpCSRMatrix **A_array_F;
//    jxf_ParVector **F_array_F;
//    jxf_ParVector **U_array_F;
//    jxf_ParVector *Vtemp_F;
//    jxf_ParVector *Ztemp_F = NULL;
//    jxf_ParCSRMatrix **P_array_F;
//    jxf_ParCSRMatrix **R_array_F;
//    jxf_ParVector *Residual_array_F;
//    JXF_Real **AI_measure_array_F;
//    JXF_Int **CF_marker_array_F;
//    JXF_Int **dof_func_array_F;
//    JXF_Int *dof_func_F;
//    JXF_Int *col_offd_S_to_A_F;
//    JXF_Int *col_offd_Sabs_to_A_F = NULL;
//    JXF_Int *AIR_maxsize_ls_F;

//    JXF_Real strong_threshold;
//    JXF_Real AIR_strong_th;
//    JXF_Real max_row_sum;
//    JXF_Real trunc_factor;
//    JXF_Real S_commpkg_switch;
//    JXF_Int max_levels;
//    JXF_Int amg_logging;
//    JXF_Int amg_print_level;
//    JXF_Int debug_flag;
//    JXF_Int local_num_vars;
//    JXF_Int P_max_elmts;
//    JXF_Int R_max_size;

//    JX_Solver *smoother = NULL;
//    JXF_Solver *smoother_F = NULL;
//    JXF_Int smooth_type = jxf_ParAMGDataSmoothType(amg_data_F);
//    JXF_Int smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data_F);
//    // JX_Int smooth_type = jx_ParAMGDataSmoothType(amg_data_D);
//    // JX_Int smooth_num_levels = jx_ParAMGDataSmoothNumLevels(amg_data_D);
//    char *euclidfile;
//    JXF_Int eu_level;
//    JXF_Int eu_bj;
//    JXF_Real eu_sparse_A_F;

//    /* Local variables */
//    /* Local variables (mainly F-side temporaries) */
//    JXF_Real *AI_measure_F;
//    JXF_Real *AI_measure2_F;
//    JXF_Int *CF_marker_F;
//    JXF_Int *CFN_marker_F = NULL;
//    jxf_ParCSRMatrix *par_S_F;
//    jxf_ParCSRMatrix *Sabs_F = NULL;
//    jxf_ParCSRMatrix *par_S2_F;
//    jxf_ParCSRMatrix *par_P_F;
//    jxf_ParCSRMatrix *R_F;
//    jxf_hpCSRMatrix *A_H_F;

//    JXF_Int wall_time_option;         /* newly added Yue Xiaoqiang 2015/09/30 */
//    JXF_Real wall_time_coarsen = 0.0; /* newly added Yue Xiaoqiang 2015/09/30 */
//    JXF_Real wall_time_rap = 0.0;     /* newly added Yue Xiaoqiang 2015/09/30 */
//    JXF_Real wall_time_interp = 0.0;  /* newly added Yue Xiaoqiang 2015/09/30 */
//    JXF_Real tmp_wall_time = 0.0;     /* newly added Yue Xiaoqiang 2015/09/30 */

//    JXF_Int old_num_levels, num_levels;
//    JXF_Int level;
//    JXF_Int local_size, i;
//    JXF_Int first_local_row;
//    JXF_Int coarse_size_last = -1;
//    JXF_Int coarse_size;
//    JXF_Int coarsen_type;
//    JXF_Int interp_type;
//    JXF_Int restri_type;
//    JXF_Int measure_type;
//    JXF_Int setup_type;
//    JXF_Int fine_size;
//    JXF_Int rest, tms, indx;
//    JXF_Int not_finished_coarsening = 1;
//    JXF_Int Setup_err_flag = 0;
//    JXF_Int coarse_threshold;   /* we don't fix it to be 9. peghoty 2010/04/14 */
//    JXF_Int coarse_threshold_2; /* by xu: 2013/11/25 */
//    JXF_Int j, k;
//    JXF_Int num_procs, my_id, num_threads;

//    JXF_Int *grid_relax_type = jxf_ParAMGDataGridRelaxType(amg_data_F);
//    JXF_Int relax_order = jxf_ParAMGDataRelaxOrder(amg_data_F);
//    JXF_Int num_functions = jxf_ParAMGDataNumFunctions(amg_data_F);
//    JXF_Int spmt_rap_type = jxf_ParAMGDataSpMtRapType(amg_data_F);
//    JXF_Int ai_measure_type = jxf_ParAMGDataAIMeasureType(amg_data_F);
//    JXF_Int num_paths = jxf_ParAMGDataNumPaths(amg_data_F);
//    JXF_Int agg_num_levels = jxf_ParAMGDataAggNumLevels(amg_data_F);
//    JXF_Int agg_interp_type = jxf_ParAMGDataAggInterpType(amg_data_F);
//    JXF_Int agg_P_max_elmts = jxf_ParAMGDataAggPMaxElmts(amg_data_F);
//    /* JXF_Int        agg_P12_max_elmts = jxf_ParAMGDataAggP12MaxElmts(amg_data_F); */
//    JXF_Real agg_trunc_factor = jxf_ParAMGDataAggTruncFactor(amg_data_F);
//    /* JXF_Real       agg_P12_trunc_factor = jxf_ParAMGDataAggP12TruncFactor(amg_data_F); */
//    JXF_Int rap2 = jxf_ParAMGDataRAP2(amg_data_F);
//    JXF_Int keepTranspose = jxf_ParAMGDataKeepTranspose(amg_data_F);
//    JXF_Int print_coarse_matrix = jxf_ParAMGDataPrintCoarseSystem(amg_data_F);

//    JXF_Int *opt_icor;
//    JXF_Int *coarse_dof_func;
//    JXF_Int *coarse_pnts_global;
//    JXF_Int *coarse_pnts_global1;
//    JXF_Real size;
//    JXF_Real coarse_ratio;          /* newly added peghoty 2010/04/14 */
//    JXF_Real wall_time = 0.0;       /* for debugging instrumentation */
//    JXF_Int measure_type_rlx;       /* newly added peghoty 2010/05/29 */
//    JXF_Int number_syn_rlx;         /* newly added peghoty 2010/05/29 */
//    JXF_Real measure_threshold_rlx; /* newly added peghoty 2010/05/29 */
//    JXF_Real **l1_norms = NULL;

//    /* coarse matrices */
//    char FileNameCoaMat[256];

//    /* AI statistic information */
//    JXF_Int num_vars_local = 0, num_vars_global;
//    JXF_Int num_ai_local = 0, num_ai_global;
//    JXF_Int num_ai_local_valid = 0, num_ai_global_valid;
//    JXF_Int num_ai_c_local = 0, num_ai_c_global;
//    JXF_Real mai_local = 0.0, mai_global;
//    JXF_Real mai_local_valid = 0.0, mai_global_valid;
//    JXF_Real mai_c_local = 0.0, mai_c_global;

//    JXF_Int num_vars_local_0 = 0;
//    JXF_Int num_ai_local_0 = 0;
//    JXF_Int num_ai_local_valid_0 = 0;
//    JXF_Int num_ai_c_local_0 = 0;
//    JXF_Real mai_local_0 = 0.0;
//    JXF_Real mai_local_valid_0 = 0.0;
//    JXF_Real mai_c_local_0 = 0.0;
//    JXF_Real mai_threshold = 0.1;
//    JXF_Real gai_threshold = 0.3;

//    jxf_MPI_Comm_size(comm, &num_procs);
//    jxf_MPI_Comm_rank(comm, &my_id);

//    if ((num_procs > 1) && (spmt_rap_type != 1)) /* Yue Xiaoqiang 2012/10/12 */
//    {
//       spmt_rap_type = 1;
//    }

//    num_threads = jxf_NumThreads();                       /* Yue Xiaoqiang 2012/10/12 */
//    opt_icor = jxf_CTAlloc(JXF_Int, 5 * num_threads + 2); /* Yue Xiaoqiang 2012/10/12 */

//    wall_time_option = jxf_ParAMGDataWallTimeOption(amg_data_F);
//    old_num_levels = jxf_ParAMGDataNumLevels(amg_data_F);
//    max_levels = jxf_ParAMGDataMaxLevels(amg_data_F);
//    amg_logging = jxf_ParAMGDataLogging(amg_data_F);
//    amg_print_level = jxf_ParAMGDataPrintLevel(amg_data_F);
//    coarsen_type = jxf_ParAMGDataCoarsenType(amg_data_F);
//    measure_type = jxf_ParAMGDataMeasureType(amg_data_F);
//    setup_type = jxf_ParAMGDataSetupType(amg_data_F);
//    debug_flag = jxf_ParAMGDataDebugFlag(amg_data_F);
//    dof_func_F = jxf_ParAMGDataDofFunc(amg_data_F);
//    interp_type = jxf_ParAMGDataInterpType(amg_data_F);
//    restri_type = jxf_ParAMGDataRestriction(amg_data_F);
//    AIR_strong_th = jxf_ParAMGDataAIRStrongTh(amg_data_F);
//    euclidfile = jxf_ParAMGDataEuclidFile(amg_data_F);
//    eu_level = jxf_ParAMGDataEuLevel(amg_data_F);
//    eu_bj = jxf_ParAMGDataEuBJ(amg_data_F);
//    eu_sparse_A_F = jxf_ParAMGDataEuSparseA(amg_data_F);
//    coarse_threshold = jxf_ParAMGDataCoarseThreshold(amg_data_F);          /* newly added peghoty 2010/04/14 */
//    coarse_ratio = jxf_ParAMGDataCoarseRatio(amg_data_F);                  /* newly added peghoty 2010/04/14 */
//    measure_type_rlx = jxf_ParAMGDataMeasureTypeRlx(amg_data_F);           /* newly added peghoty 2010/05/29 */
//    number_syn_rlx = jxf_ParAMGDataNumberSynRlx(amg_data_F);               /* newly added peghoty 2010/05/29 */
//    measure_threshold_rlx = jxf_ParAMGDataMeasureThresholdRlx(amg_data_F); /* newly added peghoty 2010/05/29 */
//    coarse_threshold_2 = -1;

//    jxf_ParCSRMatrixSetNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix_F));
//    jxf_ParCSRMatrixSetDNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix_F));
//    jxf_ParAMGDataNumVariables(amg_data_F) = jxf_hpCSRMatrixNumRows(hp_matrix_F);

//    if (setup_type == 0)
//    {
//       return Setup_err_flag;
//    }

//    par_S_F = NULL;
//    A_H_F = NULL;

//    A_array_F = jxf_hpAMGDataAArray(amg_data_F);
//    P_array_F = jxf_ParAMGDataPArray(amg_data_F);
//    R_array_F = jxf_ParAMGDataRArray(amg_data_F);
//    AIR_maxsize_ls_F = jxf_ParAMGDataAIRMaxSizeLS(amg_data_F);
//    AI_measure_array_F = jxf_ParAMGDataAIMeasureArray(amg_data_F);
//    CF_marker_array_F = jxf_ParAMGDataCFMarkerArray(amg_data_F);
//    dof_func_array_F = jxf_ParAMGDataDofFuncArray(amg_data_F);

//    local_size = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix_F));

//    grid_relax_type[3] = jxf_ParAMGDataUserCoarseRelaxType(amg_data_F);

//    if (A_array_F || P_array_F || R_array_F || AIR_maxsize_ls_F ||
//        AI_measure_array_F || CF_marker_array_F || dof_func_array_F)
//    {
//       for (j = 1; j < old_num_levels; j++)
//       {
//          if (A_array_F[j])
//          {
//             jxf_hpCSRMatrixDestroy(A_array_F[j]);
//             A_array_F[j] = NULL;
//          }

//          if (dof_func_array_F[j])
//          {
//             jxf_TFree(dof_func_array_F[j]);
//             dof_func_array_F[j] = NULL;
//          }
//       }

//       for (j = 0; j < old_num_levels - 1; j++)
//       {
//          if (P_array_F[j])
//          {
//             jxf_ParCSRMatrixDestroy(P_array_F[j]);
//             P_array_F[j] = NULL;
//          }
//          if (R_array_F[j])
//          {
//             jxf_ParCSRMatrixDestroy(R_array_F[j]);
//             R_array_F[j] = NULL;
//          }
//       }

//       jxf_TFree(AIR_maxsize_ls_F);
//       AIR_maxsize_ls_F = NULL;

//       /*-------------------------------------------------------------------
//        *  Special case use of CF_marker_array when old_num_levels = 1
//        *  requires us to attempt this deallocation every time
//        *------------------------------------------------------------------*/

//       if (CF_marker_array_F && CF_marker_array_F[0])
//       {
//          jxf_TFree(CF_marker_array_F[0]);
//          CF_marker_array_F[0] = NULL;
//       }

//       for (j = 1; CF_marker_array_F && j < old_num_levels - 1; j++)
//       {
//          if (CF_marker_array_F[j])
//          {
//             jxf_TFree(CF_marker_array_F[j]);
//             CF_marker_array_F[j] = NULL;
//          }
//       }

//       /* for AI_measure: added by xwxu */
//       if (AI_measure_array_F && AI_measure_array_F[0])
//       {
//          jxf_TFree(AI_measure_array_F[0]);
//          AI_measure_array_F[0] = NULL;
//       }

//       for (j = 1; AI_measure_array_F && j < old_num_levels - 1; j++)
//       {
//          if (AI_measure_array_F[j])
//          {
//             jxf_TFree(AI_measure_array_F[j]);
//             AI_measure_array_F[j] = NULL;
//          }
//       }
//    }

//    if (A_array_F == NULL)
//    {
//       A_array_F = jxf_CTAlloc(jxf_hpCSRMatrix *, max_levels);
//    }

//    if (P_array_F == NULL && max_levels > 1)
//    {
//       P_array_F = jxf_CTAlloc(jxf_ParCSRMatrix *, max_levels - 1);
//    }

//    /* If restri_type != 0, R != P^T, allocate R matrices */
//    if (restri_type)
//    {
//       if (R_array_F == NULL && max_levels > 1)
//       {
//          R_array_F = jxf_CTAlloc(jxf_ParCSRMatrix *, max_levels - 1);
//          AIR_maxsize_ls_F = jxf_CTAlloc(JXF_Int, max_levels - 1);
//       }
//    }

//    if (AI_measure_array_F == NULL)
//    {
//       AI_measure_array_F = jxf_CTAlloc(JXF_Real *, max_levels);
//    }
//    if (CF_marker_array_F == NULL)
//    {
//       CF_marker_array_F = jxf_CTAlloc(JXF_Int *, max_levels);
//    }
//    if (dof_func_array_F == NULL)
//    {
//       dof_func_array_F = jxf_CTAlloc(JXF_Int *, max_levels);
//    }

//    /*
//     * 这一版先保留基于 hp_matrix_F 的 dof_func_F 初始化逻辑，
//     * 后面如果改成“直接复制 D 端 dof_func / dof_func_array”，再替换这里。
//     */
//    if (num_functions > 1 && dof_func_F == NULL)
//    {
//       first_local_row = jxf_hpCSRMatrixFirstRowIndex(hp_matrix_F);
//       dof_func_F = jxf_CTAlloc(JXF_Int, local_size);

//       rest = first_local_row - ((first_local_row / num_functions) * num_functions);
//       indx = num_functions - rest;
//       if (rest == 0)
//       {
//          indx = 0;
//       }

//       k = num_functions - 1;
//       for (j = indx - 1; j > -1; j--)
//       {
//          dof_func_F[j] = k--;
//       }

//       tms = local_size / num_functions;
//       if (tms * num_functions + indx > local_size)
//       {
//          tms--;
//       }

//       for (j = 0; j < tms; j++)
//       {
//          for (k = 0; k < num_functions; k++)
//          {
//             dof_func_F[indx++] = k;
//          }
//       }

//       k = 0;
//       while (indx < local_size)
//       {
//          dof_func_F[indx++] = k++;
//       }

//       jxf_ParAMGDataDofFunc(amg_data_F) = dof_func_F;
//    }

//    /*
//     * 这里先把 level-0 入口挂到 hp_matrix_F。
//     * 后面如果改成“level-0 也从 D 端 A_array_D[0] 转成 F”，再替换这一句。
//     */
//    A_array_F[0] = hp_matrix_F;

//    dof_func_array_F[0] = dof_func_F;
//    jxf_ParAMGDataAIMeasureArray(amg_data_F) = AI_measure_array_F;
//    jxf_ParAMGDataCFMarkerArray(amg_data_F) = CF_marker_array_F;
//    jxf_ParAMGDataDofFuncArray(amg_data_F) = dof_func_array_F;
//    jxf_hpAMGDataAArray(amg_data_F) = A_array_F;
//    jxf_ParAMGDataPArray(amg_data_F) = P_array_F;

//    /* If R != P^T */
//    if (restri_type)
//    {
//       jxf_ParAMGDataRArray(amg_data_F) = R_array_F;
//       jxf_ParAMGDataAIRMaxSizeLS(amg_data_F) = AIR_maxsize_ls_F;
//    }
//    else
//    {
//       jxf_ParAMGDataRArray(amg_data_F) = P_array_F;
//    }

//    Vtemp_F = jxf_ParAMGDataVtemp(amg_data_F);

//    if (Vtemp_F != NULL)
//    {
//       jxf_ParVectorDestroy(Vtemp_F);
//       Vtemp_F = NULL;
//    }

//    Vtemp_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array_F[0]),
//                                  jxf_hpCSRMatrixGlobalNumRows(A_array_F[0]),
//                                  jxf_hpCSRMatrixRowStarts(A_array_F[0]));
//    jxf_ParVectorInitialize(Vtemp_F);
//    jxf_ParVectorSetPartitioningOwner(Vtemp_F, 0);
//    jxf_ParAMGDataVtemp(amg_data_F) = Vtemp_F;

//    if (num_threads > 1)
//    {
//       Ztemp_F = jxf_ParAMGDataZtemp(amg_data_F);

//       if (Ztemp_F != NULL)
//       {
//          jxf_ParVectorDestroy(Ztemp_F);
//          Ztemp_F = NULL;
//       }

//       for (j = 1; j < 4; j++)
//       {
//          if (grid_relax_type[j] == 3 || grid_relax_type[j] == 4 || grid_relax_type[j] == 6 ||
//              grid_relax_type[j] == 8 || grid_relax_type[j] == 13 || grid_relax_type[j] == 14)
//          {
//             Ztemp_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array_F[0]),
//                                           jxf_hpCSRMatrixGlobalNumRows(A_array_F[0]),
//                                           jxf_hpCSRMatrixRowStarts(A_array_F[0]));
//             jxf_ParVectorInitialize(Ztemp_F);
//             jxf_ParVectorSetPartitioningOwner(Ztemp_F, 0);
//             jxf_ParAMGDataZtemp(amg_data_F) = Ztemp_F;
//             break;
//          }
//       }
//    }

//    F_array_F = jxf_ParAMGDataFArray(amg_data_F);
//    U_array_F = jxf_ParAMGDataUArray(amg_data_F);

//    if (F_array_F != NULL || U_array_F != NULL)
//    {
//       for (j = 1; j < old_num_levels; j++)
//       {
//          if (F_array_F[j] != NULL)
//          {
//             jxf_ParVectorDestroy(F_array_F[j]);
//             F_array_F[j] = NULL;
//          }
//          if (U_array_F[j] != NULL)
//          {
//             jxf_ParVectorDestroy(U_array_F[j]);
//             U_array_F[j] = NULL;
//          }
//       }
//    }

//    if (F_array_F == NULL)
//    {
//       F_array_F = jxf_CTAlloc(jxf_ParVector *, max_levels);
//    }
//    if (U_array_F == NULL)
//    {
//       U_array_F = jxf_CTAlloc(jxf_ParVector *, max_levels);
//    }

//    jxf_ParAMGDataFArray(amg_data_F) = F_array_F;
//    jxf_ParAMGDataUArray(amg_data_F) = U_array_F;

//    /*----------------------------------------------------------
//     *   Initialize jxf_ParAMGData
//     *---------------------------------------------------------*/
//    not_finished_coarsening = 1;
//    level = 0;
//    strong_threshold = jxf_ParAMGDataStrongThreshold(amg_data_F);
//    max_row_sum = jxf_ParAMGDataMaxRowSum(amg_data_F);
//    trunc_factor = jxf_ParAMGDataTruncFactor(amg_data_F);
//    P_max_elmts = jxf_ParAMGDataPMaxElmts(amg_data_F);
//    S_commpkg_switch = jxf_ParAMGDataSCommPkgSwitch(amg_data_F);

//    if (smooth_num_levels > level)
//    {
//       smoother_F = jxf_CTAlloc(JXF_Solver, smooth_num_levels);
//       jxf_ParAMGDataSmoother(amg_data_F) = smoother_F;
//    }

//    /* if AI-based coarsening is used, set ai_measure_type = 1. */
//    if (coarsen_type == 990 || coarsen_type == 991 || coarsen_type == 993 ||
//        coarsen_type == 96 || coarsen_type == 98 || coarsen_type == 910 ||
//        coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 ||
//        coarsen_type == 938 || coarsen_type == 968)
//    {
//       ai_measure_type = 1;
//    }

//    /*-----------------------------------------------------
//     *   DtoF conversion starts here
//     *   DO NOT continue the original coarsening loop below.
//     *-----------------------------------------------------*/

//    if (amg_print_level > 0 && my_id == 0 && ai_measure_type)
//    {
//       jxf_printf(" \n");
//       jxf_printf("================================================================== \n");
//       jxf_printf("+++++++++++++ multi-scale/AI infomation for levels +++++++++++++++ \n");
//    }

//    /*-----------------------------------------------------
//     *   DtoF conversion starts here
//     *-----------------------------------------------------*/

//    jx_hpCSRMatrix **A_array_D;
//    jx_ParCSRMatrix **P_array_D;
//    jx_ParCSRMatrix **R_array_D;
//    jx_ParVector **F_array_D;
//    jx_ParVector **U_array_D;
//    JX_Real **AI_measure_array_D;
//    JX_Int **CF_marker_array_D;
//    JX_Int **dof_func_array_D;
//    JX_Int local_n;

//    num_levels = jx_ParAMGDataNumLevels(amg_data_D);
//    jxf_ParAMGDataNumLevels(amg_data_F) = num_levels;

//    A_array_D = jx_hpAMGDataAArray(amg_data_D);
//    P_array_D = jx_ParAMGDataPArray(amg_data_D);
//    R_array_D = jx_ParAMGDataRArray(amg_data_D);
//    F_array_D = jx_ParAMGDataFArray(amg_data_D);
//    U_array_D = jx_ParAMGDataUArray(amg_data_D);
//    AI_measure_array_D = jx_ParAMGDataAIMeasureArray(amg_data_D);
//    CF_marker_array_D = jx_ParAMGDataCFMarkerArray(amg_data_D);
//    dof_func_array_D = jx_ParAMGDataDofFuncArray(amg_data_D);

//    /* Convert A_array_D to A_array_F */
//    A_array_F[0] = hp_matrix_F;

//    for (level = 1; level < num_levels; level++)
//    {
//       jx_ParCSRMatrix *par_A_D;
//       jxf_ParCSRMatrix *par_A_F;

//       par_A_D = jx_hpCSRMatrixPar(A_array_D[level]);

//       A_array_F[level] = jxf_hpInithpCSRMatrix();

//       par_A_F = jxf_ParCSRMatrixCreate(jx_ParCSRMatrixComm(par_A_D),
//                                        jx_ParCSRMatrixGlobalNumRows(par_A_D),
//                                        jx_ParCSRMatrixGlobalNumCols(par_A_D),
//                                        jx_ParCSRMatrixRowStarts(par_A_D),
//                                        jx_ParCSRMatrixColStarts(par_A_D),
//                                        jx_CSRMatrixNumCols(jx_ParCSRMatrixOffd(par_A_D)),
//                                        jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixDiag(par_A_D)),
//                                        jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixOffd(par_A_D)));

//       jxf_ParCSRMatrixInitialize(par_A_F);

//       /* 如果 Create 默认接管 row_starts / col_starts，这里要显式设为不拥有 */
//       jxf_ParCSRMatrixOwnsRowStarts(par_A_F) = 0;
//       jxf_ParCSRMatrixOwnsColStarts(par_A_F) = 0;

//       jxmp_ParCSRMatrixDtoF(par_A_D, par_A_F, 1);

//       jxf_ParCSRMatrixSetNumNonzeros(par_A_F);
//       jxf_ParCSRMatrixSetDNumNonzeros(par_A_F);

//       jxf_hpCSRMatrixPar(A_array_F[level]) = par_A_F;
//    }

//    /* P / R 只有 0 ... num_levels-2 */
//    for (level = 0; level < num_levels - 1; level++)
//    {
//       if (P_array_D && P_array_D[level])
//       {
//          jx_ParCSRMatrix *par_P_D;
//          jxf_ParCSRMatrix *par_P_F;

//          par_P_D = P_array_D[level];

//          par_P_F = jxf_ParCSRMatrixCreate(jx_ParCSRMatrixComm(par_P_D),
//                                           jx_ParCSRMatrixGlobalNumRows(par_P_D),
//                                           jx_ParCSRMatrixGlobalNumCols(par_P_D),
//                                           jx_ParCSRMatrixRowStarts(par_P_D),
//                                           jx_ParCSRMatrixColStarts(par_P_D),
//                                           jx_CSRMatrixNumCols(jx_ParCSRMatrixOffd(par_P_D)),
//                                           jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixDiag(par_P_D)),
//                                           jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixOffd(par_P_D)));

//          jxf_ParCSRMatrixInitialize(par_P_F);

//          /* row/col starts 直接引用 D 端的，不让 F 端负责释放 */
//          jxf_ParCSRMatrixOwnsRowStarts(par_P_F) = 0;
//          jxf_ParCSRMatrixOwnsColStarts(par_P_F) = 0;

//          jxmp_ParCSRMatrixDtoF(par_P_D, par_P_F, 1);

//          jxf_ParCSRMatrixSetNumNonzeros(par_P_F);
//          jxf_ParCSRMatrixSetDNumNonzeros(par_P_F);

//          P_array_F[level] = par_P_F;
//       }

//       if (restri_type)
//       {
//          if (R_array_D && R_array_D[level])
//          {
//             jx_ParCSRMatrix *par_R_D;
//             jxf_ParCSRMatrix *par_R_F;

//             par_R_D = R_array_D[level];

//             par_R_F = jxf_ParCSRMatrixCreate(jx_ParCSRMatrixComm(par_R_D),
//                                              jx_ParCSRMatrixGlobalNumRows(par_R_D),
//                                              jx_ParCSRMatrixGlobalNumCols(par_R_D),
//                                              jx_ParCSRMatrixRowStarts(par_R_D),
//                                              jx_ParCSRMatrixColStarts(par_R_D),
//                                              jx_CSRMatrixNumCols(jx_ParCSRMatrixOffd(par_R_D)),
//                                              jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixDiag(par_R_D)),
//                                              jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixOffd(par_R_D)));

//             jxf_ParCSRMatrixInitialize(par_R_F);

//             jxf_ParCSRMatrixOwnsRowStarts(par_R_F) = 0;
//             jxf_ParCSRMatrixOwnsColStarts(par_R_F) = 0;

//             jxmp_ParCSRMatrixDtoF(par_R_D, par_R_F, 1);

//             jxf_ParCSRMatrixSetNumNonzeros(par_R_F);
//             jxf_ParCSRMatrixSetDNumNonzeros(par_R_F);

//             R_array_F[level] = par_R_F;
//          }
//       }
//    }

//    jxf_ParAMGDataPArray(amg_data_F) = P_array_F;

//    if (restri_type == 0)
//    {
//       jxf_ParAMGDataRArray(amg_data_F) = P_array_F;
//    }
//    else
//    {
//       jxf_ParAMGDataRArray(amg_data_F) = R_array_F;
//    }

//    /*-----------------------------------------------------
//     * Convert auxiliary vectors: F_array / U_array
//     * Note:
//     *   - level 0 通常不在这里处理，保持与原 setup 风格一致
//     *   - A_array_F[level] 必须已经准备好
//     *-----------------------------------------------------*/
//    for (level = 1; level < num_levels; level++)
//    {
//       if (F_array_D && F_array_D[level])
//       {
//          if (F_array_F[level] == NULL)
//          {
//             F_array_F[level] = jxf_ParVectorCreate(
//                 jxf_hpCSRMatrixComm(A_array_F[level]),
//                 jxf_hpCSRMatrixGlobalNumRows(A_array_F[level]),
//                 jxf_hpCSRMatrixRowStarts(A_array_F[level]));
//             jxf_ParVectorInitialize(F_array_F[level]);
//             jxf_ParVectorSetPartitioningOwner(F_array_F[level], 0);
//          }

//          jxmp_ParVectorDtoF(F_array_D[level], F_array_F[level]);
//       }

//       if (U_array_D && U_array_D[level])
//       {
//          if (U_array_F[level] == NULL)
//          {
//             U_array_F[level] = jxf_ParVectorCreate(
//                 jxf_hpCSRMatrixComm(A_array_F[level]),
//                 jxf_hpCSRMatrixGlobalNumRows(A_array_F[level]),
//                 jxf_hpCSRMatrixRowStarts(A_array_F[level]));
//             jxf_ParVectorInitialize(U_array_F[level]);
//             jxf_ParVectorSetPartitioningOwner(U_array_F[level], 0);
//          }

//          jxmp_ParVectorDtoF(U_array_D[level], U_array_F[level]);
//       }
//    }

//    /*-----------------------------------------------------
//     * Convert auxiliary arrays: AI_measure / CF_marker / dof_func
//     * 这些按层复制即可，level 0 也需要保留
//     *-----------------------------------------------------*/
//    for (level = 0; level < num_levels; level++)
//    {
//       JXF_Int local_n = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array_D[level]));

//       if (AI_measure_array_D && AI_measure_array_D[level])
//       {
//          if (AI_measure_array_F[level])
//          {
//             jxf_TFree(AI_measure_array_F[level]);
//             AI_measure_array_F[level] = NULL;
//          }

//          AI_measure_array_F[level] =
//              jxmp_CopyRealArrayDtoF(AI_measure_array_D[level], local_n);
//       }

//       if (CF_marker_array_D && CF_marker_array_D[level])
//       {
//          if (CF_marker_array_F[level])
//          {
//             jxf_TFree(CF_marker_array_F[level]);
//             CF_marker_array_F[level] = NULL;
//          }

//          CF_marker_array_F[level] =
//              jxmp_CopyIntArrayDtoF(CF_marker_array_D[level], local_n);
//       }

//       if (dof_func_array_D && dof_func_array_D[level])
//       {
//          if (dof_func_array_F[level])
//          {
//             jxf_TFree(dof_func_array_F[level]);
//             dof_func_array_F[level] = NULL;
//          }

//          dof_func_array_F[level] =
//              jxmp_CopyIntArrayDtoF(dof_func_array_D[level], local_n);
//       }
//    }

//    /* level-0 dof_func 与 F 端 solver 主指针保持一致 */
//    if (dof_func_array_F && dof_func_array_F[0])
//    {
//       jxf_ParAMGDataDofFunc(amg_data_F) = dof_func_array_F[0];
//    }

//    return Setup_err_flag;
// }

// JXF_Int
// JX_PAMGPrecond_MxP(JXF_Solver solver,
//                    jxf_hpCSRMatrix par_matrix,
//                    JX_ParVector par_rhs,
//                    JX_ParVector par_app)
// {

//    return (jx_PAMGPrecond_MxP((void *)solver,
//                               (jxf_hpCSRMatrix *)par_matrix,
//                               (jx_ParVector *)par_rhs,
//                               (jx_ParVector *)par_app));
// }

// JXF_Int
// jx_PAMGPrecond_MxP(void *amg_vdata,
//                    jxf_hpCSRMatrix *hp_matrix,
//                    jx_ParVector *par_rhs,
//                    jx_ParVector *par_app)
// {
//    // jx_printf("\n >>> jx_PAMGPrecond_MxP \n");
//    jxf_ParAMGData *amg_data = amg_vdata;
//    jxf_ParVector **F_array = jxf_ParAMGDataFArray(amg_data);
//    jxf_ParVector **U_array = jxf_ParAMGDataUArray(amg_data);

//    JXF_Int max_iter = jxf_ParAMGDataMaxIter(amg_data);
//    JXF_Int max_levels = jxf_ParAMGDataMaxLevels(amg_data);
//    JXF_Int cycle_count = 0;

//    jxf_ParVector *par_rhs_F = NULL;
//    jxf_ParVector *par_app_F = NULL;

//    par_rhs_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
//                                    jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
//                                    jxf_hpCSRMatrixRowStarts(hp_matrix));
//    jxf_ParVectorInitialize(par_rhs_F);
//    jxmp_ParVectorDtoF(par_rhs, par_rhs_F);
//    par_app_F = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix),
//                                    jxf_hpCSRMatrixGlobalNumRows(hp_matrix),
//                                    jxf_hpCSRMatrixRowStarts(hp_matrix));
//    jxf_ParVectorInitialize(par_app_F);
//    jxmp_ParVectorDtoF(par_app, par_app_F);

//    if (max_levels == 1)
//    {
//       jxf_CoarsestSolver(amg_data, hp_matrix, par_rhs_F, par_app_F);
//       return (0);
//    }

//    F_array[0] = par_rhs_F;
//    U_array[0] = par_app_F;

//    /*----------------------------------------------------------------
//     * Generally, max_iter should be set to 1 when AMG is used as
//     * a preconditioner, here, we keep the parameter max_iter for
//     * convenience if max_iter should be larger than 1.
//     *---------------------------------------------------------------*/
//    while (cycle_count < max_iter)
//    {
//       jxf_hpPAMGCycle(amg_data, F_array, U_array);
//       cycle_count++;
//    }

//    jxmp_ParVectorFtoD(U_array[0], par_app);

//    return 0;
// }

// ==============================================================
// typedef struct
// {
//    jx_ParAMGData *amg_data_D;  /* 双精度 AMG 数据 */
//    jxf_ParAMGData *amg_data_F; /* 单精度 AMG 数据 */

//    JX_Int num_levels;
//    JX_Int *level_prec;

// } jx_AdpAMGData;

// typedef jx_AdpAMGData *jx_AdpSolver;

// JX_Int
// JX_AdpAMGCreate(jx_AdpSolver *adp_solver)
// {
//    jx_AdpAMGData *adp_data = NULL;

//    if (adp_solver == NULL)
//    {
//       return -1;
//    }

//    adp_data = (jx_AdpAMGData *)jx_CTAlloc(1, sizeof(jx_AdpAMGData));

//    if (adp_data == NULL)
//    {
//       *adp_solver = NULL;
//       return -1;
//    }

//    adp_data->amg_data_D = NULL;
//    adp_data->amg_data_F = NULL;

//    *adp_solver = adp_data;

//    return 0;
// }

// JX_Int
// JX_AdpAMGSetData(jx_AdpSolver adp_solver,
//                  jx_ParAMGData *amg_data_D,
//                  jxf_ParAMGData *amg_data_F)
// {
//    if (adp_solver == NULL)
//    {
//       return -1;
//    }

//    adp_solver->amg_data_D = amg_data_D;
//    adp_solver->amg_data_F = amg_data_F;

//    return 0;
// }

// JX_Int
// JX_AdpAMGDestroy(jx_AdpSolver adp_solver)
// {
//    if (adp_solver == NULL)
//    {
//       return 0;
//    }

//    jx_TFree(adp_solver);

//    return 0;
// }

// JX_Int
// JX_PAMGPrecond_MxP_V2(jx_AdpSolver adp_solver,
//                       JXF_hpCSRMatrix hp_matrix,
//                       JX_ParVector par_rhs,
//                       JX_ParVector par_app)
// {

//    return (jx_PAMGPrecond_MxP_V2((void *)adp_solver,
//                                  (jxf_hpCSRMatrix *)hp_matrix,
//                                  (jx_ParVector *)par_rhs,
//                                  (jx_ParVector *)par_app));
// }

// JX_Int
// jx_PAMGPrecond_MxP_V2(void *adp_vdata,
//                       jxf_hpCSRMatrix *hp_matrix,
//                       jx_ParVector *par_rhs,
//                       jx_ParVector *par_app)
// {
//    jx_AdpAMGData *adp_data;
//    jx_ParAMGData *amg_data_D;
//    jxf_ParAMGData *amg_data_F;

//    jx_ParVector **F_array_D;
//    jx_ParVector **U_array_D;
//    jxf_ParVector **F_array_F;
//    jxf_ParVector **U_array_F;

//    JX_Int max_iter;
//    JX_Int max_levels;
//    JX_Int cycle_count = 0;

//    if (adp_vdata == NULL)
//    {
//       return -1;
//    }

//    adp_data = (jx_AdpAMGData *)adp_vdata;

//    if (adp_data->amg_data_D == NULL || adp_data->amg_data_F == NULL)
//    {
//       return -1;
//    }

//    amg_data_D = adp_data->amg_data_D;
//    amg_data_F = adp_data->amg_data_F;

//    F_array_D = jx_ParAMGDataFArray(amg_data_D);
//    U_array_D = jx_ParAMGDataUArray(amg_data_D);

//    F_array_F = jxf_ParAMGDataFArray(amg_data_F);
//    U_array_F = jxf_ParAMGDataUArray(amg_data_F);

//    max_iter = jxf_ParAMGDataMaxIter(amg_data_F);
//    max_levels = jxf_ParAMGDataMaxLevels(amg_data_F);

//    if (max_levels == 1)
//    {
//       if (level_prec[level] == JX_PRECISION_FP64)
//       {
//       }
//       else if (level_prec[level] == JX_PRECISION_FP32)
//       {
//       }
//    }

//    /* level 0 绑定双精度向量 */
//    F_array_D[0] = par_rhs;
//    U_array_D[0] = par_app;

//    while (cycle_count < max_iter)
//    {
//       jxmp_hpPAMGCycle_Recursive_adp(adp_data, 0);
//       cycle_count++;
//    }

//    return 0;
// }