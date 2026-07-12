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

#include "solver_src.h"

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

   /* adaptive precision package */
   jxmp_AdpAMGData adp_data;
   JX_Int *level_prec = NULL;

   /* JXFPAMG solver - double / float */
   JX_Solver amg_solver_D;
   JXF_Solver amg_solver_F;

   /* shared AMG parameters */
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
   JX_Real strong_threshold;
   JX_Real max_row_sum;

   JX_Real norm1, norm2;
   JX_Real relax_wt;
   JX_Real outer_wt;
   JX_Real S_commpkg_switch;
   JX_Real AIR_strong_th;

   JX_Int coarse_threshold;
   JX_Real coarse_ratio;
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
   max_iter = 50;          /* 迭代法最大迭代次数 */
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
      MatFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_A.bin";
      RhsFile = "/vol8/home/xtu_ljc/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_b.bin";
      file_base = 2;
      break;
   case 2:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/mat_csr_128X128X128.bin";
      RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/rhs_128X128X128.bin";
      file_base = 2;
      break;
      // case 3:
      //    MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e3.bin";
      //    RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      //    file_base = 2;
      //    break;
      // case 4:
      //    MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e6.bin";
      //    RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      //    file_base = 2;
      //    break;
      // case 5:
      //    MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Rand/mat_csr_128X128X128.bin";
      //    RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/Rand/rhs_128X128X128.bin";
      //    file_base = 2;
      //    break;
      // case 6:
      //    MatFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/test/Checkerboard_mat_csr_128X128X128.bin";
      //    RhsFile = "/home/spring/xingxin/jxpamg-basic-xingxin/matrix/test/Checkerboard_rhs_128X128X128.bin";
      //    file_base = 2;
      //    break;
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

   par_rhs_F = jxf_ParVectorCreate(comm, glosize, partitioning);
   jxf_ParVectorSetPartitioningOwner(par_rhs_F, 0);
   jxf_ParVectorInitialize(par_rhs_F);
   jxmp_ParVectorDtoF(par_rhs, par_rhs_F);

   par_sol_F = jxf_ParVectorCreate(comm, glosize, partitioning);
   jxf_ParVectorSetPartitioningOwner(par_sol_F, 0);
   jxf_ParVectorInitialize(par_sol_F);
   jxmp_ParVectorDtoF(par_sol, par_sol_F);
   // jxf_ParVectorSetConstantValues(par_sol_F, 0.0);

   //--------------------------------------------------------------------
   //  求解线性代数系统
   //--------------------------------------------------------------------
   starttimeT = jx_MPI_Wtime();

   switch (solver_id)
   {
   case 22: /* Adaptive-Precision PAMG-GMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: Adaptive-Precision PAMG-GMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      /* ---------------------------------------------------------
       * 1. 创建高精度 AMG solver_D
       * --------------------------------------------------------- */
      JX_PAMGCreate(&amg_solver_D);

      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver_D, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver_D, grid_relax_points);
      }

      JX_PAMGSetMaxLevels(amg_solver_D, max_levels);
      JX_PAMGSetMaxIter(amg_solver_D, 1);
      JX_PAMGSetCycleType(amg_solver_D, cycle_type);
      JX_PAMGSetMeasureType(amg_solver_D, measure_type);
      JX_PAMGSetRAP2(amg_solver_D, rap2);
      JX_PAMGSetKeepTranspose(amg_solver_D, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver_D, coarsen_type);
      JX_PAMGSetInterpType(amg_solver_D, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver_D, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver_D, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver_D, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver_D, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver_D, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver_D, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver_D, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver_D, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver_D, relax_wt);
      JX_PAMGSetOuterWt(amg_solver_D, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver_D, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver_D, AIR_strong_th);

      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_down, 1);
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_up, 2);
      JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_coarse, 3);
      JX_PAMGSetCycleRelaxType(amg_solver_D, relax_type, 1);
      JX_PAMGSetCycleRelaxType(amg_solver_D, relax_type, 2);
      JX_PAMGSetCycleRelaxType(amg_solver_D, 9, 3);

      /* ---------------------------------------------------------
       * 2. 创建低精度 AMG solver_F
       * --------------------------------------------------------- */
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
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_down, 1);
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1);
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2);
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3);

      /* ---------------------------------------------------------
       * 3. 分别做高/低精度 AMG setup
       * --------------------------------------------------------- */
      JX_PAMGSetup(amg_solver_D, (JX_hpCSRMatrix)hp_matrix);
      JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

      JX_Int num_levels_D = jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver_D);
      JX_Int num_levels_F = jxf_ParAMGDataNumLevels((jxf_ParAMGData *)amg_solver_F);

      if (num_levels_D != num_levels_F)
      {
         if (myid == 0)
            jx_printf("[ERROR] AMG hierarchy mismatch: num_levels_D=%d, num_levels_F=%d\n",
                      num_levels_D, num_levels_F);

         JXF_PAMGDestroy(amg_solver_F);
         JX_PAMGDestroy(amg_solver_D);
         return -1;
      }

      JX_Int actual_num_levels = num_levels_D;
      level_prec = jx_CTAlloc(JX_Int, actual_num_levels);

      for (i = 0; i < actual_num_levels; i++)
         level_prec[i] = JX_PRECISION_FP32;

      /* ---------------------------------------------------------
       * 6. 创建 GMRES 外层
       * --------------------------------------------------------- */
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level);

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                         (JX_PtrToSolverFcn)JXF_PAMGSetup, (JXF_Solver)amg_solver_F);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "Adaptive PAMG-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver, (JX_Matrix)hp_matrix, // preOperater
                    (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "Adaptive PAMG-GMRES Solve", starttime, endtime, 0, 2);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      jx_TFree(level_prec);
      level_prec = NULL;

      JX_GMRESDestroy(solver);
      JXF_PAMGDestroy(amg_solver_F);
      JX_PAMGDestroy(amg_solver_D);
   }
   break;

   case 2201: /* Adaptive-Precision PAMG-GMRES (FP64-only isolation test) */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: Adaptive-Precision PAMG-GMRES(%d) [FP64-only isolation] \n\n", k_dim);

      if (TTest)
         starttime = jx_MPI_Wtime();

      /* ---------------------------------------------------------
       * 1. 创建高精度 AMG solver_D
       * --------------------------------------------------------- */
      JX_PAMGCreate(&amg_solver_D);

      if (restri_type)
      {
         jx_assert(restri_type >= 0);
         JX_PAMGSetRestriction(amg_solver_D, restri_type);
         JX_PAMGSetGridRelaxPoints(amg_solver_D, grid_relax_points);
      }

      JX_PAMGSetMaxLevels(amg_solver_D, max_levels);
      JX_PAMGSetMaxIter(amg_solver_D, 1);
      JX_PAMGSetCycleType(amg_solver_D, cycle_type);
      JX_PAMGSetMeasureType(amg_solver_D, measure_type);
      JX_PAMGSetRAP2(amg_solver_D, rap2);
      JX_PAMGSetKeepTranspose(amg_solver_D, keepTranspose);
      JX_PAMGSetCoarsenType(amg_solver_D, coarsen_type);
      JX_PAMGSetInterpType(amg_solver_D, interp_type);
      JX_PAMGSetPMaxElmts(amg_solver_D, P_max_elmts);
      JX_PAMGSetAggNumLevels(amg_solver_D, agg_num_levels);
      JX_PAMGSetAIMeasureType(amg_solver_D, ai_measure_type);
      JX_PAMGSetAIRelaxType(amg_solver_D, ai_relax_type);
      JX_PAMGSetStrongThreshold(amg_solver_D, strong_threshold);
      JX_PAMGSetMaxRowSum(amg_solver_D, max_row_sum);
      JX_PAMGSetPrintLevel(amg_solver_D, amg_print_level);
      JX_PAMGSetCoarseThreshold(amg_solver_D, coarse_threshold);
      JX_PAMGSetRelaxWt(amg_solver_D, relax_wt);
      JX_PAMGSetOuterWt(amg_solver_D, outer_wt);
      JX_PAMGSetSCommPkgSwitch(amg_solver_D, S_commpkg_switch);
      JX_PAMGSetAIRStrongTh(amg_solver_D, AIR_strong_th);

      if (ns_down > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_down, 1);
      if (ns_up > -1)
         JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_up, 2);
      JX_PAMGSetCycleNumSweeps(amg_solver_D, ns_coarse, 3);
      JX_PAMGSetCycleRelaxType(amg_solver_D, relax_type, 1);
      JX_PAMGSetCycleRelaxType(amg_solver_D, relax_type, 2);
      JX_PAMGSetCycleRelaxType(amg_solver_D, 9, 3);

      /* ---------------------------------------------------------
       * 2. 只做高精度 AMG setup
       * --------------------------------------------------------- */
      JX_PAMGSetup(amg_solver_D, (JX_hpCSRMatrix)hp_matrix);

      /* ---------------------------------------------------------
       * 3. 构造 level_prec：全部 FP64
       * --------------------------------------------------------- */
      JX_Int actual_num_levels = jx_ParAMGDataNumLevels((jx_ParAMGData *)amg_solver_D);

      level_prec = jx_CTAlloc(JX_Int, actual_num_levels);
      for (i = 0; i < actual_num_levels; i++)
         level_prec[i] = JX_PRECISION_FP64;

      /* ---------------------------------------------------------
       * 4. 绑定 adp 数据包：F 侧先置空
       * --------------------------------------------------------- */
      adp_data.amg_data_D = (jx_ParAMGData *)amg_solver_D;
      adp_data.amg_data_F = NULL;
      adp_data.level_prec = level_prec;

      /* ---------------------------------------------------------
       * 5. 创建 GMRES 外层
       * --------------------------------------------------------- */
      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted);
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level);

      JX_GMRESSetPrecond(solver,
                         (JX_PtrToSolverFcn)JXMP_PAMGPrecond_adp,
                         (JX_PtrToSolverFcn)JXMP_PAMGSetup_adp,
                         (JX_Solver)&adp_data);

      JX_GMRESSetup(solver, (JX_Matrix)hp_matrix, (JX_Vector)par_rhs, (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "Adaptive PAMG-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jx_MPI_Wtime();

      JX_GMRESSolve(solver,
                    (JX_Matrix)hp_matrix,
                    (JX_Matrix)hp_matrix,
                    (JX_Vector)par_rhs,
                    (JX_Vector)par_sol);

      if (TTest)
      {
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "Adaptive PAMG-GMRES Solve", starttime, endtime, 0, 2);
      }

      JX_GMRESGetNumIterations(solver, &num_iterations);
      JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jx_printf(" >>> num_iterations = %d\n", num_iterations);
         jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JX_GMRESDestroy(solver);
      jx_TFree(level_prec);
      level_prec = NULL;
      JX_PAMGDestroy(amg_solver_D);
   }
   break;

   case 2202: /* PAMG-GMRES */
   {
      if (myid == 0)
         jx_printf("\n >>> Solver: PAMG-GMRES(%d) \n\n", k_dim);

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
         JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_up, 2);                          /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver_F, ns_coarse, 3);                         /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 1);                        /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, relax_type, 2);                        /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver_F, 9, 3); /* relax_type for "coarsest" */ /* relax_type for "coarsest" */

      JX_GMRESCreate(comm, &solver);
      JX_GMRESSetKDim(solver, k_dim);
      JX_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JX_GMRESSetMaxIter(solver, max_iter);
      JX_GMRESSetTol(solver, tol);
      JX_GMRESSetLogging(solver, 1);
      JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond_MxP,
                         (JX_PtrToSolverFcn)JXF_PAMGSetup, (JXF_Solver)amg_solver_F);

      JXF_PAMGSetup(amg_solver_F, (JXF_hpCSRMatrix)hp_matrix_F);

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

      JXF_PAMGDestroy(amg_solver_F);
      JX_GMRESDestroy(solver);
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
