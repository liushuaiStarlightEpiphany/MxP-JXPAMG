//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * solver_strong.c -- This is an example demonstrating how to solve a given
 * linear system in the form of files by calling JXFPAMG solver and
 * the JXFPAMG-based Krylov iterative methods.
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

#include "jxf_pamg.h"
#include "jxf_krylov.h"
#include "jxf_diagscale.h"
#include "jxf_blockprec.h"
#include "jxf_euclid.h"
#include "jxf_combined.h"

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

JXF_Int
jxf_ParVectorSetRandomValues(jxf_ParVector *v, JXF_Int seed);
JXF_Int
jxf_SeqVectorSetRandomValues(jxf_Vector *x, JXF_Int seed);

int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   JXF_Int myid, nprocs;
#if JXF_USING_OPENMP || defined(JXF_USING_PGCC_SMP)
   JXF_Int nthreads;
#endif

   JXF_Int arg_index = 0;
   JXF_Int print_usage = 0;

   JXF_Real starttime, endtime;
   JXF_Real starttimeT, endtimeT;

   char *MatFile = NULL;
   char *RhsFile = NULL;
   char *GusFile = NULL;
   char AppFile[120];

   jxf_hpCSRMatrix *hp_matrix = NULL;
   jxf_ParVector *par_rhs = NULL;
   jxf_ParVector *par_sol = NULL;

   JXF_Int *partitioning = NULL;

   JXF_Int **grid_relax_points = NULL;

   JXF_Int build_matrix_arg_index;
   JXF_Int build_rhs_arg_index;
   /* DiagScale Precond */
   JXF_Solver ds_solver;

   /* Euclid solver */
   JXF_Solver euclid_solver;
   JXF_Int euclid_level;
   JXF_Int euclid_bj;

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
   jxf_MPI_Init(&argc, &argv);
   jxf_MPI_Comm_rank(comm, &myid);
   jxf_MPI_Comm_size(comm, &nprocs);
#ifdef USING_HWLOC
   jxf_hpCreateHardwareInfo(comm);
#endif

   //-----------------------
   //  参数设置
   //-----------------------
   build_matrix_arg_index = argc;
   build_rhs_arg_index = argc;
   max_levels = 25;  /* 最大网格层数 */
   cycle_type = 1;   /* Cycle 类型  1: V_Cycle; 2：W_Cycle */
   relax_type = 6;   /* Relax 类型  3: hGS; 6：hSGS */
   measure_type = 0; /* 影响值的计算方式 0：局部；1：全局 */
   rap2 = 0;         /* RAP计算方式  0：RAP；1：先算Q=AP，再算RQ */
   num_functions = 1;
   ns_down = 1;
   ns_up = 1;
   ns_coarse = 1;
   restri_type = 0;         /* 限制算子类型 0: P^T, 1: AIR, 2: AIR-2 */
   keepTranspose = 0;       /* 存放限制算子  0：no；1：yes */
   coarsen_type = 8;        /* 粗化策略 */
   interp_type = 0;         /* 插值策略 */
   P_max_elmts = 0;         /* 插值算子每行最大非零元个数 */
   agg_num_levels = 0;      /* Aggressive粗化的层数 */
   ai_measure_type = 0;     /* AI-策略  0: no; 1: yes */
   ai_relax_type = 0;       /* AI-磨光  0: no; 1: yes */
   amg_print_level = 0;     /* work only when AMG as preconditioner */
   strong_threshold = 0.25; /* 强弱连通参数, 0.25 for 2D, 0.5 for 3D is recommended */
   max_row_sum = 0.9;       /* 行和参数 */
   relax_wt = 1.0;
   outer_wt = 1.0;
   S_commpkg_switch = 1.0;
   AIR_strong_th = 0.25;

   coarse_threshold = 1000; /* 最粗网格层上网格节点个数的最大值 */
   coarse_ratio = 0.75;  /* 相邻两个网格层的粗点个数超过细点个数的 coarse_ratio, 则换成 CLJP 粗化 */
   coarsestsolverid = 9; /* 最粗网格层解法器 */
   conv_criteria = 0;    /* 收敛准则类型 */

   euclid_level = 0; /* level of fill-in */
   euclid_bj = 0;    /* Select PILU (0) or Block Jacobi ILU (1) */

   tol = 1.0e-8;          /* 控制精度 */
   max_iter = 500;         /* 迭代法最大迭代次数 */
   k_dim = 50;             /* 回头数 */
   is_check_restarted = 0; /* peghoty, 2011/11/08 */
   twonorm = 1;            /* PCG 法中的范数控制类型，0: B 范数; 1: l2 范数 */
   print_level = 2;        /* 0: 关闭；1：Setup参数；2：Solve参数；3：Setup+Solve参数 */
   keepsol = 0;            /* 是否保存解向量 */
   TTest = 1;              /* 是否测试时间 */
   solver_id = 22;
   problem_id = -1; /* 默认为-1，需要通过命令行读取矩阵请勿修改此参数*/
   file_base = 0;
#if JXF_USING_OPENMP || defined(JXF_USING_PGCC_SMP)
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
#if JXF_USING_OPENMP || defined(JXF_USING_PGCC_SMP)
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
      jxf_printf("\n");
      jxf_printf("  Usage: %s [<options>]\n", argv[0]);
      jxf_printf("\n");
      jxf_printf("    -bt <val>   : bind type: 0:bind to node, 1:bind to cpu, 2:bind to core\n");
      jxf_printf("    -sid <val>   : solver id\n");
      jxf_printf("    -pid <val>   : problem id\n");
      jxf_printf("  -fromonecsrfile <filename> : ");
      jxf_printf("matrix read from a single file (CSR format)\n");
      jxf_printf("  -rhsfromonefile        : ");
      jxf_printf("rhs read from a single file (CSR format)\n");
      jxf_printf("    -nf  <val>   : number of functions\n");
      jxf_printf("    -nts <val>   : threads number\n");
      jxf_printf("    -rap2        : 2nd implementation of RAP\n");
      jxf_printf("    -kt          : keep transpose\n");
      jxf_printf("    -Pmx <val>   : maximal number of elements per row for interpolation\n");
      jxf_printf(" -agg_nl <val>   : num_levels for aggressive coarsening\n");
      jxf_printf("   -tnrm <val>   : two_norm in PCG\n");
      jxf_printf("   -kdim <val>   : krylov dimension\n");
      jxf_printf("   -dtol <val>   : Drop-tolerance for ILU(0) factorization\n");
      jxf_printf("-euc_lvl <val>   : level of fill-in for Euclid\n");
      jxf_printf(" -euc_bj <val>   : PILU or Block Jacobi ILU\n");
      jxf_printf("  -rlx <val>     : relaxation type\n");
      jxf_printf("  -ai_rlx <val>  : AI relaxation type\n");
      jxf_printf("  -ai_mt  <val>  : AI measure type\n");
      jxf_printf("     -ct <val>   : coarsening type\n");
      jxf_printf("    -ipt <val>   : interpolation type\n");
      jxf_printf("   -mxct <val>   : max. size on coarsest grid\n");
      jxf_printf("    -str <val>   : AMG strength threshold\n");
      jxf_printf("   -mxrs <val>   : maximum row sum threshold for dependency weakening\n");
      jxf_printf("-amg_ptlv <val>  : print_level of AMG when AMG as preconditioner\n");
      jxf_printf("-theta_ms <iiiv> : threshold for multi-scale judgement\n");
      jxf_printf("   -ptlv <val>   : print_level\n");
      jxf_printf("   -help         : using help message\n\n");
      exit(1);
   }

   //-----------------------
   //  提供线性系统文件
   //-----------------------
   switch (problem_id)
   {
   case 1:
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Constant/mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Constant/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 2:
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Anisotropy/mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Anisotropy/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 3:
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e3.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 4:
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Jump/mat_csr_128X128X128_1e6.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Jump/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 5: // 论文测试矩阵
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Rand/mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/Rand/rhs_128X128X128.bin";
      file_base = 2;
      break;
   case 6: // 论文测试矩阵
      MatFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/test/Checkerboard_mat_csr_128X128X128.bin";
      RhsFile = "/home/spring/xingxin/jxfpamg-basic-xingxin/matrix/test/Checkerboard_rhs_128X128X128.bin";
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
      jxf_printf("Error: No filename specified \n");
      exit(1);
   }

   if (build_rhs_arg_index < argc)
   {
      RhsFile = argv[build_rhs_arg_index];
   }
   else if (problem_id < 0)
   {
      jxf_printf("Error: No filename specified \n");
      exit(1);
   }

   if (myid == 0)
   {
      jxf_printf("\n\n+++++++++++++++++++++\n matrix:\n %s\n+++++++++++++++++++++\n", MatFile);
      jxf_printf("\n+++++++++++++++++++++\n rhs:\n %s\n+++++++++++++++++++++\n", RhsFile);
      jxf_printf("\n+++++++++++++++++++++");
#if JXF_USING_BIG_INT
      jxf_printf(" Using BIG_INT,");
#endif
#if JXF_USING_BIG_DOUBLE
      jxf_printf(" BIG_DOUBLE,");
#endif
#if JXF_USING_OPENMP || defined(JXF_USING_PGCC_SMP)
      jxf_printf(" With OpenMP using %d threads,", nthreads);
#endif
      jxf_printf(" MPI using %d processors +++++++++++++++++++++\n\n", nprocs);
   }

   //----------------------------------------------------------------
   // 设定线程数
   //----------------------------------------------------------------
#if JXF_USING_OPENMP || defined(JXF_USING_PGCC_SMP)
   omp_set_num_threads(nthreads);
#endif

   if (restri_type) /* Set Restriction to be AIR */
   {
      interp_type = 100; /* 1-pt Interp */
      relax_type = 0;
      ns_down = 0;
      ns_up = 3;
      grid_relax_points = jxf_CTAlloc(JXF_Int *, 4);
      grid_relax_points[0] = NULL;
      grid_relax_points[1] = jxf_CTAlloc(JXF_Int, ns_down);
      grid_relax_points[2] = jxf_CTAlloc(JXF_Int, ns_up);
      grid_relax_points[3] = jxf_CTAlloc(JXF_Int, ns_coarse);
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
      starttime = jxf_MPI_Wtime();

   hp_matrix = jxf_hpBuildMatParFromOneFile(MatFile, 1, file_base);
   par_rhs = jxf_hpBuildRhsParFromOneFile(RhsFile, hp_matrix, file_base);

   if (!GusFile)
   {
      glosize = jxf_ParVectorGlobalSize(par_rhs);
      partitioning = jxf_ParVectorPartitioning(par_rhs);
      par_sol = jxf_ParVectorCreate(comm, glosize, partitioning);
      jxf_ParVectorSetPartitioningOwner(par_sol, 0);
      jxf_ParVectorInitialize(par_sol);
      if (initguess == 0)
      {
         jxf_ParVectorSetConstantValues(par_sol, 0.0);
      }
      else
      {
         jxf_ParVectorSetRandomValues(par_sol, 22775);
         norm = jxf_ParVectorInnerProd(par_sol, par_sol);
         norm = 1. / sqrt(norm);
         jxf_ParVectorScale(norm, par_sol);
      }
   }
   else
   {
      par_sol = jxf_hpBuildRhsParFromOneFile(GusFile, hp_matrix, 0);
   }

   endtime = jxf_MPI_Wtime();
   jxf_GetWallTime(comm, "BuildParLinearSystem", starttime, endtime, 0, 2);

   //-----------------------------------------------------------
   //  求解线性代数系统
   //-----------------------------------------------------------
   starttimeT = jxf_MPI_Wtime();

   switch (solver_id)
   {
   case 0: /* PAMG */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: PAMG \n\n");

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
      if (restri_type)
      {
         jxf_assert(restri_type >= 0);
         JXF_PAMGSetRestriction(amg_solver, restri_type);
         JXF_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
      }
      JXF_PAMGSetMaxLevels(amg_solver, max_levels);
      JXF_PAMGSetMaxIter(amg_solver, max_iter);
      JXF_PAMGSetNumFunctions(amg_solver, num_functions);
      JXF_PAMGSetCycleType(amg_solver, cycle_type);
      JXF_PAMGSetMeasureType(amg_solver, measure_type);
      JXF_PAMGSetRAP2(amg_solver, rap2);
      JXF_PAMGSetKeepTranspose(amg_solver, keepTranspose);
      JXF_PAMGSetTol(amg_solver, tol);
      JXF_PAMGSetConvCriteria(amg_solver, conv_criteria);
      JXF_PAMGSetCoarsenType(amg_solver, coarsen_type);
      JXF_PAMGSetInterpType(amg_solver, interp_type);
      JXF_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
      JXF_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
      JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
      JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
      JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
      JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
      JXF_PAMGSetPrintLevel(amg_solver, print_level);
      JXF_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
      JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JXF_PAMGSetCoarseRatio(amg_solver, coarse_ratio);
      JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver, outer_wt);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */
      //------------------------------------------------------------
      //    JXF_PAMG Setup
      //------------------------------------------------------------
      if (max_levels != 1)
      {
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
      }

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG Setup", starttime, endtime, 0, 2);
      }

      //------------------------------------------------------------
      //    JXF_PAMG Solve
      //------------------------------------------------------------
      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_PAMGSolve(amg_solver, (JXF_hpCSRMatrix)(hp_matrix), (JXF_ParVector)par_rhs, (JXF_ParVector)par_sol);

      JXF_PAMGGetNumIterations(amg_solver, &num_iterations);
      JXF_PAMGGetFinalRelativeResidualNorm(amg_solver, &final_res_norm);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG Solve", starttime, endtime, 0, 2);
      }

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }
      JXF_PAMGDestroy(amg_solver);
   }
   break;

   case 11: /* CG */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: CG \n\n");

      JXF_PCGCreate(comm, &solver);

      JXF_PCGSetMaxIter(solver, max_iter);
      JXF_PCGSetTol(solver, tol);
      JXF_PCGSetTwoNorm(solver, twonorm); // 0: B 范数； 1：l2 范数
      JXF_PCGSetLogging(solver, 1);
      JXF_PCGSetPrintLevel(solver, print_level);

      JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      JXF_PCGGetNumIterations(solver, &num_iterations);
      JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      JXF_PCGDestroy(solver);
   }
   break;

   case 12: /* PAMG-CG */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: PAMG-CG \n\n");

      starttime = jxf_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
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
      JXF_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
      JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
      JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
      JXF_PAMGSetOuterWt(amg_solver, outer_wt);
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */
      JXF_PCGCreate(comm, &solver);

      JXF_PCGSetMaxIter(solver, max_iter);
      JXF_PCGSetTol(solver, tol);
      JXF_PCGSetTwoNorm(solver, twonorm); // 0: B 范数； 1：l2 范数
      JXF_PCGSetLogging(solver, 1);
      JXF_PCGSetPrintLevel(solver, print_level);

      JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_PAMGPrecond, (JXF_PtrToSolverFcn)JXF_PAMGSetup, amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "PAMG-CG Setup", starttime, endtime, 0, 2);

      starttime = jxf_MPI_Wtime();

      JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "PAMG-CG Solve", starttime, endtime, 0, 2);

      JXF_PCGGetNumIterations(solver, &num_iterations);
      JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      JXF_PAMGDestroy(amg_solver);
      JXF_PCGDestroy(solver);
   }
   break;

   case 13: /* DS-CG */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: DS-CG \n\n");

      starttime = jxf_MPI_Wtime();

      ds_solver = NULL;

      JXF_PCGCreate(comm, &solver);
      JXF_PCGSetMaxIter(solver, max_iter);
      JXF_PCGSetTol(solver, tol);
      JXF_PCGSetTwoNorm(solver, twonorm); // 0: B 范数； 1：l2 范数
      JXF_PCGSetLogging(solver, 1);
      JXF_PCGSetPrintLevel(solver, print_level);

      JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_DiagScalePrecond,
                       (JXF_PtrToSolverFcn)JXF_DiagScaleSetup, ds_solver);

      JXF_DiagScaleSetup(ds_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "DS-CG Setup", starttime, endtime, 0, 2);

      starttime = jxf_MPI_Wtime();

      JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "DS-CG Solve", starttime, endtime, 0, 2);

      JXF_PCGGetNumIterations(solver, &num_iterations);
      JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      JXF_PCGDestroy(solver);
   }
   break;

   case 14: /* Euclid-CG */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: Euclid-CG \n\n");

      starttime = jxf_MPI_Wtime();

      JXF_EuclidCreate(comm, &euclid_solver);
      // JXF_EuclidSetParams(euclid_solver, argc, argv);
      JXF_EuclidSetLevel(euclid_solver, euclid_level);
      JXF_EuclidSetBJ(euclid_solver, euclid_bj);

      JXF_PCGCreate(comm, &solver);
      JXF_PCGSetMaxIter(solver, max_iter);
      JXF_PCGSetTol(solver, tol);
      JXF_PCGSetTwoNorm(solver, twonorm); // 0: B 范数； 1：l2 范数
      JXF_PCGSetLogging(solver, 1);
      JXF_PCGSetPrintLevel(solver, print_level);

      JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_EuclidSolve,
                       (JXF_PtrToSolverFcn)JXF_EuclidSetup, euclid_solver);

      JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "Euclid-CG Setup", starttime, endtime, 0, 2);

      starttime = jxf_MPI_Wtime();

      JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      endtime = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "Euclid-CG Solve", starttime, endtime, 0, 2);

      JXF_PCGGetNumIterations(solver, &num_iterations);
      JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      JXF_EuclidDestroy(euclid_solver);
      JXF_PCGDestroy(solver);
   }
   break;

   case 21: /* GMRES */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: GMRES(%d) \n\n", k_dim);

      JXF_GMRESCreate(comm, &solver);
      JXF_GMRESSetKDim(solver, k_dim);
      JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JXF_GMRESSetMaxIter(solver, max_iter);
      JXF_GMRESSetTol(solver, tol);
      JXF_GMRESSetLogging(solver, 1);
      JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                    (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      JXF_GMRESGetNumIterations(solver, &num_iterations);
      JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      JXF_GMRESDestroy(solver);
   }
   break;

   case 22: /* PAMG-GMRES */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: PAMG-GMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jxf_MPI_Wtime();

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

      JXF_GMRESCreate(comm, &solver);
      JXF_GMRESSetKDim(solver, k_dim);
      JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JXF_GMRESSetMaxIter(solver, max_iter);
      JXF_GMRESSetTol(solver, tol);
      JXF_GMRESSetLogging(solver, 1);
      JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_PAMGPrecond,
                         (JXF_PtrToSolverFcn)JXF_PAMGSetup, amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                    (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-GMRES Solve", starttime, endtime, 0, 2);
      }

      JXF_GMRESGetNumIterations(solver, &num_iterations);
      JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver);
      JXF_GMRESDestroy(solver);
   }
   break;

   case 23: /* DS-GMRES */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: DS-GMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jxf_MPI_Wtime();

      ds_solver = NULL;

      JXF_GMRESCreate(comm, &solver);
      JXF_GMRESSetKDim(solver, k_dim);
      JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JXF_GMRESSetMaxIter(solver, max_iter);
      JXF_GMRESSetTol(solver, tol);
      JXF_GMRESSetLogging(solver, 1);
      JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_DiagScalePrecond,
                         (JXF_PtrToSolverFcn)JXF_DiagScaleSetup, ds_solver);

      JXF_DiagScaleSetup(ds_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                    (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-GMRES Solve", starttime, endtime, 0, 2);
      }

      JXF_GMRESGetNumIterations(solver, &num_iterations);
      JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_GMRESDestroy(solver);
   }
   break;

   case 24: /* Euclid-GMRES */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: Euclid-GMRES(%d) \n\n", k_dim);

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_EuclidCreate(comm, &euclid_solver);
      // JXF_EuclidSetParams(euclid_solver, argc, argv);
      JXF_EuclidSetLevel(euclid_solver, euclid_level);
      JXF_EuclidSetBJ(euclid_solver, euclid_bj);

      JXF_GMRESCreate(comm, &solver);
      JXF_GMRESSetKDim(solver, k_dim);
      JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
      JXF_GMRESSetMaxIter(solver, max_iter);
      JXF_GMRESSetTol(solver, tol);
      JXF_GMRESSetLogging(solver, 1);
      JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */

      JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_EuclidSolve,
                         (JXF_PtrToSolverFcn)JXF_EuclidSetup, euclid_solver);

      JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-GMRES Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                    (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-GMRES Solve", starttime, endtime, 0, 2);
      }

      JXF_GMRESGetNumIterations(solver, &num_iterations);
      JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_EuclidDestroy(euclid_solver);
      JXF_GMRESDestroy(solver);
   }
   break;

   case 31: /* BiCGSTAB */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: BiCGSTAB \n\n");

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_BiCGSTABCreate(comm, &solver);
      JXF_BiCGSTABSetMaxIter(solver, max_iter);
      JXF_BiCGSTABSetTol(solver, tol);
      JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JXF_BiCGSTABSetConvCriteria(solver, 0);
      JXF_BiCGSTABSetLogging(solver, 1);
      JXF_BiCGSTABSetPrintLevel(solver, print_level);

      JXF_BiCGSTABSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "BiCGSTAB Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                       (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "BiCGSTAB Solve", starttime, endtime, 0, 2);
      }

      JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
      JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_BiCGSTABDestroy(solver);
   }
   break;

   case 32: /* PAMG-BiCGSTAB */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: PAMG-BiCGSTAB \n\n");

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_PAMGCreate(&amg_solver);
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
      if (ns_down > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1); /* sweep for "down" */
      if (ns_up > -1)
         JXF_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);   /* sweep for "up" */
      JXF_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);  /* sweep for "coarsest" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1); /* relax_type for "down" */
      JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2); /* relax_type for "up" */
      JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);          /* relax_type for "coarsest" */

      JXF_BiCGSTABCreate(comm, &solver);
      JXF_BiCGSTABSetMaxIter(solver, max_iter);
      JXF_BiCGSTABSetTol(solver, tol);
      JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JXF_BiCGSTABSetConvCriteria(solver, 0);
      JXF_BiCGSTABSetLogging(solver, 1);
      JXF_BiCGSTABSetPrintLevel(solver, print_level);

      JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_PAMGPrecond,
                            (JXF_PtrToSolverFcn)JXF_PAMGSetup, amg_solver);

      JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_BiCGSTABSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-BiCGSTAB Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                       (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-BiCGSTAB Solve", starttime, endtime, 0, 2);
      }

      JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
      JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_PAMGDestroy(amg_solver);
      JXF_BiCGSTABDestroy(solver);
   }
   break;

   case 33: /* DS-BiCGSTAB */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: DS-BiCGSTAB \n\n");

      if (TTest)
         starttime = jxf_MPI_Wtime();

      ds_solver = NULL;

      JXF_BiCGSTABCreate(comm, &solver);
      JXF_BiCGSTABSetMaxIter(solver, max_iter);
      JXF_BiCGSTABSetTol(solver, tol);
      JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JXF_BiCGSTABSetConvCriteria(solver, 0);
      JXF_BiCGSTABSetLogging(solver, 1);
      JXF_BiCGSTABSetPrintLevel(solver, print_level);

      JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_DiagScalePrecond,
                            (JXF_PtrToSolverFcn)JXF_DiagScaleSetup, ds_solver);

      JXF_DiagScaleSetup(ds_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_BiCGSTABSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-BiCGSTAB Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                       (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-BiCGSTAB Solve", starttime, endtime, 0, 2);
      }

      JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
      JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_BiCGSTABDestroy(solver);
   }
   break;

   case 34: /* Euclid-BiCGSTAB */
   {
      if (myid == 0)
         jxf_printf("\n >>> Solver: Euclid-BiCGSTAB \n\n");

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_EuclidCreate(comm, &euclid_solver);
      // JXF_EuclidSetParams(euclid_solver, argc, argv);
      JXF_EuclidSetLevel(euclid_solver, euclid_level);
      JXF_EuclidSetBJ(euclid_solver, euclid_bj);

      JXF_BiCGSTABCreate(comm, &solver);
      JXF_BiCGSTABSetMaxIter(solver, max_iter);
      JXF_BiCGSTABSetTol(solver, tol);
      JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
      JXF_BiCGSTABSetConvCriteria(solver, 0);
      JXF_BiCGSTABSetLogging(solver, 1);
      JXF_BiCGSTABSetPrintLevel(solver, print_level);

      JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_EuclidSolve,
                            (JXF_PtrToSolverFcn)JXF_EuclidSetup, euclid_solver);

      JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_matrix);

      JXF_BiCGSTABSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-BiCGSTAB Setup", starttime, endtime, 0, 2);
      }

      if (TTest)
         starttime = jxf_MPI_Wtime();

      JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                       (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);

      if (TTest)
      {
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-BiCGSTAB Solve", starttime, endtime, 0, 2);
      }

      JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
      JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (print_level == 0 && myid == 0)
      {
         jxf_printf(" >>> num_iterations = %d\n", num_iterations);
         jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
      }

      JXF_EuclidDestroy(euclid_solver);
      JXF_BiCGSTABDestroy(solver);
   }
   break;
   }

   if (TTest)
   {
      endtimeT = jxf_MPI_Wtime();
      jxf_GetWallTime(comm, "Total Solve Time", starttimeT, endtimeT, 0, 6);
   }

   // rel-residual
   norm2 = jxf_ParVectorNorm2(par_rhs);
   jxf_hpCSRMatrixMatvec(1.0, hp_matrix, par_sol, -1.0, par_rhs);
   norm1 = jxf_ParVectorNorm2(par_rhs);
   norm1 = norm1 / norm2;

   if (myid == 0)
   {
      jxf_printf("\n true rel residual is %.4le\n", norm1);
   }

   //-----------------------------------------------------------
   //  将近似解向量保存到文件中
   //-----------------------------------------------------------
   if (keepsol)
   {
      jxf_Vector *ser_sol = NULL;
      ser_sol = jxf_ParVectorToVectorAll(par_sol);
      if (myid == 0)
      {
         jxf_sprintf(AppFile, "%s%06d%s%06d%s%03d%s%03d",
                    "app_", nprocs, "_", problem_id, "_", solver_id, "_", relax_type);
         jxf_SeqVectorPrint(ser_sol, AppFile);
      }
      jxf_SeqVectorDestroy(ser_sol);
   }

//-----------------------------------------------------------
//  销毁并行矩阵和并行向量
//-----------------------------------------------------------
#ifdef USING_HWLOC
   jxf_hpCSRhardwareDestroy();
#endif
   jxf_hpCSRMatrixDestroy(hp_matrix);
   jxf_ParVectorDestroy(par_rhs);
   jxf_ParVectorDestroy(par_sol);

   //-----------------------------------------------------------
   //  终止 MPI
   //-----------------------------------------------------------
   jxf_MPI_Finalize();

   return 0;
}

/*!
 * \fn JXF_Int jxf_ParVectorSetRandomValues
 */
JXF_Int
jxf_ParVectorSetRandomValues(jxf_ParVector *v, JXF_Int seed)
{
   JXF_Int my_id;
   jxf_Vector *v_local = jxf_ParVectorLocalVector(v);

   MPI_Comm comm = jxf_ParVectorComm(v);
   jxf_MPI_Comm_rank(comm, &my_id);

   seed *= (my_id + 1);

   return jxf_SeqVectorSetRandomValues(v_local, seed);
}

/*!
 * \fn JXF_Int jxf_SeqVectorSetRandomValues
 */
JXF_Int
jxf_SeqVectorSetRandomValues(jxf_Vector *x, JXF_Int seed)
{
   JXF_Real *vector_data = jxf_VectorData(x);
   JXF_Int size = jxf_VectorSize(x);
   JXF_Int i, ierr = 0;

   jxf_SeedRand(seed);

   size *= jxf_VectorNumVectors(x);

   /* RDF: threading this loop may cause problems because of jxf_Rand() */
   for (i = 0; i < size; i++)
   {
      vector_data[i] = 2.0 * jxf_Rand() - 1.0;
   }

   return ierr;
}
