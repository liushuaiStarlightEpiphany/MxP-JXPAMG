//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * solver_strong_ij.c -- This is an example demonstrating how to solve a given
 * linear system in the form of files by calling JXPAMG solver and
 * the JXPAMG-based Krylov iterative methods.
 *
 * One of the following solvers can be used by
 * assigning parameter 'solver_id'
 *
 *  solver_id = 0:  PAMG
 *  solver_id = 11: CG
 *  solver_id = 12: PAMG-CG
 *  solver_id = 21: GMRES
 *  solver_id = 22: PAMG-GMRES
 *  solver_id = 31: BiCGSTAB
 *  solver_id = 32: PAMG-BiCGSTAB
 *
 *  Created by Yue Xiaoqiang 2022/08/15
 *
 *  Xiangtan University
 *  yuexq@xtu.edu.cn
 *
 */

#include "jx_pamg.h"
#include "jx_krylov.h"
#include "jx_diagscale.h"
#include "jx_euclid.h"
#include "jx_combined.h"

int
main( int argc, char *argv[] )
{
   MPI_Comm  comm = MPI_COMM_WORLD;
   JX_Int       myid, nprocs;
#if JX_USING_OPENMP || defined (JX_USING_PGCC_SMP)
   JX_Int       nthreads;
#endif
   
   JX_Int arg_index   = 0;
   JX_Int print_usage = 0;
   
   JX_Real starttime, endtime;
   JX_Real starttimeT, endtimeT;
   
   char *MatFile = NULL;
   char *RhsFile = NULL;
   
   JX_IJMatrix  mat_ij;
   JX_IJVector  vec_ij;
   JX_IJVector  sol_ij;

   JX_hpCSRMatrix   mat_csr;
   JX_ParVector     vec_csr;
   JX_ParVector     sol_csr;

   jx_CSRMatrix *ser_mat_csr;
   jx_Vector    *ser_vec_csr;

   JX_Int  first_local_row;
   JX_Int  last_local_row;
   JX_Int  first_local_col;
   JX_Int  last_local_col;
   JX_Int  local_num_cols;

   void  *object;

   JX_Real  *values;

   JX_Int **grid_relax_points = NULL;
   
   /* JXPAMG solver */
   JX_Solver    amg_solver;
   JX_Int       max_levels;
   JX_Int       cycle_type;
   JX_Int       relax_type;
   JX_Int       measure_type;
   JX_Int       rap2;
   JX_Int       num_functions;
   JX_Int       ns_down;
   JX_Int       ns_up;
   JX_Int       ns_coarse;
   JX_Int       restri_type;
   JX_Int       keepTranspose;
   JX_Int       coarsen_type;
   JX_Int       interp_type;
   JX_Int       P_max_elmts;
   JX_Int       agg_num_levels;
   JX_Int       ai_measure_type;
   JX_Int       ai_relax_type;
   JX_Real      strong_threshold;
   JX_Real      max_row_sum;
   
   JX_Real    relax_wt;
   JX_Real    outer_wt;
   
   JX_Real S_commpkg_switch;
   JX_Real AIR_strong_th;

   JX_Int       coarse_threshold;
   JX_Real      coarse_ratio;
   JX_Int       coarsestsolverid;
   JX_Int       conv_criteria;
   JX_Int       amg_print_level;
   
   /* iterative method */
   JX_Solver    solver;
   JX_Real      tol;
   JX_Int       max_iter;
   JX_Int       k_dim;
   JX_Int       is_check_restarted;  /* peghoty, 2011/11/08 */
   JX_Int       twonorm;
   JX_Int       solver_id;
   JX_Int       problem_id;
   JX_Int       print_level;
   JX_Int       print_system;
   JX_Int       TTest;
   
   /* other variables */
   JX_Int       i;
   JX_Int       num_iterations;
   JX_Real      final_res_norm;
   
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
   max_levels       = 25;       /* 最大网格层数 */
   cycle_type       = 1;        /* Cycle 类型  1: V_Cycle; 2：W_Cycle */
   relax_type       = 3;        /* Relax 类型  3: hGS; 6：hSGS */
   measure_type     = 0;        /* 影响值的计算方式 0：局部；1：全局 */
   rap2             = 0;        /* RAP计算方式  0：RAP；1：先算Q=AP，再算RQ */
   num_functions    = 1;
   ns_down          = 1;
   ns_up            = 1;
   ns_coarse        = 1;
   restri_type      = 0;        /* 限制算子类型 0: P^T, 1: AIR, 2: AIR-2 */
   keepTranspose    = 0;        /* 存放限制算子  0：no；1：yes */
   coarsen_type     = 6;        /* 粗化策略 */
   interp_type      = 0;        /* 插值策略 */
   P_max_elmts      = 0;        /* 插值算子每行最大非零元个数 */
   agg_num_levels   = 0;        /* Aggressive粗化的层数 */
   ai_measure_type  = 0;        /* AI-策略  0: no; 1: yes */
   ai_relax_type    = 0;        /* AI-磨光  0: no; 1: yes */
   amg_print_level  = 3;        /* work only when AMG as preconditioner */
   strong_threshold = 0.25;     /* 强弱连通参数, 0.25 for 2D, 0.5 for 3D is recommended */
   max_row_sum      = 0.9;      /* 行和参数 */
   relax_wt         = 1.0;
   outer_wt         = 1.0;
   S_commpkg_switch = 1.0;
   AIR_strong_th    = 0.25;
   
   coarse_threshold = 9;      /* 最粗网格层上网格节点个数的最大值 */
   coarse_ratio     = 0.75;     /* 相邻两个网格层的粗点个数超过细点个数的 coarse_ratio, 则换成 CLJP 粗化 */
   coarsestsolverid = 9;        /* 最粗网格层解法器 */
   conv_criteria    = 0;        /* 收敛准则类型 */
   
   tol                 = 1.0e-8; /* 控制精度 */
   max_iter            = 200;    /* 迭代法最大迭代次数 */
   k_dim               = 30;     /* 回头数 */
   is_check_restarted  = 1;      /* peghoty, 2011/11/08 */
   twonorm             = 1;      /* PCG 法中的范数控制类型，0: B 范数; 1: l2 范数 */
   print_level         = 3;      /* 0: 关闭；1：Setup参数；2：Solve参数；3：Setup+Solve参数 */
   print_system        = 0;      /* 是否打印线性系统 */
   TTest               = 1;      /* 是否测试时间 */
   solver_id           = 12;
   problem_id          = 14;
#if JX_USING_OPENMP || defined (JX_USING_PGCC_SMP)
   nthreads            = 1;      /* 线程数 */
#endif
   
   //-----------------------
   //  命令行修改参数
   //-----------------------
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-sid") == 0 )
      {
         arg_index ++;
         solver_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-pid") == 0 )
      {
         arg_index ++;
         problem_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mxlvl") == 0 )
      {
         arg_index ++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ctsid") == 0 )
      {
         arg_index ++;
         coarsestsolverid = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nf") == 0 )
      {
         arg_index ++;
         num_functions = atoi(argv[arg_index++]);
      }
#if JX_USING_OPENMP || defined (JX_USING_PGCC_SMP)
      else if ( strcmp(argv[arg_index], "-nts") == 0 )
      {
         arg_index ++;
         nthreads = atoi(argv[arg_index++]);
      }
#endif
      else if ( strcmp(argv[arg_index], "-rap2") == 0 )
      {
         arg_index ++;
         rap2 = 1;
      }
      else if ( strcmp(argv[arg_index], "-air") == 0 )
      {
         arg_index ++;
         restri_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-kt") == 0 )
      {
         arg_index ++;
         keepTranspose = 1;
      }
      else if (strcmp(argv[arg_index], "-Pmx") == 0)
      {
         arg_index ++;
         P_max_elmts = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-agg_nl") == 0 )
      {
         arg_index ++;
         agg_num_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tnrm") == 0 )
      {
         arg_index ++;
         twonorm = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-kdim") == 0 )
      {
         arg_index ++;
         k_dim = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rlx") == 0 )
      {
         arg_index ++;
         relax_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ai_rlx") == 0 )
      {
         arg_index ++;
         ai_relax_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ai_mt") == 0 )
      {
         arg_index ++;
         ai_measure_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ct") == 0 )
      {
         arg_index ++;
         coarsen_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ipt") == 0 )
      {
         arg_index ++;
         interp_type = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mxct") == 0 )
      {
         arg_index ++;
         coarse_threshold = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-str") == 0 )
      {
         arg_index ++;
         strong_threshold = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mxrs") == 0 )
      {
         arg_index ++;
         max_row_sum = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-scs") == 0 )
      {
         arg_index ++;
         S_commpkg_switch = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-airst") == 0 )
      {
         arg_index ++;
         AIR_strong_th = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-amg_ptlv") == 0 )
      {
         arg_index ++;
         amg_print_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ptlv") == 0 )
      {
         arg_index ++;
         print_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage = 1;
         break;
      }
      else
      {
         arg_index ++;
      }
   }
   if (print_usage)
   {
      jx_printf("\n");
      jx_printf("  Usage: %s [<options>]\n", argv[0]);
      jx_printf("\n");
      jx_printf("    -sid <val>   : solver id\n");
      jx_printf("    -pid <val>   : problem id\n");
      jx_printf("    -nf  <val>   : number of functions\n");
      jx_printf("    -nts <val>   : threads number\n");
      jx_printf("    -rap2        : 2nd implementation of RAP\n");
      jx_printf("    -kt          : keep transpose\n");
      jx_printf("    -Pmx <val>   : maximal number of elements per row for interpolation\n");
      jx_printf(" -agg_nl <val>   : num_levels for aggressive coarsening\n");
      jx_printf("   -tnrm <val>   : two_norm in PCG\n");
      jx_printf("   -kdim <val>   : krylov dimension\n");
      jx_printf("  -rlx <val>     : relaxation type\n");
      jx_printf("  -ai_rlx <val>  : AI relaxation type\n");
      jx_printf("  -ai_mt  <val>  : AI measure type\n");
      jx_printf("     -ct <val>   : coarsening type\n");
      jx_printf("    -ipt <val>   : interpolation type\n");
      jx_printf("   -mxct <val>   : max. size on coarsest grid\n");
      jx_printf("    -str <val>   : AMG strength threshold\n");
      jx_printf("   -mxrs <val>   : maximum row sum threshold for dependency weakening\n");
      jx_printf("-amg_ptlv <val>  : print_level of AMG when AMG as preconditioner\n");
      jx_printf("   -ptlv <val>   : print_level\n");
      jx_printf("   -help         : using help message\n\n");
      exit(1);
   }
   
   //-----------------------
   //  提供线性系统文件
   //-----------------------
   switch (problem_id)
   {
      case 0:
      MatFile = "/home/matrix/iter00001/A/job0/A.iter00001.job0.P.mat";
      RhsFile = "/home/matrix/iter00001/b/job0/b.iter00001.job0.P.vec";
      break;
      case 1:
      MatFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/A.iter00001.job0.P.mat";
      RhsFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/b.iter00001.job0.P.vec";
      break;
      case 2:
      MatFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/A.iter00001.job1.P.mat";
      RhsFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/b.iter00001.job1.P.vec";
      break;
      case 3:
      MatFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/A.iter15000.job0.P.mat";
      RhsFile = "/home/yuexq/data/comac/Axb_LDI_struct_94w/b.iter15000.job0.P.vec";
      break;
   }
   
   if (myid == 0)
   {
      jx_printf("\n\n+++++++++++++++++++++ Problem %d,", problem_id);
#if JX_USING_BIG_INT
      jx_printf(" Using BIG_INT,");
#endif
#if JX_USING_BIG_DOUBLE
      jx_printf(" BIG_DOUBLE,");
#endif
#if JX_USING_OPENMP || defined (JX_USING_PGCC_SMP)
      jx_printf(" With OpenMP using %d threads,", nthreads);
#endif
      jx_printf(" MPI using %d processors +++++++++++++++++++++\n\n", nprocs);
   }
   
   //----------------------------------------------------------------
   // 设定线程数
   //----------------------------------------------------------------
#if JX_USING_OPENMP || defined (JX_USING_PGCC_SMP)
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
      for (i = 0; i < ns_down; i ++) grid_relax_points[1][i] = 0; /* down cycle */
      if (ns_up == 3) /* up cycle */
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
      for (i = 0; i < ns_coarse; i ++) grid_relax_points[3][i] = 0; /* coarse: all */
      coarse_threshold = 100;
      agg_num_levels = 0; /* does not support aggressive coarsening */
   }

   //--------------------------------------------------------------------
   //  利用文件创建并行矩阵和并行右端, 以及并行初始迭代向量(零向量或随机向量)
   //--------------------------------------------------------------------
   if (TTest) starttime = jx_MPI_Wtime();
   
   if (myid == 0)
   {
      jx_printf(" IJ matrix read from file %s\n", MatFile);
   }
   JX_IJMatrixRead(MatFile, comm, JX_HPCSR, &mat_ij);
   JX_IJMatrixGetLocalRange(mat_ij, &first_local_row, &last_local_row, &first_local_col, &last_local_col);
   local_num_cols = last_local_col - first_local_col + 1;
   JX_IJMatrixGetObject(mat_ij, &object);
   mat_csr = (JX_hpCSRMatrix)object;

   if (myid == 0)
   {
      jx_printf(" RHS vector read from file %s\n", RhsFile);
   }
   JX_IJVectorRead(RhsFile, comm, JX_HPCSR, &vec_ij);
   JX_IJVectorGetObject(vec_ij, &object);
   vec_csr = (JX_ParVector)object;

   JX_IJVectorCreate(comm, first_local_col, last_local_col, &sol_ij);
   JX_IJVectorSetObjectType(sol_ij, JX_HPCSR);
   JX_IJVectorInitialize(sol_ij);
   values = jx_CTAlloc(JX_Real, local_num_cols);
   for (i = 0; i < local_num_cols; i ++) values[i] = 0.0;
   JX_IJVectorSetValues(sol_ij, local_num_cols, NULL, values);
   jx_TFree(values);
   JX_IJVectorGetObject(sol_ij, &object);
   sol_csr = (JX_ParVector)object;
   
   endtime = jx_MPI_Wtime();
   jx_GetWallTime(comm, "BuildParLinearSystem", starttime, endtime, 0, 2);

   if (TTest) starttime = jx_MPI_Wtime();
   if (print_system)
   {
      ser_mat_csr = jx_ParCSRMatrixToCSRMatrixAll((jx_ParCSRMatrix *) mat_csr);
      if (myid == 0)
      {
         jx_CSRMatrixPrint(ser_mat_csr, "./matA");
      }
      jx_CSRMatrixDestroy(ser_mat_csr);
      ser_vec_csr = jx_ParVectorToVectorAll((jx_ParVector *) vec_csr);
      if (myid == 0)
      {
         jx_SeqVectorPrint(ser_vec_csr, "./vecB");
      }
      jx_SeqVectorDestroy(ser_vec_csr);
   }
   endtime = jx_MPI_Wtime();
   jx_GetWallTime(comm, "PrintLinearSystem", starttime, endtime, 0, 2);

   //-----------------------------------------------------------
   //  求解线性代数系统
   //-----------------------------------------------------------
   starttimeT = jx_MPI_Wtime();

   switch (solver_id)
   {
      case 0:  /* PAMG */
      {
         if (myid == 0) jx_printf("\n >>> Solver: PAMG \n\n");
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_PAMGCreate(&amg_solver);
         if (restri_type)
         {
            jx_assert(restri_type >= 0);
            JX_PAMGSetRestriction(amg_solver, restri_type);
            JX_PAMGSetGridRelaxPoints(amg_solver, grid_relax_points);
         }
         JX_PAMGSetMaxLevels(amg_solver, max_levels);
         JX_PAMGSetMaxIter(amg_solver, max_iter);
         JX_PAMGSetNumFunctions(amg_solver, num_functions);
         JX_PAMGSetCycleType(amg_solver, cycle_type);
         JX_PAMGSetMeasureType(amg_solver, measure_type);
         JX_PAMGSetRAP2(amg_solver, rap2);
         JX_PAMGSetKeepTranspose(amg_solver, keepTranspose);
         JX_PAMGSetTol(amg_solver, tol);
         JX_PAMGSetConvCriteria(amg_solver, conv_criteria);
         JX_PAMGSetCoarsenType(amg_solver, coarsen_type);
         JX_PAMGSetInterpType(amg_solver, interp_type);
         JX_PAMGSetPMaxElmts(amg_solver, P_max_elmts);
         JX_PAMGSetAggNumLevels(amg_solver, agg_num_levels);
         JX_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
         JX_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
         JX_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JX_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JX_PAMGSetPrintLevel(amg_solver, print_level);
         JX_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetCoarseRatio(amg_solver, coarse_ratio);
         JX_PAMGSetRelaxWt(amg_solver, relax_wt);
         JX_PAMGSetOuterWt(amg_solver, outer_wt);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);           /* sweep for "down" */
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);           /* sweep for "up" */
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);           /* sweep for "coarsest" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);           /* relax_type for "coarsest" */
         //------------------------------------------------------------
         //    JX_PAMG Setup
         //------------------------------------------------------------
         if (max_levels != 1)
         {
            JX_PAMGSetup(amg_solver, mat_csr);
         }
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG Setup", starttime, endtime, 0, 2);
         }

         //------------------------------------------------------------
         //    JX_PAMG Solve
         //------------------------------------------------------------
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_PAMGSolve(amg_solver, mat_csr, vec_csr, sol_csr);
         
         JX_PAMGGetNumIterations(amg_solver, &num_iterations);
         JX_PAMGGetFinalRelativeResidualNorm(amg_solver, &final_res_norm);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG Solve", starttime, endtime, 0, 2);
         }
         
         if (print_level == 0 && myid == 0)
         {
            jx_printf(" >>> num_iterations = %d\n", num_iterations);
            jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JX_PAMGDestroy(amg_solver);
      }
      break;

      case 11:  /* CG */
      {
         if (myid == 0) jx_printf("\n >>> Solver: CG \n\n");
         
         JX_PCGCreate(comm, &solver);
         
         JX_PCGSetMaxIter(solver, max_iter);
         JX_PCGSetTol(solver, tol);
         JX_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JX_PCGSetLogging(solver, 1);
         JX_PCGSetPrintLevel(solver, print_level);
         
         JX_PCGSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         JX_PCGSolve(solver, (JX_Matrix)mat_csr, // preOperater
                             (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         JX_PCGGetNumIterations(solver, &num_iterations);
         JX_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         JX_PCGDestroy(solver);
      }
      break;

      case 12:  /* PAMG-CG */
      {
         if (myid == 0) jx_printf("\n >>> Solver: PAMG-CG \n\n");
         
         starttime = jx_MPI_Wtime();
         
         JX_PAMGCreate(&amg_solver);
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
         JX_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetRelaxWt(amg_solver, relax_wt);
         JX_PAMGSetOuterWt(amg_solver, outer_wt);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);           /* sweep for "down" */
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);           /* sweep for "up" */
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);           /* sweep for "coarsest" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);           /* relax_type for "coarsest" */
         JX_PCGCreate(comm, &solver);
         
         JX_PCGSetMaxIter(solver, max_iter);
         JX_PCGSetTol(solver, tol);
         JX_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JX_PCGSetLogging(solver, 1);
         JX_PCGSetPrintLevel(solver, print_level);
         
         JX_PCGSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond, (JX_PtrToSolverFcn)JX_PAMGSetup, amg_solver);
         
         JX_PAMGSetup(amg_solver, mat_csr);
         
         JX_PCGSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-CG Setup", starttime, endtime, 0, 2);
         
         starttime = jx_MPI_Wtime();
         
         JX_PCGSolve(solver, (JX_Matrix)mat_csr, // preOperater
                             (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         endtime = jx_MPI_Wtime();
         jx_GetWallTime(comm, "PAMG-CG Solve", starttime, endtime, 0, 2);
         
         JX_PCGGetNumIterations(solver, &num_iterations);
         JX_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         if (print_level == 3 && myid == 0)
         {
            jx_printf(" >>> num_iterations = %d\n", num_iterations);
            jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }

         JX_PAMGDestroy(amg_solver);
         JX_PCGDestroy(solver);
      }
      break;
      
      case 21:  /* GMRES */
      {
         if (myid == 0) jx_printf("\n >>> Solver: GMRES(%d) \n\n", k_dim);
         
         JX_GMRESCreate(comm, &solver);
         JX_GMRESSetKDim(solver, k_dim);
         JX_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
         JX_GMRESSetMaxIter(solver, max_iter);
         JX_GMRESSetTol(solver, tol);
         JX_GMRESSetLogging(solver, 1);
         JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */
         
         JX_GMRESSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         JX_GMRESSolve(solver, (JX_Matrix)mat_csr, // preOperater
                               (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         JX_GMRESGetNumIterations(solver, &num_iterations);
         JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         JX_GMRESDestroy(solver);
      }
      break;

      case 22:  /* PAMG-GMRES */
      {
         if (myid == 0) jx_printf("\n >>> Solver: PAMG-GMRES(%d) \n\n", k_dim);
         
         if (TTest) starttime = jx_MPI_Wtime();
         
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
         JX_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
         JX_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JX_PAMGSetRelaxWt(amg_solver, relax_wt);
         JX_PAMGSetOuterWt(amg_solver, outer_wt);
         JX_PAMGSetSCommPkgSwitch(amg_solver, S_commpkg_switch);
         JX_PAMGSetAIRStrongTh(amg_solver, AIR_strong_th);
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);           /* sweep for "down" */
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);           /* sweep for "up" */
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);           /* sweep for "coarsest" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);           /* relax_type for "coarsest" */

         JX_GMRESCreate(comm, &solver);
         JX_GMRESSetKDim(solver, k_dim);
         JX_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* peghoty 2011/11/08 */
         JX_GMRESSetMaxIter(solver, max_iter);
         JX_GMRESSetTol(solver, tol);
         JX_GMRESSetLogging(solver, 1);
         JX_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */
         
         JX_GMRESSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond,
                                    (JX_PtrToSolverFcn)JX_PAMGSetup, amg_solver);
         
         JX_PAMGSetup(amg_solver, mat_csr);
         
         JX_GMRESSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG-GMRES Setup", starttime, endtime, 0, 2);
         }
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_GMRESSolve(solver, (JX_Matrix)mat_csr, // preOperater
                               (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG-GMRES Solve", starttime, endtime, 0, 2);
         }
         
         JX_GMRESGetNumIterations(solver, &num_iterations);
         JX_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         if (print_level == 3 && myid == 0)
         {
            jx_printf(" >>> num_iterations = %d\n", num_iterations);
            jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         
         JX_PAMGDestroy(amg_solver);
         JX_GMRESDestroy(solver);

        //jx_GetWallTime2(comm, "Setup Comm Time", jx__setup_commun_time, 0, 2);
        //jx_GetWallTime2(comm, "Solve Comm Time", jx__solve_commun_time, 0, 2);
      }
      break;
      
      case 31:  /* BiCGSTAB */
      {
         if (myid == 0) jx_printf("\n >>> Solver: BiCGSTAB \n\n");
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_BiCGSTABCreate(comm, &solver);
         JX_BiCGSTABSetMaxIter(solver, max_iter);
         JX_BiCGSTABSetTol(solver, tol);
         JX_BiCGSTABSetAbsoluteTol(solver, 0.0);
         JX_BiCGSTABSetConvCriteria(solver, 0);
         JX_BiCGSTABSetLogging(solver, 1);
         JX_BiCGSTABSetPrintLevel(solver, print_level);
         
         JX_BiCGSTABSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "BiCGSTAB Setup", starttime, endtime, 0, 2);
         }
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_BiCGSTABSolve(solver, (JX_Matrix)mat_csr, // preOperater
                                  (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "BiCGSTAB Solve", starttime, endtime, 0, 2);
         }
         
         JX_BiCGSTABGetNumIterations(solver, &num_iterations);
         JX_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         if (print_level == 0 && myid == 0)
         {
            jx_printf(" >>> num_iterations = %d\n", num_iterations);
            jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         
         JX_BiCGSTABDestroy(solver);
      }
      break;

      case 32:  /* PAMG-BiCGSTAB */
      {
         if (myid == 0) jx_printf("\n >>> Solver: PAMG-BiCGSTAB \n\n");
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_PAMGCreate(&amg_solver);
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
         if (ns_down > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_down, 1);           /* sweep for "down" */
         if (ns_up > -1) JX_PAMGSetCycleNumSweeps(amg_solver, ns_up, 2);           /* sweep for "up" */
         JX_PAMGSetCycleNumSweeps(amg_solver, ns_coarse, 3);           /* sweep for "coarsest" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);           /* relax_type for "coarsest" */

         JX_BiCGSTABCreate(comm, &solver);
         JX_BiCGSTABSetMaxIter(solver, max_iter);
         JX_BiCGSTABSetTol(solver, tol);
         JX_BiCGSTABSetAbsoluteTol(solver, 0.0);
         JX_BiCGSTABSetConvCriteria(solver, 0);
         JX_BiCGSTABSetLogging(solver, 1);
         JX_BiCGSTABSetPrintLevel(solver, print_level);
         
         JX_BiCGSTABSetPrecond(solver, (JX_PtrToSolverFcn)JX_PAMGPrecond,
                                       (JX_PtrToSolverFcn)JX_PAMGSetup, amg_solver);
         
         JX_PAMGSetup(amg_solver, mat_csr);
         
         JX_BiCGSTABSetup(solver, (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG-BiCGSTAB Setup", starttime, endtime, 0, 2);
         }
         
         if (TTest) starttime = jx_MPI_Wtime();
         
         JX_BiCGSTABSolve(solver, (JX_Matrix)mat_csr, // preOperater
                                  (JX_Matrix)mat_csr, (JX_Vector)vec_csr, (JX_Vector)sol_csr);
         
         if (TTest)
         {
            endtime = jx_MPI_Wtime();
            jx_GetWallTime(comm, "PAMG-BiCGSTAB Solve", starttime, endtime, 0, 2);
         }
         
         JX_BiCGSTABGetNumIterations(solver, &num_iterations);
         JX_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
         
         if (print_level == 0 && myid == 0)
         {
            jx_printf(" >>> num_iterations = %d\n", num_iterations);
            jx_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         
         JX_PAMGDestroy(amg_solver);
         JX_BiCGSTABDestroy(solver);
      }
      break;
   }
   
   if (TTest)
   {
      endtimeT = jx_MPI_Wtime();
      jx_GetWallTime(comm, "Total Sove Time", starttimeT, endtimeT, 0, 2);
   }
   
   //-----------------------------------------------------------
   //  销毁并行矩阵和并行向量
   //-----------------------------------------------------------
   #ifdef USING_HWLOC
   jx_hpCSRhardwareDestroy();
   #endif
   JX_IJMatrixDestroy(mat_ij);
   JX_IJVectorDestroy(vec_ij);
   JX_IJVectorDestroy(sol_ij);
   
   //-----------------------------------------------------------
   //  终止 MPI
   //-----------------------------------------------------------
   jx_MPI_Finalize();
   
   return 0;
}
