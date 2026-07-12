//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 * solver_laplace.c -- This is an example demonstrating how to solve a sparse
 * linear system arising from discretizations of Diffusion Operator
 * by calling JXFPAMG solver and the JXFPAMG-based Krylov iterative methods.
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
 * Created by Yue Xiaoqiang 2012/09/20
 *
 *  Xiangtan University
 *  yuexq1111@163.com
 *
 */

#include "jxf_pamg.h"
#include "jxf_diagscale.h"
#include "jxf_euclid.h"
#include "jxf_krylov.h"

JXF_Int
jxf_ParVectorSetRandomValues( jxf_ParVector *v, JXF_Int seed );
JXF_Int
jxf_SeqVectorSetRandomValues( jxf_Vector *x, JXF_Int seed );
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian( MPI_Comm comm,
                         JXF_Int nx,
                         JXF_Int ny,
                         JXF_Int nz,
                         JXF_Int P,
                         JXF_Int Q,
                         JXF_Int R,
                         JXF_Int p,
                         JXF_Int q,
                         JXF_Int r,
                         JXF_Real *value );
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian2d9pt( MPI_Comm comm,
                              JXF_Int nx,
                              JXF_Int ny,
                              JXF_Int P,
                              JXF_Int Q,
                              JXF_Int p,
                              JXF_Int q,
                              JXF_Real *value );
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian3d27pt( MPI_Comm comm,
                               JXF_Int nx,
                               JXF_Int ny,
                               JXF_Int nz,
                               JXF_Int P,
                               JXF_Int Q,
                               JXF_Int R,
                               JXF_Int p,
                               JXF_Int q,
                               JXF_Int r,
                               JXF_Real *value );
jxf_ParCSRMatrix *
jxf_GenerateParConvecDiff( MPI_Comm comm,
                          JXF_Int nx,
                          JXF_Int ny,
                          JXF_Int nz,
                          JXF_Int P,
                          JXF_Int Q,
                          JXF_Int R,
                          JXF_Int p,
                          JXF_Int q,
                          JXF_Int r,
                          JXF_Real  *value );
void
jxf_BuildParLaplacianMat( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr );
void
jxf_BuildParLaplacianMat2d9pt( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr );
void
jxf_BuildParLaplacianMat3d27pt( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr );
void
jxf_BuildParConvecDiff( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr );
void
jxf_BuildParLaplacianRhs( JXF_Int argc, char *argv[], jxf_hpCSRMatrix *A, jxf_ParVector **b_ptr );

int
main( int argc, char *argv[] )
{
   MPI_Comm     comm = MPI_COMM_WORLD;
   JXF_Int       myid, nprocs;
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
   JXF_Int       nthreads;
#endif
   
   JXF_Int arg_index   = 0;
   JXF_Int print_usage = 0;
   JXF_Int build_matrix_type = 0; /* Yue Xiaoqiang 2012/10/16 */
   
   JXF_Real starttime, endtime;
   JXF_Real starttimeT, endtimeT;
   
   // jxf_ParCSRMatrix *par_matrix = NULL;
   jxf_hpCSRMatrix *hp_matrix = NULL;
   jxf_ParVector   *par_rhs   = NULL;
   jxf_ParVector   *par_sol   = NULL;
   
   JXF_Int *partitioning = NULL;
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
   JXF_Int keepTranspose;
   JXF_Int coarsen_type;
   JXF_Int interp_type;
   JXF_Int P_max_elmts;
   JXF_Int agg_num_levels;
   JXF_Int ai_relax_type;
   JXF_Int ai_measure_type;
   JXF_Int coarsestsolverid;
   JXF_Int conv_criteria;
   JXF_Int coarse_threshold;
   JXF_Int amg_print_level;
   JXF_Real strong_threshold;
   JXF_Real max_row_sum;
   JXF_Real relax_wt;
   JXF_Real outer_wt;
   JXF_Real coarse_ratio;
   
   /* iterative method */
   JXF_Solver solver;
   JXF_Real tol;
   JXF_Int max_iter;
   JXF_Int k_dim; 
   JXF_Int is_check_restarted; /* Zhou Zhiyang 2011/11/08 */
   JXF_Int twonorm;
   JXF_Int solver_id;
   JXF_Int print_level;
   JXF_Int keepsol;
   JXF_Int new_mvcpu_flag;
   
   /* other variables */
   JXF_Int glosize;
   JXF_Int initguess = 0;
   JXF_Int num_iterations;
   JXF_Real final_res_norm;
   JXF_Real norm;
   JXF_Real mv_time_min = 0.0;
   JXF_Real mv_time_max = 0.0;
   JXF_Real mv_time_avg = 0.0;
   
   //----------------------------------------------------------------
   // 启动 MPI
   //----------------------------------------------------------------
   jxf_MPI_Init(&argc, &argv);
   jxf_MPI_Comm_rank(comm, &myid);
   jxf_MPI_Comm_size(comm, &nprocs);
   #ifdef USING_HWLOC
   jxf_hpCreateHardwareInfo(comm);
   #endif
   
   //----------------------------------------------------------------
   // 参数设置
   //----------------------------------------------------------------
   max_levels       = 25;       /* 最大网格层数 */
   cycle_type       = 1;        /* Cycle 类型  1: V_Cycle; 2：W_Cycle */
   relax_type       = 3;        /* Relax 类型  3: hGS; 6：hSGS */
   measure_type     = 0;        /* 影响值的计算方式 0：局部；1：全局 */
   rap2             = 0;        /* RAP计算方式  0：RAP；1：先算Q=AP，再算RQ */
   keepTranspose    = 0;        /* 存放限制算子  0：no；1：yes */
   coarsen_type     = 6;        /* 粗化策略 */
   interp_type      = 0;        /* 插值策略 */
   P_max_elmts      = 0;        /* 插值算子每行最大非零元个数 */
   agg_num_levels   = 0;        /* Aggressive粗化的层数 */
   ai_measure_type  = 0;        /* AI-策略  0: no; 1: yes, 注意调整AI-粗化和AI-磨光 */
   ai_relax_type    = 0;        /* AI-磨光  0: no; 1: yes */
   amg_print_level  = 0;        /* work only when AMG as preconditioner */
   strong_threshold = 0.25;     /* 强弱连通参数, 0.25 for 2D, 0.5 for 3D is recommended */
   max_row_sum      = 0.9;      /* 行和参数 */
   relax_wt         = 1.0;
   outer_wt         = 1.0;
   
   coarse_threshold = 9;        /* 最粗网格层上网格节点个数的最大值 */
   coarse_ratio     = 0.75;     /* 相邻两个网格层的粗点个数超过细点个数的 coarse_ratio, 则换成 CLJP 粗化 */
   coarsestsolverid = 9;        /* 最粗网格层解法器 */
   conv_criteria    = 0;        /* 收敛准则类型 */
   
   euclid_level     = 0;        /* level of fill-in */
   euclid_bj        = 0;        /* Select PILU (0) or Block Jacobi ILU (1) */
   
   tol                = 1.0e-7;  /* 控制精度 */
   max_iter           = 200;     /* 迭代法最大迭代次数 */
   k_dim              = 100;     /* 回头数 */
   is_check_restarted = 0;       /* Zhou Zhiyang 2011/11/08 */
   twonorm            = 0;       /* PCG 法中的范数控制类型，0: B 范数; 1: l2 范数 */
   print_level        = 3;       /* 0: 关闭；1：Setup参数；2：Solve参数；3：Setup+Solve参数 */
   keepsol            = 0;       /* 是否保存解向量 */
   new_mvcpu_flag     = 0;       /* 是否记录并行矩阵向量乘的CPU时间，0：不记录；1：记录 */
   solver_id          = 0;
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
   nthreads           = 1;       /* 线程数 Yue Xiaoqiang 2012/10/12 */
#endif
   
   //----------------------------------------------------------------
   // 命令行修改参数
   //----------------------------------------------------------------
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-sid") == 0)
      {
         arg_index ++;
         solver_id = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rap2") == 0 )
      {
         arg_index ++;
         rap2 = 1;
      }
      else if ( strcmp(argv[arg_index], "-kt") == 0 )
      {
         arg_index ++;
         keepTranspose = 1;
      }
      else if ( strcmp(argv[arg_index], "-mvcpu") == 0 )
      {
         arg_index ++;
         new_mvcpu_flag = 1;
      }
      else if ( strcmp(argv[arg_index], "-tnrm") == 0 )
      {
         arg_index ++;
         twonorm = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-kdim") == 0)
      {
         arg_index ++;
         k_dim = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-euc_lvl") == 0)
      {
         arg_index ++;
         euclid_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-euc_bj") == 0)
      {
         arg_index ++;
         euclid_bj = atoi(argv[arg_index++]);
      }
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
      else if (strcmp(argv[arg_index], "-nts") == 0)
      {
         arg_index ++;
         nthreads = atoi(argv[arg_index++]);
      }
#endif
      else if (strcmp(argv[arg_index], "-Pmx") == 0)
      {
         arg_index ++;
         P_max_elmts = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-agg_nl") == 0)
      {
         arg_index ++;
         agg_num_levels = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-9pt") == 0)
      {
         arg_index ++;
         build_matrix_type = 1;
      }
      else if (strcmp(argv[arg_index], "-27pt") == 0)
      {
         arg_index ++;
         build_matrix_type = 2;
      }
      else if (strcmp(argv[arg_index], "-ct") == 0)
      {
         arg_index ++;
         coarsen_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ipt") == 0)
      {
         arg_index ++;
         interp_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mxct") == 0)
      {
         arg_index ++;
         coarse_threshold = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-str") == 0)
      {
         arg_index ++;
         strong_threshold = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mxrs") == 0)
      {
         arg_index ++;
         max_row_sum = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-rlx") == 0)
      {
         arg_index ++;
         relax_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ai_rlx") == 0)
      {
         arg_index ++;
         ai_relax_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ai_mt") == 0)
      {
         arg_index ++;
         ai_measure_type = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-amg_ptlv") == 0)
      {
         arg_index ++;
         amg_print_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ptlv") == 0)
      {
         arg_index ++;
         print_level = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-help") == 0)
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
      jxf_printf("\n");
      jxf_printf("  Usage: %s [<options>]\n", argv[0]);
      jxf_printf("\n");
      jxf_printf("  -n <nx> <ny> <nz> : total problem size[10 10 10 being default]\n");
      jxf_printf("  -P <Px> <Py> <Pz> : processor topology[1 np 1 being default]\n");
      jxf_printf("  -c <cx> <cy> <cz> : diffusion coefficients[1.0 1.0 1.0 being default]\n");
      jxf_printf("     diffusion Opts : -cx Dxx - cy Dyy - cz Dzz\n");
      jxf_printf("  -one              : rhs is vector with unit components[default]\n");
      jxf_printf("  -rdm              : rhs is random vector and unit 2-norm\n");
      jxf_printf("  -sid  <val>       : solver id\n");
      jxf_printf("  -rap2             : 2nd implementation of RAP\n");
      jxf_printf("  -kt               : keep transpose\n");
      jxf_printf("  -mvcpu            : record CPU time of Matvec\n");
      jxf_printf("  -tnrm <val>       : two_norm in PCG\n");
      jxf_printf("  -kdim <val>       : krylov dimension\n");
      jxf_printf("  -euc_lvl <val>    : level of fill-in for Euclid\n");
      jxf_printf("  -euc_bj <val>     : PILU or Block Jacobi ILU\n");
      jxf_printf("  -nts <val>        : threads number\n");
      jxf_printf("  -Pmx <val>        : maximal number of elements per row for interpolation\n");
      jxf_printf("  -agg_nl <val>     : num_levels for aggressive coarsening\n");
      jxf_printf("  -9pt              : build 9pt 2D laplacian problem\n");
      jxf_printf("  -27pt             : build 27pt 3D laplacian problem\n");
      jxf_printf("  -ct <val>         : coarsening type\n");
      jxf_printf("  -ipt <val>        : interpolation type\n");
      jxf_printf("  -mxct <val>       : max. size on coarsest grid\n");
      jxf_printf("  -str <val>        : AMG strength threshold\n");
      jxf_printf("  -mxrs <val>       : maximum row sum threshold for dependency weakening\n");
      jxf_printf("  -rlx <val>        : relaxation type\n");
      jxf_printf("  -ai_rlx <val>     : AI relaxation type\n");
      jxf_printf("  -ai_mt  <val>     : AI measure type\n");
      jxf_printf("  -amg_ptlv <val>   : print_level of AMG when AMG as preconditioner\n");
      jxf_printf("  -ptlv <val>       : print_level\n");
      jxf_printf("  -help             : using help message\n\n");
      exit(1);
   }
   
   //----------------------------------------------------------------
   // 提供线性系统文件
   //----------------------------------------------------------------
   if (myid == 0)
   {
      jxf_printf("\n\n+++++++++++++++++++++");
#if JXF_USING_BIG_INT
      jxf_printf(" Using BIG_INT,");
#endif
#if JXF_USING_BIG_DOUBLE
      jxf_printf(" BIG_DOUBLE,");
#endif
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
      jxf_printf(" With OpenMP using %d threads,", nthreads);
#endif
      jxf_printf(" MPI using %d processors +++++++++++++++++++++\n\n", nprocs);
   }
   
   //----------------------------------------------------------------
   // 设定线程数
   //----------------------------------------------------------------
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
   omp_set_num_threads(nthreads);
#endif
   
   //----------------------------------------------------------------
   // 利用矩阵文件创建并行矩阵和并行右端以及并行初始迭代向量(零向量或随机向量)
   //----------------------------------------------------------------
   starttime = jxf_MPI_Wtime();
   if (build_matrix_type == 0)
   {
      jxf_BuildParLaplacianMat(argc, argv, &hp_matrix);
   }
   else if (build_matrix_type == 1)
   {
      jxf_BuildParLaplacianMat2d9pt(argc, argv, &hp_matrix);
   }
   else if (build_matrix_type == 2)
   {
      jxf_BuildParLaplacianMat3d27pt(argc, argv, &hp_matrix);
   }
   else if (build_matrix_type == 3)
   {
      jxf_BuildParConvecDiff(argc, argv, &hp_matrix);
   }
   jxf_BuildParLaplacianRhs(argc, argv, hp_matrix, &par_rhs);
   glosize      = jxf_ParVectorGlobalSize(par_rhs);
   partitioning = jxf_ParVectorPartitioning(par_rhs);
   par_sol      = jxf_ParVectorCreate(comm, glosize, partitioning);
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
      norm = 1.0 / sqrt(norm);
      jxf_ParVectorScale(norm, par_sol);
   }
   endtime = jxf_MPI_Wtime();
   jxf_GetWallTime(comm, "BuildParLinearSystem", starttime, endtime, 0, 2);
   
   jxf_set_mvcpu_handler(new_mvcpu_flag);
   
   //----------------------------------------------------------------
   // 求解线性代数系统
   //----------------------------------------------------------------
   starttimeT = jxf_MPI_Wtime();
   switch (solver_id)
   {
      case 0:  /* PAMG */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: PAMG \n\n");
         }
         starttime = jxf_MPI_Wtime();
         JXF_PAMGCreate(&amg_solver);
         JXF_PAMGSetMaxLevels(amg_solver, max_levels);
         JXF_PAMGSetMaxIter(amg_solver, max_iter);
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
         JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
         JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetPrintLevel(amg_solver, print_level);
         JXF_PAMGSetCoarsestSolverID(amg_solver, coarsestsolverid);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetCoarseRatio(amg_solver, coarse_ratio);
         JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
         JXF_PAMGSetOuterWt(amg_solver, outer_wt);
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);  /* sweep for "down" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);  /* sweep for "up" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);  /* sweep for "coarsest" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);  /* relax_type for "coarsest" */
         
         //----------------------------------------------------------------
         // JXF_PAMG Setup
         //----------------------------------------------------------------
         if (max_levels != 1)
         {
            JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) hp_matrix);
         }
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG Setup", starttime, endtime, 0, 2);
         
         //----------------------------------------------------------------
         // JXF_PAMG Solve
         //----------------------------------------------------------------
         starttime = jxf_MPI_Wtime();
         JXF_PAMGSolve(amg_solver, (JXF_hpCSRMatrix)hp_matrix, (JXF_ParVector)par_rhs, (JXF_ParVector)par_sol);
         JXF_PAMGGetNumIterations(amg_solver, &num_iterations);
         JXF_PAMGGetFinalRelativeResidualNorm(amg_solver, &final_res_norm);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG Solve", starttime, endtime, 0, 2);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PAMGDestroy(amg_solver);
      }
      break;
      
      case 11:  /* CG */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: CG \n\n");
         }
         JXF_PCGCreate(comm, &solver);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetTol(solver, tol);
         JXF_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JXF_PCGSetLogging(solver, 1);
         JXF_PCGSetPrintLevel(solver, print_level);
         
         //----------------------------------------------------------------
         // JXF_PCG Setup
         //----------------------------------------------------------------
         JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         
         //----------------------------------------------------------------
         // JXF_PCG Solve
         //----------------------------------------------------------------
         JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                             (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PCGDestroy(solver);
      }
      break;
      
      case 12:  /* PAMG-CG */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: PAMG-CG \n\n");
         }
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
         JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
         JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
         JXF_PAMGSetOuterWt(amg_solver, outer_wt);
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);  /* sweep for "down" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);  /* sweep for "up" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);  /* sweep for "coarsest" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);  /* relax_type for "coarsest" */
         JXF_PCGCreate(comm, &solver);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetTol(solver, tol);
         JXF_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JXF_PCGSetLogging(solver, 1);
         JXF_PCGSetPrintLevel(solver, print_level);
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_PAMGPrecond,
                                  (JXF_PtrToSolverFcn)JXF_PAMGSetup, amg_solver);
         
         //----------------------------------------------------------------
         // JXF_PAMG and JXF_PCG Setup
         //----------------------------------------------------------------
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
         JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-CG Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         
         //----------------------------------------------------------------
         // JXF_PCG Solve
         //----------------------------------------------------------------
         JXF_PCGSolve(solver, (JXF_Matrix) hp_matrix, // preOperater
                             (JXF_Matrix) hp_matrix, (JXF_Vector) par_rhs, (JXF_Vector) par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-CG Solve", starttime, endtime, 0, 2);
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PAMGDestroy(amg_solver);
         JXF_PCGDestroy(solver);
      }
      break;
      
      case 13:  /* DS-CG */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: DS-CG \n\n");
         }
         starttime = jxf_MPI_Wtime();
         ds_solver = NULL;
         JXF_PCGCreate(comm, &solver);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetTol(solver, tol);
         JXF_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JXF_PCGSetLogging(solver, 1);
         JXF_PCGSetPrintLevel(solver, print_level);
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_DiagScalePrecond,
                                  (JXF_PtrToSolverFcn)JXF_DiagScaleSetup, ds_solver);
         
         //----------------------------------------------------------------
         // JXF_DS and JXF_PCG Setup
         //----------------------------------------------------------------
         JXF_DiagScaleSetup(ds_solver, (JXF_hpCSRMatrix)hp_matrix);
         JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-CG Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         
         //----------------------------------------------------------------
         // JXF_PCG Solve
         //----------------------------------------------------------------
         JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                             (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-CG Solve", starttime, endtime, 0, 2);
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PCGDestroy(solver);
      }
      break;
      
      case 14:  /* Euclid-CG */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: Euclid-CG \n\n");
         }
         starttime = jxf_MPI_Wtime();
         JXF_EuclidCreate(comm, &euclid_solver);
         JXF_EuclidSetParams(euclid_solver, argc, argv);
         JXF_EuclidSetLevel(euclid_solver, euclid_level);
         JXF_EuclidSetBJ(euclid_solver, euclid_bj);
         JXF_PCGCreate(comm, &solver);
         JXF_PCGSetMaxIter(solver, max_iter);
         JXF_PCGSetTol(solver, tol);
         JXF_PCGSetTwoNorm(solver, twonorm);  // 0: B 范数； 1：l2 范数
         JXF_PCGSetLogging(solver, 1);
         JXF_PCGSetPrintLevel(solver, print_level);
         JXF_PCGSetPrecond(solver, (JXF_PtrToSolverFcn)JXF_EuclidSolve,
                                  (JXF_PtrToSolverFcn)JXF_EuclidSetup, euclid_solver);
         
         //----------------------------------------------------------------
         // JXF_Euclid and JXF_PCG Setup
         //----------------------------------------------------------------
         JXF_EuclidSetup(euclid_solver, (JXF_hpCSRMatrix)hp_matrix);
         JXF_PCGSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-CG Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         
         //----------------------------------------------------------------
         // JXF_PCG Solve
         //----------------------------------------------------------------
         JXF_PCGSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                             (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-CG Solve", starttime, endtime, 0, 2);
         JXF_PCGGetNumIterations(solver, &num_iterations);
         JXF_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_EuclidDestroy(euclid_solver);
         JXF_PCGDestroy(solver);
      }
      break;
      
      case 21:  /* GMRES */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: GMRES(%d) \n\n", k_dim);
         }
         JXF_GMRESCreate(comm, &solver);
         JXF_GMRESSetKDim(solver, k_dim);
         JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* Zhou Zhiyang 2011/11/08 */
         JXF_GMRESSetMaxIter(solver, max_iter);
         JXF_GMRESSetTol(solver, tol);
         JXF_GMRESSetLogging(solver, 1);
         JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */
         
         //----------------------------------------------------------------
         // JXF_GMRES Setup
         //----------------------------------------------------------------
         JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         
         //----------------------------------------------------------------
         // JXF_GMRES Solve
         //----------------------------------------------------------------
         JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                               (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         JXF_GMRESGetNumIterations(solver, &num_iterations);
         JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_GMRESDestroy(solver);
      }
      break;
      
      case 22:  /* PAMG-GMRES */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: PAMG-GMRES(%d) \n\n", k_dim);
         }
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
         JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
         JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
         JXF_PAMGSetOuterWt(amg_solver, outer_wt);
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);  /* sweep for "down" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);  /* sweep for "up" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);  /* sweep for "coarsest" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);  /* relax_type for "coarsest" */
         JXF_GMRESCreate(comm, &solver);
         JXF_GMRESSetKDim(solver, k_dim);
         JXF_GMRESSetIsCheckRestarted(solver, is_check_restarted); /* Zhou Zhiyang 2011/11/08 */
         JXF_GMRESSetMaxIter(solver, max_iter);
         JXF_GMRESSetTol(solver, tol);
         JXF_GMRESSetLogging(solver, 1);
         JXF_GMRESSetPrintLevel(solver, print_level); /* 是否在屏幕上打印残量等信息 */
         JXF_GMRESSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond,
                                    (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         
         //----------------------------------------------------------------
         // JXF_PAMG and JXF_GMRES Solve
         //----------------------------------------------------------------
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix)hp_matrix);
         JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-GMRES Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         
         //----------------------------------------------------------------
         // JXF_GMRES Solve
         //----------------------------------------------------------------
         JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                               (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "PAMG-GMRES Solve", starttime, endtime, 0, 2);
         JXF_GMRESGetNumIterations(solver, &num_iterations);
         JXF_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PAMGDestroy(amg_solver);
         JXF_GMRESDestroy(solver);
      }
      break;
      
      case 23:  /* DS-GMRES */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: DS-GMRES(%d) \n\n", k_dim);
         }
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
         
         //----------------------------------------------------------------
         // JXF_DS and JXF_GMRES Setup
         //----------------------------------------------------------------
         JXF_DiagScaleSetup(ds_solver, (JXF_hpCSRMatrix)hp_matrix);
         JXF_GMRESSetup(solver, (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-GMRES Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         
         //----------------------------------------------------------------
         // JXF_GMRES Solve
         //----------------------------------------------------------------
         JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                               (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-GMRES Solve", starttime, endtime, 0, 2);
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
      
      case 24:  /* Euclid-GMRES */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: Euclid-GMRES(%d) \n\n", k_dim);
         }
         starttime = jxf_MPI_Wtime();
         JXF_EuclidCreate(comm, &euclid_solver);
         JXF_EuclidSetParams(euclid_solver, argc, argv);
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
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-GMRES Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         JXF_GMRESSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                               (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-GMRES Solve", starttime, endtime, 0, 2);
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
      
      case 31:  /* BICGStab */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: BICGStab \n\n");
         }
         JXF_BiCGSTABCreate(comm, &solver);
         JXF_BiCGSTABSetMaxIter(solver, max_iter);
         JXF_BiCGSTABSetTol(solver, tol);
         JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
         JXF_BiCGSTABSetConvCriteria(solver, 0);
         JXF_BiCGSTABSetLogging(solver, 1);
         JXF_BiCGSTABSetPrintLevel(solver, print_level);
         JXF_BiCGSTABSetup(solver, (JXF_Matrix) hp_matrix, (JXF_Vector) par_rhs, (JXF_Vector) par_sol);
         JXF_BiCGSTABSolve(solver, (JXF_Matrix) hp_matrix, // preOperater
                                  (JXF_Matrix) hp_matrix, (JXF_Vector) par_rhs, (JXF_Vector) par_sol);
         JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
         JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_BiCGSTABDestroy(solver);
      }
      break;
      
      case 32:  /* PAMG-BICGStab */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: PAMG-BICGStab \n\n");
         }
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
         JXF_PAMGSetAIRelaxType(amg_solver, ai_relax_type);
         JXF_PAMGSetAIMeasureType(amg_solver, ai_measure_type);
         JXF_PAMGSetStrongThreshold(amg_solver, strong_threshold);
         JXF_PAMGSetMaxRowSum(amg_solver, max_row_sum);
         JXF_PAMGSetPrintLevel(amg_solver, amg_print_level);
         JXF_PAMGSetCoarseThreshold(amg_solver, coarse_threshold);
         JXF_PAMGSetRelaxWt(amg_solver, relax_wt);
         JXF_PAMGSetOuterWt(amg_solver, outer_wt);
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 1);  /* sweep for "down" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 2);  /* sweep for "up" */
         JXF_PAMGSetCycleNumSweeps(amg_solver, 1, 3);  /* sweep for "coarsest" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);  /* relax_type for "down" */
         JXF_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);  /* relax_type for "up" */
         JXF_PAMGSetCycleRelaxType(amg_solver, 9, 3);  /* relax_type for "coarsest" */
         JXF_BiCGSTABCreate(comm, &solver);
         JXF_BiCGSTABSetMaxIter(solver, max_iter);
         JXF_BiCGSTABSetTol(solver, tol);
         JXF_BiCGSTABSetAbsoluteTol(solver, 0.0);
         JXF_BiCGSTABSetConvCriteria(solver, 0);
         JXF_BiCGSTABSetLogging(solver, 1);
         JXF_BiCGSTABSetPrintLevel(solver, print_level);
         JXF_BiCGSTABSetPrecond(solver, (JXF_PtrToSolverFcn) JXF_PAMGPrecond,
                                       (JXF_PtrToSolverFcn) JXF_PAMGSetup, amg_solver);
         JXF_PAMGSetup(amg_solver, (JXF_hpCSRMatrix) hp_matrix);
         JXF_BiCGSTABSetup(solver, (JXF_Matrix) hp_matrix, (JXF_Vector) par_rhs, (JXF_Vector) par_sol);
         JXF_BiCGSTABSolve(solver, (JXF_Matrix) hp_matrix, // preOperater
                                  (JXF_Matrix) hp_matrix, (JXF_Vector) par_rhs, (JXF_Vector) par_sol);
         JXF_BiCGSTABGetNumIterations(solver, &num_iterations);
         JXF_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
         if ((print_level == 0) && (myid == 0))
         {
            jxf_printf(" >>> num_iterations = %d\n", num_iterations);
            jxf_printf(" >>> final_res_norm = %.4le\n", final_res_norm);
         }
         JXF_PAMGDestroy(amg_solver);
         JXF_BiCGSTABDestroy(solver);
      }
      break;
      
      case 33:  /* DS-BiCGSTAB */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: DS-BiCGSTAB \n\n");
         }
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
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-BiCGSTAB Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "DS-BiCGSTAB Solve", starttime, endtime, 0, 2);
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
      
      case 34:  /* Euclid-BiCGSTAB */
      {
         if (myid == 0)
         {
            jxf_printf("\n >>> Solver: Euclid-BiCGSTAB \n\n");
         }
         starttime = jxf_MPI_Wtime();
         JXF_EuclidCreate(comm, &euclid_solver);
         JXF_EuclidSetParams(euclid_solver, argc, argv);
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
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-BiCGSTAB Setup", starttime, endtime, 0, 2);
         starttime = jxf_MPI_Wtime();
         JXF_BiCGSTABSolve(solver, (JXF_Matrix)hp_matrix, // preOperater
                                  (JXF_Matrix)hp_matrix, (JXF_Vector)par_rhs, (JXF_Vector)par_sol);
         endtime = jxf_MPI_Wtime();
         jxf_GetWallTime(comm, "Euclid-BiCGSTAB Solve", starttime, endtime, 0, 2);
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
   endtimeT = jxf_MPI_Wtime();
   jxf_GetWallTime(comm, "Total Sove Time", starttimeT, endtimeT, 0, 2);

   if (jxf__global_mvcpu_flag)
   {
      jxf_MPI_Reduce(&jxf_total_elapsed_time_matvec, &mv_time_min, 1, JXF_MPI_REAL, MPI_MIN, 0, comm);
      jxf_MPI_Reduce(&jxf_total_elapsed_time_matvec, &mv_time_max, 1, JXF_MPI_REAL, MPI_MAX, 0, comm);
      jxf_MPI_Reduce(&jxf_total_elapsed_time_matvec, &mv_time_avg, 1, JXF_MPI_REAL, MPI_SUM, 0, comm);
      if (myid == 0) jxf_printf("\n >> %s: time(min,max,ave) = (%.1f, %.1f, %.1f) seconds\n\n",
                         "Matvec Time", mv_time_min, mv_time_max, mv_time_avg / nprocs);
   }

   //----------------------------------------------------------------
   // 将近似解向量保存到文件中
   //----------------------------------------------------------------
   if (keepsol)
   {
      jxf_Vector *ser_sol = NULL;
      ser_sol = jxf_ParVectorToVectorAll(par_sol);
      if (myid == 0)
      {
         jxf_SeqVectorPrint(ser_sol, "./app");
      }
      jxf_SeqVectorDestroy(ser_sol);
   }
   
   //----------------------------------------------------------------
   // 销毁并行矩阵和并行向量
   //----------------------------------------------------------------
   #ifdef USING_HWLOC
   jxf_hpCSRhardwareDestroy();
   #endif
   jxf_hpCSRMatrixDestroy(hp_matrix);
   jxf_ParVectorDestroy(par_rhs);
   jxf_ParVectorDestroy(par_sol);
   
   //----------------------------------------------------------------
   // 终止 MPI
   //----------------------------------------------------------------
   jxf_MPI_Finalize();
   
   return 0;
}

/*!
 * \fn JXF_Int jxf_ParVectorSetRandomValues
 */
JXF_Int
jxf_ParVectorSetRandomValues( jxf_ParVector *v, JXF_Int seed )
{
   JXF_Int my_id;
   jxf_Vector *v_local = jxf_ParVectorLocalVector(v);

   MPI_Comm 	comm = jxf_ParVectorComm(v);
   jxf_MPI_Comm_rank(comm, &my_id); 

   seed *= (my_id + 1);
           
   return jxf_SeqVectorSetRandomValues(v_local,seed);
}

/*!
 * \fn JXF_Int jxf_SeqVectorSetRandomValues
 */
JXF_Int
jxf_SeqVectorSetRandomValues( jxf_Vector *x, JXF_Int seed )
{
   JXF_Real  *vector_data = jxf_VectorData(x);
   JXF_Int      size        = jxf_VectorSize(x);
   JXF_Int      i, ierr = 0;

   jxf_SeedRand(seed);

   size *= jxf_VectorNumVectors(x);

   /* RDF: threading this loop may cause problems because of jxf_Rand() */
   for (i = 0; i < size; i ++)
   {
      vector_data[i] = 2.0 * jxf_Rand() - 1.0;
   }
   
   return ierr;
}

/*!
 * \fn void jxf_BuildParLaplacianMat
 */
void
jxf_BuildParLaplacianMat( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr )
{
    jxf_hpCSRMatrix *hp_A = NULL;

    JXF_Real *values = NULL;
    JXF_Int nx, ny, nz;
    JXF_Int P, Q, R;
    JXF_Int p, q, r;
    JXF_Int num_procs, my_id;
    JXF_Int arg_index;
    JXF_Real cx, cy, cz;
    
   /*-----------------------
    * Initialize some stuff
    *-----------------------*/
    jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    
   /*--------------
    * Set defaults
    *--------------*/
    nx = 10;
    ny = 10;
    nz = 10;
    P = 1;
    Q = num_procs;
    R = 1;
    cx = 1.0;
    cy = 1.0;
    cz = 1.0;
    
   /*--------------------
    * Parse command line
    *--------------------*/
    arg_index = 0;
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-n") == 0)
        {
            arg_index ++;
            nx = atoi(argv[arg_index++]);
            ny = atoi(argv[arg_index++]);
            nz = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-P") == 0)
        {
            arg_index ++;
            P = atoi(argv[arg_index++]);
            Q = atoi(argv[arg_index++]);
            R = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-c") == 0)
        {
            arg_index ++;
            cx = atof(argv[arg_index++]);
            cy = atof(argv[arg_index++]);
            cz = atof(argv[arg_index++]);
        }
        else
        {
            arg_index ++;
        }
    }
    
   /*--------------------
    * Check a few things
    *--------------------*/
    if ((P*Q*R) != num_procs)
    {
        jxf_printf("Error: Invalid number of processors or processor topology\n");
        exit(1);
    }
    
   /*-------------------------
    * Print driver parameters
    *-------------------------*/
    if (my_id == 0)
    {
        jxf_printf("  Laplacian:\n");
        jxf_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
        jxf_printf("    (Px, Py, Pz) = (%d, %d, %d)\n", P, Q, R);
        jxf_printf("    (cx, cy, cz) = (%f, %f, %f)\n\n", cx, cy, cz);
    }
    
   /*---------------------------
    * Set up the grid structure
    *---------------------------*/
    p = my_id % P;
    q = ((my_id - p) / P) % Q;
    r = (my_id - p - P * q) / (P * Q);
    
   /*---------------------
    * Generate the matrix
    *---------------------*/
    values = jxf_CTAlloc(JXF_Real, 4);
    values[1] = -cx;
    values[2] = -cy;
    values[3] = -cz;
    values[0] = 0.0;
    if (nx > 1)
    {
        values[0] += 2.0 * cx;
    }
    if (ny > 1)
    {
        values[0] += 2.0 * cy;
    }
    if (nz > 1)
    {
        values[0] += 2.0 * cz;
    }

    hp_A = jxf_hpInithpCSRMatrix();
    jxf_hpCSRMatrixPar(hp_A) = jxf_GenerateParLaplacian(MPI_COMM_WORLD, nx, ny, nz, P, Q, R, p, q, r, values);
    jxf_TFree(values);
   *A_ptr = hp_A;
}

/*!
 * \fn void jxf_BuildParLaplacianMat2d9pt
 */
void
jxf_BuildParLaplacianMat2d9pt( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr )
{
    jxf_hpCSRMatrix *hp_A = NULL;
    JXF_Real *values = NULL;
    JXF_Int num_procs, my_id;
    JXF_Int nx, ny;
    JXF_Int p, q;
    JXF_Int P, Q;
    JXF_Int arg_index;
    
   /*-----------------------
    * Initialize some stuff
    *-----------------------*/
    jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    
   /*--------------
    * Set defaults
    *--------------*/
    nx = 10;
    ny = 10;
    P = 1;
    Q = num_procs;
    
   /*--------------------
    * Parse command line
    *--------------------*/
    arg_index = 0;
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-n") == 0)
        {
            arg_index ++;
            nx = atoi(argv[arg_index++]);
            ny = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-P") == 0)
        {
            arg_index ++;
            P = atoi(argv[arg_index++]);
            Q = atoi(argv[arg_index++]);
        }
        else
        {
            arg_index ++;
        }
    }
    
   /*--------------------
    * Check a few things
    *--------------------*/
    if ((P*Q) != num_procs)
    {
        jxf_printf("Error: Invalid number of processors or processor topology\n");
        exit(1);
    }
    
   /*-------------------------
    * Print driver parameters
    *-------------------------*/
    if (my_id == 0)
    {
        jxf_printf("  Laplacian 9pt:\n");
        jxf_printf("    (nx, ny) = (%d, %d)\n", nx, ny);
        jxf_printf("    (Px, Py) = (%d, %d)\n\n", P, Q);
    }
    
   /*---------------------------
    * Set up the grid structure
    *---------------------------*/
    p = my_id % P;
    q = (my_id - p) / P;
    
   /*---------------------
    * Generate the matrix
    *---------------------*/
    values = jxf_CTAlloc(JXF_Real, 2);
    values[1] = -1.0;
    values[0] = 0.0;
    if (nx > 1)
    {
        values[0] += 2.0;
    }
    if (ny > 1)
    {
        values[0] += 2.0;
    }
    if ((nx > 1) && (ny > 1))
    {
        values[0] += 4.0;
    }
    hp_A = jxf_hpInithpCSRMatrix();
    jxf_hpCSRMatrixPar(hp_A) = jxf_GenerateParLaplacian2d9pt(MPI_COMM_WORLD, nx, ny, P, Q, p, q, values);
    jxf_TFree(values);
   *A_ptr = hp_A;
}

/*!
 * \fn void jxf_BuildParLaplacianMat3d27pt
 */
void
jxf_BuildParLaplacianMat3d27pt( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr )
{
    jxf_hpCSRMatrix *hp_A = NULL;
    JXF_Real *values = NULL;
    JXF_Int nx, ny, nz;
    JXF_Int P, Q, R;
    JXF_Int p, q, r;
    JXF_Int num_procs, my_id;
    JXF_Int arg_index;
    
   /*-----------------------
    * Initialize some stuff
    *-----------------------*/
    jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    
   /*--------------
    * Set defaults
    *--------------*/
    nx = 10;
    ny = 10;
    nz = 10;
    P = 1;
    Q = num_procs;
    R = 1;
    
   /*--------------------
    * Parse command line
    *--------------------*/
    arg_index = 0;
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-n") == 0)
        {
            arg_index ++;
            nx = atoi(argv[arg_index++]);
            ny = atoi(argv[arg_index++]);
            nz = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-P") == 0)
        {
            arg_index ++;
            P = atoi(argv[arg_index++]);
            Q = atoi(argv[arg_index++]);
            R = atoi(argv[arg_index++]);
        }
        else
        {
            arg_index ++;
        }
    }
    
   /*--------------------
    * Check a few things
    *--------------------*/
    if ((P*Q*R) != num_procs)
    {
        jxf_printf("Error: Invalid number of processors or processor topology\n");
        exit(1);
    }
    
   /*-------------------------
    * Print driver parameters
    *-------------------------*/
    if (my_id == 0)
    {
        jxf_printf("  Laplacian_27pt:\n");
        jxf_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
        jxf_printf("    (Px, Py, Pz) = (%d, %d, %d)\n\n", P, Q, R);
    }
    
   /*---------------------------
    * Set up the grid structure
    *---------------------------*/
    p = my_id % P;
    q = ((my_id - p) / P) % Q;
    r = (my_id - p - P * q) / (P * Q);
    
   /*---------------------
    * Generate the matrix
    *---------------------*/
    values = jxf_CTAlloc(JXF_Real, 2);
    values[0] = 26.0;
    if ((nx == 1) || (ny == 1) || (nz == 1))
    {
        values[0] = 8.0;
    }
    if ((nx*ny == 1) || (nx*nz == 1) || (ny*nz == 1))
    {
        values[0] = 2.0;
    }
    values[1] = -1.0;
    hp_A = jxf_hpInithpCSRMatrix();
    jxf_hpCSRMatrixPar(hp_A) = jxf_GenerateParLaplacian3d27pt(MPI_COMM_WORLD, nx, ny, nz, P, Q, R, p, q, r, values);
    jxf_TFree(values);
   *A_ptr = hp_A;
}

/*!
 * \fn void jxf_BuildParConvecDiff
 */
void
jxf_BuildParConvecDiff( JXF_Int argc, char *argv[], jxf_hpCSRMatrix **A_ptr )
{
    jxf_hpCSRMatrix *hp_A = NULL;
    JXF_Real *values = NULL;
    JXF_Int nx, ny, nz;
    JXF_Int P, Q, R;
    JXF_Int p, q, r;
    JXF_Int num_procs, my_id;
    JXF_Int arg_index;
    JXF_Real cx, cy, cz;
    JXF_Real ax, ay, az;
    JXF_Real hinx, hiny, hinz;
    
   /*-----------------------
    * Initialize some stuff
    *-----------------------*/
    jxf_MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    
   /*--------------
    * Set defaults
    *--------------*/
    nx = 10;
    ny = 10;
    nz = 10;
    P = 1;
    Q = num_procs;
    R = 1;
    cx = 1.0;
    cy = 1.0;
    cz = 1.0;
    ax = 1.0;
    ay = 1.0;
    az = 1.0;
   /*--------------------
    * Parse command line
    *--------------------*/
    arg_index = 0;
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-n") == 0)
        {
            arg_index ++;
            nx = atoi(argv[arg_index++]);
            ny = atoi(argv[arg_index++]);
            nz = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-P") == 0)
        {
            arg_index ++;
            P = atoi(argv[arg_index++]);
            Q = atoi(argv[arg_index++]);
            R = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-c") == 0)
        {
            arg_index ++;
            cx = atof(argv[arg_index++]);
            cy = atof(argv[arg_index++]);
            cz = atof(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-a") == 0)
        {
            arg_index ++;
            ax = atof(argv[arg_index++]);
            ay = atof(argv[arg_index++]);
            az = atof(argv[arg_index++]);
        }
        else
        {
            arg_index ++;
        }
    }
    
   /*--------------------
    * Check a few things
    *--------------------*/
    if ((P*Q*R) != num_procs)
    {
        jxf_printf("Error: Invalid number of processors or processor topology\n");
        exit(1);
    }
    
   /*-------------------------
    * Print driver parameters
    *-------------------------*/
    if (my_id == 0)
    {
        jxf_printf("  Convection-Diffusion: \n");
        jxf_printf("    -cx Dxx - cy Dyy - cz Dzz + ax Dx + ay Dy + az Dz = f\n");
        jxf_printf("    (nx, ny, nz) = (%d, %d, %d)\n", nx, ny, nz);
        jxf_printf("    (Px, Py, Pz) = (%d, %d, %d)\n", P,  Q,  R);
        jxf_printf("    (cx, cy, cz) = (%f, %f, %f)\n", cx, cy, cz);
        jxf_printf("    (ax, ay, az) = (%f, %f, %f)\n\n", ax, ay, az);
    }
    
   /*---------------------------
    * Set up the grid structure
    *---------------------------*/
    p = my_id % P;
    q = ((my_id - p) / P) % Q;
    r = (my_id - p - P * q) / (P * Q);
    hinx = 1.0 / (nx + 1);
    hiny = 1.0 / (ny + 1);
    hinz = 1.0 / (nz + 1);
    
   /*---------------------
    * Generate the matrix
    *---------------------*/
    values = jxf_CTAlloc(JXF_Real, 7);
    values[1] = -cx / (hinx * hinx);
    values[2] = -cy / (hiny * hiny);
    values[3] = -cz / (hinz * hinz);
    values[4] = -cx / (hinx * hinx) + ax / hinx;
    values[5] = -cy / (hiny * hiny) + ay / hiny;
    values[6] = -cz / (hinz * hinz) + az / hinz;
    values[0] = 0.0;
    if (nx > 1)
    {
        values[0] += 2.0 * cx / (hinx * hinx) - 1.0 * ax / hinx;
    }
    if (ny > 1)
    {
        values[0] += 2.0 * cy / (hiny * hiny) - 1.0 * ay / hiny;
    }
    if (nz > 1)
    {
        values[0] += 2.0 * cz / (hinz * hinz) - 1.0 * az / hinz;
    }
    hp_A = jxf_hpInithpCSRMatrix();
    jxf_hpCSRMatrixPar(hp_A) = jxf_GenerateParConvecDiff(MPI_COMM_WORLD, nx, ny, nz, P, Q, R, p, q, r, values);
    jxf_TFree(values);
   *A_ptr = hp_A;
}

/*!
 * \fn void jxf_BuildParLaplacianRhs
 */
void
jxf_BuildParLaplacianRhs( JXF_Int argc, char *argv[], jxf_hpCSRMatrix *hp_A, jxf_ParVector **b_ptr )
{
    JXF_Int build_rhs_type = 1;
    jxf_ParVector *b = NULL;
    jxf_Vector *loc_b = NULL;
    JXF_Int *partition = NULL;
    JXF_Int my_id, global_num_rows;
    JXF_Int arg_index;
    JXF_Real norm;
    jxf_ParCSRMatrix *A = jxf_hpCSRMatrixPar(hp_A);
    
   /*-----------------------
    * Initialize some stuff
    *-----------------------*/
    jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    
   /*-------------------
    * Get the partition
    *-------------------*/
    jxf_ParCSRMatrixGetRowPartitioning(A, &partition);
    global_num_rows = jxf_ParCSRMatrixGlobalNumRows(A);
    
   /*--------------------
    * Parse command line
    *--------------------*/
    arg_index = 0;
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-one") == 0)
        {
            arg_index ++;
            build_rhs_type = 1;
        }
        else if (strcmp(argv[arg_index], "-rdm") == 0)
        {
            arg_index ++;
            build_rhs_type = 2;
        }
        else
        {
            arg_index ++;
        }
    }
    if (build_rhs_type == 1)
    {
        if (my_id == 0)
        {
            jxf_printf("  RHS vector has unit components\n");
            loc_b = jxf_SeqVectorCreate(global_num_rows);
            jxf_SeqVectorInitialize(loc_b);
            jxf_SeqVectorSetConstantValues(loc_b, 1.0);
        }
    }
    else if (build_rhs_type == 2)
    {
        if (my_id == 0)
        {
            jxf_printf("  RHS vector has random components and unit 2-norm\n");
            loc_b = jxf_SeqVectorCreate(global_num_rows);
            jxf_SeqVectorInitialize(loc_b);
            jxf_SeqVectorSetRandomValues(loc_b, 22775);
            norm = jxf_SeqVectorInnerProd(loc_b, loc_b);
            norm = 1.0 / sqrt(norm);
            jxf_SeqVectorScale(norm, loc_b);
        }
    }
    b = jxf_VectorToParVector(MPI_COMM_WORLD, loc_b, partition);
    if (my_id == 0)
    {
        jxf_SeqVectorDestroy(loc_b);
    }
   *b_ptr = b;
}
