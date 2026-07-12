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

#include "hthread_host.h"

JX_Int myid;
JX_Int coreNums = 20;
JX_Int dot_type = 0;
JX_Int jx_spmv_type = 0;
JXF_Int jxf_spmv_type = 0;
JXF_Int jxh_spmv_type = 0;
JXF_Int mix_spmv_type = 0;

JX_Int
jx_ParVectorSetRandomValues(jx_ParVector *v, JX_Int seed);
JX_Int
jx_SeqVectorSetRandomValues(jx_Vector *x, JX_Int seed);

/* 高精度矩阵转低精度 */
JXF_Int
jxmp_CSRMatrixDtoF(jx_CSRMatrix *A, jxf_CSRMatrix *B, JX_Int copy_data);
JXF_Int
jxmp_ParCSRMatrixDtoF(jx_ParCSRMatrix *A, jxf_ParCSRMatrix *B, JX_Int copy_data);

/* 高精度向量转低精度 */
JXF_Int
jxmp_ParVectorDtoF(jx_ParVector *x, jxf_ParVector *y);
JXF_Int
jxmp_SeqVectorDtoF(jx_Vector *x, jxf_Vector *y);

static JX_Real
diff_norm2_f32_vs_fp64(jx_ParVector *ref_D,
                       jxf_ParVector *cand_F,
                       MPI_Comm comm)
{
   jx_Vector *ref_local = jx_ParVectorLocalVector(ref_D);
   jxf_Vector *cand_local = jxf_ParVectorLocalVector(cand_F);

   JX_Real *ref_data = jx_VectorData(ref_local);
   JXF_Real *cand_data = jxf_VectorData(cand_local);

   JX_Int n = jx_VectorSize(ref_local);

   JX_Real local_sum = 0.0;
   JX_Real global_sum = 0.0;

   for (JX_Int i = 0; i < n; i++)
   {
      JX_Real d = ref_data[i] - (JX_Real)cand_data[i];
      local_sum += d * d;
   }

   MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

   return sqrt(global_sum);
}

int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   JX_Int nprocs;

#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
   JX_Int nthreads = 1;
#endif

   JX_Int arg_index = 0;
   JX_Int build_matrix_arg_index;
   JX_Int build_rhs_arg_index;
   JX_Int problem_id = -1;
   JX_Int file_base = 13;
   JX_Int TTest = 1;

   char *MatFile = NULL;
   char *RhsFile = NULL;
   char *GusFile = NULL;

   JX_Int *partitioning = NULL;
   JX_Int glosize;

   JX_Real alpha = 1.0;
   JX_Real beta = 2.0;
   JXF_Real alpha_F = 1.0f;
   JXF_Real beta_F = 2.0f;

   JX_Int num_calls = 10;

   jx_hpCSRMatrix *hp_matrix = NULL;
   jx_ParCSRMatrix *par_matrix = NULL;
   jx_ParVector *par_rhs = NULL;
   jx_ParVector *par_x = NULL;
   jx_ParVector *par_rhs_ref = NULL;
   jx_ParVector *par_rhs_v1 = NULL;
   jx_ParVector *par_diff = NULL;

   jxf_hpCSRMatrix *hp_matrix_F = NULL;
   jxf_ParCSRMatrix *par_matrix_F = NULL;
   jxf_ParVector *par_rhs_F = NULL;
   jxf_ParVector *par_x_F = NULL;
   jxf_ParVector *par_rhs_ref_F = NULL;
   jxf_ParVector *par_rhs_v1_F = NULL;
   jxf_ParVector *par_diff_F = NULL;

   jxh_hpCSRMatrix *hp_matrix_H = NULL;
   jxh_ParCSRMatrix *par_matrix_H = NULL;
   jxf_ParVector *par_rhs_ref_H = NULL;
   jxf_ParVector *par_rhs_v1_H = NULL;
   jxf_ParVector *par_rhs_v2_H = NULL;
   jxf_ParVector *par_diff_H = NULL;

   JX_Real norm_ref = 0.0, norm_v1 = 0.0, norm_diff = 0.0, error = 0.0;
   JXF_Real norm_ref_F = 0.0f, norm_v1_F = 0.0f, norm_diff_F = 0.0f, error_F = 0.0f;
   JXF_Real norm_ref_H = 0.0f, norm_v1_H = 0.0f, norm_diff_H = 0.0f, error_H = 0.0f;
   JXF_Real norm_v2_H = 0.0f;

   JX_Real norm_diff_F_vs_D = 0.0;
   JX_Real error_F_vs_D = 0.0;
   JX_Real norm_diff_H_vs_D = 0.0;
   JX_Real error_H_vs_D = 0.0;
   JX_Real norm_diff_Hmul_vs_D = 0.0;
   JX_Real error_Hmul_vs_D = 0.0;

   JX_Real starttime = 0.0, endtime = 0.0;

   jx_MPI_Init(&argc, &argv);
   jx_MPI_Comm_rank(comm, &myid);
   jx_MPI_Comm_size(comm, &nprocs);

#ifdef USING_HWLOC
   jx_hpCreateHardwareInfo(comm);
#endif

   {
      JX_Int cluster_id = myid % 4;
      hthread_dev_open(cluster_id);
      hthread_dat_load(cluster_id, "kernel.dat");
   }

   build_matrix_arg_index = argc;
   build_rhs_arg_index = argc;

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
      else if (strcmp(argv[arg_index], "-pid") == 0)
      {
         arg_index++;
         problem_id = atoi(argv[arg_index++]);
      }
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
      else if (strcmp(argv[arg_index], "-nts") == 0)
      {
         arg_index++;
         nthreads = atoi(argv[arg_index++]);
      }
#endif
      else
      {
         arg_index++;
      }
   }

   //-----------------------
   //  提供线性系统文件
   //-----------------------
   switch (problem_id)
   {
   // case 1: 报错，不可以用
   //    MatFile = "/vol8/home/xtu_ljc/矩阵数据/1-ASIC_680k.bin";
   //    break;
   case 2:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/2-ASIC_680ks.bin";
      break;
   case 3:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/3-Bump_2911.bin";
      break;
   case 4:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/4-Chevron4.bin";
      break;
   case 5:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/5-Cube_Coup_dt0.bin";
      break;
   case 6:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/6-Cube_Coup_dt6.bin";
      break;
   case 7:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/7-CurlCurl_2.bin";
      break;
   case 8:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/8-CurlCurl_3.bin";
      break;
   case 9:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/9-CurlCurl_4.bin";
      break;
   case 10:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/10-Emilia_923.bin";
      break;
   case 11:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/11-Fault_639.bin";
      break;
   case 12:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/12-Flan_1565.bin";
      break;
   case 13:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/13-Freescale1.bin";
      break;
   case 14:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/14-Freescale2.bin";
      break;
   case 15:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/15-FullChip.bin";
      break;
   case 16:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/16-G3_circuit.bin";
      break;
   case 17:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/17-Geo_1438.bin";
      break;
   case 18:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/18-Hardesty1.bin";
      break;
   case 19:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/19-Hook_1498.bin";
      break;
   case 20:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/20-Long_Coup_dt0.bin";
      break;
   case 21:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/21-Long_Coup_dt6.bin";
      break;
   case 22:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/22-PFlow_742.bin";
      break;
   case 23:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/23-Serena.bin";
      break;
   case 24:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/24-StocF-1465.bin";
      break;
   case 25:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/25-adaptive.bin";
      break;
   case 26:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/26-af_shell10.bin";
      break;
   case 27:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/27-amazon-2008.bin";
      break;
   case 28:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/28-apache2.bin";
      break;
   case 29:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/29-as-Skitter.bin";
      break;
   case 30:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/30-audikw_1.bin";
      break;
   case 31:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/31-cage14.bin";
      break;
   case 32:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/32-cage15.bin";
      break;
   case 33:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/33-channel-500x100x100-b050.bin";
      break;
   case 34:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/34-circuit5M.bin";
      break;
   case 35:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/35-circuit5M_dc.bin";
      break;
   case 36:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/36-delaunay_n20.bin";
      break;
   case 37:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/37-delaunay_n21.bin";
      break;
   case 38:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/38-delaunay_n22.bin";
      break;
   case 39:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/39-delaunay_n23.bin";
      break;
   case 40:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/40-delaunay_n24.bin";
      break;
   case 41:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/41-dgreen.bin";
      break;
   case 42:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/42-ecology1.bin";
      break;
   case 43:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/43-ecology2.bin";
      break;
   case 44:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/44-eu-2005.bin";
      break;
   case 45:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/45-flickr.bin";
      break;
   case 46:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/46-hollywood-2009.bin";
      break;
   case 47:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/47-in-2004.bin";
      break;
   case 48:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/48-indochina-2004.bin";
      break;
   case 49:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/49-kkt_power.bin";
      break;
   case 50:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/50-kmer_V2a.bin";
      break;
   case 51:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/51-kron_g500-logn20.bin";
      break;
   case 52:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/52-kron_g500-logn21.bin";
      break;
   case 53:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/53-ldoor.bin";
      break;
   case 54:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/54-ljournal-2008.bin";
      break;
   case 55:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/55-mawi_201512012345.bin";
      break;
   case 56:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/56-mawi_201512020000.bin";
      break;
   case 57:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/57-memchip.bin";
      break;
   case 58:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/58-nlpkkt120.bin";
      break;
   case 59:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/59-nlpkkt160.bin";
      break;
   case 60:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/60-nlpkkt200.bin";
      break;
   case 61:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/61-nlpkkt80.bin";
      break;
   case 62:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/62-nv2.bin";
      break;
   case 63:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/63-packing-500x100x100-b050.bin";
      break;
   case 64:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/64-pre2.bin";
      break;
   case 65:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/65-rajat29.bin";
      break;
   case 66:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/66-rajat30.bin";
      break;
   case 67:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/67-rajat31.bin";
      break;
   case 68:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/68-roadNet-CA.bin";
      break;
   case 69:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/69-roadNet-PA.bin";
      break;
   case 70:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/70-roadNet-TX.bin";
      break;
   case 71:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/71-road_central.bin";
      break;
   case 72:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/72-road_usa.bin";
      break;
   case 73:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/73-soc-LiveJournal1.bin";
      break;
   case 74:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/74-ss.bin";
      break;
   case 75:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/75-t2em.bin";
      break;
   case 76:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/76-tmt_sym.bin";
      break;
   case 77:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/77-uk-2002.bin";
      break;
   case 78:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/78-vas_stokes_1M.bin";
      break;
   case 79:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/79-vas_stokes_2M.bin";
      break;
   case 80:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/80-vas_stokes_4M.bin";
      break;
   case 81:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/81-venturiLevel3.bin";
      break;
   case 82:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/82-wb-edu.bin";
      break;
   case 83:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/83-web-BerkStan.bin";
      break;
   case 84:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/84-web-Google.bin";
      break;
   case 85:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/85-webbase-1M.bin";
      break;
   case 86:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/86-wiki-Talk.bin";
      break;
   case 87:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/87-wikipedia-20051105.bin";
      break;
   case 88:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/88-wikipedia-20060925.bin";
      break;
   case 89:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/89-wikipedia-20061104.bin";
      break;
   case 90:
      MatFile = "/vol8/home/xtu_ljc/矩阵数据/90-wikipedia-20070206.bin";
      break;
   case 101:
      MatFile = "/vol8/home/xtu_ljc/XingXin/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_A.bin";
      // RhsFile = "/vol8/home/xtu_ljc/XingXin/jxpamg-basic-GAI/solverchallenge/solverchallenge23_01/solverchallenge23_01_b.bin";
      file_base = 2;
      break;
   case 102:
      MatFile = "/vol8/home/xtu_ljc/XingXin/jxpamg-basic-GAI/solverchallenge/solverchallenge23_02/solverchallenge23_02_A.bin";
      // RhsFile = "/vol8/home/xtu_ljc/XingXin/jxpamg-basic-GAI/solverchallenge/solverchallenge23_02/solverchallenge23_02_b.bin";
      file_base = 2;
      break;
   case 103:
      MatFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/mat_csr_128X128X128.bin";
      // RhsFile = "/vol8/home/xtu_ljc/XX-test/jxpamg-basic-xingxin/matrix/Constant/rhs_128X128X128.bin";
      file_base = 2;
      break;
   }

   if (build_matrix_arg_index < argc)
   {
      MatFile = argv[build_matrix_arg_index];
      file_base = 0;
   }
   else if (problem_id < 0 || MatFile == NULL)
   {
      jx_printf("Error: No matrix filename specified.\n");
      exit(1);
   }

   if (build_rhs_arg_index < argc)
   {
      RhsFile = argv[build_rhs_arg_index];
   }

   if (myid == 0)
   {
      jx_printf("\n\n+++++++++++++++++++++\n matrix:\n %s\n+++++++++++++++++++++\n", MatFile);
      jx_printf("\n+++++++++++++++++++++\n rhs:\n %s\n+++++++++++++++++++++\n", RhsFile);
      jx_printf("\n+++++++++++++++++++++");
#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
      jx_printf(" With OpenMP using %d threads,", nthreads);
#endif
      jx_printf(" MPI using %d processors +++++++++++++++++++++\n\n", nprocs);
   }

#if JX_USING_OPENMP || defined(JX_USING_PGCC_SMP)
   omp_set_num_threads(nthreads);
#endif

   if (TTest)
      starttime = jx_MPI_Wtime();

   hp_matrix = jx_hpBuildMatParFromOneFile(MatFile, 1, file_base);
   par_matrix = jx_hpCSRMatrixPar(hp_matrix);

   if (jx_ParCSRMatrixSpMVPrecondFP64Create(par_matrix, myid))
   {
      if (myid == 0)
         jx_printf("Error: jx_ParCSRMatrixSpMVPrecondFP64Create failed.\n");
      exit(1);
   }

   if (myid == 0)
      jx_printf("Info: Generating a new RHS vector...\n");

   {
      JX_Int rhs_glosize = jx_ParCSRMatrixGlobalNumRows(par_matrix);
      JX_Int *rhs_partitioning = jx_ParCSRMatrixRowStarts(par_matrix);

      par_rhs = jx_ParVectorCreate(comm, rhs_glosize, rhs_partitioning);
      jx_ParVectorSetPartitioningOwner(par_rhs, 0);
      jx_ParVectorInitialize(par_rhs);
      jx_ParVectorSetConstantValues(par_rhs, 2.0);
   }

   if (!GusFile)
   {
      glosize = jx_ParVectorGlobalSize(par_rhs);
      partitioning = jx_ParVectorPartitioning(par_rhs);

      par_x = jx_ParVectorCreate(comm, glosize, partitioning);
      jx_ParVectorSetPartitioningOwner(par_x, 0);
      jx_ParVectorInitialize(par_x);
      jx_ParVectorSetRandomValues(par_x, 42);
   }
   else
   {
      par_x = jx_hpBuildRhsParFromOneFile(GusFile, hp_matrix, 0);
   }

   glosize = jx_ParVectorGlobalSize(par_rhs);
   partitioning = jx_ParVectorPartitioning(par_rhs);

   par_rhs_ref = jx_ParVectorCreate(comm, glosize, partitioning);
   jx_ParVectorSetPartitioningOwner(par_rhs_ref, 0);
   jx_ParVectorInitialize(par_rhs_ref);
   jx_ParVectorSetConstantValues(par_rhs_ref, 2.0);

   par_rhs_v1 = jx_ParVectorCreate(comm, glosize, partitioning);
   jx_ParVectorSetPartitioningOwner(par_rhs_v1, 0);
   jx_ParVectorInitialize(par_rhs_v1);
   jx_ParVectorSetConstantValues(par_rhs_v1, 2.0);

   par_diff = jx_ParVectorCreate(comm, glosize, partitioning);
   jx_ParVectorSetPartitioningOwner(par_diff, 0);
   jx_ParVectorInitialize(par_diff);
   jx_ParVectorSetConstantValues(par_diff, 0.0);

   if (myid == 0)
      jx_printf(" 将并行矩阵和并行右端以及并行初始迭代向量转成低精度版本+++++ \n\n");

   par_matrix_F = jxf_ParCSRMatrixCreate(comm,
                                         jx_ParCSRMatrixGlobalNumRows(par_matrix),
                                         jx_ParCSRMatrixGlobalNumCols(par_matrix),
                                         jx_ParCSRMatrixRowStarts(par_matrix),
                                         jx_ParCSRMatrixColStarts(par_matrix),
                                         jx_CSRMatrixNumCols(jx_ParCSRMatrixOffd(par_matrix)),
                                         jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixDiag(par_matrix)),
                                         jx_CSRMatrixNumNonzeros(jx_ParCSRMatrixOffd(par_matrix)));

   jxf_ParCSRMatrixInitialize(par_matrix_F);
   jxmp_ParCSRMatrixDtoF(par_matrix, par_matrix_F, 1);

   if (jxf_ParCSRMatrixSpMVPrecondFP32Create(par_matrix_F, myid))
   {
      if (myid == 0)
         jx_printf("Error: jxf_ParCSRMatrixSpMVPrecondFP32Create failed.\n");
      exit(1);
   }

   hp_matrix_F = jxf_hpInithpCSRMatrix();
   jxf_hpCSRMatrixPar(hp_matrix_F) = par_matrix_F;

   par_matrix_H = jxh_ParCSRMatrixCreate(comm,
                                         jxf_ParCSRMatrixGlobalNumRows(par_matrix_F),
                                         jxf_ParCSRMatrixGlobalNumCols(par_matrix_F),
                                         jxf_ParCSRMatrixRowStarts(par_matrix_F),
                                         jxf_ParCSRMatrixColStarts(par_matrix_F),
                                         jxf_CSRMatrixNumCols(jxf_ParCSRMatrixOffd(par_matrix_F)),
                                         jxf_CSRMatrixNumNonzeros(jxf_ParCSRMatrixDiag(par_matrix_F)),
                                         jxf_CSRMatrixNumNonzeros(jxf_ParCSRMatrixOffd(par_matrix_F)));

   jxh_ParCSRMatrixOwnsRowStarts(par_matrix_H) = 0;
   jxh_ParCSRMatrixOwnsColStarts(par_matrix_H) = 0;

   jxh_ParCSRMatrixInitialize(par_matrix_H);
   jxf_ParCSRMatrixFtoH(par_matrix_F, par_matrix_H, 1);
   jxh_ParCSRMatrixSetNumNonzeros(par_matrix_H);
   jxh_ParCSRMatrixSetDNumNonzeros(par_matrix_H);

   if (jxh_ParCSRMatrixSpMVPrecondFP16Create(par_matrix_H, myid))
   {
      if (myid == 0)
         jx_printf("Error: jxh_ParCSRMatrixSpMVPrecondFP16Create failed.\n");
      exit(1);
   }

   hp_matrix_H = jxh_hpInithpCSRMatrix();
   jxh_hpCSRMatrixPar(hp_matrix_H) = par_matrix_H;

#define CREATE_F_VECTOR(vec, size_func, starts_func)              \
   do                                                             \
   {                                                              \
      vec = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(hp_matrix_F), \
                                size_func(hp_matrix_F),           \
                                starts_func(hp_matrix_F));        \
      jxf_ParVectorInitialize(vec);                               \
   } while (0)

   CREATE_F_VECTOR(par_rhs_F, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxmp_ParVectorDtoF(par_rhs, par_rhs_F);

   CREATE_F_VECTOR(par_rhs_ref_F, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxmp_ParVectorDtoF(par_rhs_ref, par_rhs_ref_F);

   CREATE_F_VECTOR(par_x_F, jxf_hpCSRMatrixGlobalNumCols, jxf_hpCSRMatrixColStarts);
   jxmp_ParVectorDtoF(par_x, par_x_F);

   CREATE_F_VECTOR(par_rhs_v1_F, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxmp_ParVectorDtoF(par_rhs_v1, par_rhs_v1_F);

   CREATE_F_VECTOR(par_diff_F, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxf_ParVectorSetConstantValues(par_diff_F, 0.0f);

   CREATE_F_VECTOR(par_rhs_ref_H, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxf_ParVectorSetConstantValues(par_rhs_ref_H, 2.0f);

   CREATE_F_VECTOR(par_rhs_v1_H, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxf_ParVectorSetConstantValues(par_rhs_v1_H, 2.0f);

   CREATE_F_VECTOR(par_rhs_v2_H, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxf_ParVectorSetConstantValues(par_rhs_v2_H, 2.0f);

   CREATE_F_VECTOR(par_diff_H, jxf_hpCSRMatrixGlobalNumRows, jxf_hpCSRMatrixRowStarts);
   jxf_ParVectorSetConstantValues(par_diff_H, 0.0f);

#undef CREATE_F_VECTOR

   endtime = jx_MPI_Wtime();
   (void)starttime;
   (void)endtime;

   /* ============================================================
    * FP64 reference
    * ============================================================ */
   if (myid == 0)
      jx_printf("\n================ FP64 reference ================\n");

   jx_spmv_type = 0;
   jx_ParVectorSetConstantValues(par_rhs_ref, 2.0);
   jx_ParCSRMatrixMatvec(alpha, par_matrix, par_x, beta, par_rhs_ref);
   norm_ref = jx_ParVectorNorm2(par_rhs_ref);

   if (myid == 0)
   {
      jx_printf("\n[FP64 reference]\n");
      jx_printf("  jx_spmv_type : %d\n", jx_spmv_type);
      jx_printf("  ||ref||_2    : %.12e\n", norm_ref);
   }

   /* FP64 timing */
   {
      JX_Real total = 0.0, begin, end, local_time, max_time;

      jx_spmv_type = 1;

      for (JX_Int call = 0; call < num_calls; call++)
      {
         jx_ParVectorSetConstantValues(par_rhs_v1, 2.0);
         MPI_Barrier(comm);

         begin = jx_MPI_Wtime();
         jx_ParCSRMatrixMatvec(alpha, par_matrix, par_x, beta, par_rhs_v1);
         end = jx_MPI_Wtime();

         local_time = end - begin;
         MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
         total += max_time;
      }

      total /= num_calls;

      if (myid == 0)
      {
         jx_printf("\n----------------------------\n");
         jx_printf("FP64 method     : block2 + dma_sg\n");
         jx_printf("jx_spmv_type    : %d\n", jx_spmv_type);
         jx_printf("FP64 SpMV_total : %e seconds\n", total);
         jx_printf("----------------------------\n");
      }
   }

   /* ============================================================
    * FP32 candidate vs FP64
    * ============================================================ */
   if (myid == 0)
      printf("\n================ FP32 vs FP64 SpMV error ================\n");

   jxf_spmv_type = 1;

   jxf_ParVectorSetConstantValues(par_rhs_v1_F, 2.0f);
   jxf_ParCSRMatrixMatvec(alpha_F, par_matrix_F, par_x_F, beta_F, par_rhs_v1_F);

   norm_v1_F = jxf_ParVectorNorm2(par_rhs_v1_F);
   norm_diff_F_vs_D = diff_norm2_f32_vs_fp64(par_rhs_ref, par_rhs_v1_F, comm);
   error_F_vs_D = norm_diff_F_vs_D / norm_ref;

   if (myid == 0)
   {
      printf("\n[FP32 block2 + dma_sg]\n");
      printf("  jxf_spmv_type          : %d\n", jxf_spmv_type);
      printf("  ||v1_F||_2             : %.12e\n", (double)norm_v1_F);
      printf("  ||FP64-ref - v1_F||_2  : %.12e\n", (double)norm_diff_F_vs_D);
      printf("  relative error vs FP64 : %.12e\n", (double)error_F_vs_D);
   }

   {
      JX_Real total = 0.0, begin, end, local_time, max_time;

      jxf_spmv_type = 1;

      for (JX_Int call = 0; call < num_calls; call++)
      {
         jxf_ParVectorSetConstantValues(par_rhs_v1_F, 2.0f);
         MPI_Barrier(comm);

         begin = jxf_MPI_Wtime();
         jxf_ParCSRMatrixMatvec(alpha_F, par_matrix_F, par_x_F, beta_F, par_rhs_v1_F);
         end = jxf_MPI_Wtime();

         local_time = end - begin;
         MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
         total += max_time;
      }

      total /= num_calls;

      if (myid == 0)
      {
         printf("\n----------------------------\n");
         printf("FP32 method       : block2 + dma_sg\n");
         printf("jxf_spmv_type     : %d\n", jxf_spmv_type);
         printf("FP32 SpMV_total   : %e seconds\n", total);
         printf("----------------------------\n");
      }
   }

   /* ============================================================
    * HF32/SPLIT candidate vs FP64
    * ============================================================ */
   if (myid == 0)
      printf("\n================ HF32/SPLIT vs FP64 SpMV error ================\n");

   mix_spmv_type = 0;
   jxh_spmv_type = 1;

   jxf_ParVectorSetConstantValues(par_rhs_v1_H, 2.0f);
   jxh_ParCSRMatrixMatvec(alpha_F, par_matrix_H, par_x_F, beta_F, par_rhs_v1_H);

   norm_v1_H = jxf_ParVectorNorm2(par_rhs_v1_H);
   norm_diff_H_vs_D = diff_norm2_f32_vs_fp64(par_rhs_ref, par_rhs_v1_H, comm);
   error_H_vs_D = norm_diff_H_vs_D / norm_ref;

   if (myid == 0)
   {
      printf("\n[HF32/SPLIT GSM block4 + dma_sg]\n");
      printf("  mix_spmv_type          : %d\n", mix_spmv_type);
      printf("  jxh_spmv_type          : %d\n", jxh_spmv_type);
      printf("  ||v1_H||_2             : %.12e\n", (double)norm_v1_H);
      printf("  ||FP64-ref - v1_H||_2  : %.12e\n", (double)norm_diff_H_vs_D);
      printf("  relative error vs FP64 : %.12e\n", (double)error_H_vs_D);
   }

   {
      JX_Real total = 0.0, begin, end, local_time, max_time;

      mix_spmv_type = 0;
      jxh_spmv_type = 1;

      for (JX_Int call = 0; call < num_calls; call++)
      {
         jxf_ParVectorSetConstantValues(par_rhs_v1_H, 2.0f);
         MPI_Barrier(comm);

         begin = jxf_MPI_Wtime();
         jxh_ParCSRMatrixMatvec(alpha_F, par_matrix_H, par_x_F, beta_F, par_rhs_v1_H);
         end = jxf_MPI_Wtime();

         local_time = end - begin;
         MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
         total += max_time;
      }

      total /= num_calls;

      if (myid == 0)
      {
         printf("\n----------------------------\n");
         printf("HF32/PAIRSUM method : FP16 storage + FP16 x + FP32 product pair-sum block4 + dma_sg\n");
         printf("mix_spmv_type     : %d\n", mix_spmv_type);
         printf("jxh_spmv_type     : %d\n", jxh_spmv_type);
         printf("HF32/SPLIT total  : %e seconds\n", total);
         printf("----------------------------\n");
      }
   }

   /* ============================================================
    * HmulFacc candidate vs FP64
    * ============================================================ */
   if (myid == 0)
      printf("\n================ HmulFacc vs FP64 SpMV error ================\n");

   mix_spmv_type = 1;
   jxh_spmv_type = 1;

   jxf_ParVectorSetConstantValues(par_rhs_v2_H, 2.0f);
   jxh_ParCSRMatrixMatvec(alpha_F, par_matrix_H, par_x_F, beta_F, par_rhs_v2_H);

   norm_v2_H = jxf_ParVectorNorm2(par_rhs_v2_H);
   norm_diff_Hmul_vs_D = diff_norm2_f32_vs_fp64(par_rhs_ref, par_rhs_v2_H, comm);
   error_Hmul_vs_D = norm_diff_Hmul_vs_D / norm_ref;

   if (myid == 0)
   {
      printf("\n[HmulFacc GSM block4 host-reduce]\n");
      printf("  mix_spmv_type          : %d\n", mix_spmv_type);
      printf("  jxh_spmv_type          : %d\n", jxh_spmv_type);
      printf("  ||v2_H||_2             : %.12e\n", (double)norm_v2_H);
      printf("  ||FP64-ref - v2_H||_2  : %.12e\n", (double)norm_diff_Hmul_vs_D);
      printf("  relative error vs FP64 : %.12e\n", (double)error_Hmul_vs_D);
   }

   {
      JX_Real total = 0.0, begin, end, local_time, max_time;

      mix_spmv_type = 1;
      jxh_spmv_type = 1;

      for (JX_Int call = 0; call < num_calls; call++)
      {
         jxf_ParVectorSetConstantValues(par_rhs_v2_H, 2.0f);
         MPI_Barrier(comm);

         begin = jxf_MPI_Wtime();
         jxh_ParCSRMatrixMatvec(alpha_F, par_matrix_H, par_x_F, beta_F, par_rhs_v2_H);
         end = jxf_MPI_Wtime();

         local_time = end - begin;
         MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
         total += max_time;
      }

      total /= num_calls;

      if (myid == 0)
      {
         printf("\n----------------------------\n");
         printf("HmulFacc method  : FP16 storage + FP16 multiply + FP32 host reduce\n");
         printf("mix_spmv_type    : %d\n", mix_spmv_type);
         printf("jxh_spmv_type    : %d\n", jxh_spmv_type);
         printf("HmulFacc total   : %e seconds\n", total);
         printf("----------------------------\n");
      }
   }

#ifdef USING_HWLOC
   jx_hpCSRhardwareDestroy();
#endif

   jx_hpCSRMatrixDestroy(hp_matrix);
   jx_ParVectorDestroy(par_rhs);
   jx_ParVectorDestroy(par_x);
   jx_ParVectorDestroy(par_rhs_ref);
   jx_ParVectorDestroy(par_rhs_v1);
   jx_ParVectorDestroy(par_diff);

   jxf_ParVectorDestroy(par_rhs_F);
   jxf_ParVectorDestroy(par_x_F);
   jxf_ParVectorDestroy(par_rhs_ref_F);
   jxf_ParVectorDestroy(par_rhs_v1_F);
   jxf_ParVectorDestroy(par_diff_F);

   jxf_ParVectorDestroy(par_rhs_ref_H);
   jxf_ParVectorDestroy(par_rhs_v1_H);
   jxf_ParVectorDestroy(par_rhs_v2_H);
   jxf_ParVectorDestroy(par_diff_H);

   /*
    * 注意：par_matrix_H 复用了 par_matrix_F 的 row/col starts，
    * owns flags 已经置 0，因此可以正常 destroy。
    */
   if (hp_matrix_H)
      jxh_hpCSRMatrixDestroy(hp_matrix_H);

   if (hp_matrix_F)
      jxf_hpCSRMatrixDestroy(hp_matrix_F);

   jx_MPI_Finalize();

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
   // jx_printf("\n >>> jxmp_CSRMatrixDtoF \n");
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
   // jx_printf("\n >>> jxmp_ParCSRMatrixDtoF \n");
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
#include "../jxpamg-basic-DGAI/include/jx_smp_forloop.h"
   for (i = 0; i < size; i++)
   {
      y_data[i] = (JXF_Real)x_data[i]; // 双精度到单精度转换
   }

   return 0;
}