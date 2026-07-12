//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_interp_clsfix.c
 *
 *  两阶段插值构建：
 *    Stage 1: 对所有 F 点使用经典修正插值 (interp_type=0)
 *             构建完整 P 矩阵 P_classical
 *    Stage 2: 对 eta_i >= eta_threshold 的 F 点，
 *             用扩展插值 (interp_type=6) 重新计算该行并覆盖
 *
 *  实现方式：
 *    1. 调用 jx_PAMGBuildInterp 得到纯经典 P
 *    2. 调用 jx_PAMGBuildExtPIInterp 得到纯扩展 P
 *    3. 对 EI 点，将扩展 P 的行复制覆盖经典 P 的对应行
 *
 *  Date: 2026/07/10
 */

#include "jx_pamg.h"

/* 在 P 矩阵中扩展一行，为其增加 extra 个非零元 */
static JX_Int
expand_row_in_P(jx_CSRMatrix *diag, jx_CSRMatrix *offd,
                JX_Int row, JX_Int extra_diag, JX_Int extra_offd)
{
   JX_Int *diag_i = jx_CSRMatrixI(diag);
   JX_Int *offd_i = offd ? jx_CSRMatrixI(offd) : NULL;
   JX_Int old_diag_nz = jx_CSRMatrixNumNonzeros(diag);
   JX_Int old_offd_nz = offd ? jx_CSRMatrixNumNonzeros(offd) : 0;
   JX_Int new_diag_nz = old_diag_nz + extra_diag;
   JX_Int new_offd_nz = old_offd_nz + extra_offd;

   if (extra_diag > 0)
   {
      JX_Int r = jx_CSRMatrixNumRows(diag);
      JX_Int shift_start = diag_i[row + 1];
      JX_Int shift_len = old_diag_nz - shift_start;

      jx_TReAlloc(jx_CSRMatrixJ(diag), JX_Int, new_diag_nz);
      jx_TReAlloc(jx_CSRMatrixData(diag), JX_Real, new_diag_nz);

      /* 将 row+1 之后的数据后移 extra_diag 位 */
      for (JX_Int k = shift_len - 1; k >= 0; k--)
      {
         jx_CSRMatrixJ(diag)[shift_start + k + extra_diag] = jx_CSRMatrixJ(diag)[shift_start + k];
         jx_CSRMatrixData(diag)[shift_start + k + extra_diag] = jx_CSRMatrixData(diag)[shift_start + k];
      }

      /* 更新行指针 */
      for (JX_Int k = row + 1; k <= r; k++)
         diag_i[k] += extra_diag;
   }

   if (extra_offd > 0 && offd)
   {
      JX_Int r = jx_CSRMatrixNumRows(offd);
      JX_Int shift_start = offd_i[row + 1];
      JX_Int shift_len = old_offd_nz - shift_start;

      jx_TReAlloc(jx_CSRMatrixJ(offd), JX_Int, new_offd_nz);
      jx_TReAlloc(jx_CSRMatrixData(offd), JX_Real, new_offd_nz);

      for (JX_Int k = shift_len - 1; k >= 0; k--)
      {
         jx_CSRMatrixJ(offd)[shift_start + k + extra_offd] = jx_CSRMatrixJ(offd)[shift_start + k];
         jx_CSRMatrixData(offd)[shift_start + k + extra_offd] = jx_CSRMatrixData(offd)[shift_start + k];
      }

      for (JX_Int k = row + 1; k <= r; k++)
         offd_i[k] += extra_offd;
   }

   return 0;
}

/*!
 * \fn JX_Int jx_PAMGBuildClassicalFix
 *    经典修正插值 + 扩展插值修复
 *
 *  Stage 1: 对所有 F 点使用经典修正插值 (interp_type=0)
 *  Stage 2: 对 eta_i >= eta_threshold 的 F 点，
 *           使用扩展插值 (interp_type=6) 修复该行
 */
JX_Int
jx_PAMGBuildClassicalFix(jx_ParCSRMatrix *par_A,
                         JX_Int *CF_marker,
                         jx_ParCSRMatrix *par_S,
                         JX_Real *eta,
                         JX_Real eta_threshold,
                         JX_Int *num_cpts_global,
                         JX_Int num_functions,
                         JX_Int *dof_func,
                         JX_Int debug_flag,
                         JX_Real trunc_factor,
                         JX_Int max_elmts,
                         JX_Int *col_offd_S_to_A,
                         jx_ParCSRMatrix **P_ptr,
                         JX_Int *num_fpts_global_ptr,
                         JX_Int *num_ei_fpts_global_ptr)
{
   MPI_Comm comm = jx_ParCSRMatrixComm(par_A);
   JX_Int my_id, num_procs;
   jx_MPI_Comm_rank(comm, &my_id);
   jx_MPI_Comm_size(comm, &num_procs);

   JX_Int n_fine = jx_CSRMatrixNumRows(jx_ParCSRMatrixDiag(par_A));
   JX_Int local_num_fpts = 0;
   JX_Int local_num_ei_fpts = 0;

   /*--------------------------------------------------------------
    *  Stage 1: 构建经典修正插值 P_classical
    *           使用 interp_type=0 的逻辑
    *--------------------------------------------------------------*/
   jx_ParCSRMatrix *P_cls = NULL;
   jx_PAMGBuildInterp(par_A, CF_marker, par_S,
                       num_cpts_global, num_functions,
                       dof_func, debug_flag,
                       trunc_factor, max_elmts,
                       col_offd_S_to_A, &P_cls);

   /*--------------------------------------------------------------
    *  Stage 2: 对 EI 点, 使用扩展插值修复
    *--------------------------------------------------------------*/
   jx_ParCSRMatrix *P_ext = NULL;
   jx_PAMGBuildExtPIInterp(par_A, CF_marker, par_S,
                            num_cpts_global, num_functions,
                            dof_func, debug_flag,
                            trunc_factor, max_elmts,
                            col_offd_S_to_A, &P_ext);

   jx_CSRMatrix *cls_diag = jx_ParCSRMatrixDiag(P_cls);
   jx_CSRMatrix *cls_offd = jx_ParCSRMatrixOffd(P_cls);
   JX_Int *cls_diag_i = jx_CSRMatrixI(cls_diag);
   JX_Int *cls_offd_i = cls_offd ? jx_CSRMatrixI(cls_offd) : NULL;
   JX_Int *cls_diag_j = jx_CSRMatrixJ(cls_diag);
   JX_Real *cls_diag_data = jx_CSRMatrixData(cls_diag);
   JX_Int *cls_offd_j = cls_offd ? jx_CSRMatrixJ(cls_offd) : NULL;
   JX_Real *cls_offd_data = cls_offd ? jx_CSRMatrixData(cls_offd) : NULL;

   jx_CSRMatrix *ext_diag = jx_ParCSRMatrixDiag(P_ext);
   jx_CSRMatrix *ext_offd = jx_ParCSRMatrixOffd(P_ext);
   JX_Int *ext_diag_i = jx_CSRMatrixI(ext_diag);
   JX_Int *ext_offd_i = ext_offd ? jx_CSRMatrixI(ext_offd) : NULL;
   JX_Int *ext_diag_j = jx_CSRMatrixJ(ext_diag);
   JX_Real *ext_diag_data = jx_CSRMatrixData(ext_diag);
   JX_Int *ext_offd_j = ext_offd ? jx_CSRMatrixJ(ext_offd) : NULL;
   JX_Real *ext_offd_data = ext_offd ? jx_CSRMatrixData(ext_offd) : NULL;

   /* 关闭扩展 P 的列映射所有权 */
   jx_ParCSRMatrixColMapOffd(P_ext) = NULL;

   for (JX_Int i = 0; i < n_fine; i++)
   {
      if (CF_marker[i] >= 0 || CF_marker[i] == -3)
         continue;

      local_num_fpts++;

      if (eta == NULL || eta[i] < eta_threshold)
         continue;

      local_num_ei_fpts++;

      /* 该点需要修复: 用扩展行覆盖经典行 */
      JX_Int cls_d_start = cls_diag_i[i];
      JX_Int cls_d_end = cls_diag_i[i + 1];
      JX_Int cls_d_len = cls_d_end - cls_d_start;

      JX_Int ext_d_start = ext_diag_i[i];
      JX_Int ext_d_end = ext_diag_i[i + 1];
      JX_Int ext_d_len = ext_d_end - ext_d_start;

      JX_Int cls_o_len = 0, ext_o_len = 0;
      JX_Int cls_o_start = 0, ext_o_start = 0;

      if (cls_offd_i && ext_offd_i)
      {
         cls_o_start = cls_offd_i[i];
         cls_o_len = cls_offd_i[i + 1] - cls_o_start;
         ext_o_start = ext_offd_i[i];
         ext_o_len = ext_offd_i[i + 1] - ext_o_start;
      }

      JX_Int extra_d = ext_d_len - cls_d_len;
      JX_Int extra_o = ext_o_len - cls_o_len;

      /* 在经典 P 中为额外非零元腾出空间 */
      if (extra_d > 0 || extra_o > 0)
         expand_row_in_P(cls_diag, cls_offd, i, extra_d > 0 ? extra_d : 0, extra_o > 0 ? extra_o : 0);

      /* 重新获取行指针（expand后可能已变化） */
      cls_d_start = cls_diag_i[i];
      cls_d_end = cls_diag_i[i + 1];
      if (cls_offd_i)
      {
         cls_o_start = cls_offd_i[i];
         cls_o_len = cls_offd_i[i + 1] - cls_o_start;
      }

      /* 复制扩展行数据到经典行 */
      for (JX_Int k = 0; k < ext_d_len; k++)
      {
         cls_diag_j[cls_d_start + k] = ext_diag_j[ext_d_start + k];
         cls_diag_data[cls_d_start + k] = ext_diag_data[ext_d_start + k];
      }

      if (cls_offd_j && ext_offd_j)
      {
         for (JX_Int k = 0; k < ext_o_len; k++)
         {
            cls_offd_j[cls_o_start + k] = ext_offd_j[ext_o_start + k];
            cls_offd_data[cls_o_start + k] = ext_offd_data[ext_o_start + k];
         }
      }
   }

   /* 释放扩展 P, 返回 P_classical */
   if (P_ext)
   {
      jx_ParCSRMatrixColMapOffd(P_ext) = NULL;
      jx_ParCSRMatrixDestroy(P_ext);
   }

   *P_ptr = P_cls;

   if (num_fpts_global_ptr)
   {
      JX_Int global_fpts = 0;
      jx_MPI_Allreduce(&local_num_fpts, &global_fpts, 1, JX_MPI_INT, MPI_SUM, comm);
      *num_fpts_global_ptr = global_fpts;
   }

   if (num_ei_fpts_global_ptr)
   {
      JX_Int global_ei = 0;
      jx_MPI_Allreduce(&local_num_ei_fpts, &global_ei, 1, JX_MPI_INT, MPI_SUM, comm);
      *num_ei_fpts_global_ptr = global_ei;
   }

   return 0;
}
