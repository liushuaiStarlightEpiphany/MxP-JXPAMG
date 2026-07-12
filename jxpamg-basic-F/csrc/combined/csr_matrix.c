//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  csr_matrix.c -- operations for CSR matrices.
 *  Date: 2014/07/05
 */

#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_CSRMatrixWeaklyDiagDominant
 * \brief Weakly Diag-Dominant
 * \author Yue Xiaoqiang
 * \date 2014/01/03
 */
JXF_Int
jxf_CSRMatrixWeaklyDiagDominant( jxf_CSRMatrix *A,
                                JXF_Int num_equns,
                                JXF_Real theta,
                                JXF_Real gamma_3,
                                JXF_Real gamma_11 )
{
   JXF_Int ierr = 1;
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int block_size = num_rows / num_equns;
   JXF_Int t_b_size = 2 * block_size;
   JXF_Int f_b_size = 3 * block_size;
   JXF_Int tt_b_size = 10 * block_size;
   JXF_Int ff_b_size = 11 * block_size;
   JXF_Int *IA = jxf_CSRMatrixI(A);
   JXF_Int *JA = jxf_CSRMatrixJ(A);
   JXF_Real *AA = jxf_CSRMatrixData(A);
   JXF_Int cnt = 0, gnt = 0;
   JXF_Int row, row_srt, row_end, col;
   JXF_Real row_sum;
   
   // block - 11
   for (row = tt_b_size; row < ff_b_size; row ++)
   {
      row_sum = 0.0;
      row_srt = IA[row];
      row_end = IA[row+1];
      for (col = row_srt; col < row_end; col ++)
      {
         if ((JA[col] >= tt_b_size) && (JA[col] < ff_b_size))
         {
            row_sum += AA[col];
         }
      }
      if (fabs(row_sum) <= (theta * fabs(AA[row_srt])))
      {
         gnt ++;
      }
   }
   if ((JXF_Real)gnt/block_size < gamma_11)
   {
      return ierr;
   }
   // block - 3
   if (gamma_3 > 0.0)
   {
      for (row = t_b_size; row < f_b_size; row ++)
      {
         row_sum = 0.0;
         row_srt = IA[row];
         row_end = IA[row+1];
         for (col = row_srt; col < row_end; col ++)
         {
            if ((JA[col] >= t_b_size) && (JA[col] < f_b_size))
            {
               row_sum += AA[col];
            }
         }
         if (fabs(row_sum) <= (theta * fabs(AA[row_srt])))
         {
            cnt ++;
         }
      }
      if ((JXF_Real)cnt/block_size >= gamma_3)
      {
          ierr = 0;
      }
   }
   else
   {
       ierr = 0;
   }
   
   return ierr;
}
