//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  seq_csr_matvec.c -- basic operations for mat-vec multiplication.
 *  Date: 2011/09/03
 */

#include "jxf_mv.h"
#include "hthread_host.h"

extern JXF_Int myid;
extern JXF_Int coreNums;
extern JXF_Int jxf_spmv_type;
extern JXF_Int jxh_spmv_type;

struct jxf_SpMVPrecondFP32Data_struct
{
   JXF_Int cluster_id;

   JXF_Int num_rows;
   JXF_Int num_cols;
   JXF_Int A_nnz;

   JXF_Int num_blocks;
   JXF_Int num_tasks;
   JXF_Int x_pair_size;

   JXF_Int *A_i_block;
   JXF_Int *A_j_block;
   JXF_Real *A_data_block;

   JXF_Int *task_row_bounds;
   JXF_Real *x_pair;
   JXF_Real *dy_data;
};

static void
jxf_SpMVPrecondFP32DataDestroy(jxf_SpMVPrecondFP32Data *pre)
{
   if (pre == NULL)
      return;

   if (pre->A_i_block != NULL)
      hthread_free(pre->A_i_block);

   if (pre->A_j_block != NULL)
      hthread_free(pre->A_j_block);

   if (pre->A_data_block != NULL)
      hthread_free(pre->A_data_block);

   if (pre->task_row_bounds != NULL)
      hthread_free(pre->task_row_bounds);

   if (pre->x_pair != NULL)
      hthread_free(pre->x_pair);

   if (pre->dy_data != NULL)
      hthread_free(pre->dy_data);

   free(pre);
}

static jxf_SpMVPrecondFP32Data *
jxf_SpMVPrecondFP32DataCreate(jxf_CSRMatrix *A,
                              JXF_Int x_size,
                              JXF_Int myid)
{
   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A);
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   JXF_Int A_nnz = jxf_CSRMatrixNumNonzeros(A);

   JXF_Int cluster_id = myid % 4;

   JXF_Int AM_size = (96000 / 3) * 2;
   JXF_Int max_blocks_per_task = (AM_size / 2) & ~1;

   jxf_SpMVPrecondFP32Data *pre =
       (jxf_SpMVPrecondFP32Data *)malloc(sizeof(jxf_SpMVPrecondFP32Data));

   if (pre == NULL)
      return NULL;

   memset(pre, 0, sizeof(jxf_SpMVPrecondFP32Data));

   pre->cluster_id = cluster_id;
   pre->num_rows = num_rows;
   pre->num_cols = num_cols;
   pre->A_nnz = A_nnz;

   pre->x_pair_size = x_size;
   if (pre->x_pair_size % 2 != 0)
      pre->x_pair_size++;

   pre->A_i_block =
       (JXF_Int *)hthread_malloc(cluster_id,
                                 sizeof(JXF_Int) * (num_rows + 1),
                                 HT_MEM_RW);

   pre->A_j_block =
       (JXF_Int *)hthread_malloc(cluster_id,
                                 sizeof(JXF_Int) * A_nnz,
                                 HT_MEM_RW);

   pre->A_data_block =
       (JXF_Real *)hthread_malloc(cluster_id,
                                  sizeof(JXF_Real) * 2 * A_nnz,
                                  HT_MEM_RW);

   pre->task_row_bounds =
       (JXF_Int *)hthread_malloc(cluster_id,
                                 sizeof(JXF_Int) * (num_rows + 1),
                                 HT_MEM_RW);

   pre->x_pair =
       (JXF_Real *)hthread_malloc(cluster_id,
                                  sizeof(JXF_Real) * pre->x_pair_size,
                                  HT_MEM_RW);

   if (pre->A_i_block == NULL || pre->A_j_block == NULL ||
       pre->A_data_block == NULL || pre->task_row_bounds == NULL ||
       pre->x_pair == NULL)
   {
      jxf_SpMVPrecondFP32DataDestroy(pre);
      return NULL;
   }

   if (x_size < pre->x_pair_size)
   {
      pre->x_pair[x_size] = 0.0f;
   }

   JXF_Int block_pos = 0;
   pre->A_i_block[0] = 0;

   for (JXF_Int i = 0; i < num_rows; i++)
   {
      JXF_Int row_block_start = block_pos;

      for (JXF_Int jj = A_i[i]; jj < A_i[i + 1]; jj++)
      {
         JXF_Int col = A_j[jj];
         JXF_Int block_col = col / 2;
         JXF_Int lane = col & 1;
         JXF_Int offset = block_col * 8;
         JXF_Int found = -1;

         for (JXF_Int b = row_block_start; b < block_pos; b++)
         {
            if (pre->A_j_block[b] == offset)
            {
               found = b;
               break;
            }
         }

         if (found < 0)
         {
            found = block_pos++;
            pre->A_j_block[found] = offset;
            pre->A_data_block[2 * found] = 0.0f;
            pre->A_data_block[2 * found + 1] = 0.0f;
         }

         pre->A_data_block[2 * found + lane] += A_data[jj];
      }

      pre->A_i_block[i + 1] = block_pos;
   }

   pre->num_blocks = block_pos;

   JXF_Int size = 0;
   JXF_Int start_row = 0;

   pre->task_row_bounds[size++] = 0;

   while (start_row < num_rows)
   {
      JXF_Int target_blocks = pre->A_i_block[start_row] + max_blocks_per_task;
      JXF_Int current_row = start_row + 1;
      JXF_Int end_row;

      if (pre->A_i_block[num_rows] <= target_blocks)
      {
         pre->task_row_bounds[size++] = num_rows;
         break;
      }

      while (current_row < num_rows &&
             pre->A_i_block[current_row] < target_blocks)
      {
         current_row++;
      }

      if (current_row <= start_row + 1)
         end_row = start_row + 1;
      else
         end_row = current_row - 1;

      if (end_row > num_rows)
         end_row = num_rows;

      while (end_row > start_row + 1 && (pre->A_i_block[end_row] & 1))
      // while (end_row > start_row + 1 && (pre->A_i_block[end_row] & 7))
      {
         end_row--;
      }

      if (end_row <= start_row)
         end_row = start_row + 1;

      pre->task_row_bounds[size++] = end_row;
      start_row = end_row;
   }

   pre->num_tasks = size - 1;

   pre->dy_data =
       (JXF_Real *)hthread_malloc(cluster_id,
                                  sizeof(JXF_Real) * 2 * pre->num_blocks,
                                  HT_MEM_RW);

   if (pre->dy_data == NULL)
   {
      jxf_SpMVPrecondFP32DataDestroy(pre);
      return NULL;
   }

   // if (myid == 0)
   // {
   //    double util2 = (double)A_nnz / (2.0 * (double)pre->num_blocks);
   //    double A_bytes_per_nnz = (double)(sizeof(JXF_Real) * 2 * pre->num_blocks) / (double)A_nnz;
   //    double idx_bytes_per_nnz = (double)(sizeof(JXF_Int) * pre->num_blocks) / (double)A_nnz;
   //    double x_bytes_per_nnz = (double)(sizeof(JXF_Real) * 2 * pre->num_blocks) / (double)A_nnz;
   //    double dy_bytes_per_nnz = (double)(sizeof(JXF_Real) * 2 * pre->num_blocks) / (double)A_nnz;

   //    jxf_printf("[FP32-BLOCK2] rows = %d, nnz = %d, blocks = %d, util = %.4f\n",
   //               num_rows, A_nnz, pre->num_blocks, util2);

   //    jxf_printf("[FP32-BLOCK2] bytes/nnz: A = %.3f, idx = %.3f, x = %.3f, dy = %.3f, total = %.3f\n",
   //               A_bytes_per_nnz,
   //               idx_bytes_per_nnz,
   //               x_bytes_per_nnz,
   //               dy_bytes_per_nnz,
   //               A_bytes_per_nnz + idx_bytes_per_nnz + x_bytes_per_nnz + dy_bytes_per_nnz);
   // }

   return pre;
}

JXF_Int
jxf_CSRMatrixSpMVPrecondFP32Create(jxf_CSRMatrix *A,
                                   JXF_Int x_size,
                                   JXF_Int myid)
{
   if (A == NULL)
      return 1;

   if (jxf_CSRMatrixSpMVPrecondFP32(A) != NULL)
      return 0;

   jxf_CSRMatrixSpMVPrecondFP32(A) =
       jxf_SpMVPrecondFP32DataCreate(A, x_size, myid);

   if (jxf_CSRMatrixSpMVPrecondFP32(A) == NULL)
      return 1;

   return 0;
}

void jxf_CSRMatrixSpMVPrecondFP32Destroy(jxf_CSRMatrix *A)
{
   if (A == NULL)
      return;

   if (jxf_CSRMatrixSpMVPrecondFP32(A) == NULL)
      return;

   jxf_SpMVPrecondFP32DataDestroy(jxf_CSRMatrixSpMVPrecondFP32(A));
   jxf_CSRMatrixSpMVPrecondFP32(A) = NULL;
}


// =============================================
JXF_Int
jxf_CSRMatrixMatvec(JXF_Real alpha,
                    jxf_CSRMatrix *A,
                    jxf_Vector *x,
                    JXF_Real beta,
                    jxf_Vector *y,
                    JXF_Int myid)
{
   // 根据全局设置的 jxf_spmv_type  来选择具体的实现
   switch (jxf_spmv_type)
   {
   case 0:
      return jxf_CSRMatrixMatvec_origin(alpha, A, x, beta, y);
   case 1: // GSM 版本
      return jxf_CSRMatrixMatvec_v1(alpha, A, x, beta, y, myid);
   case 2: // x_new 版本
      return jxf_CSRMatrixMatvec_v2(alpha, A, x, beta, y, myid);
   default:
      if (myid == 0)
         fprintf(stderr, "Warning: Unknown jxf_spmv_type  (%d). Falling back to the origin version.\n", jxf_spmv_type);
      return jxf_CSRMatrixMatvec_origin(alpha, A, x, beta, y);
   }
}

/*!
 * \fn JXF_Int jxf_CSRMatrixMatvec
 * \brief Perform y = alpha*A*x + beta*y.
 * \date 2011/09/03
 */
JXF_Int
jxf_CSRMatrixMatvec_origin(JXF_Real alpha,
                           jxf_CSRMatrix *A,
                           jxf_Vector *x,
                           JXF_Real beta,
                           jxf_Vector *y)
{
   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A);
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);

   JXF_Int *A_rownnz = jxf_CSRMatrixRownnz(A);
   JXF_Int num_rownnz = jxf_CSRMatrixNumRownnz(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);
   JXF_Int idxstride_y = jxf_VectorIndexStride(y);
   JXF_Int vecstride_y = jxf_VectorVectorStride(y);
   JXF_Int idxstride_x = jxf_VectorIndexStride(x);
   JXF_Int vecstride_x = jxf_VectorVectorStride(x);

   JXF_Real temp, tempx;

   JXF_Int i, j, jj;

   JXF_Int m;

   JXF_Real xpar = 0.7;

   JXF_Int ierr = 0;

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  Matvec returns ierr = 1 if
    *  length of X doesn't equal the number of columns of A,
    *  ierr = 2 if the length of Y doesn't equal the number of rows
    *  of A, and ierr = 3 if both are true.
    *
    *  Because temporary vectors are often used in Matvec, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_cols != x_size)
      ierr = 1;

   if (num_rows != y_size)
      ierr = 2;

   if (num_cols != x_size && num_rows != y_size)
      ierr = 3;

   /*-----------------------------------------------------------------------
    * Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
    *-----------------------------------------------------------------------*/

   if (alpha == 0.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;

      return ierr;
   }

   /*--------------------------------------------
    * y = (beta/alpha)*y
    *------------------------------------------*/

   temp = beta / alpha;

   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] = 0.0;
      }
      else
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] *= temp;
      }
   }

   /*-----------------------------------------------
    * y += A*x
    *---------------------------------------------*/

   /* use rownnz pointer to do the A*x multiplication
      when num_rownnz is smaller than num_rows */

   if (num_rownnz < xpar * (num_rows))
   {

#define JXF_SMP_PRIVATE i, jj, m, tempx
#include "../../../include/jxf_smp_forloop.h"

      for (i = 0; i < num_rownnz; i++)
      {
         m = A_rownnz[i];

         if (num_vectors == 1)
         {
            tempx = y_data[m];
            for (jj = A_i[m]; jj < A_i[m + 1]; jj++)
            {
               tempx += A_data[jj] * x_data[A_j[jj]];
            }
            y_data[m] = tempx;
         }
         else
         {
            for (j = 0; j < num_vectors; ++j)
            {
               tempx = y_data[j * vecstride_y + m * idxstride_y];
               for (jj = A_i[m]; jj < A_i[m + 1]; jj++)
               {
                  tempx += A_data[jj] * x_data[j * vecstride_x + A_j[jj] * idxstride_x];
               }
               y_data[j * vecstride_y + m * idxstride_y] = tempx;
            }
         }
      }
   }
   else
   {

#define JXF_SMP_PRIVATE i, jj, temp
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows; i++)
      {
         if (num_vectors == 1)
         {
            temp = y_data[i];
            for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
            {
               temp += A_data[jj] * x_data[A_j[jj]];
            }
            y_data[i] = temp;
         }
         else
         {
            for (j = 0; j < num_vectors; ++j)
            {
               temp = y_data[j * vecstride_y + i * idxstride_y];
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  temp += A_data[jj] * x_data[j * vecstride_x + A_j[jj] * idxstride_x];
               }
               y_data[j * vecstride_y + i * idxstride_y] = temp;
            }
         }
      }
   }

   /*----------------------------------------------------
    *  y = alpha*y
    *---------------------------------------------------*/

   if (alpha != 1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}

/* y[offset:end] = alpha*A[offset:end,:]*x + beta*b[offset:end] */
JXF_Int
jxf_CSRMatrixMatvecOutOfPlace(JXF_Real alpha,
                              jxf_CSRMatrix *A,
                              jxf_Vector *x,
                              JXF_Real beta,
                              jxf_Vector *b,
                              jxf_Vector *y,
                              JXF_Int offset)
{
   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A) + offset;
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A) - offset;
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   /*JXF_Int         num_nnz  = jxf_CSRMatrixNumNonzeros(A);*/

   JXF_Int *A_rownnz = jxf_CSRMatrixRownnz(A);
   JXF_Int num_rownnz = jxf_CSRMatrixNumRownnz(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *b_data = jxf_VectorData(b) + offset;
   JXF_Real *y_data = jxf_VectorData(y) + offset;
   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int b_size = jxf_VectorSize(b) - offset;
   JXF_Int y_size = jxf_VectorSize(y) - offset;
   JXF_Int num_vectors = jxf_VectorNumVectors(x);
   JXF_Int idxstride_y = jxf_VectorIndexStride(y);
   JXF_Int vecstride_y = jxf_VectorVectorStride(y);
   /*JXF_Int         idxstride_b = jxf_VectorIndexStride(b);
   JXF_Int         vecstride_b = jxf_VectorVectorStride(b);*/
   JXF_Int idxstride_x = jxf_VectorIndexStride(x);
   JXF_Int vecstride_x = jxf_VectorVectorStride(x);
   JXF_Real temp, tempx;
   JXF_Int i, j, jj, m, ierr = 0;
   JXF_Real xpar = 0.7;
   jxf_Vector *x_tmp = NULL;

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  Matvec returns ierr = 1 if
    *  length of X doesn't equal the number of columns of A,
    *  ierr = 2 if the length of Y doesn't equal the number of rows
    *  of A, and ierr = 3 if both are true.
    *
    *  Because temporary vectors are often used in Matvec, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));
   jxf_assert(num_vectors == jxf_VectorNumVectors(b));

   if (num_cols != x_size)
      ierr = 1;

   if (num_rows != y_size || num_rows != b_size)
      ierr = 2;

   if (num_cols != x_size && (num_rows != y_size || num_rows != b_size))
      ierr = 3;

   /*-----------------------------------------------------------------------
    * Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
    *-----------------------------------------------------------------------*/

   if (alpha == 0.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] = beta * b_data[i];

      return ierr;
   }

   if (x == y)
   {
      x_tmp = jxf_SeqVectorCloneDeep(x);
      x_data = jxf_VectorData(x_tmp);
   }

   /*-----------------------------------------------------------------------
    * y = (beta/alpha)*y
    *-----------------------------------------------------------------------*/

   temp = beta / alpha;

   /* use rownnz pointer to do the A*x multiplication  when num_rownnz is smaller than num_rows */

   if (num_rownnz < xpar * (num_rows) || num_vectors > 1)
   {
      /*-----------------------------------------------------------------------
       * y = (beta/alpha)*y
       *-----------------------------------------------------------------------*/

      if (temp != 1.0)
      {
         if (temp == 0.0)
         {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < num_rows * num_vectors; i++)
               y_data[i] = 0.0;
         }
         else
         {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < num_rows * num_vectors; i++)
               y_data[i] = b_data[i] * temp;
         }
      }
      else
      {
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] = b_data[i];
      }

      /*-----------------------------------------------------------------
       * y += A*x
       *-----------------------------------------------------------------*/

      if (num_rownnz < xpar * (num_rows))
      {
#define JXF_SMP_PRIVATE i, j, jj, m, tempx
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_rownnz; i++)
         {
            m = A_rownnz[i];

            if (num_vectors == 1)
            {
               tempx = 0;
               for (jj = A_i[m]; jj < A_i[m + 1]; jj++)
                  tempx += A_data[jj] * x_data[A_j[jj]];
               y_data[m] += tempx;
            }
            else
               for (j = 0; j < num_vectors; ++j)
               {
                  tempx = 0;
                  for (jj = A_i[m]; jj < A_i[m + 1]; jj++)
                     tempx += A_data[jj] * x_data[j * vecstride_x + A_j[jj] * idxstride_x];
                  y_data[j * vecstride_y + m * idxstride_y] += tempx;
               }
         }
      }
      else // num_vectors > 1
      {
#define JXF_SMP_PRIVATE i, j, jj, tempx
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_rows; i++)
         {
            for (j = 0; j < num_vectors; ++j)
            {
               tempx = 0;
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  tempx += A_data[jj] * x_data[j * vecstride_x + A_j[jj] * idxstride_x];
               }
               y_data[j * vecstride_y + i * idxstride_y] += tempx;
            }
         }
      }

      /*-----------------------------------------------------------------
       * y = alpha*y
       *-----------------------------------------------------------------*/

      if (alpha != 1.0)
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] *= alpha;
      }
   }
   else
   { // JSP: this is currently the only path optimized
#define JXF_SMP_PRIVATE i, jj, tempx
#define JXF_SMP_PAR_REGION
#include "../../../include/jxf_smp_forloop.h"
      {
         JXF_Int iBegin = jxf_CSRMatrixGetLoadBalancedPartitionBegin(A);
         JXF_Int iEnd = jxf_CSRMatrixGetLoadBalancedPartitionEnd(A);
         jxf_assert(iBegin <= iEnd);
         jxf_assert(iBegin >= 0 && iBegin <= num_rows);
         jxf_assert(iEnd >= 0 && iEnd <= num_rows);

         if (0 == temp)
         {
            if (1 == alpha) // JSP: a common path
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = 0.0;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = A*x
            else if (-1 == alpha)
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = 0.0;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx -= A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = -A*x
            else
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = 0.0;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = alpha * tempx;
               }
            } // y = alpha*A*x
         } // temp == 0
         else if (-1 == temp) // beta == -alpha
         {
            if (1 == alpha) // JSP: a common path
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = -b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = A*x - y
            else if (-1 == alpha) // JSP: a common path
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx -= A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = -A*x + y
            else
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = -b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = alpha * tempx;
               }
            } // y = alpha*(A*x - y)
         } // temp == -1
         else if (1 == temp)
         {
            if (1 == alpha) // JSP: a common path
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = A*x + y
            else if (-1 == alpha)
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = -b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx -= A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = -A*x - y
            else
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = b_data[i];
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = alpha * tempx;
               }
            } // y = alpha*(A*x + y)
         }
         else
         {
            if (1 == alpha) // JSP: a common path
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = b_data[i] * temp;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = A*x + temp*y
            else if (-1 == alpha)
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = -b_data[i] * temp;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx -= A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = tempx;
               }
            } // y = -A*x - temp*y
            else
            {
               for (i = iBegin; i < iEnd; i++)
               {
                  tempx = b_data[i] * temp;
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     tempx += A_data[jj] * x_data[A_j[jj]];
                  }
                  y_data[i] = alpha * tempx;
               }
            } // y = alpha*(A*x + temp*y)
         } // temp != 0 && temp != -1 && temp != 1
      } // omp parallel
   }

   if (x == y)
   {
      jxf_SeqVectorDestroy(x_tmp);
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxf_CSRMatrixMatvecT
 * \brief Perform y = alpha*A^T*x + beta*y.
 * \date 2011/09/03
 */
JXF_Int
jxf_CSRMatrixMatvecT(JXF_Real alpha,
                     jxf_CSRMatrix *A,
                     jxf_Vector *x,
                     JXF_Real beta,
                     jxf_Vector *y)
{
   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A);
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);
   JXF_Int idxstride_y = jxf_VectorIndexStride(y);
   JXF_Int vecstride_y = jxf_VectorVectorStride(y);
   JXF_Int idxstride_x = jxf_VectorIndexStride(x);
   JXF_Int vecstride_x = jxf_VectorVectorStride(x);

   JXF_Real temp;

   JXF_Int i, i1, j, jv, jj, ns, ne, size, rest;
   JXF_Int num_threads;

   JXF_Int ierr = 0;

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  MatvecT returns ierr = 1 if
    *  length of X doesn't equal the number of rows of A,
    *  ierr = 2 if the length of Y doesn't equal the number of
    *  columns of A, and ierr = 3 if both are true.
    *
    *  Because temporary vectors are often used in MatvecT, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_rows != x_size)
      ierr = 1;

   if (num_cols != y_size)
      ierr = 2;

   if (num_rows != x_size && num_cols != y_size)
      ierr = 3;

   /*--------------------------------------------------------------
    * Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
    *-------------------------------------------------------------*/

   if (alpha == 0.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= beta;
      }

      return ierr;
   }

   /*------------------------------------------------
    * y = (beta/alpha)*y
    *-----------------------------------------------*/

   temp = beta / alpha;

   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] = 0.0;
         }
      }
      else
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] *= temp;
         }
      }
   }

   /*-----------------------------------------------------------------
    * y += A^T*x
    *-----------------------------------------------------------------*/
   num_threads = jxf_NumThreads();
   if (num_threads > 1)
   {

#define JXF_SMP_PRIVATE i, i1, jj, j, ns, ne, size, rest
#include "../../../include/jxf_smp_forloop.h"
      for (i1 = 0; i1 < num_threads; i1++)
      {
         size = num_cols / num_threads;
         rest = num_cols - size * num_threads;
         if (i1 < rest)
         {
            ns = i1 * size + i1 - 1;
            ne = (i1 + 1) * size + i1 + 1;
         }
         else
         {
            ns = i1 * size + rest - 1;
            ne = (i1 + 1) * size + rest;
         }

         if (num_vectors == 1)
         {
            for (i = 0; i < num_rows; i++)
            {
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];
                  if (j > ns && j < ne)
                  {
                     y_data[j] += A_data[jj] * x_data[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < num_rows; i++)
            {
               for (jv = 0; jv < num_vectors; ++jv)
               {
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     j = A_j[jj];
                     if (j > ns && j < ne)
                     {
                        y_data[j * idxstride_y + jv * vecstride_y] +=
                            A_data[jj] * x_data[i * idxstride_x + jv * vecstride_x];
                     }
                  }
               }
            }
         }
      }
   }
   else
   {
      for (i = 0; i < num_rows; i++)
      {
         if (num_vectors == 1)
         {
            for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
            {
               j = A_j[jj];
               y_data[j] += A_data[jj] * x_data[i];
            }
         }
         else
         {
            for (jv = 0; jv < num_vectors; ++jv)
            {
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];
                  y_data[j * idxstride_y + jv * vecstride_y] +=
                      A_data[jj] * x_data[i * idxstride_x + jv * vecstride_x];
               }
            }
         }
      }
   }

   /*--------------------------------------------
    * y = alpha*y
    *-------------------------------------------*/

   if (alpha != 1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxh_CSRMatrixMatvecT
 * \brief Perform y = alpha*(scale*A_H^T*x with F-mul and F-acc) + beta*y.
 *
 * A is stored in FP16.
 * Multiplication and accumulation are performed in JXF_Real.
 *
 * Numerical model:
 *     aij_fp32 = float(A_H_ij)
 *     y_j     += scale * aij_fp32 * x_i
 */
JXF_Int
jxh_CSRMatrixMatvecT(JXF_Real alpha,
                     jxh_CSRMatrix *A,
                     jxf_Vector *x,
                     JXF_Real beta,
                     jxf_Vector *y)
{
   // static JXF_Int print_count_T_HmulFacc = 0;
   // if (print_count_T_HmulFacc < 2)
   // {
   //    jxf_printf("[DEBUG] enter jxh_CSRMatrixMatvecT\n");
   //    print_count_T_HmulFacc++;
   // }

   JXH_Real *A_data = jxh_CSRMatrixData(A);
   JXF_Int *A_i = jxh_CSRMatrixI(A);
   JXF_Int *A_j = jxh_CSRMatrixJ(A);

   JXF_Int num_rows = jxh_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxh_CSRMatrixNumCols(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);

   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);

   JXF_Int num_vectors = jxf_VectorNumVectors(x);

   JXF_Int idxstride_y = jxf_VectorIndexStride(y);
   JXF_Int vecstride_y = jxf_VectorVectorStride(y);
   JXF_Int idxstride_x = jxf_VectorIndexStride(x);
   JXF_Int vecstride_x = jxf_VectorVectorStride(x);

   JXF_Real temp;
   JXF_Real scale;
   JXF_Real aij_fp32;

   JXF_Int i, i1, j, jv, jj;
   JXF_Int ns, ne, size, rest;
   JXF_Int num_threads;

   JXF_Int ierr = 0;

   /*---------------------------------------------------------------------
    * Check for size compatibility.
    *
    * MatvecT computes:
    *
    *     y = alpha * A^T * x + beta * y
    *
    * Therefore:
    *     size(x) should equal num_rows(A),
    *     size(y) should equal num_cols(A).
    *--------------------------------------------------------------------*/

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_rows != x_size)
   {
      ierr = 1;
   }

   if (num_cols != y_size)
   {
      ierr = 2;
   }

   if (num_rows != x_size && num_cols != y_size)
   {
      ierr = 3;
   }

   /*
    * Recover scaling factor for FP16-stored matrix.
    *
    * The represented matrix is approximately:
    *
    *     A_F ~= scale * A_H
    */
   scale = jxh_CSRMatrixScale(A);

   if (scale <= 0.0f)
   {
      scale = 1.0f;
   }

   /*--------------------------------------------------------------
    * alpha == 0:
    *
    * y = beta * y
    *-------------------------------------------------------------*/

   if (alpha == 0.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= beta;
      }

      return ierr;
   }

   /*------------------------------------------------
    * Same algebraic structure as jxf_CSRMatrixMatvecT:
    *
    *     y = (beta / alpha) * y
    *     y += scale * A_H^T * x
    *     y = alpha * y
    *
    * Final result:
    *
    *     y = beta*y_old + alpha*scale*A_H^T*x
    *-----------------------------------------------*/

   temp = beta / alpha;

   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] = 0.0;
         }
      }
      else
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] *= temp;
         }
      }
   }

   /*-----------------------------------------------------------------
    * y += scale * A_H^T * x
    *
    * FP16 storage + FP32 multiplication + FP32 accumulation:
    *
    *     aij_fp32 = float(A_H_ij)
    *     y_j     += scale * aij_fp32 * x_i
    *-----------------------------------------------------------------*/

   num_threads = jxf_NumThreads();

   if (num_threads > 1)
   {
#define JXF_SMP_PRIVATE i, i1, jj, j, jv, ns, ne, size, rest, aij_fp32
#include "../../../include/jxf_smp_forloop.h"
      for (i1 = 0; i1 < num_threads; i1++)
      {
         size = num_cols / num_threads;
         rest = num_cols - size * num_threads;

         if (i1 < rest)
         {
            ns = i1 * size + i1 - 1;
            ne = (i1 + 1) * size + i1 + 1;
         }
         else
         {
            ns = i1 * size + rest - 1;
            ne = (i1 + 1) * size + rest;
         }

         if (num_vectors == 1)
         {
            for (i = 0; i < num_rows; i++)
            {
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];

                  if (j > ns && j < ne)
                  {
                     aij_fp32 = (JXF_Real)A_data[jj];

                     y_data[j] += scale * aij_fp32 * x_data[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < num_rows; i++)
            {
               for (jv = 0; jv < num_vectors; ++jv)
               {
                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     j = A_j[jj];

                     if (j > ns && j < ne)
                     {
                        aij_fp32 = (JXF_Real)A_data[jj];

                        y_data[j * idxstride_y + jv * vecstride_y] +=
                            scale * aij_fp32 *
                            x_data[i * idxstride_x + jv * vecstride_x];
                     }
                  }
               }
            }
         }
      }
   }
   else
   {
      if (num_vectors == 1)
      {
         for (i = 0; i < num_rows; i++)
         {
            for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
            {
               j = A_j[jj];

               aij_fp32 = (JXF_Real)A_data[jj];

               y_data[j] += scale * aij_fp32 * x_data[i];
            }
         }
      }
      else
      {
         for (i = 0; i < num_rows; i++)
         {
            for (jv = 0; jv < num_vectors; ++jv)
            {
               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];

                  aij_fp32 = (JXF_Real)A_data[jj];

                  y_data[j * idxstride_y + jv * vecstride_y] +=
                      scale * aij_fp32 *
                      x_data[i * idxstride_x + jv * vecstride_x];
               }
            }
         }
      }
   }

   /*--------------------------------------------
    * y = alpha * y
    *
    * Final result:
    *
    *     y = beta*y_old + alpha*scale*A_H^T*x
    *-------------------------------------------*/

   if (alpha != 1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxh_CSRMatrixMatvecT_HmulFacc
 * \brief Perform y = alpha*(scale*A_H^T*x_H with H-mul and F-acc) + beta*y.
 *
 * A is stored in half precision.
 * x/y are stored in JXF_Real.
 * Each product A_H[jj] * x[i] is rounded to JXH_Real first,
 * then accumulated in JXF_Real.
 *
 * Numerical model:
 *     x_h    = half(x_i)
 *     prod_h = half(A_H_ij * x_h)
 *     y_j   += scale * float(prod_h)
 */
JXF_Int
jxh_CSRMatrixMatvecT_HmulFacc(JXF_Real alpha,
                              jxh_CSRMatrix *A,
                              jxf_Vector *x,
                              JXF_Real beta,
                              jxf_Vector *y)
{
   // static JXF_Int print_count_T_HmulFacc = 0;
   // if (print_count_T_HmulFacc < 2)
   // {
   //    jxf_printf("[DEBUG] enter jxh_CSRMatrixMatvecT_HmulFacc\n");
   //    print_count_T_HmulFacc++;
   // }

   JXH_Real *A_data = jxh_CSRMatrixData(A);
   JXF_Int *A_i = jxh_CSRMatrixI(A);
   JXF_Int *A_j = jxh_CSRMatrixJ(A);

   JXF_Int num_rows = jxh_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxh_CSRMatrixNumCols(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);

   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);

   JXF_Int num_vectors = jxf_VectorNumVectors(x);

   JXF_Int idxstride_y = jxf_VectorIndexStride(y);
   JXF_Int vecstride_y = jxf_VectorVectorStride(y);
   JXF_Int idxstride_x = jxf_VectorIndexStride(x);
   JXF_Int vecstride_x = jxf_VectorVectorStride(x);

   JXF_Real temp;
   JXF_Real scale;

   /*
    * Half variables used to emulate FP16 multiplication.
    */
   JXH_Real x_h;
   JXH_Real prod_h;

   JXF_Int i, i1, j, jv, jj;
   JXF_Int ns, ne, size, rest;
   JXF_Int num_threads;

   JXF_Int ierr = 0;

   /*---------------------------------------------------------------------
    * Check for size compatibility.
    *
    * MatvecT computes:
    *
    *     y = alpha * A^T * x + beta * y
    *
    * Therefore:
    *     size(x) should equal num_rows(A),
    *     size(y) should equal num_cols(A).
    *--------------------------------------------------------------------*/

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_rows != x_size)
   {
      ierr = 1;
   }

   if (num_cols != y_size)
   {
      ierr = 2;
   }

   if (num_rows != x_size && num_cols != y_size)
   {
      ierr = 3;
   }

   /*
    * Recover scaling factor for FP16-stored matrix.
    *
    * The represented matrix is approximately:
    *
    *     A_F ~= scale * A_H
    */
   scale = jxh_CSRMatrixScale(A);

   if (scale <= 0.0f)
   {
      scale = 1.0f;
   }

   /*--------------------------------------------------------------
    * alpha == 0:
    *
    * y = beta * y
    *-------------------------------------------------------------*/

   if (alpha == 0.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= beta;
      }

      return ierr;
   }

   /*------------------------------------------------
    * Same algebraic structure as jxf_CSRMatrixMatvecT:
    *
    *     y = (beta / alpha) * y
    *     y += scale * A_H^T * x
    *     y = alpha * y
    *
    * Final result:
    *
    *     y = beta*y_old + alpha*scale*Hmul(A_H^T, x_H)
    *-----------------------------------------------*/

   temp = beta / alpha;

   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] = 0.0;
         }
      }
      else
      {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] *= temp;
         }
      }
   }

   /*-----------------------------------------------------------------
    * y += scale * A_H^T * x
    *
    * HmulFacc:
    *
    *     x_h    = half(x_i)
    *     prod_h = half(A_H_ij * x_h)
    *     y_j   += scale * float(prod_h)
    *
    * Accumulation remains FP32.
    *-----------------------------------------------------------------*/

   num_threads = jxf_NumThreads();

   if (num_threads > 1)
   {
#define JXF_SMP_PRIVATE i, i1, jj, j, jv, ns, ne, size, rest, x_h, prod_h
#include "../../../include/jxf_smp_forloop.h"
      for (i1 = 0; i1 < num_threads; i1++)
      {
         size = num_cols / num_threads;
         rest = num_cols - size * num_threads;

         if (i1 < rest)
         {
            ns = i1 * size + i1 - 1;
            ne = (i1 + 1) * size + i1 + 1;
         }
         else
         {
            ns = i1 * size + rest - 1;
            ne = (i1 + 1) * size + rest;
         }

         if (num_vectors == 1)
         {
            for (i = 0; i < num_rows; i++)
            {
               /*
                * x_i is rounded to FP16 before multiplication.
                */
               x_h = (JXH_Real)x_data[i];

               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];

                  if (j > ns && j < ne)
                  {
                     /*
                      * A_data[jj] is already FP16.
                      *
                      * prod_h is explicitly rounded to FP16,
                      * then converted back to FP32 for accumulation.
                      */
                     prod_h = (JXH_Real)(A_data[jj] * x_h);

                     y_data[j] += scale * ((JXF_Real)prod_h);
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < num_rows; i++)
            {
               for (jv = 0; jv < num_vectors; ++jv)
               {
                  x_h = (JXH_Real)x_data[i * idxstride_x +
                                         jv * vecstride_x];

                  for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
                  {
                     j = A_j[jj];

                     if (j > ns && j < ne)
                     {
                        prod_h = (JXH_Real)(A_data[jj] * x_h);

                        y_data[j * idxstride_y + jv * vecstride_y] +=
                            scale * ((JXF_Real)prod_h);
                     }
                  }
               }
            }
         }
      }
   }
   else
   {
      if (num_vectors == 1)
      {
         for (i = 0; i < num_rows; i++)
         {
            /*
             * x_i is rounded to FP16 before multiplication.
             */
            x_h = (JXH_Real)x_data[i];

            for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
            {
               j = A_j[jj];

               /*
                * FP16 multiplication model:
                *
                *     prod_h = half(A_H_ij * half(x_i))
                *
                * FP32 accumulation:
                *
                *     y_j += scale * float(prod_h)
                */
               prod_h = (JXH_Real)(A_data[jj] * x_h);

               y_data[j] += scale * ((JXF_Real)prod_h);
            }
         }
      }
      else
      {
         for (i = 0; i < num_rows; i++)
         {
            for (jv = 0; jv < num_vectors; ++jv)
            {
               x_h = (JXH_Real)x_data[i * idxstride_x +
                                      jv * vecstride_x];

               for (jj = A_i[i]; jj < A_i[i + 1]; jj++)
               {
                  j = A_j[jj];

                  prod_h = (JXH_Real)(A_data[jj] * x_h);

                  y_data[j * idxstride_y + jv * vecstride_y] +=
                      scale * ((JXF_Real)prod_h);
               }
            }
         }
      }
   }

   /*--------------------------------------------
    * y = alpha * y
    *
    * Keep final scaling in FP32.
    *-------------------------------------------*/

   if (alpha != 1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}

// =======================================================
// 异构版本 ===============================================
JXF_Int
jxf_CSRMatrixMatvec_baseline_FP32(JXF_Real alpha,
                                  jxf_CSRMatrix *A,
                                  jxf_Vector *x,
                                  JXF_Real beta,
                                  jxf_Vector *y,
                                  JXF_Int myid)
{
   if (myid == 0)
      printf("\n jxf_CSRMatrixMatvec_baseline_FP32 (GSM) \n");

   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A);
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   JXF_Int A_nnz = jxf_CSRMatrixNumNonzeros(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);

   JXF_Int i;
   JXF_Int ierr = 0;
   JXF_Int cluster_id = myid % 4;

   JXF_Int *core_row_bounds = NULL;
   JXF_Int *dA_j = NULL;
   JXF_Real *x_ext = NULL;
   JXF_Real *dy_data = NULL;
   JXF_Int threads_id = -1;
   JXF_Int barrier_id = -1;

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_cols != x_size)
      ierr = 1;

   if (num_rows != y_size)
      ierr = 2;

   if (num_cols != x_size && num_rows != y_size)
      ierr = 3;

   if (alpha == 0.0f)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;

      return ierr;
   }

   core_row_bounds = hthread_malloc(cluster_id, sizeof(JXF_Int) * (coreNums + 1), HT_MEM_RW);
   dA_j = hthread_malloc(cluster_id, sizeof(JXF_Int) * A_nnz, HT_MEM_RW);

   JXF_Int x_ext_size = 2 * x_size;
   x_ext = hthread_malloc(cluster_id, sizeof(JXF_Real) * x_ext_size, HT_MEM_RW);
   dy_data = hthread_malloc(cluster_id, sizeof(JXF_Real) * num_rows, HT_MEM_RW);

   JXF_Int nnz_per_thread = A_nnz / coreNums;
   core_row_bounds[0] = 0;

   for (JXF_Int tid = 0; tid < coreNums; tid++)
   {
      JXF_Int start_row = core_row_bounds[tid];

      if (start_row >= num_rows)
      {
         core_row_bounds[tid + 1] = num_rows;
         continue;
      }

      if (tid == coreNums - 1)
      {
         core_row_bounds[tid + 1] = num_rows;
         break;
      }

      JXF_Int target_nnz = A_i[start_row] + nnz_per_thread;
      JXF_Int low = start_row + 1;
      JXF_Int high = num_rows;
      JXF_Int boundary = start_row + 1;

      while (low <= high)
      {
         JXF_Int mid = low + (high - low) / 2;

         if (A_i[mid] <= target_nnz)
         {
            boundary = mid;
            low = mid + 1;
         }
         else
         {
            high = mid - 1;
         }
      }

      if (boundary <= start_row)
         boundary = start_row + 1;

      if (boundary > num_rows)
         boundary = num_rows;

      core_row_bounds[tid + 1] = boundary;
   }

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < A_nnz; i++)
   {
      dA_j[i] = A_j[i] * 8;
   }

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < x_size; i++)
   {
      x_ext[2 * i] = x_data[i];
      x_ext[2 * i + 1] = 0.0f;
   }

   threads_id = hthread_group_create(cluster_id, coreNums);
   barrier_id = hthread_barrier_create(cluster_id);

   unsigned long args[10];
   args[0] = (unsigned long)num_rows;
   args[1] = (unsigned long)x_ext_size;
   args[2] = (unsigned long)coreNums;
   args[3] = (unsigned long)barrier_id;
   args[4] = (unsigned long)A_i;
   args[5] = (unsigned long)dA_j;
   args[6] = (unsigned long)A_data;
   args[7] = (unsigned long)x_ext;
   args[8] = (unsigned long)dy_data;
   args[9] = (unsigned long)core_row_bounds;

   hthread_group_exec(threads_id, "spmv_baseline_FP32", 4, 6, args);
   hthread_group_wait(threads_id);

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < num_rows; i++)
   {
      y_data[i] = alpha * dy_data[i] + beta * y_data[i];
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);
   hthread_free(core_row_bounds);
   hthread_free(dA_j);
   hthread_free(x_ext);
   hthread_free(dy_data);

   return ierr;
}

// JXF_Int
// jxf_CSRMatrixMatvec_v1(JXF_Real alpha,
//                        jxf_CSRMatrix *A,
//                        jxf_Vector *x,
//                        JXF_Real beta,
//                        jxf_Vector *y,
//                        JXF_Int myid)
// {
//    static JXF_Int print_count_timing = 0;

//    JXF_Real t0, t1, t2, t3, t4, t5;
//    JXF_Real xcopy_time, create_time, kernel_time;
//    JXF_Real reduce_time, cleanup_time, total_time;

//    JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
//    JXF_Int num_cols = jxf_CSRMatrixNumCols(A);

//    JXF_Real *x_data = jxf_VectorData(x);
//    JXF_Real *y_data = jxf_VectorData(y);

//    JXF_Int x_size = jxf_VectorSize(x);
//    JXF_Int y_size = jxf_VectorSize(y);
//    JXF_Int num_vectors = jxf_VectorNumVectors(x);

//    JXF_Int i;
//    JXF_Int ierr = 0;
//    JXF_Int cluster_id = myid % 4;

//    jxf_SpMVPrecondFP32Data *pre = NULL;

//    jxf_assert(num_vectors == jxf_VectorNumVectors(y));

//    if (num_cols != x_size)
//       ierr = 1;

//    if (num_rows != y_size)
//       ierr = 2;

//    if (num_cols != x_size && num_rows != y_size)
//       ierr = 3;

//    if (alpha == 0.0f)
//    {
// #define JXF_SMP_PRIVATE i
// #include "../../../include/jxf_smp_forloop.h"
//       for (i = 0; i < num_rows * num_vectors; i++)
//       {
//          y_data[i] *= beta;
//       }

//       return ierr;
//    }

//    pre = jxf_CSRMatrixSpMVPrecondFP32(A);

//    if (pre == NULL)
//    {
//       ierr = jxf_CSRMatrixSpMVPrecondFP32Create(A, x_size, myid);
//       if (ierr)
//       {
//          if (myid == 0)
//          {
//             fprintf(stderr, "Fatal Error: jxf_CSRMatrixSpMVPrecondFP32Create failed.\n");
//          }
//          return ierr;
//       }

//       pre = jxf_CSRMatrixSpMVPrecondFP32(A);
//    }

//    if (pre == NULL)
//       return 1;

//    t0 = jxf_MPI_Wtime();

//    memcpy(pre->x_pair, x_data, sizeof(JXF_Real) * x_size);

//    if (x_size < pre->x_pair_size)
//    {
//       pre->x_pair[x_size] = 0.0f;
//    }

//    t1 = jxf_MPI_Wtime();

//    JXF_Int threads_id = hthread_group_create(cluster_id, coreNums);
//    JXF_Int barrier_id = hthread_barrier_create(cluster_id);

//    unsigned long args[11];

//    args[0] = (unsigned long)pre->num_rows;
//    args[1] = (unsigned long)pre->x_pair_size;
//    args[2] = (unsigned long)coreNums;
//    args[3] = (unsigned long)pre->num_tasks;
//    args[4] = (unsigned long)barrier_id;
//    args[5] = (unsigned long)pre->A_i_block;
//    args[6] = (unsigned long)pre->A_j_block;
//    args[7] = (unsigned long)pre->A_data_block;
//    args[8] = (unsigned long)pre->x_pair;
//    args[9] = (unsigned long)pre->dy_data;
//    args[10] = (unsigned long)pre->task_row_bounds;

//    t2 = jxf_MPI_Wtime();

//    hthread_group_exec(threads_id, "SpMV_GSM_FP32_BLOCK2", 5, 6, args);
//    hthread_group_wait(threads_id);

//    t3 = jxf_MPI_Wtime();

//    JXF_Real row_sum;
//    JXF_Int row_start, row_end;
//    JXF_Int p;

//    JXF_Int *Ai_block = pre->A_i_block;
//    JXF_Real *dy = pre->dy_data;

// #define JXF_SMP_PRIVATE i, p, row_start, row_end, row_sum
// #include "../../../include/jxf_smp_forloop.h"
//    for (i = 0; i < pre->num_rows; i++)
//    {
//       row_sum = 0.0f;
//       row_start = 2 * Ai_block[i];
//       row_end = 2 * Ai_block[i + 1];

//       for (p = row_start; p < row_end; p += 2)
//       {
//          row_sum += dy[p] + dy[p + 1];
//       }

//       y_data[i] = beta * y_data[i] + alpha * row_sum;
//    }

//    t4 = jxf_MPI_Wtime();

//    hthread_barrier_destroy(barrier_id);
//    hthread_group_destroy(threads_id);

//    t5 = jxf_MPI_Wtime();

//    xcopy_time = t1 - t0;
//    create_time = t2 - t1;
//    kernel_time = t3 - t2;
//    reduce_time = t4 - t3;
//    cleanup_time = t5 - t4;
//    total_time = t5 - t0;

//    if (myid == 0 && print_count_timing < 12)
//    {
//       printf("[FP32 v1 timing] "
//              "xcopy=%e, create=%e, kernel=%e, reduce=%e, cleanup=%e, total=%e\n",
//              (double)xcopy_time,
//              (double)create_time,
//              (double)kernel_time,
//              (double)reduce_time,
//              (double)cleanup_time,
//              (double)total_time);

//       print_count_timing++;
//    }

//    return ierr;
// }

/* 干净版，无计时 */
JXF_Int
jxf_CSRMatrixMatvec_v1(JXF_Real alpha,
                       jxf_CSRMatrix *A,
                       jxf_Vector *x,
                       JXF_Real beta,
                       jxf_Vector *y,
                       JXF_Int myid)
{
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);

   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);

   JXF_Int i;
   JXF_Int ierr = 0;
   JXF_Int cluster_id = myid % 4;

   jxf_SpMVPrecondFP32Data *pre = NULL;

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_cols != x_size)
      ierr = 1;

   if (num_rows != y_size)
      ierr = 2;

   if (num_cols != x_size && num_rows != y_size)
      ierr = 3;

   if (alpha == 0.0f)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
      {
         y_data[i] *= beta;
      }

      return ierr;
   }

   pre = jxf_CSRMatrixSpMVPrecondFP32(A);

   if (pre == NULL)
   {
      ierr = jxf_CSRMatrixSpMVPrecondFP32Create(A, x_size, myid);
      if (ierr)
      {
         if (myid == 0)
         {
            fprintf(stderr, "Fatal Error: jxf_CSRMatrixSpMVPrecondFP32Create failed.\n");
         }
         return ierr;
      }

      pre = jxf_CSRMatrixSpMVPrecondFP32(A);
   }

   if (pre == NULL)
      return 1;

   memcpy(pre->x_pair, x_data, sizeof(JXF_Real) * x_size);

   if (x_size < pre->x_pair_size)
   {
      pre->x_pair[x_size] = 0.0f;
   }

   JXF_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JXF_Int barrier_id = hthread_barrier_create(cluster_id);

   unsigned long args[11];

   args[0] = (unsigned long)pre->num_rows;
   args[1] = (unsigned long)pre->x_pair_size;
   args[2] = (unsigned long)coreNums;
   args[3] = (unsigned long)pre->num_tasks;
   args[4] = (unsigned long)barrier_id;
   args[5] = (unsigned long)pre->A_i_block;
   args[6] = (unsigned long)pre->A_j_block;
   args[7] = (unsigned long)pre->A_data_block;
   args[8] = (unsigned long)pre->x_pair;
   args[9] = (unsigned long)pre->dy_data;
   args[10] = (unsigned long)pre->task_row_bounds;

   hthread_group_exec(threads_id, "SpMV_GSM_FP32_BLOCK2", 5, 6, args);
   hthread_group_wait(threads_id);

   JXF_Real row_sum;
   JXF_Int row_start, row_end;
   JXF_Int p;

   JXF_Int *Ai_block = pre->A_i_block;
   JXF_Real *dy = pre->dy_data;

#define JXF_SMP_PRIVATE i, p, row_start, row_end, row_sum
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < pre->num_rows; i++)
   {
      row_sum = 0.0f;
      row_start = 2 * Ai_block[i];
      row_end = 2 * Ai_block[i + 1];

      for (p = row_start; p < row_end; p += 2)
      {
         row_sum += dy[p] + dy[p + 1];
      }

      y_data[i] = beta * y_data[i] + alpha * row_sum;
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   return ierr;
}

JXF_Int
jxf_CSRMatrixMatvec_v2(JXF_Real alpha,
                       jxf_CSRMatrix *A,
                       jxf_Vector *x,
                       JXF_Real beta,
                       jxf_Vector *y,
                       JXF_Int myid)
{
   JXF_Real *A_data = jxf_CSRMatrixData(A);
   JXF_Int *A_i = jxf_CSRMatrixI(A);
   JXF_Int *A_j = jxf_CSRMatrixJ(A);
   JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
   JXF_Int num_cols = jxf_CSRMatrixNumCols(A);
   JXF_Int A_nnz = jxf_CSRMatrixNumNonzeros(A);

   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int x_size = jxf_VectorSize(x);
   JXF_Int y_size = jxf_VectorSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);

   JXF_Int i, j, jj;
   JXF_Int ierr = 0;
   JXF_Int cluster_id = myid % 4;

   jxf_assert(num_vectors == jxf_VectorNumVectors(y));

   if (num_cols != x_size)
      ierr = 1;
   if (num_rows != y_size)
      ierr = 2;
   if (num_cols != x_size && num_rows != y_size)
      ierr = 3;

   if (alpha == 0.0f)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;

      return ierr;
   }

   const JXF_Int ALIGN_SIZE = 2;
   JXF_Int AM_size = (96000 / 3) * 2;

   JXF_Int capacity = num_rows + 1;
   JXF_Int *row_bounds =
       (JXF_Int *)hthread_malloc(cluster_id,
                                 sizeof(JXF_Int) * capacity,
                                 HT_MEM_RW);

   JXF_Int size = 0;
   JXF_Int start_row = 0;

   row_bounds[size++] = 0;

   while (start_row < num_rows)
   {
      JXF_Int target_nnz = A_i[start_row] + AM_size;
      JXF_Int current_row = start_row + 1;
      JXF_Int end_row;

      if (A_i[num_rows] <= target_nnz)
      {
         row_bounds[size++] = num_rows;
         break;
      }

      while (current_row < num_rows && A_i[current_row] < target_nnz)
      {
         current_row++;
      }

      if (current_row <= start_row + 1)
         end_row = start_row + 1;
      else
         end_row = current_row - 1;

      if (end_row > num_rows)
         end_row = num_rows;

      while (end_row > start_row && (A_i[end_row] % ALIGN_SIZE != 0))
      {
         end_row--;
      }

      if (end_row <= start_row)
         end_row = start_row + 1;

      row_bounds[size++] = end_row;
      start_row = end_row;
   }

   JXF_Int num_tasks = size - 1;

   JXF_Int padded_nnz = A_nnz;
   if (padded_nnz % ALIGN_SIZE != 0)
      padded_nnz += ALIGN_SIZE - (padded_nnz % ALIGN_SIZE);

   JXF_Real *dx_data =
       (JXF_Real *)hthread_malloc(cluster_id,
                                  sizeof(JXF_Real) * padded_nnz,
                                  HT_MEM_RW);

#pragma omp parallel for private(jj) num_threads(1)
   for (jj = 0; jj < A_nnz; jj++)
   {
      dx_data[jj] = x_data[A_j[jj]];
   }

   for (jj = A_nnz; jj < padded_nnz; jj++)
   {
      dx_data[jj] = 0.0f;
   }

   JXF_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JXF_Int barrier_id = hthread_barrier_create(cluster_id);

   unsigned long args[6];

   args[0] = (unsigned long)coreNums;
   args[1] = (unsigned long)num_tasks;
   args[2] = (unsigned long)A_i;
   args[3] = (unsigned long)A_data;
   args[4] = (unsigned long)dx_data;
   args[5] = (unsigned long)row_bounds;

   hthread_group_exec(threads_id, "kernel_SpMV_FP32B", 2, 4, args);
   hthread_group_wait(threads_id);

   JXF_Real row_sum;

   if (beta == 0.0f)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] = 0.0f;
   }
   else
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;
   }

#define JXF_SMP_PRIVATE i, j, row_sum
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < num_rows; i++)
   {
      row_sum = 0.0f;

      for (j = A_i[i]; j < A_i[i + 1]; j++)
      {
         row_sum += dx_data[j];
      }

      y_data[i] += alpha * row_sum;
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   hthread_free(dx_data);
   hthread_free(row_bounds);

   return ierr;
}

