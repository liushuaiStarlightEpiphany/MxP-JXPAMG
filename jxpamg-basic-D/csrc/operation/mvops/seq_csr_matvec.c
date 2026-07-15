//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  seq_csr_matvec.c -- basic operations for mat-vec multiplication.
 *  Date: 2011/09/03
 */

#include "jx_mv.h"
#include "hthread_host.h"

// -------------------------------------
extern JX_Int myid;
extern JX_Int coreNums;
extern JX_Int jx_spmv_type;
// -------------------------------------

struct jx_SpMVPrecondFP64Data_struct
{
   JX_Int cluster_id;

   JX_Int num_rows;
   JX_Int num_cols;
   JX_Int A_nnz;

   JX_Int num_tasks;

   JX_Int *dA_j;
   JX_Int *task_row_bounds;
   JX_Real *dy_data;
};

/* Device-side preprocessing owned by a relax-103 color cache. */
struct jx_Relax103SpMVData_struct
{
   JX_Int cluster_id;
   JX_Int num_rows;
   JX_Int num_nonzeros;
   JX_Int num_colors;

   JX_Int *dA_j;
   JX_Real *dy_data;

   JX_Int *task_offsets;
   JX_Int *task_bounds;
};

static void
jx_SpMVPrecondFP64DataDestroy(jx_SpMVPrecondFP64Data *pre)
{
   if (pre == NULL)
      return;

   if (pre->dA_j != NULL)
      hthread_free(pre->dA_j);

   if (pre->dy_data != NULL)
      hthread_free(pre->dy_data);

   if (pre->task_row_bounds != NULL)
      hthread_free(pre->task_row_bounds);

   jx_TFree(pre);
}

static void
jx_Relax103SpMVDataDestroy(void *data_ptr)
{
   struct jx_Relax103SpMVData_struct *pre =
      (struct jx_Relax103SpMVData_struct *)data_ptr;

   if (pre == NULL)
   {
      return;
   }

   if (pre->dA_j != NULL)
   {
      hthread_free(pre->dA_j);
   }
   if (pre->dy_data != NULL)
   {
      hthread_free(pre->dy_data);
   }
   if (pre->task_bounds != NULL)
   {
      hthread_free(pre->task_bounds);
   }

   jx_TFree(pre->task_offsets);
   jx_TFree(pre);
}

static struct jx_Relax103SpMVData_struct *
jx_Relax103SpMVDataCreate(jx_Relax103ColorData *color_data,
                          JX_Int myid)
{
   struct jx_Relax103SpMVData_struct *pre;
   JX_Int cluster_id = myid % 4;
   JX_Int AM_size = 96000 / 3;
   JX_Int bounds_capacity =
      color_data->num_rows + color_data->nlev;
   JX_Int nnz_capacity = color_data->num_nonzeros;
   JX_Int size = 0;
   JX_Int c, i;

   if (bounds_capacity < 1)
   {
      bounds_capacity = 1;
   }
   if (nnz_capacity < 1)
   {
      nnz_capacity = 1;
   }

   pre = jx_CTAlloc(struct jx_Relax103SpMVData_struct, 1);
   if (pre == NULL)
   {
      return NULL;
   }

   pre->cluster_id = cluster_id;
   pre->num_rows = color_data->num_rows;
   pre->num_nonzeros = color_data->num_nonzeros;
   pre->num_colors = color_data->nlev;

   pre->task_offsets = jx_CTAlloc(JX_Int, color_data->nlev + 1);
   pre->task_bounds =
      (JX_Int *)hthread_malloc(cluster_id,
                               bounds_capacity * sizeof(JX_Int),
                               HT_MEM_RW);
   pre->dA_j =
      (JX_Int *)hthread_malloc(cluster_id,
                               nnz_capacity * sizeof(JX_Int),
                               HT_MEM_RW);
   pre->dy_data =
      (JX_Real *)hthread_malloc(cluster_id,
                                nnz_capacity * sizeof(JX_Real),
                                HT_MEM_RW);

   if (pre->task_offsets == NULL || pre->task_bounds == NULL ||
       pre->dA_j == NULL || pre->dy_data == NULL)
   {
      jx_Relax103SpMVDataDestroy(pre);
      return NULL;
   }

   /*
    * Each color is a contiguous row interval in the reordered CSR.  Build
    * one flattened collection of NNZ-balanced task bounds and remember the
    * segment belonging to every color.
    */
   for (c = 0; c < color_data->nlev; c++)
   {
      JX_Int start_row = color_data->ilev[c];
      JX_Int end_color_row = color_data->ilev[c + 1];

      pre->task_offsets[c] = size;
      pre->task_bounds[size++] = start_row;

      while (start_row < end_color_row)
      {
         JX_Int target_nnz = color_data->A_i[start_row] + AM_size;
         JX_Int current_row = start_row + 1;
         JX_Int end_row;

         if (color_data->A_i[end_color_row] <= target_nnz)
         {
            end_row = end_color_row;
         }
         else
         {
            while (current_row < end_color_row &&
                   color_data->A_i[current_row] < target_nnz)
            {
               current_row++;
            }

            if (current_row <= start_row + 1)
            {
               end_row = start_row + 1;
            }
            else
            {
               end_row = current_row - 1;
            }
         }

         if (end_row > end_color_row)
         {
            end_row = end_color_row;
         }
         if (end_row <= start_row)
         {
            end_row = start_row + 1;
         }

         pre->task_bounds[size++] = end_row;
         start_row = end_row;
      }
   }
   pre->task_offsets[color_data->nlev] = size;

   for (i = 0; i < color_data->num_nonzeros; i++)
   {
      pre->dA_j[i] = color_data->A_j[i] * 8;
   }

   return pre;
}

static JX_Int
jx_Relax103ColorMatvecCPU(jx_Relax103ColorData *color_data,
                          JX_Int color,
                          jx_Vector *x,
                          jx_Vector *y)
{
   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int k, jj, row;
   JX_Real row_sum;

#pragma omp parallel for private(k, jj, row, row_sum)
   for (k = color_data->ilev[color];
        k < color_data->ilev[color + 1];
        k++)
   {
      row_sum = 0.0;
      for (jj = color_data->A_i[k]; jj < color_data->A_i[k + 1]; jj++)
      {
         row_sum += color_data->A_data[jj] *
                    x_data[color_data->A_j[jj]];
      }

      row = color_data->jlev[k];
      y_data[row] = row_sum;
   }

   return 0;
}

static JX_Int
jx_Relax103ColorMatvecFP64(jx_Relax103ColorData *color_data,
                           JX_Int color,
                           jx_Vector *x,
                           jx_Vector *y,
                           JX_Int myid)
{
   struct jx_Relax103SpMVData_struct *pre =
      (struct jx_Relax103SpMVData_struct *)color_data->spmv_data;
   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int task_start;
   JX_Int num_tasks;
   JX_Int threads_id;
   JX_Int barrier_id;
   JX_Int k, jj, row;
   JX_Real row_sum;
   unsigned long args[10];

   if (pre == NULL)
   {
      pre = jx_Relax103SpMVDataCreate(color_data, myid);
      if (pre == NULL)
      {
         if (myid == 0)
         {
            fprintf(stderr,
                    "Fatal Error: relax-103 color SpMV preprocessing failed.\n");
         }
         return 1;
      }

      color_data->spmv_data = pre;
      color_data->spmv_data_destroy = jx_Relax103SpMVDataDestroy;
   }

   task_start = pre->task_offsets[color];
   num_tasks = pre->task_offsets[color + 1] - task_start - 1;

   threads_id = hthread_group_create(pre->cluster_id, coreNums);
   barrier_id = hthread_barrier_create(pre->cluster_id);

   args[0] = (unsigned long)pre->num_rows;
   args[1] = (unsigned long)coreNums;
   args[2] = (unsigned long)num_tasks;
   args[3] = (unsigned long)barrier_id;
   args[4] = (unsigned long)color_data->A_i;
   args[5] = (unsigned long)pre->dA_j;
   args[6] = (unsigned long)color_data->A_data;
   args[7] = (unsigned long)x_data;
   args[8] = (unsigned long)pre->dy_data;
   args[9] = (unsigned long)&pre->task_bounds[task_start];

   hthread_group_exec(threads_id, "SpMV_GSM_FP64", 4, 6, args);
   hthread_group_wait(threads_id);

#pragma omp parallel for private(k, jj, row, row_sum)
   for (k = color_data->ilev[color];
        k < color_data->ilev[color + 1];
        k++)
   {
      row_sum = 0.0;
      for (jj = color_data->A_i[k]; jj < color_data->A_i[k + 1]; jj++)
      {
         row_sum += pre->dy_data[jj];
      }

      row = color_data->jlev[k];
      y_data[row] = row_sum;
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   return 0;
}

JX_Int
jx_Relax103ColorMatvec(jx_Relax103ColorData *color_data,
                       JX_Int color,
                       jx_Vector *x,
                       jx_Vector *y,
                       JX_Int myid)
{
   if (color_data == NULL || color < 0 || color >= color_data->nlev)
   {
      return 1;
   }
   if (jx_VectorSize(x) != color_data->num_rows ||
       jx_VectorSize(y) != color_data->num_rows)
   {
      return 2;
   }

   if (jx_spmv_type == 1)
   {
      return jx_Relax103ColorMatvecFP64(color_data, color, x, y, myid);
   }
   if (jx_spmv_type != 0 && myid == 0)
   {
      fprintf(stderr,
              "Warning: relax-103 only supports jx_spmv_type 0/1; "
              "falling back to CPU.\n");
   }

   return jx_Relax103ColorMatvecCPU(color_data, color, x, y);
}

static jx_SpMVPrecondFP64Data *
jx_SpMVPrecondFP64DataCreate(jx_CSRMatrix *A,
                             JX_Int myid)
{
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   JX_Int A_nnz = jx_CSRMatrixNumNonzeros(A);

   JX_Int cluster_id = myid % 4;
   JX_Int AM_size = 96000 / 3;
   JX_Int capacity = num_rows + 1;

   jx_SpMVPrecondFP64Data *pre =
       jx_CTAlloc(jx_SpMVPrecondFP64Data, 1);

   if (pre == NULL)
      return NULL;

   pre->cluster_id = cluster_id;
   pre->num_rows = num_rows;
   pre->num_cols = num_cols;
   pre->A_nnz = A_nnz;
   pre->num_tasks = 0;

   pre->dA_j = NULL;
   pre->task_row_bounds = NULL;
   pre->dy_data = NULL;

   pre->task_row_bounds =
       (JX_Int *)hthread_malloc(cluster_id,
                                capacity * sizeof(JX_Int),
                                HT_MEM_RW);

   pre->dA_j =
       (JX_Int *)hthread_malloc(cluster_id,
                                sizeof(JX_Int) * A_nnz,
                                HT_MEM_RW);

   pre->dy_data =
       (JX_Real *)hthread_malloc(cluster_id,
                                 sizeof(JX_Real) * A_nnz,
                                 HT_MEM_RW);

   if (pre->task_row_bounds == NULL ||
       pre->dA_j == NULL ||
       pre->dy_data == NULL)
   {
      jx_SpMVPrecondFP64DataDestroy(pre);
      return NULL;
   }

   JX_Int size = 0;
   JX_Int start_row = 0;

   pre->task_row_bounds[size++] = 0;

   while (start_row < num_rows)
   {
      JX_Int target_nnz = A_i[start_row] + AM_size;
      JX_Int current_row = start_row + 1;
      JX_Int end_row;

      if (A_i[num_rows] <= target_nnz)
      {
         pre->task_row_bounds[size++] = num_rows;
         break;
      }

      while (current_row < num_rows && A_i[current_row] < target_nnz)
      {
         current_row++;
      }

      if (current_row <= start_row + 1)
      {
         end_row = start_row + 1;
      }
      else
      {
         end_row = current_row - 1;
      }

      if (end_row > num_rows)
         end_row = num_rows;

      if (end_row <= start_row)
         end_row = start_row + 1;

      pre->task_row_bounds[size++] = end_row;
      start_row = end_row;
   }

   pre->num_tasks = size - 1;

   for (JX_Int i = 0; i < A_nnz; i++)
   {
      pre->dA_j[i] = A_j[i] * 8;
   }

   return pre;
}

JX_Int
jx_CSRMatrixSpMVPrecondFP64Create(jx_CSRMatrix *A,
                                  JX_Int myid)
{
   if (A == NULL)
      return 1;

   if (jx_CSRMatrixSpMVPrecondFP64(A) != NULL)
      return 0;

   jx_CSRMatrixSpMVPrecondFP64(A) =
       jx_SpMVPrecondFP64DataCreate(A, myid);

   if (jx_CSRMatrixSpMVPrecondFP64(A) == NULL)
      return 1;

   return 0;
}

void jx_CSRMatrixSpMVPrecondFP64Destroy(jx_CSRMatrix *A)
{
   if (A == NULL)
      return;

   if (jx_CSRMatrixSpMVPrecondFP64(A) == NULL)
      return;

   jx_SpMVPrecondFP64DataDestroy(jx_CSRMatrixSpMVPrecondFP64(A));
   jx_CSRMatrixSpMVPrecondFP64(A) = NULL;
}

JX_Int
jx_CSRMatrixMatvec(JX_Real alpha,
                   jx_CSRMatrix *A,
                   jx_Vector *x,
                   JX_Real beta,
                   jx_Vector *y,
                   JX_Int myid)
{
   // 根据全局设置的 jx_spmv_type 来选择具体的实现
   switch (jx_spmv_type)
   {
   case 0:
      return jx_CSRMatrixMatvec_origin(alpha, A, x, beta, y);
   case 1: // GSM 版本
      return jx_CSRMatrixMatvec_v1(alpha, A, x, beta, y, myid);
   case 2: // x_new 版本
      return jx_CSRMatrixMatvec_v2(alpha, A, x, beta, y, myid);
   case 3:
      return jx_CSRMatrixMatvec_baseline(alpha, A, x, beta, y, myid);
   default:
      if (myid == 0)
         fprintf(stderr, "Warning: Unknown jx_spmv_type (%d). Falling back to the origin version.\n", jx_spmv_type);
      return jx_CSRMatrixMatvec_origin(alpha, A, x, beta, y);
   }
}

/*!
 * \fn JX_Int jx_CSRMatrixMatvec
 * \brief Perform y = alpha*A*x + beta*y.
 * \date 2011/09/03
 */
JX_Int
jx_CSRMatrixMatvec_origin(JX_Real alpha,
                          jx_CSRMatrix *A,
                          jx_Vector *x,
                          JX_Real beta,
                          jx_Vector *y)
{
   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);

   JX_Int *A_rownnz = jx_CSRMatrixRownnz(A);
   JX_Int num_rownnz = jx_CSRMatrixNumRownnz(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int x_size = jx_VectorSize(x);
   JX_Int y_size = jx_VectorSize(y);
   JX_Int num_vectors = jx_VectorNumVectors(x);
   JX_Int idxstride_y = jx_VectorIndexStride(y);
   JX_Int vecstride_y = jx_VectorVectorStride(y);
   JX_Int idxstride_x = jx_VectorIndexStride(x);
   JX_Int vecstride_x = jx_VectorVectorStride(x);

   JX_Real temp, tempx;

   JX_Int i, j, jj;

   JX_Int m;

   JX_Real xpar = 0.7;

   JX_Int ierr = 0;

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

   jx_assert(num_vectors == jx_VectorNumVectors(y));

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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] = 0.0;
      }
      else
      {
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
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

#define JX_SMP_PRIVATE i, jj, m, tempx
#include "../../../include/jx_smp_forloop.h"

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

#define JX_SMP_PRIVATE i, jj, temp
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}

/* y[offset:end] = alpha*A[offset:end,:]*x + beta*b[offset:end] */
JX_Int
jx_CSRMatrixMatvecOutOfPlace(JX_Real alpha,
                             jx_CSRMatrix *A,
                             jx_Vector *x,
                             JX_Real beta,
                             jx_Vector *b,
                             jx_Vector *y,
                             JX_Int offset)
{
   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A) + offset;
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A) - offset;
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   /*JX_Int         num_nnz  = jx_CSRMatrixNumNonzeros(A);*/

   JX_Int *A_rownnz = jx_CSRMatrixRownnz(A);
   JX_Int num_rownnz = jx_CSRMatrixNumRownnz(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *b_data = jx_VectorData(b) + offset;
   JX_Real *y_data = jx_VectorData(y) + offset;
   JX_Int x_size = jx_VectorSize(x);
   JX_Int b_size = jx_VectorSize(b) - offset;
   JX_Int y_size = jx_VectorSize(y) - offset;
   JX_Int num_vectors = jx_VectorNumVectors(x);
   JX_Int idxstride_y = jx_VectorIndexStride(y);
   JX_Int vecstride_y = jx_VectorVectorStride(y);
   /*JX_Int         idxstride_b = jx_VectorIndexStride(b);
   JX_Int         vecstride_b = jx_VectorVectorStride(b);*/
   JX_Int idxstride_x = jx_VectorIndexStride(x);
   JX_Int vecstride_x = jx_VectorVectorStride(x);
   JX_Real temp, tempx;
   JX_Int i, j, jj, m, ierr = 0;
   JX_Real xpar = 0.7;
   jx_Vector *x_tmp = NULL;

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

   jx_assert(num_vectors == jx_VectorNumVectors(y));
   jx_assert(num_vectors == jx_VectorNumVectors(b));

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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] = beta * b_data[i];

      return ierr;
   }

   if (x == y)
   {
      x_tmp = jx_SeqVectorCloneDeep(x);
      x_data = jx_VectorData(x_tmp);
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
            for (i = 0; i < num_rows * num_vectors; i++)
               y_data[i] = 0.0;
         }
         else
         {
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i, j, jj, m, tempx
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i, j, jj, tempx
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
         for (i = 0; i < num_rows * num_vectors; i++)
            y_data[i] *= alpha;
      }
   }
   else
   { // JSP: this is currently the only path optimized
#define JX_SMP_PRIVATE i, jj, tempx
#define JX_SMP_PAR_REGION
#include "../../../include/jx_smp_forloop.h"
      {
         JX_Int iBegin = jx_CSRMatrixGetLoadBalancedPartitionBegin(A);
         JX_Int iEnd = jx_CSRMatrixGetLoadBalancedPartitionEnd(A);
         jx_assert(iBegin <= iEnd);
         jx_assert(iBegin >= 0 && iBegin <= num_rows);
         jx_assert(iEnd >= 0 && iEnd <= num_rows);

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
      jx_SeqVectorDestroy(x_tmp);
   }

   return ierr;
}

/*!
 * \fn JX_Int jx_CSRMatrixMatvecT
 * \brief Perform y = alpha*A^T*x + beta*y (dispatcher with jx_spmv_type).
 * \date 2026/07
 */
JX_Int
jx_CSRMatrixMatvecT(JX_Real alpha,
                    jx_CSRMatrix *A,
                    jx_Vector *x,
                    JX_Real beta,
                    jx_Vector *y)
{
   switch (jx_spmv_type)
   {
   case 2: return jx_CSRMatrixMatvecT_v2(alpha, A, x, beta, y);
   default: return jx_CSRMatrixMatvecT_origin(alpha, A, x, beta, y);
   }
}

/*!
 * \fn JX_Int jx_CSRMatrixMatvecT_origin
 * \brief Perform y = alpha*A^T*x + beta*y.
 * \date 2011/09/03
 */
JX_Int
jx_CSRMatrixMatvecT_origin(JX_Real alpha,
                    jx_CSRMatrix *A,
                    jx_Vector *x,
                    JX_Real beta,
                    jx_Vector *y)
{
   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int x_size = jx_VectorSize(x);
   JX_Int y_size = jx_VectorSize(y);
   JX_Int num_vectors = jx_VectorNumVectors(x);
   JX_Int idxstride_y = jx_VectorIndexStride(y);
   JX_Int vecstride_y = jx_VectorVectorStride(y);
   JX_Int idxstride_x = jx_VectorIndexStride(x);
   JX_Int vecstride_x = jx_VectorVectorStride(x);

   JX_Real temp;

   JX_Int i, i1, j, jv, jj, ns, ne, size, rest;
   JX_Int num_threads;

   JX_Int ierr = 0;

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

   jx_assert(num_vectors == jx_VectorNumVectors(y));

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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] = 0.0;
         }
      }
      else
      {
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
         for (i = 0; i < num_cols * num_vectors; i++)
         {
            y_data[i] *= temp;
         }
      }
   }

   /*-----------------------------------------------------------------
    * y += A^T*x
    *-----------------------------------------------------------------*/
   num_threads = jx_NumThreads();
   if (num_threads > 1)
   {

#define JX_SMP_PRIVATE i, i1, jj, j, ns, ne, size, rest
#include "../../../include/jx_smp_forloop.h"
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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_cols * num_vectors; i++)
      {
         y_data[i] *= alpha;
      }
   }

   return ierr;
}


/*!
 * \fn JX_Int jx_CSRMatrixMatvecT_v2
 * \brief y = alpha*A^T*x + beta*y (DOT-split DSP-accelerated transpose).
 * \date 2026/07
 */
JX_Int
jx_CSRMatrixMatvecT_v2(JX_Real alpha,
                       jx_CSRMatrix *A,
                       jx_Vector *x,
                       JX_Real beta,
                       jx_Vector *y)
{
   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   JX_Int A_nnz = jx_CSRMatrixNumNonzeros(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int i, jj;
   JX_Int ierr = 0;

   if (num_rows != jx_VectorSize(x)) ierr = 1;
   if (num_cols != jx_VectorSize(y)) ierr = 2;
   if (num_rows != jx_VectorSize(x) && num_cols != jx_VectorSize(y)) ierr = 3;

   if (alpha == 0.0)
   {
      for (i = 0; i < num_cols; i++) y_data[i] *= beta;
      return ierr;
   }

   JX_Real rtemp = beta / alpha;
   if (rtemp == 0.0) { for (i = 0; i < num_cols; i++) y_data[i] = 0.0; }
   else if (rtemp != 1.0) { for (i = 0; i < num_cols; i++) y_data[i] *= rtemp; }

   JX_Int *row_of_nz = (JX_Int *)malloc(A_nnz * sizeof(JX_Int));
   for (i = 0; i < num_rows; i++) {
      for (jj = A_i[i]; jj < A_i[i + 1]; jj++) row_of_nz[jj] = i;
   }

   JX_Int cluster_id = 0;

   JX_Real *dx_data = (JX_Real *)hthread_malloc(cluster_id, sizeof(JX_Real) * A_nnz, HT_MEM_RW);
   for (jj = 0; jj < A_nnz; jj++) dx_data[jj] = x_data[row_of_nz[jj]];

   JX_Int AM_size = 96000 / 3;
   JX_Int capacity = num_rows + 1;
   JX_Int sz = 0, start_row = 0;
   JX_Int *task_row_bounds = (JX_Int *)hthread_malloc(cluster_id, capacity * sizeof(JX_Int), HT_MEM_RW);
   task_row_bounds[sz++] = 0;
   while (start_row < num_rows)
   {
      JX_Int target_nnz = A_i[start_row] + AM_size;
      JX_Int current_row = start_row + 1;
      JX_Int end_row;
      if (A_i[num_rows] <= target_nnz) { task_row_bounds[sz++] = num_rows; break; }
      while (current_row < num_rows && A_i[current_row] < target_nnz) current_row++;
      end_row = (current_row <= start_row + 1) ? start_row + 1 : current_row - 1;
      if (end_row > num_rows) end_row = num_rows;
      task_row_bounds[sz++] = end_row;
      start_row = end_row;
   }
   JX_Int num_tasks = sz - 1;

   JX_Real *dy_data = (JX_Real *)hthread_malloc(cluster_id, sizeof(JX_Real) * A_nnz, HT_MEM_RW);

   JX_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JX_Int barrier_id = hthread_barrier_create(cluster_id);
   unsigned long args[7];
   args[0] = coreNums;
   args[1] = num_tasks;
   args[2] = (unsigned long)task_row_bounds;
   args[3] = (unsigned long)A_i;
   args[4] = (unsigned long)A_data;
   args[5] = (unsigned long)dx_data;
   args[6] = (unsigned long)dy_data;
   hthread_group_exec(threads_id, "SpMV_DOT_FP64", 2, 5, args);
   hthread_group_wait(threads_id);

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   for (jj = 0; jj < A_nnz; jj++)
   {
      y_data[A_j[jj]] += dy_data[jj];
   }

   hthread_free(dx_data);
   hthread_free(dy_data);
   hthread_free(task_row_bounds);
   free(row_of_nz);

   if (alpha != 1.0)
   {
      for (i = 0; i < num_cols; i++) y_data[i] *= alpha;
   }

   return ierr;
}
// 2025.12.18 ============================
// GSM 版本 baseline  ====================
JX_Int
jx_CSRMatrixMatvec_baseline(JX_Real alpha,
                            jx_CSRMatrix *A,
                            jx_Vector *x,
                            JX_Real beta,
                            jx_Vector *y,
                            JX_Int myid)
{
   // if (myid == 0)
   //    printf("\n jx_CSRMatrixMatvec_baseline (GSM) \n");

   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   JX_Int A_nnz = jx_CSRMatrixNumNonzeros(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int x_size = jx_VectorSize(x);
   JX_Int y_size = jx_VectorSize(y);
   JX_Int num_vectors = jx_VectorNumVectors(x);
   JX_Int idxstride_y = jx_VectorIndexStride(y);
   JX_Int vecstride_y = jx_VectorVectorStride(y);
   JX_Int idxstride_x = jx_VectorIndexStride(x);
   JX_Int vecstride_x = jx_VectorVectorStride(x);

   JX_Real temp, tempx;
   JX_Int i, j, jj;
   JX_Int m;
   JX_Int ierr = 0;

   JX_Int cluster_id = myid % 4;

   jx_assert(num_vectors == jx_VectorNumVectors(y));

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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;

      return ierr;
   }

   /*-----------------------------------------------
    * y += A*x   天河平台实现
    *---------------------------------------------*/
   JX_Int *core_row_bounds = hthread_malloc(cluster_id, sizeof(JX_Int) * (coreNums + 1), HT_MEM_RW);
   JX_Int nnz_per_thread = A_nnz / coreNums;
   core_row_bounds[0] = 0;

   for (JX_Int tid = 0; tid < coreNums; tid++)
   {
      JX_Int start_row = core_row_bounds[tid];

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

      JX_Int target_nnz = A_i[start_row] + nnz_per_thread;
      JX_Int low = start_row + 1;
      JX_Int high = num_rows;
      JX_Int boundary = start_row + 1;

      while (low <= high)
      {
         JX_Int mid = low + (high - low) / 2;

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

   // if (myid == 0)
   // {
   //    for (JX_Int tid = 0; tid <= coreNums; tid++)
   //       printf("core_row_bounds[%d] = %d\n", tid, core_row_bounds[tid]);
   // }

   // 将 A_j 转换为按字节偏移
   JX_Int *dA_j = hthread_malloc(cluster_id, sizeof(JX_Int) * A_nnz, HT_MEM_RW);

   for (JX_Int i = 0; i < A_nnz; i++)
   {
      dA_j[i] = A_j[i] * 8;
   }

   // 设备端计算结果
   JX_Real *dy_data = hthread_malloc(cluster_id, sizeof(JX_Real) * num_rows, HT_MEM_RW);
   JX_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JX_Int barrier_id = hthread_barrier_create(cluster_id);

   unsigned long args[9];
   args[0] = (unsigned long)num_rows;
   args[1] = (unsigned long)coreNums;
   args[2] = (unsigned long)barrier_id;
   args[3] = (unsigned long)A_i;
   args[4] = (unsigned long)dA_j;
   args[5] = (unsigned long)A_data;
   args[6] = (unsigned long)x_data;
   args[7] = (unsigned long)dy_data;
   args[8] = (unsigned long)core_row_bounds;

   hthread_group_exec(threads_id, "spmv_baseline", 3, 6, args);
   hthread_group_wait(threads_id);

   // #define JX_SMP_PRIVATE i
   // #include "../../../include/jx_smp_forloop.h"
   for (i = 0; i < num_rows; i++)
   {
      y_data[i] = alpha * dy_data[i] + beta * y_data[i];
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   hthread_free(core_row_bounds);
   hthread_free(dA_j);
   hthread_free(dy_data);

   return 0;
}

JX_Int
jx_CSRMatrixMatvec_v1(JX_Real alpha,
                      jx_CSRMatrix *A,
                      jx_Vector *x,
                      JX_Real beta,
                      jx_Vector *y,
                      JX_Int myid)
{
   // if (myid == 0)
   //    jxf_printf("jx_CSRMatrixMatvec_v1 \n");

   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int x_size = jx_VectorSize(x);
   JX_Int y_size = jx_VectorSize(y);
   JX_Int num_vectors = jx_VectorNumVectors(x);

   JX_Int i, j;
   JX_Int ierr = 0;

   jx_SpMVPrecondFP64Data *pre = NULL;

   jx_assert(num_vectors == jx_VectorNumVectors(y));

   if (num_cols != x_size)
      ierr = 1;

   if (num_rows != y_size)
      ierr = 2;

   if (num_cols != x_size && num_rows != y_size)
      ierr = 3;

   if (alpha == 0.0)
   {
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
      {
         y_data[i] *= beta;
      }

      return ierr;
   }

   pre = jx_CSRMatrixSpMVPrecondFP64(A);

   if (pre == NULL)
   {
      ierr = jx_CSRMatrixSpMVPrecondFP64Create(A, myid);
      if (ierr)
      {
         if (myid == 0)
         {
            fprintf(stderr, "Fatal Error: jx_CSRMatrixSpMVPrecondFP64Create failed.\n");
         }
         return ierr;
      }

      pre = jx_CSRMatrixSpMVPrecondFP64(A);
   }

   if (pre == NULL)
      return 1;

   JX_Int cluster_id = myid % 4;
   JX_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JX_Int barrier_id = hthread_barrier_create(cluster_id);

   unsigned long args[10];

   args[0] = (unsigned long)pre->num_rows;
   args[1] = (unsigned long)coreNums;
   args[2] = (unsigned long)pre->num_tasks;
   args[3] = (unsigned long)barrier_id;
   args[4] = (unsigned long)A_i;
   args[5] = (unsigned long)pre->dA_j;
   args[6] = (unsigned long)A_data;
   args[7] = (unsigned long)x_data;
   args[8] = (unsigned long)pre->dy_data;
   args[9] = (unsigned long)pre->task_row_bounds;

   hthread_group_exec(threads_id, "SpMV_GSM_FP64", 4, 6, args);
   hthread_group_wait(threads_id);

   JX_Real row_sum;

#define JX_SMP_PRIVATE i, j, row_sum
#include "../../../include/jx_smp_forloop.h"
   for (i = 0; i < pre->num_rows; i++)
   {
      row_sum = 0.0;

      for (j = A_i[i]; j < A_i[i + 1]; j++)
      {
         row_sum += pre->dy_data[j];
      }

      if (beta == 0.0)
      {
         y_data[i] = alpha * row_sum;
      }
      else
      {
         y_data[i] = beta * y_data[i] + alpha * row_sum;
      }
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   return ierr;
}

/* 乘加分裂版本, 不使用 GSM */
JX_Int
jx_CSRMatrixMatvec_v2(JX_Real alpha,
                      jx_CSRMatrix *A,
                      jx_Vector *x,
                      JX_Real beta,
                      jx_Vector *y,
                      JX_Int myid)
{
   // if (myid == 0)
   //    jxf_printf("jx_CSRMatrixMatvec_v2 \n");

   JX_Real *A_data = jx_CSRMatrixData(A);
   JX_Int *A_i = jx_CSRMatrixI(A);
   JX_Int *A_j = jx_CSRMatrixJ(A);
   JX_Int num_rows = jx_CSRMatrixNumRows(A);
   JX_Int num_cols = jx_CSRMatrixNumCols(A);
   JX_Int A_nnz = jx_CSRMatrixNumNonzeros(A);

   JX_Real *x_data = jx_VectorData(x);
   JX_Real *y_data = jx_VectorData(y);
   JX_Int x_size = jx_VectorSize(x);
   JX_Int y_size = jx_VectorSize(y);
   JX_Int num_vectors = jx_VectorNumVectors(x);
   JX_Int idxstride_y = jx_VectorIndexStride(y);
   JX_Int vecstride_y = jx_VectorVectorStride(y);
   JX_Int idxstride_x = jx_VectorIndexStride(x);
   JX_Int vecstride_x = jx_VectorVectorStride(x);

   JX_Real temp, tempx;
   JX_Int i, j, jj;
   JX_Int m;
   JX_Int ierr = 0;

   JX_Int cluster_id = myid % 4;

   jx_assert(num_vectors == jx_VectorNumVectors(y));

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
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < num_rows * num_vectors; i++)
         y_data[i] *= beta;

      return ierr;
   }
   /*-----------------------------------------------
    * y += A*x   天河平台实现
    *---------------------------------------------*/
   JX_Int AM_size = 96000 / 3;
   JX_Int capacity = num_rows + 1;
   JX_Int size = 0;
   JX_Int start_row = 0;

   JX_Int *task_row_bounds = (JX_Int *)hthread_malloc(cluster_id, capacity * sizeof(JX_Int), HT_MEM_RW);
   task_row_bounds[size++] = 0;
   while (start_row < num_rows)
   {
      JX_Int target_nnz = A_i[start_row] + AM_size;
      JX_Int current_row = start_row + 1;
      JX_Int end_row;

      if (A_i[num_rows] <= target_nnz)
      {
         task_row_bounds[size++] = num_rows;
         break;
      }

      while (current_row < num_rows && A_i[current_row] < target_nnz)
      {
         current_row++;
      }

      if (current_row <= start_row + 1)
      {
         end_row = start_row + 1;
      }
      else
      {
         end_row = current_row - 1;
      }

      if (end_row > num_rows)
         end_row = num_rows;

      task_row_bounds[size++] = end_row;
      start_row = end_row;
   }

   JX_Int num_tasks = size - 1;

   JX_Real *dx_data = hthread_malloc(cluster_id, sizeof(JX_Real) * A_nnz, HT_MEM_RW);
#pragma omp parallel for private(jj) num_threads(4)
   for (jj = 0; jj < A_nnz; jj++)
   {
      dx_data[jj] = x_data[A_j[jj]];
   }

   JX_Real *dy_data = hthread_malloc(cluster_id, sizeof(JX_Real) * A_nnz, HT_MEM_RW);
   // memset(dy_data, 0.0, sizeof(JX_Real) * A_nnz);

   JX_Int threads_id = hthread_group_create(cluster_id, coreNums);
   JX_Int barrier_id = hthread_barrier_create(cluster_id);
   unsigned long args[7];
   args[0] = (unsigned long)coreNums;
   args[1] = (unsigned long)num_tasks;
   args[2] = (unsigned long)A_i;
   args[3] = (unsigned long)A_data;
   args[4] = (unsigned long)dx_data;
   args[5] = (unsigned long)dy_data;
   args[6] = (unsigned long)task_row_bounds;

   hthread_group_exec(threads_id, "SpMV_DOT_FP64", 2, 5, args);
   hthread_group_wait(threads_id);

   // step 4 : 行规约
   JX_Real row_sum;
#define JX_SMP_PRIVATE i, j, row_sum
#include "../../../include/jx_smp_forloop.h"
   for (i = 0; i < num_rows; i++)
   {
      row_sum = 0.0;
      for (j = A_i[i]; j < A_i[i + 1]; j++)
      {
         row_sum += dy_data[j];
      }

      if (beta == 0.0)
      {
         y_data[i] = alpha * row_sum;
      }
      else
      {
         y_data[i] = beta * y_data[i] + alpha * row_sum;
      }
   }

   hthread_barrier_destroy(barrier_id);
   hthread_group_destroy(threads_id);

   hthread_free(dx_data);
   hthread_free(dy_data);
   hthread_free(task_row_bounds);

   return 0;
}
