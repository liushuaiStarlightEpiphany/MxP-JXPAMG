//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  l1norm.c
 *  Date: 2024/01/29
 */

#include "jx_mv.h"

JX_Int
jx_ParCSRComputeL1Norms( jx_ParCSRMatrix *A, JX_Int option, JX_Int *cf_marker, JX_Real **l1_norm_ptr )
{
   JX_Int i, j;
   JX_Int num_rows = jx_ParCSRMatrixNumRows(A);

   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   JX_Int *A_diag_I = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_J = jx_CSRMatrixJ(A_diag);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);

   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_offd_I = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_J = jx_CSRMatrixJ(A_offd);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int num_cols_offd = jx_CSRMatrixNumCols(A_offd);

   JX_Real diag;
   JX_Real *l1_norm = jx_TAlloc(JX_Real, num_rows);

   JX_Int *cf_marker_offd = NULL;
   JX_Int cf_diag;

   if (cf_marker != NULL)
   {
      JX_Int index;
      JX_Int num_sends;
      JX_Int start;
      JX_Int *int_buf_data = NULL;

      jx_ParCSRCommPkg  *comm_pkg = jx_ParCSRMatrixCommPkg(A);
      jx_ParCSRCommHandle *comm_handle;

      if (num_cols_offd)
         cf_marker_offd = jx_CTAlloc(JX_Int, num_cols_offd);
      num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
      if (jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends))
         int_buf_data = jx_CTAlloc(JX_Int, jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
         {
            int_buf_data[index++] = cf_marker[jx_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
         }
      }
      comm_handle = jx_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, cf_marker_offd);
      jx_ParCSRCommHandleDestroy(comm_handle);
      jx_TFree(int_buf_data);
   }

   if (option == 1)
   {
      for (i = 0; i < num_rows; i ++)
      {
         l1_norm[i] = 0.0;
         if (cf_marker == NULL)
         {
            for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               l1_norm[i] += fabs(A_diag_data[j]);
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  l1_norm[i] += fabs(A_offd_data[j]);
            }
         }
         else
         {
            cf_diag = cf_marker[i];
            for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               if (cf_diag == cf_marker[A_diag_J[j]])
                  l1_norm[i] += fabs(A_diag_data[j]);
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  if (cf_diag == cf_marker_offd[A_offd_J[j]])
                     l1_norm[i] += fabs(A_offd_data[j]);
            }
         }
      }
   }
   else if (option == 2)
   {
      for (i = 0; i < num_rows; i ++)
      {
         l1_norm[i] = fabs(A_diag_data[A_diag_I[i]]);
         if (cf_marker == NULL)
         {
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  l1_norm[i] += fabs(A_offd_data[j]);
            }
         }
         else
         {
            cf_diag = cf_marker[i];
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  if (cf_diag == cf_marker_offd[A_offd_J[j]])
                     l1_norm[i] += fabs(A_offd_data[j]);
            }
         }
      }
   }
   else if (option == 3)
   {
      for (i = 0; i < num_rows; i ++)
      {
         l1_norm[i] = 0.0;
         for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
            l1_norm[i] += A_diag_data[j] * A_diag_data[j];
         if (num_cols_offd)
            for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
               l1_norm[i] += A_offd_data[j] * A_offd_data[j];
      }
   }
   else if (option == 4)
   {
      for (i = 0; i < num_rows; i ++)
      {
         diag = l1_norm[i] = fabs(A_diag_data[A_diag_I[i]]);
         if (cf_marker == NULL)
         {
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  l1_norm[i] += 0.5*fabs(A_offd_data[j]);
            }
         }
         else
         {
            cf_diag = cf_marker[i];
            if (num_cols_offd)
            {
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  if (cf_diag == cf_marker_offd[A_offd_J[j]])
                     l1_norm[i] += 0.5*fabs(A_offd_data[j]);
            }
         }

         if (l1_norm[i] <= 4.0/3.0*diag)
            l1_norm[i] = diag;
      }
   }

   for (i = 0; i < num_rows; i ++)
      if (A_diag_data[A_diag_I[i]] < 0)
         l1_norm[i] = -l1_norm[i];

   for (i = 0; i < num_rows; i ++)
      if (fabs(l1_norm[i]) == 0.0)
      {
         jx_error_in_arg(1);
         break;
      }

   jx_TFree(cf_marker_offd);

  *l1_norm_ptr = l1_norm;

   return jx_error_flag;
}

JX_Int
jx_ParCSRComputeL1NormsThreads( jx_ParCSRMatrix *A, JX_Int option, JX_Int num_threads, JX_Int *cf_marker, JX_Real **l1_norm_ptr )
{
   JX_Int i, j, k;
   JX_Int num_rows = jx_ParCSRMatrixNumRows(A);

   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   JX_Int *A_diag_I = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_J = jx_CSRMatrixJ(A_diag);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);

   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_offd_I = jx_CSRMatrixI(A_offd);
   JX_Int *A_offd_J = jx_CSRMatrixJ(A_offd);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int num_cols_offd = jx_CSRMatrixNumCols(A_offd);

   JX_Real diag = 0.0;
   JX_Real *l1_norm = jx_CTAlloc(JX_Real, num_rows);
   JX_Int ii, ns, ne, rest, size;

   JX_Int *cf_marker_offd = NULL;
   JX_Int cf_diag;

   if (cf_marker != NULL)
   {
      JX_Int index;
      JX_Int num_sends;
      JX_Int start;
      JX_Int *int_buf_data = NULL;

      jx_ParCSRCommPkg  *comm_pkg = jx_ParCSRMatrixCommPkg(A);
      jx_ParCSRCommHandle *comm_handle;

      if (num_cols_offd)
         cf_marker_offd = jx_CTAlloc(JX_Int, num_cols_offd);
      num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
      if (jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends))
         int_buf_data = jx_CTAlloc(JX_Int, jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
         {
            int_buf_data[index++] = cf_marker[jx_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
         }
      }
      comm_handle = jx_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, cf_marker_offd);
      jx_ParCSRCommHandleDestroy(comm_handle);
      jx_TFree(int_buf_data);
   }

#define JX_SMP_PRIVATE i,ii,j,k,ns,ne,rest,size,diag,cf_diag
#include "../../include/jx_smp_forloop.h"
   for (k = 0; k < num_threads; k ++)
   {
      size = num_rows / num_threads;
      rest = num_rows - size * num_threads;
      if (k < rest)
      {
         ns = k * size + k;
         ne = (k + 1) * size + k + 1;
      }
      else
      {
         ns = k * size + rest;
         ne = (k + 1) * size + rest;
      }

      if (option == 1)
      {
         for (i = ns; i < ne; i ++)
         {
            l1_norm[i] = 0.0;
            if (cf_marker == NULL)
            {
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
                  l1_norm[i] += fabs(A_diag_data[j]);
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     l1_norm[i] += fabs(A_offd_data[j]);
               }
            }
            else
            {
               cf_diag = cf_marker[i];
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
                  if (cf_diag == cf_marker[A_diag_J[j]])
                     l1_norm[i] += fabs(A_diag_data[j]);
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     if (cf_diag == cf_marker_offd[A_offd_J[j]])
                        l1_norm[i] += fabs(A_offd_data[j]);
               }
            }
         }
      }
      else if (option == 2)
      {
         for (i = ns; i < ne; i ++)
         {
            l1_norm[i] = 0.0;
            if (cf_marker == NULL)
            {
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j++)
               {
                  ii = A_diag_J[j];
                  if (ii == i || ii < ns || ii >= ne)
                     l1_norm[i] += fabs(A_diag_data[j]);
               }
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     l1_norm[i] += fabs(A_offd_data[j]);
               }
            }
            else
            {
               cf_diag = cf_marker[i];
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               {
                  ii = A_diag_J[j];
                  if ((ii == i || ii < ns || ii >= ne) && (cf_diag == cf_marker[A_diag_J[j]]))
                     l1_norm[i] += fabs(A_diag_data[j]);
               }
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     if (cf_diag == cf_marker_offd[A_offd_J[j]])
                        l1_norm[i] += fabs(A_offd_data[j]);
               }
            }
         }
      }
      else if (option == 3)
      {
         for (i = ns; i < ne; i ++)
         {
            l1_norm[i] = 0.0;
            for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               l1_norm[i] += A_diag_data[j] * A_diag_data[j];
            if (num_cols_offd)
               for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                  l1_norm[i] += A_offd_data[j] * A_offd_data[j];
         }
      }
      else if (option == 4)
      {
         for (i = ns; i < ne; i ++)
         {
            l1_norm[i] = 0.0;
            if (cf_marker == NULL)
            {
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               {
                  ii = A_diag_J[j];
                  if (ii == i || ii < ns || ii >= ne)
                  {
                     if (ii == i)
                     {
                        diag = fabs(A_diag_data[j]);
                        l1_norm[i] += fabs(A_diag_data[j]);
                     }
                     else
                        l1_norm[i] += 0.5 * fabs(A_diag_data[j]);
                  }
               }
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     l1_norm[i] += 0.5 * fabs(A_offd_data[j]);
               }
            }
            else
            {
               cf_diag = cf_marker[i];
               for (j = A_diag_I[i]; j < A_diag_I[i+1]; j ++)
               {
                  ii = A_diag_J[j];
                  if ((ii == i || ii < ns || ii >= ne) && (cf_diag == cf_marker[A_diag_J[j]]))
                  {
                     if (ii == i)
                     {
                        diag = fabs(A_diag_data[j]);
                        l1_norm[i] += fabs(A_diag_data[j]);
                     }
                     else
                        l1_norm[i] += 0.5 * fabs(A_diag_data[j]);
                  }
               }
               if (num_cols_offd)
               {
                  for (j = A_offd_I[i]; j < A_offd_I[i+1]; j ++)
                     if (cf_diag == cf_marker_offd[A_offd_J[j]])
                        l1_norm[i] += 0.5 * fabs(A_offd_data[j]);
               }
            }

            if (l1_norm[i] <= 4.0/3.0*diag)
               l1_norm[i] = diag;
         }
      }

      for (i = ns; i < ne; i ++)
         if (A_diag_data[A_diag_I[i]] < 0)
            l1_norm[i] = -l1_norm[i];

      for (i = ns; i < ne; i ++)
         if (fabs(l1_norm[i]) == 0.0)
         {
            jx_error_in_arg(1);
            break;
         }
   }

   jx_TFree(cf_marker_offd);

  *l1_norm_ptr = l1_norm;

   return jx_error_flag;
}
