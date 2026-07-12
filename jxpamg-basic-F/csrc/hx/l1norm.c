//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  l1norm.c
 *  Date: 2024/01/29
 */

#include "jxf_mv.h"

JXF_Int
jxf_ParCSRComputeL1Norms( jxf_ParCSRMatrix *A, JXF_Int option, JXF_Int *cf_marker, JXF_Real **l1_norm_ptr )
{
   JXF_Int i, j;
   JXF_Int num_rows = jxf_ParCSRMatrixNumRows(A);

   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(A);
   JXF_Int *A_diag_I = jxf_CSRMatrixI(A_diag);
   JXF_Int *A_diag_J = jxf_CSRMatrixJ(A_diag);
   JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);

   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(A);
   JXF_Int *A_offd_I = jxf_CSRMatrixI(A_offd);
   JXF_Int *A_offd_J = jxf_CSRMatrixJ(A_offd);
   JXF_Real *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   JXF_Real diag;
   JXF_Real *l1_norm = jxf_TAlloc(JXF_Real, num_rows);

   JXF_Int *cf_marker_offd = NULL;
   JXF_Int cf_diag;

   if (cf_marker != NULL)
   {
      JXF_Int index;
      JXF_Int num_sends;
      JXF_Int start;
      JXF_Int *int_buf_data = NULL;

      jxf_ParCSRCommPkg  *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
      jxf_ParCSRCommHandle *comm_handle;

      if (num_cols_offd)
         cf_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      if (jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends))
         int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
         {
            int_buf_data[index++] = cf_marker[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
         }
      }
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, cf_marker_offd);
      jxf_ParCSRCommHandleDestroy(comm_handle);
      jxf_TFree(int_buf_data);
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
         jxf_error_in_arg(1);
         break;
      }

   jxf_TFree(cf_marker_offd);

  *l1_norm_ptr = l1_norm;

   return jxf_error_flag;
}

JXF_Int
jxf_ParCSRComputeL1NormsThreads( jxf_ParCSRMatrix *A, JXF_Int option, JXF_Int num_threads, JXF_Int *cf_marker, JXF_Real **l1_norm_ptr )
{
   JXF_Int i, j, k;
   JXF_Int num_rows = jxf_ParCSRMatrixNumRows(A);

   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(A);
   JXF_Int *A_diag_I = jxf_CSRMatrixI(A_diag);
   JXF_Int *A_diag_J = jxf_CSRMatrixJ(A_diag);
   JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);

   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(A);
   JXF_Int *A_offd_I = jxf_CSRMatrixI(A_offd);
   JXF_Int *A_offd_J = jxf_CSRMatrixJ(A_offd);
   JXF_Real *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   JXF_Real diag = 0.0;
   JXF_Real *l1_norm = jxf_CTAlloc(JXF_Real, num_rows);
   JXF_Int ii, ns, ne, rest, size;

   JXF_Int *cf_marker_offd = NULL;
   JXF_Int cf_diag;

   if (cf_marker != NULL)
   {
      JXF_Int index;
      JXF_Int num_sends;
      JXF_Int start;
      JXF_Int *int_buf_data = NULL;

      jxf_ParCSRCommPkg  *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
      jxf_ParCSRCommHandle *comm_handle;

      if (num_cols_offd)
         cf_marker_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      if (jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends))
         int_buf_data = jxf_CTAlloc(JXF_Int, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
         {
            int_buf_data[index++] = cf_marker[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
         }
      }
      comm_handle = jxf_ParCSRCommHandleCreate(11, comm_pkg, int_buf_data, cf_marker_offd);
      jxf_ParCSRCommHandleDestroy(comm_handle);
      jxf_TFree(int_buf_data);
   }

#define JXF_SMP_PRIVATE i,ii,j,k,ns,ne,rest,size,diag,cf_diag
#include "../../include/jxf_smp_forloop.h"
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
            jxf_error_in_arg(1);
            break;
         }
   }

   jxf_TFree(cf_marker_offd);

  *l1_norm_ptr = l1_norm;

   return jxf_error_flag;
}
