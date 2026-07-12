//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_8.c
 *  Date: 2024/01/31
 */ 

#include "jx_hpcsr.h"

JX_Int
jx_PAMGRelax8( jx_ParCSRMatrix *A,
               jx_ParVector *f,
               JX_Int *cf_marker,
               JX_Int relax_points,
               JX_Real relax_weight,
               JX_Real omega,
               JX_Real *l1_norms,
               jx_ParVector *u,
               jx_ParVector *Vtemp,
               jx_ParVector *Ztemp )
{ 
   MPI_Comm	comm = jx_ParCSRMatrixComm(A);
   jx_CSRMatrix *A_diag = jx_ParCSRMatrixDiag(A);
   JX_Real *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int *A_diag_i = jx_CSRMatrixI(A_diag);
   JX_Int *A_diag_j = jx_CSRMatrixJ(A_diag);
   jx_CSRMatrix *A_offd = jx_ParCSRMatrixOffd(A);
   JX_Int *A_offd_i = jx_CSRMatrixI(A_offd);
   JX_Real *A_offd_data = jx_CSRMatrixData(A_offd);
   JX_Int *A_offd_j = jx_CSRMatrixJ(A_offd);
   jx_ParCSRCommPkg *comm_pkg = jx_ParCSRMatrixCommPkg(A);
   jx_ParCSRCommHandle *comm_handle;

   JX_Int n = jx_CSRMatrixNumRows(A_diag);
   JX_Int num_cols_offd = jx_CSRMatrixNumCols(A_offd);

   jx_Vector *u_local = jx_ParVectorLocalVector(u);
   JX_Real *u_data  = jx_VectorData(u_local);

   jx_Vector *f_local = jx_ParVectorLocalVector(f);
   JX_Real *f_data  = jx_VectorData(f_local);

   jx_Vector *Vtemp_local = jx_ParVectorLocalVector(Vtemp);
   JX_Real *Vtemp_data = jx_VectorData(Vtemp_local);
   JX_Real *Vext_data = NULL;
   JX_Real *v_buf_data = NULL;
   JX_Real *tmp_data;

   jx_Vector *Ztemp_local;
   JX_Real *Ztemp_data;

   JX_Int i, j;
   JX_Int ii, jj;
   JX_Int ns, ne, size, rest;
   JX_Int relax_error = 0;
   JX_Int num_sends;
   JX_Int index, start;
   JX_Int num_procs, num_threads, my_id;

   JX_Real zero = 0.0;
   JX_Real res, res0, res2;
   JX_Real one_minus_omega;
   JX_Real prod;

   one_minus_omega = 1.0 - omega;

   jx_MPI_Comm_size(comm, &num_procs);
   jx_MPI_Comm_rank(comm, &my_id);
   num_threads = jx_NumThreads();

   if (num_threads > 1)
   {
      Ztemp_local = jx_ParVectorLocalVector(Ztemp);
      Ztemp_data = jx_VectorData(Ztemp_local);
   }

   if (num_procs > 1)
   {
      num_sends = jx_ParCSRCommPkgNumSends(comm_pkg);
      v_buf_data = jx_CTAlloc(JX_Real, jx_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      Vext_data = jx_CTAlloc(JX_Real, num_cols_offd);

      if (num_cols_offd)
      {
         A_offd_j = jx_CSRMatrixJ(A_offd);
         A_offd_data = jx_CSRMatrixData(A_offd);
      }

      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
            v_buf_data[index++] = u_data[jx_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
      }

      comm_handle = jx_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);

      jx_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;
   }

   if (relax_weight == 1 && omega == 1)
   {
      if (relax_points == 0)
      {
         if (num_threads > 1)
         {
            tmp_data = Ztemp_data;
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JX_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jx_smp_forloop.h"
            for (j = 0; j < num_threads; j ++)
            {
               size = n / num_threads;
               rest = n - size * num_threads;
               if (j < rest)
               {
                  ns = j * size + j;
                  ne = (j + 1) * size + j + 1;
               }
               else
               {
                  ns = j * size + rest;
                  ne = (j + 1) * size + rest;
               }
               for (i = ns; i < ne; i ++)
               {
                  if (l1_norms[i] != zero)
                  {
                     res = f_data[i];
                     for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] += res / l1_norms[i];
                  }
               }
               for (i = ne-1; i > ns-1; i --)
               {
                  if (l1_norms[i] != zero)
                  {
                     res = f_data[i];
                     for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] += res / l1_norms[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < n; i ++)
            {
               if ( l1_norms[i] != zero)
               {
                  res = f_data[i];
                  for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res -= A_diag_data[jj] * u_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] += res / l1_norms[i];
               }
            }
            for (i = n-1; i > -1; i --)
            {
               if ( l1_norms[i] != zero)
               {
                  res = f_data[i];
                  for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res -= A_diag_data[jj] * u_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] += res / l1_norms[i];
               }
            }
         }
      }
      else
      {
         if (num_threads > 1)
         {
            tmp_data = Ztemp_data;
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JX_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jx_smp_forloop.h"
            for (j = 0; j < num_threads; j ++)
            {
               size = n / num_threads;
               rest = n - size * num_threads;
               if (j < rest)
               {
                  ns = j * size + j;
                  ne = (j + 1) * size + j + 1;
               }
               else
               {
                  ns = j * size + rest;
                  ne = (j + 1) * size + rest;
               }
               for (i = ns; i < ne; i ++)
               {
                  if (cf_marker[i] == relax_points && l1_norms[i] != zero)
                  {
                     res = f_data[i];
                     for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] += res / l1_norms[i];
                  }
               }
               for (i = ne-1; i > ns-1; i --)
               {
                  if (cf_marker[i] == relax_points && l1_norms[i] != zero)
                  {
                     res = f_data[i];
                     for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] += res / l1_norms[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < n; i ++)
            {
               if (cf_marker[i] == relax_points && l1_norms[i] != zero)
               {
                  res = f_data[i];
                  for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res -= A_diag_data[jj] * u_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] += res / l1_norms[i];
               }
            }
            for (i = n-1; i > -1; i --)
            {
               if (cf_marker[i] == relax_points && l1_norms[i] != zero)
               {
                  res = f_data[i];
                  for (jj = A_diag_i[i]; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res -= A_diag_data[jj] * u_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] += res / l1_norms[i];
               }
            }
         }
      }
   }
   else
   {
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
      for (i = 0; i < n; i ++)
      {
         Vtemp_data[i] = u_data[i];
      }
      prod = (1.0 - relax_weight * omega);
      if (relax_points == 0)
      {
         if (num_threads > 1)
         {
            tmp_data = Ztemp_data;
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JX_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jx_smp_forloop.h"
            for (j = 0; j < num_threads; j ++)
            {
               size = n / num_threads;
               rest = n - size * num_threads;
               if (j < rest)
               {
                  ns = j * size + j;
                  ne = (j + 1) * size + j + 1;
               }
               else
               {
                  ns = j * size + rest;
                  ne = (j + 1) * size + rest;
               }
               for (i = ns; i < ne; i ++)
               {
                  if (l1_norms[i] != zero)
                  {
                     res0 = 0.0;
                     res2 = 0.0;
                     res = f_data[i];
                     for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res0 -= A_diag_data[jj] * u_data[ii];
                           res2 += A_diag_data[jj] * Vtemp_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] *= prod;
                     u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
                  }
               }
               for (i = ne-1; i > ns-1; i --)
               {
                  if (l1_norms[i] != zero)
                  {
                     res0 = 0.0;
                     res2 = 0.0;
                     res = f_data[i];
                     for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res0 -= A_diag_data[jj] * u_data[ii];
                           res2 += A_diag_data[jj] * Vtemp_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] *= prod;
                     u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < n; i ++)
            {
               if (l1_norms[i] != zero)
               {
                  res0 = 0.0;
                  res = f_data[i];
                  res2 = 0.0;
                  for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res0 -= A_diag_data[jj] * u_data[ii];
                     res2 += A_diag_data[jj] * Vtemp_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] *= prod;
                  u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
               }
            }
            for (i = n-1; i > -1; i --)
            {
               if (l1_norms[i] != zero)
               {
                  res0 = 0.0;
                  res = f_data[i];
                  res2 = 0.0;
                  for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res0 -= A_diag_data[jj] * u_data[ii];
                     res2 += A_diag_data[jj] * Vtemp_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] *= prod;
                  u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
               }
            }
         }
      }
      else
      {
         if (num_threads > 1)
         {
            tmp_data = Ztemp_data;
#define JX_SMP_PRIVATE i
#include "../../../include/jx_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JX_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jx_smp_forloop.h"
            for (j = 0; j < num_threads; j ++)
            {
               size = n / num_threads;
               rest = n - size * num_threads;
               if (j < rest)
               {
                  ns = j * size + j;
                  ne = (j + 1) * size + j + 1;
               }
               else
               {
                  ns = j * size + rest;
                  ne = (j + 1) * size + rest;
               }
               for (i = ns; i < ne; i ++)
               {
                  if (cf_marker[i] == relax_points && l1_norms[i] != zero)
                  {
                     res0 = 0.0;
                     res2 = 0.0;
                     res = f_data[i];
                     for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res2 += A_diag_data[jj] * Vtemp_data[ii];
                           res0 -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] *= prod;
                     u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
                  }
               }
               for (i = ne-1; i > ns-1; i --)
               {
                  if (cf_marker[i] == relax_points && l1_norms[i] != zero)
                  {
                     res0 = 0.0;
                     res2 = 0.0;
                     res = f_data[i];
                     for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                     {
                        ii = A_diag_j[jj];
                        if (ii >= ns && ii < ne)
                        {
                           res2 += A_diag_data[jj] * Vtemp_data[ii];
                           res0 -= A_diag_data[jj] * u_data[ii];
                        }
                        else
                           res -= A_diag_data[jj] * tmp_data[ii];
                     }
                     for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                     {
                        ii = A_offd_j[jj];
                        res -= A_offd_data[jj] * Vext_data[ii];
                     }
                     u_data[i] *= prod;
                     u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
                  }
               }
            }
         }
         else
         {
            for (i = 0; i < n; i ++)
            {
               if (cf_marker[i] == relax_points && l1_norms[i] != zero)
               {
                  res = f_data[i];
                  res0 = 0.0;
                  res2 = 0.0;
                  for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res0 -= A_diag_data[jj] * u_data[ii];
                     res2 += A_diag_data[jj] * Vtemp_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] *= prod;
                  u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
               }
            }
            for (i = n-1; i > -1; i --)
            {
               if (cf_marker[i] == relax_points && l1_norms[i] != zero)
               {
                  res = f_data[i];
                  res0 = 0.0;
                  res2 = 0.0;
                  for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                  {
                     ii = A_diag_j[jj];
                     res0 -= A_diag_data[jj] * u_data[ii];
                     res2 += A_diag_data[jj] * Vtemp_data[ii];
                  }
                  for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                  {
                     ii = A_offd_j[jj];
                     res -= A_offd_data[jj] * Vext_data[ii];
                  }
                  u_data[i] *= prod;
                  u_data[i] += relax_weight * (omega * res + res0 + one_minus_omega * res2) / l1_norms[i];
               }
            }
         }
      }
   }

   if (num_procs > 1)
   {
      jx_TFree(Vext_data);
      jx_TFree(v_buf_data);
   }

   return(relax_error);
}

JX_Int
jx_hpPAMGRelax8( jx_hpCSRMatrix *A,
               jx_ParVector *f,
               JX_Int *cf_marker,
               JX_Int relax_points,
               JX_Real relax_weight,
               JX_Real omega,
               JX_Real *l1_norms,
               jx_ParVector *u,
               jx_ParVector *Vtemp,
               jx_ParVector *Ztemp )
{
   return jx_PAMGRelax8(jx_hpCSRMatrixPar(A), f, cf_marker, relax_points, relax_weight, omega, l1_norms, u, Vtemp, Ztemp);
}