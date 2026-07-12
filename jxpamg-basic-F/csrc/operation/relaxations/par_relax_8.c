//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_8.c
 *  Date: 2024/01/31
 */ 

#include "jxf_hpcsr.h"

JXF_Int
jxf_PAMGRelax8( jxf_ParCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp )
{ 
   MPI_Comm	comm = jxf_ParCSRMatrixComm(A);
   jxf_CSRMatrix *A_diag = jxf_ParCSRMatrixDiag(A);
   JXF_Real *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int *A_diag_i = jxf_CSRMatrixI(A_diag);
   JXF_Int *A_diag_j = jxf_CSRMatrixJ(A_diag);
   jxf_CSRMatrix *A_offd = jxf_ParCSRMatrixOffd(A);
   JXF_Int *A_offd_i = jxf_CSRMatrixI(A_offd);
   JXF_Real *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int *A_offd_j = jxf_CSRMatrixJ(A_offd);
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_ParCSRCommHandle *comm_handle;

   JXF_Int n = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   jxf_Vector *u_local = jxf_ParVectorLocalVector(u);
   JXF_Real *u_data  = jxf_VectorData(u_local);

   jxf_Vector *f_local = jxf_ParVectorLocalVector(f);
   JXF_Real *f_data  = jxf_VectorData(f_local);

   jxf_Vector *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real *Vtemp_data = jxf_VectorData(Vtemp_local);
   JXF_Real *Vext_data = NULL;
   JXF_Real *v_buf_data = NULL;
   JXF_Real *tmp_data;

   jxf_Vector *Ztemp_local;
   JXF_Real *Ztemp_data;

   JXF_Int i, j;
   JXF_Int ii, jj;
   JXF_Int ns, ne, size, rest;
   JXF_Int relax_error = 0;
   JXF_Int num_sends;
   JXF_Int index, start;
   JXF_Int num_procs, num_threads, my_id;

   JXF_Real zero = 0.0;
   JXF_Real res, res0, res2;
   JXF_Real one_minus_omega;
   JXF_Real prod;

   one_minus_omega = 1.0 - omega;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   num_threads = jxf_NumThreads();

   if (num_threads > 1)
   {
      Ztemp_local = jxf_ParVectorLocalVector(Ztemp);
      Ztemp_data = jxf_VectorData(Ztemp_local);
   }

   if (num_procs > 1)
   {
      num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
      v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
      Vext_data = jxf_CTAlloc(JXF_Real, num_cols_offd);

      if (num_cols_offd)
      {
         A_offd_j = jxf_CSRMatrixJ(A_offd);
         A_offd_data = jxf_CSRMatrixData(A_offd);
      }

      index = 0;
      for (i = 0; i < num_sends; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1); j ++)
            v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
      }

      comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);

      jxf_ParCSRCommHandleDestroy(comm_handle);
      comm_handle = NULL;
   }

   if (relax_weight == 1 && omega == 1)
   {
      if (relax_points == 0)
      {
         if (num_threads > 1)
         {
            tmp_data = Ztemp_data;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
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
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
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
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
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
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
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
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++)
               tmp_data[i] = u_data[i];
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
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
      jxf_TFree(Vext_data);
      jxf_TFree(v_buf_data);
   }

   return(relax_error);
}

JXF_Int
jxf_hpPAMGRelax8( jxf_hpCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp )
{
   return jxf_PAMGRelax8(jxf_hpCSRMatrixPar(A), f, cf_marker, relax_points, relax_weight, omega, l1_norms, u, Vtemp, Ztemp);
}