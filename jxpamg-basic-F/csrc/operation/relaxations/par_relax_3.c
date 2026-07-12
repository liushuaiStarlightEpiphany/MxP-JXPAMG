//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_3.c
 *  Date: 2011/09/03
 */ 
#include "jxf_hpcsr.h"

/*!
 * \fn JXF_Int jxf_PAMGRelax3
 * \brief hybrid: SOR-J mix off-processor, SOR on-processor
          with outer relaxation parameters (forward solve)
 * \date 2011/09/03
 */
JXF_Int  
jxf_PAMGRelax3( jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   MPI_Comm         comm         = jxf_ParCSRMatrixComm(par_matrix);
   jxf_CSRMatrix    *A_diag       = jxf_ParCSRMatrixDiag(par_matrix);
   JXF_Real          *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int             *A_diag_i     = jxf_CSRMatrixI(A_diag);
   JXF_Int             *A_diag_j     = jxf_CSRMatrixJ(A_diag);
   jxf_CSRMatrix    *A_offd       = jxf_ParCSRMatrixOffd(par_matrix);
   JXF_Int             *A_offd_i     = jxf_CSRMatrixI(A_offd);
   JXF_Real          *A_offd_data  = jxf_CSRMatrixData(A_offd);
   JXF_Int             *A_offd_j     = jxf_CSRMatrixJ(A_offd);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix);
   jxf_ParCSRCommHandle *comm_handle = NULL;

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *f_local = jxf_ParVectorLocalVector(par_rhs);
   JXF_Real         *f_data  = jxf_VectorData(f_local);

   jxf_Vector      *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real         *Vtemp_data  = jxf_VectorData(Vtemp_local);
   JXF_Real 	  *Vext_data  = NULL;
   JXF_Real 	  *v_buf_data = NULL;
   JXF_Real 	  *tmp_data   = NULL;

   JXF_Int             i, j;
   JXF_Int             ii, jj;
   JXF_Int             ns, ne, size, rest;
   JXF_Int             relax_error = 0;
   JXF_Int             num_sends = 0;
   JXF_Int             index, start;
   JXF_Int             num_procs, num_threads, my_id;

   JXF_Real          zero = 0.0;
   JXF_Real          res, res0, res2;
   JXF_Real          one_minus_omega;
   JXF_Real          prod;

   one_minus_omega  = 1.0 - omega;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  
   num_threads = jxf_NumThreads();

  /*-------------------------------------------------------------------------
   * added by peghoty, if the comm_pkg of par_matrix is not created  
   * previously, something will be wrong when the "relax" function  
   * is called on multi-processors occasions.  2009/07/24
   *----------------------------------------------------------------------- */
   
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_matrix);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix); 
   }

      //------------------------------------------------------------------------------//
      //   Hybrid: Jacobi off-processor, Gauss-Seidel on-processor (forward loop)     //
      //------------------------------------------------------------------------------//
      
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
               for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg,i+1); j ++)
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
            }
 
            comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);

            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
         }

        /*-----------------------------------------------------------------
         * Relax all points.
         *-----------------------------------------------------------------*/

         if (relax_weight == 1 && omega == 1)
         {
            if (relax_points == 0)
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)	/* interior points first */
                     {

                       /*--------------------------------------------------------------------
                        * If diagonal is nonzero, relax point i; otherwise, skip it.
                        *-------------------------------------------------------------------*/
             
                        if ( A_diag_data[A_diag_i[i]] != zero)
                        {
                           res = f_data[i];
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res -= A_diag_data[jj] * u_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                  
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] = res / A_diag_data[A_diag_i[i]];
                        }
                     }
                  }
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++) /* interior points first */
                  {

                    /*-----------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *----------------------------------------------------------------*/
             
                     if ( A_diag_data[A_diag_i[i]] != zero)
                     {
                        res = f_data[i];
                        for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                        {
                           ii = A_diag_j[jj];
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                        {
                           ii = A_offd_j[jj];
                           res -= A_offd_data[jj] * Vext_data[ii];
                        }
                        u_data[i] = res / A_diag_data[A_diag_i[i]];
                     }
                  }
               }
            }

           /*-----------------------------------------------------------------
            * Relax only C or F points as determined by relax_points.
            *-----------------------------------------------------------------*/

            else
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)  /* relax interior points */
                     {

                       /*-----------------------------------------------------------
                        * If i is of the right type ( C or F ) and diagonal is
                        * nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero)
                        {
                           res = f_data[i];
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res -= A_diag_data[jj] * u_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] = res / A_diag_data[A_diag_i[i]];
                        }
                     }     
                  }     
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* relax interior points */
                  {

                    /*-----------------------------------------------------------
                     * If i is of the right type ( C or F ) and diagonal is
                     * nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------*/
             
                     if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero)
                     {
                        res = f_data[i];
                        for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                        {
                           ii = A_diag_j[jj];
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                        {
                           ii = A_offd_j[jj];
                           res -= A_offd_data[jj] * Vext_data[ii];
                        }
                        u_data[i] = res / A_diag_data[A_diag_i[i]];
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
            prod = (1.0 - relax_weight*omega);
            
            if (relax_points == 0)
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)  /* interior points first */
                     {

                       /*-----------------------------------------------------------
                        * If diagonal is nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if ( A_diag_data[A_diag_i[i]] != zero)
                        {
                           res = f_data[i];
                           res0 = 0.0;
                           res2 = 0.0;
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res0 -= A_diag_data[jj] * u_data[ii];
                                 res2 += A_diag_data[jj] * Vtemp_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] *= prod;
                           u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) 
                                        / A_diag_data[A_diag_i[i]];
                        }
                     }
                  }
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* interior points first */
                  {

                    /*------------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------------*/
             
                     if ( A_diag_data[A_diag_i[i]] != zero)
                     {
                        res0 = 0.0;
                        res2 = 0.0;
                        res = f_data[i];
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
                        u_data[i] += relax_weight*(omega*res+res0+one_minus_omega*res2)/A_diag_data[A_diag_i[i]];
                     }
                  }
               }
            }

           /*-----------------------------------------------------------------
            * Relax only C or F points as determined by relax_points.
            *-----------------------------------------------------------------*/

            else
           {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     for (i = ns; i < ne; i ++)  /* relax interior points */
                     {

                       /*-----------------------------------------------------------
                        * If i is of the right type ( C or F ) and diagonal is
                        * nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero)
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
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] *= prod;
                           u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) 
                                        / A_diag_data[A_diag_i[i]];
                        }
                     }     
                  }     
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* relax interior points */
                  {

                    /*-----------------------------------------------------------
                     * If i is of the right type ( C or F ) and diagonal is
                     * nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------*/
             
                     if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero)
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
                        u_data[i] += relax_weight*(omega*res+res0+one_minus_omega*res2)/A_diag_data[A_diag_i[i]];
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
jxf_hpPAMGRelax3( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelax3(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}
/*!
 * \fn JXF_Int jxf_PAMGRelaxAI3
 */
JXF_Int  
jxf_PAMGRelaxAI3( jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp )
{
   MPI_Comm         comm         = jxf_ParCSRMatrixComm(par_matrix);
   jxf_CSRMatrix    *A_diag       = jxf_ParCSRMatrixDiag(par_matrix);
   JXF_Real          *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int             *A_diag_i     = jxf_CSRMatrixI(A_diag);
   JXF_Int             *A_diag_j     = jxf_CSRMatrixJ(A_diag);
   jxf_CSRMatrix    *A_offd       = jxf_ParCSRMatrixOffd(par_matrix);
   JXF_Int             *A_offd_i     = jxf_CSRMatrixI(A_offd);
   JXF_Real          *A_offd_data  = jxf_CSRMatrixData(A_offd);
   JXF_Int             *A_offd_j     = jxf_CSRMatrixJ(A_offd);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix);
   jxf_ParCSRCommHandle *comm_handle = NULL;

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *f_local = jxf_ParVectorLocalVector(par_rhs);
   JXF_Real         *f_data  = jxf_VectorData(f_local);

   jxf_Vector      *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real         *Vtemp_data  = jxf_VectorData(Vtemp_local);
   JXF_Real 	  *Vext_data  = NULL;
   JXF_Real 	  *v_buf_data = NULL;
   JXF_Real 	  *tmp_data   = NULL;
   
   JXF_Int             i, j;
   JXF_Int             ii, jj;
   JXF_Int             ns, ne, size, rest;
   JXF_Int             relax_error = 0;
   JXF_Int             num_sends = 0;
   JXF_Int             index, start;
   JXF_Int             num_procs, num_threads, my_id;

   JXF_Real          zero = 0.0;
   JXF_Real          res, res0, res2;
   JXF_Real          one_minus_omega;
   JXF_Real          prod;

   one_minus_omega  = 1.0 - omega;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  
   num_threads = jxf_NumThreads();

  /*-------------------------------------------------------------------------
   * added by peghoty, if the comm_pkg of par_matrix is not created  
   * previously, something will be wrong when the "relax" function  
   * is called on multi-processors occasions.  2009/07/24
   *----------------------------------------------------------------------- */
   
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_matrix);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix); 
   }

   
      //------------------------------------------------------------------------------//
      //   Hybrid: Jacobi off-processor, Gauss-Seidel on-processor (forward loop)     //
      //------------------------------------------------------------------------------//

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
               for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg,i+1); j ++)
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
            }
 
            comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);

            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
         }

        /*-----------------------------------------------------------------
         * Relax all points.
         *-----------------------------------------------------------------*/

         if (relax_weight == 1 && omega == 1)
         {
            if (relax_points == 0)
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)	/* interior points first */
                     {

                       /*--------------------------------------------------------------------
                        * If diagonal is nonzero, relax point i; otherwise, skip it.
                        *-------------------------------------------------------------------*/
             
                        if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                        {
                           res = f_data[i];
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res -= A_diag_data[jj] * u_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                  
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] = res / A_diag_data[A_diag_i[i]];
                        }
                     }
                  }
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++) /* interior points first */
                  {

                    /*-----------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *----------------------------------------------------------------*/
             
                     if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                     {
                        res = f_data[i];
                        for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                        {
                           ii = A_diag_j[jj];
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                        {
                           ii = A_offd_j[jj];
                           res -= A_offd_data[jj] * Vext_data[ii];
                        }
                        u_data[i] = res / A_diag_data[A_diag_i[i]];
                     }
                  }
               }
            }

           /*-----------------------------------------------------------------
            * Relax only C or F points as determined by relax_points.
            *-----------------------------------------------------------------*/

            else
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)  /* relax interior points */
                     {

                       /*-----------------------------------------------------------
                        * If i is of the right type ( C or F ) and diagonal is
                        * nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                        {
                           res = f_data[i];
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res -= A_diag_data[jj] * u_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] = res / A_diag_data[A_diag_i[i]];
                        }
                     }     
                  }     
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* relax interior points */
                  {

                    /*-----------------------------------------------------------
                     * If i is of the right type ( C or F ) and diagonal is
                     * nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------*/
             
                     if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                     {
                        res = f_data[i];
                        for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                        {
                           ii = A_diag_j[jj];
                           res -= A_diag_data[jj] * u_data[ii];
                        }
                        for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                        {
                           ii = A_offd_j[jj];
                           res -= A_offd_data[jj] * Vext_data[ii];
                        }
                        u_data[i] = res / A_diag_data[A_diag_i[i]];
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
            prod = (1.0 - relax_weight*omega);
            
            if (relax_points == 0)
            {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     
                     for (i = ns; i < ne; i ++)  /* interior points first */
                     {

                       /*-----------------------------------------------------------
                        * If diagonal is nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                        {
                           res = f_data[i];
                           res0 = 0.0;
                           res2 = 0.0;
                           for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
                           {
                              ii = A_diag_j[jj];
                              if (ii >= ns && ii < ne)
                              {
                                 res0 -= A_diag_data[jj] * u_data[ii];
                                 res2 += A_diag_data[jj] * Vtemp_data[ii];
                              }
                              else
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] *= prod;
                           u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) 
                                        / A_diag_data[A_diag_i[i]];
                        }
                     }
                  }
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* interior points first */
                  {

                    /*------------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------------*/
             
                     if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
                     {
                        res0 = 0.0;
                        res2 = 0.0;
                        res = f_data[i];
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
                        u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
                     }
                  }
               }
            }

           /*-----------------------------------------------------------------
            * Relax only C or F points as determined by relax_points.
            *-----------------------------------------------------------------*/

            else
           {
               if (num_threads > 1)
               {
                  tmp_data = jxf_CTAlloc(JXF_Real,n);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
                  for (i = 0; i < n; i ++)
                  {
                     tmp_data[i] = u_data[i];
                  }
#define JXF_SMP_PRIVATE i,ii,j,jj,ns,ne,res,rest,size
#include "../../../include/jxf_smp_forloop.h"
                  for (j = 0; j < num_threads; j ++)
                  {
                     size = n/num_threads;
                     rest = n - size*num_threads;
                     if (j < rest)
                     {
                        ns = j*size+j;
                        ne = (j+1)*size+j+1;
                     }
                     else
                     {
                        ns = j*size+rest;
                        ne = (j+1)*size+rest;
                     }
                     for (i = ns; i < ne; i ++)  /* relax interior points */
                     {

                       /*-----------------------------------------------------------
                        * If i is of the right type ( C or F ) and diagonal is
                        * nonzero, relax point i; otherwise, skip it.
                        *-----------------------------------------------------------*/
             
                        if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
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
                              {
                                 res -= A_diag_data[jj] * tmp_data[ii];
                              }
                           }
                           for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
                           {
                              ii = A_offd_j[jj];
                              res -= A_offd_data[jj] * Vext_data[ii];
                           }
                           u_data[i] *= prod;
                           u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) 
                                        / A_diag_data[A_diag_i[i]];
                        }
                        //}
                     }     
                  }     
                  jxf_TFree(tmp_data);
               }
               else
               {
                  for (i = 0; i < n; i ++)  /* relax interior points */
                  {

                    /*-----------------------------------------------------------
                     * If i is of the right type ( C or F ) and diagonal is
                     * nonzero, relax point i; otherwise, skip it.
                     *-----------------------------------------------------------*/
             
                     if (cf_marker[i] == relax_points && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
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
                        u_data[i] += relax_weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
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
jxf_hpPAMGRelaxAI3( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelaxAI3(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}