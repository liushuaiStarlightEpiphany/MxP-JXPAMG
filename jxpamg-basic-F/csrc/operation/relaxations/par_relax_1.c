//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_1.c
 *  Date: 2011/09/03
 */ 

#include "jxf_hpcsr.h"

/*!
 * \fn JXF_Int jxf_PAMGRelax1
 * \brief Gauss-Seidel <--- very slow, sequential
 * \date 2011/09/03
 */
JXF_Int  
jxf_PAMGRelax1( jxf_ParCSRMatrix *par_matrix,
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

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *f_local = jxf_ParVectorLocalVector(par_rhs);
   JXF_Real         *f_data  = jxf_VectorData(f_local);

   JXF_Real 	  *Vext_data  = NULL;
   JXF_Real 	  *v_buf_data = NULL;

   JXF_Int             i, j, jr;
   JXF_Int             ii, jj;
   JXF_Int             relax_error = 0;
   JXF_Int             num_sends = 0;
   JXF_Int             num_recvs = 0;
   JXF_Int             num_procs, my_id, ip, p;
   JXF_Int             vec_start, vec_len;
   MPI_Status     *status = NULL;
   MPI_Request    *requests = NULL;

   JXF_Real          zero = 0.0;
   JXF_Real          res;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  

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
      //                             Gauss-Seidel VERY SLOW                           //
      //------------------------------------------------------------------------------//

         if (num_procs > 1)
         {
            num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
            num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);

            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

            Vext_data = jxf_CTAlloc(JXF_Real,num_cols_offd);
        
            status   = jxf_CTAlloc(MPI_Status,num_recvs+num_sends);
            requests = jxf_CTAlloc(MPI_Request, num_recvs+num_sends);

            if (num_cols_offd)
            {
               A_offd_j = jxf_CSRMatrixJ(A_offd);
               A_offd_data = jxf_CSRMatrixData(A_offd);
            }
         } 

        /*-----------------------------------------------------------------
         * Relax all points.
         *-----------------------------------------------------------------*/
         for (p = 0; p < num_procs; p ++)
         {
            jr = 0;
	
            if (p != my_id)
            {
               for (i = 0; i < num_sends; i ++)
               {
                  ip = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
                  if (ip == p)
                  {
                     vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
                     vec_len = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
                     for (j = vec_start; j < vec_start+vec_len; j ++)
                     {
                        v_buf_data[j] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
                     }
                     jxf_MPI_Isend(&v_buf_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
               }
               jxf_MPI_Waitall(jr,requests,status);
               jxf_MPI_Barrier(comm);
            }
            else
            {
               if (num_procs > 1)
               {
                  for (i = 0; i < num_recvs; i ++)
                  {
                     ip = jxf_ParCSRCommPkgRecvProc(comm_pkg, i);
                     vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
                     vec_len = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1)-vec_start;
                     jxf_MPI_Irecv(&Vext_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
                  jxf_MPI_Waitall(jr, requests, status);
               }
               
               if (relax_points == 0)
               {
               
                  for (i = 0; i < n; i ++)	
                  {

                    /*--------------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *------------------------------------------------------------------*/
             
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
                  } // end for i
                  
               }

              /*-----------------------------------------------------------------
               * Relax only C or F points as determined by relax_points.
               *-----------------------------------------------------------------*/

               else
               {
                  for (i = 0; i < n; i ++)
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
                  }  // end for i   
               }
               
               if (num_procs > 1)
               {
                  jxf_MPI_Barrier(comm);
               }
            }
            
         } // end for p
         
         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
            jxf_TFree(status);
            jxf_TFree(requests);
         }

   return(relax_error); 
}

JXF_Int  
jxf_hpPAMGRelax1( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelax1(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}

/*!
 * \fn JXF_Int jxf_PAMGRelax1
 */
JXF_Int  
jxf_PAMGRelaxAI1( jxf_ParCSRMatrix *par_matrix,
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

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *f_local = jxf_ParVectorLocalVector(par_rhs);
   JXF_Real         *f_data  = jxf_VectorData(f_local);

   JXF_Real 	  *Vext_data  = NULL;
   JXF_Real 	  *v_buf_data = NULL;
   
   JXF_Int             i, j, jr;
   JXF_Int             ii, jj;
   JXF_Int             relax_error = 0;
   JXF_Int             num_sends = 0;
   JXF_Int             num_recvs = 0;
   JXF_Int             num_procs, my_id, ip, p;
   JXF_Int             vec_start, vec_len;
   MPI_Status     *status = NULL;
   MPI_Request    *requests = NULL;

   JXF_Real          zero = 0.0;
   JXF_Real          res;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  

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
      //                             Gauss-Seidel VERY SLOW                           //
      //------------------------------------------------------------------------------//
      
         if (num_procs > 1)
         {
            num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
            num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);

            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));

            Vext_data = jxf_CTAlloc(JXF_Real,num_cols_offd);
        
            status   = jxf_CTAlloc(MPI_Status,num_recvs+num_sends);
            requests = jxf_CTAlloc(MPI_Request, num_recvs+num_sends);

            if (num_cols_offd)
            {
               A_offd_j = jxf_CSRMatrixJ(A_offd);
               A_offd_data = jxf_CSRMatrixData(A_offd);
            }
         } 

        /*-----------------------------------------------------------------
         * Relax all points.
         *-----------------------------------------------------------------*/
         for (p = 0; p < num_procs; p ++)
         {
            jr = 0;
	
            if (p != my_id)
            {
               for (i = 0; i < num_sends; i ++)
               {
                  ip = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
                  if (ip == p)
                  {
                     vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
                     vec_len = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
                     for (j = vec_start; j < vec_start+vec_len; j ++)
                     {
                        v_buf_data[j] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
                     }
                     jxf_MPI_Isend(&v_buf_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
               }
               jxf_MPI_Waitall(jr,requests,status);
               jxf_MPI_Barrier(comm);
            }
            else
            {
               if (num_procs > 1)
               {
                  for (i = 0; i < num_recvs; i ++)
                  {
                     ip = jxf_ParCSRCommPkgRecvProc(comm_pkg, i);
                     vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i);
                     vec_len = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1)-vec_start;
                     jxf_MPI_Irecv(&Vext_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
                  jxf_MPI_Waitall(jr, requests, status);
               }
               
               if (relax_points == 0)
               {
               
                  for (i = 0; i < n; i ++)	
                  {

                    /*--------------------------------------------------------------------
                     * If diagonal is nonzero, relax point i; otherwise, skip it.
                     *------------------------------------------------------------------*/
             
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
                  } // end for i
                  
               }

              /*-----------------------------------------------------------------
               * Relax only C or F points as determined by relax_points.
               *-----------------------------------------------------------------*/

               else
               {
                  for (i = 0; i < n; i ++)
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
                  }  // end for i   
               }
               
               if (num_procs > 1)
               {
                  jxf_MPI_Barrier(comm);
               }
            }
            
         } // end for p
         
         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
            jxf_TFree(status);
            jxf_TFree(requests);
         }

   return(relax_error); 
}

JXF_Int  
jxf_hpPAMGRelaxAI1( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelaxAI1(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}