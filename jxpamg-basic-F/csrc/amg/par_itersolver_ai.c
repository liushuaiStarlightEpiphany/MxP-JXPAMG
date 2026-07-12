//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_itersolver.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_CoarsestIterativeMethod
 * \brief One iteration for some iterative method.
 * \date 2011/09/03
 */
JXF_Int  
jxf_CoarsestIterativeMethodAI( jxf_hpCSRMatrix *hp_matrix,
                            jxf_ParVector    *par_rhs,
                            JXF_Int             *relax_marker,
                            JXF_Int              solverid,
                            JXF_Real           weight,
                            JXF_Real           omega,
                            jxf_ParVector    *par_app,
                            jxf_ParVector    *Vtemp )
{
   MPI_Comm             comm        = jxf_hpCSRMatrixComm(hp_matrix);
   jxf_CSRMatrix        *A_diag      = jxf_hpCSRMatrixDiag(hp_matrix);
   JXF_Real              *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int                 *A_diag_i    = jxf_CSRMatrixI(A_diag);
   JXF_Int                 *A_diag_j    = jxf_CSRMatrixJ(A_diag);
   jxf_CSRMatrix        *A_offd      = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int                 *A_offd_i    = jxf_CSRMatrixI(A_offd);
   JXF_Real              *A_offd_data = jxf_CSRMatrixData(A_offd);
   JXF_Int                 *A_offd_j    = jxf_CSRMatrixJ(A_offd);
   jxf_ParCSRCommPkg    *comm_pkg    = jxf_hpCSRMatrixCommPkg(hp_matrix);
   jxf_ParCSRCommHandle *comm_handle = NULL;

   JXF_Int             n             = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             num_cols_offd = jxf_CSRMatrixNumCols(A_offd);
 
   jxf_Vector   *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real      *u_data  = jxf_VectorData(u_local);

   jxf_Vector   *f_local = jxf_ParVectorLocalVector(par_rhs);
   JXF_Real      *f_data  = jxf_VectorData(f_local);

   jxf_Vector   *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real      *Vtemp_data  = jxf_VectorData(Vtemp_local);
   JXF_Real      *Vext_data   = NULL;
   JXF_Real      *v_buf_data  = NULL;
   
   JXF_Int             i, j, jr;
   JXF_Int             ii, jj;
   JXF_Int             num_sends = 0;
   JXF_Int             num_recvs = 0;
   JXF_Int             index, start;
   JXF_Int             num_procs, my_id, ip, p;
   JXF_Int             vec_start, vec_len;
   MPI_Status     *status   = NULL;
   MPI_Request    *requests = NULL;

   JXF_Real          zero = 0.0;
   JXF_Real          res, res0, res2;
   JXF_Real          one_minus_weight;
   JXF_Real          one_minus_omega;
   JXF_Real          prod;

   one_minus_weight = 1.0 - weight;
   one_minus_omega  = 1.0 - omega;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  

  /*-------------------------------------------------------------------------
   * added by peghoty, if the comm_pkg of par_matrix is not created  
   * previously, something will be wrong when the "relax" function  
   * is called on multi-processors occasions.  2009/07/24
   *----------------------------------------------------------------------- */
   
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(jxf_hpCSRMatrixPar(hp_matrix));
      comm_pkg = jxf_hpCSRMatrixCommPkg(hp_matrix); 
   }


  /*-----------------------------------------------------------------------------
   * Switch statement to direct control based on relax_type:
   *     relax_type = 0 -> Jacobi or CF-Jacobi
   *     relax_type = 1 -> Gauss-Seidel <--- very slow, sequential
   *     relax_type = 2 -> Gauss_Seidel: interior points in parallel ,
   *			 	   	  boundary sequential 
   *     relax_type = 3 -> hybrid: SOR-J mix off-processor, SOR on-processor
   *     		    with outer relaxation parameters (forward solve)
   *     relax_type = 4 -> hybrid: SOR-J mix off-processor, SOR on-processor
   *     		    with outer relaxation parameters (backward solve)
   *     relax_type = 5 -> hybrid: GS-J mix off-processor, chaotic GS on-node
   *     relax_type = 6 -> hybrid: SSOR-J mix off-processor, SSOR on-processor
   *     		    with outer relaxation parameters
   *-----------------------------------------------------------------------------*/
   
   switch (solverid)
   {        

      //------------------------------------------------------------------------------//
      //-----------------------------   Weighted Jacobi  -----------------------------//
      //------------------------------------------------------------------------------//
      
      case 0: 
      {
	 if (num_procs > 1)
	 {
   	    num_sends  = jxf_ParCSRCommPkgNumSends(comm_pkg);
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
   	       {
   	          v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
   	       }
   	    }
 
   	    comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);
	 }

        /*-----------------------------------------------------------------
         * Copy current approximation into temporary vector.
         *-----------------------------------------------------------------*/
#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
         for (i = 0; i < n; i ++)
         {
            Vtemp_data[i] = u_data[i];
         }
         
	 if (num_procs > 1)
	 { 
            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
	 } 

#define JXF_SMP_PRIVATE i,ii,jj,res
#include "../../include/jxf_smp_forloop.h"
         for (i = 0; i < n; i ++)
         {
           /*----------------------------------------------------------------
            * If diagonal is nonzero, relax point i; otherwise, skip it.
            *---------------------------------------------------------------*/
            if (A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
            {
               res = f_data[i];
               for (jj = A_diag_i[i]+1; jj < A_diag_i[i+1]; jj ++)
               {
                  ii = A_diag_j[jj];
                  res -= A_diag_data[jj] * Vtemp_data[ii];
               }
               for (jj = A_offd_i[i]; jj < A_offd_i[i+1]; jj ++)
               {
                  ii = A_offd_j[jj];
                  res -= A_offd_data[jj] * Vext_data[ii];
               }
               u_data[i] *= one_minus_weight; 
               u_data[i] += weight * res / A_diag_data[A_diag_i[i]];
            }
         }

         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
         }
      }
      break;


      //------------------------------------------------------------------------------//
      //                             Gauss-Seidel VERY SLOW                           //
      //------------------------------------------------------------------------------//

      case 1:  
      {
         if (num_procs > 1)
         {
   	    num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   	    num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);

   	    v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   	    Vext_data  = jxf_CTAlloc(JXF_Real, num_cols_offd);
        
   	    status   = jxf_CTAlloc(MPI_Status, num_recvs + num_sends);
   	    requests = jxf_CTAlloc(MPI_Request, num_recvs + num_sends);

   	    if (num_cols_offd)
   	    {
   	       A_offd_j = jxf_CSRMatrixJ(A_offd);
   	       A_offd_data = jxf_CSRMatrixData(A_offd);
   	    }
   	 } 

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
                     vec_len   = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
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
                    vec_len = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
                    jxf_MPI_Irecv(&Vext_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
                  jxf_MPI_Waitall(jr, requests, status);
               }

               for (i = 0; i < n; i ++)	
               {
                  //if ( A_diag_data[A_diag_i[i]] != zero)
                  if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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

               if (num_procs > 1) jxf_MPI_Barrier(comm);
	    }
	 } // end for

	 if (num_procs > 1)
	 {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
            jxf_TFree(status);
            jxf_TFree(requests);
	 }
      }
      break;
      

      //------------------------------------------------------------------------------//
      //   Gauss-Seidel: relax interior points in parallel, boundary sequentially     //
      //------------------------------------------------------------------------------//
       
      case 2: 
      {
	 if (num_procs > 1)
	 {
            num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
            num_recvs = jxf_ParCSRCommPkgNumRecvs(comm_pkg);
            
            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
            Vext_data = jxf_CTAlloc(JXF_Real, num_cols_offd);
        
            status  = jxf_CTAlloc(MPI_Status, num_recvs + num_sends);
            requests= jxf_CTAlloc(MPI_Request, num_recvs + num_sends);

            if (num_cols_offd)
            {
		A_offd_j = jxf_CSRMatrixJ(A_offd);
		A_offd_data = jxf_CSRMatrixData(A_offd);
            }
         }

        /*-----------------------------------------------------------------
         * Relax interior points first
         *-----------------------------------------------------------------*/
      
         for (i = 0; i < n; i ++)	
         {
            if ((A_offd_i[i+1]-A_offd_i[i]) == zero && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
            {
               res = f_data[i];
               for (jj = A_diag_i[i] + 1; jj < A_diag_i[i+1]; jj ++)
               {
                  ii = A_diag_j[jj];
                  res -= A_diag_data[jj] * u_data[ii];
               }
               u_data[i] = res / A_diag_data[A_diag_i[i]];
            }
         }
          

         for (p = 0; p < num_procs; p++)
         {
            jr = 0;
            if (p != my_id)
            {
               for (i = 0; i < num_sends; i++)
               {
                  ip = jxf_ParCSRCommPkgSendProc(comm_pkg, i);
                  if (ip == p)
                  {
                     vec_start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
                     vec_len = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i+1) - vec_start;
                     for (j = vec_start; j < vec_start + vec_len; j ++)
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
                     vec_start = jxf_ParCSRCommPkgRecvVecStart(comm_pkg, i);
                     vec_len = jxf_ParCSRCommPkgRecvVecStart(comm_pkg,i+1) - vec_start;
                     jxf_MPI_Irecv(&Vext_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[jr++]);
                  }
                  jxf_MPI_Waitall(jr, requests, status);
               }

               for (i = 0; i < n; i++)	
               {
                  if ((A_offd_i[i+1]-A_offd_i[i]) != zero && A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
         
               if (num_procs > 1) 
               {
                  jxf_MPI_Barrier(comm);
               }
            }
         }

         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
            jxf_TFree(status);
            jxf_TFree(requests);
         }
      }
      break;


      //------------------------------------------------------------------------------//
      //   Hybrid: Jacobi off-processor, Gauss-Seidel on-processor (forward loop)     //
      //------------------------------------------------------------------------------//

      case 3: 
      {
         if (num_procs > 1)
         {
            num_sends  = jxf_ParCSRCommPkgNumSends(comm_pkg);
            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
            Vext_data  = jxf_CTAlloc(JXF_Real, num_cols_offd);

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
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
            }
 
            comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);
            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
         }

         if (weight == 1 && omega == 1)
         {

            for (i = 0; i < n; i ++)   /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
         else
         {
#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++) 
            {
               Vtemp_data[i] = u_data[i];
            }
            
            prod = 1.0 - weight*omega;

            for (i = 0; i < n; i ++)  /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
                  u_data[i] += weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
               }
            }
         }

         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
         }

      }
      break;


      //------------------------------------------------------------------------------//
      //  Hybrid: Jacobi off-processor, Gauss-Seidel/SOR on-processor (backward loop) //
      //------------------------------------------------------------------------------//

      case 4: 
      {
         if (num_procs > 1)
         {
            num_sends  = jxf_ParCSRCommPkgNumSends(comm_pkg);
            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
            Vext_data  = jxf_CTAlloc(JXF_Real, num_cols_offd);
        
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
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
            }
 
            comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);
            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
         }

         if (weight == 1 && omega == 1)
         {
            for (i = n-1; i > -1; i --)  /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
         else
         {
#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++) 
            {
               Vtemp_data[i] = u_data[i];
            }
            
            prod = 1.0 - weight*omega;

            for (i = n-1; i > -1; i --)   /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
                  u_data[i] += weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
               }
            }
         }

         if (num_procs > 1)
         {
	    jxf_TFree(Vext_data);
	    jxf_TFree(v_buf_data);
         }
      }
      break;


      //------------------------------------------------------------------------------//
      //       Hybrid: Jacobi off-processor, chaotic Gauss-Seidel on-processor        //
      //------------------------------------------------------------------------------//

      case 5:
      {
	  if (num_procs > 1)
	  {
   	     num_sends  = jxf_ParCSRCommPkgNumSends(comm_pkg);
   	     v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
	     Vext_data  = jxf_CTAlloc(JXF_Real, num_cols_offd);
        
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
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
   	    }
 
   	    comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);
   	    jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
	 }

#define JXF_SMP_PRIVATE i,ii,jj,res
#include "../../include/jxf_smp_forloop.h"
         for (i = 0; i < n; i ++)	  /* interior points first */
         {
            if (A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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

         if (num_procs > 1)
         {
	    jxf_TFree(Vext_data);
	    jxf_TFree(v_buf_data);
         }

      }
      break;


      //------------------------------------------------------------------------------//
      //  Hybrid: Jacobi off-processor, Symm. Gauss-Seidel/ SSOR on-processor         //
      //  with outer relaxation parameter                                             //
      //------------------------------------------------------------------------------//
      
      case 6:
      {
         if (num_procs > 1)
         {
            num_sends  = jxf_ParCSRCommPkgNumSends(comm_pkg);
            v_buf_data = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
            Vext_data  = jxf_CTAlloc(JXF_Real, num_cols_offd);
        
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
               {
                  v_buf_data[index++] = u_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg,j)];
               }
            }
 
            comm_handle = jxf_ParCSRCommHandleCreate(1, comm_pkg, v_buf_data, Vext_data);
            jxf_ParCSRCommHandleDestroy(comm_handle);
            comm_handle = NULL;
	 }

         if (weight == 1 && omega == 1)
         {

            for (i = 0; i < n; i ++)   /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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

            for (i = n-1; i > -1; i --)   /* interior points first */
            {
               if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
         else
         {
#define JXF_SMP_PRIVATE i
#include "../../include/jxf_smp_forloop.h"
            for (i = 0; i < n; i ++) 
            {
               Vtemp_data[i] = u_data[i];
            }
            
            prod = 1.0 - weight*omega;

            for (i = 0; i < n; i ++)   /* interior points first */
            {
                if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
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
                   u_data[i] += weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
                }
            }

            for (i = n-1; i > -1; i --)   /* interior points first */
            {
                if ( A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1)
                {
                   res0 = 0.0;
                   res  = f_data[i];
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
                   u_data[i] += weight*(omega*res + res0 + one_minus_omega*res2) / A_diag_data[A_diag_i[i]];
                }
            }
         }

         if (num_procs > 1)
         {
            jxf_TFree(Vext_data);
            jxf_TFree(v_buf_data);
         }

      }
      break;

   } // end switch

   return(0); 
}
