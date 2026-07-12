//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_9.c
 *  Date: 2011/09/03
 */ 

#include "jxf_hpcsr.h"

/*!
 * \fn JXF_Int jxf_PAMGRelax9
 * \brief Direct Solve
 * \date 2011/09/03
 */
JXF_Int  
jxf_PAMGRelax9( jxf_ParCSRMatrix *par_matrix,
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
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix);

   JXF_Int             n_global = jxf_ParCSRMatrixGlobalNumRows(par_matrix);
   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int             first_index   = jxf_ParVectorFirstIndex(par_app);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_CSRMatrix    *A_CSR;
   JXF_Int             *A_CSR_i;   
   JXF_Int             *A_CSR_j;
   JXF_Real          *A_CSR_data;
   
   jxf_Vector       *f_vector;
   JXF_Real          *f_vector_data;

   JXF_Int             i;
   JXF_Int             jj;
   JXF_Int             column;
   JXF_Int             relax_error = 0;
   JXF_Int             num_procs, my_id;

   JXF_Real         *A_mat;
   JXF_Real         *b_vec;

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
      //                   Direct solve: use gaussian elimination                     //
      //------------------------------------------------------------------------------//

        /*-----------------------------------------------------------------
         *  Generate CSR matrix from ParCSRMatrix par_matrix
         *-----------------------------------------------------------------*/
#ifdef JXF_NO_GLOBAL_PARTITION
         /* all processors are needed for these routines */
         A_CSR = jxf_ParCSRMatrixToCSRMatrixAll(par_matrix);
         f_vector = jxf_ParVectorToVectorAll(par_rhs);
	 if (n)
	 {
#else
	 if (n)
	 {
	    A_CSR = jxf_ParCSRMatrixToCSRMatrixAll(par_matrix);
	    f_vector = jxf_ParVectorToVectorAll(par_rhs);
#endif
 	    A_CSR_i = jxf_CSRMatrixI(A_CSR);
 	    A_CSR_j = jxf_CSRMatrixJ(A_CSR);
 	    A_CSR_data = jxf_CSRMatrixData(A_CSR);
   	    f_vector_data = jxf_VectorData(f_vector);
            
            A_mat = jxf_CTAlloc(JXF_Real, n_global*n_global);
            
            b_vec = jxf_CTAlloc(JXF_Real, n_global);    

           /*---------------------------------------------------------------
            *  Load CSR matrix into A_mat.
            *---------------------------------------------------------------*/

            for (i = 0; i < n_global; i ++)
            {
               for (jj = A_CSR_i[i]; jj < A_CSR_i[i+1]; jj ++)
               {
                  column = A_CSR_j[jj];
                  A_mat[i*n_global+column] = A_CSR_data[jj];
               }
               b_vec[i] = f_vector_data[i];
            }

            relax_error = jxf_gselim(A_mat,b_vec,n_global);
            /* use version with pivoting */
            /* relax_error = jxf_gselim_piv(A_mat,b_vec,n_global);*/

            for (i = 0; i < n; i ++)
            {
               u_data[i] = b_vec[first_index+i];
            }

	    jxf_TFree(A_mat); 
            jxf_TFree(b_vec);
            jxf_CSRMatrixDestroy(A_CSR);
            A_CSR = NULL;
            jxf_SeqVectorDestroy(f_vector);
            f_vector = NULL;
         
         }
#ifdef JXF_NO_GLOBAL_PARTITION
         else
         {
            jxf_CSRMatrixDestroy(A_CSR);
            A_CSR = NULL;
            jxf_SeqVectorDestroy(f_vector);
            f_vector = NULL;
         }
#endif

   return(relax_error); 
}

JXF_Int  
jxf_hpPAMGRelax9( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelax9(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}