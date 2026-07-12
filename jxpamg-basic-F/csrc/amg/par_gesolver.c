//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_gesolver.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_GaussElimination
 * \brief Gauss Elimination solver.
 * \author peghoty
 * \date 2009/07/25
 */
JXF_Int 
jxf_GaussElimination( JXF_Int              solverid,
                     jxf_ParCSRMatrix *par_matrix, 
                     jxf_ParVector    *par_rhs,
                     jxf_ParVector    *par_app )
{
   jxf_CSRMatrix    *A_diag      = jxf_ParCSRMatrixDiag(par_matrix);
   JXF_Int              n_global    = jxf_ParCSRMatrixGlobalNumRows(par_matrix);
   JXF_Int              n           = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int              first_index = jxf_ParVectorFirstIndex(par_app);

   jxf_CSRMatrix    *A_CSR;
   JXF_Int             *A_CSR_i;  
   JXF_Int             *A_CSR_j;
   JXF_Real          *A_CSR_data;
   
   jxf_Vector       *f_vector;
   JXF_Real          *f_vector_data;

   jxf_Vector       *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real          *u_data  = jxf_VectorData(u_local);

   JXF_Real          *A_mat;
   JXF_Real          *b_vec;

   JXF_Int              err_flag = 0;
   JXF_Int              i, jj;
   JXF_Int              column;

   if (n)   
   {
      A_CSR    = jxf_ParCSRMatrixToCSRMatrixAll(par_matrix);
      f_vector = jxf_ParVectorToVectorAll(par_rhs);

      A_CSR_i = jxf_CSRMatrixI(A_CSR);
      A_CSR_j = jxf_CSRMatrixJ(A_CSR);
      A_CSR_data    = jxf_CSRMatrixData(A_CSR);
      f_vector_data = jxf_VectorData(f_vector);

      A_mat = jxf_CTAlloc(JXF_Real, n_global*n_global);
      b_vec = jxf_CTAlloc(JXF_Real, n_global);    

     /*---------------------------------------------------------------
      *  Load CSR matrix into A_mat and b_vec.
      *--------------------------------------------------------------*/

      for (i = 0; i < n_global; i ++)
      {
         for (jj = A_CSR_i[i]; jj < A_CSR_i[i+1]; jj ++)
         {
            column = A_CSR_j[jj];
            A_mat[i*n_global+column] = A_CSR_data[jj];
         }
         b_vec[i] = f_vector_data[i];
      }

      if (solverid == 9)
      {
         err_flag = jxf_gselim(A_mat, b_vec, n_global);
      }
      if (solverid == 10)  /* newly-added 2009/08/08 */
      {
         err_flag = jxf_gselim_piv(A_mat, b_vec, n_global);
      }

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
   
   return err_flag; 
}
