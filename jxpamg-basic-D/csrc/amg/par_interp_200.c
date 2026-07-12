//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_interp_200.c
 *  Date: 2022/04/11
 */ 

#include "jx_pamg.h"

/*!
 * \fn JX_Int jx_PAMGBuildInterpDinvRATrans
 * \date 2022/04/11
 */
JX_Int
jx_PAMGBuildInterpDinvRATrans( jx_ParCSRMatrix  *A,
                               jx_ParCSRMatrix  *R,
                               jx_ParCSRMatrix **P_ptr )
{
   MPI_Comm              comm     = jx_ParCSRMatrixComm(A);

   jx_CSRMatrix         *A_diag      = jx_ParCSRMatrixDiag(A);
   JX_Real              *A_diag_data = jx_CSRMatrixData(A_diag);
   JX_Int               *A_diag_i    = jx_CSRMatrixI(A_diag);

   /* Interpolation matrix P */
   jx_ParCSRMatrix      *AR;
   /* csr's */
   jx_CSRMatrix    *AR_diag, *P_diag;
   /* arrays */
   JX_Real         *P_diag_data;
   JX_Int          *P_diag_i;
   /* local size */
   JX_Int           n_fine = jx_CSRMatrixNumRows(A_diag);
   JX_Int           my_id;

   JX_Int           i, j;
   JX_Real          factor;

   jx_MPI_Comm_rank(comm, &my_id);

   AR = jx_ParMatmul(R, A);

   if (my_id == 0)
   {
      AR_diag = jx_ParCSRMatrixDiag(AR);

      jx_CSRMatrixTranspose(AR_diag, &P_diag, 1);

      P_diag_data = jx_CSRMatrixData(P_diag);
      P_diag_i    = jx_CSRMatrixI(P_diag);

      for (i = 0; i < n_fine; i++)
      {
         factor = 1.0 / A_diag_data[A_diag_i[i]];
         for (j = P_diag_i[i]; j < P_diag_i[i+1]; j ++)
         {
            P_diag_data[j] *= factor;
         }
      }
   }

   *P_ptr = jx_CSRMatrixToParCSRMatrix(comm, P_diag, NULL, NULL);

   jx_ParCSRMatrixDestroy(AR);

   if (my_id == 0)
   {
      jx_CSRMatrixDestroy(P_diag);
   }

   return 0;
}
