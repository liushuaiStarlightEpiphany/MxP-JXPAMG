//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_interp_200.c
 *  Date: 2022/04/11
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGBuildInterpDinvRATrans
 * \date 2022/04/11
 */
JXF_Int
jxf_PAMGBuildInterpDinvRATrans( jxf_ParCSRMatrix  *A,
                               jxf_ParCSRMatrix  *R,
                               jxf_ParCSRMatrix **P_ptr )
{
   MPI_Comm              comm     = jxf_ParCSRMatrixComm(A);

   jxf_CSRMatrix         *A_diag      = jxf_ParCSRMatrixDiag(A);
   JXF_Real              *A_diag_data = jxf_CSRMatrixData(A_diag);
   JXF_Int               *A_diag_i    = jxf_CSRMatrixI(A_diag);

   /* Interpolation matrix P */
   jxf_ParCSRMatrix      *AR;
   /* csr's */
   jxf_CSRMatrix    *AR_diag, *P_diag;
   /* arrays */
   JXF_Real         *P_diag_data;
   JXF_Int          *P_diag_i;
   /* local size */
   JXF_Int           n_fine = jxf_CSRMatrixNumRows(A_diag);
   JXF_Int           my_id;

   JXF_Int           i, j;
   JXF_Real          factor;

   jxf_MPI_Comm_rank(comm, &my_id);

   AR = jxf_ParMatmul(R, A);

   if (my_id == 0)
   {
      AR_diag = jxf_ParCSRMatrixDiag(AR);

      jxf_CSRMatrixTranspose(AR_diag, &P_diag, 1);

      P_diag_data = jxf_CSRMatrixData(P_diag);
      P_diag_i    = jxf_CSRMatrixI(P_diag);

      for (i = 0; i < n_fine; i++)
      {
         factor = 1.0 / A_diag_data[A_diag_i[i]];
         for (j = P_diag_i[i]; j < P_diag_i[i+1]; j ++)
         {
            P_diag_data[j] *= factor;
         }
      }
   }

   *P_ptr = jxf_CSRMatrixToParCSRMatrix(comm, P_diag, NULL, NULL);

   jxf_ParCSRMatrixDestroy(AR);

   if (my_id == 0)
   {
      jxf_CSRMatrixDestroy(P_diag);
   }

   return 0;
}
