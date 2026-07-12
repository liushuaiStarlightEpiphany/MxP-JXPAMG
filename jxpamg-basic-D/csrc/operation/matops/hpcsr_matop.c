#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

JX_Int 
jx_hpCSRMatrixCopy( jx_hpCSRMatrix *A, 
                     jx_hpCSRMatrix *B, 
                     JX_Int              copy_data )
{
   jx_CSRMatrix *A_diag;
   jx_CSRMatrix *A_offd;
   JX_Int *col_map_offd_A;
   jx_CSRMatrix *B_diag;
   jx_CSRMatrix *B_offd;
   JX_Int *col_map_offd_B;
   JX_Int num_cols_offd;
   JX_Int i;

   if (!A)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   if (!B)
   {
      jx_error_in_arg(1);
      return jx_error_flag;
   }
   A_diag = jx_hpCSRMatrixDiag(A);
   A_offd = jx_hpCSRMatrixOffd(A);
   col_map_offd_A = jx_hpCSRMatrixColMapOffd(A);
   B_diag = jx_hpCSRMatrixDiag(B);
   B_offd = jx_hpCSRMatrixOffd(B);
   col_map_offd_B = jx_hpCSRMatrixColMapOffd(B);
   num_cols_offd = jx_CSRMatrixNumCols(A_offd);

   jx_CSRMatrixCopy(A_diag, B_diag, copy_data);
   jx_CSRMatrixCopy(A_offd, B_offd, copy_data);
   if (num_cols_offd && col_map_offd_B == NULL)
   {
      col_map_offd_B = jx_CTAlloc(JX_Int,num_cols_offd);
      jx_hpCSRMatrixColMapOffd(B) = col_map_offd_B;
   }
   for (i = 0; i < num_cols_offd; i ++)
   {
      col_map_offd_B[i] = col_map_offd_A[i];
   }
        
   return jx_error_flag;
}

// jx_hpCSRMatrix *
// jx_hpMatmulL( jx_hpCSRMatrix *A, jx_ParCSRMatrix *B )
// {
//    jx_hpCSRMatrix *C = jx_CTAlloc(jx_hpCSRMatrix, 1);
//    jx_hpCSRMatrixPar(C) = jx_ParMatmul(jx_hpCSRMatrixPar(A), B);
//    jx_hpCSRMatrixOffdCpu(C) = NULL;
//    jx_hpCSRMatrixOffdNode(C) = NULL;
//    jx_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jx_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }

// jx_hpCSRMatrix *
// jx_hpMatmulR( jx_ParCSRMatrix *A, jx_hpCSRMatrix *B )
// {
//    jx_hpCSRMatrix *C = jx_CTAlloc(jx_hpCSRMatrix, 1);
//    jx_hpCSRMatrixPar(C) = jx_ParMatmul(A, jx_hpCSRMatrixPar(B));
//    jx_hpCSRMatrixOffdCpu(C) = NULL;
//    jx_hpCSRMatrixOffdNode(C) = NULL;
//    jx_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jx_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }


// jx_hpCSRMatrix *
// jx_hpTMatmulL( jx_hpCSRMatrix *A, jx_ParCSRMatrix *B )
// {
//    jx_hpCSRMatrix *C = jx_CTAlloc(jx_hpCSRMatrix, 1);
//    jx_hpCSRMatrixPar(C) = jx_ParTMatmul(jx_hpCSRMatrixPar(A), B);
//    jx_hpCSRMatrixOffdCpu(C) = NULL;
//    jx_hpCSRMatrixOffdNode(C) = NULL;
//    jx_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jx_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }

// jx_hpCSRMatrix *
// jx_hpTMatmulR( jx_ParCSRMatrix *A, jx_hpCSRMatrix *B  )
// {
//    jx_hpCSRMatrix *C = jx_CTAlloc(jx_hpCSRMatrix, 1);
//    jx_hpCSRMatrixPar(C) = jx_ParTMatmul(A, jx_hpCSRMatrixPar(B));
//    jx_hpCSRMatrixOffdCpu(C) = NULL;
//    jx_hpCSRMatrixOffdNode(C) = NULL;
//    jx_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jx_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }