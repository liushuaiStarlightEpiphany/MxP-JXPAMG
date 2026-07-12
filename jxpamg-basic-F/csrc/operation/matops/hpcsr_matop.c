#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

JXF_Int 
jxf_hpCSRMatrixCopy( jxf_hpCSRMatrix *A, 
                     jxf_hpCSRMatrix *B, 
                     JXF_Int              copy_data )
{
   jxf_CSRMatrix *A_diag;
   jxf_CSRMatrix *A_offd;
   JXF_Int *col_map_offd_A;
   jxf_CSRMatrix *B_diag;
   jxf_CSRMatrix *B_offd;
   JXF_Int *col_map_offd_B;
   JXF_Int num_cols_offd;
   JXF_Int i;

   if (!A)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (!B)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   A_diag = jxf_hpCSRMatrixDiag(A);
   A_offd = jxf_hpCSRMatrixOffd(A);
   col_map_offd_A = jxf_hpCSRMatrixColMapOffd(A);
   B_diag = jxf_hpCSRMatrixDiag(B);
   B_offd = jxf_hpCSRMatrixOffd(B);
   col_map_offd_B = jxf_hpCSRMatrixColMapOffd(B);
   num_cols_offd = jxf_CSRMatrixNumCols(A_offd);

   jxf_CSRMatrixCopy(A_diag, B_diag, copy_data);
   jxf_CSRMatrixCopy(A_offd, B_offd, copy_data);
   if (num_cols_offd && col_map_offd_B == NULL)
   {
      col_map_offd_B = jxf_CTAlloc(JXF_Int,num_cols_offd);
      jxf_hpCSRMatrixColMapOffd(B) = col_map_offd_B;
   }
   for (i = 0; i < num_cols_offd; i ++)
   {
      col_map_offd_B[i] = col_map_offd_A[i];
   }
        
   return jxf_error_flag;
}

// jxf_hpCSRMatrix *
// jxf_hpMatmulL( jxf_hpCSRMatrix *A, jxf_ParCSRMatrix *B )
// {
//    jxf_hpCSRMatrix *C = jxf_CTAlloc(jxf_hpCSRMatrix, 1);
//    jxf_hpCSRMatrixPar(C) = jxf_ParMatmul(jxf_hpCSRMatrixPar(A), B);
//    jxf_hpCSRMatrixOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixOffdNode(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }

// jxf_hpCSRMatrix *
// jxf_hpMatmulR( jxf_ParCSRMatrix *A, jxf_hpCSRMatrix *B )
// {
//    jxf_hpCSRMatrix *C = jxf_CTAlloc(jxf_hpCSRMatrix, 1);
//    jxf_hpCSRMatrixPar(C) = jxf_ParMatmul(A, jxf_hpCSRMatrixPar(B));
//    jxf_hpCSRMatrixOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixOffdNode(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }


// jxf_hpCSRMatrix *
// jxf_hpTMatmulL( jxf_hpCSRMatrix *A, jxf_ParCSRMatrix *B )
// {
//    jxf_hpCSRMatrix *C = jxf_CTAlloc(jxf_hpCSRMatrix, 1);
//    jxf_hpCSRMatrixPar(C) = jxf_ParTMatmul(jxf_hpCSRMatrixPar(A), B);
//    jxf_hpCSRMatrixOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixOffdNode(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }

// jxf_hpCSRMatrix *
// jxf_hpTMatmulR( jxf_ParCSRMatrix *A, jxf_hpCSRMatrix *B  )
// {
//    jxf_hpCSRMatrix *C = jxf_CTAlloc(jxf_hpCSRMatrix, 1);
//    jxf_hpCSRMatrixPar(C) = jxf_ParTMatmul(A, jxf_hpCSRMatrixPar(B));
//    jxf_hpCSRMatrixOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixOffdNode(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdCpu(C) = NULL;
//    jxf_hpCSRMatrixColMapOffdNode(C) = NULL;
//    return C;
// }