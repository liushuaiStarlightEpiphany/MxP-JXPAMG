//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_matrix.c -- operations for parallel matrices.
 *  Date: 2014/07/05
 */

#include "jxf_mv.h"

/*!
 * \fn jxf_ParCSRMatrix *jxf_BuildMatParFromOneFile4
 * \brief Build a symmetric parallel matrix by reading data from a given file.
 * \param filename pointer to the matrix file.
 * \param row_part pointer to the row-start partitioning.
 * \param row_part pointer to the col-start partitioning.
 * \note The storage format of the matrix is expected to be CSR, i.e.
 *         the 1st row:                           n                          // JXF_Int
 *         the 2 to (n+2)-th row:                ia[k]     k=0,1,...,n       // JXF_Int
 *         the (n+3) to (n+nz+2)-th row:         ja[k]     k=0,1,...,nz-1    // JXF_Int
 *         the (n+nz+3) to ((n+2nz+2))-th row:    a[k]     k=0,1,...,nz-1    // JXF_Real 
 *       where n is the order of the matrix, nz is the number of non-zeros of the matrix.
 * \author Yue Xiaoqiang
 * \date 2014/04/14
 */
jxf_ParCSRMatrix *
jxf_BuildMatParFromOneFile4( char *filename, JXF_Int file_base, JXF_Int  *row_part, JXF_Int  *col_part )
{
   jxf_ParCSRMatrix  *A     = NULL;
   jxf_CSRMatrix     *A_CSR = NULL;
   jxf_CSRMatrix     *B_CSR = NULL;

   JXF_Int my_id;
   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0)
   {
      B_CSR = jxf_CSRMatrixRead(filename, file_base);
      jxf_CSRMatrixReorder(B_CSR);
      A_CSR = jxf_CSRMatrixSymmetrization(B_CSR);
      jxf_CSRMatrixDestroy(B_CSR);
      jxf_CSRMatrixReorder(A_CSR);
   }
   
   /* ser -> par transfering */
   A = jxf_CSRMatrixToParCSRMatrix(MPI_COMM_WORLD, A_CSR, row_part, col_part);

   if (my_id == 0) 
   {
      jxf_CSRMatrixDestroy(A_CSR);
   }

   return (A);
}
