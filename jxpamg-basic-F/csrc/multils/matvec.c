/*========================================================================*
 *  FSLS (Fast Solvers for Linear System)  (c) 2010-2012                  *
 *  School of Mathematics and Computational Science                       *
 *  Xiangtan University                                                   *
 *  Email: peghoty@163.com                                                *
 *========================================================================*/
 
/*!
 *  matvec.c
 *
 *  Created by peghoty 2010/08/27
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 */ 
 
#include "jxf_multils.h" 
 
/*!
 * \fn JXF_Int fsls_BuildCSRMatFromFile
 * \brief Build a fsls_CSRMatrix matrix from a given file.
 * \param *filename pointer to the matrix file
 * \param **A_ptr pointer to the resulting matrix
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_BuildCSRMatFromFile( char *filename, fsls_CSRMatrix **A_ptr )
{
   fsls_CSRMatrix *A = NULL;

   A = fsls_CSRMatrixRead(filename);
   fsls_CSRMatrixReorder(A); // diagonal elements first
   
   *A_ptr = A;

   return (0);
}

/*!
 * \fn fsls_CSRMatrix * fsls_CSRMatrixRead
 * \brief Read a CSR matrix from a file and Create a 'fsls_CSRMatrix' type matrix.
 * \param *file_name the source matrix name
 * \note the matrix file should be written in the following format:
 *     n
 *     ia[i],i=0(1)n
 *     ja[i],i=0(1)nz-1 (nz=ia[n])
 *      a[i],i=0(1)nz-1 (nz=ia[n])
 *
 * \note my question: according to the code, we know it can deal with either  
 *       square matrix or rectangular matrix as long as the file is written 
 *       in the format refered above, but as for the latter, does 'max_col' 
 *       get the true number of colum if the matrix has a colum with all zero 
 *       entries? In fact, the same problem also exists for a square matrix.
 * \author peghoty
 * \date 2009/12/05
 */
fsls_CSRMatrix *
fsls_CSRMatrixRead( char *file_name )
{
   fsls_CSRMatrix  *matrix      = NULL;
   FILE            *fp          = NULL;
   JXF_Real          *matrix_data = NULL;
   JXF_Int             *matrix_i    = NULL;
   JXF_Int             *matrix_j    = NULL;
   JXF_Int              num_rows;
   JXF_Int              num_nonzeros;
   JXF_Int              max_col     = 0;
   JXF_Int              file_base   = 1;
   JXF_Int              j;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &num_rows);

   matrix_i = fsls_CTAlloc(JXF_Int, num_rows + 1);
   for (j = 0; j < num_rows+1; j ++)
   {
      jxf_fscanf(fp, "%d", &matrix_i[j]);
      matrix_i[j] -= file_base;
   }

   num_nonzeros = matrix_i[num_rows];

   matrix = fsls_CSRMatrixCreate(num_rows, num_rows, matrix_i[num_rows]);
   fsls_CSRMatrixI(matrix) = matrix_i;
   fsls_CSRMatrixInitialize(matrix);

   matrix_j = fsls_CSRMatrixJ(matrix);
   for (j = 0; j < num_nonzeros; j ++)
   {
      jxf_fscanf(fp, "%d", &matrix_j[j]);
      matrix_j[j] -= file_base;

      if (matrix_j[j] > max_col)
      {
         max_col = matrix_j[j];
      }
   }

   matrix_data = fsls_CSRMatrixData(matrix);
   for (j = 0; j < matrix_i[num_rows]; j ++)
   {
      jxf_fscanf(fp, "%le", &matrix_data[j]);
   }

   fclose(fp);

   fsls_CSRMatrixNumNonzeros(matrix) = num_nonzeros;
   fsls_CSRMatrixNumCols(matrix) = ++ max_col;

   return matrix;
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixCreate
 * \brief Create a fsls_CSRMatrix matrix.
 * \param num_rows number of rows
 * \param num_cols number of cols
 * \param num_nonzeros number of nonzeros
 * \return pointer to a fsls_CSRMatrix matrix
 * \author peghoty
 * \date 2010/01/07
 */
fsls_CSRMatrix *
fsls_CSRMatrixCreate( JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros )
{
   fsls_CSRMatrix *matrix = NULL;

   matrix = fsls_CTAlloc(fsls_CSRMatrix, 1);

   fsls_CSRMatrixData(matrix)        = NULL;
   fsls_CSRMatrixI(matrix)           = NULL;
   fsls_CSRMatrixJ(matrix)           = NULL;
   fsls_CSRMatrixNumRows(matrix)     = num_rows;
   fsls_CSRMatrixNumCols(matrix)     = num_cols;
   fsls_CSRMatrixNumNonzeros(matrix) = num_nonzeros;

   return matrix;
}

/*!
 * \fn fsls_ICSRMatrix * fsls_ICSRMatrixCreate
 * \brief Create a fsls_ICSRMatrix matrix.
 * \param num_rows number of rows
 * \param num_cols number of cols
 * \param num_nonzeros number of nonzeros
 * \return pointer to a fsls_ICSRMatrix matrix
 * \author peghoty
 * \date 2010/11/11
 */
fsls_ICSRMatrix *
fsls_ICSRMatrixCreate( JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros )
{
   fsls_ICSRMatrix *matrix = NULL;

   matrix = fsls_CTAlloc(fsls_ICSRMatrix, 1);

   fsls_ICSRMatrixData(matrix)        = NULL;
   fsls_ICSRMatrixI(matrix)           = NULL;
   fsls_ICSRMatrixJ(matrix)           = NULL;
   fsls_ICSRMatrixNumRows(matrix)     = num_rows;
   fsls_ICSRMatrixNumCols(matrix)     = num_cols;
   fsls_ICSRMatrixNumNonzeros(matrix) = num_nonzeros;

   return matrix;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixInitialize
 * \brief Initialize a fsls_CSRMatrix matrix.
 * \param *matrix pointer to the matrix to be initialized
 * \date 2010/01/07
 */
JXF_Int 
fsls_CSRMatrixInitialize( fsls_CSRMatrix *matrix )
{
   JXF_Int  num_rows     = fsls_CSRMatrixNumRows(matrix);
   JXF_Int  num_nonzeros = fsls_CSRMatrixNumNonzeros(matrix);
   JXF_Int  ierr = 0;

   if (!fsls_CSRMatrixData(matrix) && num_nonzeros)
   {
      fsls_CSRMatrixData(matrix) = fsls_CTAlloc(JXF_Real, num_nonzeros);
   }
   if (!fsls_CSRMatrixI(matrix))
   {
      fsls_CSRMatrixI(matrix) = fsls_CTAlloc(JXF_Int, num_rows + 1);
   }
   if (!fsls_CSRMatrixJ(matrix) && num_nonzeros)
   {
      fsls_CSRMatrixJ(matrix) = fsls_CTAlloc(JXF_Int, num_nonzeros);
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_ICSRMatrixInitialize
 * \brief Initialize a fsls_ICSRMatrix matrix.
 * \param *matrix pointer to the matrix to be initialized
 * \date 2010/11/11
 */
JXF_Int 
fsls_ICSRMatrixInitialize( fsls_ICSRMatrix *matrix )
{
   JXF_Int  num_rows     = fsls_ICSRMatrixNumRows(matrix);
   JXF_Int  num_nonzeros = fsls_ICSRMatrixNumNonzeros(matrix);
   JXF_Int  ierr = 0;

   if (!fsls_ICSRMatrixData(matrix) && num_nonzeros)
   {
      fsls_ICSRMatrixData(matrix) = fsls_CTAlloc(JXF_Int, num_nonzeros);
   }
   if (!fsls_ICSRMatrixI(matrix))
   {
      fsls_ICSRMatrixI(matrix) = fsls_CTAlloc(JXF_Int, num_rows + 1);
   }
   if (!fsls_ICSRMatrixJ(matrix) && num_nonzeros)
   {
      fsls_ICSRMatrixJ(matrix) = fsls_CTAlloc(JXF_Int, num_nonzeros);
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixReorder
 * \brief Reorder the column and data arrays of a square CSR matrix, 
 *        so that the first entry in each row is the diagonal one.
 * \param *A pointer to the matrix to be reordered
 * \date 2010/01/07
 */
JXF_Int 
fsls_CSRMatrixReorder( fsls_CSRMatrix *A )
{
   JXF_Int    i, j, tempi, row_size;
   JXF_Real tempd;

   JXF_Real *A_data = fsls_CSRMatrixData(A);
   JXF_Int    *A_i    = fsls_CSRMatrixI(A);
   JXF_Int    *A_j    = fsls_CSRMatrixJ(A);
   JXF_Int     num_rowsA = fsls_CSRMatrixNumRows(A);
   JXF_Int     num_colsA = fsls_CSRMatrixNumCols(A);

   /* the matrix should be square */
   if (num_rowsA != num_colsA) return -1;

   for (i = 0; i < num_rowsA; i ++)
   {
      row_size = A_i[i+1] - A_i[i];

      for (j = 0; j < row_size; j ++)
      {
         if (A_j[j] == i)
         {
            if (j != 0)
            {
               tempi  = A_j[0];
               A_j[0] = A_j[j];
               A_j[j] = tempi;

               tempd     = A_data[0];
               A_data[0] = A_data[j];
               A_data[j] = tempd;
            }
            break;
         }

         /* diagonal element is missing */
         if (j == row_size-1) return -2;
      }

      A_j    += row_size;
      A_data += row_size;
   }

   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixPrint
 * \brief Write a fsls_CSRMatrix into a file.
 * \param *A pointer to the matrix to be written
 * \param *file_name pointer to the file where the matrix to be kept
 * \date 2010/01/07
 */
JXF_Int
fsls_CSRMatrixPrint( fsls_CSRMatrix *matrix, char *file_name )
{
   FILE    *fp = NULL;

   JXF_Real  *matrix_data = fsls_CSRMatrixData(matrix);
   JXF_Int     *matrix_i    = fsls_CSRMatrixI(matrix);
   JXF_Int     *matrix_j    = fsls_CSRMatrixJ(matrix);
   JXF_Int      num_rows    = fsls_CSRMatrixNumRows(matrix);
   
   JXF_Int      file_base = 1; 
   JXF_Int      j, ierr = 0;

   fp = fopen(file_name, "w");

   jxf_fprintf(fp, "%d\n", num_rows);

   for (j = 0; j <= num_rows; j ++)
   {
      jxf_fprintf(fp, "%d\n", matrix_i[j] + file_base);
   }

   for (j = 0; j < matrix_i[num_rows]; j ++)
   {
      jxf_fprintf(fp, "%d\n", matrix_j[j] + file_base);
   }

   if (matrix_data)
   {
      for (j = 0; j < matrix_i[num_rows]; j ++)
      {
         jxf_fprintf(fp, "%.15le\n", matrix_data[j]); // we always use "%.15le\n" 
      }
   }
   else
   {
      jxf_fprintf(fp, "Warning: No matrix data!\n");
   }

   fclose(fp);

   return ierr;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixDestroy
 * \brief Destroy a fsls_CSRMatrix.
 * \param *matrix pointer to the matrix to be destroyed
 * \date 2010/01/07
 */
JXF_Int 
fsls_CSRMatrixDestroy( fsls_CSRMatrix *matrix )
{
   JXF_Int  ierr = 0;
   if (matrix)
   {  
      if (fsls_CSRMatrixI(matrix))    fsls_TFree(fsls_CSRMatrixI(matrix));
      if (fsls_CSRMatrixData(matrix)) fsls_TFree(fsls_CSRMatrixData(matrix));
      if (fsls_CSRMatrixJ(matrix))    fsls_TFree(fsls_CSRMatrixJ(matrix));
      fsls_TFree(matrix);
   }
   return ierr;
}

/*!
 * \fn JXF_Int fsls_ICSRMatrixDestroy
 * \brief Destroy a fsls_ICSRMatrix.
 * \param *matrix pointer to the matrix to be destroyed
 * \date 2010/11/11
 */
JXF_Int 
fsls_ICSRMatrixDestroy( fsls_ICSRMatrix *matrix )
{
   JXF_Int  ierr = 0;
   if (matrix)
   {  
      if (fsls_ICSRMatrixI(matrix))    fsls_TFree(fsls_ICSRMatrixI(matrix));
      if (fsls_ICSRMatrixData(matrix)) fsls_TFree(fsls_ICSRMatrixData(matrix));
      if (fsls_ICSRMatrixJ(matrix))    fsls_TFree(fsls_ICSRMatrixJ(matrix));
      fsls_TFree(matrix);
   }
   return ierr;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixTranspose
 * \brief get the transpose of a given CSR matrix
 * \param *A  the pointer to the given CSR matrix
 * \param *AT the pointer to the result CSR matrix after transposing
 * \param data integer, data = 0 means just transpoing the 
 *        nonzeros-structure while data > 0 means transpoing all 
 *        including the nonzeros entries.
 * \note After transposing, diagonal entries are not firsly-stored 
 *       in each row any more.  2011/03/29
 * \date 2009/12/05
 */
JXF_Int 
fsls_CSRMatrixTranspose( fsls_CSRMatrix *A, fsls_CSRMatrix **AT, JXF_Int data )
{
   JXF_Real *A_data        = fsls_CSRMatrixData(A);
   JXF_Int    *A_i           = fsls_CSRMatrixI(A);
   JXF_Int    *A_j           = fsls_CSRMatrixJ(A);
   JXF_Int     num_rowsA     = fsls_CSRMatrixNumRows(A);
   JXF_Int     num_colsA     = fsls_CSRMatrixNumCols(A);
   JXF_Int     num_nonzerosA = fsls_CSRMatrixNumNonzeros(A);

   JXF_Real *AT_data = NULL;
   JXF_Int    *AT_i    = NULL;
   JXF_Int    *AT_j    = NULL;
   JXF_Int     num_rowsAT;
   JXF_Int     num_colsAT;
   JXF_Int     num_nonzerosAT;

   JXF_Int     max_col;
   JXF_Int     index;
   JXF_Int     i,j;

  /*------------------------------------------------------- 
   * First, ascertain that num_cols and num_nonzeros 
   * has been set. If not, set them.
   *------------------------------------------------------*/
   
   if (!num_nonzerosA)
   {
      num_nonzerosA = A_i[num_rowsA];
   }

   if (num_rowsA && !num_colsA)
   {
      max_col = -1;
      for (i = 0; i < num_rowsA; ++ i)
      {
          for (j = A_i[i]; j < A_i[i+1]; j++)
          {
             if (A_j[j] > max_col) max_col = A_j[j];
          }
      }
      num_colsA = max_col + 1;
   }

   num_rowsAT = num_colsA;
   num_colsAT = num_rowsA;
   num_nonzerosAT = num_nonzerosA;

   *AT = fsls_CSRMatrixCreate(num_rowsAT, num_colsAT, num_nonzerosAT);

   AT_i = fsls_CTAlloc(JXF_Int, num_rowsAT+1);
   AT_j = fsls_CTAlloc(JXF_Int, num_nonzerosAT);
   fsls_CSRMatrixI(*AT) = AT_i;
   fsls_CSRMatrixJ(*AT) = AT_j;
   if (data) 
   {
      AT_data = fsls_CTAlloc(JXF_Real, num_nonzerosAT);
      fsls_CSRMatrixData(*AT) = AT_data;
   }

  /*-------------------------------------------------------
   * Count the number of entries in each column of A 
   * (row of AT) and fill the AT_i array.
   *------------------------------------------------------*/

   for (i = 0; i < num_nonzerosA; i ++)
   {
       ++ AT_i[A_j[i]+1];
   }

   for (i = 2; i <= num_rowsAT; i ++)
   {
       AT_i[i] += AT_i[i-1];
   }

  /*------------------------------------------------------
   * Load the data and column numbers of AT
   *----------------------------------------------------*/

   for (i = 0; i < num_rowsA; i ++)
   {
      for (j = A_i[i]; j < A_i[i+1]; j ++)
      {
         index = A_j[j];
         AT_j[AT_i[index]] = i;
         if (data) AT_data[AT_i[index]] = A_data[j];
         AT_i[index] ++;
      }
   }

  /*------------------------------------------------------------
   * AT_i[j] now points to the *end* of the jth row of entries
   * instead of the beginning.  Restore AT_i to front of row.
   *------------------------------------------------------------*/

   for (i = num_rowsAT; i > 0; i --)
   {
       AT_i[i] = AT_i[i-1];
   }

   AT_i[0] = 0;

   return(0);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixScaledNorm
 * \brief Compute the 'inf-norm' of the matrix D^{-1/2}*A*D^{-1/2}.
 * \param *A the pointer to the matrix(fsls_CSRMatrix)
 * \param *scnorm pointer to the 'inf-norm'
 * \date 2009/12/29
 * \author peghoty
 */
JXF_Int
fsls_CSRMatrixScaledNorm( fsls_CSRMatrix *A, JXF_Real *scnorm )
{
   JXF_Int      *diag_i    = fsls_CSRMatrixI(A);
   JXF_Int      *diag_j    = fsls_CSRMatrixJ(A);
   JXF_Real   *diag_data = fsls_CSRMatrixData(A);
   JXF_Int       num_rows  = fsls_CSRMatrixNumRows(A);

   fsls_Vector   *dinvsqrt = NULL;
   fsls_Vector   *sum      = NULL;
   JXF_Real        *dis_data = NULL;   
   JXF_Real        *sum_data = NULL;
  
   JXF_Int	      i, j;
   JXF_Real     mat_norm;

   dinvsqrt = fsls_SeqVectorCreate(num_rows);
   fsls_SeqVectorInitialize(dinvsqrt);
   dis_data = fsls_VectorData(dinvsqrt);
   sum = fsls_SeqVectorCreate(num_rows);
   fsls_SeqVectorInitialize(sum);
   sum_data = fsls_VectorData(sum);

   /* generate dinvsqrt */
   for (i = 0; i < num_rows; i ++)
   {
      dis_data[i] = 1.0 / sqrt(fabs(diag_data[diag_i[i]]));
   }
   
   for (i = 0; i < num_rows; i ++)
   {
      for (j = diag_i[i]; j < diag_i[i+1]; j ++)
      {
	 sum_data[i] += fabs(diag_data[j])*dis_data[i]*dis_data[diag_j[j]];
      }
   }   

   mat_norm = 0.0;
   for (i = 0; i < num_rows; i ++)
   {
      if (mat_norm < sum_data[i]) mat_norm = sum_data[i];
   }	

   fsls_SeqVectorDestroy(dinvsqrt);
   fsls_SeqVectorDestroy(sum);

   *scnorm = mat_norm;  

   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixCopy
 * \brief Copy A to B (A -> B).
 * \param *A pointer to the matrix being copied.
 * \param *B pointer to the resulting matrix.
 * \note if copy_data = 0 only the structure of A is copied to B.
 *       the routine does not check if the dimensions of A and B match !!!
 * \date 2010/01/07
 */
JXF_Int 
fsls_CSRMatrixCopy( fsls_CSRMatrix *A, fsls_CSRMatrix *B, JXF_Int copy_data )
{
   JXF_Int     ierr     = 0;
   JXF_Int     num_rows = fsls_CSRMatrixNumRows(A);
   JXF_Int    *A_i      = fsls_CSRMatrixI(A);
   JXF_Int    *A_j      = fsls_CSRMatrixJ(A);
   JXF_Real *A_data;
   JXF_Int    *B_i      = fsls_CSRMatrixI(B);
   JXF_Int    *B_j      = fsls_CSRMatrixJ(B);
   JXF_Real *B_data;

   JXF_Int i,j;

   for (i = 0; i < num_rows; i ++)
   {
	B_i[i] = A_i[i];
	for (j = A_i[i]; j < A_i[i+1]; j ++)
	{
           B_j[j] = A_j[j];
	}
   }
   B_i[num_rows] = A_i[num_rows];
   if (copy_data)
   {
	A_data = fsls_CSRMatrixData(A);
	B_data = fsls_CSRMatrixData(B);
   	for (i = 0; i < num_rows; i ++)
   	{
	   for (j = A_i[i]; j < A_i[i+1]; j ++)
	   {
               B_data[j] = A_data[j];
	   }
	}
   }
   return ierr;
}

/*!
 * \fn fsls_CSRMatrix * fsls_CSRMatrixAdd
 * \brief adds two CSR Matrices A and B and returns a CSR Matrix C
 * \param *A the pointer to the first matrix
 * \param *B the pointer to the second matrix
 * \note The routine does not check for 0-elements which might be generated
 *       through cancellation of elements in A and B or already contained
 *	 in A and B. To remove those, use fsls_CSRMatrixDeleteZeros
 * \date 2009/12/05
 */
fsls_CSRMatrix *
fsls_CSRMatrixAdd( fsls_CSRMatrix *A, fsls_CSRMatrix *B )
{
   JXF_Real     *A_data   = fsls_CSRMatrixData(A);
   JXF_Int        *A_i      = fsls_CSRMatrixI(A);
   JXF_Int        *A_j      = fsls_CSRMatrixJ(A);
   JXF_Int         nrows_A  = fsls_CSRMatrixNumRows(A);
   JXF_Int         ncols_A  = fsls_CSRMatrixNumCols(A);
   
   JXF_Real     *B_data   = fsls_CSRMatrixData(B);
   JXF_Int        *B_i      = fsls_CSRMatrixI(B);
   JXF_Int        *B_j      = fsls_CSRMatrixJ(B);
   JXF_Int         nrows_B  = fsls_CSRMatrixNumRows(B);
   JXF_Int         ncols_B  = fsls_CSRMatrixNumCols(B);
   
   fsls_CSRMatrix  *C      = NULL;
   JXF_Real          *C_data = NULL;
   JXF_Int	           *C_i    = NULL;
   JXF_Int             *C_j    = NULL;

   JXF_Int         ia, ib, ic, jcol, num_nonzeros;
   JXF_Int	       pos;
   JXF_Int        *marker;

   /* check the compatibility of A and B */
   if (nrows_A != nrows_B || ncols_A != ncols_B)
   {
      jxf_printf(" >>> Warning: Incompatible matrix dimensions!\n");
      return NULL;
   }

   marker = fsls_CTAlloc(JXF_Int, ncols_A);
   C_i = fsls_CTAlloc(JXF_Int, nrows_A+1);

   /* initialize the marker array */
   for (ia = 0; ia < ncols_A; ia ++)
   {
      marker[ia] = -1;
   }

   num_nonzeros = 0;
   C_i[0] = 0;
   for (ic = 0; ic < nrows_A; ic ++)
   {
	for (ia = A_i[ic]; ia < A_i[ic+1]; ia ++)
	{
            jcol = A_j[ia];
            marker[jcol] = ic;
            num_nonzeros ++;
	}
	for (ib = B_i[ic]; ib < B_i[ic+1]; ib ++)
	{
            jcol = B_j[ib];
            if (marker[jcol] != ic)
            {
               marker[jcol] = ic;
               num_nonzeros ++;
            }
   	}
	C_i[ic+1] = num_nonzeros;
   }

   C = fsls_CSRMatrixCreate(nrows_A, ncols_A, num_nonzeros);
   fsls_CSRMatrixI(C) = C_i;
   fsls_CSRMatrixInitialize(C);
   C_j = fsls_CSRMatrixJ(C);
   C_data = fsls_CSRMatrixData(C);

   /* initialize the marker array again */
   for (ia = 0; ia < ncols_A; ia ++)
   {
      marker[ia] = -1;
   }

   pos = 0;
   for (ic = 0; ic < nrows_A; ic ++)
   {
	for (ia = A_i[ic]; ia < A_i[ic+1]; ia ++)
	{
            jcol = A_j[ia];
            C_j[pos] = jcol;
            C_data[pos] = A_data[ia];
            marker[jcol] = pos;
            pos ++;
	}
	for (ib = B_i[ic]; ib < B_i[ic+1]; ib ++)
	{
            jcol = B_j[ib];
            if (marker[jcol] < C_i[ic])  /* very skillful */
            {  /* the col number 'jcol' doesn't appear in the ic-th row of A */
               C_j[pos] = jcol;
               C_data[pos] = B_data[ib];
               marker[jcol] = pos;
               pos ++;
            }
            else
            {  /* the col number 'jcol' has appeared in the ic-th row of A */
               C_data[marker[jcol]] += B_data[ib];
            }
   	}
   }

   fsls_TFree(marker);
   return C;
}	

/*!
 * \fn fsls_CSRMatrix * fsls_CSRMatrixMultiply
 * \brief multiplies two CSR Matrices A and B and returns a CSR Matrix C
 * \param *A the pointer to the first matrix
 * \param *B the pointer to the second matrix
 * \note The routine does not check for 0-elements which might be generated
 *       through cancellation of elements in A and B or already contained
 *	 in A and B. To remove those, use fsls_CSRMatrixDeleteZeros
 * \note The similar strategy is applied in both 'fsls_CSRMatrixAdd'
 *       and 'fsls_CSRMatrixMultiply'
 * \date 2009/12/05
 */
fsls_CSRMatrix *
fsls_CSRMatrixMultiply( fsls_CSRMatrix *A, fsls_CSRMatrix *B )
{
   JXF_Real     *A_data   = fsls_CSRMatrixData(A);
   JXF_Int        *A_i      = fsls_CSRMatrixI(A);
   JXF_Int        *A_j      = fsls_CSRMatrixJ(A);
   JXF_Int         nrows_A  = fsls_CSRMatrixNumRows(A);
   JXF_Int         ncols_A  = fsls_CSRMatrixNumCols(A);
   
   JXF_Real     *B_data   = fsls_CSRMatrixData(B);
   JXF_Int        *B_i      = fsls_CSRMatrixI(B);
   JXF_Int        *B_j      = fsls_CSRMatrixJ(B);
   JXF_Int         nrows_B  = fsls_CSRMatrixNumRows(B);
   JXF_Int         ncols_B  = fsls_CSRMatrixNumCols(B);
   
   fsls_CSRMatrix *C      = NULL;
   JXF_Real         *C_data = NULL;
   JXF_Int	          *C_i    = NULL;
   JXF_Int            *C_j    = NULL;

   JXF_Int         ia, ib, ic, ja, jb;
   JXF_Int         num_nonzeros = 0;
   JXF_Int	       row_start, counter;
   JXF_Real      a_entry, b_entry;
   JXF_Int        *B_marker = NULL;

   /* check the compatibility of A and B */
   if (ncols_A != nrows_B)
   {
      jxf_printf(" >>> Warning: Incompatible matrix dimensions!\n");
      return NULL;
   }

   B_marker = fsls_CTAlloc(JXF_Int, ncols_B);
   C_i = fsls_CTAlloc(JXF_Int, nrows_A+1);

   /* initialize the marker array */
   for (ib = 0; ib < ncols_B; ib ++)
   {
      B_marker[ib] = -1;
   }

   /* obtain the nonzero-structure of C */
   for (ic = 0; ic < nrows_A; ic ++)
   {
	for (ia = A_i[ic]; ia < A_i[ic+1]; ia ++)
	{
            ja = A_j[ia];
            for (ib = B_i[ja]; ib < B_i[ja+1]; ib ++)
            {
               jb = B_j[ib];
               if (B_marker[jb] != ic)
               {
                  B_marker[jb] = ic;
                  num_nonzeros ++;
               }
            }
   	}
	C_i[ic+1] = num_nonzeros;
   }

   C = fsls_CSRMatrixCreate(nrows_A, ncols_B, num_nonzeros);
   fsls_CSRMatrixI(C) = C_i;
   fsls_CSRMatrixInitialize(C);
   C_j = fsls_CSRMatrixJ(C);
   C_data = fsls_CSRMatrixData(C);

   /* initialize the marker array again */
   for (ib = 0; ib < ncols_B; ib ++)
   {
      B_marker[ib] = -1;
   }

   /* fill in the nonzero entries of C */
   counter = 0;
   for (ic = 0; ic < nrows_A; ic ++)
   {
	row_start = C_i[ic];
	for (ia = A_i[ic]; ia < A_i[ic+1]; ia ++)
	{
             ja = A_j[ia];
             a_entry = A_data[ia];
             for (ib = B_i[ja]; ib < B_i[ja+1]; ib ++)
             {
                 jb = B_j[ib];
                 b_entry = B_data[ib];
                 if (B_marker[jb] < row_start)
                 {
                    B_marker[jb] = counter;
                    C_j[B_marker[jb]] = jb;
                    C_data[B_marker[jb]] = a_entry*b_entry;
                    counter ++;
                 }
                 else
                 {
                    C_data[B_marker[jb]] += a_entry*b_entry;
                 }	 
             }
	}
   }
   fsls_TFree(B_marker);
   return C;
}	

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixTriMultiply
 * \brief Compute RAP = R*A*P.
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix 
 * \note RAP satisfies the "diagonal-entries-firstly-stored-in-each-row" requirement. 2011/03/29 
 * \author peghoty
 * \date 2010/12/09
 */
fsls_CSRMatrix *
fsls_CSRMatrixTriMultiply( fsls_CSRMatrix   *R,
                           fsls_CSRMatrix   *A,
                           fsls_CSRMatrix   *P  )
{
   fsls_CSRMatrix  *RAP;
   
   JXF_Real          *A_data;
   JXF_Int             *A_i;
   JXF_Int             *A_j;

   JXF_Real          *P_data;
   JXF_Int             *P_i;
   JXF_Int             *P_j;

   JXF_Real          *RAP_data;
   JXF_Int             *RAP_i;
   JXF_Int             *RAP_j;

   JXF_Int              RAP_size;
   
   JXF_Real          *R_data;
   JXF_Int             *R_i;
   JXF_Int             *R_j;

   JXF_Int             *P_marker;
   JXF_Int             *A_marker;

   JXF_Int              n_coarse;
   JXF_Int              n_fine;
   
   JXF_Int              ic, i;
   JXF_Int              i1, i2, i3;
   JXF_Int              jj1, jj2, jj3;
   
   JXF_Int              jj_counter;
   JXF_Int              jj_row_begining;

   JXF_Real           r_entry;
   JXF_Real           r_a_product;
   JXF_Real           r_a_p_product;
   
   JXF_Real           zero = 0.0;

  /*-------------------------------------------------
   *  Access the CSR vectors for R, A, P. 
   *  Also get sizes of fine and coarse grids.
   *------------------------------------------------*/

   R_data = fsls_CSRMatrixData(R);
   R_i    = fsls_CSRMatrixI(R);
   R_j    = fsls_CSRMatrixJ(R);

   A_data = fsls_CSRMatrixData(A);
   A_i    = fsls_CSRMatrixI(A);
   A_j    = fsls_CSRMatrixJ(A);

   P_data = fsls_CSRMatrixData(P);
   P_i    = fsls_CSRMatrixI(P);
   P_j    = fsls_CSRMatrixJ(P);

   n_fine   = fsls_CSRMatrixNumRows(A);
   n_coarse = fsls_CSRMatrixNumRows(R);

  /*-----------------------------------------------------------------------
   *  Allocate RAP_i and marker arrays.
   *----------------------------------------------------------------------*/

   RAP_i    = fsls_CTAlloc(JXF_Int, n_coarse+1);
   P_marker = fsls_CTAlloc(JXF_Int, n_coarse);
   A_marker = fsls_CTAlloc(JXF_Int, n_fine);


  /*-----------------------------------------------------------------------
   *  First Pass: Determine size of RAP and set up RAP_i
   *----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------
   *  Intialize some stuff.
   *----------------------------------------------------------------------*/

   jj_counter = 0;
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -1;
   }
   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -1;
   }   

  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < n_coarse; ic ++)
   {
      
     /*--------------------------------------------------------------------
      *  Set marker for diagonal entry, RAP_{ic,ic}.
      *-------------------------------------------------------------------*/

      P_marker[ic] = jj_counter;
      jj_row_begining = jj_counter;
      jj_counter ++;

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/
   
      for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
      {
         i1 = R_j[jj1];

        /*-----------------------------------------------------------------
         *  Loop over entries in row i1 of A.
         *----------------------------------------------------------------*/
         
         for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
         {
            i2 = A_j[jj2];

           /*----------------------------------------------------------------
            *  Check A_marker to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/

            if (A_marker[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_marker[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
               {
                  i3 = P_j[jj3];
                  
                 /*------------------------------------------------------------
                  *  Check P_marker to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, mark it and increment
                  *  counter.
                  *-----------------------------------------------------------*/

                  if (P_marker[i3] < jj_row_begining)
                  {
                     P_marker[i3] = jj_counter;
                     jj_counter ++;
                  }
               }
            }
         }
      }
            
     /*--------------------------------------------------------------------
      * Set RAP_i for this row.
      *--------------------------------------------------------------------*/

      RAP_i[ic] = jj_row_begining;
      
   }
  
   RAP_i[n_coarse] = jj_counter;
 
  /*-----------------------------------------------------------------------
   *  Allocate RAP_data and RAP_j arrays.
   *----------------------------------------------------------------------*/

   RAP_size = jj_counter;
   RAP_data = fsls_CTAlloc(JXF_Real, RAP_size);
   RAP_j    = fsls_CTAlloc(JXF_Int, RAP_size);

  /*-----------------------------------------------------------------------
   *  Second Pass: Fill in RAP_data and RAP_j.
   *----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------
   *  Intialize some stuff.
   *----------------------------------------------------------------------*/

   jj_counter = 0;
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -1;
   }
   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -1;
   }   
   
  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < n_coarse; ic ++)
   {
      
     /*--------------------------------------------------------------------
      *  Create diagonal entry, RAP_{ic,ic}.
      *-------------------------------------------------------------------*/

      P_marker[ic] = jj_counter;
      jj_row_begining = jj_counter;
      RAP_data[jj_counter] = zero;
      RAP_j[jj_counter] = ic;
      jj_counter ++;

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/
   
      for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
      {
         i1 = R_j[jj1];
         r_entry = R_data[jj1];

        /*-----------------------------------------------------------------
         *  Loop over entries in row i1 of A.
         *-----------------------------------------------------------------*/
         
         for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
         {
            i2 = A_j[jj2];
            r_a_product = r_entry * A_data[jj2];
            
           /*----------------------------------------------------------------
            *  Check A_marker to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/

            if (A_marker[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_marker[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product * P_data[jj3];
                  
                 /*-----------------------------------------------------------
                  *  Check P_marker to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, create a new entry.
                  *  If it has, add new contribution.
                  *----------------------------------------------------------*/

                  if (P_marker[i3] < jj_row_begining)
                  {
                     P_marker[i3] = jj_counter;
                     RAP_data[jj_counter] = r_a_p_product;
                     RAP_j[jj_counter] = i3;
                     jj_counter ++;
                  }
                  else
                  {
                     RAP_data[P_marker[i3]] += r_a_p_product;
                  }
               }
            }

           /*--------------------------------------------------------------
            *  If i2 is previously visted ( A_marker[12]=ic ) it yields
            *  no new entries in RAP and can just add new contributions.
            *--------------------------------------------------------------*/

            else
            {
               for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product*P_data[jj3];
                  RAP_data[P_marker[i3]] += r_a_p_product;
               }
            }
         }
      }
   }

   RAP = fsls_CSRMatrixCreate(n_coarse, n_coarse, RAP_size);
   fsls_CSRMatrixData(RAP) = RAP_data; 
   fsls_CSRMatrixI(RAP) = RAP_i; 
   fsls_CSRMatrixJ(RAP) = RAP_j;
   
  /*--------------------------------------------
   *  Free R and marker arrays.
   *-------------------------------------------*/

   fsls_TFree(P_marker);   
   fsls_TFree(A_marker);

   return (RAP);
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_opt
 * \brief Compute RAP = R*A*P.
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix 
 * \note optimize the "ia[k+1]" in the for-loops.
 * \note RAP satisfies the "diagonal-entries-firstly-stored-in-each-row" requirement. 2011/03/29
 * \author peghoty
 * \date 2011/03/21
 */
fsls_CSRMatrix *
fsls_CSRMatrixTriMultiply_opt( fsls_CSRMatrix   *R,
                               fsls_CSRMatrix   *A,
                               fsls_CSRMatrix   *P  )
{
   //struct timeval tStart,tEnd;

   fsls_CSRMatrix  *RAP;
   
   JXF_Real          *A_data;
   JXF_Int             *A_i;
   JXF_Int             *A_j;

   JXF_Real          *P_data;
   JXF_Int             *P_i;
   JXF_Int             *P_j;

   JXF_Real          *RAP_data;
   JXF_Int             *RAP_i;
   JXF_Int             *RAP_j;

   JXF_Int              RAP_size;
   
   JXF_Real          *R_data;
   JXF_Int             *R_i;
   JXF_Int             *R_j;

   JXF_Int             *P_marker;
   JXF_Int             *A_marker;

   JXF_Int              n_coarse;
   JXF_Int              n_fine;
   
   JXF_Int              ic, i;
   JXF_Int              i1, i2, i3;
   JXF_Int              jj1, jj2, jj3;
   
   JXF_Int              jj_counter;
   JXF_Int              jj_row_begining;
   
   JXF_Int              beginR, beginA, beginP;
   JXF_Int              endR, endA, endP;
   
   JXF_Real           r_entry;
   JXF_Real           r_a_product;
   JXF_Real           r_a_p_product;
   
   JXF_Real           zero = 0.0;

  /*--------------------------------------------------------
   *  Access the CSR vectors for R, A, P. 
   *  Also get sizes of fine and coarse grids.
   *-------------------------------------------------------*/

   R_data = fsls_CSRMatrixData(R);
   R_i    = fsls_CSRMatrixI(R);
   R_j    = fsls_CSRMatrixJ(R);

   A_data = fsls_CSRMatrixData(A);
   A_i    = fsls_CSRMatrixI(A);
   A_j    = fsls_CSRMatrixJ(A);

   P_data = fsls_CSRMatrixData(P);
   P_i    = fsls_CSRMatrixI(P);
   P_j    = fsls_CSRMatrixJ(P);

   n_fine   = fsls_CSRMatrixNumRows(A);
   n_coarse = fsls_CSRMatrixNumRows(R);

  /*-----------------------------------------------------------------------
   *  Allocate RAP_i and marker arrays.
   *----------------------------------------------------------------------*/

   RAP_i    = fsls_CTAlloc(JXF_Int, n_coarse+1);
   P_marker = fsls_CTAlloc(JXF_Int, n_coarse);
   A_marker = fsls_CTAlloc(JXF_Int, n_fine);
   
   
  /*-----------------------------------------------------------------------
   *  First Pass: Determine size of RAP and set up RAP_i
   *----------------------------------------------------------------------*/

   //GetTime(tStart);
   
  /*-----------------------------------------------------------------------
   *  Intialize some stuff.
   *----------------------------------------------------------------------*/

   jj_counter = 0;
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -1;
   }
   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -1;
   }   

  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < n_coarse; ic ++)
   {
      
     /*--------------------------------------------------------------------
      *  Set marker for diagonal entry, RAP_{ic,ic}.
      *-------------------------------------------------------------------*/

      P_marker[ic]    = jj_counter;
      jj_row_begining = jj_counter;
      jj_counter ++;

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/
      
      beginR = R_i[ic];
      endR = R_i[ic+1];
      for (jj1 = beginR; jj1 < endR; jj1 ++)
      {
         i1 = R_j[jj1];

        /*-----------------------------------------------------------------
         *  Loop over entries in row i1 of A.
         *----------------------------------------------------------------*/

         beginA = A_i[i1];
         endA = A_i[i1+1];
         for (jj2 = beginA; jj2 < endA; jj2 ++)
         {
            i2 = A_j[jj2];

           /*----------------------------------------------------------------
            *  Check A_marker to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/

            if (A_marker[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_marker[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               beginP = P_i[i2];
               endP = P_i[i2+1];
               for (jj3 = beginP; jj3 < endP; jj3 ++)
               {
                  i3 = P_j[jj3];
                  
                 /*------------------------------------------------------------
                  *  Check P_marker to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, mark it and increment
                  *  counter.
                  *-----------------------------------------------------------*/

                  if (P_marker[i3] < jj_row_begining)
                  {
                     P_marker[i3] = jj_counter;
                     jj_counter ++;
                  }
               }
            }
         }
      }
            
     /*--------------------------------------------------------------------
      * Set RAP_i for this row.
      *--------------------------------------------------------------------*/

      RAP_i[ic] = jj_row_begining;
      
   }
  
   RAP_i[n_coarse] = jj_counter;
   
   //GetTime(tEnd);
   //jxf_printf(" >>> \033[31mthe 1st pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd));     
 
  /*-----------------------------------------------------------------------
   *  Allocate RAP_data and RAP_j arrays.
   *----------------------------------------------------------------------*/

   RAP_size = jj_counter;
   RAP_data = fsls_CTAlloc(JXF_Real, RAP_size);
   RAP_j    = fsls_CTAlloc(JXF_Int, RAP_size);


  /*-----------------------------------------------------------------------
   *  Second Pass: Fill in RAP_data and RAP_j.
   *----------------------------------------------------------------------*/
   
   //GetTime(tStart);
 
  /*-----------------------------------------------------------------------
   *  Intialize some stuff.
   *----------------------------------------------------------------------*/

   jj_counter = 0;
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -1;
   }
   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -1;
   }   
   
  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < n_coarse; ic ++)
   {
      
     /*--------------------------------------------------------------------
      *  Create diagonal entry, RAP_{ic,ic}.
      *-------------------------------------------------------------------*/

      P_marker[ic]         = jj_counter;
      jj_row_begining      = jj_counter;
      RAP_data[jj_counter] = zero;
      RAP_j[jj_counter]    = ic;
      jj_counter ++;

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/
 
      beginR = R_i[ic];
      endR = R_i[ic+1];
      for (jj1 = beginR; jj1 < endR; jj1 ++)
      {
         i1 = R_j[jj1];
         r_entry = R_data[jj1];

        /*-----------------------------------------------------------------
         *  Loop over entries in row i1 of A.
         *-----------------------------------------------------------------*/

         beginA = A_i[i1];
         endA = A_i[i1+1];
         for (jj2 = beginA; jj2 < endA; jj2 ++)
         {
            i2 = A_j[jj2];
            r_a_product = r_entry * A_data[jj2];
            
           /*----------------------------------------------------------------
            *  Check A_marker to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/
             
            beginP = P_i[i2];  
            endP = P_i[i2+1];
            
            if (A_marker[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_marker[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               for (jj3 = beginP; jj3 < endP; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product * P_data[jj3];
                  
                 /*-----------------------------------------------------------
                  *  Check P_marker to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, create a new entry.
                  *  If it has, add new contribution.
                  *----------------------------------------------------------*/

                  if (P_marker[i3] < jj_row_begining)
                  {
                     P_marker[i3]         = jj_counter;
                     RAP_data[jj_counter] = r_a_p_product;
                     RAP_j[jj_counter]    = i3;
                     jj_counter ++;
                  }
                  else
                  {
                     RAP_data[P_marker[i3]] += r_a_p_product;
                  }
               }
            }

           /*--------------------------------------------------------------
            *  If i2 is previously visted ( A_marker[12]=ic ) it yields
            *  no new entries in RAP and can just add new contributions.
            *--------------------------------------------------------------*/

            else
            {
               for (jj3 = beginP; jj3 < endP; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product*P_data[jj3];
                  RAP_data[P_marker[i3]] += r_a_p_product;
               }
            }
         }
      }
   }

   //GetTime(tEnd);
   //jxf_printf(" >>> \033[31mthe 2nd pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd)); 
   
   RAP = fsls_CSRMatrixCreate(n_coarse, n_coarse, RAP_size);
   fsls_CSRMatrixData(RAP) = RAP_data; 
   fsls_CSRMatrixI(RAP)    = RAP_i; 
   fsls_CSRMatrixJ(RAP)    = RAP_j;
   
  /*--------------------------------------------
   *  Free R and marker arrays.
   *-------------------------------------------*/

   fsls_TFree(P_marker);   
   fsls_TFree(A_marker);

   return (RAP);
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_01
 * \brief Compute RAP = R*A*P.
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix 
 * \note Compute 'iac','jac' and 'ac' in only one pass, 
 *       this may cause too many 'if' judgements.
 * \note RAP satisfies the "diagonal-entries-firstly-stored-in-each-row" requirement. 2011/03/29 
 * \author peghoty
 * \date 2011/03/21
 */
fsls_CSRMatrix *
fsls_CSRMatrixTriMultiply_01( fsls_CSRMatrix   *R,
                              fsls_CSRMatrix   *A,
                              fsls_CSRMatrix   *P  )
{
   //struct timeval tStart,tEnd;
   //JXF_Real allocate_time = 0.0;
   
   fsls_CSRMatrix  *RAP;
   
   JXF_Real          *A_data;
   JXF_Int             *A_i;
   JXF_Int             *A_j;

   JXF_Real          *P_data;
   JXF_Int             *P_i;
   JXF_Int             *P_j;

   JXF_Real          *RAP_data;
   JXF_Int             *RAP_i;
   JXF_Int             *RAP_j;

   JXF_Int              RAP_size;

   JXF_Real          *R_data;
   JXF_Int             *R_i;
   JXF_Int             *R_j;

   JXF_Int             *P_marker;
   JXF_Int             *A_marker;

   JXF_Int              n_coarse;
   JXF_Int              n_fine;
   
   JXF_Int              ic, i;
   JXF_Int              i1, i2, i3;
   JXF_Int              jj1, jj2, jj3;
   
   JXF_Int              jj_counter;
   JXF_Int              jj_row_begining;

   JXF_Int              endR, endA, endP;
   
   JXF_Real           r_entry;
   JXF_Real           r_a_product;
   JXF_Real           r_a_p_product;
   
   JXF_Real           zero = 0.0;

  /*----------------------------------------------------
   *  Access the CSR vectors for R, A, P. 
   *  Also get sizes of fine and coarse grids.
   *---------------------------------------------------*/

   R_data = fsls_CSRMatrixData(R);
   R_i    = fsls_CSRMatrixI(R);
   R_j    = fsls_CSRMatrixJ(R);

   A_data = fsls_CSRMatrixData(A);
   A_i    = fsls_CSRMatrixI(A);
   A_j    = fsls_CSRMatrixJ(A);

   P_data = fsls_CSRMatrixData(P);
   P_i    = fsls_CSRMatrixI(P);
   P_j    = fsls_CSRMatrixJ(P);

   n_fine   = fsls_CSRMatrixNumRows(A);
   n_coarse = fsls_CSRMatrixNumRows(R);

  /*-----------------------------------------------------------------------
   *  Allocate RAP_i and marker arrays.
   *----------------------------------------------------------------------*/

   RAP_i    = fsls_CTAlloc(JXF_Int, n_coarse + 1);
   P_marker = fsls_CTAlloc(JXF_Int, n_coarse);
   A_marker = fsls_CTAlloc(JXF_Int, n_fine);
   
  /*-----------------------------------------------------------------------
   *  Pre-allocate RAP_data and RAP_j arrays, here we take 
   *  the number of nonzeros in A for RAP. This usually
   *  proves to be not enough. Therefore, it will be adjusted
   *  later.   peghoty 2011/03/22
   *----------------------------------------------------------------------*/   
   
   RAP_size = fsls_CSRMatrixNumNonzeros(A);   
   RAP_data = fsls_CTAlloc(JXF_Real, RAP_size);
   RAP_j    = fsls_CTAlloc(JXF_Int, RAP_size);

  /*----------------------------------------------
   *  Intialize some stuff.
   *---------------------------------------------*/

   jj_counter = 0;
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -1;
   }
   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -1;
   }   
   
  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < n_coarse; ic ++)
   {
      
     /*----------------------------------------------------------
      *  Create diagonal entry, RAP_{ic,ic}. This ensures 
      *  the diagonal element be first stored. peghoty
      *---------------------------------------------------------*/

      P_marker[ic]    = jj_counter;
      jj_row_begining = jj_counter;
      RAP_data[jj_counter] = zero;
      RAP_j[jj_counter] = ic;
      jj_counter ++;
      
      /* added by peghoty 2011/03/20 */
      if (jj_counter > RAP_size)
      {
         //jxf_printf(" realloc 01!\n");
         RAP_size = 2*RAP_size;
         
         //GetTime(tStart);
         RAP_j    = (JXF_Int*)realloc(RAP_j, RAP_size*sizeof(JXF_Int));
         /* keep in mind that reallocating JXF_Real type will be time-consuming. peghoty */
         RAP_data = (JXF_Real*)realloc(RAP_data, RAP_size*sizeof(JXF_Real)); 
         //GetTime(tEnd);
         //allocate_time += mytime(tStart,tEnd);
      }

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/

      endR = R_i[ic+1];
      for (jj1 = R_i[ic]; jj1 < endR; jj1 ++)
      {
         i1 = R_j[jj1];
         r_entry = R_data[jj1];

        /*--------------------------------------------------
         *  Loop over entries in row i1 of A.
         *-------------------------------------------------*/

         endA = A_i[i1+1];
         for (jj2 = A_i[i1]; jj2 < endA; jj2 ++)
         {
            i2 = A_j[jj2];
            r_a_product = r_entry * A_data[jj2];
            
           /*----------------------------------------------------------------
            *  Check A_marker to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/

            endP = P_i[i2+1];
               
            if (A_marker[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_marker[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               for (jj3 = P_i[i2]; jj3 < endP; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product * P_data[jj3];
                  
                 /*-----------------------------------------------------------
                  *  Check P_marker to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, create a new entry.
                  *  If it has, add new contribution.
                  *----------------------------------------------------------*/

                  if (P_marker[i3] < jj_row_begining)
                  {
                     P_marker[i3] = jj_counter;
                     RAP_data[jj_counter] = r_a_p_product;
                     RAP_j[jj_counter] = i3;
                     jj_counter ++;
                     
                     /* added by peghoty 2011/03/20 */
                     if (jj_counter > RAP_size)
                     {
                        //jxf_printf(" realloc 02!\n");
                        RAP_size = 2*RAP_size;
                        
                        //GetTime(tStart);
                        RAP_j    = (JXF_Int *)realloc(RAP_j, RAP_size*sizeof(JXF_Int));
                        RAP_data = (JXF_Real *)realloc(RAP_data, RAP_size*sizeof(JXF_Real)); 
                        //GetTime(tEnd);
                        //allocate_time += mytime(tStart,tEnd);
                     }
                  }
                  else
                  {
                     RAP_data[P_marker[i3]] += r_a_p_product;
                  }
               }
            }

           /*--------------------------------------------------------------
            *  If i2 is previously visted ( A_marker[12]=ic ) it yields
            *  no new entries in RAP and can just add new contributions.
            *--------------------------------------------------------------*/

            else
            {
               for (jj3 = P_i[i2]; jj3 < endP; jj3 ++)
               {
                  i3 = P_j[jj3];
                  r_a_p_product = r_a_product*P_data[jj3];
                  RAP_data[P_marker[i3]] += r_a_p_product;
               }
            }
         }
      }
      
     /*--------------------------------------------------------------------
      * Set RAP_i for this row.
      *--------------------------------------------------------------------*/

      RAP_i[ic] = jj_row_begining;
   }
   
   RAP_i[n_coarse] = jj_counter; 
   
   /* added by peghoty 2011/03/20 */
   if (RAP_size > jj_counter)
   {
      RAP_size = jj_counter;
      //jxf_printf(" realloc 03!\n");
      
      //GetTime(tStart);
      RAP_j    = (JXF_Int*)realloc(RAP_j, RAP_size*sizeof(JXF_Int));
      RAP_data = (JXF_Real*)realloc(RAP_data, RAP_size*sizeof(JXF_Real)); 
      //GetTime(tEnd);
      //allocate_time += mytime(tStart,tEnd);
   }    
   
   //jxf_printf(" Total reallocate time = %.6lf\n\n", allocate_time); 

   RAP = fsls_CSRMatrixCreate(n_coarse, n_coarse, RAP_size);
   fsls_CSRMatrixData(RAP) = RAP_data; 
   fsls_CSRMatrixI(RAP)    = RAP_i; 
   fsls_CSRMatrixJ(RAP)    = RAP_j;
   
  /*--------------------------------------------
   *  Free marker arrays.
   *-------------------------------------------*/
   fsls_TFree(P_marker);   
   fsls_TFree(A_marker);

   return (RAP);
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_02
 * \brief Compute RAP = R*A*P.
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix 
 * \note Compute 'iac' and 'jac' in the first pass(in FASP fashion), 
 *       then 'ac' in the second pass.
 * \note Numerical results show that the 2nd pass is slower than the one in FASP.
 * \note RAP doesn't necessarily satisfy the "diagonal-entries-firstly-stored-in-each-row"
 *       requirement. 2011/03/29 
 * \author peghoty
 * \date 2011/03/21
 */
fsls_CSRMatrix *
fsls_CSRMatrixTriMultiply_02( fsls_CSRMatrix   *R,
                              fsls_CSRMatrix   *A,
                              fsls_CSRMatrix   *P  )
{
   struct timeval tStart,tEnd;     
   
   fsls_CSRMatrix  *RAP;
   
   JXF_Real          *A_data;
   JXF_Int             *A_i;
   JXF_Int             *A_j;

   JXF_Real          *P_data;
   JXF_Int             *P_i;
   JXF_Int             *P_j;

   JXF_Real          *RAP_data;
   JXF_Int             *RAP_i;
   JXF_Int             *RAP_j;
   
   JXF_Real          *R_data;
   JXF_Int             *R_i;
   JXF_Int             *R_j;

   JXF_Int             *P_marker;
   JXF_Int             *A_marker;
   JXF_Int             *BTindex;

   JXF_Int              n_coarse;
   JXF_Int              n_fine;
   JXF_Int              nzrap;
   
   JXF_Int              ic, i;
   JXF_Int              i1, i2, i3;
   JXF_Int              jj1, jj2, jj3;
   
   JXF_Int              j, jj, jc, k, length;	 
   JXF_Int              istart, iistart, count;
   
   JXF_Int              endB, endR, endA, endP;  

   JXF_Real           r_entry;
   JXF_Real           r_a_product;
   JXF_Real           r_a_p_product;

  /*------------------------------------------------
   *  Access the CSR vectors for R, A, P. 
   *  Also get sizes of fine and coarse grids.
   *-----------------------------------------------*/

   R_data = fsls_CSRMatrixData(R);
   R_i    = fsls_CSRMatrixI(R);
   R_j    = fsls_CSRMatrixJ(R);

   A_data = fsls_CSRMatrixData(A);
   A_i    = fsls_CSRMatrixI(A);
   A_j    = fsls_CSRMatrixJ(A);
   nzrap  = fsls_CSRMatrixNumNonzeros(A);

   P_data = fsls_CSRMatrixData(P);
   P_i    = fsls_CSRMatrixI(P);
   P_j    = fsls_CSRMatrixJ(P);

   n_fine   = fsls_CSRMatrixNumRows(A);
   n_coarse = fsls_CSRMatrixNumRows(R);

  /*--------------------------------------------------
   *  Allocate RAP_i, RAP_j, and marker arrays.
   *-------------------------------------------------*/

   RAP_i    = fsls_CTAlloc(JXF_Int, n_coarse+1);
   RAP_j    = fsls_CTAlloc(JXF_Int, nzrap);
   P_marker = fsls_CTAlloc(JXF_Int, n_coarse);
   A_marker = fsls_CTAlloc(JXF_Int, n_fine);

   GetTime(tStart);     

   //---------------------------------------------------------------------//
   //  First Pass: Determine size of RAP and set up RAP_i and RAP_j       //
   //---------------------------------------------------------------------//

  /*-----------------------------------
   *  Intialize some stuff.
   *----------------------------------*/

   for (i = 0; i < n_fine; i ++)
   {      
      A_marker[i] = -2;
   }
   
   //memcpy(P_marker, A_marker, n_coarse*sizeof(JXF_Int));
   for (ic = 0; ic < n_coarse; ic ++)
   {      
      P_marker[ic] = -2;
   } 


  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
   
   RAP_i[0] = 0;
   
   for (i = 0; i < n_coarse; i ++)
   {
      // reset istart and length at the begining of each loop
      istart = -1;
      length = 0;
      i1     = i + 1;
      
      // go across the rows in R
      endR = R_i[i1];
      for (jj = R_i[i]; jj < endR; ++ jj)
      {
         j    = R_j[jj];
         endA = A_i[j+1];			
         for (k = A_i[j]; k < endA; ++ k) // for each column in A
         {
            jc = A_j[k];
            if (A_marker[jc] == -2)
            {
               A_marker[jc] = istart;
               istart = jc;
               ++ length;
            }
         }
      }    
		
      // book-keeping [reseting length and setting iistart]
      count   = length; 
      iistart = -1; 
      length  = 0;
		
      // use each column that would have resulted from R*A
      for (j = 0; j < count; ++ j)
      {
         jj           = istart;
         istart       = A_marker[istart];
         A_marker[jj] = -2;
         
         endP = P_i[jj+1];
         for (k = P_i[jj]; k < endP; ++ k) // go across the row of P
         {
            jc = P_j[k];
            // pull out the appropriate columns of P
            if (P_marker[jc] == -2)
            {
               P_marker[jc] = iistart;
               iistart = jc;
               ++ length;
            }
         } // end for k
      } // end for j
		
      // set RAP->IA
      RAP_i[i1] = RAP_i[i] + length;
		
      if (RAP_i[i1] > nzrap) 
      {
         nzrap = nzrap*2;
         RAP_j = (JXF_Int*)realloc(RAP_j, nzrap*sizeof(JXF_Int));
         //jxf_printf(" realloc 02!\n");
      }
		
      // put the correct columns of p into the column list of the products
      endB = RAP_i[i1];
      for (j = RAP_i[i]; j < endB; ++ j)
      {
         RAP_j[j] = iistart;           // put the value in B->JA
         iistart = P_marker[iistart];  // set istart to the next value
         P_marker[RAP_j[j]] = -2;      // set the P_marker spot to 0
      }  // end j
   }
 
   GetTime(tEnd);
   jxf_printf(" >>> \033[31mthe 1st pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd));    
 
  /*-----------------------------------------------------------------------
   *  Reallocate RAP_j, and allocate RAP_data arrays.
   *----------------------------------------------------------------------*/
   nzrap    = RAP_i[n_coarse];
   RAP_j    = (JXF_Int*)realloc(RAP_j, nzrap*sizeof(JXF_Int));
   //jxf_printf(" realloc 03!\n");
   RAP_data = fsls_CTAlloc(JXF_Real, nzrap);
   BTindex  = fsls_CTAlloc(JXF_Int, n_coarse);


   GetTime(tStart);

   //---------------------------------------------------------------------//
   //  Second Pass: Fill in RAP_data and RAP_j.                           //
   //---------------------------------------------------------------------//
   
  /*--------------------------------------------------
   *  Loop over c-points.
   *-------------------------------------------------*/

   for (ic = 0; ic < n_coarse; ic ++)
   {
      endB = RAP_i[ic+1];		
      for (k = RAP_i[ic]; k < endB; ++ k) 
      {
         BTindex[RAP_j[k]] = k;
      }
		 
     /*----------------------------------------------
      *  Loop over entries in row ic of R.
      *---------------------------------------------*/
      
      endR = R_i[ic+1];	
      for (jj1 = R_i[ic]; jj1 < endR; jj1 ++)
      {
         i1 = R_j[jj1];
         r_entry = R_data[jj1];

        /*----------------------------------------------
         *  Loop over entries in row i1 of A.
         *---------------------------------------------*/

         endA = A_i[i1+1];
         for (jj2 = A_i[i1]; jj2 < endA; jj2 ++)
         {
            i2 = A_j[jj2];
            r_a_product = r_entry*A_data[jj2];

           /*----------------------------------------------
            *  Loop over entries in row i2 of P.
            *---------------------------------------------*/

            endP = P_i[i2+1];
            for (jj3 = P_i[i2]; jj3 < endP; jj3 ++)
            {
               i3 = P_j[jj3];
               r_a_p_product = r_a_product * P_data[jj3];
               RAP_data[BTindex[i3]] += r_a_p_product;
            }
         }
      }
   }

   GetTime(tEnd);
   jxf_printf(" >>> \033[31mthe 2nd pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd));  

   RAP = fsls_CSRMatrixCreate(n_coarse, n_coarse, nzrap);
   fsls_CSRMatrixData(RAP) = RAP_data; 
   fsls_CSRMatrixI(RAP)    = RAP_i; 
   fsls_CSRMatrixJ(RAP)    = RAP_j;
   
  /*--------------------------------------------
   *  Free R and marker arrays.
   *-------------------------------------------*/

   fsls_TFree(P_marker);   
   fsls_TFree(A_marker);
   fsls_TFree(BTindex);

   return (RAP);
}

/**
 * \fn fsls_CSRMatrix *fasp_CSRMatrixTriMultiply
 * \brief Triple sparse matrix multiplication Ac = R*A*P
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \note RAP doesn't necessarily satisfy the "diagonal-entries-firstly-stored-in-each-row"
 *       requirement. 2011/03/29 
 * \author Xuehai Huang, Chensong Zhang, modified by peghoty
 * \date 2011/03/17
 */
fsls_CSRMatrix *
fasp_CSRMatrixTriMultiply( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P )
{
   //struct timeval tStart,tEnd;
   
   JXF_Int  row, col, acol, nB;
   JXF_Int  i, i1, j, jj, k, length;	
   JXF_Int  begin_row, end_row;
   JXF_Int  begin_rowA, end_rowA;
   JXF_Int  begin_rowR, end_rowR;
   JXF_Int  istart, iistart, count;

   JXF_Real *rj, *aj, *pj, *acj;
   JXF_Int    *ir, *ia, *ip, *iac;
   JXF_Int    *jr, *ja, *jp, *jac;

   JXF_Int    *index   = NULL;
   JXF_Int    *iindex  = NULL;	
   JXF_Real *temp    = NULL;
   JXF_Int    *BTindex = NULL;

   fsls_CSRMatrix *RAP = NULL; 

   row  = fsls_CSRMatrixNumRows(R);
   col  = fsls_CSRMatrixNumCols(P);
   acol = fsls_CSRMatrixNumCols(A);
   nB   = fsls_CSRMatrixNumNonzeros(A);

   rj = fsls_CSRMatrixData(R);
   ir = fsls_CSRMatrixI(R);
   jr = fsls_CSRMatrixJ(R);

   aj = fsls_CSRMatrixData(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);

   pj = fsls_CSRMatrixData(P);
   ip = fsls_CSRMatrixI(P);
   jp = fsls_CSRMatrixJ(P);

   jac  = fsls_CTAlloc(JXF_Int, nB);	
   iac  = fsls_CTAlloc(JXF_Int, row+1);

   index  = fsls_CTAlloc(JXF_Int, acol);
   iindex = fsls_CTAlloc(JXF_Int, col);
   temp = fsls_CTAlloc(JXF_Real, acol);	   
    
    
   //GetTime(tStart); 
    
   for (i = 0; i < acol; ++ i) 
   {
      index[i] = -2;
   }
	
   memcpy(iindex, index, col*sizeof(JXF_Int));

   //----------------------------------------------//
   //  First loop: form sparsity partern of R*A*P  //
   //----------------------------------------------//
   
   iac[0] = 0;
   
   for (i = 0; i < row; ++ i) 
   {		
      // reset istart and length at the begining of each loop
      istart = -1; 
      length = 0; 
      i1     = i + 1;
		
      // go across the rows in R
      begin_rowR = ir[i]; 
      end_rowR   = ir[i1];
      for (jj = begin_rowR; jj < end_rowR; ++ jj)
      {
         j = jr[jj];
         // for each column in A
         begin_rowA = ia[j]; 
         end_rowA   = ia[j+1];			
         for (k = begin_rowA; k < end_rowA; ++ k)
         {
            if (index[ja[k]] == -2)
            {
               index[ja[k]] = istart;
               istart = ja[k];
               ++ length;
            }
         }
      }    
		
      // book-keeping [reseting length and setting iistart]
      count   = length; 
      iistart = -1; 
      length  = 0;
		
      // use each column that would have resulted from R*A
      for (j = 0; j < count; ++ j)
      {
         jj        = istart;
         istart    = index[istart];
         index[jj] = -2;
			
         // go across the row of P
         begin_row = ip[jj]; 
         end_row   = ip[jj+1];
         for (k = begin_row; k < end_row; ++ k)
         {
            // pull out the appropriate columns of P
            if (iindex[jp[k]] == -2)
            {
               iindex[jp[k]] = iistart;
               iistart = jp[k];
               ++ length;
            }
         } // end for k
      } // end for j
		
      // set B->IA
      iac[i1] = iac[i] + length;
		
      if (iac[i1] > nB) 
      {
         nB  = nB*2;
         jac = (JXF_Int*)realloc(jac, nB*sizeof(JXF_Int));
      }
		
      // put the correct columns of p into the column list of the products
      begin_row = iac[i]; 
      end_row   = iac[i1];
      for (j = begin_row; j < end_row; ++ j)
      {
         jac[j] = iistart;           // put the value in B->JA
         iistart = iindex[iistart];  // set istart to the next value
         iindex[jac[j]] = -2;        // set the iindex spot to 0
      } // end j
		
   } // end i: First loop
   
   //GetTime(tEnd);
   //jxf_printf(" >>> \033[31mthe 1st pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd));   
   
		
   jac = (JXF_Int*)realloc(jac, (iac[row])*sizeof(JXF_Int));	
   acj = fsls_CTAlloc(JXF_Real, iac[row]);
   BTindex = fsls_CTAlloc(JXF_Int, col);

   //GetTime(tStart); 

   //----------------------------------------------//
   //  Second loop: compute entries of R*A*P       //
   //----------------------------------------------//	

   for (i = 0; i < row; ++ i) 
   {
      i1 = i + 1;
		
      // each col of B
      begin_row = iac[i]; 
      end_row   = iac[i1];		
      for (j = begin_row; j < end_row; ++ j) 
      {
         BTindex[jac[j]] = j;
      }
		
      // reset istart and length at the begining of each loop
      istart = -1; 
      length = 0;
		
      // go across the rows in R
      begin_rowR = ir[i]; 
      end_rowR   = ir[i1];		
      for (jj = begin_rowR; jj < end_rowR; ++ jj)
      {
         j = jr[jj];
			
         // for each column in A
         begin_rowA = ia[j]; 
         end_rowA   = ia[j+1];					
         for (k = begin_rowA; k < end_rowA; ++ k)
         {
            if (index[ja[k]] == -2)
            {
               index[ja[k]] = istart;
               istart = ja[k];
               ++ length;
            }
            temp[ja[k]] += rj[jj]*aj[k];
         }
      } 
		
      // book-keeping [reseting length and setting iistart]
      // use each column that would have resulted from R*A
      for (j = 0; j < length; ++ j) 
      {
         jj = istart;
         istart = index[istart];
         index[jj] = -2;
			
         // go across the row of P
         begin_row = ip[jj]; 
         end_row   = ip[jj+1];		
         for (k = begin_row; k < end_row; ++ k) 
         {
            acj[BTindex[jp[k]]] += temp[jj]*pj[k];  // pull out the appropriate columns of P
         }
         temp[jj] = 0.0;
      }
		
   } // end for i: Second loop

   //GetTime(tEnd);
   //jxf_printf(" >>> \033[31mthe 2nd pass = \033[00m %.3lf seconds\n", mytime(tStart,tEnd)); 

   // setup coarse matrix RAP
   RAP = fsls_CSRMatrixCreate(row, col, iac[row] - iac[0]);
   fsls_CSRMatrixData(RAP) = acj; 
   fsls_CSRMatrixI(RAP) = iac; 
   fsls_CSRMatrixJ(RAP) = jac;
   
   // Free R and marker arrays 
   fsls_TFree(temp);
   fsls_TFree(index);
   fsls_TFree(iindex);
   fsls_TFree(BTindex);
   
   return (RAP);
}

/**
 * \fn fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_improve
 * \brief Triple sparse matrix multiplication Ac = R*A*P
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Xuehai Huang, Chensong Zhang, modified by peghoty
 * \note Compute 'iac','jac' and 'ac' in one pass. Compared with the integer type 
 *       auxiliary array 'BTindex' of size nc, a JXF_Real type auxiliary array 'P_temp' of 
 *       size nc will be needed. 
 *
 * \note RAP doesn't necessarily satisfy the "diagonal-entries-firstly-stored-in-each-row"
 *       requirement. 2011/03/29  
 * \date 2011/03/23
 */
fsls_CSRMatrix *
fasp_CSRMatrixTriMultiply_improve( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P )
{   
   JXF_Int  row, col, acol, RAP_size;
   JXF_Int  i, i1, j, jj, k, length;
   JXF_Int  begin_row, begin_rowA, begin_rowR;	
   JXF_Int  end_row, end_rowA, end_rowR;
   JXF_Int  istart, iistart, count;

   JXF_Real *rj, *aj, *pj, *acj;
   JXF_Int    *ir, *ia, *ip, *iac;
   JXF_Int    *jr, *ja, *jp, *jac;

   JXF_Int    *A_mark = NULL;
   JXF_Int    *P_mark = NULL;	
   JXF_Real *A_temp = NULL;
   JXF_Real *P_temp = NULL;

   fsls_CSRMatrix *RAP = NULL; 

   row  = fsls_CSRMatrixNumRows(R);
   col  = fsls_CSRMatrixNumCols(P);
   acol = fsls_CSRMatrixNumCols(A);
   RAP_size = fsls_CSRMatrixNumNonzeros(A);

   rj = fsls_CSRMatrixData(R);
   ir = fsls_CSRMatrixI(R);
   jr = fsls_CSRMatrixJ(R);

   aj = fsls_CSRMatrixData(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);

   pj = fsls_CSRMatrixData(P);
   ip = fsls_CSRMatrixI(P);
   jp = fsls_CSRMatrixJ(P);
	
   iac  = fsls_CTAlloc(JXF_Int, row+1);
   jac  = fsls_CTAlloc(JXF_Int, RAP_size);    // pre-allocate with size of A, will be changed later if not enough
   acj  = fsls_CTAlloc(JXF_Real, RAP_size); // pre-allocate with size of A, will be changed later if not enough

   A_mark = fsls_CTAlloc(JXF_Int, acol);
   P_mark = fsls_CTAlloc(JXF_Int, col);
   A_temp = fsls_CTAlloc(JXF_Real, acol);	 
   P_temp = fsls_CTAlloc(JXF_Real, col);  
    
   for (i = 0; i < acol; ++ i) 
   {
      A_mark[i] = -2;
   }
	
   memcpy(P_mark, A_mark, col*sizeof(JXF_Int));

   //----------------------------------------------//
   //  Form sparsity partern of R*A*P  and         //
   //  compute all nonzeros in ONE LOOP.           //
   //----------------------------------------------//
   
   iac[0] = 0;
   
   for (i = 0; i < row; ++ i) 
   {		
      // reset istart and length at the begining of each loop
      istart = -1; 
      length = 0; 
      i1     = i + 1;
		
      // go across the rows in R
      begin_rowR = ir[i];
      end_rowR   = ir[i1];
      for (jj = begin_rowR; jj < end_rowR; ++ jj)
      {
         j = jr[jj];
			
         // for each column in A
         begin_rowA = ia[j];
         end_rowA   = ia[j+1];			
         for (k = begin_rowA; k < end_rowA; ++ k)
         {
            if (A_mark[ja[k]] == -2)
            {
               A_mark[ja[k]] = istart;
               istart = ja[k];
               ++ length;
            }
            A_temp[ja[k]] += rj[jj]*aj[k];
         }
      }    
		
      // book-keeping [reseting length and setting iistart]
      count   = length; 
      iistart = -1; 
      length  = 0;
		
      // use each column that would have resulted from R*A
      for (j = 0; j < count; ++ j)
      {
         jj         = istart;
         istart     = A_mark[istart];
         A_mark[jj] = -2;
			
         // go across the row of P
         begin_row = ip[jj];
         end_row   = ip[jj+1];
         for (k = begin_row; k < end_row; ++ k)
         {
            // pull out the appropriate columns of P
            if (P_mark[jp[k]] == -2)
            {
               P_mark[jp[k]] = iistart;
               iistart = jp[k];
               ++ length;
            }
            P_temp[jp[k]] += A_temp[jj]*pj[k];
         } // end for k
         A_temp[jj] = 0.0;
      } // end for j
		
      // set B->IA
      iac[i1] = iac[i] + length;
		
      if (iac[i1] > RAP_size) 
      {
         RAP_size  = RAP_size*2;
         jac = (JXF_Int *)realloc(jac, RAP_size*sizeof(JXF_Int));
         acj = (JXF_Real *)realloc(acj, RAP_size*sizeof(JXF_Real));
      }
		
      // put the correct columns of p into the column list of the products
      begin_row = iac[i];
      end_row   = iac[i1];
      for (j = begin_row; j < end_row; ++ j)
      {
         jac[j] = iistart;           // put the value in B->JA
         acj[j] = P_temp[iistart];   
         iistart = P_mark[iistart];  // set iistart to the next value
         P_mark[jac[j]] = -2;        // set the P_mark spot to 0
         P_temp[jac[j]] = 0.0;      // set the P_temp spot to 0         	
      } // end j
		
   } // end i
		
   jac = (JXF_Int *)realloc(jac, (iac[row])*sizeof(JXF_Int));	
   acj = (JXF_Real *)realloc(acj, (iac[row])*sizeof(JXF_Real));

   // setup coarse matrix RAP
   RAP = fsls_CSRMatrixCreate(row, col, iac[row] - iac[0]);
   fsls_CSRMatrixData(RAP) = acj; 
   fsls_CSRMatrixI(RAP)    = iac; 
   fsls_CSRMatrixJ(RAP)    = jac;
   
   // Free R and marker arrays 
   fsls_TFree(A_temp);
   fsls_TFree(P_temp);
   fsls_TFree(A_mark);
   fsls_TFree(P_mark);
   
   return (RAP);
}

/**
 * \fn fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_improve2
 * \brief Triple sparse matrix multiplication Ac = R*A*P
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Xuehai Huang, Chensong Zhang, modified by peghoty
 * \note Compute 'iac','jac' and 'ac' in one pass. Compared with the integer type 
 *       auxiliary array 'BTindex' of size nc, a JXF_Real type auxiliary array 'P_temp' of 
 *       size nc will be needed. 
 * \note RAP doesn't necessarily satisfy the "diagonal-entries-firstly-stored-in-each-row"
 *       requirement. 2011/03/29  
 * \note Pre-allocate 2*nz(A) for 'ac' and 'jac'.
 * \date 2011/03/24
 */
fsls_CSRMatrix *
fasp_CSRMatrixTriMultiply_improve2( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P )
{   
   JXF_Int  row, col, acol, RAP_size;
   JXF_Int  i, i1, j, jj, k, length;
   JXF_Int  begin_row, begin_rowA, begin_rowR;	
   JXF_Int  end_row, end_rowA, end_rowR;
   JXF_Int  istart, iistart, count;

   JXF_Real *rj, *aj, *pj, *acj;
   JXF_Int    *ir, *ia, *ip, *iac;
   JXF_Int    *jr, *ja, *jp, *jac;

   JXF_Int    *A_mark = NULL;
   JXF_Int    *P_mark = NULL;	
   JXF_Real *A_temp = NULL;
   JXF_Real *P_temp = NULL;

   fsls_CSRMatrix *RAP = NULL; 

   row  = fsls_CSRMatrixNumRows(R);
   col  = fsls_CSRMatrixNumCols(P);
   acol = fsls_CSRMatrixNumCols(A);
   RAP_size = fsls_CSRMatrixNumNonzeros(A);

   rj = fsls_CSRMatrixData(R);
   ir = fsls_CSRMatrixI(R);
   jr = fsls_CSRMatrixJ(R);

   aj = fsls_CSRMatrixData(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);

   pj = fsls_CSRMatrixData(P);
   ip = fsls_CSRMatrixI(P);
   jp = fsls_CSRMatrixJ(P);
	
   iac  = fsls_CTAlloc(JXF_Int, row+1);
   jac  = fsls_CTAlloc(JXF_Int, 2*RAP_size);    
   acj  = fsls_CTAlloc(JXF_Real, 2*RAP_size);

   A_mark = fsls_CTAlloc(JXF_Int, acol);
   P_mark = fsls_CTAlloc(JXF_Int, col);
   A_temp = fsls_CTAlloc(JXF_Real, acol);	 
   P_temp = fsls_CTAlloc(JXF_Real, col);  
    
   for (i = 0; i < acol; ++ i) 
   {
      A_mark[i] = -2;
   }
	
   memcpy(P_mark, A_mark, col*sizeof(JXF_Int));

   //----------------------------------------------//
   //  Form sparsity partern of R*A*P  and         //
   //  compute all nonzeros in ONE LOOP.           //
   //----------------------------------------------//
   
   iac[0] = 0;
   
   for (i = 0; i < row; ++ i) 
   {		
      // reset istart and length at the begining of each loop
      istart = -1; 
      length = 0; 
      i1     = i + 1;
		
      // go across the rows in R
      begin_rowR = ir[i];
      end_rowR   = ir[i1];
      for (jj = begin_rowR; jj < end_rowR; ++ jj)
      {
         j = jr[jj];
			
         // for each column in A
         begin_rowA = ia[j];
         end_rowA   = ia[j+1];			
         for (k = begin_rowA; k < end_rowA; ++ k)
         {
            if (A_mark[ja[k]] == -2)
            {
               A_mark[ja[k]] = istart;
               istart = ja[k];
               ++ length;
            }
            A_temp[ja[k]] += rj[jj]*aj[k];
         }
      }    
		
      // book-keeping [reseting length and setting iistart]
      count   = length; 
      iistart = -1; 
      length  = 0;
		
      // use each column that would have resulted from R*A
      for (j = 0; j < count; ++ j)
      {
         jj         = istart;
         istart     = A_mark[istart];
         A_mark[jj] = -2;
			
         // go across the row of P
         begin_row = ip[jj];
         end_row   = ip[jj+1];
         for (k = begin_row; k < end_row; ++ k)
         {
            // pull out the appropriate columns of P
            if (P_mark[jp[k]] == -2)
            {
               P_mark[jp[k]] = iistart;
               iistart = jp[k];
               ++ length;
            }
            P_temp[jp[k]] += A_temp[jj]*pj[k];
         } // end for k
         A_temp[jj] = 0.0;
      } // end for j
		
      // set B->IA
      iac[i1] = iac[i] + length;
		
      if (iac[i1] > RAP_size) 
      {
         RAP_size  = RAP_size*2;
         jac = (JXF_Int *)realloc(jac, RAP_size*sizeof(JXF_Int));
         acj = (JXF_Real *)realloc(acj, RAP_size*sizeof(JXF_Real));
      }
		
      // put the correct columns of p into the column list of the products
      begin_row = iac[i];
      end_row   = iac[i1];
      for (j = begin_row; j < end_row; ++ j)
      {
         jac[j] = iistart;           // put the value in B->JA
         acj[j] = P_temp[iistart];   
         iistart = P_mark[iistart];  // set iistart to the next value
         P_mark[jac[j]] = -2;        // set the P_mark spot to 0
         P_temp[jac[j]] = 0.0;      // set the P_temp spot to 0         	
      } // end j
		
   } // end i
		
   jac = (JXF_Int *)realloc(jac, (iac[row])*sizeof(JXF_Int));	
   acj = (JXF_Real *)realloc(acj, (iac[row])*sizeof(JXF_Real));

   // setup coarse matrix RAP
   RAP = fsls_CSRMatrixCreate(row, col, iac[row] - iac[0]);
   fsls_CSRMatrixData(RAP) = acj; 
   fsls_CSRMatrixI(RAP)    = iac; 
   fsls_CSRMatrixJ(RAP)    = jac;
   
   // Free R and marker arrays 
   fsls_TFree(A_temp);
   fsls_TFree(P_temp);
   fsls_TFree(A_mark);
   fsls_TFree(P_mark);
   
   return (RAP);
}

/**
 * \fn fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_mixed
 * \brief Triple sparse matrix multiplication Ac = R*A*P
 * \param *R pointer to the fsls_CSRMatrix matrix
 * \param *A pointer to the fsls_CSRMatrix matrix
 * \param *P pointer to the fsls_CSRMatrix matrix
 * \return pointer to the resulting matrix
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Xuehai Huang, Chensong Zhang, modified by peghoty
 * \note Compute 'iac' in HYPRE fasion, and 'jac','ac' in FASP fasion.
 * \note RAP doesn't necessarily satisfy the "diagonal-entries-firstly-stored-in-each-row"
 *       requirement. 2011/03/29  
 * \date 2011/03/23
 */
fsls_CSRMatrix *
fasp_CSRMatrixTriMultiply_mixed( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P )
{   
   JXF_Int  row, col, acol, RAP_size;
   JXF_Int  i, i1, j, jj, k, length;
   JXF_Int  begin_row, begin_rowA, begin_rowR;	
   JXF_Int  end_row, end_rowA, end_rowR;
   JXF_Int  istart, iistart, count;
   
   JXF_Int  ic, i2, i3, jj1, jj2, jj3;
   JXF_Int  jj_row_begining, jj_counter;

   JXF_Real *rj, *aj, *pj, *acj;
   JXF_Int    *ir, *ia, *ip, *iac;
   JXF_Int    *jr, *ja, *jp, *jac;

   JXF_Int    *A_mark = NULL;
   JXF_Int    *P_mark = NULL;	
   JXF_Real *A_temp = NULL;
   JXF_Real *P_temp = NULL;

   fsls_CSRMatrix *RAP = NULL; 

   row  = fsls_CSRMatrixNumRows(R);
   col  = fsls_CSRMatrixNumCols(P);
   acol = fsls_CSRMatrixNumCols(A);
   RAP_size = fsls_CSRMatrixNumNonzeros(A);

   rj = fsls_CSRMatrixData(R);
   ir = fsls_CSRMatrixI(R);
   jr = fsls_CSRMatrixJ(R);

   aj = fsls_CSRMatrixData(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);

   pj = fsls_CSRMatrixData(P);
   ip = fsls_CSRMatrixI(P);
   jp = fsls_CSRMatrixJ(P);
	
   iac  = fsls_CTAlloc(JXF_Int, row+1);

   A_mark = fsls_CTAlloc(JXF_Int, acol);
   P_mark = fsls_CTAlloc(JXF_Int, col);
   A_temp = fsls_CTAlloc(JXF_Real, acol);	 
   P_temp = fsls_CTAlloc(JXF_Real, col);  
    

   //----------------------------------------------//
   //  First Loop: Form sparsity partern of R*A*P  //
   //----------------------------------------------// 
   
   jj_counter = 0;
   for (i = 0; i < acol; ++ i) 
   {
      A_mark[i] = -2;
   }
   //for (i = 0; i < col; ++ i) 
   //{
   //   P_mark[i] = -2;
   //}
   memcpy(P_mark, A_mark, col*sizeof(JXF_Int)); // Though fast, this will be dangerous if acol < col. peghoty 

  /*-----------------------------------------------------------------------
   *  Loop over c-points.
   *----------------------------------------------------------------------*/
    
   for (ic = 0; ic < row; ic ++)
   {
      
     /*--------------------------------------------------------------------
      *  Set marker for diagonal entry, RAP_{ic,ic}.
      *-------------------------------------------------------------------*/

      P_mark[ic]      = jj_counter;
      jj_row_begining = jj_counter;
      jj_counter ++;

     /*--------------------------------------------------------------------
      *  Loop over entries in row ic of R.
      *-------------------------------------------------------------------*/
      
      begin_rowR = ir[ic];
      end_rowR = ir[ic+1];
      for (jj1 = begin_rowR; jj1 < end_rowR; jj1 ++)
      {
         i1 = jr[jj1];

        /*-----------------------------------------------------------------
         *  Loop over entries in row i1 of A.
         *----------------------------------------------------------------*/

         begin_rowA = ia[i1];
         end_rowA = ia[i1+1];
         for (jj2 = begin_rowA; jj2 < end_rowA; jj2 ++)
         {
            i2 = ja[jj2];

           /*----------------------------------------------------------------
            *  Check A_mark to see if point i2 has been previously
            *  visited. New entries in RAP only occur from unmarked points.
            *---------------------------------------------------------------*/

            if (A_mark[i2] != ic)
            {

              /*-----------------------------------------------------------
               *  Mark i2 as visited.
               *----------------------------------------------------------*/

               A_mark[i2] = ic;
               
              /*-----------------------------------------------------------
               *  Loop over entries in row i2 of P.
               *----------------------------------------------------------*/

               begin_row = ip[i2];
               end_row = ip[i2+1];
               for (jj3 = begin_row; jj3 < end_row; jj3 ++)
               {
                  i3 = jp[jj3];
                  
                 /*------------------------------------------------------------
                  *  Check P_mark to see that RAP_{ic,i3} has not already
                  *  been accounted for. If it has not, mark it and increment
                  *  counter.
                  *-----------------------------------------------------------*/

                  if (P_mark[i3] < jj_row_begining)
                  {
                     P_mark[i3] = jj_counter;
                     jj_counter ++;
                  }
               }
            }
         }
      }
            
     /*--------------------------------------------------------------------
      * Set RAP_i for this row.
      *--------------------------------------------------------------------*/

      iac[ic] = jj_row_begining;
      
   }
  
   iac[row] = jj_counter; 
   
   RAP_size = jj_counter;
   jac = fsls_CTAlloc(JXF_Int, RAP_size);    
   acj = fsls_CTAlloc(JXF_Real, RAP_size);

   //----------------------------------------------------------------//
   //  Second pass:  compute all nonzeros and their column numbers   //
   //----------------------------------------------------------------//

   for (i = 0; i < acol; ++ i) 
   {
      A_mark[i] = -2;
   }
   //for (i = 0; i < col; ++ i) 
   //{
   //   P_mark[i] = -2;
   //}
   memcpy(P_mark, A_mark, col*sizeof(JXF_Int)); // Though fast, this will be dangerous if acol < col. peghoty 
   
   
   for (i = 0; i < row; ++ i) 
   {		
      // reset istart and length at the begining of each loop
      istart = -1; 
      length = 0; 
      i1     = i + 1;
		
      // go across the rows in R
      begin_rowR = ir[i];
      end_rowR   = ir[i1];
      for (jj = begin_rowR; jj < end_rowR; ++ jj)
      {
         j = jr[jj];
			
         // for each column in A
         begin_rowA = ia[j];
         end_rowA   = ia[j+1];			
         for (k = begin_rowA; k < end_rowA; ++ k)
         {
            if (A_mark[ja[k]] == -2)
            {
               A_mark[ja[k]] = istart;
               istart = ja[k];
               ++ length;
            }
            A_temp[ja[k]] += rj[jj]*aj[k];
         }
      }    
		
      // book-keeping [reseting length and setting iistart]
      count   = length; 
      iistart = -1; 
      length  = 0;
		
      // use each column that would have resulted from R*A
      for (j = 0; j < count; ++ j)
      {
         jj         = istart;
         istart     = A_mark[istart];
         A_mark[jj] = -2;
			
         // go across the row of P
         begin_row = ip[jj];
         end_row   = ip[jj+1];
         for (k = begin_row; k < end_row; ++ k)
         {
            // pull out the appropriate columns of P
            if (P_mark[jp[k]] == -2)
            {
               P_mark[jp[k]] = iistart;
               iistart = jp[k];
               ++ length;
            }
            P_temp[jp[k]] += A_temp[jj]*pj[k];
         } // end for k
         A_temp[jj] = 0.0;
      } // end for j
		
      // put the correct columns of p into the column list of the products
      begin_row = iac[i];
      end_row   = iac[i1];
      for (j = begin_row; j < end_row; ++ j)
      {
         jac[j] = iistart;           // put the value in B->JA
         acj[j] = P_temp[iistart];   
         iistart = P_mark[iistart];  // set iistart to the next value
         P_mark[jac[j]] = -2;        // set the P_mark spot to 0
         P_temp[jac[j]] = 0.0;       // set the P_temp spot to 0         	
      } // end j
		
   } // end i

   // setup coarse matrix RAP
   RAP = fsls_CSRMatrixCreate(row, col, RAP_size);
   fsls_CSRMatrixData(RAP) = acj; 
   fsls_CSRMatrixI(RAP)    = iac; 
   fsls_CSRMatrixJ(RAP)    = jac;
   
   // Free R and marker arrays 
   fsls_TFree(A_temp);
   fsls_TFree(P_temp);
   fsls_TFree(A_mark);
   fsls_TFree(P_mark);
   
   return (RAP);
}

/*!
 * \fn fsls_CSRMatrix * fsls_CSRMatrixDeleteZeros
 * \brief get rid of zeros (here, zero means the element which 
 *        satisfies "|a_{ij}|<=tol") from the matrix A.
 * \param *A the pointer to the matrix which possibly contains zeros.
 * \param tol the tolerance to define zero
 * \return the pointer to the matrix after processing
 * \author peghoty
 * \date 2009/12/05
 */
fsls_CSRMatrix *
fsls_CSRMatrixDeleteZeros( fsls_CSRMatrix *A, JXF_Real tol )
{
   JXF_Real *A_data       = fsls_CSRMatrixData(A);
   JXF_Int    *A_i          = fsls_CSRMatrixI(A);
   JXF_Int    *A_j          = fsls_CSRMatrixJ(A);
   JXF_Int     nrows_A      = fsls_CSRMatrixNumRows(A);
   JXF_Int     ncols_A      = fsls_CSRMatrixNumCols(A);
   JXF_Int     num_nonzeros = fsls_CSRMatrixNumNonzeros(A);

   fsls_CSRMatrix *B      = NULL;
   JXF_Real         *B_data = NULL; 
   JXF_Int            *B_i    = NULL;
   JXF_Int            *B_j    = NULL;

   JXF_Int zeros;
   JXF_Int i,j;
   JXF_Int nzB;

   /* get the total number of zeros in matrix A */
   zeros = 0;
   for (i = 0; i < num_nonzeros; i ++)
   {
      if (fabs(A_data[i]) <= tol) zeros ++;
   }

   /* there exists zeros in the matrix A */
   if (zeros)
   {
      B = fsls_CSRMatrixCreate(nrows_A,ncols_A,num_nonzeros-zeros);
      fsls_CSRMatrixInitialize(B);
      B_i = fsls_CSRMatrixI(B);
      B_j = fsls_CSRMatrixJ(B);
      B_data = fsls_CSRMatrixData(B);
      B_i[0] = 0;
      nzB = 0;
      for (i = 0; i < nrows_A; i ++)
      {  
         for (j = A_i[i]; j < A_i[i+1]; j ++)
         {  /* modified by peghoty 2009/12/06 */
            if (fabs(A_data[j]) > tol)
            {
               B_data[nzB] = A_data[j];
               B_j[nzB]    = A_j[j];
               nzB ++;
            }
         }
         B_i[i+1] = nzB;
      }
      fsls_CSRMatrixDestroy(A); /* added by peghoty 2009/12/06 */
      A = NULL;
      return B;
   }
   else /* there doesn't exist zeros in the matrix A */
   {
      return A; /* modified by peghoty 2009/12/06 */
   }
}

/*!
 * \fn fsls_CSRMatrix * fsls_CSRMatrixSymmetrization
 * \brief Symmetrize a given matrix A using "(A+A^T)/2".
 * \param *A pointer to the matrix to be symmetrized.
 * \note We still keep the original matrix A after symmetrizing.
 * \author peghoty
 * \date 2010/01/11  
 */
fsls_CSRMatrix *
fsls_CSRMatrixSymmetrization( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix *AT = NULL;
   fsls_CSRMatrix *B  = NULL;

   fsls_CSRMatrixTranspose(A, &AT, 1);
   B = fsls_CSRMatrixAdd(A, AT);
   fsls_CSRMatrixScale(B, 0.5);
   
   fsls_CSRMatrixDestroy(AT);

   return(B);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixScale
 * \brief Scale a given matrix by multiplying each entries with alpha.
 * \param *A pointer to the matrix to be scaled.
 * \param alpha the scaler
 * \author peghoty 
 * \date 2010/01/11 
 */
JXF_Int 
fsls_CSRMatrixScale( fsls_CSRMatrix *A, JXF_Real alpha )
{
   JXF_Int i,nz;

   nz = fsls_CSRMatrixNumNonzeros(A);
   for (i = 0; i < nz; i ++)
   {
      fsls_CSRMatrixData(A)[i] *= alpha; 
   } 

   return(0);
}

/*!
 * \fn JXF_Real fsls_CSRMatrixFrobeniusNorm
 * \brief Compute the Frobenius norm of a given matrix
 *        Frobenius norm of A = (\sum\sum_a_{ij}^2)^{1/2}
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/01/15  
 */
JXF_Real 
fsls_CSRMatrixFrobeniusNorm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real tmp;
   JXF_Real fnorm;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   fnorm = 0.0;
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         tmp = a[j];
         fnorm += tmp*tmp;
      }
   }
   fnorm = sqrt(fnorm);
   
   return(fnorm);
}

/*!
 * \fn JXF_Real fsls_CSRMatrixInfiniteNorm
 * \brief Compute the Infinite norm of a given matrix
 *        ||A||_{infinite} = max_{i}{\sum_{j}|a_{ij}|}
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/09/11  
 */
JXF_Real 
fsls_CSRMatrixInfiniteNorm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real sum  = 0.0;
   JXF_Real norm = 0.0;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   for (i = 0; i < n; i ++)
   {
      sum = 0.0;
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         sum += fabs(a[j]);
      }
      if (sum > norm) norm = sum;
   }

   return (norm);
}

/*!
 * \fn fsls_CSRMatrix * fsls_GetSubCSRMatrix
 * \brief Get a submatrix B from A.
 * \param *A pointer to the matrix
 * \param *ROW pointer to selected rows
 * \param *COL pointer to selected cols
 * \param mr number of selected rows
 * \param mc number of selected cols
 * \author peghoty 
 * \date 2010/01/19
 * \note Optimized a bit by peghoty on 2010/10/27
 */
fsls_CSRMatrix *
fsls_GetSubCSRMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc )
{
   JXF_Int    i,j;
   JXF_Int    row,colL;
   JXF_Int    nzB;

   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   fsls_CSRMatrix *B = NULL;
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   
   JXF_Int    *G2LCol = fsls_CTAlloc(JXF_Int, n);
   
   /* Initialize colglo2sub */
   for (i = 0; i < n; i ++) G2LCol[i] = -1;
   
   /* Build the map from glo to loc */
   for (i = 0; i < mc; i ++)
   {
      G2LCol[COL[i]] = i;
   }
   
   /* Compute the number of nonzeros in B */
   nzB = 0;
   for (i = 0; i < mr; i ++)
   {
      row = ROW[i];
      for (j = ia[row]; j < ia[row+1]; j ++)
      {
         if (G2LCol[ja[j]] >= 0) nzB ++;
      }
   }
   
   /* Create a submatrix */
   B = fsls_CSRMatrixCreate(mr,mc,nzB);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);
   b  = fsls_CSRMatrixData(B);
   
   /* Generate the submatrix B */
   nzB = 0;
   for (i = 0; i < mr; i ++)
   {
      ib[i] = nzB;
      row = ROW[i];
      for (j = ia[row]; j < ia[row+1]; j ++)
      {
          colL = G2LCol[ja[j]];
          if (colL > -1) 
          {
             b[nzB]  = a[j];
             jb[nzB] = colL;
             nzB ++;
          }
      }
   }
   ib[mr] = nzB;
   
   /* free some stuff */
   fsls_TFree(G2LCol); 

   return(B);
}


/*!
 * \fn JXF_Int fsls_GetSubDenseMatrix
 * \brief Abstract a dense submatrix B from a CSR matrix A, B is stored
 *        in a one-dimension array 'dense' row by row with the column 
 *        numbers ascendingly ordered.
 * \param *A pointer to the matrix
 * \param *ROW pointer to selected rows
 * \param *COL pointer to selected cols
 * \param mr number of selected rows
 * \param mc number of selected cols
 * \param *dense pointer to the dense submatrix, 'dense' should be previously zero-initialized. 
 * \param *G2LCol a working array of the size 'fsls_CSRMatrixNumRows(A)',
 *        which should be previously allocated and will serve as the column number
 *        mapping from global to local.
 * \author peghoty 
 * \date 2010/10/27 
 */
JXF_Int
fsls_GetSubDenseMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc, JXF_Real *dense, JXF_Int *G2LCol )
{
   /* members of A */
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   JXF_Int i,j;
   JXF_Int row,colL,colG;
    
  /*-------------------------------------------
   * Initialize 'G2LCol' 
   *------------------------------------------*/
   for (i = 0; i < n; i ++)
   { 
      G2LCol[i] = -1;
   }
   
  /*---------------------------------------------------
   * Build the column mapping from global to local 
   *-------------------------------------------------*/   
   for (i = 0; i < mc; i ++)
   {
      G2LCol[COL[i]] = i;
   }
   
  /*---------------------------------------------------
   * Generate the dense submatrix 'B' 
   *-------------------------------------------------*/   
   for (i = 0; i < mr; i ++)
   {
      row = ROW[i];            // global row number
      for (j = ia[row]; j < ia[row+1]; j ++)
      {
          colG = ja[j];        // global column number
          colL = G2LCol[colG]; // local  column number
          if (colL > -1) 
          {
             dense[i*mc+colL] = a[j];
          }
      }
   }

   return (0);
}

/*!
 * \fn fsls_DenseMatrix *fsls_SubDenseMatrix
 * \brief Abstract a dense submatrix B from a CSR matrix A, B is stored
 *        in a one-dimension array 'dense' row by row with the column 
 *        numbers ascendingly ordered.
 * \param *A pointer to the matrix
 * \param *ROW pointer to selected rows
 * \param *COL pointer to selected cols
 * \param mr number of selected rows
 * \param mc number of selected cols
 * \author peghoty 
 * \date 2010/11/19 
 */
fsls_DenseMatrix *
fsls_SubDenseMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc )
{
   /* members of A */
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   fsls_DenseMatrix *B = NULL;
   
   JXF_Int    *G2LCol = fsls_CTAlloc(JXF_Int, n); // serve as the column number mapping from global to local.
   JXF_Real *dense  = NULL;
   
   JXF_Int i,j;
   JXF_Int row,colL,colG;
   
   B = fsls_DenseMatrixCreate(mr,mc);
   fsls_DenseMatrixInitialize(B);
   dense = fsls_DenseMatrixData(B);
    
  /*-------------------------------------------
   * Initialize 'G2LCol' 
   *------------------------------------------*/
   for (i = 0; i < n; i ++)
   { 
      G2LCol[i] = -1;
   }
   
  /*---------------------------------------------------
   * Build the column mapping from global to local 
   *-------------------------------------------------*/   
   for (i = 0; i < mc; i ++)
   {
      G2LCol[COL[i]] = i;
   }
   
  /*---------------------------------------------------
   * Generate the dense submatrix 'B' 
   *-------------------------------------------------*/   
   for (i = 0; i < mr; i ++)
   {
      row = ROW[i];            // global row number
      for (j = ia[row]; j < ia[row+1]; j ++)
      {
          colG = ja[j];        // global column number
          colL = G2LCol[colG]; // local  column number
          if (colL > -1) 
          {
             dense[i*mc+colL] = a[j];
          }
      }
   }
   
   fsls_TFree(G2LCol);

   return (B);
}

/*!
 * \fn JXF_Int fsls_BuildCSRMatFromArrays
 * \brief Build a CSR Matrix from 'IA,JA,A' arrays
 * \param *ia IA of CSR format of matrix 'A'
 * \param *ja JA of CSR format of matrix 'A'
 * \param *a  DATA of CSR format of matrix 'A'
 * \param n size of matrix 'A'
 * \param nz number of nonzeros in 'A'
 * \param **A_ptr pointer to the resulting matrix
 * \author peghoty
 * \date 2010/01/20
 */
JXF_Int
fsls_BuildCSRMatFromArrays( JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n, JXF_Int nz, fsls_CSRMatrix **A_ptr )
{
   JXF_Int i;
   fsls_CSRMatrix *A = NULL;

   /* Generate the matrix */
   A = fsls_CSRMatrixCreate(n, n, nz);
   fsls_CSRMatrixInitialize(A);
   for (i = 0; i < n+1; i ++)
     fsls_CSRMatrixI(A)[i] = ia[i];
   for (i = 0; i < nz; i ++)
     fsls_CSRMatrixJ(A)[i] = ja[i];
   for (i = 0; i < nz; i ++)
     fsls_CSRMatrixData(A)[i] = a[i];
   fsls_CSRMatrixReorder(A); 

   *A_ptr = A;

   return (0);
}

/*!
 * \fn fsls_CSRMatrixReorderByPoint
 * \brief Reorder a matrix which is ordered by directions by points.
 * \param *A pointer to the matrix to be ordered.
 * \return a reordered matrix.
 * \note A is ordered in such way:
 *          (u1,u2,...,un; v1,v2,...,vn; w1,w2,...,wn)   [by direction]
 *       after reordering, B is ordered as fallows:
 *          (u1,v1,w1; u2,v2,w2; ..., un,vn,wn)          [by points]  
 * \author peghoty
 * \date 2011/08/21 
 */
fsls_CSRMatrix *
fsls_CSRMatrixReorderByPoint( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix  *B = NULL;

   /* CSR information of A */
   JXF_Int     N   = fsls_CSRMatrixNumRows(A);
   JXF_Int     nza = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   /* CSR information of B */
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   JXF_Int     nzb;
   
   /* local variables */
   JXF_Int    j,k,col;
   JXF_Int    pt;
   JXF_Int    RowA = 0;
   JXF_Int    RowB = 0;
   JXF_Int    n = N / 3;
   JXF_Int    nplusn = 2*n;
   JXF_Real val;
   

  /*-----------------------------------------------------
   * Create a new CSR matrix
   *----------------------------------------------------*/
      
   B = fsls_CSRMatrixCreate(N,N,nza);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);
   
  /*-----------------------------------------------------
   * Reorder A by n different points
   *----------------------------------------------------*/
   
   /* Initialize */
   ib[0] = 0;
   nzb = 0;
   
   /* Reordering */
   for (pt = 0; pt < n; pt ++)
   {
      for (k = 0; k < 3; k ++)
      {
         RowB = 3*pt + k;
         RowA = pt + k*n;
         
         /* dealing with the current row in matrix A */
         for (j = ia[RowA]; j < ia[RowA+1]; j ++)
         {
            val =  a[j];
            col = ja[j];
            if (col < n)
            {
               col = 3*col;
            }
            else if (col >= n && col < nplusn)
            {
               col = 3*(col - n) + 1;
            }
            else if (col >= nplusn)
            {
               col = 3*(col - nplusn) + 2;
            }         
            jb[nzb] = col;
             b[nzb] = val;
            nzb ++;
         }
         ib[RowB+1] = nzb;
      }
   }
  
#if 0 /* Unnecessary if A has been reordered! */
   fsls_CSRMatrixReorder(B);
#endif
   
   return (B);
}

/*!
 * \fn fsls_CSRMatrixReorderByDirection
 * \brief Reorder a matrix which is ordered by points by directions.
 * \param *A pointer to the matrix to be ordered.
 * \return a reordered matrix.
 * \note A is ordered in such way:
 *          (u1,v1,w1; u2,v2,w2; ..., un,vn,wn)          [by points]
 *       after reordering, B is ordered as fallows:
 *          (u1,u2,...,un; v1,v2,...,vn; w1,w2,...,wn)   [by direction]
 * \author peghoty
 * \date 2010/04/06 
 */
fsls_CSRMatrix *
fsls_CSRMatrixReorderByDirection( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix  *B = NULL;

   /* CSR information of A */
   JXF_Int     N   = fsls_CSRMatrixNumRows(A);
   JXF_Int     nza = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   /* CSR information of B */
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   JXF_Int     nzb;
   
   /* local variables */
   JXF_Int    j,col;
   JXF_Int    RowA = 0;
   JXF_Int    RowB = 0;
   JXF_Int    n = N / 3;
   JXF_Int    nplusn = 2*n;
   JXF_Real val;
   

  /*-----------------------------------------------------
   * Create a new CSR matrix
   *----------------------------------------------------*/
      
   B = fsls_CSRMatrixCreate(N,N,nza);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);
   
  /*-----------------------------------------------------
   * Reorder A by three different directions
   *----------------------------------------------------*/
   
   /* Initialize */
   ib[0] = 0;
   nzb = 0;
   
   /* Reordering */
   for (RowB = 0; RowB < N; RowB ++)
   {
      /* find the row number in matrix A */
      if (RowB < n)
      {
         RowA = 3*RowB;
      }
      else if (RowB >= n && RowB < nplusn)
      {
         RowA = (RowB-n)*3 + 1;
      }
      else if (RowB >= nplusn)
      {
         RowA = (RowB-nplusn)*3 + 2;
      }
      
      /* dealing with the current row in matrix A */
      for (j = ia[RowA]; j < ia[RowA+1]; j ++)
      {
         val =  a[j];
         col = ja[j];
         if (col % 3 == 0)
         {
            col = col / 3;
         }
         else if (col % 3 == 1)
         {
            col = n + (col-1) / 3;
         }
         else if (col % 3 == 2)
         {
            col = nplusn + (col-2) / 3;
         }         
         jb[nzb] = col;
          b[nzb] = val;
         nzb ++;
      }
      ib[RowB+1] = nzb;
   }

#if 0 /* Unnecessary if A has been reordered! */
   fsls_CSRMatrixReorder(B);
#endif
   
   return (B);
}

/*!
 * \fn fsls_CSRMatrixReorderByDirection2
 * \brief Reorder a matrix which is ordered by points by directions.
 * \param *A pointer to the matrix to be ordered.
 * \return a reordered matrix.
 * \note A is ordered in such way:
 *          (u1,v1,w1; u2,v2,w2; ..., un,vn,wn)          [by points]
 *       after reordering, B is ordered as fallows:
 *          (u1,u2,...,un; v1,v2,...,vn; w1,w2,...,wn)   [by direction]
 * \note This function is another implementation of fsls_CSRMatrixReorderByDirection,
 *       but much more longthy.
 * \author peghoty
 * \date 2010/04/06
 */
fsls_CSRMatrix *
fsls_CSRMatrixReorderByDirection2( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix  *B = NULL;

   /* CSR information of A */
   JXF_Int     N   = fsls_CSRMatrixNumRows(A);
   JXF_Int     nza = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   /* CSR information of B */
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   
   /* local variables */
   JXF_Int i,j,col;
   JXF_Int start,cnt;
   JXF_Int rowA;
   JXF_Int n = N / 3;
   JXF_Int nplusone = n + 1;
   JXF_Int nplusn   = 2*n;
   JXF_Int nplusnplusone = 2*n + 1;  
   JXF_Real val;
   

  /*-----------------------------------------------------
   * Create a new CSR matrix
   *----------------------------------------------------*/
      
   B = fsls_CSRMatrixCreate(N,N,nza);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);
   
   
  /*-----------------------------------------------------
   * Cenerate 'ib' array
   *----------------------------------------------------*/
   
   ib[0] = ia[0];
   
   for (i = 0; i < n; i ++)
   {
      j = 3*i;     // row number in A
      ib[i+1] = ib[i] + (ia[j+1] - ia[j]);
   }
   for (i = 0; i < n; i ++)
   {
      j = 3*i + 1; // row number in A
      ib[nplusone+i] = ib[n+i] + (ia[j+1] - ia[j]);
   }   
   for (i = 0; i < n; i ++)
   {
      j = 3*i + 2; // row number in A
      ib[nplusnplusone+i] = ib[nplusn+i] + (ia[j+1] - ia[j]);
   }   
   
   
  /*-----------------------------------------------------
   * Reorder A by three different directions
   *----------------------------------------------------*/

   /* Loop for the points */
   for (i = 0; i < n; i ++)
   {
       /* fist direction */
       start = ib[i];
       cnt = 0;
       rowA = 3*i;
       for (j = ia[rowA]; j < ia[rowA+1]; j ++)
       {
          col = ja[j];
          val =  a[j];
          
          if (col % 3 == 0)
          {
             col = col / 3;
          }
          else if (col % 3 == 1)
          {
             col = n + (col-1) / 3;
          }
          else if (col % 3 == 2)
          {
             col = nplusn + (col-2) / 3;
          }
          jb[start+cnt] = col;
           b[start+cnt] = val;
          cnt ++;
       }

       /* second direction */
       start = ib[n+i];
       cnt = 0;
       rowA = 3*i + 1;
       for (j = ia[rowA]; j < ia[rowA+1]; j ++)
       {
          col = ja[j];
          val =  a[j];
          
          if (col % 3 == 0)
          {
             col = col / 3;
          }
          else if (col % 3 == 1)
          {
             col = n + (col-1) / 3;
          }
          else if (col % 3 == 2)
          {
             col = nplusn + (col-2) / 3;
          }
          jb[start+cnt] = col;
           b[start+cnt] = val;
          cnt ++;
       }

       /* third direction */
       start = ib[2*n+i];
       cnt = 0;
       rowA = 3*i + 2;
       for (j = ia[rowA]; j < ia[rowA+1]; j ++)
       {
          col = ja[j];
          val =  a[j];
          
          if (col % 3 == 0)
          {
             col = col / 3;
          }
          else if (col % 3 == 1)
          {
             col = n + (col-1) / 3;
          }
          else if (col % 3 == 2)
          {
             col = nplusn + (col-2) / 3;
          }
          jb[start+cnt] = col;
           b[start+cnt] = val;
          cnt ++;
       }

   }
   
#if 0 /* Unnecessary if A has been reordered! */
   fsls_CSRMatrixReorder(B);
#endif

   return (B);
}  
 
/*!
 * \fn fsls_CSRMatrixReorderByDirectionRect
 * \brief Reorder a rectangle matrix which is ordered 
 *        by points by directions.
 * \param *A pointer to the matrix to be ordered.
 * \return a reordered matrix.
 * \note Suppose the rows of A are ordered in such way:
 *          (u1,v1,w1; u2,v2,w2; ..., un,vn,wn)          [by points]
 *       after reordering, B is ordered as fallows:
 *          (u1,u2,...,un; v1,v2,...,vn; w1,w2,...,wn)   [by direction]
 *       and the cols of A should be reordered in the same fasion.
 * \author peghoty
 * \date 2010/04/07
 */
fsls_CSRMatrix *
fsls_CSRMatrixReorderByDirectionRect( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix  *B = NULL;

   /* CSR information of A */
   JXF_Int     NR  = fsls_CSRMatrixNumRows(A);
   JXF_Int     NC  = fsls_CSRMatrixNumCols(A);
   JXF_Int     nza = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   /* CSR information of B */
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   JXF_Int     nzb;
   
   /* local variables */
   JXF_Int    j,col;
   JXF_Int    RowA = 0;
   JXF_Int    RowB = 0;
   JXF_Int    n = NR / 3;  // number of nodes in row
   JXF_Int    m = NC / 3;  // number of nodes in col 
     
   JXF_Int    nplusn = 2*n;
   JXF_Int    mplusm = 2*m;
   JXF_Real val;
   

  /*-----------------------------------------------------
   * Create a new CSR matrix
   *----------------------------------------------------*/
      
   B = fsls_CSRMatrixCreate(NR,NC,nza);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);
   
  /*-----------------------------------------------------
   * Reorder A by three different directions
   *----------------------------------------------------*/
   
   /* Initialize */
   ib[0] = 0;
   nzb = 0;
   
   /* Reordering */
   for (RowB = 0; RowB < NR; RowB ++)
   {
      /* find the row number in matrix A */
      if (RowB < n)
      {
         RowA = 3*RowB;
      }
      else if (RowB >= n && RowB < nplusn)
      {
         RowA = (RowB-n)*3 + 1;
      }
      else if (RowB >= nplusn)
      {
         RowA = (RowB-nplusn)*3 + 2;
      }
      
      /* dealing with the current row in matrix A */
      for (j = ia[RowA]; j < ia[RowA+1]; j ++)
      {
         val =  a[j];
         col = ja[j];
         if (col % 3 == 0)
         {
            col = col / 3;
         }
         else if (col % 3 == 1)
         {
            col = m + (col-1) / 3;
         }
         else if (col % 3 == 2)
         {
            col = mplusm + (col-2) / 3;
         }         
         jb[nzb] = col;
          b[nzb] = val;
         nzb ++;
      }
      ib[RowB+1] = nzb;
   }

#if 0 /* Unnecessary if A has been reordered! */
   fsls_CSRMatrixReorder(B);
#endif
   
   return (B);
}

/*!
 * \fn fsls_CSRMatrixGet3Block
 * \brief Abstract the three diagonal blocks of A to 
 *        form a new matrix.
 * \param fsls_CSRMatrix *A pointer to the matrix to be abstracted.
 * \return the new matrix B
 * \note If A is partitioned as follows:
 *
 *            / A11  A12  A13 \
 *       A =  | A21  A22  A23 |
 *            \ A31  A32  A33 /
 *
 *       then B equals to
 * 
 *            / A11   0    0  \
 *       B =  |  0   A22   0  |
 *            \  0    0   A33 /
 * \author peghoty
 * \date 2010/04/06 
 *
 */
fsls_CSRMatrix *
fsls_CSRMatrixGet3Block( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix  *B = NULL;

   /* CSR information of A */
   JXF_Int     N   = fsls_CSRMatrixNumRows(A);
   JXF_Int     nza = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   /* CSR information of B */
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;
   
   /* local variables */
   JXF_Int row,col;
   JXF_Int j,n = N / 3;
   JXF_Int cnt;
   
  /*-----------------------------------------------------
   * Create a new CSR matrix
   *----------------------------------------------------*/
      
   B = fsls_CSRMatrixCreate(N,N,nza);
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);
   
  /*-----------------------------------------------------
   * Reorder A by three different directions
   *----------------------------------------------------*/

   ib[0] = 0;
   cnt = 0;
   for (row = 0; row < N; row ++)
   {
      if (row < n)
      {
          for (j = ia[row]; j < ia[row+1]; j ++)
          {
             col = ja[j];
             if (col < n)
             {
                jb[cnt] = col;
                 b[cnt] = a[j];
                cnt ++;
             }
          }
      }
      else if (row >= n && row < 2*n)
      {
          for (j = ia[row]; j < ia[row+1]; j ++)
          {
             col = ja[j];
             if (col >= n && col < 2*n)
             {
                jb[cnt] = col;
                 b[cnt] = a[j];
                cnt ++;
             }
          }
      }
      else
      {
          for (j = ia[row]; j < ia[row+1]; j ++)
          {
             col = ja[j];
             if (col >= 2*n)
             {
                jb[cnt] = col;
                 b[cnt] = a[j];
                cnt ++;
             }
          }
      }
      ib[row+1] = cnt;
   }
   
   /* free redundant memory in advance and fill 
      the real number of nonzeros for B */
   jb = (JXF_Int *)realloc(jb, cnt*sizeof(JXF_Int));
   b  = (JXF_Real *)realloc(b, cnt*sizeof(JXF_Real));   
   fsls_CSRMatrixNumNonzeros(B) = cnt;

   fsls_CSRMatrixReorder(B);

   return (B);
}


/*!
 * \fn fsls_CSRMatrixReorderColumnNumber
 * \brief Reorder the column number of each row ascendingly while keeping
 *        all the diagonal entries firstly-stored. 
 * \note Both the column number and data part are reordered. 
 * \author peghoty
 * \date 2010/10/20
 */
void
fsls_CSRMatrixReorderColumnNumber( fsls_CSRMatrix *A )
{
   JXF_Int i,begin,end;
   JXF_Int n = fsls_CSRMatrixNumRows(A);
   JXF_Int *ia = fsls_CSRMatrixI(A);
   JXF_Int *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a = fsls_CSRMatrixData(A);

   for (i = 0; i < n; i ++)
   {
      begin = ia[i] + 1;
      end   = ia[i+1] - 1;
      fsls_diQuickSort12(a, ja, begin, end);       
   }
}

/*!
 * \fn fsls_CSRMatrixReorderColumnNumber2
 * \brief Reorder the column number of each row ascendingly 
 *        and the diagonal entries are included. 
 * \note Both the column number and data part are reordered.  
 * \author peghoty
 * \date 2010/10/29
 */
void
fsls_CSRMatrixReorderColumnNumber2( fsls_CSRMatrix *A )
{
   JXF_Int i,begin,end;
   JXF_Int n = fsls_CSRMatrixNumRows(A);
   JXF_Int *ia = fsls_CSRMatrixI(A);
   JXF_Int *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a = fsls_CSRMatrixData(A);

   for (i = 0; i < n; i ++)
   {
      begin = ia[i];
      end   = ia[i+1] - 1;
      fsls_diQuickSort12(a, ja, begin, end);      
   }
}

/*!
 * \fn fsls_CSRMatrixReorderColumnNumber3
 * \brief Reorder the column number of each row ascendingly 
 *        and the diagonal entries are included.
 * \note only the column number are reordered. 
 * \author peghoty
 * \date 2010/10/29
 */
void
fsls_CSRMatrixReorderColumnNumber3( fsls_CSRMatrix *A )
{
   JXF_Int i,begin,end;
   JXF_Int n = fsls_CSRMatrixNumRows(A);
   JXF_Int *ia = fsls_CSRMatrixI(A);
   JXF_Int *ja = fsls_CSRMatrixJ(A);

   for (i = 0; i < n; i ++)
   {
      begin = ia[i];
      end   = ia[i+1] - 1;
      fsls_iQuickSort12(ja, begin, end);     
   }
}

/*!
 * \fn JXF_Real fsls_CSRMatrixFindMinElm
 * \brief Find the minimal entry of the CSRMatrix
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/10/20  
 */
JXF_Real 
fsls_CSRMatrixFindMinElm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real minelm;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   minelm = a[0];
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         if (a[j] < minelm) minelm = a[j];
      }
   }
   
   return(minelm);
}

/*!
 * \fn JXF_Real fsls_CSRMatrixFindMinFabsElm
 * \brief Find the minimal absolute entry of the CSRMatrix
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/10/20  
 */
JXF_Real 
fsls_CSRMatrixFindMinFabsElm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real minelm;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   minelm = fabs(a[0]);
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         if (fabs(a[j]) < minelm) minelm = fabs(a[j]);
      }
   }
   
   return(minelm);
}

/*!
 * \fn JXF_Real fsls_CSRMatrixFindMaxElm
 * \brief Find the maximal entry of the CSRMatrix
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/10/20  
 */
JXF_Real 
fsls_CSRMatrixFindMaxElm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real maxelm = 0.0;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   maxelm = a[0];
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         if (a[j] > maxelm) maxelm = a[j];
      }
   }
   
   return(maxelm);
}

/*!
 * \fn JXF_Real fsls_CSRMatrixFindFabsMaxElm
 * \brief Find the maximal absolute entry of the CSRMatrix
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/10/20  
 */
JXF_Real 
fsls_CSRMatrixFindMaxFabsElm( fsls_CSRMatrix *A )
{
   JXF_Int    i,j;
   JXF_Real maxelm = 0.0;
   
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   maxelm = fabs(a[0]);
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         if (fabs(a[j]) > maxelm) maxelm = fabs(a[j]);
      }
   }
   
   return(maxelm);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixSparsityPatternCheck
 * \brief Compare two CSR matrices to make sure whether they have the same
 *        pattern of nonzeros.
 * \param col_reordered integer to indicate whether the column numbers in 
 *        each row have been ordered ascendingly.
 * \param *A the pointer to the first matrix
 * \param *B the pointer to the second matrix
 * \author peghoty
 * \date 2009/12/05
 */
JXF_Int
fsls_CSRMatrixSparsityPatternCheck( JXF_Int col_reordered, fsls_CSRMatrix *A, fsls_CSRMatrix *B )
{
   JXF_Int  *A_i      = fsls_CSRMatrixI(A);
   JXF_Int  *A_j      = fsls_CSRMatrixJ(A);
   JXF_Int   nrows_A  = fsls_CSRMatrixNumRows(A);
   JXF_Int   ncols_A  = fsls_CSRMatrixNumCols(A);
   JXF_Int   nz_A     = fsls_CSRMatrixNumNonzeros(A);
   
   JXF_Int  *B_i      = fsls_CSRMatrixI(B);
   JXF_Int  *B_j      = fsls_CSRMatrixJ(B);
   JXF_Int   nrows_B  = fsls_CSRMatrixNumRows(B);
   JXF_Int   ncols_B  = fsls_CSRMatrixNumCols(B);
   JXF_Int   nz_B     = fsls_CSRMatrixNumNonzeros(B);
   
   JXF_Int   i,same_sparsity_pattern = 1;
 
   if (nrows_A == nrows_B && ncols_A == ncols_B && nz_A == nz_B)
   {
      if (col_reordered == 0)
      {
         fsls_CSRMatrixReorderColumnNumber(A);
         fsls_CSRMatrixReorderColumnNumber(B);
      }
      
      for (i = 0; i < nrows_A; i ++)
      {
         if (A_i[i] != B_i[i]) 
         {
            same_sparsity_pattern = 0;
            jxf_printf("\n >> 1: A and B are of different sparsity pattern!\n\n");
            return (same_sparsity_pattern);
         }
      }

      for (i = 0; i < nz_A; i ++)
      {
         if (A_j[i] != B_j[i]) 
         {
            same_sparsity_pattern = 0;
            jxf_printf("\n >> 2: A and B are of different sparsity pattern!\n\n");
            return (same_sparsity_pattern);
         }
      }
      jxf_printf("\n >> A and B are of the same sparsity pattern!\n\n");
   }
   else
   {
      if (nrows_A != nrows_B || ncols_A != ncols_B)
      {
         jxf_printf("\n >> A and B are of different size!\n\n");
      }
      if (nz_A != nz_B)
      {
         jxf_printf("\n >> The numbers of nonzeros of A and B are different!\n\n");
      }      
      same_sparsity_pattern = 0;
   }
 
   return (same_sparsity_pattern);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixDiffFromFiles
 * \brief Find the difference between A and B
 * \note A and B must be of the same size.
 * \author peghoty
 * \date 2010/10/21
 */
JXF_Int 
fsls_CSRMatrixDiffFromFiles( char *MatFile1, char *MatFile2 )
{
   fsls_CSRMatrix  *A = NULL;
   fsls_CSRMatrix  *B = NULL;
   fsls_CSRMatrix  *C = NULL;

   JXF_Int nrows_A, ncols_A;
   JXF_Int nrows_B, ncols_B;
   
   fsls_BuildCSRMatFromFile(MatFile1, &A);
   fsls_BuildCSRMatFromFile(MatFile2, &B);
   
   nrows_A = fsls_CSRMatrixNumRows(A);
   ncols_A = fsls_CSRMatrixNumCols(A);
   nrows_B = fsls_CSRMatrixNumRows(B);
   ncols_B = fsls_CSRMatrixNumCols(B);
   
   if (nrows_A != nrows_B || ncols_A != ncols_B)
   {
      jxf_printf("\n >> A and B are of different size!\n\n");
      return (-1);
   }
   
   fsls_CSRMatrixScale(B, -1.0);
   C = fsls_CSRMatrixAdd(A, B);
   
   jxf_printf("\n >>> ||A-B||_F = %.16le", fsls_CSRMatrixFrobeniusNorm(C));
   jxf_printf("\n >>> min(A-B)  = %.16le", fsls_CSRMatrixFindMinElm(C));
   jxf_printf("\n >>> max(A-B)  = %.16le", fsls_CSRMatrixFindMaxElm(C));
   jxf_printf("\n >>> min|A-B|  = %.16le", fsls_CSRMatrixFindMinFabsElm(C));
   jxf_printf("\n >>> max|A-B|  = %.16le\n\n", fsls_CSRMatrixFindMaxFabsElm(C));

   fsls_CSRMatrixDestroy(A);
   fsls_CSRMatrixDestroy(B);
   fsls_CSRMatrixDestroy(C);
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixDiff
 * \brief Find the difference between A and B
 * \note A and B must be of the same size.
 * \author peghoty
 * \date 2010/12/12
 */
JXF_Int 
fsls_CSRMatrixDiff( fsls_CSRMatrix *A, fsls_CSRMatrix *B )
{
   fsls_CSRMatrix  *C = NULL;
   fsls_CSRMatrix  *B1 = NULL;

   JXF_Int nrows_A, ncols_A;
   JXF_Int nrows_B, ncols_B;
   
   nrows_A = fsls_CSRMatrixNumRows(A);
   ncols_A = fsls_CSRMatrixNumCols(A);
   nrows_B = fsls_CSRMatrixNumRows(B);
   ncols_B = fsls_CSRMatrixNumCols(B);
   
   if (nrows_A != nrows_B || ncols_A != ncols_B)
   {
      jxf_printf("\n >> A and B are of different size!\n\n");
      return (-1);
   }
   
   /* modified by peghoty 2010/12/28 */
   B1 = fsls_CSRMatrixCreate(nrows_B, ncols_B, fsls_CSRMatrixNumNonzeros(B));
   fsls_CSRMatrixInitialize(B1);
   fsls_CSRMatrixCopy(B, B1, 1);   
   
   fsls_CSRMatrixScale(B1, -1.0);
   C = fsls_CSRMatrixAdd(A, B1);
   
   jxf_printf("\n >>> ||A-B||_F = %.16le", fsls_CSRMatrixFrobeniusNorm(C));
   jxf_printf("\n >>> min(A-B)  = %.16le", fsls_CSRMatrixFindMinElm(C));
   jxf_printf("\n >>> max(A-B)  = %.16le", fsls_CSRMatrixFindMaxElm(C));
   jxf_printf("\n >>> min|A-B|  = %.16le", fsls_CSRMatrixFindMinFabsElm(C));
   jxf_printf("\n >>> max|A-B|  = %.16le\n\n", fsls_CSRMatrixFindMaxFabsElm(C));

   fsls_CSRMatrixDestroy(B1);
   fsls_CSRMatrixDestroy(C);
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixRelativeDiff
 * \brief Find the difference between A and B with the same nonzero pattern.
 * \note A and B must be of the same size and of the same nonzero pattern.
 * \author peghoty
 * \date 2011/08/27
 */
JXF_Int 
fsls_CSRMatrixRelativeDiff( fsls_CSRMatrix *A, fsls_CSRMatrix *B, char *filename )
{
   JXF_Int nrows_A, ncols_A;
   JXF_Int nrows_B, ncols_B;
   
   nrows_A = fsls_CSRMatrixNumRows(A);
   ncols_A = fsls_CSRMatrixNumCols(A);
   nrows_B = fsls_CSRMatrixNumRows(B);
   ncols_B = fsls_CSRMatrixNumCols(B);
   
   fsls_CSRMatrix  *C = NULL;
    
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *aa = fsls_CSRMatrixData(A); 
   JXF_Real *bb = fsls_CSRMatrixData(B);
   
   JXF_Real *cc = NULL; 
   
   JXF_Int i, j; 
   JXF_Real diff;       
   
   if (nrows_A != nrows_B || ncols_A != ncols_B)
   {
      jxf_printf("\n >> A and B are of different size!\n\n");
      return (-1);
   }
   
   /* modified by peghoty 2010/12/28 */
   C = fsls_CSRMatrixCreate(nrows_A, ncols_A, fsls_CSRMatrixNumNonzeros(A));
   fsls_CSRMatrixInitialize(C);
   cc = fsls_CSRMatrixData(C); 
   fsls_CSRMatrixCopy(A, C, 0);
   
   for (i = 0; i < nrows_A; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         diff = fabs(aa[j] - bb[j]);
         if (aa[j])
         {
            cc[j] = diff / fabs(aa[j]);
         }
         else
         {
            jxf_printf(" A[%d][%d] = 0.0", i, ja[j]);
            cc[j] = diff;
         }
      }
   }
   
   if (filename)
   {
      fsls_CSRMatrixPrint(C, filename);
   }

   jxf_printf("\n >>> ||A-B||_F = %.16le", fsls_CSRMatrixFrobeniusNorm(C));
   jxf_printf("\n >>> min(A-B)  = %.16le", fsls_CSRMatrixFindMinElm(C));
   jxf_printf("\n >>> max(A-B)  = %.16le", fsls_CSRMatrixFindMaxElm(C));
   jxf_printf("\n >>> min|A-B|  = %.16le", fsls_CSRMatrixFindMinFabsElm(C));
   jxf_printf("\n >>> max|A-B|  = %.16le\n\n", fsls_CSRMatrixFindMaxFabsElm(C));

   fsls_CSRMatrixDestroy(C);
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRGroupCreate
 * \brief Create a fsls_CSRGroup.  
 * \author peghoty
 * \date 2010/11/16
 */
fsls_CSRGroup *
fsls_CSRGroupCreate( JXF_Int n )
{
   fsls_CSRGroup *csrgroup = fsls_CTAlloc(fsls_CSRGroup, 1);
   
   fsls_CSRGroupN(csrgroup)      = n;
   fsls_CSRGroupAArray(csrgroup) = NULL;
   return csrgroup;
} 
 
/*!
 * \fn JXF_Int fsls_CSRGroupInitialize
 * \brief Initialize a fsls_CSRGroup.  
 * \author peghoty
 * \date 2010/11/16
 */
JXF_Int
fsls_CSRGroupInitialize( fsls_CSRGroup *csrgroup )
{
   JXF_Int  ierr = 0;
   JXF_Int  n = fsls_CSRGroupN(csrgroup);
   
   if ( n && !fsls_CSRGroupAArray(csrgroup) )
   { 
      fsls_CSRGroupAArray(csrgroup) = fsls_CTAlloc(fsls_CSRMatrix *, n);
   }
   
   return ierr;
} 

/*!
 * \fn JXF_Int fsls_CSRGroupDiagCombine
 * \brief Diagonally Combine a series of CSR matrices.  
 *         |A1                  |
 *         |   A2               |
 *  A  =   |     ...            |
 *         |         A_{n-1}    | 
 *         |             A_{n}  |
 * \author peghoty
 * \date 2010/11/16
 */
fsls_CSRMatrix *
fsls_CSRGroupDiagCombine( fsls_CSRGroup *csrgroup )
{
   /* members of csrgoup */
   JXF_Int               n       = fsls_CSRGroupN(csrgroup);
   fsls_CSRMatrix  **A_array = fsls_CSRGroupAArray(csrgroup);

   fsls_CSRMatrix   *A  = NULL;
   fsls_CSRMatrix   *Ag = NULL;
   JXF_Int    *iag = NULL;
   JXF_Int    *jag = NULL;
   JXF_Real *ag  = NULL;   
   
   JXF_Int num_rows = 0;
   JXF_Int num_cols = 0;
   JXF_Int num_nonzeros = 0;
   JXF_Int    *ia = NULL;
   JXF_Int    *ja = NULL;
   JXF_Real *a  = NULL;
   
   JXF_Int i,j,m;
   JXF_Int row,rowg,colg,cnt,dis;
   
   if (n == 1)
   {
      num_rows = fsls_CSRMatrixNumRows(A_array[0]);
      num_cols = fsls_CSRMatrixNumCols(A_array[0]);
      num_nonzeros = fsls_CSRMatrixNumNonzeros(A_array[0]);
      A = fsls_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
      fsls_CSRMatrixInitialize(A);
      fsls_CSRMatrixCopy(A_array[0], A, 1);
      return A;
   }
   
   //---------------------------------------------------------
   //  Caculate 'num_rows, num_cols, num_nonzeros' for A
   //---------------------------------------------------------
   for (i = 0; i < n; i ++)
   {
      num_rows += fsls_CSRMatrixNumRows(A_array[i]);
      num_cols += fsls_CSRMatrixNumCols(A_array[i]);
      num_nonzeros += fsls_CSRMatrixNumNonzeros(A_array[i]);
   }

   //---------------------------------------------------------
   //  Create and Initialize a CSR matrix A
   //---------------------------------------------------------
   A = fsls_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   fsls_CSRMatrixInitialize(A);
   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);
   a  = fsls_CSRMatrixData(A);
   
   //---------------------------------------------------------
   //  Generate the 'ia' for A
   //---------------------------------------------------------
   cnt = 0; dis = 0; row = 0;
   ia[0] = 0;
   for (m = 0; m < n; m ++)
   {
      Ag   = A_array[m];
      rowg = fsls_CSRMatrixNumRows(Ag);
      colg = fsls_CSRMatrixNumCols(Ag);
      iag  = fsls_CSRMatrixI(Ag);
      jag  = fsls_CSRMatrixJ(Ag);
      ag   = fsls_CSRMatrixData(Ag);
      for (i = 0; i < rowg; i ++)
      {
         row ++;
         ia[row] = ia[row-1] + iag[i+1] - iag[i];
         for (j = iag[i]; j < iag[i+1]; j ++)
         {
            ja[cnt] = jag[j] + dis;
             a[cnt] = ag[j];
            cnt ++;
         }
      }
      dis += colg;
   }
   
   return A;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixPartitionByRow0
 * \brief Split a CSR matrix into 'np' rectangular matrices by row-partition.
 * \param *A pointer to the CSR matrix to be splitted
 * \param np number of rectangular matrices after splitting
 * \param *partition the row-partition(size is 'np+1')
 * \param ***recA_array_ptr pointer to all the rectangular matrices
 * \note This function is not efficient enough and it has been replaced 
 *       by 'fsls_CSRMatrixPartitionByRow'.
 * \author peghoty
 * \date 2010/11/22
 */
JXF_Int 
fsls_CSRMatrixPartitionByRow0( fsls_CSRMatrix *A, JXF_Int np, JXF_Int *partition, fsls_CSRMatrix ***recA_array_ptr )
{
   fsls_CSRMatrix **recA_array = fsls_CTAlloc(fsls_CSRMatrix *, np);
   
   JXF_Int num_rows = fsls_CSRMatrixNumRows(A);
   JXF_Int num_cols = fsls_CSRMatrixNumCols(A);
   
   JXF_Int *ROW = fsls_CTAlloc(JXF_Int, num_rows);
   JXF_Int *COL = fsls_CTAlloc(JXF_Int, num_cols);
   
   JXF_Int i,begin;
   JXF_Int local_num_rows;
   
   for (i = 0; i < num_rows; i ++)
   {
      ROW[i] = i;
   }   
   for (i = 0; i < num_cols; i ++)
   {
      COL[i] = i;
   }
   
   for (i = 0; i < np; i ++)
   {
      local_num_rows = partition[i+1] - partition[i];
      begin = partition[i];
      recA_array[i] = fsls_GetSubCSRMatrix(A, &ROW[begin], COL, local_num_rows, num_cols);
   }
   
   fsls_TFree(ROW);
   fsls_TFree(COL); 
   
   *recA_array_ptr = recA_array; 
   
   return 0; 
}

/*!
 * \fn JXF_Int fsls_CSRMatrixPartitionByRow
 * \brief Split a CSR matrix into 'np' rectangular matrices by row-partition.
 * \param *A pointer to the CSR matrix to be splitted
 * \param np number of rectangular matrices after splitting
 * \param *partition the row-partition(size is 'np+1')
 * \param ***recA_array_ptr pointer to all the rectangular matrices
 * \author peghoty
 * \date 2010/11/23
 */
JXF_Int 
fsls_CSRMatrixPartitionByRow( fsls_CSRMatrix *A, JXF_Int np, JXF_Int *partition, fsls_CSRMatrix ***recA_array_ptr )
{
   JXF_Int num_cols = fsls_CSRMatrixNumCols(A);
   JXF_Int    *ia   = fsls_CSRMatrixI(A);
   JXF_Int    *ja   = fsls_CSRMatrixJ(A);
   JXF_Real *a    = fsls_CSRMatrixData(A);
   
   fsls_CSRMatrix **recA_array = fsls_CTAlloc(fsls_CSRMatrix *, np);
   
   JXF_Int    *ia_rec = NULL;
   JXF_Int    *ja_rec = NULL;
   JXF_Real *a_rec  = NULL;
     
   JXF_Int i,j,begin,end;
   JXF_Int num_nonzeros;
   JXF_Int local_num_rows, local_num_rows1;
   
   for (i = 0; i < np; i ++)
   {
      begin = partition[i];
      end   = partition[i+1];
      local_num_rows = end - begin;
      num_nonzeros = ia[end] - ia[begin];
      recA_array[i] = fsls_CSRMatrixCreate(local_num_rows, num_cols, num_nonzeros);
      fsls_CSRMatrixInitialize(recA_array[i]);
      ia_rec = fsls_CSRMatrixI(recA_array[i]);
      ja_rec = fsls_CSRMatrixJ(recA_array[i]);
      a_rec  = fsls_CSRMatrixData(recA_array[i]);
      memcpy(ia_rec, &ia[begin],     (local_num_rows+1)*sizeof(JXF_Int));
      memcpy(ja_rec, &ja[ia[begin]], num_nonzeros*sizeof(JXF_Int));
      memcpy(a_rec,  &a[ia[begin]],  num_nonzeros*sizeof(JXF_Real));
      local_num_rows1 = local_num_rows + 1;
      if (i > 0) // Generraly speaking, this is not necessary.
      {
         for (j = 0; j < local_num_rows1; j ++)
         {
            ia_rec[j] -= ia[begin];
         }
      }
   }
   
   *recA_array_ptr = recA_array; 
   
   return 0; 
}

/*!
 * \fn JXF_Int *fsls_CSRMatrixNumNZPerRow
 * \brief Compute the number of nonzeros in each row of A.
 * \param *A pointer to the CSR matrix
 * \author peghoty
 * \date 2010/12/12
 */
JXF_Int * 
fsls_CSRMatrixNumNZPerRow( fsls_CSRMatrix *A )
{
   /* information of the matrix */
   JXF_Int  n  = fsls_CSRMatrixNumRows(A);
   JXF_Int *ia = fsls_CSRMatrixI(A);
   
   JXF_Int *nzperrow = fsls_CTAlloc(JXF_Int, n);
   JXF_Int  i;

   for (i = 0; i < n; i ++)
   {
      nzperrow[i] = ia[i+1] - ia[i];
   }
 
   return (nzperrow);
}

/*!
 * \fn JXF_Int *fsls_CSRMatrixReOrderByNZPerRow
 * \brief Get the ordering in which the number of nonzeros 
 *        in each row of A is ascending.
 * \param *A pointer to the CSR matrix
 * \param order ascendingly or descendingly? 12/default: ascending; 21: descending 
 * \author peghoty
 * \date 2010/12/12
 */
JXF_Int * 
fsls_CSRMatrixReOrderByNZPerRow( fsls_CSRMatrix *A, JXF_Int order )
{
   JXF_Int  n = fsls_CSRMatrixNumRows(A);
   JXF_Int *p = fsls_CTAlloc(JXF_Int, n);
   JXF_Int *q = NULL;
   JXF_Int  i;
   
   for (i = 0; i < n; i ++) p[i] = i;
   
   q = fsls_CSRMatrixNumNZPerRow(A);
   
   switch (order)
   {
      case 12: // ascendingly
      fsls_iQuickSortIndex12(q, 0, n-1, p);
      break;
      
      case 21: // descendingly
      fsls_iQuickSortIndex21(q, 0, n-1, p);
      break;
      
      default: // ascendingly
      fsls_iQuickSortIndex12(q, 0, n-1, p);
      break;      
   }
   
   fsls_TFree(q);
   return (p);
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixPermuteTest
 * \brief Reorder a CSR matrix A by computing P*A*P^T, where P is defined by p, an ordering array.
 * \param *A pointer to the CSR matrix to be reordered
 * \param *p pointer to the reordering(same size as A)
 * \note This is a test version for checking the correctness of 'fsls_CSRMatrixPAPT'
 *       and it's of low efficiency.
 * \author peghoty
 * \date 2010/12/09
 */
fsls_CSRMatrix * 
fsls_CSRMatrixPermuteTest( fsls_CSRMatrix *A, JXF_Int *p )
{
   JXF_Int i;
   JXF_Int n = fsls_CSRMatrixNumRows(A);
   fsls_CSRMatrix *P = NULL;
   fsls_CSRMatrix *R = NULL;
   fsls_CSRMatrix *B = NULL;
   
   P = fsls_CSRMatrixCreate(n, n, n);
   fsls_CSRMatrixInitialize(P);
   
   for (i = 0; i < n+1; i ++)
   {
      fsls_CSRMatrixI(P)[i] = i;
   }
   memcpy(fsls_CSRMatrixJ(P), p, n*sizeof(JXF_Int));
   for (i = 0; i < n; i ++)
   {
      fsls_CSRMatrixData(P)[i] = 1.0;
   }
   
   fsls_CSRMatrixTranspose(P, &R, 1);

   B = fsls_CSRMatrixTriMultiply(P, A, R);
   
   fsls_CSRMatrixDestroy(R);
   fsls_CSRMatrixDestroy(P);

   return B;
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixPermute
 * \brief Reorder a CSR matrix A by computing P*A*P^T, where P is defined by p, an ordering array.
 * \param *A pointer to the CSR matrix to be reordered
 * \param *p pointer to the reordering(same size as A)
 * \author XuSenlin,modified by peghoty
 * \date 2010/12/12
 */
fsls_CSRMatrix * 
fsls_CSRMatrixPermute( fsls_CSRMatrix *A, JXF_Int *p )
{
   /* information of the matrix */
   JXF_Int     n  = fsls_CSRMatrixNumRows(A);
   JXF_Int     nz = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);

   /* the resulting CSR matrix B */
   fsls_CSRMatrix *B  = NULL;
   JXF_Int            *ib = NULL;
   JXF_Int            *jb = NULL;
   JXF_Real         *b  = NULL;
   
   /* local veriables */
   JXF_Int  i,j,k,m;
   JXF_Int *q = fsls_CTAlloc(JXF_Int, n);   
   
   //----------------------------------------//
   //       Create a CSR matrix B            //
   //----------------------------------------//
   
   B = fsls_CSRMatrixCreate(n, n, nz);
   fsls_CSRMatrixInitialize(B);
   b  = fsls_CSRMatrixData(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);  

   //-----------------------------------------------//
   // Generate the CSR components of B: (ib,jb,b)   //
   //-----------------------------------------------//
   	
  /*------------------------------------------------
   * get 'q' which is the inverse mapping of 'p' 
   *-----------------------------------------------*/
   for (i = 0; i < n; i ++)
   {
      q[p[i]] = i;
   }

  /*------------------------------------------------
   * get 'ib' through 'ia' and 'p' 
   *-----------------------------------------------*/
   ib[0] = 0;
   for (i = 0; i < n; i ++)
   {
      ib[i+1] = ib[i] + ia[p[i]+1] - ia[p[i]];
   }

  /*------------------------------------------------
   * get 'b' and 'jb' through 'ja', 'a' and 'p' 
   *-----------------------------------------------*/
   for (i = 0; i < n; i ++)
   {
      for (k = ia[i]; k < ia[i+1]; k ++)
      {
         j = ja[k];     // i,j: row and column number in A
         m = ib[q[i]];  // q[i]: row number in B
         jb[m] = q[j];  // q[j]: column number in B
          b[m] = a[k];
         ib[q[i]] ++;
      }
   }

  /*------------------------------------------------
   * recover the 'ib' 
   *-----------------------------------------------*/
   for (i = n-1; i > 0; i --)
   {
      ib[i] = ib[i-1];
   }
   ib[0] = 0;

  /*---------------------------
   * free 'q' and return
   *--------------------------*/
   fsls_TFree(q);
   return (B);
}

/*!
 * \fn void fsls_CSRMatrixIndexC2Fortran
 * \brief Index displacing from the C manner to Fortran manner for a CSR matrix.
 * \param n size of the matrix
 * \param nz number of nonzeros in the matrix
 * \param *ia IA of the CSR format
 * \param *ja JA of the CSR format
 * \author peghoty
 * \date 2010/12/28
 */
void
fsls_CSRMatrixIndexC2Fortran( JXF_Int n, JXF_Int nz, JXF_Int *ia, JXF_Int *ja )
{
   JXF_Int i;
   JXF_Int n1 = n + 1;
   
   for (i = 0; i < n1; i ++) 
   {
      ia[i] ++;
   }
   for (i = 0; i < nz; i ++) 
   {
      ja[i] ++;
   }
}

/*!
 * \fn void fsls_CSRMatrixIndexFortran2C
 * \brief Index displacing from the Fortran manner to C manner for a CSR matrix.
 * \param n size of the matrix
 * \param nz number of nonzeros in the matrix
 * \param *ia IA of the CSR format
 * \param *ja JA of the CSR format
 * \author peghoty
 * \date 2010/12/28
 */
void
fsls_CSRMatrixIndexFortran2C( JXF_Int n, JXF_Int nz, JXF_Int *ia, JXF_Int *ja )
{
   JXF_Int i;
   JXF_Int n1 = n + 1;
   
   for (i = 0; i < n1; i ++) 
   {
      ia[i] --;
   }
   for (i = 0; i < nz; i ++) 
   {
      ja[i] --;
   }
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixFilter
 * \brief Form a new CSR matrix after filtering the 'small' entries 
 *        which are defined by 'strong_threshold'.
 * \param *A pointer to the CSR matrix to be filtered
 * \param strong_threshold the parameter to define the strong-dependence.
 * \note The original matrix will be kept.
 * \author peghoty
 * \date 2011/01/10
 */
fsls_CSRMatrix * 
fsls_CSRMatrixFilter( fsls_CSRMatrix *A, JXF_Real strong_threshold )
{
   /* information of the matrix */
   JXF_Int num_rows     = fsls_CSRMatrixNumRows(A);
   JXF_Int num_cols     = fsls_CSRMatrixNumCols(A);
   JXF_Int num_nonzeros = fsls_CSRMatrixNumNonzeros(A);

   fsls_CSRMatrix *B = NULL;
   JXF_Int    *ib = NULL;
   JXF_Real *b  = NULL;
   
   JXF_Int    i,j;   
   JXF_Real offdmax; 
   JXF_Real tmp; 

   B = fsls_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
   fsls_CSRMatrixInitialize(B);
   fsls_CSRMatrixCopy(A, B, 1);
   fsls_CSRMatrixReorder(B); // reorder B
   ib = fsls_CSRMatrixI(B);
   b  = fsls_CSRMatrixData(B);
   
   for (i = 0; i < num_rows; i ++)
   {
      offdmax = 0.0;
      for (j = ib[i]+1; j < ib[i+1]; j ++)
      {
         tmp = fabs(b[j]);
         if (tmp > offdmax)
         {
            offdmax = tmp;
         }
      }
      
      for (j = ib[i]+1; j < ib[i+1]; j ++)
      {
         tmp = fabs(b[j]);
         if ( tmp <= strong_threshold*offdmax )
         {
            b[j] = 0.0;
         }
      } 
   }
   
   B = fsls_CSRMatrixDeleteZeros(B, 0.0);
   
   return (B);
}

/*!
 * \fn fsls_CSRMatrix *fsls_CSRMatrixIdentity
 * \brief Generate a Identity matrix of the given size.
 * \param n size of the matrix
 * \author peghoty
 * \date 2011/02/25
 */
fsls_CSRMatrix * 
fsls_CSRMatrixIdentity( JXF_Int n )
{
   fsls_CSRMatrix *matrix = NULL;
   
   JXF_Int *ia = NULL;
   JXF_Int *ja = NULL;
   JXF_Real *a = NULL;
   JXF_Int i;
   JXF_Int nplus1 = n + 1;

   matrix = fsls_CSRMatrixCreate(n, n, n);
   fsls_CSRMatrixInitialize(matrix);

   a  = fsls_CSRMatrixData(matrix);
   ia = fsls_CSRMatrixI(matrix);
   ja = fsls_CSRMatrixJ(matrix);
   
   for (i = 0; i < nplus1; i ++) ia[i] = i;
   for (i = 0; i < n; i ++)
   {
      ja[i] = i;
      a[i]  = 1.0;
   }

   return(matrix);   
}

/*!
 * \fn JXF_Int fsls_CSRMatrixLowerPart
 * \brief Obtain the lower part of a CSR matrix.
 * \param *A pointer to the CSR matrix
 * \param reorder flag to indicate whether the diagonal entries have been
 *        stored firstly, if reorder = 1, the function will reorder the
 *        matrix in the very beginning. if reorder = 0, the reordering is 
 *        unneccessary.
 * \param **L_ptr pointer to the pointer to the resulting Lower part matrix
 * \note we always assume all the diagonal entries are nonzeros. 
 * \note The diagonal entries of 'L' are stored firstly in each row.  
 * \author peghoty
 * \date 2011/04/13 
 */
JXF_Int 
fsls_CSRMatrixLowerPart( fsls_CSRMatrix *A, JXF_Int reorder, fsls_CSRMatrix **L_ptr )
{
   JXF_Int     n  = fsls_CSRMatrixNumRows(A); 
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A); 
   
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;    
   
   fsls_CSRMatrix *L = NULL; 

   JXF_Int i,j,k;
   JXF_Int begin,end;
   JXF_Int Lcnt;
   
   /* Ensure the diagonal entries firstly stored */
   if (reorder == 1)
   {
      fsls_CSRMatrixReorder(A);
   }
   
   /* Get the number of nonzeros of matrix 'L' */
   Lcnt = 0;
   for (i = 0; i < n; i ++)
   {
      begin = ia[i];
      end   = ia[i+1];
      for (k = begin; k < end; k ++)
      {
         if (ja[k] <= i) Lcnt ++;
      }
   }
   
   /* Create a CSR matrix 'L' */
   L = fsls_CSRMatrixCreate(n, n, Lcnt);
   fsls_CSRMatrixInitialize(L);
   ib = fsls_CSRMatrixI(L);
   jb = fsls_CSRMatrixJ(L);
   b  = fsls_CSRMatrixData(L); 
   
   /* Generate matrix 'L' */
   Lcnt = 0;
   ib[0] = Lcnt;
   for (i = 0; i < n; i ++)
   {
      b[Lcnt]  = a[ia[i]];
      jb[Lcnt] = i;
      Lcnt ++;
       
      begin = ia[i] + 1;
      end   = ia[i+1];
      for (k = begin; k < end; k ++)
      {
         j = ja[k];
         if (j <= i) 
         {
            b[Lcnt]  = a[k];
            jb[Lcnt] = j;
            Lcnt ++;
         }
      }
      ib[i+1] = Lcnt;
   }
   
   *L_ptr = L;
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixUpperPart
 * \brief Obtain the upper part of a CSR matrix.
 * \param *A pointer to the CSR matrix
 * \param reorder flag to indicate whether the diagonal entries have been
 *        stored firstly, if reorder = 1, the function will reorder the
 *        matrix in the very beginning. if reorder = 0, the reordering is 
 *        unneccessary.
 * \param **U_ptr pointer to the pointer to the resulting upper part matrix
 * \note we always assume all the diagonal entries are nonzeros. 
 * \note The diagonal entries of 'U' are stored firstly in each row. 
 * \author peghoty
 * \date 2011/04/13 
 */
JXF_Int 
fsls_CSRMatrixUpperPart( fsls_CSRMatrix *A, JXF_Int reorder, fsls_CSRMatrix **U_ptr )
{
   JXF_Int     n  = fsls_CSRMatrixNumRows(A); 
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A); 
   
   JXF_Int    *ib = NULL;
   JXF_Int    *jb = NULL;
   JXF_Real *b  = NULL;    
   
   fsls_CSRMatrix *U = NULL; 

   JXF_Int i,j,k;
   JXF_Int begin,end;
   JXF_Int Ucnt;
   
   /* Ensure the diagonal entries firstly stored */
   if (reorder == 1)
   {
      fsls_CSRMatrixReorder(A);
   }
   
   /* Get the number of nonzeros of matrix 'U' */
   Ucnt = 0;
   for (i = 0; i < n; i ++)
   {
      begin = ia[i];
      end   = ia[i+1];
      for (k = begin; k < end; k ++)
      {
         if (ja[k] >= i) Ucnt ++;
      }
   }
   
   /* Create a CSR matrix 'U' */
   U = fsls_CSRMatrixCreate(n, n, Ucnt);
   fsls_CSRMatrixInitialize(U);
   ib = fsls_CSRMatrixI(U);
   jb = fsls_CSRMatrixJ(U);
   b  = fsls_CSRMatrixData(U); 
   
   /* Generate matrix 'U' */
   Ucnt = 0;
   ib[0] = Ucnt;
   for (i = 0; i < n; i ++)
   {
      b[Ucnt]  = a[ia[i]];
      jb[Ucnt] = i;
      Ucnt ++;
       
      begin = ia[i] + 1;
      end   = ia[i+1];
      for (k = begin; k < end; k ++)
      {
         j = ja[k];
         if (j >= i) 
         {
            b[Ucnt]  = a[k];
            jb[Ucnt] = j;
            Ucnt ++;
         }
      }
      ib[i+1] = Ucnt;
   }
   
   *U_ptr = U;
   
   return 0;
}

/*!
 * \fn JXF_Real *fsls_CSRMatrixDiagonalEntries
 * \brief Obtain the diagonal part of a CSR matrix.
 * \param *A pointer to the CSR matrix
 * \author peghoty
 * \date 2011/04/13 
 */
JXF_Real *
fsls_CSRMatrixDiagonalEntries( fsls_CSRMatrix *A )
{
   JXF_Int     n  = fsls_CSRMatrixNumRows(A); 
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Int    *ja = fsls_CSRMatrixJ(A);
   JXF_Real *a  = fsls_CSRMatrixData(A); 
   
   JXF_Int i, k, begin, end;
   JXF_Real *diag = fsls_CTAlloc(JXF_Real, n);
   
   for (i = 0; i < n; i ++)
   {
      begin = ia[i];
      end   = ia[i+1];
      for (k = begin; k < end; k ++)
      {
         if (ja[k] == i) 
         {
            diag[i] = a[k];
            break;
         }
      }
   }
   
   return (diag);
}

/*!
 * \fn fsls_CSRMatrix *fsls_MatfromREItoEIR
 * \brief Transfer a R-E-I ordered CSRMatrix into an E-I-R ordered CSRMatrix.
 *    /           \      /           \  
 *   | Arr Are  0  |    | Aee Aei Aer |  
 *   | Aer Aee Aei | => | Aie Aii  0  |
 *   |  0  Aie Aii |    | Are  0  Arr |
 *    \           /      \           / 
 * \author peghoty
 * \date 2011/06/26 
 */
fsls_CSRMatrix *
fsls_MatfromREItoEIR( fsls_CSRMatrix *A )
{
   fsls_CSRMatrix *A11 = NULL;
   fsls_CSRMatrix *A22 = NULL;
   fsls_CSRMatrix *A33 = NULL;
  
   fsls_Vector *V12 = NULL; 
   fsls_Vector *V21 = NULL;
   fsls_Vector *V23 = NULL;
   fsls_Vector *V32 = NULL;
   
   /* CSR information of A */
   JXF_Int     N   = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);

   /* CSR information of B */
   fsls_CSRMatrix *B = NULL;
   JXF_Int    *ib  = NULL;
   JXF_Int    *jb  = NULL;     
   JXF_Real *b   = NULL;
   
   /* CSR information of submatrices */
   JXF_Int    *ia11 = NULL;
   JXF_Int    *ia22 = NULL;
   JXF_Int    *ia33 = NULL;
   
   JXF_Int    *ja11 = NULL;
   JXF_Int    *ja22 = NULL;
   JXF_Int    *ja33 = NULL;
   
   JXF_Real *a11 = NULL;
   JXF_Real *a22 = NULL;
   JXF_Real *a33 = NULL;
   
   /* Vector information */
   JXF_Real *a12 = NULL;
   JXF_Real *a21 = NULL;
   JXF_Real *a23 = NULL;
   JXF_Real *a32 = NULL;
   
   /* local variables */
   JXF_Int row,col;
   JXF_Int j,n = N / 3;
   JXF_Int nza11,nza22,nza33;
   JXF_Int nzv12,nzv23,nzv21,nzv32;
   JXF_Int nzb;
   JXF_Int nplusn = 2*n;
   
  /*--------------------------------------------------------------------
   * Counter the numbers of nonzero entries in each diagonal block
   *------------------------------------------------------------------*/
   nza11 = 0;
   nza22 = 0;
   nza33 = 0;

   //=================================================================
   // Remark: Since Arr, Aee and Aii have the same nonzero-structure,
   //         only nza11 should be counted, nza22 = nza33 = nza11.
   //                                        peghoty 2011/06/26
   //=================================================================
      
   for (row = 0; row < n; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col < n)
          {
             nza11 ++;
          }
       }     
   }
 
   for (row = n; row < nplusn; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col >= n && col < nplusn)
          {
             nza22 ++;
          }
       }     
   }
 
   for (row = nplusn; row < N; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col >= nplusn)
          {
             nza33 ++;
          }
       }
   }   

#if 1
   jxf_printf(" nza11 = %d\n", nza11);
   jxf_printf(" nza22 = %d\n", nza22);
   jxf_printf(" nza33 = %d\n", nza33);
   jxf_printf(" nztot = %d\n", nza11+nza22+nza33+4*n);
   jxf_printf(" nz(A) = %d\n", fsls_CSRMatrixNumNonzeros(A));
#endif

  /*-----------------------------------------------------
   * Create 3 CSR submatrices
   *----------------------------------------------------*/
   A11 = fsls_CSRMatrixCreate(n,n,nza11);
   fsls_CSRMatrixInitialize(A11);
   ia11 = fsls_CSRMatrixI(A11);
   ja11 = fsls_CSRMatrixJ(A11);     
   a11  = fsls_CSRMatrixData(A11);

   A22 = fsls_CSRMatrixCreate(n,n,nza22);
   fsls_CSRMatrixInitialize(A22);
   ia22 = fsls_CSRMatrixI(A22);
   ja22 = fsls_CSRMatrixJ(A22);     
   a22  = fsls_CSRMatrixData(A22);

   A33 = fsls_CSRMatrixCreate(n,n,nza33);
   fsls_CSRMatrixInitialize(A33);
   ia33 = fsls_CSRMatrixI(A33);
   ja33 = fsls_CSRMatrixJ(A33);     
   a33  = fsls_CSRMatrixData(A33);
  
  /*-----------------------------------------------------
   * Create 4 Vectors
   *----------------------------------------------------*/   
   V12 = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(V12);
   a12 = fsls_VectorData(V12);

   V21 = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(V21);
   a21 = fsls_VectorData(V21);
 
   V23 = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(V23);
   a23 = fsls_VectorData(V23);
  
   V32 = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(V32);
   a32 = fsls_VectorData(V32);
            
  /*-----------------------------------------------------
   * Initialization
   *----------------------------------------------------*/
   nza11 = 0; nza22 = 0; nza33 = 0;
   nzv12 = 0; nzv21 = 0; nzv23 = 0;  nzv32 = 0; 
   ia11[0] = 0; ia22[0] = 0; ia33[0] = 0;

  /*-----------------------------------------------------
   * Get the 7 submatrices
   *----------------------------------------------------*/
  
   for (row = 0; row < n; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col < n)
          {
             ja11[nza11] = col;
              a11[nza11] = a[j];
             nza11 ++;
          }
          else if (col >= n && col < nplusn)
          {
             a12[nzv12++] = a[j];
          }
          else
          {
             jxf_printf("\n\nA13 != 0\n\n");
          }
       }
       ia11[row+1] = nza11;
   }
 
   for (row = n; row < nplusn; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col < n)
          {
             a21[nzv21++] = a[j];
          }
          else if (col >= n && col < nplusn)
          {
             ja22[nza22] = col - n;
              a22[nza22] = a[j];
             nza22 ++;
          }
          else
          {
             a23[nzv23++] = a[j];
          }
       }
       ia22[row-n+1] = nza22;   
   }   
  
   for (row = nplusn; row < N; row ++)
   {
       for (j = ia[row]; j < ia[row+1]; j ++)
       {
          col = ja[j];
          if (col < n)
          {
             jxf_printf("\n\nA31 != 0\n\n");
          }
          else if (col >= nplusn)
          {
             ja33[nza33] = col - nplusn;
              a33[nza33] = a[j];
             nza33 ++;
          }
          else
          {
             a32[nzv32++] = a[j];
          }
       }
       ia33[row-nplusn+1] = nza33;
   }  
#if 0
   fsls_CSRMatrixPrint(A11, "Arr");
   fsls_CSRMatrixPrint(A22, "Aee");
   fsls_CSRMatrixPrint(A33, "Aii");
   fsls_SeqVectorPrint(V12, "Vre");     
   fsls_SeqVectorPrint(V21, "Ver");
   fsls_SeqVectorPrint(V23, "Vei");
   fsls_SeqVectorPrint(V32, "Vie");  
#endif
  
  /*------------------------------------------------------------
   * Generate B
   *-----------------------------------------------------------*/
   B = fsls_CSRMatrixCreate(N, N, fsls_CSRMatrixNumNonzeros(A));
   fsls_CSRMatrixInitialize(B);
   ib = fsls_CSRMatrixI(B);
   jb = fsls_CSRMatrixJ(B);     
   b  = fsls_CSRMatrixData(B);

   nzb   = 0;
   ib[0] = 0;
   
   for (row = 0; row < n; row ++)
   {
      // Aee
      for (j = ia22[row]; j < ia22[row+1]; j ++)
      {
         jb[nzb] = ja22[j];
          b[nzb] =  a22[j];
         nzb ++;
      }
      // Aei
       b[nzb] = a23[row];
      jb[nzb] = n + row;
      nzb ++;
      // Aer
       b[nzb] = a21[row];
      jb[nzb] = nplusn + row;
      nzb ++;
      
      ib[row+1] = nzb;
   }
   
   for (row = 0; row < n; row ++)
   {
      // Aii
      for (j = ia33[row]; j < ia33[row+1]; j ++)
      {
         jb[nzb] = ja33[j] + n;
          b[nzb] =  a33[j];
         nzb ++;
      }
      // Aie
      jb[nzb] = row;
       b[nzb] = a32[row];
      nzb ++;

      ib[row+n+1] = nzb;
   }
   
   
   for (row = 0; row < n; row ++)
   {
      // Arr
      for (j = ia11[row]; j < ia11[row+1]; j ++)
      {
         jb[nzb] = ja11[j] + nplusn;
          b[nzb] =  a11[j];
         nzb ++;
      }
      // Are
      jb[nzb] = row;
       b[nzb] = a12[row];
      nzb ++;

      ib[row+nplusn+1] = nzb;
   }
   
   fsls_CSRMatrixDestroy(A11);
   fsls_CSRMatrixDestroy(A22);
   fsls_CSRMatrixDestroy(A33);
   
   fsls_SeqVectorDestroy(V12);
   fsls_SeqVectorDestroy(V21);
   fsls_SeqVectorDestroy(V23);
   fsls_SeqVectorDestroy(V32);
   
   return (B); 
}

/*!
 * \fn JXF_Int fsls_CSRMatrixMatvec00
 * \brief Compute y := A*x, here, A is n*n matrix,
 *        so we don't check the size compatibility of A,x,y 
 * \param *A pointer to the matrix
 * \param *x pointer to the vector
 * \param *y pointer to the vector
 * \note This function has been called by PEV modul, so please
 *       modify the function names at the same time. 2011/05/26 
 * \author peghoty
 * \date 2010/05/01 
 */
JXF_Int
fsls_CSRMatrixMatvec00( fsls_CSRMatrix *A, fsls_Vector *x, fsls_Vector *y )
{
   JXF_Int      n  = fsls_CSRMatrixNumRows(A);
   JXF_Int     *ia = fsls_CSRMatrixI(A);
   JXF_Int     *ja = fsls_CSRMatrixJ(A);
   JXF_Real  *a  = fsls_CSRMatrixData(A);
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Real   temp;
   JXF_Int      i,j;

   // y = A*x
   for (i = 0; i < n; i ++)
   {
       temp = 0.0;
       for (j = ia[i]; j < ia[i+1]; j ++)
       {
          temp += a[j]*x_data[ja[j]];
       }
       y_data[i] = temp;
   }

   return (0);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixMatvec01
 * \brief Compute y := alpha*A*x + beta*y, here, A is n*n matrix, 
 *        so we don't check the size compatibility of A,x,y
 * \param alpha a real number
 * \param *A pointer to the matrix
 * \param *x pointer to the vector
 * \param beta a real number
 * \param *y pointer to the vector
 * \note x and y are vectors. 
 * \author peghoty
 * \date 2010/05/01 
 */
JXF_Int
fsls_CSRMatrixMatvec01( JXF_Real          alpha,
                        fsls_CSRMatrix *A,
                        fsls_Vector    *x,
                        JXF_Real          beta,
                        fsls_Vector    *y )
{
   JXF_Int      n  = fsls_CSRMatrixNumRows(A);
   JXF_Int     *ia = fsls_CSRMatrixI(A);
   JXF_Int     *ja = fsls_CSRMatrixJ(A);
   JXF_Real  *a  = fsls_CSRMatrixData(A);
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Real   temp;
   JXF_Int      i,j,k;

   JXF_Real *py   = NULL;
   JXF_Real  gama = 0.0;
	
   if (alpha == 0)
   {
      if (beta == 0)
      {
         memset(y_data, 0X0, n*sizeof(JXF_Real));
      }
      else
      {
         py = &y_data[0];
         for (k = 0; k < n; k++, py++)
         {
            *py *= beta;
         }
      }
      return (0);
   }
	
   gama = beta / alpha;
	   
   // y: = (beta / alpha)*y
   if (gama != 1.0)
   {
      if (gama == 0)
      {
         memset(y_data, 0X0, n*sizeof(JXF_Real));
      }
      else
      {
         py = &y_data[0];
         for (k = 0; k < n; k ++, py ++)
         {
            *py *= gama;
         }
      }
   }    

   // y += A*x
   for (i = 0; i < n; i ++)
   {
       temp = y_data[i];
       for (j = ia[i]; j < ia[i+1]; j ++)
       {
          temp += a[j]*x_data[ja[j]];
       }
       y_data[i] = temp;
   }
   
   // y: = alpha*y
   if (alpha != 1.0)
   {
      py = &y_data[0];
      for (k = 0; k < n; k ++, py ++)
      {
         *py *= alpha;
      }  
   }
      
   return (0);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixMatvec02
 * \brief Compute y := alpha*A*x + beta*y, here, A is n*n matrix, 
 *        so we don't check the size compatibility of A,x,y
 * \param alpha a real number
 * \param *A pointer to the matrix
 * \param *x pointer to the vector
 * \param beta a real number
 * \param *y pointer to the vector
 * \note x and y are JXF_Real arrays.
 * \author peghoty
 * \date 2010/12/07
 */
JXF_Int
fsls_CSRMatrixMatvec02( JXF_Real          alpha,
                        fsls_CSRMatrix *A,
                        JXF_Real         *x,
                        JXF_Real          beta,
                        JXF_Real         *y )
{
   JXF_Int      n  = fsls_CSRMatrixNumRows(A);
   JXF_Int     *ia = fsls_CSRMatrixI(A);
   JXF_Int     *ja = fsls_CSRMatrixJ(A);
   JXF_Real  *a  = fsls_CSRMatrixData(A);
   JXF_Real   temp;
   JXF_Int      i,j,k;

   JXF_Real *py   = NULL;
   JXF_Real  gama = 0.0;
	
   if (alpha == 0)
   {
      if (beta == 0)
      {
         memset(y, 0X0, n*sizeof(JXF_Real));
      }
      else
      {
         py = &y[0];
         for (k = 0; k < n; k++, py++)
         {
            *py *= beta;
         }
      }
      return (0);
   }
	
   gama = beta / alpha;
	   
   // y: = (beta / alpha)*y
   if (gama != 1.0)
   {
      if (gama == 0)
      {
         memset(y, 0X0, n*sizeof(JXF_Real));
      }
      else
      {
         py = &y[0];
         for (k = 0; k < n; k ++, py ++)
         {
            *py *= gama;
         }
      }
   }    

   // y += A*x
   for (i = 0; i < n; i ++)
   {
       temp = y[i];
       for (j = ia[i]; j < ia[i+1]; j ++)
       {
          temp += a[j]*x[ja[j]];
       }
       y[i] = temp;
   }
   
   // y: = alpha*y
   if (alpha != 1.0)
   {
      py = &y[0];
      for (k = 0; k < n; k ++, py ++)
      {
         *py *= alpha;
      }  
   }
      
   return (0);
}

/*!
 * \fn JXF_Int fsls_CSRMatrixMatvec
 * \brief Compute y := alpha*A*x + beta*y
 * \param alpha a real number
 * \param *A pointer to the matrix
 * \param *x pointer to the vector
 * \param beta a real number
 * \param *y pointer to the vector
 * \author peghoty
 * \date 2010/01/03  
 */
JXF_Int
fsls_CSRMatrixMatvec( JXF_Real          alpha,
                      fsls_CSRMatrix *A,
                      fsls_Vector    *x,
                      JXF_Real          beta,
                      fsls_Vector    *y     )
{
   JXF_Real  *A_data   = fsls_CSRMatrixData(A);
   JXF_Int     *A_i      = fsls_CSRMatrixI(A);
   JXF_Int     *A_j      = fsls_CSRMatrixJ(A);
   JXF_Int      num_rows = fsls_CSRMatrixNumRows(A);
   JXF_Int      num_cols = fsls_CSRMatrixNumCols(A);

   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      x_size = fsls_VectorSize(x);
   JXF_Int      y_size = fsls_VectorSize(y);

   JXF_Real   temp;
   JXF_Int      i, jj, ierr = 0;

  /*-----------------------------------------------------------------
   *  Check for size compatibility.  Matvec returns ierr = 1 if
   *  length of X doesn't equal the number of columns of A,
   *  ierr = 2 if the length of Y doesn't equal the number of rows
   *  of A, and ierr = 3 if both are true.
   *
   *  Because temporary vectors are often used in Matvec, none of 
   *  these conditions terminates processing, and the ierr flag
   *  is informational only.
   *----------------------------------------------------------------*/
   if (num_cols != x_size) ierr = 1;
   if (num_rows != y_size) ierr = 2;
   if (num_cols != x_size && num_rows != y_size) ierr = 3;

  /*----------------------------------------------------------
   * Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
   *---------------------------------------------------------*/
   if (alpha == 0.0)
   {
      for (i = 0; i < num_rows; i ++) y_data[i] *= beta;
      return ierr;
   }

  /*-------------------------------------------------
   * y = (beta/alpha)*y
   *------------------------------------------------*/
   temp = beta / alpha;
   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
	 for (i = 0; i < num_rows; i ++) y_data[i] = 0.0;
      }
      else
      {
	 for (i = 0; i < num_rows; i ++) y_data[i] *= temp;
      }
   }

  /*-----------------------------------------------------------------
   * y += A*x
   *-----------------------------------------------------------------*/
   for (i = 0; i < num_rows; i ++)
   {
      temp = y_data[i];
      for (jj = A_i[i]; jj < A_i[i+1]; jj ++)
      {
         temp += A_data[jj]*x_data[A_j[jj]];
      }
      y_data[i] = temp;
   }

  /*------------------------------------------
   * y = alpha*y
   *-----------------------------------------*/
   if (alpha != 1.0)
   {
      for (i = 0; i < num_rows; i ++)
      {
	 y_data[i] *= alpha;
      }
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_CSRMatrixMatvecT
 * \brief Compute y := alpha * A^T*x + beta*y.
 * \note From Van Henson's modification of fsls_CSRMatrixMatvec.
 * \date 2010/01/07
 */
JXF_Int
fsls_CSRMatrixMatvecT( JXF_Real          alpha,
                       fsls_CSRMatrix *A,
                       fsls_Vector    *x,
                       JXF_Real          beta,
                       fsls_Vector    *y     )
{
   JXF_Real  *A_data    = fsls_CSRMatrixData(A);
   JXF_Int     *A_i       = fsls_CSRMatrixI(A);
   JXF_Int     *A_j       = fsls_CSRMatrixJ(A);
   JXF_Int      num_rows  = fsls_CSRMatrixNumRows(A);
   JXF_Int      num_cols  = fsls_CSRMatrixNumCols(A);

   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      x_size = fsls_VectorSize(x);
   JXF_Int      y_size = fsls_VectorSize(y);

   JXF_Real   temp;
   JXF_Int      i, j, jj, ierr = 0;

  /*-----------------------------------------------------------------
   *  Check for size compatibility.  Matvec returns ierr = 1 if
   *  length of X doesn't equal the number of columns of A,
   *  ierr = 2 if the length of Y doesn't equal the number of rows
   *  of A, and ierr = 3 if both are true.
   *
   *  Because temporary vectors are often used in Matvec, none of 
   *  these conditions terminates processing, and the ierr flag
   *  is informational only.
   *----------------------------------------------------------------*/
   if (num_rows != x_size) ierr = 1;
   if (num_cols != y_size) ierr = 2;
   if (num_rows != x_size && num_cols != y_size) ierr = 3;

  /*---------------------------------------------------------------
   * Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
   *--------------------------------------------------------------*/
   if (alpha == 0.0)
   {
      for (i = 0; i < num_cols; i ++) y_data[i] *= beta;
      return ierr;
   }

  /*-----------------------------------
   * y = (beta/alpha)*y
   *---------------------------------*/
   temp = beta / alpha;
   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
	 for (i = 0; i < num_cols; i ++) y_data[i] = 0.0;
      }
      else
      {
	 for (i = 0; i < num_cols; i ++) y_data[i] *= temp;
      }
   }

  /*-----------------------------------------------------------------
   * y += A^T*x
   *-----------------------------------------------------------------*/
   for (i = 0; i < num_rows; i++)
   {
      for (jj = A_i[i]; jj < A_i[i+1]; jj ++)
      {
         j = A_j[jj];
         y_data[j] += A_data[jj]*x_data[i];
      }
   }

  /*-----------------------------------------------------------------
   * y = alpha*y
   *-----------------------------------------------------------------*/
   if (alpha != 1.0)
   {
      for (i = 0; i < num_cols; i ++) y_data[i] *= alpha;
   }

   return ierr;
}

/*!
 * \fn JXF_Real fsls_CSRMatrixVecTMatVec
 * \param *y pointer to the vector 
 * \brief Compute value := y^{T}*A*x
 * \param *A pointer to the matrix
 * \param *x pointer to the vector
 * \return value the resulting value 
 * \author peghoty
 * \date 2010/11/14 
 */
JXF_Real
fsls_CSRMatrixVecTMatVec( fsls_Vector    *y,
                          fsls_CSRMatrix *A,
                          fsls_Vector    *x  )
{
   JXF_Int      n  = fsls_CSRMatrixNumRows(A);
   JXF_Int     *ia = fsls_CSRMatrixI(A);
   JXF_Int     *ja = fsls_CSRMatrixJ(A);
   JXF_Real  *a  = fsls_CSRMatrixData(A);
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Real   temp, value = 0.0;
   JXF_Int      i,j;
	
   // value = y^{T}*A*x
   for (i = 0; i < n; i ++)
   {
       temp = 0.0;
       for (j = ia[i]; j < ia[i+1]; j ++)
       {
          temp += a[j]*x_data[ja[j]];
       }
       value += temp*y_data[i];
   }
      
   return value;
}

/*!
 * \fn JXF_Int fsls_CSRMatRhsRead
 * \brief Read the matrix and rhs files and obtain the 
 *        corresponding variales and arrays.
 * \param *MatFile pointer to the matrix file
 * \param *RhsFile pointer to the right hand side vector file
 * \param **ia_ptr IA of CSR format of matrix 'A'
 * \param **ja_ptr JA of CSR format of matrix 'A'
 * \param **a_ptr  DATA of CSR format of matrix 'A'
 * \param n_ptr size of matrix 'A'
 * \param nz_ptr number of nonzeros in 'A'
 * \param **f_ptr pointer to the right hand side vector
 *
 * \date 2009/11/29
 * \author peghoty
 */
JXF_Int 
fsls_CSRMatRhsRead( char    *MatFile,
                    char    *RhsFile,
                    JXF_Int    **ia_ptr,
                    JXF_Int    **ja_ptr,
                    JXF_Real **a_ptr,
                    JXF_Int     *n_ptr,
                    JXF_Int     *nz_ptr,
                    JXF_Real **f_ptr   )
{
    FILE   *fp = NULL;
    JXF_Int    *ia = NULL;
    JXF_Int    *ja = NULL;
    JXF_Real *a  = NULL;
    JXF_Real *f  = NULL;
    JXF_Int     j,m,n,nz;
    
    /* Read the matrix */
    if ( (fp = fopen(MatFile, "r")) == NULL )
    {
       jxf_printf("Can't open the file %s!\n",MatFile);
       exit(-1);
    }    
    jxf_fscanf(fp, "%d", &n);
    ia = fsls_CTAlloc(JXF_Int, n+1);
    for (j = 0; j < n+1; j ++)
    {
       jxf_fscanf(fp, "%d", &ia[j]);
       ia[j] --;
    }
    nz = ia[n];
    ja = fsls_CTAlloc(JXF_Int, nz);
    a  = fsls_CTAlloc(JXF_Real, nz);
    for (j = 0; j < nz; j ++)
    {
       jxf_fscanf(fp, "%d", &ja[j]);
       ja[j] --;
    }
    for (j = 0; j < nz; j ++)
    {
       jxf_fscanf(fp, "%le", &a[j]);
    }
    fclose(fp);
    
    /* Read the rhs vector */
    if ( (fp = fopen(RhsFile, "r")) == NULL )
    {
       jxf_printf("Can't open the file %s!\n",RhsFile);
       exit(-1);
    }  
    jxf_fscanf(fp, "%d", &m);
    if (m != n) 
    {
       jxf_printf(" >>> Warning: The sizes of matrix and rhs are not compatible!\n");
       exit(0);
    }
    f = fsls_CTAlloc(JXF_Real, m);
    for (j = 0; j < m; j ++)
    {
       jxf_fscanf(fp, "%le", &f[j]);
    }
    fclose(fp);
    
    /* return */
    *n_ptr  = n;
    *nz_ptr = nz;
    *ia_ptr = ia;
    *ja_ptr = ja;
    *a_ptr  = a;
    *f_ptr  = f;

    return (0);
}


/*!
 * \fn JXF_Int fsls_Arrays2CSRMatVec
 * \brief From arrays to CSR matrix and rhs vector, initial guess
 * \param *ia IA of CSR format of matrix 'A'
 * \param *ja JA of CSR format of matrix 'A'
 * \param *a  DATA of CSR format of matrix 'A'
 * \param n size of matrix 'A'
 * \param nz number of nonzeros in 'A'
 * \param *f pointer to the right hand side vector
 * \param *u pointer to the initial guess vector
 * \param **A_ptr pointer to the resulting matrix
 * \param **b_ptr pointer to the resulting right hand side
 * \param **x_ptr pointer to the resulting initial guess
 * \remark we allocate memories in this subroutine.
 * \author peghoty
 * \date 2009/11/29 
 */
JXF_Int
fsls_Arrays2CSRMatVec( JXF_Int    *ia,
                       JXF_Int    *ja,
                       JXF_Real *a,
                       JXF_Int     n,
                       JXF_Int     nz,
                       JXF_Real *f,
                       JXF_Real *u,
                       fsls_CSRMatrix **A_ptr,
                       fsls_Vector    **b_ptr,
                       fsls_Vector    **x_ptr )
{
    JXF_Int i;
    fsls_CSRMatrix *A = NULL;
    fsls_Vector    *b = NULL;
    fsls_Vector    *x = NULL;

    /* Generate the matrix */
    A = fsls_CSRMatrixCreate(n, n, nz);
    fsls_CSRMatrixInitialize(A);
    for (i = 0; i < n+1; i ++)
      fsls_CSRMatrixI(A)[i] = ia[i];
    for (i = 0; i < nz; i ++)
      fsls_CSRMatrixJ(A)[i] = ja[i];
    for (i = 0; i < nz; i ++)
      fsls_CSRMatrixData(A)[i] = a[i];
    fsls_CSRMatrixReorder(A); 
      
    /* Generate the rhs */
    b = fsls_SeqVectorCreate(n);
    fsls_SeqVectorInitialize(b);
    for (i = 0; i < n; i ++)
      fsls_VectorData(b)[i] = f[i]; 
      
    /* Generate the initial guess */
    x = fsls_SeqVectorCreate(n);
    fsls_SeqVectorInitialize(x);
    for (i = 0; i < n; i ++)
      fsls_VectorData(x)[i] = u[i];
    
   *A_ptr = A;
   *b_ptr = b;
   *x_ptr = x;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_Arrays2CSRMatVecNoAllocate
 * \brief From arrays to CSR matrix and rhs vector, initial guess
 * \param *ia IA of CSR format of matrix 'A'
 * \param *ja JA of CSR format of matrix 'A'
 * \param *a  DATA of CSR format of matrix 'A'
 * \param n size of matrix 'A'
 * \param nz number of nonzeros in 'A'
 * \param *f pointer to the right hand side vector
 * \param *u pointer to the initial guess vector
 * \param **A_ptr pointer to the resulting matrix
 * \param **b_ptr pointer to the resulting right hand side
 * \param **x_ptr pointer to the resulting initial guess
 * \remark we don't allocate memory in this subroutine.
 * \author peghoty
 * \date 2011/01/24 
 */
JXF_Int
fsls_Arrays2CSRMatVecNoAllocate( JXF_Int    *ia,
                                 JXF_Int    *ja,
                                 JXF_Real *a,
                                 JXF_Int     n,
                                 JXF_Int     nz,
                                 JXF_Real *f,
                                 JXF_Real *u,
                                 fsls_CSRMatrix **A_ptr,
                                 fsls_Vector    **b_ptr,
                                 fsls_Vector    **x_ptr )
{
    fsls_CSRMatrix *A = NULL;
    fsls_Vector    *b = NULL;
    fsls_Vector    *x = NULL;

    /* Generate the matrix */
    A = fsls_CSRMatrixCreate(n, n, nz);
    fsls_CSRMatrixI(A) = ia;
    fsls_CSRMatrixJ(A) = ja;
    fsls_CSRMatrixData(A) = a;
    
    fsls_CSRMatrixReorder(A);
          
    /* Generate the rhs */
    b = fsls_SeqVectorCreate(n);
    fsls_VectorData(b) = f; 
      
    /* Generate the initial guess */
    x = fsls_SeqVectorCreate(n);
    fsls_VectorData(x) = u;
    
   *A_ptr = A;
   *b_ptr = b;
   *x_ptr = x;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_ComputeResidual
 * \brief Compute the residual y := b-Ax
 * \param *ia IA of CSR format of matrix 'A'
 * \param *ja JA of CSR format of matrix 'A'
 * \param *a  DATA of CSR format of matrix 'A'
 * \param n size of matrix 'A'
 * \param *x the vector
 * \param *b the right hand side
 * \param *y y = b - Ax
 * \author peghoty
 * \date 2009/12/30
 */
JXF_Int
fsls_ComputeResidual( JXF_Int     *ia,
                      JXF_Int     *ja,
                      JXF_Real  *a,
                      JXF_Int      n,
                      JXF_Real  *x,
                      JXF_Real  *b,
                      JXF_Real  *y )
{
   JXF_Int i;
   
   fsls_MatVecMultiply(ia,ja,a,n,x,y); // y = A*x
   
   for (i = 0; i < n; i ++)
   {
      y[i] = b[i] - y[i];
   }
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_MatVecMultiply
 * \brief Compute y := Ax
 * \param *ia IA of CSR format of matrix 'A'
 * \param *ja JA of CSR format of matrix 'A'
 * \param *a  DATA of CSR format of matrix 'A'
 * \param n size of matrix 'A'
 * \param *x the vector
 * \param *y (y = Ax)
 * \author peghoty
 * \date 2009/12/30 
 */
JXF_Int
fsls_MatVecMultiply( JXF_Int     *ia,
                     JXF_Int     *ja,
                     JXF_Real  *a,
                     JXF_Int      n,
                     JXF_Real  *x,
                     JXF_Real  *y )
{
   JXF_Int i,j;
      
   for (i = 0; i < n; i ++)
   {
       y[i] = 0.0;
       for (j = ia[i]; j < ia[i+1]; j ++)
       {
          y[i] += a[j]*x[ja[j]];
       }  
   }
      
   return 0;
}

/*!
 * \fn fsls_Vector *fsls_CreateRhsByMatrix
 * \brief Construct a vector y = A*v, where v = (1,1,...,1)^T.
 * \param *A pointer to the matrix.
 * \author peghoty
 * \date 2010/01/11  
 */
fsls_Vector * 
fsls_CreateRhsByMatrix( fsls_CSRMatrix *A )
{
   JXF_Int i,j;
   JXF_Real tmp;

   JXF_Int     n  = fsls_CSRMatrixNumRows(A);   
   JXF_Int    *ia = fsls_CSRMatrixI(A);
   JXF_Real *a  = fsls_CSRMatrixData(A);
   
   fsls_Vector *x = NULL;

   x = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(x);

   for (i = 0; i < n; i ++)
   {
      tmp = 0.0;
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         tmp += a[j];
      }
      fsls_VectorData(x)[i] = tmp;
   } 

   return(x);
}

/*!
 * \fn JXF_Int fsls_LinearSystemSymDiagScale0
 * \brief Diag-Scaling of a linear system symmetrically.
 * \note diag = D^{-1/2}, where D = diag(A).
 * \note The original LS is kept.
 * \date 2011/07/03
 * \author peghoty
 */
JXF_Int
fsls_LinearSystemSymDiagScale0( fsls_CSRMatrix   *A,
                                fsls_Vector      *f,
                                fsls_CSRMatrix  **AA_ptr,
                                fsls_Vector     **ff_ptr,
                                JXF_Real          **diag_ptr )
{
   JXF_Int      n     = fsls_CSRMatrixNumRows(A);
   JXF_Int      nz    = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia     = NULL;
   JXF_Int    *ja     = NULL;
   JXF_Real *a      = NULL;   
   JXF_Real *diag   = NULL;
   JXF_Real *f_data = NULL;

   fsls_CSRMatrix  *AA = NULL;
   fsls_Vector     *ff = NULL;
                               
   JXF_Int i,j,k;

   // AA
   AA = fsls_CSRMatrixCreate(n, n, nz);
   fsls_CSRMatrixInitialize(AA);
   fsls_CSRMatrixCopy(A, AA, 1);
   ia = fsls_CSRMatrixI(AA);
   ja = fsls_CSRMatrixJ(AA);
   a  = fsls_CSRMatrixData(AA); 
   
   // ff
   ff = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(ff);
   fsls_SeqVectorCopy(f, ff);   
   f_data = fsls_VectorData(ff);  
   
   // diag
   diag = fsls_CSRMatrixDiagonalEntries(AA);
   
   // check the diagonal elements, and set diag := D^{-1/2}
   for (i = 0; i < n; i ++)
   {
      if (diag[i] > 0.0)
      {
         diag[i] = 1.0 / sqrt(diag[i]);
      }
      else
      {
         jxf_printf("\n Warning: diag[%d] = %le <= 0\n\n", i, diag[i]);
         exit(0);
      }
   }
    
   // AA = D^{-1/2}*A*D^{-1/2}: a_{ij} = a_{ij}*d_i*d_j 
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         k = ja[j];
         a[j] = a[j]*diag[i]*diag[k];
      }
   } 
   
   // ff = D^{-1/2}*f
   for (i = 0; i < n; i ++)
   {
      f_data[i] = f_data[i]*diag[i];
   }
   
   *AA_ptr   = AA;
   *ff_ptr   = ff;
   *diag_ptr = diag;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_LinearSystemSymDiagScale
 * \brief Diag-Scaling of a linear system symmetrically.
 * \note diag = D^{-1/2}, where D = diag(A).
 * \note The original LS is changed. 
 * \date 2011/07/03
 * \author peghoty
 */
JXF_Int
fsls_LinearSystemSymDiagScale( fsls_CSRMatrix   *A,
                               fsls_Vector      *f,
                               JXF_Real          **diag_ptr )
{
   JXF_Int      n     = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia     = NULL;
   JXF_Int    *ja     = NULL;
   JXF_Real *a      = NULL;   
   JXF_Real *diag   = NULL;
   JXF_Real *f_data = NULL;
                               
   JXF_Int i,j,k;

   ia = fsls_CSRMatrixI(A);
   ja = fsls_CSRMatrixJ(A);
   a  = fsls_CSRMatrixData(A);  
   f_data = fsls_VectorData(f);  
   
   // diag
   diag = fsls_CSRMatrixDiagonalEntries(A);
   
   // check the diagonal elements, and set diag := D^{-1/2}
   for (i = 0; i < n; i ++)
   {
      if (diag[i] > 0.0)
      {
         diag[i] = 1.0 / sqrt(diag[i]);
      }
      else
      {
         jxf_printf("\n Warning: diag[%d] = %le <= 0\n\n", i, diag[i]);
         exit(0);
      }
   }
    
   // A = D^{-1/2}*A*D^{-1/2}: a_{ij} = a_{ij}*d_i*d_j 
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         k = ja[j];
         a[j] = a[j]*diag[i]*diag[k];
      }
   } 
   
   // f = D^{-1/2}*f
   for (i = 0; i < n; i ++)
   {
      f_data[i] = f_data[i]*diag[i];
   }

   *diag_ptr = diag;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_LinearSystemUnSymDiagScale0
 * \brief Diag-Scaling of a linear system unsymmetrically.
 * \note diag = D^{-1}, where D = diag(A). 
 * \note The original LS is kept. 
 * \date 2011/07/03
 * \author peghoty
 */
JXF_Int
fsls_LinearSystemUnSymDiagScale0( fsls_CSRMatrix   *A,
                                  fsls_Vector      *f,
                                  fsls_CSRMatrix  **AA_ptr,
                                  fsls_Vector     **ff_ptr )
{
   JXF_Int      n     = fsls_CSRMatrixNumRows(A);
   JXF_Int      nz    = fsls_CSRMatrixNumNonzeros(A);
   JXF_Int    *ia     = NULL;
   JXF_Real *a      = NULL;   
   JXF_Real *diag   = NULL;
   JXF_Real *f_data = NULL;

   fsls_CSRMatrix  *AA = NULL;
   fsls_Vector     *ff = NULL;
                               
   JXF_Int i,j;

   // AA
   AA = fsls_CSRMatrixCreate(n, n, nz);
   fsls_CSRMatrixInitialize(AA);
   fsls_CSRMatrixCopy(A, AA, 1);
   ia = fsls_CSRMatrixI(AA);
   a  = fsls_CSRMatrixData(AA); 
   
   // ff
   ff = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(ff);
   fsls_SeqVectorCopy(f, ff);   
   f_data = fsls_VectorData(ff);  
   
   // diag
   diag = fsls_CSRMatrixDiagonalEntries(AA);
   
   // check the diagonal elements, and set diag := D^{-1}
   for (i = 0; i < n; i ++)
   {
      if (diag[i] != 0.0)
      {
         diag[i] = 1.0 / diag[i];
      }
      else
      {
         jxf_printf("\n Warning: diag[%d] = 0\n\n", i);
         exit(0);
      }
   }
    
   // AA = D^{-1}*A: a_{ij} = d_i*a_{ij} 
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         a[j] = a[j]*diag[i];
      }
   } 
   
   // ff = D^{-1}*f
   for (i = 0; i < n; i ++)
   {
      f_data[i] = f_data[i]*diag[i];
   }
   
   fsls_TFree(diag);
   
   *AA_ptr   = AA;
   *ff_ptr   = ff;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_LinearSystemUnSymDiagScale
 * \brief Diag-Scaling of a linear system unsymmetrically.
 * \note diag = D^{-1}, where D = diag(A).
 * \note The original LS is changed.   
 * \date 2011/07/03
 * \author peghoty
 */
JXF_Int
fsls_LinearSystemUnSymDiagScale( fsls_CSRMatrix   *A,
                                 fsls_Vector      *f  )
{
   JXF_Int      n     = fsls_CSRMatrixNumRows(A);
   JXF_Int    *ia     = NULL;
   JXF_Real *a      = NULL;   
   JXF_Real *diag   = NULL;
   JXF_Real *f_data = NULL;
                               
   JXF_Int i,j;

   ia = fsls_CSRMatrixI(A);
   a  = fsls_CSRMatrixData(A);   
   f_data = fsls_VectorData(f);  
   
   // diag
   diag = fsls_CSRMatrixDiagonalEntries(A);
   
   // check the diagonal elements, and set diag := D^{-1}
   for (i = 0; i < n; i ++)
   {
      if (diag[i] != 0.0)
      {
         diag[i] = 1.0 / diag[i];
      }
      else
      {
         jxf_printf("\n Warning: diag[%d] = 0\n\n", i);
         exit(0);
      }
   }
    
   // A = D^{-1}*A: a_{ij} = d_i*a_{ij} 
   for (i = 0; i < n; i ++)
   {
      for (j = ia[i]; j < ia[i+1]; j ++)
      {
         a[j] = a[j]*diag[i];
      }
   } 
   
   // f = D^{-1}*f
   for (i = 0; i < n; i ++)
   {
      f_data[i] = f_data[i]*diag[i];
   }
   
   fsls_TFree(diag);
   
   return (0);
}

/**
 * \fn fsls_DenseMatrix *fsls_DenseMatrixCreate
 * \brief Create a dense matrix 
 * \param size size of the dense matrix
 * \author peghoty
 * \date 2010/11/19
 */
fsls_DenseMatrix * 
fsls_DenseMatrixCreate( JXF_Int num_rows, JXF_Int num_cols )
{
   fsls_DenseMatrix *A = fsls_CTAlloc(fsls_DenseMatrix, 1);
   fsls_DenseMatrixNumRows(A) = num_rows;
   fsls_DenseMatrixNumCols(A) = num_cols;
   fsls_DenseMatrixData(A)    = NULL;
   return (A);
}

/**
 * \fn JXF_Int fsls_DenseMatrixInitialize
 * \brief Initialize a dense matrix 
 * \param size size of the dense matrix
 * \author peghoty
 * \date 2010/11/19
 */
JXF_Int 
fsls_DenseMatrixInitialize( fsls_DenseMatrix *A )
{
   JXF_Int num_rows = fsls_DenseMatrixNumRows(A);
   JXF_Int num_cols = fsls_DenseMatrixNumCols(A);
   JXF_Int ierr = 0;
   
   if (num_rows && num_cols && !fsls_DenseMatrixData(A))
   {
      fsls_DenseMatrixData(A) = fsls_CTAlloc(JXF_Real, num_rows*num_cols);
   }
   
   return (ierr);
}

/**
 * \fn JXF_Int fsls_DenseMatrixDestroy
 * \brief Destroy a dense matrix 
 * \param *A the dense matrix to be destroyed
 * \author peghoty
 * \date 2010/11/19
 */
JXF_Int 
fsls_DenseMatrixDestroy( fsls_DenseMatrix *A )
{
   JXF_Int ierr = 0;
   
   if (A)
   {
      if ( fsls_DenseMatrixData(A) )
      {
         fsls_TFree(fsls_DenseMatrixData(A));
      }
      fsls_TFree(A);
   }
   
   return (ierr);
}

/*!
 * \fn JXF_Int fsls_DenseMatrixPrint
 * \brief Write a full matrix into a given file 
 * \param JXF_Int row number of rows 
 * \param JXF_Int col number of colums
 * \param JXF_Real *full the row*col entries 
 * \param char *filename file name 
 * \note the matrix file is written in the following format:
 *       num_row
 *       num_col
 *       full[0]
 *       full[1]
 *        ...
 *       full[num_row*num_col-1]
 * \author peghoty
 * \date 2010/08/01  
 */
JXF_Int
fsls_DenseMatrixPrint( JXF_Int row, JXF_Int col, JXF_Real *full, char *filename )
{
   FILE *fp = NULL;
   JXF_Int i;

   fp = fopen(filename, "w");
   jxf_fprintf(fp, "%d\n", row);
   jxf_fprintf(fp, "%d\n", col);
   for (i = 0; i < row*col; i ++)
   {
      jxf_fprintf(fp, "%.15le\n", full[i]);
   }
   fclose(fp);

   return (0);   
}

/**
 * \fn void fsls_DenseymAx(JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y)
 * \brief Compute y := y - Ax, where 'A' is a n*n dense matrix 
 * \param n the dimension of the dense matrix
 * \param *A pointer to the n*n dense matrix
 * \param *x pointer to the JXF_Real vector with length n
 * \param *y pointer to the JXF_Real vector with length n
 * \author peghoty
 * \date 2010/10/25
 */
void 
fsls_DenseymAx( JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y )
{
   JXF_Int i,j,k;

   for (i = 0; i < n; i ++) 
   {
      k = i*n;
      for (j = 0; j < n; j ++) 
      {
         y[i] -= A[k+j]*x[j];
      }
   }
} 

/**
 * \fn void fsls_DenseaAxpby(JXF_Real alpha, JXF_Real beta, JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y)  
 * \brief Compute y:=alpha*A*x + beta*y, here x,y are vectors, A is a n*n full matrix.   
 * \param alpha a real number
 * \param beta  a real number
 * \param n the length of vector x and y
 * \param *A pointer to the JXF_Real vector which stands for a n*n full matrix 
 * \param *x pointer to the JXF_Real vector with length n
 * \param *y pointer to the JXF_Real vector with length n
 * \author peghoty
 * \date 2010/10/25
 */
void 
fsls_DenseaAxpby( JXF_Real alpha, JXF_Real beta, JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y )
{
   JXF_Int    i,j,k;
   JXF_Real tmp = 0.0;
   
   if (alpha == 0)
   {
      for (i = 0; i < n; i ++) y[i] *= beta;
      return;
   }
   
   //! y := (beta/alpha)y
   tmp = beta / alpha;
   if (tmp != 1.0)
   {
      for (i = 0; i < n; i ++) y[i] *= tmp;
   }
   
   //! y := y + A*x
   for (i = 0; i < n; i ++)
   {
      k = i*n;
      for (j = 0; j < n; j ++)
      {
         y[i] += A[k+j]*x[j];
      }
   }
   
   //! y := alpha*y
   if (alpha != 1.0)
   {
      for (i = 0; i < n; i ++) y[i] *= alpha;
   }
}

/**
 * \fn void fsls_DenseMatMul
 * \comput the matrix product of two small full matrices a and b, stored in c
 * \param *a pointer to the JXF_Real vector which stands a n*n matrix
 * \param *b pointer to the JXF_Real vector which stands a n*n matrix
 * \param *c pointer to the JXF_Real vector which stands a n*n matrix
 * \param n the dimension of the matrix
 * \note this function comes from FASP package 
 * \date 04/29/2010    
 */
void 
fsls_DenseMatMul( JXF_Real *a, JXF_Real *b, JXF_Int n, JXF_Real *c )
{ 
   JXF_Int i,j,k,l,m;
    
   for (i = 0; i < n; i ++)
   {
      k = i*n;
      for (j = 0; j < n; j ++)
      { 
         m = k + j; 
         c[m] = 0.0;
         for (l = 0; l < n; l ++)
         {
            c[m] += a[k+l]*b[l*n+j];
         }
      }
   }
}

/**
 * \fn void fsls_DenseMatvec
 * \brief Compute the product of a small full matrix 'a' and a vector 'b', stored in 'c'.
 * \param *a pointer to the JXF_Real vector which stands a n*n matrix
 * \param *b pointer to the JXF_Real vector with length n
 * \param n the dimension of the matrix
 * \param *c pointer to the JXF_Real vector with length n
 * \note this function comes from FASP package  
 * \date 2010/04/21
 */
void 
fsls_DenseMatvec( JXF_Real *a, JXF_Real *b, JXF_Int n, JXF_Real *c )
{ 
   JXF_Int i,j;
    
   for (i = 0; i < n; i ++)
   {
      c[i] = 0.0;
      for (j = 0; j < n; j ++)
      {
         c[i] += a[i*n+j]*b[j];
      }
   }
}

/**
 * \fn JXF_Int fsls_DenseInverse(JXF_Int pivot_type, JXF_Real *a, JXF_Int n)
 * \brief Compute the inverse of a small dense matrix which is stored in 
 *        a one-dimension array in row-major order.
 * \param pivot_type pivoting type 
 * \param *a pointer to the JXF_Real vector which stands a n*n matrix
 * \param n dimension of the matrix
 * \param *is pointer to the auxiliary array (when pivot_type = 1)
 * \param *js pointer to the auxiliary array (when pivot_type = 1,2)
 * \return *a as the inverse of the original matrix
 * 
 *  Description of Gauss-Jordan Complete Pivoting Method for Solving the Inverse of A Given Matrix
 *
 *  For k, from 0 to n-1, do the following:
 *     1. Find the max abs. entry of the (n-k+1) order matrix, mark its row num in IS(k), 
 *        and col num in JS(k), and exchange it with a(k,k).
 *     2. a(k,k) = 1.0/a(k,k)
 *     3. a(k,j) = a(k,j)*a(k,k),            j   = 0, 1, ..., n-1;   j != k
 *     4. a(i,j) = a(i,j)-a(i,k)*a(k,j),     i,j = 0, 1, ..., n-1; i,j != k
 *     5. a(i,k) = -a(i,k)*a(k,k),           i   = 0, 1, ..., n-1;   i != k
 *  For k, from n-1 to 0, do the following:
 *     1. exchange the k-th row with the JS(k)-th row
 *     2. exchange the k-th col with the IS(k)-th col
 * \note This function comes from FASP package, options of 'NOPIVOT' 
 *       and 'COLUMNPIVOT' are added by Yue Xiaoqiang
 * \date 2010/05/05  
 */
JXF_Int 
fsls_DenseInverse( JXF_Int pivot_type, JXF_Real *a, JXF_Int n, JXF_Int *is, JXF_Int *js )
{ 
    JXF_Int    i, j, k, l, u, v, in, kn, isn, ink;
    JXF_Real d, p;
    
    if (pivot_type == NOPIVOT)
    {
        for (k = 0; k < n; k ++)
        {
            kn = k * n;
            l  = kn + k;
            
            a[l] = 1.0 / a[l];
                        
            for (j = 0; j < n; j ++)
            {
                if (j != k)
                {
                    u    = kn + j;
                    a[u] = a[u] * a[l];
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                {
                    in  = i * n;
                    ink = in + k;
                    
                    for (j = 0; j < n; j ++)
                    {
                        if (j != k)
                        {
                            u    = in + j;
                            a[u] = a[u] - a[ink] * a[kn+j];
                        }
                    }
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                {
                    u    = i * n + k;
                    a[u] = - a[u] * a[l];
                }
            }
        }
    }
    else if (pivot_type == COLUMNPIVOT)
    {
        for (k = 0; k < n; k ++)
        {
            d  = 0.0;
            kn = k * n;
            
            for (i = k; i < n; i ++)
            {
                in = i * n;
                
                for (j = k; j < n; j ++)
                {
                    l = in + j;
                    p = fabs(a[l]);
                    
                    if (p > d)
                    {
                        d     = p;
                        is[k] = i;
                    }
                }
            }
            
            if ((d+1.0) == 1.0)
            {
                jxf_printf("\n warning: The matrix is singular\n\n");
                return 0;
            }
            
            if (is[k] != k)
            {
                isn = is[k] * n;
                
                for (j = 0; j < n; j ++)
                {
                    u = kn + j;
                    v = isn + j;
                    
                    p    = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            }
            
            l = kn + k;
            
            a[l] = 1.0 / a[l];
            
            for (j = 0; j < n; j++)
            {
                if (j != k)
                {
                    u    = kn + j;
                    a[u] = a[u] * a[l];
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                {
                    in  = i * n;
                    ink = in + k;
                    
                    for (j = 0; j < n; j ++)
                    {
                        if (j != k)
                        {
                            u    = in + j;
                            a[u] = a[u] - a[ink]*a[kn+j];
                        }
                    }
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                {
                    u    = i * n + k;
                    a[u] = -a[u] * a[l];
                }
            }
        }
        
        for (k = n-1; k >= 0; k --)
        {
            if (is[k] != k)
            {
                for (i = 0; i < n; i ++)
                {
                    in = i * n;
                    u  = in + k;
                    v  = in + is[k];
                    
                    p    = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            }
        }
    }
    else if (pivot_type == COMPLETEPIVOT)
    {
        for (k = 0; k < n; k ++)
        { 
            d  = 0.0;
            kn = k * n;
            
            for (i = k; i < n; i ++)
            {
                in = i * n;
                
                for (j = k; j < n; j ++)
                { 
                    l = in + j; 
                    p = fabs(a[l]);
                    
                    if (p > d) 
                    { 
                        d     = p; 
                        is[k] = i; 
                        js[k] = j;
                    }
                }
            }
            
            if (d + 1.0 == 1.0)
            {
                jxf_printf("\n warning: The matrix is singular\n\n");
                return 0;
            }
            
            if (is[k] != k)
            {
                isn = is[k] * n;
                
                for (j = 0; j < n; j ++)
                { 
                    u = kn + j;
                    v = isn + j;
                    
                    p    = a[u];
                    a[u] = a[v]; 
                    a[v] = p;
                }
            }
            
            if (js[k] !=k )
            {
                for (i = 0; i < n; i ++)
                { 
                    in = i * n;
                    u  = in + k; 
                    v  = in + js[k];
                    
                    p    = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            }
            
            l = kn + k;
            
            a[l] = 1.0 / a[l];
            
            for (j = 0; j < n; j ++)
            {
                if (j != k)
                { 
                    u    = kn + j;
                    a[u] = a[u] * a[l];
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                {
                    in  = i * n;
                    ink = in + k;
                    
                    for (j = 0; j < n; j ++)
                    {
                        if (j != k)
                        { 
                            u    = in + j;
                            a[u] = a[u] - a[ink] * a[kn+j];
                        }
                    }
                }
            }
            
            for (i = 0; i < n; i ++)
            {
                if (i != k)
                { 
                    u    = i * n + k; 
                    a[u] = -a[u] * a[l];
                }
            }
            
        }
        
        for (k = n-1; k >= 0; k --)
        { 
            if (js[k] != k)
            {
                kn  = k * n;
                isn = js[k] * n;
                
                for (j = 0; j < n; j ++)
                { 
                    u = kn + j; 
                    v = isn + j;
                    
                    p    = a[u]; 
                    a[u] = a[v]; 
                    a[v] = p;
                }
            }
            
            if (is[k] != k)
            {
                for (i = 0; i < n; i ++)
                { 
                    in = i * n;
                    u  = in + k;
                    v  = in + is[k];
                    
                    p    = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            }     
        }
    }

    return 1;
}

/*!
 * \fn JXF_Int fsls_CSR2FullMatrix
 * \brief Transfer a CSR matrix to a full matrix(stored in one-dimensional array row by row)
 * \param fsls_CSRMatrix *A pointer to the CSR matrix
 * \param JXF_Real **full_ptr pointer to the Full matrix
 * \author peghoty
 * \date 2010/05/05  
 */
JXF_Int 
fsls_CSR2FullMatrix( fsls_CSRMatrix *A, JXF_Real **full_ptr )
{   
   /* CSR information of A */
   JXF_Int     n   = fsls_CSRMatrixNumRows(A);
   JXF_Int     m   = fsls_CSRMatrixNumCols(A);
   JXF_Int    *ia  = fsls_CSRMatrixI(A);
   JXF_Int    *ja  = fsls_CSRMatrixJ(A);     
   JXF_Real *a   = fsls_CSRMatrixData(A);
   
   JXF_Real *full = NULL;
   JXF_Int size = n*m;
   JXF_Int row,start;
   JXF_Int k;
      
   full = fsls_CTAlloc(JXF_Real, size);
  
   for (row = 0; row < n; row ++)
   {
      start = m*row;
      for (k = ia[row]; k < ia[row+1]; k ++)
      {
         full[start+ja[k]] = a[k];
      }
   }
   
   *full_ptr = full; 
    
   return (0);
}


/*!
 * \fn JXF_Int fsls_Full2CSRMatrix
 * \brief Transfer a full matrix(stored in one-dimensional array row by row) to a CSR matrix 
 * \param JXF_Real *full pointer to the Full matrix 
 * \param fsls_CSRMatrix **A_ptr pointer to pointer to the result CSR matrix
 * \note work for nXn full matrix
 * \date 2010/06/21
 * \author peghoty 
 */
JXF_Int 
fsls_Full2CSRMatrix( JXF_Int n, JXF_Real *full, fsls_CSRMatrix **A_ptr )
{   
   fsls_CSRMatrix *A = NULL;
   JXF_Int     *ia = NULL;
   JXF_Int     *ja = NULL;
   JXF_Real  *a  = NULL; 
   
   JXF_Int i,j,k;
   JXF_Int nz = 0;  

   ia = fsls_CTAlloc(JXF_Int, n+1);
    
   ia[0] = 0; 
   for (i = 0; i < n; i ++)
   {
      for (j = 0; j < n; j ++)
      {
         k = i*n + j;
         if (full[k] != 0.0)
         {
            nz ++;
         }
      }
      ia[i+1] = nz;
   }

   ja = fsls_CTAlloc(JXF_Int, nz);
   a  = fsls_CTAlloc(JXF_Real, nz);
   
   nz = 0;
   for (i = 0; i < n; i ++)
   {
      for (j = 0; j < n; j ++)
      {
         k = i*n + j;
         if (full[k] != 0.0)
         {
            ja[nz] = j;
            a[nz]  = full[k];
            nz ++;
         }
      }
   } 

   A = fsls_CSRMatrixCreate(n, n, nz);
   fsls_CSRMatrixI(A) = ia;
   fsls_CSRMatrixJ(A) = ja;
   fsls_CSRMatrixData(A) = a;
   
   *A_ptr = A;
   
   return (0);
}

/*!
 * \fn JXF_Int fsls_BuildVecFromFile
 * \brief Build a fsls_Vector from a given file.
 * \param *filename pointer to the vector file
 * \param **b_ptr pointer to the resulting vector
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_BuildVecFromFile( char *filename, fsls_Vector **b_ptr )
{
   fsls_Vector *b;
   
   b = fsls_SeqVectorRead(filename);
  *b_ptr = b;
 
   return (0);
}

/*!
 * \fn fsls_Vector * fsls_SeqVectorRead
 * \brief Build a fsls_Vector from a given file.
 * \param *filename pointer to the vector file
 * \note File format:
 *	 nrow
 *       val(j), j=0:nrow-1 Long E
 * \return a fsls_Vector
 * \author peghoty
 * \date 2010/01/07 
 */
fsls_Vector *
fsls_SeqVectorRead( char *file_name )
{
   fsls_Vector  *vector;
   FILE    *fp = NULL;
   JXF_Real  *data;
   JXF_Int      size;
   JXF_Int      j;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &size);

   vector = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(vector);

   data = fsls_VectorData(vector);
   for (j = 0; j < size; j ++)
   {
      jxf_fscanf(fp, "%le", &data[j]);
   }

   fclose(fp);

   return vector;
}

/*!
 * \fn fsls_IVector * fsls_IVectorRead
 * \brief Build a fsls_IVector from a given file.
 * \param *filename pointer to the vector file
 * \note File format:
 *	 nrow
 *       val(j), j=0:nrow-1 JXF_Int
 * \return a fsls_IVector
 * \author peghoty
 * \date 2010/12/09 
 */
fsls_IVector *
fsls_IVectorRead( char *file_name )
{
   fsls_IVector  *vector;
   FILE    *fp = NULL;
   JXF_Int     *data;
   JXF_Int      size;
   JXF_Int      j;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &size);

   vector = fsls_IVectorCreate(size);
   fsls_IVectorInitialize(vector);

   data = fsls_IVectorData(vector);
   for (j = 0; j < size; j ++)
   {
      jxf_fscanf(fp, "%d", &data[j]);
   }

   fclose(fp);

   return vector;
}

/*!
 * \fn fsls_Vector *fsls_IJVectorRead
 * \brief Build a fsls_Vector from a given file.
 * \param *filename pointer to the vector file
 * \note File format:
 *	 nrow
 *       ind(j) val(j), j=0:nrow-1 Long E 
 * \return a fsls_Vector
 * \author peghoty
 * \date 2010/11/14 
 */
fsls_Vector *
fsls_IJVectorRead( char *file_name )
{
   fsls_Vector  *vector;
   FILE    *fp = NULL;
   JXF_Real  *data;
   JXF_Int      size;
   JXF_Int      j;
   JXF_Int      index;
   JXF_Real   value;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &size);

   vector = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(vector);

   data = fsls_VectorData(vector);
   for (j = 0; j < size; j ++)
   {
      jxf_fscanf(fp, "%d %le", &index, &value);
      data[index] = value;
   }

   fclose(fp);

   return vector;
}

/*!
 * \fn JXF_Int fsls_SeqVectorPrint
 * \brief Write a fsls_Vector into a given file.
 * \param *vector pointer to the vector to be kept
 * \param *file_name pointer to the vector file
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_SeqVectorPrint( fsls_Vector *vector, char *file_name )
{
   FILE    *fp = NULL;

   JXF_Real  *data = fsls_VectorData(vector);
   JXF_Int      size = fsls_VectorSize(vector);
   
   JXF_Int      i,ierr = 0;

   fp = fopen(file_name, "w");

   jxf_fprintf(fp, "%d\n", size);

   for (i = 0; i < size; i ++)
   {
      jxf_fprintf(fp, "%.15le\n", data[i]);
   }

   fclose(fp);

   return ierr;
}

/*!
 * \fn JXF_Int fsls_IVectorPrint
 * \brief Write a fsls_IVector into a given file.
 * \param *vector pointer to the vector to be kept
 * \param *file_name pointer to the vector file
 * \author peghoty
 * \date 2010/11/16 
 */
JXF_Int
fsls_IVectorPrint( fsls_IVector *vector, char *file_name )
{
   FILE    *fp = NULL;

   JXF_Int     *data = fsls_IVectorData(vector);
   JXF_Int      size = fsls_IVectorSize(vector);
   
   JXF_Int      i,ierr = 0;

   fp = fopen(file_name, "w");

   jxf_fprintf(fp, "%d\n", size);

   for (i = 0; i < size; i ++)
   {
      jxf_fprintf(fp, "%d\n", data[i]);
   }

   fclose(fp);

   return ierr;
}

/*!
 * \fn JXF_Int fsls_SeqVectorSetConstantValues
 * \brief Set a fsls_Vector to be a constant vector.
 * \param *x pointer to the vector to be set
 * \param value the constant value
 * \author peghoty
 * \date 2010/01/07  
 */
JXF_Int
fsls_SeqVectorSetConstantValues( fsls_Vector *x, JXF_Real value )
{
   JXF_Real  *vector_data = fsls_VectorData(x);
   JXF_Int      size        = fsls_VectorSize(x);
           
   JXF_Int      i, ierr = 0;

   for (i = 0; i < size; i ++)
   {
      vector_data[i] = value;
   }
   
   return ierr;
}

/*!
 * \fn JXF_Int fsls_IVectorSetConstantValues
 * \brief Set a fsls_IVector to be a constant vector.
 * \param *x pointer to the vector to be set
 * \param value the constant value
 * \author peghoty
 * \date 2010/11/14  
 */
JXF_Int
fsls_IVectorSetConstantValues( fsls_IVector *x, JXF_Int value )
{
   JXF_Int *vector_data = fsls_IVectorData(x);
   JXF_Int  size        = fsls_IVectorSize(x);
           
   JXF_Int  i, ierr = 0;

   for (i = 0; i < size; i ++)
   {
      vector_data[i] = value;
   }
   
   return ierr;
}

/*!
 * \fn JXF_Int fsls_IVectorSetNatureOrder
 * \brief Set a fsls_IVector to be a nature order.
 * \param *x pointer to the vector to be set
 * \author peghoty
 * \date 2010/11/14  
 */
JXF_Int
fsls_IVectorSetNatureOrder( fsls_IVector *x )
{
   JXF_Int *vector_data = fsls_IVectorData(x);
   JXF_Int  size        = fsls_IVectorSize(x);
           
   JXF_Int  i, ierr = 0;

   for (i = 0; i < size; i ++)
   {
      vector_data[i] = i;
   }
   
   return ierr;
}

/*!
 * \fn fsls_Vector * fsls_SeqVectorCreate
 * \brief Create a fsls_Vector with 'size' as its length.
 * \param size length of the fsls_Vector
 * \return a fsls_Vector
 * \author peghoty
 * \date 2010/01/07 
 */
fsls_Vector *
fsls_SeqVectorCreate( JXF_Int size )
{
   fsls_Vector  *vector = NULL;

   vector = fsls_CTAlloc(fsls_Vector, 1);

   fsls_VectorData(vector) = NULL;
   fsls_VectorSize(vector) = size;

   return vector;
}

/*!
 * \fn fsls_IVector * fsls_IVectorCreate
 * \brief Create a fsls_IVector with 'size' as its length.
 * \param size length of the fsls_IVector
 * \return a fsls_IVector
 * \author peghoty
 * \date 2010/11/14 
 */
fsls_IVector *
fsls_IVectorCreate( JXF_Int size )
{
   fsls_IVector  *vector = NULL;

   vector = fsls_CTAlloc(fsls_IVector, 1);

   fsls_IVectorData(vector) = NULL;
   fsls_IVectorSize(vector) = size;

   return vector;
}

/*!
 * \fn JXF_Int fsls_SeqVectorDestroy
 * \brief Destroy a fsls_Vector.
 * \param *vector pointer to the fsls_Vector
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int 
fsls_SeqVectorDestroy( fsls_Vector *vector )
{
   JXF_Int  ierr = 0;

   if (vector)
   {
      if (fsls_VectorData(vector))
      {
         fsls_TFree(fsls_VectorData(vector));
      }
      fsls_TFree(vector);
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_IVectorDestroy
 * \brief Destroy a fsls_IVector.
 * \param *vector pointer to the fsls_IVector
 * \author peghoty
 * \date 2010/11/14 
 */
JXF_Int 
fsls_IVectorDestroy( fsls_IVector *vector )
{
   JXF_Int  ierr = 0;

   if (vector)
   {
      if (fsls_IVectorData(vector))
      {
         fsls_TFree(fsls_IVectorData(vector));
      }
      fsls_TFree(vector);
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_SeqVectorInitialize
 * \brief Initialize a fsls_Vector.
 * \param *vector pointer to the fsls_Vector to be initialized
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int 
fsls_SeqVectorInitialize( fsls_Vector *vector )
{
   JXF_Int  size = fsls_VectorSize(vector);
   JXF_Int  ierr = 0;

   if ( ! fsls_VectorData(vector) )
   {
      fsls_VectorData(vector) = fsls_CTAlloc(JXF_Real, size);
   }
   
   return ierr;
}

/*!
 * \fn JXF_Int fsls_IVectorInitialize
 * \brief Initialize a fsls_IVector.
 * \param *vector pointer to the fsls_IVector to be initialized
 * \author peghoty
 * \date 2010/11/14 
 */
JXF_Int 
fsls_IVectorInitialize( fsls_IVector *vector )
{
   JXF_Int  size = fsls_IVectorSize(vector);
   JXF_Int  ierr = 0;

   if ( ! fsls_IVectorData(vector) )
   {
      fsls_IVectorData(vector) = fsls_CTAlloc(JXF_Int, size);
   }
   
   return ierr;
}

/*!
 * \fn JXF_Real fsls_SeqVectorInnerProd
 * \brief Compute the inner product of x and y.
 * \param *x pointer to the first vector
 * \param *y pointer to the second vector
 * \return the inner product (x,y)
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Real   
fsls_SeqVectorInnerProd( fsls_Vector *x, fsls_Vector *y )
{
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      size   = fsls_VectorSize(x);
           
   JXF_Int      i;
   JXF_Real   result = 0.0;

   for (i = 0; i < size; i ++)
   {
      result += y_data[i]*x_data[i];
   }
   
   return result;
}

/*!
 * \fn JXF_Int fsls_SeqVectorCopy
 * \brief Copy data from x to y (x -> y)
 * \note y should have already been initialized
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_SeqVectorCopy( fsls_Vector *x, fsls_Vector *y )
{
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      size   = fsls_VectorSize(x);  /* here, it's dangerous to use fsls_VectorSize(y) 
                                             since sometimes y is longer than x, e.g. in AMG */
           
   JXF_Int      i,ierr = 0;

   for (i = 0; i < size; i ++)
   {
      y_data[i] = x_data[i];
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_SeqVectorSetRandomValues
 * \brief returns vector of values randomly distributed between 0.0 and 1.0
 * \note v should have already been initialized
 * \author peghoty
 * \date 2010/09/05 
 */
JXF_Int
fsls_SeqVectorSetRandomValues( fsls_Vector *v, JXF_Int seed )
{
   JXF_Real  *vector_data = fsls_VectorData(v);
   JXF_Int      size        = fsls_VectorSize(v);
           
   JXF_Int      i, ierr  = 0;
   
   fsls_SeedRand(seed);

   for (i = 0; i < size; i++)
   {
      vector_data[i] = fsls_Rand();
      //vector_data[i] = 2.0*fsls_Rand() - 1.0;
   }

   return ierr;
}

/*!
 * \fn JXF_Int fsls_SeqVectorScale
 * \brief Compute y := alpha*y
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_SeqVectorScale( JXF_Real alpha, fsls_Vector *y )
{
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      size   = fsls_VectorSize(y);
           
   JXF_Int      i,ierr = 0;

   for (i = 0; i < size; i ++)
   {
      y_data[i] *= alpha;
   }
   
   return ierr;
}

/*!
 * \fn fsls_Vector * fsls_SeqVectorAllOne
 * \brief Create a vector x = (1,1,...,1)^T.
 * \param n length of the vector
 * \author peghoty
 * \date 2010/01/11 
 */
fsls_Vector *
fsls_SeqVectorAllOne( JXF_Int n )
{
   fsls_Vector *x = NULL;

   x = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(x);
   fsls_SeqVectorSetConstantValues(x,1.0);
   
   return x;
}

/*!
 * \fn fsls_Vector *fsls_SeqVectorOne2N
 * \brief Create a vector x = (1,2,...,n)^T.
 * \param n length of the vector
 * \author peghoty
 * \date 2011/01/02 
 */
fsls_Vector *
fsls_SeqVectorOne2N( JXF_Int n )
{
   JXF_Int i;
   fsls_Vector *x = NULL;
   JXF_Real *x_data = NULL;

   x = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(x);
   x_data = fsls_VectorData(x);
   
   for (i = 0; i < n; i ++)
   {
      x_data[i] = (JXF_Real)(i + 1.0);
   }
   
   return x;
}

/*!
 * \fn JXF_Int fsls_SeqVectorAxpy
 * \brief Compute y := y + alpha*x
 * \author peghoty
 * \date 2010/01/07 
 */
JXF_Int
fsls_SeqVectorAxpy( JXF_Real alpha, fsls_Vector *x, fsls_Vector *y )
{
   JXF_Real  *x_data = fsls_VectorData(x);
   JXF_Real  *y_data = fsls_VectorData(y);
   JXF_Int      size   = fsls_VectorSize(y); /* here, it's dangerous to use fsls_VectorSize(x) 
                                            since sometimes x is longer than y, e.g. in AMG */
           
   JXF_Int      i,ierr = 0;

   for (i = 0; i < size; i ++)
   {
      y_data[i] += alpha*x_data[i];
   }

   return ierr;
}


/*!
 * \fn JXF_Int fsls_SeqVectorsDiff
 * \brief Compare two vectors
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_SeqVectorsDiff( fsls_Vector *x, fsls_Vector *y )
{
   JXF_Int     n1     = fsls_VectorSize(x);
   JXF_Int     n2     = fsls_VectorSize(y);
   JXF_Real *x_data = fsls_VectorData(x);
   JXF_Real *y_data = fsls_VectorData(y);
   
   JXF_Real *dif    = NULL;  // dif[i]    = |x[i] - y[i]|
   JXF_Real *reldif = NULL;  // reldif[i] = dif[i] / |x[i]|,  if x[i]!= 0.0
                           //           = dif[i] / |y[i]|,  if x[i] = 0.0 and y[i] != 0.0
                           //           = 0.0            ,  if x[i] = y[i] = 0.0 
   JXF_Real difmax,reldifmax,difl2norm_common,difl2norm_dividp;
   JXF_Int    difmax_index,reldifmax_index;
   
   JXF_Int i,n;
   
   if (n1 != n2)
   {
      jxf_printf("\n >>> Warning: Size Uncompatible !\n\n");
      return (-1);
   }
   else
   {
      n = n1;
   }
   
   dif    = fsls_CTAlloc(JXF_Real, n);
   reldif = fsls_CTAlloc(JXF_Real, n);
   
   /* generate dif and reldif */
   for (i = 0; i < n; i ++)
   {
      dif[i] = fabs(x_data[i] - y_data[i]);
      if (x_data[i])
      {
         reldif[i] = dif[i] / fabs(x_data[i]);
      }
      else if (y_data[i])
      {
         reldif[i] = dif[i] / (y_data[i]);
      }
      else
      {
         reldif[i] = 0.0;
      }
   }
   
   /* get difmax, reldifmax, difmax_index, reldifmax_index */
   difmax_index = fsls_ArrayDoubleAbsMax(dif, n, &difmax);
   reldifmax_index = fsls_ArrayDoubleAbsMax(reldif, n, &reldifmax);
   
   /* get difl2norm_common and difl2norm_dividp */
   difl2norm_common = fsls_Arrayl2Norm(dif, n);
   difl2norm_dividp = fsls_ArrayAveL2Norm(dif, n);
   
   /* output the results */
   jxf_printf("\n++++++++++++++++++++++++++++++++++\n");
   jxf_printf(" >>> Vector Difference Testing     \n");
   jxf_printf("++++++++++++++++++++++++++++++++++\n");
   jxf_printf(" difmax = %.2le\n", difmax);
   jxf_printf(" index  = %d\n", difmax_index);
   jxf_printf(" x[%d] = %.20le\n",difmax_index,x_data[difmax_index]);
   jxf_printf(" y[%d] = %.20le\n\n",difmax_index,y_data[difmax_index]);
   jxf_printf(" difl2norm_common = %.4le\n",difl2norm_common);
   jxf_printf(" difl2norm_dividp = %.4le\n\n",difl2norm_dividp);
   jxf_printf(" reldifmax = %.2le\n", reldifmax);
   jxf_printf(" index     = %d\n", reldifmax_index);
   jxf_printf(" x[%d] = %.20le\n",reldifmax_index,x_data[reldifmax_index]);
   jxf_printf(" y[%d] = %.20le\n\n",reldifmax_index,y_data[reldifmax_index]);
   
   /* free some staff */
   fsls_TFree(dif);
   fsls_TFree(reldif);
   
   return 1;
}

/*!
 * \fn JXF_Int fsls_SeqVectorFilesDiff
 * \brief compare two vectors from files
 * \param *file1 pointer to the first vector to be compared
 * \param *file2 pointer to the second vector to be compared
 * \author peghoty
 * \date 2010/08/26
 */
JXF_Int
fsls_SeqVectorFilesDiff( char *file1, char *file2 )
{
   FILE *fp1 = NULL;
   FILE *fp2 = NULL;

   fsls_Vector *x = NULL;
   fsls_Vector *y = NULL; 
   JXF_Real *x_data = NULL;
   JXF_Real *y_data = NULL; 
      
   JXF_Int i,n,n1,n2;
   
   /* compatibility checking */
   if ((fp1 = fopen(file1,"r")) == NULL)
   {
      jxf_printf("Can't open the %s!\n",file1);
      exit(0);
   }
   jxf_fscanf(fp1,"%d\n",&n1);
   
   if ((fp2 = fopen(file2,"r")) == NULL)
   {
      jxf_printf("Can't open the %s!\n",file2);
      exit(0);
   }
   jxf_fscanf(fp2,"%d\n",&n2);
   
   if (n1 != n2) 
   {
      jxf_printf("\n >>> Warning: Size Uncompatible !\n\n");
      fclose(fp1); fclose(fp2);
      return (-1);
   }
   else
   {
      n = n1;
   }  
   
   /* read the two vectors from file */
   x = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(x);
   y = fsls_SeqVectorCreate(n);
   fsls_SeqVectorInitialize(y);
   x_data = fsls_VectorData(x);  
   y_data = fsls_VectorData(y);    
   
   for (i = 0; i < n; i ++)
   {   
      jxf_fscanf(fp1, "%le\n", &x_data[i]);
      jxf_fscanf(fp2, "%le\n", &y_data[i]);
   }
   fclose(fp1);
   fclose(fp2);
   
   /* compare the two vectors */
   fsls_SeqVectorsDiff(x,y);
   
   fsls_SeqVectorDestroy(x);
   fsls_SeqVectorDestroy(y);
   
   return 1;
}

/*!
 * \fn JXF_Real fsls_SeqVectorAveL2Norm
 * \brief Compute the average l^2 norm of x 
 * \param x pointer to the vector to be computed
 * \author peghoty
 * \date 2010/01/05 
 */
JXF_Real 
fsls_SeqVectorAveL2Norm( fsls_Vector *x )
{
   JXF_Int     n = fsls_VectorSize(x);
   JXF_Real  avel2norm = 0.0;
   avel2norm = sqrt(fsls_SeqVectorInnerProd(x,x)/(JXF_Real)n);
   return(avel2norm);
}

/*!
 * \fn JXF_Real fsls_SeqVectorL2Norm
 * \brief Compute the l^2 norm of a given vector
 * \param *x pointer to the vector
 * \return l2norm, the l^2 norm of x
 * \author peghoty
 * \date 2009/12/30 
 */
JXF_Real   
fsls_SeqVectorL2Norm( fsls_Vector *x )
{
   JXF_Real l2norm;
   l2norm = sqrt(fsls_SeqVectorInnerProd(x,x));
   return l2norm;
}

/*!
 * \fn JXF_Real fsls_SeqVectorDiagInvNorm
 * \brief Compute the D^{-1} norm of a given vector, where D=diag(A).
 * \param *diag pointer to the diagonal element in D
 * \param *x pointer to the vector
 * \return diaginvnorm, the D^{-1} norm norm of x
 * \author peghoty
 * \date 2011/07/08 
 */
JXF_Real   
fsls_SeqVectorDiagInvNorm( JXF_Real *diag, fsls_Vector *x )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x); 
   
   JXF_Real diaginvnorm = 0.0;
   JXF_Int    i;   
       
   for (i = 0; i < size; i ++)
   {
      if (diag[i])
      {
         diaginvnorm += x_data[i]*x_data[i]/diag[i];
      }
   }
   
   if (diaginvnorm >= 0.0)
   {
      diaginvnorm = sqrt(diaginvnorm);
   }
   else
   {
      jxf_printf("\n >> Warning: diaginvnorm < 0 in fsls_SeqVectorDiagInvNorm!!\n\n");
      exit(0);
   }
   return diaginvnorm;
}

/*!
 * \fn fsls_Vector * fsls_SeqVectorReorderByDirection
 * \brief Reorder a vector which is ordered by points by directions.
 * \author peghoty
 * \date 2010/04/06
 */
fsls_Vector *
fsls_SeqVectorReorderByDirection( fsls_Vector *x )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
  
   fsls_Vector *y = NULL;
   JXF_Real *y_data = NULL;
   
   /* local variables */
   JXF_Int i, n = size / 3, n2 = 2*n;
   
   /* Create a new vector */   
   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);
   y_data = fsls_VectorData(y);

   /* Reordering */
   for (i = 0; i < n; i ++)
   {
      y_data[i] = x_data[3*i];  
   }
   for (i = n; i < n2; i ++)
   {
      y_data[i] = x_data[(i-n)*3 + 1];  
   }   
   for (i = n2; i < size; i ++)
   {
      y_data[i] = x_data[(i-n2)*3 + 2];  
   }  
   
   return (y);
}

/*!
 * \fn JXF_Int fsls_SeqVectorReorderByPoint
 * \brief Reorder a vector which is ordered by directions by points.
 * \author peghoty
 * \date 2010/04/06
 */
JXF_Int
fsls_SeqVectorReorderByPoint( fsls_Vector *x )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
  
   fsls_Vector *y = NULL;
   JXF_Real *y_data = NULL;
   
   /* local variables */
   JXF_Int i;
   JXF_Int n = size / 3;
   JXF_Int n2 = 2*n;
   
   /* Create a new vector */   
   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);
   y_data = fsls_VectorData(y);

   /* Looping for points */
   for (i = 0; i < n; i ++)
   {
      y_data[3*i]   = x_data[i];
      y_data[3*i+1] = x_data[i+n];
      y_data[3*i+2] = x_data[i+n2];
   }  
   
   for (i = 0; i < size; i ++)
   {   
      x_data[i] = y_data[i];
   }

   fsls_SeqVectorDestroy(y);

   return (0);
}

/*!
 * \fn JXF_Int fsls_SeqVectorRecoverByPoint
 * \brief Recover x which is ordered by directions to y which is ordered by points.
 * \author peghoty
 * \date 2011/05/15
 */
JXF_Int
fsls_SeqVectorRecoverByPoint( fsls_Vector *x, fsls_Vector *y )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
   JXF_Real *y_data = fsls_VectorData(y);
   
   /* local variables */
   JXF_Int i;
   JXF_Int n = size / 3;
   JXF_Int n2 = 2*n;

   /* Looping for points */
   for (i = 0; i < n; i ++)
   {
      y_data[3*i]   = x_data[i];
      y_data[3*i+1] = x_data[i+n];
      y_data[3*i+2] = x_data[i+n2];
   }  

   return (0);
}

/*!
 * \fn JXF_Int fsls_VecGroupCreate
 * \brief Create a fsls_VecGroup.  
 * \author peghoty
 * \date 2010/11/16
 */
fsls_VecGroup *
fsls_VecGroupCreate( JXF_Int n )
{
   fsls_VecGroup *vecgroup = fsls_CTAlloc(fsls_VecGroup, 1);
   fsls_VecGroupN(vecgroup)       = n;
   fsls_VecGroupXArray(vecgroup)  = NULL;
   return vecgroup;
} 
 
/*!
 * \fn JXF_Int fsls_VecGroupInitialize
 * \brief Initialize a fsls_VecGroup.  
 * \author peghoty
 * \date 2010/11/16
 */
JXF_Int
fsls_VecGroupInitialize( fsls_VecGroup *vecgroup )
{
   JXF_Int  ierr = 0;
   JXF_Int  n = fsls_VecGroupN(vecgroup);
   
   if ( n && !fsls_VecGroupXArray(vecgroup) )
   { 
      fsls_VecGroupXArray(vecgroup) = fsls_CTAlloc(fsls_Vector *, n);
   }
   
   return ierr;
} 

/*!
 * \fn JXF_Int fsls_VecGroupCombine
 * \brief Diagonally Combine a series of Vector.  
 *         |  x1   |
 *         |  x2   |
 *  x  =   |  ...  |
 *         |x_{n-1}| 
 *         | x_{n} |
 * \author peghoty
 * \date 2010/11/16
 */
fsls_Vector *
fsls_VecGroupCombine( fsls_VecGroup *vecgroup )
{
   JXF_Int           n       = fsls_VecGroupN(vecgroup);
   fsls_Vector **x_array = fsls_VecGroupXArray(vecgroup);
   
   fsls_Vector  *x  = NULL;
   fsls_Vector  *xg = NULL;
   
   JXF_Real *x_data  = NULL;
   JXF_Real *xg_data = NULL;
   
   JXF_Int i,j;
   JXF_Int cnt;
   JXF_Int sizeg,size = 0;
   
   for (i = 0; i < n; i ++)
   {
      size += fsls_VectorSize(x_array[i]);
   }

   x = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(x);
   x_data = fsls_VectorData(x);
   
   cnt = 0;
   for (i = 0; i < n; i ++)
   {
      xg = x_array[i];
      xg_data = fsls_VectorData(xg);
      sizeg   = fsls_VectorSize(xg);
      for (j = 0; j < sizeg; j ++)
      {
         x_data[cnt++] = xg_data[j]; 
      }
   }

   return x;
}

/*!
 * \fn JXF_Int fsls_SeqVectorPartitionByRow
 * \brief Split a vector into 'np' sub vectors by row-partition.
 * \param *x pointer to the vector to be splitted
 * \param np number of sub-vectors after splitting
 * \param *partition the row-partition(size is 'np+1')
 * \param ***x_array_ptr pointer to all the sub-vectors
 * \author peghoty
 * \date 2010/11/22
 */
JXF_Int 
fsls_SeqVectorPartitionByRow( fsls_Vector *x, JXF_Int np, JXF_Int *partition, fsls_Vector ***x_array_ptr )
{
   JXF_Real *x_data = fsls_VectorData(x);
   JXF_Int     i, begin, local_size;
   
   fsls_Vector **x_array = fsls_CTAlloc(fsls_Vector *, np);
     
   for (i = 0; i < np; i ++)
   {
      begin = partition[i];
      local_size = partition[i+1] - begin;
      x_array[i] = fsls_SeqVectorCreate(local_size);
      fsls_SeqVectorInitialize(x_array[i]);
      memcpy(fsls_VectorData(x_array[i]), &x_data[begin], local_size*sizeof(JXF_Real) );
   }

   *x_array_ptr = x_array;
   
   return 0;
}

/*!
 * \fn JXF_Int fsls_SeqVectorPermute
 * \brief Permute a vector according to the ordering 'p', i.e., y = P*x.
 * \param *x pointer to the vector to be permuted
 * \param *p pointer to the ordering
 * \author peghoty
 * \date 2010/12/12
 */
fsls_Vector * 
fsls_SeqVectorPermute( fsls_Vector *x, JXF_Int *p )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
   
   fsls_Vector *y      = NULL;
   JXF_Real      *y_data = NULL;
   
   JXF_Int i;
   
   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);
   y_data = fsls_VectorData(y);
   
   for (i = 0; i < size; i ++)
   {
      y_data[i] = x_data[p[i]];
   }
   
   return y;
}

/*!
 * \fn JXF_Int fsls_SeqVectorPermuteT
 * \brief Permute a vector according to the ordering 'p', i.e., y = P^T*x.
 * \param *x pointer to the vector to be permuted
 * \param *p pointer to the ordering
 * \author peghoty
 * \date 2010/12/13
 */
fsls_Vector * 
fsls_SeqVectorPermuteT( fsls_Vector *x, JXF_Int *p )
{
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
   
   fsls_Vector *y      = NULL;
   JXF_Real      *y_data = NULL;
   
   JXF_Int i;
   
   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);
   y_data = fsls_VectorData(y);
   
   for (i = 0; i < size; i ++)
   {
      y_data[p[i]] = x_data[i];
   }
   
   return y;
}

/*!
 * \fn JXF_Int fsls_SeqVectorScaleBack
 * \brief Diag-Scaling back of a vector.
 * \date 2011/07/03
 * \author peghoty
 */
JXF_Int
fsls_SeqVectorScaleBack( fsls_Vector *f, JXF_Real *v )
{
   JXF_Int     size   = fsls_VectorSize(f);
   JXF_Real *f_data = fsls_VectorData(f); 
                             
   JXF_Int i;
      
   // f = V*f
   for (i = 0; i < size; i ++)
   {
      f_data[i] = f_data[i]*v[i];
   }
   
   return (0);
}

/*!
 * \fn fsls_Vector *fsls_VecfromREItoEIR
 * \brief Transfer a R-E-I ordered vector into an E-I-R ordered vector.
 *    /  \      /  \  
 *   | Xr |    | Xe |  
 *   | Xe | => | Xi |
 *   | Xi |    | Xr |
 *    \  /      \  / 
 * \author peghoty
 * \date 2011/06/27 
 */
fsls_Vector *
fsls_VecfromREItoEIR( fsls_Vector *x )
{
   fsls_Vector *y = NULL;
   
   JXF_Int     size   = fsls_VectorSize(x);
   JXF_Real *x_data = fsls_VectorData(x);
   JXF_Real *y_data = NULL;
   
   JXF_Int i;
   JXF_Int n = size / 3;
   JXF_Int nplusn = 2*n;
   
   y = fsls_SeqVectorCreate(size);
   fsls_SeqVectorInitialize(y);
   y_data = fsls_VectorData(y);
   
   for (i = 0; i < n; i ++)
   {
      y_data[i] = x_data[i+n];
   }
   for (i = n; i < nplusn; i ++)
   {
      y_data[i] = x_data[i+n];
   }
   for (i = nplusn; i < size; i ++)
   {
      y_data[i] = x_data[i-nplusn];
   }
   
   return (y);
} 
