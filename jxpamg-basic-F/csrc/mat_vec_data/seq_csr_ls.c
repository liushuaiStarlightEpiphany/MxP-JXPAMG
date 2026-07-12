//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_mv.h"

/*!
 * \fn void jxf_ReadLinearSystemData
 * \brief Read the data of a linear system from given two files.
 * \date 2011/09/08
 */ 
void
jxf_ReadLinearSystemData( char    *MatFile, 
                         char    *RhsFile,
                         JXF_Int    **IA_ptr,
                         JXF_Int    **JA_ptr,
                         JXF_Real **AA_ptr,
                         JXF_Real **F_ptr,
                         JXF_Int     *n_ptr   )
{
   FILE   *fp = NULL;
   JXF_Int     n, nz, j;

   JXF_Int    *IA = NULL;
   JXF_Int    *JA = NULL;
   JXF_Real *AA = NULL;
   JXF_Real *F  = NULL;
                            
   /* Read the matrix from file */
   
   fp = fopen(MatFile, "r");
   
   jxf_fscanf(fp, "%d", &n);
   
   IA = jxf_CTAlloc(JXF_Int, n + 1);
   
   for (j = 0; j < n + 1; j ++)
   {
      jxf_fscanf(fp, "%d", &IA[j]);
      IA[j] -= 1; // index-displacing
   }
   
   nz = IA[n];
   JA = jxf_CTAlloc(JXF_Int, nz);
   
   for (j = 0; j < nz; j ++)
   {
      jxf_fscanf(fp, "%d", &JA[j]);
      JA[j] -= 1; // index-displacing
   }
   
   AA = jxf_CTAlloc(JXF_Real, nz);
   
   for (j = 0; j < nz; j ++)
   {
      jxf_fscanf(fp, "%le", &AA[j]);
   }
   fclose(fp);
      
   /* Read the right hand side vector from file */
   
   fp = fopen(RhsFile, "r");
   
   jxf_fscanf(fp, "%d", &n); 
   
   F = jxf_CTAlloc(JXF_Real, n);
   
   for (j = 0; j < n; j ++)
   {
      jxf_fscanf(fp, "%le", &F[j]);
   }
   fclose(fp);
   
   *IA_ptr = IA;
   *JA_ptr = JA;
   *AA_ptr = AA;
   *F_ptr  = F;
   *n_ptr  = n;
}      

/*!
 * \fn void jxf_InitializeLinearSystemStruct
 * \brief Initialize the structs of a sequential matrix and a vector for a linear system.
 * \date 2011/09/08
 */ 
void
jxf_InitializeLinearSystemStruct( jxf_CSRMatrix **matrix_ptr, 
                                 jxf_Vector    **vectorF_ptr, 
                                 jxf_Vector    **vectorX_ptr )
{
   jxf_CSRMatrix *matrix  = NULL; 
   jxf_Vector    *vectorF = NULL; 
   jxf_Vector    *vectorX = NULL;
                                 
   matrix  = jxf_CTAlloc(jxf_CSRMatrix, 1);
   jxf_CSRMatrixData(matrix)     = NULL;
   jxf_CSRMatrixI(matrix)        = NULL;
   jxf_CSRMatrixJ(matrix)        = NULL;
   jxf_CSRMatrixRownnz(matrix)   = NULL;
   jxf_CSRMatrixOwnsData(matrix) = 1;
   
   vectorF = jxf_CTAlloc(jxf_Vector, 1);
   jxf_VectorData(vectorF)                  = NULL;
   jxf_VectorNumVectors(vectorF)            = 1;
   jxf_VectorMultiVecStorageMethod(vectorF) = 0;
   jxf_VectorOwnsData(vectorF)              = 1;   

   vectorX = jxf_CTAlloc(jxf_Vector, 1);
   jxf_VectorData(vectorX)                  = NULL;
   jxf_VectorNumVectors(vectorX)            = 1;
   jxf_VectorMultiVecStorageMethod(vectorX) = 0;
   jxf_VectorOwnsData(vectorX)              = 1; 
   
   *matrix_ptr  = matrix;
   *vectorF_ptr = vectorF;  
   *vectorX_ptr = vectorX; 
}

/*!
 * \fn void jxf_FillData2Struct
 * \brief Arrays -> structs for a linear system.
 * \note The zeros entries are all deleted during the scanning. 
 * \date 2011/09/08
 */ 
void
jxf_FillData2Struct( JXF_Int              n, 
                    JXF_Int             *ia, 
                    JXF_Int             *ja, 
                    JXF_Real          *aa, 
                    JXF_Real          *f, 
                    JXF_Real          *x, 
                    jxf_CSRMatrix    *matrix, 
                    jxf_Vector       *vectorF, 
                    jxf_Vector       *vectorX  )
{
   JXF_Int nz;
   JXF_Int i, j;
   JXF_Int start;
  
   /* get rid of the zeros in the matrix */ 
   nz = 0;
   for (i = 0; i < n; i ++)
   {
      start = ia[i];
      ia[i] = nz;
      for (j = start; j < ia[i+1]; j ++)
      {
         if (aa[j])
         {
            ja[nz] = ja[j];
            aa[nz] = aa[j];
            nz ++;
         }
      }
   } 
   ia[n] = nz;

   jxf_CSRMatrixNumRows(matrix)     = n;  
   jxf_CSRMatrixNumCols(matrix)     = n;
   jxf_CSRMatrixNumNonzeros(matrix) = nz;
   jxf_CSRMatrixData(matrix)        = aa;
   jxf_CSRMatrixI(matrix)           = ia;
   jxf_CSRMatrixJ(matrix)           = ja;
   jxf_CSRMatrixOwnsData(matrix)    = 1;
   jxf_CSRMatrixRownnz(matrix)      = NULL; 
   jxf_CSRMatrixNumRownnz(matrix)   = n; 
   
   jxf_CSRMatrixReorder(matrix); 
   
   /* for the vectors */
   jxf_VectorData(vectorF) = f;
   jxf_VectorSize(vectorF) = n;
   jxf_VectorNumVectors(vectorF)            = 1;
   jxf_VectorMultiVecStorageMethod(vectorF) = 0;
   jxf_VectorOwnsData(vectorF)              = 1; 
   
   jxf_VectorData(vectorX) = x;
   jxf_VectorSize(vectorX) = n;
   jxf_VectorNumVectors(vectorX)            = 1;
   jxf_VectorMultiVecStorageMethod(vectorX) = 0;
   jxf_VectorOwnsData(vectorX)              = 1; 
} 
