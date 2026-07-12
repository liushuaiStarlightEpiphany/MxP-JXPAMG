//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  seq_vector.c -- basic operations for sequential vectors
 *  Date: 2011/05/13
 */

#include "jxf_mv.h"

/*!
 * \fn jxf_Vector *jxf_SeqVectorCreate
 * \brief Create a sequential vector with size of 'size'.
 * \note No memories for the vector are allocated here.
 * \date 2011/05/13
 */
jxf_Vector *
jxf_SeqVectorCreate(JXF_Int size)
{
   jxf_Vector *vector;

   vector = jxf_CTAlloc(jxf_Vector, 1);

   jxf_VectorData(vector) = NULL;
   jxf_VectorSize(vector) = size;

   jxf_VectorNumVectors(vector) = 1;
   jxf_VectorMultiVecStorageMethod(vector) = 0;

   /* set defaults */
   jxf_VectorOwnsData(vector) = 1;

   return vector;
}

/*!
 * \fn JXF_Int jxf_SeqVectorInitialize
 * \brief Initialize a sequential vector.
 * \note Allocate memory for the vector if needed.
 * \date 2011/05/13
 */
JXF_Int
jxf_SeqVectorInitialize(jxf_Vector *vector)
{
   JXF_Int size = jxf_VectorSize(vector);
   JXF_Int num_vectors = jxf_VectorNumVectors(vector);
   JXF_Int multivec_storage_method = jxf_VectorMultiVecStorageMethod(vector);
   JXF_Int ierr = 0;

   if (!jxf_VectorData(vector))
   {
      jxf_VectorData(vector) = jxf_CTAlloc(JXF_Real, num_vectors * size);
   }

   if (multivec_storage_method == 0)
   {
      jxf_VectorVectorStride(vector) = size;
      jxf_VectorIndexStride(vector) = 1;
   }
   else if (multivec_storage_method == 1)
   {
      jxf_VectorVectorStride(vector) = 1;
      jxf_VectorIndexStride(vector) = num_vectors;
   }
   else
   {
      ++ierr;
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxf_SeqVectorSetDataOwner
 * \brief Set the member data_owner of a sequential vector.
 * \date 2011/09/07
 */
JXF_Int
jxf_SeqVectorSetDataOwner(jxf_Vector *vector,
                         JXF_Int owns_data)
{
   JXF_Int ierr = 0;

   jxf_VectorOwnsData(vector) = owns_data;

   return ierr;
}

/*!
 * \fn JXF_Int jxf_SeqVectorDestroy
 * \brief Destroy a sequential vector.
 * \date 2011/05/13
 */
JXF_Int
jxf_SeqVectorDestroy(jxf_Vector *vector)
{
   JXF_Int ierr = 0;

   if (vector)
   {
      if (jxf_VectorOwnsData(vector))
      {
         jxf_TFree(jxf_VectorData(vector));
      }
      jxf_TFree(vector);
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxf_SeqVectorPrint
 * \brief Print a sequential vector into a given file.
 * \date 2011/05/13
 */
JXF_Int
jxf_SeqVectorPrint(jxf_Vector *vector, char *file_name)
{
   FILE *fp = NULL;

   JXF_Real *data;
   JXF_Int size, num_vectors, vecstride, idxstride;

   JXF_Int i, j, ierr = 0;

   num_vectors = jxf_VectorNumVectors(vector);
   vecstride = jxf_VectorVectorStride(vector);
   idxstride = jxf_VectorIndexStride(vector);
   data = jxf_VectorData(vector);
   size = jxf_VectorSize(vector);

   fp = fopen(file_name, "w");

   if (num_vectors == 1)
   {
      jxf_fprintf(fp, "%d\n", size);
   }
   else
   {
      jxf_fprintf(fp, "%d vectors of size %d\n", num_vectors, size);
   }

   if (num_vectors > 1)
   {
      for (j = 0; j < num_vectors; j++)
      {
         jxf_fprintf(fp, "vector %d\n", j);
         for (i = 0; i < size; i++)
         {
            jxf_fprintf(fp, "%.15e\n", data[j * vecstride + i * idxstride]);
         }
      }
   }
   else // if (num_vectors == 1)
   {
      for (i = 0; i < size; i++)
      {
         jxf_fprintf(fp, "%.15e\n", data[i]);
      }
   }

   fclose(fp);

   return ierr;
}

/*!
 * \fn jxf_Vector *jxf_SeqVectorRead
 * \brief Read a sequential vector from a given file.
 * \date 2011/05/13
 */
jxf_Vector *
jxf_SeqVectorRead(char *file_name)
{
   jxf_Vector *vector;

   FILE *fp;

   JXF_Real *data;
   JXF_Int size;

   JXF_Int j;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &size);

   vector = jxf_SeqVectorCreate(size);
   jxf_SeqVectorInitialize(vector);

   data = jxf_VectorData(vector);
   for (j = 0; j < size; j++)
   {
      jxf_fscanf(fp, "%le", &data[j]);
   }

   fclose(fp);

   /* multivector code not written yet */
   jxf_assert(jxf_VectorNumVectors(vector) == 1);

   return vector;
}

// jxf_Vector *
// jxf_SeqVectorBinaryRead( char *file_name )
// {
//    jxf_Vector *vector;

//    FILE      *fp;

//    JXF_Real   *data;
//    JXF_Int     size;

//    JXF_Int j;

//    fp = fopen(file_name, "rb");

//    fread(&size, sizeof(JXF_Int), 1, fp);

//    vector = jxf_SeqVectorCreate(size);
//    jxf_SeqVectorInitialize(vector);

//    data = jxf_VectorData(vector);
//    for (j = 0; j < size; j ++)
//    {
//       fread(&data[j], sizeof(JXF_Real), 1, fp);
//    }

//    fclose(fp);

//    /* multivector code not written yet */
//    jxf_assert( jxf_VectorNumVectors(vector) == 1);

//    return vector;
// }

jxf_Vector *
jxf_SeqVectorBinaryRead(char *file_name)
{
   jxf_Vector *vector;
   FILE *fp;
   JXF_Real *data;
   JXF_Int size;

   fp = fopen(file_name, "rb");
   if (!fp)
      return NULL;

   // 1. 读取大小
   fread(&size, sizeof(JXF_Int), 1, fp);

   // 2. 创建向量并初始化
   vector = jxf_SeqVectorCreate(size);
   jxf_SeqVectorInitialize(vector);

   // 3. 块读取数据
   data = jxf_VectorData(vector);
   fread(data, sizeof(JXF_Real), size, fp);

   fclose(fp);

   /* multivector code not written yet */
   jxf_assert(jxf_VectorNumVectors(vector) == 1);

   printf("Vector loaded successfully.\n\n");

   return vector;
}

jxf_Vector *
jxf_SeqVectorRead2(char *file_name)
{
   jxf_Vector *vector;

   FILE *fp;

   JXF_Real *data;
   JXF_Int size;

   JXF_Int j, tmp;

   fp = fopen(file_name, "r");

   jxf_fscanf(fp, "%d", &size);

   vector = jxf_SeqVectorCreate(size);
   jxf_SeqVectorInitialize(vector);

   data = jxf_VectorData(vector);
   for (j = 0; j < size; j++)
   {
      jxf_fscanf(fp, "%d%le", &tmp, &data[j]);
   }

   fclose(fp);

   /* multivector code not written yet */
   jxf_assert(jxf_VectorNumVectors(vector) == 1);

   return vector;
}

/*!
 * \fn jxf_Vector *jxf_SeqMultiVectorCreate
 * \brief Create a multi vector.
 * \date 2011/05/13
 */
jxf_Vector *
jxf_SeqMultiVectorCreate(JXF_Int size, JXF_Int num_vectors)
{
   jxf_Vector *vector = jxf_SeqVectorCreate(size);
   jxf_VectorNumVectors(vector) = num_vectors;
   return vector;
}
