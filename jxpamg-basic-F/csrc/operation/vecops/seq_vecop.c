//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  seq_vecop.c -- basic operations for sequential vectors
 *  Date: 2011/05/13
 */ 

#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_SeqVectorCopy
 * \brief Copy data from x to y.
 * \note y should have already been initialized with the same size as x.
 * \date 2011/05/13
 */ 
JXF_Int
jxf_SeqVectorCopy( jxf_Vector *x, jxf_Vector *y )
{
   JXF_Real  *x_data = jxf_VectorData(x);
   JXF_Real  *y_data = jxf_VectorData(y);
   JXF_Int      size   = jxf_VectorSize(x);
           
   JXF_Int      i, ierr = 0;

   size *= jxf_VectorNumVectors(x);
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      y_data[i] = x_data[i];
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxf_SeqVectorSetConstantValues
 * \brief Set all components of the vector to be a given constant.
 * \date 2011/05/13
 */  
JXF_Int
jxf_SeqVectorSetConstantValues( jxf_Vector *x, JXF_Real value )
{
   JXF_Real  *vector_data = jxf_VectorData(x);
   JXF_Int      size        = jxf_VectorSize(x);
           
   JXF_Int      i, ierr  = 0;

   size *= jxf_VectorNumVectors(x);

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      vector_data[i] = value;
   }

   return ierr;
}

/*!
 * \fn JXF_Real jxf_SeqVectorNorm1
 * \brief Compute the 1-norm of vector x.
 *        norm1 = max{ |x_i|, i = 0,1,...,n-1 }
 * \date 2011/05/13
 */ 
JXF_Real 
jxf_SeqVectorNorm1( jxf_Vector *x )
{
   JXF_Real  *x_data = jxf_VectorData(x);
   JXF_Real   norm1  = 0.0;
   JXF_Real   tmp    = 0.0;
   JXF_Int      size   = jxf_VectorSize(x);
   JXF_Int      i;
   
   norm1 = fabs(x_data[0]);
   for (i = 1; i < size; i ++)
   {
      tmp = fabs(x_data[i]);
      if( tmp > norm1 ) 
      {
         norm1 = tmp;
      }
   }

   return norm1;
}

/*!
 * \fn JXF_Real jxf_SeqVectorPointWiseRelNorm1
 * \brief subroutine for Relative Point-wise 1-norm.
 *        result = max{ |x_i| / |y_i|, i = 0,1,...,n-1 }
 * \date 2011/05/13
 */ 
JXF_Real 
jxf_SeqVectorPointWiseRelNorm1( jxf_Vector *x, jxf_Vector *y )
{
   JXF_Real  *x_data = jxf_VectorData(x);
   JXF_Real  *y_data = jxf_VectorData(y);
   JXF_Int      size   = jxf_VectorSize(x);
           
   JXF_Int      i;
   JXF_Real   tmp    = 0.0;
   JXF_Real   result = 0.0;

   for (i = 0; i < size; i ++)
   {
      if (y_data[i] < EPS)
      {
         tmp = fabs(x_data[i]);
      }
      else
      {
         tmp = fabs(x_data[i] / y_data[i]);
      }
      
      if (tmp > result) result = tmp;
   }

   return result;
}

/*!
 * \fn JXF_Int jxf_SeqVectorScale
 * \brief Scale a vector, i.e., y := alpha*y.
 * \date 2011/05/13
 */
JXF_Int
jxf_SeqVectorScale( JXF_Real alpha, jxf_Vector *y )
{
   JXF_Real  *y_data = jxf_VectorData(y);
   JXF_Int      size   = jxf_VectorSize(y);    
   JXF_Int      i, ierr = 0;

   size *= jxf_VectorNumVectors(y);

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      y_data[i] *= alpha;
   }

   return ierr;
}

/*!
 * \fn JXF_Int jxf_SeqVecMul
 * \brief Performs z = x*y, where z_i is defiend by x_i*y_i.
 * \note x, y, z should be allocated and of the same size.
 * \author peghoty
 * \date 2011/09/09 
 */ 
JXF_Int
jxf_SeqVecMul( jxf_Vector *x,
              jxf_Vector *y,
              jxf_Vector *z )
{
   JXF_Int i, size = jxf_VectorSize(x);
   
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Real *z_data = jxf_VectorData(z);
   
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      z_data[i] = x_data[i] * y_data[i];
   }
   
   return 0;
}

/*!
 * \fn JXF_Int jxf_SeqVecZXY
 * \brief Performs z = z + alpha*x*y, where xy_i is defiend by x_i*y_i.
 * \note x, y, z should be allocated and be of the same size.
 * \author peghoty
 * \date 2011/09/09 
 */ 
JXF_Int
jxf_SeqVecZXY( jxf_Vector  *z,
              JXF_Real      alpha,
              jxf_Vector  *x,
              jxf_Vector  *y )
{
   JXF_Int i, size = jxf_VectorSize(x);

   JXF_Real *z_data = jxf_VectorData(z);
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   
   if (alpha == 0.0)
   {
      return 0; 
   }
   else if (alpha == 1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i ++)
      {
         z_data[i] += x_data[i] * y_data[i];
      }
   }
   else if (alpha == -1.0)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i ++)
      {
         z_data[i] -= x_data[i] * y_data[i];
      }
   }
   else
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i ++)
      {
         z_data[i] += alpha * x_data[i] * y_data[i];
      }
   }

   return 0;   
}

/*!
 * \fn JXF_Real jxf_SeqVectorInnerProd
 * \brief Compute the inner product of vectors x and y.
 * \date 2011/05/13
 */   
JXF_Real 
jxf_SeqVectorInnerProd( jxf_Vector *x, jxf_Vector *y )
{
   JXF_Real  *x_data = jxf_VectorData(x);
   JXF_Real  *y_data = jxf_VectorData(y);
   JXF_Int      size   = jxf_VectorSize(x);
           
   JXF_Int      i;
   JXF_Real   result = 0.0;

   size *= jxf_VectorNumVectors(x);

#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS result
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      result += y_data[i] * x_data[i];
   }

   return result;
}

/*!
 * \fn JXF_Int jxf_SeqVectorAxpy
 * \brief Compute y: = alpha*x + y.
 * \date 2011/05/13
 */ 
JXF_Int
jxf_SeqVectorAxpy( JXF_Real alpha, jxf_Vector *x, jxf_Vector *y )
{
   JXF_Real  *x_data = jxf_VectorData(x);
   JXF_Real  *y_data = jxf_VectorData(y);
   JXF_Int      size   = jxf_VectorSize(x);
           
   JXF_Int      i, ierr = 0;

   size *= jxf_VectorNumVectors(x);

#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
   for (i = 0; i < size; i ++)
   {
      y_data[i] += alpha * x_data[i];
   }

   return ierr;
}

JXF_Int
jxf_SeqVectorMassDotpTwo( jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Int unroll, JXF_Real *result_x, JXF_Real *result_y )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Real *z_data = jxf_VectorData(z[0]);
   JXF_Real  res_x, res_y;
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int      i, j, jstart;

   if (unroll == 8)
   {
      jxf_SeqVectorMassDotpTwo8(x, y, z, k, result_x, result_y);
      return jxf_error_flag;
   }
   else if (unroll == 4)
   {
      jxf_SeqVectorMassDotpTwo4(x, y, z, k, result_x, result_y);
      return jxf_error_flag;
   }
   else
   {
      for (j = 0; j < k; j++)
      {
         res_x = result_x[j];
         res_y = result_y[j];
         jstart = j*size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x,res_y
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            res_x += z_data[jstart+i] * x_data[i];
            res_y += z_data[jstart+i] * y_data[i];
         }
         result_x[j] = res_x;
         result_y[j] = res_y;
      }
   }
   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassDotpTwo8( jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Real *result_x, JXF_Real *result_y )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Real *z_data = jxf_VectorData(z[0]);
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int   i, j, restk;
   JXF_Real  res_x1, res_x2, res_x3, res_x4, res_x5, res_x6, res_x7, res_x8;
   JXF_Real  res_y1, res_y2, res_y3, res_y4, res_y5, res_y6, res_y7, res_y8;
   JXF_Int   jstart, jstart1, jstart2, jstart3, jstart4, jstart5, jstart6, jstart7;

   restk = (k-(k/8*8));
   if (k > 7)
   {
      for (j = 0; j < k-7; j += 8)
      {
         res_x1 = 0;
         res_x2 = 0;
         res_x3 = 0;
         res_x4 = 0;
         res_x5 = 0;
         res_x6 = 0;
         res_x7 = 0;
         res_x8 = 0;
         res_y1 = 0;
         res_y2 = 0;
         res_y3 = 0;
         res_y4 = 0;
         res_y5 = 0;
         res_y6 = 0;
         res_y7 = 0;
         res_y8 = 0;
         jstart = j*size;
         jstart1 = jstart+size;
         jstart2 = jstart1+size;
         jstart3 = jstart2+size;
         jstart4 = jstart3+size;
         jstart5 = jstart4+size;
         jstart6 = jstart5+size;
         jstart7 = jstart6+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_x5,res_x6,res_x7,res_x8,res_y1,res_y2,res_y3,res_y4,res_y5,res_y6,res_y7,res_y8
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            res_x1 += z_data[jstart+i] * x_data[i];
            res_y1 += z_data[jstart+i] * y_data[i];
            res_x2 += z_data[jstart1+i] * x_data[i];
            res_y2 += z_data[jstart1+i] * y_data[i];
            res_x3 += z_data[jstart2+i] * x_data[i];
            res_y3 += z_data[jstart2+i] * y_data[i];
            res_x4 += z_data[jstart3+i] * x_data[i];
            res_y4 += z_data[jstart3+i] * y_data[i];
            res_x5 += z_data[jstart4+i] * x_data[i];
            res_y5 += z_data[jstart4+i] * y_data[i];
            res_x6 += z_data[jstart5+i] * x_data[i];
            res_y6 += z_data[jstart5+i] * y_data[i];
            res_x7 += z_data[jstart6+i] * x_data[i];
            res_y7 += z_data[jstart6+i] * y_data[i];
            res_x8 += z_data[jstart7+i] * x_data[i];
            res_y8 += z_data[jstart7+i] * y_data[i];
         }
         result_x[j] = res_x1;
         result_x[j+1] = res_x2;
         result_x[j+2] = res_x3;
         result_x[j+3] = res_x4;
         result_x[j+4] = res_x5;
         result_x[j+5] = res_x6;
         result_x[j+6] = res_x7;
         result_x[j+7] = res_x8;
         result_y[j] = res_y1;
         result_y[j+1] = res_y2;
         result_y[j+2] = res_y3;
         result_y[j+3] = res_y4;
         result_y[j+4] = res_y5;
         result_y[j+5] = res_y6;
         result_y[j+6] = res_y7;
         result_y[j+7] = res_y8;
      }
   }
   if (restk == 1)
   {
      res_x1 = 0;
      res_y1 = 0;
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_y1
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
      }
      result_x[k-1] = res_x1;
      result_y[k-1] = res_y1;
   }
   else if (restk == 2)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_y1 = 0;
      res_y2 = 0;
      jstart = (k-2)*size;
      jstart1 = jstart+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_y1,res_y2
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
      }
      result_x[k-2] = res_x1;
      result_x[k-1] = res_x2;
      result_y[k-2] = res_y1;
      result_y[k-1] = res_y2;
   }
   else if (restk == 3)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      jstart = (k-3)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_y1,res_y2,res_y3
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
      }
      result_x[k-3] = res_x1;
      result_x[k-2] = res_x2;
      result_x[k-1] = res_x3;
      result_y[k-3] = res_y1;
      result_y[k-2] = res_y2;
      result_y[k-1] = res_y3;
   }
   else if (restk == 4)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_x4 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      res_y4 = 0;
      jstart = (k-4)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_y1,res_y2,res_y3,res_y4
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
         res_x4 += z_data[jstart3+i] * x_data[i];
         res_y4 += z_data[jstart3+i] * y_data[i];
      }
      result_x[k-4] = res_x1;
      result_x[k-3] = res_x2;
      result_x[k-2] = res_x3;
      result_x[k-1] = res_x4;
      result_y[k-4] = res_y1;
      result_y[k-3] = res_y2;
      result_y[k-2] = res_y3;
      result_y[k-1] = res_y4;
   }
   else if (restk == 5)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_x4 = 0;
      res_x5 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      res_y4 = 0;
      res_y5 = 0;
      jstart = (k-5)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_x5,res_y1,res_y2,res_y3,res_y4,res_y5
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
         res_x4 += z_data[jstart3+i] * x_data[i];
         res_y4 += z_data[jstart3+i] * y_data[i];
         res_x5 += z_data[jstart4+i] * x_data[i];
         res_y5 += z_data[jstart4+i] * y_data[i];
      }
      result_x[k-5] = res_x1;
      result_x[k-4] = res_x2;
      result_x[k-3] = res_x3;
      result_x[k-2] = res_x4;
      result_x[k-1] = res_x5;
      result_y[k-5] = res_y1;
      result_y[k-4] = res_y2;
      result_y[k-3] = res_y3;
      result_y[k-2] = res_y4;
      result_y[k-1] = res_y5;
   }
   else if (restk == 6)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_x4 = 0;
      res_x5 = 0;
      res_x6 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      res_y4 = 0;
      res_y5 = 0;
      res_y6 = 0;
      jstart = (k-6)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
      jstart5 = jstart4+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_x5,res_x6,res_y1,res_y2,res_y3,res_y4,res_y5,res_y6
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
         res_x4 += z_data[jstart3+i] * x_data[i];
         res_y4 += z_data[jstart3+i] * y_data[i];
         res_x5 += z_data[jstart4+i] * x_data[i];
         res_y5 += z_data[jstart4+i] * y_data[i];
         res_x6 += z_data[jstart5+i] * x_data[i];
         res_y6 += z_data[jstart5+i] * y_data[i];
      }
      result_x[k-6] = res_x1;
      result_x[k-5] = res_x2;
      result_x[k-4] = res_x3;
      result_x[k-3] = res_x4;
      result_x[k-2] = res_x5;
      result_x[k-1] = res_x6;
      result_y[k-6] = res_y1;
      result_y[k-5] = res_y2;
      result_y[k-4] = res_y3;
      result_y[k-3] = res_y4;
      result_y[k-2] = res_y5;
      result_y[k-1] = res_y6;
   }
   else if (restk == 7)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_x4 = 0;
      res_x5 = 0;
      res_x6 = 0;
      res_x7 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      res_y4 = 0;
      res_y5 = 0;
      res_y6 = 0;
      res_y7 = 0;
      jstart = (k-7)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
      jstart5 = jstart4+size;
      jstart6 = jstart5+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_x5,res_x6,res_x7,res_y1,res_y2,res_y3,res_y4,res_y5,res_y6,res_y7
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
         res_x4 += z_data[jstart3+i] * x_data[i];
         res_y4 += z_data[jstart3+i] * y_data[i];
         res_x5 += z_data[jstart4+i] * x_data[i];
         res_y5 += z_data[jstart4+i] * y_data[i];
         res_x6 += z_data[jstart5+i] * x_data[i];
         res_y6 += z_data[jstart5+i] * y_data[i];
         res_x7 += z_data[jstart6+i] * x_data[i];
         res_y7 += z_data[jstart6+i] * y_data[i];
      }
      result_x[k-7] = res_x1;
      result_x[k-6] = res_x2;
      result_x[k-5] = res_x3;
      result_x[k-4] = res_x4;
      result_x[k-3] = res_x5;
      result_x[k-2] = res_x6;
      result_x[k-1] = res_x7;
      result_y[k-7] = res_y1;
      result_y[k-6] = res_y2;
      result_y[k-5] = res_y3;
      result_y[k-4] = res_y4;
      result_y[k-3] = res_y5;
      result_y[k-2] = res_y6;
      result_y[k-1] = res_y7;
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassDotpTwo4( jxf_Vector *x, jxf_Vector *y, jxf_Vector **z, JXF_Int k, JXF_Real *result_x, JXF_Real *result_y )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Real *z_data = jxf_VectorData(z[0]);
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int   i, j, restk;
   JXF_Real  res_x1, res_x2, res_x3, res_x4, res_y1, res_y2, res_y3, res_y4;
   JXF_Int   jstart, jstart1, jstart2, jstart3;

   restk = (k-(k/4*4));
   if (k > 3)
   {
      for (j = 0; j < k-3; j += 4)
      {
         res_x1 = 0;
         res_x2 = 0;
         res_x3 = 0;
         res_x4 = 0;
         res_y1 = 0;
         res_y2 = 0;
         res_y3 = 0;
         res_y4 = 0;
         jstart = j*size;
         jstart1 = jstart+size;
         jstart2 = jstart1+size;
         jstart3 = jstart2+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_x4,res_y1,res_y2,res_y3,res_y4
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            res_x1 += z_data[jstart+i] * x_data[i];
            res_y1 += z_data[jstart+i] * y_data[i];
            res_x2 += z_data[jstart1+i] * x_data[i];
            res_y2 += z_data[jstart1+i] * y_data[i];
            res_x3 += z_data[jstart2+i] * x_data[i];
            res_y3 += z_data[jstart2+i] * y_data[i];
            res_x4 += z_data[jstart3+i] * x_data[i];
            res_y4 += z_data[jstart3+i] * y_data[i];
         }
         result_x[j] = res_x1;
         result_x[j+1] = res_x2;
         result_x[j+2] = res_x3;
         result_x[j+3] = res_x4;
         result_y[j] = res_y1;
         result_y[j+1] = res_y2;
         result_y[j+2] = res_y3;
         result_y[j+3] = res_y4;
      }
   }
   if (restk == 1)
   {
      res_x1 = 0;
      res_y1 = 0;
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_y1
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
      }
      result_x[k-1] = res_x1;
      result_y[k-1] = res_y1;
   }
   else if (restk == 2)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_y1 = 0;
      res_y2 = 0;
      jstart = (k-2)*size;
      jstart1 = jstart+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_y1,res_y2
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
      }
      result_x[k-2] = res_x1;
      result_x[k-1] = res_x2;
      result_y[k-2] = res_y1;
      result_y[k-1] = res_y2;
   }
   else if (restk == 3)
   {
      res_x1 = 0;
      res_x2 = 0;
      res_x3 = 0;
      res_y1 = 0;
      res_y2 = 0;
      res_y3 = 0;
      jstart = (k-3)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res_x1,res_x2,res_x3,res_y1,res_y2,res_y3
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res_x1 += z_data[jstart+i] * x_data[i];
         res_y1 += z_data[jstart+i] * y_data[i];
         res_x2 += z_data[jstart1+i] * x_data[i];
         res_y2 += z_data[jstart1+i] * y_data[i];
         res_x3 += z_data[jstart2+i] * x_data[i];
         res_y3 += z_data[jstart2+i] * y_data[i];
      }
      result_x[k-3] = res_x1;
      result_x[k-2] = res_x2;
      result_x[k-1] = res_x3;
      result_y[k-3] = res_y1;
      result_y[k-2] = res_y2;
      result_y[k-1] = res_y3;
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassAxpy( JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k, JXF_Int unroll )
{
   JXF_Real *x_data = jxf_VectorData(x[0]);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int   size   = jxf_VectorSize(x[0]);
   JXF_Int   i, j, jstart;

   if (unroll == 8)
   {
      jxf_SeqVectorMassAxpy8(alpha, x, y, k);
      return jxf_error_flag;
   }
   else if (unroll == 4)
   {
      jxf_SeqVectorMassAxpy4(alpha, x, y, k);
      return jxf_error_flag;
   }
   else
   {
      for (j = 0; j < k; j++)
      {
         jstart = j*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            y_data[i] += alpha[j]*x_data[jstart+i];
         }
      }
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassAxpy8( JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k )
{
   JXF_Real *x_data = jxf_VectorData(x[0]);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int   size   = jxf_VectorSize(x[0]);
   JXF_Int   i, j, jstart, restk;

   restk = (k-(k/8*8));
   if (k > 7)
   {
      for (j = 0; j < k-7; j += 8)
      {
         jstart = j*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            y_data[i] += alpha[j]*x_data[jstart+i] + alpha[j+1]*x_data[jstart+i+size]
            + alpha[j+2]*x_data[(j+2)*size+i] + alpha[j+3]*x_data[(j+3)*size+i]
            + alpha[j+4]*x_data[(j+4)*size+i] + alpha[j+5]*x_data[(j+5)*size+i]
            + alpha[j+6]*x_data[(j+6)*size+i] + alpha[j+7]*x_data[(j+7)*size+i];
         }
      }
   }
   if (restk == 1)
   {
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-1] * x_data[jstart+i];
      }
   }
   else if (restk == 2)
   {
      jstart = (k-2)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-2] * x_data[jstart+i] + alpha[k-1] * x_data[jstart+size+i];
      }
   }
   else if (restk == 3)
   {
      jstart = (k-3)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-3] * x_data[jstart+i] + alpha[k-2] * x_data[jstart+size+i] + alpha[k-1] * x_data[(k-1)*size+i];
      }
   }
   else if (restk == 4)
   {
      jstart = (k-4)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
            y_data[i] += alpha[k-4]*x_data[(k-4)*size+i] + alpha[k-3]*x_data[(k-3)*size+i]
            + alpha[k-2]*x_data[(k-2)*size+i] + alpha[k-1]*x_data[(k-1)*size+i];
      }
   }
   else if (restk == 5)
   {
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
            y_data[i] += + alpha[k-5]*x_data[(k-5)*size+i] + alpha[k-4]*x_data[(k-4)*size+i]
            + alpha[k-3]*x_data[(k-3)*size+i] + alpha[k-2]*x_data[(k-2)*size+i]
            + alpha[k-1]*x_data[(k-1)*size+i];
      }
   }
   else if (restk == 6)
   {
      jstart = (k-6)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
            y_data[i] += alpha[k-6]*x_data[jstart+i] + alpha[k-5]*x_data[jstart+i+size]
            + alpha[k-4]*x_data[(k-4)*size+i] + alpha[k-3]*x_data[(k-3)*size+i]
            + alpha[k-2]*x_data[(k-2)*size+i] + alpha[k-1]*x_data[(k-1)*size+i];
      }
   }
   else if (restk == 7)
   {
      jstart = (k-7)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
            y_data[i] += alpha[k-7]*x_data[jstart+i] + alpha[k-6]*x_data[jstart+i+size]
            + alpha[k-5]*x_data[(k-5)*size+i] + alpha[k-4]*x_data[(k-4)*size+i]
            + alpha[k-3]*x_data[(k-3)*size+i] + alpha[k-2]*x_data[(k-2)*size+i]
            + alpha[k-1]*x_data[(k-1)*size+i];
      }
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassAxpy4( JXF_Real *alpha, jxf_Vector **x, jxf_Vector *y, JXF_Int k )
{
   JXF_Real *x_data = jxf_VectorData(x[0]);
   JXF_Real *y_data = jxf_VectorData(y);
   JXF_Int   size   = jxf_VectorSize(x[0]);
   JXF_Int   i, j, jstart, restk;

   restk = (k-(k/4*4));
   if (k > 3)
   {
      for (j = 0; j < k-3; j += 4)
      {
         jstart = j*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            y_data[i] += alpha[j]*x_data[jstart+i] + alpha[j+1]*x_data[jstart+i+size]
            + alpha[j+2]*x_data[(j+2)*size+i] + alpha[j+3]*x_data[(j+3)*size+i];
         }
      }
   }
   if (restk == 1)
   {
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-1] * x_data[jstart+i];
      }
   }
   else if (restk == 2)
   {
      jstart = (k-2)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-2] * x_data[jstart+i] + alpha[k-1] * x_data[jstart+size+i];
      }
   }
   else if (restk == 3)
   {
      jstart = (k-3)*size;
#define JXF_SMP_PRIVATE i
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         y_data[i] += alpha[k-3] * x_data[jstart+i] + alpha[k-2] * x_data[jstart+size+i] + alpha[k-1] * x_data[(k-1)*size+i];
      }
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassInnerProd( jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Int unroll, JXF_Real *result )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y[0]);
   JXF_Real  res;
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int   i, j, jstart;

   if (unroll == 8)
   {
      jxf_SeqVectorMassInnerProd8(x, y, k, result);
      return jxf_error_flag;
   }
   else if (unroll == 4)
   {
      jxf_SeqVectorMassInnerProd4(x, y, k, result);
      return jxf_error_flag;
   }
   else
   {
      for (j = 0; j < k; j ++)
      {
         res = 0;
         jstart = j * size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i ++)
         {
            res += y_data[jstart+i] * x_data[i];
         }
         result[j] = res;
      }
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassInnerProd8( jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Real *result )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y[0]);
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int   i, j, restk;
   JXF_Real  res1, res2, res3, res4, res5, res6, res7, res8;
   JXF_Int   jstart, jstart1, jstart2, jstart3, jstart4, jstart5, jstart6, jstart7;

   restk = (k-(k/8*8));
   if (k > 7)
   {
      for (j = 0; j < k-7; j += 8)
      {
         res1 = 0;
         res2 = 0;
         res3 = 0;
         res4 = 0;
         res5 = 0;
         res6 = 0;
         res7 = 0;
         res8 = 0;
         jstart = j*size;
         jstart1 = jstart+size;
         jstart2 = jstart1+size;
         jstart3 = jstart2+size;
         jstart4 = jstart3+size;
         jstart5 = jstart4+size;
         jstart6 = jstart5+size;
         jstart7 = jstart6+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4,res5,res6,res7,res8
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            res1 += y_data[jstart+i] * x_data[i];
            res2 += y_data[jstart1+i] * x_data[i];
            res3 += y_data[jstart2+i] * x_data[i];
            res4 += y_data[jstart3+i] * x_data[i];
            res5 += y_data[jstart4+i] * x_data[i];
            res6 += y_data[jstart5+i] * x_data[i];
            res7 += y_data[jstart6+i] * x_data[i];
            res8 += y_data[jstart7+i] * x_data[i];
         }
         result[j] = res1;
         result[j+1] = res2;
         result[j+2] = res3;
         result[j+3] = res4;
         result[j+4] = res5;
         result[j+5] = res6;
         result[j+6] = res7;
         result[j+7] = res8;
      }
   }
   if (restk == 1)
   {
      res1 = 0;
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
      }
      result[k-1] = res1;
   }
   else if (restk == 2)
   {
      res1 = 0;
      res2 = 0;
      jstart = (k-2)*size;
      jstart1 = jstart+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
      }
      result[k-2] = res1;
      result[k-1] = res2;
   }
   else if (restk == 3)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      jstart = (k-3)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
      }
      result[k-3] = res1;
      result[k-2] = res2;
      result[k-1] = res3;
   }
   else if (restk == 4)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      res4 = 0;
      jstart = (k-4)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
         res4 += y_data[jstart3+i] * x_data[i];
      }
      result[k-4] = res1;
      result[k-3] = res2;
      result[k-2] = res3;
      result[k-1] = res4;
   }
   else if (restk == 5)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      res4 = 0;
      res5 = 0;
      jstart = (k-5)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4,res5
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
         res4 += y_data[jstart3+i] * x_data[i];
         res5 += y_data[jstart4+i] * x_data[i];
      }
      result[k-5] = res1;
      result[k-4] = res2;
      result[k-3] = res3;
      result[k-2] = res4;
      result[k-1] = res5;
   }
   else if (restk == 6)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      res4 = 0;
      res5 = 0;
      res6 = 0;
      jstart = (k-6)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
      jstart5 = jstart4+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4,res5,res6
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
         res4 += y_data[jstart3+i] * x_data[i];
         res5 += y_data[jstart4+i] * x_data[i];
         res6 += y_data[jstart5+i] * x_data[i];
      }
      result[k-6] = res1;
      result[k-5] = res2;
      result[k-4] = res3;
      result[k-3] = res4;
      result[k-2] = res5;
      result[k-1] = res6;
   }
   else if (restk == 7)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      res4 = 0;
      res5 = 0;
      res6 = 0;
      res7 = 0;
      jstart = (k-7)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
      jstart3 = jstart2+size;
      jstart4 = jstart3+size;
      jstart5 = jstart4+size;
      jstart6 = jstart5+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4,res5,res6,res7
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
         res4 += y_data[jstart3+i] * x_data[i];
         res5 += y_data[jstart4+i] * x_data[i];
         res6 += y_data[jstart5+i] * x_data[i];
         res7 += y_data[jstart6+i] * x_data[i];
      }
      result[k-7] = res1;
      result[k-6] = res2;
      result[k-5] = res3;
      result[k-4] = res4;
      result[k-3] = res5;
      result[k-2] = res6;
      result[k-1] = res7;
   }

   return jxf_error_flag;
}

JXF_Int
jxf_SeqVectorMassInnerProd4( jxf_Vector *x, jxf_Vector **y, JXF_Int k, JXF_Real *result )
{
   JXF_Real *x_data = jxf_VectorData(x);
   JXF_Real *y_data = jxf_VectorData(y[0]);
   JXF_Int   size   = jxf_VectorSize(x);
   JXF_Int   i, j, restk;
   JXF_Real  res1, res2, res3, res4;
   JXF_Int   jstart, jstart1, jstart2, jstart3;

   restk = (k-(k/4*4));
   if (k > 3)
   {
      for (j = 0; j < k-3; j += 4)
      {
         res1 = 0;
         res2 = 0;
         res3 = 0;
         res4 = 0;
         jstart = j*size;
         jstart1 = jstart+size;
         jstart2 = jstart1+size;
         jstart3 = jstart2+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3,res4
#include "../../../include/jxf_smp_forloop.h"
         for (i = 0; i < size; i++)
         {
            res1 += y_data[jstart+i] * x_data[i];
            res2 += y_data[jstart1+i] * x_data[i];
            res3 += y_data[jstart2+i] * x_data[i];
            res4 += y_data[jstart3+i] * x_data[i];
         }
         result[j] = res1;
         result[j+1] = res2;
         result[j+2] = res3;
         result[j+3] = res4;
      }
   }
   if (restk == 1)
   {
      res1 = 0;
      jstart = (k-1)*size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
      }
      result[k-1] = res1;
   }
   else if (restk == 2)
   {
      res1 = 0;
      res2 = 0;
      jstart = (k-2)*size;
      jstart1 = jstart+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
      }
      result[k-2] = res1;
      result[k-1] = res2;
   }
   else if (restk == 3)
   {
      res1 = 0;
      res2 = 0;
      res3 = 0;
      jstart = (k-3)*size;
      jstart1 = jstart+size;
      jstart2 = jstart1+size;
#define JXF_SMP_PRIVATE i
#define JXF_SMP_REDUCTION_OP +
#define JXF_SMP_REDUCTION_VARS res1,res2,res3
#include "../../../include/jxf_smp_forloop.h"
      for (i = 0; i < size; i++)
      {
         res1 += y_data[jstart+i] * x_data[i];
         res2 += y_data[jstart1+i] * x_data[i];
         res3 += y_data[jstart2+i] * x_data[i];
      }
      result[k-3] = res1;
      result[k-2] = res2;
      result[k-1] = res3;
   }

   return jxf_error_flag;
}

jxf_Vector *
jxf_SeqVectorMGReorderByVariables( jxf_Vector *x, JXF_Int num_groups )
{
    JXF_Int size = jxf_VectorSize(x);
    JXF_Real *x_data = jxf_VectorData(x);
    jxf_Vector *y = jxf_SeqVectorCreate(size);
    JXF_Real *y_data = NULL;
    JXF_Int ng_p_two = num_groups + 2;
    JXF_Int sub_size = size / ng_p_two;
    JXF_Int Rowx = 0, Rowy, mdo;
    
    jxf_SeqVectorInitialize(y);
    y_data = jxf_VectorData(y);
    for (Rowy = 0; Rowy < size; Rowy ++)
    {
        mdo = Rowy / sub_size;
        Rowx = ng_p_two * (Rowy - mdo * sub_size) + mdo;
        y_data[Rowy] = x_data[Rowx];
    }
    
    return y;
}

jxf_Vector *
jxf_SeqVectorCloneDeep( jxf_Vector *x )
{
   JXF_Int size = jxf_VectorSize(x);
   JXF_Int num_vectors = jxf_VectorNumVectors(x);
   jxf_Vector *y = jxf_SeqMultiVectorCreate(size, num_vectors);
   jxf_VectorMultiVecStorageMethod(y) = jxf_VectorMultiVecStorageMethod(x);
   jxf_VectorVectorStride(y) = jxf_VectorVectorStride(x);
   jxf_VectorIndexStride(y) = jxf_VectorIndexStride(x);
   jxf_SeqVectorInitialize(y);
   jxf_SeqVectorCopy(x, y);
   return y;
}
