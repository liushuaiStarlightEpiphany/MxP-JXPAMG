//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_vecop.c -- basic operations for parallel vectors
 *  Date: 2011/09/03
 */ 

#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_ParVectorCopy
 * \brief Copy x -> y.
 * \date 2011/09/03
 */ 
JXF_Int
jxf_ParVectorCopy( jxf_ParVector *x, jxf_ParVector *y )
{
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   return jxf_SeqVectorCopy(x_local, y_local);
}

/*!
 * \fn JXF_Int jxf_ParVectorSetConstantValues
 * \brief Set a parallel vector to be a given constant.
 * \date 2011/09/03
 */   
JXF_Int
jxf_ParVectorSetConstantValues( jxf_ParVector *v, JXF_Real value )
{
   jxf_Vector *v_local = jxf_ParVectorLocalVector(v);
           
   return jxf_SeqVectorSetConstantValues(v_local,value);
}

/*!
 * \fn JXF_Real jxf_ParVectorNorm1
 * \brief Compute the 1-norm of vector x.
 * \author peghoty 
 * \date 2009/07/23
 */ 
JXF_Real
jxf_ParVectorNorm1( jxf_ParVector *x )
{
   MPI_Comm      comm = jxf_ParVectorComm(x);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);

   JXF_Real result = 0.0;
   JXF_Real local_result = jxf_SeqVectorNorm1(x_local);

   jxf_MPI_Allreduce(&local_result, &result, 1, JXF_MPI_REAL, MPI_MAX, comm);

   return result;
}

/*!
 * \fn JXF_Real jxf_ParVectorPointWiseRelNorm1
 * \brief Compute the pointwise relative 1-norm of vector x and y.
 * \author peghoty 
 * \date 2009/07/23
 */ 
JXF_Real
jxf_ParVectorPointWiseRelNorm1(jxf_ParVector *x, jxf_ParVector *y)
{
   MPI_Comm      comm = jxf_ParVectorComm(x);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   
   JXF_Real result = 0.0;
   JXF_Real local_result = jxf_SeqVectorPointWiseRelNorm1(x_local, y_local);
   
   jxf_MPI_Allreduce(&local_result, &result, 1, JXF_MPI_REAL, MPI_MAX, comm);

   return result;
}

/*!
 * \fn JXF_Int jxf_ParVectorScale
 * \brief Perform y := alpha*y. 
 * \date 2011/09/03
 */ 
JXF_Int
jxf_ParVectorScale( JXF_Real alpha, jxf_ParVector *y )
{
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);

   return jxf_SeqVectorScale( alpha, y_local);
}

/*!
 * \fn JXF_Int jxf_ParVecMul
 * \brief Performs z = x*y, where z_i is defiend by x_i*y_i.
 * \note x, y, z should be allocated and have the same partitioning.
 * \author peghoty
 * \date 2011/09/09 
 */ 
JXF_Int
jxf_ParVecMul( jxf_ParVector *x,
              jxf_ParVector *y,
              jxf_ParVector *z )
{
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   jxf_Vector *z_local = jxf_ParVectorLocalVector(z);
   
   jxf_SeqVecMul(x_local, y_local, z_local);
   
   return 0;
}

/*!
 * \fn JXF_Int jxf_ParVecZXY
 * \brief Performs z = z + alpha*x*y, where xy_i is defiend by x_i*y_i.
 * \note x, y, z should be allocated and have the same partitioning.
 * \author peghoty
 * \date 2011/09/09 
 */ 
JXF_Int
jxf_ParVecZXY( jxf_ParVector  *z,
              JXF_Real         alpha,
              jxf_ParVector  *x,
              jxf_ParVector  *y )
{
   jxf_Vector *z_local = jxf_ParVectorLocalVector(z);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
     
   if (alpha == 0) 
   {
      return 0;
   }
     
   jxf_SeqVecZXY(z_local, alpha, x_local, y_local);
   
   return 0;
}

/*!
 * \fn JXF_Real jxf_ParVectorInnerProd
 * \brief Compute the inner product of two parallel vectors.
 * \date 2011/09/03
 */ 
JXF_Real
jxf_ParVectorInnerProd( jxf_ParVector *x, jxf_ParVector *y )
{
   MPI_Comm      comm = jxf_ParVectorComm(x);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
           
   JXF_Real result = 0.0;
   JXF_Real local_result = jxf_SeqVectorInnerProd(x_local, y_local);
   
   jxf_MPI_Allreduce(&local_result, &result, 1, JXF_MPI_REAL, MPI_SUM, comm);
   
   return result;
}

JXF_Int
jxf_ParVectorMassInnerProd( jxf_ParVector *x, jxf_ParVector **y, JXF_Int k, JXF_Int unroll, JXF_Real *result )
{
   MPI_Comm   comm    = jxf_ParVectorComm(x);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   JXF_Real *local_result;
   JXF_Int i;
   jxf_Vector **y_local;

   y_local = jxf_TAlloc(jxf_Vector *, k);
   for (i = 0; i < k; i ++)
   {
      y_local[i] = (jxf_Vector *)jxf_ParVectorLocalVector(y[i]);
   }
   local_result = jxf_CTAlloc(JXF_Real, k);
   jxf_SeqVectorMassInnerProd(x_local, y_local, k, unroll, local_result);
   jxf_MPI_Allreduce(local_result, result, k, JXF_MPI_REAL, MPI_SUM, comm);
   jxf_TFree(y_local);
   jxf_TFree(local_result);
 
   return jxf_error_flag;
}

/*!
 * \fn JXF_Real jxf_ParVectorNorm2
 * \brief Compute the 2-norm of vector x.
 * \author peghoty 
 * \date 2009/07/23
 */ 
JXF_Real
jxf_ParVectorNorm2( jxf_ParVector *x )
{
   JXF_Real result = 0.0;
   result = sqrt( jxf_ParVectorInnerProd(x,x) );
   return result;
}

/*!
 * \fn JXF_Int jxf_ParVectorAxpy
 * \brief Perform y = alpha*x + y 
 * \date 2011/09/03
 */ 
JXF_Int
jxf_ParVectorAxpy( JXF_Real        alpha,
                  jxf_ParVector *x,
                  jxf_ParVector *y  )
{
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   
   return jxf_SeqVectorAxpy( alpha, x_local, y_local);
}

JXF_Int
jxf_ParVectorMassDotpTwo( jxf_ParVector *x, jxf_ParVector *y, jxf_ParVector **z, JXF_Int k, JXF_Int unroll, JXF_Real *result_x, JXF_Real *result_y )
{
   MPI_Comm   comm    = jxf_ParVectorComm(x);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   JXF_Real *local_result, *result;
   JXF_Int i; 
   jxf_Vector **z_local;
   z_local = jxf_TAlloc(jxf_Vector *, k);

   for (i = 0; i < k; i ++)
      z_local[i] = (jxf_Vector *) jxf_ParVectorLocalVector(z[i]);

   local_result = jxf_CTAlloc(JXF_Real, 2*k);
   result = jxf_CTAlloc(JXF_Real, 2*k);

   jxf_SeqVectorMassDotpTwo(x_local, y_local, z_local, k, unroll, &local_result[0], &local_result[k]);

   jxf_MPI_Allreduce(local_result, result, 2*k, JXF_MPI_REAL, MPI_SUM, comm);

   for (i=0; i < k; i++)
   {
      result_x[i] = result[i];
      result_y[i] = result[k+i];
   }
   jxf_TFree(z_local);
   jxf_TFree(local_result);
   jxf_TFree(result);
 
   return jxf_error_flag;
}

JXF_Int
jxf_ParVectorMassAxpy( JXF_Real *alpha, jxf_ParVector **x, jxf_ParVector *y, JXF_Int k, JXF_Int unroll )
{
   JXF_Int i;
   jxf_Vector **x_local;
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   x_local = jxf_TAlloc(jxf_Vector *, k);

   for (i=0; i < k; i++)
      x_local[i] = jxf_ParVectorLocalVector(x[i]);

   jxf_SeqVectorMassAxpy(alpha, x_local, y_local, k, unroll);

   jxf_TFree(x_local);

   return jxf_error_flag;
}
