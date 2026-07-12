//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  krylov.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"
#include "jxf_krylov.h"

/*!
 * \fn char *jxf_ParKrylovCAlloc
 */ 
char *
jxf_ParKrylovCAlloc( JXF_Int count, JXF_Int elt_size )
{
   return ( jxf_CAlloc( count, elt_size ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovFree
 */ 
JXF_Int
jxf_ParKrylovFree( char *ptr )
{
   JXF_Int ierr = 0;
   jxf_Free( ptr );
   return ierr;
}

/*!
 * \fn JXF_Int jxf_ParKrylovCommInfo
 */ 
JXF_Int
jxf_ParKrylovCommInfo( void *A, JXF_Int *my_id, JXF_Int *num_procs )
{
   MPI_Comm comm = jxf_hpCSRMatrixComm ((jxf_hpCSRMatrix *) A);
   jxf_MPI_Comm_size(comm, num_procs);
   jxf_MPI_Comm_rank(comm, my_id);
   return 0;
}

/*!
 * \fn void *jxf_ParKrylovCreateVector
 */ 
void *
jxf_ParKrylovCreateVector( void *vvector )
{
   jxf_ParVector *vector = vvector;
   jxf_ParVector *new_vector;

   new_vector = jxf_ParVectorCreate( jxf_ParVectorComm(vector),
				    jxf_ParVectorGlobalSize(vector),	
                                    jxf_ParVectorPartitioning(vector) );
   jxf_ParVectorSetPartitioningOwner(new_vector,0);
   jxf_ParVectorInitialize(new_vector);

   return ( (void *) new_vector );
}

/*!
 * \fn void *jxf_ParKrylovCreateVectorArray
 */ 
void *
jxf_ParKrylovCreateVectorArray( JXF_Int n, void *vvector )
{
   jxf_ParVector  *vector = vvector;
   jxf_ParVector **new_vector;
   JXF_Int i, j;

   new_vector = jxf_CTAlloc(jxf_ParVector*, n);
   for (i = 0; i < n; i ++)
   {
      new_vector[i] = jxf_ParVectorCreate( jxf_ParVectorComm(vector),
				          jxf_ParVectorGlobalSize(vector),	
                                          jxf_ParVectorPartitioning(vector) );
      jxf_ParVectorSetPartitioningOwner(new_vector[i], 0);
      if (jxf_ParVectorInitialize(new_vector[i]) == JXF_ERROR_MEMORY)  // Feng Chunsheng & Yue Xiaoqiang 2012/10/26
      {
         jxf__kdim_memory_error = i - 1;
         jxf_error_flag = 0;
         for (j = 0; j < i; j ++)
         {
            jxf_ParVectorDestroy(new_vector[j]);
         }
         jxf_TFree(new_vector);
         new_vector = NULL;
         break;
      }
   }

   return ( (void *) new_vector );
}

/*!
 * \fn JXF_Int jxf_ParKrylovDestroyVector
 */ 
JXF_Int
jxf_ParKrylovDestroyVector( void *vvector )
{
   jxf_ParVector *vector = vvector;
   return( jxf_ParVectorDestroy( vector ) );
}

/*!
 * \fn void *jxf_ParKrylovMatvecCreate
 */ 
void *
jxf_ParKrylovMatvecCreate( void *A, void *x )
{
   void *matvec_data;
   matvec_data = NULL;
   return ( matvec_data );
}

/*!
 * \fn JXF_Int jxf_ParKrylovMatvec
 */
JXF_Int
jxf_ParKrylovMatvec( void   *matvec_data,
                    JXF_Real  alpha,
                    void   *A,
                    void   *x,
                    JXF_Real  beta,
                    void   *y   )
{
   return ( jxf_ParCSRMatrixMatvec ( alpha,
                                    jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) A),
                                    (jxf_ParVector *) x,
                                    beta,
                                    (jxf_ParVector *) y ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovMatvecDestroy
 */
JXF_Int
jxf_ParKrylovMatvecDestroy( void *matvec_data )
{
   return 0;
}

/*!
 * \fn JXF_Real jxf_ParKrylovInnerProd
 */
JXF_Real
jxf_ParKrylovInnerProd( void *x, void *y )
{
   return ( jxf_ParVectorInnerProd( (jxf_ParVector *) x, (jxf_ParVector *) y ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovMassInnerProd
 */
JXF_Int
jxf_ParKrylovMassInnerProd( void *x, void **y, JXF_Int k, JXF_Int unroll, void  *result )
{
   return ( jxf_ParVectorMassInnerProd( (jxf_ParVector *) x, (jxf_ParVector **) y, k, unroll, (JXF_Real *)result ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovMassDotpTwo
 */
JXF_Int
jxf_ParKrylovMassDotpTwo( void *x, void *y, void **z, JXF_Int k, JXF_Int unroll, void *result_x, void *result_y )
{
   return ( jxf_ParVectorMassDotpTwo( (jxf_ParVector *) x, (jxf_ParVector *) y, (jxf_ParVector **) z,
              k, unroll, (JXF_Real *) result_x, (JXF_Real *) result_y ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovCopyVector
 */
JXF_Int
jxf_ParKrylovCopyVector( void *x, void *y )
{
   return ( jxf_ParVectorCopy( (jxf_ParVector *) x, (jxf_ParVector *) y ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovClearVector
 */
JXF_Int
jxf_ParKrylovClearVector( void *x )
{
   return ( jxf_ParVectorSetConstantValues( (jxf_ParVector *) x, 0.0 ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovScaleVector
 */
JXF_Int
jxf_ParKrylovScaleVector( JXF_Real alpha, void *x )
{
   return ( jxf_ParVectorScale( alpha, (jxf_ParVector *) x ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovAxpy
 */
JXF_Int
jxf_ParKrylovAxpy( JXF_Real alpha, void *x, void *y )
{
   return ( jxf_ParVectorAxpy( alpha, (jxf_ParVector *) x, (jxf_ParVector *) y ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovMassAxpy
 */
JXF_Int
jxf_ParKrylovMassAxpy( JXF_Real *alpha, void **x, void *y, JXF_Int k, JXF_Int unroll )
{
   return ( jxf_ParVectorMassAxpy( alpha, (jxf_ParVector **) x, (jxf_ParVector *) y , k, unroll ) );
}

/*!
 * \fn JXF_Int jxf_ParKrylovIdentitySetup
 */
JXF_Int
jxf_ParKrylovIdentitySetup( void *vdata, void *A )
{
   return 0;
}

/*!
 * \fn JXF_Int jxf_ParKrylovIdentity
 */
JXF_Int
jxf_ParKrylovIdentity( void *vdata, void *A, void *b, void *x )
{
   return( jxf_ParKrylovCopyVector( b, x ) );
}
