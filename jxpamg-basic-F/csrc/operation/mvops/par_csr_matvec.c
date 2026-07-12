//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_csr_matvec.c -- basic operations for parallel matrix-vector multiplication.
 *  Date: 2011/09/05
 */

#include "jxf_mv.h"

extern JXF_Int myid;
extern JXF_Int coreNums;
extern JXF_Int jxf_spmv_type;

JXF_Real jxf_total_elapsed_time_matvec = 0.0;

JXF_Int
jxf_ParCSRMatrixSpMVPrecondFP32Create(jxf_ParCSRMatrix *A,
                                      JXF_Int myid)
{
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(A);
   JXF_Int x_size = jxf_CSRMatrixNumCols(diag);

   return jxf_CSRMatrixSpMVPrecondFP32Create(diag, x_size, myid);
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixMatvec
 * \brief Perform y := alpha*A*x + beta*y.
 * \date 2011/09/05
 */
JXF_Int
jxf_ParCSRMatrixMatvec(JXF_Real alpha,
                       jxf_ParCSRMatrix *A,
                       jxf_ParVector *x,
                       JXF_Real beta,
                       jxf_ParVector *y)
{
   jxf_ParCSRCommHandle **comm_handle;
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(A);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   JXF_Int num_rows = jxf_ParCSRMatrixGlobalNumRows(A);
   JXF_Int num_cols = jxf_ParCSRMatrixGlobalNumCols(A);

   jxf_Vector *x_tmp;
   JXF_Int x_size = jxf_ParVectorGlobalSize(x);
   JXF_Int y_size = jxf_ParVectorGlobalSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(x_local);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(offd);
   JXF_Int ierr = 0;
   JXF_Int num_sends, i, j, jv, index, start;

   JXF_Int vecstride = jxf_VectorVectorStride(x_local);
   JXF_Int idxstride = jxf_VectorIndexStride(x_local);

   JXF_Real *x_tmp_data, **x_buf_data;
   JXF_Real *x_local_data = jxf_VectorData(x_local);

   JXF_Real wall_time = 0.0; /* for debugging instrumentation  */

   if (jxf__global_mvcpu_flag)
      wall_time = jxf_time_getWallclockSeconds();

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  ParMatvec returns ierr = 11 if
    *  length of X doesn't equal the number of columns of A,
    *  ierr = 12 if the length of Y doesn't equal the number of rows
    *  of A, and ierr = 13 if both are true.
    *
    *  Because temporary vectors are often used in ParMatvec, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/

   jxf_assert(idxstride > 0);

   if (num_cols != x_size)
      ierr = 11;

   if (num_rows != y_size)
      ierr = 12;

   if (num_cols != x_size && num_rows != y_size)
      ierr = 13;

   jxf_assert(jxf_VectorNumVectors(y_local) == num_vectors);

   if (num_vectors == 1)
   {
      x_tmp = jxf_SeqVectorCreate(num_cols_offd);
   }
   else
   {
      jxf_assert(num_vectors > 1);
      x_tmp = jxf_SeqMultiVectorCreate(num_cols_offd, num_vectors);
   }
   jxf_SeqVectorInitialize(x_tmp);
   x_tmp_data = jxf_VectorData(x_tmp);

   comm_handle = jxf_CTAlloc(jxf_ParCSRCommHandle *, num_vectors);

   /*---------------------------------------------------------------------
    * If there exists no CommPkg for A, a CommPkg is generated using
    * equally load balanced partitionings
    *--------------------------------------------------------------------*/

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   x_buf_data = jxf_CTAlloc(JXF_Real *, num_vectors);
   for (jv = 0; jv < num_vectors; ++jv)
   {
      x_buf_data[jv] = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   }

   if (num_vectors == 1)
   {
      index = 0;
      for (i = 0; i < num_sends; i++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
         {
            x_buf_data[0][index++] = x_local_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
         }
      }
   }
   else
   {
      for (jv = 0; jv < num_vectors; ++jv)
      {
         index = 0;
         for (i = 0; i < num_sends; i++)
         {
            start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
            for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
            {
               x_buf_data[jv][index++] = x_local_data[jv * vecstride + idxstride * jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)];
            }
         }
      }
   }

   jxf_assert(idxstride == 1);

   /* >>> ... The assert is because the following loop only works for 'column' storage of a multivector <<<
      >>> This needs to be fixed to work more generally, at least for 'row' storage. <<<
      >>> This in turn, means either change CommPkg so num_sends is no.zones*no.vectors (not no.zones)
      >>> or, less dangerously, put a stride in the logic of CommHandleCreate (stride either from a
      >>> new arg or a new variable inside CommPkg).  Or put the num_vector iteration inside
      >>> CommHandleCreate (perhaps a new multivector variant of it). */

   for (jv = 0; jv < num_vectors; ++jv)
   {
      comm_handle[jv] = jxf_ParCSRCommHandleCreate(1, comm_pkg, x_buf_data[jv], &(x_tmp_data[jv * num_cols_offd]));
   }

   jxf_CSRMatrixMatvec(alpha, diag, x_local, beta, y_local, myid);

   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle[jv]);
      comm_handle[jv] = NULL;
   }
   jxf_TFree(comm_handle);

   if (num_cols_offd)
   {
      jxf_CSRMatrixMatvec_origin(alpha, offd, x_tmp, 1.0, y_local);
   }

   jxf_SeqVectorDestroy(x_tmp);
   x_tmp = NULL;
   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_TFree(x_buf_data[jv]);
   }
   jxf_TFree(x_buf_data);

   if (jxf__global_mvcpu_flag)
      jxf_total_elapsed_time_matvec += (jxf_time_getWallclockSeconds() - wall_time);

   return ierr;
}

// y = alpha*A*x + beta*b
JXF_Int
jxf_ParCSRMatrixMatvecOutOfPlace(JXF_Real alpha,
                                 jxf_ParCSRMatrix *A,
                                 jxf_ParVector *x,
                                 JXF_Real beta,
                                 jxf_ParVector *b,
                                 jxf_ParVector *y)
{
   jxf_ParCSRCommHandle **comm_handle;
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(A);

   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(A);

   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *b_local = jxf_ParVectorLocalVector(b);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   jxf_Vector *x_tmp;

   JXF_Int num_rows = jxf_ParCSRMatrixGlobalNumRows(A);
   JXF_Int num_cols = jxf_ParCSRMatrixGlobalNumCols(A);
   JXF_Int x_size = jxf_ParVectorGlobalSize(x);
   JXF_Int b_size = jxf_ParVectorGlobalSize(b);
   JXF_Int y_size = jxf_ParVectorGlobalSize(y);

   JXF_Int num_vectors = jxf_VectorNumVectors(x_local);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(offd);
   JXF_Int ierr = 0;
   JXF_Int num_sends, jv;

   JXF_Int vecstride = jxf_VectorVectorStride(x_local);
   JXF_Int idxstride = jxf_VectorIndexStride(x_local);

   JXF_Real *x_tmp_data, **x_buf_data;
   JXF_Real *x_local_data = jxf_VectorData(x_local);

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  ParMatvec returns ierr = 11 if
    *  length of X doesn't equal the number of columns of A,
    *  ierr = 12 if the length of Y doesn't equal the number of rows
    *  of A, and ierr = 13 if both are true.
    *
    *  Because temporary vectors are often used in ParMatvec, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/
   jxf_assert(idxstride > 0);

   if (num_cols != x_size)
   {
      ierr = 11;
   }

   if (num_rows != y_size || num_rows != b_size)
   {
      ierr = 12;
   }

   if (num_cols != x_size && (num_rows != y_size || num_rows != b_size))
   {
      ierr = 13;
   }

   jxf_assert(jxf_VectorNumVectors(b_local) == num_vectors);
   jxf_assert(jxf_VectorNumVectors(y_local) == num_vectors);

   if (num_vectors == 1)
   {
      x_tmp = jxf_SeqVectorCreate(num_cols_offd);
   }
   else
   {
      jxf_assert(num_vectors > 1);
      x_tmp = jxf_SeqMultiVectorCreate(num_cols_offd, num_vectors);
   }

   /*---------------------------------------------------------------------
    * If there exists no CommPkg for A, a CommPkg is generated using
    * equally load balanced partitionings
    *--------------------------------------------------------------------*/
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   jxf_assert(num_cols_offd == jxf_ParCSRCommPkgRecvVecStart(comm_pkg, jxf_ParCSRCommPkgNumRecvs(comm_pkg)));
   jxf_assert(jxf_ParCSRCommPkgSendMapStart(comm_pkg, 0) == 0);

   comm_handle = jxf_CTAlloc(jxf_ParCSRCommHandle *, num_vectors);

   jxf_SeqVectorInitialize(x_tmp);
   x_tmp_data = jxf_VectorData(x_tmp);

   /* x_buff_data */
   x_buf_data = jxf_CTAlloc(JXF_Real *, num_vectors);

   for (jv = 0; jv < num_vectors; ++jv)
   {
      x_buf_data[jv] = jxf_TAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   }

   /* The assert is because the following loop only works for 'column'
      storage of a multivector. This needs to be fixed to work more generally,
      at least for 'row' storage. This in turn, means either change CommPkg so
      num_sends is no.zones*no.vectors (not no.zones) or, less dangerously, put
      a stride in the logic of CommHandleCreate (stride either from a new arg or
      a new variable inside CommPkg).  Or put the num_vector iteration inside
      CommHandleCreate (perhaps a new multivector variant of it).
   */

   jxf_assert(idxstride == 1);

   for (jv = 0; jv < num_vectors; ++jv)
   {
      JXF_Real *send_data = (JXF_Real *)x_buf_data[jv];
      JXF_Real *locl_data = x_local_data + jv * vecstride;

      JXF_Int i;
      /* pack send data on host */
#include "../../../include/jxf_smp_forloop.h"
      for (i = jxf_ParCSRCommPkgSendMapStart(comm_pkg, 0); i < jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends); i++)
      {
         send_data[i] = locl_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, i)];
      }
   }

   for (jv = 0; jv < num_vectors; ++jv)
   {
      comm_handle[jv] = jxf_ParCSRCommHandleCreate(1, comm_pkg, x_buf_data[jv], &x_tmp_data[jv * num_cols_offd]);
   }

   /* overlapped local computation */
   jxf_CSRMatrixMatvecOutOfPlace(alpha, diag, x_local, beta, b_local, y_local, 0);

   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle[jv]);
      comm_handle[jv] = NULL;
   }
   jxf_TFree(comm_handle);

   /* computation offd part */
   if (num_cols_offd)
   {
      jxf_CSRMatrixMatvec_origin(alpha, offd, x_tmp, 1.0, y_local);
   }

   jxf_SeqVectorDestroy(x_tmp);
   x_tmp = NULL;

   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_TFree(x_buf_data[jv]);
   }
   jxf_TFree(x_buf_data);

   return ierr;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixMatvecT
 * \brief Perform y := alpha*A^T*x + beta*y.
 * \date 2011/09/05
 */
JXF_Int
jxf_ParCSRMatrixMatvecT(JXF_Real alpha,
                        jxf_ParCSRMatrix *A,
                        jxf_ParVector *x,
                        JXF_Real beta,
                        jxf_ParVector *y)
{
   jxf_ParCSRCommHandle **comm_handle;
   jxf_ParCSRCommPkg *comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   jxf_CSRMatrix *diag = jxf_ParCSRMatrixDiag(A);
   jxf_CSRMatrix *offd = jxf_ParCSRMatrixOffd(A);
   jxf_Vector *x_local = jxf_ParVectorLocalVector(x);
   jxf_Vector *y_local = jxf_ParVectorLocalVector(y);
   jxf_Vector *y_tmp;

   JXF_Int vecstride = jxf_VectorVectorStride(y_local);
   JXF_Int idxstride = jxf_VectorIndexStride(y_local);
   JXF_Real *y_tmp_data, **y_buf_data;
   JXF_Real *y_local_data = jxf_VectorData(y_local);

   JXF_Int num_rows = jxf_ParCSRMatrixGlobalNumRows(A);
   JXF_Int num_cols = jxf_ParCSRMatrixGlobalNumCols(A);
   JXF_Int num_cols_offd = jxf_CSRMatrixNumCols(offd);
   JXF_Int x_size = jxf_ParVectorGlobalSize(x);
   JXF_Int y_size = jxf_ParVectorGlobalSize(y);
   JXF_Int num_vectors = jxf_VectorNumVectors(y_local);

   JXF_Int i, j, jv, index, start, num_sends;

   JXF_Int ierr = 0;

   JXF_Real wall_time = 0.0; /* for debugging instrumentation  */

   if (jxf__global_mvcpu_flag)
      wall_time = jxf_time_getWallclockSeconds();

   /*---------------------------------------------------------------------
    *  Check for size compatibility.  MatvecT returns ierr = 1 if
    *  length of X doesn't equal the number of rows of A,
    *  ierr = 2 if the length of Y doesn't equal the number of
    *  columns of A, and ierr = 3 if both are true.
    *
    *  Because temporary vectors are often used in MatvecT, none of
    *  these conditions terminates processing, and the ierr flag
    *  is informational only.
    *--------------------------------------------------------------------*/

   if (num_rows != x_size)
      ierr = 1;

   if (num_cols != y_size)
      ierr = 2;

   if (num_rows != x_size && num_cols != y_size)
      ierr = 3;

   comm_handle = jxf_CTAlloc(jxf_ParCSRCommHandle *, num_vectors);

   if (num_vectors == 1)
   {
      y_tmp = jxf_SeqVectorCreate(num_cols_offd);
   }
   else
   {
      y_tmp = jxf_SeqMultiVectorCreate(num_cols_offd, num_vectors);
   }
   jxf_SeqVectorInitialize(y_tmp);

   /*---------------------------------------------------------------------
    * If there exists no CommPkg for A, a CommPkg is generated using
    * equally load balanced partitionings
    *--------------------------------------------------------------------*/

   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(A);
      comm_pkg = jxf_ParCSRMatrixCommPkg(A);
   }

   num_sends = jxf_ParCSRCommPkgNumSends(comm_pkg);
   y_buf_data = jxf_CTAlloc(JXF_Real *, num_vectors);
   for (jv = 0; jv < num_vectors; ++jv)
   {
      y_buf_data[jv] = jxf_CTAlloc(JXF_Real, jxf_ParCSRCommPkgSendMapStart(comm_pkg, num_sends));
   }

   y_tmp_data = jxf_VectorData(y_tmp);
   y_local_data = jxf_VectorData(y_local);

   jxf_assert(idxstride == 1); /* >>> only 'column' storage of multivectors implemented so far */

   if (num_cols_offd)
   {
      // offdT is optional. Used only if it's present.
      if (A->offdT)
      {
         jxf_CSRMatrixMatvec_origin(alpha, A->offdT, x_local, 0.0, y_tmp);
      }
      else
      {
         jxf_CSRMatrixMatvecT(alpha, offd, x_local, 0.0, y_tmp);
      }
   }

   for (jv = 0; jv < num_vectors; ++jv)
   {
      /* >>> this is where we assume multivectors are 'column' storage */
      comm_handle[jv] = jxf_ParCSRCommHandleCreate(2, comm_pkg, &(y_tmp_data[jv * num_cols_offd]), y_buf_data[jv]);
   }

   if (A->diagT)
   {
      // diagT is optional. Used only if it's present.
      jxf_CSRMatrixMatvec(alpha, A->diagT, x_local, beta, y_local, myid);
   }
   else
   {
      jxf_CSRMatrixMatvecT(alpha, diag, x_local, beta, y_local);
   }

   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle[jv]);
      comm_handle[jv] = NULL;
   }
   jxf_TFree(comm_handle);

   if (num_vectors == 1)
   {
      index = 0;
      for (i = 0; i < num_sends; i++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
         {
            y_local_data[jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)] += y_buf_data[0][index++];
         }
      }
   }
   else
   {
      for (jv = 0; jv < num_vectors; ++jv)
      {
         index = 0;
         for (i = 0; i < num_sends; i++)
         {
            start = jxf_ParCSRCommPkgSendMapStart(comm_pkg, i);
            for (j = start; j < jxf_ParCSRCommPkgSendMapStart(comm_pkg, i + 1); j++)
            {
               y_local_data[jv * vecstride + idxstride * jxf_ParCSRCommPkgSendMapElmt(comm_pkg, j)] += y_buf_data[jv][index++];
            }
         }
      }
   }

   jxf_SeqVectorDestroy(y_tmp);
   y_tmp = NULL;
   for (jv = 0; jv < num_vectors; ++jv)
   {
      jxf_TFree(y_buf_data[jv]);
   }
   jxf_TFree(y_buf_data);

   if (jxf__global_mvcpu_flag)
      jxf_total_elapsed_time_matvec += (jxf_time_getWallclockSeconds() - wall_time);

   return ierr;
}
