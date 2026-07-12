//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hp_csr_matvec.c -- basic operations for hierarchy parallel matrix-vector multiplication.
 *  Date: 2011/09/05
 */ 

#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

JX_Real jx_hp_total_elapsed_time_matvec = 0.0;

/*!
 * \fn JX_Int jx_ParCSRMatrixMatvec
 * \brief Perform y := alpha*A*x + beta*y.
 * \date 2011/09/05
 */
JX_Int
jx_hpCSRMatrixMatvec( JX_Real           alpha,
              	       jx_hpCSRMatrix *hp_A,
                       jx_ParVector    *x,
                       JX_Real           beta,
                       jx_ParVector    *y)
{
   return jx_ParCSRMatrixMatvec(alpha, jx_hpCSRMatrixPar(hp_A), x, beta, y);
}

#if 0
/*!
// code backup:拆分为三级计算与通信重叠
 * \fn JX_Int jx_hpCSRMatrixMatvecMultiOverlaping
 * \brief Perform y := alpha*A*x + beta*y.
 * \date 2022/03/07
 * \author mrz
 */
JX_Int
jx_hpCSRMatrixMatvecMultiOverlaping( JX_Real           alpha,
              	                  jx_hpCSRMatrix *hp_A,
                                 jx_ParVector    *x,
                                 JX_Real           beta,
                                 jx_ParVector    *y)
{
   //TODO创建结点外，结点内CPU外，CPU内核外矩阵
   jx_hpCSRMatrixLevel *node_outside_matrix  = jx_hpCSRMatrixNodeOutside(hp_A);
   jx_hpCSRMatrixLevel *cpu_outside_matrix   = jx_hpCSRMatrixCpuOutside(hp_A);
   jx_hpCSRMatrixLevel *core_outside_matrix  = jx_hpCSRMatrixCpu(hp_A);

   jx_CSRMatrix         *diag      = jx_hpCSRMatrixDiag(hp_A);
   jx_CSRMatrix *offd_node_outside = jx_hpCSRMatrixOffdlevel(node_outside_matrix);
   jx_CSRMatrix *offd_cpu_outside  = jx_hpCSRMatrixOffdlevel(cpu_outside_matrix);
   jx_CSRMatrix *offd_core_outside = jx_hpCSRMatrixOffdlevel(core_outside_matrix);

   jx_Vector *x_node_outside, *x_cpu_outside, *x_core_outside;
   JX_Real   *x_node_outside_data, *x_cpu_outside_data, *x_core_outside_data;
   JX_Real  **x_node_outside_buf_data, **x_cpu_outside_buf_data, **x_core_outside_buf_data;

   jx_ParCSRCommHandle **comm_handle_node_outside, **comm_handle_cpu_outside, **comm_handle_core_outside;
   
   jx_ParCSRCommPkg	*node_outside_commpkg = jx_hpCSRMatrixlevelCommPkg(node_outside_matrix);
   jx_ParCSRCommPkg	*cpu_outside_commpkg  = jx_hpCSRMatrixlevelCommPkg(cpu_outside_matrix);
   jx_ParCSRCommPkg	*core_outside_commpkg = jx_hpCSRMatrixlevelCommPkg(core_outside_matrix);
   
   jx_Vector            *x_local  = jx_ParVectorLocalVector(x);   
   jx_Vector            *y_local  = jx_ParVectorLocalVector(y);   
   JX_Int         num_rows = jx_hpCSRMatrixGlobalNumRows(hp_A);
   JX_Int         num_cols = jx_hpCSRMatrixGlobalNumCols(hp_A);

   JX_Int        x_size = jx_ParVectorGlobalSize(x);
   JX_Int        y_size = jx_ParVectorGlobalSize(y);
   JX_Int        num_vectors   = jx_VectorNumVectors(x_local);

   JX_Int	      num_cols_node_outside = jx_CSRMatrixNumCols(offd_node_outside);
   JX_Int	      num_cols_cpu_outside = jx_CSRMatrixNumCols(offd_cpu_outside);
   JX_Int	      num_cols_core_outside = jx_CSRMatrixNumCols(offd_core_outside);

   JX_Int         ierr = 0;
   JX_Int	      i, j, jv, index, start;
   JX_Int         num_sends_node_outside, num_sends_cpu_outside, num_sends_core_outside;

   JX_Int        vecstride = jx_VectorVectorStride( x_local );
   JX_Int        idxstride = jx_VectorIndexStride( x_local );

   JX_Real     *x_local_data = jx_VectorData(x_local);

   JX_Real      wall_time = 0.0;  /* for debugging instrumentation  */

   if (jx__global_mvcpu_flag) wall_time = jx_time_getWallclockSeconds();
 
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
 
   jx_assert( idxstride > 0 );

    if (num_cols != x_size)
              ierr = 11;

    if (num_rows != y_size)
              ierr = 12;

    if (num_cols != x_size && num_rows != y_size)
              ierr = 13;

    jx_assert( jx_VectorNumVectors(y_local)==num_vectors );

    if ( num_vectors == 1 )
    {
      x_node_outside = jx_SeqVectorCreate( num_cols_node_outside );
    }
    else
    {
      jx_assert( num_vectors > 1 );
      x_node_outside = jx_SeqMultiVectorCreate(num_cols_node_outside, num_vectors);
    }

    jx_SeqVectorInitialize(x_node_outside);

   x_node_outside_data = jx_VectorData(x_node_outside);
   
   comm_handle_node_outside = jx_CTAlloc(jx_ParCSRCommHandle*, num_vectors);

  /*---------------------------------------------------------------------
   * If there exists no CommPkg for A, a CommPkg is generated using
   * equally load balanced partitionings
   *--------------------------------------------------------------------*/
   //node外块的初始化
   if (!node_outside_commpkg)
   {  
      jx_hpMatvecNodeOutsideCommPkgCreate(hp_A);
      node_outside_commpkg = jx_hpCSRMatrixlevelCommPkg(node_outside_matrix);
   }

   num_sends_node_outside = jx_ParCSRCommPkgNumSends(node_outside_commpkg);
   x_node_outside_buf_data = jx_CTAlloc(JX_Real *, num_vectors);
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_node_outside_buf_data[jv] = jx_CTAlloc( JX_Real, jx_ParCSRCommPkgSendMapStart(node_outside_commpkg, num_sends_node_outside) );
   }
   if ( num_vectors == 1 )
   {
      index = 0;
      for (i = 0; i < num_sends_node_outside; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(node_outside_commpkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(node_outside_commpkg, i+1); j ++)
         {
            x_node_outside_buf_data[0][index++] = x_local_data[jx_ParCSRCommPkgSendMapElmt(node_outside_commpkg,j)];
         }
      }
   }
   else
   {
      for (jv = 0; jv < num_vectors; ++ jv)
      {
         index = 0;
         for (i = 0; i < num_cols_node_outside; i ++)
         {
            start = jx_ParCSRCommPkgSendMapStart(node_outside_commpkg, i);
            for (j = start; j < jx_ParCSRCommPkgSendMapStart(node_outside_commpkg, i+1); j ++)
            {
               x_node_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jx_ParCSRCommPkgSendMapElmt(node_outside_commpkg,j) ];
            }
         }
      }
   }

   jx_assert( idxstride == 1 );
   for (jv = 0; jv < num_vectors; ++ jv)
   {
        comm_handle_node_outside[jv] 
      = jx_ParCSRCommHandleCreate( 1, node_outside_commpkg, x_node_outside_buf_data[jv], &(x_node_outside_data[jv*num_cols_node_outside]) );
   }

//TODO完成对角块计算
   jx_CSRMatrixMatvec(alpha, diag, x_local, beta, y_local);

//TODO完成CPU外块通信
    if ( num_vectors == 1 )
    {
      x_cpu_outside = jx_SeqVectorCreate( num_cols_cpu_outside );
    }
    else
    {
      jx_assert( num_vectors > 1 );
      x_cpu_outside = jx_SeqMultiVectorCreate(num_cols_cpu_outside, num_vectors);
    }
    jx_SeqVectorInitialize(x_cpu_outside);
    //TODO获得每个向量的数据值
   x_cpu_outside_data = jx_VectorData(x_cpu_outside);
   
   //为多个handle创建空间
   comm_handle_cpu_outside = jx_CTAlloc(jx_ParCSRCommHandle*, num_vectors);

   if (!cpu_outside_commpkg)
   {
      jx_hpMatvecCPUOutsideCommPkgCreate(hp_A);
      cpu_outside_commpkg = jx_hpCSRMatrixlevelCommPkg(cpu_outside_matrix);
   }  

   num_sends_cpu_outside  = jx_ParCSRCommPkgNumSends(cpu_outside_commpkg);
   x_cpu_outside_buf_data  = jx_CTAlloc(JX_Real *, num_vectors);

   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_cpu_outside_buf_data[jv] = jx_CTAlloc( JX_Real, jx_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, num_sends_cpu_outside) );
   }

   if ( num_vectors == 1 )
   {

      index = 0;
      for (i = 0; i < num_sends_cpu_outside; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i+1); j ++)
         {
            x_cpu_outside_buf_data[0][index++] = x_local_data[jx_ParCSRCommPkgSendMapElmt(cpu_outside_commpkg,j)];
         }
      }
   }
   else
   {

      for (jv = 0; jv < num_vectors; ++ jv)
      {
         index = 0;
         for (i = 0; i < num_cols_cpu_outside; i ++)
         {
            start = jx_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i);
            for (j = start; j < jx_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i+1); j ++)
            {
               x_cpu_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jx_ParCSRCommPkgSendMapElmt(cpu_outside_commpkg,j) ];
            }
         }
      }

   }
   
   jx_assert( idxstride == 1 );
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      comm_handle_cpu_outside[jv] 
      = jx_ParCSRCommHandleCreate( 1, cpu_outside_commpkg, x_cpu_outside_buf_data[jv], &(x_cpu_outside_data[jv*num_cols_cpu_outside]) );
   }

//TODO完成结点外通信
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jx_ParCSRCommHandleDestroy(comm_handle_node_outside[jv]);
      comm_handle_node_outside[jv] = NULL;
   }
   jx_TFree(comm_handle_node_outside);
//TODO计算结点外
   if (num_cols_node_outside) 
   {
      jx_CSRMatrixMatvec(alpha, offd_node_outside, x_node_outside, 1.0, y_local);
   }

   jx_SeqVectorDestroy(x_node_outside);
   x_node_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jx_TFree(x_node_outside_buf_data[jv]);
   }
   jx_TFree(x_node_outside_buf_data);

   //TODO初始化核外的通信
   if ( num_vectors == 1 )
   {
      x_core_outside = jx_SeqVectorCreate( num_cols_core_outside );
   }
   else
   {
      jx_assert( num_vectors > 1 );
      x_core_outside = jx_SeqMultiVectorCreate(num_cols_core_outside, num_vectors);
      //x_tmp = jx_SeqMultiVectorCreate( num_cols_offd, num_vectors );
   }
   jx_SeqVectorInitialize(x_core_outside);
    //TODO获得每个向量的数据值
   x_core_outside_data = jx_VectorData(x_core_outside);
    //x_tmp_data = jx_VectorData(x_tmp);
   
   //为多个handle创建空间
   comm_handle_core_outside = jx_CTAlloc(jx_ParCSRCommHandle*, num_vectors);

   if (!core_outside_commpkg)
   {
      jx_hpMatvecCoreOutsideCommPkgCreate(hp_A);
      core_outside_commpkg = jx_hpCSRMatrixlevelCommPkg(core_outside_matrix);
   }
   num_sends_core_outside = jx_ParCSRCommPkgNumSends(core_outside_commpkg);
   x_core_outside_buf_data = jx_CTAlloc(JX_Real *, num_vectors);

   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_core_outside_buf_data[jv] = jx_CTAlloc( JX_Real, jx_ParCSRCommPkgSendMapStart(core_outside_commpkg, num_sends_core_outside) );
   }

   if ( num_vectors == 1 )
   {

      index = 0;
      for (i = 0; i < num_sends_core_outside; i ++)
      {
         start = jx_ParCSRCommPkgSendMapStart(core_outside_commpkg, i);
         for (j = start; j < jx_ParCSRCommPkgSendMapStart(core_outside_commpkg, i+1); j ++)
         {
            x_core_outside_buf_data[0][index++] = x_local_data[jx_ParCSRCommPkgSendMapElmt(core_outside_commpkg,j)];
         }
      }

   }
   else
   {

      for (jv = 0; jv < num_vectors; ++ jv)
      {
         index = 0;
         for (i = 0; i < num_cols_core_outside; i ++)
         {
            start = jx_ParCSRCommPkgSendMapStart(core_outside_commpkg, i);
            for (j = start; j < jx_ParCSRCommPkgSendMapStart(core_outside_commpkg, i+1); j ++)
            {
               x_core_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jx_ParCSRCommPkgSendMapElmt(core_outside_commpkg,j) ];
            }
         }
      }

   }

   jx_assert( idxstride == 1 );
   
   /* >>> ... The assert is because the following loop only works for 'column' storage of a multivector <<<
      >>> This needs to be fixed to work more generally, at least for 'row' storage. <<<
      >>> This in turn, means either change CommPkg so num_sends is no.zones*no.vectors (not no.zones)
      >>> or, less dangerously, put a stride in the logic of CommHandleCreate (stride either from a
      >>> new arg or a new variable inside CommPkg).  Or put the num_vector iteration inside
      >>> CommHandleCreate (perhaps a new multivector variant of it). */
   
   //TODO先通信所有块

   for (jv = 0; jv < num_vectors; ++ jv)
   {
      comm_handle_core_outside[jv] 
      = jx_ParCSRCommHandleCreate( 1, core_outside_commpkg, x_core_outside_buf_data[jv], &(x_core_outside_data[jv*num_cols_core_outside]) );
   }

//TODO完成结点内CPU外通信
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jx_ParCSRCommHandleDestroy(comm_handle_cpu_outside[jv]);
      comm_handle_cpu_outside[jv] = NULL;
   }
   jx_TFree(comm_handle_cpu_outside);

//TODO计算结点内CPU外
   if (num_cols_cpu_outside) 
   {
      jx_CSRMatrixMatvec(alpha, offd_cpu_outside, x_cpu_outside, 1.0, y_local);
   }

   jx_SeqVectorDestroy(x_cpu_outside);
   x_cpu_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jx_TFree(x_cpu_outside_buf_data[jv]);
   }
   jx_TFree(x_cpu_outside_buf_data);

//TODO完成CPU内核外通信 
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jx_ParCSRCommHandleDestroy(comm_handle_core_outside[jv]);
      comm_handle_core_outside[jv] = NULL;
   }
   jx_TFree(comm_handle_core_outside);

//TODO计算CPU内核外
   if (num_cols_core_outside) 
   {
      jx_CSRMatrixMatvec(alpha, offd_core_outside, x_core_outside, 1.0, y_local);
   }

   jx_SeqVectorDestroy(x_core_outside);
   x_core_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jx_TFree(x_core_outside_buf_data[jv]);
   }
   jx_TFree(x_core_outside_buf_data);   

   if (jx__global_mvcpu_flag) jx_total_elapsed_time_matvec += (jx_time_getWallclockSeconds() - wall_time);

   return ierr;
}
#endif
JX_Int
jx_hpCSRMatrixMatvecLevel( JX_Real           alpha,
              	       jx_hpCSRMatrix *hp_A,
                       jx_ParVector    *x,
                       JX_Real           beta,
                       jx_ParVector    *y,
                       JX_Int Level)
{
   JX_Int ierr = 0;
   jx_ParCSRMatrix *par_matrix = jx_CTAlloc(jx_ParCSRMatrix, 1);

   par_matrix = jx_hpMatrixLevelToPar(hp_A, Level);
   jx_ParCSRMatrixMatvec(alpha, par_matrix, x, beta, y);

   jx_TFree(par_matrix);
   return ierr;
}

/*!
 * \fn JX_Int jx_ParCSRMatrixMatvecT
 * \brief Perform y := alpha*A^T*x + beta*y.
 * \date 2011/09/05
 */
JX_Int
jx_hpCSRMatrixMatvecT( JX_Real           alpha,
                        jx_hpCSRMatrix *hp_A,
                        jx_ParVector    *x,
                        JX_Real         beta,
                        jx_ParVector    *y)
{
   return jx_ParCSRMatrixMatvecT(alpha, jx_hpCSRMatrixPar(hp_A), x, beta, y);
}
