//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hp_csr_matvec.c -- basic operations for hierarchy parallel matrix-vector multiplication.
 *  Date: 2011/09/05
 */ 

#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

JXF_Real jxf_hp_total_elapsed_time_matvec = 0.0;

/*!
 * \fn JXF_Int jxf_ParCSRMatrixMatvec
 * \brief Perform y := alpha*A*x + beta*y.
 * \date 2011/09/05
 */
JXF_Int
jxf_hpCSRMatrixMatvec( JXF_Real           alpha,
              	       jxf_hpCSRMatrix *hp_A,
                       jxf_ParVector    *x,
                       JXF_Real           beta,
                       jxf_ParVector    *y)
{
   return jxf_ParCSRMatrixMatvec(alpha, jxf_hpCSRMatrixPar(hp_A), x, beta, y);
}

#if 0
/*!
// code backup:拆分为三级计算与通信重叠
 * \fn JXF_Int jxf_hpCSRMatrixMatvecMultiOverlaping
 * \brief Perform y := alpha*A*x + beta*y.
 * \date 2022/03/07
 * \author mrz
 */
JXF_Int
jxf_hpCSRMatrixMatvecMultiOverlaping( JXF_Real           alpha,
              	                  jxf_hpCSRMatrix *hp_A,
                                 jxf_ParVector    *x,
                                 JXF_Real           beta,
                                 jxf_ParVector    *y)
{
   //TODO创建结点外，结点内CPU外，CPU内核外矩阵
   jxf_hpCSRMatrixLevel *node_outside_matrix  = jxf_hpCSRMatrixNodeOutside(hp_A);
   jxf_hpCSRMatrixLevel *cpu_outside_matrix   = jxf_hpCSRMatrixCpuOutside(hp_A);
   jxf_hpCSRMatrixLevel *core_outside_matrix  = jxf_hpCSRMatrixCpu(hp_A);

   jxf_CSRMatrix         *diag      = jxf_hpCSRMatrixDiag(hp_A);
   jxf_CSRMatrix *offd_node_outside = jxf_hpCSRMatrixOffdlevel(node_outside_matrix);
   jxf_CSRMatrix *offd_cpu_outside  = jxf_hpCSRMatrixOffdlevel(cpu_outside_matrix);
   jxf_CSRMatrix *offd_core_outside = jxf_hpCSRMatrixOffdlevel(core_outside_matrix);

   jxf_Vector *x_node_outside, *x_cpu_outside, *x_core_outside;
   JXF_Real   *x_node_outside_data, *x_cpu_outside_data, *x_core_outside_data;
   JXF_Real  **x_node_outside_buf_data, **x_cpu_outside_buf_data, **x_core_outside_buf_data;

   jxf_ParCSRCommHandle **comm_handle_node_outside, **comm_handle_cpu_outside, **comm_handle_core_outside;
   
   jxf_ParCSRCommPkg	*node_outside_commpkg = jxf_hpCSRMatrixlevelCommPkg(node_outside_matrix);
   jxf_ParCSRCommPkg	*cpu_outside_commpkg  = jxf_hpCSRMatrixlevelCommPkg(cpu_outside_matrix);
   jxf_ParCSRCommPkg	*core_outside_commpkg = jxf_hpCSRMatrixlevelCommPkg(core_outside_matrix);
   
   jxf_Vector            *x_local  = jxf_ParVectorLocalVector(x);   
   jxf_Vector            *y_local  = jxf_ParVectorLocalVector(y);   
   JXF_Int         num_rows = jxf_hpCSRMatrixGlobalNumRows(hp_A);
   JXF_Int         num_cols = jxf_hpCSRMatrixGlobalNumCols(hp_A);

   JXF_Int        x_size = jxf_ParVectorGlobalSize(x);
   JXF_Int        y_size = jxf_ParVectorGlobalSize(y);
   JXF_Int        num_vectors   = jxf_VectorNumVectors(x_local);

   JXF_Int	      num_cols_node_outside = jxf_CSRMatrixNumCols(offd_node_outside);
   JXF_Int	      num_cols_cpu_outside = jxf_CSRMatrixNumCols(offd_cpu_outside);
   JXF_Int	      num_cols_core_outside = jxf_CSRMatrixNumCols(offd_core_outside);

   JXF_Int         ierr = 0;
   JXF_Int	      i, j, jv, index, start;
   JXF_Int         num_sends_node_outside, num_sends_cpu_outside, num_sends_core_outside;

   JXF_Int        vecstride = jxf_VectorVectorStride( x_local );
   JXF_Int        idxstride = jxf_VectorIndexStride( x_local );

   JXF_Real     *x_local_data = jxf_VectorData(x_local);

   JXF_Real      wall_time = 0.0;  /* for debugging instrumentation  */

   if (jxf__global_mvcpu_flag) wall_time = jxf_time_getWallclockSeconds();
 
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
 
   jxf_assert( idxstride > 0 );

    if (num_cols != x_size)
              ierr = 11;

    if (num_rows != y_size)
              ierr = 12;

    if (num_cols != x_size && num_rows != y_size)
              ierr = 13;

    jxf_assert( jxf_VectorNumVectors(y_local)==num_vectors );

    if ( num_vectors == 1 )
    {
      x_node_outside = jxf_SeqVectorCreate( num_cols_node_outside );
    }
    else
    {
      jxf_assert( num_vectors > 1 );
      x_node_outside = jxf_SeqMultiVectorCreate(num_cols_node_outside, num_vectors);
    }

    jxf_SeqVectorInitialize(x_node_outside);

   x_node_outside_data = jxf_VectorData(x_node_outside);
   
   comm_handle_node_outside = jxf_CTAlloc(jxf_ParCSRCommHandle*, num_vectors);

  /*---------------------------------------------------------------------
   * If there exists no CommPkg for A, a CommPkg is generated using
   * equally load balanced partitionings
   *--------------------------------------------------------------------*/
   //node外块的初始化
   if (!node_outside_commpkg)
   {  
      jxf_hpMatvecNodeOutsideCommPkgCreate(hp_A);
      node_outside_commpkg = jxf_hpCSRMatrixlevelCommPkg(node_outside_matrix);
   }

   num_sends_node_outside = jxf_ParCSRCommPkgNumSends(node_outside_commpkg);
   x_node_outside_buf_data = jxf_CTAlloc(JXF_Real *, num_vectors);
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_node_outside_buf_data[jv] = jxf_CTAlloc( JXF_Real, jxf_ParCSRCommPkgSendMapStart(node_outside_commpkg, num_sends_node_outside) );
   }
   if ( num_vectors == 1 )
   {
      index = 0;
      for (i = 0; i < num_sends_node_outside; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(node_outside_commpkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(node_outside_commpkg, i+1); j ++)
         {
            x_node_outside_buf_data[0][index++] = x_local_data[jxf_ParCSRCommPkgSendMapElmt(node_outside_commpkg,j)];
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
            start = jxf_ParCSRCommPkgSendMapStart(node_outside_commpkg, i);
            for (j = start; j < jxf_ParCSRCommPkgSendMapStart(node_outside_commpkg, i+1); j ++)
            {
               x_node_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jxf_ParCSRCommPkgSendMapElmt(node_outside_commpkg,j) ];
            }
         }
      }
   }

   jxf_assert( idxstride == 1 );
   for (jv = 0; jv < num_vectors; ++ jv)
   {
        comm_handle_node_outside[jv] 
      = jxf_ParCSRCommHandleCreate( 1, node_outside_commpkg, x_node_outside_buf_data[jv], &(x_node_outside_data[jv*num_cols_node_outside]) );
   }

//TODO完成对角块计算
   jxf_CSRMatrixMatvec(alpha, diag, x_local, beta, y_local);

//TODO完成CPU外块通信
    if ( num_vectors == 1 )
    {
      x_cpu_outside = jxf_SeqVectorCreate( num_cols_cpu_outside );
    }
    else
    {
      jxf_assert( num_vectors > 1 );
      x_cpu_outside = jxf_SeqMultiVectorCreate(num_cols_cpu_outside, num_vectors);
    }
    jxf_SeqVectorInitialize(x_cpu_outside);
    //TODO获得每个向量的数据值
   x_cpu_outside_data = jxf_VectorData(x_cpu_outside);
   
   //为多个handle创建空间
   comm_handle_cpu_outside = jxf_CTAlloc(jxf_ParCSRCommHandle*, num_vectors);

   if (!cpu_outside_commpkg)
   {
      jxf_hpMatvecCPUOutsideCommPkgCreate(hp_A);
      cpu_outside_commpkg = jxf_hpCSRMatrixlevelCommPkg(cpu_outside_matrix);
   }  

   num_sends_cpu_outside  = jxf_ParCSRCommPkgNumSends(cpu_outside_commpkg);
   x_cpu_outside_buf_data  = jxf_CTAlloc(JXF_Real *, num_vectors);

   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_cpu_outside_buf_data[jv] = jxf_CTAlloc( JXF_Real, jxf_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, num_sends_cpu_outside) );
   }

   if ( num_vectors == 1 )
   {

      index = 0;
      for (i = 0; i < num_sends_cpu_outside; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i+1); j ++)
         {
            x_cpu_outside_buf_data[0][index++] = x_local_data[jxf_ParCSRCommPkgSendMapElmt(cpu_outside_commpkg,j)];
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
            start = jxf_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i);
            for (j = start; j < jxf_ParCSRCommPkgSendMapStart(cpu_outside_commpkg, i+1); j ++)
            {
               x_cpu_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jxf_ParCSRCommPkgSendMapElmt(cpu_outside_commpkg,j) ];
            }
         }
      }

   }
   
   jxf_assert( idxstride == 1 );
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      comm_handle_cpu_outside[jv] 
      = jxf_ParCSRCommHandleCreate( 1, cpu_outside_commpkg, x_cpu_outside_buf_data[jv], &(x_cpu_outside_data[jv*num_cols_cpu_outside]) );
   }

//TODO完成结点外通信
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle_node_outside[jv]);
      comm_handle_node_outside[jv] = NULL;
   }
   jxf_TFree(comm_handle_node_outside);
//TODO计算结点外
   if (num_cols_node_outside) 
   {
      jxf_CSRMatrixMatvec(alpha, offd_node_outside, x_node_outside, 1.0, y_local);
   }

   jxf_SeqVectorDestroy(x_node_outside);
   x_node_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jxf_TFree(x_node_outside_buf_data[jv]);
   }
   jxf_TFree(x_node_outside_buf_data);

   //TODO初始化核外的通信
   if ( num_vectors == 1 )
   {
      x_core_outside = jxf_SeqVectorCreate( num_cols_core_outside );
   }
   else
   {
      jxf_assert( num_vectors > 1 );
      x_core_outside = jxf_SeqMultiVectorCreate(num_cols_core_outside, num_vectors);
      //x_tmp = jxf_SeqMultiVectorCreate( num_cols_offd, num_vectors );
   }
   jxf_SeqVectorInitialize(x_core_outside);
    //TODO获得每个向量的数据值
   x_core_outside_data = jxf_VectorData(x_core_outside);
    //x_tmp_data = jxf_VectorData(x_tmp);
   
   //为多个handle创建空间
   comm_handle_core_outside = jxf_CTAlloc(jxf_ParCSRCommHandle*, num_vectors);

   if (!core_outside_commpkg)
   {
      jxf_hpMatvecCoreOutsideCommPkgCreate(hp_A);
      core_outside_commpkg = jxf_hpCSRMatrixlevelCommPkg(core_outside_matrix);
   }
   num_sends_core_outside = jxf_ParCSRCommPkgNumSends(core_outside_commpkg);
   x_core_outside_buf_data = jxf_CTAlloc(JXF_Real *, num_vectors);

   for (jv = 0; jv < num_vectors; ++ jv)
   {
      x_core_outside_buf_data[jv] = jxf_CTAlloc( JXF_Real, jxf_ParCSRCommPkgSendMapStart(core_outside_commpkg, num_sends_core_outside) );
   }

   if ( num_vectors == 1 )
   {

      index = 0;
      for (i = 0; i < num_sends_core_outside; i ++)
      {
         start = jxf_ParCSRCommPkgSendMapStart(core_outside_commpkg, i);
         for (j = start; j < jxf_ParCSRCommPkgSendMapStart(core_outside_commpkg, i+1); j ++)
         {
            x_core_outside_buf_data[0][index++] = x_local_data[jxf_ParCSRCommPkgSendMapElmt(core_outside_commpkg,j)];
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
            start = jxf_ParCSRCommPkgSendMapStart(core_outside_commpkg, i);
            for (j = start; j < jxf_ParCSRCommPkgSendMapStart(core_outside_commpkg, i+1); j ++)
            {
               x_core_outside_buf_data[jv][index++] 
               = x_local_data[ jv*vecstride + idxstride*jxf_ParCSRCommPkgSendMapElmt(core_outside_commpkg,j) ];
            }
         }
      }

   }

   jxf_assert( idxstride == 1 );
   
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
      = jxf_ParCSRCommHandleCreate( 1, core_outside_commpkg, x_core_outside_buf_data[jv], &(x_core_outside_data[jv*num_cols_core_outside]) );
   }

//TODO完成结点内CPU外通信
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle_cpu_outside[jv]);
      comm_handle_cpu_outside[jv] = NULL;
   }
   jxf_TFree(comm_handle_cpu_outside);

//TODO计算结点内CPU外
   if (num_cols_cpu_outside) 
   {
      jxf_CSRMatrixMatvec(alpha, offd_cpu_outside, x_cpu_outside, 1.0, y_local);
   }

   jxf_SeqVectorDestroy(x_cpu_outside);
   x_cpu_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jxf_TFree(x_cpu_outside_buf_data[jv]);
   }
   jxf_TFree(x_cpu_outside_buf_data);

//TODO完成CPU内核外通信 
   for (jv = 0; jv < num_vectors; ++ jv)
   {
      jxf_ParCSRCommHandleDestroy(comm_handle_core_outside[jv]);
      comm_handle_core_outside[jv] = NULL;
   }
   jxf_TFree(comm_handle_core_outside);

//TODO计算CPU内核外
   if (num_cols_core_outside) 
   {
      jxf_CSRMatrixMatvec(alpha, offd_core_outside, x_core_outside, 1.0, y_local);
   }

   jxf_SeqVectorDestroy(x_core_outside);
   x_core_outside = NULL;
   for (jv = 0; jv < num_vectors; ++ jv) 
   {
      jxf_TFree(x_core_outside_buf_data[jv]);
   }
   jxf_TFree(x_core_outside_buf_data);   

   if (jxf__global_mvcpu_flag) jxf_total_elapsed_time_matvec += (jxf_time_getWallclockSeconds() - wall_time);

   return ierr;
}
#endif
JXF_Int
jxf_hpCSRMatrixMatvecLevel( JXF_Real           alpha,
              	       jxf_hpCSRMatrix *hp_A,
                       jxf_ParVector    *x,
                       JXF_Real           beta,
                       jxf_ParVector    *y,
                       JXF_Int Level)
{
   JXF_Int ierr = 0;
   jxf_ParCSRMatrix *par_matrix = jxf_CTAlloc(jxf_ParCSRMatrix, 1);

   par_matrix = jxf_hpMatrixLevelToPar(hp_A, Level);
   jxf_ParCSRMatrixMatvec(alpha, par_matrix, x, beta, y);

   jxf_TFree(par_matrix);
   return ierr;
}

/*!
 * \fn JXF_Int jxf_ParCSRMatrixMatvecT
 * \brief Perform y := alpha*A^T*x + beta*y.
 * \date 2011/09/05
 */
JXF_Int
jxf_hpCSRMatrixMatvecT( JXF_Real           alpha,
                        jxf_hpCSRMatrix *hp_A,
                        jxf_ParVector    *x,
                        JXF_Real         beta,
                        jxf_ParVector    *y)
{
   return jxf_ParCSRMatrixMatvecT(alpha, jxf_hpCSRMatrixPar(hp_A), x, beta, y);
}
