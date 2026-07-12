//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ij_matrix.c -- basic operations for IJ matrices.
 *  Date: 2015/11/28
 */ 

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

JXF_Int
JXF_IJMatrixCreate( MPI_Comm     comm,
                   JXF_Int          ilower,
                   JXF_Int          iupper,
                   JXF_Int          jlower,
                   JXF_Int          jupper,
                   JXF_IJMatrix *matrix )
{
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   JXF_Int *info;
   JXF_Int num_procs;
   JXF_Int myid;

   jxf_IJMatrix *ijmatrix;

#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int  row0, col0, rowN, colN;
#else
   JXF_Int *recv_buf;
   JXF_Int i, i4;
   JXF_Int square;
#endif

   ijmatrix = jxf_CTAlloc(jxf_IJMatrix, 1);

   jxf_IJMatrixComm(ijmatrix)         = comm;
   jxf_IJMatrixObject(ijmatrix)       = NULL;
   jxf_IJMatrixTranslator(ijmatrix)   = NULL;
   jxf_IJMatrixAssumedPart(ijmatrix)  = NULL;
   jxf_IJMatrixObjectType(ijmatrix)   = JXF_UNITIALIZED;
   jxf_IJMatrixAssembleFlag(ijmatrix) = 0;
   jxf_IJMatrixPrintLevel(ijmatrix)   = 0;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &myid);

   if (ilower > iupper+1 || ilower < 0)
   {
      jxf_error_in_arg(2);
      jxf_TFree(ijmatrix);
      return jxf_error_flag;
   }

   if (iupper < -1)
   {
      jxf_error_in_arg(3);
      jxf_TFree(ijmatrix);
      return jxf_error_flag;
   }

   if (jlower > jupper+1 || jlower < 0)
   {
      jxf_error_in_arg(4);
      jxf_TFree(ijmatrix);
      return jxf_error_flag;
   }

   if (jupper < -1)
   {
      jxf_error_in_arg(5);
      jxf_TFree(ijmatrix);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION

   info = jxf_CTAlloc(JXF_Int, 2);

   row_partitioning = jxf_CTAlloc(JXF_Int, 2);
   col_partitioning = jxf_CTAlloc(JXF_Int, 2);

   row_partitioning[0] = ilower;
   row_partitioning[1] = iupper + 1;
   col_partitioning[0] = jlower;
   col_partitioning[1] = jupper + 1;

   /* now we need the global number of rows and columns as well
      as the global first row and column index */

   /* proc 0 has the first row and col */
   if (myid == 0)
   {
      info[0] = ilower;
      info[1] = jlower;
   }
   jxf_MPI_Bcast(info, 2, JXF_MPI_INT, 0, comm);
   row0 = info[0];
   col0 = info[1];

   /* proc (num_procs-1) has the last row and col */
   if (myid == (num_procs-1))
   {
      info[0] = iupper;
      info[1] = jupper;
   }
   jxf_MPI_Bcast(info, 2, JXF_MPI_INT, num_procs-1, comm);

   rowN = info[0];
   colN = info[1];

   jxf_IJMatrixGlobalFirstRow(ijmatrix) = row0;
   jxf_IJMatrixGlobalFirstCol(ijmatrix) = col0;
   jxf_IJMatrixGlobalNumRows(ijmatrix) = rowN - row0 + 1;
   jxf_IJMatrixGlobalNumCols(ijmatrix) = colN - col0 + 1;

   jxf_TFree(info);

#else

   info = jxf_CTAlloc(JXF_Int, 4);
   recv_buf = jxf_CTAlloc(JXF_Int, 4*num_procs);
   row_partitioning = jxf_CTAlloc(JXF_Int, num_procs+1);

   info[0] = ilower;
   info[1] = iupper;
   info[2] = jlower;
   info[3] = jupper;

   /* Generate row- and column-partitioning through information exchange
      across all processors, check whether the matrix is square, and
      if the partitionings match. i.e. no overlaps or gaps,
      if there are overlaps or gaps in the row partitioning or column
      partitioning , ierr will be set to -9 or -10, respectively */

   jxf_MPI_Allgather(info, 4, JXF_MPI_INT, recv_buf, 4, JXF_MPI_INT, comm);

   row_partitioning[0] = recv_buf[0];
   square = 1;
   for (i=0; i < num_procs-1; i++)
   {
      i4 = 4*i;
      if ( recv_buf[i4+1] != (recv_buf[i4+4]-1) )
      {
         jxf_error(JXF_ERROR_GENERIC);
         jxf_TFree(ijmatrix);
         jxf_TFree(info);
         jxf_TFree(recv_buf);
         jxf_TFree(row_partitioning);
   	 return jxf_error_flag;
      }
      else
      {
         row_partitioning[i+1] = recv_buf[i4+4];
      }

      if ((square && (recv_buf[i4] != recv_buf[i4+2])) || (recv_buf[i4+1] != recv_buf[i4+3]))
      {
         square = 0;
      }
   }
   i4 = (num_procs-1)*4;
   row_partitioning[num_procs] = recv_buf[i4+1]+1;

   if ((recv_buf[i4] != recv_buf[i4+2]) || (recv_buf[i4+1] != recv_buf[i4+3]))
   {
      square = 0;
   }

   if (square) col_partitioning = row_partitioning;
   else
   {
      col_partitioning = jxf_CTAlloc(JXF_Int, num_procs+1);
      col_partitioning[0] = recv_buf[2];
      for (i=0; i < num_procs-1; i++)
      {
         i4 = 4*i;
         if (recv_buf[i4+3] != recv_buf[i4+6]-1)
         {
           jxf_error(JXF_ERROR_GENERIC);
           jxf_TFree(ijmatrix);
           jxf_TFree(info);
           jxf_TFree(recv_buf);
           jxf_TFree(row_partitioning);
           jxf_TFree(col_partitioning);
   	   return jxf_error_flag;
         }
         else
   	   col_partitioning[i+1] = recv_buf[i4+6];
      }
      col_partitioning[num_procs] = recv_buf[num_procs*4-1]+1;
   }

   jxf_IJMatrixGlobalFirstRow(ijmatrix) = row_partitioning[0];
   jxf_IJMatrixGlobalFirstCol(ijmatrix) = col_partitioning[0];
   jxf_IJMatrixGlobalNumRows(ijmatrix) = row_partitioning[num_procs] - row_partitioning[0];
   jxf_IJMatrixGlobalNumCols(ijmatrix) = col_partitioning[num_procs] - col_partitioning[0];

   jxf_TFree(info);
   jxf_TFree(recv_buf);

#endif

   jxf_IJMatrixRowPartitioning(ijmatrix) = row_partitioning;
   jxf_IJMatrixColPartitioning(ijmatrix) = col_partitioning;

  *matrix = (JXF_IJMatrix) ijmatrix;

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixCreatehpCSR( jxf_IJMatrix *matrix)
{
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int *row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   JXF_Int *col_partitioning = jxf_IJMatrixColPartitioning(matrix);
   jxf_hpCSRMatrix *hp_matrix;
   JXF_Int *row_starts;
   JXF_Int *col_starts;
   JXF_Int num_procs;
   JXF_Int i;


   jxf_MPI_Comm_size(comm, &num_procs);
   

#ifdef JXF_NO_GLOBAL_PARTITION
   row_starts = jxf_CTAlloc(JXF_Int, 2);
   if (jxf_IJMatrixGlobalFirstRow(matrix))
      for (i=0; i < 2; i++)
	 row_starts[i] = row_partitioning[i] - jxf_IJMatrixGlobalFirstRow(matrix);
   else
      for (i=0; i < 2; i++)
	 row_starts[i] = row_partitioning[i];
   if (row_partitioning != col_partitioning)
   {
      col_starts = jxf_CTAlloc(JXF_Int, 2);
      if (jxf_IJMatrixGlobalFirstCol(matrix))
	 for (i=0; i < 2; i++)
	    col_starts[i] = col_partitioning[i] - jxf_IJMatrixGlobalFirstCol(matrix);
      else
	 for (i=0; i < 2; i++)
	    col_starts[i] = col_partitioning[i];
   }
   else
      col_starts = row_starts;

   hp_matrix = jxf_hpCSRMatrixCreate(comm, jxf_IJMatrixGlobalNumRows(matrix),
                                         jxf_IJMatrixGlobalNumCols(matrix),
                                         row_starts, col_starts, 0, 0, 0);

#else
   row_starts = jxf_CTAlloc(JXF_Int, num_procs+1);
   if (row_partitioning[0])
      for (i=0; i < num_procs+1; i++)
	 row_starts[i] = row_partitioning[i] - row_partitioning[0];
   else
      for (i=0; i < num_procs+1; i++)
	 row_starts[i] = row_partitioning[i];
   if (row_partitioning != col_partitioning)
   {
      col_starts = jxf_CTAlloc(JXF_Int, num_procs+1);
      if (col_partitioning[0])
	 for (i=0; i < num_procs+1; i++)
	    col_starts[i] = col_partitioning[i] - col_partitioning[0];
      else
	 for (i=0; i < num_procs+1; i++)
	    col_starts[i] = col_partitioning[i];
   }
   else
      col_starts = row_starts;
   hp_matrix = jxf_hpCSRMatrixCreate(comm, row_starts[num_procs],
                                         col_starts[num_procs],
                                         row_starts, col_starts, 0, 0, 0);
#endif

   jxf_IJMatrixObject(matrix) = hp_matrix;

   return jxf_error_flag;
}

JXF_Int
jxf_AuxParCSRMatrixCreate( jxf_AuxParCSRMatrix **aux_matrix,
			  JXF_Int  local_num_rows,
                       	  JXF_Int  local_num_cols,
			  JXF_Int *sizes )
{
   jxf_AuxParCSRMatrix  *matrix;

   matrix = jxf_CTAlloc(jxf_AuxParCSRMatrix, 1);

   jxf_AuxParCSRMatrixLocalNumRows(matrix) = local_num_rows;
   jxf_AuxParCSRMatrixLocalNumCols(matrix) = local_num_cols;

   if (sizes)
   {
      jxf_AuxParCSRMatrixRowSpace(matrix) = sizes;
   }
   else
   {
      jxf_AuxParCSRMatrixRowSpace(matrix) = NULL;
   }

   /* set defaults */
   jxf_AuxParCSRMatrixNeedAux(matrix) = 1;
   jxf_AuxParCSRMatrixMaxOffProcElmts(matrix) = 0;
   jxf_AuxParCSRMatrixCurrentNumElmts(matrix) = 0;
   jxf_AuxParCSRMatrixOffProcIIndx(matrix) = 0;
   jxf_AuxParCSRMatrixRowLength(matrix) = NULL;
   jxf_AuxParCSRMatrixAuxJ(matrix) = NULL;
   jxf_AuxParCSRMatrixAuxData(matrix) = NULL;
   jxf_AuxParCSRMatrixIndxDiag(matrix) = NULL;
   jxf_AuxParCSRMatrixIndxOffd(matrix) = NULL;
   /* stash for setting or adding off processor values */
   jxf_AuxParCSRMatrixOffProcI(matrix) = NULL;
   jxf_AuxParCSRMatrixOffProcJ(matrix) = NULL;
   jxf_AuxParCSRMatrixOffProcData(matrix) = NULL;
  *aux_matrix = matrix;
   
   return 0;
}

JXF_Int
JXF_IJMatrixInitialize( JXF_IJMatrix matrix)
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;

   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   if ( jxf_IJMatrixObjectType(ijmatrix) == JXF_HPCSR )
      jxf_IJMatrixInitializehpCSR(ijmatrix);
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixInitializehpCSR( jxf_IJMatrix *matrix)
{
   jxf_hpCSRMatrix *hp_matrix = jxf_IJMatrixObject(matrix);
   jxf_AuxParCSRMatrix *aux_matrix = jxf_IJMatrixTranslator(matrix);
   JXF_Int local_num_rows;

   if (jxf_IJMatrixAssembleFlag(matrix) == 0)
   {
      if (!hp_matrix)
      {
         jxf_IJMatrixCreatehpCSR(matrix);
         hp_matrix = jxf_IJMatrixObject(matrix);
      }
      local_num_rows = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix));
      if (!aux_matrix)
      {
         jxf_AuxParCSRMatrixCreate(&aux_matrix, local_num_rows, 
            jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_matrix)), NULL);
         jxf_IJMatrixTranslator(matrix) = aux_matrix;
      }

      jxf_hpCSRMatrixInitialize(hp_matrix);
      jxf_AuxParCSRMatrixInitialize(aux_matrix);
      if (!jxf_AuxParCSRMatrixNeedAux(aux_matrix))
      {
         JXF_Int i, *indx_diag, *indx_offd, *diag_i, *offd_i;
         diag_i = jxf_CSRMatrixI(jxf_hpCSRMatrixDiag(hp_matrix));
         offd_i = jxf_CSRMatrixI(jxf_hpCSRMatrixOffd(hp_matrix));
         indx_diag = jxf_AuxParCSRMatrixIndxDiag(aux_matrix);
         indx_offd = jxf_AuxParCSRMatrixIndxOffd(aux_matrix);
#if JXF_USING_OPENMP
#pragma omp parallel for private(i) schedule(static)
#endif
         for (i=0; i < local_num_rows; i++)
         {
	    indx_diag[i] = diag_i[i];
	    indx_offd[i] = offd_i[i];
         }
      }
   }
   else
   {
      if (!aux_matrix)
      {
         local_num_rows = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix));
         jxf_AuxParCSRMatrixCreate(&aux_matrix, local_num_rows,
            jxf_CSRMatrixNumCols(jxf_hpCSRMatrixDiag(hp_matrix)), NULL);
         jxf_AuxParCSRMatrixNeedAux(aux_matrix) = 0;
         jxf_IJMatrixTranslator(matrix) = aux_matrix;
      }
   }

   return jxf_error_flag;
}

JXF_Int
jxf_AuxParCSRMatrixInitialize( jxf_AuxParCSRMatrix *matrix )
{
   JXF_Int local_num_rows = jxf_AuxParCSRMatrixLocalNumRows(matrix);
   JXF_Int *row_space = jxf_AuxParCSRMatrixRowSpace(matrix);
   JXF_Int max_off_proc_elmts = jxf_AuxParCSRMatrixMaxOffProcElmts(matrix);
   JXF_Int **aux_j;
   JXF_Real **aux_data;
   JXF_Int i;

   if (local_num_rows < 0)
      return -1;
   if (local_num_rows == 0)
      return 0;
   /* allocate stash for setting or adding off processor values */
   if (max_off_proc_elmts > 0)
   {
      jxf_AuxParCSRMatrixOffProcI(matrix) = jxf_CTAlloc(JXF_Int, 2*max_off_proc_elmts);
      jxf_AuxParCSRMatrixOffProcJ(matrix) = jxf_CTAlloc(JXF_Int, max_off_proc_elmts);
      jxf_AuxParCSRMatrixOffProcData(matrix) = jxf_CTAlloc(JXF_Real, max_off_proc_elmts);
   }
   if (jxf_AuxParCSRMatrixNeedAux(matrix))
   {
      aux_j = jxf_CTAlloc(JXF_Int *, local_num_rows);
      aux_data = jxf_CTAlloc(JXF_Real *, local_num_rows);
      if (!jxf_AuxParCSRMatrixRowLength(matrix))
         jxf_AuxParCSRMatrixRowLength(matrix) = jxf_CTAlloc(JXF_Int, local_num_rows);
      if (row_space)
      {
         for (i=0; i < local_num_rows; i++)
         {
            aux_j[i] = jxf_CTAlloc(JXF_Int, row_space[i]);
            aux_data[i] = jxf_CTAlloc(JXF_Real, row_space[i]);
         }
      }
      else
      {
         row_space = jxf_CTAlloc(JXF_Int, local_num_rows);
         for (i=0; i < local_num_rows; i++)
         {
            row_space[i] = 30;
            aux_j[i] = jxf_CTAlloc(JXF_Int, 30);
            aux_data[i] = jxf_CTAlloc(JXF_Real, 30);
         }
         jxf_AuxParCSRMatrixRowSpace(matrix) = row_space;
      }
      jxf_AuxParCSRMatrixAuxJ(matrix) = aux_j;
      jxf_AuxParCSRMatrixAuxData(matrix) = aux_data;
   }
   else
   {
      jxf_AuxParCSRMatrixIndxDiag(matrix) = jxf_CTAlloc(JXF_Int, local_num_rows);
      jxf_AuxParCSRMatrixIndxOffd(matrix) = jxf_CTAlloc(JXF_Int, local_num_rows);
   }

   return 0;
}

JXF_Int
JXF_IJMatrixDestroy( JXF_IJMatrix matrix )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;

   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   if (ijmatrix)
   {
      if (jxf_IJMatrixRowPartitioning(ijmatrix) == jxf_IJMatrixColPartitioning(ijmatrix))
         jxf_TFree(jxf_IJMatrixRowPartitioning(ijmatrix));
      else
      {
         jxf_TFree(jxf_IJMatrixRowPartitioning(ijmatrix));
         jxf_TFree(jxf_IJMatrixColPartitioning(ijmatrix));
      }
      if (jxf_IJMatrixAssumedPart(ijmatrix))
         jxf_AssumedPartitionDestroy(jxf_IJMatrixAssumedPart(ijmatrix));
      if (jxf_IJMatrixObjectType(ijmatrix) == JXF_HPCSR)
         jxf_IJMatrixDestroyhpCSR(ijmatrix);
      else if (jxf_IJMatrixObjectType(ijmatrix) != -1)
      {
         jxf_error_in_arg(1);
         return jxf_error_flag;
      }
   }

   jxf_TFree(ijmatrix);

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixDestroyhpCSR( jxf_IJMatrix *matrix )
{
   jxf_hpCSRMatrixDestroy(jxf_IJMatrixObject(matrix));
   jxf_AuxParCSRMatrixDestroy(jxf_IJMatrixTranslator(matrix));
   return jxf_error_flag;
}

JXF_Int
jxf_AuxParCSRMatrixDestroy( jxf_AuxParCSRMatrix *matrix )
{
   JXF_Int ierr=0;
   JXF_Int i;
   JXF_Int num_rows;

   if (matrix)
   {
      num_rows = jxf_AuxParCSRMatrixLocalNumRows(matrix);
      if (jxf_AuxParCSRMatrixRowLength(matrix))
         jxf_TFree(jxf_AuxParCSRMatrixRowLength(matrix));
      if (jxf_AuxParCSRMatrixRowSpace(matrix))
         jxf_TFree(jxf_AuxParCSRMatrixRowSpace(matrix));
      if (jxf_AuxParCSRMatrixAuxJ(matrix))
      {
         for (i=0; i < num_rows; i++)
	    jxf_TFree(jxf_AuxParCSRMatrixAuxJ(matrix)[i]);
	 jxf_TFree(jxf_AuxParCSRMatrixAuxJ(matrix));
      }
      if (jxf_AuxParCSRMatrixAuxData(matrix))
      {
         for (i=0; i < num_rows; i++)
            jxf_TFree(jxf_AuxParCSRMatrixAuxData(matrix)[i]);
	 jxf_TFree(jxf_AuxParCSRMatrixAuxData(matrix));
      }
      if (jxf_AuxParCSRMatrixIndxDiag(matrix))
            jxf_TFree(jxf_AuxParCSRMatrixIndxDiag(matrix));
      if (jxf_AuxParCSRMatrixIndxOffd(matrix))
            jxf_TFree(jxf_AuxParCSRMatrixIndxOffd(matrix));
      if (jxf_AuxParCSRMatrixOffProcI(matrix))
      	    jxf_TFree(jxf_AuxParCSRMatrixOffProcI(matrix));
      if (jxf_AuxParCSRMatrixOffProcJ(matrix))
      	    jxf_TFree(jxf_AuxParCSRMatrixOffProcJ(matrix));
      if (jxf_AuxParCSRMatrixOffProcData(matrix))
      	    jxf_TFree(jxf_AuxParCSRMatrixOffProcData(matrix));
      jxf_TFree(matrix);
   }

   return ierr;
}

JXF_Int
jxf_AssumedPartitionDestroy( jxf_IJAssumedPart *apart )
{
   if (apart->storage_length > 0)
   {
      jxf_TFree(apart->proc_list);
      jxf_TFree(apart->row_start_list);
      jxf_TFree(apart->row_end_list);
      jxf_TFree(apart->sort_index);
   }
   jxf_TFree(apart);

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixSetObjectType( JXF_IJMatrix matrix, JXF_Int type )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;
   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_IJMatrixObjectType(ijmatrix) = type;

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixSetValues( JXF_IJMatrix    matrix,
                      JXF_Int            nrows,
                      JXF_Int           *ncols,
                      const JXF_Int     *rows,
                      const JXF_Int     *cols,
                      const JXF_Real  *values )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;

   if (nrows == 0)
      return jxf_error_flag;

   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   if (nrows < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (!ncols)
   {
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

   if (!rows)
   {
      jxf_error_in_arg(4);
      return jxf_error_flag;
   }

   if (!cols)
   {
      jxf_error_in_arg(5);
      return jxf_error_flag;
   }

   if (!values)
   {
      jxf_error_in_arg(6);
      return jxf_error_flag;
   }

   if (jxf_IJMatrixObjectType(ijmatrix) == JXF_HPCSR)
   {
      if (jxf_IJMatrixOMPFlag(ijmatrix))
	 return(jxf_IJMatrixSetValuesOMPhpCSR(ijmatrix, nrows, ncols, rows, cols, values));
      else
         return(jxf_IJMatrixSetValueshpCSR(ijmatrix, nrows, ncols, rows, cols, values));
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixSetValuesOMPhpCSR( jxf_IJMatrix  *matrix,
                               JXF_Int           nrows,
                               JXF_Int          *ncols,
                               const JXF_Int    *rows,
                               const JXF_Int    *cols,
                               const JXF_Real *values )
{
   jxf_hpCSRMatrix *hp_matrix;
   jxf_CSRMatrix *diag, *offd;
   jxf_AuxParCSRMatrix *aux_matrix;
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int num_procs, my_id;
   JXF_Int col_0, col_n;
   JXF_Int cancel_indx = 0;
   JXF_Int **aux_j = NULL;
   JXF_Real **aux_data = NULL;
   JXF_Int *row_length, *row_space;
   JXF_Int need_aux;
   JXF_Int *diag_i = NULL;
   JXF_Int *diag_j = NULL;
   JXF_Real *diag_data = NULL;
   JXF_Int *offd_i = NULL;
   JXF_Int *offd_j = NULL;
   JXF_Real *offd_data = NULL;
   JXF_Int first, pstart;
   //JXF_Int current_num_elmts;
   /*JXF_Int max_off_proc_elmts;*/
   JXF_Int off_proc_i_indx = 0;
   JXF_Int *off_proc_i = NULL;
   JXF_Int *off_proc_j = NULL;
   JXF_Int *value_start, *offproc_cnt;

   JXF_Int print_level = jxf_IJMatrixPrintLevel(matrix);
   JXF_Int max_num_threads;
   JXF_Int error_flag = 0;
   JXF_Int i1;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   max_num_threads = jxf_NumThreads();
   hp_matrix = jxf_IJMatrixObject(matrix);
   row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   col_partitioning = jxf_IJMatrixColPartitioning(matrix);

   value_start = jxf_CTAlloc(JXF_Int, max_num_threads+1);
   offproc_cnt = jxf_CTAlloc(JXF_Int, max_num_threads);

#ifdef JXF_NO_GLOBAL_PARTITION
   col_0 = col_partitioning[0];
   col_n = col_partitioning[1]-1;
   first =  jxf_IJMatrixGlobalFirstCol(matrix);
   pstart = 0;
#else
   col_0 = col_partitioning[my_id];
   col_n = col_partitioning[my_id+1]-1;
   first = col_partitioning[0];
   pstart = my_id;
#endif
   if (nrows < 0)
   {
      jxf_error_in_arg(2);
      if (print_level)
         jxf_printf("Error! nrows negative! JXF_IJMatrixSetValues\n");
      return jxf_error_flag;
   }

   if (jxf_IJMatrixAssembleFlag(matrix))  /* matrix already assembled*/
   {
      JXF_Int *col_map_offd = NULL;
      JXF_Int num_cols_offd;

      diag = jxf_hpCSRMatrixDiag(hp_matrix);
      diag_i = jxf_CSRMatrixI(diag);
      diag_j = jxf_CSRMatrixJ(diag);
      diag_data = jxf_CSRMatrixData(diag);
      offd = jxf_hpCSRMatrixOffd(hp_matrix);
      offd_i = jxf_CSRMatrixI(offd);
      num_cols_offd = jxf_CSRMatrixNumCols(offd);
      if (num_cols_offd)
      {
          col_map_offd = jxf_hpCSRMatrixColMapOffd(hp_matrix);
          offd_j = jxf_CSRMatrixJ(offd);
          offd_data = jxf_CSRMatrixData(offd);
      }
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      if (aux_matrix)
      {
//         current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
         off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
         cancel_indx = jxf_AuxParCSRMatrixCancelIndx(aux_matrix);
      }

#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
      {
         JXF_Int j_offd;
         JXF_Int num_threads, my_thread_num;
         JXF_Int len, rest, ns, ne;
         JXF_Int pos_diag, pos_offd;
         JXF_Int len_diag, len_offd;
         JXF_Int row_len;
         JXF_Int row_local;
         JXF_Int i, j, k, ii, n, row;
         JXF_Int not_found, size, indx, cnt1, col_indx;

         num_threads = jxf_NumActiveThreads();
         my_thread_num = jxf_GetThreadNum();

         len = nrows/num_threads; 
         rest = nrows - len*num_threads;

         if (my_thread_num < rest)
         {
            ns = my_thread_num*(len+1);
            ne = (my_thread_num+1)*(len+1);
         }
         else
         {
            ns = my_thread_num*len+rest;
            ne = (my_thread_num+1)*len+rest;
         }

         value_start[my_thread_num] = 0;
         for (ii=ns; ii < ne; ii++)
            value_start[my_thread_num] += ncols[ii];

#if JXF_USING_OPENMP
#pragma omp barrier
#endif  
         if (my_thread_num == 0)
         {
            for (i=0; i < max_num_threads; i++)
               value_start[i+1] += value_start[i];
         }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif  
         indx = 0;
	 if (my_thread_num)
	   indx = value_start[my_thread_num-1];
         for (ii=ns; ii < ne; ii++)
         {
            row = rows[ii];
            n = ncols[ii];
            /* processor owns the row */ 
            if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
            {
               row_local = row - row_partitioning[pstart];

               /* compute local row number */
               size = diag_i[row_local+1] - diag_i[row_local] + offd_i[row_local+1] - offd_i[row_local];

               if (n > size)
               {
                  jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                  error_flag++;
      	          if (print_level) jxf_printf (" row %d too long! \n", row);
      	          break;
               }
       
               pos_diag = diag_i[row_local];
               pos_offd = offd_i[row_local];
               len_diag = diag_i[row_local+1];
               len_offd = offd_i[row_local+1];
               not_found = 1;
      	
               for (i=0; i < n; i++)
               {
                  if (cols[indx] < col_0 || cols[indx] > col_n)
                  /* insert into offd */	
                  {
      	             j_offd = jxf_BinarySearch(col_map_offd,cols[indx]-first,
                                              num_cols_offd);
      	             if (j_offd == -1)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	                break;
      	             }
      	             for (j=pos_offd; j < len_offd; j++)
      	             {
      	                if (offd_j[j] == j_offd)
      	                {
                           offd_data[j] = values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	                break;
      	             }
      	             not_found = 1;
                  }
                  /* diagonal element */
      	          else if (cols[indx] == row)
      	          {
      	             if (diag_j[pos_diag] != row_local)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	                break;
      	             }
      	             diag_data[pos_diag] = values[indx];
      	          }
                  else  /* insert into diag */
                  {
      	             for (j=pos_diag; j < len_diag; j++)
      	             {
      	                if (diag_j[j] == (cols[indx]-col_0))
      	                {
                           diag_data[j] = values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	                break;
      	             }
                  }
                  indx++;
               }
            }
         /* processor does not own the row */
            else /*search for previous occurrences and cancel them */
	    {

   	       if (aux_matrix)
               {
                  col_indx = 0;
                  for (i=0; i < off_proc_i_indx; i=i+2)
                  {
	             row_len = off_proc_i[i+1];
	             if (off_proc_i[i] == row)
		     {
		        for (j=0; j < n; j++)
		        {
			   cnt1 = col_indx;
			   for (k=0; k < row_len; k++)
			   {
			      if (off_proc_j[cnt1] == cols[j])
			      {
                                 off_proc_j[cnt1++] = -1;
	                         offproc_cnt[my_thread_num]++;
               		         /*cancel_indx++;*/
			         /* if no repetition allowed */
                                 /* off_proc_j[col_indx] = -1;
                                    col_indx -= k;
                                 break; */
			      }
			      else
			      {
			         cnt1++;
			      }
			   }
		        }
		        col_indx += row_len;
                     }
                     else
                     {
		        col_indx += row_len;
                     }
                  }
	       }
	    } 
         } 
      } /*end parallel region */
   }
   else  /* matrix not assembled */
   {
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      if (aux_matrix)
      {
//         current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
         off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
         cancel_indx = jxf_AuxParCSRMatrixCancelIndx(aux_matrix);
      }
      row_space = jxf_AuxParCSRMatrixRowSpace(aux_matrix);
      row_length = jxf_AuxParCSRMatrixRowLength(aux_matrix);
      need_aux = jxf_AuxParCSRMatrixNeedAux(aux_matrix);
      if (need_aux)
      {
         aux_j = jxf_AuxParCSRMatrixAuxJ(aux_matrix);
         aux_data = jxf_AuxParCSRMatrixAuxData(aux_matrix);
      }
      else
      {
         diag = jxf_hpCSRMatrixDiag(hp_matrix);
         diag_i = jxf_CSRMatrixI(diag);
         diag_j = jxf_CSRMatrixJ(diag);
         diag_data = jxf_CSRMatrixData(diag);
         offd = jxf_hpCSRMatrixOffd(hp_matrix);
         offd_i = jxf_CSRMatrixI(offd);
         if (num_procs > 1)
	 {
	    offd_j = jxf_CSRMatrixJ(offd);
            offd_data = jxf_CSRMatrixData(offd);
         }
      }
#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
      {
         JXF_Int num_threads, my_thread_num;
         JXF_Int len, rest, ns, ne;
         JXF_Int *tmp_j = NULL;
         JXF_Int *local_j = NULL;
         JXF_Real *tmp_data = NULL;
         JXF_Real *local_data = NULL;
         JXF_Int tmp_indx;
         JXF_Int row_len;
         JXF_Int row_local;
         JXF_Int i, j, k, ii, n, row;
         JXF_Int not_found, size, indx, cnt1, col_indx;
         JXF_Int old_size, space, cnt;

         num_threads = jxf_NumActiveThreads();
         my_thread_num = jxf_GetThreadNum();

         len = nrows / num_threads; 
         rest = nrows - len*num_threads;

         if (my_thread_num < rest)
         {
            ns = my_thread_num*(len+1);
            ne = (my_thread_num+1)*(len+1);
         }
         else
         {
            ns = my_thread_num*len+rest;
            ne = (my_thread_num+1)*len+rest;
         }

         value_start[my_thread_num] = 0;
         for (ii=ns; ii < ne; ii++)
            value_start[my_thread_num] += ncols[ii];

#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         if (my_thread_num == 0)
         {
            for (i=0; i < max_num_threads; i++)
               value_start[i+1] += value_start[i];
         }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         indx = 0;
         if (my_thread_num)
            indx = value_start[my_thread_num-1];
         for (ii=ns; ii < ne; ii++)
         {
            row = rows[ii];
            n = ncols[ii];
            /* processor owns the row */ 
            if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
            {
               row_local = row - row_partitioning[pstart]; 
               /* compute local row number */
               if (need_aux)
               {
                  local_j = aux_j[row_local];
                  local_data = aux_data[row_local];
   	          space = row_space[row_local]; 
   	          old_size = row_length[row_local]; 
   	          size = space - old_size;
   	          if (size < n)
      	          {
      	             size = n - size;
      	             tmp_j = jxf_CTAlloc(JXF_Int, size);
      	             tmp_data = jxf_CTAlloc(JXF_Real, size);
      	          }
      	          tmp_indx = 0;
      	          not_found = 1;
      	          size = old_size;
                  for (i=0; i < n; i++)
      	          {
      	             for (j=0; j < old_size; j++)
      	             {
      	                if (local_j[j] == cols[indx])
      	                {
                           local_data[j] = values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
      	                if (size < space)
      	                {
      	                   local_j[size] = cols[indx];
      	                   local_data[size++] = values[indx];
      	                }
      	                else
      	                {
      	                   tmp_j[tmp_indx] = cols[indx];
      	                   tmp_data[tmp_indx++] = values[indx];
      	                }
      	             }
      	             not_found = 1;
        	     indx++;
      	          }
      	    
                  row_length[row_local] = size+tmp_indx;
                  
                  if (tmp_indx)
                  {
   	             aux_j[row_local] = jxf_TReAlloc(aux_j[row_local], JXF_Int, size+tmp_indx);
   	             aux_data[row_local] = jxf_TReAlloc(aux_data[row_local], JXF_Real, size+tmp_indx);
                     row_space[row_local] = size+tmp_indx;
                     local_j = aux_j[row_local];
                     local_data = aux_data[row_local];
                  }
   
   	          cnt = size; 
   
   	          for (i=0; i < tmp_indx; i++)
   	          {
   	             local_j[cnt] = tmp_j[i];
   	             local_data[cnt++] = tmp_data[i];
	          }
  
	          if (tmp_j)
	          { 
	             jxf_TFree(tmp_j); 
	             jxf_TFree(tmp_data); 
	          } 
               }
               else /* insert immediately into data in ParCSRMatrix structure */
               {
                  JXF_Int offd_indx, diag_indx;
                  JXF_Int offd_space, diag_space;
                  JXF_Int cnt_diag, cnt_offd;
	          offd_indx = jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local];
	          diag_indx = jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local];
	          cnt_diag = diag_indx;
	          cnt_offd = offd_indx;
	          diag_space = diag_i[row_local+1];
	          offd_space = offd_i[row_local+1];
	          not_found = 1;
  	          for (i=0; i < n; i++)
	          {
	             if (cols[indx] < col_0 || cols[indx] > col_n)
                     /* insert into offd */	
	             {
	                for (j=offd_i[row_local]; j < offd_indx; j++)
	                {
		           if (offd_j[j] == cols[indx])
		           {
                              offd_data[j] = values[indx];
		              not_found = 0;
		              break;
		           }
	                }
	                if (not_found)
	                { 
	                   if (cnt_offd < offd_space) 
	                   { 
	                      offd_j[cnt_offd] = cols[indx];
	                      offd_data[cnt_offd++] = values[indx];
	                   } 
	                   else 
	 	           {
                              jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                              error_flag++;
	    	              if (print_level)
                                 jxf_printf("Error in row %d ! Too many elements!\n", row);
                              break;
	 	           }
	                }  
	                not_found = 1;
	             }
	             else  /* insert into diag */
	             {
	                for (j=diag_i[row_local]; j < diag_indx; j++)
	                {
		           if (diag_j[j] == cols[indx])
		           {
                              diag_data[j] = values[indx];
		              not_found = 0;
		              break;
		           }
	                } 
	                if (not_found)
	                { 
	                   if (cnt_diag < diag_space) 
	                   { 
	                      diag_j[cnt_diag] = cols[indx];
	                      diag_data[cnt_diag++] = values[indx];
	                   } 
	                   else 
	 	           {
                              jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                              error_flag++;
	    	              if (print_level)
                                 jxf_printf("Error in row %d ! Too many elements !\n", row);
                              break;
	 	           }
	                } 
	                not_found = 1;
	             }
	             indx++;
	          }

                  jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local] = cnt_diag;
                  jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local] = cnt_offd;

               }
            }

            /* processor does not own the row */
            else
	    {
               indx +=  n;
   	       if (aux_matrix)
               {
                  col_indx = 0;
                  for (i=0; i < off_proc_i_indx; i=i+2)
                  {
	             row_len = off_proc_i[i+1];
	             if (off_proc_i[i] == row)
		     {
		        for (j=0; j < n; j++)
		        {
			   cnt1 = col_indx;
			   for (k=0; k < row_len; k++)
			   {
			      if (off_proc_j[cnt1] == cols[j])
			      {
                                 off_proc_j[cnt1++] = -1;
               		         /*cancel_indx++;*/
	                         offproc_cnt[my_thread_num]++;
			         /* if no repetition allowed */
                                 /* off_proc_j[col_indx] = -1;
                                 col_indx -= k;
                                 break; */
			      }
			      else
			      {
			         cnt1++;
			      }
			   }
		        }
		        col_indx += row_len;
                     }
                     else
                     {
		        col_indx += row_len;
                     }
                  }
               }
	    }
         }
      } /* end parallel region */
   }
   if (error_flag) return jxf_error_flag;
   if (aux_matrix)
   {
      for (i1=0; i1 < max_num_threads; i1++)
         cancel_indx += offproc_cnt[i1];
      jxf_AuxParCSRMatrixCancelIndx(aux_matrix) = cancel_indx;
   }
   jxf_TFree(value_start);
   jxf_TFree(offproc_cnt);
   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixSetValueshpCSR( jxf_IJMatrix  *matrix,
                            JXF_Int           nrows,
                            JXF_Int          *ncols,
                            const JXF_Int    *rows,
                            const JXF_Int    *cols,
                            const JXF_Real *values )
{
   jxf_hpCSRMatrix *hp_matrix;
   jxf_CSRMatrix *diag, *offd;
   jxf_AuxParCSRMatrix *aux_matrix;
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int num_procs, my_id;
   JXF_Int row_local, row;
   JXF_Int row_len;
   JXF_Int col_0, col_n;
   JXF_Int i, ii, j, k, n, not_found;
   JXF_Int col_indx, cancel_indx, cnt1;
   JXF_Int **aux_j;
   JXF_Int *local_j;
   JXF_Int *tmp_j;
   JXF_Real **aux_data;
   JXF_Real  *local_data;
   JXF_Real  *tmp_data = NULL;
   JXF_Int diag_space, offd_space;
   JXF_Int *row_length, *row_space;
   JXF_Int need_aux;
   JXF_Int tmp_indx, indx;
   JXF_Int space, size, old_size;
   JXF_Int cnt, cnt_diag, cnt_offd;
   JXF_Int pos_diag, pos_offd;
   JXF_Int len_diag, len_offd;
   JXF_Int offd_indx, diag_indx;
   JXF_Int *diag_i;
   JXF_Int *diag_j;
   JXF_Real *diag_data;
   JXF_Int *offd_i;
   JXF_Int *offd_j = NULL;
   JXF_Real *offd_data = NULL;
   JXF_Int first, pstart;
   //JXF_Int current_num_elmts;
   /*JXF_Int max_off_proc_elmts;*/
   JXF_Int off_proc_i_indx;
   JXF_Int *off_proc_i;
   JXF_Int *off_proc_j;
   JXF_Int print_level = jxf_IJMatrixPrintLevel(matrix);
   /*JXF_Real *off_proc_data;*/
   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   hp_matrix = jxf_IJMatrixObject( matrix );
   row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   col_partitioning = jxf_IJMatrixColPartitioning(matrix);

#ifdef JXF_NO_GLOBAL_PARTITION
   col_0 = col_partitioning[0];
   col_n = col_partitioning[1]-1;
   first =  jxf_IJMatrixGlobalFirstCol(matrix);
   pstart = 0;
#else
   col_0 = col_partitioning[my_id];
   col_n = col_partitioning[my_id+1]-1;
   first = col_partitioning[0];
   pstart = my_id;
#endif
   if (nrows < 0)
   {
      jxf_error_in_arg(2);
      if (print_level)
         jxf_printf("Error! nrows negative! JXF_IJMatrixSetValues\n");
   }

   if (jxf_IJMatrixAssembleFlag(matrix))  /* matrix already assembled*/
   {
      JXF_Int *col_map_offd = NULL;
      JXF_Int num_cols_offd;
      JXF_Int j_offd;
      indx = 0;   
      for (ii=0; ii < nrows; ii++)
      {
         row = rows[ii];
         n = ncols[ii];

         /* processor owns the row */ 
         if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
         {
            row_local = row - row_partitioning[pstart];

            /* compute local row number */
            diag = jxf_hpCSRMatrixDiag(hp_matrix);
            diag_i = jxf_CSRMatrixI(diag);
            diag_j = jxf_CSRMatrixJ(diag);
            diag_data = jxf_CSRMatrixData(diag);
            offd = jxf_hpCSRMatrixOffd(hp_matrix);
            offd_i = jxf_CSRMatrixI(offd);
            num_cols_offd = jxf_CSRMatrixNumCols(offd);
            if (num_cols_offd)
            {
               col_map_offd = jxf_hpCSRMatrixColMapOffd(hp_matrix);
               offd_j = jxf_CSRMatrixJ(offd);
               offd_data = jxf_CSRMatrixData(offd);
            }
            size = diag_i[row_local+1] - diag_i[row_local] + offd_i[row_local+1] - offd_i[row_local];
      
            if (n > size)
            {
               jxf_error(JXF_ERROR_GENERIC);
      	       if (print_level) jxf_printf (" row %d too long! \n", row);
               return jxf_error_flag;
            }

            pos_diag = diag_i[row_local];
            pos_offd = offd_i[row_local];
            len_diag = diag_i[row_local+1];
            len_offd = offd_i[row_local+1];
            not_found = 1;
      	
            for (i=0; i < n; i++)
            {
               if (cols[indx] < col_0 || cols[indx] > col_n)
                  /* insert into offd */	
               {
      	          j_offd = jxf_BinarySearch(col_map_offd,cols[indx]-first,
                                              num_cols_offd);
      	          if (j_offd == -1)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
                     return jxf_error_flag;
      	          }
      	          for (j=pos_offd; j < len_offd; j++)
      	          {
      	             if (offd_j[j] == j_offd)
      	             {
                        offd_data[j] = values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
                     return jxf_error_flag;
      	          }
      	          not_found = 1;
               }
               /* diagonal element */
      	       else if (cols[indx] == row)
      	       {
      	          if (diag_j[pos_diag] != row_local)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	             /* return -1;*/
                     return jxf_error_flag;
      	          }
      	          diag_data[pos_diag] = values[indx];
      	       }
               else  /* insert into diag */
               {
      	          for (j=pos_diag; j < len_diag; j++)
      	          {
      	             if (diag_j[j] == (cols[indx]-col_0))
      	             {
                        diag_data[j] = values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf (" Error, element %d %d does not exist\n", row, cols[indx]);
      	             /* return -1; */
                     return jxf_error_flag;
      	          }
               }
               indx++;
            }
         }
         
         /* processor does not own the row */  
        
         else /*search for previous occurrences and cancel them */
	 {
            aux_matrix = jxf_IJMatrixTranslator(matrix);
   	    if (aux_matrix)
            {
   	       //current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
   	       off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
   	       off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	       off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
               col_indx = 0;
               cancel_indx = jxf_AuxParCSRMatrixCancelIndx(aux_matrix);
               for (i=0; i < off_proc_i_indx; i=i+2)
               {
	          row_len = off_proc_i[i+1];
	          if (off_proc_i[i] == row)
		  {
		     for (j=0; j < n; j++)
		     {
			cnt1 = col_indx;
			for (k=0; k < row_len; k++)
			{
			   if (off_proc_j[cnt1] == cols[j])
			   {
                              off_proc_j[cnt1++] = -1;
               		      cancel_indx++;
			      /* if no repetition allowed */
                              /* off_proc_j[col_indx] = -1;
                                 col_indx -= k;
                                 break; */
			   }
			   else
			   {
			      cnt1++;
			   }
			}
		     }
		     col_indx += row_len;
                  }
                  else
                  {
		     col_indx += row_len;
                  }
               }
               jxf_AuxParCSRMatrixCancelIndx(aux_matrix) = cancel_indx;
	    }
	 } 
      } 
   }
   else
   {
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      row_space = jxf_AuxParCSRMatrixRowSpace(aux_matrix);
      row_length = jxf_AuxParCSRMatrixRowLength(aux_matrix);
      need_aux = jxf_AuxParCSRMatrixNeedAux(aux_matrix);
      indx = 0;   
      for (ii=0; ii < nrows; ii++)
      {
         row = rows[ii];
         n = ncols[ii];
         /* processor owns the row */ 
         if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
         {
            row_local = row - row_partitioning[pstart]; 
            /* compute local row number */
            if (need_aux)
            {
               aux_j = jxf_AuxParCSRMatrixAuxJ(aux_matrix);
               aux_data = jxf_AuxParCSRMatrixAuxData(aux_matrix);
               local_j = aux_j[row_local];
               local_data = aux_data[row_local];
   	       space = row_space[row_local]; 
   	       old_size = row_length[row_local]; 
   	       size = space - old_size;
   	       if (size < n)
      	       {
      	          size = n - size;
      	          tmp_j = jxf_CTAlloc(JXF_Int, size);
      	          tmp_data = jxf_CTAlloc(JXF_Real, size);
      	       }
      	       else
      	       {
      	          tmp_j = NULL;
      	       }
      	       tmp_indx = 0;
      	       not_found = 1;
      	       size = old_size;
               for (i=0; i < n; i++)
      	       {
      	          for (j=0; j < old_size; j++)
      	          {
      	             if (local_j[j] == cols[indx])
      	             {
                        local_data[j] = values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
      	             if (size < space)
      	             {
      	                local_j[size] = cols[indx];
      	                local_data[size++] = values[indx];
      	             }
      	             else
      	             {
      	                tmp_j[tmp_indx] = cols[indx];
      	                tmp_data[tmp_indx++] = values[indx];
      	             }
      	          }
      	          not_found = 1;
        	  indx++;
      	       }
      	    
               row_length[row_local] = size+tmp_indx;
               
               if (tmp_indx)
               {
   	          aux_j[row_local] = jxf_TReAlloc(aux_j[row_local], JXF_Int, size+tmp_indx);
   	          aux_data[row_local] = jxf_TReAlloc(aux_data[row_local], JXF_Real, size+tmp_indx);
                  row_space[row_local] = size+tmp_indx;
                  local_j = aux_j[row_local];
                  local_data = aux_data[row_local];
               }
   
   	       cnt = size; 
   
   	       for (i=0; i < tmp_indx; i++)
   	       {
   	          local_j[cnt] = tmp_j[i];
   	          local_data[cnt++] = tmp_data[i];
	       }
  
	       if (tmp_j)
	       { 
	          jxf_TFree(tmp_j); 
	          jxf_TFree(tmp_data); 
	       } 
            }
            else /* insert immediately into data in ParCSRMatrix structure */
            {
	       offd_indx = jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local];
	       diag_indx = jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local];
               diag = jxf_hpCSRMatrixDiag(hp_matrix);
               diag_i = jxf_CSRMatrixI(diag);
               diag_j = jxf_CSRMatrixJ(diag);
               diag_data = jxf_CSRMatrixData(diag);
               offd = jxf_hpCSRMatrixOffd(hp_matrix);
               offd_i = jxf_CSRMatrixI(offd);
               if (num_procs > 1)
	       {
	          offd_j = jxf_CSRMatrixJ(offd);
                  offd_data = jxf_CSRMatrixData(offd);
               }
	       cnt_diag = diag_indx;
	       cnt_offd = offd_indx;
	       diag_space = diag_i[row_local+1];
	       offd_space = offd_i[row_local+1];
	       not_found = 1;
  	       for (i=0; i < n; i++)
	       {
	          if (cols[indx] < col_0 || cols[indx] > col_n)
                     /* insert into offd */	
	          {
	             for (j=offd_i[row_local]; j < offd_indx; j++)
	             {
		        if (offd_j[j] == cols[indx])
		        {
                           offd_data[j] = values[indx];
		           not_found = 0;
		           break;
		        }
	             }
	             if (not_found)
	             { 
	                if (cnt_offd < offd_space) 
	                { 
	                   offd_j[cnt_offd] = cols[indx];
	                   offd_data[cnt_offd++] = values[indx];
	                } 
	                else 
	 	        {
                           jxf_error(JXF_ERROR_GENERIC);
	    	           if (print_level)
                              jxf_printf("Error in row %d ! Too many elements!\n", row);
	    	           /* return 1; */
                           return jxf_error_flag;
	 	        }
	             }  
	             not_found = 1;
	          }
	          else  /* insert into diag */
	          {
	             for (j=diag_i[row_local]; j < diag_indx; j++)
	             {
		        if (diag_j[j] == cols[indx])
		        {
                           diag_data[j] = values[indx];
		           not_found = 0;
		           break;
		        }
	             } 
	             if (not_found)
	             { 
	                if (cnt_diag < diag_space) 
	                { 
	                   diag_j[cnt_diag] = cols[indx];
	                   diag_data[cnt_diag++] = values[indx];
	                } 
	                else 
	 	        {
                           jxf_error(JXF_ERROR_GENERIC);
	    	           if (print_level)
                              jxf_printf("Error in row %d ! Too many elements !\n", row);
	    	           /* return 1; */
                           return jxf_error_flag;
	 	        }
	             } 
	             not_found = 1;
	          }
	          indx++;
	       }

               jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local] = cnt_diag;
               jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local] = cnt_offd;

            }
         }

         /* processor does not own the row */
         else
	 {
            indx +=  n;
	    aux_matrix = jxf_IJMatrixTranslator(matrix);
   	    if (aux_matrix)
            {
   	       //current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
   	       off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
   	       off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	       off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
               col_indx = 0;
               cancel_indx = jxf_AuxParCSRMatrixCancelIndx(aux_matrix);
               for (i=0; i < off_proc_i_indx; i=i+2)
               {
	          row_len = off_proc_i[i+1];
	          if (off_proc_i[i] == row)
		  {
		     for (j=0; j < n; j++)
		     {
			cnt1 = col_indx;
			for (k=0; k < row_len; k++)
			{
			   if (off_proc_j[cnt1] == cols[j])
			   {
                              off_proc_j[cnt1++] = -1;
               		      cancel_indx++;
			      /* if no repetition allowed */
                              /* off_proc_j[col_indx] = -1;
                                 col_indx -= k;
                                 break; */
			   }
			   else
			   {
			      cnt1++;
			   }
			}
		     }
		     col_indx += row_len;
                  }
                  else
                  {
		     col_indx += row_len;
                  }
               }
               jxf_AuxParCSRMatrixCancelIndx(aux_matrix) = cancel_indx;
            }
	 }
      }
   }

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixAssemble( JXF_IJMatrix matrix )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;

   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   if ( jxf_IJMatrixObjectType(ijmatrix) == JXF_HPCSR )
   {
      return(jxf_IJMatrixAssemblehpCSR(ijmatrix));
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixAssemblehpCSR( jxf_IJMatrix *matrix )
{
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   jxf_hpCSRMatrix *hp_matrix = jxf_IJMatrixObject(matrix);
   jxf_AuxParCSRMatrix *aux_matrix = jxf_IJMatrixTranslator(matrix);
   JXF_Int *row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   JXF_Int *col_partitioning = jxf_IJMatrixColPartitioning(matrix);

   jxf_CSRMatrix *diag = jxf_hpCSRMatrixDiag(hp_matrix);
   jxf_CSRMatrix *offd = jxf_hpCSRMatrixOffd(hp_matrix);
   JXF_Int *diag_i = jxf_CSRMatrixI(diag);
   JXF_Int *offd_i = jxf_CSRMatrixI(offd);
   JXF_Int *diag_j;
   JXF_Int *offd_j = NULL;
   JXF_Real *diag_data;
   JXF_Real *offd_data = NULL;
   JXF_Int i, j, j0;
   JXF_Int num_cols_offd;
   JXF_Int *diag_pos;
   JXF_Int *col_map_offd;
   JXF_Int *row_length;
   JXF_Int **aux_j;
   JXF_Real **aux_data;
   JXF_Int my_id, num_procs;
   JXF_Int num_rows;
   JXF_Int i_diag, i_offd;
   JXF_Int col_0, col_n;
   JXF_Int nnz_offd;
   JXF_Int *aux_offd_j;
   JXF_Real temp; 
#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int base = jxf_IJMatrixGlobalFirstCol(matrix);
#else
   JXF_Int base = col_partitioning[0];
#endif
   JXF_Int off_proc_i_indx;
   JXF_Int max_off_proc_elmts;
   JXF_Int current_num_elmts;
   JXF_Int *off_proc_i;
   JXF_Int *off_proc_j;
   JXF_Real *off_proc_data;
   JXF_Int offd_proc_elmts;
   JXF_Int new_off_proc_i_indx;
   JXF_Int cancel_indx;
   JXF_Int col_indx;
   JXF_Int current_indx;
   JXF_Int current_i;
   JXF_Int row_len;
   JXF_Int max_num_threads;
   JXF_Int aux_flag, aux_flag_global;

   max_num_threads = jxf_NumThreads();

   /* first find out if anyone has an aux_matrix, and create one if you don't
    * have one, but other procs do */
   aux_flag = 0;
   aux_flag_global = 0;
   if(aux_matrix)
   {
      aux_flag = 1;
   }
   jxf_MPI_Allreduce(&aux_flag, &aux_flag_global, 1, JXF_MPI_INT, MPI_SUM, comm);
   if(aux_flag_global && (!aux_flag))
   {
      jxf_MPI_Comm_rank(comm, &my_id);
      num_rows = row_partitioning[my_id+1] - row_partitioning[my_id];
      jxf_AuxParCSRMatrixCreate(&aux_matrix, num_rows, num_rows, NULL);
      jxf_AuxParCSRMatrixNeedAux(aux_matrix) = 0;
      jxf_IJMatrixTranslator(matrix) = aux_matrix;
   }

   if (aux_matrix)
   {
      /* first delete all cancelled elements */
      cancel_indx = jxf_AuxParCSRMatrixCancelIndx(aux_matrix);
      if (cancel_indx)
      {
         current_num_elmts=jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i=jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j=jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
         off_proc_data=jxf_AuxParCSRMatrixOffProcData(aux_matrix);
         off_proc_i_indx=jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
         col_indx = 0;
         current_i = 0;
         current_indx = 0;
	 new_off_proc_i_indx = off_proc_i_indx;
         for (i=0; i < off_proc_i_indx; i= i+2)
         {
            row_len = off_proc_i[i+1];
            for (j=0; j < off_proc_i[i+1]; j++)
            {
 	       if (off_proc_j[col_indx] == -1)
	       {
		  col_indx++;
		  row_len--;
		  current_num_elmts--;
               }
               else
	       {
		  off_proc_j[current_indx] = off_proc_j[col_indx];
		  off_proc_data[current_indx++] = off_proc_data[col_indx++];
	       }
            }
            if (row_len)
	    {
	       off_proc_i[current_i] = off_proc_i[i];
	       off_proc_i[current_i+1] = row_len;
	       current_i += 2;
            }
	    else
            {
	       new_off_proc_i_indx -= 2;
            }
         }
         jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix) = new_off_proc_i_indx;
         jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix) = current_num_elmts;
      }
      off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
      jxf_MPI_Allreduce(&off_proc_i_indx, &offd_proc_elmts, 1, JXF_MPI_INT, MPI_SUM, comm);
      if (offd_proc_elmts)
      {
         max_off_proc_elmts=jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix);
         current_num_elmts=jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i=jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j=jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
         off_proc_data=jxf_AuxParCSRMatrixOffProcData(aux_matrix);
         jxf_IJMatrixAssembleOffProcValshpCSR(
            matrix,off_proc_i_indx, max_off_proc_elmts, current_num_elmts,
            off_proc_i, off_proc_j, off_proc_data);
      }
   }

   if (jxf_IJMatrixAssembleFlag(matrix) == 0)
   {
      jxf_MPI_Comm_size(comm, &num_procs); 
      jxf_MPI_Comm_rank(comm, &my_id);
#ifdef JXF_NO_GLOBAL_PARTITION
      num_rows = row_partitioning[1] - row_partitioning[0]; 
      col_0 = col_partitioning[0];
      col_n = col_partitioning[1]-1;
#else
      num_rows = row_partitioning[my_id+1] - row_partitioning[my_id]; 
      col_0 = col_partitioning[my_id];
      col_n = col_partitioning[my_id+1]-1;
#endif
      /* move data into ParCSRMatrix if not there already */ 
      if (jxf_AuxParCSRMatrixNeedAux(aux_matrix))
      {
         JXF_Int *diag_array, *offd_array;
         diag_array = jxf_CTAlloc(JXF_Int, max_num_threads);
         offd_array = jxf_CTAlloc(JXF_Int, max_num_threads);
         aux_j = jxf_AuxParCSRMatrixAuxJ(aux_matrix);
         aux_data = jxf_AuxParCSRMatrixAuxData(aux_matrix);
         row_length = jxf_AuxParCSRMatrixRowLength(aux_matrix);
         diag_pos = jxf_CTAlloc(JXF_Int, num_rows);
         i_diag = 0;
         i_offd = 0;
#if JXF_USING_OPENMP
#pragma omp parallel private(i,j,i_diag,i_offd)
#endif
         {
          JXF_Int *local_j;
          JXF_Real *local_data;
          JXF_Int rest, size, ns, ne;
          JXF_Int num_threads, my_thread_num;
          num_threads = jxf_NumActiveThreads();
          my_thread_num = jxf_GetThreadNum();

          size = num_rows/num_threads;
          rest = num_rows - size*num_threads;
        
          if (my_thread_num < rest)
          {
             ns = my_thread_num*(size + 1); 
             ne = (my_thread_num+1)*(size + 1); 
          } 
          else
          {
             ns = my_thread_num*size + rest; 
             ne = (my_thread_num+1)*size + rest; 
          } 

          i_diag = 0;
          i_offd = 0;
          for (i=ns; i < ne; i++)
          {
   	    local_j = aux_j[i];
   	    local_data = aux_data[i];
   	    diag_pos[i] = -1;
   	    for (j=0; j < row_length[i]; j++)
   	    {
   	       if (local_j[j] < col_0 || local_j[j] > col_n)
   	          i_offd++;
   	       else
   	       {
   	          i_diag++;
   	          if (local_j[j]-col_0 == i) diag_pos[i] = j;
   	       }
   	    }
          }
   	  diag_array[my_thread_num] = i_diag;
   	  offd_array[my_thread_num] = i_offd;
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
          if (my_thread_num == 0)
          {
            i_diag = 0;
            i_offd = 0;
            for (i = 0; i < num_threads; i++)
	    {
	       i_diag += diag_array[i];
	       i_offd += offd_array[i];
	       diag_array[i] = i_diag;
	       offd_array[i] = i_offd;
	    }
            diag_i[num_rows] = i_diag;         
            offd_i[num_rows] = i_offd;         
            if (jxf_CSRMatrixJ(diag))
               jxf_TFree(jxf_CSRMatrixJ(diag));
            if (jxf_CSRMatrixData(diag))
               jxf_TFree(jxf_CSRMatrixData(diag));
            if (jxf_CSRMatrixJ(offd))
               jxf_TFree(jxf_CSRMatrixJ(offd));
            if (jxf_CSRMatrixData(offd))
               jxf_TFree(jxf_CSRMatrixData(offd));
            diag_j = jxf_CTAlloc(JXF_Int, i_diag);
            diag_data = jxf_CTAlloc(JXF_Real, i_diag);
            if (i_offd > 0)
            {
    	       offd_j = jxf_CTAlloc(JXF_Int,i_offd);
               offd_data = jxf_CTAlloc(JXF_Real,i_offd);
            }
          }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
          if (my_thread_num)
          {
             i_diag = diag_array[my_thread_num-1];
             i_offd = offd_array[my_thread_num-1];
          }
          else
          {
             i_diag = 0;
             i_offd = 0;
          }
          for (i=ns; i < ne; i++)
          {
   	    diag_i[i] = i_diag;
   	    offd_i[i] = i_offd;
   	    local_j = aux_j[i];
   	    local_data = aux_data[i];
            if (diag_pos[i] > -1)
            {
   	       diag_j[i_diag] = local_j[diag_pos[i]] - col_0;
               diag_data[i_diag++] = local_data[diag_pos[i]];
            }
   	    for (j=0; j < row_length[i]; j++)
   	    {
   	       if (local_j[j] < col_0 || local_j[j] > col_n)
   	       {
   	          offd_j[i_offd] = local_j[j];
   	          offd_data[i_offd++] = local_data[j];
   	       }
   	       else if (j != diag_pos[i])
   	       {
   	          diag_j[i_diag] = local_j[j] - col_0;
   	          diag_data[i_diag++] = local_data[j];
   	       }
   	    }
          }
         } /* end parallel region */

         jxf_TFree(diag_array);
         jxf_TFree(offd_array);

         jxf_CSRMatrixJ(diag) = diag_j;      
         jxf_CSRMatrixData(diag) = diag_data;      
         jxf_CSRMatrixNumNonzeros(diag) = diag_i[num_rows];      
         if (offd_i[num_rows] > 0)
         {
            jxf_CSRMatrixJ(offd) = offd_j;      
            jxf_CSRMatrixData(offd) = offd_data;      
         }
         jxf_CSRMatrixNumNonzeros(offd) = offd_i[num_rows];      
         jxf_TFree(diag_pos);
      }
      else
      {
         /* move diagonal element into first space */
         diag_j = jxf_CSRMatrixJ(diag);
         diag_data = jxf_CSRMatrixData(diag);
#if JXF_USING_OPENMP
#pragma omp parallel for private (i,j,j0,temp)
#endif
         for (i=0; i < num_rows; i++)
         {
   	    j0 = diag_i[i];
   	    for (j=j0; j < diag_i[i+1]; j++)
   	    {
   	       diag_j[j] -= col_0;
   	       if (diag_j[j] == i)
   	       {
   	          temp = diag_data[j0];
   	          diag_data[j0] = diag_data[j];
   	          diag_data[j] = temp;
   	          diag_j[j] = diag_j[j0];
   	          diag_j[j0] = i;
   	       }
   	    }
         }
         offd_j = jxf_CSRMatrixJ(offd);
      }

      /*  generate the nonzero rows inside offd and diag by calling */

      jxf_CSRMatrixSetRownnz(diag);
      jxf_CSRMatrixSetRownnz(offd);

      /*  generate col_map_offd */
      nnz_offd = offd_i[num_rows];
      if (nnz_offd)
      {
         aux_offd_j = jxf_CTAlloc(JXF_Int, nnz_offd);
         for (i=0; i < nnz_offd; i++)
            aux_offd_j[i] = offd_j[i];
         jxf_qsort0(aux_offd_j,0,nnz_offd-1);
         num_cols_offd = 1;
         for (i=0; i < nnz_offd-1; i++)
         {
            if (aux_offd_j[i+1] > aux_offd_j[i])
               aux_offd_j[num_cols_offd++] = aux_offd_j[i+1];
         }
         col_map_offd = jxf_CTAlloc(JXF_Int,num_cols_offd);
         for (i=0; i < num_cols_offd; i++)
   	    col_map_offd[i] = aux_offd_j[i];
#if JXF_USING_OPENMP
#pragma omp parallel for private(i)
#endif
         for (i=0; i < nnz_offd; i++)
            offd_j[i] = jxf_BinarySearch(col_map_offd,offd_j[i],num_cols_offd);

 	 if (base)
 	 {
	    for (i=0; i < num_cols_offd; i++)
	       col_map_offd[i] -= base;
	 } 
         jxf_hpCSRMatrixColMapOffd(hp_matrix) = col_map_offd;    
         jxf_CSRMatrixNumCols(offd) = num_cols_offd;    
         jxf_TFree(aux_offd_j);
      }
      jxf_AuxParCSRMatrixDestroy(aux_matrix);
      jxf_IJMatrixTranslator(matrix) = NULL;
      jxf_IJMatrixAssembleFlag(matrix) = 1;
   }
   #ifdef USING_HWLOC
   hp_matrix = jxf_hpCreateMatrixLevelBlock(hp_matrix);
   #endif

   return jxf_error_flag;
}

#ifndef JXF_NO_GLOBAL_PARTITION

JXF_Int
jxf_IJMatrixAssembleOffProcValshpCSR( jxf_IJMatrix *matrix,
                                      JXF_Int          off_proc_i_indx,
                                      JXF_Int          max_off_proc_elmts,
                                      JXF_Int          current_num_elmts,
                                      JXF_Int         *off_proc_i,
                                      JXF_Int         *off_proc_j,
                                      JXF_Real      *off_proc_data )
{
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   MPI_Request *requests = NULL;
   MPI_Status *status = NULL;
   JXF_Int i, ii, j, j2, jj, n, row;
   JXF_Int iii, iid, indx, ip;
   JXF_Int proc_id, num_procs, my_id;
   JXF_Int num_sends, num_sends3;
   JXF_Int num_recvs;
   JXF_Int num_requests;
   JXF_Int vec_start, vec_len;
   JXF_Int *send_procs;
   JXF_Int *chunks;
   JXF_Int *send_i;
   JXF_Int *send_map_starts;
   JXF_Int *dbl_send_map_starts;
   JXF_Int *recv_procs;
   JXF_Int *recv_chunks;
   JXF_Int *recv_i;
   JXF_Int *recv_vec_starts;
   JXF_Int *dbl_recv_vec_starts;
   JXF_Int *info;
   JXF_Int *int_buffer;
   JXF_Int *proc_id_mem;
   JXF_Int *partitioning;
   JXF_Int *displs;
   JXF_Int *recv_buf;
   JXF_Real *send_data;
   JXF_Real *recv_data;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   partitioning = jxf_IJMatrixRowPartitioning(matrix);

   info = jxf_CTAlloc(JXF_Int,num_procs);
   chunks = jxf_CTAlloc(JXF_Int,num_procs);
   proc_id_mem = jxf_CTAlloc(JXF_Int, off_proc_i_indx/2);
   j=0;
   for (i=0; i < off_proc_i_indx; i++)
   {
      row = off_proc_i[i++];
      if (row < 0) row = -row-1; 
      n = off_proc_i[i];
      proc_id = jxf_FindProc(partitioning,row,num_procs);
      proc_id_mem[j++] = proc_id; 
      info[proc_id] += n;
      chunks[proc_id]++;
   }

   /* determine send_procs and amount of data to be sent */
   num_sends = 0;
   for (i=0; i < num_procs; i++)
   {
      if (info[i])
      {
         num_sends++;
      }
   }
   send_procs =  jxf_CTAlloc(JXF_Int,num_sends);
   send_map_starts =  jxf_CTAlloc(JXF_Int,num_sends+1);
   dbl_send_map_starts =  jxf_CTAlloc(JXF_Int,num_sends+1);
   num_sends3 = 3*num_sends;
   int_buffer =  jxf_CTAlloc(JXF_Int,3*num_sends);
   j = 0;
   j2 = 0;
   send_map_starts[0] = 0;
   dbl_send_map_starts[0] = 0;
   for (i=0; i < num_procs; i++)
   {
      if (info[i])
      {
         send_procs[j++] = i;
         send_map_starts[j] = send_map_starts[j-1]+2*chunks[i]+info[i];
         dbl_send_map_starts[j] = dbl_send_map_starts[j-1]+info[i];
         int_buffer[j2++] = i;
	 int_buffer[j2++] = chunks[i];
	 int_buffer[j2++] = info[i];
      }
   }

   jxf_TFree(chunks);

   jxf_MPI_Allgather(&num_sends3,1,JXF_MPI_INT,info,1,JXF_MPI_INT,comm);

   displs = jxf_CTAlloc(JXF_Int, num_procs+1);
   displs[0] = 0;
   for (i=1; i < num_procs+1; i++)
      displs[i] = displs[i-1]+info[i-1];
   recv_buf = jxf_CTAlloc(JXF_Int, displs[num_procs]);

   jxf_MPI_Allgatherv(int_buffer,num_sends3,JXF_MPI_INT,recv_buf,info,displs,JXF_MPI_INT,comm);

   jxf_TFree(int_buffer);
   jxf_TFree(info);

   /* determine recv procs and amount of data to be received */
   num_recvs = 0;
   for (j=0; j < displs[num_procs]; j+=3)
   {
      if (recv_buf[j] == my_id)
	 num_recvs++;
   }

   recv_procs = jxf_CTAlloc(JXF_Int,num_recvs);
   recv_chunks = jxf_CTAlloc(JXF_Int,num_recvs);
   recv_vec_starts = jxf_CTAlloc(JXF_Int,num_recvs+1);
   dbl_recv_vec_starts = jxf_CTAlloc(JXF_Int,num_recvs+1);

   j2 = 0;
   recv_vec_starts[0] = 0;
   dbl_recv_vec_starts[0] = 0;
   for (i=0; i < num_procs; i++)
   {
      for (j=displs[i]; j < displs[i+1]; j+=3)
      {
         if (recv_buf[j] == my_id)
         {
	    recv_procs[j2] = i;
	    recv_chunks[j2++] = recv_buf[j+1];
	    recv_vec_starts[j2] = recv_vec_starts[j2-1]+2*recv_buf[j+1]+recv_buf[j+2];
	    dbl_recv_vec_starts[j2] = dbl_recv_vec_starts[j2-1]+recv_buf[j+2];
         }
         if (j2 == num_recvs) break;
      }
   }
   jxf_TFree(recv_buf);
   jxf_TFree(displs);

   /* set up data to be sent to send procs */
   /* send_i contains for each send proc : row no., no. of elmts and column
      indices, send_data contains corresponding values */
      
   send_i = jxf_CTAlloc(JXF_Int,send_map_starts[num_sends]);
   send_data = jxf_CTAlloc(JXF_Real,dbl_send_map_starts[num_sends]);
   recv_i = jxf_CTAlloc(JXF_Int,recv_vec_starts[num_recvs]);
   recv_data = jxf_CTAlloc(JXF_Real,dbl_recv_vec_starts[num_recvs]);
    
   j=0;
   jj=0;
   for (i=0; i < off_proc_i_indx; i++)
   {
      row = off_proc_i[i++]; 
      n = off_proc_i[i];
      proc_id = proc_id_mem[i/2];
      indx = jxf_BinarySearch(send_procs,proc_id,num_sends);
      iii = send_map_starts[indx];
      iid = dbl_send_map_starts[indx];
      send_i[iii++] = row;
      send_i[iii++] = n;
      for (ii = 0; ii < n; ii++)
      {
         send_i[iii++] = off_proc_j[jj];
         send_data[iid++] = off_proc_data[jj++];
      }
      send_map_starts[indx] = iii;
      dbl_send_map_starts[indx] = iid;
   }

   jxf_TFree(proc_id_mem);

   for (i=num_sends; i > 0; i--)
   {
      send_map_starts[i] = send_map_starts[i-1];
      dbl_send_map_starts[i] = dbl_send_map_starts[i-1];
   }
   send_map_starts[0] = 0;
   dbl_send_map_starts[0] = 0;

   num_requests = num_recvs+num_sends;

   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   j=0; 
   for (i=0; i < num_recvs; i++)
   {
      vec_start = recv_vec_starts[i];
      vec_len = recv_vec_starts[i+1] - vec_start;
      ip = recv_procs[i];
      jxf_MPI_Irecv(&recv_i[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
   }

   for (i=0; i < num_sends; i++)
   {
      vec_start = send_map_starts[i];
      vec_len = send_map_starts[i+1] - vec_start;
      ip = send_procs[i];
      jxf_MPI_Isend(&send_i[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
   }

   if (num_requests)
   {
      jxf_MPI_Waitall(num_requests, requests, status);
   }

   j=0;
   for (i=0; i < num_recvs; i++)
   {
      vec_start = dbl_recv_vec_starts[i];
      vec_len = dbl_recv_vec_starts[i+1] - vec_start;
      ip = recv_procs[i];
      jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
   }

   for (i=0; i < num_sends; i++)
   {
      vec_start = dbl_send_map_starts[i];
      vec_len = dbl_send_map_starts[i+1] - vec_start;
      ip = send_procs[i];
      jxf_MPI_Isend(&send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
   }

   if (num_requests)
   {
      jxf_MPI_Waitall(num_requests, requests, status);
   }

   jxf_TFree(requests);
   jxf_TFree(status);
   jxf_TFree(send_i);
   jxf_TFree(send_data);
   jxf_TFree(send_procs);
   jxf_TFree(send_map_starts);
   jxf_TFree(dbl_send_map_starts);
   jxf_TFree(recv_procs);
   jxf_TFree(recv_vec_starts);
   jxf_TFree(dbl_recv_vec_starts);

   j = 0;
   j2 = 0;
   for (i=0; i < num_recvs; i++)
   {
      for (ii=0; ii < recv_chunks[i]; ii++)
      {
         row = recv_i[j];
 	 jxf_IJMatrixAddToValueshpCSR(matrix,1,&recv_i[j+1],&row,&recv_i[j+2],&recv_data[j2]);
	 j2 += recv_i[j+1]; 
	 j += recv_i[j+1]+2; 
      }
   }
   jxf_TFree(recv_chunks);
   jxf_TFree(recv_i);
   jxf_TFree(recv_data);

   return jxf_error_flag;
}

#else

/* assumed partition version */

JXF_Int
jxf_IJMatrixAssembleOffProcValshpCSR( jxf_IJMatrix *matrix,
                                      JXF_Int          off_proc_i_indx,
                                      JXF_Int          max_off_proc_elmts,
                                      JXF_Int          current_num_elmts,
                                      JXF_Int         *off_proc_i,
                                      JXF_Int         *off_proc_j,
                                      JXF_Real      *off_proc_data )
{
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int i, j, k, in_i;
   JXF_Int myid;
   JXF_Int proc_id, last_proc, prev_id, tmp_id;
   JXF_Int max_response_size;
   JXF_Int global_num_cols;
   JXF_Int global_first_col;
   JXF_Int global_first_row;
   JXF_Int ex_num_contacts = 0, num_rows = 0;
   JXF_Int range_start, range_end;
   JXF_Int num_elements;
   JXF_Int storage;
   JXF_Int indx;
   JXF_Int row, num_ranges;
   JXF_Int num_recvs;
   JXF_Int counter, upper_bound;
   JXF_Int num_real_procs;
   JXF_Int current_proc, original_proc_indx;
   JXF_Int *row_list=NULL, *row_list_num_elements=NULL;
   JXF_Int *a_proc_id=NULL, *orig_order=NULL;
   JXF_Int *real_proc_id = NULL, *us_real_proc_id = NULL;
   JXF_Int *ex_contact_procs = NULL, *ex_contact_vec_starts = NULL, *ex_contact_buf = NULL;
   JXF_Int *recv_starts=NULL;
   JXF_Int *response_buf = NULL, *response_buf_starts=NULL;
   JXF_Int *num_rows_per_proc = NULL, *num_elements_total = NULL;
   JXF_Int *argsort_contact_procs = NULL;
   JXF_Int  obj_size_bytes, int_size, complex_size;
   JXF_Int  tmp_int;
   JXF_Int *col_ptr;
   JXF_Int *int_data = NULL;
   JXF_Int int_data_size = 0, complex_data_size = 0;
   void *void_contact_buf = NULL;
   void *index_ptr;
   void *recv_data_ptr;
   JXF_Real  tmp_complex;
   JXF_Real *col_data_ptr;
   JXF_Real *complex_data = NULL;
   jxf_DataExchangeResponse  response_obj1, response_obj2;
   jxf_ProcListElements      send_proc_obj;
   jxf_IJAssumedPart   *apart;

   jxf_MPI_Comm_rank(comm, &myid);
   global_num_cols = jxf_IJMatrixGlobalNumCols(matrix);
   global_first_col = jxf_IJMatrixGlobalFirstCol(matrix);
   global_first_row = jxf_IJMatrixGlobalFirstRow(matrix);

   num_rows = off_proc_i_indx/2;
   
   /* verify that we have created the assumed partition */
   if  (jxf_IJMatrixAssumedPart(matrix) == NULL)
   {
      jxf_IJMatrixCreateAssumedPartition(matrix);
   }

   apart = jxf_IJMatrixAssumedPart(matrix);

   /*if  (jxf_ParCSRMatrixAssumedPartition(par_matrix) == NULL)
   {
      jxf_ParCSRMatrixCreateAssumedPartition(par_matrix);
   }

   apart = jxf_ParCSRMatrixAssumedPartition(par_matrix);*/

   row_list = jxf_CTAlloc(JXF_Int, num_rows);
   row_list_num_elements = jxf_CTAlloc(JXF_Int, num_rows);
   a_proc_id = jxf_CTAlloc(JXF_Int, num_rows);
   orig_order =  jxf_CTAlloc(JXF_Int, num_rows);
   real_proc_id = jxf_CTAlloc(JXF_Int, num_rows);

   /* get the assumed processor id for each row */
   if (num_rows > 0 )
   {
      for (i=0; i < num_rows; i++)
      {
         row = off_proc_i[i*2];
         if (row < 0) row = -row-1; 
         row_list[i] = row;
         row_list_num_elements[i] = off_proc_i[i*2+1];
         jxf_GetAssumedPartitionProcFromRow(comm, row, global_first_row, global_num_cols, &proc_id);
         a_proc_id[i] = proc_id;
         orig_order[i] = i;
      }

      /* now we need to find the actual order of each row  - sort on row -
         this will result in proc ids sorted also...*/

      jxf_qsort3i(row_list, a_proc_id, orig_order, 0, num_rows -1);

      /* calculate the number of contacts */
      ex_num_contacts = 1;
      last_proc = a_proc_id[0];
      for (i=1; i < num_rows; i++)
      {
         if (a_proc_id[i] > last_proc)      
         {
            ex_num_contacts++;
            last_proc = a_proc_id[i];
         }
      }
      
   } 
   
   /* now we will go through a create a contact list - need to contact assumed
      processors and find out who the actual row owner is - we will contact with
      a range (2 numbers) */

   ex_contact_procs = jxf_CTAlloc(JXF_Int, ex_num_contacts);
   ex_contact_vec_starts =  jxf_CTAlloc(JXF_Int, ex_num_contacts+1);
   ex_contact_buf =  jxf_CTAlloc(JXF_Int, ex_num_contacts*2);

   counter = 0;
   range_end = -1;
   for (i=0; i< num_rows; i++) 
   {
      if (row_list[i] > range_end)
      {
         /* assumed proc */
         proc_id = a_proc_id[i];

         /* end of prev. range */
         if (counter > 0)  ex_contact_buf[counter*2 - 1] = row_list[i-1];
         
         /*start new range*/
    	 ex_contact_procs[counter] = proc_id;
         ex_contact_vec_starts[counter] = counter*2;
         ex_contact_buf[counter*2] =  row_list[i];
         counter++;
         
         jxf_GetAssumedPartitionRowRange(comm, proc_id, global_first_col, global_num_cols, 
                                           &range_start, &range_end); 
      }
   }
   /*finish the starts*/
   ex_contact_vec_starts[counter] =  counter*2;
   /*finish the last range*/
   if (counter > 0)  
      ex_contact_buf[counter*2 - 1] = row_list[num_rows - 1];

   /*don't allocate space for responses */
    
   /* create response object - can use same fill response as used in the commpkg
      routine */
   response_obj1.fill_response = jxf_RangeFillResponseIJDetermineRecvProcs;
   response_obj1.data1 =  apart; /* this is necessary so we can fill responses*/ 
   response_obj1.data2 = NULL;
   
   max_response_size = 6;  /* 6 means we can fit 3 ranges*/
   
   jxf_DataExchangeList(ex_num_contacts, ex_contact_procs, 
                          ex_contact_buf, ex_contact_vec_starts, sizeof(JXF_Int), 
                          sizeof(JXF_Int), &response_obj1, max_response_size, 1, 
                          comm, (void**) &response_buf, &response_buf_starts);

   /* now response_buf contains a proc_id followed by a range upper bound */

   jxf_TFree(ex_contact_procs);
   jxf_TFree(ex_contact_buf);
   jxf_TFree(ex_contact_vec_starts);

   jxf_TFree(a_proc_id);

   /*how many ranges were returned?*/
   num_ranges = response_buf_starts[ex_num_contacts];   
   num_ranges = num_ranges/2;

   prev_id = -1;
   j = 0;
   counter = 0;
   num_real_procs = 0;
  
   /* loop through ranges - create a list of actual processor ids*/
   for (i=0; i<num_ranges; i++)
   {
      upper_bound = response_buf[i*2+1];
      counter = 0;
      tmp_id = response_buf[i*2];

      /* loop through row_list entries - counting how many are in the range */
      while (j < num_rows && row_list[j] <= upper_bound)    
      {
         real_proc_id[j] = tmp_id;
         j++;
         counter++;       
      }
      if (counter > 0 && tmp_id != prev_id)        
      {
         num_real_procs++;
      }

      prev_id = tmp_id;
   }

   /* now we have the list of real procesors ids (real_proc_id) - and the number
      of distinct ones - so now we can set up data to be sent - we have
      JXF_Int data and JXF_Real data.  that we will need to pack
      together */
   
   /* first find out how many rows and elements we need to send per proc - so we
      can do storage */
   
   ex_contact_procs = jxf_CTAlloc(JXF_Int, num_real_procs);
   num_rows_per_proc = jxf_CTAlloc(JXF_Int, num_real_procs);
   num_elements_total  =  jxf_CTAlloc(JXF_Int, num_real_procs); 
   
   counter = 0;
   
   if (num_real_procs > 0 )
   {
      ex_contact_procs[0] = real_proc_id[0];
      num_rows_per_proc[0] = 1;
      num_elements_total[0] = row_list_num_elements[orig_order[0]];

      /* loop through real procs - these are sorted (row_list is sorted also)*/
      for (i=1; i < num_rows; i++)
      {
         if (real_proc_id[i] == ex_contact_procs[counter]) /* same processor */
         {
            num_rows_per_proc[counter] += 1; /*another row */
            num_elements_total[counter] += row_list_num_elements[orig_order[i]];
         }
         else /* new processor */
         {
            counter++;
            ex_contact_procs[counter] = real_proc_id[i];
            num_rows_per_proc[counter] = 1;
            num_elements_total[counter] = row_list_num_elements[orig_order[i]];
         }
      }
   }
   
   /* to pack together, we need to use the largest obj. size of
      (JXF_Int) and (JXF_Real) - if these are much different, then we are
      wasting some storage, but I do not think that it will be a
      large amount since this function should not be used on really
      large amounts of data anyway*/
   int_size = sizeof(JXF_Int);
   complex_size = sizeof(JXF_Real);
   
   obj_size_bytes = jxf_max(int_size, complex_size);

   /* set up data to be sent to send procs */
   /* for each proc, ex_contact_buf contains #rows, row #,
      no. elements, col indicies, col data, row #, no. elements, col
      indicies, col data, etc. */
      
   /* first calculate total storage and make vec_starts arrays */
   storage = 0;
   ex_contact_vec_starts = jxf_CTAlloc(JXF_Int, num_real_procs + 1);
   ex_contact_vec_starts[0] = -1;
   
   for (i=0; i < num_real_procs; i++)
   {
      storage += 1 + 2 * num_rows_per_proc[i] + 2* num_elements_total[i];
      ex_contact_vec_starts[i+1] = -storage-1; /* need negative for next loop */
   }      

   jxf_TFree(num_elements_total);

   void_contact_buf = jxf_MAlloc(storage*obj_size_bytes);
   index_ptr = void_contact_buf; /* step through with this index */

   /* for each proc: #rows, row #, no. elements, 
      col indicies, col data, row #, no. elements, col indicies, col data, etc. */
      
   /* un-sort real_proc_id - we want to access data arrays in order, so 
      cheaper to do this*/
   us_real_proc_id =  jxf_CTAlloc(JXF_Int, num_rows);
   for (i=0; i < num_rows; i++)
   {
      us_real_proc_id[orig_order[i]] = real_proc_id[i];
   }
   jxf_TFree(real_proc_id);

   counter = 0; /* index into data arrays */
   prev_id = -1;
   for (i=0; i < num_rows; i++)
   {
      proc_id = us_real_proc_id[i];
      /* can't use row list[i] - you loose the negative signs that differentiate
         add/set values */
      row = off_proc_i[i*2];
      num_elements = row_list_num_elements[i];
      /* find position of this processor */
      indx = jxf_BinarySearch(ex_contact_procs, proc_id, num_real_procs);
      in_i = ex_contact_vec_starts[indx];

      index_ptr = (void *) ((char *) void_contact_buf + in_i*obj_size_bytes);

      /* first time for this processor - add the number of rows to the buffer */
      if (in_i < 0)
      {
         in_i = -in_i - 1;
         /* re-calc. index_ptr since in_i was negative */
         index_ptr = (void *) ((char *) void_contact_buf + in_i*obj_size_bytes);

         tmp_int =  num_rows_per_proc[indx];
         memcpy( index_ptr, &tmp_int, int_size);
         index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);

         in_i++;
      }
      /* add row # */   
      memcpy( index_ptr, &row, int_size);
      index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
      in_i++;
            
      /* add number of elements */   
      memcpy( index_ptr, &num_elements, int_size);
      index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
      in_i++;

      /* now add col indices  */
      for (j=0; j< num_elements; j++)
      {
         tmp_int = off_proc_j[counter+j]; /* col number */

         memcpy( index_ptr, &tmp_int, int_size);
         index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
         in_i ++;
      }

      /* now add data */
      for (j=0; j< num_elements; j++)
      {
         tmp_complex = off_proc_data[counter++]; /* value */

         memcpy(index_ptr, &tmp_complex, complex_size);
         index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
         in_i++;
      }

      /* increment the indexes to keep track of where we are - we
       * adjust below to be actual starts*/
      ex_contact_vec_starts[indx] = in_i;
   }
   
   /* some clean up */
 
   jxf_TFree(response_buf);
   jxf_TFree(response_buf_starts);

   jxf_TFree(us_real_proc_id);
   jxf_TFree(orig_order);
   jxf_TFree(row_list);
   jxf_TFree(row_list_num_elements);
   jxf_TFree(num_rows_per_proc);
   
   for (i=num_real_procs; i > 0; i--)
   {
      ex_contact_vec_starts[i] = ex_contact_vec_starts[i-1];
   }

   ex_contact_vec_starts[0] = 0;

   /* now send the data */

   /***********************************/

   /* first get the interger info in send_proc_obj */

   /* the response we expect is just a confirmation*/
   response_buf = NULL;
   response_buf_starts = NULL;

   /*build the response object*/

   /* use the send_proc_obj for the info kept from contacts */
   /*estimate inital storage allocation */
   send_proc_obj.length = 0;
   send_proc_obj.storage_length = num_real_procs + 5;
   send_proc_obj.id = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length+1);
   send_proc_obj.vec_starts = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length+1);
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = storage + 20;
   send_proc_obj.v_elements = jxf_MAlloc(obj_size_bytes*send_proc_obj.element_storage_length);

   response_obj2.fill_response = jxf_FillResponseIJOffProcVals;
   response_obj2.data1 = NULL;
   response_obj2.data2 = &send_proc_obj;

   max_response_size = 0;

   jxf_DataExchangeList(num_real_procs, ex_contact_procs, 
                          void_contact_buf, ex_contact_vec_starts, obj_size_bytes,
                          0, &response_obj2, max_response_size, 2, 
                          comm,  (void **) &response_buf, &response_buf_starts);

   jxf_TFree(response_buf);
   jxf_TFree(response_buf_starts);

   jxf_TFree(ex_contact_procs);
   jxf_TFree(void_contact_buf);
   jxf_TFree(ex_contact_vec_starts);

   /* Now we can unpack the send_proc_objects and call set 
      and add to values functions.  We unpack messages in a 
      deterministic order, using processor rank */
   
   num_recvs = send_proc_obj.length; 
   argsort_contact_procs = jxf_CTAlloc(JXF_Int, num_recvs);
   for(i=0; i < num_recvs; i++)
   {
      argsort_contact_procs[i] = i;
   }
   /* This sort's the id array, but the original indices are stored in
    * argsort_contact_procs */
   jxf_qsort2i(send_proc_obj.id, argsort_contact_procs, 0, num_recvs-1);

   /* alias */
   recv_data_ptr = send_proc_obj.v_elements;
   recv_starts = send_proc_obj.vec_starts;

   for (i=0; i < num_recvs; i++)
   {
      
      /* Find the current processor in order, and reset recv_data_ptr to that processor's message */
      original_proc_indx = argsort_contact_procs[i];
      current_proc = send_proc_obj.id[i];
      indx = recv_starts[original_proc_indx];
      recv_data_ptr = (void *) ((char *) send_proc_obj.v_elements + indx*obj_size_bytes);

      /* get the number of rows for this recv */
      memcpy( &num_rows, recv_data_ptr, int_size);
      recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
      indx++;

      for (j=0; j < num_rows; j++) /* for each row: unpack info */
      {
         /* row # */
         memcpy( &row, recv_data_ptr, int_size);
         recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
         indx++;

         /* num elements for this row */
         memcpy( &num_elements, recv_data_ptr, int_size);
         recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
         indx++;

         /* col indices */
         if (int_size == obj_size_bytes)
         {
            col_ptr = (JXF_Int *) recv_data_ptr;
            recv_data_ptr =
               (void *) ((char *)recv_data_ptr + num_elements*obj_size_bytes);
         }
         else /* copy data */
         {
            if (int_data_size < num_elements)
            {
               int_data = jxf_TReAlloc(int_data, JXF_Int, num_elements + 10);
            }
            for (k=0; k< num_elements; k++)
            { 
               memcpy( &int_data[k], recv_data_ptr, int_size);
               recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
            }
            col_ptr = int_data;
         }
         
         /* col data */
         if (complex_size == obj_size_bytes)
         {
            col_data_ptr = (JXF_Real *) recv_data_ptr;
            recv_data_ptr =
               (void *) ((char *)recv_data_ptr + num_elements*obj_size_bytes);
         }
         else /* copy data */
         {
            if (complex_data_size < num_elements)
            {
               complex_data =
                  jxf_TReAlloc(complex_data, JXF_Real, num_elements + 10);
            }
            for (k=0; k< num_elements; k++)
            { 
               memcpy( &complex_data[k], recv_data_ptr, complex_size);
               recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
            }
            col_data_ptr = complex_data;
         }

 	 jxf_IJMatrixAddToValueshpCSR(matrix,1,&num_elements,&row,col_ptr,col_data_ptr);
         indx += (num_elements*2); 

      }
   }
   jxf_TFree(send_proc_obj.v_elements);
   jxf_TFree(send_proc_obj.vec_starts);
   jxf_TFree(send_proc_obj.id);
   jxf_TFree(argsort_contact_procs);
 
   if (int_data) jxf_TFree(int_data);
   if (complex_data) jxf_TFree(complex_data);

   return jxf_error_flag;
}

#endif

JXF_Int
jxf_IJMatrixCreateAssumedPartition( jxf_IJMatrix *matrix )
{
   JXF_Int global_num_rows;
   JXF_Int global_first_row;
   JXF_Int myid;
   JXF_Int  row_start = 0, row_end = 0;
   JXF_Int *row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   MPI_Comm comm;
   jxf_IJAssumedPart *apart;

   global_num_rows = jxf_IJMatrixGlobalNumRows(matrix);
   global_first_row = jxf_IJMatrixGlobalFirstRow(matrix);
   comm = jxf_IJMatrixComm(matrix);

   /* find out my actual range of rows and rowumns */
   row_start = row_partitioning[0];
   row_end = row_partitioning[1]-1;
   jxf_MPI_Comm_rank(comm, &myid);

   /* allocate space */
   apart = jxf_CTAlloc(jxf_IJAssumedPart, 1);

  /* get my assumed partitioning  - we want row partitioning of the matrix
      for off processor values - so we use the row start and end 
     Note that this is different from the assumed partitioning for the parcsr matrix
     which needs it for matvec multiplications and therefore needs to do it for
     the col partitioning */
   jxf_GetAssumedPartitionRowRange( comm, myid, global_first_row, 
			global_num_rows, &(apart->row_start), &(apart->row_end));

  /*allocate some space for the partition of the assumed partition */
   apart->length = 0;

  /*room for 10 owners of the assumed partition*/ 
   apart->storage_length = 10; /*need to be >=1 */ 
   apart->proc_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_start_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_end_list = jxf_TAlloc(JXF_Int, apart->storage_length);

  /* now we want to reconcile our actual partition with the assumed partition */
   jxf_LocateAssummedPartition(comm, row_start, row_end, global_first_row,
			global_num_rows, apart, myid);

  /* this partition will be saved in the matrix data structure until the matrix is destroyed */
   jxf_IJMatrixAssumedPart(matrix) = apart;

   return jxf_error_flag;
}

JXF_Int
jxf_FillResponseIJOffProcVals( void      *p_recv_contact_buf, 
                              JXF_Int        contact_size,
                              JXF_Int        contact_proc,
                              void      *ro,
                              MPI_Comm   comm,
                              void     **p_send_response_buf,
                              JXF_Int       *response_message_size )
{
   JXF_Int    myid;
   JXF_Int    index, count, elength;
   JXF_Int object_size;
   void *index_ptr;

   jxf_DataExchangeResponse *response_obj = ro;

   jxf_ProcListElements *send_proc_obj = response_obj->data2;

   object_size = jxf_max(sizeof(JXF_Int), sizeof(JXF_Real));

   jxf_MPI_Comm_rank(comm, &myid);

   /*check to see if we need to allocate more space in send_proc_obj for vec starts
    * and id */
   if (send_proc_obj->length == send_proc_obj->storage_length)
   {
      send_proc_obj->storage_length += 20; /*add space for 20 more contact*/
      send_proc_obj->vec_starts = jxf_TReAlloc(send_proc_obj->vec_starts,JXF_Int,
                                                 send_proc_obj->storage_length + 1);
      if(send_proc_obj->id != NULL)
      {
         send_proc_obj->id = jxf_TReAlloc(send_proc_obj->id, JXF_Int,
                                         send_proc_obj->storage_length + 1);
      }
   }

   /*initialize*/ 
   count = send_proc_obj->length;
   index = send_proc_obj->vec_starts[count]; /* current number of elements */
   if( send_proc_obj->id != NULL)
   {
      send_proc_obj->id[count] = contact_proc;
   }

   /*do we need more storage for the elements?*/
   if (send_proc_obj->element_storage_length < index + contact_size)
   {
      elength = jxf_max(contact_size, 100);
      elength += index;
      send_proc_obj->v_elements = jxf_ReAlloc(send_proc_obj->v_elements, elength*object_size);
      send_proc_obj->element_storage_length = elength;
   }
   /*populate send_proc_obj*/
   index_ptr = (void *) ((char *) send_proc_obj->v_elements + index*object_size);

   memcpy(index_ptr, p_recv_contact_buf , object_size*contact_size);

   send_proc_obj->vec_starts[count+1] = index + contact_size;
   send_proc_obj->length++;

   /* output - no message to return (confirmation) */
  *response_message_size = 0; 

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixAddToValueshpCSR( jxf_IJMatrix   *matrix,
                              JXF_Int            nrows,
                              JXF_Int           *ncols,
                              const JXF_Int     *rows,
                              const JXF_Int     *cols,
                              const JXF_Real  *values )
{
   jxf_hpCSRMatrix *hp_matrix;
   jxf_CSRMatrix *diag, *offd;
   jxf_AuxParCSRMatrix *aux_matrix;
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int num_procs, my_id;
   JXF_Int row_local, row;
   JXF_Int col_0, col_n;
   JXF_Int i, ii, j, n, not_found;
   JXF_Int **aux_j;
   JXF_Int *local_j;
   JXF_Int *tmp_j;
   JXF_Real **aux_data;
   JXF_Real  *local_data;
   JXF_Real  *tmp_data = NULL;
   JXF_Int diag_space, offd_space;
   JXF_Int *row_length, *row_space;
   JXF_Int need_aux;
   JXF_Int tmp_indx, indx;
   JXF_Int space, size, old_size;
   JXF_Int cnt, cnt_diag, cnt_offd;
   JXF_Int pos_diag, pos_offd;
   JXF_Int len_diag, len_offd;
   JXF_Int offd_indx, diag_indx;
   JXF_Int first, pstart;
   JXF_Int *diag_i;
   JXF_Int *diag_j;
   JXF_Real *diag_data;
   JXF_Int *offd_i;
   JXF_Int *offd_j = NULL;
   JXF_Real *offd_data = NULL;
   JXF_Int current_num_elmts;
   JXF_Int max_off_proc_elmts;
   JXF_Int off_proc_i_indx;
   JXF_Int *off_proc_i;
   JXF_Int *off_proc_j;
   JXF_Real *off_proc_data;
   JXF_Int print_level = jxf_IJMatrixPrintLevel(matrix);

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   hp_matrix = jxf_IJMatrixObject( matrix );
   row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   col_partitioning = jxf_IJMatrixColPartitioning(matrix);
#ifdef JXF_NO_GLOBAL_PARTITION
   col_0 = col_partitioning[0];
   col_n = col_partitioning[1]-1;
   first = jxf_IJMatrixGlobalFirstCol(matrix);
   pstart = 0;
#else
   col_0 = col_partitioning[my_id];
   col_n = col_partitioning[my_id+1]-1;
   first = col_partitioning[0];
   pstart = my_id;
#endif
   if (jxf_IJMatrixAssembleFlag(matrix))
   {
      JXF_Int num_cols_offd;
      JXF_Int *col_map_offd = NULL;
      JXF_Int j_offd;
      indx = 0;

      /* AB - 4/06 - need to get this object*/
      aux_matrix = jxf_IJMatrixTranslator(matrix);

      for (ii=0; ii < nrows; ii++)
      {
         row = rows[ii];
         n = ncols[ii];
         if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
         {
            row_local = row - row_partitioning[pstart]; 
            /* compute local row number */
            diag = jxf_hpCSRMatrixDiag(hp_matrix);
            diag_i = jxf_CSRMatrixI(diag);
            diag_j = jxf_CSRMatrixJ(diag);
            diag_data = jxf_CSRMatrixData(diag);
            offd = jxf_hpCSRMatrixOffd(hp_matrix);
            offd_i = jxf_CSRMatrixI(offd);
            num_cols_offd = jxf_CSRMatrixNumCols(offd);
            if (num_cols_offd)
            {
               col_map_offd = jxf_hpCSRMatrixColMapOffd(hp_matrix);
               offd_j = jxf_CSRMatrixJ(offd);
               offd_data = jxf_CSRMatrixData(offd);
            }
            size = diag_i[row_local+1] - diag_i[row_local]
               + offd_i[row_local+1] - offd_i[row_local];
      
            if (n > size)
            {
               jxf_error(JXF_ERROR_GENERIC);
      	       if (print_level) jxf_printf (" row %d too long! \n", row);
      	       /* return -1; */
               return jxf_error_flag;
            }

            pos_diag = diag_i[row_local];
            pos_offd = offd_i[row_local];
            len_diag = diag_i[row_local+1];
            len_offd = offd_i[row_local+1];
            not_found = 1;
      	
            for (i=0; i < n; i++)
            {
               if (cols[indx] < col_0 || cols[indx] > col_n)
                  /* insert into offd */	
               {
      	          j_offd = jxf_BinarySearch(col_map_offd,cols[indx]-first,
                                              num_cols_offd);
      	          if (j_offd == -1)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
                     return jxf_error_flag;
      	             /* return -1; */
      	          }
      	          for (j=pos_offd; j < len_offd; j++)
      	          {
      	             if (offd_j[j] == j_offd)
      	             {
                        offd_data[j] += values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
      	             /* return -1;*/
                     return jxf_error_flag;
      	          }
      	          not_found = 1;
               }
               /* diagonal element */
      	       else if (cols[indx] == row)
      	       {
      	          if (diag_j[pos_diag] != row_local)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
      	             /* return -1; */
                     return jxf_error_flag;
      	          }
      	          diag_data[pos_diag] += values[indx];
      	       }
               else  /* insert into diag */
               {
      	          for (j=pos_diag; j < len_diag; j++)
      	          {
      	             if (diag_j[j] == (cols[indx]-col_0))
      	             {
                        diag_data[j] += values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
                     jxf_error(JXF_ERROR_GENERIC);
      	             if (print_level)
			jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
      	             /* return -1;*/
                     return jxf_error_flag;
      	          }
               }
               indx++;
            }
         }
         /* not my row */
         else
	 {
   	    if (!aux_matrix)
            {
               size = row_partitioning[pstart+1]-row_partitioning[pstart];
	       jxf_AuxParCSRMatrixCreate(&aux_matrix, size, size, NULL);
      	       jxf_AuxParCSRMatrixNeedAux(aux_matrix) = 0;
      	       jxf_IJMatrixTranslator(matrix) = aux_matrix;
            }
   	    current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
   	    max_off_proc_elmts = jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix);
   	    off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
   	    off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	    off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	    off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);

	    if (!max_off_proc_elmts)
	    {
	       max_off_proc_elmts = jxf_max(n,1000);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) =
                  max_off_proc_elmts;
   	       jxf_AuxParCSRMatrixOffProcI(aux_matrix)
                  = jxf_CTAlloc(JXF_Int,2*max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix)
                  = jxf_CTAlloc(JXF_Int,max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcData(aux_matrix)
                  = jxf_CTAlloc(JXF_Real,max_off_proc_elmts);
   	       off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	       off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	       off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);
	    }
            else if (current_num_elmts + n > max_off_proc_elmts)
            {
               max_off_proc_elmts += 3*n;
               off_proc_i = jxf_TReAlloc(off_proc_i,JXF_Int,2*max_off_proc_elmts);
               off_proc_j = jxf_TReAlloc(off_proc_j,JXF_Int,max_off_proc_elmts);
               off_proc_data = jxf_TReAlloc(off_proc_data,JXF_Real,max_off_proc_elmts);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) = max_off_proc_elmts;
	       jxf_AuxParCSRMatrixOffProcI(aux_matrix) = off_proc_i;
	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix) = off_proc_j;
	       jxf_AuxParCSRMatrixOffProcData(aux_matrix) = off_proc_data;
	    }

            /* AB - 4/6 - the row should be negative to indicate an add */
            /* UMY - 12/28/09 - now positive since we eliminated the feature of
               setting on other processors */
            /* off_proc_i[off_proc_i_indx++] = row; */
            off_proc_i[off_proc_i_indx++] = row;
            
            off_proc_i[off_proc_i_indx++] = n; 
	    for (i=0; i < n; i++)
	    {
	       off_proc_j[current_num_elmts] = cols[indx];
	       off_proc_data[current_num_elmts++] = values[indx++];
	    }
	    jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix) = off_proc_i_indx; 
	    jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix) = current_num_elmts; 
	 }
      }
   }
   
   /* not assembled */
   else
   {
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      row_space = jxf_AuxParCSRMatrixRowSpace(aux_matrix);
      row_length = jxf_AuxParCSRMatrixRowLength(aux_matrix);
      need_aux = jxf_AuxParCSRMatrixNeedAux(aux_matrix);
      indx = 0;   
      for (ii=0; ii < nrows; ii++)
      {
         row = rows[ii];
         n = ncols[ii];
         if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
         {
            row_local = row - row_partitioning[pstart]; 
            /* compute local row number */
            if (need_aux)
            {
               aux_j = jxf_AuxParCSRMatrixAuxJ(aux_matrix);
               aux_data = jxf_AuxParCSRMatrixAuxData(aux_matrix);
               local_j = aux_j[row_local];
               local_data = aux_data[row_local];
   	       space = row_space[row_local]; 
   	       old_size = row_length[row_local]; 
   	       size = space - old_size;
   	       if (size < n)
      	       {
      	          size = n - size;
      	          tmp_j = jxf_CTAlloc(JXF_Int,size);
      	          tmp_data = jxf_CTAlloc(JXF_Real,size);
      	       }
      	       else
      	       {
      	          tmp_j = NULL;
      	       }
      	       tmp_indx = 0;
      	       not_found = 1;
      	       size = old_size;
               for (i=0; i < n; i++)
      	       {
      	          for (j=0; j < old_size; j++)
      	          {
      	             if (local_j[j] == cols[indx])
      	             {
                        local_data[j] += values[indx];
      		        not_found = 0;
      		        break;
      	             }
      	          }
      	          if (not_found)
      	          {
      	             if (size < space)
      	             {
      	                local_j[size] = cols[indx];
      	                local_data[size++] = values[indx];
      	             }
      	             else
      	             {
      	                tmp_j[tmp_indx] = cols[indx];
      	                tmp_data[tmp_indx++] = values[indx];
      	             }
      	          }
      	          not_found = 1;
        	  indx++;
      	       }

               row_length[row_local] = size+tmp_indx;

               if (tmp_indx)
               {
   	          aux_j[row_local] = jxf_TReAlloc(aux_j[row_local],JXF_Int,size+tmp_indx);
   	          aux_data[row_local] = jxf_TReAlloc(aux_data[row_local],JXF_Real,size+tmp_indx);
                  row_space[row_local] = size+tmp_indx;
                  local_j = aux_j[row_local];
                  local_data = aux_data[row_local];
               }

   	       cnt = size; 

   	       for (i=0; i < tmp_indx; i++)
   	       {
   	          local_j[cnt] = tmp_j[i];
   	          local_data[cnt++] = tmp_data[i];
	       }
  
	       if (tmp_j)
	       { 
	          jxf_TFree(tmp_j); 
	          jxf_TFree(tmp_data); 
	       } 
            }
            else /* insert immediately into data in ParCSRMatrix structure */
            {
	       offd_indx = jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local];
	       diag_indx = jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local];
               diag = jxf_hpCSRMatrixDiag(hp_matrix);
               diag_i = jxf_CSRMatrixI(diag);
               diag_j = jxf_CSRMatrixJ(diag);
               diag_data = jxf_CSRMatrixData(diag);
               offd = jxf_hpCSRMatrixOffd(hp_matrix);
               offd_i = jxf_CSRMatrixI(offd);
               if (num_procs > 1)
	       {
	          offd_j = jxf_CSRMatrixJ(offd);
                  offd_data = jxf_CSRMatrixData(offd);
               }
	       cnt_diag = diag_indx;
	       cnt_offd = offd_indx;
	       diag_space = diag_i[row_local+1];
	       offd_space = offd_i[row_local+1];
	       not_found = 1;
  	       for (i=0; i < n; i++)
	       {
	          if (cols[indx] < col_0 || cols[indx] > col_n)
                     /* insert into offd */	
	          {
	             for (j=offd_i[row_local]; j < offd_indx; j++)
	             {
		        if (offd_j[j] == cols[indx])
		        {
                           offd_data[j] += values[indx];
		           not_found = 0;
		           break;
		        }
	             }
	             if (not_found)
	             { 
	                if (cnt_offd < offd_space) 
	                { 
	                   offd_j[cnt_offd] = cols[indx];
	                   offd_data[cnt_offd++] = values[indx];
	                } 
	                else 
	 	        {
                           jxf_error(JXF_ERROR_GENERIC);
	    	           if (print_level)
                              jxf_printf("Error in row %d ! Too many elements!\n", row);
	    	           /* return 1;*/
                           return jxf_error_flag;
	 	        }
	             }  
	             not_found = 1;
	          }
	          else  /* insert into diag */
	          {
	             for (j=diag_i[row_local]; j < diag_indx; j++)
	             {
		        if (diag_j[j] == cols[indx])
		        {
                           diag_data[j] += values[indx];
		           not_found = 0;
		           break;
		        }
	             } 
	             if (not_found)
	             { 
	                if (cnt_diag < diag_space) 
	                { 
	                   diag_j[cnt_diag] = cols[indx];
	                   diag_data[cnt_diag++] = values[indx];
	                } 
	                else 
	 	        {
                           jxf_error(JXF_ERROR_GENERIC);
	    	           if (print_level)
                              jxf_printf("Error in row %d ! Too many elements !\n", row);
	    	           /* return 1; */
                           return jxf_error_flag;
	 	        }
	             } 
	             not_found = 1;
	          }
	          indx++;
	       }

               jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local] = cnt_diag;
               jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local] = cnt_offd;

            }
         }
         /* not my row */
         else
         {
   	    current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
   	    max_off_proc_elmts = jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix);
   	    off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
   	    off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	    off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	    off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);

	    if (!max_off_proc_elmts)
	    {
	       max_off_proc_elmts = jxf_max(n,1000);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) = max_off_proc_elmts;
   	       jxf_AuxParCSRMatrixOffProcI(aux_matrix) = jxf_CTAlloc(JXF_Int,2*max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix) = jxf_CTAlloc(JXF_Int,max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcData(aux_matrix) = jxf_CTAlloc(JXF_Real,max_off_proc_elmts);
   	       off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	       off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	       off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);
	    }
            else if (current_num_elmts + n > max_off_proc_elmts)
            {
               max_off_proc_elmts += 3*n;
               off_proc_i = jxf_TReAlloc(off_proc_i,JXF_Int,2*max_off_proc_elmts);
               off_proc_j = jxf_TReAlloc(off_proc_j,JXF_Int,max_off_proc_elmts);
               off_proc_data = jxf_TReAlloc(off_proc_data,JXF_Real,max_off_proc_elmts);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) = max_off_proc_elmts;
	       jxf_AuxParCSRMatrixOffProcI(aux_matrix) = off_proc_i;
	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix) = off_proc_j;
	       jxf_AuxParCSRMatrixOffProcData(aux_matrix) = off_proc_data;
	    }
            off_proc_i[off_proc_i_indx++] = row; 
            off_proc_i[off_proc_i_indx++] = n; 
	    for (i=0; i < n; i++)
	    {
	       off_proc_j[current_num_elmts] = cols[indx];
	       off_proc_data[current_num_elmts++] = values[indx++];
	    }
	    jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix) = off_proc_i_indx; 
	    jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix) = current_num_elmts; 
         }
      }
   }

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixAddToValues( JXF_IJMatrix   matrix,
                        JXF_Int           nrows,
                        JXF_Int          *ncols,
                        const JXF_Int    *rows,
                        const JXF_Int    *cols,
                        const JXF_Real *values )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;
   if (nrows == 0)
      return jxf_error_flag;
   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (nrows < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   if (!ncols)
   {
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (!rows)
   {
      jxf_error_in_arg(4);
      return jxf_error_flag;
   }
   if (!cols)
   {
      jxf_error_in_arg(5);
      return jxf_error_flag;
   }
   if (!values)
   {
      jxf_error_in_arg(6);
      return jxf_error_flag;
   }
   if (jxf_IJMatrixObjectType(ijmatrix) == JXF_HPCSR)
   {
      if (jxf_IJMatrixOMPFlag(ijmatrix))
	 return( jxf_IJMatrixAddToValuesOMPhpCSR( ijmatrix, nrows, ncols, rows, cols, values ) );
      else
         return( jxf_IJMatrixAddToValueshpCSR( ijmatrix, nrows, ncols, rows, cols, values ) );
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixRead(  const char  *filename,
                  MPI_Comm     comm,
                  JXF_Int          type,
		            JXF_IJMatrix *matrix_ptr )
{
   JXF_IJMatrix matrix;
   JXF_Int       ilower, iupper, jlower, jupper;
   JXF_Int       ncols, I, J;
   JXF_Real    value;
   JXF_Int       myid, ret;
   char      new_filename[255];
   FILE     *file;

   jxf_MPI_Comm_rank(comm, &myid);
   jxf_sprintf(new_filename, "%s.%05d", filename, myid);
   if ((file = fopen(new_filename, "r")) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_fscanf(file, "%d %d %d %d", &ilower, &iupper, &jlower, &jupper);
   JXF_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &matrix);
   JXF_IJMatrixSetObjectType(matrix, type);
   JXF_IJMatrixInitialize(matrix);

   /* It is important to ensure that whitespace follows the index value to help
    * catch mistakes in the input file.  See comments in IJVectorRead(). */
   ncols = 1;
   while ( (ret = jxf_fscanf(file, "%d %d%*[ \t]%le", &I, &J, &value)) != EOF )
   {
      if (ret != 3)
      {
         jxf_error_w_msg(JXF_ERROR_GENERIC, "Error in IJ matrix input file.");
         return jxf_error_flag;
      }
      if (I < ilower || I > iupper)
         JXF_IJMatrixAddToValues(matrix, 1, &ncols, &I, &J, &value);
      else
         JXF_IJMatrixSetValues(matrix, 1, &ncols, &I, &J, &value);
   }
   JXF_IJMatrixAssemble(matrix);
   fclose(file);
  *matrix_ptr = matrix;

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixPrint( JXF_IJMatrix matrix, const char *filename )
{
   MPI_Comm  comm;
   JXF_Int      *row_partitioning;
   JXF_Int      *col_partitioning;
   JXF_Int       ilower, iupper, jlower, jupper;
   JXF_Int       i, j, ii = 0;
   JXF_Int       ncols, *cols;
   JXF_Real   *values;
   JXF_Int       myid;
   char      new_filename[255];
   FILE     *file;
   void     *object;

   if (!matrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if ((jxf_IJMatrixObjectType(matrix) != JXF_HPCSR))
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   comm = jxf_IJMatrixComm(matrix);
   jxf_MPI_Comm_rank(comm, &myid);
   jxf_sprintf(new_filename,"%s.%05d", filename, myid);
   if ((file = fopen(new_filename, "w")) == NULL)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   col_partitioning = jxf_IJMatrixColPartitioning(matrix);
#ifdef JXF_NO_GLOBAL_PARTITION
   ilower = row_partitioning[0];
   iupper = row_partitioning[1] - 1;
   jlower = col_partitioning[0];
   jupper = col_partitioning[1] - 1;
#else
   ilower = row_partitioning[myid];
   iupper = row_partitioning[myid+1] - 1;
   jlower = col_partitioning[myid];
   jupper = col_partitioning[myid+1] - 1;
#endif
   jxf_fprintf(file, "%d %d %d %d\n", ilower, iupper, jlower, jupper);
   JXF_IJMatrixGetObject(matrix, &object);
   for (i = ilower; i <= iupper; i++)
   {
      if (jxf_IJMatrixObjectType(matrix) == JXF_HPCSR)
      {
#ifdef JXF_NO_GLOBAL_PARTITION
         ii = i -  jxf_IJMatrixGlobalFirstRow(matrix);
#else
         ii = i - row_partitioning[0];
#endif
         JXF_ParCSRMatrixGetRow((JXF_ParCSRMatrix) jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)object), ii, &ncols, &cols, &values);
         for (j = 0; j < ncols; j++)
         {
#ifdef JXF_NO_GLOBAL_PARTITION
            cols[j] +=  jxf_IJMatrixGlobalFirstCol(matrix);
#else
            cols[j] += col_partitioning[0];
#endif
         }
      }
      for (j = 0; j < ncols; j++)
      {
         jxf_fprintf(file, "%d %d %.14e\n", i, cols[j], values[j]);
      }
      if (jxf_IJMatrixObjectType(matrix) == JXF_HPCSR)
      {
         for (j = 0; j < ncols; j++)
         {
#ifdef JXF_NO_GLOBAL_PARTITION
            cols[j] -=  jxf_IJMatrixGlobalFirstCol(matrix);
#else
            cols[j] -= col_partitioning[0];
#endif
         }
         JXF_ParCSRMatrixRestoreRow((JXF_ParCSRMatrix) jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *)object), ii, &ncols, &cols, &values);
      }
   }
   fclose(file);

   return jxf_error_flag;
}

JXF_Int
jxf_IJMatrixAddToValuesOMPhpCSR( jxf_IJMatrix  *matrix,
                                 JXF_Int           nrows,
                                 JXF_Int          *ncols,
                                 const JXF_Int    *rows,
                                 const JXF_Int    *cols,
                                 const JXF_Real *values )
{
   jxf_hpCSRMatrix *hp_matrix;
   jxf_CSRMatrix *diag, *offd;
   jxf_AuxParCSRMatrix *aux_matrix;
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   MPI_Comm comm = jxf_IJMatrixComm(matrix);
   JXF_Int num_procs, my_id;
   JXF_Int col_0, col_n;
   JXF_Int **aux_j = NULL;
   JXF_Real **aux_data = NULL;
   JXF_Int *row_length, *row_space;
   JXF_Int need_aux;
   JXF_Int first, pstart;
   JXF_Int *diag_i = NULL;
   JXF_Int *diag_j = NULL;
   JXF_Real *diag_data = NULL;
   JXF_Int *offd_i = NULL;
   JXF_Int *offd_j = NULL;
   JXF_Real *offd_data = NULL;
   JXF_Int current_num_elmts;
   JXF_Int max_off_proc_elmts;
   JXF_Int off_proc_i_indx;
   JXF_Int *off_proc_i;
   JXF_Int *off_proc_j;
   JXF_Real *off_proc_data;
   JXF_Int *value_start, **offproc_cnt;
   JXF_Int print_level = jxf_IJMatrixPrintLevel(matrix);
   JXF_Int max_num_threads;
   JXF_Int error_flag = 0;
   JXF_Int i1;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   max_num_threads = jxf_NumThreads();
   hp_matrix = jxf_IJMatrixObject( matrix );
   row_partitioning = jxf_IJMatrixRowPartitioning(matrix);
   col_partitioning = jxf_IJMatrixColPartitioning(matrix);
   value_start = jxf_CTAlloc(JXF_Int, max_num_threads+1);
   offproc_cnt = jxf_CTAlloc(JXF_Int *, max_num_threads);
   for (i1=0; i1 < max_num_threads; i1++)
      offproc_cnt[i1] = NULL;

#ifdef JXF_NO_GLOBAL_PARTITION
   col_0 = col_partitioning[0];
   col_n = col_partitioning[1]-1;
   first = jxf_IJMatrixGlobalFirstCol(matrix);
   pstart = 0;
#else
   col_0 = col_partitioning[my_id];
   col_n = col_partitioning[my_id+1]-1;
   first = col_partitioning[0];
   pstart = my_id;
#endif
   if (jxf_IJMatrixAssembleFlag(matrix)) /* matrix already assembled */
   {
      JXF_Int num_cols_offd;
      JXF_Int *col_map_offd = NULL;

      diag = jxf_hpCSRMatrixDiag(hp_matrix);
      diag_i = jxf_CSRMatrixI(diag);
      diag_j = jxf_CSRMatrixJ(diag);
      diag_data = jxf_CSRMatrixData(diag);
      offd = jxf_hpCSRMatrixOffd(hp_matrix);
      offd_i = jxf_CSRMatrixI(offd);
      num_cols_offd = jxf_CSRMatrixNumCols(offd);
      if (num_cols_offd)
      {
          col_map_offd = jxf_hpCSRMatrixColMapOffd(hp_matrix);
          offd_j = jxf_CSRMatrixJ(offd);
          offd_data = jxf_CSRMatrixData(offd);
      }
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      if (aux_matrix)
      {
         current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
         off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
      }
#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
      {
         JXF_Int j_offd;
         JXF_Int num_threads, my_thread_num;
         JXF_Int len, rest, ns, ne;
         JXF_Int pos_diag, pos_offd;
         JXF_Int len_diag, len_offd;
         JXF_Int row_local;
         JXF_Int i, j, ii, n, row;
         JXF_Int not_found, size, indx;
         JXF_Int *my_offproc_cnt = NULL;

         num_threads = jxf_NumActiveThreads();
         my_thread_num = jxf_GetThreadNum();

         len = nrows/num_threads; 
         rest = nrows - len*num_threads;

         if (my_thread_num < rest)
         {
            ns = my_thread_num*(len+1);
            ne = (my_thread_num+1)*(len+1);
         }
         else
         {
            ns = my_thread_num*len+rest;
            ne = (my_thread_num+1)*len+rest;
         }

         value_start[my_thread_num] = 0;
         for (ii=ns; ii < ne; ii++)
            value_start[my_thread_num] += ncols[ii];

#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         if (my_thread_num == 0)
         {
            for (i=0; i < max_num_threads; i++)
               value_start[i+1] += value_start[i];
         }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         indx = 0;
         if (my_thread_num) indx = value_start[my_thread_num-1];
         for (ii=ns; ii < ne; ii++)
         {
            row = rows[ii];
            n = ncols[ii];
            if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
            {
               row_local = row - row_partitioning[pstart]; 
               /* compute local row number */
               size = diag_i[row_local+1] - diag_i[row_local] + offd_i[row_local+1] - offd_i[row_local];
               if (n > size)
               {
                  jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                  error_flag++;
      	          if (print_level) jxf_printf (" row %d too long! \n", row);
                  break;
               }
       
               pos_diag = diag_i[row_local];
               pos_offd = offd_i[row_local];
               len_diag = diag_i[row_local+1];
               len_offd = offd_i[row_local+1];
               not_found = 1;
      	
               for (i=0; i < n; i++)
               {
                  if (cols[indx] < col_0 || cols[indx] > col_n)
                  /* insert into offd */	
                  {
      	             j_offd = jxf_BinarySearch(col_map_offd,cols[indx]-first,num_cols_offd);
      	             if (j_offd == -1)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
                        break;
      	             }
      	             for (j=pos_offd; j < len_offd; j++)
      	             {
      	                if (offd_j[j] == j_offd)
      	                {
                           offd_data[j] += values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
                        break;
      	             }
      	             not_found = 1;
                  }
                  /* diagonal element */
      	          else if (cols[indx] == row)
      	          {
      	             if (diag_j[pos_diag] != row_local)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
                        break;
      	             }
      	             diag_data[pos_diag] += values[indx];
      	          }
                  else  /* insert into diag */
                  {
      	             for (j=pos_diag; j < len_diag; j++)
      	             {
      	                if (diag_j[j] == (cols[indx]-col_0))
      	                {
                           diag_data[j] += values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
                        jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                        error_flag++;
      	                if (print_level)
			   jxf_printf(" Error, element %d %d does not exist\n", row, cols[indx]);
                        break;
      	             }
                  }
                  indx++;
               }
            }
            /* not my row */
            /* need to find solution for threaded version!!!! */
            /* could save row number and process later .... */
            else
	    {
               if (!my_offproc_cnt)
               {
                  my_offproc_cnt = jxf_CTAlloc(JXF_Int, 200);
                  offproc_cnt[my_thread_num] = my_offproc_cnt;
                  my_offproc_cnt[0] = 200;
	          my_offproc_cnt[1] = 2;
               }
               i = my_offproc_cnt[1];
               if (i+2 < my_offproc_cnt[0])
               {
                  my_offproc_cnt[i] = ii;
                  my_offproc_cnt[i+1] = indx;
                  my_offproc_cnt[1] += 2;
               }
               else
               {
                  size = my_offproc_cnt[0];
                  my_offproc_cnt = jxf_TReAlloc(my_offproc_cnt,JXF_Int,size+200);
                  my_offproc_cnt[0] += 200;
                  my_offproc_cnt[i] = ii;
                  my_offproc_cnt[i+1] = indx;
                  my_offproc_cnt[1] += 2;
               }
               indx +=n;
	    }
         }
      } /* end parallel region */
   }
   
   /* not assembled */
   else
   {
      aux_matrix = jxf_IJMatrixTranslator(matrix);
      if (aux_matrix)
      {
         current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
         off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
         off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
         off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
      }
      row_space = jxf_AuxParCSRMatrixRowSpace(aux_matrix);
      row_length = jxf_AuxParCSRMatrixRowLength(aux_matrix);
      need_aux = jxf_AuxParCSRMatrixNeedAux(aux_matrix);
      if (need_aux)
      {
         aux_j = jxf_AuxParCSRMatrixAuxJ(aux_matrix);
         aux_data = jxf_AuxParCSRMatrixAuxData(aux_matrix);
      }
      else
      {
         diag = jxf_hpCSRMatrixDiag(hp_matrix);
         diag_i = jxf_CSRMatrixI(diag);
         diag_j = jxf_CSRMatrixJ(diag);
         diag_data = jxf_CSRMatrixData(diag);
         offd = jxf_hpCSRMatrixOffd(hp_matrix);
         offd_i = jxf_CSRMatrixI(offd);
         if (num_procs > 1)
         {
            offd_j = jxf_CSRMatrixJ(offd);
            offd_data = jxf_CSRMatrixData(offd);
         }
      }
#if JXF_USING_OPENMP
#pragma omp parallel 
#endif
      {
         JXF_Int num_threads, my_thread_num;
         JXF_Int len, rest, ns, ne;
         JXF_Int *tmp_j = NULL;
         JXF_Int *local_j = NULL;
         JXF_Real *tmp_data = NULL;
         JXF_Real *local_data = NULL;
         JXF_Int tmp_indx;
         JXF_Int row_local;
         JXF_Int i, j, ii, n, row;
         JXF_Int not_found, size, indx;
         JXF_Int old_size, space, cnt;
         JXF_Int *my_offproc_cnt = NULL;

         num_threads = jxf_NumActiveThreads();
         my_thread_num = jxf_GetThreadNum();

         len = nrows/num_threads; 
         rest = nrows - len*num_threads;

         if (my_thread_num < rest)
         {
            ns = my_thread_num*(len+1);
            ne = (my_thread_num+1)*(len+1);
         }
         else
         {
            ns = my_thread_num*len+rest;
            ne = (my_thread_num+1)*len+rest;
         }

         value_start[my_thread_num] = 0;
         for (ii=ns; ii < ne; ii++)
            value_start[my_thread_num] += ncols[ii];

#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         if (my_thread_num == 0)
         {
            for (i=0; i < max_num_threads; i++)
               value_start[i+1] += value_start[i];
         }
#if JXF_USING_OPENMP
#pragma omp barrier
#endif
         indx = 0;
         if (my_thread_num) indx = value_start[my_thread_num-1];
         for (ii=ns; ii < ne; ii++)
         {
            row = rows[ii];
            n = ncols[ii];
            if (row >= row_partitioning[pstart] && row < row_partitioning[pstart+1])
            {
               row_local = row - row_partitioning[pstart]; 
               /* compute local row number */
               if (need_aux)
               {
                  local_j = aux_j[row_local];
                  local_data = aux_data[row_local];
   	          space = row_space[row_local]; 
   	          old_size = row_length[row_local]; 
   	          size = space - old_size;
   	          if (size < n)
      	          {
      	             size = n - size;
      	             tmp_j = jxf_CTAlloc(JXF_Int,size);
      	             tmp_data = jxf_CTAlloc(JXF_Real,size);
      	          }
      	          tmp_indx = 0;
      	          not_found = 1;
      	          size = old_size;
                  for (i=0; i < n; i++)
      	          {
      	             for (j=0; j < old_size; j++)
      	             {
      	                if (local_j[j] == cols[indx])
      	                {
                           local_data[j] += values[indx];
      		           not_found = 0;
      		           break;
      	                }
      	             }
      	             if (not_found)
      	             {
      	                if (size < space)
      	                {
      	                   local_j[size] = cols[indx];
      	                   local_data[size++] = values[indx];
      	                }
      	                else
      	                {
      	                   tmp_j[tmp_indx] = cols[indx];
      	                   tmp_data[tmp_indx++] = values[indx];
      	                }
      	             }
      	             not_found = 1;
        	     indx++;
      	          }
      	    
                  row_length[row_local] = size+tmp_indx;
               
                  if (tmp_indx)
                  {
   	             aux_j[row_local] = jxf_TReAlloc(aux_j[row_local],JXF_Int,size+tmp_indx);
   	             aux_data[row_local] = jxf_TReAlloc(aux_data[row_local],JXF_Real,size+tmp_indx);
                     row_space[row_local] = size+tmp_indx;
                     local_j = aux_j[row_local];
                     local_data = aux_data[row_local];
                  }
   
   	          cnt = size; 
   
   	          for (i=0; i < tmp_indx; i++)
   	          {
   	             local_j[cnt] = tmp_j[i];
   	             local_data[cnt++] = tmp_data[i];
	          }
  
	          if (tmp_j)
	          { 
	             jxf_TFree(tmp_j); 
	             jxf_TFree(tmp_data); 
	          } 
               }
               else /* insert immediately into data in ParCSRMatrix structure */
               {
                  JXF_Int offd_indx, diag_indx;
                  JXF_Int offd_space, diag_space;
                  JXF_Int cnt_diag, cnt_offd;
	          offd_indx = jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local];
	          diag_indx = jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local];
	          cnt_diag = diag_indx;
	          cnt_offd = offd_indx;
	          diag_space = diag_i[row_local+1];
	          offd_space = offd_i[row_local+1];
	          not_found = 1;
  	          for (i=0; i < n; i++)
	          {
	             if (cols[indx] < col_0 || cols[indx] > col_n)
                     /* insert into offd */	
	             {
	                for (j=offd_i[row_local]; j < offd_indx; j++)
	                {
		           if (offd_j[j] == cols[indx])
		           {
                              offd_data[j] += values[indx];
		              not_found = 0;
		              break;
		           }
	                }
	                if (not_found)
	                { 
	                   if (cnt_offd < offd_space) 
	                   { 
	                      offd_j[cnt_offd] = cols[indx];
	                      offd_data[cnt_offd++] = values[indx];
	                   } 
	                   else 
	 	           {
                              jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                              error_flag++;
	    	              if (print_level)
                                 jxf_printf("Error in row %d ! Too many elements!\n", row);
                              break;
	 	           }
	                }  
	                not_found = 1;
	             }
	             else  /* insert into diag */
	             {
	                for (j=diag_i[row_local]; j < diag_indx; j++)
	                {
		           if (diag_j[j] == cols[indx])
		           {
                              diag_data[j] += values[indx];
		              not_found = 0;
		              break;
		           }
	                } 
	                if (not_found)
	                { 
	                   if (cnt_diag < diag_space) 
	                   { 
	                      diag_j[cnt_diag] = cols[indx];
	                      diag_data[cnt_diag++] = values[indx];
	                   } 
	                   else 
	 	           {
                              jxf_error(JXF_ERROR_GENERIC);
#if JXF_USING_OPENMP
#pragma omp atomic
#endif
                              error_flag++;
	    	              if (print_level)
                                 jxf_printf("Error in row %d ! Too many elements !\n", row);
                              break;
	 	           }
	                } 
	                not_found = 1;
	             }
	             indx++;
	          }

                  jxf_AuxParCSRMatrixIndxDiag(aux_matrix)[row_local] = cnt_diag;
                  jxf_AuxParCSRMatrixIndxOffd(aux_matrix)[row_local] = cnt_offd;

               }
            }
            /* not my row */
            else
            {
//               JXF_Int offproc_indx = 0;
               if (!my_offproc_cnt)
               {
                  my_offproc_cnt = jxf_CTAlloc(JXF_Int, 200);
                  offproc_cnt[my_thread_num] = my_offproc_cnt;
                  my_offproc_cnt[0] = 200;
	          my_offproc_cnt[1] = 2;
//	          offproc_indx = 2;
               }
               i = my_offproc_cnt[1];
               if (i+2 < my_offproc_cnt[0])
               {
                  my_offproc_cnt[i] = ii;
                  my_offproc_cnt[i+1] = indx;
                  my_offproc_cnt[1] += 2;
               }
               else
               {
                  size = my_offproc_cnt[0];
                  my_offproc_cnt = jxf_TReAlloc(my_offproc_cnt,JXF_Int,size+200);
                  my_offproc_cnt[0] += 200;
                  my_offproc_cnt[i] = ii;
                  my_offproc_cnt[i+1] = indx;
                  my_offproc_cnt[1] += 2;
               }
               indx +=n;
            }
         }
      } /*end parallel region */
   }
   if (error_flag) return jxf_error_flag;
   jxf_TFree(value_start);
   if (!aux_matrix)
   {
       JXF_Int size = row_partitioning[pstart+1]-row_partitioning[pstart];
       jxf_AuxParCSRMatrixCreate(&aux_matrix, size, size, NULL);
       jxf_AuxParCSRMatrixNeedAux(aux_matrix) = 0;
       jxf_IJMatrixTranslator(matrix) = aux_matrix;
   }
   for (i1 = 0; i1 < max_num_threads; i1++)
   {
      if (offproc_cnt[i1])
      {
         JXF_Int *my_offproc_cnt = offproc_cnt[i1];
         JXF_Int i, i2, ii, row, n, indx;
         for (i2 = 2; i2 < my_offproc_cnt[1]; i2+=2)
         {
            ii = my_offproc_cnt[i2];
            row = rows[ii];
            n = ncols[ii];
            indx = my_offproc_cnt[i2+1];
   	    current_num_elmts = jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix);
   	    max_off_proc_elmts = jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix);
   	    off_proc_i_indx = jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix);
   	    off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	    off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	    off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);
   	    
	    if (!max_off_proc_elmts)
	    {
	       max_off_proc_elmts = jxf_max(n,1000);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) = max_off_proc_elmts;
   	       jxf_AuxParCSRMatrixOffProcI(aux_matrix) = jxf_CTAlloc(JXF_Int,2*max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix) = jxf_CTAlloc(JXF_Int,max_off_proc_elmts);
   	       jxf_AuxParCSRMatrixOffProcData(aux_matrix) = jxf_CTAlloc(JXF_Real,max_off_proc_elmts);
   	       off_proc_i = jxf_AuxParCSRMatrixOffProcI(aux_matrix);
   	       off_proc_j = jxf_AuxParCSRMatrixOffProcJ(aux_matrix);
   	       off_proc_data = jxf_AuxParCSRMatrixOffProcData(aux_matrix);
	    }
            else if (current_num_elmts + n > max_off_proc_elmts)
            {
               max_off_proc_elmts += 3*n;
               off_proc_i = jxf_TReAlloc(off_proc_i,JXF_Int,2*max_off_proc_elmts);
               off_proc_j = jxf_TReAlloc(off_proc_j,JXF_Int,max_off_proc_elmts);
               off_proc_data = jxf_TReAlloc(off_proc_data,JXF_Real,max_off_proc_elmts);
	       jxf_AuxParCSRMatrixMaxOffProcElmts(aux_matrix) = max_off_proc_elmts;
	       jxf_AuxParCSRMatrixOffProcI(aux_matrix) = off_proc_i;
	       jxf_AuxParCSRMatrixOffProcJ(aux_matrix) = off_proc_j;
	       jxf_AuxParCSRMatrixOffProcData(aux_matrix) = off_proc_data;
	    }
            off_proc_i[off_proc_i_indx++] = row; 
            off_proc_i[off_proc_i_indx++] = n; 
	    for (i=0; i < n; i++)
	    {
	       off_proc_j[current_num_elmts] = cols[indx];
	       off_proc_data[current_num_elmts++] = values[indx++];
	    }
	    jxf_AuxParCSRMatrixOffProcIIndx(aux_matrix) = off_proc_i_indx; 
	    jxf_AuxParCSRMatrixCurrentNumElmts(aux_matrix) = current_num_elmts; 
	 }
	 jxf_TFree (offproc_cnt[i1]);
      }
   }
   jxf_TFree(offproc_cnt);
   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixGetObject( JXF_IJMatrix matrix, void **object )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;
   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *object = jxf_IJMatrixObject(ijmatrix);

   return jxf_error_flag;
}

JXF_Int
JXF_IJMatrixGetLocalRange( JXF_IJMatrix  matrix,
                          JXF_Int      *ilower,
                          JXF_Int      *iupper,
                          JXF_Int      *jlower,
                          JXF_Int      *jupper )
{
   jxf_IJMatrix *ijmatrix = (jxf_IJMatrix *) matrix;
   MPI_Comm comm;
   JXF_Int *row_partitioning;
   JXF_Int *col_partitioning;
   JXF_Int my_id;

   if (!ijmatrix)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   comm = jxf_IJMatrixComm(ijmatrix);
   row_partitioning = jxf_IJMatrixRowPartitioning(ijmatrix);
   col_partitioning = jxf_IJMatrixColPartitioning(ijmatrix);

   jxf_MPI_Comm_rank(comm, &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION
   *ilower = row_partitioning[0];
   *iupper = row_partitioning[1]-1;
   *jlower = col_partitioning[0];
   *jupper = col_partitioning[1]-1;
#else
   *ilower = row_partitioning[my_id];
   *iupper = row_partitioning[my_id+1]-1;
   *jlower = col_partitioning[my_id];
   *jupper = col_partitioning[my_id+1]-1;
#endif

   return jxf_error_flag;
}
