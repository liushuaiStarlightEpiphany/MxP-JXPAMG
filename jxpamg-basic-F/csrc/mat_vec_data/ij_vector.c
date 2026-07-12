//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  ij_vector.c -- basic operations for IJ vectors
 *  Date: 2015/11/28
 */

#include "jxf_mv.h"

JXF_Int
JXF_IJVectorCreate( MPI_Comm     comm,
                   JXF_Int          jlower,
                   JXF_Int          jupper,
                   JXF_IJVector *vector )
{
   jxf_IJVector *vec;
   JXF_Int num_procs, my_id, *partitioning;
#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int  row0, rowN;
#else
   JXF_Int *recv_buf;
   JXF_Int *info;
   JXF_Int i, i2;
#endif

   vec = jxf_CTAlloc(jxf_IJVector, 1);
   if (!vec)
   {
      jxf_error(JXF_ERROR_MEMORY);
      return jxf_error_flag;
   }

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (jlower > jupper+1 || jlower < 0)
   {
      jxf_error_in_arg(2);
      jxf_TFree(vec);
      return jxf_error_flag;
   }
   if (jupper < -1)
   {
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION

   partitioning = jxf_CTAlloc(JXF_Int, 2);

   partitioning[0] = jlower;
   partitioning[1] = jupper+1;

   /* now we need the global number of rows as well
      as the global first row index */

   /* proc 0 has the first row  */
   if (my_id==0) 
   {
      row0 = jlower;
   }
   jxf_MPI_Bcast(&row0, 1, JXF_MPI_INT, 0, comm);
   /* proc (num_procs-1) has the last row */
   if (my_id == (num_procs-1))
   {
      rowN = jupper;
   }
   jxf_MPI_Bcast(&rowN, 1, JXF_MPI_INT, num_procs-1, comm);

   jxf_IJVectorGlobalFirstRow(vec) = row0;
   jxf_IJVectorGlobalNumRows(vec) = rowN - row0 + 1;

#else

   info = jxf_CTAlloc(JXF_Int,2);
   recv_buf = jxf_CTAlloc(JXF_Int, 2*num_procs);
   partitioning = jxf_CTAlloc(JXF_Int, num_procs+1);

   info[0] = jlower;
   info[1] = jupper;

   jxf_MPI_Allgather(info, 2, JXF_MPI_INT, recv_buf, 2, JXF_MPI_INT, comm);

   partitioning[0] = recv_buf[0];
   for (i=0; i < num_procs-1; i++)
   {
      i2 = i+i;
      if (recv_buf[i2+1] != (recv_buf[i2+2]-1))
      {
	 jxf_error(JXF_ERROR_GENERIC);
         jxf_TFree(info);
         jxf_TFree(recv_buf);
         jxf_TFree(partitioning);
         jxf_TFree(vec);
         return jxf_error_flag;
      }
      else
	 partitioning[i+1] = recv_buf[i2+2];
   }
   i2 = (num_procs-1)*2;
   partitioning[num_procs] = recv_buf[i2+1]+1;

   jxf_TFree(info);
   jxf_TFree(recv_buf);

   jxf_IJVectorGlobalFirstRow(vec) = partitioning[0];
   jxf_IJVectorGlobalNumRows(vec)= partitioning[num_procs]-partitioning[0];

#endif

   jxf_IJVectorComm(vec)         = comm;
   jxf_IJVectorPartitioning(vec) = partitioning;
   jxf_IJVectorObjectType(vec)   = JXF_UNITIALIZED;
   jxf_IJVectorObject(vec)       = NULL;
   jxf_IJVectorTranslator(vec)   = NULL;
   jxf_IJVectorAssumedPart(vec)  = NULL;
   jxf_IJVectorPrintLevel(vec)   = 0;

  *vector = (JXF_IJVector) vec;
   
   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorCreatePar( jxf_IJVector *vector, JXF_Int *IJpartitioning )
{
   MPI_Comm comm = jxf_IJVectorComm(vector);
   JXF_Int num_procs, jmin, global_n, *partitioning, j;
   jxf_MPI_Comm_size(comm, &num_procs);

#ifdef JXF_NO_GLOBAL_PARTITION
   jmin = jxf_IJVectorGlobalFirstRow(vector);
   global_n = jxf_IJVectorGlobalNumRows(vector);

   partitioning = jxf_CTAlloc(JXF_Int, 2); 

   /* Shift to zero-based partitioning for ParVector object */
   for (j = 0; j < 2; j++) 
      partitioning[j] = IJpartitioning[j] - jmin;

#else
   jmin = IJpartitioning[0];
   global_n = IJpartitioning[num_procs] - jmin;

   partitioning = jxf_CTAlloc(JXF_Int, num_procs+1); 

   /* Shift to zero-based partitioning for ParVector object */
   for (j = 0; j < num_procs+1; j++) 
      partitioning[j] = IJpartitioning[j] - jmin;

#endif

   jxf_IJVectorObject(vector) = jxf_ParVectorCreate(comm, global_n, (JXF_Int *) partitioning);

   return jxf_error_flag;
}

JXF_Int
jxf_AuxParVectorCreate( jxf_AuxParVector **aux_vector )
{
   jxf_AuxParVector  *vector;

   vector = jxf_CTAlloc(jxf_AuxParVector, 1);

   /* set defaults */
   jxf_AuxParVectorMaxOffProcElmts(vector) = 0;
   jxf_AuxParVectorCurrentNumElmts(vector) = 0;
   /* stash for setting or adding off processor values */
   jxf_AuxParVectorOffProcI(vector) = NULL;
   jxf_AuxParVectorOffProcData(vector) = NULL;

  *aux_vector = vector;
   return 0;
}

JXF_Int
JXF_IJVectorInitialize( JXF_IJVector vector )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (jxf_IJVectorObjectType(vec) == JXF_HPCSR)
   {
      if (!jxf_IJVectorObject(vec)) jxf_IJVectorCreatePar(vec, jxf_IJVectorPartitioning(vec));
      jxf_IJVectorInitializePar(vec);
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorInitializePar( jxf_IJVector *vector )
{
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   jxf_AuxParVector *aux_vector = jxf_IJVectorTranslator(vector);
   JXF_Int *partitioning = jxf_ParVectorPartitioning(par_vector);
   jxf_Vector *local_vector = jxf_ParVectorLocalVector(par_vector);
   JXF_Int my_id;
   JXF_Int print_level = jxf_IJVectorPrintLevel(vector);

   MPI_Comm  comm = jxf_IJVectorComm(vector);
   jxf_MPI_Comm_rank(comm, &my_id);
   if (!partitioning)
   {
      if (print_level)
      {
         jxf_printf("No ParVector partitioning for initialization -- ");
         jxf_printf("jxf_IJVectorInitializePar\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   jxf_VectorSize(local_vector) = partitioning[1] - partitioning[0];
#else
   jxf_VectorSize(local_vector) = partitioning[my_id+1] - partitioning[my_id];
#endif

   jxf_ParVectorInitialize(par_vector);
   if (!aux_vector)
   {
      jxf_AuxParVectorCreate(&aux_vector);
      jxf_IJVectorTranslator(vector) = aux_vector;
   }
   jxf_AuxParVectorInitialize(aux_vector);

   return jxf_error_flag;
}

JXF_Int
jxf_AuxParVectorInitialize( jxf_AuxParVector *vector )
{
   JXF_Int max_off_proc_elmts = jxf_AuxParVectorMaxOffProcElmts(vector);
   /* allocate stash for setting or adding off processor values */
   if (max_off_proc_elmts > 0)
   {
      jxf_AuxParVectorOffProcI(vector) = jxf_CTAlloc(JXF_Int, max_off_proc_elmts);
      jxf_AuxParVectorOffProcData(vector) = jxf_CTAlloc(JXF_Real, max_off_proc_elmts);
   }
   return 0;
}

JXF_Int
JXF_IJVectorDestroy( JXF_IJVector vector )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (jxf_IJVectorPartitioning(vec))
      jxf_TFree(jxf_IJVectorPartitioning(vec));
   if (jxf_IJVectorAssumedPart(vec))
      jxf_AssumedPartitionDestroy(jxf_IJVectorAssumedPart(vec));
   if ( jxf_IJVectorObjectType(vec) == JXF_HPCSR)
   {
      jxf_IJVectorDestroyPar(vec);
      if (jxf_IJVectorTranslator(vec))
      {
         jxf_AuxParVectorDestroy((jxf_AuxParVector *) (jxf_IJVectorTranslator(vec)));
      }
   }
   else if (jxf_IJVectorObjectType(vec) != -1)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_TFree(vec);

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorDestroyPar( jxf_IJVector *vector )
{
   return jxf_ParVectorDestroy(jxf_IJVectorObject(vector));
}

JXF_Int
jxf_AuxParVectorDestroy( jxf_AuxParVector *vector )
{
   JXF_Int ierr=0;
   if (vector)
   {
      if (jxf_AuxParVectorOffProcI(vector))
         jxf_TFree(jxf_AuxParVectorOffProcI(vector));
      if (jxf_AuxParVectorOffProcData(vector))
         jxf_TFree(jxf_AuxParVectorOffProcData(vector));
      jxf_TFree(vector);
   }

   return ierr;
}

JXF_Int
JXF_IJVectorSetObjectType( JXF_IJVector vector, JXF_Int type )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_IJVectorObjectType(vec) = type;

   return jxf_error_flag;
}

JXF_Int
JXF_IJVectorSetValues( JXF_IJVector   vector,
                      JXF_Int           nvalues,
                      const JXF_Int    *indices,
                      const JXF_Real *values )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;

   if (nvalues == 0) return jxf_error_flag;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (nvalues < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   if (!values)
   {
      jxf_error_in_arg(4);
      return jxf_error_flag;
   }
   if (jxf_IJVectorObjectType(vec) == JXF_HPCSR)
   {
      return(jxf_IJVectorSetValuesPar(vec, nvalues, indices, values));
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorSetValuesPar( jxf_IJVector  *vector,
                         JXF_Int           num_values,
                         const JXF_Int    *indices,
                         const JXF_Real *values )
{
   JXF_Int my_id;
   JXF_Int i, j, vec_start, vec_stop;
   JXF_Real *data;
   JXF_Int print_level = jxf_IJVectorPrintLevel(vector);
   JXF_Int *IJpartitioning = jxf_IJVectorPartitioning(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   jxf_AuxParVector *aux_vector = jxf_IJVectorTranslator(vector);
   MPI_Comm comm = jxf_IJVectorComm(vector);
   jxf_Vector *local_vector;

   /* If no components are to be set, perform no checking and return */
   if (num_values < 1) return 0;

   jxf_MPI_Comm_rank(comm, &my_id);

   /* If par_vector == NULL or partitioning == NULL or local_vector == NULL 
      let user know of catastrophe and exit */
   if (!par_vector)
   {
      if (print_level)
      {
         jxf_printf("par_vector == NULL -- ");
         jxf_printf("jxf_IJVectorSetValuesPar\n");
         jxf_printf("**** Vector storage is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   local_vector = jxf_ParVectorLocalVector(par_vector);
   if (!IJpartitioning)
   {
      if (print_level)
      {
         jxf_printf("IJpartitioning == NULL -- ");
         jxf_printf("jxf_IJVectorSetValuesPar\n");
         jxf_printf("**** IJVector partitioning is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (!local_vector)
   {
      if (print_level)
      {
         jxf_printf("local_vector == NULL -- ");
         jxf_printf("jxf_IJVectorSetValuesPar\n");
         jxf_printf("**** Vector local data is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   vec_start = IJpartitioning[0];
   vec_stop  = IJpartitioning[1]-1;
#else
   vec_start = IJpartitioning[my_id];
   vec_stop  = IJpartitioning[my_id+1]-1;
#endif

   if (vec_start > vec_stop) 
   {
      if (print_level)
      {
         jxf_printf("vec_start > vec_stop -- ");
         jxf_printf("jxf_IJVectorSetValuesPar\n");
         jxf_printf("**** This vector partitioning should not occur ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   /* Determine whether indices points to local indices only, and if not, store
      indices and values in auxiliary vector structure.  If indices == NULL,
      assume that num_values components are to be set in a block starting at
      vec_start.  NOTE: If indices == NULL off proc values are ignored!!! */

   data = jxf_VectorData(local_vector);

   if (indices)
   {
      JXF_Int current_num_elmts = jxf_AuxParVectorCurrentNumElmts(aux_vector);
      JXF_Int *off_proc_i = jxf_AuxParVectorOffProcI(aux_vector);
      JXF_Int cancel_indx = jxf_AuxParVectorCancelIndx(aux_vector);
      JXF_Int ii;
      for (j = 0; j < num_values; j++)
      {
	 i = indices[j];
	 if (i < vec_start || i > vec_stop)
         {
            for (ii = 0; ii < current_num_elmts; ii++)
	    {
	       if (i == off_proc_i[ii])
	       {
		  off_proc_i[ii] = -1;
		  cancel_indx++;
               }
            }
            jxf_AuxParVectorCancelIndx(aux_vector) = cancel_indx;
         }
         else /* local values are inserted into the vector */
         {
            i -= vec_start;
            data[i] = values[j];
         }
      } 
   }
   else
   {
      if (num_values > vec_stop - vec_start + 1)
      {
         if (print_level)
         {
            jxf_printf("Warning! Indices beyond local range  not identified!\n ");
            jxf_printf("Off processor values have been ignored!\n");
         }
	 num_values = vec_stop - vec_start +1;
      }
#if JXF_USING_OPENMP
#pragma omp parallel for private(j) schedule(static)
#endif
      for (j = 0; j < num_values; j++) data[j] = values[j];
   } 
  
   return jxf_error_flag;
}

JXF_Int
JXF_IJVectorAssemble( JXF_IJVector vector )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (jxf_IJVectorObjectType(vec) == JXF_HPCSR)
   {
      return(jxf_IJVectorAssemblePar(vec));
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorAssemblePar( jxf_IJVector *vector )
{
   JXF_Int *IJpartitioning = jxf_IJVectorPartitioning(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   jxf_AuxParVector *aux_vector = jxf_IJVectorTranslator(vector);
   JXF_Int *partitioning;
   MPI_Comm comm = jxf_IJVectorComm(vector);
   JXF_Int print_level = jxf_IJVectorPrintLevel(vector);

   if (!par_vector)
   {
      if (print_level)
      {
         jxf_printf("par_vector == NULL -- ");
         jxf_printf("jxf_IJVectorAssemblePar\n");
         jxf_printf("**** Vector storage is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
   } 
   partitioning = jxf_ParVectorPartitioning(par_vector);
   if (!IJpartitioning)
   { 
      if (print_level)
      {
         jxf_printf("IJpartitioning == NULL -- ");
         jxf_printf("jxf_IJVectorAssemblePar\n");
         jxf_printf("**** IJVector partitioning is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
   }
   if (!partitioning)
   { 
      if (print_level)
      {
         jxf_printf("partitioning == NULL -- ");
         jxf_printf("jxf_IJVectorAssemblePar\n");
         jxf_printf("**** ParVector partitioning is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
   }

   if (aux_vector)
   {
      JXF_Int off_proc_elmts, current_num_elmts;
      JXF_Int max_off_proc_elmts;
      JXF_Int *off_proc_i;
      JXF_Real *off_proc_data;
      JXF_Int cancel_indx = jxf_AuxParVectorCancelIndx(aux_vector);
      JXF_Int current_i, ii;
      current_num_elmts = jxf_AuxParVectorCurrentNumElmts(aux_vector);
      if (cancel_indx)
      {
         off_proc_i=jxf_AuxParVectorOffProcI(aux_vector);
         off_proc_data=jxf_AuxParVectorOffProcData(aux_vector);
         current_i = 0;
	 for (ii=0; ii < current_num_elmts; ii++) 
         {
            if (off_proc_i[ii] != -1)
	    {
	       off_proc_i[current_i] = off_proc_i[ii];
	       off_proc_data[current_i++] = off_proc_data[ii];
	    }
         }
         jxf_AuxParVectorCurrentNumElmts(aux_vector) = current_i;
         current_num_elmts = current_i;
      }
      jxf_MPI_Allreduce(&current_num_elmts,&off_proc_elmts,1,JXF_MPI_INT,MPI_SUM,comm);
      if (off_proc_elmts)
      {
         max_off_proc_elmts=jxf_AuxParVectorMaxOffProcElmts(aux_vector);
         off_proc_i=jxf_AuxParVectorOffProcI(aux_vector);
         off_proc_data=jxf_AuxParVectorOffProcData(aux_vector);
         jxf_IJVectorAssembleOffProcValsPar(vector, max_off_proc_elmts, 
                                  current_num_elmts, off_proc_i, off_proc_data);
	 jxf_TFree(jxf_AuxParVectorOffProcI(aux_vector));
	 jxf_TFree(jxf_AuxParVectorOffProcData(aux_vector));
	 jxf_AuxParVectorMaxOffProcElmts(aux_vector) = 0;
	 jxf_AuxParVectorCurrentNumElmts(aux_vector) = 0;
      }
   }

   return jxf_error_flag;
}

#ifndef JXF_NO_GLOBAL_PARTITION

JXF_Int
jxf_IJVectorAssembleOffProcValsPar( jxf_IJVector *vector, 
   				   JXF_Int          max_off_proc_elmts,
   				   JXF_Int          current_num_elmts,
   				   JXF_Int         *off_proc_i,
   			     	   JXF_Real      *off_proc_data )
{
   MPI_Comm comm = jxf_IJVectorComm(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   MPI_Request *requests = NULL;
   MPI_Status *status = NULL;
   JXF_Int i, j, j2, row;
   JXF_Int iii, indx, ip, first_index;
   JXF_Int proc_id, num_procs, my_id;
   JXF_Int num_sends, num_sends2;
   JXF_Int num_recvs;
   JXF_Int num_requests;
   JXF_Int vec_start, vec_len;
   JXF_Int *send_procs;
   JXF_Int *send_i;
   JXF_Int *send_map_starts;
   JXF_Int *recv_procs;
   JXF_Int *recv_i;
   JXF_Int *recv_vec_starts;
   JXF_Int *info;
   JXF_Int *int_buffer;
   JXF_Int *proc_id_mem;
   JXF_Int *partitioning;
   JXF_Int *displs;
   JXF_Int *recv_buf;
   JXF_Real *send_data;
   JXF_Real *recv_data;
   JXF_Real *data = jxf_VectorData(jxf_ParVectorLocalVector(par_vector));

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   partitioning = jxf_IJVectorPartitioning(vector);
   first_index = partitioning[my_id];
   info = jxf_CTAlloc(JXF_Int,num_procs);  
   proc_id_mem = jxf_CTAlloc(JXF_Int,current_num_elmts);
   for (i=0; i < current_num_elmts; i++)
   {
      row = off_proc_i[i];
      proc_id = jxf_FindProc(partitioning,row,num_procs);
      proc_id_mem[i] = proc_id; 
      info[proc_id]++;
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
   num_sends2 = 2*num_sends;
   send_procs = jxf_CTAlloc(JXF_Int,num_sends);
   send_map_starts = jxf_CTAlloc(JXF_Int,num_sends+1);
   int_buffer = jxf_CTAlloc(JXF_Int,num_sends2);
   j = 0;
   j2 = 0;
   send_map_starts[0] = 0;
   for (i=0; i < num_procs; i++)
   {
      if (info[i])
      {
         send_procs[j++] = i;
         send_map_starts[j] = send_map_starts[j-1]+info[i];
         int_buffer[j2++] = i;
	 int_buffer[j2++] = info[i];
      }
   }

   jxf_MPI_Allgather(&num_sends2,1,JXF_MPI_INT,info,1,JXF_MPI_INT,comm);

   displs = jxf_CTAlloc(JXF_Int, num_procs+1);
   displs[0] = 0;
   for (i=1; i < num_procs+1; i++)
      displs[i] = displs[i-1]+info[i-1];
   recv_buf = jxf_CTAlloc(JXF_Int, displs[num_procs]);

   jxf_MPI_Allgatherv(int_buffer,num_sends2,JXF_MPI_INT,recv_buf,info,displs,JXF_MPI_INT,comm);

   jxf_TFree(int_buffer);
   jxf_TFree(info);

   /* determine recv procs and amount of data to be received */
   num_recvs = 0;
   for (j=0; j < displs[num_procs]; j+=2)
   {
      if (recv_buf[j] == my_id)
	 num_recvs++;
   }

   recv_procs = jxf_CTAlloc(JXF_Int,num_recvs);
   recv_vec_starts = jxf_CTAlloc(JXF_Int,num_recvs+1);

   j2 = 0;
   recv_vec_starts[0] = 0;
   for (i=0; i < num_procs; i++)
   {
      for (j=displs[i]; j < displs[i+1]; j+=2)
      {
         if (recv_buf[j] == my_id)
         {
	    recv_procs[j2++] = i;
	    recv_vec_starts[j2] = recv_vec_starts[j2-1]+recv_buf[j+1];
         }
         if (j2 == num_recvs) break;
      }
   }
   jxf_TFree(recv_buf);
   jxf_TFree(displs);

   /* set up data to be sent to send procs */
   /* send_i contains for each send proc 
      indices, send_data contains corresponding values */
   send_i = jxf_CTAlloc(JXF_Int,send_map_starts[num_sends]);
   send_data = jxf_CTAlloc(JXF_Real,send_map_starts[num_sends]);
   recv_i = jxf_CTAlloc(JXF_Int,recv_vec_starts[num_recvs]);
   recv_data = jxf_CTAlloc(JXF_Real,recv_vec_starts[num_recvs]);

   for (i=0; i < current_num_elmts; i++)
   {
      proc_id = proc_id_mem[i];
      indx = jxf_BinarySearch(send_procs,proc_id,num_sends);
      iii = send_map_starts[indx];
      send_i[iii] = off_proc_i[i]; 
      send_data[iii] = off_proc_data[i];
      send_map_starts[indx]++;
   }

   jxf_TFree(proc_id_mem);

   for (i=num_sends; i > 0; i--)
   {
      send_map_starts[i] = send_map_starts[i-1];
   }
   send_map_starts[0] = 0;

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
      vec_start = recv_vec_starts[i];
      vec_len = recv_vec_starts[i+1] - vec_start;
      ip = recv_procs[i];
      jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
   }

   for (i=0; i < num_sends; i++)
   {
      vec_start = send_map_starts[i];
      vec_len = send_map_starts[i+1] - vec_start;
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
   jxf_TFree(recv_procs);

   for (i=0; i < recv_vec_starts[num_recvs]; i++)
   {
      row = recv_i[i];
      j = row - first_index;
      data[j] += recv_data[i];
   }

   jxf_TFree(recv_vec_starts);
   jxf_TFree(recv_i);
   jxf_TFree(recv_data);

   return jxf_error_flag;
}

#else

/* assumed partition version */

JXF_Int
jxf_IJVectorAssembleOffProcValsPar( jxf_IJVector *vector, 
   				   JXF_Int          max_off_proc_elmts,
   				   JXF_Int          current_num_elmts,
   				   JXF_Int         *off_proc_i,
   			     	   JXF_Real      *off_proc_data )
{
   JXF_Int myid, global_num_rows;
   JXF_Int global_first_row;
   JXF_Int i, j, in, k;
   JXF_Int proc_id, last_proc, prev_id, tmp_id;
   JXF_Int max_response_size;
   JXF_Int ex_num_contacts = 0;
   JXF_Int range_start, range_end;
   JXF_Int storage;
   JXF_Int indx;
   JXF_Int row, num_ranges, row_count;
   JXF_Int num_recvs;
   JXF_Int counter, upper_bound;
   JXF_Int num_real_procs;
   JXF_Int *row_list=NULL;
   JXF_Int *a_proc_id=NULL, *orig_order=NULL;
   JXF_Int *real_proc_id = NULL, *us_real_proc_id = NULL;
   JXF_Int *ex_contact_procs = NULL, *ex_contact_vec_starts = NULL;
   JXF_Int *recv_starts=NULL;
   JXF_Int *response_buf = NULL, *response_buf_starts=NULL;
   JXF_Int *num_rows_per_proc = NULL;
   JXF_Int  tmp_int;
   JXF_Int  obj_size_bytes, int_size, complex_size;
   JXF_Int  first_index;
   void *void_contact_buf = NULL;
   void *index_ptr;
   void *recv_data_ptr;
   JXF_Real tmp_complex;
   JXF_Int *ex_contact_buf=NULL;
   JXF_Real *vector_data;
   JXF_Real value;
   jxf_DataExchangeResponse response_obj1, response_obj2;
   jxf_ProcListElements send_proc_obj; 

   MPI_Comm comm = jxf_IJVectorComm(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);

   jxf_IJAssumedPart *apart;

   jxf_MPI_Comm_rank(comm, &myid);
   
   global_num_rows = jxf_IJVectorGlobalNumRows(vector);
   global_first_row = jxf_IJVectorGlobalFirstRow(vector);
 
   /* verify that we have created the assumed partition */

   if  (jxf_IJVectorAssumedPart(vector) == NULL)
   {
      jxf_IJVectorCreateAssumedPartition(vector);
   }
   apart = jxf_IJVectorAssumedPart(vector);

   /* get the assumed processor id for each row */
   a_proc_id = jxf_CTAlloc(JXF_Int, current_num_elmts);
   orig_order =  jxf_CTAlloc(JXF_Int, current_num_elmts);
   real_proc_id = jxf_CTAlloc(JXF_Int, current_num_elmts);
   row_list =   jxf_CTAlloc(JXF_Int, current_num_elmts);

   if (current_num_elmts > 0)
   {
      for (i=0; i < current_num_elmts; i++)
      {
         row = off_proc_i[i]; 
         row_list[i] = row;
         jxf_GetAssumedPartitionProcFromRow(comm, row, global_first_row, global_num_rows, &proc_id);
         a_proc_id[i] = proc_id;
         orig_order[i] = i;
      }

      /* now we need to find the actual order of each row  - sort on row -
         this will result in proc ids sorted also...*/
      
      jxf_qsort3i(row_list, a_proc_id, orig_order, 0, current_num_elmts -1);

      /* calculate the number of contacts */
      ex_num_contacts = 1;
      last_proc = a_proc_id[0];
      for (i=1; i < current_num_elmts; i++)
      {
         if (a_proc_id[i] > last_proc)      
         {
            ex_num_contacts++;
            last_proc = a_proc_id[i];
         }
      }
      
   }
   
   /* now we will go through a create a contact list - need to contact
      assumed processors and find out who the actual row owner is - we
      will contact with a range (2 numbers) */

   ex_contact_procs = jxf_CTAlloc(JXF_Int, ex_num_contacts);
   ex_contact_vec_starts =  jxf_CTAlloc(JXF_Int, ex_num_contacts+1);
   ex_contact_buf =  jxf_CTAlloc(JXF_Int, ex_num_contacts*2);

   counter = 0;
   range_end = -1;
   for (i=0; i< current_num_elmts; i++) 
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
         
         jxf_GetAssumedPartitionRowRange(comm, proc_id, global_first_row,
					global_num_rows, &range_start, &range_end);
      }
   }

   /*finish the starts*/
   ex_contact_vec_starts[counter] =  counter*2;
   /*finish the last range*/
   if (counter > 0)  
      ex_contact_buf[counter*2 - 1] = row_list[current_num_elmts - 1];

   /* create response object - can use same fill response as used in the commpkg
      routine */
   response_obj1.fill_response = jxf_RangeFillResponseIJDetermineRecvProcs;
   response_obj1.data1 =  apart; /* this is necessary so we can fill responses*/ 
   response_obj1.data2 = NULL;
   
   max_response_size = 6;  /* 6 means we can fit 3 ranges*/
   
   jxf_DataExchangeList(ex_num_contacts, ex_contact_procs, 
                          ex_contact_buf, ex_contact_vec_starts, sizeof(JXF_Int), 
                          sizeof(JXF_Int), &response_obj1, max_response_size, 4, 
                          comm, (void**) &response_buf, &response_buf_starts);

   /* now response_buf contains a proc_id followed by an upper bound for the
      range.  */

   jxf_TFree(ex_contact_procs);
   jxf_TFree(ex_contact_buf);
   jxf_TFree(ex_contact_vec_starts);

   jxf_TFree(a_proc_id);
   a_proc_id = NULL;

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
      while (j < current_num_elmts && row_list[j] <= upper_bound)     
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
      JXF_Int and JXF_Real data.  (row number and value) - we will send
      everything as a void since we may not know the rel sizes of ints and
      doubles */
 
   /* first find out how many elements to send per proc - so we can do
      storage */
 
   int_size = sizeof(JXF_Int);
   complex_size = sizeof(JXF_Real);
   
   obj_size_bytes = jxf_max(int_size, complex_size);
    
   ex_contact_procs = jxf_CTAlloc(JXF_Int, num_real_procs);
   num_rows_per_proc = jxf_CTAlloc(JXF_Int, num_real_procs);
   
   counter = 0;
   
   if (num_real_procs > 0 )
   {
      ex_contact_procs[0] = real_proc_id[0];
      num_rows_per_proc[0] = 1;

      /* loop through real procs - these are sorted (row_list is sorted also)*/
      for (i=1; i < current_num_elmts; i++)
      {
         if (real_proc_id[i] == ex_contact_procs[counter]) /* same processor */
         {
            num_rows_per_proc[counter] += 1; /*another row */
         }
         else /* new processor */
         {
            counter++;
            ex_contact_procs[counter] = real_proc_id[i];
            num_rows_per_proc[counter] = 1;
         }
      }
   }
     
   /* calculate total storage and make vec_starts arrays */
   storage = 0;
   ex_contact_vec_starts = jxf_CTAlloc(JXF_Int, num_real_procs + 1);
   ex_contact_vec_starts[0] = -1;

   for (i=0; i < num_real_procs; i++)
   {
      storage += 1 + 2*  num_rows_per_proc[i];
      ex_contact_vec_starts[i+1] = -storage-1; /* need negative for next loop */
   }      

   void_contact_buf = jxf_MAlloc(storage*obj_size_bytes);
   index_ptr = void_contact_buf; /* step through with this index */

   /* set up data to be sent to send procs */
   /* for each proc, ex_contact_buf_d contains #rows, row #, data, etc. */
      
   /* un-sort real_proc_id  - we want to access data arrays in order */

   us_real_proc_id =  jxf_CTAlloc(JXF_Int, current_num_elmts);
   for (i=0; i < current_num_elmts; i++)
   {
      us_real_proc_id[orig_order[i]] = real_proc_id[i];
   }
   jxf_TFree(real_proc_id);

   prev_id = -1;
   for (i=0; i < current_num_elmts; i++)
   {
      proc_id = us_real_proc_id[i];
      /* can't use row list[i] - you loose the negative signs that differentiate
         add/set values */
      row = off_proc_i[i];
      /* find position of this processor */
      indx = jxf_BinarySearch(ex_contact_procs, proc_id, num_real_procs);
      in =  ex_contact_vec_starts[indx];

      index_ptr = (void *) ((char *) void_contact_buf + in*obj_size_bytes);

      /* first time for this processor - add the number of rows to the buffer */
      if (in < 0)
      {
         in = -in - 1;
         /* re-calc. index_ptr since in_i was negative */
         index_ptr = (void *) ((char *) void_contact_buf + in*obj_size_bytes);

         tmp_int =  num_rows_per_proc[indx];
         memcpy( index_ptr, &tmp_int, int_size);
         index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);

         in++;
      }
      /* add row # */   
      memcpy( index_ptr, &row, int_size);
      index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
      in++;

      /* add value */
      tmp_complex = off_proc_data[i];
      memcpy( index_ptr, &tmp_complex, complex_size);
      index_ptr = (void *) ((char *) index_ptr + obj_size_bytes);
      in++;
      
      /* increment the indexes to keep track of where we are - fix later */
      ex_contact_vec_starts[indx] = in;
   }
   
   /* some clean up */
 
   jxf_TFree(response_buf);
   jxf_TFree(response_buf_starts);

   jxf_TFree(us_real_proc_id);
   jxf_TFree(orig_order);
   jxf_TFree(row_list);
   jxf_TFree(num_rows_per_proc);

   for (i=num_real_procs; i > 0; i--)
   {
      ex_contact_vec_starts[i] =   ex_contact_vec_starts[i-1];
   }

   ex_contact_vec_starts[0] = 0;

   /* now send the data */

   /***********************************/
   /* now get the info in send_proc_obj_d */

   /* the response we expect is just a confirmation*/
   response_buf = NULL;
   response_buf_starts = NULL;

   /*build the response object*/

   /* use the send_proc_obj for the info kept from contacts */
   /*estimate inital storage allocation */

   send_proc_obj.length = 0;
   send_proc_obj.storage_length = num_real_procs + 5;
   send_proc_obj.id = NULL; /* don't care who sent it to us */
   send_proc_obj.vec_starts =
      jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length + 1); 
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = storage + 20;
   send_proc_obj.v_elements =
      jxf_MAlloc(obj_size_bytes*send_proc_obj.element_storage_length);

   response_obj2.fill_response = jxf_FillResponseIJOffProcVals;
   response_obj2.data1 = NULL;
   response_obj2.data2 = &send_proc_obj;

   max_response_size = 0;

   jxf_DataExchangeList(num_real_procs, ex_contact_procs, 
                          void_contact_buf, ex_contact_vec_starts, obj_size_bytes,
                          0, &response_obj2, max_response_size, 5, 
                          comm,  (void **) &response_buf, &response_buf_starts);

   /***********************************/

   jxf_TFree(response_buf);
   jxf_TFree(response_buf_starts);

   jxf_TFree(ex_contact_procs);
   jxf_TFree(void_contact_buf);
   jxf_TFree(ex_contact_vec_starts);

   /* Now we can unpack the send_proc_objects and either set or add to the
      vector data */

   num_recvs = send_proc_obj.length; 

   /* alias */
   recv_data_ptr = send_proc_obj.v_elements;
   recv_starts = send_proc_obj.vec_starts;
   
   vector_data = jxf_VectorData(jxf_ParVectorLocalVector(par_vector));
   first_index =  jxf_ParVectorFirstIndex(par_vector);

   for (i=0; i < num_recvs; i++)
   {
      indx = recv_starts[i];

      /* get the number of rows for  this recv */
      memcpy( &row_count, recv_data_ptr, int_size);
      recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
      indx++;

      for (j=0; j < row_count; j++) /* for each row: unpack info */
      {
         /* row # */
         memcpy( &row, recv_data_ptr, int_size);
         recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
         indx++;

         /* value */
         memcpy( &value, recv_data_ptr, complex_size);
         recv_data_ptr = (void *) ((char *)recv_data_ptr + obj_size_bytes);
         indx++;

         k = row - first_index - global_first_row;
         vector_data[k] += value;
      }
   }
   
   jxf_TFree(send_proc_obj.v_elements);
   jxf_TFree(send_proc_obj.vec_starts);
 
   return jxf_error_flag;
}

#endif

JXF_Int
jxf_IJVectorCreateAssumedPartition( jxf_IJVector *vector )
{
   JXF_Int global_num, global_first_row;
   JXF_Int myid;
   JXF_Int start=0, end=0;
   JXF_Int *partitioning = jxf_IJVectorPartitioning(vector);
   MPI_Comm   comm;
   jxf_IJAssumedPart *apart;

   global_num = jxf_IJVectorGlobalNumRows(vector);
   global_first_row = jxf_IJVectorGlobalFirstRow(vector);
   comm = jxf_ParVectorComm(vector);

   /* find out my actualy range of rows */
   start =  partitioning[0];
   end = partitioning[1]-1;

   jxf_MPI_Comm_rank(comm, &myid);

   /* allocate space */
   apart = jxf_CTAlloc(jxf_IJAssumedPart, 1);

  /* get my assumed partitioning  - we want partitioning of the vector that the
      matrix multiplies - so we use the col start and end */
   jxf_GetAssumedPartitionRowRange(comm, myid, global_first_row, 
				global_num, &(apart->row_start), &(apart->row_end));

  /*allocate some space for the partition of the assumed partition */
   apart->length = 0;
  /*room for 10 owners of the assumed partition*/ 
   apart->storage_length = 10; /*need to be >=1 */ 
   apart->proc_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_start_list = jxf_TAlloc(JXF_Int, apart->storage_length);
   apart->row_end_list = jxf_TAlloc(JXF_Int, apart->storage_length);

  /* now we want to reconcile our actual partition with the assumed partition */
   jxf_LocateAssummedPartition(comm, start, end, global_first_row, global_num, apart, myid);

  /* this partition will be saved in the vector data structure until the vector is destroyed */
   jxf_IJVectorAssumedPart(vector) = apart;

   return jxf_error_flag;
}

JXF_Int
jxf_FindProc( JXF_Int *list, JXF_Int value, JXF_Int list_length )
{
   JXF_Int low, high, m;

   low = 0;
   high = list_length;
   if (value >= list[high] || value < list[low])
      return -1;
   else
   {
      while (low+1 < high)
      {
         m = (low + high) / 2;
         if (value < list[m])
         {
            high = m;
         }
         else if (value >= list[m])
         {
            low = m;
         }
      }
      return low;
   }
}

JXF_Int
JXF_IJVectorGetValues( JXF_IJVector  vector,
                      JXF_Int          nvalues,
                      const JXF_Int   *indices,
                      JXF_Real      *values )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (nvalues == 0) return jxf_error_flag;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (nvalues < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   if (!values)
   {
      jxf_error_in_arg(4);
      return jxf_error_flag;
   }
   if (jxf_IJVectorObjectType(vec) == JXF_HPCSR)
   {
      return(jxf_IJVectorGetValuesPar(vec, nvalues, indices, values));
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorGetValuesPar( jxf_IJVector *vector,
                         JXF_Int          num_values,
                         const JXF_Int   *indices,
                         JXF_Real      *values )
{
   JXF_Int my_id;
   JXF_Int i, j, vec_start, vec_stop;
   JXF_Real *data;
   JXF_Int ierr = 0;
   JXF_Int *IJpartitioning = jxf_IJVectorPartitioning(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   MPI_Comm comm = jxf_IJVectorComm(vector);
   jxf_Vector *local_vector;
   JXF_Int print_level = jxf_IJVectorPrintLevel(vector);

   /* If no components are to be retrieved, perform no checking and return */
   if (num_values < 1) return 0;

   jxf_MPI_Comm_rank(comm, &my_id);

   /* If par_vector == NULL or partitioning == NULL or local_vector == NULL 
      let user know of catastrophe and exit */

   if (!par_vector)
   {
      if (print_level)
      {
         jxf_printf("par_vector == NULL -- ");
         jxf_printf("jxf_IJVectorGetValuesPar\n");
         jxf_printf("**** Vector storage is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   local_vector = jxf_ParVectorLocalVector(par_vector);
   if (!IJpartitioning)
   {
      if (print_level)
      {
         jxf_printf("IJpartitioning == NULL -- ");
         jxf_printf("jxf_IJVectorGetValuesPar\n");
         jxf_printf("**** IJVector partitioning is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (!local_vector)
   {
      if (print_level)
      {
         jxf_printf("local_vector == NULL -- ");
         jxf_printf("jxf_IJVectorGetValuesPar\n");
         jxf_printf("**** Vector local data is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   vec_start = IJpartitioning[0];
   vec_stop  = IJpartitioning[1];
#else
   vec_start = IJpartitioning[my_id];
   vec_stop  = IJpartitioning[my_id+1];
#endif

   if (vec_start > vec_stop) 
   {
      if (print_level)
      {
         jxf_printf("vec_start > vec_stop -- ");
         jxf_printf("jxf_IJVectorGetValuesPar\n");
         jxf_printf("**** This vector partitioning should not occur ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   /* Determine whether indices points to local indices only, and if not, let
      user know of catastrophe and exit.  If indices == NULL, assume that
      num_values components are to be retrieved from block starting at
      vec_start */

   if (indices)
   {
      for (i = 0; i < num_values; i++)
      { 	
         ierr += (indices[i] <  vec_start);
         ierr += (indices[i] >= vec_stop);
      }
   }

   if (ierr)
   {
      if (print_level)
      {
         jxf_printf("indices beyond local range -- ");
         jxf_printf("jxf_IJVectorGetValuesPar\n");
         jxf_printf("**** Indices specified are unusable ****\n");
      }
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

   data = jxf_VectorData(local_vector);

   if (indices)
   {
#if JXF_USING_OPENMP
#pragma omp parallel for private(i,j) schedule(static)
#endif
      for (j = 0; j < num_values; j++)
      {
         i = indices[j] - vec_start;
         values[j] = data[i];
      }
   }
   else
   {
     if (num_values > (vec_stop-vec_start))
     {
        jxf_error_in_arg(2);
        return jxf_error_flag;
     }
#if JXF_USING_OPENMP
#pragma omp parallel for private(j) schedule(static)
#endif
      for (j = 0; j < num_values; j++)
         values[j] = data[j];
   }

   return jxf_error_flag;
}

JXF_Int
JXF_IJVectorAddToValues( JXF_IJVector    vector,
                        JXF_Int            nvalues,
                        const JXF_Int     *indices,
                        const JXF_Real  *values )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;
   if (nvalues == 0) return jxf_error_flag;
   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (nvalues < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   if (!values)
   {
      jxf_error_in_arg(4);
      return jxf_error_flag;
   }
   if ( jxf_IJVectorObjectType(vec) == JXF_HPCSR )
   {
      return( jxf_IJVectorAddToValuesPar(vec, nvalues, indices, values) );
   }
   else
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

JXF_Int
JXF_IJVectorRead( const char  *filename,
                 MPI_Comm     comm,
                 JXF_Int          type,
                 JXF_IJVector *vector_ptr )
{
   JXF_IJVector  vector;
   JXF_Int       jlower, jupper, j;
   JXF_Real    value;
   JXF_Int       myid, ret;
   char      new_filename[255];
   FILE     *file;

   jxf_MPI_Comm_rank(comm, &myid);
   jxf_sprintf(new_filename,"%s.%05d", filename, myid);
   if ((file = fopen(new_filename, "r")) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_fscanf(file, "%d %d", &jlower, &jupper);
   JXF_IJVectorCreate(comm, jlower, jupper, &vector);
   JXF_IJVectorSetObjectType(vector, type);
   JXF_IJVectorInitialize(vector);
   /* It is important to ensure that whitespace follows the index value to help
    * catch mistakes in the input file.  This is done with %*[ \t].  Using a
    * space here causes an input line with a single decimal value on it to be
    * read as if it were an integer followed by a decimal value. */
   while ( (ret = jxf_fscanf(file, "%d%*[ \t]%le", &j, &value)) != EOF )
   {
      if (ret != 2)
      {
         jxf_error_w_msg(JXF_ERROR_GENERIC, "Error in IJ vector input file.");
         return jxf_error_flag;
      }
      if (j < jlower || j > jupper)
	 JXF_IJVectorAddToValues(vector, 1, &j, &value);
      else
	 JXF_IJVectorSetValues(vector, 1, &j, &value);
   }
   JXF_IJVectorAssemble(vector);
   fclose(file);
  *vector_ptr = vector;

   return jxf_error_flag;
}

JXF_Int
JXF_IJVectorPrint( JXF_IJVector vector, const char *filename )
{
   MPI_Comm  comm;
   JXF_Int      *partitioning;
   JXF_Int       jlower, jupper, j;
   JXF_Real    value;
   JXF_Int       myid;
   char      new_filename[255];
   FILE     *file;

   if (!vector)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   comm = jxf_IJVectorComm(vector);
   jxf_MPI_Comm_rank(comm, &myid);
   jxf_sprintf(new_filename,"%s.%05d", filename, myid);
   if ((file = fopen(new_filename, "w")) == NULL)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   partitioning = jxf_IJVectorPartitioning(vector);
#ifdef JXF_NO_GLOBAL_PARTITION
   jlower = partitioning[0];
   jupper = partitioning[1] - 1;
#else
   jlower = partitioning[myid];
   jupper = partitioning[myid+1] - 1;
#endif
   jxf_fprintf(file, "%d %d\n", jlower, jupper);
   for (j = jlower; j <= jupper; j++)
   {
      JXF_IJVectorGetValues(vector, 1, &j, &value);
      jxf_fprintf(file, "%d %.14e\n", j, value);
   }
   fclose(file);

   return jxf_error_flag;
}

JXF_Int
jxf_IJVectorAddToValuesPar( jxf_IJVector   *vector,
                           JXF_Int            num_values,
                           const JXF_Int     *indices,
                           const JXF_Real  *values )
{
   JXF_Int my_id;
   JXF_Int i, j, vec_start, vec_stop;
   JXF_Real *data;
   JXF_Int print_level = jxf_IJVectorPrintLevel(vector);
   JXF_Int *IJpartitioning = jxf_IJVectorPartitioning(vector);
   jxf_ParVector *par_vector = jxf_IJVectorObject(vector);
   jxf_AuxParVector *aux_vector = jxf_IJVectorTranslator(vector);
   MPI_Comm comm = jxf_IJVectorComm(vector);
   jxf_Vector *local_vector;

   /* If no components are to be retrieved, perform no checking and return */
   if (num_values < 1) return 0;

   jxf_MPI_Comm_rank(comm, &my_id);

   /* If par_vector == NULL or partitioning == NULL or local_vector == NULL 
      let user know of catastrophe and exit */

   if (!par_vector)
   {
      if (print_level)
      {
         jxf_printf("par_vector == NULL -- ");
         jxf_printf("jxf_IJVectorAddToValuesPar\n");
         jxf_printf("**** Vector storage is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   local_vector = jxf_ParVectorLocalVector(par_vector);
   if (!IJpartitioning)
   {
      if (print_level)
      {
         jxf_printf("IJpartitioning == NULL -- ");
         jxf_printf("jxf_IJVectorAddToValuesPar\n");
         jxf_printf("**** IJVector partitioning is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (!local_vector)
   {
      if (print_level)
      {
         jxf_printf("local_vector == NULL -- ");
         jxf_printf("jxf_IJVectorAddToValuesPar\n");
         jxf_printf("**** Vector local data is either unallocated or orphaned ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

#ifdef JXF_NO_GLOBAL_PARTITION
   vec_start = IJpartitioning[0];
   vec_stop  = IJpartitioning[1]-1;
#else
   vec_start = IJpartitioning[my_id];
   vec_stop  = IJpartitioning[my_id+1]-1;
#endif

   if (vec_start > vec_stop) 
   {
      if (print_level)
      {
         jxf_printf("vec_start > vec_stop -- ");
         jxf_printf("jxf_IJVectorAddToValuesPar\n");
         jxf_printf("**** This vector partitioning should not occur ****\n");
      }
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   data = jxf_VectorData(local_vector);

   if (indices)
   {
      JXF_Int current_num_elmts = jxf_AuxParVectorCurrentNumElmts(aux_vector);
      JXF_Int max_off_proc_elmts = jxf_AuxParVectorMaxOffProcElmts(aux_vector);
      JXF_Int *off_proc_i = jxf_AuxParVectorOffProcI(aux_vector);
      JXF_Real *off_proc_data = jxf_AuxParVectorOffProcData(aux_vector);

      for (j = 0; j < num_values; j++)
      {
	 i = indices[j];
	 if (i < vec_start || i > vec_stop)
         {
            /* if elements outside processor boundaries, store in off processor
               stash */
	    if (!max_off_proc_elmts)
            {
               max_off_proc_elmts = 100;
               jxf_AuxParVectorMaxOffProcElmts(aux_vector) = max_off_proc_elmts;
               jxf_AuxParVectorOffProcI(aux_vector) = jxf_CTAlloc(JXF_Int,max_off_proc_elmts);
               jxf_AuxParVectorOffProcData(aux_vector) = jxf_CTAlloc(JXF_Real,max_off_proc_elmts);
               off_proc_i = jxf_AuxParVectorOffProcI(aux_vector);
               off_proc_data = jxf_AuxParVectorOffProcData(aux_vector);
            }
            else if (current_num_elmts + 1 > max_off_proc_elmts)
            {
               max_off_proc_elmts += 10;
               off_proc_i = jxf_TReAlloc(off_proc_i,JXF_Int,max_off_proc_elmts);
               off_proc_data = jxf_TReAlloc(off_proc_data,JXF_Real,max_off_proc_elmts);
               jxf_AuxParVectorMaxOffProcElmts(aux_vector) = max_off_proc_elmts;
               jxf_AuxParVectorOffProcI(aux_vector) = off_proc_i;
               jxf_AuxParVectorOffProcData(aux_vector) = off_proc_data;
            }
            off_proc_i[current_num_elmts] = i;
            off_proc_data[current_num_elmts++] = values[j];
            jxf_AuxParVectorCurrentNumElmts(aux_vector)=current_num_elmts;
         }
         else /* local values are added to the vector */
         {
            i -= vec_start;
            data[i] += values[j];
         }
      } 
   }
   else 
   {
      if (num_values > vec_stop - vec_start + 1)
      {
         if (print_level)
         {
            jxf_printf("Warning! Indices beyond local range  not identified!\n ");
            jxf_printf("Off processor values have been ignored!\n");
         }
	 num_values = vec_stop - vec_start +1;
      }
#if JXF_USING_OPENMP
#pragma omp parallel for private(j) schedule(static)
#endif
      for (j = 0; j < num_values; j++)
         data[j] += values[j];
   } 
  
   return jxf_error_flag;
}

/*--------------------------------------------------------------------------
 * JXF_IJVectorGetObject
 *--------------------------------------------------------------------------*/

JXF_Int
JXF_IJVectorGetObject( JXF_IJVector   vector,
                         void           **object )
{
   jxf_IJVector *vec = (jxf_IJVector *) vector;

   if (!vec)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   *object = jxf_IJVectorObject(vec);

   return jxf_error_flag;
}
