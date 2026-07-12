//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_vector.c -- basic operations for parallel vectors
 *  Date: 2011/09/03
 */

#include "jxf_mv.h"

/*!
 * \fn jxf_ParVector *jxf_BuildRhsParFromOneFile
 * \brief Build a parallel vector from a given file.
 * \param *filename pointer to the file name.
 * \param *A pointer to the matrix which will provide the row partitioning.
 * \note The storage format of the vector is expected to be as follows:
 *         the 1st row:                   n                       // JXF_Int-type
 *         the 2 to (n+1)-th row:         f[i]    i=0(1)n-1       // JXF_Real-type
 *      where n is the total size of the vector, f[i] is the i-th entry of the vector.
 * \author peghoty
 * \date 2009/07/10
 */
jxf_ParVector *
jxf_BuildRhsParFromOneFile(char *filename, jxf_ParCSRMatrix *A, JXF_Int file_type)
{
   jxf_ParVector *b = NULL;
   jxf_Vector *b_CSR = NULL;
   JXF_Int my_id;
   JXF_Int *partitioning;

   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0)
   {
      if (file_type == 2 || file_type == 12)
      {
         b_CSR = jxf_SeqVectorBinaryRead(filename);
      }
      else
      {
         b_CSR = jxf_SeqVectorRead(filename);
      }
   }

   jxf_ParCSRMatrixGetRowPartitioning(A, &partitioning);

   b = jxf_VectorToParVector(MPI_COMM_WORLD, b_CSR, partitioning);

   if (my_id == 0)
      jxf_SeqVectorDestroy(b_CSR);

   return (b);
}

jxf_ParVector *
jxf_BuildRhsParFromOneFile3(char *filename, jxf_ParCSRMatrix *A)
{
   jxf_ParVector *b = NULL;
   jxf_Vector *b_CSR = NULL;
   JXF_Int my_id;
   JXF_Int *partitioning;

   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0)
   {
      b_CSR = jxf_SeqVectorRead2(filename);
   }

   jxf_ParCSRMatrixGetRowPartitioning(A, &partitioning);

   b = jxf_VectorToParVector(MPI_COMM_WORLD, b_CSR, partitioning);

   if (my_id == 0)
      jxf_SeqVectorDestroy(b_CSR);

   return (b);
}

/*!
 * \fn jxf_ParVector *jxf_ParMultiVectorCreate
 * \brief Create a parallel Multivector.
 * \param comm communicator.
 * \param global_size global size of a single vector.
 * \param *partitioning partitioning of the parallel vector.
 * \param num_vector number of vectors.
 * \date 2011/09/03
 */
jxf_ParVector *
jxf_ParMultiVectorCreate(MPI_Comm comm,
                        JXF_Int global_size,
                        JXF_Int *partitioning,
                        JXF_Int num_vectors)
{
   /* note that global_size is the global length of a single vector */
   jxf_ParVector *vector = jxf_ParVectorCreate(comm, global_size, partitioning);
   jxf_ParVectorNumVectors(vector) = num_vectors;

   return vector;
}

/*!
 * \fn jxf_ParVector *jxf_ParVectorCreate
 * \brief Create a parallel vector.
 * \param comm communicator.
 * \param global_size global size of the vector.
 * \param *partitioning partitioning of the parallel vector.
 * \note If create is called for JXF_NO_GLOBAL_PARTITION and partitioning is NOT null,
 *       then it is assumed that it is array of length 2 containing the start row of
 *       the calling processor followed by the start row of the next processor - AHB 6/05
 * \date 2011/09/03
 */
jxf_ParVector *
jxf_ParVectorCreate(MPI_Comm comm,
                   JXF_Int global_size,
                   JXF_Int *partitioning)
{
   jxf_ParVector *vector;
   JXF_Int num_procs, my_id;

   if (global_size < 0)
   {
      jxf_error_in_arg(2);
      return NULL;
   }
   vector = jxf_CTAlloc(jxf_ParVector, 1);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (!partitioning)
   {
      jxf_MPI_Comm_size(comm, &num_procs);
#ifdef JXF_NO_GLOBAL_PARTITION
      jxf_GenerateLocalPartitioning(global_size, num_procs, my_id, &partitioning);
#else
      jxf_GeneratePartitioning(global_size, num_procs, &partitioning);
#endif
   }

   jxf_ParVectorAssumedPartition(vector) = NULL;

   jxf_ParVectorComm(vector) = comm;
   jxf_ParVectorGlobalSize(vector) = global_size;

#ifdef JXF_NO_GLOBAL_PARTITION
   jxf_ParVectorFirstIndex(vector) = partitioning[0];
   jxf_ParVectorLastIndex(vector) = partitioning[1] - 1;
   jxf_ParVectorPartitioning(vector) = partitioning;
   jxf_ParVectorLocalVector(vector) = jxf_SeqVectorCreate(partitioning[1] - partitioning[0]);
#else
   jxf_ParVectorFirstIndex(vector) = partitioning[my_id];
   jxf_ParVectorLastIndex(vector) = partitioning[my_id + 1] - 1;
   jxf_ParVectorPartitioning(vector) = partitioning;
   jxf_ParVectorLocalVector(vector) = jxf_SeqVectorCreate(partitioning[my_id + 1] - partitioning[my_id]);
#endif

   /* set defaults */
   jxf_ParVectorOwnsData(vector) = 1;
   jxf_ParVectorOwnsPartitioning(vector) = 1;

   return vector;
}

/*!
 * \fn JXF_Int jxf_ParVectorInitialize
 * \brief Initialize a parallel vector.
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorInitialize(jxf_ParVector *vector)
{
   if (!vector)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_SeqVectorInitialize(jxf_ParVectorLocalVector(vector));

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParVectorDestroy
 * \brief Destroy a parallel vector.
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorDestroy(jxf_ParVector *vector)
{
   if (vector)
   {
      if (jxf_ParVectorOwnsData(vector))
      {
         jxf_SeqVectorDestroy(jxf_ParVectorLocalVector(vector));
      }

      if (jxf_ParVectorOwnsPartitioning(vector))
      {
         jxf_TFree(jxf_ParVectorPartitioning(vector));
      }

      if (jxf_ParVectorAssumedPartition(vector))
      {
         jxf_ParVectorDestroyAssumedPartition(vector);
      }

      jxf_TFree(vector);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParVectorDestroyAssumedPartition
 * \brief Destroy the assumedpartition of a parallel vector.
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorDestroyAssumedPartition(jxf_ParVector *vector)
{
   jxf_IJAssumedPart *apart;

   apart = jxf_ParVectorAssumedPartition(vector);

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

/*!
 * \fn jxf_Vector *jxf_ParVectorToVectorAll
 * \brief Generates a Vector on every proc which has a piece of the data
 *        from a ParVector on several procs in comm, vec_starts needs to
 *        contain the partitioning across all procs in comm
 * \date 2011/09/03
 */
jxf_Vector *
jxf_ParVectorToVectorAll(jxf_ParVector *par_v)
{
   MPI_Comm comm = jxf_ParVectorComm(par_v);
   JXF_Int global_size = jxf_ParVectorGlobalSize(par_v);
#ifndef JXF_NO_GLOBAL_PARTITION
   JXF_Int *vec_starts = jxf_ParVectorPartitioning(par_v);
#endif
   jxf_Vector *local_vector = jxf_ParVectorLocalVector(par_v);
   JXF_Int num_procs, my_id;
   JXF_Int num_vectors = jxf_ParVectorNumVectors(par_v);
   jxf_Vector *vector;
   JXF_Real *vector_data;
   JXF_Real *local_data;
   JXF_Int local_size;
   MPI_Request *requests;
   MPI_Status *status;
   JXF_Int i, j;
   JXF_Int *used_procs;
   JXF_Int num_types, num_requests;
   JXF_Int vec_len, proc_id;

#ifdef JXF_NO_GLOBAL_PARTITION
   JXF_Int *new_vec_starts;
   JXF_Int num_contacts;
   JXF_Int contact_proc_list[1];
   JXF_Int contact_send_buf[1];
   JXF_Int contact_send_buf_starts[2];
   JXF_Int max_response_size;
   JXF_Int *response_recv_buf = NULL;
   JXF_Int *response_recv_buf_starts = NULL;
   jxf_DataExchangeResponse response_obj;
   jxf_ProcListElements send_proc_obj;

   JXF_Int *send_info = NULL;
   MPI_Status status1;
   JXF_Int count, tag1 = 112, tag2 = 223;
   JXF_Int start;
#endif

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

#ifdef JXF_NO_GLOBAL_PARTITION

   local_size = jxf_ParVectorLastIndex(par_v) - jxf_ParVectorFirstIndex(par_v) + 1;

   /* determine procs which hold data of par_v and store ids in used_procs */
   /* we need to do an exchange data for this.  If I own row then I will contact
      processor 0 with the endpoint of my local range */

   if (local_size > 0)
   {
      num_contacts = 1;
      contact_proc_list[0] = 0;
      contact_send_buf[0] = jxf_ParVectorLastIndex(par_v);
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 1;
   }
   else
   {
      num_contacts = 0;
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 0;
   }

   /*build the response object*/
   /*send_proc_obj will  be for saving info from contacts */
   send_proc_obj.length = 0;
   send_proc_obj.storage_length = 10;
   send_proc_obj.id = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length);
   send_proc_obj.vec_starts = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length + 1);
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = 10;
   send_proc_obj.elements = jxf_CTAlloc(JXF_Int, send_proc_obj.element_storage_length);

   max_response_size = 0; /* each response is null */
   response_obj.fill_response = jxf_FillResponseParToVectorAll;
   response_obj.data1 = NULL;
   response_obj.data2 = &send_proc_obj; /*this is where we keep info from contacts*/

   jxf_DataExchangeList(num_contacts,
                       contact_proc_list, contact_send_buf,
                       contact_send_buf_starts, sizeof(JXF_Int),
                       sizeof(JXF_Int), &response_obj,
                       max_response_size, 1,
                       comm, (void **)&response_recv_buf,
                       &response_recv_buf_starts);

   /* now processor 0 should have a list of ranges for processors that have rows -
      these are in send_proc_obj - it needs to create the new list of processors
      and also an array of vec starts - and send to those who own row */
   if (my_id)
   {
      if (local_size)
      {
         /* look for a message from processor 0 */
         jxf_MPI_Probe(0, tag1, comm, &status1);
         jxf_MPI_Get_count(&status1, JXF_MPI_INT, &count);

         send_info = jxf_CTAlloc(JXF_Int, count);
         jxf_MPI_Recv(send_info, count, JXF_MPI_INT, 0, tag1, comm, &status1);

         /* now unpack */
         num_types = send_info[0];
         used_procs = jxf_CTAlloc(JXF_Int, num_types);
         new_vec_starts = jxf_CTAlloc(JXF_Int, num_types + 1);

         for (i = 1; i <= num_types; i++)
         {
            used_procs[i - 1] = send_info[i];
         }
         for (i = num_types + 1; i < count; i++)
         {
            new_vec_starts[i - num_types - 1] = send_info[i];
         }
      }
      else /* clean up and exit */
      {
         jxf_TFree(send_proc_obj.vec_starts);
         jxf_TFree(send_proc_obj.id);
         jxf_TFree(send_proc_obj.elements);
         if (response_recv_buf)
            jxf_TFree(response_recv_buf);
         if (response_recv_buf_starts)
            jxf_TFree(response_recv_buf_starts);
         return NULL;
      }
   }
   else /* my_id ==0 */
   {
      num_types = send_proc_obj.length;
      used_procs = jxf_CTAlloc(JXF_Int, num_types);
      new_vec_starts = jxf_CTAlloc(JXF_Int, num_types + 1);

      new_vec_starts[0] = 0;
      for (i = 0; i < num_types; i++)
      {
         used_procs[i] = send_proc_obj.id[i];
         new_vec_starts[i + 1] = send_proc_obj.elements[i] + 1;
      }
      jxf_qsort0(used_procs, 0, num_types - 1);
      jxf_qsort0(new_vec_starts, 0, num_types);

      /*now we need to put into an array to send */
      count = 2 * num_types + 2;
      send_info = jxf_CTAlloc(JXF_Int, count);
      send_info[0] = num_types;
      for (i = 1; i <= num_types; i++)
      {
         send_info[i] = used_procs[i - 1];
      }
      for (i = num_types + 1; i < count; i++)
      {
         send_info[i] = new_vec_starts[i - num_types - 1];
      }
      requests = jxf_CTAlloc(MPI_Request, num_types);
      status = jxf_CTAlloc(MPI_Status, num_types);

      /* don't send to myself - these are sorted so my id would be first */
      start = 0;
      if (used_procs[0] == 0)
      {
         start = 1;
      }

      for (i = start; i < num_types; i++)
      {
         jxf_MPI_Isend(send_info, count, JXF_MPI_INT, used_procs[i], tag1, comm, &requests[i - start]);
      }
      jxf_MPI_Waitall(num_types - start, requests, status);

      jxf_TFree(status);
      jxf_TFree(requests);
   }

   /* clean up */
   jxf_TFree(send_proc_obj.vec_starts);
   jxf_TFree(send_proc_obj.id);
   jxf_TFree(send_proc_obj.elements);
   jxf_TFree(send_info);
   if (response_recv_buf)
      jxf_TFree(response_recv_buf);
   if (response_recv_buf_starts)
      jxf_TFree(response_recv_buf_starts);

   /* now proc 0 can exit if it has no rows */
   if (!local_size)
   {
      jxf_TFree(used_procs);
      jxf_TFree(new_vec_starts);
      return NULL;
   }

   /* everyone left has rows and knows: new_vec_starts, num_types, and used_procs */

   /* this vector should be rather small */

   local_data = jxf_VectorData(local_vector);
   vector = jxf_SeqVectorCreate(global_size);
   jxf_VectorNumVectors(vector) = num_vectors;
   jxf_SeqVectorInitialize(vector);
   vector_data = jxf_VectorData(vector);

   num_requests = 2 * num_types;

   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   /* initialize data exchange among used_procs and generate vector
      - here we send to ourself also */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = new_vec_starts[i + 1] - new_vec_starts[i];
      jxf_MPI_Irecv(&vector_data[new_vec_starts[i]], num_vectors * vec_len, JXF_MPI_REAL,
                   proc_id, tag2, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      jxf_MPI_Isend(local_data, num_vectors * local_size, JXF_MPI_REAL, used_procs[i],
                   tag2, comm, &requests[j++]);
   }

   jxf_MPI_Waitall(num_requests, requests, status);

   if (num_requests)
   {
      jxf_TFree(requests);
      jxf_TFree(status);
      jxf_TFree(used_procs);
   }

   jxf_TFree(new_vec_starts);

#else

   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   /* if my_id contains no data, return NULL  */

   if (!local_size)
   {
      return NULL;
   }

   local_data = jxf_VectorData(local_vector);
   vector = jxf_SeqVectorCreate(global_size);
   jxf_VectorNumVectors(vector) = num_vectors;
   jxf_SeqVectorInitialize(vector);
   vector_data = jxf_VectorData(vector);

   /* determine procs which hold data of par_v and store ids in used_procs */

   num_types = -1;
   for (i = 0; i < num_procs; i++)
   {
      if (vec_starts[i + 1] - vec_starts[i])
      {
         num_types++;
      }
   }
   num_requests = 2 * num_types;

   used_procs = jxf_CTAlloc(JXF_Int, num_types);
   j = 0;
   for (i = 0; i < num_procs; i++)
   {
      if (vec_starts[i + 1] - vec_starts[i] && i - my_id)
      {
         used_procs[j++] = i;
      }
   }

   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   /* initialize data exchange among used_procs and generate vector */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = vec_starts[proc_id + 1] - vec_starts[proc_id];
      jxf_MPI_Irecv(&vector_data[vec_starts[proc_id]], num_vectors * vec_len, JXF_MPI_REAL,
                   proc_id, 0, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      jxf_MPI_Isend(local_data, num_vectors * local_size, JXF_MPI_REAL, used_procs[i],
                   0, comm, &requests[j++]);
   }

   for (i = 0; i < num_vectors * local_size; i++)
   {
      vector_data[vec_starts[my_id] + i] = local_data[i];
   }

   jxf_MPI_Waitall(num_requests, requests, status);

   if (num_requests)
   {
      jxf_TFree(used_procs);
      jxf_TFree(requests);
      jxf_TFree(status);
   }

#endif

   return vector;
}

/*!
 * \fn JXF_Int jxf_ParVectorToVector_Alloctaed
 * \brief Generates a Vector on every proc which has a piece of the data
 *        from a ParVector on several procs in comm, vec_starts needs
 *        to contain the partitioning across all procs in comm.
 * \note The memory for the data of 'vector' has been already allocated.
 * \author peghoty
 * \date 2011/08/31
 */
JXF_Int
jxf_ParVectorToVector_Alloctaed(jxf_ParVector *par_v,
                               jxf_Vector *vector)
{
   MPI_Comm comm = jxf_ParVectorComm(par_v);
   JXF_Int *vec_starts = jxf_ParVectorPartitioning(par_v);
   jxf_Vector *local_vector = jxf_ParVectorLocalVector(par_v);
   JXF_Int num_procs, my_id;
   JXF_Int num_vectors = jxf_ParVectorNumVectors(par_v);
   JXF_Real *vector_data;
   JXF_Real *local_data;
   JXF_Int local_size;
   MPI_Request *requests;
   MPI_Status *status;
   JXF_Int *used_procs;
   JXF_Int num_types, num_requests;
   JXF_Int i, j, vec_len, proc_id;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   /* if my_id contains no data, return NULL */

   if (!local_size)
   {
      return -1;
   }

   local_data = jxf_VectorData(local_vector);
   vector_data = jxf_VectorData(vector);

   /*------------------------------------------------------------------------
    * determine procs which hold data of par_v and store ids in used_procs
    *-----------------------------------------------------------------------*/

   num_types = -1;
   for (i = 0; i < num_procs; i++)
   {
      if (vec_starts[i + 1] - vec_starts[i])

      {
         num_types++;
      }
   }
   num_requests = 2 * num_types;

   used_procs = jxf_CTAlloc(JXF_Int, num_types);
   j = 0;
   for (i = 0; i < num_procs; i++)
   {
      if (vec_starts[i + 1] - vec_starts[i] && i - my_id)
      {
         used_procs[j++] = i;
      }
   }

   requests = jxf_CTAlloc(MPI_Request, num_requests);
   status = jxf_CTAlloc(MPI_Status, num_requests);

   /*------------------------------------------------------------------------
    * initialize data exchange among used_procs and generate vector
    *-----------------------------------------------------------------------*/

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = vec_starts[proc_id + 1] - vec_starts[proc_id];
      jxf_MPI_Irecv(&vector_data[vec_starts[proc_id]], num_vectors * vec_len,
                   JXF_MPI_REAL, proc_id, 0, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      jxf_MPI_Isend(local_data, num_vectors * local_size, JXF_MPI_REAL,
                   used_procs[i], 0, comm, &requests[j++]);
   }

   for (i = 0; i < num_vectors * local_size; i++)
   {
      vector_data[vec_starts[my_id] + i] = local_data[i];
   }

   jxf_MPI_Waitall(num_requests, requests, status);

   if (num_requests)
   {
      jxf_TFree(used_procs);
      jxf_TFree(requests);
      jxf_TFree(status);
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_GenerateLocalPartitioning
 * \brief It only returns  the portion of the partition
 *        belonging to the individual process
 *        - to do this it requires the processor id as well. AHB 6/05
 * \date 2011/09/03
 */
JXF_Int
jxf_GenerateLocalPartitioning(JXF_Int length, JXF_Int num_procs, JXF_Int myid, JXF_Int **part_ptr)
{
   JXF_Int ierr = 0;
   JXF_Int *part;
   JXF_Int size, rest;

   part = jxf_CTAlloc(JXF_Int, 2);
   size = length / num_procs;
   rest = length - size * num_procs;

   /* first row I own */
   part[0] = size * myid;
   part[0] += jxf_min(myid, rest);

   /* last row I own */
   part[1] = size * (myid + 1);
   part[1] += jxf_min(myid + 1, rest);
   part[1] = part[1] - 1;

   /* add 1 to last row since this is for "starts" vector */
   part[1] = part[1] + 1;

   *part_ptr = part;

   return ierr;
}

/*!
 * \fn JXF_Int jxf_ParVectorPrint
 * \brief Print a parallel vector into given file(s).
 * \author peghoty
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorPrint(jxf_ParVector *vector, const char *file_name)
{
   char new_file_name[80];
   jxf_Vector *local_vector;
   MPI_Comm comm;
   JXF_Int my_id, num_procs, i;
   JXF_Int *partitioning;
   JXF_Int global_size;
   FILE *fp;

   if (!vector)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   local_vector = jxf_ParVectorLocalVector(vector);
   comm = jxf_ParVectorComm(vector);
   partitioning = jxf_ParVectorPartitioning(vector);
   global_size = jxf_ParVectorGlobalSize(vector);

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   jxf_sprintf(new_file_name, "%s.%d", file_name, my_id);
   jxf_SeqVectorPrint(local_vector, new_file_name);
   jxf_sprintf(new_file_name, "%s.INFO.%d", file_name, my_id);

   fp = fopen(new_file_name, "w");
   jxf_fprintf(fp, "%d\n", global_size);

#ifdef JXF_NO_GLOBAL_PARTITION
   for (i = 0; i < 2; i++)
   {
      jxf_fprintf(fp, "%d\n", partitioning[i]);
   }
#else
   for (i = 0; i < num_procs; i++)
   {
      jxf_fprintf(fp, "%d\n", partitioning[i]);
   }
#endif

   fclose(fp);
   return jxf_error_flag;
}

/*!
 * \fn jxf_ParVector *jxf_ParVectorRead
 * \brief Generate a parallel vector by reading data from given file(s)
 * \date 2011/09/03
 */
jxf_ParVector *
jxf_ParVectorRead(MPI_Comm comm, const char *file_name)
{
   char new_file_name[80];
   jxf_ParVector *par_vector;
   JXF_Int my_id, num_procs;
   JXF_Int *partitioning;
   JXF_Int global_size, i;
   FILE *fp;

   jxf_MPI_Comm_rank(comm, &my_id);
   jxf_MPI_Comm_size(comm, &num_procs);

   partitioning = jxf_CTAlloc(JXF_Int, num_procs + 1);

   jxf_sprintf(new_file_name, "%s.INFO.%d", file_name, my_id);
   fp = fopen(new_file_name, "r");
   jxf_fscanf(fp, "%d\n", &global_size);

#ifdef JXF_NO_GLOBAL_PARTITION
   for (i = 0; i < 2; i++)
   {
      jxf_fscanf(fp, "%d\n", &partitioning[i]);
   }
   fclose(fp);
#else
   for (i = 0; i < num_procs; i++)
   {
      jxf_fscanf(fp, "%d\n", &partitioning[i]);
   }
   fclose(fp);
   partitioning[num_procs] = global_size;
#endif

   par_vector = jxf_CTAlloc(jxf_ParVector, 1);

   jxf_ParVectorComm(par_vector) = comm;
   jxf_ParVectorGlobalSize(par_vector) = global_size;

#ifdef JXF_NO_GLOBAL_PARTITION
   jxf_ParVectorFirstIndex(par_vector) = partitioning[0];
   jxf_ParVectorLastIndex(par_vector) = partitioning[1] - 1;
#else
   jxf_ParVectorFirstIndex(par_vector) = partitioning[my_id];
   jxf_ParVectorLastIndex(par_vector) = partitioning[my_id + 1] - 1;
#endif

   jxf_ParVectorPartitioning(par_vector) = partitioning;

   jxf_ParVectorOwnsData(par_vector) = 1;
   jxf_ParVectorOwnsPartitioning(par_vector) = 1;

   jxf_sprintf(new_file_name, "%s.%d", file_name, my_id);
   jxf_ParVectorLocalVector(par_vector) = jxf_SeqVectorRead(new_file_name);

   /* multivector code not written yet >>> */
   jxf_assert(jxf_ParVectorNumVectors(par_vector) == 1);

   return par_vector;
}

/*!
 * \fn JXF_Int jxf_GeneratePartitioning
 * \brief Generates load balanced partitioning of a 1-d array.
 * \note For multivectors, length should be the (global) length of a single vector.
 *       Thus each of the vectors of the multivector will get the same data distribution.
 * \date 2011/09/03
 */
JXF_Int
jxf_GeneratePartitioning(JXF_Int length, JXF_Int num_procs, JXF_Int **part_ptr)
{
   JXF_Int ierr = 0;
   JXF_Int *part;
   JXF_Int size, rest;
   JXF_Int i;

   part = jxf_CTAlloc(JXF_Int, num_procs + 1);
   size = length / num_procs;
   rest = length - size * num_procs;
   part[0] = 0;
   for (i = 0; i < num_procs; i++)
   {
      part[i + 1] = part[i] + size;
      if (i < rest)
         part[i + 1]++;
   }

   *part_ptr = part;
   return ierr;
}

/*!
 * \fn JXF_Int jxf_FillResponseParToVectorAll
 * \brief Fill response function for determining the
 *        send processors data exchange.
 * \date 2011/09/03
 */
JXF_Int
jxf_FillResponseParToVectorAll(void *p_recv_contact_buf,
                              JXF_Int contact_size,
                              JXF_Int contact_proc,
                              void *ro,
                              MPI_Comm comm,
                              void **p_send_response_buf,
                              JXF_Int *response_message_size)
{
   JXF_Int myid;
   JXF_Int i, index, count, elength;

   JXF_Int *recv_contact_buf = (JXF_Int *)p_recv_contact_buf;

   jxf_DataExchangeResponse *response_obj = ro;

   jxf_ProcListElements *send_proc_obj = response_obj->data2;

   jxf_MPI_Comm_rank(comm, &myid);

   /* check to see if we need to allocate more space in send_proc_obj for ids */
   if (send_proc_obj->length == send_proc_obj->storage_length)
   {
      send_proc_obj->storage_length += 10; /* add space for 10 more processors */
      send_proc_obj->id = jxf_TReAlloc(send_proc_obj->id, JXF_Int, send_proc_obj->storage_length);
      send_proc_obj->vec_starts = jxf_TReAlloc(send_proc_obj->vec_starts, JXF_Int, send_proc_obj->storage_length + 1);
   }

   /* initialize */
   count = send_proc_obj->length;
   index = send_proc_obj->vec_starts[count]; /* this is the number of elements */

   /* send proc */
   send_proc_obj->id[count] = contact_proc;

   /* do we need more storage for the elements? */
   if (send_proc_obj->element_storage_length < index + contact_size)
   {
      elength = jxf_max(contact_size, 10);
      elength += index;
      send_proc_obj->elements = jxf_TReAlloc(send_proc_obj->elements, JXF_Int, elength);
      send_proc_obj->element_storage_length = elength;
   }
   /* populate send_proc_obj */
   for (i = 0; i < contact_size; i++)
   {
      send_proc_obj->elements[index++] = recv_contact_buf[i];
   }
   send_proc_obj->vec_starts[count + 1] = index;
   send_proc_obj->length++;

   /* output - no message to return (confirmation) */
   *response_message_size = 0;

   return jxf_error_flag;
}

/*!
 * \fn jxf_ParVector *jxf_VectorToParVector
 * \brief Generates a ParVector from a Vector on proc 0 and
 *        distributes the pieces to the other procs in comm.
 * \note This is not being optimized to use JXF_NO_GLOBAL_PARTITION.
 * \date 2011/09/03
 */
jxf_ParVector *
jxf_VectorToParVector(MPI_Comm comm, jxf_Vector *v, JXF_Int *vec_starts)
{
   JXF_Int global_size;
   JXF_Int local_size;
   JXF_Int num_vectors;
   JXF_Int num_procs, my_id;
   JXF_Int global_vecstride, vecstride, idxstride;
   jxf_ParVector *par_vector;
   jxf_Vector *local_vector;
   JXF_Real *v_data = NULL;
   JXF_Real *local_data;
   MPI_Request *requests;
   MPI_Status *status, status0;
   JXF_Int i, j, k, p;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == 0)
   {
      global_size = jxf_VectorSize(v);
      v_data = jxf_VectorData(v);
      num_vectors = jxf_VectorNumVectors(v); /* for multivectors */
      global_vecstride = jxf_VectorVectorStride(v);
   }

   jxf_MPI_Bcast(&global_size, 1, JXF_MPI_INT, 0, comm);
   jxf_MPI_Bcast(&num_vectors, 1, JXF_MPI_INT, 0, comm);
   jxf_MPI_Bcast(&global_vecstride, 1, JXF_MPI_INT, 0, comm);

   if (num_vectors == 1)
   {
      par_vector = jxf_ParVectorCreate(comm, global_size, vec_starts);
   }
   else
   {
      par_vector = jxf_ParMultiVectorCreate(comm, global_size, vec_starts, num_vectors);
   }

   vec_starts = jxf_ParVectorPartitioning(par_vector);

   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   jxf_ParVectorInitialize(par_vector);
   local_vector = jxf_ParVectorLocalVector(par_vector);
   local_data = jxf_VectorData(local_vector);
   vecstride = jxf_VectorVectorStride(local_vector);
   idxstride = jxf_VectorIndexStride(local_vector);
   jxf_assert(idxstride == 1); /* <<< so far only the only implemented multivector StorageMethod is 0 <<< */

   if (my_id == 0)
   {
      requests = jxf_CTAlloc(MPI_Request, num_vectors * (num_procs - 1));
      status = jxf_CTAlloc(MPI_Status, num_vectors * (num_procs - 1));
      k = 0;
      for (p = 1; p < num_procs; p++)
      {
         for (j = 0; j < num_vectors; ++j)
         {
            jxf_MPI_Isend(&v_data[vec_starts[p]] + j * global_vecstride,
                         (vec_starts[p + 1] - vec_starts[p]),
                         JXF_MPI_REAL, p, 0, comm, &requests[k++]);
         }
      }

      if (num_vectors == 1)
      {
         for (i = 0; i < local_size; i++)
         {
            local_data[i] = v_data[i];
         }
      }
      else
      {
         for (j = 0; j < num_vectors; ++j)
         {
            for (i = 0; i < local_size; i++)
            {
               local_data[i + j * vecstride] = v_data[i + j * global_vecstride];
            }
         }
      }
      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      for (j = 0; j < num_vectors; ++j)
      {
         jxf_MPI_Recv(local_data + j * vecstride, local_size, JXF_MPI_REAL, 0, 0, comm, &status0);
      }
   }

   return par_vector;
}

/*!
 * \fn jxf_ParVector *jxf_VectorToParVector_FromGivenPro
 * \brief Generates a ParVector from a Vector on the given proc (srcid#) and
 *        distributes the pieces to the other procs in comm.
 *        利用指定进程(srcid 号进程)上的一个串行向量生成一个并行向量. 这里只考虑 num_vectors = 1 的情形!
 * \note This is not being optimized to use JXF_NO_GLOBAL_PARTITION.
 * \author peghoty
 * \date 2012/02/29
 */
jxf_ParVector *
jxf_VectorToParVector_FromGivenPro(MPI_Comm comm, JXF_Int srcid, jxf_Vector *v, JXF_Int *vec_starts)
{
   JXF_Int global_size;
   JXF_Int local_size;
   JXF_Int num_procs, my_id;
   jxf_ParVector *par_vector;
   jxf_Vector *local_vector;
   JXF_Real *v_data = NULL;
   JXF_Real *local_data;
   MPI_Request *requests;
   MPI_Status *status, status0;
   JXF_Int i, k, p;
   JXF_Int start = 0;
   JXF_Int length = 0;
   JXF_Int partitioning_owner;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == srcid)
   {
      global_size = jxf_VectorSize(v);
      v_data = jxf_VectorData(v);
   }

   jxf_MPI_Bcast(&global_size, 1, JXF_MPI_INT, srcid, comm); // 从 srcid 号进程广播出去. peghoty, 2012/03/05

   /* set partitioning_owner to be zero is quite important here, 2012/03/05 */
   if (vec_starts)
   {
      partitioning_owner = 0;
   }
   else
   {
      partitioning_owner = 1;
   }

   par_vector = jxf_ParVectorCreate(comm, global_size, vec_starts);
   jxf_ParVectorSetPartitioningOwner(par_vector, partitioning_owner);
   vec_starts = jxf_ParVectorPartitioning(par_vector);
   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   jxf_ParVectorInitialize(par_vector);
   local_vector = jxf_ParVectorLocalVector(par_vector);
   local_data = jxf_VectorData(local_vector);

   if (my_id == srcid)
   {
      requests = jxf_CTAlloc(MPI_Request, num_procs - 1);
      status = jxf_CTAlloc(MPI_Status, num_procs - 1);
      k = 0;

      for (p = 0; p < num_procs; p++)
      {
         if (p != srcid)
         {
            start = vec_starts[p];
            length = vec_starts[p + 1] - start;
            jxf_MPI_Isend(&v_data[start], length, JXF_MPI_REAL, p, 0, comm, &requests[k++]);
         }
      }

      for (i = 0; i < local_size; i++)
      {
         local_data[i] = v_data[i];
      }

      start = vec_starts[my_id];
      for (i = 0; i < local_size; i++)
      {
         local_data[i] = v_data[start + i]; // be careful here, DON't forget the start. peghoty,  2012/03/07
      }

      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      jxf_MPI_Recv(local_data, local_size, JXF_MPI_REAL, srcid, 0, comm, &status0);
   }

   return par_vector;
}

/*!
 * \fn jxf_ParVector *jxf_VectorToParVector_sp
 * \brief Transfer a serial vector to a parallel one (only for single-processor case).
 * \note 该模块的特点是: 转换后的并行向量和串行向量共享数据部分。因此释放内存时要特别小心，
 *       不要轻易将串行向量 Destroy，一般采用 Free 即可.
 * \date 2011/09/27
 */
jxf_ParVector *
jxf_VectorToParVector_sp(MPI_Comm comm, jxf_Vector *vector)
{
   jxf_ParVector *par_vector = NULL;
   JXF_Int global_size = jxf_VectorSize(vector);

#if 0  
   jxf_Vector *local_vector = NULL; 
   par_vector = jxf_ParVectorCreate(comm, global_size, NULL);
   local_vector = jxf_ParVectorLocalVector(par_vector);
   jxf_VectorData(local_vector) = jxf_VectorData(vector);
#else // modified by peghoty, 2012/02/25
   JXF_Real temp_adrress = 0.0;
   par_vector = jxf_ParVectorCreate(comm, global_size, NULL);
   par_vector->local_vector->data = &temp_adrress; // 保证调用 jxf_ParVectorInitialize 时不另外开设空间
   jxf_ParVectorInitialize(par_vector);             // 保证成员 multivec_storage_method, vecstride, idxstride 被赋上适当的默认值，
                                                   // 因为在调用 jxf_ParCSRMatrixMatvec 时要检查，如果赋值不对，会有警告信息.
   jxf_ParVectorSetDataOwner(par_vector, 1);        // 保证 Destroy 并行向量时，能释放内存.
   par_vector->local_vector->data = vector->data;  // 建立数据共享关系
#endif

   return (par_vector);
}

/*!
 * \fn JXF_Int jxf_VectorToParVector_Allocated
 * \brief Generates a ParVector from a Vector on proc 0 and distributes the pieces
 *        to the other procs in comm.
 * \note (1) The memory for the data of 'par_vector' has been already allocated.
 *       (2) vec_starts should not be NULL!
 *       (3) This is not being optimized to use JXF_NO_GLOBAL_PARTITION.
 * \author peghoty
 * \date 2011/08/31
 */
JXF_Int
jxf_VectorToParVector_Allocated(MPI_Comm comm,
                               jxf_Vector *v,
                               JXF_Int *vec_starts,
                               jxf_ParVector *par_vector)
{
   JXF_Int global_size;
   JXF_Int local_size;
   JXF_Int num_vectors;
   JXF_Int num_procs, my_id;
   JXF_Int global_vecstride, vecstride, idxstride;
   jxf_Vector *local_vector;
   JXF_Real *v_data = NULL;
   JXF_Real *local_data;
   MPI_Request *requests;
   MPI_Status *status, status0;
   JXF_Int i, j, k, p;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == 0)
   {
      global_size = jxf_VectorSize(v);
      v_data = jxf_VectorData(v);
      num_vectors = jxf_VectorNumVectors(v); /* for multivectors */
      global_vecstride = jxf_VectorVectorStride(v);
   }

   jxf_MPI_Bcast(&global_size, 1, JXF_MPI_INT, 0, comm);
   jxf_MPI_Bcast(&num_vectors, 1, JXF_MPI_INT, 0, comm);
   jxf_MPI_Bcast(&global_vecstride, 1, JXF_MPI_INT, 0, comm);

   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   local_vector = jxf_ParVectorLocalVector(par_vector);
   local_data = jxf_VectorData(local_vector);

   vecstride = jxf_VectorVectorStride(local_vector);
   idxstride = jxf_VectorIndexStride(local_vector);
   jxf_assert(idxstride == 1); /* <<< so far only the only implemented multivector StorageMethod is 0 <<< */

   if (my_id == 0)
   {
      requests = jxf_CTAlloc(MPI_Request, num_vectors * (num_procs - 1));
      status = jxf_CTAlloc(MPI_Status, num_vectors * (num_procs - 1));
      k = 0;
      for (p = 1; p < num_procs; p++)
      {
         for (j = 0; j < num_vectors; ++j)
         {
            jxf_MPI_Isend(&v_data[vec_starts[p]] + j * global_vecstride,
                         (vec_starts[p + 1] - vec_starts[p]),
                         JXF_MPI_REAL, p, 0, comm, &requests[k++]);
         }
      }
      if (num_vectors == 1)
      {
         for (i = 0; i < local_size; i++)
         {
            local_data[i] = v_data[i];
         }
      }
      else
      {
         for (j = 0; j < num_vectors; ++j)
         {
            for (i = 0; i < local_size; i++)
            {
               local_data[i + j * vecstride] = v_data[i + j * global_vecstride];
            }
         }
      }
      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      for (j = 0; j < num_vectors; ++j)
      {
         jxf_MPI_Recv(local_data + j * vecstride, local_size, JXF_MPI_REAL, 0, 0, comm, &status0);
      }
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_VectorToParVector_Allocated2
 * \brief 利用 0 号进程上的一个串行向量生成一个并行向量，其中并行向量具有并行分划且已开设数据部分.
 * \note 在 jxf_VectorToParVector_Allocated 的基础上进行了简化，取 num_vectors = 1.
 * \author peghoty
 * \date 2012/02/27
 */
JXF_Int
jxf_VectorToParVector_Allocated2(MPI_Comm comm,
                                jxf_Vector *v,
                                JXF_Int *vec_starts,
                                jxf_ParVector *par_vector)
{
   JXF_Int local_size;
   JXF_Int num_procs, my_id;
   jxf_Vector *local_vector;
   JXF_Real *v_data = NULL;
   JXF_Real *local_data;
   MPI_Request *requests;
   MPI_Status *status, status0;
   JXF_Int i, k, p;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == 0)
   {
      v_data = jxf_VectorData(v);
   }

   local_size = vec_starts[my_id + 1] - vec_starts[my_id];

   local_vector = jxf_ParVectorLocalVector(par_vector);
   local_data = jxf_VectorData(local_vector);

   if (my_id == 0)
   {
      requests = jxf_CTAlloc(MPI_Request, num_procs - 1);
      status = jxf_CTAlloc(MPI_Status, num_procs - 1);
      k = 0;
      for (p = 1; p < num_procs; p++)
      {
         jxf_MPI_Isend(&v_data[vec_starts[p]],
                      (vec_starts[p + 1] - vec_starts[p]),
                      JXF_MPI_REAL, p, 0, comm, &requests[k++]);
      }

      for (i = 0; i < local_size; i++)
      {
         local_data[i] = v_data[i];
      }

      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      jxf_MPI_Recv(local_data, local_size, JXF_MPI_REAL, 0, 0, comm, &status0);
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_VectorToParVector_Allocated_FromGivenPro
 * \brief 利用指定进程(srcid 号进程)上的一个串行向量生成一个并行向量，其中并行向量具有并行分划且已开设数据部分.
 * \note 该模块基于 jxf_VectorToParVector_Allocated2 来改写的，只适合于 num_vectors = 1 的情形.
 * \author peghoty
 * \date 2012/02/27
 */
JXF_Int
jxf_VectorToParVector_Allocated_FromGivenPro(MPI_Comm comm,
                                            JXF_Int srcid,
                                            jxf_Vector *v,
                                            JXF_Int *vec_starts,
                                            jxf_ParVector *par_vector)
{
   JXF_Int local_size;
   JXF_Int num_procs, my_id;
   jxf_Vector *local_vector = NULL;
   JXF_Real *v_data = NULL;
   JXF_Real *local_data = NULL;
   MPI_Request *requests = NULL;
   MPI_Status *status = NULL;
   MPI_Status status0;
   JXF_Int i, k, p;
   JXF_Int start = 0;
   JXF_Int length = 0;
   JXF_Int tag = 111;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   if (my_id == srcid)
   {
      v_data = jxf_VectorData(v);
   }

   local_size = vec_starts[my_id + 1] - vec_starts[my_id]; // 本地进程上向量的长度

   local_vector = jxf_ParVectorLocalVector(par_vector);
   local_data = jxf_VectorData(local_vector);

   if (my_id == srcid)
   {
      requests = jxf_CTAlloc(MPI_Request, num_procs - 1);
      status = jxf_CTAlloc(MPI_Status, num_procs - 1);
      k = 0;

      for (p = 0; p < num_procs; p++)
      {
         if (p != srcid)
         {
            start = vec_starts[p];
            length = vec_starts[p + 1] - start;
            jxf_MPI_Isend(&v_data[start], length, JXF_MPI_REAL, p, tag, comm, &requests[k++]);
         }
      }

      start = vec_starts[my_id];
      for (i = 0; i < local_size; i++)
      {
         local_data[i] = v_data[start + i]; // be careful here, DON't forget the start. peghoty,  2012/03/07
      }

      jxf_MPI_Waitall(num_procs - 1, requests, status);
      jxf_TFree(requests);
      jxf_TFree(status);
   }
   else
   {
      jxf_MPI_Recv(local_data, local_size, JXF_MPI_REAL, srcid, tag, comm, &status0);
   }

   return 0;
}

/*!
 * \fn JXF_Int jxf_ParVectorSetDataOwner
 * \brief Set owns_data of a parallel vector.
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorSetDataOwner(jxf_ParVector *vector, JXF_Int owns_data)
{
   if (!vector)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParVectorOwnsData(vector) = owns_data;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParVectorSetPartitioningOwner
 * \brief Set owns_partitioning of a parallel vector.
 * \date 2011/09/03
 */
JXF_Int
jxf_ParVectorSetPartitioningOwner(jxf_ParVector *vector, JXF_Int owns_partitioning)
{
   if (!vector)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }

   jxf_ParVectorOwnsPartitioning(vector) = owns_partitioning;

   return jxf_error_flag;
}

/*!
 * \fn jxf_ParVecCommPkg *jxf_ParVecCommPkgCreate
 * \brief Create a communication package for the switch between two
 *        parallel vectors with different partitionings.
 * \param comm communication
 * \param *px pointer to the partioning of the source vector.
 * \param *py pointer to the partioning of the target vector.
 * \author peghoty
 * \date 2011/09/07
 */
jxf_ParVecCommPkg *
jxf_ParVecCommPkgCreate(MPI_Comm comm, JXF_Int *px, JXF_Int *py)
{
   jxf_ParVecCommPkg *comm_pkg = NULL;

   JXF_Int num_sends;
   JXF_Int *send_procs = NULL;
   JXF_Int *send_starts = NULL;

   JXF_Int num_recvs;
   JXF_Int *recv_procs = NULL;
   JXF_Int *recv_starts = NULL;

   JXF_Int i, j, k, m;
   JXF_Int last = 0;
   JXF_Int start, end;
   JXF_Int local_size;
   JXF_Int myid, np;

   jxf_MPI_Comm_rank(comm, &myid);
   jxf_MPI_Comm_size(comm, &np);

   //-------------------------------------------
   //  统计发送进程个数 num_sends
   //-------------------------------------------

   start = px[myid];
   end = px[myid + 1];
   local_size = end - start;

   k = 0;
   for (i = 0; i < np; i++)
   {
      j = py[i];
      if (j > start && j < end)
      {
         k++;
      }
      else if (j >= end)
      {
         break;
      }
   }
   num_sends = k + 1; //

   //-------------------------------------------
   //  提取发送进程， 生成数组 send_procs
   //  同时生成发送数据管理数组 send_starts
   //-------------------------------------------

   send_procs = jxf_CTAlloc(JXF_Int, num_sends);
   send_starts = jxf_CTAlloc(JXF_Int, num_sends + 1);

   k = 0;
   m = start;
   send_starts[0] = 0;
   for (i = 0; i < np + 1; i++)
   {
      j = py[i];
      if (j > start && j < end)
      {
         send_procs[k] = i - 1;
         k++;

         send_starts[k] = send_starts[k - 1] + (j - m);
         m = j;
      }
      else if (j >= end)
      {
         last = i - 1;
         break;
      }
   }
   if (num_sends > 1)
   {
      send_procs[num_sends - 1] = send_procs[num_sends - 2] + 1; // 利用发送进程编号的单调性和连续性
   }
   else
   {
      send_procs[num_sends - 1] = last;
   }
   send_starts[num_sends] = local_size;

   //-------------------------------------------
   //  统计接收进程个数 num_recvs
   //-------------------------------------------

   start = py[myid];
   end = py[myid + 1];
   local_size = end - start;

   k = 0;
   for (i = 0; i < np; i++)
   {
      j = px[i];
      if (j > start && j < end)
      {
         k++;
      }
      else if (j >= end)
      {
         break;
      }
   }
   num_recvs = k + 1; //

   //-------------------------------------------
   //  提取接收进程， 生成数组 recv_procs
   //  同时生成接收数据管理数组 recv_starts
   //-------------------------------------------

   recv_procs = jxf_CTAlloc(JXF_Int, num_recvs);
   recv_starts = jxf_CTAlloc(JXF_Int, num_recvs + 1);

   k = 0;
   m = start;
   recv_starts[0] = 0;
   for (i = 0; i < np + 1; i++)
   {
      j = px[i];
      if (j > start && j < end)
      {
         recv_procs[k] = i - 1;
         k++;

         recv_starts[k] = recv_starts[k - 1] + (j - m);
         m = j;
      }
      else if (j >= end)
      {
         last = i - 1;
         break;
      }
   }
   if (num_recvs > 1)
   {
      recv_procs[num_recvs - 1] = recv_procs[num_recvs - 2] + 1; // 利用发送进程编号的单调性和连续性
   }
   else
   {
      recv_procs[num_recvs - 1] = last;
   }
   recv_starts[num_recvs] = local_size;

   comm_pkg = jxf_CTAlloc(jxf_ParVecCommPkg, 1);

   jxf_ParVecCommPkgComm(comm_pkg) = comm;
   jxf_ParVecCommPkgNumSends(comm_pkg) = num_sends;
   jxf_ParVecCommPkgSendProcs(comm_pkg) = send_procs;
   jxf_ParVecCommPkgSendStarts(comm_pkg) = send_starts;
   jxf_ParVecCommPkgNumRecvs(comm_pkg) = num_recvs;
   jxf_ParVecCommPkgRecvProcs(comm_pkg) = recv_procs;
   jxf_ParVecCommPkgRecvStarts(comm_pkg) = recv_starts;

   return (comm_pkg);
}

/*!
 * \fn jxf_ParVecCommHandle *jxf_ParVecCommHandleCreate
 * \brief Create a communication handle for the switch between two
 *        parallel vectors with different partitionings.
 *        Data exchanging takes place in this function.
 * \author peghoty
 * \date 2011/09/07
 */
jxf_ParVecCommHandle *
jxf_ParVecCommHandleCreate(JXF_Int job,
                          jxf_ParVecCommPkg *comm_pkg,
                          void *ssend_data,
                          void *rrecv_data)
{
   jxf_ParVecCommHandle *comm_handle = NULL;

   /* members of the object comm_pkg */
   MPI_Comm comm = jxf_ParVecCommPkgComm(comm_pkg);
   JXF_Int num_sends = jxf_ParVecCommPkgNumSends(comm_pkg);
   JXF_Int num_recvs = jxf_ParVecCommPkgNumRecvs(comm_pkg);
   JXF_Int *send_procs = jxf_ParVecCommPkgSendProcs(comm_pkg);
   JXF_Int *send_starts = jxf_ParVecCommPkgSendStarts(comm_pkg);
   JXF_Int *recv_procs = jxf_ParVecCommPkgRecvProcs(comm_pkg);
   JXF_Int *recv_starts = jxf_ParVecCommPkgRecvStarts(comm_pkg);

   JXF_Int num_requests;
   MPI_Request *requests = NULL;

   JXF_Int i, j;
   JXF_Int my_id, num_procs;
   JXF_Int ip, vec_start, vec_len;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);

   num_requests = num_sends + num_recvs;
   requests = jxf_CTAlloc(MPI_Request, num_requests);

   switch (job)
   {
   case 1: // JXF_Real
   {
      JXF_Real *send_data = (JXF_Real *)ssend_data;
      JXF_Real *recv_data = (JXF_Real *)rrecv_data;
      j = 0;
      for (i = 0; i < num_recvs; i++)
      {
         ip = recv_procs[i];
         vec_start = recv_starts[i];
         vec_len = recv_starts[i + 1] - vec_start;
         jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
      }
      for (i = 0; i < num_sends; i++)
      {
         ip = send_procs[i];
         vec_start = send_starts[i];
         vec_len = send_starts[i + 1] - vec_start;
         jxf_MPI_Isend(&send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
      }
   }
   break;

   case 2: // JXF_Real
   {
      JXF_Real *send_data = (JXF_Real *)ssend_data;
      JXF_Real *recv_data = (JXF_Real *)rrecv_data;
      j = 0;
      for (i = 0; i < num_sends; i++)
      {
         ip = send_procs[i];
         vec_start = send_starts[i];
         vec_len = send_starts[i + 1] - vec_start;
         jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
      }
      for (i = 0; i < num_recvs; i++)
      {
         ip = recv_procs[i];
         vec_start = recv_starts[i];
         vec_len = recv_starts[i + 1] - vec_start;
         jxf_MPI_Isend(&send_data[vec_start], vec_len, JXF_MPI_REAL, ip, 0, comm, &requests[j++]);
      }
   }
   break;

   case 11: // JXF_Int
   {
      JXF_Int *send_data = (JXF_Int *)ssend_data;
      JXF_Int *recv_data = (JXF_Int *)rrecv_data;
      j = 0;
      for (i = 0; i < num_recvs; i++)
      {
         ip = recv_procs[i];
         vec_start = recv_starts[i];
         vec_len = recv_starts[i + 1] - vec_start;
         jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
      }
      for (i = 0; i < num_sends; i++)
      {
         ip = send_procs[i];
         vec_start = send_starts[i];
         vec_len = send_starts[i + 1] - vec_start;
         jxf_MPI_Isend(&send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
      }
   }
   break;

   case 12: // JXF_Int
   {
      JXF_Int *send_data = (JXF_Int *)ssend_data;
      JXF_Int *recv_data = (JXF_Int *)rrecv_data;
      j = 0;
      for (i = 0; i < num_sends; i++)
      {
         ip = send_procs[i];
         vec_start = send_starts[i];
         vec_len = send_starts[i + 1] - vec_start;
         jxf_MPI_Irecv(&recv_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
      }
      for (i = 0; i < num_recvs; i++)
      {
         ip = recv_procs[i];
         vec_start = recv_starts[i];
         vec_len = recv_starts[i + 1] - vec_start;
         jxf_MPI_Isend(&send_data[vec_start], vec_len, JXF_MPI_INT, ip, 0, comm, &requests[j++]);
      }
   }
   break;
   }

   /*--------------------------------------------------------------------
    * set up comm_handle and return
    *--------------------------------------------------------------------*/

   comm_handle = jxf_CTAlloc(jxf_ParVecCommHandle, 1);

   /* 这三个成员目前还没有用到. peghoty, 2012/02/29 */
   // jxf_ParVecCommHandleCommPkg(comm_handle)     = comm_pkg;
   // jxf_ParVecCommHandleSendData(comm_handle)    = ssend_data;
   // jxf_ParVecCommHandleRecvData(comm_handle)    = rrecv_data;
   jxf_ParVecCommHandleNumRequests(comm_handle) = num_requests;
   jxf_ParVecCommHandleRequests(comm_handle) = requests;

   return (comm_handle);
}

/*!
 * \fn JXF_Int jxf_ParVecCommHandleDestroy
 * \brief Destroy a communication handle for the switch between two
 *        parallel vectors with different partitionings.
 * \author peghoty
 * \date 2011/09/07
 */
JXF_Int
jxf_ParVecCommHandleDestroy(jxf_ParVecCommHandle *comm_handle)
{
   MPI_Status *status = NULL;

   MPI_Request *requests = jxf_ParVecCommHandleRequests(comm_handle);
   JXF_Int num_requests = jxf_ParVecCommHandleNumRequests(comm_handle);

   if (comm_handle == NULL)
   {
      return jxf_error_flag;
   }

   if (num_requests)
   {
      status = jxf_CTAlloc(MPI_Status, num_requests);
      jxf_MPI_Waitall(num_requests, requests, status);
      jxf_TFree(status);
   }

   jxf_TFree(requests);
   jxf_TFree(comm_handle);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParVecCommPkgDestroy
 * \brief Destroy a communication package for the switch between two
 *        parallel vectors with different partitionings.
 * \author peghoty
 * \date 2011/09/07
 */
JXF_Int
jxf_ParVecCommPkgDestroy(jxf_ParVecCommPkg *comm_pkg)
{
   if (jxf_ParVecCommPkgNumSends(comm_pkg))
   {
      jxf_TFree(jxf_ParVecCommPkgSendProcs(comm_pkg));
   }
   jxf_TFree(jxf_ParVecCommPkgSendStarts(comm_pkg));

   if (jxf_ParVecCommPkgNumRecvs(comm_pkg))
   {
      jxf_TFree(jxf_ParVecCommPkgRecvProcs(comm_pkg));
   }
   jxf_TFree(jxf_ParVecCommPkgRecvStarts(comm_pkg));

   jxf_TFree(comm_pkg);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_ParVecSwitch
 * \brief Switch between two parallel vectors with different partitionings.
 * \param comm_pkg pointer to the communication package.
 * \param job flag to indicate what to do.
 *        job = 1: x -> y
 *        job = 2: y -> x
 * \param x pointer to the first vector
 * \param y pointer to the second vector
 * \author peghoty
 * \date 2011/09/07
 */
JXF_Int
jxf_ParVecSwitch(jxf_ParVecCommPkg *comm_pkg,
                JXF_Int job,
                jxf_ParVector *x,
                jxf_ParVector *y)
{
   JXF_Real *send_buf = NULL;
   JXF_Real *recv_buf = NULL;
   jxf_ParVecCommHandle *comm_handle = NULL;

   switch (job)
   {
   case 1: // x -> y
   {
      send_buf = jxf_VectorData(jxf_ParVectorLocalVector(x));
      recv_buf = jxf_VectorData(jxf_ParVectorLocalVector(y));

      comm_handle = jxf_ParVecCommHandleCreate(1, comm_pkg, send_buf, recv_buf);
      jxf_ParVecCommHandleDestroy(comm_handle);
   }
   break;

   case 2: // y -> x
   {
      send_buf = jxf_VectorData(jxf_ParVectorLocalVector(y));
      recv_buf = jxf_VectorData(jxf_ParVectorLocalVector(x));
      comm_handle = jxf_ParVecCommHandleCreate(2, comm_pkg, send_buf, recv_buf);
      jxf_ParVecCommHandleDestroy(comm_handle);
   }
   break;
   }

   return 0;
}
