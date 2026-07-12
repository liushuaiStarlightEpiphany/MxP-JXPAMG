//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  new_commpkg.c -- new commpkg subroutines.
 *  Date: 2011/09/07
 */ 

#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_NewCommPkgCreate_core
 * \brief This does the work for jxf_NewCommPkgCreate - we have to split it 
 *        off so that it can also be used for block matrices.
 * \date 2011/09/07
 */
JXF_Int 
jxf_NewCommPkgCreate_core( /* input args: */
                          MPI_Comm            comm, 
                          JXF_Int                *col_map_off_d, 
                          JXF_Int                 first_col_diag,
                          JXF_Int                 col_start, 
                          JXF_Int                 col_end, 
                          JXF_Int                 num_cols_off_d, 
                          JXF_Int                 global_num_cols,
                          /* pointers to output args: */
                          JXF_Int                *p_num_recvs, 
                          JXF_Int               **p_recv_procs, 
                          JXF_Int               **p_recv_vec_starts,
                          JXF_Int                *p_num_sends, 
                          JXF_Int               **p_send_procs, 
                          JXF_Int               **p_send_map_starts,
                          JXF_Int               **p_send_map_elements, 
                          jxf_IJAssumedPart   *apart)
{
   JXF_Int        num_procs, myid;
   JXF_Int        j, i;
   JXF_Int        range_start, range_end; 

   JXF_Int        size;
   JXF_Int        count;  

   JXF_Int        num_recvs, *recv_procs = NULL, *recv_vec_starts = NULL;
   JXF_Int        tmp_id, prev_id;

   JXF_Int        num_sends;

   JXF_Int        ex_num_contacts, *ex_contact_procs = NULL, *ex_contact_vec_starts = NULL;
   JXF_Int       *ex_contact_buf = NULL;
    
   JXF_Int        num_ranges, upper_bound;
   JXF_Int       *response_buf = NULL, *response_buf_starts = NULL;

   JXF_Int        max_response_size;
   
   jxf_DataExchangeResponse  response_obj1, response_obj2;
   jxf_ProcListElements      send_proc_obj; 

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &myid);

   /*-------------------------------------------------------------------
    *  Everyone knows where their assumed range is located
    * (because of the assumed partition object (apart).
    *  For the comm. package, each proc must know it's receive
    *  procs (who it will receive data from and how much data) 
    *  and its send procs 
    *  (who it will send data to) and the indices of the elements
    *  to be sent.  This is based on the non-zero
    *  entries in its rows. Each proc should know this from the user. 
    *------------------------------------------------------------------*/

   /*------------------------------------------------------------
    *  First, get the receive processors each par_csr
    *  matrix will have a certain number of columns
    *  (num_cols_off_d) given in col_map_offd[] for which 
    *  it needs data from another processor. 
    *------------------------------------------------------------*/

   /* calculate the assumed receive processors */

   /* need to populate num_recvs, *recv_procs, and *recv_vec_starts (correlates 
      to starts in col_map_off_d for recv_procs) for the comm. package */

   /* create contact information */

   ex_num_contacts = 0;

   /* estimate the storage needed */
   
   if ( num_cols_off_d > 0 && (apart->row_end - apart->row_start) > 0 )
   {
      size = col_map_off_d[num_cols_off_d-1] - col_map_off_d[0];
      size = (size / (apart->row_end - apart->row_start)) + 2;
   }
   else
   {
      size = 0;
   }
   

   /*we will contact each with a range of cols that we need */
   /* it is ok to contact yourself - because then there doesn't need to be separate code */

   ex_contact_procs = jxf_CTAlloc(JXF_Int, size);
   ex_contact_vec_starts = jxf_CTAlloc(JXF_Int, size+1);
   ex_contact_buf = jxf_CTAlloc(JXF_Int, size*2);

   range_end = -1;
   for (i = 0; i < num_cols_off_d; i ++) 
   { 
      if (col_map_off_d[i] > range_end)
      {
         jxf_GetAssumedPartitionProcFromRow(comm, col_map_off_d[i], 0, global_num_cols, &tmp_id);

         if (ex_num_contacts == size) /* need more space? */ 
         {
            size += 20;
            ex_contact_procs = jxf_TReAlloc(ex_contact_procs, JXF_Int, size);
            ex_contact_vec_starts = jxf_TReAlloc(ex_contact_vec_starts, JXF_Int, size + 1);
            ex_contact_buf = jxf_TReAlloc(ex_contact_buf, JXF_Int, size*2);
         }

         /* end of prev. range */
         if (ex_num_contacts > 0)  
         {
            ex_contact_buf[ex_num_contacts*2 - 1] = col_map_off_d[i-1];
         }
         
         /* start new range */
    	 ex_contact_procs[ex_num_contacts] = tmp_id;
         ex_contact_vec_starts[ex_num_contacts] = ex_num_contacts*2;
         ex_contact_buf[ex_num_contacts*2] = col_map_off_d[i];
         
         ex_num_contacts ++;

         jxf_GetAssumedPartitionRowRange(comm, tmp_id, 0, global_num_cols, &range_start, &range_end);
      }
   }

   /* finish the starts */
   ex_contact_vec_starts[ex_num_contacts] =  ex_num_contacts*2;
   
   /* finish the last range */
   if (ex_num_contacts > 0)  
   {
      ex_contact_buf[ex_num_contacts*2 - 1] = col_map_off_d[num_cols_off_d-1];
   }


   /* don't allocate space for responses */
    

   /* create response object */
   response_obj1.fill_response = jxf_RangeFillResponseIJDetermineRecvProcs;
   response_obj1.data1 = apart; /* this is necessary so we can fill responses */ 
   response_obj1.data2 = NULL;
   
   max_response_size = 6;  /* 6 means we can fit 3 ranges */
   
   jxf_DataExchangeList( ex_num_contacts, ex_contact_procs, 
                        ex_contact_buf, ex_contact_vec_starts, sizeof(JXF_Int), 
                        sizeof(JXF_Int), &response_obj1, max_response_size, 1, 
                        comm, (void**) &response_buf, &response_buf_starts );

   /* now create recv_procs[] and recv_vec_starts[] and num_recvs 
      from the complete data in response_buf - this array contains
      a proc_id followed by an upper bound for the range. */

   /* initialize */ 
   num_recvs = 0;
   size = ex_num_contacts + 20;  /* num of recv procs should be roughly similar size 
                                    to number of contacts - add a buffer of 20 */

   recv_procs = jxf_CTAlloc(JXF_Int, size);
   recv_vec_starts = jxf_CTAlloc(JXF_Int, size + 1);
   recv_vec_starts[0] = 0;
   
   /* how many ranges were returned? */
   num_ranges = response_buf_starts[ex_num_contacts];   
   num_ranges = num_ranges / 2;
   
   prev_id = -1;
   j = 0;
   count = 0;
   
   /* loop through ranges */
   for (i = 0; i < num_ranges; i ++)
   {
      upper_bound = response_buf[i*2+1];
      count = 0;
      
      /* loop through off_d entries - counting how many are in the range */
      while (j < num_cols_off_d && col_map_off_d[j] <= upper_bound)
      {
         j ++;
         count ++;       
      }
      if (count > 0)        
      {
         /* add the range if the proc id != myid */    
         tmp_id = response_buf[i*2];
         if (tmp_id != myid)
         {
            if (tmp_id != prev_id) /* increment the number of recvs */
            {
               /* check size of recv buffers */
               if (num_recvs == size) 
               {
                  size += 20;
                  recv_procs = jxf_TReAlloc(recv_procs, JXF_Int, size);
                  recv_vec_starts = jxf_TReAlloc(recv_vec_starts, JXF_Int, size + 1);
               }
            
               recv_vec_starts[num_recvs+1] = j; /* the new start is at this element */
               recv_procs[num_recvs] = tmp_id;   /* add the new processor */
               num_recvs++;

            }
            else
            {
               /* same processor - just change the vec starts */
               recv_vec_starts[num_recvs] = j;  /* the new start is at this element */
            }
         }
         prev_id = tmp_id;
      }
      
   }


   /*------------------------------------------------------------
    *  determine the send processors
    *  each processor contacts its recv procs to let them
    *  know they are a send processor
    *-----------------------------------------------------------*/

   /* the contact information is the recv_processor infomation
      - so nothing more to do to generate contact info */

   /* the response we expect is just a confirmation */
   jxf_TFree(response_buf);
   jxf_TFree(response_buf_starts);
   response_buf = NULL;
   response_buf_starts = NULL;

   /* build the response object */
   /* estimate for inital storage allocation that we send to as many procs as we recv from + pad by 5 */
   send_proc_obj.length = 0;
   send_proc_obj.storage_length = num_recvs + 5;
   send_proc_obj.id = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length);
   send_proc_obj.vec_starts = jxf_CTAlloc(JXF_Int, send_proc_obj.storage_length + 1); 
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = num_cols_off_d;
   send_proc_obj.elements = jxf_CTAlloc(JXF_Int, send_proc_obj.element_storage_length);

   response_obj2.fill_response = jxf_FillResponseIJDetermineSendProcs;
   response_obj2.data1 = NULL;
   response_obj2.data2 = &send_proc_obj; /* this is where we keep info from contacts */
  
   max_response_size = 0;
      
   jxf_DataExchangeList( num_recvs, recv_procs, 
                        col_map_off_d, recv_vec_starts, sizeof(JXF_Int),
                        sizeof(JXF_Int), &response_obj2, max_response_size, 2, 
                        comm, (void **) &response_buf, &response_buf_starts );

   num_sends = send_proc_obj.length; 

   /* send procs are in send_proc_object.id */
   /* send proc starts are in send_proc_obj.vec_starts */


  /*-----------------------------------------------------------------
   *  We need to sort the send procs and send elements (to produce
   *  the same result as with the standard comm package) 11/07/05
   *----------------------------------------------------------------*/
   
   {
      JXF_Int *orig_order;
      JXF_Int *orig_send_map_starts;
      JXF_Int *orig_send_elements;
      JXF_Int  ct, sz, pos;
      
      orig_order = jxf_CTAlloc(JXF_Int, num_sends);
      orig_send_map_starts = jxf_CTAlloc(JXF_Int, num_sends + 1);
      orig_send_elements = jxf_CTAlloc(JXF_Int, send_proc_obj.vec_starts[num_sends]);
      
      orig_send_map_starts[0] = 0;
      
      /* copy send map starts and elements */ 
      for (i = 0; i < num_sends; i ++)
      {
         orig_order[i] = i;
         orig_send_map_starts[i+1] = send_proc_obj.vec_starts[i+1];
      }
      for (i = 0; i < send_proc_obj.vec_starts[num_sends]; i ++)
      {
         orig_send_elements[i] = send_proc_obj.elements[i];
      }
      
      /* sort processor ids - keep track of original order */
      jxf_qsort2i( send_proc_obj.id, orig_order, 0, num_sends-1 );
      
      /* now rearrange vec starts and send elements to correspond to proc ids */ 
      ct = 0;
      for (i = 0; i < num_sends; i ++)
      {
         pos = orig_order[i];
         sz = orig_send_map_starts[pos + 1] - orig_send_map_starts[pos];
         send_proc_obj.vec_starts[i+1] = ct + sz;
         for (j = 0; j < sz; j ++)
         {
            send_proc_obj.elements[ct+j] = orig_send_elements[orig_send_map_starts[pos]+j];
         }
         ct += sz;
      }
      
      /* clean up */
      jxf_TFree(orig_order);
      jxf_TFree(orig_send_elements);
      jxf_TFree(orig_send_map_starts);
   }
      
  /*-----------------------------------------------------------
   *  Return output info for setting up the comm package
   *----------------------------------------------------------*/
  
   if (!num_recvs)
   {
      jxf_TFree(recv_procs);
      recv_procs = NULL;
   }
   if (!num_sends)
   {
      jxf_TFree(send_proc_obj.id);
      send_proc_obj.id = NULL;
   }

   *p_num_recvs = num_recvs;
   *p_recv_procs = recv_procs;
   *p_recv_vec_starts = recv_vec_starts;
   *p_num_sends = num_sends;
   *p_send_procs = send_proc_obj.id;
   *p_send_map_starts = send_proc_obj.vec_starts;


   /* send map elements have global index - need local instead */

   if (num_sends)
   {
      for (i = 0; i < send_proc_obj.vec_starts[num_sends]; i ++)
      {   
         send_proc_obj.elements[i] -= first_col_diag;
      }
   }
   else
   {
      jxf_TFree(send_proc_obj.elements);
      send_proc_obj.elements = NULL;
   }
   
   *p_send_map_elements = send_proc_obj.elements;


  /*-----------------------------------------------------------
   *  Clean up
   *----------------------------------------------------------*/

   if (ex_contact_procs)      jxf_TFree(ex_contact_procs);
   if (ex_contact_vec_starts) jxf_TFree(ex_contact_vec_starts);
   jxf_TFree(ex_contact_buf);
   
   if (response_buf)        jxf_TFree(response_buf);
   if (response_buf_starts) jxf_TFree(response_buf_starts);


  /*---------------------------------------------------------------------------------
   *  Don't free send_proc_obj.id, send_proc_obj.vec_starts, send_proc_obj.elements;
   *  recv_procs, recv_vec_starts. These are aliased to the comm package and
   *  will be destroyed there. 
   *--------------------------------------------------------------------------------*/

   return jxf_error_flag;
}


/*!
 * \fn JXF_Int jxf_RangeFillResponseIJDetermineRecvProcs
 * \brief Fill response function for determining the 
 *        recv. processors data exchange.
 * \date 2011/09/07
 */
JXF_Int
jxf_RangeFillResponseIJDetermineRecvProcs( void       *p_recv_contact_buf,
                                          JXF_Int         contact_size, 
                                          JXF_Int         contact_proc, 
                                          void       *ro, 
                                          MPI_Comm   comm, 
                                          void     **p_send_response_buf, 
                                          JXF_Int       *response_message_size )
{
   JXF_Int    myid, tmp_id, row_end;
   JXF_Int    j, row_val, index, size;
   
   JXF_Int   *send_response_buf = (JXF_Int *) *p_send_response_buf;
   JXF_Int   *recv_contact_buf  = (JXF_Int *)  p_recv_contact_buf;


   jxf_DataExchangeResponse  *response_obj = ro; 
   jxf_IJAssumedPart         *part = response_obj->data1;
   
   JXF_Int overhead = response_obj->send_response_overhead;

  /*------------------------------------------------------------------------
   * we are getting a range of off_d entries - need to see if we own them
   * or how many ranges to send back  - send back with format
   * [proc_id end_row  proc_id #end_row  proc_id #end_row etc...].
   *-----------------------------------------------------------------------*/ 

   jxf_MPI_Comm_rank(comm, &myid);


   /* populate send_response_buf */
      
   index = 0; /* count entries in send_response_buf */
   
   j = 0; /* marks which partition of the assumed partition we are in */
   row_val = recv_contact_buf[0]; /* beginning of range */
   row_end = part->row_end_list[part->sort_index[j]];
   tmp_id  = part->proc_list[part->sort_index[j]];

   /* check storage in send_buf for adding the ranges */
   size = 2*(part->length);
         
   if (response_obj->send_response_storage < size)
   {
      response_obj->send_response_storage = jxf_max(size, 20); 
      send_response_buf = jxf_TReAlloc( send_response_buf, JXF_Int, 
                                       response_obj->send_response_storage + overhead );
      *p_send_response_buf = send_response_buf; /* needed when using ReAlloc */
   }

   while (row_val > row_end) /* which partition to start in */
   {
      j ++;
      row_end = part->row_end_list[part->sort_index[j]];   
      tmp_id  = part->proc_list[part->sort_index[j]];
   }

   /* add this range */
   send_response_buf[index++] = tmp_id;
   send_response_buf[index++] = row_end; 

   j ++; /* increase j to look in next partition */
      
   /* any more? - now compare with end of range value */
   row_val = recv_contact_buf[1]; /* end of range */
   while ( j < part->length && row_val > row_end )
   {
      row_end = part->row_end_list[part->sort_index[j]];  
      tmp_id  = part->proc_list[part->sort_index[j]];

      send_response_buf[index++] = tmp_id;
      send_response_buf[index++] = row_end; 

      j++;
   }

   *response_message_size = index;
   *p_send_response_buf   = send_response_buf;

   return jxf_error_flag;
}


/*!
 * \fn JXF_Int jxf_FillResponseIJDetermineSendProcs
 * \brief Fill response function for determining the send 
 *        processors data exchange.
 * \date 2011/09/07
 */
JXF_Int
jxf_FillResponseIJDetermineSendProcs( void      *p_recv_contact_buf, 
                                     JXF_Int        contact_size, 
                                     JXF_Int        contact_proc, 
                                     void      *ro, 
                                     MPI_Comm   comm, 
                                     void     **p_send_response_buf, 
                                     JXF_Int       *response_message_size )
{
   JXF_Int    myid;
   JXF_Int    i, index, count, elength;

   JXF_Int    *recv_contact_buf = (JXF_Int *) p_recv_contact_buf;

   jxf_DataExchangeResponse  *response_obj = ro;  

   jxf_ProcListElements      *send_proc_obj = response_obj->data2;   


   jxf_MPI_Comm_rank(comm, &myid);

   /* check to see if we need to allocate more space in send_proc_obj for ids */
   if (send_proc_obj->length == send_proc_obj->storage_length)
   {
      send_proc_obj->storage_length += 20; /* add space for 20 more processors */
      send_proc_obj->id 
      = jxf_TReAlloc(send_proc_obj->id, JXF_Int, send_proc_obj->storage_length);
      send_proc_obj->vec_starts 
      = jxf_TReAlloc(send_proc_obj->vec_starts, JXF_Int, send_proc_obj->storage_length + 1);
   }
  
   /* initialize */ 
   count = send_proc_obj->length;
   index = send_proc_obj->vec_starts[count]; /* this is the number of elements */

   /* send proc */ 
   send_proc_obj->id[count] = contact_proc; 

   /* do we need more storage for the elements? */
   if (send_proc_obj->element_storage_length < index + contact_size)
   {
      elength = jxf_max(contact_size, 50);   
      elength += index;
      send_proc_obj->elements = jxf_TReAlloc(send_proc_obj->elements, JXF_Int, elength);
      send_proc_obj->element_storage_length = elength; 
   }
   
   /* populate send_proc_obj */
   for (i = 0; i < contact_size; i ++) 
   { 
      send_proc_obj->elements[index++] = recv_contact_buf[i];
   }
   send_proc_obj->vec_starts[count+1] = index;
   send_proc_obj->length ++;
   
   /* output - no message to return (confirmation) */
   *response_message_size = 0; 
  
   return jxf_error_flag;
}


/*!
 * \fn JXF_Int jxf_GetAssumedPartitionRowRange
 * \brief Assumed partition for IJ case. Given a particular processor id,
 *        return the assumed range of rows ([row_start, row_end]) for that processor.
 * \date 2011/09/03
 */
JXF_Int
jxf_GetAssumedPartitionRowRange( MPI_Comm comm,
                                JXF_Int  proc_id,
                                JXF_Int  global_first_row,
                                JXF_Int  global_num_rows,
                                JXF_Int *row_start,
                                JXF_Int *row_end )
{
   JXF_Int num_procs;
   JXF_Int size, extra;
   
   jxf_MPI_Comm_size(comm, &num_procs);

  /*----------------------------------------------------------------------- 
   *  this may look non-intuitive, but we have to be very careful that
   *  this function and the next are inverses - and avoiding overflow and
   *  rounding errors makes this difficult! 
   *----------------------------------------------------------------------*/

   size  = global_num_rows / num_procs;
   extra = global_num_rows - size*num_procs;

   *row_start = global_first_row + size*proc_id;
   *row_start += jxf_min(proc_id, extra);
   
   *row_end  =  global_first_row + size*(proc_id+1);
   *row_end +=  jxf_min(proc_id+1, extra);
   *row_end  = *row_end - 1;

   return jxf_error_flag;
}
