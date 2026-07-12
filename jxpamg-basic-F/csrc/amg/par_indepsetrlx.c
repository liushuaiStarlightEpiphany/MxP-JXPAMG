//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_indepsetrlx.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGIndepSetRelaxInit
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGIndepSetRelaxInit( jxf_ParCSRMatrix *S,
                          JXF_Real          *measure_array )
{
   jxf_CSRMatrix    *S_diag = jxf_ParCSRMatrixDiag(S);
   JXF_Int              S_num_nodes = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int              first_index = jxf_ParCSRMatrixFirstRowIndex(S);
   JXF_Int              i;
   JXF_Int              ierr = 0;

   jxf_SeedRand(2747);
   
   for (i = 0; i < first_index; i ++)
   {
      jxf_Rand();
   }
   
   for (i = 0; i < S_num_nodes; i ++)
   {
      measure_array[i] += jxf_Rand();
   }

   return (ierr);
}


/*!
 * \fn JXF_Int jxf_PAMGIndepSetRelax
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGIndepSetRelax( jxf_ParCSRMatrix *S,
		      jxf_CSRMatrix    *S_ext,
                      JXF_Real          *measure_array,
                      JXF_Int             *graph_array,
                      JXF_Int              graph_array_size,
                      JXF_Int             *graph_array_offd,
                      JXF_Int              graph_array_offd_size,
                      JXF_Int             *IS_marker,
                      JXF_Int             *IS_marker_offd,
                      JXF_Int              nc_per,
                      JXF_Real           measure_th )
{
   jxf_CSRMatrix *S_diag   = jxf_ParCSRMatrixDiag(S);
   JXF_Int          *S_diag_i = jxf_CSRMatrixI(S_diag);
   JXF_Int          *S_diag_j = jxf_CSRMatrixJ(S_diag);
   jxf_CSRMatrix *S_offd   = jxf_ParCSRMatrixOffd(S);
   JXF_Int          *S_offd_i = jxf_CSRMatrixI(S_offd);
   JXF_Int          *S_offd_j = NULL;
   JXF_Int          *S_ext_i  = NULL;
   JXF_Int          *S_ext_j  = NULL;

   JXF_Int		 local_num_vars = jxf_CSRMatrixNumRows(S_diag);
   JXF_Int		 jc;
   JXF_Int           i, j, ig, jS, jj; 
   JXF_Int           num_per;
                   
   JXF_Int           ierr = 0;

  /*-------------------------------------------------------
   * Initialize IS_marker by putting all nodes in
   * the independent set.
   *-------------------------------------------------------*/

   if (jxf_CSRMatrixNumCols(S_offd))
   {
      S_offd_j = jxf_CSRMatrixJ(S_offd);
      S_ext_i  = jxf_CSRMatrixI(S_ext);
      S_ext_j  = jxf_CSRMatrixJ(S_ext);
   }

   for (ig = 0; ig < graph_array_size; ig ++)
   {
      i = graph_array[ig];
      if (measure_array[i] > 1)
      {
         IS_marker[i] = 1;
      }
   }
   
   for (ig = 0; ig < graph_array_offd_size; ig ++)
   {
      i = graph_array_offd[ig];
      if (measure_array[i+local_num_vars] > 1)
      {
         IS_marker_offd[i] = 1;
      }
   }
   
   
  /*-------------------------------------------------------
   * Remove nodes from the initial independent set
   *-------------------------------------------------------*/
   
   for (ig = 0; ig < graph_array_size; ig ++)
   {
      i = graph_array[ig];
      if (measure_array[i] < measure_th) 
      {
      	 IS_marker[i] = 0;      	
      }
      else if (measure_array[i] > 1)
      {
         for (jS = S_diag_i[i]; jS < S_diag_i[i+1]; jS ++)
         {
            j = S_diag_j[jS];
	    if (j < 0) 
	    {
	       j = - j - 1;
	    }
            
            /* only consider valid graph edges */
            if (measure_array[j] > 1) 
            {
               if (measure_array[i] > measure_array[j])
               {
                  IS_marker[j] = 0;                  
               }
               else if (measure_array[j] > measure_array[i])
               {
                  IS_marker[i] = 0;
               }
            }
         }
         
         for (jS = S_offd_i[i]; jS < S_offd_i[i+1]; jS ++)
         {
            jj = S_offd_j[jS];
            if (jj < 0) 
            {
               jj = - jj - 1;
            }
            j = local_num_vars + jj;
            
            /* only consider valid graph edges */
            if (measure_array[j] > 1) 
            {
               if (measure_array[i] > measure_array[j])
               {
                  IS_marker_offd[jj] = 0;
               }
               else if (measure_array[j] > measure_array[i])
               {
                  IS_marker[i] = 0;
               }
            }
         } // end for jS
      }
   } // end for ig
   
   
   for (ig = 0; ig < graph_array_offd_size; ig ++)
   {
      i = graph_array_offd[ig];
      if (measure_array[local_num_vars+i] < measure_th) 
      {
      	 IS_marker_offd[i] = 0;
      }
      else if (measure_array[local_num_vars+i] > 1)
      {
         for (jS = S_ext_i[i]; jS < S_ext_i[i+1]; jS ++)
         {
            j = S_ext_j[jS];
            
            if (j >= 0)
	    {
               /* only consider valid graph edges */
               if (measure_array[j] > 1) 
               {
                  if (measure_array[i+local_num_vars] > measure_array[j])
                  {
                     IS_marker[j] = 0;
                  }
                  else if (measure_array[j] > measure_array[i+local_num_vars])
                  {
                     IS_marker_offd[i] = 0;
                  }
               }
	    }
	    else
	    {
	       jc = -j - 1 + local_num_vars;
	     
               /* only consider valid graph edges */
               if (measure_array[jc] > 1)
               {
                  if (measure_array[i+local_num_vars] > measure_array[jc])
                  {
                     IS_marker_offd[-j-1] = 0;
                  }
                  else if (measure_array[jc] > measure_array[i+local_num_vars])
                  {
                     IS_marker_offd[i] = 0;
                  }
               }
	    }
         }
      }
   }

   num_per = 0;
   for (ig = 0; ig < graph_array_size; ig ++)
   {
      i = graph_array[ig];
      if (IS_marker[i] == 1) 
      {
         num_per ++;
      }
   }
   
   while (num_per > nc_per)
   {       
       for (ig = graph_array_size-1; ig >= 0; ig --)
       {
       	  i = graph_array[ig];
       	  if (IS_marker[i] == 1) 
       	  {
       	      IS_marker[i] = 0;
       	      num_per --;
       	      if (num_per <= nc_per) 
       	      {
       	         break;
       	      }
       	  }
       }
       break;
   }
 
   return (ierr);
}     
