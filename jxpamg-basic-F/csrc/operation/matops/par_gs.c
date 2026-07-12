//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_gs.c
 *  Date: 2011/09/03
 */ 

#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_gselim
 * \brief Gaussian Elimination
 * \date 2011/09/03
 */
JXF_Int 
jxf_gselim( JXF_Real *A_matrix, JXF_Real *x_vector, JXF_Int size )
{
   JXF_Int    err_flag = 0;
   JXF_Int    j,k,m;
   JXF_Real factor;
   
   if (size == 1) /* A_matrix is 1x1 */  
   {
      if (A_matrix[0] != 0.0)
      {
         x_vector[0] = x_vector[0] / A_matrix[0];
         return(err_flag);
      }
      else
      {
         err_flag = 1;
         return(err_flag);
      }
   }
   else /* A_matrix is nxn.  Forward elimination */ 
   {
      for (k = 0; k < size-1; k ++)
      {
          if (A_matrix[k*size+k] != 0.0)
          {          
             for (j = k+1; j < size; j ++)
             {
                 if (A_matrix[j*size+k] != 0.0)
                 {
                    factor = A_matrix[j*size+k] / A_matrix[k*size+k];
                    for (m = k+1; m < size; m ++)
                    {
                        A_matrix[j*size+m] -= factor * A_matrix[k*size+m];
                    }             
                    x_vector[j] -= factor * x_vector[k];   /* Elimination step for rhs */          
                 }
             }
          }
       }
       /* Back Substitution */
       for (k = size-1; k > 0; -- k)
       {
           x_vector[k] /= A_matrix[k*size+k];
           for (j = 0; j < k; j ++)
           {
               if (A_matrix[j*size+k] != 0.0)
               {
                  x_vector[j] -= x_vector[k] * A_matrix[j*size+k];
               }
           }
       }
       x_vector[0] /= A_matrix[0];
       return(err_flag);
    }
}
