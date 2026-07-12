/*========================================================================*
 *  FSLS (Fast Solvers for Linear System)  (c) 2010-2012                  *
 *  School of Mathematics and Computational Science                       *
 *  Xiangtan University                                                   *
 *  Email: peghoty@163.com                                                *
 *========================================================================*/
 
/*!
 *  tool.c
 *
 *  Created by peghoty 2010/08/22
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 */ 

#include "jxf_multils.h"

/*!
 * \fn JXF_Int fsls_ArrayInitialize
 * \brief Zero out a JXF_Real array.
 * \param *x pointer to the JXF_Real array.
 * \param n size of the 'x' array.
 * \author peghoty
 * \date 2010/09/16
 */
void 
fsls_ArrayInitialize( JXF_Real *x, JXF_Int n )
{
   memset(x, 0X0, n*sizeof(JXF_Real));
}

/*!
 * \fn JXF_Int fsls_ArrayPrint
 * \brief write a vector to a given file
 * \param *u pointer to the vector array
 * \param n size of the 'u' array
 * \param *filename where the vector be saved
 * \author peghoty
 * \date 2009/11/29 
 */
JXF_Int 
fsls_ArrayPrint( JXF_Real *u, JXF_Int n, char *filename )
{
    JXF_Int j;
    FILE *fp = NULL;
   
    /* Save the vector */
    fp = fopen(filename, "w");
    jxf_fprintf(fp, "%d\n", n);
    for (j = 0; j < n; j ++)
    {
       jxf_fprintf(fp, "%.15le\n", u[j]);
    }
    fclose(fp);
    
    return (0);
}

/*!
 * \fn JXF_Int fsls_IntArrayPrint
 * \brief write a vector of integer type to a given file
 * \param *u pointer to the vector array
 * \param n size of the 'u' array
 * \param *filename where the vector be saved
 * \author peghoty
 * \date 2011/05/14 
 */
JXF_Int 
fsls_IntArrayPrint( JXF_Int *u, JXF_Int n, char *filename )
{
    JXF_Int j;
    FILE *fp = NULL;
   
    /* Save the vector */
    fp = fopen(filename, "w");
    jxf_fprintf(fp, "%d\n", n);
    for (j = 0; j < n; j ++)
    {
       jxf_fprintf(fp, "%d\n", u[j]);
    }
    fclose(fp);
    
    return (0);
}

/*!
 * \fn JXF_Int fsls_IntArrayIJPrint
 * \brief write a vector of integer type to a given file in IJ manner.
 * \param *u pointer to the vector array
 * \param n size of the 'u' array
 * \param *filename where the vector be saved
 * \author peghoty
 * \date 2011/05/18 
 */
JXF_Int 
fsls_IntArrayIJPrint( JXF_Int *u, JXF_Int n, char *filename )
{
    JXF_Int j;
    FILE *fp = NULL;
   
    /* Save the vector */
    fp = fopen(filename, "w");
    for (j = 0; j < n; j ++)
    {
       jxf_fprintf(fp, "%d %d\n", j+1, u[j]);
    }
    fclose(fp);
    
    return (0);
}

/*!
 * \fn void fsls_ArrayCopy( JXF_Real *x, JXF_Real *y, JXF_Int size )
 * \brief copy x to y ( x -> y ).
 * \author peghoty
 * \date 2010/04/29 
 */
void 
fsls_ArrayCopy( JXF_Int size, JXF_Real *x, JXF_Real *y )
{
   memcpy( y, x, size*sizeof(JXF_Real) );
}

/*!
 * \fn JXF_Real fsls_Arrayl2Norm
 * \brief Compute the l^2 norm of a given vector
 * \param *x pointer to the vector
 * \return l2norm, the l^2 norm of x
 * \author peghoty
 * \date 2009/12/30 
 */
JXF_Real   
fsls_Arrayl2Norm( JXF_Real *x, JXF_Int n )
{
   JXF_Int i;
   JXF_Real l2norm = 0.0;
   
   for (i = 0; i < n; i ++)
   {
      l2norm += x[i]*x[i];
   }
   l2norm = sqrt(l2norm); 
     
   return l2norm;
}

/*!
 * \fn JXF_Real fsls_ArrayAveL2Norm
 * \brief Compute the l^2 norm of a given vector
 * \param *x pointer to the vector
 * \return l2norm, the l^2 norm of x
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Real   
fsls_ArrayAveL2Norm( JXF_Real *x, JXF_Int n )
{
   JXF_Int i;
   JXF_Real l2norm = 0.0;
   
   for (i = 0; i < n; i ++)
   {
     l2norm += x[i]*x[i];
   }
   l2norm = sqrt(l2norm/(JXF_Real)n); 
     
   return l2norm;
}

/*!
 * \fn JXF_Int fsls_ArrayDoubleAbsMax
 * \brief Find the abs. maximal value of 'x'
 * \return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayDoubleAbsMax( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr )
{
   JXF_Int i,index;
   JXF_Real elm = 0.0;
   
   elm = fabs(x[0]); index = 0;
   for (i = 1; i < n; i ++)
   {
      if (fabs(x[i]) > elm)
      {
         elm = fabs(x[i]);
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayDoubleAbsMin
 * \brief Find the abs. minimal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayDoubleAbsMin( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr )
{
   JXF_Int i,index;
   JXF_Real elm = 0.0;
   
   elm = fabs(x[0]); index = 0;
   for (i = 1; i < n; i ++)
   {
      if (fabs(x[i]) < elm)
      {
         elm = fabs(x[i]);
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayDoubleMax
 * \brief Find the maximal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayDoubleMax( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr )
{
   JXF_Int i,index;
   JXF_Real elm = 0.0;
   
   elm = x[0]; index = 0;
   for (i = 1; i < n; i ++)
   {
      if (x[i] > elm)
      {
         elm = x[i];
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayDoubleMin
 * \brief Find the minimal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayDoubleMin( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr )
{
   JXF_Int i,index;
   JXF_Real elm = 0.0;
   
   elm = x[0]; index = 0;
   for (i = 1; i < n; i ++)
   {
      if (x[i] < elm)
      {
         elm = x[i];
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayIntAbsMax
 * \brief Find the abs. maximal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayIntAbsMax( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr )
{
   JXF_Int i,index;
   JXF_Int elm = 0;
   
   elm = abs(x[0]); index = 0;
   for (i = 1; i < n; i ++)
   {
      if (abs(x[i]) > elm)
      {
         elm = abs(x[i]);
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayIntAbsMin
 * \brief Find the abs. minimal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayIntAbsMin( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr )
{
   JXF_Int i,index;
   JXF_Int elm = 0;
   
   elm = abs(x[0]); index = 0;
   for (i = 1; i < n; i ++)
   {
      if (abs(x[i]) < elm)
      {
         elm = abs(x[i]);
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayIntMax
 * \brief Find the maximal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayIntMax( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr )
{
   JXF_Int i,index;
   JXF_Int elm = 0;
   
   elm = x[0]; index = 0;
   for (i = 1; i < n; i ++)
   {
      if (x[i] > elm)
      {
         elm = x[i];
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_ArrayIntMin
 * \brief Find the minimal value of 'x', and return the corresponding index
 * \author peghoty
 * \date 2010/08/26 
 */
JXF_Int
fsls_ArrayIntMin( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr )
{
   JXF_Int i,index;
   JXF_Int elm = 0;
   
   elm = x[0]; index = 0;
   for (i = 1; i < n; i ++)
   {
      if (x[i] < elm)
      {
         elm = x[i];
         index = i;
      } 
   }
   
   *elm_ptr = elm;
   return index;
}

/*!
 * \fn JXF_Int fsls_daxpy
 * \brief Compute y := y + alpha*x
 * \author peghoty
 * \date 2010/12/30 
 */
JXF_Int
fsls_daxpy( JXF_Int n, JXF_Real alpha, JXF_Real *x, JXF_Real *y )
{           
   JXF_Int i,ierr = 0;

   for (i = 0; i < n; i ++)
   {
      y[i] += alpha*x[i];
   }

   return ierr;
}

/*!
 * \fn fsls_gselim
 * \brief Gaussian Elimination
 * \param *A the pointer to the matrix to be solved
 * \param *x the pointer to the rhs vector at first and the solution at last
 * \param n the size of the matrix
 * \return 0 if success, 1 if fail
 * \date 2009/12/21
 */
JXF_Int 
fsls_gselim( JXF_Real *A, JXF_Real *x, JXF_Int n )
{
   JXF_Int    err_flag = 0;
   JXF_Int    j,k,m;
   JXF_Real factor = 0.0;
   
   if (n == 1)   /* A is 1x1 */  
   {
      if (A[0] != 0.0)
      {
         x[0] = x[0] / A[0];
         return(err_flag);
      }
      else
      {
         err_flag = 1;
         return(err_flag);
      }
   }
   else   /* A is nxn */   
   {
   
      /* Forward elimination */
      for (k = 0; k < n-1; k ++)
      {
          if (A[k*n+k] != 0.0)
          {          
             for (j = k+1; j < n; j ++)
             {
                 if (A[j*n+k] != 0.0)
                 {
                    factor = A[j*n+k] / A[k*n+k];
                    for (m = k+1; m < n; m ++)
                    {
                        A[j*n+m] -= factor * A[k*n+m];
                    }
                    /* Elimination step for rhs */ 
                    x[j] -= factor * x[k];              
                 }
             }
          }
       }
       
       /* Back Substitution */
       for (k = n-1; k > 0; -- k)
       {
           x[k] /= A[k*n+k];
           for (j = 0; j < k; j ++)
           {
               if (A[j*n+k] != 0.0)
               {
                  x[j] -= x[k] * A[j*n+k];
               }
           }
       }
       
       x[0] /= A[0];
       
       return(err_flag);
    }
}

/*!
 * \fn fsls_gselim_piv
 * \brief Gaussian Elimination with pivoting
 * \param *A the pointer to the matrix to be solved
 * \param *x the pointer to the rhs vector at first and the solution at last
 * \param n the size of the matrix
 * \return 0 if success, 1,-1 if fail
 * \date 2010/03/03
 */
JXF_Int 
fsls_gselim_piv( JXF_Real *A, JXF_Real *x, JXF_Int n )
{
   JXF_Int    err_flag = 0;
   JXF_Int    j,k,m, piv_row;
   JXF_Real factor, piv, tmp;
   JXF_Real eps = 1e-8;
   
   if (n == 1)  /* A is 1x1 */  
   {
      if (fabs(A[0]) >  1e-10)
      {
         x[0] = x[0] / A[0];
         return(err_flag);
      }
      else
      {
         err_flag = 1;
         return(err_flag);
      }
   }
   
   else  /* A is nxn.  Forward elimination */ 
   {
      for (k = 0; k < n-1; k ++)
      {
         /* we do partial pivoting for size */
         piv = A[k*n+k];
         piv_row = k;
         /* find the largest pivot in position k*/
         for (j = k+1; j < n; j ++)         
         {
            if (fabs(A[j*n+k]) > fabs(piv))
            {
               piv =  A[j*n+k];
               piv_row = j;
            }
         }
         if (piv_row != k)   /* do a row exchange  - rows k and piv_row*/
         {
            for (j = 0; j < n; j ++)
            {
               tmp = A[k*n + j];
               A[k*n + j] = A[piv_row*n + j];
               A[piv_row*n + j] = tmp;
            }
            tmp = x[k];
            x[k] = x[piv_row];
            x[piv_row] = tmp;
         }

         if (fabs(piv) > eps)
         {          
             for (j = k+1; j < n; j ++)
             {
                 if (A[j*n+k] != 0.0)
                 {
                    factor = A[j*n+k]/A[k*n+k];
                    for (m = k+1; m < n; m ++)
                    {
                        A[j*n+m] -= factor * A[k*n+m];
                    }
                    /* Elimination step for rhs */ 
                    x[j] -= factor * x[k];              
                 }
             }
         }
         else
         {
            jxf_printf("Matrix is nearly singular: zero pivot error\n");
            return(-1);
         }
      }
      
      /* we also need to check the pivot in the last row to see if it is zero */  
      k = n - 1; /* last row */
      if ( fabs(A[k*n+k]) < eps)
      {
         jxf_printf("Block of matrix is nearly singular: zero pivot error\n");
         return(-1);
      }

      /* Back Substitution  */
      for (k = n-1; k > 0; -- k)
      {
           x[k] /= A[k*n+k];
           for (j = 0; j < k; j ++)
           {
               if (A[j*n+k] != 0.0)
               {
                  x[j] -= x[k] * A[j*n+k];
               }
           }
      }
      x[0] /= A[0];
      return(err_flag);
   }
}

/*!
 * \fn fsls_Free
 * \brief free memory
 * \date 2010/01/03
 */
void
fsls_Free( char *ptr )
{
   if (ptr) free(ptr);
}

/*!
 * \fn fsls_CAlloc
 * \brief allocate memory
 * \date 2010/01/03
 */
char *
fsls_CAlloc( size_t count, size_t elt_size )
{
   char *ptr  = NULL;
   JXF_Int   size = count*elt_size;

   if (size > 0)
   {
      ptr = calloc(count, elt_size);
      if (ptr == NULL)
      {
         fsls_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

/*!
 * \fn fsls_OutOfMemory
 * \brief print information when out of memory
 * \date 2010/01/03
 */
JXF_Int
fsls_OutOfMemory( size_t size )
{
   jxf_printf("\n \033[31mOut of memory\033[00m trying to allocate \033[31m%d\033[00m bytes!\n", (JXF_Int) size);
   fflush(stdout);
   return 0;
}

/*!
 * \fn void fsls_ISwap
 * \brief swap the i-th and j-th element in the array 'w'(JXF_Int type)
 * \param *w pointer to the array
 * \param i one position in w
 * \param j the other position in w
 * \author peghoty
 * \date 2009/11/28 
 */
void 
fsls_ISwap(JXF_Int *w, JXF_Int i, JXF_Int j)
{
   JXF_Int temp;
   temp = w[i];
   w[i] = w[j];
   w[j] = temp;
}

/*!
 * \fn fsls_DSwap
 * \brief swap the i-th and j-th element in the array 'w'(JXF_Real type)
 * \param *w pointer to the array
 * \param i one position in w
 * \param j the other position in w  
 * \author peghoty 
 * \date 2009/11/28 
 */
void 
fsls_DSwap(JXF_Real *w, JXF_Int i, JXF_Int j)
{
   JXF_Real temp;
   temp = w[i];
   w[i] = w[j];
   w[j] = temp;
}

/*!
 * \fn void fsls_IDSwap
 * \brief swap the i-th and j-th element both in the array 
 *  'v'(JXF_Int type) and 'w'(JXF_Real type)
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Real array
 * \param i one position in the array
 * \param j the other position in in the array  
 * \author peghoty 
 * \date 2009/12/04
 */
void 
fsls_IDSwap( JXF_Int *v, JXF_Real *w, JXF_Int i, JXF_Int j )
{
   JXF_Int    temp;
   JXF_Real temp2;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp2 = w[i];
   w[i] = w[j];
   w[j] = temp2;
}

/*!
 * \fn void fsls_IISwap
 * \brief swap the i-th and j-th element both in the array 
 *  'v'(JXF_Int type) and 'w'(JXF_Int type)
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param i one position in the array
 * \param j the other position in in the array  
 * \author peghoty 
 * \date 2010/05/06 
 */
void 
fsls_IISwap( JXF_Int *v, JXF_Int *w, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;
   JXF_Int temp2;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp2 = w[i];
   w[i] = w[j];
   w[j] = temp2;
}

/*!
 * \fn void fsls_IQuickSort
 * \brief sort the array 'a'(JXF_Int type) with the quick sorting algorithm
 * \param order integer, 12 means ascendingly; 21 means descendingly
 * \param *a pointer to the array needed to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the 
 * length of 'a'.
 * \author peghoty
 * \date 2009/12/05 
 */
void
fsls_IQuickSort(JXF_Int order, JXF_Int *a, JXF_Int left, JXF_Int right)
{
   JXF_Int i, last;

   if (left >= right) return;
   
   fsls_ISwap(a, left, (left+right)/2);
   
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (a[i] < a[left]) fsls_ISwap(a, ++last, i);
      }
      else if(order == 21)
      {
         if (a[i] > a[left]) fsls_ISwap(a, ++last, i);
      }
   }
   
   fsls_ISwap(a, left, last);
   
   fsls_IQuickSort(order, a, left, last-1);
   fsls_IQuickSort(order, a, last+1, right);
}

/*!
 * \fn void fsls_DQuickSort
 * \brief sort the array 'a'(JXF_Real type) with the quick sorting algorithm
 * \param order integer, 12 means ascendingly; 21 means descendingly
 * \param *a pointer to the array needed to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the 
 * length of 'a'.  
 * \author peghoty
 * \date 2009/12/05  
 */
void
fsls_DQuickSort(JXF_Int order, JXF_Real *a, JXF_Int left, JXF_Int right)
{
   JXF_Int i, last;

   if (left >= right) return;
   
   fsls_DSwap(a, left, (left+right)/2);
   
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (a[i] < a[left]) fsls_DSwap(a, ++last, i);
      }
      else if (order == 21)
      {
         if (a[i] > a[left]) fsls_DSwap(a, ++last, i);
      }
   }
   
   fsls_DSwap(a, left, last);
   
   fsls_DQuickSort(order, a, left, last-1);
   fsls_DQuickSort(order, a, last+1, right);
}

/*!
 * \fn void fsls_IQuickSortIndex
 * \brief reorder the index of 'a'(JXF_Int type) so that 'a' is ascending 
 *  or descending in such order
 * \param order integer, 12 means ascendingly; 21 means descendingly  
 * \param *a pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'a'
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the 
 * length of 'a'. 'index' should be initialized in the nature order and it has the
 * same length as 'a'.   
 * \author peghoty
 * \date 2009/12/05  
 */
void
fsls_IQuickSortIndex(JXF_Int order, JXF_Int *a, JXF_Int left, JXF_Int right, JXF_Int *index)
{
   JXF_Int i, last;

   if (left >= right) return;
   
   fsls_ISwap(index, left, (left+right)/2);
   
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (a[index[i]] < a[index[left]]) 
           fsls_ISwap(index, ++last, i);
      }
      else if (order == 21)
      {
         if (a[index[i]] > a[index[left]]) 
           fsls_ISwap(index, ++last, i);
      }
   }
   
   fsls_ISwap(index, left, last);
   
   fsls_IQuickSortIndex(order, a, left, last-1, index);
   fsls_IQuickSortIndex(order, a, last+1, right, index);
}

/*!
 * \fn void fsls_DQuickSortIndex
 * \brief reorder the index of 'a'(JXF_Real type) so that 'a' is ascending 
 *  or descending in such order
 * \param order integer, 12 means ascendingly; 21 means descendingly   
 * \param *a pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'a'
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the 
 * length of 'a'. 'index' should be initialized in the nature order and it has the
 * same length as 'a'.   
 * \author peghoty
 * \date 2009/12/05  
 */
void
fsls_DQuickSortIndex(JXF_Int order, JXF_Real *a, JXF_Int left, JXF_Int right, JXF_Int *index)
{
   JXF_Int i, last;

   if (left >= right) return;
   
   fsls_ISwap(index, left, (left+right)/2);
   
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (a[index[i]] < a[index[left]])
           fsls_ISwap(index, ++last, i);
      }
      else if (order == 21)
      {
         if (a[index[i]] > a[index[left]])
           fsls_ISwap(index, ++last, i);
      }
   }
   
   fsls_ISwap(index, left, last);
   
   fsls_DQuickSortIndex(order, a, left, last-1, index);
   fsls_DQuickSortIndex(order, a, last+1, right, index);
}

/*!
 * \fn void fsls_IIQuickSort
 * \brief sort the array 'w' with the quick sorting 
 *  algorithm, and meanwhile 'v'(usually be the indexes array of 'w')  
 *  is sorted in the same way.
 * \param order integer, 12 means ascendingly; 21 means descendingly  
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *  respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty  
 * \date 2010/05/06
 */
void 
fsls_IIQuickSort(JXF_Int order, JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right)
{
   JXF_Int i, last;

   if (left >= right) return;
   fsls_IISwap(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (w[i] < w[left])
           fsls_IISwap(v, w, ++last, i);
      }
      else if (order == 21)
      {
         if (w[i] > w[left])
           fsls_IISwap(v, w, ++last, i);
      }
   }
   fsls_IISwap(v, w, left, last);
   fsls_IIQuickSort(order, v, w, left, last-1);
   fsls_IIQuickSort(order, v, w, last+1, right);
}

/*!
 * \fn void fsls_IDQuickSort
 * \brief sort the array 'w' with the quick sorting 
 *  algorithm, and meanwhile 'v'(usually be the indexes array of 'w')  
 *  is sorted in the same way.
 * \param order integer, 12 means ascendingly; 21 means descendingly  
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Real array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *  respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2009/12/05  
 */
void 
fsls_IDQuickSort( JXF_Int order, JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right) return;
   fsls_IDSwap(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (w[i] < w[left])
           fsls_IDSwap(v, w, ++last, i);
      }
      else if (order == 21)
      {
         if (w[i] > w[left])
           fsls_IDSwap(v, w, ++last, i);
      }
   }
   fsls_IDSwap(v, w, left, last);
   fsls_IDQuickSort(order, v, w, left, last-1);
   fsls_IDQuickSort(order, v, w, last+1, right);
}

/*!
 * \fn void fsls_DIQuickSort
 * \brief sort the array 'w'(JXF_Int type) with the quick sorting 
 *  algorithm, and meanwhile 'v'(JXF_Real type)  
 *  is sorted in the same way.
 * \param order integer, 12 means ascendingly; 21 means descendingly  
 * \param *v pointer to the JXF_Real array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *  respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2009/12/05  
 */
void 
fsls_DIQuickSort( JXF_Int order, JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right) return;
   fsls_IDSwap(w, v, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (w[i] < w[left])
           fsls_IDSwap(w, v, ++last, i);
      }
      else if (order == 21)
      {
         if (w[i] > w[left])
           fsls_IDSwap(w, v, ++last, i);
      }
   }
   fsls_IDSwap(w, v, left, last);
   fsls_DIQuickSort(order, v, w, left, last-1);
   fsls_DIQuickSort(order, v, w, last+1, right);
}

/*!
 * \fn void fsls_IDfabsQuickSort
 * \brief sort the array '|w|' with the quick sorting 
 *  algorithm, and meanwhile 'v'(usually be the indexes array of 'w')  
 *  is sorted in the same way.
 * \param order integer, 12 means ascendingly; 21 means descendingly  
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Real array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *  respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2009/12/05  
 */
void 
fsls_IDfabsQuickSort( JXF_Int order, JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right) return;
   fsls_IDSwap(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (order == 12)
      {
         if (fabs(w[i]) < fabs(w[left]))
           fsls_IDSwap(v, w, ++last, i);
      }
      else if (order == 21)
      {
         if (fabs(w[i]) > fabs(w[left]))
           fsls_IDSwap(v, w, ++last, i);
      }
   }
   fsls_IDSwap(v, w, left, last);
   fsls_IDfabsQuickSort(order, v, w, left, last-1);
   fsls_IDfabsQuickSort(order, v, w, last+1, right);
}

/*!
 * \fn void fsls_iQuickSort12
 * \brief Sort the array 'data'(JXF_Int type) with quick sorting algorithm
 * \param *data pointer to the array to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively,
 *       where n is the length of 'data'.  
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_iQuickSort12( JXF_Int *data, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[j] <= key) break;
         }
         if (i < j)
         {
            tmp = data[i]; data[i] = data[j]; data[j] = tmp;
         }
      }
      tmp = data[left]; data[left] = data[j]; data[j] = tmp;
      fsls_iQuickSort12(data, left, j-1);
      fsls_iQuickSort12(data, i, right);
   }
}

/*!
 * \fn void fsls_iQuickSortIndex12
 * \brief Reorder the index of 'data'(JXF_Int type) so that 'data' is  
 *        ascending in such order
 * \param *data pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'data'
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively, 
 *       where n is the length of 'data'. 'index', with the same length with
 *       'data', should be initialized in the nature order.   
 * \date 2010/12/10 
 * \author peghoty
 */
void 
fsls_iQuickSortIndex12( JXF_Int *data, JXF_Int left, JXF_Int right, JXF_Int *index )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[index[left]];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[index[i]] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[index[j]] <= key) break;
         }
         if (i < j)
         {
            tmp = index[i]; index[i] = index[j]; index[j] = tmp;
         }
      }
      tmp = index[left]; index[left] = index[j]; index[j] = tmp;
      fsls_iQuickSortIndex12(data, left, j-1, index);
      fsls_iQuickSortIndex12(data, i, right, index);
   }
}

/*!
 * \fn void fsls_dQuickSort12
 * \brief Sort the array 'data'(JXF_Real type) with quick sorting algorithm
 * \param *a pointer to the array to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively,
 *       where n is the length of 'data'.  
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_dQuickSort12( JXF_Real *data, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Real key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[j] <= key) break;
         }
         if (i < j)
         {
            tmp = data[i]; data[i] = data[j]; data[j] = tmp;
         }
      }
      tmp = data[left]; data[left] = data[j]; data[j] = tmp;
      fsls_dQuickSort12(data, left, j-1);
      fsls_dQuickSort12(data, i, right);
   }
}

/*!
 * \fn void fsls_dQuickSortIndex12
 * \brief Reorder the index of 'data'(JXF_Real type) so that 'data' is  
 *        ascending in such order
 * \param *data pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'data'
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively, 
 *       where n is the length of 'data'. 'index', with the same length with
 *       'data', should be initialized in the nature order.   
 * \date 2010/12/10 
 * \author peghoty
 */
void 
fsls_dQuickSortIndex12( JXF_Real *data, JXF_Int left, JXF_Int right, JXF_Int *index )
{
   JXF_Int i,j;
   JXF_Real key;
   JXF_Int tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[index[left]];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[index[i]] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[index[j]] <= key) break;
         }
         if (i < j)
         {
            tmp = index[i]; index[i] = index[j]; index[j] = tmp;
         }
      }
      tmp = index[left]; index[left] = index[j]; index[j] = tmp;
      fsls_dQuickSortIndex12(data, left, j-1, index);
      fsls_dQuickSortIndex12(data, i, right, index);
   }
}

/*!
 * \fn void fsls_iiQuickSort12
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v'(usually be the indexes array of 'w')  
 *        is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_iiQuickSort12( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] <= key) break;
         }
         if (i < j)
         {
            tmp = w[i]; w[i] = w[j]; w[j] = tmp;
            tmp = v[i]; v[i] = v[j]; v[j] = tmp;
         }
      }
      tmp = w[left]; w[left] = w[j]; w[j] = tmp;
      tmp = v[left]; v[left] = v[j]; v[j] = tmp;
      fsls_iiQuickSort12(v, w, left, j-1);
      fsls_iiQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void fsls_ddQuickSort12
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_ddQuickSort12( JXF_Real *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Real key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] <= key) break;
         }
         if (i < j)
         {
            tmp = w[i]; w[i] = w[j]; w[j] = tmp;
            tmp = v[i]; v[i] = v[j]; v[j] = tmp;
         }
      }
      tmp = w[left]; w[left] = w[j]; w[j] = tmp;
      tmp = v[left]; v[left] = v[j]; v[j] = tmp;
      fsls_ddQuickSort12(v, w, left, j-1);
      fsls_ddQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void fsls_idQuickSort12
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_idQuickSort12( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int tmpv;
   JXF_Real key,tmpw;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] <= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_idQuickSort12(v, w, left, j-1);
      fsls_idQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void fsls_idFabsQuickSort12
 * \brief Sort the array '|w|' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/20  
 */
void 
fsls_idFabsQuickSort12( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int tmpv;
   JXF_Real key,tmpw;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = fabs(w[left]);
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (fabs(w[i]) >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (fabs(w[j]) <= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_idFabsQuickSort12(v, w, left, j-1);
      fsls_idFabsQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void fsls_diQuickSort12
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_diQuickSort12( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmpw;
   JXF_Real tmpv;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] >= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] <= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_diQuickSort12(v, w, left, j-1);
      fsls_diQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void fsls_iQuickSort21
 * \brief Sort the array 'data'(JXF_Int type) with quick sorting algorithm
 * \param *data pointer to the array to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively,
 *       where n is the length of 'data'.  
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_iQuickSort21( JXF_Int *data, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[j] >= key) break;
         }
         if (i < j)
         {
            tmp = data[i]; data[i] = data[j]; data[j] = tmp;
         }
      }
      tmp = data[left]; data[left] = data[j]; data[j] = tmp;
      fsls_iQuickSort21(data, left, j-1);
      fsls_iQuickSort21(data, i, right);
   }
}

/*!
 * \fn void fsls_iQuickSortIndex21
 * \brief Reorder the index of 'data'(JXF_Int type) so that 'data' is  
 *        descending in such order
 * \param *data pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'data'
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively, 
 *       where n is the length of 'data'. 'index', with the same length with
 *       'data', should be initialized in the nature order.   
 * \date 2010/12/10 
 * \author peghoty
 */
void 
fsls_iQuickSortIndex21( JXF_Int *data, JXF_Int left, JXF_Int right, JXF_Int *index )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[index[left]];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[index[i]] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[index[j]] >= key) break;
         }
         if (i < j)
         {
            tmp = index[i]; index[i] = index[j]; index[j] = tmp;
         }
      }
      tmp = index[left]; index[left] = index[j]; index[j] = tmp;
      fsls_iQuickSortIndex21(data, left, j-1, index);
      fsls_iQuickSortIndex21(data, i, right, index);
   }
}

/*!
 * \fn void fsls_dQuickSort21
 * \brief Sort the array 'data'(JXF_Real type) with quick sorting algorithm
 * \param *a pointer to the array to be sorted
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively,
 *       where n is the length of 'data'.  
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_dQuickSort21( JXF_Real *data, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Real key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[j] >= key) break;
         }
         if (i < j)
         {
            tmp = data[i]; data[i] = data[j]; data[j] = tmp;
         }
      }
      tmp = data[left]; data[left] = data[j]; data[j] = tmp;
      fsls_dQuickSort21(data, left, j-1);
      fsls_dQuickSort21(data, i, right);
   }
}

/*!
 * \fn void fsls_dQuickSortIndex21
 * \brief Reorder the index of 'data'(JXF_Real type) so that 'data' is  
 *        descending in such order
 * \param *data pointer to the array 
 * \param left the starting index
 * \param right the ending index
 * \param *index the index of 'data'
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively, 
 *       where n is the length of 'data'. 'index', with the same length with
 *       'data', should be initialized in the nature order.   
 * \date 2010/12/10 
 * \author peghoty
 */
void 
fsls_dQuickSortIndex21( JXF_Real *data, JXF_Int left, JXF_Int right, JXF_Int *index )
{
   JXF_Int i,j;
   JXF_Real key;
   JXF_Int tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = data[index[left]];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (data[index[i]] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (data[index[j]] >= key) break;
         }
         if (i < j)
         {
            tmp = index[i]; index[i] = index[j]; index[j] = tmp;
         }
      }
      tmp = index[left]; index[left] = index[j]; index[j] = tmp;
      fsls_dQuickSortIndex21(data, left, j-1, index);
      fsls_dQuickSortIndex21(data, i, right, index);
   }
}

/*!
 * \fn void fsls_iiQuickSort21
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v'(usually be the indexes array of 'w')  
 *        is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_iiQuickSort21( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] >= key) break;
         }
         if (i < j)
         {
            tmp = w[i]; w[i] = w[j]; w[j] = tmp;
            tmp = v[i]; v[i] = v[j]; v[j] = tmp;
         }
      }
      tmp = w[left]; w[left] = w[j]; w[j] = tmp;
      tmp = v[left]; v[left] = v[j]; v[j] = tmp;
      fsls_iiQuickSort21(v, w, left, j-1);
      fsls_iiQuickSort21(v, w, i, right);
   }
}

/*!
 * \fn void fsls_ddQuickSort21
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_ddQuickSort21( JXF_Real *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Real key,tmp;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] >= key) break;
         }
         if (i < j)
         {
            tmp = w[i]; w[i] = w[j]; w[j] = tmp;
            tmp = v[i]; v[i] = v[j]; v[j] = tmp;
         }
      }
      tmp = w[left]; w[left] = w[j]; w[j] = tmp;
      tmp = v[left]; v[left] = v[j]; v[j] = tmp;
      fsls_ddQuickSort21(v, w, left, j-1);
      fsls_ddQuickSort21(v, w, i, right);
   }
}

/*!
 * \fn void fsls_idQuickSort21
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_idQuickSort21( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int tmpv;
   JXF_Real key,tmpw;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] >= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_idQuickSort21(v, w, left, j-1);
      fsls_idQuickSort21(v, w, i, right);
   }
}

/*!
 * \fn void fsls_idFabsQuickSort21
 * \brief Sort the array '|w|' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/20  
 */
void 
fsls_idFabsQuickSort21( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int tmpv;
   JXF_Real key,tmpw;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = fabs(w[left]);
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (fabs(w[i]) <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (fabs(w[j]) >= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_idFabsQuickSort21(v, w, left, j-1);
      fsls_idFabsQuickSort21(v, w, i, right);
   }
}

/*!
 * \fn void fsls_diQuickSort21
 * \brief Sort the array 'w' with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int array
 * \param *w pointer to the JXF_Int array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \author peghoty
 * \date 2010/12/10  
 */
void 
fsls_diQuickSort21( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i,j;
   JXF_Int key,tmpw;
   JXF_Real tmpv;
   if (left < right)
   {
      i = left;
      j = right + 1;
      key = w[left];
      while (i < j)
      {
         for (i = i+1; i < right; i ++)
         {
            if (w[i] <= key) break;
         }
         for (j = j-1; j > left; j --)
         {
            if (w[j] >= key) break;
         }
         if (i < j)
         {
            tmpw = w[i]; w[i] = w[j]; w[j] = tmpw;
            tmpv = v[i]; v[i] = v[j]; v[j] = tmpv;
         }
      }
      tmpw = w[left]; w[left] = w[j]; w[j] = tmpw;
      tmpv = v[left]; v[left] = v[j]; v[j] = tmpv;
      fsls_diQuickSort21(v, w, left, j-1);
      fsls_diQuickSort21(v, w, i, right);
   }
}

/*-----------------------------------
 * Static variables
 *----------------------------------*/
 
static JXF_Int Seed = 13579;

#define a  16807
#define m  2147483647
#define q  127773
#define r  2836

/*!
 * \fn void fsls_SeedRand
 * \brief Initializes the pseudo-random number 
 *        generator to a place in the sequence.
 * \param seed an JXF_Int containing the seed for the RNG.
 * \date 2010/08/29
 */
void  
fsls_SeedRand( JXF_Int seed )
{
   Seed = seed;
}


/*!
 * \fn JXF_Real fsls_Rand
 * \brief Computes the next pseudo-random number in 
 *        the sequence using the global variable Seed.
 * \return a JXF_Real containing the next number 
 *         in the sequence divided by 2147483647
 *         so that the numbers are in (0, 1]
 * \date 2010/08/29
 */
JXF_Real  
fsls_Rand()
{
   JXF_Int low, high, test;

   high = Seed / q;
   low  = Seed % q;
   test = a * low - r * high;
   
   if (test > 0)
   {
      Seed = test;
   }
   else
   {
      Seed = test + m;
   }
   
   return ((JXF_Real)(Seed) / m);
}

#undef a  
#undef m  
#undef q  
#undef r 

