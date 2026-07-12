//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

/*!
 * \fn void jxf_swap
 * \brief Swap two components of an integer-array.
 * \param *v pointer to the array to be swapped. 
 * \param i one of the index to be swapped.
 * \param j the other index to be swapped.
 * \date 2011/09/01 
 */  
void 
jxf_swap( JXF_Int *v, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}

/*!
 * \fn void jxf_swap1
 * \brief Swap two components of a JXF_Real-array.
 * \param *v pointer to the array to be swapped. 
 * \param i one of the index to be swapped.
 * \param j the other index to be swapped.
 * \date 2011/09/01 
 */ 
void 
jxf_swap1( JXF_Real *v, JXF_Int i, JXF_Int j )
{
   JXF_Real temp;
   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}

/*!
 * \fn void jxf_swap2
 * \brief Swap two components of an integer-array and a 
 *        JXF_Real-array at the same time.
 * \param *v pointer to the integer-array to be swapped. 
 * \param *w pointer to the JXF_Real-array to be swapped.  
 * \param i one of the index to be swapped.
 * \param j the other index to be swapped.
 * \date 2011/09/01 
 */   
void 
jxf_swap2( JXF_Int *v, JXF_Real *w, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;
   JXF_Real temp2;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp2 = w[i];
   w[i] = w[j];
   w[j] = temp2;
}

/*!
 * \fn void jxf_swap2i
 * \brief Swap two components of two integer-arrays at the same time.
 * \param *v pointer to the integer-array to be swapped. 
 * \param *w pointer to the other integer-array to be swapped.  
 * \param i one of the index to be swapped.
 * \param j the other index to be swapped.
 * \date 2011/09/01 
 */   
void 
jxf_swap2i( JXF_Int *v, JXF_Int *w, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp = w[i];
   w[i] = w[j];
   w[j] = temp;
}

void
jxf_swap3_d_perm( JXF_Int *v, JXF_Real *w, JXF_Int *z, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;
   JXF_Real temp_d;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp_d = w[i];
   w[i] = w[j];
   w[j] = temp_d;
   temp = z[v[i]];
   z[v[i]] = z[v[j]];
   z[v[j]] = temp;
}

/*!
 * \fn void jxf_swap3i
 * \brief Swap two components of three integer-arrays at the same time.
 * \param *v pointer to the first integer-array to be swapped. 
 * \param *w pointer to the second integer-array to be swapped. 
 * \param *z pointer to the third integer-array to be swapped.  
 * \param i one of the index to be swapped.
 * \param j the other index to be swapped.
 * \date 2011/09/01 
 */   
void 
jxf_swap3i( JXF_Int *v, JXF_Int *w, JXF_Int *z, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;

   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
   temp = w[i];
   w[i] = w[j];
   w[j] = temp;
   temp = z[i];
   z[i] = z[j];
   z[j] = temp;
}

/*!
 * \fn void jxf_qsort0
 * \brief Sort parts components of an integer-array ascendingly.
 * \param *v pointer to the integer-array to be sorted. 
 * \param left starting index of the components to be sorted.
 * \param right ending index of the components to be sorted.
 * \date 2011/09/01 
 */   
void 
jxf_qsort0( JXF_Int *v, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap(v, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (v[i] < v[left])
      {
         jxf_swap(v, ++last, i);
      }
   }
   jxf_swap(v, left, last);
   jxf_qsort0(v, left, last-1);
   jxf_qsort0(v, last+1, right);
}

/*!
 * \fn void jxf_qsortd
 * \brief Sort parts components of a JXF_Real-array ascendingly.
 * \param *v pointer to the JXF_Real-array to be sorted. 
 * \param left starting index of the components to be sorted.
 * \param right ending index of the components to be sorted.
 * \date 2011/09/01 
 */ 
void 
jxf_qsortd( JXF_Real *v, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap1(v, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (v[i] < v[left])
      {
         jxf_swap1(v, ++last, i);
      }
   }
   jxf_swap1(v, left, last);
   jxf_qsortd(v, left, last-1);
   jxf_qsortd(v, last+1, right);
}

/*!
 * \fn void jxf_qsort1 (sort on w, move v)
 * \brief Sort the array 'v' ascendingly with quick sorting algorithm, 
 *        and meanwhile 'w' is sorted in the same way.
 * \param *v pointer to the JXF_Int-array
 * \param *w pointer to the JXF_Real-array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \date 2011/09/01  
 */  
void 
jxf_qsort1( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap2(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (v[i] < v[left])
      {
         jxf_swap2(v, w, ++last, i);
      }
   }   
   jxf_swap2(v, w, left, last);
   jxf_qsort1(v, w, left, last-1);
   jxf_qsort1(v, w, last+1, right);
}

/*!
 * \fn void jxf_qsort2i (sort on v, move w)
 * \brief Sort the array 'v' ascendingly with quick sorting algorithm, 
 *        and meanwhile 'w' is sorted in the same way.
 * \param *v pointer to the JXF_Int-array
 * \param *w pointer to the JXF_Int-array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \date 2011/09/01  
 */   
void 
jxf_qsort2i( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap2i(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (v[i] < v[left])
      {
         jxf_swap2i(v, w, ++last, i);
      }
   }
   jxf_swap2i(v, w, left, last);
   jxf_qsort2i(v, w, left, last-1);
   jxf_qsort2i(v, w, last+1, right);
}

/*!
 * \fn void jxf_qsort2 (sort on w, move v)
 * \brief Sort the array 'w' ascendingly with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int-array
 * \param *w pointer to the JXF_Real-array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \date 2011/09/01  
 */  
void 
jxf_qsort2( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap2(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (w[i] < w[left])
      {
         jxf_swap2(v, w, ++last, i);
      }
   }
   jxf_swap2(v, w, left, last);
   jxf_qsort2(v, w, left, last-1);
   jxf_qsort2(v, w, last+1, right);
}

/*!
 * \fn void jxf_qsort3i (sort on v, move w and z)
 * \brief Sort the array 'v' ascendingly with quick sorting algorithm, 
 *        and meanwhile 'w', 'z' are sorted in the same way.
 * \param *v pointer to the first JXF_Int-array
 * \param *w pointer to the second JXF_Int-array
 * \param *z pointer to the third JXF_Int-array 
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \date 2011/09/01  
 */  
void 
jxf_qsort3i( JXF_Int *v, JXF_Int *w, JXF_Int *z, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap3i(v, w, z, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (v[i] < v[left])
      {
         jxf_swap3i(v, w, z, ++last, i);
      }
   }
   jxf_swap3i(v, w, z, left, last);
   jxf_qsort3i(v, w, z, left, last-1);
   jxf_qsort3i(v, w, z, last+1, right);
}

void
jxf_qsort3ir( JXF_Int *v, JXF_Real *w, JXF_Int *z, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;

   if (left >= right)
   {
      return;
   }
   jxf_swap3_d_perm( v, w, z, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
   {
      if (v[i] < v[left])
      {
         jxf_swap3_d_perm(v, w, z, ++last, i);
      }
   }
   jxf_swap3_d_perm(v, w, z, left, last);
   jxf_qsort3ir(v, w, z, left, last-1);
   jxf_qsort3ir(v, w, z, last+1, right);
}

/*!
 * \fn void jxf_qsort2abs (sort on w, move v)
 * \brief Sort the array '|w|' ascendingly with quick sorting algorithm, 
 *        and meanwhile 'v' is sorted in the same way.
 * \param *v pointer to the JXF_Int-array
 * \param *w pointer to the JXF_Real-array
 * \param left the starting index
 * \param right the ending index
 * \note 'left' and 'right' are usually set to be 0 and n-1, 
 *        respectively, where n is the length of 'v' or 'w'. 
 * \date 2011/09/01  
 */ 
void 
jxf_qsort2abs( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;
   if (left >= right)
   {
      return;
   }
   jxf_swap2(v, w, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i ++)
   {
      if (fabs(w[i]) > fabs(w[left]))
      {
         jxf_swap2(v, w, ++last, i);
      }
   }
   jxf_swap2(v, w, left, last);
   jxf_qsort2abs(v, w, left, last-1);
   jxf_qsort2abs(v, w, last+1, right);
}

/*!
 * \fn void jxf_diQuickSort12
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
jxf_diQuickSort12( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right )
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
      jxf_diQuickSort12(v, w, left, j-1);
      jxf_diQuickSort12(v, w, i, right);
   }
}

/*!
 * \fn void jxf_diQuickSort21
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
jxf_diQuickSort21( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right )
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
      jxf_diQuickSort21(v, w, left, j-1);
      jxf_diQuickSort21(v, w, i, right);
   }
}

void
jxf_union2( JXF_Int n1, JXF_Int *arr1, JXF_Int n2, JXF_Int *arr2, JXF_Int *n3, JXF_Int *arr3, JXF_Int *map1, JXF_Int *map2 )
{
   JXF_Int i = 0, j = 0, k = 0;
   while (i < n1 && j < n2)
   {
      if (arr1[i] < arr2[j])
      {
         if (map1) { map1[i] = k; }
         arr3[k++] = arr1[i++];
      }
      else if (arr1[i] > arr2[j])
      {
         if (map2) { map2[j] = k; }
         arr3[k++] = arr2[j++];
      }
      else /* == */
      {
         if (map1) { map1[i] = k; }
         if (map2) { map2[j] = k; }
         arr3[k++] = arr1[i++];
         j++;
      }
   }
   while (i < n1)
   {
      if (map1) { map1[i] = k; }
      arr3[k++] = arr1[i++];
   }
   while (j < n2)
   {
      if (map2) { map2[j] = k; }
      arr3[k++] = arr2[j++];
   }
   *n3 = k;
}

void
jxf_swap3_d( JXF_Real *v, JXF_Int *w, JXF_Int *z, JXF_Int i, JXF_Int j )
{
   JXF_Int temp;
   JXF_Real temp_d;
   temp_d = v[i];
   v[i] = v[j];
   v[j] = temp_d;
   temp = w[i];
   w[i] = w[j];
   w[j] = temp;
   temp = z[i];
   z[i] = z[j];
   z[j] = temp;
}

void
jxf_qsort3( JXF_Real *v, JXF_Int *w, JXF_Int *z, JXF_Int left, JXF_Int right )
{
   JXF_Int i, last;
   if (left >= right)
      return;
   jxf_swap3_d( v, w, z, left, (left+right)/2);
   last = left;
   for (i = left+1; i <= right; i++)
      if (v[i] < v[left])
      {
         jxf_swap3_d(v,w, z, ++last, i);
      }
   jxf_swap3_d(v, w, z, left, last);
   jxf_qsort3(v, w, z, left, last-1);
   jxf_qsort3(v, w, z, last+1, right);
}
