//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jx_util.h"

/*!
 * \fn JX_Int jx_BinarySearch
 * \brief performs a binary search for value on array list where 
 *        list needs to contain ordered nonnegative numbers.
 * \param *list pointer to the array to be searched. 
 * \param value the component to be searched.
 * \param list_length size of the array.
 * \note the routine returns the location of the value or -1
 * \date 2011/09/01 
 */ 
JX_Int 
jx_BinarySearch( JX_Int *list, JX_Int value, JX_Int list_length )
{
   JX_Int low, high, m;
   JX_Int not_found = 1;

   low = 0;
   high = list_length - 1; 
   while (not_found && low <= high)
   {
      m = (low + high) / 2;
      if (value < list[m])
      {
         high = m - 1;
      }
      else if (value > list[m])
      {
         low = m + 1;
      }
      else
      {
         not_found = 0;
         return m;
      }
   }
   return -1;
}

/*
 * Equivalent to C++ std::lower_bound
 */
JX_Int *
jx_LowerBound( JX_Int *first, JX_Int *last, JX_Int value )
{
   JX_Int *it;
   JX_Int count = last - first, step;
   while (count > 0) {
      it = first; step = count/2; it += step;
      if (*it < value) {
         first = ++it;
         count -= step + 1;
      }
      else count = step;
   }
   return first;
}

JX_Int
jx_BinarySearch2( JX_Int *list, JX_Int value, JX_Int low, JX_Int high, JX_Int *spot )
{

   JX_Int m;
   while (low <= high)
   {
      m = low + (high - low)/2;
      if (value < list[m])
         high = m - 1;
      else if (value > list[m])
         low = m + 1;
      else
      {
         *spot = m;
         return m;
      }
   }
   /* not found (high = low-1) - so insert at low */
   *spot = low;
   return -1;
}
