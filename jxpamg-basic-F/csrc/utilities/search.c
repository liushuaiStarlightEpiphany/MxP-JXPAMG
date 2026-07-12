//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

/*!
 * \fn JXF_Int jxf_BinarySearch
 * \brief performs a binary search for value on array list where 
 *        list needs to contain ordered nonnegative numbers.
 * \param *list pointer to the array to be searched. 
 * \param value the component to be searched.
 * \param list_length size of the array.
 * \note the routine returns the location of the value or -1
 * \date 2011/09/01 
 */ 
JXF_Int 
jxf_BinarySearch( JXF_Int *list, JXF_Int value, JXF_Int list_length )
{
   JXF_Int low, high, m;
   JXF_Int not_found = 1;

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
JXF_Int *
jxf_LowerBound( JXF_Int *first, JXF_Int *last, JXF_Int value )
{
   JXF_Int *it;
   JXF_Int count = last - first, step;
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

JXF_Int
jxf_BinarySearch2( JXF_Int *list, JXF_Int value, JXF_Int low, JXF_Int high, JXF_Int *spot )
{

   JXF_Int m;
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
