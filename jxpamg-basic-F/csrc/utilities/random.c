//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

/* Static variables */
static JXF_Int Seed = 13579;

#define a  16807
#define m  2147483647
#define q  127773
#define r  2836

/*!
 * \fn void jxf_SeedRand
 * \brief Initializes the pseudo-random number generator 
 *        to a place in the sequence.
 * \param seed an JXF_Int containing the seed for the RNG.
 * \date 2011/09/01 
 */ 
void  
jxf_SeedRand( JXF_Int seed )
{
   Seed = seed;
}

/*!
 * \fn JXF_Real jxf_Rand
 * \brief Computes the next pseudo-random number in the sequence 
 *        using the global variable Seed.
 * \return a JXF_Real containing the next number in the sequence divided by
 *         2147483647 so that the numbers are in (0, 1].
 * \date 2011/09/01 
 */
JXF_Real 
jxf_Rand()
{
   JXF_Int  low, high, test;

   high = Seed / q;
   low  = Seed % q;
   test = a * low - r * high;
   if(test > 0)
   {
      Seed = test;
   }
   else
   {
      Seed = test + m;
   }

   return ((JXF_Real)(Seed) / m);
}
