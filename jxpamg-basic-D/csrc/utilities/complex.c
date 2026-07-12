//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jx_util.h"

#if JX_USING_COMPLEX

#include <complex.h>

JX_Complex
jx_conj( JX_Complex value )
{
   return conj(value);
}

JX_Real
jx_cabs( JX_Complex value )
{
   return cabs(value);
}

JX_Real
jx_creal( JX_Complex value )
{
   return creal(value);
}

JX_Real
jx_cimag( JX_Complex value )
{
   return cimag(value);
}

#endif
