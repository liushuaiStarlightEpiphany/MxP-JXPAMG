//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

#if JXF_USING_COMPLEX

#include <complex.h>

JXF_Complex
jxf_conj( JXF_Complex value )
{
   return conj(value);
}

JXF_Real
jxf_cabs( JXF_Complex value )
{
   return cabs(value);
}

JXF_Real
jxf_creal( JXF_Complex value )
{
   return creal(value);
}

JXF_Real
jxf_cimag( JXF_Complex value )
{
   return cimag(value);
}

#endif
