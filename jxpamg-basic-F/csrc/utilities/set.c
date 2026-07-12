//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

JXF_Int jxf__global_mvcpu_flag = 0;

/*!
 * \fn void jxf_set_mvcpu_handler
 * \brief Reset jxf__global_mvcpu_flag.
 * \date 2018/01/16 
 */  
void
jxf_set_mvcpu_handler( JXF_Int new_mvcpu_flag )
{
   jxf__global_mvcpu_flag = new_mvcpu_flag;
}
