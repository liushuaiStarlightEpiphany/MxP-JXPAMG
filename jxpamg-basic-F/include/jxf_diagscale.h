//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_diagscale.h -- head files for diagscale preconditioner
 *  Date: 2013/10/27
 */ 

#ifndef JXF_DIAGSCALE_HEADER
#define JXF_DIAGSCALE_HEADER
 
#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

#ifndef JXF_MV_HEADER 
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER 
#include "jxf_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/diagscale/diagscale.c */
JXF_Int 
JXF_DiagScaleSetup( JXF_Solver        solver, 
                   JXF_hpCSRMatrix  hp_matrix );
JXF_Int
jxf_DiagScaleSetup( void             *solver, 
                   jxf_ParCSRMatrix  *par_matrix );
JXF_Int 
JXF_DiagScalePrecond( JXF_Solver       solver,
                     JXF_hpCSRMatrix hp_matrix,
                     JXF_ParVector    par_rhs,
                     JXF_ParVector    par_app  );
JXF_Int
jxf_DiagScalePrecond( void            *solver,
                     jxf_ParCSRMatrix *par_matrix,
                     jxf_ParVector    *par_rhs,
                     jxf_ParVector    *par_app  );

#endif
