//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_coarsen_falgout.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGCoarsenFalgout
 * \brief Falgout(RS + CLJP) coarsening routine.
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGCoarsenFalgout( jxf_ParCSRMatrix  *par_S,
                       jxf_ParCSRMatrix  *par_matrix,
                       JXF_Int               measure_type,
                       JXF_Int               debug_flag,
                       JXF_Int             **CF_marker_ptr )
{
   JXF_Int ierr = 0;

  /*-------------------------------------------------------
   * Perform Ruge coarsening followed by CLJP coarsening
   *-------------------------------------------------------*/

   ierr += jxf_PAMGCoarsenRuge (par_S, par_matrix, measure_type, 6, debug_flag, CF_marker_ptr);

   ierr += jxf_PAMGCoarsen (par_S, par_matrix, 1, debug_flag, CF_marker_ptr);

   return (ierr);
}
