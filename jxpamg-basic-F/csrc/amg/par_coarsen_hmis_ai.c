//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_coarsen_hmis.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGCoarsenHMIS
 * \brief HMIS (RS + PMIS) coarsening routine.
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGCoarsenHMISAI( jxf_ParCSRMatrix  *par_S,
                    jxf_ParCSRMatrix  *par_matrix,
                    JXF_Real           *AI_measure,
                    JXF_Int               CF_init,
                    JXF_Int               measure_type,
                    JXF_Int               debug_flag,
                    JXF_Int             **CF_marker_ptr )
{
   JXF_Int ierr = 0;

  /*--------------------------------------------------------
   * Perform Ruge coarsening followed by PMIS coarsening
   *--------------------------------------------------------*/

   ierr += jxf_PAMGCoarsenRugeAI (par_S, par_matrix, AI_measure, CF_init, measure_type, 10, debug_flag, CF_marker_ptr);

   ierr += jxf_PAMGCoarsenPMISAI (par_S, par_matrix, AI_measure, 1, debug_flag, CF_marker_ptr);

   return (ierr);
}
