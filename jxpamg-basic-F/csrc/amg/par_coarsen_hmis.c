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
jxf_PAMGCoarsenHMIS( jxf_ParCSRMatrix  *par_S,
                    jxf_ParCSRMatrix  *par_matrix,
                    JXF_Int               measure_type,
                    JXF_Int               debug_flag,
                    JXF_Int             **CF_marker_ptr )
{
   JXF_Int ierr = 0;

  /*--------------------------------------------------------
   * Perform Ruge coarsening followed by PMIS coarsening
   *--------------------------------------------------------*/

   ierr += jxf_PAMGCoarsenRuge (par_S, par_matrix, measure_type, 10, debug_flag, CF_marker_ptr);

   ierr += jxf_PAMGCoarsenPMIS (par_S, par_matrix, 1, debug_flag, CF_marker_ptr);

   return (ierr);
}
