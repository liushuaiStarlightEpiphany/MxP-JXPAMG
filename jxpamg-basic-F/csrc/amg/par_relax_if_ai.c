//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"

/*!
 * \fn JXF_Int jxf_PAMGRelaxIFAI
 */
JXF_Int  
jxf_hpPAMGRelaxIFAI( jxf_hpCSRMatrix *hp_matrix,
                  jxf_ParVector    *par_rhs,
                  JXF_Int             *cf_marker,
                  JXF_Int             *relax_marker,
                  JXF_Int              relax_type,
                  JXF_Int              relax_order,
                  JXF_Int              cycle_type,
                  JXF_Real           relax_weight,
                  JXF_Real           omega,
                  jxf_ParVector    *par_app,
                  jxf_ParVector    *Vtemp )
{
   JXF_Int i, Solve_err_flag = 0;
   JXF_Int relax_points[2];

   if (relax_order == 1 && cycle_type < 3)
   {
       if (cycle_type < 2)
       {
          relax_points[0] =  1;
          relax_points[1] = -1;
       }
       else
       {
	  relax_points[0] = -1;
	  relax_points[1] =  1;
       }

       for (i = 0; i < 2; i ++)
       {
          Solve_err_flag = jxf_hpPAMGRelaxAI( hp_matrix,
                                           par_rhs,
                                           cf_marker,
                                           relax_marker,
                                           relax_type,
                                           relax_points[i],
                                           relax_weight,
                                           omega,
                                           par_app,
                                           Vtemp ); 
       }
   }
   else
   {
       Solve_err_flag = jxf_hpPAMGRelaxAI( hp_matrix,
                                        par_rhs,
                                        cf_marker,
                                        relax_marker,
                                        relax_type,
                                        0,
                                        relax_weight,
                                        omega,
                                        par_app,
                                        Vtemp ); 
   }

   return Solve_err_flag;
}
