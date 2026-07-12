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
JXF_Int  
jxf_hpPAMGRelax( jxf_hpCSRMatrix *hp_matrix,
              jxf_ParVector    *par_rhs,
              JXF_Int             *cf_marker,
              JXF_Int              relax_type,
              JXF_Int              relax_points,
              JXF_Real           relax_weight,
              JXF_Real           omega,
              JXF_Real          *l1_norms,
              jxf_ParVector    *par_app,
              jxf_ParVector    *Vtemp,
              jxf_ParVector    *Ztemp )
{
   JXF_Int relax_error = 0;
   
  /*-----------------------------------------------------------------------------
   * Switch statement to direct control based on relax_type:
   *     relax_type = 0 -> Jacobi or CF-Jacobi
   *     relax_type = 1 -> Gauss-Seidel <--- very slow, sequential
   *     relax_type = 2 -> Gauss_Seidel: interior points in parallel ,
   *			 	   	  boundary sequential 
   *     relax_type = 3 -> hybrid: SOR-J mix off-processor, SOR on-processor
   *     		    with outer relaxation parameters (forward solve)
   *     relax_type = 4 -> hybrid: SOR-J mix off-processor, SOR on-processor
   *     		    with outer relaxation parameters (backward solve)
   *     relax_type = 5 -> hybrid: GS-J mix off-processor, chaotic GS on-node
   *     relax_type = 6 -> hybrid: SSOR-J mix off-processor, SSOR on-processor
   *     		    with outer relaxation parameters 
   *     relax_type = 7 -> Jacobi (uses Matvec), only needed in CGNR
   *     relax_type = 9 -> Direct Solve
   *-----------------------------------------------------------------------------*/
  
   if (relax_type == 0)
   {
      relax_error = jxf_hpPAMGRelax0(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 1)
   {
      relax_error = jxf_hpPAMGRelax1(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 2)
   {
      relax_error = jxf_hpPAMGRelax2(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 3)
   {
      relax_error = jxf_hpPAMGRelax3(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 4)
   {
      relax_error = jxf_hpPAMGRelax4(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 5)
   {
      relax_error = jxf_hpPAMGRelax5(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 6)
   {
      relax_error = jxf_hpPAMGRelax6(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 7)
   {
      relax_error = jxf_hpPAMGRelax7(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }
   else if (relax_type == 8)
   {
      relax_error = jxf_hpPAMGRelax8(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, l1_norms, par_app, Vtemp, Ztemp);   
   }
   else if (relax_type == 9)
   {
      relax_error = jxf_hpPAMGRelax9(hp_matrix, par_rhs, cf_marker, relax_points,
                                  relax_weight, omega, par_app, Vtemp);
   }

   else if (relax_type == 13)
   {
      relax_error = jxf_hpPAMGRelax13(hp_matrix, par_rhs, cf_marker, relax_points,
                                   relax_weight, omega, l1_norms, par_app, Vtemp, Ztemp);
   }
   else if (relax_type == 14)
   {
      relax_error = jxf_hpPAMGRelax14(hp_matrix, par_rhs, cf_marker, relax_points,
                                   relax_weight, omega, l1_norms, par_app, Vtemp, Ztemp);
   }
   else if (relax_type == 109)
   {
      relax_error = jxf_Mumps(hp_matrix, par_rhs, par_app);
   }
   return(relax_error); 
}
