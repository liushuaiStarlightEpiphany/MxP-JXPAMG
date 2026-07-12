//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  diagscale.c
 *  Date: 2013/10/27
 */ 
#ifndef JXF_HPCSRMV_HEADER 
#include "jxf_hpcsr.h"
#endif
#include "jxf_diagscale.h"


/*!
 * \fn JXF_Int JXF_DiagScaleSetup
 * \date 2013/10/27
 */
JXF_Int 
JXF_DiagScaleSetup( JXF_Solver        solver, 
                   JXF_hpCSRMatrix  hp_matrix )
{
   return( jxf_DiagScaleSetup( (void *) solver, jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) hp_matrix )) );
}
/*!
 * \fn JXF_Int jxf_DiagScaleSetup
 * \brief Setup phase of DiagScale preconditioner.
 * \date 2013/10/27
 */
JXF_Int
jxf_DiagScaleSetup( void             *solver, 
                   jxf_ParCSRMatrix  *par_matrix )
{
   return 0;
}

/*!
 * \fn JXF_Int JXF_DiagScalePrecond
 * \brief Solve phase of DiagScale preconditioner.
 * \date 2013/10/27
 */
JXF_Int 
JXF_DiagScalePrecond( JXF_Solver       solver,
                     JXF_hpCSRMatrix hp_matrix,
                     JXF_ParVector    par_rhs,
                     JXF_ParVector    par_app  )
{
   return( jxf_DiagScalePrecond( (void *) solver,
                              jxf_hpCSRMatrixPar((jxf_hpCSRMatrix *) hp_matrix),
                                (jxf_ParVector *) par_rhs,
                                (jxf_ParVector *) par_app ) );
}

/*!
 * \fn JXF_Int jxf_DiagScalePrecond
 * \brief DiagScale preconditioner.
 * \param solver pointer to NULL
 * \param par_matrix pointer to the coefficient matrix 
 * \param par_rhs pointer to the right hand side vector
 * \param par_app pointer to the approxamation vector
 * \date 2013/10/27
 */
JXF_Int
jxf_DiagScalePrecond( void            *solver,
                     jxf_ParCSRMatrix *par_matrix,
                     jxf_ParVector    *par_rhs,
                     jxf_ParVector    *par_app  )
{
   JXF_Int    *A_i     = jxf_CSRMatrixI(jxf_ParCSRMatrixDiag(par_matrix));
   JXF_Real *A_data  = jxf_CSRMatrixData(jxf_ParCSRMatrixDiag(par_matrix));
   JXF_Real *x_data  = jxf_VectorData(jxf_ParVectorLocalVector(par_app));
   JXF_Real *y_data  = jxf_VectorData(jxf_ParVectorLocalVector(par_rhs));
   JXF_Int  local_size = jxf_VectorSize(jxf_ParVectorLocalVector(par_app));
   JXF_Int  i, ierr = 0;
   
   for (i = 0; i < local_size; i ++)
   {
      x_data[i] = y_data[i] / A_data[A_i[i]];
   }
   
   return ierr;
}
