//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_relax_7.c
 *  Date: 2011/09/03
 */ 

#include "jxf_hpcsr.h"
/*!
 * \fn JXF_Int jxf_PAMGRelax7
 * \brief Jacobi (uses Matvec), only needed in CGNR
 * \date 2011/09/03
 */
JXF_Int  
jxf_PAMGRelax7( jxf_ParCSRMatrix *par_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   MPI_Comm         comm         = jxf_ParCSRMatrixComm(par_matrix);
   jxf_CSRMatrix    *A_diag       = jxf_ParCSRMatrixDiag(par_matrix);
   JXF_Real          *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int             *A_diag_i     = jxf_CSRMatrixI(A_diag);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix);

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real         *Vtemp_data  = jxf_VectorData(Vtemp_local);

   JXF_Int             i;
   JXF_Int             relax_error = 0;
   JXF_Int             num_procs, my_id;

   JXF_Real          zero = 0.0;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  

  /*-------------------------------------------------------------------------
   * added by peghoty, if the comm_pkg of par_matrix is not created  
   * previously, something will be wrong when the "relax" function  
   * is called on multi-processors occasions.  2009/07/24
   *----------------------------------------------------------------------- */
   
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_matrix);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix); 
   }

      //------------------------------------------------------------------------------//
      //                       Jacobi (uses ParMatvec)                                //
      //------------------------------------------------------------------------------//

        /*-----------------------------------------------------------------
         * Copy par_rhs into temporary vector.
         *-----------------------------------------------------------------*/

         jxf_ParVectorCopy(par_rhs,Vtemp);

        /*-----------------------------------------------------------------
         * Perform Matvec Vtemp = par_rhs - A u
         *-----------------------------------------------------------------*/

         jxf_ParCSRMatrixMatvec(-1.0,par_matrix, par_app, 1.0, Vtemp);
         for (i = 0; i < n; i ++)
         {
           /*-----------------------------------------------------------
            * If diagonal is nonzero, relax point i; otherwise, skip it.
            *-----------------------------------------------------------*/

            if (A_diag_data[A_diag_i[i]] != zero)
            {
               u_data[i] += relax_weight * Vtemp_data[i] / A_diag_data[A_diag_i[i]];
            }
         }

   return(relax_error); 
}
JXF_Int  
jxf_hpPAMGRelax7( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelax7(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}
/*!
 * \fn JXF_Int jxf_PAMGRelaxAI7
 */
JXF_Int  
jxf_PAMGRelaxAI7( jxf_ParCSRMatrix *par_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp )
{
   MPI_Comm         comm         = jxf_ParCSRMatrixComm(par_matrix);
   jxf_CSRMatrix    *A_diag       = jxf_ParCSRMatrixDiag(par_matrix);
   JXF_Real          *A_diag_data  = jxf_CSRMatrixData(A_diag);
   JXF_Int             *A_diag_i     = jxf_CSRMatrixI(A_diag);
   jxf_ParCSRCommPkg    *comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix);

   JXF_Int             n        = jxf_CSRMatrixNumRows(A_diag);
   
   jxf_Vector      *u_local = jxf_ParVectorLocalVector(par_app);
   JXF_Real         *u_data  = jxf_VectorData(u_local);

   jxf_Vector      *Vtemp_local = jxf_ParVectorLocalVector(Vtemp);
   JXF_Real         *Vtemp_data  = jxf_VectorData(Vtemp_local);
   
   JXF_Int             i;
   JXF_Int             relax_error = 0;
   JXF_Int             num_procs, my_id;

   JXF_Real          zero = 0.0;

   jxf_MPI_Comm_size(comm, &num_procs);  
   jxf_MPI_Comm_rank(comm, &my_id);  

  /*-------------------------------------------------------------------------
   * added by peghoty, if the comm_pkg of par_matrix is not created  
   * previously, something will be wrong when the "relax" function  
   * is called on multi-processors occasions.  2009/07/24
   *----------------------------------------------------------------------- */
   
   if (!comm_pkg)
   {
      jxf_MatvecCommPkgCreate(par_matrix);
      comm_pkg = jxf_ParCSRMatrixCommPkg(par_matrix); 
   }

   
      //------------------------------------------------------------------------------//
      //                       Jacobi (uses ParMatvec)                                //
      //------------------------------------------------------------------------------//
      
        /*-----------------------------------------------------------------
         * Copy par_rhs into temporary vector.
         *-----------------------------------------------------------------*/

         jxf_ParVectorCopy(par_rhs,Vtemp);

        /*-----------------------------------------------------------------
         * Perform Matvec Vtemp = par_rhs - A u
         *-----------------------------------------------------------------*/

         jxf_ParCSRMatrixMatvec(-1.0, par_matrix, par_app, 1.0, Vtemp);
         for (i = 0; i < n; i ++)
         {
           /*-----------------------------------------------------------
            * If diagonal is nonzero, relax point i; otherwise, skip it.
            *-----------------------------------------------------------*/

            if (A_diag_data[A_diag_i[i]] != zero && relax_marker[i] == 1 )
            {
               u_data[i] += relax_weight * Vtemp_data[i] / A_diag_data[A_diag_i[i]];
            }
         }
         
   return(relax_error); 
}

JXF_Int  
jxf_hpPAMGRelaxAI7( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp )
{
   return jxf_PAMGRelaxAI7(jxf_hpCSRMatrixPar(hp_matrix), par_rhs, cf_marker, relax_marker, relax_points, relax_weight, omega, par_app, Vtemp);
}