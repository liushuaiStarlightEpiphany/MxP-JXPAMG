//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_amg_setup.c
 *  Date: 2011/09/03
 */ 

#include "jx_pamg.h"
#include "jx_euclid.h"

/* forward declaration from spgemm_adapter.h */
int spgemm_adapter_multiply(jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C);
int spgemm_adapter_numeric_multiply(jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix *C);
jx_CSRMatrix *jx_CSRMatrixMultiplySymbolic(jx_CSRMatrix *A, jx_CSRMatrix *B);
#include "jx_ilu.h"

JX_Int jx_error_flag = 0;

JX_UInt jx_long_size_length_rap = 0;         /* Yue Xiaoqiang 2012/10/12 */
JX_UInt jx_long_size_length_interp = 0;      /* Yue Xiaoqiang 2012/10/12 */
JX_UInt jx_real_long_size_length_rap = 0;    /* Yue Xiaoqiang 2012/10/12 */
JX_UInt jx_real_long_size_length_interp = 0; /* Yue Xiaoqiang 2012/10/12 */

extern JX_Int *Pmarkers_global_rap;    /* Yue Xiaoqiang 2012/10/17 */
extern JX_Int *Pmarkers_global_interp; /* Yue Xiaoqiang 2012/10/17 */

/*!
 * \fn JX_Int jx_PAMGSetup
 * \brief Setup phase of PAMG solver.
 * \date 2011/09/03
 */

JX_Int
jx_hpPAMGSetup( void             *amg_vdata, 
              jx_hpCSRMatrix  *hp_matrix)
{
   MPI_Comm         comm     = jx_hpCSRMatrixComm(hp_matrix); 
   jx_ParAMGData   *amg_data = amg_vdata;
   /* Data Structure variables */
   jx_hpCSRMatrix **A_array;
   jx_ParVector    **F_array;
   jx_ParVector    **U_array;
   jx_ParVector     *Vtemp;
   jx_ParVector     *Ztemp = NULL;
   jx_ParCSRMatrix **P_array;
   jx_ParCSRMatrix **R_array;
   jx_ParVector     *Residual_array;
   JX_Real          **AI_measure_array;   
   JX_Int             **CF_marker_array;   
   JX_Int             **dof_func_array;   
   JX_Int              *dof_func;
   JX_Int              *col_offd_S_to_A;
   JX_Int              *col_offd_Sabs_to_A = NULL;
   JX_Int            *AIR_maxsize_ls;
   JX_Real            strong_threshold;
   JX_Real            AIR_strong_th;
   JX_Real            max_row_sum;
   JX_Real            trunc_factor;
   JX_Real            S_commpkg_switch;

   JX_Int      max_levels; 
   JX_Int      amg_logging;
   JX_Int      amg_print_level;
   JX_Int      debug_flag;
   JX_Int      local_num_vars;
   JX_Int      P_max_elmts;
   JX_Int      R_max_size; 
   
   JX_Solver *smoother = NULL;
   JX_Int        smooth_type = jx_ParAMGDataSmoothType(amg_data);
   JX_Int        smooth_num_levels = jx_ParAMGDataSmoothNumLevels(amg_data);
   char      *euclidfile;
   JX_Int	      eu_level;
   JX_Int	      eu_bj;
   JX_Real     eu_sparse_A;

   /* Local variables */
   JX_Real           *AI_measure;
   JX_Real           *AI_measure2;
   JX_Int              *CF_marker;
   JX_Int              *CFN_marker = NULL;
   jx_ParCSRMatrix  *par_S;
   jx_ParCSRMatrix  *Sabs = NULL;
   jx_ParCSRMatrix  *par_S2;
   jx_ParCSRMatrix  *par_P;
   jx_ParCSRMatrix  *R;
   jx_hpCSRMatrix  *A_H;

   JX_Int       wall_time_option;         /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_coarsen = 0.0;  /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_rap = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    wall_time_interp = 0.0;   /* newly added Yue Xiaoqiang 2015/09/30 */
   JX_Real    tmp_wall_time = 0.0;      /* newly added Yue Xiaoqiang 2015/09/30 */

   JX_Int       old_num_levels, num_levels;
   JX_Int       level;
   JX_Int       local_size, i;
   JX_Int       first_local_row;
   JX_Int       coarse_size_last = -1;
   JX_Int       coarse_size;
   JX_Int       coarsen_type;
   JX_Int       interp_type;
   JX_Int       restri_type;
   JX_Int       measure_type;
   JX_Int       setup_type;
   JX_Int       fine_size;
   JX_Int       rest, tms, indx;
   JX_Int       not_finished_coarsening = 1;
   JX_Int       Setup_err_flag = 0;
   JX_Int       coarse_threshold; /* we don't fix it to be 9.  peghoty 2010/04/14 */
   JX_Int       coarse_threshold_2; /* by xu: 2013/11/25 */
   JX_Int       j, k;
   JX_Int       num_procs, my_id, num_threads;
   JX_Int      *grid_relax_type = jx_ParAMGDataGridRelaxType(amg_data);
   JX_Int       relax_order = jx_ParAMGDataRelaxOrder(amg_data);
   JX_Int       num_functions = jx_ParAMGDataNumFunctions(amg_data);
   JX_Int       spmt_rap_type = jx_ParAMGDataSpMtRapType(amg_data);
   JX_Int       ai_measure_type = jx_ParAMGDataAIMeasureType(amg_data);
   JX_Int       num_paths = jx_ParAMGDataNumPaths(amg_data);
   JX_Int       agg_num_levels = jx_ParAMGDataAggNumLevels(amg_data);
   JX_Int       agg_interp_type = jx_ParAMGDataAggInterpType(amg_data);
   JX_Int       agg_P_max_elmts = jx_ParAMGDataAggPMaxElmts(amg_data);
   //JX_Int       agg_P12_max_elmts = jx_ParAMGDataAggP12MaxElmts(amg_data);
   JX_Real    agg_trunc_factor = jx_ParAMGDataAggTruncFactor(amg_data);
   //JX_Real    agg_P12_trunc_factor = jx_ParAMGDataAggP12TruncFactor(amg_data);
   JX_Int       rap2 = jx_ParAMGDataRAP2(amg_data);
   JX_Int       keepTranspose = jx_ParAMGDataKeepTranspose(amg_data);
   JX_Int       print_coarse_matrix = jx_ParAMGDataPrintCoarseSystem(amg_data);
   JX_Int      *opt_icor;
   JX_Int      *coarse_dof_func;
   JX_Int      *coarse_pnts_global;
   JX_Int	    *coarse_pnts_global1;
   JX_Real    size;
   JX_Real    coarse_ratio;           /* newly added peghoty 2010/04/14 */
   JX_Real    wall_time = 0.0;        /* for debugging instrumentation */
   JX_Int       measure_type_rlx;       // newly added peghoty 2010/05/29
   JX_Int       number_syn_rlx;         // newly added peghoty 2010/05/29
   JX_Real    measure_threshold_rlx;  // newly added peghoty 2010/05/29
   JX_Real         **l1_norms = NULL;

   /* coarse matrices */
   char FileNameCoaMat[256];

   /* ai statistic information*/
   JX_Int       num_vars_local = 0, num_vars_global;
   JX_Int       num_ai_local = 0, num_ai_global;
   JX_Int       num_ai_local_valid = 0, num_ai_global_valid;
   JX_Int       num_ai_c_local = 0, num_ai_c_global;
   JX_Real    mai_local = 0.0, mai_global;
   JX_Real    mai_local_valid = 0.0, mai_global_valid;
   JX_Real    mai_c_local = 0.0, mai_c_global;

   JX_Int       num_vars_local_0 = 0;
   JX_Int       num_ai_local_0 = 0;
   JX_Int       num_ai_local_valid_0 = 0;
   JX_Int       num_ai_c_local_0 = 0;
   JX_Real    mai_local_0 = 0.0;
   JX_Real    mai_local_valid_0 = 0.0;
   JX_Real    mai_c_local_0 = 0.0;
   JX_Real    mai_threshold = 0.1;

   jx_MPI_Comm_size(comm, &num_procs);   
   jx_MPI_Comm_rank(comm, &my_id);
   if ((num_procs > 1) && (spmt_rap_type != 1) && (spmt_rap_type != 3) && (spmt_rap_type != 7))  /* Yue Xiaoqiang 2012/10/12 */
   {
      spmt_rap_type = 1;
   }
   num_threads = jx_NumThreads();               /* Yue Xiaoqiang 2012/10/12 */
   opt_icor = jx_CTAlloc(JX_Int, 5*num_threads+2); /* Yue Xiaoqiang 2012/10/12 */

   wall_time_option = jx_ParAMGDataWallTimeOption(amg_data);
   old_num_levels   = jx_ParAMGDataNumLevels(amg_data);
   max_levels       = jx_ParAMGDataMaxLevels(amg_data);
   amg_logging      = jx_ParAMGDataLogging(amg_data);
   amg_print_level  = jx_ParAMGDataPrintLevel(amg_data);
   coarsen_type     = jx_ParAMGDataCoarsenType(amg_data);
   measure_type     = jx_ParAMGDataMeasureType(amg_data);
   setup_type       = jx_ParAMGDataSetupType(amg_data);
   debug_flag       = jx_ParAMGDataDebugFlag(amg_data);
   dof_func         = jx_ParAMGDataDofFunc(amg_data);
   interp_type      = jx_ParAMGDataInterpType(amg_data);
   restri_type      = jx_ParAMGDataRestriction(amg_data);
   AIR_strong_th    = jx_ParAMGDataAIRStrongTh(amg_data);
   euclidfile       = jx_ParAMGDataEuclidFile(amg_data);
   eu_level         = jx_ParAMGDataEuLevel(amg_data);
   eu_bj            = jx_ParAMGDataEuBJ(amg_data);
   eu_sparse_A      = jx_ParAMGDataEuSparseA(amg_data);
   coarse_threshold = jx_ParAMGDataCoarseThreshold(amg_data);           /* newly added peghoty 2010/04/14 */
   coarse_ratio     = jx_ParAMGDataCoarseRatio(amg_data);               /* newly added peghoty 2010/04/14 */
   measure_type_rlx = jx_ParAMGDataMeasureTypeRlx(amg_data);            /* newly added peghoty 2010/05/29 */
   number_syn_rlx   = jx_ParAMGDataNumberSynRlx(amg_data);              /* newly added peghoty 2010/05/29 */
   measure_threshold_rlx = jx_ParAMGDataMeasureThresholdRlx(amg_data);  /* newly added peghoty 2010/05/29 */

   coarse_threshold_2 = -1;
   
   jx_ParCSRMatrixSetNumNonzeros(jx_hpCSRMatrixPar(hp_matrix));
   jx_ParCSRMatrixSetDNumNonzeros(jx_hpCSRMatrixPar(hp_matrix));
   jx_ParAMGDataNumVariables(amg_data) = jx_hpCSRMatrixNumRows(hp_matrix);
   
   if (setup_type == 0) 
   {
      return(Setup_err_flag);
   }

   par_S = NULL;
   A_H = NULL;

   A_array = jx_hpAMGDataAArray(amg_data);
   P_array = jx_ParAMGDataPArray(amg_data);
   R_array = jx_ParAMGDataRArray(amg_data);
   AIR_maxsize_ls = jx_ParAMGDataAIRMaxSizeLS(amg_data);
   AI_measure_array = jx_ParAMGDataAIMeasureArray(amg_data);
   CF_marker_array = jx_ParAMGDataCFMarkerArray(amg_data);
   dof_func_array  = jx_ParAMGDataDofFuncArray(amg_data);
   local_size = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(hp_matrix));

   grid_relax_type[3] = jx_ParAMGDataUserCoarseRelaxType(amg_data); 

   if (A_array || P_array || R_array || AIR_maxsize_ls || AI_measure_array || CF_marker_array || dof_func_array)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (A_array[j])
         {
            jx_hpCSRMatrixDestroy(A_array[j]);
            A_array[j] = NULL;
         }
       
         if (dof_func_array[j])
         {
            jx_TFree(dof_func_array[j]);
            dof_func_array[j] = NULL;
         }
      }

      for (j = 0; j < old_num_levels - 1; j ++)
      {
         if (P_array[j])
         {
            jx_ParCSRMatrixDestroy(P_array[j]);
            P_array[j] = NULL;
         }
         if (R_array[j])
         {
            jx_ParCSRMatrixDestroy(R_array[j]);
            R_array[j] = NULL;
         }
      }

      jx_TFree(AIR_maxsize_ls);
      AIR_maxsize_ls = NULL;

     /*-------------------------------------------------------------------
      *  Special case use of CF_marker_array when old_num_levels = 1
      *  requires us to attempt this deallocation every time
      *------------------------------------------------------------------*/
      
      if (CF_marker_array[0])
      {
         jx_TFree(CF_marker_array[0]);
         CF_marker_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (CF_marker_array[j])
         {
            jx_TFree(CF_marker_array[j]);
            CF_marker_array[j] = NULL;
         }
      }
      
      /* for AI_measure: added by xwxu */
      if (AI_measure_array[0])
      {
         jx_TFree(AI_measure_array[0]);
         AI_measure_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels-1; j ++)
      {
         if (AI_measure_array[j])
         {
            jx_TFree(AI_measure_array[j]);
            AI_measure_array[j] = NULL;
         }
      }
   }

   if (A_array == NULL)
   {
      A_array = jx_CTAlloc(jx_hpCSRMatrix*, max_levels);
   }

   if (P_array == NULL && max_levels > 1)
   {
      P_array = jx_CTAlloc(jx_ParCSRMatrix*, max_levels - 1);
   }

   /* If retri_type != 0, R != P^T, allocate R matrices */
   if (restri_type)
   {
      if (R_array == NULL && max_levels > 1)
      {
         R_array = jx_CTAlloc(jx_ParCSRMatrix*, max_levels-1);
         AIR_maxsize_ls = jx_CTAlloc(JX_Int, max_levels-1);
      }
   }

   if (AI_measure_array == NULL)
   {
      AI_measure_array = jx_CTAlloc(JX_Real*, max_levels);
   }
   if (CF_marker_array == NULL)
   {
      CF_marker_array = jx_CTAlloc(JX_Int*, max_levels);
   }
   if (dof_func_array == NULL)
   {
      dof_func_array = jx_CTAlloc(JX_Int*, max_levels);
   }

   if (num_functions > 1 && dof_func == NULL)
   {
      first_local_row = jx_hpCSRMatrixFirstRowIndex(hp_matrix);
      dof_func = jx_CTAlloc(JX_Int, local_size);
      rest = first_local_row-((first_local_row / num_functions)*num_functions);
      indx = num_functions - rest;
      if (rest == 0) 
      {
         indx = 0;
      }
      k = num_functions - 1;
      for (j = indx - 1; j > -1; j --)
      {
         dof_func[j] = k --;
      }
      tms = local_size / num_functions;
      if (tms*num_functions + indx > local_size) 
      {
         tms --;
      }
      for (j = 0; j < tms; j ++)
      {
         for (k = 0; k < num_functions; k ++)
         {
            dof_func[indx++] = k;
         }
      }
      k = 0;
      while (indx < local_size)
      {
         dof_func[indx++] = k ++;
      }
      jx_ParAMGDataDofFunc(amg_data) = dof_func;
   }

   A_array[0] = hp_matrix;

   dof_func_array[0] = dof_func;
   jx_ParAMGDataAIMeasureArray(amg_data) = AI_measure_array;
   jx_ParAMGDataCFMarkerArray(amg_data) = CF_marker_array;
   jx_ParAMGDataDofFuncArray(amg_data) = dof_func_array;
   jx_hpAMGDataAArray(amg_data) = A_array;
   jx_ParAMGDataPArray(amg_data) = P_array;
   /* If R != P^T */
   if (restri_type)
   {
      jx_ParAMGDataRArray(amg_data) = R_array;
      jx_ParAMGDataAIRMaxSizeLS(amg_data) = AIR_maxsize_ls;
   }
   else
   {
      jx_ParAMGDataRArray(amg_data) = P_array;
   }

   Vtemp = jx_ParAMGDataVtemp(amg_data);

   if (Vtemp != NULL)
   {
      jx_ParVectorDestroy(Vtemp);
      Vtemp = NULL;
   }

   Vtemp = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[0]),
                              jx_hpCSRMatrixGlobalNumRows(A_array[0]),
                              jx_hpCSRMatrixRowStarts(A_array[0]));
   jx_ParVectorInitialize(Vtemp);
   jx_ParVectorSetPartitioningOwner(Vtemp,0);
   jx_ParAMGDataVtemp(amg_data) = Vtemp;

   if (num_threads > 1)
   {
      Ztemp = jx_ParAMGDataZtemp(amg_data);

      if (Ztemp != NULL)
      {
         jx_ParVectorDestroy(Ztemp);
         Ztemp = NULL;
      }

      for (j = 1; j < 4; j ++)
      {
         if (grid_relax_type[j] == 3 || grid_relax_type[j] == 4 || grid_relax_type[j] == 6 ||
             grid_relax_type[j] == 8 || grid_relax_type[j] == 13 || grid_relax_type[j] == 14)
         {
            Ztemp = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[0]), jx_hpCSRMatrixGlobalNumRows(A_array[0]), jx_hpCSRMatrixRowStarts(A_array[0]));
            jx_ParVectorInitialize(Ztemp);
            jx_ParVectorSetPartitioningOwner(Ztemp, 0);
            jx_ParAMGDataZtemp(amg_data) = Ztemp;
            break;
         }
      }
   }

   F_array = jx_ParAMGDataFArray(amg_data);
   U_array = jx_ParAMGDataUArray(amg_data);

   if (F_array != NULL || U_array != NULL)
   {
      for (j = 1; j < old_num_levels; j ++)
      {
         if (F_array[j] != NULL)
         {
            jx_ParVectorDestroy(F_array[j]);
            F_array[j] = NULL;
         }
         if (U_array[j] != NULL)
         {
            jx_ParVectorDestroy(U_array[j]);
            U_array[j] = NULL;
         }
      }
   }

   if (F_array == NULL)
   {
      F_array = jx_CTAlloc(jx_ParVector*, max_levels);
   }
   if (U_array == NULL)
   {
      U_array = jx_CTAlloc(jx_ParVector*, max_levels);
   }

  /*---------------------------------------------------------
   *  They have been moved to "amg_solve" since they are   
   *  not used here.  peghoty 2009/07/27 
   *--------------------------------------------------------*/
   //F_array[0] = f;
   //U_array[0] = u;

   jx_ParAMGDataFArray(amg_data) = F_array;
   jx_ParAMGDataUArray(amg_data) = U_array;

  /*----------------------------------------------------------
   *   Initialize jx_ParAMGData
   *---------------------------------------------------------*/
   not_finished_coarsening = 1;
   level = 0;
   strong_threshold = jx_ParAMGDataStrongThreshold(amg_data);
   max_row_sum = jx_ParAMGDataMaxRowSum(amg_data);
   trunc_factor = jx_ParAMGDataTruncFactor(amg_data);
   P_max_elmts = jx_ParAMGDataPMaxElmts(amg_data);
   S_commpkg_switch = jx_ParAMGDataSCommPkgSwitch(amg_data);
   if (smooth_num_levels > level)
   {
      smoother = jx_CTAlloc(JX_Solver, smooth_num_levels);
      jx_ParAMGDataSmoother(amg_data) = smoother;
   }

   /* if AI-based coarsening is used, set ai_measure_type = 1. */
   if (coarsen_type == 990 || coarsen_type == 991 || coarsen_type == 993 || 
       coarsen_type == 96 || coarsen_type == 98 || coarsen_type == 910 || 
       coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || 
       coarsen_type == 938 || coarsen_type == 968)
   {
      ai_measure_type = 1; 

   }

  /*-----------------------------------------------------
   *   Enter Coarsening Loop
   *-----------------------------------------------------*/

   if (amg_print_level > 0 && my_id ==0 && ai_measure_type)
   {
      jx_printf(" \n");
      jx_printf("================================================================== \n");
      jx_printf("+++++++++++++ multi-scale/AI infomation for levels +++++++++++++++ \n");
   }
   while (not_finished_coarsening)
   {  
      fine_size = jx_hpCSRMatrixGlobalNumRows(A_array[level]);

      if (level > 0)
      {   
         F_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                             jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jx_hpCSRMatrixRowStarts(A_array[level]));
         jx_ParVectorInitialize(F_array[level]);
         jx_ParVectorSetPartitioningOwner(F_array[level],0);
            
         U_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                             jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jx_hpCSRMatrixRowStarts(A_array[level]));
         jx_ParVectorInitialize(U_array[level]);
         jx_ParVectorSetPartitioningOwner(U_array[level],0);
      }


     /*--------------------------------------------------------------
      *  Select coarse-grid points on 'level' : returns CF_marker
      *  for the level.  Returns strength matrix, par_S 
      *  Returns AI_Measure. 
      *--------------------------------------------------------------*/
      
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();
      if (debug_flag == 3)
      {
         jx_printf("\n ===== Proc = %d     Level = %d  =====\n",my_id, level);
         fflush(NULL);
      }

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      if (max_levels > 1)
      {
         local_num_vars = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level]));

        /*--------------------------------------------------------------------------
         *  Get the Strength Matrix 
         *-------------------------------------------------------------------------*/ 
         jx_PAMGCreateS(jx_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);      
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]),par_S,&col_offd_S_to_A);
         }

         /* for AIR, need absolute value SOC */
         if ((restri_type == 1) || (restri_type == 2))
         {
            jx_PAMGCreateSabs(jx_hpCSRMatrixPar(A_array[level]),AIR_strong_th,1.0,num_functions,dof_func_array[level],&Sabs);
            
            
            col_offd_Sabs_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]), Sabs, &col_offd_Sabs_to_A);
            }
         }
         
         if (ai_measure_type == 1)
         {
            
           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/   
            jx_PAMGMeasureAI(par_S, jx_hpCSRMatrixPar(A_array[level]),debug_flag, &AI_measure);
           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

        /*--------------------------------------------------------------------------
         *    Do the appropriate coarsening as follows:
         *
         *      coarsen_type =  0: CLJP 
         *      coarsen_type =  1: Ruge
         *      coarsen_type = 11: Ruge 1st pass only
         *      coarsen_type =  2: Ruge2B
         *      coarsen_type =  3: Ruge3
         *      coarsen_type =  4: Ruge3c
         *      coarsen_type =  5: Ruge relax special points
         *      coarsen_type =  6: Falgout
         *      coarsen_type =  8: PMIS
         *      coarsen_type = 10: HMIS
         *      coarsen_type = 90: RCLJP 
         *      coarsen_type = 91: RRS0
         *      coarsen_type = 990: CLJP_AI 
         *      coarsen_type = 991: Ruge_AI
         *      coarsen_type = 993: Ruge3_AI
         *      coarsen_type = 96: Falgout_AI
         *      coarsen_type = 98: PMIS_AI
         *      coarsen_type = 910: HMIS_AI
         *      coarsen_type = 908, 918, 928, 938, 968: AI-TYPE.
         *                                          peghoty 2010/05/29
         *-------------------------------------------------------------------------*/ 

         //jx_printf("=========== coarsen_type = %d \n", coarsen_type); 
         //jx_printf("myid is %d,coarsen_type is %d,before Coarsen run success\n",my_id,coarsen_type);
         if (coarsen_type == 6)  
         {
            jx_PAMGCoarsenFalgout(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 96) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenFalgoutAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 8) 
         {
            jx_PAMGCoarsenPMIS(par_S, jx_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 10) 
         {
            jx_PAMGCoarsenHMIS(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }      
         else if ((coarsen_type == 910) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenHMISAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, debug_flag, &CF_marker);
         }      
         else if (coarsen_type == 1 || coarsen_type == 2 ||
                  coarsen_type == 3 || coarsen_type == 4 ||
                  coarsen_type == 5 || coarsen_type == 11  ) 
         {
            jx_PAMGCoarsenRuge(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, coarsen_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 991) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 11, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 992) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 1, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 993) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenRugeAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 3, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 0)    
         {
            jx_PAMGCoarsen(par_S, jx_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 990) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 90)
         {
            jx_PAMGCoarsenRCLJP(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type_rlx, number_syn_rlx,
                                measure_threshold_rlx, 0, debug_flag, &CF_marker);
         }         
         else if (coarsen_type == 91)
         {
            jx_PAMGCoarsenRRS0(par_S, jx_hpCSRMatrixPar(A_array[level]), measure_type, measure_type_rlx,
                               number_syn_rlx, measure_threshold_rlx, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 98) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenPMISAI(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || coarsen_type == 938 || 
                  coarsen_type == 968) && (ai_measure_type == 1))
         {
            jx_PAMGCoarsenXML(par_S, jx_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, coarsen_type, debug_flag, &CF_marker);
         }

         if (level < agg_num_levels)
         {  
            jx_PAMGCoarseParms(comm, local_num_vars, 1,
                               dof_func_array[level], CF_marker,
                               &coarse_dof_func, &coarse_pnts_global1);
            jx_PAMGCreate2ndS(par_S, CF_marker, num_paths, coarse_pnts_global1, &par_S2);
           /*--------------------------------------------------------------------------
            *  Get the AI Measure
            *-------------------------------------------------------------------------*/
            if (ai_measure_type == 1)
            {
               jx_PAMGMeasureAI(par_S2, par_S2, debug_flag, &AI_measure2);
            }

            if (coarsen_type == 10)
            {
               jx_PAMGCoarsenHMIS(par_S2, par_S2, measure_type+3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 910) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenHMISAI(par_S2, par_S2, AI_measure2, 0, measure_type+3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 8)
            {
               jx_PAMGCoarsenPMIS(par_S2, par_S2, 3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 98) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenPMISAI(par_S2, par_S2, AI_measure2, 3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 6)
            {
               jx_PAMGCoarsenFalgout(par_S2, par_S2, measure_type, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 96) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenFalgoutAI(par_S2, par_S2, AI_measure2, measure_type, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 0)
            {
               jx_PAMGCoarsen(par_S2, par_S2, 0, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 990) && (ai_measure_type == 1))
            {
               jx_PAMGCoarsenAI(par_S2, par_S2, AI_measure2, 0, debug_flag, &CFN_marker);
            }
            jx_ParCSRMatrixDestroy(par_S2);
            if (ai_measure_type == 1)
            {
               jx_TFree(AI_measure2);
            }
         }

         /* Here for changes of min_coarse_size */

         if (level < agg_num_levels)
         {
            if (agg_interp_type == 4)
            {
               jx_PAMGCorrectCFMarker(CF_marker, local_num_vars, CFN_marker);
               jx_TFree(coarse_pnts_global1);
               jx_TFree(CFN_marker);
            }
         }

        /*--------------------------------------------
         *  store the CF array
         *-------------------------------------------*/ 
         CF_marker_array[level] = CF_marker;

         if (debug_flag == 1)
         {
            wall_time = jx_time_getWallclockSeconds() - wall_time;
            jx_printf("Proc = %d    Level = %d    Coarsen Time = %f\n",my_id,level, wall_time); 
            fflush(NULL);
         }

         if (wall_time_option == 1)
         {
            wall_time_coarsen += (jx_time_getWallclockSeconds() - tmp_wall_time);
         }

        /*-------------------------------------------------------------------------
         *  Get the coarse parameters
         *-------------------------------------------------------------------------*/ 
         jx_PAMGCoarseParms(comm, local_num_vars, num_functions, 
                            dof_func_array[level], CF_marker,
                            &coarse_dof_func, &coarse_pnts_global);

         dof_func_array[level+1] = NULL;
         if (num_functions > 1) dof_func_array[level+1] = coarse_dof_func;

#ifdef JX_NO_GLOBAL_PARTITION
         if (my_id == (num_procs -1)) 
         {
            coarse_size = coarse_pnts_global;
         }
         jx_MPI_Bcast(&coarse_size, 1, JX_MPI_INT, num_procs-1, comm);
#else
         coarse_size = coarse_pnts_global[num_procs];
#endif
      
        if ( coarse_size <= coarse_threshold && coarse_size_last >= coarse_size ) {
           break;  
        } 

      }
      else  /* max_levels = 1 */
      {

        /*--------------------------------------------------------------------------
         *  Get the Strength Matrix 
         *-------------------------------------------------------------------------*/        
         jx_PAMGCreateS(jx_hpCSRMatrixPar(A_array[level]),strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jx_PAMGCreateSCommPkg(jx_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }

         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jx_PAMGMeasureAI(par_S, jx_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }
         par_S = NULL;
         coarse_pnts_global = NULL;
         CF_marker = jx_CTAlloc(JX_Int, local_size );
         for (i = 0; i < local_size ; i ++) 
         {
            CF_marker[i] = 1;
         }
         CF_marker_array[level] = CF_marker;
         coarse_size = fine_size;

         local_num_vars = jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level]));
      }
      
     /*-------------------------------------------------------------------------
      *  compute ai-information for all-pnts and C-pnts!
      *-------------------------------------------------------------------------*/ 
     
      if ((max_levels == 1 || level != max_levels-1) && (ai_measure_type == 1))
      { 
         for (i = 0; i < local_num_vars; i++) {
             num_vars_local++ ;
             //if (AI_measure[i] > 0.0 && AI_measure[i] < 1.0)
             if (AI_measure[i] > 1.0e-10) 
             {
                mai_local += AI_measure[i];
                num_ai_local++ ;
                if (AI_measure[i] >= mai_threshold) {
                   mai_local_valid += AI_measure[i];
                   num_ai_local_valid++ ;
                   if (CF_marker[i] == 1) {
                      //if (level == 0 ) jx_printf("i_c = %d, ai_measure= %f \n", i, AI_measure[i]);
                      mai_c_local += AI_measure[i];
                      num_ai_c_local++;
                   }
                }
             }
         }
      }

      if (level == 0) {
        num_ai_local_0 = num_ai_local;
        num_ai_local_valid_0 = num_ai_local_valid;
        num_vars_local_0 = num_vars_local;
        num_ai_c_local_0 = num_ai_c_local;
        mai_local_0 = mai_local;
        mai_local_valid_0 = mai_local_valid;
        mai_c_local_0 = mai_c_local;
      }


     /*-------------------------------------------------------------------------
      *  if no coarse-grid, stop coarsening, and set the
      *  coarsest solve to be a single sweep of Jacobi !
      *-------------------------------------------------------------------------*/ 

      if ( (coarse_size == 0) || (coarse_size == fine_size) )
      {
         JX_Int  *num_grid_sweeps   = jx_ParAMGDataNumGridSweeps(amg_data);
         JX_Int **grid_relax_points = jx_ParAMGDataGridRelaxPoints(amg_data);
         if (grid_relax_type[3] == 9)
	 {
	    grid_relax_type[3] = grid_relax_type[0];
	    num_grid_sweeps[3] = 1;
	    if (grid_relax_points)
	    {
	       grid_relax_points[3][0] = 0; 
	    }
	 }
         if (par_S) 
         {
            jx_ParCSRMatrixDestroy(par_S);
         }
         jx_TFree(coarse_pnts_global);
         if (level > 0)
         {
           /*-------------------------------------------------------------
            * Note special case treatment of CF_marker is necessary
            * to do CF relaxation correctly when num_levels = 1
            *------------------------------------------------------------*/ 
            jx_TFree(CF_marker_array[level]);
            jx_ParVectorDestroy(F_array[level]);
            jx_ParVectorDestroy(U_array[level]);
         }
         break; 
      }

      /* Build restriction */
      if (restri_type)
      {
         /* !!! Ensure that CF_marker contains -1 or 1 !!! */
         for (i = 0; i < jx_CSRMatrixNumRows(jx_hpCSRMatrixDiag(A_array[level])); i ++)
         {
            CF_marker[i] = CF_marker[i] > 0 ? 1 : -1;
         }
         if (restri_type == 1) /* distance-1 AIR */
         {
            jx_PAMGBuildRestrAIR(jx_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                 Sabs, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], 
                                 debug_flag, trunc_factor, P_max_elmts, 
                                 col_offd_Sabs_to_A, &R, &R_max_size );
         }
         else if (restri_type == 2) /* distance-2 AIR */
         {
            jx_PAMGBuildRestrDist2AIR(jx_hpCSRMatrixPar(A_array[level]), CF_marker, 
                                      Sabs, coarse_pnts_global, num_functions, 
                                      dof_func_array[level], 
                                      debug_flag, trunc_factor, P_max_elmts, 
                                      col_offd_Sabs_to_A, &R, &R_max_size );
         }
         else /* R=[0 I_c] */
         {
            jx_PAMGBuildOPntRestr(jx_hpCSRMatrixPar(A_array[level]), CF_marker, coarse_pnts_global, &R, &R_max_size );
         }
         if (Sabs)
         {
            jx_ParCSRMatrixDestroy(Sabs);
            Sabs = NULL;
         }
         jx_TFree(col_offd_Sabs_to_A);
      }


     /*-----------------------------------------------------------------------------------
      *   Build prolongation matrix, P, and place in P_array[level]
      *
      *     interp_type =  0: modified classical interpolation
      *     interp_type =  3: direct interpolation (with separation of weights)
      *     interp_type =  4: multipass interpolation
      *     interp_type =  5: multipass interpolation (with separation of weights)
      *     interp_type =  6: extended classical modified interpolation
      *     interp_type =  7: extended (if no common C neighbor) classical 
      *                       modified interpolation
      *     interp_type =  8: standard interpolation
      *     interp_type =  9: standard interpolation (with separation of weights)
      *-----------------------------------------------------------------------------------*/
      
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      if (level < agg_num_levels)
      {
         if (agg_interp_type == 4)
         {
            jx_PAMGBuildMultipass(jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S, coarse_pnts_global,
                                  num_functions, dof_func_array[level], debug_flag,
                                  agg_trunc_factor, agg_P_max_elmts, 0, col_offd_S_to_A, &par_P);
         }
      }
      else
      {  
         if (interp_type == 0)  // classical modified interpolation
         {
            if (spmt_rap_type == 1 || spmt_rap_type == 3 || spmt_rap_type == 7)
            {
               jx_PAMGBuildInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, trunc_factor,
                                P_max_elmts, col_offd_S_to_A, &par_P ); 
            }
            else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
            {
               jx_PAMGBuildInterp1( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                 coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag, trunc_factor,
                                 P_max_elmts, col_offd_S_to_A, &par_P, opt_icor );
            }
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 3)  // direct interpolation (with separation of weights)
         {
            jx_PAMGBuildDirInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag, 
                                trunc_factor, P_max_elmts, 
                                col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         } 
         else if (interp_type == 4)  // multipass interpolation
         {
            jx_PAMGBuildMultipass( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 5) // multipass interpolation (with separation of weights)
         {
            jx_PAMGBuildMultipass( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                par_S, coarse_pnts_global, num_functions, 
                                dof_func_array[level], debug_flag,
                                trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 6)  // extended classical modified interpolation
         {
             jx_PAMGBuildExtPIInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                   par_S, coarse_pnts_global, num_functions, 
                                   dof_func_array[level], debug_flag,
                                   trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
             jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 7) // extended (if no common C neighbor) classical modified interpolation
         {
            jx_PAMGBuildExtPICCInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                     par_S, coarse_pnts_global, num_functions, 
                                     dof_func_array[level], debug_flag,
                                     trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 8) // standard interpolation
         {
            jx_PAMGBuildStdInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 9) // standard interpolation (with separation of weights)
         {
            jx_PAMGBuildStdInterp( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], 
                                 par_S, coarse_pnts_global, num_functions, 
                                 dof_func_array[level], debug_flag,
                                 trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 100) /* 1pt interpolation */
         {
            jx_PAMGBuildInterpOnePnt( jx_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                    par_S, coarse_pnts_global, num_functions,
                                    dof_func_array[level], debug_flag, col_offd_S_to_A, &par_P );
            jx_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 20)
         {
            if (restri_type)
            {
               jx_PAMGBuildInterpDinvRATrans( jx_hpCSRMatrixPar(A_array[level]), R, &par_P );
            }
            jx_TFree(col_offd_S_to_A);
         }
      }

      if (debug_flag == 1)
      {
         wall_time = jx_time_getWallclockSeconds() - wall_time;
         jx_printf("Proc = %d    Level = %d    Build Interp Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_interp += (jx_time_getWallclockSeconds() - tmp_wall_time);
      }

      P_array[level] = par_P;
      
      if (restri_type)
      {
         R_array[level] = R;
         AIR_maxsize_ls[level] = R_max_size;
      }

      if (par_S) 
      {
         jx_ParCSRMatrixDestroy(par_S);
      }
      par_S = NULL;

     /*--------------------------------------------------------------------------------
      *   Build coarse-grid operator, A_array[level+1] by R*A*P
      *-------------------------------------------------------------------------------*/
       
      if (debug_flag == 1) wall_time = jx_time_getWallclockSeconds();

      if (wall_time_option == 1) tmp_wall_time = jx_time_getWallclockSeconds();

      A_H = jx_hpInithpCSRMatrix();

      if (restri_type)
      {  
         /* Use two matrix products to generate A_H */
         jx_ParCSRMatrix *Q = jx_ParMatmul(jx_hpCSRMatrixPar(A_array[level]), P_array[level]);
         jx_hpCSRMatrixPar(A_H) = jx_ParMatmul(R_array[level], Q);
         jx_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jx_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jx_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         jx_ParCSRMatrixOwnsRowStarts(R_array[level]) = 0;
         if (num_procs > 1) jx_MatvecCommPkgCreate(jx_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jx_ParCSRMatrixDestroy(Q);
      }
      else if (rap2)
      {
         /* Use two matrix products to generate A_H */
         jx_ParCSRMatrix *Q = jx_ParMatmul(jx_hpCSRMatrixPar(A_array[level]), P_array[level]);
         jx_hpCSRMatrixPar(A_H) = jx_ParTMatmul(P_array[level], Q);
         jx_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jx_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jx_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         if (num_procs > 1) jx_MatvecCommPkgCreate(jx_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jx_ParCSRMatrixDestroy(Q);
      }
      else if (spmt_rap_type == 1)
      {
         jx_PAMGBuildCoarseOperatorKT( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &(jx_hpCSRMatrixPar(A_H)) );
      }
      else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
      {
         jx_PAMGBuildCoarseOperatorOMP( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], &(jx_hpCSRMatrixPar(A_H)), opt_icor );
      }
      else if (spmt_rap_type == 3) /* Liushuai 2026/07/12: BeidouBLAS SpGEMM RAP */
      {
         if (num_procs == 1)
         {
            jx_ParCSRMatrix *P_par = P_array[level];
            jx_ParCSRMatrix *A_par = jx_hpCSRMatrixPar(A_array[level]);
            jx_CSRMatrix *P_csr = jx_ParCSRMatrixDiag(P_par);
            jx_CSRMatrix *A_csr = jx_ParCSRMatrixDiag(A_par);
            jx_CSRMatrix *R_csr = NULL;
            jx_CSRMatrix *T_csr = NULL;
            jx_CSRMatrix *AH_csr = NULL;

            jx_CSRMatrixTranspose(P_csr, &R_csr, 1);

            if (spgemm_adapter_multiply(R_csr, A_csr, &T_csr) != 0)
            {
               jx_printf("Proc %d: T3 T=R*A failed at level %d\n", my_id, level);
               jx_CSRMatrixDestroy(R_csr);
               break;
            }
            if (spgemm_adapter_multiply(T_csr, P_csr, &AH_csr) != 0)
            {
               jx_printf("Proc %d: T3 AH=T*P failed at level %d\n", my_id, level);
               jx_CSRMatrixDestroy(R_csr);
               jx_CSRMatrixDestroy(T_csr);
               break;
            }
            A_H = jx_hpInithpCSRMatrix();
            jx_ParCSRMatrix *A_H_par = jx_CSRMatrixToParCSRMatrix_sp(
               jx_ParCSRMatrixComm(A_par), AH_csr);
            jx_hpCSRMatrixPar(A_H) = A_H_par;
            jx_CSRMatrixI(AH_csr) = NULL;
            jx_CSRMatrixJ(AH_csr) = NULL;
            jx_CSRMatrixData(AH_csr) = NULL;
            jx_CSRMatrixDestroy(AH_csr);
            jx_CSRMatrixDestroy(T_csr);
            jx_CSRMatrixDestroy(R_csr);
         }
         else
         {
            jx_PAMGBuildCoarseOperatorKT( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &(jx_hpCSRMatrixPar(A_H)) );
         }
      }
      else if (spmt_rap_type == 7) /* Liushuai 2026/07/12: JXPAMG sym + BeidouBLAS num RAP */
      {
         if (num_procs == 1)
         {
            printf("1111");
            jx_ParCSRMatrix *P_par = P_array[level];
            jx_ParCSRMatrix *A_par = jx_hpCSRMatrixPar(A_array[level]);
            jx_CSRMatrix *P_csr = jx_ParCSRMatrixDiag(P_par);
            jx_CSRMatrix *A_csr = jx_ParCSRMatrixDiag(A_par);
            jx_CSRMatrix *R_csr = NULL;
            jx_CSRMatrix *T_sym = NULL;
            jx_CSRMatrix *T_csr = NULL;
            jx_CSRMatrix *AH_sym = NULL;
            jx_CSRMatrix *AH_csr = NULL;

            jx_CSRMatrixTranspose(P_csr, &R_csr, 1);

            T_sym = jx_CSRMatrixMultiplySymbolic(R_csr, A_csr);
            if (!T_sym) { jx_CSRMatrixDestroy(R_csr); break; }
            T_csr = jx_CSRMatrixCreate(T_sym->num_rows, T_sym->num_cols, T_sym->num_nonzeros);
            jx_CSRMatrixInitialize(T_csr);
            jx_CSRMatrixCopy(T_sym, T_csr, 0);
            jx_CSRMatrixDestroy(T_sym);
            printf("2222");
            if (spgemm_adapter_numeric_multiply(R_csr, A_csr, T_csr) != 0)
            {
               jx_printf("Proc %d: T7 num T=R*A failed at level %d\n", my_id, level);
               jx_CSRMatrixDestroy(R_csr); jx_CSRMatrixDestroy(T_csr); break;
            }

            AH_sym = jx_CSRMatrixMultiplySymbolic(T_csr, P_csr);
            if (!AH_sym) { jx_CSRMatrixDestroy(R_csr); jx_CSRMatrixDestroy(T_csr); break; }
            AH_csr = jx_CSRMatrixCreate(AH_sym->num_rows, AH_sym->num_cols, AH_sym->num_nonzeros);
            jx_CSRMatrixInitialize(AH_csr);
            jx_CSRMatrixCopy(AH_sym, AH_csr, 0);
            jx_CSRMatrixDestroy(AH_sym);
            if (spgemm_adapter_numeric_multiply(T_csr, P_csr, AH_csr) != 0)
            {
               jx_printf("Proc %d: T7 num AH=T*P failed at level %d\n", my_id, level);
               jx_CSRMatrixDestroy(R_csr); jx_CSRMatrixDestroy(T_csr);
               jx_CSRMatrixDestroy(AH_csr); break;
            }

            A_H = jx_hpInithpCSRMatrix();
            jx_ParCSRMatrix *A_H_par = jx_CSRMatrixToParCSRMatrix_sp(jx_ParCSRMatrixComm(A_par), AH_csr);
            jx_hpCSRMatrixPar(A_H) = A_H_par;
            jx_CSRMatrixI(AH_csr) = NULL; jx_CSRMatrixJ(AH_csr) = NULL;
            jx_CSRMatrixData(AH_csr) = NULL;
            jx_CSRMatrixDestroy(AH_csr); jx_CSRMatrixDestroy(T_csr); jx_CSRMatrixDestroy(R_csr);
         }
         else
         {
            jx_PAMGBuildCoarseOperatorKT( P_array[level], jx_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &(jx_hpCSRMatrixPar(A_H)) );
         }
      }

      /*Create hierarchy block of coarse matrix*/
      #ifdef USING_HWLOC
      A_H = jx_hpCreateMatrixLevelBlock(A_H);
      #endif

      if (debug_flag == 1)
      {
         wall_time = jx_time_getWallclockSeconds() - wall_time;
         jx_printf("Proc = %d    Level = %d    Build Coarse Operator Time = %f\n", my_id,level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_rap += (jx_time_getWallclockSeconds() - tmp_wall_time);
      }

      ++ level;

      jx_ParCSRMatrixSetNumNonzeros(jx_hpCSRMatrixPar(A_H));
      jx_ParCSRMatrixSetDNumNonzeros(jx_hpCSRMatrixPar(A_H));
      A_array[level] = A_H;

      if (coarse_size <= coarse_threshold_2) {

         coarse_size_last = coarse_size;

      }

      /* print coarser operator. */
	  if (print_coarse_matrix == 1) {

         jx_sprintf(FileNameCoaMat, "cmat_%d", level);
         jx_hpCSRMatrixPrint(A_array[level], FileNameCoaMat);

      }
      
     /*-------------------------------------------------------------------------------
      *   Switch to CLJP when coarsening slows
      *                                            peghoty  2009/07/09
      *------------------------------------------------------------------------------*/ 
      
      size = ((JX_Real) fine_size )*coarse_ratio;   /* peghoty 2010/04/14 */
      if (coarsen_type > 0 && coarse_size >= (JX_Int) size)
      {
	 coarsen_type = 0;      
      }

     /*--------------------------------------------------------------------------------
      *   How to stop the loop "while"
      *                                            peghoty  2009/07/09
      *--------------------------------------------------------------------------------*/ 
      if ( (level == max_levels-1) || (coarse_size <= coarse_threshold) )
      {
         not_finished_coarsening = 0;

#if 0
         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the Strength Matrix 
            *-------------------------------------------------------------------------*/        
            jx_PAMGCreateS(A_array[level],strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
            col_offd_S_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jx_PAMGCreateSCommPkg(A_array[level], par_S, &col_offd_S_to_A);
            }

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jx_PAMGMeasureAI(par_S, A_array[level], debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         /*-------------------------------------------------------------------------
          *  compute ai-information for all-pnts and C-pnts!
          *-------------------------------------------------------------------------*/ 
         local_num_vars = jx_CSRMatrixNumRows(jx_ParCSRMatrixDiag(A_array[level]));

          for (i = 0; i < local_num_vars; i++) {
              num_vars_local++ ;
              //if (AI_measure[i] > 0.0 && AI_measure[i] < 1.0)
              if (AI_measure[i] > 1.0e-10) 
              {
                 mai_local += AI_measure[i];
                 num_ai_local++ ;
                 if (AI_measure[i] > 0.1) {
                    mai_local_valid += AI_measure[i];
                    num_ai_local_valid++ ;
                    mai_c_local += AI_measure[i];
                    num_ai_c_local++;
                 }
              }
          }
#endif
      }

      /*jx_printf("level = %d, not_finished_coarsening = %d \n", level, not_finished_coarsening);
      if (level == 1) {
        num_ai_local_0 = num_ai_local;
        num_ai_local_valid_0 = num_ai_local_valid;
        num_vars_local_0 = num_vars_local;
        num_ai_c_local_0 = num_ai_c_local;
        mai_local_0 = mai_local;
        mai_local_valid_0 = mai_local_valid;
        mai_c_local_0 = mai_c_local;
      }*/
   } // end while loop

   if (amg_print_level > 0 && ai_measure_type)
   {

      jx_MPI_Allreduce(&num_ai_local_0,&num_ai_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_local_valid_0,&num_ai_global_valid,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_vars_local_0,&num_vars_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_c_local_0,&num_ai_c_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_0,&mai_global,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_valid_0,&mai_global_valid,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_c_local_0,&mai_c_global,1,JX_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JX_Real num_ai = num_ai_global;
         JX_Real num_ai_valid = num_ai_global_valid;
         JX_Real num_var = num_vars_global;
         JX_Real num_ai_c = num_ai_c_global;
         jx_printf(" \n");
         jx_printf("Level 0: \n");
         jx_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jx_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jx_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai) jx_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid) jx_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jx_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jx_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jx_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
      }

      jx_MPI_Allreduce(&num_ai_local,&num_ai_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_local_valid,&num_ai_global_valid,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_vars_local,&num_vars_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&num_ai_c_local,&num_ai_c_global,1,JX_MPI_INT,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local,&mai_global,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_local_valid,&mai_global_valid,1,JX_MPI_REAL,MPI_SUM,comm);
      jx_MPI_Allreduce(&mai_c_local,&mai_c_global,1,JX_MPI_REAL,MPI_SUM,comm);

      if (my_id == 0)
      {
         JX_Real num_ai = num_ai_global;
         JX_Real num_ai_valid = num_ai_global_valid;
         JX_Real num_var = num_vars_global;
         JX_Real num_ai_c = num_ai_c_global;
         jx_printf(" \n");
         jx_printf("All levels (except coarsest level): \n");
         jx_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jx_printf(" - num_ai_ratio = %f\n", num_ai/num_var);
         jx_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid/num_var);
         if (num_ai>0) jx_printf(" - num_ai_c_ratio = %f\n", num_ai_c/num_ai);
         if (num_ai_valid>0) jx_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c/num_ai_valid);
         if (mai_global > 0.0) jx_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid/mai_global);
         if (mai_global > 0.0) jx_printf(" - measure_ai_c_ratio = %f\n", mai_c_global/mai_global);
         if (mai_global_valid > 0.0) jx_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global/mai_global_valid);
         jx_printf("================================================================== \n");
      }
   }

   if (level > 0)
   {
      F_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jx_hpCSRMatrixRowStarts(A_array[level]));
      jx_ParVectorInitialize(F_array[level]);
      jx_ParVectorSetPartitioningOwner(F_array[level],0);
         
      U_array[level] = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[level]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jx_hpCSRMatrixRowStarts(A_array[level]));
      jx_ParVectorInitialize(U_array[level]);
      jx_ParVectorSetPartitioningOwner(U_array[level],0);
   }
   


  /*-----------------------------------------------------------------------
   * enter all the stuff created, A[level], P[level], CF_marker[level],
   * for levels 1 through coarsest, into amg_data data structure
   *-----------------------------------------------------------------------*/

   num_levels = level + 1;
   jx_ParAMGDataNumLevels(amg_data) = num_levels;
   if (jx_ParAMGDataSmoothNumLevels(amg_data) > level)
   {
      jx_ParAMGDataSmoothNumLevels(amg_data) = level;
   }
   smooth_num_levels = jx_ParAMGDataSmoothNumLevels(amg_data);

   if (grid_relax_type[1] == 8 || grid_relax_type[2] == 8 || grid_relax_type[3] == 8 ||
       grid_relax_type[1] == 13 || grid_relax_type[2] == 13 || grid_relax_type[3] == 13 ||
       grid_relax_type[1] == 14 || grid_relax_type[2] == 14 || grid_relax_type[3] == 14)
   {
      l1_norms = jx_CTAlloc(JX_Real *, num_levels);
      jx_ParAMGDataL1Norms(amg_data) = l1_norms;
   }

   for (j = 0; j < num_levels; j ++)
   {
      if (num_threads == 1)
      {
         if (j < num_levels-1 && (grid_relax_type[1] == 8 || grid_relax_type[1] == 13 || grid_relax_type[1] == 14 ||
                                  grid_relax_type[2] == 8 || grid_relax_type[2] == 13 || grid_relax_type[2] == 14))
         {
            if (relax_order)
               jx_ParCSRComputeL1Norms(jx_hpCSRMatrixPar(A_array[j]), 4, CF_marker_array[j], &l1_norms[j]);
            else
               jx_ParCSRComputeL1Norms(jx_hpCSRMatrixPar(A_array[j]), 4, NULL, &l1_norms[j]);
         }
         else if ((grid_relax_type[3] == 8 || grid_relax_type[3] == 13 || grid_relax_type[3] == 14) && j == num_levels-1)
         {
            jx_ParCSRComputeL1Norms(jx_hpCSRMatrixPar(A_array[j]), 4, NULL, &l1_norms[j]);
         }
      }
      else
      {
         if (j < num_levels-1 && (grid_relax_type[1] == 8 || grid_relax_type[1] == 13 || grid_relax_type[1] == 14 ||
                                  grid_relax_type[2] == 8 || grid_relax_type[2] == 13 || grid_relax_type[2] == 14))
         {
            if (relax_order)
               jx_ParCSRComputeL1NormsThreads(jx_hpCSRMatrixPar(A_array[j]), 4, num_threads, CF_marker_array[j] , &l1_norms[j]);
            else
               jx_ParCSRComputeL1NormsThreads(jx_hpCSRMatrixPar(A_array[j]), 4, num_threads, NULL, &l1_norms[j]);
         }
         else if ((grid_relax_type[3] == 8 || grid_relax_type[3] == 13 || grid_relax_type[3] == 14) && j == num_levels-1)
         {
            jx_ParCSRComputeL1NormsThreads(jx_hpCSRMatrixPar(A_array[j]), 4, num_threads, NULL, &l1_norms[j]);
         }
      }
   }

   for (j = 0; j < smooth_num_levels; j++)
   {
      if (smooth_type == 9 || smooth_type == 19) // Euclid smoothers
      {
         JX_EuclidCreate(comm, &smoother[j]);
         if (euclidfile)
         {
            JX_EuclidSetParamsFromFile(smoother[j], euclidfile);
         }
         JX_EuclidSetLevel(smoother[j], eu_level);
         if (eu_bj)
         {
            JX_EuclidSetBJ(smoother[j], eu_bj);
         }
         if (eu_sparse_A)
         {
            JX_EuclidSetSparseA(smoother[j], eu_sparse_A);
         }
         JX_EuclidSetup(smoother[j], (JX_hpCSRMatrix)A_array[j]);
      }
      else if (smooth_type == 5 || smooth_type == 15) // ILU smoothers
      {
         JX_ILUCreate(&smoother[j]);
         if (eu_bj)
         {
            JX_ILUSetType(smoother[j], 0);
         }
         else
         {
            JX_ILUSetType(smoother[j], 30);
         }
         JX_ILUSetMaxIter(smoother[j], 1);
         /* set tol to zero since we are doing just 1 iteration */
         JX_ILUSetTol(smoother[j], 0.);
         JX_ILUSetLogging(smoother[j], 0);
         JX_ILUSetPrintLevel(smoother[j], 0);
         JX_ILUSetLevelOfFill(smoother[j], eu_level);
         JX_ILUSetup(smoother[j], (jx_hpCSRMatrix *)A_array[j], NULL, NULL);
      }
   }

   if (amg_logging > 1) 
   {
      Residual_array = jx_ParVectorCreate(jx_hpCSRMatrixComm(A_array[0]),
                                          jx_hpCSRMatrixGlobalNumRows(A_array[0]),
                                          jx_hpCSRMatrixRowStarts(A_array[0]));
      jx_ParVectorInitialize(Residual_array);
      jx_ParVectorSetPartitioningOwner(Residual_array,0);
      jx_ParAMGDataResidual(amg_data) = Residual_array;
   }
   else
   {
      jx_ParAMGDataResidual(amg_data) = NULL;
   }

   jx_TFree(opt_icor); /* Yue Xiaoqiang 2012/10/12 */
   if (Pmarkers_global_rap)
   {
      jx_TFree(Pmarkers_global_rap); /* Yue Xiaoqiang 2012/10/17 */
   }
   if (Pmarkers_global_interp)
   {
      jx_TFree(Pmarkers_global_interp); /* Yue Xiaoqiang 2012/10/17 */
   }

  /*--------------------------------------------------------------
   *    Print some stuff
   *-------------------------------------------------------------*/

   //if (amg_print_level == 1 || amg_print_level == 3)
   /*change by mrz 2022/10/26*/
   if (amg_print_level > 1)  
   {
      /* Write the SETUP parameters */
      jx_hpPAMGSetupStatus(amg_data, hp_matrix);
   }

   if (wall_time_option == 1)
   {
      jx_printf("\n\nProc = %d, Coarsen Time = %f\n", my_id, wall_time_coarsen);
      jx_printf("Proc = %d, Build Coarse Operator Time = %f\n", my_id, wall_time_rap);
      jx_printf("Proc = %d, Build Interp Time = %f\n\n", my_id, wall_time_interp);
   }

   if (my_id == 0 && amg_print_level > 1)
   {
      /* Write the SOLVE parameters */
      jx_PAMGSolveStatus(amg_data); /* Yue Xiaoqiang 2014/04/12 */
   }

   return(Setup_err_flag);
}
