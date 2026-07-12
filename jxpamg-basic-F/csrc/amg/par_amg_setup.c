//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_amg_setup.c
 *  Date: 2011/09/03
 */

#include "jxf_pamg.h"
#include "jxf_euclid.h"
#include "jxf_ilu.h"

JXF_Int jxf_error_flag = 0;

JXF_UInt jxf_long_size_length_rap = 0;         /* Yue Xiaoqiang 2012/10/12 */
JXF_UInt jxf_long_size_length_interp = 0;      /* Yue Xiaoqiang 2012/10/12 */
JXF_UInt jxf_real_long_size_length_rap = 0;    /* Yue Xiaoqiang 2012/10/12 */
JXF_UInt jxf_real_long_size_length_interp = 0; /* Yue Xiaoqiang 2012/10/12 */

extern JXF_Int *Jxf_Pmarkers_global_rap;    /* Yue Xiaoqiang 2012/10/17 */
extern JXF_Int *Jxf_Pmarkers_global_interp; /* Yue Xiaoqiang 2012/10/17 */

/*!
 * \fn JXF_Int jxf_PAMGSetup
 * \brief Setup phase of PAMG solver.
 * \date 2011/09/03
 */

JXF_Int
jxf_hpPAMGSetup(void *amg_vdata,
               jxf_hpCSRMatrix *hp_matrix)
{
   MPI_Comm comm = jxf_hpCSRMatrixComm(hp_matrix);
   jxf_ParAMGData *amg_data = amg_vdata;
   /* Data Structure variables */
   jxf_hpCSRMatrix **A_array;
   jxf_ParVector **F_array;
   jxf_ParVector **U_array;
   jxf_ParVector *Vtemp;
   jxf_ParVector *Ztemp = NULL;
   jxf_ParCSRMatrix **P_array;
   jxf_ParCSRMatrix **R_array;
   jxf_ParVector *Residual_array;
   JXF_Real **AI_measure_array;
   JXF_Int **CF_marker_array;
   JXF_Int **dof_func_array;
   JXF_Int *dof_func;
   JXF_Int *col_offd_S_to_A;
   JXF_Int *col_offd_Sabs_to_A = NULL;
   JXF_Int *AIR_maxsize_ls;
   JXF_Real strong_threshold;
   JXF_Real AIR_strong_th;
   JXF_Real max_row_sum;
   JXF_Real trunc_factor;
   JXF_Real S_commpkg_switch;

   JXF_Int max_levels;
   JXF_Int amg_logging;
   JXF_Int amg_print_level;
   JXF_Int debug_flag;
   JXF_Int local_num_vars;
   JXF_Int P_max_elmts;
   JXF_Int R_max_size;

   JXF_Solver *smoother = NULL;
   JXF_Int smooth_type = jxf_ParAMGDataSmoothType(amg_data);
   JXF_Int smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data);
   char *euclidfile;
   JXF_Int eu_level;
   JXF_Int eu_bj;
   JXF_Real eu_sparse_A;

   /* Local variables */
   JXF_Real *AI_measure;
   JXF_Real *AI_measure2;
   JXF_Int *CF_marker;
   JXF_Int *CFN_marker = NULL;
   jxf_ParCSRMatrix *par_S;
   jxf_ParCSRMatrix *Sabs = NULL;
   jxf_ParCSRMatrix *par_S2;
   jxf_ParCSRMatrix *par_P;
   jxf_ParCSRMatrix *R;
   jxf_hpCSRMatrix *A_H;

   JXF_Int wall_time_option;         /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real wall_time_coarsen = 0.0; /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real wall_time_rap = 0.0;     /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real wall_time_interp = 0.0;  /* newly added Yue Xiaoqiang 2015/09/30 */
   JXF_Real tmp_wall_time = 0.0;     /* newly added Yue Xiaoqiang 2015/09/30 */

   JXF_Int old_num_levels, num_levels;
   JXF_Int level;
   JXF_Int local_size, i;
   JXF_Int first_local_row;
   JXF_Int coarse_size_last = -1;
   JXF_Int coarse_size;
   JXF_Int coarsen_type;
   JXF_Int interp_type;
   JXF_Int restri_type;
   JXF_Int measure_type;
   JXF_Int setup_type;
   JXF_Int fine_size;
   JXF_Int rest, tms, indx;
   JXF_Int not_finished_coarsening = 1;
   JXF_Int Setup_err_flag = 0;
   JXF_Int coarse_threshold;   /* we don't fix it to be 9.  peghoty 2010/04/14 */
   JXF_Int coarse_threshold_2; /* by xu: 2013/11/25 */
   JXF_Int j, k;
   JXF_Int num_procs, my_id, num_threads;
   JXF_Int *grid_relax_type = jxf_ParAMGDataGridRelaxType(amg_data);
   JXF_Int relax_order = jxf_ParAMGDataRelaxOrder(amg_data);
   JXF_Int num_functions = jxf_ParAMGDataNumFunctions(amg_data);
   JXF_Int spmt_rap_type = jxf_ParAMGDataSpMtRapType(amg_data);
   JXF_Int ai_measure_type = jxf_ParAMGDataAIMeasureType(amg_data);
   JXF_Int num_paths = jxf_ParAMGDataNumPaths(amg_data);
   JXF_Int agg_num_levels = jxf_ParAMGDataAggNumLevels(amg_data);
   JXF_Int agg_interp_type = jxf_ParAMGDataAggInterpType(amg_data);
   JXF_Int agg_P_max_elmts = jxf_ParAMGDataAggPMaxElmts(amg_data);
   // JXF_Int       agg_P12_max_elmts = jxf_ParAMGDataAggP12MaxElmts(amg_data);
   JXF_Real agg_trunc_factor = jxf_ParAMGDataAggTruncFactor(amg_data);
   // JXF_Real    agg_P12_trunc_factor = jxf_ParAMGDataAggP12TruncFactor(amg_data);
   JXF_Int rap2 = jxf_ParAMGDataRAP2(amg_data);
   JXF_Int keepTranspose = jxf_ParAMGDataKeepTranspose(amg_data);
   JXF_Int print_coarse_matrix = jxf_ParAMGDataPrintCoarseSystem(amg_data);
   JXF_Int *opt_icor;
   JXF_Int *coarse_dof_func;
   JXF_Int *coarse_pnts_global;
   JXF_Int *coarse_pnts_global1;
   JXF_Real size;
   JXF_Real coarse_ratio;          /* newly added peghoty 2010/04/14 */
   JXF_Real wall_time = 0.0;       /* for debugging instrumentation */
   JXF_Int measure_type_rlx;       // newly added peghoty 2010/05/29
   JXF_Int number_syn_rlx;         // newly added peghoty 2010/05/29
   JXF_Real measure_threshold_rlx; // newly added peghoty 2010/05/29
   JXF_Real **l1_norms = NULL;

   /* coarse matrices */
   char FileNameCoaMat[256];

   /* ai statistic information*/
   JXF_Int num_vars_local = 0, num_vars_global;
   JXF_Int num_ai_local = 0, num_ai_global;
   JXF_Int num_ai_local_valid = 0, num_ai_global_valid;
   JXF_Int num_ai_c_local = 0, num_ai_c_global;
   JXF_Real mai_local = 0.0, mai_global;
   JXF_Real mai_local_valid = 0.0, mai_global_valid;
   JXF_Real mai_c_local = 0.0, mai_c_global;

   JXF_Int num_vars_local_0 = 0;
   JXF_Int num_ai_local_0 = 0;
   JXF_Int num_ai_local_valid_0 = 0;
   JXF_Int num_ai_c_local_0 = 0;
   JXF_Real mai_local_0 = 0.0;
   JXF_Real mai_local_valid_0 = 0.0;
   JXF_Real mai_c_local_0 = 0.0;
   JXF_Real mai_threshold = 0.1;
   JXF_Real gai_threshold = 0.3;

   jxf_MPI_Comm_size(comm, &num_procs);
   jxf_MPI_Comm_rank(comm, &my_id);
   if ((num_procs > 1) && (spmt_rap_type != 1)) /* Yue Xiaoqiang 2012/10/12 */
   {
      spmt_rap_type = 1;
   }
   num_threads = jxf_NumThreads();                      /* Yue Xiaoqiang 2012/10/12 */
   opt_icor = jxf_CTAlloc(JXF_Int, 5 * num_threads + 2); /* Yue Xiaoqiang 2012/10/12 */

   wall_time_option = jxf_ParAMGDataWallTimeOption(amg_data);
   old_num_levels = jxf_ParAMGDataNumLevels(amg_data);
   max_levels = jxf_ParAMGDataMaxLevels(amg_data);
   amg_logging = jxf_ParAMGDataLogging(amg_data);
   amg_print_level = jxf_ParAMGDataPrintLevel(amg_data);
   coarsen_type = jxf_ParAMGDataCoarsenType(amg_data);
   measure_type = jxf_ParAMGDataMeasureType(amg_data);
   setup_type = jxf_ParAMGDataSetupType(amg_data);
   debug_flag = jxf_ParAMGDataDebugFlag(amg_data);
   dof_func = jxf_ParAMGDataDofFunc(amg_data);
   interp_type = jxf_ParAMGDataInterpType(amg_data);
   restri_type = jxf_ParAMGDataRestriction(amg_data);
   AIR_strong_th = jxf_ParAMGDataAIRStrongTh(amg_data);
   euclidfile = jxf_ParAMGDataEuclidFile(amg_data);
   eu_level = jxf_ParAMGDataEuLevel(amg_data);
   eu_bj = jxf_ParAMGDataEuBJ(amg_data);
   eu_sparse_A = jxf_ParAMGDataEuSparseA(amg_data);
   coarse_threshold = jxf_ParAMGDataCoarseThreshold(amg_data);          /* newly added peghoty 2010/04/14 */
   coarse_ratio = jxf_ParAMGDataCoarseRatio(amg_data);                  /* newly added peghoty 2010/04/14 */
   measure_type_rlx = jxf_ParAMGDataMeasureTypeRlx(amg_data);           /* newly added peghoty 2010/05/29 */
   number_syn_rlx = jxf_ParAMGDataNumberSynRlx(amg_data);               /* newly added peghoty 2010/05/29 */
   measure_threshold_rlx = jxf_ParAMGDataMeasureThresholdRlx(amg_data); /* newly added peghoty 2010/05/29 */

   coarse_threshold_2 = -1;

   jxf_ParCSRMatrixSetNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix));
   jxf_ParCSRMatrixSetDNumNonzeros(jxf_hpCSRMatrixPar(hp_matrix));
   jxf_ParAMGDataNumVariables(amg_data) = jxf_hpCSRMatrixNumRows(hp_matrix);

   if (setup_type == 0)
   {
      return (Setup_err_flag);
   }

   par_S = NULL;
   A_H = NULL;

   A_array = jxf_hpAMGDataAArray(amg_data);
   P_array = jxf_ParAMGDataPArray(amg_data);
   R_array = jxf_ParAMGDataRArray(amg_data);
   AIR_maxsize_ls = jxf_ParAMGDataAIRMaxSizeLS(amg_data);
   AI_measure_array = jxf_ParAMGDataAIMeasureArray(amg_data);
   CF_marker_array = jxf_ParAMGDataCFMarkerArray(amg_data);
   dof_func_array = jxf_ParAMGDataDofFuncArray(amg_data);
   local_size = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(hp_matrix));

   grid_relax_type[3] = jxf_ParAMGDataUserCoarseRelaxType(amg_data);

   if (A_array || P_array || R_array || AIR_maxsize_ls || AI_measure_array || CF_marker_array || dof_func_array)
   {
      for (j = 1; j < old_num_levels; j++)
      {
         if (A_array[j])
         {
            jxf_hpCSRMatrixDestroy(A_array[j]);
            A_array[j] = NULL;
         }

         if (dof_func_array[j])
         {
            jxf_TFree(dof_func_array[j]);
            dof_func_array[j] = NULL;
         }
      }

      for (j = 0; j < old_num_levels - 1; j++)
      {
         if (P_array[j])
         {
            jxf_ParCSRMatrixDestroy(P_array[j]);
            P_array[j] = NULL;
         }
         if (R_array[j])
         {
            jxf_ParCSRMatrixDestroy(R_array[j]);
            R_array[j] = NULL;
         }
      }

      jxf_TFree(AIR_maxsize_ls);
      AIR_maxsize_ls = NULL;

      /*-------------------------------------------------------------------
       *  Special case use of CF_marker_array when old_num_levels = 1
       *  requires us to attempt this deallocation every time
       *------------------------------------------------------------------*/

      if (CF_marker_array[0])
      {
         jxf_TFree(CF_marker_array[0]);
         CF_marker_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels - 1; j++)
      {
         if (CF_marker_array[j])
         {
            jxf_TFree(CF_marker_array[j]);
            CF_marker_array[j] = NULL;
         }
      }

      /* for AI_measure: added by xwxu */
      if (AI_measure_array[0])
      {
         jxf_TFree(AI_measure_array[0]);
         AI_measure_array[0] = NULL;
      }

      for (j = 1; j < old_num_levels - 1; j++)
      {
         if (AI_measure_array[j])
         {
            jxf_TFree(AI_measure_array[j]);
            AI_measure_array[j] = NULL;
         }
      }
   }

   if (A_array == NULL)
   {
      A_array = jxf_CTAlloc(jxf_hpCSRMatrix *, max_levels);
   }

   if (P_array == NULL && max_levels > 1)
   {
      P_array = jxf_CTAlloc(jxf_ParCSRMatrix *, max_levels - 1);
   }

   /* If retri_type != 0, R != P^T, allocate R matrices */
   if (restri_type)
   {
      if (R_array == NULL && max_levels > 1)
      {
         R_array = jxf_CTAlloc(jxf_ParCSRMatrix *, max_levels - 1);
         AIR_maxsize_ls = jxf_CTAlloc(JXF_Int, max_levels - 1);
      }
   }

   if (AI_measure_array == NULL)
   {
      AI_measure_array = jxf_CTAlloc(JXF_Real *, max_levels);
   }
   if (CF_marker_array == NULL)
   {
      CF_marker_array = jxf_CTAlloc(JXF_Int *, max_levels);
   }
   if (dof_func_array == NULL)
   {
      dof_func_array = jxf_CTAlloc(JXF_Int *, max_levels);
   }

   if (num_functions > 1 && dof_func == NULL)
   {
      first_local_row = jxf_hpCSRMatrixFirstRowIndex(hp_matrix);
      dof_func = jxf_CTAlloc(JXF_Int, local_size);
      rest = first_local_row - ((first_local_row / num_functions) * num_functions);
      indx = num_functions - rest;
      if (rest == 0)
      {
         indx = 0;
      }
      k = num_functions - 1;
      for (j = indx - 1; j > -1; j--)
      {
         dof_func[j] = k--;
      }
      tms = local_size / num_functions;
      if (tms * num_functions + indx > local_size)
      {
         tms--;
      }
      for (j = 0; j < tms; j++)
      {
         for (k = 0; k < num_functions; k++)
         {
            dof_func[indx++] = k;
         }
      }
      k = 0;
      while (indx < local_size)
      {
         dof_func[indx++] = k++;
      }
      jxf_ParAMGDataDofFunc(amg_data) = dof_func;
   }

   A_array[0] = hp_matrix;

   dof_func_array[0] = dof_func;
   jxf_ParAMGDataAIMeasureArray(amg_data) = AI_measure_array;
   jxf_ParAMGDataCFMarkerArray(amg_data) = CF_marker_array;
   jxf_ParAMGDataDofFuncArray(amg_data) = dof_func_array;
   jxf_hpAMGDataAArray(amg_data) = A_array;
   jxf_ParAMGDataPArray(amg_data) = P_array;
   /* If R != P^T */
   if (restri_type)
   {
      jxf_ParAMGDataRArray(amg_data) = R_array;
      jxf_ParAMGDataAIRMaxSizeLS(amg_data) = AIR_maxsize_ls;
   }
   else
   {
      jxf_ParAMGDataRArray(amg_data) = P_array;
   }

   Vtemp = jxf_ParAMGDataVtemp(amg_data);

   if (Vtemp != NULL)
   {
      jxf_ParVectorDestroy(Vtemp);
      Vtemp = NULL;
   }

   Vtemp = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[0]),
                              jxf_hpCSRMatrixGlobalNumRows(A_array[0]),
                              jxf_hpCSRMatrixRowStarts(A_array[0]));
   jxf_ParVectorInitialize(Vtemp);
   jxf_ParVectorSetPartitioningOwner(Vtemp, 0);
   jxf_ParAMGDataVtemp(amg_data) = Vtemp;

   if (num_threads > 1)
   {
      Ztemp = jxf_ParAMGDataZtemp(amg_data);

      if (Ztemp != NULL)
      {
         jxf_ParVectorDestroy(Ztemp);
         Ztemp = NULL;
      }

      for (j = 1; j < 4; j++)
      {
         if (grid_relax_type[j] == 3 || grid_relax_type[j] == 4 || grid_relax_type[j] == 6 ||
             grid_relax_type[j] == 8 || grid_relax_type[j] == 13 || grid_relax_type[j] == 14)
         {
            Ztemp = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[0]), jxf_hpCSRMatrixGlobalNumRows(A_array[0]), jxf_hpCSRMatrixRowStarts(A_array[0]));
            jxf_ParVectorInitialize(Ztemp);
            jxf_ParVectorSetPartitioningOwner(Ztemp, 0);
            jxf_ParAMGDataZtemp(amg_data) = Ztemp;
            break;
         }
      }
   }

   F_array = jxf_ParAMGDataFArray(amg_data);
   U_array = jxf_ParAMGDataUArray(amg_data);

   if (F_array != NULL || U_array != NULL)
   {
      for (j = 1; j < old_num_levels; j++)
      {
         if (F_array[j] != NULL)
         {
            jxf_ParVectorDestroy(F_array[j]);
            F_array[j] = NULL;
         }
         if (U_array[j] != NULL)
         {
            jxf_ParVectorDestroy(U_array[j]);
            U_array[j] = NULL;
         }
      }
   }

   if (F_array == NULL)
   {
      F_array = jxf_CTAlloc(jxf_ParVector *, max_levels);
   }
   if (U_array == NULL)
   {
      U_array = jxf_CTAlloc(jxf_ParVector *, max_levels);
   }

   /*---------------------------------------------------------
    *  They have been moved to "amg_solve" since they are
    *  not used here.  peghoty 2009/07/27
    *--------------------------------------------------------*/
   // F_array[0] = f;
   // U_array[0] = u;

   jxf_ParAMGDataFArray(amg_data) = F_array;
   jxf_ParAMGDataUArray(amg_data) = U_array;

   /*----------------------------------------------------------
    *   Initialize jxf_ParAMGData
    *---------------------------------------------------------*/
   not_finished_coarsening = 1;
   level = 0;
   strong_threshold = jxf_ParAMGDataStrongThreshold(amg_data);
   max_row_sum = jxf_ParAMGDataMaxRowSum(amg_data);
   trunc_factor = jxf_ParAMGDataTruncFactor(amg_data);
   P_max_elmts = jxf_ParAMGDataPMaxElmts(amg_data);
   S_commpkg_switch = jxf_ParAMGDataSCommPkgSwitch(amg_data);
   if (smooth_num_levels > level)
   {
      smoother = jxf_CTAlloc(JXF_Solver, smooth_num_levels);
      jxf_ParAMGDataSmoother(amg_data) = smoother;
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

   if (amg_print_level > 0 && my_id == 0 && ai_measure_type)
   {
      jxf_printf(" \n");
      jxf_printf("================================================================== \n");
      jxf_printf("+++++++++++++ multi-scale/AI infomation for levels +++++++++++++++ \n");
   }
   while (not_finished_coarsening)
   {
      /* 建议在循环内部定义临时指针，并存入数组 */
      JXF_Real *measure_gai = NULL;
      JXF_Int *marker_gai = NULL;
      JXF_Int *side_gai = NULL;
      JXF_Real tau_actual = 0.0;

      fine_size = jxf_hpCSRMatrixGlobalNumRows(A_array[level]);

      if (level > 0)
      {
         F_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                             jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jxf_hpCSRMatrixRowStarts(A_array[level]));
         jxf_ParVectorInitialize(F_array[level]);
         jxf_ParVectorSetPartitioningOwner(F_array[level], 0);

         U_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                             jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                             jxf_hpCSRMatrixRowStarts(A_array[level]));
         jxf_ParVectorInitialize(U_array[level]);
         jxf_ParVectorSetPartitioningOwner(U_array[level], 0);
      }

      /*--------------------------------------------------------------
       *  Select coarse-grid points on 'level' : returns CF_marker
       *  for the level.  Returns strength matrix, par_S
       *  Returns AI_Measure.
       *--------------------------------------------------------------*/

      if (debug_flag == 1)
         wall_time = jxf_time_getWallclockSeconds();
      if (debug_flag == 3)
      {
         jxf_printf("\n ===== Proc = %d     Level = %d  =====\n", my_id, level);
         fflush(NULL);
      }

      if (wall_time_option == 1)
         tmp_wall_time = jxf_time_getWallclockSeconds();

      if (max_levels > 1)
      {
         local_num_vars = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level]));

         /*--------------------------------------------------------------------------
          *  Get the Strength Matrix
          *-------------------------------------------------------------------------*/
         jxf_PAMGCreateS(jxf_hpCSRMatrixPar(A_array[level]), strong_threshold, max_row_sum, num_functions, dof_func_array[level], &par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }

         /* for AIR, need absolute value SOC */
         if ((restri_type == 1) || (restri_type == 2))
         {
            jxf_PAMGCreateSabs(jxf_hpCSRMatrixPar(A_array[level]), AIR_strong_th, 1.0, num_functions, dof_func_array[level], &Sabs);

            col_offd_Sabs_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), Sabs, &col_offd_Sabs_to_A);
            }
         }

         if (ai_measure_type == 1)
         {

            /*--------------------------------------------------------------------------
             *  Get the AI Measure
             *-------------------------------------------------------------------------*/
            jxf_PAMGMeasureAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);
            /*--------------------------------------------
             *  store the AI array
             *-------------------------------------------*/
            AI_measure_array[level] = AI_measure;
         }

         jxf_PAMGMeasureGAI_Base(par_S,
                                jxf_hpCSRMatrixPar(A_array[level]),
                                level,
                                debug_flag,
                                gai_threshold, // 这里传入你的阈值 tau
                                &measure_gai,
                                &marker_gai,
                                &side_gai,
                                &tau_actual);

         {
            JXF_Int local_gai_count = 0;
            for (i = 0; i < local_num_vars; i++)
            {
               if (marker_gai[i])
                  local_gai_count++;
            }
            printf("[Level %d] GAI points: %d / %d (%.2f%%)\n",
                   level, local_gai_count, local_num_vars,
                   (float)local_gai_count / local_num_vars * 100.0);
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

         // jxf_printf("=========== coarsen_type = %d \n", coarsen_type);
         // jxf_printf("myid is %d,coarsen_type is %d,before Coarsen run success\n",my_id,coarsen_type);
         if (coarsen_type == 6)
         {
            jxf_PAMGCoarsenFalgout(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 96) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenFalgoutAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 8)
         {
            jxf_PAMGCoarsenPMIS(par_S, jxf_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 10)
         {
            jxf_PAMGCoarsenHMIS(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 910) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenHMISAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 1 || coarsen_type == 2 ||
                  coarsen_type == 3 || coarsen_type == 4 ||
                  coarsen_type == 5 || coarsen_type == 11)
         {
            jxf_PAMGCoarsenRuge(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, coarsen_type, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 991) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 11, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 992) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 1, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 993) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenRugeAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, measure_type, 3, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 0)
         {
            jxf_PAMGCoarsen(par_S, jxf_hpCSRMatrixPar(A_array[level]), 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 990) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 90)
         {
            jxf_PAMGCoarsenRCLJP(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type_rlx, number_syn_rlx,
                                measure_threshold_rlx, 0, debug_flag, &CF_marker);
         }
         else if (coarsen_type == 91)
         {
            jxf_PAMGCoarsenRRS0(par_S, jxf_hpCSRMatrixPar(A_array[level]), measure_type, measure_type_rlx,
                               number_syn_rlx, measure_threshold_rlx, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 98) && (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenPMISAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], 0, debug_flag, &CF_marker);
         }
         else if ((coarsen_type == 908 || coarsen_type == 918 || coarsen_type == 928 || coarsen_type == 938 ||
                   coarsen_type == 968) &&
                  (ai_measure_type == 1))
         {
            jxf_PAMGCoarsenXML(par_S, jxf_hpCSRMatrixPar(A_array[level]), AI_measure_array[level], measure_type, coarsen_type, debug_flag, &CF_marker);
         }

         if (level < agg_num_levels)
         {
            jxf_PAMGCoarseParms(comm, local_num_vars, 1,
                               dof_func_array[level], CF_marker,
                               &coarse_dof_func, &coarse_pnts_global1);
            jxf_PAMGCreate2ndS(par_S, CF_marker, num_paths, coarse_pnts_global1, &par_S2);
            /*--------------------------------------------------------------------------
             *  Get the AI Measure
             *-------------------------------------------------------------------------*/
            if (ai_measure_type == 1)
            {
               jxf_PAMGMeasureAI(par_S2, par_S2, debug_flag, &AI_measure2);
            }

            if (coarsen_type == 10)
            {
               jxf_PAMGCoarsenHMIS(par_S2, par_S2, measure_type + 3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 910) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenHMISAI(par_S2, par_S2, AI_measure2, 0, measure_type + 3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 8)
            {
               jxf_PAMGCoarsenPMIS(par_S2, par_S2, 3, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 98) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenPMISAI(par_S2, par_S2, AI_measure2, 3, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 6)
            {
               jxf_PAMGCoarsenFalgout(par_S2, par_S2, measure_type, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 96) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenFalgoutAI(par_S2, par_S2, AI_measure2, measure_type, debug_flag, &CFN_marker);
            }
            else if (coarsen_type == 0)
            {
               jxf_PAMGCoarsen(par_S2, par_S2, 0, debug_flag, &CFN_marker);
            }
            else if ((coarsen_type == 990) && (ai_measure_type == 1))
            {
               jxf_PAMGCoarsenAI(par_S2, par_S2, AI_measure2, 0, debug_flag, &CFN_marker);
            }
            jxf_ParCSRMatrixDestroy(par_S2);
            if (ai_measure_type == 1)
            {
               jxf_TFree(AI_measure2);
            }
         }

         /* Here for changes of min_coarse_size */

         if (level < agg_num_levels)
         {
            if (agg_interp_type == 4)
            {
               jxf_PAMGCorrectCFMarker(CF_marker, local_num_vars, CFN_marker);
               jxf_TFree(coarse_pnts_global1);
               jxf_TFree(CFN_marker);
            }
         }

         /*--------------------------------------------
          *  store the CF array
          *-------------------------------------------*/
         CF_marker_array[level] = CF_marker;

         if (debug_flag == 1)
         {
            wall_time = jxf_time_getWallclockSeconds() - wall_time;
            jxf_printf("Proc = %d    Level = %d    Coarsen Time = %f\n", my_id, level, wall_time);
            fflush(NULL);
         }

         if (wall_time_option == 1)
         {
            wall_time_coarsen += (jxf_time_getWallclockSeconds() - tmp_wall_time);
         }

         /*-------------------------------------------------------------------------
          *  Get the coarse parameters
          *-------------------------------------------------------------------------*/
         jxf_PAMGCoarseParms(comm, local_num_vars, num_functions,
                            dof_func_array[level], CF_marker,
                            &coarse_dof_func, &coarse_pnts_global);

         dof_func_array[level + 1] = NULL;
         if (num_functions > 1)
            dof_func_array[level + 1] = coarse_dof_func;

#ifdef JXF_NO_GLOBAL_PARTITION
         if (my_id == (num_procs - 1))
         {
            coarse_size = coarse_pnts_global[1];
         }
         jxf_MPI_Bcast(&coarse_size, 1, JXF_MPI_INT, num_procs - 1, comm);
#else
         coarse_size = coarse_pnts_global[num_procs];
#endif

         if (coarse_size <= coarse_threshold && coarse_size_last >= coarse_size)
         {
            break;
         }
      }
      else /* max_levels = 1 */
      {

         /*--------------------------------------------------------------------------
          *  Get the Strength Matrix
          *-------------------------------------------------------------------------*/
         jxf_PAMGCreateS(jxf_hpCSRMatrixPar(A_array[level]), strong_threshold, max_row_sum, num_functions, dof_func_array[level], &par_S);
         col_offd_S_to_A = NULL;
         if (strong_threshold > S_commpkg_switch)
         {
            jxf_PAMGCreateSCommPkg(jxf_hpCSRMatrixPar(A_array[level]), par_S, &col_offd_S_to_A);
         }

         if (ai_measure_type == 1)
         {

            /*--------------------------------------------------------------------------
             *  Get the AI Measure
             *-------------------------------------------------------------------------*/
            jxf_PAMGMeasureAI(par_S, jxf_hpCSRMatrixPar(A_array[level]), debug_flag, &AI_measure);

            /*--------------------------------------------
             *  store the AI array
             *-------------------------------------------*/
            AI_measure_array[level] = AI_measure;
         }
         par_S = NULL;
         coarse_pnts_global = NULL;
         CF_marker = jxf_CTAlloc(JXF_Int, local_size);
         for (i = 0; i < local_size; i++)
         {
            CF_marker[i] = 1;
         }
         CF_marker_array[level] = CF_marker;
         coarse_size = fine_size;

         local_num_vars = jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level]));
      }

      /*-------------------------------------------------------------------------
       *  compute ai-information for all-pnts and C-pnts!
       *-------------------------------------------------------------------------*/

      if ((max_levels == 1 || level != max_levels - 1) && (ai_measure_type == 1))
      {
         for (i = 0; i < local_num_vars; i++)
         {
            num_vars_local++;
            // if (AI_measure[i] > 0.0 && AI_measure[i] < 1.0)
            if (AI_measure[i] > 1.0e-10)
            {
               mai_local += AI_measure[i];
               num_ai_local++;
               if (AI_measure[i] >= mai_threshold)
               {
                  mai_local_valid += AI_measure[i];
                  num_ai_local_valid++;
                  if (CF_marker[i] == 1)
                  {
                     // if (level == 0 ) jxf_printf("i_c = %d, ai_measure= %f \n", i, AI_measure[i]);
                     mai_c_local += AI_measure[i];
                     num_ai_c_local++;
                  }
               }
            }
         }
      }

      if (level == 0)
      {
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

      if ((coarse_size == 0) || (coarse_size == fine_size))
      {
         JXF_Int *num_grid_sweeps = jxf_ParAMGDataNumGridSweeps(amg_data);
         JXF_Int **grid_relax_points = jxf_ParAMGDataGridRelaxPoints(amg_data);
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
            jxf_ParCSRMatrixDestroy(par_S);
         }
         jxf_TFree(coarse_pnts_global);
         if (level > 0)
         {
            /*-------------------------------------------------------------
             * Note special case treatment of CF_marker is necessary
             * to do CF relaxation correctly when num_levels = 1
             *------------------------------------------------------------*/
            jxf_TFree(CF_marker_array[level]);
            jxf_ParVectorDestroy(F_array[level]);
            jxf_ParVectorDestroy(U_array[level]);
         }
         break;
      }

      /* Build restriction */
      if (restri_type)
      {
         /* !!! Ensure that CF_marker contains -1 or 1 !!! */
         for (i = 0; i < jxf_CSRMatrixNumRows(jxf_hpCSRMatrixDiag(A_array[level])); i++)
         {
            CF_marker[i] = CF_marker[i] > 0 ? 1 : -1;
         }
         if (restri_type == 1) /* distance-1 AIR */
         {
            jxf_PAMGBuildRestrAIR(jxf_hpCSRMatrixPar(A_array[level]), CF_marker,
                                 Sabs, coarse_pnts_global, num_functions,
                                 dof_func_array[level],
                                 debug_flag, trunc_factor, P_max_elmts,
                                 col_offd_Sabs_to_A, &R, &R_max_size);
         }
         else if (restri_type == 2) /* distance-2 AIR */
         {
            jxf_PAMGBuildRestrDist2AIR(jxf_hpCSRMatrixPar(A_array[level]), CF_marker,
                                      Sabs, coarse_pnts_global, num_functions,
                                      dof_func_array[level],
                                      debug_flag, trunc_factor, P_max_elmts,
                                      col_offd_Sabs_to_A, &R, &R_max_size);
         }
         else /* R=[0 I_c] */
         {
            jxf_PAMGBuildOPntRestr(jxf_hpCSRMatrixPar(A_array[level]), CF_marker, coarse_pnts_global, &R, &R_max_size);
         }
         if (Sabs)
         {
            jxf_ParCSRMatrixDestroy(Sabs);
            Sabs = NULL;
         }
         jxf_TFree(col_offd_Sabs_to_A);
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

      if (debug_flag == 1)
         wall_time = jxf_time_getWallclockSeconds();

      if (wall_time_option == 1)
         tmp_wall_time = jxf_time_getWallclockSeconds();

      if (level < agg_num_levels)
      {
         if (agg_interp_type == 4)
         {
            jxf_PAMGBuildMultipass(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S, coarse_pnts_global,
                                  num_functions, dof_func_array[level], debug_flag,
                                  agg_trunc_factor, agg_P_max_elmts, 0, col_offd_S_to_A, &par_P);
         }
      }
      else
      {
         if (interp_type == 0) // classical modified interpolation
         {
            if (spmt_rap_type == 1)
            {
               jxf_PAMGBuildInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                  coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag, trunc_factor,
                                  P_max_elmts, col_offd_S_to_A, &par_P);
            }
            else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
            {
               jxf_PAMGBuildInterp1(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level], par_S,
                                   coarse_pnts_global, num_functions,
                                   dof_func_array[level], debug_flag, trunc_factor,
                                   P_max_elmts, col_offd_S_to_A, &par_P, opt_icor);
            }
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 3) // direct interpolation (with separation of weights)
         {
            jxf_PAMGBuildDirInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                  par_S, coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag,
                                  trunc_factor, P_max_elmts,
                                  col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 4) // multipass interpolation
         {
            jxf_PAMGBuildMultipass(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                  par_S, coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag,
                                  trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 5) // multipass interpolation (with separation of weights)
         {
            jxf_PAMGBuildMultipass(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                  par_S, coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag,
                                  trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 6) // extended classical modified interpolation
         {
            jxf_PAMGBuildExtPIInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                    par_S, coarse_pnts_global, num_functions,
                                    dof_func_array[level], debug_flag,
                                    trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 7) // extended (if no common C neighbor) classical modified interpolation
         {
            jxf_PAMGBuildExtPICCInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                      par_S, coarse_pnts_global, num_functions,
                                      dof_func_array[level], debug_flag,
                                      trunc_factor, P_max_elmts, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 8) // standard interpolation
         {
            jxf_PAMGBuildStdInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                  par_S, coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag,
                                  trunc_factor, P_max_elmts, 0, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 9) // standard interpolation (with separation of weights)
         {
            jxf_PAMGBuildStdInterp(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                  par_S, coarse_pnts_global, num_functions,
                                  dof_func_array[level], debug_flag,
                                  trunc_factor, P_max_elmts, 1, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 100) /* 1pt interpolation */
         {
            jxf_PAMGBuildInterpOnePnt(jxf_hpCSRMatrixPar(A_array[level]), CF_marker_array[level],
                                     par_S, coarse_pnts_global, num_functions,
                                     dof_func_array[level], debug_flag, col_offd_S_to_A, &par_P);
            jxf_TFree(col_offd_S_to_A);
         }
         else if (interp_type == 20)
         {
            if (restri_type)
            {
               jxf_PAMGBuildInterpDinvRATrans(jxf_hpCSRMatrixPar(A_array[level]), R, &par_P);
            }
            jxf_TFree(col_offd_S_to_A);
         }
      }

      if (debug_flag == 1)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         jxf_printf("Proc = %d    Level = %d    Build Interp Time = %f\n", my_id, level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_interp += (jxf_time_getWallclockSeconds() - tmp_wall_time);
      }

      P_array[level] = par_P;

      if (restri_type)
      {
         R_array[level] = R;
         AIR_maxsize_ls[level] = R_max_size;
      }

      if (par_S)
      {
         jxf_ParCSRMatrixDestroy(par_S);
      }
      par_S = NULL;

      /*--------------------------------------------------------------------------------
       *   Build coarse-grid operator, A_array[level+1] by R*A*P
       *-------------------------------------------------------------------------------*/

      if (debug_flag == 1)
         wall_time = jxf_time_getWallclockSeconds();

      if (wall_time_option == 1)
         tmp_wall_time = jxf_time_getWallclockSeconds();

      A_H = jxf_hpInithpCSRMatrix();

      if (restri_type)
      {
         /* Use two matrix products to generate A_H */
         jxf_ParCSRMatrix *Q = jxf_ParMatmul(jxf_hpCSRMatrixPar(A_array[level]), P_array[level]);
         jxf_hpCSRMatrixPar(A_H) = jxf_ParMatmul(R_array[level], Q);
         jxf_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jxf_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jxf_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         jxf_ParCSRMatrixOwnsRowStarts(R_array[level]) = 0;
         if (num_procs > 1)
            jxf_MatvecCommPkgCreate(jxf_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jxf_ParCSRMatrixDestroy(Q);
      }
      else if (rap2)
      {
         /* Use two matrix products to generate A_H */
         jxf_ParCSRMatrix *Q = jxf_ParMatmul(jxf_hpCSRMatrixPar(A_array[level]), P_array[level]);
         jxf_hpCSRMatrixPar(A_H) = jxf_ParTMatmul(P_array[level], Q);
         jxf_hpCSRMatrixOwnsRowStarts(A_H) = 1;
         jxf_hpCSRMatrixOwnsColStarts(A_H) = 0;
         jxf_ParCSRMatrixOwnsColStarts(P_array[level]) = 0;
         if (num_procs > 1)
            jxf_MatvecCommPkgCreate(jxf_hpCSRMatrixPar(A_H));
         /* Delete AP */
         jxf_ParCSRMatrixDestroy(Q);
      }
      else if (spmt_rap_type == 1)
      {
         jxf_PAMGBuildCoarseOperatorKT(P_array[level], jxf_hpCSRMatrixPar(A_array[level]), P_array[level], keepTranspose, &(jxf_hpCSRMatrixPar(A_H)));
      }
      else if (spmt_rap_type == 2) /* Yue Xiaoqiang 2012/10/12 */
      {
         jxf_PAMGBuildCoarseOperatorOMP(P_array[level], jxf_hpCSRMatrixPar(A_array[level]), P_array[level], &(jxf_hpCSRMatrixPar(A_H)), opt_icor);
      }

/*Create hierarchy block of coarse matrix*/
#ifdef USING_HWLOC
      A_H = jxf_hpCreateMatrixLevelBlock(A_H);
#endif

      if (debug_flag == 1)
      {
         wall_time = jxf_time_getWallclockSeconds() - wall_time;
         jxf_printf("Proc = %d    Level = %d    Build Coarse Operator Time = %f\n", my_id, level, wall_time);
         fflush(NULL);
      }

      if (wall_time_option == 1)
      {
         wall_time_rap += (jxf_time_getWallclockSeconds() - tmp_wall_time);
      }

      ++level;

      jxf_ParCSRMatrixSetNumNonzeros(jxf_hpCSRMatrixPar(A_H));
      jxf_ParCSRMatrixSetDNumNonzeros(jxf_hpCSRMatrixPar(A_H));
      A_array[level] = A_H;

      if (coarse_size <= coarse_threshold_2)
      {

         coarse_size_last = coarse_size;
      }

      /* print coarser operator. */
      if (print_coarse_matrix == 1)
      {

         jxf_sprintf(FileNameCoaMat, "cmat_%d", level);
         jxf_hpCSRMatrixPrint(A_array[level], FileNameCoaMat);
      }

      /*-------------------------------------------------------------------------------
       *   Switch to CLJP when coarsening slows
       *                                            peghoty  2009/07/09
       *------------------------------------------------------------------------------*/

      size = ((JXF_Real)fine_size) * coarse_ratio; /* peghoty 2010/04/14 */
      if (coarsen_type > 0 && coarse_size >= (JXF_Int)size)
      {
         coarsen_type = 0;
      }

      /*--------------------------------------------------------------------------------
       *   How to stop the loop "while"
       *                                            peghoty  2009/07/09
       *--------------------------------------------------------------------------------*/
      if ((level == max_levels - 1) || (coarse_size <= coarse_threshold))
      {
         not_finished_coarsening = 0;

#if 0
         if (ai_measure_type == 1)
         {

           /*--------------------------------------------------------------------------
            *  Get the Strength Matrix 
            *-------------------------------------------------------------------------*/        
            jxf_PAMGCreateS(A_array[level],strong_threshold,max_row_sum,num_functions,dof_func_array[level],&par_S);
            col_offd_S_to_A = NULL;
            if (strong_threshold > S_commpkg_switch)
            {
               jxf_PAMGCreateSCommPkg(A_array[level], par_S, &col_offd_S_to_A);
            }

           /*--------------------------------------------------------------------------
            *  Get the AI Measure 
            *-------------------------------------------------------------------------*/        
            jxf_PAMGMeasureAI(par_S, A_array[level], debug_flag, &AI_measure);

           /*--------------------------------------------
            *  store the AI array
            *-------------------------------------------*/ 
            AI_measure_array[level] = AI_measure;

         }

         /*-------------------------------------------------------------------------
          *  compute ai-information for all-pnts and C-pnts!
          *-------------------------------------------------------------------------*/ 
         local_num_vars = jxf_CSRMatrixNumRows(jxf_ParCSRMatrixDiag(A_array[level]));

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

      /*jxf_printf("level = %d, not_finished_coarsening = %d \n", level, not_finished_coarsening);
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

      jxf_MPI_Allreduce(&num_ai_local_0, &num_ai_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_ai_local_valid_0, &num_ai_global_valid, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_vars_local_0, &num_vars_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_ai_c_local_0, &num_ai_c_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_local_0, &mai_global, 1, JXF_MPI_REAL, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_local_valid_0, &mai_global_valid, 1, JXF_MPI_REAL, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_c_local_0, &mai_c_global, 1, JXF_MPI_REAL, MPI_SUM, comm);

      if (my_id == 0)
      {
         JXF_Real num_ai = num_ai_global;
         JXF_Real num_ai_valid = num_ai_global_valid;
         JXF_Real num_var = num_vars_global;
         JXF_Real num_ai_c = num_ai_c_global;
         jxf_printf(" \n");
         jxf_printf("Level 0: \n");
         jxf_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jxf_printf(" - num_ai_ratio = %f\n", num_ai / num_var);
         jxf_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid / num_var);
         if (num_ai)
            jxf_printf(" - num_ai_c_ratio = %f\n", num_ai_c / num_ai);
         if (num_ai_valid)
            jxf_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c / num_ai_valid);
         if (mai_global > 0.0)
            jxf_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid / mai_global);
         if (mai_global > 0.0)
            jxf_printf(" - measure_ai_c_ratio = %f\n", mai_c_global / mai_global);
         if (mai_global_valid > 0.0)
            jxf_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global / mai_global_valid);
      }

      jxf_MPI_Allreduce(&num_ai_local, &num_ai_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_ai_local_valid, &num_ai_global_valid, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_vars_local, &num_vars_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&num_ai_c_local, &num_ai_c_global, 1, JXF_MPI_INT, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_local, &mai_global, 1, JXF_MPI_REAL, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_local_valid, &mai_global_valid, 1, JXF_MPI_REAL, MPI_SUM, comm);
      jxf_MPI_Allreduce(&mai_c_local, &mai_c_global, 1, JXF_MPI_REAL, MPI_SUM, comm);

      if (my_id == 0)
      {
         JXF_Real num_ai = num_ai_global;
         JXF_Real num_ai_valid = num_ai_global_valid;
         JXF_Real num_var = num_vars_global;
         JXF_Real num_ai_c = num_ai_c_global;
         jxf_printf(" \n");
         jxf_printf("All levels (except coarsest level): \n");
         jxf_printf(" - num_vars = %d, num_ai = %d, num_ai_valid = %d \n", num_vars_global, num_ai_global, num_ai_global_valid);
         jxf_printf(" - num_ai_ratio = %f\n", num_ai / num_var);
         jxf_printf(" - num_ai_ratio_valid = %f\n", num_ai_valid / num_var);
         if (num_ai > 0)
            jxf_printf(" - num_ai_c_ratio = %f\n", num_ai_c / num_ai);
         if (num_ai_valid > 0)
            jxf_printf(" - num_ai_c_ratio_valid = %f\n", num_ai_c / num_ai_valid);
         if (mai_global > 0.0)
            jxf_printf(" - measure_ai_ratio_valid = %f\n", mai_global_valid / mai_global);
         if (mai_global > 0.0)
            jxf_printf(" - measure_ai_c_ratio = %f\n", mai_c_global / mai_global);
         if (mai_global_valid > 0.0)
            jxf_printf(" - measure_ai_c_ratio_valid = %f\n", mai_c_global / mai_global_valid);
         jxf_printf("================================================================== \n");
      }
   }

   if (level > 0)
   {
      F_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jxf_hpCSRMatrixRowStarts(A_array[level]));
      jxf_ParVectorInitialize(F_array[level]);
      jxf_ParVectorSetPartitioningOwner(F_array[level], 0);

      U_array[level] = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[level]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[level]),
                                          jxf_hpCSRMatrixRowStarts(A_array[level]));
      jxf_ParVectorInitialize(U_array[level]);
      jxf_ParVectorSetPartitioningOwner(U_array[level], 0);
   }

   /*-----------------------------------------------------------------------
    * enter all the stuff created, A[level], P[level], CF_marker[level],
    * for levels 1 through coarsest, into amg_data data structure
    *-----------------------------------------------------------------------*/

   num_levels = level + 1;
   jxf_ParAMGDataNumLevels(amg_data) = num_levels;
   if (jxf_ParAMGDataSmoothNumLevels(amg_data) > level)
   {
      jxf_ParAMGDataSmoothNumLevels(amg_data) = level;
   }
   smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data);

   if (grid_relax_type[1] == 8 || grid_relax_type[2] == 8 || grid_relax_type[3] == 8 ||
       grid_relax_type[1] == 13 || grid_relax_type[2] == 13 || grid_relax_type[3] == 13 ||
       grid_relax_type[1] == 14 || grid_relax_type[2] == 14 || grid_relax_type[3] == 14)
   {
      l1_norms = jxf_CTAlloc(JXF_Real *, num_levels);
      jxf_ParAMGDataL1Norms(amg_data) = l1_norms;
   }

   for (j = 0; j < num_levels; j++)
   {
      if (num_threads == 1)
      {
         if (j < num_levels - 1 && (grid_relax_type[1] == 8 || grid_relax_type[1] == 13 || grid_relax_type[1] == 14 ||
                                    grid_relax_type[2] == 8 || grid_relax_type[2] == 13 || grid_relax_type[2] == 14))
         {
            if (relax_order)
               jxf_ParCSRComputeL1Norms(jxf_hpCSRMatrixPar(A_array[j]), 4, CF_marker_array[j], &l1_norms[j]);
            else
               jxf_ParCSRComputeL1Norms(jxf_hpCSRMatrixPar(A_array[j]), 4, NULL, &l1_norms[j]);
         }
         else if ((grid_relax_type[3] == 8 || grid_relax_type[3] == 13 || grid_relax_type[3] == 14) && j == num_levels - 1)
         {
            jxf_ParCSRComputeL1Norms(jxf_hpCSRMatrixPar(A_array[j]), 4, NULL, &l1_norms[j]);
         }
      }
      else
      {
         if (j < num_levels - 1 && (grid_relax_type[1] == 8 || grid_relax_type[1] == 13 || grid_relax_type[1] == 14 ||
                                    grid_relax_type[2] == 8 || grid_relax_type[2] == 13 || grid_relax_type[2] == 14))
         {
            if (relax_order)
               jxf_ParCSRComputeL1NormsThreads(jxf_hpCSRMatrixPar(A_array[j]), 4, num_threads, CF_marker_array[j], &l1_norms[j]);
            else
               jxf_ParCSRComputeL1NormsThreads(jxf_hpCSRMatrixPar(A_array[j]), 4, num_threads, NULL, &l1_norms[j]);
         }
         else if ((grid_relax_type[3] == 8 || grid_relax_type[3] == 13 || grid_relax_type[3] == 14) && j == num_levels - 1)
         {
            jxf_ParCSRComputeL1NormsThreads(jxf_hpCSRMatrixPar(A_array[j]), 4, num_threads, NULL, &l1_norms[j]);
         }
      }
   }

   for (j = 0; j < smooth_num_levels; j++)
   {
      if (smooth_type == 9 || smooth_type == 19) // Euclid smoothers
      {
         JXF_EuclidCreate(comm, &smoother[j]);
         if (euclidfile)
         {
            JXF_EuclidSetParamsFromFile(smoother[j], euclidfile);
         }
         JXF_EuclidSetLevel(smoother[j], eu_level);
         if (eu_bj)
         {
            JXF_EuclidSetBJ(smoother[j], eu_bj);
         }
         if (eu_sparse_A)
         {
            JXF_EuclidSetSparseA(smoother[j], eu_sparse_A);
         }
         JXF_EuclidSetup(smoother[j], (JXF_hpCSRMatrix)A_array[j]);
      }
      else if (smooth_type == 5 || smooth_type == 15) // ILU smoothers
      {
         JXF_ILUCreate(&smoother[j]);
         if (eu_bj)
         {
            JXF_ILUSetType(smoother[j], 0);
         }
         else
         {
            JXF_ILUSetType(smoother[j], 30);
         }
         JXF_ILUSetMaxIter(smoother[j], 1);
         /* set tol to zero since we are doing just 1 iteration */
         JXF_ILUSetTol(smoother[j], 0.);
         JXF_ILUSetLogging(smoother[j], 0);
         JXF_ILUSetPrintLevel(smoother[j], 0);
         JXF_ILUSetLevelOfFill(smoother[j], eu_level);
         JXF_ILUSetup(smoother[j], (JXF_ParCSRMatrix)A_array[j]);
      }
   }

   if (amg_logging > 1)
   {
      Residual_array = jxf_ParVectorCreate(jxf_hpCSRMatrixComm(A_array[0]),
                                          jxf_hpCSRMatrixGlobalNumRows(A_array[0]),
                                          jxf_hpCSRMatrixRowStarts(A_array[0]));
      jxf_ParVectorInitialize(Residual_array);
      jxf_ParVectorSetPartitioningOwner(Residual_array, 0);
      jxf_ParAMGDataResidual(amg_data) = Residual_array;
   }
   else
   {
      jxf_ParAMGDataResidual(amg_data) = NULL;
   }

   jxf_TFree(opt_icor); /* Yue Xiaoqiang 2012/10/12 */
   if (Jxf_Pmarkers_global_rap)
   {
      jxf_TFree(Jxf_Pmarkers_global_rap); /* Yue Xiaoqiang 2012/10/17 */
   }
   if (Jxf_Pmarkers_global_interp)
   {
      jxf_TFree(Jxf_Pmarkers_global_interp); /* Yue Xiaoqiang 2012/10/17 */
   }

   /*--------------------------------------------------------------
    *    Print some stuff
    *-------------------------------------------------------------*/

   // if (amg_print_level == 1 || amg_print_level == 3)
   /*change by mrz 2022/10/26*/
   if (amg_print_level > 1)
   {
      /* Write the SETUP parameters */
      jxf_hpPAMGSetupStatus(amg_data, hp_matrix);
   }

   if (wall_time_option == 1)
   {
      jxf_printf("\n\nProc = %d, Coarsen Time = %f\n", my_id, wall_time_coarsen);
      jxf_printf("Proc = %d, Build Coarse Operator Time = %f\n", my_id, wall_time_rap);
      jxf_printf("Proc = %d, Build Interp Time = %f\n\n", my_id, wall_time_interp);
   }

   if (my_id == 0 && amg_print_level > 1)
   {
      /* Write the SOLVE parameters */
      jxf_PAMGSolveStatus(amg_data); /* Yue Xiaoqiang 2014/04/12 */
   }

   return (Setup_err_flag);
}
