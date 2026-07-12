//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_amg.c
 *  Date: 2011/09/03
 */ 

#include "jxf_pamg.h"


/*!
 * \fn JXF_Int JXF_PAMGCreate
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGCreate( JXF_Solver *solver)
{
   *solver = (JXF_Solver) jxf_PAMGCreate( ) ;
   if (!solver)
   {
      jxf_error_in_arg(1);
   }

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int JXF_PAMGDestroy
 * \date 2011/09/03
 */
JXF_Int 
JXF_PAMGDestroy( JXF_Solver solver)
{
   return( jxf_hpPAMGDestroy( (void *) solver) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetup
 * \date 2011/09/03
 */
JXF_Int 
JXF_PAMGSetup( JXF_Solver        solver, 
              JXF_hpCSRMatrix  hp_matrix)
{
   return( jxf_hpPAMGSetup( (void *) solver, (jxf_hpCSRMatrix *) hp_matrix));
}

/*!
 * \fn JXF_Int JXF_PAMGSolve
 * \date 2011/09/03
 */
JXF_Int 
JXF_PAMGSolve( JXF_Solver       solver,
              JXF_hpCSRMatrix hp_matrix,
              JXF_ParVector    par_rhs,
              JXF_ParVector    par_app)
{
   return( jxf_hpPAMGSolve( (void *) solver,
                         (jxf_hpCSRMatrix *) hp_matrix,
                         (jxf_ParVector *) par_rhs,
                         (jxf_ParVector *) par_app));
}

/*!
 * \fn JXF_Int JXF_PAMGPrecond
 * \author peghoty
 * \date 2011/09/04
 */
JXF_Int 
JXF_PAMGPrecond( JXF_Solver       solver,
                JXF_hpCSRMatrix par_matrix,
                JXF_ParVector    par_rhs,
                JXF_ParVector    par_app)
{
   return( jxf_PAMGPrecond( (void *) solver,
                           (jxf_hpCSRMatrix *) par_matrix,
                           (jxf_ParVector *) par_rhs,
                           (jxf_ParVector *) par_app));
}

/*!
 * \fn void *jxf_PAMGCreate
 * \date 2011/09/03
 */

void *
jxf_PAMGCreate()
{
   jxf_ParAMGData  *amg_data;

   /* setup params */
   JXF_Int      max_levels;
   JXF_Real   strong_threshold;
   JXF_Real   max_row_sum;
   JXF_Real   trunc_factor;
   JXF_Real   jacobi_trunc_threshold;
   JXF_Real   S_commpkg_switch;
   JXF_Real   CR_rate;
   JXF_Real   CR_strong_th;
   JXF_Real   AIR_strong_th;
   JXF_Int      interp_type;
   JXF_Int      coarsen_type;
   JXF_Int      measure_type;
   JXF_Int      setup_type;
   JXF_Int      P_max_elmts;
   JXF_Int      num_functions;
   JXF_Int      nodal, nodal_diag;
   JXF_Int      num_paths;
   JXF_Int      agg_num_levels;
   JXF_Int      agg_interp_type;
   JXF_Int      agg_P_max_elmts;
   JXF_Int      agg_P12_max_elmts;
   JXF_Real   agg_trunc_factor;
   JXF_Real   agg_P12_trunc_factor;
   JXF_Int      post_interp_type;
   JXF_Int      num_CR_relax_steps;
   JXF_Int      IS_type;
   JXF_Int      CR_use_CG;
   JXF_Int      cgc_its;
   JXF_Int      spmt_rap_type;          // Yue Xiaoqiang 2012/10/22
   JXF_Int      wall_time_option;       // Yue Xiaoqiang 2015/09/30
   JXF_Int      ai_measure_type;        // Yue Xiaoqiang 2014/02/26
   JXF_Int      ai_relax_type;          // Yue Xiaoqiang 2014/07/06
   JXF_Int      measure_type_rlx;       // peghoty 2010/05/29
   JXF_Int      number_syn_rlx;         // peghoty 2010/05/29
   JXF_Real   measure_threshold_rlx;  // peghoty 2010/05/29   

   /* solve params */
   JXF_Int      min_iter;
   JXF_Int      max_iter;
   JXF_Int      cycle_type;    
 
   JXF_Real   tol;
   JXF_Real   rhsnrm_threshold; /* peghoty, 2011/09/28 */

   JXF_Int      num_sweeps;  
   JXF_Int      relax_type;   
   JXF_Int      relax_order;   
   JXF_Real   relax_wt;
   JXF_Real   outer_wt;
   JXF_Int      smooth_type;
   JXF_Int      smooth_num_levels;
   JXF_Int      smooth_num_sweeps;

   JXF_Int      variant, overlap, domain_type, schwarz_use_nonsymm;
   JXF_Real   schwarz_rlx_weight;
   JXF_Int      level, sym;
   JXF_Int      eu_level, eu_bj;
   JXF_Int      max_nz_per_row;
   JXF_Real   thresh, filter;
   JXF_Real   drop_tol;
   JXF_Real   eu_sparse_A;
   char    *euclidfile;

   JXF_Int block_mode;
   
   /* log info */
   JXF_Int      num_iterations;
   JXF_Int      cum_num_iterations;

   /* output params */
   JXF_Int      print_level;
   JXF_Int      print_coarse_matrix;
   JXF_Int      logging;
   char     log_file_name[256];
   char     plot_file_name[251];
   JXF_Int      debug_flag;
   
   JXF_Int      conv_criteria;     /* peghoty  2009/07/25 */
   JXF_Int      coarsestsolverid;  /* peghoty  2009/07/27 */
   JXF_Int      coarse_threshold;  /* peghoty  2010/04/14 */
   JXF_Real   coarse_ratio;      /* peghoty  2010/04/14 */
   JXF_Real   convfac_threshold;      /* peghoty  2010/04/14 */

  /*-----------------------------------------------------------------------
   * Setup default values for parameters
   *-----------------------------------------------------------------------*/

   /* setup params */
   max_levels = 25;
   strong_threshold = 0.25;
   max_row_sum = 0.9;
   trunc_factor = 0.0;
   jacobi_trunc_threshold = 0.01;
   S_commpkg_switch = 1.0;
   AIR_strong_th = 0.1;
   interp_type = 0;
   coarsen_type = 6;
   measure_type = 0;
   setup_type = 1;
   P_max_elmts = 0;
   num_functions = 1;
   nodal = 0;
   nodal_diag = 0;
   num_paths = 1;
   agg_num_levels = 0;
   agg_interp_type = 4;
   agg_P_max_elmts = 0;
   agg_P12_max_elmts = 0;
   agg_trunc_factor = 0.0;
   agg_P12_trunc_factor = 0.0;
   post_interp_type = 0;
   num_CR_relax_steps = 2;
   CR_rate = 0.7;
   CR_strong_th = 0;
   IS_type = 1;
   CR_use_CG = 0;
   cgc_its = 1;
   spmt_rap_type = 1;            // Yue Xiaoqiang 2012/10/22
   wall_time_option = 0;         // Yue Xiaoqiang 2015/09/30
   ai_measure_type = 0;          // Yue Xiaoqiang 2014/02/26
   ai_relax_type = 0;            // Yue Xiaoqiang 2014/07/06
   measure_type_rlx = 1;         // peghoty 2010/05/29
   number_syn_rlx = 2;           // peghoty 2010/05/29
   measure_threshold_rlx = 0.8;  // peghoty 2010/05/29

   variant = 0;
   overlap = 1;
   domain_type = 2;
   schwarz_rlx_weight = 1.0;
   smooth_num_sweeps = 1;
   smooth_num_levels = 0;
   smooth_type = 6;
   schwarz_use_nonsymm = 0;
   
   level = 1;
   sym = 0;
   thresh = 0.1;
   filter = 0.05;
   drop_tol = 0.0001;
   max_nz_per_row = 20;
   euclidfile = NULL;
   eu_level = 0;
   eu_sparse_A = 0.0;
   eu_bj = 0;

   /* solve params */
   min_iter  = 0;
   max_iter  = 20;
   cycle_type = 1;
   tol = 1.0e-7;
   rhsnrm_threshold = 0.0;

   num_sweeps = 1;
   relax_type = 3;
   relax_order = 1;
   relax_wt = 1.0;
   outer_wt = 1.0;

   block_mode = 0;

   /* log info */
   num_iterations = 0;
   cum_num_iterations = 0;

   /* output params */
   print_level = 0;
   print_coarse_matrix = 0;
   logging = 0;
   jxf_sprintf(log_file_name, "%s", "amg.out.log");
   jxf_sprintf(plot_file_name, "%s", "amg.plot.log"); /* peghoty 2009/08/25 */
   debug_flag = 0;

   
   conv_criteria    = 0;    /* peghoty 2009/07/25 */
   coarsestsolverid = 9;    /* peghoty 2009/07/27 */
   coarse_threshold = 9;    /* peghoty 2010/04/14 */
   coarse_ratio     = 0.75; /* peghoty 2010/04/14 */    

   convfac_threshold = 1.0e4;

  /*-----------------------------------------------------------------------
   * Create the jxf_ParAMGData structure and return
   *-----------------------------------------------------------------------*/

   amg_data = jxf_CTAlloc(jxf_ParAMGData, 1);

   jxf_ParAMGDataUserCoarseRelaxType(amg_data) = 9;
   jxf_ParAMGDataUserRelaxType(amg_data) = -1;
   jxf_PAMGSetMaxLevels(amg_data, max_levels);
   jxf_PAMGSetStrongThreshold(amg_data, strong_threshold);
   jxf_PAMGSetMaxRowSum(amg_data, max_row_sum);
   jxf_PAMGSetTruncFactor(amg_data, trunc_factor);
   jxf_PAMGSetJacobiTruncThreshold(amg_data, jacobi_trunc_threshold);
   jxf_PAMGSetSCommPkgSwitch(amg_data, S_commpkg_switch);
   jxf_PAMGSetInterpType(amg_data, interp_type);
   jxf_PAMGSetMeasureType(amg_data, measure_type);
   jxf_PAMGSetCoarsenType(amg_data, coarsen_type);
   jxf_PAMGSetSetupType(amg_data, setup_type);
   jxf_PAMGSetPMaxElmts(amg_data, P_max_elmts);
   jxf_PAMGSetNumFunctions(amg_data, num_functions);
   jxf_PAMGSetNodal(amg_data, nodal);
   jxf_PAMGSetNodal(amg_data, nodal_diag);
   jxf_PAMGSetNumPaths(amg_data, num_paths);
   jxf_PAMGSetAggNumLevels(amg_data, agg_num_levels);
   jxf_PAMGSetAggInterpType(amg_data, agg_interp_type);
   jxf_PAMGSetAggPMaxElmts(amg_data, agg_P_max_elmts);
   jxf_PAMGSetAggP12MaxElmts(amg_data, agg_P12_max_elmts);
   jxf_PAMGSetAggTruncFactor(amg_data, agg_trunc_factor);
   jxf_PAMGSetAggP12TruncFactor(amg_data, agg_P12_trunc_factor);
   jxf_PAMGSetPostInterpType(amg_data, post_interp_type);
   jxf_PAMGSetNumCRRelaxSteps(amg_data, num_CR_relax_steps);
   jxf_PAMGSetCRRate(amg_data, CR_rate);
   jxf_PAMGSetCRStrongTh(amg_data, CR_strong_th);
   jxf_PAMGSetISType(amg_data, IS_type);
   jxf_PAMGSetCRUseCG(amg_data, CR_use_CG);
   jxf_PAMGSetCGCIts(amg_data, cgc_its);
   
   jxf_PAMGSetSpMtRapType(amg_data, spmt_rap_type); // Yue Xiaoqiang 2012/10/13
   jxf_PAMGSetWallTimeOption(amg_data, wall_time_option); // Yue Xiaoqiang 2015/09/30
   jxf_PAMGSetAIMeasureType(amg_data, ai_measure_type); // Yue Xiaoqiang 2014/02/26
   jxf_PAMGSetAIRelaxType(amg_data, ai_relax_type); // Yue Xiaoqiang 2014/07/06
   
   jxf_PAMGSetRelaxedCoarsenMeasureType(amg_data,measure_type_rlx);           // peghoty 2010/05/29
   jxf_PAMGSetRelaxedCoarsenNumberSyn(amg_data,number_syn_rlx);               // peghoty 2010/05/29
   jxf_PAMGSetRelaxedCoarsenMeasureThreshold(amg_data,measure_threshold_rlx); // peghoty 2010/05/29

   jxf_PAMGSetVariant(amg_data, variant);
   jxf_PAMGSetOverlap(amg_data, overlap);
   jxf_PAMGSetSchwarzRlxWeight(amg_data, schwarz_rlx_weight);
   jxf_PAMGSetSchwarzUseNonSymm(amg_data, schwarz_use_nonsymm);
   jxf_PAMGSetDomainType(amg_data, domain_type);
   jxf_PAMGSetSym(amg_data, sym);
   jxf_PAMGSetLevel(amg_data, level);
   jxf_PAMGSetThreshold(amg_data, thresh);
   jxf_PAMGSetFilter(amg_data, filter);
   jxf_PAMGSetDropTol(amg_data, drop_tol);
   jxf_PAMGSetMaxNzPerRow(amg_data, max_nz_per_row);
   jxf_PAMGSetEuclidFile(amg_data, euclidfile);
   jxf_PAMGSetEuLevel(amg_data, eu_level);
   jxf_PAMGSetEuSparseA(amg_data, eu_sparse_A);
   jxf_PAMGSetEuBJ(amg_data, eu_bj);

   jxf_PAMGSetMinIter(amg_data, min_iter);
   jxf_PAMGSetMaxIter(amg_data, max_iter);
   jxf_PAMGSetCycleType(amg_data, cycle_type);
   jxf_PAMGSetTol(amg_data, tol); 
   jxf_PAMGSetRhsNrmThreshold(amg_data, rhsnrm_threshold);
   jxf_PAMGSetNumSweeps(amg_data, num_sweeps);
   jxf_PAMGSetRelaxType(amg_data, relax_type);
   jxf_PAMGSetRelaxOrder(amg_data, relax_order);
   jxf_PAMGSetRelaxWt(amg_data, relax_wt);
   jxf_PAMGSetOuterWt(amg_data, outer_wt);
   jxf_PAMGSetSmoothType(amg_data, smooth_type);
   jxf_PAMGSetSmoothNumLevels(amg_data, smooth_num_levels);
   jxf_PAMGSetSmoothNumSweeps(amg_data, smooth_num_sweeps);

   jxf_PAMGSetNumIterations(amg_data, num_iterations);
#ifdef CUMNUMIT
   jxf_ParAMGDataCumNumIterations(amg_data) = cum_num_iterations;
#endif
   jxf_PAMGSetPrintLevel(amg_data, print_level);
   jxf_PAMGSetPrintCoarseMatrix(amg_data, print_coarse_matrix);
   jxf_PAMGSetLogging(amg_data, logging);
   jxf_PAMGSetPrintFileName(amg_data, log_file_name); 
   jxf_PAMGSetDebugFlag(amg_data, debug_flag);
   
   jxf_PAMGSetConvCriteria(amg_data, conv_criteria);        /* peghoty 2009/07/25 */
   jxf_PAMGSetCoarsestSolverID(amg_data, coarsestsolverid); /* peghoty 2009/07/27 */
   jxf_PAMGSetCoarseThreshold(amg_data, coarse_threshold ); /* peghoty 2010/04/14 */
   jxf_PAMGSetCoarseRatio(amg_data, coarse_ratio );         /* peghoty 2010/04/14 */

   jxf_PAMGSetConvFacThreshold(amg_data, convfac_threshold);

   jxf_PAMGSetRestriction(amg_data, 0);
   jxf_PAMGSetAIRStrongTh(amg_data, AIR_strong_th);


   jxf_PAMGSetGSMG(amg_data, 0);
   jxf_PAMGSetNumSamples(amg_data, 0);
   
   jxf_hpAMGDataAArray(amg_data) = NULL;
   jxf_ParAMGDataPArray(amg_data) = NULL;
   jxf_ParAMGDataRArray(amg_data) = NULL;
   jxf_ParAMGDataAIMeasureArray(amg_data) = NULL;
   jxf_ParAMGDataCFMarkerArray(amg_data) = NULL;
   jxf_ParAMGDataVtemp(amg_data)  = NULL;
   jxf_ParAMGDataRtemp(amg_data)  = NULL;
   jxf_ParAMGDataPtemp(amg_data)  = NULL;
   jxf_ParAMGDataZtemp(amg_data)  = NULL;
   jxf_ParAMGDataFArray(amg_data) = NULL;
   jxf_ParAMGDataUArray(amg_data) = NULL;
   jxf_ParAMGDataDofFunc(amg_data) = NULL;
   jxf_ParAMGDataDofFuncArray(amg_data) = NULL;
   jxf_ParAMGDataDofPointArray(amg_data) = NULL;
   jxf_ParAMGDataDofPointArray(amg_data) = NULL;
   jxf_ParAMGDataPointDofMapArray(amg_data) = NULL;
   jxf_ParAMGDataSmoother(amg_data) = NULL;
   jxf_ParAMGDataL1Norms(amg_data) = NULL;
  
   jxf_ParAMGDataABlockArray(amg_data) = NULL;
   jxf_ParAMGDataPBlockArray(amg_data) = NULL;
   jxf_ParAMGDataRBlockArray(amg_data) = NULL;

   /* this can not be set by the user currently */
   jxf_ParAMGDataBlockMode(amg_data) = block_mode;

   /* BM Oct 22, 2006 */
   jxf_ParAMGDataPlotGrids(amg_data) = 0;
   jxf_PAMGSetPlotFileName(amg_data, plot_file_name);

   /* BM Oct 17, 2006 */
   jxf_ParAMGDataCoordDim(amg_data) = 0;
   jxf_ParAMGDataCoordinates(amg_data) = NULL;

   jxf_ParAMGDataRAP2(amg_data) = 0;
   jxf_ParAMGDataKeepTranspose(amg_data) = 0;

     
   return (void *) amg_data;
}

JXF_Int
jxf_hpPAMGDestroy( void *data )
{
   jxf_ParAMGData  *amg_data = data;
   JXF_Int num_levels = jxf_ParAMGDataNumLevels(amg_data);
   JXF_Int i;

   if (jxf_ParAMGDataNumGridSweeps(amg_data))
   {
      jxf_TFree (jxf_ParAMGDataNumGridSweeps(amg_data));
      jxf_ParAMGDataNumGridSweeps(amg_data) = NULL; 
   }
   if (jxf_ParAMGDataGridRelaxType(amg_data))
   {
      jxf_TFree (jxf_ParAMGDataGridRelaxType(amg_data));
      jxf_ParAMGDataGridRelaxType(amg_data) = NULL; 
   }
   if (jxf_ParAMGDataRelaxWeight(amg_data))
   {
      jxf_TFree (jxf_ParAMGDataRelaxWeight(amg_data));
      jxf_ParAMGDataRelaxWeight(amg_data) = NULL; 
   }
   if (jxf_ParAMGDataOmega(amg_data))
   {
      jxf_TFree (jxf_ParAMGDataOmega(amg_data));
      jxf_ParAMGDataOmega(amg_data) = NULL; 
   }
   if (jxf_ParAMGDataDofFunc(amg_data))
   {
      jxf_TFree (jxf_ParAMGDataDofFunc(amg_data));
      jxf_ParAMGDataDofFunc(amg_data) = NULL; 
   }
   if (jxf_ParAMGDataGridRelaxPoints(amg_data))
   {
      for (i = 0; i < 4; i ++)
      {
   	 jxf_TFree (jxf_ParAMGDataGridRelaxPoints(amg_data)[i]);
      }
      jxf_TFree (jxf_ParAMGDataGridRelaxPoints(amg_data));
      jxf_ParAMGDataGridRelaxPoints(amg_data) = NULL; 
   }
   for (i = 1; i < num_levels; i ++)
   {
	jxf_ParVectorDestroy(jxf_ParAMGDataFArray(amg_data)[i]);
	jxf_ParVectorDestroy(jxf_ParAMGDataUArray(amg_data)[i]);

        if (jxf_hpAMGDataAArray(amg_data)[i])
        {
           jxf_hpCSRMatrixDestroy(jxf_hpAMGDataAArray(amg_data)[i]);
        }

        if (jxf_ParAMGDataPArray(amg_data)[i-1])
        {
            jxf_ParCSRMatrixDestroy(jxf_ParAMGDataPArray(amg_data)[i-1]);
        }

        if (jxf_ParAMGDataRestriction(amg_data))
        {
           if (jxf_ParAMGDataRArray(amg_data)[i-1])
           {
              jxf_ParCSRMatrixDestroy(jxf_ParAMGDataRArray(amg_data)[i-1]);
           }
           jxf_TFree(jxf_ParAMGDataAIRMaxSizeLS(amg_data));
        }

        if (jxf_ParAMGDataAIMeasureArray(amg_data)[i-1])
        {
            jxf_TFree(jxf_ParAMGDataAIMeasureArray(amg_data)[i-1]);
        }

	jxf_TFree(jxf_ParAMGDataCFMarkerArray(amg_data)[i-1]);
   }

   if (jxf_ParAMGDataL1Norms(amg_data))
   {
      for (i = 0; i < num_levels; i ++)
      {
         if (jxf_ParAMGDataL1Norms(amg_data)[i])
         {
            jxf_TFree(jxf_ParAMGDataL1Norms(amg_data)[i]);
         }
      }
      jxf_TFree(jxf_ParAMGDataL1Norms(amg_data));
   }

   /* see comments in par_coarsen.c regarding special case for CF_marker */
   if (num_levels == 1)
   {
      jxf_TFree(jxf_ParAMGDataAIMeasureArray(amg_data)[0]);
      jxf_TFree(jxf_ParAMGDataCFMarkerArray(amg_data)[0]);
   }
   jxf_ParVectorDestroy(jxf_ParAMGDataVtemp(amg_data));
   jxf_TFree(jxf_ParAMGDataFArray(amg_data));
   jxf_TFree(jxf_ParAMGDataUArray(amg_data));
   jxf_TFree(jxf_hpAMGDataAArray(amg_data));
   jxf_TFree(jxf_ParAMGDataPArray(amg_data));
   if (jxf_ParAMGDataRestriction(amg_data)) jxf_TFree(jxf_ParAMGDataRArray(amg_data));
   jxf_TFree(jxf_ParAMGDataAIMeasureArray(amg_data));
   jxf_TFree(jxf_ParAMGDataCFMarkerArray(amg_data));
   if (jxf_ParAMGDataRtemp(amg_data))
   {
      jxf_ParVectorDestroy(jxf_ParAMGDataRtemp(amg_data));
   }
   if (jxf_ParAMGDataPtemp(amg_data))
   {
      jxf_ParVectorDestroy(jxf_ParAMGDataPtemp(amg_data));
   }
   if (jxf_ParAMGDataZtemp(amg_data))
   {
      jxf_ParVectorDestroy(jxf_ParAMGDataZtemp(amg_data));
   }

   if (jxf_ParAMGDataDofFuncArray(amg_data))
   {
      for (i = 1; i < num_levels; i ++)
      {
         jxf_TFree(jxf_ParAMGDataDofFuncArray(amg_data)[i]);
      }
      jxf_TFree(jxf_ParAMGDataDofFuncArray(amg_data));
      jxf_ParAMGDataDofFuncArray(amg_data) = NULL;
   }
   if (jxf_ParAMGDataRestriction(amg_data))
   {
      jxf_TFree(jxf_ParAMGDataRArray(amg_data));
      jxf_ParAMGDataRArray(amg_data) = NULL;
   }
   if (jxf_ParAMGDataDofPointArray(amg_data))
   {
      for (i = 0; i < num_levels; i ++)
      {
         jxf_TFree(jxf_ParAMGDataDofPointArray(amg_data)[i]);
      }
      jxf_TFree(jxf_ParAMGDataDofPointArray(amg_data));
      jxf_ParAMGDataDofPointArray(amg_data) = NULL;
   }
   if (jxf_ParAMGDataPointDofMapArray(amg_data))
   {
      for (i = 0; i < num_levels; i ++)
      {
         jxf_TFree(jxf_ParAMGDataPointDofMapArray(amg_data)[i]);
      }
      jxf_TFree(jxf_ParAMGDataPointDofMapArray(amg_data));
      jxf_ParAMGDataPointDofMapArray(amg_data) = NULL;
   }
   if ( jxf_ParAMGDataResidual(amg_data) ) 
   {
      jxf_ParVectorDestroy( jxf_ParAMGDataResidual(amg_data) );
      jxf_ParAMGDataResidual(amg_data) = NULL;
   }

   jxf_TFree(amg_data);
   return jxf_error_flag;
}


/*!
 * \fn JXF_Int JXF_PAMGSetXXXX
 * \brief  Routines to set the setup phase parameters begin with "JXF_"
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetRestriction( JXF_Solver solver, JXF_Int restr_par )
{
   return( jxf_PAMGSetRestriction( (void *) solver, restr_par ) );
}


JXF_Int
JXF_PAMGSetAIRStrongTh( JXF_Solver solver, JXF_Real AIR_strong_th )
{
   return( jxf_PAMGSetAIRStrongTh( (void *) solver, AIR_strong_th ) );
}


JXF_Int
JXF_PAMGSetMaxLevels( JXF_Solver solver, JXF_Int  max_levels )
{
   return( jxf_PAMGSetMaxLevels( (void *) solver, max_levels ) );
}

JXF_Int
JXF_PAMGGetMaxLevels( JXF_Solver solver, JXF_Int *max_levels )
{
   return( jxf_PAMGGetMaxLevels( (void *) solver, max_levels ) );
}

JXF_Int
JXF_PAMGSetStrongThreshold( JXF_Solver solver, JXF_Real strong_threshold  )
{
   return( jxf_PAMGSetStrongThreshold( (void *) solver, strong_threshold ) );
}

JXF_Int
JXF_PAMGGetStrongThreshold( JXF_Solver solver, JXF_Real *strong_threshold )
{
   return( jxf_PAMGGetStrongThreshold( (void *) solver, strong_threshold ) );
}

JXF_Int
JXF_PAMGSetMaxRowSum( JXF_Solver solver, JXF_Real max_row_sum )
{
   return( jxf_PAMGSetMaxRowSum( (void *) solver, max_row_sum ) );
}

JXF_Int
JXF_PAMGGetMaxRowSum( JXF_Solver solver, JXF_Real *max_row_sum )
{
   return( jxf_PAMGGetMaxRowSum( (void *) solver, max_row_sum ) );
}

JXF_Int
JXF_PAMGSetTruncFactor( JXF_Solver solver, JXF_Real trunc_factor )
{
   return( jxf_PAMGSetTruncFactor( (void *) solver, trunc_factor ) );
}

JXF_Int
JXF_PAMGGetTruncFactor( JXF_Solver solver, JXF_Real *trunc_factor )
{
   return( jxf_PAMGGetTruncFactor( (void *) solver, trunc_factor ) );
}

JXF_Int
JXF_PAMGSetPMaxElmts( JXF_Solver solver, JXF_Int P_max_elmts )
{
   return( jxf_PAMGSetPMaxElmts( (void *) solver, P_max_elmts ) );
}

JXF_Int
JXF_PAMGGetPMaxElmts( JXF_Solver solver, JXF_Int *P_max_elmts )
{
   return( jxf_PAMGGetPMaxElmts( (void *) solver, P_max_elmts ) );
}

JXF_Int
JXF_PAMGSetJacobiTruncThreshold( JXF_Solver solver, JXF_Real jacobi_trunc_threshold )
{
   return( jxf_PAMGSetJacobiTruncThreshold( (void *) solver, jacobi_trunc_threshold ) );
}

JXF_Int
JXF_PAMGGetJacobiTruncThreshold( JXF_Solver solver, JXF_Real *jacobi_trunc_threshold )
{
   return( jxf_PAMGGetJacobiTruncThreshold( (void *) solver, jacobi_trunc_threshold ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetPostInterpType
 * \brief If > 0, specifies something to do to improve a computed 
 *        interpolation matrix. defaults to 0, for nothing.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetPostInterpType( JXF_Solver solver, JXF_Int post_interp_type )
{
   return( jxf_PAMGSetPostInterpType( (void *) solver, post_interp_type ) );
}

JXF_Int
JXF_PAMGGetPostInterpType( JXF_Solver solver, JXF_Int *post_interp_type )
{
   return( jxf_PAMGGetPostInterpType( (void *) solver, post_interp_type ) );
}

JXF_Int
JXF_PAMGSetSCommPkgSwitch( JXF_Solver solver, JXF_Real S_commpkg_switch )
{
   return( jxf_PAMGSetSCommPkgSwitch( (void *) solver, S_commpkg_switch ) );
}

JXF_Int
JXF_PAMGSetInterpType( JXF_Solver solver, JXF_Int interp_type )
{
   return( jxf_PAMGSetInterpType( (void *) solver, interp_type ) );
}

JXF_Int
JXF_PAMGSetMinIter( JXF_Solver solver, JXF_Int min_iter )
{
   return( jxf_PAMGSetMinIter( (void *) solver, min_iter ) );
}

JXF_Int
JXF_PAMGSetMaxIter( JXF_Solver solver, JXF_Int max_iter )
{
   return( jxf_PAMGSetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_PAMGGetMaxIter( JXF_Solver solver, JXF_Int *max_iter )
{
   return( jxf_PAMGGetMaxIter( (void *) solver, max_iter ) );
}

JXF_Int
JXF_PAMGSetCoarsenType( JXF_Solver solver, JXF_Int coarsen_type )
{
   return( jxf_PAMGSetCoarsenType( (void *) solver, coarsen_type ) );
}

JXF_Int
JXF_PAMGGetCoarsenType( JXF_Solver solver, JXF_Int *coarsen_type )
{
   return( jxf_PAMGGetCoarsenType( (void *) solver, coarsen_type ) );
}

JXF_Int
JXF_PAMGSetMeasureType( JXF_Solver solver, JXF_Int measure_type )
{
   return( jxf_PAMGSetMeasureType( (void *) solver, measure_type ) );
}

JXF_Int
JXF_PAMGGetMeasureType( JXF_Solver solver, JXF_Int *measure_type )
{
   return( jxf_PAMGGetMeasureType( (void *) solver, measure_type ) );
}

JXF_Int
JXF_PAMGSetSetupType( JXF_Solver solver, JXF_Int setup_type )
{
   return( jxf_PAMGSetSetupType( (void *) solver, setup_type ) );
}

JXF_Int
JXF_PAMGSetCycleType( JXF_Solver solver,JXF_Int cycle_type )
{
   return( jxf_PAMGSetCycleType( (void *) solver, cycle_type ) );
}

JXF_Int
JXF_PAMGGetCycleType( JXF_Solver solver, JXF_Int *cycle_type )
{
   return( jxf_PAMGGetCycleType( (void *) solver, cycle_type ) );
}

JXF_Int
JXF_PAMGSetTol( JXF_Solver solver, JXF_Real tol )
{
   return( jxf_PAMGSetTol( (void *) solver, tol ) );
}

JXF_Int
JXF_PAMGSetRhsNrmThreshold( JXF_Solver solver, JXF_Real rhsnrm_threshold )
{
   return( jxf_PAMGSetRhsNrmThreshold( (void *) solver, rhsnrm_threshold ) );
}

JXF_Int
JXF_PAMGGetTol( JXF_Solver solver, JXF_Real *tol )
{
   return( jxf_PAMGGetTol( (void *) solver, tol ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetNumGridSweeps
 * \note DEPRECATED.  There are memory management problems associated with the
 *       use of a user-supplied array (who releases it?).
 *       Use SetNumSweeps and SetCycleNumSweeps instead.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetNumGridSweeps( JXF_Solver  solver, JXF_Int *num_grid_sweeps )
{
   return( jxf_PAMGSetNumGridSweeps( (void *) solver, num_grid_sweeps ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetNumSweeps
 * \note There is no corresponding Get function.  Use GetCycleNumSweeps.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetNumSweeps( JXF_Solver  solver, JXF_Int num_sweeps )
{
   return( jxf_PAMGSetNumSweeps( (void *) solver, num_sweeps ) );
}

JXF_Int
JXF_PAMGSetCycleNumSweeps( JXF_Solver solver, JXF_Int num_sweeps, JXF_Int k )
{
   return( jxf_PAMGSetCycleNumSweeps( (void *) solver, num_sweeps, k ) );
}

JXF_Int
JXF_PAMGGetCycleNumSweeps( JXF_Solver solver, JXF_Int *num_sweeps, JXF_Int k )
{
   return( jxf_PAMGGetCycleNumSweeps( (void *) solver, num_sweeps, k ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetGridRelaxType
 * \note DEPRECATED.  There are memory management problems associated with the
 *       use of a user-supplied array (who releases it?).
 *       Use SetRelaxType and SetCycleRelaxType instead.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetGridRelaxType( JXF_Solver  solver, JXF_Int *grid_relax_type )
{
   return( jxf_PAMGSetGridRelaxType( (void *) solver, grid_relax_type ) );
}

JXF_Int
JXF_PAMGSetRelaxType( JXF_Solver solver, JXF_Int relax_type  )
{
   return( jxf_PAMGSetRelaxType( (void *) solver, relax_type ) );
}

JXF_Int
JXF_PAMGSetCycleRelaxType( JXF_Solver solver, JXF_Int relax_type, JXF_Int k )
{
   return( jxf_PAMGSetCycleRelaxType( (void *) solver, relax_type, k ) );
}

JXF_Int
JXF_PAMGGetCycleRelaxType( JXF_Solver solver, JXF_Int * relax_type, JXF_Int k )
{
   return( jxf_PAMGGetCycleRelaxType( (void *) solver, relax_type, k ) );
}

JXF_Int
JXF_PAMGSetRelaxOrder( JXF_Solver solver, JXF_Int relax_order )
{
   return( jxf_PAMGSetRelaxOrder( (void *) solver, relax_order ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetGridRelaxPoints
 * \note DEPRECATED.  There are memory management problems associated with the
 *       use of a user-supplied array (who releases it?).
 *       Ulrike Yang suspects that nobody uses this function.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetGridRelaxPoints( JXF_Solver solver, JXF_Int **grid_relax_points )
{
   return( jxf_PAMGSetGridRelaxPoints( (void *) solver, grid_relax_points ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetRelaxWeight
 * \note DEPRECATED.  There are memory management problems associated with the
 *       use of a user-supplied array (who releases it?).
 *       Use SetRelaxWt and SetLevelRelaxWt instead.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetRelaxWeight( JXF_Solver solver, JXF_Real *relax_weight )
{
   return( jxf_PAMGSetRelaxWeight( (void *) solver, relax_weight ) );
}

JXF_Int
JXF_PAMGSetRelaxWt( JXF_Solver solver, JXF_Real relax_wt )
{
   return( jxf_PAMGSetRelaxWt( (void *) solver, relax_wt ) );
}

JXF_Int
JXF_PAMGSetLevelRelaxWt( JXF_Solver solver, JXF_Real relax_wt, JXF_Int level )
{
   return( jxf_PAMGSetLevelRelaxWt( (void *) solver, relax_wt, level ) );
}

JXF_Int
JXF_PAMGSetOmega( JXF_Solver solver, JXF_Real *omega )
{
   return( jxf_PAMGSetOmega( (void *) solver, omega ) );
}

JXF_Int
JXF_PAMGSetOuterWt( JXF_Solver solver, JXF_Real outer_wt )
{
   return( jxf_PAMGSetOuterWt( (void *) solver, outer_wt ) );
}

JXF_Int
JXF_PAMGSetLevelOuterWt( JXF_Solver solver, JXF_Real outer_wt, JXF_Int level )
{
   return( jxf_PAMGSetLevelOuterWt( (void *) solver, outer_wt, level ) );
}

JXF_Int
JXF_PAMGSetSmoothType( JXF_Solver solver, JXF_Int smooth_type )
{
   return( jxf_PAMGSetSmoothType( (void *) solver, smooth_type ) );
}

JXF_Int
JXF_PAMGGetSmoothType( JXF_Solver solver, JXF_Int * smooth_type )
{
   return( jxf_PAMGGetSmoothType( (void *) solver, smooth_type ) );
}

JXF_Int
JXF_PAMGSetSmoothNumLevels( JXF_Solver solver, JXF_Int smooth_num_levels )
{
   return( jxf_PAMGSetSmoothNumLevels((void *)solver,smooth_num_levels ));
}

JXF_Int
JXF_PAMGGetSmoothNumLevels( JXF_Solver solver, JXF_Int *smooth_num_levels )
{
   return( jxf_PAMGGetSmoothNumLevels((void *)solver,smooth_num_levels ));
}

JXF_Int
JXF_PAMGSetSmoothNumSweeps( JXF_Solver solver, JXF_Int smooth_num_sweeps )
{
   return( jxf_PAMGSetSmoothNumSweeps((void *)solver,smooth_num_sweeps ));
}

JXF_Int
JXF_PAMGGetSmoothNumSweeps( JXF_Solver solver, JXF_Int *smooth_num_sweeps )
{
   return( jxf_PAMGGetSmoothNumSweeps((void *)solver,smooth_num_sweeps ));
}

JXF_Int
JXF_PAMGSetLogging( JXF_Solver solver, JXF_Int logging )
{
  /*-------------------------------------------------------------------------
   * This function should be called before Setup.  Logging changes
   * may require allocation or freeing of arrays, which is presently
   * only done there.
   * It may be possible to support logging changes at other times,
   * but there is little need.
   *------------------------------------------------------------------------*/
   return( jxf_PAMGSetLogging( (void *) solver, logging ) );
}

JXF_Int
JXF_PAMGGetLogging( JXF_Solver solver, JXF_Int * logging )
{
   return( jxf_PAMGGetLogging( (void *) solver, logging ) );
}

JXF_Int
JXF_PAMGSetPrintLevel( JXF_Solver solver, JXF_Int print_level )
{
   return( jxf_PAMGSetPrintLevel( (void *) solver, print_level ) );
}

JXF_Int
JXF_PAMGSetPrintCoarseMatrix( JXF_Solver solver, JXF_Int print_coarse_matrix )
{
   return( jxf_PAMGSetPrintCoarseMatrix( (void *) solver, print_coarse_matrix ) );
}

JXF_Int
JXF_PAMGGetPrintLevel( JXF_Solver solver, JXF_Int *print_level )
{
   return( jxf_PAMGGetPrintLevel( (void *) solver, print_level ) );
}

JXF_Int
JXF_PAMGSetPrintFileName( JXF_Solver solver, const char *print_file_name )
{
   return( jxf_PAMGSetPrintFileName( (void *) solver, print_file_name ) );
}

JXF_Int
JXF_PAMGSetDebugFlag( JXF_Solver solver, JXF_Int debug_flag )
{
   return( jxf_PAMGSetDebugFlag( (void *) solver, debug_flag ) );
}

JXF_Int
JXF_PAMGGetDebugFlag( JXF_Solver solver, JXF_Int *debug_flag )
{
   return( jxf_PAMGGetDebugFlag( (void *) solver, debug_flag ) );
}

JXF_Int
JXF_PAMGGetNumLevels( JXF_Solver solver, JXF_Int *num_levels )
{
   return( jxf_PAMGGetNumLevels( (void *) solver, num_levels ) );
}

JXF_Int
JXF_PAMGGetNumIterations( JXF_Solver solver, JXF_Int *num_iterations )
{
   return( jxf_PAMGGetNumIterations( (void *) solver, num_iterations ) );
}

JXF_Int
JXF_PAMGGetCumNumIterations( JXF_Solver solver, JXF_Int *cum_num_iterations )
{
   return( jxf_PAMGGetCumNumIterations( (void *) solver, cum_num_iterations ) );
}

JXF_Int
JXF_PAMGGetResidual( JXF_Solver solver, JXF_ParVector *residual )
{
   return jxf_PAMGGetResidual( (void *) solver, (jxf_ParVector **)residual );
}
                            
JXF_Int
JXF_PAMGGetFinalRelativeResidualNorm( JXF_Solver solver, JXF_Real *rel_resid_norm )
{
   return( jxf_PAMGGetRelResidualNorm( (void *) solver, rel_resid_norm ) );
}

JXF_Int
JXF_PAMGSetVariant( JXF_Solver solver,  JXF_Int variant )
{
   return( jxf_PAMGSetVariant( (void *) solver, variant ) );
}

JXF_Int
JXF_PAMGGetVariant( JXF_Solver solver, JXF_Int *variant )
{
   return( jxf_PAMGGetVariant( (void *) solver, variant ) );
}

JXF_Int
JXF_PAMGSetOverlap( JXF_Solver solver, JXF_Int overlap )
{
   return( jxf_PAMGSetOverlap( (void *) solver, overlap ) );
}

JXF_Int
JXF_PAMGGetOverlap( JXF_Solver solver, JXF_Int *overlap )
{
   return( jxf_PAMGGetOverlap( (void *) solver, overlap ) );
}

JXF_Int
JXF_PAMGSetDomainType( JXF_Solver solver, JXF_Int domain_type )
{
   return( jxf_PAMGSetDomainType( (void *) solver, domain_type ) );
}

JXF_Int
JXF_PAMGGetDomainType( JXF_Solver solver, JXF_Int *domain_type )
{
   return( jxf_PAMGGetDomainType( (void *) solver, domain_type ) );
}

JXF_Int
JXF_PAMGSetSchwarzRlxWeight( JXF_Solver solver, JXF_Real schwarz_rlx_weight )
{
   return( jxf_PAMGSetSchwarzRlxWeight( (void *)solver, schwarz_rlx_weight ) );
}

JXF_Int
JXF_PAMGGetSchwarzRlxWeight( JXF_Solver solver, JXF_Real *schwarz_rlx_weight )
{
   return( jxf_PAMGGetSchwarzRlxWeight( (void *)solver, schwarz_rlx_weight ) );
}

JXF_Int
JXF_PAMGSetSchwarzUseNonSymm( JXF_Solver solver, JXF_Int use_nonsymm )
{
   return( jxf_PAMGSetSchwarzUseNonSymm( (void *) solver, use_nonsymm ) );
}

JXF_Int
JXF_PAMGSetSym( JXF_Solver solver, JXF_Int sym )
{
   return( jxf_PAMGSetSym( (void *) solver, sym ) );
}

JXF_Int
JXF_PAMGSetLevel( JXF_Solver solver, JXF_Int level )
{
   return( jxf_PAMGSetLevel( (void *) solver, level ) );
}

JXF_Int
JXF_PAMGSetThreshold( JXF_Solver solver, JXF_Real threshold )
{
   return( jxf_PAMGSetThreshold( (void *) solver, threshold ) );
}

JXF_Int
JXF_PAMGSetFilter( JXF_Solver solver, JXF_Real filter )
{
   return( jxf_PAMGSetFilter( (void *) solver, filter ) );
}

JXF_Int
JXF_PAMGSetDropTol( JXF_Solver solver, JXF_Real drop_tol )
{
   return( jxf_PAMGSetDropTol( (void *) solver, drop_tol ) );
}

JXF_Int
JXF_PAMGSetMaxNzPerRow( JXF_Solver solver, JXF_Int max_nz_per_row )
{
   return( jxf_PAMGSetMaxNzPerRow( (void *) solver, max_nz_per_row ) );
}

JXF_Int
JXF_PAMGSetEuclidFile( JXF_Solver solver, char *euclidfile )
{
   return( jxf_PAMGSetEuclidFile( (void *) solver, euclidfile ) );
}

JXF_Int
JXF_PAMGSetEuLevel( JXF_Solver solver, JXF_Int eu_level )
{
   return( jxf_PAMGSetEuLevel( (void *) solver, eu_level ) );
}

JXF_Int
JXF_PAMGSetEuSparseA( JXF_Solver solver, JXF_Real eu_sparse_A )
{
   return( jxf_PAMGSetEuSparseA( (void *) solver, eu_sparse_A ) );
}

JXF_Int
JXF_PAMGSetEuBJ( JXF_Solver solver, JXF_Int eu_bj )
{
   return( jxf_PAMGSetEuBJ( (void *) solver, eu_bj ) );
}

JXF_Int
JXF_PAMGSetNumFunctions( JXF_Solver solver, JXF_Int num_functions )
{
   return( jxf_PAMGSetNumFunctions( (void *) solver, num_functions ) );
}

JXF_Int
JXF_PAMGGetNumFunctions( JXF_Solver solver, JXF_Int *num_functions )
{
   return( jxf_PAMGGetNumFunctions( (void *) solver, num_functions ) );
}

JXF_Int
JXF_PAMGSetNodal( JXF_Solver solver, JXF_Int nodal )
{
   return( jxf_PAMGSetNodal( (void *) solver, nodal ) );
}

JXF_Int
JXF_PAMGSetNodalDiag( JXF_Solver solver, JXF_Int nodal )
{
   return( jxf_PAMGSetNodalDiag( (void *) solver, nodal ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetDofFunc
 * \note Warning about a possible memory problem: When the PAMG object is 
 *       destroyed in jxf_PAMGDestroy, dof_func aka DofFunc will be destroyed
 *       Normally this is what we want.  But if the user provided dof_func
 *       by calling JXF_PAMGSetDofFunc, this could be an unwanted surprise.
 *       As jxf is currently commonly used, this situation is likely to be rare.
 * \date 2011/09/03
 */
JXF_Int
JXF_PAMGSetDofFunc( JXF_Solver solver, JXF_Int  *dof_func )
{
   return( jxf_PAMGSetDofFunc( (void *) solver, dof_func ) );
}

JXF_Int
JXF_PAMGSetNumPaths( JXF_Solver solver, JXF_Int num_paths )
{
   return( jxf_PAMGSetNumPaths( (void *) solver, num_paths ) );
}

JXF_Int
JXF_PAMGSetAggNumLevels( JXF_Solver solver, JXF_Int agg_num_levels )
{
   return( jxf_PAMGSetAggNumLevels( (void *) solver, agg_num_levels ) );
}

JXF_Int
JXF_PAMGSetAggInterpType( JXF_Solver solver, JXF_Int agg_interp_type )
{
   return( jxf_PAMGSetAggInterpType( (void *) solver, agg_interp_type ) );
}

JXF_Int
JXF_PAMGSetAggPMaxElmts( JXF_Solver solver, JXF_Int agg_P_max_elmts )
{
   return( jxf_PAMGSetAggPMaxElmts( (void *) solver, agg_P_max_elmts ) );
}

JXF_Int
JXF_PAMGSetAggP12MaxElmts( JXF_Solver solver, JXF_Int agg_P12_max_elmts )
{
   return( jxf_PAMGSetAggP12MaxElmts( (void *) solver, agg_P12_max_elmts ) );
}

JXF_Int
JXF_PAMGSetAggTruncFactor( JXF_Solver solver, JXF_Real agg_trunc_factor )
{
   return( jxf_PAMGSetAggTruncFactor( (void *) solver, agg_trunc_factor ) );
}

JXF_Int
JXF_PAMGSetAggP12TruncFactor( JXF_Solver solver, JXF_Real agg_P12_trunc_factor )
{
   return( jxf_PAMGSetAggP12TruncFactor( (void *) solver, agg_P12_trunc_factor ) );
}

JXF_Int
JXF_PAMGSetNumCRRelaxSteps( JXF_Solver solver, JXF_Int num_CR_relax_steps )
{
   return( jxf_PAMGSetNumCRRelaxSteps( (void *) solver, num_CR_relax_steps ) );
}

JXF_Int
JXF_PAMGSetCRRate( JXF_Solver solver, JXF_Real CR_rate )
{
   return( jxf_PAMGSetCRRate( (void *) solver, CR_rate ) );
}

JXF_Int
JXF_PAMGSetCRStrongTh( JXF_Solver solver, JXF_Real CR_strong_th )
{
   return( jxf_PAMGSetCRStrongTh( (void *) solver, CR_strong_th ) );
}

JXF_Int
JXF_PAMGSetISType( JXF_Solver solver, JXF_Int IS_type )
{
   return( jxf_PAMGSetISType( (void *) solver, IS_type ) );
}

JXF_Int
JXF_PAMGSetCRUseCG( JXF_Solver solver, JXF_Int CR_use_CG )
{
   return( jxf_PAMGSetCRUseCG( (void *) solver, CR_use_CG ) );
}

JXF_Int
JXF_PAMGSetGSMG( JXF_Solver solver, JXF_Int gsmg )
{
   return( jxf_PAMGSetGSMG( (void *) solver, gsmg ) );
}

JXF_Int
JXF_PAMGSetNumSamples( JXF_Solver solver, JXF_Int gsmg )
{
   return( jxf_PAMGSetNumSamples( (void *) solver, gsmg ) );
}
                                                                                                                                                                                                            
JXF_Int
JXF_PAMGSetCGCIts ( JXF_Solver solver, JXF_Int its )
{
  return (jxf_PAMGSetCGCIts ( (void *) solver, its ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetSpMtRapType
 * \author Yue Xiaoqiang
 * \date 2012/10/13
 */
JXF_Int
JXF_PAMGSetSpMtRapType( JXF_Solver solver, JXF_Int spmt_rap_type )
{
   return (jxf_PAMGSetSpMtRapType ( (void *) solver, spmt_rap_type ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetWallTimeOption
 * \author Yue Xiaoqiang
 * \date 2015/09/30
 */
JXF_Int
JXF_PAMGSetWallTimeOption( JXF_Solver solver, JXF_Int wall_time_option )
{
   return (jxf_PAMGSetWallTimeOption ( (void *) solver, wall_time_option ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetAIMeasureType
 * \author Yue Xiaoqiang
 * \date 2014/02/26
 */
JXF_Int
JXF_PAMGSetAIMeasureType( JXF_Solver solver, JXF_Int ai_measure_type )
{
   return (jxf_PAMGSetAIMeasureType ( (void *) solver, ai_measure_type ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetAIRelaxType
 * \author Yue Xiaoqiang
 * \date 2014/07/06
 */
JXF_Int
JXF_PAMGSetAIRelaxType( JXF_Solver solver, JXF_Int ai_relax_type )
{
   return (jxf_PAMGSetAIRelaxType ( (void *) solver, ai_relax_type ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetRelaxedCoarsenMeasureType
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
JXF_PAMGSetRelaxedCoarsenMeasureType( JXF_Solver solver, JXF_Int measure_type_rlx )
{
   return (jxf_PAMGSetRelaxedCoarsenMeasureType ( (void *) solver, measure_type_rlx ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetRelaxedCoarsenNumberSyn
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
JXF_PAMGSetRelaxedCoarsenNumberSyn( JXF_Solver solver, JXF_Int number_syn_rlx )
{
   return (jxf_PAMGSetRelaxedCoarsenNumberSyn ( (void *) solver, number_syn_rlx ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetRelaxedCoarsenMeasureThreshold
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
JXF_PAMGSetRelaxedCoarsenMeasureThreshold( JXF_Solver solver, JXF_Real measure_threshold_rlx )
{
   return (jxf_PAMGSetRelaxedCoarsenMeasureThreshold ( (void *) solver, measure_threshold_rlx ) );
}
                                                                                                       
JXF_Int
JXF_PAMGSetPlotGrids ( JXF_Solver solver, JXF_Int plotgrids )
{
   return (jxf_PAMGSetPlotGrids ( (void *) solver, plotgrids ) );
}
                                                                                                                                                                                                            
JXF_Int
JXF_PAMGSetPlotFileName ( JXF_Solver solver, const char *plotfilename )
{
   return (jxf_PAMGSetPlotFileName ( (void *) solver, plotfilename ) );
}
                                                                                                                                                                                                              
JXF_Int
JXF_PAMGSetCoordDim ( JXF_Solver solver, JXF_Int coorddim )
{
   return (jxf_PAMGSetCoordDim ( (void *) solver, coorddim ) );
}
                                                                                                                                                                                                           
JXF_Int
JXF_PAMGSetCoordinates ( JXF_Solver solver, float *coordinates )
{
   return (jxf_PAMGSetCoordinates ( (void *) solver, coordinates ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetConvCriteria
 * \author peghoty
 * \date 2009/07/25
 */
JXF_Int
JXF_PAMGSetConvCriteria( JXF_Solver solver, JXF_Int conv_criteria )
{
   return (jxf_PAMGSetConvCriteria ( (void *) solver, conv_criteria ) );
}

JXF_Int
JXF_PAMGSetCoarsestSolverID( JXF_Solver solver, JXF_Int coarsestsolverid )
{
   return (jxf_PAMGSetCoarsestSolverID ( (void *) solver, coarsestsolverid ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetCoarseThreshold
 * \author peghoty
 * \date 2010/04/14
 */
JXF_Int
JXF_PAMGSetCoarseThreshold( JXF_Solver solver, JXF_Int coarse_threshold )
{
   return (jxf_PAMGSetCoarseThreshold ( (void *) solver, coarse_threshold ) );
}

/*!
 * \fn JXF_Int JXF_PAMGSetCoarseRatio
 * \author peghoty
 * \date 2010/04/14
 */
JXF_Int
JXF_PAMGSetCoarseRatio( JXF_Solver solver, JXF_Real coarse_ratio )
{
   return (jxf_PAMGSetCoarseRatio ( (void *) solver, coarse_ratio ) );
}

JXF_Int
JXF_PAMGSetConvFacThreshold( JXF_Solver solver, JXF_Real convfac_threshold )
{
   return (jxf_PAMGSetConvFacThreshold ( (void *) solver, convfac_threshold ) );
}

JXF_Int
JXF_PAMGSetRAP2( JXF_Solver solver, JXF_Int rap2 )
{
   return (jxf_PAMGSetRAP2 ( (void *) solver, rap2 ) );
}

JXF_Int
JXF_PAMGSetKeepTranspose( JXF_Solver solver, JXF_Int keepTranspose )
{
   return (jxf_PAMGSetKeepTranspose ( (void *) solver, keepTranspose ) );
}


/*!
 * \fn JXF_Int jxf_PAMGSet*
 * \brief  Routines to set the setup phase parameters begin with "jxf_"
 * \date 2010/04/14
 */
JXF_Int
jxf_PAMGSetRestriction( void *data, JXF_Int restr_par )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   /* RL: currently, 0: R = P^T, 1: AIR, >1: AIR-2 */
   if (restr_par < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataRestriction(amg_data) = restr_par;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetAIRStrongTh( void *data, JXF_Real AIR_strong_th )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataAIRStrongTh(amg_data) = AIR_strong_th;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetMaxLevels( void *data, JXF_Int max_levels )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (max_levels < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataMaxLevels(amg_data) = max_levels;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetMaxLevels( void *data, JXF_Int *max_levels )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *max_levels = jxf_ParAMGDataMaxLevels(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetStrongThreshold( void *data, JXF_Real strong_threshold )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (strong_threshold < 0 || strong_threshold > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataStrongThreshold(amg_data) = strong_threshold;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetStrongThreshold( void *data, JXF_Real *strong_threshold )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *strong_threshold = jxf_ParAMGDataStrongThreshold(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetMaxRowSum( void *data, JXF_Real max_row_sum )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (max_row_sum <= 0 || max_row_sum > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataMaxRowSum(amg_data) = max_row_sum;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetMaxRowSum( void *data, JXF_Real *max_row_sum )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *max_row_sum = jxf_ParAMGDataMaxRowSum(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetTruncFactor( void *data, JXF_Real trunc_factor )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (trunc_factor < 0 || trunc_factor >= 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataTruncFactor(amg_data) = trunc_factor;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetTruncFactor( void *data, JXF_Real *trunc_factor )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *trunc_factor = jxf_ParAMGDataTruncFactor(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPMaxElmts( void *data, JXF_Int P_max_elmts )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (P_max_elmts < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataPMaxElmts(amg_data) = P_max_elmts;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetPMaxElmts( void *data, JXF_Int *P_max_elmts )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *P_max_elmts = jxf_ParAMGDataPMaxElmts(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetJacobiTruncThreshold( void *data, JXF_Real jacobi_trunc_threshold )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (jacobi_trunc_threshold < 0 || jacobi_trunc_threshold >= 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataJacobiTruncThreshold(amg_data) = jacobi_trunc_threshold;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetJacobiTruncThreshold( void *data, JXF_Real *jacobi_trunc_threshold )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *jacobi_trunc_threshold = jxf_ParAMGDataJacobiTruncThreshold(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPostInterpType( void *data, JXF_Int post_interp_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (post_interp_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataPostInterpType(amg_data) = post_interp_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetPostInterpType( void *data, JXF_Int *post_interp_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *post_interp_type = jxf_ParAMGDataPostInterpType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSCommPkgSwitch( void *data, JXF_Real S_commpkg_switch )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataSCommPkgSwitch(amg_data) = S_commpkg_switch;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSCommPkgSwitch( void *data, JXF_Real *S_commpkg_switch )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *S_commpkg_switch = jxf_ParAMGDataSCommPkgSwitch(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetInterpType( void *data, JXF_Int interp_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

  /*------------------------------------------------------------------------
   *  “interp_type > 25” has been modified to “interp_type > 100” since
   *  we may add other interpolations. peghoty  2009/08/08  
   *------------------------------------------------------------------------*/
   if (interp_type < 0 || interp_type > 100)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataInterpType(amg_data) = interp_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetInterpType( void *data, JXF_Int *interp_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *interp_type = jxf_ParAMGDataInterpType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetMinIter( void *data, JXF_Int min_iter )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataMinIter(amg_data) = min_iter;

   return jxf_error_flag;
} 

JXF_Int
jxf_PAMGGetMinIter( void *data, JXF_Int *min_iter )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *min_iter = jxf_ParAMGDataMinIter(amg_data);

   return jxf_error_flag;
} 

JXF_Int
jxf_PAMGSetMaxIter( void *data, JXF_Int max_iter )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (max_iter < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataMaxIter(amg_data) = max_iter;

   return jxf_error_flag;
} 

JXF_Int
jxf_PAMGGetMaxIter( void *data, JXF_Int *max_iter )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *max_iter = jxf_ParAMGDataMaxIter(amg_data);

   return jxf_error_flag;
} 

JXF_Int
jxf_PAMGSetCoarsenType( void *data, JXF_Int coarsen_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataCoarsenType(amg_data) = coarsen_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetCoarsenType( void *data, JXF_Int *coarsen_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *coarsen_type = jxf_ParAMGDataCoarsenType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetMeasureType( void *data, JXF_Int measure_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataMeasureType(amg_data) = measure_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetMeasureType( void *data, JXF_Int *measure_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *measure_type = jxf_ParAMGDataMeasureType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSetupType( void *data, JXF_Int setup_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   jxf_ParAMGDataSetupType(amg_data) = setup_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSetupType( void *data, JXF_Int *setup_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *setup_type = jxf_ParAMGDataSetupType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetCycleType( void *data, JXF_Int cycle_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (cycle_type < 0 || cycle_type > 2)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataCycleType(amg_data) = cycle_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetCycleType( void *data, JXF_Int *cycle_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *cycle_type = jxf_ParAMGDataCycleType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetTol( void *data, JXF_Real tol  )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (tol < 0 || tol > 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataTol(amg_data) = tol;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetRhsNrmThreshold( void *data, JXF_Real rhsnrm_threshold )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (rhsnrm_threshold < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   jxf_ParAMGDataRhsNrmThreshold(amg_data) = rhsnrm_threshold;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetTol( void *data, JXF_Real *tol )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   *tol = jxf_ParAMGDataTol(amg_data);

   return jxf_error_flag;
}

/* The "Get" function for SetNumSweeps is GetCycleNumSweeps. */
JXF_Int
jxf_PAMGSetNumSweeps( void *data, JXF_Int num_sweeps )  /* default: num_sweeps = 1 */
{
   JXF_Int i;
   JXF_Int *num_grid_sweeps;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (num_sweeps < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataNumGridSweeps(amg_data) == NULL)
   {
      jxf_ParAMGDataNumGridSweeps(amg_data) = jxf_CTAlloc(JXF_Int, 4);
   }
       
   num_grid_sweeps = jxf_ParAMGDataNumGridSweeps(amg_data);

   for (i = 0; i < 3; i ++)
   {
      num_grid_sweeps[i] = num_sweeps;
   }
   num_grid_sweeps[3] = 1;

   return jxf_error_flag;
}
 
JXF_Int
jxf_PAMGSetCycleNumSweeps( void *data, JXF_Int num_sweeps, JXF_Int k )
{
   JXF_Int i;
   JXF_Int *num_grid_sweeps;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 

   if (num_sweeps < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (k < 1 || k > 3)
   {
      jxf_printf (" Warning! Invalid cycle! num_sweeps not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataNumGridSweeps(amg_data) == NULL)
   {
       num_grid_sweeps = jxf_CTAlloc(JXF_Int, 4);
       for (i = 0; i < 4; i ++)
       {
	  num_grid_sweeps[i] = 1;
       }
       jxf_ParAMGDataNumGridSweeps(amg_data) = num_grid_sweeps;
   }
       
   jxf_ParAMGDataNumGridSweeps(amg_data)[k] = num_sweeps;

   return jxf_error_flag;
}
 
JXF_Int
jxf_PAMGGetCycleNumSweeps( void *data, JXF_Int *num_sweeps, JXF_Int k )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (k < 1 || k > 3)
   {
      jxf_printf (" Warning! Invalid cycle! No num_sweeps to get!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataNumGridSweeps(amg_data) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
       
   *num_sweeps = jxf_ParAMGDataNumGridSweeps(amg_data)[k];

   return jxf_error_flag;
}
 
JXF_Int
jxf_PAMGSetNumGridSweeps( void *data, JXF_Int *num_grid_sweeps )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (!num_grid_sweeps)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataNumGridSweeps(amg_data))
   {
      jxf_TFree(jxf_ParAMGDataNumGridSweeps(amg_data));
   }
   jxf_ParAMGDataNumGridSweeps(amg_data) = num_grid_sweeps;

   return jxf_error_flag;
}
 
JXF_Int
jxf_PAMGGetNumGridSweeps( void *data, JXF_Int **num_grid_sweeps )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *num_grid_sweeps = jxf_ParAMGDataNumGridSweeps(amg_data);

   return jxf_error_flag;
}

/*!
 * \fn void jxf_dispose_elt
 * \note The "Get" function for SetRelaxType is GetCycleRelaxType.
 * \date 2011/09/03
 */
JXF_Int
jxf_PAMGSetRelaxType( void *data, JXF_Int relax_type ) /* default: relax_type = 3 */
{
   JXF_Int i;
   JXF_Int *grid_relax_type;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (relax_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataGridRelaxType(amg_data) == NULL)
   {
      jxf_ParAMGDataGridRelaxType(amg_data) = jxf_CTAlloc(JXF_Int, 4);
   }
   grid_relax_type = jxf_ParAMGDataGridRelaxType(amg_data);

   for (i = 0; i < 3; i ++)
   {
      grid_relax_type[i] = relax_type;
   }
   grid_relax_type[3] = 9;
   jxf_ParAMGDataUserCoarseRelaxType(amg_data) = 9;
   jxf_ParAMGDataUserRelaxType(amg_data) = relax_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetCycleRelaxType( void *data, JXF_Int relax_type, JXF_Int k )
{
   JXF_Int i;
   JXF_Int *grid_relax_type;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (k < 1 || k > 3)
   {
      jxf_printf (" Warning! Invalid cycle! relax_type not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (relax_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataGridRelaxType(amg_data) == NULL)
   {
      grid_relax_type = jxf_CTAlloc(JXF_Int, 4);
      for (i = 0; i < 3; i ++)
      {
         grid_relax_type[i] = 3;
      }
      grid_relax_type[3] = 9;
      jxf_ParAMGDataGridRelaxType(amg_data) = grid_relax_type;
   }
      
   jxf_ParAMGDataGridRelaxType(amg_data)[k] = relax_type;
   if (k == 3)
   {
      jxf_ParAMGDataUserCoarseRelaxType(amg_data) = relax_type;
   }
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetCycleRelaxType( void *data, JXF_Int *relax_type, JXF_Int k )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (k < 1 || k > 3)
   {
      jxf_printf (" Warning! Invalid cycle! relax_type not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataGridRelaxType(amg_data) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
      
   *relax_type = jxf_ParAMGDataGridRelaxType(amg_data)[k];

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetRelaxOrder( void *data, JXF_Int relax_order )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataRelaxOrder(amg_data) = relax_order;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetRelaxOrder( void *data, JXF_Int *relax_order )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *relax_order = jxf_ParAMGDataRelaxOrder(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetGridRelaxType( void *data, JXF_Int  *grid_relax_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (!grid_relax_type)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataGridRelaxType(amg_data))
   {
      jxf_TFree(jxf_ParAMGDataGridRelaxType(amg_data));
   }
   jxf_ParAMGDataGridRelaxType(amg_data) = grid_relax_type;
   jxf_ParAMGDataUserCoarseRelaxType(amg_data) = grid_relax_type[3];

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetGridRelaxType( void *data, JXF_Int **grid_relax_type )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *grid_relax_type = jxf_ParAMGDataGridRelaxType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetGridRelaxPoints( void *data, JXF_Int **grid_relax_points )
{
   JXF_Int i;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (!grid_relax_points)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataGridRelaxPoints(amg_data))
   {
      for (i = 0; i < 4; i ++)
      {
   	 jxf_TFree (jxf_ParAMGDataGridRelaxPoints(amg_data)[i]);
      }
      jxf_TFree(jxf_ParAMGDataGridRelaxPoints(amg_data));
   }
   jxf_ParAMGDataGridRelaxPoints(amg_data) = grid_relax_points; 

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetGridRelaxPoints( void *data, JXF_Int ***grid_relax_points )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *grid_relax_points = jxf_ParAMGDataGridRelaxPoints(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetRelaxWeight( void *data, JXF_Real *relax_weight )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (!relax_weight)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }

   if (jxf_ParAMGDataRelaxWeight(amg_data))
   {
      jxf_TFree(jxf_ParAMGDataRelaxWeight(amg_data));
   }
   jxf_ParAMGDataRelaxWeight(amg_data) = relax_weight;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetRelaxWeight( void *data, JXF_Real **relax_weight )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *relax_weight = jxf_ParAMGDataRelaxWeight(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetRelaxWt( void *data, JXF_Real relax_weight )
{
   JXF_Int i, num_levels;
   JXF_Real *relax_weight_array;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (jxf_ParAMGDataRelaxWeight(amg_data) == NULL)
   {
      jxf_ParAMGDataRelaxWeight(amg_data) = jxf_CTAlloc(JXF_Real, num_levels);
   }
                     
   relax_weight_array = jxf_ParAMGDataRelaxWeight(amg_data);
   for (i = 0; i < num_levels; i ++)
   {
      relax_weight_array[i] = relax_weight;
   }
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetLevelRelaxWt( void *data, JXF_Real relax_weight, JXF_Int level )
{
   JXF_Int i, num_levels;
   jxf_ParAMGData  *amg_data = data;
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (level > num_levels-1 || level < 0) 
   {
      jxf_printf (" Warning! Invalid level! Relax weight not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (jxf_ParAMGDataRelaxWeight(amg_data) == NULL)
   {
      jxf_ParAMGDataRelaxWeight(amg_data) = jxf_CTAlloc(JXF_Real, num_levels);
      for (i = 0; i < num_levels; i ++)
      {
         jxf_ParAMGDataRelaxWeight(amg_data)[i] = 1.0;
      }
   }
               
   jxf_ParAMGDataRelaxWeight(amg_data)[level] = relax_weight;
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetLevelRelaxWt( void *data, JXF_Real *relax_weight, JXF_Int level )
{
   JXF_Int num_levels;
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (level > num_levels-1 || level < 0) 
   {
      jxf_printf (" Warning! Invalid level! Relax weight not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (jxf_ParAMGDataRelaxWeight(amg_data) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
               
   *relax_weight = jxf_ParAMGDataRelaxWeight(amg_data)[level];
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetOmega( void *data, JXF_Real *omega )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (!omega)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   if (jxf_ParAMGDataOmega(amg_data))
   {
      jxf_TFree(jxf_ParAMGDataOmega(amg_data));
   }
   jxf_ParAMGDataOmega(amg_data) = omega;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetOmega( void *data, JXF_Real **omega )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *omega = jxf_ParAMGDataOmega(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetOuterWt( void *data, JXF_Real omega )
{
   JXF_Int i, num_levels;
   JXF_Real *omega_array;
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (jxf_ParAMGDataOmega(amg_data) == NULL)
   {
      jxf_ParAMGDataOmega(amg_data) = jxf_CTAlloc(JXF_Real, num_levels);
   }
                     
   omega_array = jxf_ParAMGDataOmega(amg_data);
   for (i = 0; i < num_levels; i ++)
   {      
      omega_array[i] = omega;
   }
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetLevelOuterWt( void *data, JXF_Real omega, JXF_Int level )
{
   JXF_Int i, num_levels;
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (level > num_levels-1) 
   {
      jxf_printf (" Warning! Invalid level! Outer weight not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (jxf_ParAMGDataOmega(amg_data) == NULL)
   {
      jxf_ParAMGDataOmega(amg_data) = jxf_CTAlloc(JXF_Real, num_levels);
      for (i = 0; i < num_levels; i ++)
      {
         jxf_ParAMGDataOmega(amg_data)[i] = 1.0;
      }
   }
               
   jxf_ParAMGDataOmega(amg_data)[level] = omega;
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetLevelOuterWt( void *data, JXF_Real  *omega, JXF_Int level )
{
   JXF_Int num_levels;
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   num_levels = jxf_ParAMGDataMaxLevels(amg_data);
   if (level > num_levels-1) 
   {
      jxf_printf (" Warning! Invalid level! Outer weight not set!\n");
      jxf_error_in_arg(3);
      return jxf_error_flag;
   }
   if (jxf_ParAMGDataOmega(amg_data) == NULL)
   {
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
               
   *omega = jxf_ParAMGDataOmega(amg_data)[level];
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSmoothType( void *data, JXF_Int smooth_type )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
               
   jxf_ParAMGDataSmoothType(amg_data) = smooth_type;
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSmoothType( void *data, JXF_Int *smooth_type )
{
   jxf_ParAMGData  *amg_data = data;
               
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *smooth_type = jxf_ParAMGDataSmoothType(amg_data);
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSmoothNumLevels( void *data, JXF_Int smooth_num_levels )
{
   jxf_ParAMGData  *amg_data = data;
               
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (smooth_num_levels < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataSmoothNumLevels(amg_data) = smooth_num_levels;
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSmoothNumLevels( void *data, JXF_Int *smooth_num_levels )
{
   jxf_ParAMGData  *amg_data = data;
               
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *smooth_num_levels = jxf_ParAMGDataSmoothNumLevels(amg_data);
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSmoothNumSweeps( void *data, JXF_Int smooth_num_sweeps )
{
   jxf_ParAMGData  *amg_data = data;
               
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (smooth_num_sweeps < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataSmoothNumSweeps(amg_data) = smooth_num_sweeps;
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSmoothNumSweeps( void *data, JXF_Int *smooth_num_sweeps )
{
   jxf_ParAMGData  *amg_data = data;
               
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *smooth_num_sweeps = jxf_ParAMGDataSmoothNumSweeps(amg_data);
   
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetLogging( void *data, JXF_Int logging )
{
  /*-------------------------------------------------------------------------
   * This function should be called before Setup.  Logging changes
   * may require allocation or freeing of arrays, which is presently
   * only done there.
   * It may be possible to support logging changes at other times,
   * but there is little need.
   *------------------------------------------------------------------------*/      
      
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataLogging(amg_data) = logging;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetLogging( void *data, JXF_Int *logging )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *logging = jxf_ParAMGDataLogging(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPrintLevel( void *data, JXF_Int print_level )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataPrintLevel(amg_data) = print_level;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPrintCoarseMatrix( void *data, JXF_Int print_coarse_matrix )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataPrintCoarseSystem(amg_data) = print_coarse_matrix;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetPrintLevel( void *data, JXF_Int *print_level )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *print_level = jxf_ParAMGDataPrintLevel(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPrintFileName( void *data, const char *print_file_name )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if( strlen(print_file_name) > 256 )
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 

   jxf_sprintf(jxf_ParAMGDataLogFileName(amg_data), "%s", print_file_name);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetPrintFileName( void *data, char **print_file_name )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_sprintf( *print_file_name, "%s", jxf_ParAMGDataLogFileName(amg_data) );

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetNumIterations( void *data, JXF_Int num_iterations )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNumIterations(amg_data) = num_iterations;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetDebugFlag( void *data, JXF_Int debug_flag )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataDebugFlag(amg_data) = debug_flag;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetDebugFlag( void *data, JXF_Int  *debug_flag )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *debug_flag = jxf_ParAMGDataDebugFlag(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetGSMG( void *data, JXF_Int par )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   amg_data->gsmg = par;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetNumSamples( void *data, JXF_Int par )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   amg_data->num_samples = par;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetCGCIts( void *data, JXF_Int its )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataCGCIts(amg_data) = its;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetSpMtRapType
 * \author Yue Xiaoqiang
 * \date 2012/10/13
 */
JXF_Int
jxf_PAMGSetSpMtRapType( void *data, JXF_Int spmt_rap_type )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataSpMtRapType(amg_data) = spmt_rap_type;
  return (ierr);
}


/*!
 * \fn JXF_Int jxf_PAMGSetWallTimeOption
 * \author Yue Xiaoqiang
 * \date 2015/09/30
 */
JXF_Int
jxf_PAMGSetWallTimeOption( void *data, JXF_Int wall_time_option )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataWallTimeOption(amg_data) = wall_time_option;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetAIMeasureType
 * \author Yue Xiaoqiang
 * \date 2014/02/26
 */
JXF_Int
jxf_PAMGSetAIMeasureType( void *data, JXF_Int ai_measure_type )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataAIMeasureType(amg_data) = ai_measure_type;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetAIRelaxType
 * \author Yue Xiaoqiang
 * \date 2014/07/06
 */
JXF_Int
jxf_PAMGSetAIRelaxType( void *data, JXF_Int ai_relax_type )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataAIRelaxType(amg_data) = ai_relax_type;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetRelaxedCoarsenMeasureType
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetRelaxedCoarsenMeasureType( void *data, JXF_Int measure_type_rlx )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataMeasureTypeRlx(amg_data) = measure_type_rlx;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetRelaxedCoarsenNumberSyn
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetRelaxedCoarsenNumberSyn( void *data, JXF_Int number_syn_rlx )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataNumberSynRlx(amg_data) = number_syn_rlx;
  return (ierr);
}

/*!
 * \fn JXF_Int jxf_PAMGSetRelaxedCoarsenMeasureThreshold
 * \author peghoty
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetRelaxedCoarsenMeasureThreshold( void *data, JXF_Real measure_threshold_rlx )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataMeasureThresholdRlx(amg_data) = measure_threshold_rlx;
  return (ierr);
}

JXF_Int
jxf_PAMGSetPlotGrids( void *data, JXF_Int plotgrids )
{
  JXF_Int ierr = 0;
  jxf_ParAMGData *amg_data = data;

  jxf_ParAMGDataPlotGrids(amg_data) = plotgrids;
  return (ierr);
}

JXF_Int
jxf_PAMGSetPlotFileName( void *data, const char *plot_file_name )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if( strlen(plot_file_name) > 251 )
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   if (strlen(plot_file_name) == 0 )
   {
      jxf_sprintf(jxf_ParAMGDataPlotFileName(amg_data), "%s", "AMGgrids.CF.dat");
   }
   else
   {
      jxf_sprintf(jxf_ParAMGDataPlotFileName(amg_data), "%s", plot_file_name);
   }

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetCoordDim( void *data, JXF_Int coorddim )
{
   JXF_Int ierr = 0;
   jxf_ParAMGData *amg_data = data;

   jxf_ParAMGDataCoordDim(amg_data) = coorddim;
   return (ierr);
}

JXF_Int
jxf_PAMGSetCoordinates( void *data, float *coordinates )
{
   JXF_Int ierr = 0;
   jxf_ParAMGData *amg_data = data;

   jxf_ParAMGDataCoordinates(amg_data) = coordinates;
   return (ierr);
}



/*--------------------------------------------------------------------------
 * Routines to set the problem data parameters
 *--------------------------------------------------------------------------*/

JXF_Int
jxf_PAMGSetNumFunctions( void *data, JXF_Int num_functions )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (num_functions < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNumFunctions(amg_data) = num_functions;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetNumFunctions( void *data, JXF_Int *num_functions )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *num_functions = jxf_ParAMGDataNumFunctions(amg_data);

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetNodal
 * \note Indicate whether to use nodal systems function.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetNodal( void *data, JXF_Int nodal )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNodal(amg_data) = nodal;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetNodalDiag
 * \note Indicate how to treat diag for primary matrix 
 *       with nodal systems function.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetNodalDiag( void *data, JXF_Int nodal )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNodalDiag(amg_data) = nodal;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetNumPaths
 * \note Indicate the degree of aggressive coarsening.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetNumPaths( void *data, JXF_Int num_paths )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (num_paths < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNumPaths(amg_data) = num_paths;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggNumLevels
 * \note Indicates the number of levels of aggressive coarsening.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetAggNumLevels( void *data, JXF_Int agg_num_levels )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (agg_num_levels < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataAggNumLevels(amg_data) = agg_num_levels;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggInterpType
 * \note Indicates the interpolation type of aggressive coarsening.
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGSetAggInterpType( void *data, JXF_Int agg_interp_type )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAMGDataAggInterpType(amg_data) = agg_interp_type;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggPMaxElmts
 * \note Indicates the max. number of entries in interpolation of aggressive coarsening.
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGSetAggPMaxElmts( void *data, JXF_Int agg_P_max_elmts )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (agg_P_max_elmts < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAMGDataAggPMaxElmts(amg_data) = agg_P_max_elmts;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggP12MaxElmts
 * \note Indicates the max. number of entries in 12 interpolation of aggressive coarsening.
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGSetAggP12MaxElmts( void *data, JXF_Int agg_P12_max_elmts )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   if (agg_P12_max_elmts < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   }
   jxf_ParAMGDataAggP12MaxElmts(amg_data) = agg_P12_max_elmts;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggTruncFactor
 * \note Indicates the truncated factor in interpolation of aggressive coarsening.
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGSetAggTruncFactor( void *data, JXF_Real agg_trunc_factor )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAMGDataAggTruncFactor(amg_data) = agg_trunc_factor;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetAggP12TruncFactor
 * \note Indicates the truncated factor in 12 interpolation of aggressive coarsening.
 * \date 2015/08/14
 */
JXF_Int
jxf_PAMGSetAggP12TruncFactor( void *data, JXF_Real agg_P12_trunc_factor )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAMGDataAggP12TruncFactor(amg_data) = agg_P12_trunc_factor;
   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetNumCRRelaxSteps
 * \note Indicates the number of relaxation steps for Compatible relaxation.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetNumCRRelaxSteps( void *data, JXF_Int num_CR_relax_steps )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (num_CR_relax_steps < 1)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNumCRRelaxSteps(amg_data) = num_CR_relax_steps;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetCRRate
 * \note Indicates the desired convergence rate for Compatible relaxation.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetCRRate( void *data, JXF_Real CR_rate )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCRRate(amg_data) = CR_rate;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetCRStrongTh
 * \note Indicates the desired convergence rate for Compatible relaxation.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetCRStrongTh( void *data, JXF_Real CR_strong_th )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCRStrongTh(amg_data) = CR_strong_th;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetISType
 * \note Indicates which independent set algorithm is used for CR.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetISType( void *data, JXF_Int IS_type )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (IS_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataISType(amg_data) = IS_type;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetCRUseCG 
 * \note Indicates whether to use CG for compatible relaxation.
 * \date 2010/05/29
 */
JXF_Int
jxf_PAMGSetCRUseCG( void *data, JXF_Int CR_use_CG )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCRUseCG(amg_data) = CR_use_CG;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetNumPoints( void *data, JXF_Int num_points )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataNumPoints(amg_data) = num_points;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetDofFunc( void *data, JXF_Int *dof_func )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_TFree(jxf_ParAMGDataDofFunc(amg_data));
   jxf_ParAMGDataDofFunc(amg_data) = dof_func;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetPointDofMap( void *data, JXF_Int *point_dof_map )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_TFree(jxf_ParAMGDataPointDofMap(amg_data));
   jxf_ParAMGDataPointDofMap(amg_data) = point_dof_map;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetDofPoint( void *data, JXF_Int *dof_point )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_TFree(jxf_ParAMGDataDofPoint(amg_data));
   jxf_ParAMGDataDofPoint(amg_data) = dof_point;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetNumLevels( void *data, JXF_Int *num_levels )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
  *num_levels = jxf_ParAMGDataNumLevels(amg_data);
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetNumIterations( void *data, JXF_Int *num_iterations )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *num_iterations = jxf_ParAMGDataNumIterations(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetCumNumIterations( void *data, JXF_Int *cum_num_iterations )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
#ifdef CUMNUMIT
   *cum_num_iterations = jxf_ParAMGDataCumNumIterations(amg_data);
#endif

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetResidual( void *data, jxf_ParVector **resid )
{
   jxf_ParAMGData  *amg_data = data;
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *resid = jxf_ParAMGDataResidual( amg_data );
   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetRelResidualNorm( void *data, JXF_Real *rel_resid_norm )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *rel_resid_norm = jxf_ParAMGDataRelativeResidualNorm(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetVariant( void *data, JXF_Int variant )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (variant < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataVariant(amg_data) = variant;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetVariant( void *data, JXF_Int *variant )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *variant = jxf_ParAMGDataVariant(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetOverlap( void *data, JXF_Int overlap )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (overlap < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataOverlap(amg_data) = overlap;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetOverlap( void *data, JXF_Int *overlap )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *overlap = jxf_ParAMGDataOverlap(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetDomainType( void *data, JXF_Int domain_type )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (domain_type < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataDomainType(amg_data) = domain_type;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetDomainType( void *data, JXF_Int *domain_type )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *domain_type = jxf_ParAMGDataDomainType(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSchwarzRlxWeight( void *data, JXF_Real schwarz_rlx_weight )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataSchwarzRlxWeight(amg_data) = schwarz_rlx_weight;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGGetSchwarzRlxWeight( void *data, JXF_Real *schwarz_rlx_weight )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   *schwarz_rlx_weight = jxf_ParAMGDataSchwarzRlxWeight(amg_data);

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSchwarzUseNonSymm( void *data, JXF_Int use_nonsymm )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataSchwarzUseNonSymm(amg_data) = use_nonsymm;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetSym( void *data, JXF_Int sym )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataSym(amg_data) = sym;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetLevel( void *data, JXF_Int level )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataLevel(amg_data) = level;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetThreshold( void *data, JXF_Real thresh )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataThreshold(amg_data) = thresh;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetFilter( void *data, JXF_Real filter )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataFilter(amg_data) = filter;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetDropTol( void *data, JXF_Real drop_tol )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataDropTol(amg_data) = drop_tol;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetMaxNzPerRow( void *data, JXF_Int max_nz_per_row )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   if (max_nz_per_row < 0)
   {
      jxf_error_in_arg(2);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataMaxNzPerRow(amg_data) = max_nz_per_row;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetEuclidFile( void *data, char *euclidfile )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataEuclidFile(amg_data) = euclidfile;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetEuLevel( void *data, JXF_Int eu_level )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataEuLevel(amg_data) = eu_level;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetEuSparseA( void *data, JXF_Real eu_sparse_A )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataEuSparseA(amg_data) = eu_sparse_A;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetEuBJ( void *data, JXF_Int eu_bj )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataEuBJ(amg_data) = eu_bj;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetConvCriteria
 * \author peghoty
 * \date 2009/07/25
 */
JXF_Int
jxf_PAMGSetConvCriteria( void *data, JXF_Int conv_criteria )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataConvCriteria(amg_data) = conv_criteria;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetConvCriteria
 * \note The option of "coarsestsolverid" is similar to "relax_type", 
 *       except that coarsestsolverid could be set to 10 (gselim_piv)
 * \author peghoty
 * \date 2009/07/27
 */
JXF_Int
jxf_PAMGSetCoarsestSolverID( void *data, JXF_Int coarsestsolverid )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCoarsestSolverID(amg_data) = coarsestsolverid;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetCoarseThreshold
 * \author peghoty
 * \date 2010/04/14
 */
JXF_Int
jxf_PAMGSetCoarseThreshold( void *data, JXF_Int coarse_threshold )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCoarseThreshold(amg_data) = coarse_threshold;

   return jxf_error_flag;
}

/*!
 * \fn JXF_Int jxf_PAMGSetCoarseRatio
 * \author peghoty
 * \date 2010/04/14
 */ 
JXF_Int
jxf_PAMGSetCoarseRatio( void *data, JXF_Real coarse_ratio )
{
   jxf_ParAMGData  *amg_data = data;
 
   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   } 
   jxf_ParAMGDataCoarseRatio(amg_data) = coarse_ratio;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetConvFacThreshold( void *data, JXF_Real convfac_threshold )
{
   jxf_ParAMGData  *amg_data = data;

   if (!amg_data)
   {
      jxf_printf("Warning! PAMG object empty!\n");
      jxf_error_in_arg(1);
      return jxf_error_flag;
   }
   jxf_ParAMGDataConvFacThreshold(amg_data) = convfac_threshold;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetRAP2( void *data, JXF_Int rap2 )
{
   jxf_ParAMGData *amg_data = data;

   jxf_ParAMGDataRAP2(amg_data) = rap2;

   return jxf_error_flag;
}

JXF_Int
jxf_PAMGSetKeepTranspose( void *data, JXF_Int keepTranspose )
{
   jxf_ParAMGData *amg_data = data;

   jxf_ParAMGDataKeepTranspose(amg_data) = keepTranspose;

   return jxf_error_flag;
}
