//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  euclid_dh.c
 *  Date: 2013/01/21
 */

#include "jxf_euclid.h"

static void jxf_get_runtime_params_private( jxf_Euclid_dh ctx );
static void jxf_invert_diagonals_private( jxf_Euclid_dh ctx );
static void jxf_compute_rho_private( jxf_Euclid_dh ctx );
static void jxf_factor_private( jxf_Euclid_dh ctx );
static void jxf_reduce_timings_private( jxf_Euclid_dh ctx );

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhCreate"
void jxf_Euclid_dhCreate( jxf_Euclid_dh *ctxOUT )
{
    JXF_START_FUNC_DH
    struct _jxf_mpi_interface_dh *ctx =
         (struct _jxf_mpi_interface_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_mpi_interface_dh)); JXF_CHECK_V_ERROR;
    
   *ctxOUT = ctx;
    ctx->isSetup = jxf_false;
    ctx->rho_init = 2.0;
    ctx->rho_final = 0.0;
    ctx->m = 0;
    ctx->n = 0;
    ctx->rhs = NULL;
    ctx->A = NULL;
    ctx->F = NULL;
    ctx->sg = NULL;
    ctx->scale = NULL;
    ctx->isScaled = jxf_false;
    ctx->work = NULL;
    ctx->work2 = NULL;
    ctx->from = 0;
    ctx->to = 0;
    strcpy(ctx->algo_par, "pilu");
    strcpy(ctx->algo_ilu, "iluk");
    ctx->level = 1;
    ctx->droptol = JXF_DEFAULT_DROP_TOL;
    ctx->sparseTolA = 0.0;
    ctx->sparseTolF = 0.0;
    ctx->pivotMin = 0.0;
    ctx->pivotFix = JXF_PIVOT_FIX_DEFAULT;
    ctx->maxVal = 0.0;
    ctx->slist = NULL;
    ctx->extRows = NULL;
    strcpy(ctx->krylovMethod, "bicgstab");
    ctx->maxIts = 200;
    ctx->rtol = 1e-5;
    ctx->atol = 1e-50;
    ctx->its = 0;
    ctx->itsTotal = 0;
    ctx->setupCount = 0;
    ctx->logging = 0;
    ctx->printStats = (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-printStats"));
    JXF_Int i;
    for (i = 0; i < JXF_TIMING_BINS; ++ i) ctx->timing[i] = 0.0;
    for (i = 0; i < JXF_STATS_BINS; ++ i) ctx->stats[i] = 0.0;
    ctx->timingsWereReduced = jxf_false;
    ++ jxf_ref_counter;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhDestroy"
void jxf_Euclid_dhDestroy( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-eu_stats") || ctx->logging)
    {
        jxf_Parser_dhInsert(jxf_parser_dh, "-eu_mem", "1"); JXF_CHECK_V_ERROR;
        jxf_Euclid_dhPrintJxpamgReport(ctx, stdout); JXF_CHECK_V_ERROR;
    }
    if (ctx->setupCount > 1 && ctx->printStats)
    {
        jxf_Euclid_dhPrintStatsShorter(ctx, stdout); JXF_CHECK_V_ERROR;
    }
    if (ctx->F != NULL)
    {
        jxf_Factor_dhDestroy(ctx->F); JXF_CHECK_V_ERROR;
    }
    if (ctx->sg != NULL)
    {
        jxf_SubdomainGraph_dhDestroy(ctx->sg); JXF_CHECK_V_ERROR;
    }
    if (ctx->scale != NULL)
    {
        JXF_FREE_DH(ctx->scale); JXF_CHECK_V_ERROR;
    }
    if (ctx->work != NULL)
    {
        JXF_FREE_DH(ctx->work); JXF_CHECK_V_ERROR;
    }
    if (ctx->work2 != NULL)
    {
        JXF_FREE_DH(ctx->work2); JXF_CHECK_V_ERROR;
    }
    if (ctx->slist != NULL)
    {
        jxf_SortedList_dhDestroy(ctx->slist); JXF_CHECK_V_ERROR;
    }
    if (ctx->extRows != NULL)
    {
        jxf_ExternalRows_dhDestroy(ctx->extRows); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(ctx); JXF_CHECK_V_ERROR;
    -- jxf_ref_counter;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhSetup"
void jxf_Euclid_dhSetup( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_Int m, n, beg_row;
    JXF_Real t1;
    jxf_bool isSetup = ctx->isSetup;
    jxf_bool bj = jxf_false;
    
    if (ctx->setupCount && ctx->printStats)
    {
        jxf_Euclid_dhPrintStatsShorter(ctx, stdout); JXF_CHECK_V_ERROR;
        ctx->its = 0;
    }
    JXF_Int i;
    for (i = 0; i < JXF_STATS_BINS; ++ i) ctx->stats[i] = 0.0;
    ctx->timing[JXF_SOLVE_START_T] = jxf_MPI_Wtime();
    ctx->timing[JXF_TOTAL_SOLVE_T] += ctx->timing[JXF_TOTAL_SOLVE_TEMP_T];
    ctx->timing[JXF_TOTAL_SOLVE_TEMP_T] = 0.0;
    if (ctx->F != NULL)
    {
        jxf_Factor_dhDestroy(ctx->F); JXF_CHECK_V_ERROR;
        ctx->F = NULL;
    }
    if (ctx->A == NULL)
    {
        JXF_SET_V_ERROR("must set ctx->A before calling init");
    }
    jxf_EuclidGetDimensions(ctx->A, &beg_row, &m, &n); JXF_CHECK_V_ERROR;
    ctx->m = m;
    ctx->n = n;
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-print_size"))
    {
        jxf_printf_dh("setting up linear system; global rows: %i  local rows: %i (on P_0)\n", n,m);
    }
    jxf_sprintf(jxf_msgBuf_dh, "localRow= %i;  globalRows= %i;  beg_row= %i", m, n, beg_row);
    JXF_SET_INFO(jxf_msgBuf_dh);
    bj = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-bj");
    if (ctx->sg == NULL)
    {
        JXF_Int blocks = jxf_np_dh;
        t1 = jxf_MPI_Wtime();
        if (jxf_np_dh == 1)
        {
            jxf_Parser_dhReadInt(jxf_parser_dh, "-blocks", &blocks); JXF_CHECK_V_ERROR;
            jxf_SubdomainGraph_dhCreate(&(ctx->sg)); JXF_CHECK_V_ERROR;
            jxf_SubdomainGraph_dhInit(ctx->sg, blocks, bj, ctx->A); JXF_CHECK_V_ERROR;
        }
        else
        {
            jxf_SubdomainGraph_dhCreate(&(ctx->sg)); JXF_CHECK_V_ERROR;
            jxf_SubdomainGraph_dhInit(ctx->sg, -1, bj, ctx->A); JXF_CHECK_V_ERROR;
        }
        ctx->timing[JXF_SUB_GRAPH_T] += (jxf_MPI_Wtime() - t1);
    }
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-doNotFactor"))
    {
        goto END_OF_FUNCTION;
    }
    if (!isSetup)
    {
        jxf_get_runtime_params_private(ctx); JXF_CHECK_V_ERROR;
    }
    if (!strcmp(ctx->algo_par, "bj")) bj = jxf_false;
    if (ctx->scale == NULL)
    {
        ctx->scale = (JXF_REAL_DH*)JXF_MALLOC_DH(m*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
    }
    for (i = 0; i < m; ++ i) ctx->scale[i] = 1.0;
    if (ctx->work == NULL)
    {
        ctx->work = (JXF_REAL_DH*)JXF_MALLOC_DH(m*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
    }
    if (ctx->work2 == NULL)
    {
        ctx->work2 = (JXF_REAL_DH*)JXF_MALLOC_DH(m*sizeof(JXF_REAL_DH)); JXF_CHECK_V_ERROR;
    }
    t1 = jxf_MPI_Wtime();
    jxf_factor_private(ctx); JXF_CHECK_V_ERROR;
    ctx->timing[JXF_FACTOR_T] += (jxf_MPI_Wtime() - t1);
    if (strcmp(ctx->algo_par, "none"))
    {
        jxf_invert_diagonals_private(ctx); JXF_CHECK_V_ERROR;
    }
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-computeRho") || jxf_np_dh == 1)
    {
        if (strcmp(ctx->algo_par, "none"))
        {
            t1 = jxf_MPI_Wtime();
            jxf_compute_rho_private(ctx); JXF_CHECK_V_ERROR;
            ctx->timing[JXF_COMPUTE_RHO_T] += (jxf_MPI_Wtime() - t1);
        }
    }
    if (!strcmp(ctx->algo_par, "pilu") && jxf_np_dh > 1)
    {
        t1 = jxf_MPI_Wtime();
        jxf_Factor_dhSolveSetup(ctx->F, ctx->sg); JXF_CHECK_V_ERROR;
        ctx->timing[JXF_SOLVE_SETUP_T] += (jxf_MPI_Wtime() - t1);
    }
    
END_OF_FUNCTION: ;
    
    ctx->timing[JXF_SETUP_T] += (jxf_MPI_Wtime() - ctx->timing[JXF_SOLVE_START_T]);
    ctx->setupCount += 1;
    ctx->isSetup = jxf_true;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_get_runtime_params_private"
void jxf_get_runtime_params_private( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    char *tmp;
    
    jxf_Parser_dhReadInt(jxf_parser_dh, "-maxIts",&(ctx->maxIts));
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-rtol", &(ctx->rtol));
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-atol", &(ctx->atol));
    tmp = NULL;
    jxf_Parser_dhReadString(jxf_parser_dh, "-par", &tmp);
    if (tmp != NULL)
    {
        strcpy(ctx->algo_par, tmp);
    }
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-bj"))
    {
        strcpy(ctx->algo_par, "bj");
    }
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-rho", &(ctx->rho_init));
    jxf_Parser_dhReadInt(jxf_parser_dh, "-level", &ctx->level);
    jxf_Parser_dhReadInt(jxf_parser_dh, "-pc_ilu_levels", &ctx->level);
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-ilut"))
    {
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-ilut", &ctx->droptol);
        ctx->isScaled = jxf_true;
        strcpy(ctx->algo_ilu, "ilut");
    }
    if (!strcmp(ctx->algo_par, "none"))
    {
        strcpy(ctx->algo_ilu, "none");
    }
    else if (!strcmp(ctx->algo_ilu, "none"))
    {
        strcpy(ctx->algo_par, "none");
    }
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-sparseA",&(ctx->sparseTolA));
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-sparseF",&(ctx->sparseTolF));
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-pivotMin", &(ctx->pivotMin));
    jxf_Parser_dhReadDouble(jxf_parser_dh, "-pivotFix", &(ctx->pivotFix));
    if (ctx->sparseTolA || ! strcmp(ctx->algo_ilu, "ilut"))
    {
        ctx->isScaled = jxf_true;
    }
    tmp = NULL;
    jxf_Parser_dhReadString(jxf_parser_dh, "-ksp_type", &tmp);
    if (tmp != NULL)
    {
        strcpy(ctx->krylovMethod, tmp);
        if (!strcmp(ctx->krylovMethod, "bcgs"))
        {
            strcpy(ctx->krylovMethod, "bicgstab");
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_invert_diagonals_private"
void jxf_invert_diagonals_private( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    JXF_REAL_DH *aval = ctx->F->aval;
    JXF_Int *diag = ctx->F->diag;
    
    if (aval == NULL || diag == NULL)
    {
        JXF_SET_INFO("can't invert diags; either F->aval or F->diag is NULL");
    }
    else
    {
        JXF_Int i, m = ctx->F->m;
        for (i = 0; i < m; ++ i)
        {
            aval[diag[i]] = 1.0 / aval[diag[i]];
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_compute_rho_private"
void jxf_compute_rho_private( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    if (ctx->F != NULL)
    {
        JXF_Real bufLocal[3], bufGlobal[3];
        JXF_Int m = ctx->m;
        
        ctx->stats[JXF_NZF_STATS] = (JXF_Real)ctx->F->rp[m];
        bufLocal[0] = ctx->stats[JXF_NZA_STATS];
        bufLocal[1] = ctx->stats[JXF_NZF_STATS];
        bufLocal[2] = ctx->stats[JXF_NZA_USED_STATS];
        if (jxf_np_dh == 1)
        {
            bufGlobal[0] = bufLocal[0];
            bufGlobal[1] = bufLocal[1];
            bufGlobal[2] = bufLocal[2];
        }
        else
        {
            jxf_MPI_Reduce(bufLocal, bufGlobal, 3, JXF_MPI_REAL, MPI_SUM, 0, jxf_comm_dh);
        }
        if (jxf_myid_dh == 0)
        {
            if (bufGlobal[0] && bufGlobal[1])
            {
                ctx->rho_final = bufGlobal[1] / bufGlobal[0];
            }
            else
            {
                ctx->rho_final = -1;
            }
            if (bufGlobal[0] && bufGlobal[2])
            {
                ctx->stats[JXF_NZA_RATIO_STATS] = 100.0 * bufGlobal[2] / bufGlobal[0];
            }
            else
            {
                ctx->stats[JXF_NZA_RATIO_STATS] = 100.0;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_factor_private"
void jxf_factor_private( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    if (!strcmp(ctx->algo_par, "none"))
    {
        goto DO_NOTHING;
    }
    JXF_Int br = 0;
    JXF_Int id = jxf_np_dh;
    if (ctx->sg != NULL)
    {
        br = ctx->sg->beg_rowP[jxf_myid_dh];
        id = ctx->sg->o2n_sub[jxf_myid_dh];
    }
    jxf_Factor_dhInit(ctx->A, jxf_true, jxf_true, ctx->rho_init, id, br, &(ctx->F)); JXF_CHECK_V_ERROR;
    ctx->F->bdry_count = ctx->sg->bdry_count[jxf_myid_dh];
    ctx->F->first_bdry = ctx->F->m - ctx->F->bdry_count;
    if (!strcmp(ctx->algo_par, "bj")) ctx->F->blockJacobi = jxf_true;
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-bj")) ctx->F->blockJacobi = jxf_true;
    if (jxf_np_dh == 1)
    {
        if (!strcmp(ctx->algo_ilu, "iluk"))
        {
            ctx->from = 0;
            ctx->to = ctx->m;
            if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-mpi"))
            {
                if (ctx->sg != NULL && ctx->sg->blocks > 1)
                {
                    JXF_SET_V_ERROR("only use -mpi, which invokes jxf_ilu_mpi_pilu(), for np = 1 and -blocks 1");
                }
                jxf_iluk_mpi_pilu(ctx); JXF_CHECK_V_ERROR;
            }
            else
            {
                jxf_iluk_seq_block(ctx); JXF_CHECK_V_ERROR;
            }
        }
        else if (!strcmp(ctx->algo_ilu, "ilut"))
        {
            ctx->from = 0;
            ctx->to = ctx->m;
            jxf_ilut_seq(ctx); JXF_CHECK_V_ERROR;
        }
        else
        {
            jxf_sprintf(jxf_msgBuf_dh, "factorization method: %s is not implemented", ctx->algo_ilu);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    else
    {
        if (!strcmp(ctx->algo_par, "bj"))
        {
            ctx->from = 0;
            ctx->to = ctx->m;
            jxf_iluk_mpi_bj(ctx); JXF_CHECK_V_ERROR;
        }
        else if (!strcmp(ctx->algo_ilu, "iluk"))
        {
            jxf_bool bj = ctx->F->blockJacobi;
            jxf_SortedList_dhCreate(&(ctx->slist)); JXF_CHECK_V_ERROR;
            jxf_SortedList_dhInit(ctx->slist, ctx->sg); JXF_CHECK_V_ERROR;
            jxf_ExternalRows_dhCreate(&(ctx->extRows)); JXF_CHECK_V_ERROR;
            jxf_ExternalRows_dhInit(ctx->extRows, ctx); JXF_CHECK_V_ERROR;
            ctx->from = 0;
            ctx->to = ctx->F->first_bdry;
            jxf_iluk_seq(ctx); JXF_CHECK_V_ERROR;
            if (!bj)
            {
                jxf_ExternalRows_dhRecvRows(ctx->extRows); JXF_CHECK_V_ERROR;
            }
            ctx->from = ctx->F->first_bdry;
            ctx->to = ctx->F->m;
            jxf_iluk_mpi_pilu(ctx); JXF_CHECK_V_ERROR;
            if (!bj)
            {
                jxf_ExternalRows_dhSendRows(ctx->extRows); JXF_CHECK_V_ERROR;
            }
            jxf_SortedList_dhDestroy(ctx->slist); JXF_CHECK_V_ERROR;
            ctx->slist = NULL;
            jxf_ExternalRows_dhDestroy(ctx->extRows); JXF_CHECK_V_ERROR;
            ctx->extRows = NULL;
        }
        else
        {
            jxf_sprintf(jxf_msgBuf_dh, "factorization method: %s is not implemented", ctx->algo_ilu);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    
DO_NOTHING: ;
    
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhSolve"
void jxf_Euclid_dhSolve( jxf_Euclid_dh ctx, jxf_Vec_dh x, jxf_Vec_dh b, JXF_Int *its )
{
    JXF_START_FUNC_DH
    JXF_Int itsOUT;
    jxf_Mat_dh A = (jxf_Mat_dh)ctx->A;
    
    if (!strcmp(ctx->krylovMethod, "cg"))
    {
        jxf_cg_euclid(A, ctx, x->vals, b->vals, &itsOUT); JXF_ERRCHKA;
    }
    else if (!strcmp(ctx->krylovMethod, "bicgstab"))
    {
        jxf_bicgstab_euclid(A, ctx, x->vals, b->vals, &itsOUT); JXF_ERRCHKA;
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "unknown krylov solver: %s", ctx->krylovMethod);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
   *its = itsOUT;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintStats"
void jxf_Euclid_dhPrintStats( jxf_Euclid_dh ctx, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Real *timing;
    JXF_Int nz;
    
    nz = jxf_Factor_dhReadNz(ctx->F); JXF_CHECK_V_ERROR;
    timing = ctx->timing;
    ctx->timing[JXF_TOTAL_SOLVE_T] += ctx->timing[JXF_TOTAL_SOLVE_TEMP_T];
    ctx->timing[JXF_TOTAL_SOLVE_TEMP_T] = 0.0;
    jxf_reduce_timings_private(ctx); JXF_CHECK_V_ERROR;
    jxf_fprintf_dh(fp, "\n==================== Euclid report (start) ====================\n");
    jxf_fprintf_dh(fp, "\nruntime parameters\n");
    jxf_fprintf_dh(fp, "------------------\n");
    jxf_fprintf_dh(fp, "   setups:                 %i\n", ctx->setupCount);
    jxf_fprintf_dh(fp, "   tri solves:             %i\n", ctx->itsTotal);
    jxf_fprintf_dh(fp, "   parallelization method: %s\n", ctx->algo_par);
    jxf_fprintf_dh(fp, "   factorization method:   %s\n", ctx->algo_ilu);
    jxf_fprintf_dh(fp, "   matrix was row scaled:  %i\n", ctx->isScaled);
    jxf_fprintf_dh(fp, "   matrix row count:       %i\n", ctx->n);
    jxf_fprintf_dh(fp, "   nzF:                    %i\n", nz);
    jxf_fprintf_dh(fp, "   rho:                    %g\n", ctx->rho_final);
    jxf_fprintf_dh(fp, "   level:                  %i\n", ctx->level);
    jxf_fprintf_dh(fp, "   sparseA:                %g\n", ctx->sparseTolA);
    jxf_fprintf_dh(fp, "\nEuclid timing report\n");
    jxf_fprintf_dh(fp, "--------------------\n");
    jxf_fprintf_dh(fp, "   solves total:  %0.2f (see docs)\n", timing[JXF_TOTAL_SOLVE_T]);
    jxf_fprintf_dh(fp, "   tri solves:    %0.2f\n", timing[JXF_TRI_SOLVE_T]);
    jxf_fprintf_dh(fp, "   setups:        %0.2f\n", timing[JXF_SETUP_T]);
    jxf_fprintf_dh(fp, "      subdomain graph setup:  %0.2f\n", timing[JXF_SUB_GRAPH_T]);
    jxf_fprintf_dh(fp, "      factorization:          %0.2f\n", timing[JXF_FACTOR_T]);
    jxf_fprintf_dh(fp, "      solve setup:            %0.2f\n", timing[JXF_SOLVE_SETUP_T]);
    jxf_fprintf_dh(fp, "      rho:                    %0.2f\n", ctx->timing[JXF_COMPUTE_RHO_T]);
    jxf_fprintf_dh(fp, "      misc (should be small): %0.2f\n", 
           timing[JXF_SETUP_T]-(timing[JXF_SUB_GRAPH_T]+timing[JXF_FACTOR_T]+timing[JXF_SOLVE_SETUP_T]+timing[JXF_COMPUTE_RHO_T]));
    if (ctx->sg != NULL)
    {
        jxf_SubdomainGraph_dhPrintStats(ctx->sg, fp); JXF_CHECK_V_ERROR;
        jxf_SubdomainGraph_dhPrintRatios(ctx->sg, fp); JXF_CHECK_V_ERROR;
    }
    jxf_fprintf_dh(fp, "\nApplicable if Euclid's internal solvers were used:\n");
    jxf_fprintf_dh(fp, "---------------------------------------------------\n");
    jxf_fprintf_dh(fp, "   solve method: %s\n", ctx->krylovMethod);
    jxf_fprintf_dh(fp, "   maxIts:       %i\n", ctx->maxIts);
    jxf_fprintf_dh(fp, "   rtol:         %g\n", ctx->rtol);
    jxf_fprintf_dh(fp, "   atol:         %g\n", ctx->atol);
    jxf_fprintf_dh(fp, "\n==================== Euclid report (end) ======================\n");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintStatsShort"
void jxf_Euclid_dhPrintStatsShort( jxf_Euclid_dh ctx, JXF_Real setup, JXF_Real solve, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Real *timing = ctx->timing;
    JXF_Real apply_total;
    JXF_Real apply_per_it;
    JXF_Real perIt;
    JXF_Int blocks = jxf_np_dh;
    
    if (jxf_np_dh == 1) blocks = ctx->sg->blocks;
    jxf_reduce_timings_private(ctx); JXF_CHECK_V_ERROR;
    apply_total = timing[JXF_TRI_SOLVE_T];
    apply_per_it = apply_total/(JXF_Real)ctx->its;
    perIt = solve / (JXF_Real)ctx->its;
    jxf_fprintf_dh(fp, "\n");
    jxf_fprintf_dh(fp, "%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s XX\n",
                   "method", "subdms", "level", "its", "setup", "solve", "total", "perIt", "perIt", "rows");
    jxf_fprintf_dh(fp, "------  -----  -----  -----  -----  -----  -----  -----  -----  -----  XX\n");
    jxf_fprintf_dh(fp, "%6s %6i %6i %6i %6.2f %6.2f %6.2f %6.4f %6.5f %6g  XXX\n",
      ctx->algo_par, blocks, ctx->level, ctx->its, setup, solve, setup+solve, perIt, apply_per_it, (JXF_Real)ctx->n);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintStatsShorter"
void jxf_Euclid_dhPrintStatsShorter( jxf_Euclid_dh ctx, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Real *stats = ctx->stats;
    JXF_Int its = ctx->its;
    JXF_Real rho = ctx->rho_final;
    JXF_Real nzUsedRatio = stats[JXF_NZA_RATIO_STATS];
    
    jxf_fprintf_dh(fp, "\nStats from last linear solve: YY\n");
    jxf_fprintf_dh(fp, "%6s %6s %6s     YY\n", "its", "rho","A_%");
    jxf_fprintf_dh(fp, " -----  -----  -----     YY\n");
    jxf_fprintf_dh(fp, "%6i %6.2f %6.2f     YYY\n", its, rho, nzUsedRatio);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintScaling"
void jxf_Euclid_dhPrintScaling( jxf_Euclid_dh ctx, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = ctx->m;
    
    if (m > 10) m = 10;
    if (ctx->scale == NULL)
    {
        JXF_SET_V_ERROR("ctx->scale is NULL; was jxf_Euclid_dhSetup() called?");
    }
    jxf_fprintf(fp, "\n---------- 1st %i row scaling values:\n", m);
    for (i = 0; i < m; ++ i)
    {
        jxf_fprintf(fp, "   %i  %g  \n", i+1, ctx->scale[i]);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_reduce_timings_private"
void jxf_reduce_timings_private( jxf_Euclid_dh ctx )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh > 1)
    {
        JXF_Real bufOUT[JXF_TIMING_BINS];
        
        memcpy(bufOUT, ctx->timing, JXF_TIMING_BINS*sizeof(JXF_Real));
        jxf_MPI_Reduce(bufOUT, ctx->timing, JXF_TIMING_BINS, JXF_MPI_REAL, MPI_MAX, 0, jxf_comm_dh);
    }
    ctx->timingsWereReduced = jxf_true;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintJxpamgReport"
void jxf_Euclid_dhPrintJxpamgReport( jxf_Euclid_dh ctx, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Real *timing;
    JXF_Int nz;
    
    nz = jxf_Factor_dhReadNz(ctx->F); JXF_CHECK_V_ERROR;
    timing = ctx->timing;
    ctx->timing[JXF_TOTAL_SOLVE_T] += ctx->timing[JXF_TOTAL_SOLVE_TEMP_T];
    ctx->timing[JXF_TOTAL_SOLVE_TEMP_T] = 0.0;
    jxf_reduce_timings_private(ctx); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        jxf_fprintf(fp, "@@@@@@@@@@@@@@@@@@@@@@ Euclid statistical report (start)\n");
        jxf_fprintf_dh(fp, "\nruntime parameters\n");
        jxf_fprintf_dh(fp, "------------------\n");
        jxf_fprintf_dh(fp, "   setups:                 %i\n", ctx->setupCount);
        jxf_fprintf_dh(fp, "   tri solves:             %i\n", ctx->itsTotal);
        jxf_fprintf_dh(fp, "   parallelization method: %s\n", ctx->algo_par);
        jxf_fprintf_dh(fp, "   factorization method:   %s\n", ctx->algo_ilu);
        if (!strcmp(ctx->algo_ilu, "iluk"))
        {
            jxf_fprintf_dh(fp, "      level:               %i\n", ctx->level);
        }
        if (ctx->isScaled)
        {
            jxf_fprintf_dh(fp, "   matrix was row scaled\n");
        }
        jxf_fprintf_dh(fp, "   global matrix row count: %i\n", ctx->n);
        jxf_fprintf_dh(fp, "   nzF:                     %i\n", nz);
        jxf_fprintf_dh(fp, "   rho:                     %g\n", ctx->rho_final);
        jxf_fprintf_dh(fp, "   sparseA:                 %g\n", ctx->sparseTolA);
        jxf_fprintf_dh(fp, "\nEuclid timing report\n");
        jxf_fprintf_dh(fp, "--------------------\n");
        jxf_fprintf_dh(fp, "   solves total:  %0.2f (see docs)\n", timing[JXF_TOTAL_SOLVE_T]);
        jxf_fprintf_dh(fp, "   tri solves:    %0.2f\n", timing[JXF_TRI_SOLVE_T]);
        jxf_fprintf_dh(fp, "   setups:        %0.2f\n", timing[JXF_SETUP_T]);
        jxf_fprintf_dh(fp, "      subdomain graph setup:  %0.2f\n", timing[JXF_SUB_GRAPH_T]);
        jxf_fprintf_dh(fp, "      factorization:          %0.2f\n", timing[JXF_FACTOR_T]);
        jxf_fprintf_dh(fp, "      solve setup:            %0.2f\n", timing[JXF_SOLVE_SETUP_T]);
        jxf_fprintf_dh(fp, "      rho:                    %0.2f\n", ctx->timing[JXF_COMPUTE_RHO_T]);
        jxf_fprintf_dh(fp, "      misc (should be small): %0.2f\n", 
           timing[JXF_SETUP_T]-(timing[JXF_SUB_GRAPH_T]+timing[JXF_FACTOR_T]+timing[JXF_SOLVE_SETUP_T]+timing[JXF_COMPUTE_RHO_T]));
        if (ctx->sg != NULL)
        {
            jxf_SubdomainGraph_dhPrintStats(ctx->sg, fp); JXF_CHECK_V_ERROR;
            jxf_SubdomainGraph_dhPrintRatios(ctx->sg, fp); JXF_CHECK_V_ERROR;
        }
        jxf_fprintf(fp, "@@@@@@@@@@@@@@@@@@@@@@ Euclid statistical report (end)\n");
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Euclid_dhPrintTestData"
void jxf_Euclid_dhPrintTestData( jxf_Euclid_dh ctx, FILE *fp )
{
    JXF_START_FUNC_DH
    if (jxf_myid_dh == 0)
    {
        jxf_fprintf(fp, "   setups:                 %i\n", ctx->setupCount);
        jxf_fprintf(fp, "   tri solves:             %i\n", ctx->its);
        jxf_fprintf(fp, "   parallelization method: %s\n", ctx->algo_par);
        jxf_fprintf(fp, "   factorization method:   %s\n", ctx->algo_ilu);
        jxf_fprintf(fp, "   level:                  %i\n", ctx->level);
        jxf_fprintf(fp, "   row scaling:            %i\n", ctx->isScaled);
    }
    jxf_SubdomainGraph_dhPrintRatios(ctx->sg, fp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}
