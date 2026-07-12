//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  parcsr_euclid.c
 *  Date: 2013/01/23
 */

#include "jxf_pamg.h"
#include "jxf_euclid.h"

#define JXF_EUCLID_ERRCHKA \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("", __FUNC__, __FILE__, __LINE__); \
        jxf_printErrorMsg(stderr); \
        jxf_MPI_Abort(jxf_comm_dh, -1); \
      }

#undef JXF_START_FUNC_DH
#undef JXF_END_FUNC_VAL
#undef JXF_END_FUNC_DH
#define JXF_START_FUNC_DH
#define JXF_END_FUNC_DH
#define JXF_END_FUNC_VAL(a) return(a);

#undef __FUNC__
#define __FUNC__ "JXF_EuclidCreate"
JXF_Int 
JXF_EuclidCreate( MPI_Comm comm, JXF_Solver *solver )
{
    JXF_START_FUNC_DH
    jxf_Euclid_dh eu;
    
    jxf_comm_dh = comm;
    jxf_MPI_Comm_size(jxf_comm_dh, &jxf_np_dh); JXF_EUCLID_ERRCHKA;
    jxf_MPI_Comm_rank(jxf_comm_dh, &jxf_myid_dh); JXF_EUCLID_ERRCHKA;
    
#ifdef ENABLE_EUCLID_LOGGING
    jxf_openLogfile_dh(0, NULL); JXF_EUCLID_ERRCHKA;
#endif
    
    if (jxf_mem_dh == NULL)
    {
        jxf_Mem_dhCreate(&jxf_mem_dh); JXF_EUCLID_ERRCHKA;
    }
    if (jxf_tlog_dh == NULL)
    {
        jxf_TimeLog_dhCreate(&jxf_tlog_dh); JXF_EUCLID_ERRCHKA;
    }
    if (jxf_parser_dh == NULL)
    {
        jxf_Parser_dhCreate(&jxf_parser_dh); JXF_EUCLID_ERRCHKA;
    }
    jxf_Parser_dhInit(jxf_parser_dh, 0, NULL); JXF_EUCLID_ERRCHKA;
    jxf_Euclid_dhCreate(&eu); JXF_EUCLID_ERRCHKA;
   *solver = (JXF_Solver) eu;
    JXF_END_FUNC_VAL(0)
}

#undef __FUNC__
#define __FUNC__ "JXF_EuclidDestroy"
JXF_Int 
JXF_EuclidDestroy( JXF_Solver solver )
{
    JXF_START_FUNC_DH
    jxf_Euclid_dh eu = (jxf_Euclid_dh)solver;
    jxf_bool printMemReport = jxf_false;
    jxf_bool printStats = jxf_false;
    jxf_bool logging = eu->logging;
    
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-printTestData"))
    {
        FILE *fp;
        char fname[] = "test_data_dh.temp", *fnamePtr = fname;
        
        jxf_Parser_dhReadString(jxf_parser_dh, "-printTestData", &fnamePtr); JXF_EUCLID_ERRCHKA;
        if (!strcmp(fnamePtr, "1"))
        {
            fnamePtr = fname;
        }
        fp = jxf_openFile_dh(fnamePtr, "w"); JXF_EUCLID_ERRCHKA;
        jxf_Euclid_dhPrintTestData(eu, fp); JXF_EUCLID_ERRCHKA;
        jxf_closeFile_dh(fp); JXF_EUCLID_ERRCHKA;
        jxf_printf_dh("\n@@@@@ Euclid test data was printed to file: %s\n\n", fnamePtr);
    }
    if (logging)
    {
        printStats = jxf_true;
        printMemReport = jxf_true;
    }
    if (jxf_parser_dh != NULL)
    {
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-eu_stats"))
        {
            printStats = jxf_true;
        }
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-eu_mem"))
        {
            printMemReport = jxf_true;
        }
    }
    if (printStats)
    {
        jxf_Euclid_dhPrintJxpamgReport(eu, stdout); JXF_EUCLID_ERRCHKA;
    }
    jxf_Euclid_dhDestroy(eu); JXF_EUCLID_ERRCHKA;
    if (jxf_parser_dh != NULL && jxf_ref_counter == 0)
    {
        jxf_Parser_dhDestroy(jxf_parser_dh); JXF_EUCLID_ERRCHKA;
        jxf_parser_dh = NULL;
    }
    if (jxf_tlog_dh != NULL && jxf_ref_counter == 0)
    {
        jxf_TimeLog_dhDestroy(jxf_tlog_dh); JXF_EUCLID_ERRCHKA;
        jxf_tlog_dh = NULL;
    }
    if (jxf_mem_dh != NULL && jxf_ref_counter == 0)
    {
        if (printMemReport)
        {
            jxf_Mem_dhPrint(jxf_mem_dh, stdout, jxf_false); JXF_EUCLID_ERRCHKA;
        }
        jxf_Mem_dhDestroy(jxf_mem_dh); JXF_EUCLID_ERRCHKA;
        jxf_mem_dh = NULL;
    }
    
#ifdef ENABLE_EUCLID_LOGGING
    jxf_closeLogfile_dh(); JXF_EUCLID_ERRCHKA;
#endif
    
    JXF_END_FUNC_VAL(0)
}

#undef __FUNC__
#define __FUNC__ "JXF_EuclidSetup"
JXF_Int 
JXF_EuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix A )
{
    JXF_START_FUNC_DH
    jxf_Euclid_dh eu = (jxf_Euclid_dh)solver;
    jxf_Euclid_dhInputJXFpamgMat(eu, (JXF_ParCSRMatrix)(jxf_hpCSRMatrixPar((jxf_hpCSRMatrix*) A))); 
    JXF_EUCLID_ERRCHKA;
    jxf_Euclid_dhSetup(eu); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

#undef __FUNC__
#define __FUNC__ "JXF_EuclidSolve"
JXF_Int 
JXF_EuclidSolve( JXF_Solver solver,
                JXF_hpCSRMatrix A,
                JXF_ParVector bb,
                JXF_ParVector xx )
{
        JXF_START_FUNC_DH
    jxf_Euclid_dh eu = (jxf_Euclid_dh)solver;
    JXF_Real *b, *x;
    
    x = jxf_VectorData(jxf_ParVectorLocalVector((jxf_ParVector *)bb));
    b = jxf_VectorData(jxf_ParVectorLocalVector((jxf_ParVector *)xx));
    jxf_Euclid_dhApply(eu, x, b); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

#undef __FUNC__
#define __FUNC__ "JXF_EuclidSetParams"
JXF_Int 
JXF_EuclidSetParams( JXF_Solver solver, JXF_Int argc, char *argv[] )
{
    JXF_START_FUNC_DH
    jxf_Parser_dhInit(jxf_parser_dh, argc, argv); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

#undef __FUNC__
#define __FUNC__ "JXF_EuclidSetParamsFromFile"
JXF_Int 
JXF_EuclidSetParamsFromFile( JXF_Solver solver, char *filename )
{
    JXF_START_FUNC_DH
    jxf_Parser_dhUpdateFromFile(jxf_parser_dh, filename); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetLevel( JXF_Solver solver, JXF_Int level )
{
    char str_level[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_level,"%d",level);
    jxf_Parser_dhInsert(jxf_parser_dh, "-level", str_level); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetBJ( JXF_Solver solver, JXF_Int bj )
{
    char str_bj[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_bj,"%d",bj);
    jxf_Parser_dhInsert(jxf_parser_dh, "-bj", str_bj); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetStats( JXF_Solver solver, JXF_Int eu_stats )
{
    char str_eu_stats[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_eu_stats,"%d",eu_stats);
    jxf_Parser_dhInsert(jxf_parser_dh, "-eu_stats", str_eu_stats); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetMem( JXF_Solver solver, JXF_Int eu_mem )
{
    char str_eu_mem[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_eu_mem,"%d",eu_mem);
    jxf_Parser_dhInsert(jxf_parser_dh, "-eu_mem", str_eu_mem); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetILUT( JXF_Solver solver, JXF_Real ilut )
{
    char str_ilut[256];
    JXF_START_FUNC_DH
    jxf_sprintf(str_ilut,"%f",ilut);
    jxf_Parser_dhInsert(jxf_parser_dh, "-ilut", str_ilut); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetSparseA( JXF_Solver solver, JXF_Real sparse_A )
{
    char str_sparse_A[256];
    JXF_START_FUNC_DH
    jxf_sprintf(str_sparse_A,"%f",sparse_A);
    jxf_Parser_dhInsert(jxf_parser_dh, "-sparseA", str_sparse_A); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}

JXF_Int 
JXF_EuclidSetRowScale( JXF_Solver solver, JXF_Int row_scale )
{
    char str_row_scale[8];
    JXF_START_FUNC_DH
    jxf_sprintf(str_row_scale,"%d",row_scale);
    jxf_Parser_dhInsert(jxf_parser_dh, "-rowScale", str_row_scale); JXF_EUCLID_ERRCHKA;
    JXF_END_FUNC_VAL(0)
}
