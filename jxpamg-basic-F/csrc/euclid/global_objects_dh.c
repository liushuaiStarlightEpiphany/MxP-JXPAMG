//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  global_objects_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

extern void jxf_sigRegister_dh();

jxf_bool jxf_errFlag_dh = jxf_false;
jxf_Parser_dh jxf_parser_dh = NULL;
jxf_TimeLog_dh jxf_tlog_dh = NULL;
jxf_Mem_dh jxf_mem_dh = NULL;
FILE *jxf_logFile = NULL;
char jxf_msgBuf_dh[JXF_MSG_BUF_SIZE_DH];
JXF_Int jxf_np_dh = 1;
JXF_Int jxf_myid_dh = 0;
MPI_Comm jxf_comm_dh = 0;
FILE *jxf_logFile;
void jxf_openLogfile_dh( JXF_Int argc, char *argv[] );
void jxf_closeLogfile_dh();
jxf_bool jxf_logInfoToStderr = jxf_false;
jxf_bool jxf_logInfoToFile = jxf_true;
jxf_bool jxf_logFuncsToStderr = jxf_false;
jxf_bool jxf_logFuncsToFile = jxf_false;
jxf_bool jxf_ignoreMe = jxf_true;
JXF_Int jxf_ref_counter = 0;

#define JXF_MAX_MSG_SIZE 1024
#define JXF_MAX_STACK_SIZE 20

static char jxf_errMsg_private[JXF_MAX_STACK_SIZE][JXF_MAX_MSG_SIZE];
static JXF_Int jxf_errCount_private = 0;
static char jxf_calling_stack[JXF_MAX_STACK_SIZE][JXF_MAX_MSG_SIZE];
static JXF_Int jxf_jxf_calling_stack_count = 0;

void jxf_openLogfile_dh( JXF_Int argc, char *argv[] )
{
    char buf[1024];
    
    if (jxf_logFile != NULL) return;
    jxf_sprintf(buf, "jxf_logFile");
    if (argc && argv != NULL)
    {
        JXF_Int j;
        
        for (j = 1; j < argc; ++ j)
        {
            if (strcmp(argv[j],"-jxf_logFile") == 0)
            {
                if (j+1 < argc)
                {
                    jxf_sprintf(buf, "%s", argv[j+1]);
                    break;
                }
            }
        }
    }
    if (strcmp(buf, "none"))
    {
        char a[5];
        
        jxf_sprintf(a, ".%i", jxf_myid_dh);
        strcat(buf, a);
        if ((jxf_logFile = fopen(buf, "w")) == NULL )
        {
            jxf_fprintf(stderr, "can't open >%s< for writing; continuing anyway\n", buf);
        }
    }
}

void jxf_closeLogfile_dh()
{
    if (jxf_logFile != NULL)
    {
        if (fclose(jxf_logFile))
        {
            jxf_fprintf(stderr, "Error closing jxf_logFile\n");
        }
        jxf_logFile = NULL;
    }
}

void jxf_setInfo_dh( char *msg, char *function, char *file, JXF_Int line )
{
    if (jxf_logInfoToFile && jxf_logFile != NULL)
    {
        jxf_fprintf(jxf_logFile, "INFO: %s;\n       function= %s  file=%s  line=%i\n", msg, function, file, line);
        fflush(jxf_logFile);
    }
    if (jxf_logInfoToStderr)
    {
        jxf_fprintf(stderr, "INFO: %s;\n       function= %s  file=%s  line=%i\n", msg, function, file, line);
    }
}

void jxf_dh_StartFunc( char *function, char *file, JXF_Int line, JXF_Int priority )
{
    if (priority == 1)
    {
        jxf_sprintf(jxf_calling_stack[jxf_jxf_calling_stack_count], "[%i]   %s  file= %s  line= %i", jxf_myid_dh, function, file, line);
        ++ jxf_jxf_calling_stack_count;
        if (jxf_jxf_calling_stack_count == JXF_MAX_STACK_SIZE)
        {
            jxf_fprintf(stderr, "_____________ jxf_dh_StartFunc: OVERFLOW _____________________\n");
            if (jxf_logFile != NULL)
            {
                jxf_fprintf(jxf_logFile, "_____________ jxf_dh_StartFunc: OVERFLOW _____________________\n");
            }
            -- jxf_jxf_calling_stack_count;
        }
    }
}

void jxf_dh_EndFunc( char *function, JXF_Int priority )
{
    if (priority == 1)
    {
        -- jxf_jxf_calling_stack_count;
        if (jxf_jxf_calling_stack_count < 0)
        {
            jxf_jxf_calling_stack_count = 0;
            jxf_fprintf(stderr, "_____________ jxf_dh_EndFunc: UNDERFLOW _____________________\n");
            if (jxf_logFile != NULL)
            {
                jxf_fprintf(jxf_logFile, "_____________ jxf_dh_EndFunc: UNDERFLOW _____________________\n");
            }
        }
    }
}

void jxf_setError_dh( char *msg, char *function, char *file, JXF_Int line )
{
    jxf_errFlag_dh = jxf_true;
    if (!strcmp(msg, ""))
    {
        jxf_sprintf(jxf_errMsg_private[jxf_errCount_private],
            "[%i] called from: %s  file= %s  line= %i", jxf_myid_dh, function, file, line);
    }
    else
    {
        jxf_sprintf(jxf_errMsg_private[jxf_errCount_private],
            "[%i] ERROR: %s\n       %s  file= %s  line= %i\n", jxf_myid_dh, msg, function, file, line);
    }
    ++ jxf_errCount_private;
    if (jxf_errCount_private == JXF_MAX_STACK_SIZE) -- jxf_errCount_private;
}

void jxf_printErrorMsg( FILE *fp )
{
    if (!jxf_errFlag_dh)
    {
        jxf_fprintf(fp, "jxf_errFlag_dh is not set; nothing to print!\n");
        fflush(fp);
    }
    else
    {
        JXF_Int i;
        
        jxf_fprintf(fp, "\n============= error stack trace ====================\n");
        for (i = 0; i < jxf_errCount_private; ++ i)
        {
            jxf_fprintf(fp, "%s\n", jxf_errMsg_private[i]);
        }
        jxf_fprintf(fp, "\n");
        fflush(fp);
    }
}

void jxf_printFunctionStack( FILE *fp )
{
    JXF_Int i;
    
    for (i = 0; i < jxf_jxf_calling_stack_count; ++ i)
    {
        jxf_fprintf(fp, "%s\n", jxf_calling_stack[i]);
    }
    jxf_fprintf(fp, "\n");
    fflush(fp);
}

#define JXF_MAX_ERROR_SPACES 200
static char jxf_spaces[JXF_MAX_ERROR_SPACES];
static JXF_Int jxf_nesting = 0;
static jxf_bool jxf_initSpaces = jxf_true;
#define JXF_INDENT_DH 3

void jxf_Error_dhStartFunc( char *function, char *file, JXF_Int line )
{
    if (jxf_initSpaces)
    {
        memset(jxf_spaces, ' ', JXF_MAX_ERROR_SPACES*sizeof(char));
        jxf_initSpaces = jxf_false;
    }
    jxf_spaces[JXF_INDENT_DH*jxf_nesting] = ' ';
    ++ jxf_nesting;
    if (jxf_nesting > JXF_MAX_ERROR_SPACES-1) jxf_nesting = JXF_MAX_ERROR_SPACES - 1;
    jxf_spaces[JXF_INDENT_DH*jxf_nesting] = '\0';
    if (jxf_logFuncsToStderr)
    {
        jxf_fprintf(stderr, "%s(%i) %s  [file= %s  line= %i]\n", jxf_spaces, jxf_nesting, function, file, line);
    }
    if (jxf_logFuncsToFile && jxf_logFile != NULL)
    {
        jxf_fprintf(jxf_logFile, "%s(%i) %s  [file= %s  line= %i]\n", jxf_spaces, jxf_nesting, function, file, line);
        fflush(jxf_logFile);
    }
}

void jxf_Error_dhEndFunc( char *function )
{
    jxf_nesting -= 1;
    if (jxf_nesting < 0) jxf_nesting = 0;
    jxf_spaces[JXF_INDENT_DH*jxf_nesting] = '\0';
}

static jxf_bool jxf_EuclidIsActive = jxf_false;

#undef __FUNC__
#define __FUNC__ "jxf_EuclidIsInitialized"
jxf_bool jxf_EuclidIsInitialized()
{
    return jxf_EuclidIsActive;
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidInitialize"
void jxf_EuclidInitialize( JXF_Int argc, char *argv[], char *help )
{
    if (!jxf_EuclidIsActive)
    {
        jxf_MPI_Comm_size(jxf_comm_dh, &jxf_np_dh);
        jxf_MPI_Comm_rank(jxf_comm_dh, &jxf_myid_dh);
        jxf_openLogfile_dh(argc, argv);
        if (jxf_mem_dh == NULL)
        {
            jxf_Mem_dhCreate(&jxf_mem_dh); JXF_CHECK_V_ERROR;
        }
        if (jxf_tlog_dh == NULL)
        {
            jxf_TimeLog_dhCreate(&jxf_tlog_dh); JXF_CHECK_V_ERROR;
        }
        if (jxf_parser_dh == NULL)
        {
            jxf_Parser_dhCreate(&jxf_parser_dh); JXF_CHECK_V_ERROR;
        }
        jxf_Parser_dhInit(jxf_parser_dh, argc, argv); JXF_CHECK_V_ERROR;
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-sig_dh"))
        {
            jxf_sigRegister_dh(); JXF_CHECK_V_ERROR;
        }
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-help"))
        {
            if (jxf_myid_dh == 0) jxf_printf("%s\n\n", help);
            JXF_EUCLID_EXIT;
        }
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-jxf_logFuncsToFile"))
        {
            jxf_logFuncsToFile = jxf_true;
        }
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-jxf_logFuncsToStderr"))
        {
            jxf_logFuncsToStderr = jxf_true;
        }
        jxf_EuclidIsActive = jxf_true;
    }
}

#undef __FUNC__
#define __FUNC__ "jxf_EuclidFinalize"
void jxf_EuclidFinalize()
{
    if (jxf_ref_counter) return;
    if (jxf_EuclidIsActive)
    {
        if (jxf_parser_dh != NULL)
        {
            jxf_Parser_dhDestroy(jxf_parser_dh); JXF_CHECK_V_ERROR;
        }
        if (jxf_tlog_dh != NULL)
        {
            jxf_TimeLog_dhDestroy(jxf_tlog_dh); JXF_CHECK_V_ERROR;
        }
        if (jxf_logFile != NULL)
        {
            jxf_Mem_dhPrint(jxf_mem_dh, jxf_logFile, jxf_true); JXF_CHECK_V_ERROR;
        }
        if (jxf_mem_dh != NULL)
        {
            jxf_Mem_dhDestroy(jxf_mem_dh); JXF_CHECK_V_ERROR;
        }
        if (jxf_logFile != NULL)
        {
            jxf_closeLogfile_dh(); JXF_CHECK_V_ERROR;
        }
        jxf_EuclidIsActive = jxf_false;
    }
}

#undef __FUNC__
#define __FUNC__ "jxf_printf_dh"
void jxf_printf_dh( char *fmt, ... )
{
    JXF_START_FUNC_DH
    va_list args;
    char *buf = jxf_msgBuf_dh;
    
    va_start(args, fmt);
    vsprintf(buf, fmt, args);
    if (jxf_myid_dh == 0)
    {
        jxf_fprintf(stdout, "%s", buf);
    }
    va_end(args);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_fprintf_dh"
void jxf_fprintf_dh( FILE *fp, char *fmt, ... )
{
    JXF_START_FUNC_DH
    va_list args;
    char *buf = jxf_msgBuf_dh;
    
    va_start(args, fmt);
    vsprintf(buf, fmt, args);
    if (jxf_myid_dh == 0)
    {
        jxf_fprintf(fp, "%s", buf);
    }
    va_end(args);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_echoInvocation_dh"
void jxf_echoInvocation_dh( MPI_Comm comm, char *prefix, JXF_Int argc, char *argv[] )
{
    JXF_START_FUNC_DH
    JXF_Int i, id;
    
    jxf_MPI_Comm_rank(comm, &id);
    if (prefix != NULL)
    {
        jxf_printf_dh("\n%s ", prefix);
    }
    else
    {
        jxf_printf_dh("\n");
    }
    jxf_printf_dh("program invocation: ");
    for (i = 0; i < argc; ++ i)
    {
        jxf_printf_dh("%s ", argv[i]);
    }
    jxf_printf_dh("\n");
    JXF_END_FUNC_DH
}
