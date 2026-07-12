//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  sig_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"
#include <signal.h>

extern void jxf_sigRegister_dh();
extern void jxf_sigHandler_dh( int sig );

#ifdef WIN32
int jxf_euclid_signals_len = 2;
int jxf_euclid_signals[] = {SIGSEGV, SIGFPE};
#else
int jxf_euclid_signals_len = 3;
int jxf_euclid_signals[] = {SIGSEGV, SIGFPE, SIGBUS};
#endif

static char *JXF_SIGNAME[] = {
    "Unknown signal",
    "HUP (Hangup detected on controlling terminal or death of controlling process)",
    "INT: Interrupt from keyboard",
    "QUIT: Quit from keyboard",
    "ILL: Illegal Instruction",
    "TRAP",
    "ABRT: Abort signal",
    "EMT",
    "FPE (Floating Point Exception)",
    "KILL: Kill signal",
    "BUS (Bus Error, possibly illegal memory access)",
    "SEGV (Segmentation Violation (memory access out of range?))",
    "SYS",
    "PIPE: Broken pipe: write to pipe with no readers",
    "ALRM: Timer signal",
    "TERM: Termination signal",
    "URG",
    "STOP",
    "TSTP",
    "CONT",
    "CHLD"
};

#undef __FUNC__
#define __FUNC__ "jxf_sigHandler_dh"
void jxf_sigHandler_dh( int sig )
{
    jxf_fprintf(stderr, "\n[%i] Euclid Signal Handler got: %s\n", jxf_myid_dh, JXF_SIGNAME[sig]);
    jxf_fprintf(stderr, "[%i] ========================================================\n", jxf_myid_dh);
    jxf_fprintf(stderr, "[%i] function calling sequence that led to the exception:\n", jxf_myid_dh);
    jxf_fprintf(stderr, "[%i] ========================================================\n", jxf_myid_dh);
    jxf_printFunctionStack(stderr);
    jxf_fprintf(stderr, "\n\n");
    if (jxf_logFile != NULL)
    {
        jxf_fprintf(jxf_logFile, "\n[%i] Euclid Signal Handler got: %s\n", jxf_myid_dh, JXF_SIGNAME[sig]);
        jxf_fprintf(jxf_logFile, "[%i] ========================================================\n", jxf_myid_dh);
        jxf_fprintf(jxf_logFile, "[%i] function calling sequence that led to the exception:\n", jxf_myid_dh);
        jxf_fprintf(jxf_logFile, "[%i] ========================================================\n", jxf_myid_dh);
        jxf_printFunctionStack(jxf_logFile);
        jxf_fprintf(jxf_logFile, "\n\n");
    }
    JXF_EUCLID_EXIT;
}

#undef __FUNC__
#define __FUNC__ "jxf_sigRegister_dh"
void jxf_sigRegister_dh()
{
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-sig_dh"))
    {
        int i;
        
        for (i = 0; i < jxf_euclid_signals_len; ++ i)
        {
            signal(jxf_euclid_signals[i], jxf_sigHandler_dh);
        }
    }
}
