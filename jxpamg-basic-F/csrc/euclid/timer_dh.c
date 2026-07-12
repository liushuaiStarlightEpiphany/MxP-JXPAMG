//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  timer_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhCreate"
void jxf_Timer_dhCreate( jxf_Timer_dh *t )
{
    JXF_START_FUNC_DH
    struct _jxf_timer_dh *tmp = (struct _jxf_timer_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_timer_dh)); JXF_CHECK_V_ERROR;
    
   *t = tmp;
    tmp->isRunning = jxf_false;
    tmp->begin_wall = 0.0;
    tmp->end_wall = 0.0;
#ifdef EUCLID_TIMING
    tmp->sc_clk_tck = sysconf(_SC_CLK_TCK);
#else
    tmp->sc_clk_tck = CLOCKS_PER_SEC;
#endif
#if defined(EUCLID_TIMING)
    jxf_sprintf(jxf_msgBuf_dh, "using EUCLID_TIMING; _SC_CLK_TCK = %i", (JXF_Int)tmp->sc_clk_tck);
    JXF_SET_INFO(jxf_msgBuf_dh);
#elif defined(jxf_MPI_TIMING)
    JXF_SET_INFO("using jxf_MPI timing")
#else
    JXF_SET_INFO("using JUNK timing")
#endif
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhDestroy"
void jxf_Timer_dhDestroy( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_FREE_DH(t);
    JXF_END_FUNC_DH
}

#ifdef EUCLID_TIMING

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStart"
void jxf_Timer_dhStart( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    t->begin_wall = times(&(t->begin_cpu));
    t->isRunning = jxf_true;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStop"
void jxf_Timer_dhStop( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    t->end_wall = times(&(t->end_cpu));
    t->isRunning = jxf_false;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadWall"
JXF_Real jxf_Timer_dhReadWall( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_Real retval = 0.0;
    long JXF_Int sc_clk_tck = t->sc_clk_tck;
    
    if (t->isRunning) t->end_wall = times(&(t->end_cpu));
    retval = (JXF_Real)(t->end_wall - t->begin_wall) / (JXF_Real)sc_clk_tck;
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadCPU"
JXF_Real jxf_Timer_dhReadCPU( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_Real retval;
    long JXF_Int sc_clk_tck = t->sc_clk_tck;
    
    if (t->isRunning) t->end_wall = times(&(t->end_cpu));
    retval = (JXF_Real)(t->end_cpu.tms_utime - t->begin_cpu.tms_utime
                    + t->end_cpu.tms_stime -  t->begin_cpu.tms_stime
                    + t->end_cpu.tms_cutime - t->begin_cpu.tms_cutime
                    + t->end_cpu.tms_cstime -  t->begin_cpu.tms_cstime) / (JXF_Real)sc_clk_tck;
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadUsage"
JXF_Real jxf_Timer_dhReadUsage( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_Real cpu = jxf_Timer_dhReadCPU(t);
    JXF_Real wall = jxf_Timer_dhReadWall(t);
    JXF_Real retval = 100.0 * cpu / wall;
    JXF_END_FUNC_VAL(retval);
}

#elif defined(jxf_MPI_TIMING)

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStart"
void jxf_Timer_dhStart( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    t->begin_wall = jxf_MPI_Wtime();
    t->isRunning = jxf_true;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStop"
void jxf_Timer_dhStop( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    t->end_wall = jxf_MPI_Wtime();
    t->isRunning = jxf_false;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadWall"
JXF_Real jxf_Timer_dhReadWall( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_Real retval;
    
    if (t->isRunning) t->end_wall = jxf_MPI_Wtime();
    retval = t->end_wall - t->begin_wall;
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadCPU"
JXF_Real jxf_Timer_dhReadCPU( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(-1.0)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadUsage"
JXF_Real jxf_Timer_dhReadUsage( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(-1.0);
}

#else

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStart"
void jxf_Timer_dhStart( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhStop"
void jxf_Timer_dhStop( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadWall"
JXF_Real jxf_Timer_dhReadWall( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(-1.0)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadCPU"
JXF_Real jxf_Timer_dhReadCPU( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(-1.0)
}

#undef __FUNC__
#define __FUNC__ "jxf_Timer_dhReadUsage"
JXF_Real jxf_Timer_dhReadUsage( jxf_Timer_dh t )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(-1.0);
}

#endif
