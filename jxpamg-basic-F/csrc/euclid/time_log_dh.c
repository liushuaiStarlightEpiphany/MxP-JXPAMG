//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  time_log_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#define JXF_MAX_TIME_MARKS 100
#define JXF_MAX_DESC_LENGTH 60

struct _jxf_timeLog_dh
{
    JXF_Int first;
    JXF_Int last;
    JXF_Real time[JXF_MAX_TIME_MARKS];
    char desc[JXF_MAX_TIME_MARKS][JXF_MAX_DESC_LENGTH];
    jxf_Timer_dh timer;
};

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhCreate"
void jxf_TimeLog_dhCreate( jxf_TimeLog_dh *t )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    struct _jxf_timeLog_dh *tmp = (struct _jxf_timeLog_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_timeLog_dh)); JXF_CHECK_V_ERROR;
    
   *t = tmp;
    tmp->first = tmp->last = 0;
    jxf_Timer_dhCreate(&tmp->timer);
    for (i = 0; i < JXF_MAX_TIME_MARKS; ++ i) strcpy(tmp->desc[i], "X");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhDestroy"
void jxf_TimeLog_dhDestroy( jxf_TimeLog_dh t )
{
    JXF_START_FUNC_DH
    jxf_Timer_dhDestroy(t->timer);
    JXF_FREE_DH(t);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhStart"
void jxf_TimeLog_dhStart( jxf_TimeLog_dh t )
{
    JXF_START_FUNC_DH
    jxf_Timer_dhStart(t->timer);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhStop"
void jxf_TimeLog_dhStop( jxf_TimeLog_dh t )
{
    JXF_START_FUNC_DH
    jxf_Timer_dhStop(t->timer);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhMark"
void jxf_TimeLog_dhMark( jxf_TimeLog_dh t, const char *desc )
{
    JXF_START_FUNC_DH
    if (t->last < JXF_MAX_TIME_MARKS - 3)
    {
        jxf_Timer_dhStop(t->timer);
        t->time[t->last] = jxf_Timer_dhReadWall(t->timer);
        jxf_Timer_dhStart(t->timer);
        jxf_sprintf(t->desc[t->last], "%s", desc);
        t->last += 1;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhReset"
void jxf_TimeLog_dhReset( jxf_TimeLog_dh t )
{
    JXF_START_FUNC_DH
    if (t->last < JXF_MAX_TIME_MARKS - 2)
    {
        JXF_Real total = 0.0;
        JXF_Int i, first = t->first, last = t->last;
        
        for (i = first; i < last; ++ i) total += t->time[i];
        t->time[last] = total;
        jxf_sprintf(t->desc[last], "========== totals, and reset ==========\n");
        t->last += 1;
        t->first = t->last;
        jxf_Timer_dhStart(t->timer);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_TimeLog_dhPrint"
void jxf_TimeLog_dhPrint( jxf_TimeLog_dh t, FILE *fp, jxf_bool allPrint )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Real total = 0.0;
    JXF_Real timeMax[JXF_MAX_TIME_MARKS];
    JXF_Real timeMin[JXF_MAX_TIME_MARKS];
    static jxf_bool wasSummed = jxf_false;
    
    if (!wasSummed)
    {
        for (i = t->first; i < t->last; ++ i) total += t->time[i];
        t->time[t->last] = total;
        jxf_sprintf(t->desc[t->last], "========== totals, and reset ==========\n");
        t->last += 1;
        jxf_MPI_Allreduce(t->time, timeMax, t->last, JXF_MPI_REAL, MPI_MAX, jxf_comm_dh);
        jxf_MPI_Allreduce(t->time, timeMin, t->last, JXF_MPI_REAL, MPI_MIN, jxf_comm_dh);
        wasSummed = jxf_true;
    }
    if (fp != NULL)
    {
        if (jxf_myid_dh == 0 || allPrint)
        {
            jxf_fprintf(fp,"\n----------------------------------------- timing report\n");
            jxf_fprintf(fp, "\n   self     max     min\n");
            for (i = 0; i < t->last; ++ i)
            {
                jxf_fprintf(fp, "%7.3f %7.3f %7.3f   #%s\n", t->time[i], timeMax[i], timeMin[i], t->desc[i]);
            }
            fflush(fp);
        }
    }
    JXF_END_FUNC_DH
}
