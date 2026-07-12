//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_mem_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

typedef struct
{
    JXF_Real size;
    JXF_Real cookie;
} jxf_memRecord_dh;

struct _jxf_mem_dh
{
    JXF_Real maxMem;
    JXF_Real curMem;
    JXF_Real totalMem;
    JXF_Real mallocCount;
    JXF_Real freeCount;
};

#undef __FUNC__
#define __FUNC__ "jxf_Mem_dhCreate"
void jxf_Mem_dhCreate( jxf_Mem_dh *m )
{
    JXF_START_FUNC_DH
    struct _jxf_mem_dh *tmp = (struct _jxf_mem_dh *)JXF_PRIVATE_MALLOC(sizeof(struct _jxf_mem_dh)); JXF_CHECK_V_ERROR;
   *m = tmp;
    tmp->maxMem = 0.0;
    tmp->curMem = 0.0;
    tmp->totalMem = 0.0;
    tmp->mallocCount = 0.0;
    tmp->freeCount = 0.0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mem_dhDestroy"
void jxf_Mem_dhDestroy( jxf_Mem_dh m )
{
    JXF_START_FUNC_DH
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-eu_mem"))
    {
        jxf_Mem_dhPrint(m, stdout, jxf_false); JXF_CHECK_V_ERROR;
    }
    JXF_PRIVATE_FREE(m);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Mem_dhMalloc"
void *jxf_Mem_dhMalloc( jxf_Mem_dh m, size_t size )
{
    JXF_START_FUNC_DH_2
    void *retval;
    jxf_memRecord_dh *tmp;
    size_t s = size + 2 * sizeof(jxf_memRecord_dh);
    void *address;
    
    address = JXF_PRIVATE_MALLOC(s);
    if (address == NULL)
    {
        jxf_sprintf(jxf_msgBuf_dh, "JXF_PRIVATE_MALLOC failed; totalMem = %g; requested additional = %i", m->totalMem, (JXF_Int)s);
        JXF_SET_ERROR(NULL, jxf_msgBuf_dh);
    }
    retval = (char *)address + sizeof(jxf_memRecord_dh);
    tmp = (jxf_memRecord_dh *)address;
    tmp->size = (JXF_Real)s;
    m->mallocCount += 1;
    m->totalMem += (JXF_Real)s;
    m->curMem += (JXF_Real)s;
    m->maxMem = JXF_MAX(m->maxMem, m->curMem);
    JXF_END_FUNC_VAL_2(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Mem_dhFree"
void jxf_Mem_dhFree( jxf_Mem_dh m, void *ptr )
{
    JXF_START_FUNC_DH_2
    JXF_Real size;
    char *tmp = (char *)ptr;
    jxf_memRecord_dh *rec;
    tmp -= sizeof(jxf_memRecord_dh);
    rec = (jxf_memRecord_dh *)tmp;
    size = rec->size;
    jxf_mem_dh->curMem -= size;
    jxf_mem_dh->freeCount += 1;
    JXF_PRIVATE_FREE(tmp);
    JXF_END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "jxf_Mem_dhPrint"
void jxf_Mem_dhPrint( jxf_Mem_dh m, FILE *fp, jxf_bool allPrint )
{
    JXF_START_FUNC_DH_2
    if (fp == NULL) JXF_SET_V_ERROR("fp == NULL");
    if (jxf_myid_dh == 0 || allPrint)
    {
        JXF_Real tmp;
        
        jxf_fprintf(fp, "---------------------- Euclid memory report (start)\n");
        jxf_fprintf(fp, "malloc calls = %g\n", m->mallocCount);
        jxf_fprintf(fp, "free   calls = %g\n", m->freeCount);
        jxf_fprintf(fp, "curMem          = %g Mbytes (should be zero)\n", m->curMem/1000000);
        tmp = m->totalMem / 1000000;
        jxf_fprintf(fp, "total allocated = %g Mbytes\n", tmp);
        jxf_fprintf(fp, "max malloc      = %g Mbytes (max allocated at any point in time)\n", m->maxMem/1000000);
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "---------------------- Euclid memory report (end)\n");
    }
    JXF_END_FUNC_DH_2
}
