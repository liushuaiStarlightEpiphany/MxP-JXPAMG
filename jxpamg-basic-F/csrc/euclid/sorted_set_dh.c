//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  sorted_set_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_SortedSet_dhCreate"
void jxf_SortedSet_dhCreate( jxf_SortedSet_dh *ss, JXF_Int size )
{
    JXF_START_FUNC_DH
    struct _jxf_sortedset_dh *tmp = (struct _jxf_sortedset_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_sortedset_dh)); JXF_CHECK_V_ERROR;
    
   *ss= tmp;
    tmp->n = size;
    tmp->list = (JXF_Int *)JXF_MALLOC_DH(size*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    tmp->count = 0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedSet_dhDestroy"
void jxf_SortedSet_dhDestroy( jxf_SortedSet_dh ss )
{
    JXF_START_FUNC_DH
    if (ss->list != NULL)
    {
        JXF_FREE_DH(ss->list); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(ss); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedSet_dhInsert"
void jxf_SortedSet_dhInsert( jxf_SortedSet_dh ss, JXF_Int idx )
{
    JXF_START_FUNC_DH
    jxf_bool isInserted = jxf_false;
    JXF_Int ct = ss->count;
    JXF_Int *list = ss->list;
    JXF_Int i, n = ss->n;
    
    for (i = 0; i < ct; ++ i)
    {
        if (list[i] == idx)
        {
            isInserted = jxf_true;
            break;
        }
    }
    if (!isInserted)
    {
        if (ct == n)
        {
            JXF_Int *tmp = (JXF_Int *)JXF_MALLOC_DH(n*2*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
            
            memcpy(tmp, list, n*sizeof(JXF_Int));
            JXF_FREE_DH(list); JXF_CHECK_V_ERROR;
            list = ss->list = tmp;
            ss->n *= 2;
        }
        list[ct] = idx;
        ss->count += 1;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedSet_dhGetList"
void jxf_SortedSet_dhGetList( jxf_SortedSet_dh ss, JXF_Int **list, JXF_Int *count )
{
    JXF_START_FUNC_DH
    jxf_shellSort_int(ss->count, ss->list);
   *list = ss->list;
   *count = ss->count;
    JXF_END_FUNC_DH
}
