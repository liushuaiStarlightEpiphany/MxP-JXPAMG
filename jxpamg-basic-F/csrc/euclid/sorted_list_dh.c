//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  sorted_list_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

struct _jxf_sortedList_dh
{
    JXF_Int m;
    JXF_Int row;
    JXF_Int beg_row;
    JXF_Int beg_rowP;
    JXF_Int count;
    JXF_Int countMax;
    JXF_Int *o2n_local;
    jxf_Hash_i_dh o2n_external;
    
    jxf_SRecord *list;
    JXF_Int alloc;
    JXF_Int getLower;
    JXF_Int get;
    
    jxf_bool debug;
};

static void jxf_lengthen_list_private( jxf_SortedList_dh sList );

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhCreate"
void jxf_SortedList_dhCreate( jxf_SortedList_dh *sList )
{
    JXF_START_FUNC_DH
    struct _jxf_sortedList_dh *tmp = (struct _jxf_sortedList_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_sortedList_dh)); JXF_CHECK_V_ERROR;
    
   *sList = tmp;
    tmp->m = 0;
    tmp->row = -1;
    tmp->beg_row = 0;
    tmp->count = 1;
    tmp->countMax = 1;
    tmp->o2n_external = NULL;
    tmp->o2n_local = NULL;
    tmp->get = 0;
    tmp->getLower = 0;
    tmp->alloc = 0;
    tmp->list = NULL;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_SortedList");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhDestroy"
void jxf_SortedList_dhDestroy( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    if (sList->list != NULL)
    {
        JXF_FREE_DH(sList->list); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(sList); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhInit"
void jxf_SortedList_dhInit( jxf_SortedList_dh sList, jxf_SubdomainGraph_dh sg )
{
    JXF_START_FUNC_DH
    sList->o2n_local = sg->o2n_col;
    sList->m = sg->m;
    sList->beg_row = sg->beg_row[jxf_myid_dh];
    sList->beg_rowP = sg->beg_rowP[jxf_myid_dh];
    sList->count = 1;
    sList->countMax = 1;
    sList->o2n_external = sg->o2n_ext;
    sList->alloc = sList->m + 5;
    sList->list = (jxf_SRecord *)JXF_MALLOC_DH(sList->alloc*sizeof(jxf_SRecord));
    sList->list[0].col = INT_MAX;
    sList->list[0].next = 0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhReset"
void jxf_SortedList_dhReset( jxf_SortedList_dh sList, JXF_Int row )
{
    JXF_START_FUNC_DH
    sList->row = row;
    sList->count = 1;
    sList->countMax = 1;
    sList->get = 0;
    sList->getLower = 0;
    sList->list[0].next = 0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhReadCount"
JXF_Int jxf_SortedList_dhReadCount( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    JXF_END_FUNC_VAL(sList->count-1)
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhResetGetSmallest"
void jxf_SortedList_dhResetGetSmallest( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    sList->getLower = 0;
    sList->get = 0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhGetSmallest"
jxf_SRecord *jxf_SortedList_dhGetSmallest( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    jxf_SRecord *node = NULL;
    jxf_SRecord *list = sList->list;
    JXF_Int get = sList->get;
    
    get = list[get].next;
    if (list[get].col < INT_MAX)
    {
        node = &(list[get]);
        sList->get = get;
    }
    JXF_END_FUNC_VAL(node)
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhGetSmallestLowerTri"
jxf_SRecord *jxf_SortedList_dhGetSmallestLowerTri( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    jxf_SRecord *node = NULL;
    jxf_SRecord *list = sList->list;
    JXF_Int getLower = sList->getLower;
    JXF_Int globalRow = sList->row + sList->beg_rowP;
    
    getLower = list[getLower].next;
    if (list[getLower].col < globalRow)
    {
        node = &(list[getLower]);
        sList->getLower = getLower;
    }
    JXF_END_FUNC_VAL(node)
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhPermuteAndInsert"
jxf_bool jxf_SortedList_dhPermuteAndInsert( jxf_SortedList_dh sList, jxf_SRecord *sr, JXF_Real thresh )
{
    JXF_START_FUNC_DH
    jxf_bool wasInserted = jxf_false;
    JXF_Int col = sr->col;
    JXF_Real testVal = fabs(sr->val);
    JXF_Int beg_row = sList->beg_row, end_row = beg_row + sList->m;
    JXF_Int beg_rowP = sList->beg_rowP;
    
    if (col >= beg_row && col < end_row)
    {
        col -= beg_row;
        col = sList->o2n_local[col];
        if (testVal > thresh || col == sList->row)
        {
            col += beg_rowP;
        }
        else
        {
            col = -1;
        }
    }
    else
    {
        if (testVal < thresh) goto END_OF_FUNCTION;
        if (sList->o2n_external == NULL)
        {
            col = -1;
        }
        else
        {
            JXF_Int tmp = jxf_Hash_i_dhLookup(sList->o2n_external, col); JXF_CHECK_ERROR(-1);
            
            if (tmp == -1)
            {
                col = -1;
            }
            else
            {
                col = tmp;
            }
        }
    }
    if (col != -1)
    {
        sr->col = col;
        jxf_SortedList_dhInsert(sList, sr); JXF_CHECK_ERROR(-1);
        wasInserted = jxf_true;
    }
    
END_OF_FUNCTION: ;
    
    JXF_END_FUNC_VAL(wasInserted)
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhInsertOrUpdate"
void jxf_SortedList_dhInsertOrUpdate( jxf_SortedList_dh sList, jxf_SRecord *sr )
{
    JXF_START_FUNC_DH
    jxf_SRecord *node = jxf_SortedList_dhFind(sList, sr); JXF_CHECK_V_ERROR;
    
    if (node == NULL)
    {
        jxf_SortedList_dhInsert(sList, sr); JXF_CHECK_V_ERROR;
    }
    else
    {
        node->level = JXF_MIN(sr->level, node->level);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhInsert"
void jxf_SortedList_dhInsert( jxf_SortedList_dh sList, jxf_SRecord *sr )
{
    JXF_START_FUNC_DH
    JXF_Int prev, next;
    JXF_Int ct, col = sr->col;
    jxf_SRecord *list = sList->list;
    
    if (sList->countMax == sList->alloc)
    {
        jxf_lengthen_list_private(sList); JXF_CHECK_V_ERROR;
        list = sList->list;
    }
    ct = sList->countMax;
    sList->countMax += 1;
    sList->count += 1;
    list[ct].col = col;
    list[ct].level = sr->level;
    list[ct].val = sr->val;
    prev = 0;
    next = list[0].next;
    while (col > list[next].col)
    {
        prev = next;
        next = list[next].next;
    }
    list[prev].next = ct;
    list[ct].next = next;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhFind"
jxf_SRecord *jxf_SortedList_dhFind( jxf_SortedList_dh sList, jxf_SRecord *sr )
{
    JXF_START_FUNC_DH
    JXF_Int i, count = sList->countMax;
    JXF_Int c = sr->col;
    jxf_SRecord *s = sList->list;
    jxf_SRecord *node = NULL;
    
    for (i = 1; i < count; ++ i)
    {
        if (s[i].col == c)
        {
            node = &(s[i]);
            break;
        }
    }
    JXF_END_FUNC_VAL(node)
}

#undef __FUNC__
#define __FUNC__ "jxf_lengthen_list_private"
void jxf_lengthen_list_private( jxf_SortedList_dh sList )
{
    JXF_START_FUNC_DH
    jxf_SRecord *tmp = sList->list;
    JXF_Int size = sList->alloc = 2 * sList->alloc;
    
    JXF_SET_INFO("lengthening list");
    sList->list = (jxf_SRecord *)JXF_MALLOC_DH(size*sizeof(jxf_SRecord));
    memcpy(sList->list, tmp, sList->countMax*sizeof(jxf_SRecord));
    JXF_SET_INFO("doubling size of sList->list");
    JXF_FREE_DH(tmp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

static jxf_bool jxf_check_constraint_private( jxf_SubdomainGraph_dh sg, JXF_Int thisSubdomain, JXF_Int col );
void jxf_delete_private( jxf_SortedList_dh sList, JXF_Int col );

#undef __FUNC__
#define __FUNC__ "jxf_SortedList_dhEnforceConstraint"
void jxf_SortedList_dhEnforceConstraint( jxf_SortedList_dh sList, jxf_SubdomainGraph_dh sg )
{
    JXF_START_FUNC_DH
    JXF_Int thisSubdomain = jxf_myid_dh;
    JXF_Int col, count;
    JXF_Int beg_rowP = sList->beg_rowP;
    JXF_Int end_rowP = beg_rowP + sList->m;
    jxf_bool debug = jxf_false;
    
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_SortedList")) debug = jxf_true;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "SLIST ======= enforcing constraint for row= %i\n", 1+sList->row);
        jxf_fprintf(jxf_logFile, "\nSLIST ---- before checking: ");
        count = jxf_SortedList_dhReadCount(sList); JXF_CHECK_V_ERROR;
        while (count --)
        {
            jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(sList); JXF_CHECK_V_ERROR;
            jxf_fprintf(jxf_logFile, "%i ", sr->col+1);
        }
        jxf_fprintf(jxf_logFile, "\n");
        sList->get = 0;
    }
    count = jxf_SortedList_dhReadCount(sList); JXF_CHECK_V_ERROR;
    while (count --)
    {
        jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(sList); JXF_CHECK_V_ERROR;
        col = sr->col;
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "SLIST  next col= %i\n", col+1);
        }
        if (col < beg_rowP || col >= end_rowP)
        {
            if (debug)
            {
                jxf_fprintf(jxf_logFile, "SLIST     external col: %i ; ", 1+col);
            }
            if (jxf_check_constraint_private(sg, thisSubdomain, col))
            {
                jxf_delete_private(sList, col); JXF_CHECK_V_ERROR;
                sList->count -= 1;
                if (debug)
                {
                    jxf_fprintf(jxf_logFile, " deleted\n");
                }
            }
            else
            {
                if (debug)
                {
                    jxf_fprintf(jxf_logFile, " kept\n");
                }
            }
        }
    }
    sList->get = 0;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "SLIST---- after checking: ");
        count = jxf_SortedList_dhReadCount(sList); JXF_CHECK_V_ERROR;
        while (count --)
        {
            jxf_SRecord *sr = jxf_SortedList_dhGetSmallest(sList); JXF_CHECK_V_ERROR;
            jxf_fprintf(jxf_logFile, "%i ", sr->col+1);
        }
        jxf_fprintf(jxf_logFile, "\n");
        fflush(jxf_logFile);
        sList->get = 0;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_check_constraint_private"
jxf_bool jxf_check_constraint_private( jxf_SubdomainGraph_dh sg, JXF_Int p1, JXF_Int j )
{
    JXF_START_FUNC_DH
    jxf_bool retval = jxf_false;
    JXF_Int i, p2;
    JXF_Int *nabors, count;
    
    p2 = jxf_SubdomainGraph_dhFindOwner(sg, j, jxf_true);
    nabors = sg->adj + sg->ptrs[p1];
    count = sg->ptrs[p1+1]  - sg->ptrs[p1];
    for (i = 0; i < count; ++ i)
    {
        if (nabors[i] == p2)
        {
            retval = jxf_true;
            break;
        }
    }
    JXF_END_FUNC_VAL(! retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_delete_private"
void jxf_delete_private( jxf_SortedList_dh sList, JXF_Int col )
{
    JXF_START_FUNC_DH
    JXF_Int curNode = 0;
    jxf_SRecord *list = sList->list;
    JXF_Int next;
    
    while (list[list[curNode].next].col != col)
    {
        curNode = list[curNode].next;
    }
    next = list[curNode].next;
    list[next].col = -1;
    next = list[next].next;
    list[curNode].next = next;
    JXF_END_FUNC_DH
}
