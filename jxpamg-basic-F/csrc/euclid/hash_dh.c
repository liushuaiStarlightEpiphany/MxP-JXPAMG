//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hash_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

static void jxf_Hash_dhInit_private( jxf_Hash_dh h, JXF_Int s );

#define JXF_CUR_MARK_INIT -1

struct _jxf_hash_node_private
{
    JXF_Int key;
    JXF_Int mark;
    jxf_HashData data;
};

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhCreate"
void jxf_Hash_dhCreate( jxf_Hash_dh *h, JXF_Int size )
{
    JXF_START_FUNC_DH
    struct _jxf_hash_dh *tmp = (struct _jxf_hash_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_hash_dh)); JXF_CHECK_V_ERROR;
   *h = tmp;
    tmp->size = 0;
    tmp->count = 0;
    tmp->curMark = JXF_CUR_MARK_INIT + 1;
    tmp->data = NULL;
    jxf_Hash_dhInit_private(*h,size); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhDestroy"
void jxf_Hash_dhDestroy( jxf_Hash_dh h )
{
    JXF_START_FUNC_DH
    if (h->data != NULL)
    {
        JXF_FREE_DH(h->data); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(h); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhReset"
void jxf_Hash_dhReset( jxf_Hash_dh h )
{
    JXF_START_FUNC_DH
    h->count = 0;
    h->curMark += 1;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhInit_private"
void jxf_Hash_dhInit_private( jxf_Hash_dh h, JXF_Int s )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int size = 16;
    jxf_HashRecord *data;
    
    while (size < s) size *= 2;
    if ((size-s) < (.1*size))
    {
        size *= 2.0;
    }
    h->size = size;
    data = h->data = (jxf_HashRecord *)JXF_MALLOC_DH(size*sizeof(jxf_HashRecord)); JXF_CHECK_V_ERROR;
    for (i = 0; i < size; ++ i)
    {
        data[i].key = -1;
        data[i].mark = -1;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhLookup"
jxf_HashData *jxf_Hash_dhLookup( jxf_Hash_dh h, JXF_Int key )
{
    JXF_START_FUNC_DH
    JXF_Int i, start;
    JXF_Int curMark = h->curMark;
    JXF_Int size = h->size;
    jxf_HashData *retval = NULL;
    jxf_HashRecord *data = h->data;
    
    JXF_HASH_1(key, size, &start)
    for (i = 0; i < size; ++ i)
    {
        JXF_Int tmp, idx;
        
        JXF_HASH_2(key, size, &tmp)
        idx = (start + i * tmp) % size;
        if (data[idx].mark != curMark)
        {
            break;
        }
        else
        {
            if (data[idx].key == key)
            {
                retval = &(data[idx].data);
                break;
            }
        }
    }
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhInsert"
void jxf_Hash_dhInsert( jxf_Hash_dh h, JXF_Int key, jxf_HashData *dataIN )
{
    JXF_START_FUNC_DH
    JXF_Int i, start, size = h->size;
    JXF_Int curMark = h->curMark;
    jxf_HashRecord *data;
    
    data = h->data;
    h->count += 1;
    if (h->count == h->size)
    {
        JXF_SET_V_ERROR("hash table overflow; rehash need implementing!");
    }
    JXF_HASH_1(key, size, &start)
    for (i = 0; i < size; ++ i)
    {
        JXF_Int tmp, idx;
        
        JXF_HASH_2(key, size, &tmp)
        idx = (start + i * tmp) % size;
        if (data[idx].mark < curMark)
        {
            data[idx].key = key;
            data[idx].mark = curMark;
            memcpy(&(data[idx].data), dataIN, sizeof(jxf_HashData));
            break;
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_dhPrint"
void jxf_Hash_dhPrint( jxf_Hash_dh h, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, size = h->size;
    JXF_Int curMark = h->curMark;
    jxf_HashRecord *data = h->data;
    
    jxf_fprintf(fp, "\n--------------------------- hash table \n");
    for (i = 0; i < size; ++ i)
    {
        if (data[i].mark == curMark)
        {
            jxf_fprintf(fp, "key = %2i;  iData = %3i;  fData = %g\n",
                               data[i].key, data[i].data.iData, data[i].data.fData);
        }
    }
    jxf_fprintf(fp, "\n");
    JXF_END_FUNC_DH
}
