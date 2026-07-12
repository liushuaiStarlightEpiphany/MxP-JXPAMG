//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  hash_i_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#define JXF_DEFAULT_TABLE_SIZE 16

static void jxf_rehash_private( jxf_Hash_i_dh h );

#define JXF_HASH_1(k,size,idxOut) \
      { \
       *idxOut = k % size; \
      }

#define JXF_HASH_2(k,size,idxOut) \
      { \
        JXF_Int r = k % (size-13); \
        r = (r % 2) ? r : r+1; \
       *idxOut = r; \
      }

typedef struct _jxf_hash_i_node_private jxf_Hash_i_Record;

struct _jxf_hash_i_node_private
{
    JXF_Int key;
    JXF_Int mark;
    JXF_Int data;
};

struct _jxf_hash_i_dh
{
    JXF_Int size;
    JXF_Int count;
    JXF_Int curMark;
    jxf_Hash_i_Record *data;
};

#undef __FUNC__
#define __FUNC__ "jxf_Hash_i_dhCreate"
void jxf_Hash_i_dhCreate( jxf_Hash_i_dh *h, JXF_Int sizeIN )
{
    JXF_START_FUNC_DH
    JXF_Int i, size;
    jxf_Hash_i_Record *tmp2;
    struct _jxf_hash_i_dh *tmp;
    
    size = JXF_DEFAULT_TABLE_SIZE;
    if (sizeIN == -1)
    {
        sizeIN = size = JXF_DEFAULT_TABLE_SIZE;
    }
    tmp = (struct _jxf_hash_i_dh *)JXF_MALLOC_DH( sizeof(struct _jxf_hash_i_dh)); JXF_CHECK_V_ERROR;
   *h = tmp;
    tmp->size = 0;
    tmp->count = 0;
    tmp->curMark = 0;
    tmp->data = NULL;
    while (size < sizeIN) size *= 2;
    if ((size-sizeIN) < (.1*size))
    {
        size *= 2.0;
    }
    tmp->size = size;
    tmp2 = tmp->data = (jxf_Hash_i_Record *)JXF_MALLOC_DH(size*sizeof(jxf_Hash_i_Record)); JXF_CHECK_V_ERROR;
    for (i = 0; i < size; ++ i)
    {
        tmp2[i].key = -1;
        tmp2[i].mark = -1;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_i_dhDestroy"
void jxf_Hash_i_dhDestroy( jxf_Hash_i_dh h )
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
#define __FUNC__ "jxf_Hash_i_dhReset"
void jxf_Hash_i_dhReset( jxf_Hash_i_dh h )
{
    JXF_START_FUNC_DH
    h->count = 0;
    h->curMark += 1;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_i_dhLookup"
JXF_Int jxf_Hash_i_dhLookup( jxf_Hash_i_dh h, JXF_Int key )
{
    JXF_START_FUNC_DH
    JXF_Int idx, inc, i, start;
    JXF_Int curMark = h->curMark;
    JXF_Int size = h->size;
    JXF_Int retval = -1;
    jxf_Hash_i_Record *data = h->data;
    
    JXF_HASH_1(key, size, &start)
    JXF_HASH_2(key, size, &inc)
    for (i = 0; i < size; ++ i)
    {
        idx = (start + i * inc) % size;
        if (data[idx].mark != curMark)
        {
            break;
        }
        else
        {
            if (data[idx].key == key)
            {
                retval = data[idx].data;
                break;
            }
        }
    }
    JXF_END_FUNC_VAL(retval)
}

#undef __FUNC__
#define __FUNC__ "jxf_Hash_i_dhInsert"
void jxf_Hash_i_dhInsert( jxf_Hash_i_dh h, JXF_Int key, JXF_Int dataIN )
{
    JXF_START_FUNC_DH
    JXF_Int i, idx, inc, start, size;
    JXF_Int curMark = h->curMark;
    jxf_Hash_i_Record *data;
    jxf_bool success = jxf_false;
    
    if (dataIN < 0)
    {
        jxf_sprintf(jxf_msgBuf_dh, "data = %i must be >= 0", dataIN);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    if (h->count >= 0.9 * h->size)
    {
        jxf_rehash_private(h); JXF_CHECK_V_ERROR;
    }
    size = h->size;
    data = h->data;
    h->count += 1;
    JXF_HASH_1(key, size, &start)
    JXF_HASH_2(key, size, &inc)
    for (i = 0; i < size; ++ i)
    {
        idx = (start + i * inc) % size;
        if (data[idx].mark == curMark && data[idx].key == key)
        {
            jxf_sprintf(jxf_msgBuf_dh, "key,data= <%i, %i> already inserted", key, dataIN);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
        if (data[idx].mark < curMark)
        {
            data[idx].key = key;
            data[idx].mark = curMark;
            data[idx].data = dataIN;
            success = jxf_true;
            break;
        }
    }
    if (!success)
    {
        jxf_sprintf(jxf_msgBuf_dh, "Failed to insert key= %i, data= %i", key, dataIN);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_rehash_private"
void jxf_rehash_private( jxf_Hash_i_dh h )
{
    JXF_START_FUNC_DH
    JXF_Int i, old_size = h->size, new_size = old_size * 2, oldCurMark = h->curMark;
    jxf_Hash_i_Record *oldData = h->data, *newData;
    
    jxf_sprintf(jxf_msgBuf_dh, "rehashing; old_size= %i, new_size= %i", old_size, new_size);
    JXF_SET_INFO(jxf_msgBuf_dh);
    newData = (jxf_Hash_i_Record *)JXF_MALLOC_DH(new_size*sizeof(jxf_Hash_i_Record)); JXF_CHECK_V_ERROR;
    for (i = 0; i < new_size; ++ i)
    {
        newData[i].key = -1;
        newData[i].mark = -1;
    }
    h->size = new_size;
    h->data = newData;
    h->count = 0;
    h->curMark = 0;
    for (i = h->count; i < new_size; ++ i)
    {
        newData[i].key = -1;
        newData[i].mark = -1;
    }
    for (i = 0; i < old_size; ++ i)
    {
        if (oldData[i].mark == oldCurMark)
        {
            jxf_Hash_i_dhInsert(h, oldData[i].key, oldData[i].data); JXF_CHECK_V_ERROR;
        }
    }
    JXF_FREE_DH(oldData); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}
