//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  numbering_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_Numbering_dhCreate"
void jxf_Numbering_dhCreate( jxf_Numbering_dh *numb )
{
    JXF_START_FUNC_DH
    struct _jxf_numbering_dh *tmp = (struct _jxf_numbering_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_numbering_dh)); JXF_CHECK_V_ERROR;
   *numb = tmp;
    tmp->size = 0;
    tmp->first = 0;
    tmp->m = 0;
    tmp->num_ext = 0;
    tmp->num_extLo = 0;
    tmp->num_extHi = 0;
    tmp->idx_ext = NULL;
    tmp->idx_extLo = NULL;
    tmp->idx_extHi = NULL;
    tmp->idx_ext = NULL;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_Numbering");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Numbering_dhDestroy"
void jxf_Numbering_dhDestroy( jxf_Numbering_dh numb )
{
    JXF_START_FUNC_DH
    if (numb->global_to_local != NULL)
    {
        jxf_Hash_i_dhDestroy(numb->global_to_local); JXF_CHECK_V_ERROR;
    }
    if (numb->idx_ext != NULL)
    {
        JXF_FREE_DH(numb->idx_ext); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(numb); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Numbering_dhSetup"
void jxf_Numbering_dhSetup( jxf_Numbering_dh numb, jxf_Mat_dh mat )
{
    JXF_START_FUNC_DH
    JXF_Int i, len, *cval = mat->cval;
    JXF_Int num_ext, num_extLo, num_extHi;
    JXF_Int m = mat->m, size;
    jxf_Hash_i_dh global_to_local_hash;
    JXF_Int first = mat->beg_row, last  = first + m;
    JXF_Int *idx_ext;
    JXF_Int data;
    
    numb->first = first;
    numb->m = m;
    numb->size = size = m;
    jxf_Hash_i_dhCreate(&(numb->global_to_local), m); JXF_CHECK_V_ERROR;
    global_to_local_hash = numb->global_to_local;
    idx_ext = numb->idx_ext = (JXF_Int*)JXF_MALLOC_DH(size*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    len = mat->rp[m];
    num_ext = num_extLo = num_extHi = 0;
    for (i = 0; i < len; i ++)
    {
        JXF_Int index = cval[i];
        
        if (index < first || index >= last)
        {
            data = jxf_Hash_i_dhLookup(global_to_local_hash, cval[i]); JXF_CHECK_V_ERROR;
            if (data == -1)
            {
                if (m+num_ext >= size)
                {
                    JXF_Int newSize = size * 1.5;
                    JXF_Int *tmp = (JXF_Int *)JXF_MALLOC_DH(newSize*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
                    
                    memcpy(tmp, idx_ext, size*sizeof(size));
                    JXF_FREE_DH(idx_ext); JXF_CHECK_V_ERROR;
                    size = numb->size = newSize;
                    numb->idx_ext = idx_ext = tmp;
                    JXF_SET_INFO("reallocated ext_idx[]");
                }
                jxf_Hash_i_dhInsert(global_to_local_hash, index, num_ext); JXF_CHECK_V_ERROR;
                idx_ext[num_ext] = index;
                num_ext ++;
                if (index < first)
                {
                    num_extLo ++;
                }
                else
                {
                    num_extHi ++;
                }
            }
        }
    }
    numb->num_ext = num_ext;
    numb->num_extLo = num_extLo;
    numb->num_extHi = num_extHi;
    numb->idx_extLo = idx_ext;
    numb->idx_extHi = idx_ext + num_extLo;
    jxf_shellSort_int(num_ext, idx_ext);
    jxf_Hash_i_dhReset(global_to_local_hash); JXF_CHECK_V_ERROR;
    for (i = 0; i < num_ext; i ++)
    {
        jxf_Hash_i_dhInsert(global_to_local_hash, idx_ext[i], i+m); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Numbering_dhGlobalToLocal"
void jxf_Numbering_dhGlobalToLocal( jxf_Numbering_dh numb, JXF_Int len, JXF_Int *global, JXF_Int *local )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int first = numb->first;
    JXF_Int last = first + numb->m;
    JXF_Int data;
    jxf_Hash_i_dh global_to_local = numb->global_to_local;
    
    for (i = 0; i < len; i ++)
    {
        JXF_Int idxGlobal = global[i];
        
        if (idxGlobal >= first && idxGlobal < last)
        {
            local[i] = idxGlobal - first;
        }
        else
        {
            data = jxf_Hash_i_dhLookup(global_to_local, idxGlobal); JXF_CHECK_V_ERROR;
            if (data == -1)
            {
                jxf_sprintf(jxf_msgBuf_dh, "global index %i not found in map\n", idxGlobal);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
            else
            {
                local[i] = data;
            }
        }
    }
    JXF_END_FUNC_DH
}
