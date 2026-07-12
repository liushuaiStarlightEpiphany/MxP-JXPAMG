//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  vec_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"
#include <stdlib.h>

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhCreate"
void jxf_Vec_dhCreate( jxf_Vec_dh *v )
{
    JXF_START_FUNC_DH
    struct _jxf_vec_dh *tmp = (struct _jxf_vec_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_vec_dh)); JXF_CHECK_V_ERROR;
   *v = tmp;
    tmp->n = 0;
    tmp->vals = NULL;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhDestroy"
void jxf_Vec_dhDestroy( jxf_Vec_dh v )
{
    JXF_START_FUNC_DH
    if (v->vals != NULL) JXF_FREE_DH(v->vals); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(v); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhInit"
void jxf_Vec_dhInit( jxf_Vec_dh v, JXF_Int size )
{
    JXF_START_FUNC_DH
    v->n = size;
    v->vals = (JXF_Real *)JXF_MALLOC_DH(size*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhCopy"
void jxf_Vec_dhCopy( jxf_Vec_dh x, jxf_Vec_dh y )
{
    JXF_START_FUNC_DH
    if (x->vals == NULL) JXF_SET_V_ERROR("x->vals is NULL");
    if (y->vals == NULL) JXF_SET_V_ERROR("y->vals is NULL");
    if (x->n != y->n) JXF_SET_V_ERROR("x and y are different lengths");
    memcpy(y->vals, x->vals, x->n*sizeof(JXF_Real));
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhDuplicate"
void jxf_Vec_dhDuplicate( jxf_Vec_dh v, jxf_Vec_dh *out )
{
    JXF_START_FUNC_DH
    jxf_Vec_dh tmp;
    JXF_Int size = v->n;
    
    if (v->vals == NULL) JXF_SET_V_ERROR("v->vals is NULL");
    jxf_Vec_dhCreate(out); JXF_CHECK_V_ERROR;
    tmp = *out;
    tmp->n = size;
    tmp->vals = (JXF_Real *)JXF_MALLOC_DH(size*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhSet"
void jxf_Vec_dhSet( jxf_Vec_dh v, JXF_Real value )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = v->n;
    JXF_Real *vals = v->vals;
    
    if (v->vals == NULL) JXF_SET_V_ERROR("v->vals is NULL");
    for (i = 0; i < m; ++ i) vals[i] = value;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhSetRand"
void jxf_Vec_dhSetRand( jxf_Vec_dh v )
{
    JXF_START_FUNC_DH
    JXF_Int i, m = v->n;
    JXF_Real max = 0.0;
    JXF_Real *vals = v->vals;
    
    if (v->vals == NULL) JXF_SET_V_ERROR("v->vals is NULL");
    for (i = 0; i < m; ++ i) vals[i] = rand();
    for (i = 0; i < m; ++ i) max = JXF_MAX(max, vals[i]);
    for (i = 0; i < m; ++ i) vals[i] = vals[i]/max;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhPrint"
void jxf_Vec_dhPrint( jxf_Vec_dh v, jxf_SubdomainGraph_dh sg, char *filename )
{
    JXF_START_FUNC_DH
    JXF_Real *vals = v->vals;
    JXF_Int pe, i, m = v->n;
    FILE *fp;
    
    if (v->vals == NULL) JXF_SET_V_ERROR("v->vals is NULL");
    if (sg == NULL)
    {
        for (pe = 0; pe < jxf_np_dh; ++ pe)
        {
            jxf_MPI_Barrier(jxf_comm_dh);
            if (pe == jxf_myid_dh)
            {
                if (pe == 0)
                {
                    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
                }
                else
                {
                    fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
                }
                for (i = 0; i < m; ++ i) jxf_fprintf(fp, "%g\n", vals[i]);
                jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    else if (jxf_np_dh == 1)
    {
        JXF_Int i, j;
        
        fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
        for (i = 0; i < sg->blocks; ++ i)
        {
            JXF_Int oldBlock = sg->n2o_sub[i];
            JXF_Int beg_row = sg->beg_rowP[oldBlock];
            JXF_Int end_row = beg_row + sg->row_count[oldBlock];
            
            jxf_printf("seq: block= %i  beg= %i  end= %i\n", oldBlock, beg_row, end_row);
            for (j = beg_row; j < end_row; ++ j)
            {
                jxf_fprintf(fp, "%g\n", vals[j]);
            }
        }
    }
    else
    {
        JXF_Int id = sg->o2n_sub[jxf_myid_dh];
        
        for (pe = 0; pe < jxf_np_dh; ++ pe)
        {
            jxf_MPI_Barrier(jxf_comm_dh);
            if (id == pe)
            {
                if (pe == 0)
                {
                    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
                }
                else
                {
                    fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
                }
                jxf_fprintf(stderr, "par: block= %i\n", id);
                for (i = 0; i < m; ++ i)
                {
                    jxf_fprintf(fp, "%g\n", vals[i]);
                }
                jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhPrintBIN"
void jxf_Vec_dhPrintBIN( jxf_Vec_dh v, jxf_SubdomainGraph_dh sg, char *filename )
{
    JXF_START_FUNC_DH
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    if (sg != NULL)
    {
        JXF_SET_V_ERROR("not implemented for reordered vector; ensure sg=NULL");
    }
    jxf_io_dh_print_ebin_vec_private(v->n, 0, v->vals, NULL, NULL, NULL, filename); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#define JXF_MAX_JUNK 200

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhRead"
void jxf_Vec_dhRead( jxf_Vec_dh *vout, JXF_Int ignore, char *filename )
{
    JXF_START_FUNC_DH
    jxf_Vec_dh tmp;
    FILE *fp;
    JXF_Int items, n, i;
    JXF_Real *v, w;
    char junk[JXF_MAX_JUNK];
    
    jxf_Vec_dhCreate(&tmp); JXF_CHECK_V_ERROR;
   *vout = tmp;
    if (jxf_np_dh > 1)
    {
        JXF_SET_V_ERROR("only implemented for a single jxf_MPI task");
    }
    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
    if (ignore)
    {
        jxf_printf("jxf_Vec_dhRead:: ignoring following header lines:\n");
        jxf_printf("--------------------------------------------------------------\n");
        for (i = 0; i < ignore; ++ i)
        {
            fgets(junk, JXF_MAX_JUNK, fp);
            jxf_printf("%s", junk);
        }
        jxf_printf("--------------------------------------------------------------\n");
    }
    n = 0;
    while (!feof(fp))
    {
        items = jxf_fscanf(fp,"%lg", &w);
        if (items != 1)
        {
            break;
        }
        ++ n;
    }
    jxf_printf("jxf_Vec_dhRead:: n= %i\n", n);
    tmp->n = n;
    v = tmp->vals =  (JXF_Real *)JXF_MALLOC_DH(n*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    rewind(fp);
    rewind(fp);
    for (i = 0; i < ignore; ++ i)
    {
        fgets(junk, JXF_MAX_JUNK, fp);
    }
    for (i = 0; i < n;  ++ i)
    {
        items = jxf_fscanf(fp,"%lg", v+i);
        if (items != 1)
        {
            jxf_sprintf(jxf_msgBuf_dh, "failed to read value %i of %i", i+1, n);
        }
    }
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Vec_dhReadBIN"
extern void jxf_Vec_dhReadBIN( jxf_Vec_dh *vout, char *filename )
{
    JXF_START_FUNC_DH
    jxf_Vec_dh tmp;
    
    jxf_Vec_dhCreate(&tmp); JXF_CHECK_V_ERROR;
   *vout = tmp;
    jxf_io_dh_read_ebin_vec_private(&tmp->n, &tmp->vals, filename); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}
