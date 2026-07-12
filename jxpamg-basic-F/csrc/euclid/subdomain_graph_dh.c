//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  subdomain_graph_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#ifndef WIN32
#include <unistd.h>
#endif

static void jxf_init_seq_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A );
static void jxf_init_mpi_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A );
static void jxf_allocate_storage_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, JXF_Int m, jxf_bool bj );
static void jxf_form_subdomaingraph_mpi_private( jxf_SubdomainGraph_dh s );
static void jxf_form_subdomaingraph_seq_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A );
static void jxf_find_all_neighbors_sym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A );
static void jxf_find_all_neighbors_unsym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A );
static void jxf_find_bdry_nodes_sym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A,
                JXF_Int *interiorNodes, JXF_Int *bdryNodes, JXF_Int *interiorCount, JXF_Int *bdryCount );
static void jxf_find_bdry_nodes_unsym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A,
                   JXF_Int *interiorNodes, JXF_Int *bdryNodes, JXF_Int *interiorCount, JXF_Int *bdryCount );
static void jxf_find_bdry_nodes_seq_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A );
static void jxf_find_ordered_neighbors_private( jxf_SubdomainGraph_dh s );
static void jxf_color_subdomain_graph_private( jxf_SubdomainGraph_dh s );
static void jxf_adjust_matrix_perms_private( jxf_SubdomainGraph_dh s, JXF_Int m );

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhCreate"
void jxf_SubdomainGraph_dhCreate( jxf_SubdomainGraph_dh *s )
{
    JXF_START_FUNC_DH
    struct _jxf_subdomain_dh *tmp = (struct _jxf_subdomain_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_subdomain_dh)); JXF_CHECK_V_ERROR;
   *s = tmp;
    tmp->blocks = 1;
    tmp->ptrs = tmp->adj = NULL;
    tmp->colors = 1;
    tmp->colorVec = NULL;
    tmp->o2n_sub = tmp->n2o_sub = NULL;
    tmp->beg_row = tmp->beg_rowP = NULL;
    tmp->bdry_count = tmp->row_count = NULL;
    tmp->loNabors = tmp->hiNabors = tmp->allNabors = NULL;
    tmp->loCount = tmp->hiCount = tmp->allCount = 0;
    tmp->m = 0;
    tmp->n2o_row = tmp->o2n_col = NULL;
    tmp->o2n_ext = tmp->n2o_ext = NULL;
    tmp->doNotColor = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-doNotColor");
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_SubGraph");
    JXF_Int i;
    for (i = 0; i < JXF_JXF_TIMING_BINS_SG; ++ i) tmp->timing[i] = 0.0;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhDestroy"
void jxf_SubdomainGraph_dhDestroy( jxf_SubdomainGraph_dh s )
{
    JXF_START_FUNC_DH
    if (s->ptrs != NULL)
    {
        JXF_FREE_DH(s->ptrs); JXF_CHECK_V_ERROR;
    }
    if (s->adj != NULL)
    {
        JXF_FREE_DH(s->adj); JXF_CHECK_V_ERROR;
    }
    if (s->colorVec != NULL)
    {
        JXF_FREE_DH(s->colorVec); JXF_CHECK_V_ERROR;
    }
    if (s->o2n_sub != NULL)
    {
        JXF_FREE_DH(s->o2n_sub); JXF_CHECK_V_ERROR;
    }
    if (s->n2o_sub != NULL)
    {
        JXF_FREE_DH(s->n2o_sub); JXF_CHECK_V_ERROR;
    }
    if (s->beg_row != NULL)
    {
        JXF_FREE_DH(s->beg_row); JXF_CHECK_V_ERROR;
    }
    if (s->beg_rowP != NULL)
    {
        JXF_FREE_DH(s->beg_rowP); JXF_CHECK_V_ERROR;
    }
    if (s->row_count != NULL)
    {
        JXF_FREE_DH(s->row_count); JXF_CHECK_V_ERROR;
    }
    if (s->bdry_count != NULL)
    {
        JXF_FREE_DH(s->bdry_count); JXF_CHECK_V_ERROR;
    }
    if (s->loNabors != NULL)
    {
        JXF_FREE_DH(s->loNabors); JXF_CHECK_V_ERROR;
    }
    if (s->hiNabors != NULL)
    {
        JXF_FREE_DH(s->hiNabors); JXF_CHECK_V_ERROR;
    }
    if (s->allNabors != NULL)
    {
        JXF_FREE_DH(s->allNabors); JXF_CHECK_V_ERROR;
    }
    if (s->n2o_row != NULL)
    {
        JXF_FREE_DH(s->n2o_row); JXF_CHECK_V_ERROR;
    }
    if (s->o2n_col != NULL)
    {
        JXF_FREE_DH(s->o2n_col); JXF_CHECK_V_ERROR;
    }
    if (s->o2n_ext != NULL)
    {
        jxf_Hash_i_dhDestroy(s->o2n_ext); JXF_CHECK_V_ERROR;
    }
    if (s->n2o_ext != NULL)
    {
        jxf_Hash_i_dhDestroy(s->n2o_ext); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(s); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhInit"
void jxf_SubdomainGraph_dhInit( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A )
{
    JXF_START_FUNC_DH
    JXF_Real t1 = jxf_MPI_Wtime();
    
    if (blocks < 1) blocks = 1;
    if (jxf_np_dh == 1 || blocks > 1)
    {
        s->blocks = blocks;
        jxf_init_seq_private(s, blocks, bj, A); JXF_CHECK_V_ERROR;
    }
    else
    {
        s->blocks = jxf_np_dh;
        jxf_init_mpi_private(s, jxf_np_dh, bj, A); JXF_CHECK_V_ERROR;
    }
    s->timing[JXF_TOTAL_SGT] += (jxf_MPI_Wtime() - t1);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhFindOwner"
JXF_Int jxf_SubdomainGraph_dhFindOwner( jxf_SubdomainGraph_dh s, JXF_Int idx, jxf_bool permuted )
{
    JXF_START_FUNC_DH
    JXF_Int sd;
    JXF_Int *beg_row = s->beg_row;
    JXF_Int *row_count = s->row_count;
    JXF_Int owner = -1, blocks = s->blocks;
    
    if (permuted) beg_row = s->beg_rowP;
    for (sd = 0; sd < blocks; ++ sd)
    {
        if (idx >= beg_row[sd] && idx < beg_row[sd]+row_count[sd])
        {
            owner = sd;
            break;
        }
    }
    if (owner == -1)
    {
        jxf_fprintf(stderr, "@@@ failed to jxf_find owner for idx = %i @@@\n", idx);
        jxf_fprintf(stderr, "blocks= %i\n", blocks);
        jxf_sprintf(jxf_msgBuf_dh, "failed to jxf_find owner for idx = %i", idx);
        JXF_SET_ERROR(-1, jxf_msgBuf_dh);
    }
    JXF_END_FUNC_VAL(owner)
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhPrintStatsLong"
void jxf_SubdomainGraph_dhPrintStatsLong( jxf_SubdomainGraph_dh s, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, k;
    JXF_Real max = 0, min = INT_MAX;
    
    jxf_fprintf(fp, "\n------------- jxf_SubdomainGraph_dhPrintStatsLong -----------\n");
    jxf_fprintf(fp, "colors used     = %i\n", s->colors);
    jxf_fprintf(fp, "subdomain count = %i\n", s->blocks);
    jxf_fprintf(fp, "\ninterior/boundary node ratios:\n");
    for (i = 0; i < s->blocks; ++ i)
    {
        JXF_Int inNodes = s->row_count[i] - s->bdry_count[i];
        JXF_Int bdNodes = s->bdry_count[i];
        JXF_Real ratio;
        
        if (bdNodes == 0)
        {
            ratio = -1;
        }
        else
        {
            ratio = (JXF_Real)inNodes / (JXF_Real)bdNodes;
        }
        max = JXF_MAX(max, ratio);
        min = JXF_MIN(min, ratio);
        jxf_fprintf(fp, "   P_%i: first= %3i  rowCount= %3i  interior= %3i  bdry= %3i  ratio= %0.1f\n",
                                          i, 1+s->beg_row[i], s->row_count[i], inNodes, bdNodes, ratio);
    }
    jxf_fprintf(fp, "\nmax interior/bdry ratio = %.1f\n", max);
    jxf_fprintf(fp, "min interior/bdry ratio = %.1f\n", min);
    if (s->adj != NULL)
    {
        jxf_fprintf(fp, "\nunpermuted subdomain graph: \n");
        for (i = 0; i < s->blocks; ++ i)
        {
            jxf_fprintf(fp, "%i :: ", i);
            for (j = s->ptrs[i]; j < s->ptrs[i+1]; ++ j)
            {
                jxf_fprintf(fp, "%i  ", s->adj[j]);
            }
            jxf_fprintf(fp, "\n");
        }
    }
    jxf_fprintf(fp, "\no2n subdomain permutation:\n");
    for (i = 0; i < s->blocks; ++ i)
    {
        jxf_fprintf(fp, "  %i %i\n", i, s->o2n_sub[i]);
    }
    jxf_fprintf(fp, "\n");
    if (jxf_np_dh > 1)
    {
        jxf_fprintf(fp, "\nlocal n2o_row permutation:\n   ");
        for (i = 0; i < s->row_count[jxf_myid_dh]; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->n2o_row[i]);
        }
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "\nlocal o2n_col permutation:\n   ");
        for (i = 0; i < s->row_count[jxf_myid_dh]; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->o2n_col[i]);
        }
        jxf_fprintf(fp, "\n");
    }
    else
    {
        jxf_fprintf(fp, "\nlocal n2o_row permutation:\n");
        jxf_fprintf(fp, "--------------------------\n");
        for (k = 0; k < s->blocks; ++ k)
        {
            JXF_Int beg_row = s->beg_row[k];
            JXF_Int end_row = beg_row + s->row_count[k];
            
            for (i = beg_row; i < end_row; ++ i)
            {
                jxf_fprintf(fp, "%i ", s->n2o_row[i]);
            }
            jxf_fprintf(fp, "\n");
        }
        jxf_fprintf(fp, "\nlocal o2n_col permutation:\n");
        jxf_fprintf(fp, "--------------------------\n");
        for (k = 0; k < s->blocks; ++ k)
        {
            JXF_Int beg_row = s->beg_row[k];
            JXF_Int end_row = beg_row + s->row_count[k];
            
            for (i = beg_row; i < end_row; ++ i)
            {
                jxf_fprintf(fp, "%i ", s->o2n_col[i]);
            }
            jxf_fprintf(fp, "\n");
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_init_seq_private"
void jxf_init_seq_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int m, n, beg_row;
    JXF_Real t1;
    
    jxf_EuclidGetDimensions(A, &beg_row, &m, &n); JXF_CHECK_V_ERROR;
    s->m = n;
    jxf_allocate_storage_private(s,blocks,m, bj); JXF_CHECK_V_ERROR;
    JXF_Int i;
    JXF_Int rpp = m / blocks;
    if (rpp*blocks < m) ++ rpp;
    s->beg_row[0] = 0;
    for (i = 1; i < blocks; ++ i) s->beg_row[i] = rpp + s->beg_row[i-1];
    for (i = 0; i < blocks; ++ i) s->row_count[i] = rpp;
    s->row_count[blocks-1] = m - rpp * (blocks - 1);
    memcpy(s->beg_rowP, s->beg_row, blocks*sizeof(JXF_Int));
    t1 = jxf_MPI_Wtime();
    if (!bj)
    {
        jxf_find_bdry_nodes_seq_private(s, m, A); JXF_CHECK_V_ERROR;
    }
    else
    {
        JXF_Int i;
        
        for (i = 0; i < m; ++ i)
        {
            s->n2o_row[i] = i;
            s->o2n_col[i] = i;
        }
    }
    s->timing[JXF_ORDER_BDRY_SGT] += (jxf_MPI_Wtime() - t1);
    t1 = jxf_MPI_Wtime();
    if (!bj)
    {
        jxf_form_subdomaingraph_seq_private(s, m, A); JXF_CHECK_V_ERROR;
        if (s->doNotColor)
        {
            JXF_Int i;
            
            jxf_printf_dh("subdomain coloring and reordering is OFF\n");
            for (i = 0; i < blocks; ++ i)
            {
                s->o2n_sub[i] = i;
                s->n2o_sub[i] = i;
                s->colorVec[i] = 0;
            }
        }
        else
        {
            JXF_SET_INFO("subdomain coloring and reordering is ON");
            jxf_color_subdomain_graph_private(s); JXF_CHECK_V_ERROR;
        }
    }
    else
    {
        JXF_Int i;
        
        for (i = 0; i < blocks; ++ i)
        {
            s->o2n_sub[i] = i;
            s->n2o_sub[i] = i;
        }
    }
    s->timing[JXF_FORM_GRAPH_SGT] += (jxf_MPI_Wtime() - t1);
    if (!bj)
    {
        jxf_adjust_matrix_perms_private(s, m); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_allocate_storage_private"
void jxf_allocate_storage_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, JXF_Int m, jxf_bool bj )
{
    JXF_START_FUNC_DH
    if (!bj)
    {
        s->ptrs = (JXF_Int *)JXF_MALLOC_DH((blocks+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        s->ptrs[0] = 0;
        s->colorVec = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        s->loNabors = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        s->hiNabors = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        s->allNabors = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    }
    s->n2o_row = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->o2n_col = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->beg_row = (JXF_Int *)JXF_MALLOC_DH((blocks)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->beg_rowP = (JXF_Int *)JXF_MALLOC_DH((blocks)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->row_count = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->bdry_count = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->o2n_sub = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    s->n2o_sub = (JXF_Int *)JXF_MALLOC_DH(blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_init_mpi_private"
void jxf_init_mpi_private( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int m, n, beg_row;
    jxf_bool symmetric;
    JXF_Real t1;
    
    symmetric = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-sym"); JXF_CHECK_V_ERROR;
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-makeSymmetric"))
    {
        symmetric = jxf_true;
    }
    jxf_EuclidGetDimensions(A, &beg_row, &m, &n); JXF_CHECK_V_ERROR;
    s->m = m;
    jxf_allocate_storage_private(s, blocks, m, bj); JXF_CHECK_V_ERROR;
    if (!bj)
    {
        jxf_MPI_Allgather(&beg_row, 1, JXF_MPI_INT, s->beg_row, 1, JXF_MPI_INT, jxf_comm_dh);
        jxf_MPI_Allgather(&m, 1, JXF_MPI_INT, s->row_count, 1, JXF_MPI_INT, jxf_comm_dh);
        memcpy(s->beg_rowP, s->beg_row, jxf_np_dh*sizeof(JXF_Int));
    }
    else
    {
        s->beg_row[jxf_myid_dh] = beg_row;
        s->beg_rowP[jxf_myid_dh] = beg_row;
        s->row_count[jxf_myid_dh] = m;
    }
    if (!bj)
    {
        t1 = jxf_MPI_Wtime();
        if (symmetric)
        {
            jxf_find_all_neighbors_sym_private(s, m, A); JXF_CHECK_V_ERROR;
        }
        else
        {
            jxf_find_all_neighbors_unsym_private(s, m, A); JXF_CHECK_V_ERROR;
        }
        s->timing[JXF_FIND_NABORS_SGT] += (jxf_MPI_Wtime() - t1);
    }
    t1 = jxf_MPI_Wtime();
    if (!bj)
    {
        JXF_Int *interiorNodes, *bdryNodes;
        JXF_Int interiorCount = 0, bdryCount;
        JXF_Int *o2n = s->o2n_col, idx;
        JXF_Int i;
        
        interiorNodes = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        bdryNodes = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        if (symmetric)
        {
            jxf_find_bdry_nodes_sym_private(s, m, A, interiorNodes,
                           bdryNodes, &interiorCount, &bdryCount); JXF_CHECK_V_ERROR;
        }
        else
        {
            jxf_find_bdry_nodes_unsym_private(s, m, A, interiorNodes,
                            bdryNodes, &interiorCount, &bdryCount); JXF_CHECK_V_ERROR;
        }
        jxf_MPI_Allgather(&bdryCount, 1, JXF_MPI_INT, s->bdry_count, 1, JXF_MPI_INT, jxf_comm_dh);
        idx = 0;
        for (i = 0; i < interiorCount; ++ i)
        {
            o2n[interiorNodes[i]] = idx ++;
        }
        for (i = 0; i < bdryCount; ++ i)
        {
            o2n[bdryNodes[i]] = idx ++;
        }
        jxf_invert_perm(m, o2n, s->n2o_row); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(interiorNodes); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(bdryNodes); JXF_CHECK_V_ERROR;
    }
    else
    {
        JXF_Int *o2n = s->o2n_col, *n2o = s->n2o_row;
        JXF_Int i, m = s->m;
        
        for (i = 0; i < m; ++ i)
        {
            o2n[i] = i;
            n2o[i] = i;
        }
    }
    s->timing[JXF_ORDER_BDRY_SGT] += (jxf_MPI_Wtime() - t1);
    if (!bj)
    {
        t1 = jxf_MPI_Wtime();
        jxf_form_subdomaingraph_mpi_private(s); JXF_CHECK_V_ERROR;
        if (s->doNotColor)
        {
            JXF_Int i;
            
            jxf_printf_dh("subdomain coloring and reordering is OFF\n");
            for (i = 0; i < blocks; ++ i)
            {
                s->o2n_sub[i] = i;
                s->n2o_sub[i] = i;
                s->colorVec[i] = 0;
            }
        }
        else
        {
            JXF_SET_INFO("subdomain coloring and reordering is ON");
            jxf_color_subdomain_graph_private(s); JXF_CHECK_V_ERROR;
        }
        s->timing[JXF_FORM_GRAPH_SGT] += (jxf_MPI_Wtime() - t1);
    }
    if (!bj)
    {
        jxf_find_ordered_neighbors_private(s); JXF_CHECK_V_ERROR;
    }
    if (!bj)
    {
        t1 = jxf_MPI_Wtime();
        jxf_SubdomainGraph_dhExchangePerms(s); JXF_CHECK_V_ERROR;
        s->timing[JXF_EXCHANGE_PERMS_SGT] += (jxf_MPI_Wtime() - t1);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhExchangePerms"
void jxf_SubdomainGraph_dhExchangePerms( jxf_SubdomainGraph_dh s )
{
    JXF_START_FUNC_DH
    MPI_Request *recv_req = NULL, *send_req = NULL;
    MPI_Status *status = NULL;
    JXF_Int *nabors = s->allNabors, naborCount = s->allCount;
    JXF_Int i, j, *sendBuf = NULL, *recvBuf = NULL, *naborIdx = NULL, nz;
    JXF_Int m = s->row_count[jxf_myid_dh];
    JXF_Int beg_row = s->beg_row[jxf_myid_dh];
    JXF_Int beg_rowP = s->beg_rowP[jxf_myid_dh];
    JXF_Int *bdryNodeCounts = s->bdry_count;
    JXF_Int myBdryCount = s->bdry_count[jxf_myid_dh];
    jxf_bool debug = jxf_false;
    JXF_Int myFirstBdry = m - myBdryCount;
    JXF_Int *n2o_row = s->n2o_row;
    jxf_Hash_i_dh n2o_table, o2n_table;
    
    if (jxf_logFile != NULL && s->debug) debug = jxf_true;
    sendBuf = (JXF_Int *)JXF_MALLOC_DH(2*myBdryCount*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nSUBG myFirstBdry= %i  myBdryCount= %i  m= %i  beg_rowP= %i\n",
                                                   1+ myFirstBdry, myBdryCount, m, 1+beg_rowP);
        fflush(jxf_logFile);
    }
    for (i = myFirstBdry, j = 0; j < myBdryCount; ++ i, ++ j)
    {
        sendBuf[2*j] = n2o_row[i] + beg_row;
        sendBuf[2*j+1] = i + beg_rowP;
    }
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "\nSUBG SEND_BUF:\n");
        for (i = myFirstBdry, j = 0; j < myBdryCount; ++ i, ++ j)
        {
            jxf_fprintf(jxf_logFile, "SUBG  %i, %i\n", 1+sendBuf[2*j], 1+sendBuf[2*j+1]);
        }
        fflush(jxf_logFile);
    }
    naborIdx = (JXF_Int *)JXF_MALLOC_DH((1+naborCount)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    naborIdx[0] = 0;
    nz = 0;
    for (i = 0; i < naborCount; ++ i)
    {
        nz += (2 * bdryNodeCounts[nabors[i]]);
        naborIdx[i+1] = nz;
    }
    recvBuf = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    recv_req = (MPI_Request*)JXF_MALLOC_DH(naborCount*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    send_req = (MPI_Request*)JXF_MALLOC_DH(naborCount*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    status = (MPI_Status*)JXF_MALLOC_DH(naborCount*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
    for (i = 0; i < naborCount; ++ i)
    {
        JXF_Int nabr = nabors[i];
        JXF_Int *buf = recvBuf + naborIdx[i];
        JXF_Int ct = 2*bdryNodeCounts[nabr];
        
        jxf_MPI_Isend(sendBuf, 2*myBdryCount, JXF_MPI_INT, nabr, 444, jxf_comm_dh, &(send_req[i]));
        if (debug)
        {
            jxf_fprintf(jxf_logFile , "SUBG   sending %i elts to %i\n", 2*myBdryCount, nabr);
            fflush(jxf_logFile);
        }
        jxf_MPI_Irecv(buf, ct, JXF_MPI_INT, nabr, 444, jxf_comm_dh, &(recv_req[i]));
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "SUBG  receiving %i elts from %i\n", ct, nabr);
            fflush(jxf_logFile);
        }
    }
    jxf_MPI_Waitall(naborCount, send_req, status);
    jxf_MPI_Waitall(naborCount, recv_req, status);
    jxf_Hash_i_dhCreate(&n2o_table, nz/2); JXF_CHECK_V_ERROR;
    jxf_Hash_i_dhCreate(&o2n_table, nz/2); JXF_CHECK_V_ERROR;
    s->n2o_ext = n2o_table;
    s->o2n_ext = o2n_table;
    for (i = 0; i < nz; i += 2)
    {
        JXF_Int old = recvBuf[i];
        JXF_Int new = recvBuf[i+1];
        
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "SUBG  i= %i  old= %i  new= %i\n", i, old+1, new+1);
            fflush(jxf_logFile);
        }
        jxf_Hash_i_dhInsert(o2n_table, old, new); JXF_CHECK_V_ERROR;
        jxf_Hash_i_dhInsert(n2o_table, new, old); JXF_CHECK_V_ERROR;
    }
    if (recvBuf != NULL)
    {
        JXF_FREE_DH(recvBuf); JXF_CHECK_V_ERROR;
    }
    if (naborIdx != NULL)
    {
        JXF_FREE_DH(naborIdx); JXF_CHECK_V_ERROR;
    }
    if (sendBuf != NULL)
    {
        JXF_FREE_DH(sendBuf); JXF_CHECK_V_ERROR;
    }
    if (recv_req != NULL)
    {
        JXF_FREE_DH(recv_req); JXF_CHECK_V_ERROR;
    }
    if (send_req != NULL)
    {
        JXF_FREE_DH(send_req); JXF_CHECK_V_ERROR;
    }
    if (status != NULL)
    {
        JXF_FREE_DH(status); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_form_subdomaingraph_mpi_private"
void jxf_form_subdomaingraph_mpi_private( jxf_SubdomainGraph_dh s )
{
    JXF_START_FUNC_DH
    JXF_Int *nabors = s->allNabors, nct = s->allCount;
    JXF_Int *idxAll = NULL;
    JXF_Int i, j, nz, *adj, *ptrs = s->ptrs;
    MPI_Request *recvReqs = NULL, sendReq;
    MPI_Status *statuses = NULL, status;
    
    if (jxf_myid_dh == 0)
    {
        idxAll = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    }
    jxf_MPI_Gather(&nct, 1, JXF_MPI_INT, idxAll, 1, JXF_MPI_INT, 0, jxf_comm_dh);
    if (jxf_myid_dh == 0)
    {
        nz = 0;
        for (i = 0; i < jxf_np_dh; ++ i) nz += idxAll[i];
    }
    jxf_MPI_Bcast(&nz, 1, JXF_MPI_INT, 0, jxf_comm_dh);
    adj = s->adj = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        recvReqs = (MPI_Request*)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
        statuses = (MPI_Status*)JXF_MALLOC_DH(jxf_np_dh*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
        ptrs[0] = 0;
        for (j = 0; j < jxf_np_dh; ++ j) ptrs[j+1] = ptrs[j] + idxAll[j];
        for (j = 0; j < jxf_np_dh; ++ j)
        {
            JXF_Int ct = idxAll[j];
            
            jxf_MPI_Irecv(adj+ptrs[j], ct, JXF_MPI_INT, j, 42, jxf_comm_dh, recvReqs+j);
        }
    }
    jxf_MPI_Isend(nabors, nct, JXF_MPI_INT, 0, 42, jxf_comm_dh, &sendReq);
    if (jxf_myid_dh == 0)
    {
        jxf_MPI_Waitall(jxf_np_dh, recvReqs, statuses);
    }
    jxf_MPI_Wait(&sendReq, &status);
    jxf_MPI_Bcast(ptrs, 1+jxf_np_dh, JXF_MPI_INT, 0, jxf_comm_dh);
    jxf_MPI_Bcast(adj, nz, JXF_MPI_INT, 0, jxf_comm_dh);
    if (idxAll != NULL)
    {
        JXF_FREE_DH(idxAll); JXF_CHECK_V_ERROR;
    }
    if (recvReqs != NULL)
    {
        JXF_FREE_DH(recvReqs); JXF_CHECK_V_ERROR;
    }
    if (statuses != NULL)
    {
        JXF_FREE_DH(statuses); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_form_subdomaingraph_seq_private"
void jxf_form_subdomaingraph_seq_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int *dense, i, j, row, blocks = s->blocks;
    JXF_Int *cval, len, *adj;
    JXF_Int idx = 0, *ptrs = s->ptrs;
    
    adj = s->adj = (JXF_Int *)JXF_MALLOC_DH(blocks*blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    dense = (JXF_Int *)JXF_MALLOC_DH(blocks*blocks*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < blocks*blocks; ++ i) dense[i] = 0;
    for (i = 0; i < blocks; ++ i)
    {
        JXF_Int beg_row = s->beg_row[i];
        JXF_Int end_row = beg_row + s->row_count[i];
        
        for (row = beg_row; row < end_row; ++ row)
        {
            jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
            for (j = 0; j < len; ++ j)
            {
                JXF_Int col = cval[j];
                
                if (col < beg_row || col >= end_row)
                {
                    JXF_Int owner = jxf_SubdomainGraph_dhFindOwner(s, col, jxf_false); JXF_CHECK_V_ERROR;
                    
                    dense[i*blocks+owner] = 1;
                    dense[owner*blocks+i] = 1;
                }
            }
            jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        }
    }
    ptrs[0] = 0;
    for (i = 0; i < blocks; ++ i)
    {
        for (j = 0; j < blocks; ++ j)
        {
            if (dense[i*blocks+j])
            {
                adj[idx++] = j;
            }
        }
        ptrs[i+1] = idx;
    }
    JXF_FREE_DH(dense); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_find_all_neighbors_sym_private"
void jxf_find_all_neighbors_sym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int *marker, i, j, beg_row, end_row;
    JXF_Int row, len, *cval, ct = 0;
    JXF_Int *nabors = s->allNabors;
    
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) marker[i] = 0;
    JXF_SET_INFO("jxf_finding nabors in subdomain graph for structurally symmetric matrix");
    JXF_SET_INFO("(if this isn't what you want, use '-sym 0' switch)");
    beg_row = s->beg_row[jxf_myid_dh];
    end_row = beg_row + s->row_count[jxf_myid_dh];
    for (row = beg_row; row < end_row; ++ row)
    {
        jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        for (j = 0; j < len; ++ j)
        {
            JXF_Int col = cval[j];
            
            if (col < beg_row || col >= end_row)
            {
                JXF_Int owner = jxf_SubdomainGraph_dhFindOwner(s, col, jxf_false); JXF_CHECK_V_ERROR;
                
                if (!marker[owner])
                {
                    marker[owner] = 1;
                    nabors[ct++] = owner;
                }
            }
        }
        jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
    }
    s->allCount = ct;
    if (marker != NULL)
    {
        JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_find_all_neighbors_unsym_private"
void jxf_find_all_neighbors_unsym_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, row, beg_row, end_row;
    JXF_Int *marker;
    JXF_Int *cval, len, idx = 0;
    JXF_Int nz, *nabors = s->allNabors, *myNabors;
    
    myNabors = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    marker = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < jxf_np_dh; ++ i) marker[i] = 0;
    JXF_SET_INFO("jxf_finding nabors in subdomain graph for structurally unsymmetric matrix");
    beg_row = s->beg_row[jxf_myid_dh];
    end_row = beg_row + s->row_count[jxf_myid_dh];
    for (row = beg_row; row < end_row; ++ row)
    {
        jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        for (j = 0; j < len; ++ j)
        {
            JXF_Int col = cval[j];
            if (col < beg_row || col >= end_row)
            {
                JXF_Int owner = jxf_SubdomainGraph_dhFindOwner(s, col, jxf_false); JXF_CHECK_V_ERROR;
                
                if (!marker[owner])
                {
                    marker[owner] = 1;
                    myNabors[idx++] = owner;
                }
            }
        }
        jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
    }
    jxf_MPI_Alltoall(marker, 1, JXF_MPI_INT, nabors, 1, JXF_MPI_INT, jxf_comm_dh); JXF_CHECK_V_ERROR;
    for (i = 0; i < idx; ++ i) nabors[myNabors[i]] = 1;
    nabors[jxf_myid_dh] = 0;
    nz = 0;
    for (i = 0; i < jxf_np_dh; ++ i)
    {
        if (nabors[i]) myNabors[nz++] = i;
    }
    s->allCount = nz;
    memcpy(nabors, myNabors, nz*sizeof(JXF_Int));
    if (marker != NULL)
    {
        JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    }
    if (myNabors != NULL)
    {
        JXF_FREE_DH(myNabors); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_find_bdry_nodes_sym_private"
void jxf_find_bdry_nodes_sym_private( jxf_SubdomainGraph_dh s,
                                  JXF_Int m,
                                  void *A,
                                  JXF_Int *interiorNodes,
                                  JXF_Int *bdryNodes,
                                  JXF_Int *interiorCount,
                                  JXF_Int *bdryCount )
{
    JXF_START_FUNC_DH
    JXF_Int beg_row = s->beg_row[jxf_myid_dh];
    JXF_Int end_row = beg_row + s->row_count[jxf_myid_dh];
    JXF_Int row, inCt = 0, bdCt = 0;
    JXF_Int j;
    JXF_Int *cval;
    
    for (row = beg_row; row < end_row; ++ row)
    {
        jxf_bool isBdry = jxf_false;
        JXF_Int len;
        
        jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        for (j = 0; j < len; ++ j)
        {
            JXF_Int col = cval[j];
            
            if (col < beg_row || col >= end_row)
            {
                isBdry = jxf_true;
                break;
            }
        }
        jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        if (isBdry)
        {
            bdryNodes[bdCt++] = row - beg_row;
        }
        else
        {
            interiorNodes[inCt++] = row - beg_row;
        }
    }
   *interiorCount = inCt;
   *bdryCount = bdCt;
    JXF_END_FUNC_DH
}

#define JXF_BDRY_NODE_TAG 42

#undef __FUNC__
#define __FUNC__ "jxf_find_bdry_nodes_unsym_private"
void jxf_find_bdry_nodes_unsym_private( jxf_SubdomainGraph_dh s,
                                    JXF_Int m,
                                    void *A,
                                    JXF_Int *interiorNodes,
                                    JXF_Int *boundaryNodes,
                                    JXF_Int *interiorCount,
                                    JXF_Int *bdryCount )
{
    JXF_START_FUNC_DH
    JXF_Int beg_row = s->beg_row[jxf_myid_dh];
    JXF_Int end_row = beg_row + s->row_count[jxf_myid_dh];
    JXF_Int i, j, row, max;
    JXF_Int *cval;
    JXF_Int *list, count;
    JXF_Int *rpIN = NULL, *rpOUT = NULL;
    JXF_Int *sendBuf, *recvBuf;
    JXF_Int *marker, inCt, bdCt;
    JXF_Int *bdryNodes, nz;
    JXF_Int sendCt, recvCt;
    MPI_Request *sendReq, *recvReq;
    MPI_Status *status;
    jxf_SortedSet_dh ss;
    
    jxf_SortedSet_dhCreate(&ss, m); JXF_CHECK_V_ERROR;
    for (row = beg_row; row < end_row; ++ row)
    {
        jxf_bool isBdry = jxf_false;
        JXF_Int len;
        
        jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        for (j = 0; j < len; ++ j)
        {
            JXF_Int col = cval[j];
            
            if (col < beg_row || col >= end_row)
            {
                isBdry = jxf_true;
                jxf_SortedSet_dhInsert(ss, col); JXF_CHECK_V_ERROR;
            }
        }
        jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        if (isBdry)
        {
            jxf_SortedSet_dhInsert(ss, row); JXF_CHECK_V_ERROR;
        }
    }
    sendBuf = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    recvBuf = (JXF_Int *)JXF_MALLOC_DH(jxf_np_dh*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rpOUT = (JXF_Int *)JXF_MALLOC_DH((jxf_np_dh+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rpOUT[0] = 0;
    for (i = 0; i < jxf_np_dh; ++ i) sendBuf[i] = 0;
    sendCt = 0;
    jxf_SortedSet_dhGetList(ss, &list, &count); JXF_CHECK_V_ERROR;
    for (i = 0; i < count;)
    {
        JXF_Int node = list[i];
        JXF_Int owner;
        JXF_Int last;
        
        owner = jxf_SubdomainGraph_dhFindOwner(s, node, jxf_false); JXF_CHECK_V_ERROR;
        last = s->beg_row[owner] + s->row_count[owner];
        while ((i < count)  && (list[i] < last)) ++ i;
        ++ sendCt;
        rpOUT[sendCt] = i;
        sendBuf[owner] = rpOUT[sendCt] - rpOUT[sendCt-1];
    }
    jxf_MPI_Alltoall(sendBuf, 1, JXF_MPI_INT, recvBuf, 1, JXF_MPI_INT, jxf_comm_dh); JXF_CHECK_V_ERROR;
    rpIN = (JXF_Int *)JXF_MALLOC_DH((jxf_np_dh+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rpIN[0] = 0;
    nz = 0;
    recvCt = 0;
    for (i = 0; i < jxf_np_dh; ++ i)
    {
        if (recvBuf[i])
        {
            ++ recvCt;
            nz += recvBuf[i];
            rpIN[recvCt] = nz;
        }
    }
    bdryNodes = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    sendReq = (MPI_Request *)JXF_MALLOC_DH(sendCt*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    recvReq = (MPI_Request *)JXF_MALLOC_DH(recvCt*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    max = JXF_MAX(sendCt, recvCt);
    status = (MPI_Status *)JXF_MALLOC_DH(max*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
    j = 0;
    for (i = 0; i < jxf_np_dh; ++ i)
    {
        if (recvBuf[i])
        {
            jxf_MPI_Irecv(bdryNodes+rpIN[j], recvBuf[i], JXF_MPI_INT, i, JXF_BDRY_NODE_TAG, jxf_comm_dh, recvReq+j);
            ++ j;
        }
    }
    j = 0;
    for (i = 0; i < jxf_np_dh; ++ i)
    {
        if (sendBuf[i])
        {
            jxf_MPI_Isend(list+rpOUT[j], sendBuf[i], JXF_MPI_INT, i, JXF_BDRY_NODE_TAG, jxf_comm_dh, sendReq+j);
            ++ j;
        }
    }
    jxf_MPI_Waitall(sendCt, sendReq, status);
    jxf_MPI_Waitall(recvCt, recvReq, status);
    for (i = 0; i < nz; ++ i) bdryNodes[i] -= beg_row;
    marker = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) marker[i] = 0;
    for (i = 0; i < nz; ++ i) marker[bdryNodes[i]] = 1;
    inCt = bdCt = 0;
    for (i = 0; i < m; ++ i)
    {
        if (marker[i])
        {
            boundaryNodes[bdCt++] = i;
        }
        else
        {
            interiorNodes[inCt++] = i;
        }
    }
   *interiorCount = inCt;
   *bdryCount = bdCt;
    jxf_SortedSet_dhDestroy(ss); JXF_CHECK_V_ERROR;
    if (rpIN != NULL)
    {
        JXF_FREE_DH(rpIN); JXF_CHECK_V_ERROR;
    }
    if (rpOUT != NULL)
    {
        JXF_FREE_DH(rpOUT); JXF_CHECK_V_ERROR;
    }
    if (sendBuf != NULL)
    {
        JXF_FREE_DH(sendBuf); JXF_CHECK_V_ERROR;
    }
    if (recvBuf != NULL)
    {
        JXF_FREE_DH(recvBuf); JXF_CHECK_V_ERROR;
    }
    if (bdryNodes != NULL)
    {
        JXF_FREE_DH(bdryNodes); JXF_CHECK_V_ERROR;
    }
    if (marker != NULL)
    {
        JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    }
    if (sendReq!= NULL)
    {
        JXF_FREE_DH(sendReq); JXF_CHECK_V_ERROR;
    }
    if (recvReq!= NULL)
    {
        JXF_FREE_DH(recvReq); JXF_CHECK_V_ERROR;
    }
    if (status!= NULL)
    {
        JXF_FREE_DH(status); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_find_ordered_neighbors_private"
void jxf_find_ordered_neighbors_private( jxf_SubdomainGraph_dh s )
{
    JXF_START_FUNC_DH
    JXF_Int *loNabors = s->loNabors;
    JXF_Int *hiNabors = s->hiNabors;
    JXF_Int *allNabors = s->allNabors, allCount = s->allCount;
    JXF_Int loCt = 0, hiCt = 0;
    JXF_Int *o2n = s->o2n_sub;
    JXF_Int i, myNewId = o2n[jxf_myid_dh];
    
    for (i = 0; i < allCount; ++ i)
    {
        JXF_Int nabor = allNabors[i];
        
        if (o2n[nabor] < myNewId)
        {
            loNabors[loCt++] = nabor;
        }
        else
        {
            hiNabors[hiCt++] = nabor;
        }
    }
    s->loCount = loCt;
    s->hiCount = hiCt;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_color_subdomain_graph_private"
void jxf_color_subdomain_graph_private( jxf_SubdomainGraph_dh s )
{
    JXF_START_FUNC_DH
    JXF_Int i, n = jxf_np_dh;
    JXF_Int *rp = s->ptrs, *cval = s->adj;
    JXF_Int j, *marker, thisNodesColor, *colorCounter;
    JXF_Int *o2n = s->o2n_sub;
    JXF_Int *color = s->colorVec;
    
    if (jxf_np_dh == 1) n = s->blocks;
    marker = (JXF_Int *)JXF_MALLOC_DH((n+1)*sizeof(JXF_Int));
    colorCounter = (JXF_Int *)JXF_MALLOC_DH((n+1)*sizeof(JXF_Int));
    for (i = 0; i <= n; ++ i)
    {
        marker[i] = -1;
        colorCounter[i] = 0;
    }
    for (i = 0; i < n; ++ i)
    {
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int nabor = cval[j];
            if (nabor < i)
            {
                JXF_Int naborsColor = color[nabor];
                marker[naborsColor] = i;
            }
        }
        thisNodesColor = -1;
        for (j = 0; j < n; ++ j)
        {
            if (marker[j] != i)
            {
                thisNodesColor = j;
                break;
            }
        }
        color[i] = thisNodesColor;
        colorCounter[1+thisNodesColor] += 1;
    }
    for (i = 1; i < n; ++ i)
    {
        if (colorCounter[i] == 0) break;
        colorCounter[i] += colorCounter[i-1];
    }
    for (i = 0; i < n; ++ i)
    {
        o2n[i] = colorCounter[color[i]];
        colorCounter[color[i]] += 1;
    }
    jxf_invert_perm(n, s->o2n_sub, s->n2o_sub); JXF_CHECK_V_ERROR;
    JXF_Int ct = 0;
    for (j = 0; j < n; ++ j)
    {
        if (marker[j] == -1) break;
        ++ ct;
    }
    s->colors = ct;
    JXF_Int sum = 0;
    for (i = 0; i < n; ++ i)
    {
        JXF_Int old = s->n2o_sub[i];
        
        s->beg_rowP[old] = sum;
        sum += s->row_count[old];
    }
    JXF_FREE_DH(marker); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(colorCounter); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhDump"
void jxf_SubdomainGraph_dhDump( jxf_SubdomainGraph_dh s, char *filename )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int sCt = jxf_np_dh;
    FILE *fp;
    
    if (jxf_np_dh == 1) sCt = s->blocks;
    fp = jxf_openFile_dh(filename, "w"); JXF_CHECK_V_ERROR;
    jxf_fprintf(fp, "----- colors used\n");
    jxf_fprintf(fp, "%i\n", s->colors);
    if (s->colorVec == NULL)
    {
        jxf_fprintf(fp, "s->colorVec == NULL\n");
    }
    else
    {
        jxf_fprintf(fp, "----- colorVec\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->colorVec[i]);
        }
        jxf_fprintf(fp, "\n");
    }
    if (s->o2n_sub == NULL || s->o2n_sub == NULL)
    {
        jxf_fprintf(fp, "s->o2n_sub == NULL || s->o2n_sub == NULL\n");
    }
    else
    {
        jxf_fprintf(fp, "----- o2n_sub\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->o2n_sub[i]);
        }
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "----- n2o_sub\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->n2o_sub[i]);
        }
        jxf_fprintf(fp, "\n");
    }
    if (s->beg_row == NULL || s->beg_rowP == NULL)
    {
        jxf_fprintf(fp, "s->beg_row == NULL || s->beg_rowP == NULL\n");
    }
    else
    {
        jxf_fprintf(fp, "----- beg_row\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", 1+s->beg_row[i]);
        }
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "----- beg_rowP\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", 1+s->beg_rowP[i]);
        }
        jxf_fprintf(fp, "\n");
    }
    if (s->row_count == NULL || s->bdry_count == NULL)
    {
        jxf_fprintf(fp, "s->row_count == NULL || s->bdry_count == NULL\n");
    }
    else
    {
        jxf_fprintf(fp, "----- row_count\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->row_count[i]);
        }
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "----- bdry_count\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i ", s->bdry_count[i]);
        }
        jxf_fprintf(fp, "\n");
    }
    if (s->ptrs == NULL || s->adj == NULL)
    {
        jxf_fprintf(fp, "s->ptrs == NULL || s->adj == NULL\n");
    }
    else
    {
        JXF_Int j;
        JXF_Int ct;
        
        jxf_fprintf(fp, "----- subdomain graph\n");
        for (i = 0; i < sCt; ++ i)
        {
            jxf_fprintf(fp, "%i :: ", i);
            ct = s->ptrs[i+1] - s->ptrs[i];
            if (ct)
            {
                jxf_shellSort_int(ct, s->adj+s->ptrs[i]); JXF_CHECK_V_ERROR;
            }
            for (j = s->ptrs[i]; j < s->ptrs[i+1]; ++ j)
            {
                jxf_fprintf(fp, "%i ", s->adj[j]);
            }
            jxf_fprintf(fp, "\n");
        }
    }
    jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    if (s->beg_rowP == NULL)
    {
        JXF_SET_V_ERROR("s->beg_rowP == NULL; can't continue");
    }
    if (s->row_count == NULL)
    {
        JXF_SET_V_ERROR("s->row_count == NULL; can't continue");
    }
    if (s->o2n_sub == NULL)
    {
        JXF_SET_V_ERROR("s->o2n_sub == NULL; can't continue");
    }
    if (jxf_np_dh == 1)
    {
        fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
        if (s->n2o_row == NULL|| s->o2n_col == NULL)
        {
            jxf_fprintf(fp, "s->n2o_row == NULL|| s->o2n_col == NULL\n");
        }
        else
        {
            jxf_fprintf(fp, "----- n2o_row\n");
            for (i = 0; i < s->m; ++ i)
            {
                jxf_fprintf(fp, "%i ", 1+s->n2o_row[i]);
            }
            jxf_fprintf(fp, "\n");
        }
        jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
    }
    else
    {
        JXF_Int id = s->n2o_sub[jxf_myid_dh];
        JXF_Int m = s->m;
        JXF_Int pe;
        JXF_Int beg_row = 0;
        
        if (s->beg_row != 0) beg_row = s->beg_row[jxf_myid_dh];
        for (pe = 0; pe < jxf_np_dh; ++ pe)
        {
            jxf_MPI_Barrier(jxf_comm_dh);
            if (id == pe)
            {
                fp = jxf_openFile_dh(filename, "a"); JXF_CHECK_V_ERROR;
                if (id == 0) jxf_fprintf(fp, "----- n2o_row\n");
                for (i = 0; i < m; ++ i)
                {
                    jxf_fprintf(fp, "%i ", 1+s->n2o_row[i]+beg_row);
                }
                if (id == jxf_np_dh - 1) jxf_fprintf(fp, "\n");
                jxf_closeFile_dh(fp); JXF_CHECK_V_ERROR;
            }
        }
    }
    JXF_END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "jxf_find_bdry_nodes_seq_private"
void jxf_find_bdry_nodes_seq_private( jxf_SubdomainGraph_dh s, JXF_Int m, void *A )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, row, blocks = s->blocks;
    JXF_Int *cval, *tmp;
    
    tmp = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) tmp[i] = 0;
    for (i = 0; i < blocks; ++ i)
    {
        JXF_Int beg_row = s->beg_row[i];
        JXF_Int end_row = beg_row + s->row_count[i];
        
        for (row = beg_row; row < end_row; ++ row)
        {
            jxf_bool isBdry = jxf_false;
            JXF_Int len;
            
            jxf_EuclidGetRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
            for (j = 0; j < len; ++ j)
            {
                JXF_Int col = cval[j];
                
                if (col < beg_row || col >= end_row)
                {
                    tmp[col] = 1;
                    isBdry = jxf_true;
                }
            }
            if (isBdry) tmp[row] = 1;
            jxf_EuclidRestoreRow(A, row, &len, &cval, NULL); JXF_CHECK_V_ERROR;
        }
    }
    for (i = 0; i < blocks; ++ i)
    {
        JXF_Int beg_row = s->beg_row[i];
        JXF_Int end_row = beg_row + s->row_count[i];
        JXF_Int ct = 0;
        
        for (row = beg_row; row < end_row; ++ row)
        {
            if (tmp[row]) ++ ct;
        }
        s->bdry_count[i] = ct;
    }
    for (i = 0; i < blocks; ++ i)
    {
        JXF_Int beg_row = s->beg_row[i];
        JXF_Int end_row = beg_row + s->row_count[i];
        JXF_Int interiorIDX = beg_row;
        JXF_Int bdryIDX = end_row - s->bdry_count[i];
        
        for (row = beg_row; row < end_row; ++ row)
        {
            if (tmp[row])
            {
                s->o2n_col[row] = bdryIDX ++;
            }
            else
            {
                s->o2n_col[row] = interiorIDX ++;
            }
        }
    }
    jxf_invert_perm(m, s->o2n_col, s->n2o_row); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(tmp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhPrintSubdomainGraph"
void jxf_SubdomainGraph_dhPrintSubdomainGraph( jxf_SubdomainGraph_dh s, FILE *fp )
{
    JXF_START_FUNC_DH
    if (jxf_myid_dh == 0)
    {
        JXF_Int i, j;
        
        jxf_fprintf(fp, "\n-----------------------------------------------------\n");
        jxf_fprintf(fp, "SubdomainGraph, and coloring and ordering information\n");
        jxf_fprintf(fp, "-----------------------------------------------------\n");
        jxf_fprintf(fp, "colors used: %i\n", s->colors);
        jxf_fprintf(fp, "o2n ordering vector: ");
        for (i = 0; i < s->blocks; ++ i) jxf_fprintf(fp, "%i ", s->o2n_sub[i]);
        jxf_fprintf(fp, "\ncoloring vector (node, color): \n");
        for (i = 0; i < s->blocks; ++ i) jxf_fprintf(fp, "  %i, %i\n", i, s->colorVec[i]);
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "Adjacency lists:\n");
        for (i = 0; i < s->blocks; ++ i)
        {
            jxf_fprintf(fp, "   P_%i :: ", i);
            for (j = s->ptrs[i]; j < s->ptrs[i+1]; ++ j)
            {
                jxf_fprintf(fp, "%i ", s->adj[j]);
            }
            jxf_fprintf(fp, "\n");
        }
        jxf_fprintf(fp, "-----------------------------------------------------\n");
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_adjust_matrix_perms_private"
void jxf_adjust_matrix_perms_private( jxf_SubdomainGraph_dh s, JXF_Int m )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, blocks = s->blocks;
    JXF_Int *o2n = s->o2n_col;
    
    for (i = 0; i < blocks; ++ i)
    {
        JXF_Int beg_row = s->beg_row[i];
        JXF_Int end_row = beg_row + s->row_count[i];
        JXF_Int adjust = s->beg_rowP[i] - s->beg_row[i];
        
        for (j = beg_row; j < end_row; ++ j) o2n[j] += adjust;
    }
    jxf_invert_perm(m, s->o2n_col, s->n2o_row); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhPrintRatios"
void jxf_SubdomainGraph_dhPrintRatios( jxf_SubdomainGraph_dh s, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int blocks = jxf_np_dh;
    JXF_Real ratio[25];
    
    if (jxf_myid_dh == 0)
    {
        if (jxf_np_dh == 1) blocks = s->blocks;
        if (blocks > 25) blocks = 25;
        jxf_fprintf(fp, "\n");
        jxf_fprintf(fp, "Subdomain interior/boundary node ratios\n");
        jxf_fprintf(fp, "---------------------------------------\n");
        for (i = 0; i < blocks; ++ i)
        {
            if (s->bdry_count[i] == 0)
            {
                ratio[i] = -1;
            }
            else
            {
                ratio[i] = (JXF_Real)(s->row_count[i] - s->bdry_count[i]) / (JXF_Real)s->bdry_count[i];
            }
        }
        jxf_shellSort_float(blocks, ratio);
        if (blocks <= 20)
        {
            JXF_Int j = 0;
            
            for (i = 0; i < blocks; ++ i)
            {
                jxf_fprintf(fp, "%0.2g  ", ratio[i]);
                ++ j;
                if (j == 10)
                {
                    jxf_fprintf(fp, "\n");
                }
            }
            jxf_fprintf(fp, "\n");
        }
        else
        {
            jxf_fprintf(fp, "10 smallest ratios: ");
            for (i = 0; i < 10; ++ i)
            {
                jxf_fprintf(fp, "%0.2g  ", ratio[i]);
            }
            jxf_fprintf(fp, "\n");
            jxf_fprintf(fp, "10 largest ratios:  ");
            JXF_Int start = blocks - 6, stop = blocks - 1;
            for (i = start; i < stop; ++ i)
            {
                jxf_fprintf(fp, "%0.2g  ", ratio[i]);
            }
            jxf_fprintf(fp, "\n");
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_SubdomainGraph_dhPrintStats"
void jxf_SubdomainGraph_dhPrintStats( jxf_SubdomainGraph_dh sg, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Real *timing = sg->timing;
    
    jxf_fprintf_dh(fp, "\nSubdomainGraph timing report\n");
    jxf_fprintf_dh(fp, "-----------------------------\n");
    jxf_fprintf_dh(fp, "total setup time: %0.2f\n", timing[JXF_TOTAL_SGT]);
    jxf_fprintf_dh(fp, "  jxf_find neighbors in subdomain graph: %0.2f\n", timing[JXF_FIND_NABORS_SGT]);
    jxf_fprintf_dh(fp, "  locally order interiors and bdry:  %0.2f\n", timing[JXF_ORDER_BDRY_SGT]);
    jxf_fprintf_dh(fp, "  form and color subdomain graph:    %0.2f\n", timing[JXF_FORM_GRAPH_SGT]);
    jxf_fprintf_dh(fp, "  exchange bdry permutations:        %0.2f\n", timing[JXF_EXCHANGE_PERMS_SGT]);
    jxf_fprintf_dh(fp, "  everything else (should be small): %0.2f\n", timing[JXF_TOTAL_SGT]-(timing[JXF_FIND_NABORS_SGT]+
                                timing[JXF_ORDER_BDRY_SGT]+timing[JXF_FORM_GRAPH_SGT]+timing[JXF_EXCHANGE_PERMS_SGT]));
    JXF_END_FUNC_DH
}
