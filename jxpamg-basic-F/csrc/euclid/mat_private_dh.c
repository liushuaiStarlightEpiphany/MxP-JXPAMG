//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  mat_private_dh.c
 *  Date: 2013/01/23
 */

#include "jxf_euclid.h"

#define JXF_IS_UPPER_TRI 97
#define JXF_IS_LOWER_TRI 98
#define JXF_IS_FULL 99

static JXF_Int jxf_isTriangular( JXF_Int m, JXF_Int *rp, JXF_Int *cval );
static void jxf_mat_par_read_allocate_private( jxf_Mat_dh *Aout, JXF_Int n, JXF_Int *rowLengths, JXF_Int *rowToBlock );
void jxf_mat_partition_private( jxf_Mat_dh A, JXF_Int blocks, JXF_Int *o2n_row, JXF_Int *rowToBlock );
static void jxf_convert_triples_to_scr_private( JXF_Int m, JXF_Int nz, JXF_Int *I, JXF_Int *J,
                                 JXF_Real *A, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval );

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_print_graph_private"
void jxf_mat_dh_print_graph_private( JXF_Int m,
                                 JXF_Int beg_row,
                                 JXF_Int *rp,
                                 JXF_Int *cval,
                                 JXF_Real *aval,
                                 JXF_Int *n2o,
                                 JXF_Int *o2n,
                                 jxf_Hash_i_dh hash,
                                 FILE* fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, row, col;
    jxf_bool private_n2o = jxf_false;
    jxf_bool private_hash = jxf_false;
    JXF_Int *work = NULL;
    
    work = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (n2o == NULL)
    {
        private_n2o = jxf_true;
        jxf_create_nat_ordering_private(m, &n2o); JXF_CHECK_V_ERROR;
        jxf_create_nat_ordering_private(m, &o2n); JXF_CHECK_V_ERROR;
    }
    if (hash == NULL)
    {
        private_hash = jxf_true;
        jxf_Hash_i_dhCreate(&hash, -1); JXF_CHECK_V_ERROR;
    }
    for (i = 0; i < m; ++ i)
    {
        for (j = 0; j < m; ++ j) work[j] = 0;
        row = n2o[i];
        for (j = rp[row]; j < rp[row+1]; ++ j)
        {
            col = cval[j];
            if (col >= beg_row || col < beg_row+m)
            {
                col = o2n[col];
            }
            else
            {
                JXF_Int tmp = col;
                
                tmp = jxf_Hash_i_dhLookup(hash, col); JXF_CHECK_V_ERROR;
                if (tmp == -1)
                {
                    jxf_sprintf(jxf_msgBuf_dh, "beg_row= %i  m= %i; nonlocal column= %i not in hash table",
                                                     beg_row, m, col);
                    JXF_SET_V_ERROR(jxf_msgBuf_dh);
                }
                else
                {
                    col = tmp;
                }
            }
            work[col] = 1;
        }
        for (j = 0; j < m; ++ j)
        {
            if (work[j])
            {
                jxf_fprintf(fp, " x ");
            }
            else
            {
                jxf_fprintf(fp, "   ");
            }
        }
        jxf_fprintf(fp, "\n");
    }
    if (private_n2o)
    {
        jxf_destroy_nat_ordering_private(n2o); JXF_CHECK_V_ERROR;
        jxf_destroy_nat_ordering_private(o2n); JXF_CHECK_V_ERROR;
    }
    if (private_hash)
    {
        jxf_Hash_i_dhDestroy(hash); JXF_CHECK_V_ERROR;
    }
    if (work != NULL)
    {
        JXF_FREE_DH(work); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_create_nat_ordering_private"
void jxf_create_nat_ordering_private( JXF_Int m, JXF_Int **p )
{
    JXF_START_FUNC_DH
    JXF_Int *tmp, i;
    
    tmp = *p = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        tmp[i] = i;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_destroy_nat_ordering_private"
void jxf_destroy_nat_ordering_private( JXF_Int *p )
{
    JXF_START_FUNC_DH
    JXF_FREE_DH(p); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_invert_perm"
void jxf_invert_perm( JXF_Int m, JXF_Int *pIN, JXF_Int *pOUT )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
    for (i = 0; i < m; ++ i) pOUT[pIN[i]] = i;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_print_csr_private"
void jxf_mat_dh_print_csr_private( JXF_Int m, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, nz = rp[m];
    
    jxf_fprintf(fp, "%i %i\n", m, rp[m]);
    for (i = 0; i <= m; ++ i) jxf_fprintf(fp, "%i ", rp[i]);
    jxf_fprintf(fp, "\n");
    for (i = 0; i < nz; ++ i) jxf_fprintf(fp, "%i ", cval[i]);
    jxf_fprintf(fp, "\n");
    for (i = 0; i < nz; ++ i) jxf_fprintf(fp, "%1.19e ", aval[i]);
    jxf_fprintf(fp, "\n");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_read_csr_private"
void jxf_mat_dh_read_csr_private( JXF_Int *mOUT, JXF_Int **rpOUT, JXF_Int **cvalOUT, JXF_Real **avalOUT, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int i, m, nz, items;
    JXF_Int *rp, *cval;
    JXF_Real *aval;
    
    items = jxf_fscanf(fp,"%d %d",&m, &nz);
    if (items != 2)
    {
        JXF_SET_V_ERROR("failed to read header");
    }
    else
    {
        jxf_printf("jxf_mat_dh_read_csr_private:: m= %i  nz= %i\n", m, nz);
    }
   *mOUT = m;
    rp = *rpOUT = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = *cvalOUT = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    aval = *avalOUT = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    for (i = 0; i <= m; ++ i)
    {
        items = jxf_fscanf(fp,"%d", &(rp[i]));
        if (items != 1)
        {
            jxf_sprintf(jxf_msgBuf_dh, "failed item %i of %i in rp block", i, m+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    for (i = 0; i < nz; ++ i)
    {
        items = jxf_fscanf(fp,"%d", &(cval[i]));
        if (items != 1)
        {
            jxf_sprintf(jxf_msgBuf_dh, "failed item %i of %i in cval block", i, m+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    for (i = 0; i < nz; ++ i)
    {
        items = jxf_fscanf(fp,"%lg", &(aval[i]));
        if (items != 1)
        {
            jxf_sprintf(jxf_msgBuf_dh, "failed item %i of %i in aval block", i, m+1);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
    }
    JXF_END_FUNC_DH
}

#define JXF_MAX_JUNK 200

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_read_triples_private"
void jxf_mat_dh_read_triples_private( JXF_Int ignore, JXF_Int *mOUT, JXF_Int **rpOUT, JXF_Int **cvalOUT, JXF_Real **avalOUT, FILE *fp )
{
    JXF_START_FUNC_DH
    JXF_Int m, n, nz, items, i, j;
    JXF_Int idx = 0;
    JXF_Int *cval, *rp, *I, *J;
    JXF_Real *aval, *A, v;
    char junk[JXF_MAX_JUNK];
    fpos_t fpos;
    
    if (ignore && jxf_myid_dh == 0)
    {
        jxf_printf("jxf_mat_dh_read_triples_private:: ignoring following header lines:\n");
        jxf_printf("--------------------------------------------------------------\n");
        for (i = 0; i < ignore; ++ i)
        {
            fgets(junk, JXF_MAX_JUNK, fp);
            jxf_printf("%s", junk);
        }
        jxf_printf("--------------------------------------------------------------\n");
        if (fgetpos(fp, &fpos)) JXF_SET_V_ERROR("fgetpos failed!");
        jxf_printf("\njxf_mat_dh_read_triples_private::1st two non-ignored lines:\n");
        jxf_printf("--------------------------------------------------------------\n");
        for (i = 0; i < 2; ++ i)
        {
            fgets(junk, JXF_MAX_JUNK, fp);
            jxf_printf("%s", junk);
        }
        jxf_printf("--------------------------------------------------------------\n");
        if (fsetpos(fp, &fpos)) JXF_SET_V_ERROR("fsetpos failed!");
    }
    if (feof(fp)) jxf_printf("trouble!");
    m = n = nz = 0;
    while (!feof(fp))
    {
        items = jxf_fscanf(fp,"%d %d %lg", &i, &j, &v);
        if (items != 3)
        {
            break;
        }
        ++ nz;
        if (i > m) m = i;
        if (j > n) n = j;
    }
    if (jxf_myid_dh == 0)
    {
        jxf_printf("jxf_mat_dh_read_triples_private: m= %i  nz= %i\n", m, nz);
    }
    rewind(fp);
    for (i = 0; i < ignore; ++ i)
    {
        fgets(junk, JXF_MAX_JUNK, fp);
    }
    if (m != n)
    {
        jxf_sprintf(jxf_msgBuf_dh, "matrix is not square; row= %i, cols= %i", m, n);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
   *mOUT = m;
    rp = *rpOUT = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = *cvalOUT = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    aval = *avalOUT = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    I = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    J = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    A = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    while (!feof(fp))
    {
        items = jxf_fscanf(fp,"%d %d %lg", &i, &j, &v);
        if (items < 3) break;
        j --;
        i --;
        I[idx] = i;
        J[idx] = j;
        A[idx] = v;
        ++ idx;
    }
    jxf_convert_triples_to_scr_private(m, nz, I, J, A, rp, cval, aval); JXF_CHECK_V_ERROR;
    JXF_Int type;
    type = jxf_isTriangular(m, rp, cval); JXF_CHECK_V_ERROR;
    if (type == JXF_IS_UPPER_TRI)
    {
        jxf_printf("CAUTION: matrix is upper triangular; converting to full\n");
    }
    else if (type == JXF_IS_LOWER_TRI)
    {
        jxf_printf("CAUTION: matrix is lower triangular; converting to full\n");
    }
    if (type == JXF_IS_UPPER_TRI || type == JXF_IS_LOWER_TRI)
    {
        jxf_make_full_private(m, &rp, &cval, &aval); JXF_CHECK_V_ERROR;
    }
   *rpOUT = rp;
   *cvalOUT = cval;
   *avalOUT = aval;
    JXF_FREE_DH(I); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(J); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(A); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_convert_triples_to_scr_private"
void jxf_convert_triples_to_scr_private( JXF_Int m, JXF_Int nz, JXF_Int *I, JXF_Int *J, JXF_Real *A, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    JXF_Int *rowCounts;
    
    rowCounts = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) rowCounts[i] = 0;
    for (i = 0; i < nz; ++ i)
    {
        JXF_Int row = I[i];
        
        rowCounts[row] += 1;
    }
    rp[0] = 0;
    for (i = 1; i <= m; ++ i)
    {
        rp[i] = rp[i-1] + rowCounts[i-1];
    }
    memcpy(rowCounts, rp, (m+1)*sizeof(JXF_Int));
    for (i = 0; i < nz; ++ i)
    {
        JXF_Int row = I[i];
        JXF_Int col = J[i];
        JXF_Real val = A[i];
        JXF_Int idx = rowCounts[row];
        
        rowCounts[row] += 1;
        cval[idx] = col;
        aval[idx] = val;
    }
    JXF_FREE_DH(rowCounts); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

void jxf_fix_diags_private( jxf_Mat_dh A );
void jxf_insert_missing_diags_private( jxf_Mat_dh A );

#undef __FUNC__
#define __FUNC__ "jxf_readMat"
void jxf_readMat( jxf_Mat_dh *Aout, char *ft, char *fn, JXF_Int ignore )
{
    JXF_START_FUNC_DH
    jxf_bool makeStructurallySymmetric;
    jxf_bool fixDiags;
    
   *Aout = NULL;
    makeStructurallySymmetric = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-makeSymmetric");
    fixDiags = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-fixDiags");
    if (fn == NULL)
    {
        JXF_SET_V_ERROR("passed NULL filename; can't open for reading!");
    }
    if (!strcmp(ft, "csr"))
    {
        jxf_Mat_dhReadCSR(Aout, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "trip"))
    {
        jxf_Mat_dhReadTriples(Aout, ignore, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "ebin"))
    {
        jxf_Mat_dhReadBIN(Aout, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "petsc"))
    {
        jxf_sprintf(jxf_msgBuf_dh, "must recompile Euclid using petsc mode!");
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "unknown filetype: -ftin %s", ft);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    if (makeStructurallySymmetric)
    {
        jxf_printf("\npadding with zeros to make structurally symmetric\n");
        jxf_Mat_dhMakeStructurallySymmetric(*Aout); JXF_CHECK_V_ERROR;
    }
    if ( (*Aout)->m == 0)
    {
        JXF_SET_V_ERROR("row count = 0; something's wrong!");
    }
    if (fixDiags)
    {
        jxf_fix_diags_private(*Aout); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_fix_diags_private"
void jxf_fix_diags_private( jxf_Mat_dh A )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, m = A->m, *rp = A->rp, *cval = A->cval;
    JXF_Real *aval = A->aval;
    jxf_bool insertDiags = jxf_false;
    
    for (i = 0; i < m; ++ i)
    {
        jxf_bool isMissing = jxf_true;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            if (cval[j] == i)
            {
                isMissing = jxf_false;
                break;
            }
        }
        if (isMissing)
        {
            insertDiags = jxf_true;
            break;
        }
    }
    if (insertDiags)
    {
        jxf_insert_missing_diags_private(A); JXF_CHECK_V_ERROR;
        rp = A->rp;
        cval = A->cval;
        aval = A->aval;
    }
    for (i = 0; i < m; ++ i)
    {
        JXF_Real sum = 0;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            sum = JXF_MAX(sum, fabs(aval[j]));
        }
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            if (cval[j] == i)
            {
                aval[j] = sum;
                break;
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_insert_missing_diags_private"
void jxf_insert_missing_diags_private( jxf_Mat_dh A )
{
    JXF_START_FUNC_DH
    JXF_Int *RP = A->rp, *CVAL = A->cval, m = A->m;
    JXF_Int *rp, *cval;
    JXF_Real *AVAL = A->aval, *aval;
    JXF_Int i, j, nz = RP[m] + m;
    JXF_Int idx = 0;
    
    rp = A->rp = (JXF_Int *)JXF_MALLOC_DH((1+m)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    cval = A->cval = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    aval = A->aval = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    rp[0] = 0;
    for (i = 0; i < m; ++ i)
    {
        jxf_bool isMissing = jxf_true;
        
        for (j = RP[i]; j < RP[i+1]; ++ j)
        {
            cval[idx] = CVAL[j];
            aval[idx] = AVAL[j];
            ++ idx;
            if (CVAL[j] == i) isMissing = jxf_false;
        }
        if (isMissing)
        {
            cval[idx] = i;
            aval[idx] = 0.0;
            ++ idx;
        }
        rp[i+1] = idx;
    }
    JXF_FREE_DH(RP); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(CVAL); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(AVAL); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_readVec"
void jxf_readVec( jxf_Vec_dh *bout, char *ft, char *fn, JXF_Int ignore )
{
    JXF_START_FUNC_DH
   *bout = NULL;
    if (fn == NULL)
    {
        JXF_SET_V_ERROR("passed NULL filename; can't open for reading!");
    }
    if (!strcmp(ft, "csr")  ||  !strcmp(ft, "trip"))
    {
        jxf_Vec_dhRead(bout, ignore, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "ebin"))
    {
        jxf_Vec_dhReadBIN(bout, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "petsc"))
    {
        jxf_sprintf(jxf_msgBuf_dh, "must recompile Euclid using petsc mode!");
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "unknown filetype: -ftin %s", ft);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_writeMat"
void jxf_writeMat( jxf_Mat_dh Ain, char *ft, char *fn )
{
    JXF_START_FUNC_DH
    if (fn == NULL)
    {
        JXF_SET_V_ERROR("passed NULL filename; can't open for writing!");
    }
    if (!strcmp(ft, "csr"))
    {
        jxf_Mat_dhPrintCSR(Ain, NULL, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "trip"))
    {
        jxf_Mat_dhPrintTriples(Ain, NULL, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "ebin"))
    {
        jxf_Mat_dhPrintBIN(Ain, NULL, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "petsc"))
    {
        jxf_sprintf(jxf_msgBuf_dh, "must recompile Euclid using petsc mode!");
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "unknown filetype: -ftout %s", ft);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_writeVec"
void jxf_writeVec( jxf_Vec_dh bin, char *ft, char *fn )
{
    JXF_START_FUNC_DH
    if (fn == NULL)
    {
        JXF_SET_V_ERROR("passed NULL filename; can't open for writing!");
    }
    if (!strcmp(ft, "csr")  ||  !strcmp(ft, "trip"))
    {
        jxf_Vec_dhPrint(bin, NULL, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "ebin"))
    {
        jxf_Vec_dhPrintBIN(bin, NULL, fn); JXF_CHECK_V_ERROR;
    }
    else if (!strcmp(ft, "petsc"))
    {
        jxf_sprintf(jxf_msgBuf_dh, "must recompile Euclid using petsc mode!");
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "unknown filetype: -ftout %s", ft);
        JXF_SET_V_ERROR(jxf_msgBuf_dh);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_isTriangular"
JXF_Int jxf_isTriangular( JXF_Int m, JXF_Int *rp, JXF_Int *cval )
{
    JXF_START_FUNC_DH
    JXF_Int row, j;
    JXF_Int type;
    jxf_bool type_lower = jxf_false, type_upper = jxf_false;
    
    if (jxf_np_dh > 1)
    {
        JXF_SET_ERROR(-1, "only implemented for a single cpu");
    }
    for (row = 0; row < m; ++ row)
    {
        for (j = rp[row]; j < rp[row+1]; ++ j)
        {
            JXF_Int col = cval[j];
            
            if (col < row) type_lower = jxf_true;
            if (col > row) type_upper = jxf_true;
        }
        if (type_lower && type_upper) break;
    }
    if (type_lower && type_upper)
    {
        type = JXF_IS_FULL;
    }
    else if (type_lower)
    {
        type = JXF_IS_LOWER_TRI;
    }
    else
    {
        type = JXF_IS_UPPER_TRI;
    }
    JXF_END_FUNC_VAL(type)
}

static void jxf_mat_dh_transpose_reuse_private_private( jxf_bool allocateMem,
             JXF_Int m, JXF_Int *rpIN, JXF_Int *cvalIN, JXF_Real *avalIN, JXF_Int **rpOUT, JXF_Int **cvalOUT, JXF_Real **avalOUT );

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_transpose_reuse_private"
void jxf_mat_dh_transpose_reuse_private( JXF_Int m,
                                     JXF_Int *rpIN,
                                     JXF_Int *cvalIN,
                                     JXF_Real *avalIN,
                                     JXF_Int *rpOUT,
                                     JXF_Int *cvalOUT,
                                     JXF_Real *avalOUT )
{
    JXF_START_FUNC_DH
    jxf_mat_dh_transpose_reuse_private_private(jxf_false, m,
            rpIN, cvalIN, avalIN, &rpOUT, &cvalOUT, &avalOUT); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_transpose_private"
void jxf_mat_dh_transpose_private( JXF_Int m,
                               JXF_Int *RP,
                               JXF_Int **rpOUT,
                               JXF_Int *CVAL,
                               JXF_Int **cvalOUT,
                               JXF_Real *AVAL,
                               JXF_Real **avalOUT )
{
    JXF_START_FUNC_DH
    jxf_mat_dh_transpose_reuse_private_private(jxf_true, m, RP, CVAL, AVAL, rpOUT, cvalOUT, avalOUT); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_dh_transpose_private_private"
void jxf_mat_dh_transpose_reuse_private_private( jxf_bool allocateMem,
                                             JXF_Int m,
                                             JXF_Int *RP,
                                             JXF_Int *CVAL,
                                             JXF_Real *AVAL,
                                             JXF_Int **rpOUT,
                                             JXF_Int **cvalOUT,
                                             JXF_Real **avalOUT )
{
    JXF_START_FUNC_DH
    JXF_Int *rp, *cval, *tmp;
    JXF_Int i, j, nz = RP[m];
    JXF_Real *aval = NULL;
    
    if (allocateMem)
    {
        rp = *rpOUT = (JXF_Int *)JXF_MALLOC_DH((1+m)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        cval = *cvalOUT = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        if (avalOUT != NULL)
        {
            aval = *avalOUT = (JXF_Real*)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
        }
    }
    else
    {
        rp = *rpOUT;
        cval = *cvalOUT;
        if (avalOUT != NULL) aval = *avalOUT;
    }
    tmp = (JXF_Int *)JXF_MALLOC_DH((1+m)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i <= m; ++ i) tmp[i] = 0;
    for (i = 0; i < m; ++ i)
    {
        for (j = RP[i]; j < RP[i+1]; ++ j)
        {
            JXF_Int col = CVAL[j];
            
            tmp[col+1] += 1;
        }
    }
    for (i = 1; i <= m; ++ i) tmp[i] += tmp[i-1];
    memcpy(rp, tmp, (m+1)*sizeof(JXF_Int));
    if (avalOUT != NULL)
    {
        for (i = 0; i < m; ++ i)
        {
            for (j = RP[i]; j < RP[i+1]; ++ j)
            {
                JXF_Int col = CVAL[j];
                JXF_Int idx = tmp[col];
                
                cval[idx] = i;
                aval[idx] = AVAL[j];
                tmp[col] += 1;
            }
        }
    }
    else
    {
        for (i = 0; i < m; ++ i)
        {
            for (j = RP[i]; j < RP[i+1]; ++ j)
            {
                JXF_Int col = CVAL[j];
                JXF_Int idx = tmp[col];
                
                cval[idx] = i;
                tmp[col] += 1;
            }
        }
    }
    JXF_FREE_DH(tmp); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_jxf_find_owner"
JXF_Int jxf_mat_jxf_find_owner( JXF_Int *beg_rows, JXF_Int *end_rows, JXF_Int index )
{
    JXF_START_FUNC_DH
    JXF_Int pe, owner = -1;
    
    for (pe = 0; pe < jxf_np_dh; ++ pe)
    {
        if (index >= beg_rows[pe] && index < end_rows[pe])
        {
            owner = pe;
            break;
        }
    }
    if (owner == -1)
    {
        jxf_sprintf(jxf_msgBuf_dh, "failed to jxf_find owner for index= %i", index);
        JXF_SET_ERROR(-1, jxf_msgBuf_dh);
    }
    JXF_END_FUNC_VAL(owner)
}

#define JXF_AVAL_TAG 2
#define JXF_CVAL_TAG 3

void jxf_partition_and_distribute_private( jxf_Mat_dh A, jxf_Mat_dh *Bout );
void jxf_partition_and_distribute_metis_private( jxf_Mat_dh A, jxf_Mat_dh *Bout );

#undef __FUNC__
#define __FUNC__ "jxf_readMat_par"
void jxf_readMat_par( jxf_Mat_dh *Aout, char *fileType, char *fileName, JXF_Int ignore )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh A = NULL;
    
    if (jxf_myid_dh == 0)
    {
        JXF_Int tmp = jxf_np_dh;
        
        jxf_np_dh = 1;
        jxf_readMat(&A, fileType, fileName, ignore); JXF_CHECK_V_ERROR;
        jxf_np_dh = tmp;
    }
    if (jxf_np_dh == 1)
    {
       *Aout = A;
    }
    else
    {
        if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-metis"))
        {
            jxf_partition_and_distribute_metis_private(A, Aout); JXF_CHECK_V_ERROR;
        }
        else
        {
            jxf_partition_and_distribute_private(A, Aout); JXF_CHECK_V_ERROR;
        }
    }
    if (jxf_np_dh > 1 && A != NULL)
    {
        jxf_Mat_dhDestroy(A); JXF_CHECK_V_ERROR;
    }
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh, "-printMAT"))
    {
        char xname[] = "A", *name = xname;
        
        jxf_Parser_dhReadString(jxf_parser_dh, "-printMat", &name);
        jxf_Mat_dhPrintTriples(*Aout, NULL, name); JXF_CHECK_V_ERROR;
        jxf_printf_dh("\n@@@ jxf_readMat_par: printed mat to %s\n\n", xname);
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_partition_and_distribute_metis_private"
void jxf_partition_and_distribute_metis_private( jxf_Mat_dh A, jxf_Mat_dh *Bout )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh B = NULL;
    jxf_Mat_dh C = NULL;
    JXF_Int i, m;
    JXF_Int *rowLengths = NULL;
    JXF_Int *o2n_row = NULL, *n2o_col = NULL, *rowToBlock = NULL;
    JXF_Int *beg_row = NULL, *row_count = NULL;
    MPI_Request *send_req = NULL;
    MPI_Request *rcv_req = NULL;
    MPI_Status *send_status = NULL;
    MPI_Status *rcv_status = NULL;
    
    jxf_MPI_Barrier(jxf_comm_dh);
    jxf_printf_dh("@@@ partitioning with metis\n");
    if (jxf_myid_dh == 0) m = A->m;
    jxf_MPI_Bcast(&m, 1, JXF_MPI_INT, 0, MPI_COMM_WORLD);
    rowLengths = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rowToBlock = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        JXF_Int *tmp = A->rp;
        
        for (i = 0; i < m; ++ i)
        {
            rowLengths[i] = tmp[i+1] - tmp[i];
        }
    }
    jxf_MPI_Bcast(rowLengths, m, JXF_MPI_INT, 0, jxf_comm_dh);
    if (jxf_myid_dh == 0)
    {
        JXF_Int idx = 0;
        JXF_Int j;
        
        jxf_Mat_dhPartition(A, jxf_np_dh, &beg_row, &row_count, &n2o_col, &o2n_row); JXF_ERRCHKA;
        jxf_Mat_dhPermute(A, n2o_col, &C); JXF_ERRCHKA;
        for (i = 0; i < jxf_np_dh; ++ i)
        {
            for (j = beg_row[i]; j < beg_row[i]+row_count[i]; ++ j)
            {
                rowToBlock[idx++] = i;
            }
        }
    }
    jxf_MPI_Bcast(rowToBlock, m, JXF_MPI_INT, 0, jxf_comm_dh);
    jxf_mat_par_read_allocate_private(&B, m, rowLengths, rowToBlock); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        JXF_Int *cval = C->cval, *rp = C->rp;
        JXF_Real *aval = C->aval;
        
        send_req = (MPI_Request *)JXF_MALLOC_DH(2*m*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
        send_status = (MPI_Status*)JXF_MALLOC_DH(2*m*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
        for (i = 0; i < m; ++ i)
        {
            JXF_Int owner = rowToBlock[i];
            JXF_Int count = rp[i+1] - rp[i];
            
            if (!count)
            {
                jxf_sprintf(jxf_msgBuf_dh, "row %i of %i is empty!", i+1, m);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
            jxf_MPI_Isend(cval+rp[i], count, JXF_MPI_INT, owner, JXF_CVAL_TAG, jxf_comm_dh, send_req+2*i);
            jxf_MPI_Isend(aval+rp[i], count, JXF_MPI_REAL, owner, JXF_AVAL_TAG, jxf_comm_dh, send_req+2*i+1);
        }
    }
    JXF_Int *cval = B->cval;
    JXF_Int *rp = B->rp;
    JXF_Real *aval = B->aval;
    m = B->m;
    rcv_req = (MPI_Request *)JXF_MALLOC_DH(2*m*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    rcv_status = (MPI_Status *)JXF_MALLOC_DH(2*m*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        JXF_Int count = rp[i+1] - rp[i];
        
        if (!count)
        {
            jxf_sprintf(jxf_msgBuf_dh, "local row %i of %i is empty!", i+1, m);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
        jxf_MPI_Irecv(cval+rp[i], count, JXF_MPI_INT, 0, JXF_CVAL_TAG, jxf_comm_dh, rcv_req+2*i);
        jxf_MPI_Irecv(aval+rp[i], count, JXF_MPI_REAL, 0, JXF_AVAL_TAG, jxf_comm_dh, rcv_req+2*i+1);
    }
    if (jxf_myid_dh == 0)
    {
        jxf_MPI_Waitall(m*2, send_req, send_status);
    }
    jxf_MPI_Waitall(2*B->m, rcv_req, rcv_status);
    if (rowLengths != NULL)
    {
        JXF_FREE_DH(rowLengths); JXF_CHECK_V_ERROR;
    }
    if (o2n_row != NULL)
    {
        JXF_FREE_DH(o2n_row); JXF_CHECK_V_ERROR;
    }
    if (n2o_col != NULL)
    {
        JXF_FREE_DH(n2o_col); JXF_CHECK_V_ERROR;
    }
    if (rowToBlock != NULL)
    {
        JXF_FREE_DH(rowToBlock); JXF_CHECK_V_ERROR;
    }
    if (send_req != NULL)
    {
        JXF_FREE_DH(send_req); JXF_CHECK_V_ERROR;
    }
    if (rcv_req != NULL)
    {
        JXF_FREE_DH(rcv_req); JXF_CHECK_V_ERROR;
    }
    if (send_status != NULL)
    {
        JXF_FREE_DH(send_status); JXF_CHECK_V_ERROR;
    }
    if (rcv_status != NULL)
    {
        JXF_FREE_DH(rcv_status); JXF_CHECK_V_ERROR;
    }
    if (beg_row != NULL)
    {
        JXF_FREE_DH(beg_row); JXF_CHECK_V_ERROR;
    }
    if (row_count != NULL)
    {
        JXF_FREE_DH(row_count); JXF_CHECK_V_ERROR;
    }
    if (C != NULL)
    {
        jxf_Mat_dhDestroy(C); JXF_ERRCHKA;
    }
   *Bout = B;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_partition_and_distribute_private"
void jxf_partition_and_distribute_private( jxf_Mat_dh A, jxf_Mat_dh *Bout )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh B = NULL;
    JXF_Int i, m;
    JXF_Int *rowLengths = NULL;
    JXF_Int *o2n_row = NULL, *n2o_col = NULL, *rowToBlock = NULL;
    MPI_Request *send_req = NULL;
    MPI_Request *rcv_req = NULL;
    MPI_Status *send_status = NULL;
    MPI_Status *rcv_status = NULL;
    
    jxf_MPI_Barrier(jxf_comm_dh);
    if (jxf_myid_dh == 0) m = A->m;
    jxf_MPI_Bcast(&m, 1, JXF_MPI_INT, 0, MPI_COMM_WORLD);
    rowLengths = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        JXF_Int *tmp = A->rp;
        for (i = 0; i < m; ++ i)
        {
            rowLengths[i] = tmp[i+1] - tmp[i];
        }
    }
    jxf_MPI_Bcast(rowLengths, m, JXF_MPI_INT, 0, jxf_comm_dh);
    rowToBlock = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        o2n_row = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        jxf_mat_partition_private(A, jxf_np_dh, o2n_row, rowToBlock); JXF_CHECK_V_ERROR;
    }
    jxf_MPI_Bcast(rowToBlock, m, JXF_MPI_INT, 0, jxf_comm_dh);
    jxf_mat_par_read_allocate_private(&B, m, rowLengths, rowToBlock); JXF_CHECK_V_ERROR;
    if (jxf_myid_dh == 0)
    {
        JXF_Int *cval = A->cval, *rp = A->rp;
        JXF_Real *aval = A->aval;
        
        send_req = (MPI_Request *)JXF_MALLOC_DH(2*m*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
        send_status = (MPI_Status *)JXF_MALLOC_DH(2*m*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
        for (i = 0; i < m; ++ i)
        {
            JXF_Int owner = rowToBlock[i];
            JXF_Int count = rp[i+1]-rp[i];
            
            if (!count)
            {
                jxf_sprintf(jxf_msgBuf_dh, "row %i of %i is empty!", i+1, m);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
            jxf_MPI_Isend(cval+rp[i], count, JXF_MPI_INT, owner, JXF_CVAL_TAG, jxf_comm_dh, send_req+2*i);
            jxf_MPI_Isend(aval+rp[i], count, JXF_MPI_REAL, owner, JXF_AVAL_TAG, jxf_comm_dh, send_req+2*i+1);
        }
    }
    JXF_Int *cval = B->cval;
    JXF_Int *rp = B->rp;
    JXF_Real *aval = B->aval;
    m = B->m;
    rcv_req = (MPI_Request *)JXF_MALLOC_DH(2*m*sizeof(MPI_Request)); JXF_CHECK_V_ERROR;
    rcv_status = (MPI_Status *)JXF_MALLOC_DH(2*m*sizeof(MPI_Status)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        JXF_Int count = rp[i+1] - rp[i];
        
        if (!count)
        {
            jxf_sprintf(jxf_msgBuf_dh, "local row %i of %i is empty!", i+1, m);
            JXF_SET_V_ERROR(jxf_msgBuf_dh);
        }
        jxf_MPI_Irecv(cval+rp[i], count, JXF_MPI_INT, 0, JXF_CVAL_TAG, jxf_comm_dh, rcv_req+2*i);
        jxf_MPI_Irecv(aval+rp[i], count, JXF_MPI_REAL, 0, JXF_AVAL_TAG, jxf_comm_dh, rcv_req+2*i+1);
    }
    if (jxf_myid_dh == 0)
    {
        jxf_MPI_Waitall(m*2, send_req, send_status);
    }
    jxf_MPI_Waitall(2*B->m, rcv_req, rcv_status);
    if (rowLengths != NULL)
    {
        JXF_FREE_DH(rowLengths); JXF_CHECK_V_ERROR;
    }
    if (o2n_row != NULL)
    {
        JXF_FREE_DH(o2n_row); JXF_CHECK_V_ERROR;
    }
    if (n2o_col != NULL)
    {
        JXF_FREE_DH(n2o_col); JXF_CHECK_V_ERROR;
    }
    if (rowToBlock != NULL)
    {
        JXF_FREE_DH(rowToBlock); JXF_CHECK_V_ERROR;
    }
    if (send_req != NULL)
    {
        JXF_FREE_DH(send_req); JXF_CHECK_V_ERROR;
    }
    if (rcv_req != NULL)
    {
        JXF_FREE_DH(rcv_req); JXF_CHECK_V_ERROR;
    }
    if (send_status != NULL)
    {
        JXF_FREE_DH(send_status); JXF_CHECK_V_ERROR;
    }
    if (rcv_status != NULL)
    {
        JXF_FREE_DH(rcv_status); JXF_CHECK_V_ERROR;
    }
   *Bout = B;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_par_read_allocate_private"
void jxf_mat_par_read_allocate_private( jxf_Mat_dh *Aout, JXF_Int n, JXF_Int *rowLengths, JXF_Int *rowToBlock )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh A;
    JXF_Int i, m, nz, beg_row, *rp, idx;
    
    jxf_Mat_dhCreate(&A); JXF_CHECK_V_ERROR;
   *Aout = A;
    A->n = n;
    m = 0;
    for (i = 0; i < n; ++ i)
    {
        if (rowToBlock[i] == jxf_myid_dh) ++ m;
    }
    A->m = m;
    beg_row = 0;
    for (i = 0; i < n; ++ i)
    {
        if (rowToBlock[i] < jxf_myid_dh) ++ beg_row;
    }
    A->beg_row = beg_row;
    A->rp = rp = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    rp[0] = 0;
    nz = 0;
    idx = 1;
    for (i = 0; i < n; ++ i)
    {
        if (rowToBlock[i] == jxf_myid_dh)
        {
            nz += rowLengths[i];
            rp[idx++] = nz;
        }
    }
    A->cval = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    A->aval = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_mat_partition_private"
void jxf_mat_partition_private( jxf_Mat_dh A, JXF_Int blocks, JXF_Int *o2n_row, JXF_Int *rowToBlock )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, n = A->n;
    JXF_Int rpb = n / blocks;
    JXF_Int idx = 0;
    
    while (rpb*blocks < n) ++ rpb;
    if (rpb*(blocks-1) == n)
    {
        -- rpb;
        jxf_printf_dh("adjusted rpb to: %i\n", rpb);
    }
    for (i = 0; i < n; ++ i) o2n_row[i] = i;
    for (i = 0; i < blocks-1; ++ i)
    {
        for (j = 0; j < rpb; ++ j)
        {
            rowToBlock[idx++] = i;
        }
    }
    i = blocks - 1;
    while (idx < n) rowToBlock[idx++] = i;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_make_full_private"
void jxf_make_full_private( JXF_Int m, JXF_Int **rpIN, JXF_Int **cvalIN, JXF_Real **avalIN )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, *rpNew, *cvalNew, *rp = *rpIN, *cval = *cvalIN;
    JXF_Real *avalNew, *aval = *avalIN;
    JXF_Int nz, *rowCounts = NULL;
    
    rowCounts = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i <= m; ++ i) rowCounts[i] = 0;
    for (i = 0; i < m; ++ i)
    {
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int col = cval[j];
            
            rowCounts[i+1] += 1;
            if (col != i) rowCounts[col+1] += 1;
        }
    }
    rpNew = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 1; i <= m; ++ i) rowCounts[i] += rowCounts[i-1];
    memcpy(rpNew, rowCounts, (m+1)*sizeof(JXF_Int));
    nz = rpNew[m];
    cvalNew = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    avalNew = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i)
    {
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int col = cval[j];
            JXF_Real val  = aval[j];
            
            cvalNew[rowCounts[i]] = col;
            avalNew[rowCounts[i]] = val;
            rowCounts[i] += 1;
            if (col != i)
            {
                cvalNew[rowCounts[col]] = i;
                avalNew[rowCounts[col]] = val;
                rowCounts[col] += 1;
            }
        }
    }
    if (rowCounts != NULL)
    {
        JXF_FREE_DH(rowCounts); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(cval); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(rp); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(aval); JXF_CHECK_V_ERROR;
   *rpIN = rpNew;
   *cvalIN = cvalNew;
   *avalIN = avalNew;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_make_symmetric_private"
void jxf_make_symmetric_private( JXF_Int m, JXF_Int **rpIN, JXF_Int **cvalIN, JXF_Real **avalIN )
{
    JXF_START_FUNC_DH
    JXF_Int i, j, *rpNew, *cvalNew, *rp = *rpIN, *cval = *cvalIN;
    JXF_Real *avalNew, *aval = *avalIN;
    JXF_Int nz, *rowCounts = NULL;
    JXF_Int *rpTrans, *cvalTrans;
    JXF_Int *work;
    JXF_Real *avalTrans;
    JXF_Int nzCount = 0, transCount = 0;
    
    jxf_mat_dh_transpose_private(m, rp, &rpTrans, cval, &cvalTrans, aval, &avalTrans); JXF_CHECK_V_ERROR;
    work = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) work[i] = -1;
    rowCounts = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 0; i <= m; ++ i) rowCounts[i] = 0;
    for (i = 0; i < m; ++ i)
    {
        JXF_Int ct = 0;
        
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int col = cval[j];
            
            work[col] = i;
            ++ ct;
            ++ nzCount;
        }
        for (j = rpTrans[i]; j < rpTrans[i+1]; ++ j)
        {
            JXF_Int col = cvalTrans[j];
            
            if (work[col] != i)
            {
                ++ ct;
                ++ transCount;
            }
        }
        rowCounts[i+1] = ct;
    }
    if (transCount == 0)
    {
        jxf_printf("jxf_make_symmetric_private: matrix is already structurally symmetric!\n");
        JXF_FREE_DH(rpTrans); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(cvalTrans); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(avalTrans); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(work); JXF_CHECK_V_ERROR;
        JXF_FREE_DH(rowCounts); JXF_CHECK_V_ERROR;
        goto END_OF_FUNCTION;
    }
    else
    {
        jxf_printf("original nz= %i\n", rp[m]);
        jxf_printf("zeros added= %i\n", transCount);
        jxf_printf("ratio of added zeros to nonzeros = %0.2f (assumes all original entries were nonzero!)\n",
                                                                       (JXF_Real)transCount/(JXF_Real)(nzCount));
    }
    rpNew = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    for (i = 1; i <= m; ++ i) rowCounts[i] += rowCounts[i-1];
    memcpy(rpNew, rowCounts, (m+1)*sizeof(JXF_Int));
    for (i = 0; i < m; ++ i) work[i] = -1;
    nz = rpNew[m];
    cvalNew = (JXF_Int *)JXF_MALLOC_DH(nz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    avalNew = (JXF_Real *)JXF_MALLOC_DH(nz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) work[i] = -1;
    for (i = 0; i < m; ++ i)
    {
        for (j = rp[i]; j < rp[i+1]; ++ j)
        {
            JXF_Int col = cval[j];
            JXF_Real val  = aval[j];
            
            work[col] = i;
            cvalNew[rowCounts[i]] = col;
            avalNew[rowCounts[i]] = val;
            rowCounts[i] += 1;
        }
        for (j = rpTrans[i]; j < rpTrans[i+1]; ++ j)
        {
            JXF_Int col = cvalTrans[j];
            
            if (work[col] != i)
            {
                cvalNew[rowCounts[i]] = col;
                avalNew[rowCounts[i]] = 0.0;
                rowCounts[i] += 1;
            }
        }
    }
    if (rowCounts != NULL)
    {
        JXF_FREE_DH(rowCounts); JXF_CHECK_V_ERROR;
    }
    JXF_FREE_DH(work); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(cval); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(rp); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(aval); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(cvalTrans); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(rpTrans); JXF_CHECK_V_ERROR;
    JXF_FREE_DH(avalTrans); JXF_CHECK_V_ERROR;
   *rpIN = rpNew;
   *cvalIN = cvalNew;
   *avalIN = avalNew;
    
END_OF_FUNCTION: ;
    
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_profileMat"
void jxf_profileMat( jxf_Mat_dh A )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh B = NULL;
    JXF_Int type;
    JXF_Int m;
    JXF_Int i, j;
    JXF_Int *work1 = NULL;
    JXF_Real *work2 = NULL;
    jxf_bool isStructurallySymmetric = jxf_true;
    jxf_bool isNumericallySymmetric = jxf_true;
    jxf_bool is_Triangular = jxf_false;
    JXF_Int zeroCount = 0, nz;
    
    if (jxf_myid_dh > 0)
    {
        JXF_SET_V_ERROR("only for a single jxf_MPI task!");
    }
    m = A->m;
    jxf_printf("\nYY----------------------------------------------------\n");
    nz = A->rp[m];
    for (i = 0; i < nz; ++ i)
    {
        if (A->aval[i] == 0) ++ zeroCount;
    }
    jxf_printf("YY  row count:      %i\n", m);
    jxf_printf("YY  nz count:       %i\n", nz);
    jxf_printf("YY  explicit zeros: %i (entire matrix)\n", zeroCount);
    JXF_Int m_diag = 0, z_diag = 0;
    for (i = 0; i < m; ++ i)
    {
        jxf_bool flag = jxf_true;
        
        for (j = A->rp[i]; j < A->rp[i+1]; ++ j)
        {
            JXF_Int col = A->cval[j];
            
            if (col == i)
            {
                JXF_Real val = A->aval[j];
                
                flag = jxf_false;
                if (val == 0.0) ++ z_diag;
                break;
            }
        }
        if (flag) ++ m_diag;
    }
    jxf_printf("YY  missing diagonals:   %i\n", m_diag);
    jxf_printf("YY  explicit zero diags: %i\n", z_diag);
    type = jxf_isTriangular(m, A->rp, A->cval); JXF_CHECK_V_ERROR;
    if (type == JXF_IS_UPPER_TRI)
    {
        jxf_printf("YY  matrix is upper triangular\n");
        is_Triangular = jxf_true;
        goto END_OF_FUNCTION;
    }
    else if (type == JXF_IS_LOWER_TRI)
    {
        jxf_printf("YY  matrix is lower triangular\n");
        is_Triangular = jxf_true;
        goto END_OF_FUNCTION;
    }
    JXF_Int unz = 0, lnz = 0;
    for (i = 0; i < m; ++ i)
    {
        for (j = A->rp[i]; j < A->rp[i+1]; ++ j)
        {
            JXF_Int col = A->cval[j];
            
            if (col < i) ++ lnz;
            if (col > i) ++ unz;
        }
    }
    jxf_printf("YY  strict upper triangular nonzeros: %i\n", unz);
    jxf_printf("YY  strict lower triangular nonzeros: %i\n", lnz);
    jxf_Mat_dhTranspose(A, &B); JXF_CHECK_V_ERROR;
    work1 = (JXF_Int *)JXF_MALLOC_DH(m*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
    work2 = (JXF_Real *)JXF_MALLOC_DH(m*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    for (i = 0; i < m; ++ i) work1[i] = -1;
    for (i = 0; i < m; ++ i) work2[i] = 0.0;
    for (i = 0; i < m; ++ i)
    {
        for (j = A->rp[i]; j < A->rp[i+1]; ++ j)
        {
            JXF_Int col = A->cval[j];
            JXF_Real val = A->aval[j];
            
            work1[col] = i;
            work2[col] = val;
        }
        for (j = B->rp[i]; j < B->rp[i+1]; ++ j)
        {
            JXF_Int col = B->cval[j];
            JXF_Real val = B->aval[j];
            
            if (work1[col] != i)
            {
                isStructurallySymmetric = jxf_false;
                isNumericallySymmetric = jxf_false;
                goto END_OF_FUNCTION;
            }
            if (work2[col] != val)
            {
                isNumericallySymmetric = jxf_false;
                work2[col] = 0.0;
            }
        }
    }
    
END_OF_FUNCTION: ;
    
    if (!is_Triangular)
    {
        jxf_printf("YY  matrix is NOT triangular\n");
        if (isStructurallySymmetric)
        {
            jxf_printf("YY  matrix IS structurally symmetric\n");
        }
        else
        {
            jxf_printf("YY  matrix is NOT structurally symmetric\n");
        }
        if (isNumericallySymmetric)
        {
            jxf_printf("YY  matrix IS numerically symmetric\n");
        }
        else
        {
            jxf_printf("YY  matrix is NOT numerically symmetric\n");
        }
    }
    if (work1 != NULL)
    {
        JXF_FREE_DH(work1); JXF_CHECK_V_ERROR;
    }
    if (work2 != NULL)
    {
        JXF_FREE_DH(work2); JXF_CHECK_V_ERROR;
    }
    if (B != NULL)
    {
        jxf_Mat_dhDestroy(B); JXF_CHECK_V_ERROR;
    }
    jxf_printf("YY----------------------------------------------------\n");
    JXF_END_FUNC_DH
}
