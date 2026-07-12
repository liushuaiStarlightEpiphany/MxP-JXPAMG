//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  mat_gen_fd_dh.c
 *  Date: 2013/01/23
 */

#include "jxf_euclid.h"

static jxf_bool jxf_isThreeD;

#define JXF_FRONT(a) a[5]
#define JXF_SOUTH(a) a[3]
#define JXF_WEST(a) a[1]
#define JXF_CENTER(a) a[0]
#define JXF_EAST(a) a[2]
#define JXF_NORTH(a) a[4]
#define JXF_BACK(a) a[6]
#define JXF_RHS(a) a[7]

static void jxf_setBoundary_private( JXF_Int node, JXF_Int *cval, JXF_Real *aval,
                 JXF_Int len, JXF_Real *rhs, JXF_Real bc, JXF_Real coeff, JXF_Real ctr, JXF_Int nabor );
static void jxf_generateStriped( jxf_MatGenFD mg, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, jxf_Mat_dh A, jxf_Vec_dh b );
static void jxf_generateBlocked( jxf_MatGenFD mg, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, jxf_Mat_dh A, jxf_Vec_dh b );
static void jxf_getstencil( jxf_MatGenFD g, JXF_Int ix, JXF_Int iy, JXF_Int iz );

#undef __FUNC__
#define __FUNC__ "jxf_MatGenFDCreate"
void jxf_MatGenFD_Create( jxf_MatGenFD *mg )
{
    JXF_START_FUNC_DH
    struct _jxf_matgenfd *tmp =(struct _jxf_matgenfd *)JXF_MALLOC_DH(sizeof(struct _jxf_matgenfd)); JXF_CHECK_V_ERROR;
    
   *mg = tmp;
    tmp->debug = jxf_Parser_dhHasSwitch(jxf_parser_dh, "-debug_matgen");
    tmp->m = 9;
    tmp->px = tmp->py = 1;
    tmp->pz = 0;
    jxf_Parser_dhReadInt(jxf_parser_dh,"-m",&tmp->m);
    jxf_Parser_dhReadInt(jxf_parser_dh,"-px",&tmp->px);
    jxf_Parser_dhReadInt(jxf_parser_dh,"-py",&tmp->py);
    jxf_Parser_dhReadInt(jxf_parser_dh,"-pz",&tmp->pz);
    if (tmp->px < 1) tmp->px = 1;
    if (tmp->py < 1) tmp->py = 1;
    if (tmp->pz < 0) tmp->pz = 0;
    tmp->threeD = jxf_false;
    if (tmp->pz)
    {
        tmp->threeD = jxf_true;
    }
    else
    {
        tmp->pz = 1;
    }
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh,"-threeD")) tmp->threeD = jxf_true;
    tmp->a = tmp->b = tmp->c = 1.0;
    tmp->d = tmp->e = tmp->f = 0.0;
    tmp->g = tmp->h = 0.0;
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-dx",&tmp->a);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-dy",&tmp->b);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-dz",&tmp->c);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-cx",&tmp->d);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-cy",&tmp->e);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-cz",&tmp->f);
    tmp->a = -1 * fabs(tmp->a);
    tmp->b = -1 * fabs(tmp->b);
    tmp->c = -1 * fabs(tmp->c);
    tmp->allocateMem = jxf_true;
    tmp->A = tmp->B = tmp->C = tmp->D = tmp->E = tmp->F = tmp->G = tmp->H = jxf_konstant;
    tmp->bcX1 = tmp->bcX2 = tmp->bcY1 = tmp->bcY2 = tmp->bcZ1 = tmp->bcZ2 = 0.0;
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcx1",&tmp->bcX1);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcx2",&tmp->bcX2);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcy1",&tmp->bcY1);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcy2",&tmp->bcY2);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcz1",&tmp->bcZ1);
    jxf_Parser_dhReadDouble(jxf_parser_dh,"-bcz2",&tmp->bcZ2);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_MatGenFD_Destroy"
void jxf_MatGenFD_Destroy( jxf_MatGenFD mg )
{
    JXF_START_FUNC_DH
    JXF_FREE_DH(mg); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_MatGenFD_Run"
void jxf_MatGenFD_Run( jxf_MatGenFD mg, JXF_Int id, JXF_Int np, jxf_Mat_dh *AOut, jxf_Vec_dh *rhsOut )
{
    JXF_START_FUNC_DH
    jxf_Mat_dh A;
    jxf_Vec_dh rhs;
    jxf_bool threeD = mg->threeD;
    JXF_Int nnz;
    JXF_Int m = mg->m;
    jxf_bool debug = jxf_false, striped;
    
    if (mg->debug && jxf_logFile != NULL) debug = jxf_true;
    striped = jxf_Parser_dhHasSwitch(jxf_parser_dh,"-striped");
    jxf_Mat_dhCreate(AOut); JXF_CHECK_V_ERROR;
    jxf_Vec_dhCreate(rhsOut); JXF_CHECK_V_ERROR;
    A = *AOut;
    rhs = *rhsOut;
    if (!jxf_Parser_dhHasSwitch(jxf_parser_dh, "-noChecks"))
    {
        if (!striped)
        {
            JXF_Int npTest = mg->px * mg->py;
            
            if (threeD) npTest *= mg->pz;
            if (npTest != np)
            {
                jxf_sprintf(jxf_msgBuf_dh, "numbers don't match: jxf_np_dh = %i, px*py*pz = %i", np, npTest);
                JXF_SET_V_ERROR(jxf_msgBuf_dh);
            }
        }
    }
    mg->cc = m;
    if (threeD)
    {
        m = mg->m = m * m * m;
    }
    else
    {
        m = mg->m = m * m;
    }
    mg->first = id * m;
    mg->hh = 1.0 / (mg->px * mg->cc - 1);
    if (debug)
    {
        jxf_sprintf(jxf_msgBuf_dh, "cc (local grid dimension) = %i", mg->cc);
        JXF_SET_INFO(jxf_msgBuf_dh);
        if (threeD)
        {
            jxf_sprintf(jxf_msgBuf_dh, "threeD = jxf_true");
        }
        else
        {
            jxf_sprintf(jxf_msgBuf_dh, "threeD = jxf_false");
        }
        JXF_SET_INFO(jxf_msgBuf_dh);
        jxf_sprintf(jxf_msgBuf_dh, "np= %i  id= %i", np, id);
        JXF_SET_INFO(jxf_msgBuf_dh);
    }
    mg->id = id;
    mg->np = np;
    nnz = threeD ? m*7 : m*5;
    if (mg->allocateMem)
    {
        A->rp = (JXF_Int *)JXF_MALLOC_DH((m+1)*sizeof(JXF_Int)); JXF_CHECK_V_ERROR;
        A->rp[0] = 0;
        A->cval = (JXF_Int *)JXF_MALLOC_DH(nnz*sizeof(JXF_Int)); JXF_CHECK_V_ERROR
        A->aval = (JXF_Real *)JXF_MALLOC_DH(nnz*sizeof(JXF_Real)); JXF_CHECK_V_ERROR;
    }
    rhs->n = m;
    A->m = m;
    A->n = m * mg->np;
    A->beg_row = mg->first;
    jxf_isThreeD = threeD;
    if (jxf_Parser_dhHasSwitch(jxf_parser_dh,"-striped"))
    {
        jxf_generateStriped(mg, A->rp, A->cval, A->aval, A, rhs); JXF_CHECK_V_ERROR;
    }
    else
    {
        jxf_generateBlocked(mg, A->rp, A->cval, A->aval, A, rhs); JXF_CHECK_V_ERROR;
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_generateStriped"
void jxf_generateStriped( jxf_MatGenFD mg, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, jxf_Mat_dh A, jxf_Vec_dh b )
{
    JXF_START_FUNC_DH
    JXF_Int mGlobal;
    JXF_Int m = mg->m;
    JXF_Int beg_row, end_row;
    JXF_Int i, j, k, row;
    jxf_bool threeD = mg->threeD;
    JXF_Int idx = 0;
    JXF_Real *stencil = mg->stencil;
    jxf_bool debug = jxf_false;
    JXF_Int plane, nodeRemainder;
    JXF_Int naborx1=0, naborx2=0, nabory1=0, nabory2=0;
    JXF_Real *rhs;
    jxf_bool applyBdry = jxf_true;
    JXF_Real hhalf;
    JXF_Real bcx1 = mg->bcX1;
    JXF_Real bcx2 = mg->bcX2;
    JXF_Real bcy1 = mg->bcY1;
    JXF_Real bcy2 = mg->bcY2;
    JXF_Int nx, ny;
    
    jxf_printf_dh("@@@ using striped partitioning\n");
    if (mg->debug && jxf_logFile != NULL) debug = jxf_true;
    m = 9;
    jxf_Parser_dhReadInt(jxf_parser_dh,"-m", &m);
    mGlobal = m * m;
    if (threeD) mGlobal *= m;
    i = mGlobal / mg->np;
    beg_row = i * mg->id;
    end_row = beg_row + i;
    if (mg->id == mg->np-1) end_row = mGlobal;
    nx = ny = m;
    mg->hh = 1.0 / (m - 1);
    hhalf = 0.5 * mg->hh;
    A->n = m * m;
    A->m = end_row - beg_row;
    A->beg_row = beg_row;
    jxf_Vec_dhInit(b, A->m); JXF_CHECK_V_ERROR;
    rhs = b->vals;
    plane = m * m;
    if (debug)
    {
        jxf_fprintf(jxf_logFile, "jxf_generateStriped: beg_row= %i; end_row= %i; m= %i\n", beg_row+1, end_row+1, m);
    }
    for (row = beg_row; row < end_row; ++ row)
    {
        JXF_Int localRow = row - beg_row;
        
        k = (row / plane);
        nodeRemainder = row - (k * plane);
        j = nodeRemainder / m;
        i = nodeRemainder % m;
        if (debug)
        {
            jxf_fprintf(jxf_logFile, "row= %i  x= %i  y= %i  z= %i\n", row+1, i,j,k);
        }
        jxf_getstencil(mg,i,j,k);
        if (threeD)
        {
            if (k > 0)
            {
                cval[idx] = row - plane;
                aval[idx++] = JXF_BACK(stencil);
            }
        }
        if (j > 0)
        {
            nabory1 = cval[idx] = row - m;
            aval[idx++] = JXF_SOUTH(stencil);
        }
        if (i > 0)
        {
            naborx1 = cval[idx] = row - 1;
            aval[idx++] = JXF_WEST(stencil);
        }
        cval[idx] = row;
        aval[idx++] = JXF_CENTER(stencil);
        if (i < m-1)
        {
            naborx2 = cval[idx] = row + 1;
            aval[idx++] = JXF_EAST(stencil);
        }
        if (j < m-1)
        {
            nabory2 = cval[idx] = row + m;
            aval[idx++] = JXF_NORTH(stencil);
        }
        if (threeD)
        {
            if (k < m-1)
            {
                cval[idx] = row + plane;
                aval[idx++] = JXF_FRONT(stencil);
            }
        }
        rhs[localRow] = 0.0;
        ++ localRow;
        rp[localRow] = idx;
        if (!threeD && applyBdry)
        {
            JXF_Int offset = rp[localRow-1];
            JXF_Int len = rp[localRow] - rp[localRow-1];
            JXF_Real ctr, coeff;
            
            if (i == 0)
            {
                coeff = mg->A(mg->a, i+hhalf, j, k);
                ctr = mg->A(mg->a, i-hhalf, j, k);
                jxf_setBoundary_private(row, cval+offset, aval+offset,
                            len, &(rhs[localRow-1]), bcx1, coeff, ctr, naborx2);
            }
            else if (i == nx-1)
            {
                coeff = mg->A(mg->a, i-hhalf, j, k);
                ctr = mg->A(mg->a, i+hhalf, j, k);
                jxf_setBoundary_private(row, cval+offset, aval+offset,
                            len, &(rhs[localRow-1]), bcx2, coeff, ctr, naborx1);
            }
            else if (j == 0)
            {
                coeff = mg->B(mg->b, i, j+hhalf, k);
                ctr = mg->B(mg->b, i, j-hhalf, k);
                jxf_setBoundary_private(row, cval+offset, aval+offset,
                            len, &(rhs[localRow-1]), bcy1, coeff, ctr, nabory2);
            }
            else if (j == ny-1)
            {
                coeff = mg->B(mg->b, i, j-hhalf, k);
                ctr = mg->B(mg->b, i, j+hhalf, k);
                jxf_setBoundary_private(row, cval+offset, aval+offset,
                            len, &(rhs[localRow-1]), bcy2, coeff, ctr, nabory1);
            }
        }
    }
    JXF_END_FUNC_DH
}

JXF_Int jxf_rownum( const jxf_bool threeD,
            const JXF_Int x,
            const JXF_Int y,
            const JXF_Int z,
            const JXF_Int nx,
            const JXF_Int ny,
            const JXF_Int nz,
            JXF_Int P,
            JXF_Int Q )
{
    JXF_Int p, q, r;
    JXF_Int lowerx, lowery, lowerz;
    JXF_Int id, startrow;
    
    p = x / nx;
    q = y / ny;
    r = z / nz;
    if (threeD)
    {
        id = r * P * Q + q * P + p;
    }
    else
    {
        id = q * P + p;
    }
    startrow = id * (nx * ny * nz);
    lowerx = nx * p;
    lowery = ny * q;
    lowerz = nz * r;
    if (threeD)
    {
        return (startrow + nx * ny * (z - lowerz) + nx * (y - lowery) + (x - lowerx));
    }
    else
    {
        return (startrow + nx * (y - lowery) + (x - lowerx));
    }
}

void jxf_getstencil( jxf_MatGenFD g, JXF_Int ix, JXF_Int iy, JXF_Int iz )
{
    JXF_Int k;
    JXF_Real h = g->hh;
    JXF_Real hhalf = h * 0.5;
    JXF_Real x = h * ix;
    JXF_Real y = h * iy;
    JXF_Real z = h * iz;
    JXF_Real cntr = 0.0;
    JXF_Real *stencil = g->stencil;
    JXF_Real coeff;
    jxf_bool threeD = g->threeD;
    
    for (k = 0; k < 8; ++ k) stencil[k] = 0.0;
    coeff = g->A(g->a, x+hhalf, y, z);
    JXF_EAST(stencil) += coeff;
    cntr += coeff;
    coeff = g->A(g->a, x-hhalf, y, z);
    JXF_WEST(stencil) += coeff;
    cntr += coeff;
    coeff = g->D(g->d, x, y, z) * hhalf;
    JXF_EAST(stencil) += coeff;
    JXF_WEST(stencil) -= coeff;
    coeff = g->B(g->b, x, y+hhalf, z);
    JXF_NORTH(stencil) += coeff;
    cntr += coeff;
    coeff = g->B(g->b, x, y-hhalf, z);
    JXF_SOUTH(stencil) += coeff;
    cntr += coeff;
    coeff = g->E(g->e, x, y, z) * hhalf;
    JXF_NORTH(stencil) += coeff;
    JXF_SOUTH(stencil) -= coeff;
    if (threeD)
    {
        coeff = g->C(g->c, x, y, z+hhalf);
        JXF_BACK(stencil) += coeff;
        cntr += coeff;
        coeff = g->C(g->c, x, y, z-hhalf);
        JXF_FRONT(stencil) += coeff;
        cntr += coeff;
        coeff = g->F(g->f, x, y, z) * hhalf;
        JXF_BACK(stencil) += coeff;
        JXF_FRONT(stencil) -= coeff;
    }
    coeff = g->G(g->g, x, y, z);
    JXF_CENTER(stencil) = h * h * coeff - cntr;
    JXF_RHS(stencil) = h * h * g->H(g->h, x, y, z);
}

JXF_Real jxf_konstant( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z )
{
    return coeff;
}

JXF_Real jxf_e2_xy( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z )
{
    return exp(coeff * x * y);
}

JXF_Real jxf_boxThreeD( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );

JXF_Real jxf_box_1( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z )
{
    static jxf_bool setup = jxf_false;
    JXF_Real retval = coeff;
    static JXF_Real dd1 = JXF_BOX1_DD;
    static JXF_Real dd2 = JXF_BOX2_DD;
    static JXF_Real dd3 = JXF_BOX3_DD;
    static JXF_Real ax1 = JXF_BOX1_X1, ay1 = JXF_BOX1_Y1;
    static JXF_Real ax2 = JXF_BOX1_X2, ay2 = JXF_BOX1_Y2;
    static JXF_Real bx1 = JXF_BOX2_X1, by1 = JXF_BOX2_Y1;
    static JXF_Real bx2 = JXF_BOX2_X2, by2 = JXF_BOX2_Y2;
    static JXF_Real cx1 = JXF_BOX3_X1, cy1 = JXF_BOX3_Y1;
    static JXF_Real cx2 = JXF_BOX3_X2, cy2 = JXF_BOX3_Y2;
    
    if (jxf_isThreeD)
    {
        return jxf_boxThreeD(coeff, x, y, z);
    }
    if (!setup)
    {
        dd1 = 0.1;
        dd2 = 0.1;
        dd3 = 10;
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-dd1", &dd1);
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-dd2", &dd2);
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-dd3", &dd3);
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-box1x1", &cx1);
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-box1x2", &cx2);
        setup = jxf_true;
    }
    if (x > ax1 && x < ax2 && y > ay1 && y < ay2)
    {
        retval = dd1 * coeff;
    }
    if (x > bx1 && x < bx2 && y > by1 && y < by2)
    {
        retval = dd2 * coeff;
    }
    if (x > cx1 && x < cx2 && y > cy1 && y < cy2)
    {
        retval = dd3 * coeff;
    }
    
    return retval;
}

JXF_Real jxf_boxThreeD( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z )
{
    static jxf_bool setup = jxf_false;
    JXF_Real retval = coeff;
    static JXF_Real dd1 = 100;
    static JXF_Real x1 = 0.2, x2 = 0.8;
    static JXF_Real y1 = 0.3, y2 = 0.7;
    static JXF_Real z1 = 0.4, z2 = 0.6;
    
    if (!setup)
    {
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-dd1", &dd1);
        setup = jxf_true;
    }
    if (x > x1 && x < x2 && y > y1 && y < y2 && z > z1 && z < z2)
    {
        retval = dd1 * coeff;
    }
    
    return retval;
}

JXF_Real jxf_box_2( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z )
{
    jxf_bool setup = jxf_false;
    static JXF_Real d1, d2;
    JXF_Real retval;
    
    if (!setup)
    {
        d1 = 1;
        d2 = 2;
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-bd1", &d1);
        jxf_Parser_dhReadDouble(jxf_parser_dh, "-bd2", &d2);
    }
    retval = d2;
    if (x < 0.5 && y < 0.5) retval = d1;
    if (x > 0.5 && y > 0.5) retval = d1;
    
    return (-1 * retval);
}

#undef __FUNC__
#define __FUNC__ "jxf_generateBlocked"
void jxf_generateBlocked( jxf_MatGenFD mg, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, jxf_Mat_dh A, jxf_Vec_dh b )
{
    JXF_START_FUNC_DH
    jxf_bool applyBdry = jxf_true;
    JXF_Real *stencil = mg->stencil;
    JXF_Int id = mg->id;
    jxf_bool threeD = mg->threeD;
    JXF_Int px = mg->px, py = mg->py, pz = mg->pz;
    JXF_Int p, q, r;
    JXF_Int cc = mg->cc;
    JXF_Int nx = cc, ny = cc, nz = cc;
    JXF_Int lowerx, upperx, lowery, uppery, lowerz, upperz;
    JXF_Int startRow;
    JXF_Int x, y, z;
    jxf_bool debug = jxf_false;
    JXF_Int idx = 0, localRow = 0;
    JXF_Int naborx1=0, naborx2=0, nabory1=0, nabory2=0, naborz1, naborz2;
    JXF_Real *rhs;
    JXF_Real hhalf = 0.5 * mg->hh;
    JXF_Real bcx1 = mg->bcX1;
    JXF_Real bcx2 = mg->bcX2;
    JXF_Real bcy1 = mg->bcY1;
    JXF_Real bcy2 = mg->bcY2;
    
    jxf_Vec_dhInit(b, A->m); JXF_CHECK_V_ERROR;
    rhs = b->vals;
    if (mg->debug && jxf_logFile != NULL) debug = jxf_true;
    if (!threeD) nz = 1;
    p = id % px;
    q = ((id - p) / px) % py;
    r = (id - p - px * q) / (px * py);
    if (debug)
    {
        jxf_sprintf(jxf_msgBuf_dh, "this proc's position in subdomain grid: p= %i  q= %i  r= %i", p, q, r);
        JXF_SET_INFO(jxf_msgBuf_dh);
    }
    lowerx = nx * p;
    upperx = lowerx + nx;
    lowery = ny * q;
    uppery = lowery + ny;
    lowerz = nz * r;
    upperz = lowerz + nz;
    if (debug)
    {
        jxf_sprintf(jxf_msgBuf_dh, "local grid parameters: lowerx= %i  upperx= %i", lowerx, upperx);
        JXF_SET_INFO(jxf_msgBuf_dh);
        jxf_sprintf(jxf_msgBuf_dh, "local grid parameters: lowery= %i  uppery= %i", lowery, uppery);
        JXF_SET_INFO(jxf_msgBuf_dh);
        jxf_sprintf(jxf_msgBuf_dh, "local grid parameters: lowerz= %i  upperz= %i", lowerz, upperz);
        JXF_SET_INFO(jxf_msgBuf_dh);
    }
    startRow = mg->first;
    rp[0] = 0;
    for (z = lowerz; z < upperz; z ++)
    {
        for (y = lowery; y < uppery; y ++)
        {
            for (x = lowerx; x < upperx; x ++)
            {
                if (debug)
                {
                    jxf_fprintf(jxf_logFile, "row= %i  x= %i  y= %i  z= %i\n", localRow+startRow+1, x, y, z);
                }
                jxf_getstencil(mg, x, y, z);
                if (threeD)
                {
                    if (z > 0)
                    {
                        naborz1 = jxf_rownum(threeD, x, y, z-1, nx, ny, nz, px, py);
                        cval[idx] = naborz1;
                        aval[idx++] = JXF_FRONT(stencil);
                    }
                }
                if (y > 0)
                {
                    nabory1 = jxf_rownum(threeD, x, y-1, z, nx, ny, nz, px, py);
                    cval[idx] = nabory1;
                    aval[idx++] = JXF_SOUTH(stencil);
                }
                if (x > 0)
                {
                    naborx1 = jxf_rownum(threeD, x-1, y, z, nx, ny, nz, px, py);
                    cval[idx] = naborx1;
                    aval[idx++] = JXF_WEST(stencil);
                }
                cval[idx] = localRow + startRow;
                aval[idx++] = JXF_CENTER(stencil);
                if (x < nx*px-1)
                {
                    naborx2 = jxf_rownum(threeD, x+1, y, z, nx, ny, nz, px, py);
                    cval[idx] = naborx2;
                    aval[idx++] = JXF_EAST(stencil);
                }
                if (y < ny*py-1)
                {
                    nabory2 = jxf_rownum(threeD, x, y+1, z, nx, ny, nz, px, py);
                    cval[idx] = nabory2;
                    aval[idx++] = JXF_NORTH(stencil);
                }
                if (threeD)
                {
                    if (z < nz*pz-1)
                    {
                        naborz2 = jxf_rownum(threeD, x, y, z+1, nx, ny, nz, px, py);
                        cval[idx] = naborz2;
                        aval[idx++] = JXF_BACK(stencil);
                    }
                }
                rhs[localRow] = 0.0;
                ++ localRow;
                rp[localRow] = idx;
                if (!threeD && applyBdry)
                {
                    JXF_Int globalRow = localRow + startRow - 1;
                    JXF_Int offset = rp[localRow-1];
                    JXF_Int len = rp[localRow] - rp[localRow-1];
                    JXF_Real ctr, coeff;
                    
                    if (x == 0)
                    {
                        coeff = mg->A(mg->a, x+hhalf, y, z);
                        ctr = mg->A(mg->a, x-hhalf, y, z);
                        jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                    len, &(rhs[localRow-1]), bcx1, coeff, ctr, naborx2);
                    }
                    else if (x == nx*px-1)
                    {
                        coeff = mg->A(mg->a, x-hhalf, y, z);
                        ctr = mg->A(mg->a, x+hhalf, y, z);
                        jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                    len, &(rhs[localRow-1]), bcx2, coeff, ctr, naborx1);
                    }
                    else if (y == 0)
                    {
                        coeff = mg->B(mg->b, x, y+hhalf, z);
                        ctr = mg->B(mg->b, x, y-hhalf, z);
                        jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                    len, &(rhs[localRow-1]), bcy1, coeff, ctr, nabory2);
                    }
                    else if (y == ny*py-1)
                    {
                        coeff = mg->B(mg->b, x, y-hhalf, z);
                        ctr = mg->B(mg->b, x, y+hhalf, z);
                        jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                    len, &(rhs[localRow-1]), bcy2, coeff, ctr, nabory1);
                    }
                    else if (threeD)
                    {
                        if (z == 0)
                        {
                            coeff = mg->B(mg->b, x, y, z+hhalf);
                            ctr = mg->B(mg->b, x, y, z-hhalf);
                            jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                        len, &(rhs[localRow-1]), bcy1, coeff, ctr, naborz2);
                        }
                        else if (z == nz*nx-1)
                        {
                            coeff = mg->B(mg->b, x, y, z-hhalf);
                            ctr = mg->B(mg->b, x, y, z+hhalf);
                            jxf_setBoundary_private(globalRow, cval+offset, aval+offset,
                                        len, &(rhs[localRow-1]), bcy1, coeff, ctr, naborz1);
                        }
                    }
                }
            }
        }
    }
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_setBoundary_private"
void jxf_setBoundary_private( JXF_Int node,
                          JXF_Int *cval,
                          JXF_Real *aval,
                          JXF_Int len,
                          JXF_Real *rhs,
                          JXF_Real bc,
                          JXF_Real coeff,
                          JXF_Real ctr,
                          JXF_Int nabor )
{
    JXF_START_FUNC_DH
    JXF_Int i;
    
    if (bc >= 0)
    {
       *rhs = bc;
        for (i = 0; i < len; ++ i)
        {
            if (cval[i] == node)
            {
                aval[i] = 1.0;
            }
            else
            {
                aval[i] = 0;
            }
        }
    }
    else
    {
        for (i = 0; i < len; ++ i)
        {
            if (cval[i] == node)
            {
                aval[i] += (ctr - coeff);
            }
            else if (cval[i] == nabor)
            {
                aval[i] = 2.0 * coeff;
            }
        }
    }
    JXF_END_FUNC_DH
}
