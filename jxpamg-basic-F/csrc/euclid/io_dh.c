//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  io_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

#undef __FUNC__
#define __FUNC__ "jxf_openFile_dh"
FILE *jxf_openFile_dh( const char *filenameIN, const char *modeIN )
{
    JXF_START_FUNC_DH
    FILE *fp = NULL;
    
    if ((fp = fopen(filenameIN, modeIN)) == NULL)
    {
        jxf_sprintf(jxf_msgBuf_dh, "can't open file: %s for mode %s\n", filenameIN, modeIN);
        JXF_SET_ERROR(NULL, jxf_msgBuf_dh);
    }
    JXF_END_FUNC_VAL(fp)
}

#undef __FUNC__
#define __FUNC__ "jxf_closeFile_dh"
void jxf_closeFile_dh( FILE *fpIN )
{
    if (fclose(fpIN))
    {
        JXF_SET_V_ERROR("attempt to close file failed");
    }
}

void jxf_io_dh_print_ebin_mat_private( JXF_Int m,
                                   JXF_Int beg_row,
                                   JXF_Int *rp,
                                   JXF_Int *cval,
                                   JXF_Real *aval,
                                   JXF_Int *n2o,
                                   JXF_Int *o2n,
                                   jxf_Hash_i_dh hash,
                                   char *filename )
{
}

extern void jxf_io_dh_read_ebin_mat_private( JXF_Int *m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval, char *filename )
{
}

void jxf_io_dh_print_ebin_vec_private( JXF_Int n,
                                   JXF_Int beg_row,
                                   JXF_Real *vals,
                                   JXF_Int *n2o,
                                   JXF_Int *o2n,
                                   jxf_Hash_i_dh hash,
                                   char *filename )
{
}

void jxf_io_dh_read_ebin_vec_private( JXF_Int *n, JXF_Real **vals, char *filename )
{
}
