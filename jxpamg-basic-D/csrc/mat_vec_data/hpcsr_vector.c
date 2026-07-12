#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER
#include "jx_hpcsr.h"
#endif

/*!
 * \fn jx_hpBuildRhsParFromOneFile
 * \author mrz
 * \date 2021/10/26
 */
jx_ParVector *
jx_hpBuildRhsParFromOneFile( char *filename, jx_hpCSRMatrix *A, JX_Int file_type )
{
    return (jx_BuildRhsParFromOneFile(filename, jx_hpCSRMatrixPar(A), file_type));
}