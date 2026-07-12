#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER
#include "jxf_hpcsr.h"
#endif

/*!
 * \fn jxf_hpBuildRhsParFromOneFile
 * \author mrz
 * \date 2021/10/26
 */
jxf_ParVector *
jxf_hpBuildRhsParFromOneFile( char *filename, jxf_hpCSRMatrix *A, JXF_Int file_type )
{
    return (jxf_BuildRhsParFromOneFile(filename, jxf_hpCSRMatrixPar(A), file_type));
}