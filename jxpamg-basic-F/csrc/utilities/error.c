//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

#define JXF_PRINT_ERRORS 1

/*!
 * \fn void jxf_error_handler
 * \brief Process the error with code ierr raised in the given 
 *        line of the given source file.
 * \date 2011/09/01 
 */  
void 
jxf_error_handler( char *filename, JXF_Int line, JXF_Int ierr, const char *msg )
{
   jxf_error_flag |= ierr;

#ifdef JXF_PRINT_ERRORS
   if (msg)
   {
      jxf_fprintf(stderr, "jxf-error in file \"%s\", line %d, error code = %d - %s\n", filename, line, ierr, msg);
   }
   else
   {
      jxf_fprintf(stderr, "jxf-error in file \"%s\", line %d, error code = %d\n", filename, line, ierr);
   }
#endif
}
