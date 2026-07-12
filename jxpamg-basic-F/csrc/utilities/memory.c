//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

/*!
 * \fn JXF_Int jxf_OutOfMemory
 * \brief Out-of-memory warning.
 * \date 2011/09/01 
 */  
JXF_Int
jxf_OutOfMemory( size_t size )
{
   jxf_printf(" Out of memory trying to allocate %d bytes\n", (JXF_Int) size);
   fflush(stdout);

   jxf_error(JXF_ERROR_MEMORY);

   return 0;
}

/*!
 * \fn char *jxf_MAlloc
 * \brief Allocating memory.
 * \date 2011/09/01 
 */  
char *
jxf_MAlloc( size_t size )
{
   char *ptr;

   if (size > 0)
   {
      ptr = malloc(size);
      if (ptr == NULL)
      {
         jxf_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

/*!
 * \fn char *jxf_CAlloc
 * \brief Allocating memory.
 * \date 2011/09/01 
 */ 
char *
jxf_CAlloc( size_t count, size_t elt_size )
{
   char  *ptr;
   size_t size = count*elt_size;  // Zhou Zhiyang & Yue Xiaoqiang 2012/11/03

   if (size > 0)
   {
      ptr = calloc(count, elt_size);
      if (ptr == NULL)
      {
         jxf_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

/*!
 * \fn char *jxf_ReAlloc
 * \brief Allocating memory.
 * \date 2011/09/01 
 */ 
char *
jxf_ReAlloc( char *ptr, size_t size )
{
   if (ptr == NULL)
      ptr = malloc(size);
   else
      ptr = realloc(ptr, size);

   if ((ptr == NULL) && (size > 0))
   {
      jxf_OutOfMemory(size);
   }

   return ptr;
}

void *
jxf_ReAlloc_v2( void *ptr, size_t old_size, size_t new_size )
{
   if (new_size == 0)
   {
      jxf_Free(ptr);
      return NULL;
   }

   if (ptr == NULL)
   {
      return jxf_MAlloc(new_size);
   }

   void *new_ptr = jxf_MAlloc(new_size);
   size_t smaller_size = new_size > old_size ? old_size : new_size;
   memcpy(new_ptr, ptr, smaller_size);
   jxf_Free(ptr);
   ptr = new_ptr;

   if (!ptr)
   {
      jxf_OutOfMemory(new_size);
   }

   return ptr;
}

/*!
 * \fn void jxf_Free
 * \brief Free memory.
 * \date 2011/09/01 
 */ 
void
jxf_Free( char *ptr )
{
   if (ptr) free(ptr);
}
