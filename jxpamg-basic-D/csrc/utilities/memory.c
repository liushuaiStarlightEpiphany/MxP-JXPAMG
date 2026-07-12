//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jx_util.h"
#include "hthread_host.h"

JX_Int jx_malloc_type = 1;
extern int myid;

/*!
 * \fn JX_Int jx_OutOfMemory
 * \brief Out-of-memory warning.
 * \date 2011/09/01
 */
JX_Int
jx_OutOfMemory(size_t size)
{
   jx_printf(" Out of memory trying to allocate %d bytes\n", (JX_Int)size);
   fflush(stdout);

   jx_error(JX_ERROR_MEMORY);

   return 0;
}

/*!
 * \fn char *jx_MAlloc
 * \brief Allocating memory.
 * \date 2011/09/01
 */
char *
jx_MAlloc(size_t size)
{
   char *ptr;
   switch (jx_malloc_type)
   {
   case 1:
      ptr = jx_MAlloc_MT(size);
      break;
   default:
      ptr = jx_MAlloc_origin(size);
      break;
   }
   return ptr;
}

char *
jx_MAlloc_origin(size_t size)
{
   // if (myid == 0)
   //    printf("jx_MAlloc_origin \n");

   char *ptr;

   if (size > 0)
   {
      ptr = malloc(size);
      if (ptr == NULL)
      {
         jx_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

char *
jx_MAlloc_MT(size_t size)
{
   // if (myid == 0)
   //    printf("jx_MAlloc_MT \n");

   char *ptr;
   int cluster_id = myid % 4;

   if (size > 0)
   {
      ptr = hthread_malloc(cluster_id, size, HT_MEM_RW);
      if (ptr == NULL)
      {
         jx_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }
   return ptr;
}

/*!
 * \fn char *jx_CAlloc
 * \brief Allocating memory.
 * \date 2011/09/01
 */
char *
jx_CAlloc(size_t count, size_t elt_size)
{
   // JX_Int jx_malloc_type = 1;

   char *ptr;
   switch (jx_malloc_type)
   {
   case 1:
      ptr = jx_CAlloc_MT(count, elt_size);
      break;
   default:
      ptr = jx_CAlloc_origin(count, elt_size);
      break;
   }
   return ptr;
}

char *
jx_CAlloc_origin(size_t count, size_t elt_size)
{
   char *ptr;
   size_t size = count * elt_size; // Zhou Zhiyang & Yue Xiaoqiang 2012/11/03

   if (size > 0)
   {
      ptr = calloc(count, elt_size);
      if (ptr == NULL)
      {
         jx_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

char *
jx_CAlloc_MT(size_t count, size_t elt_size)
{
   char *ptr;
   int cluster_id = myid % 4;
   size_t size = count * elt_size;

   if (size > 0)
   {
      ptr = hthread_malloc(cluster_id, size, HT_MEM_RW);
      if (ptr == NULL)
      {
         jx_OutOfMemory(size); // 内存分配失败，调用内存不足处理函数
      }
      else
      {
         memset(ptr, 0, size); // 初始化
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

/*!
 * \fn char *jx_ReAlloc
 * \brief Allocating memory.
 * \date 2011/09/01
 */
char *
jx_ReAlloc(char *ptr, size_t size)
{
   // JX_Int jx_malloc_type = 1;

   switch (jx_malloc_type)
   {
   case 1:
      ptr = jx_ReAlloc_MT(ptr, size);
      break;
   default:
      ptr = jx_ReAlloc_origin(ptr, size);
      break;
   }
   return ptr;
}

char *
jx_ReAlloc_origin(char *ptr, size_t size)
{
   if (ptr == NULL)
      ptr = malloc(size);
   else
      ptr = realloc(ptr, size);

   if ((ptr == NULL) && (size > 0))
   {
      jx_OutOfMemory(size);
   }

   return ptr;
}

char *
jx_ReAlloc_MT(char *ptr, size_t size)
{
   int cluster_id = myid % 4;

   if (ptr == NULL)
   {
      ptr = hthread_malloc(cluster_id, size, HT_MEM_RW);
   }
   else
   {
      // 如果原指针不为 NULL，则重新分配内存
      char *new_ptr = hthread_malloc(cluster_id, size, HT_MEM_RW);
      if (new_ptr != NULL)
      {
         memcpy(new_ptr, ptr, size); // 拷贝原内存数据到新内存
         hthread_free(ptr);          // 释放原内存
         ptr = new_ptr;              // 更新指针
      }
      else
      {
         // 内存分配失败时
         jx_OutOfMemory(size);
      }
   }

   return ptr;
}

void *
jx_ReAlloc_v2(void *ptr, size_t old_size, size_t new_size)
{
   if (new_size == 0)
   {
      jx_Free(ptr);
      return NULL;
   }

   if (ptr == NULL)
   {
      return jx_MAlloc(new_size);
   }

   void *new_ptr = jx_MAlloc(new_size);
   size_t smaller_size = new_size > old_size ? old_size : new_size;
   memcpy(new_ptr, ptr, smaller_size);
   jx_Free(ptr);
   ptr = new_ptr;

   if (!ptr)
   {
      jx_OutOfMemory(new_size);
   }

   return ptr;
}

/*!
 * \fn void jx_Free
 * \brief Free memory.
 * \date 2011/09/01
 */
void jx_Free(char *ptr)
{
   switch (jx_malloc_type)
   {
   case 1:
      jx_Free_MT(ptr);
      break;
   default:
      jx_Free_origin(ptr);
      break;
   }
}

void jx_Free_origin(char *ptr)
{
   if (ptr)
      free(ptr);
}

void jx_Free_MT(char *ptr)
{
   if (ptr)
      hthread_free(ptr);
}
