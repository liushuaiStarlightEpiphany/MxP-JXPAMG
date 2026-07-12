//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)

/*!
 * \fn JXF_Int jxf_NumThreads
 * \brief Get the number of all the threads
 * \date 2012/09/04
 */
JXF_Int
jxf_NumThreads()
{
   JXF_Int num_threads;
   
#if JXF_USING_OPENMP
#pragma omp parallel
   num_threads = omp_get_num_threads();
#endif
#ifdef JXF_USING_PGCC_SMP
   num_threads = 2;
#endif
   
   return num_threads;
}

/* This next function must be called from within a parallel region! */

JXF_Int
jxf_NumActiveThreads()
{
   JXF_Int num_threads;
   
   num_threads = omp_get_num_threads();
   
   return num_threads;
}

/*!
 * \fn JXF_Int jxf_GetThreadNum
 * \brief Get the id of the current thread
 * \date 2012/09/04
 */
JXF_Int
jxf_GetThreadNum()
{
   JXF_Int my_thread_num;
   
#if JXF_USING_OPENMP
   my_thread_num = omp_get_thread_num();
#endif
#ifdef JXF_USING_PGCC_SMP
   /* THIS NEEDS TO BE FIXED */
   my_thread_num = 0;
#endif
   
   return my_thread_num;
}

#endif

void
jxf_GetSimpleThreadPartition( JXF_Int *begin, JXF_Int *end, JXF_Int n )
{
   JXF_Int num_threads = jxf_NumActiveThreads();
   JXF_Int my_thread_num = jxf_GetThreadNum();

   JXF_Int n_per_thread = (n + num_threads - 1)/num_threads;

   *begin = jxf_min(n_per_thread*my_thread_num, n);
   *end = jxf_min(*begin + n_per_thread, n);
}
