//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

void
jxf_prefix_sum_pair(JXF_Int *in_out1, JXF_Int *sum1, JXF_Int *in_out2, JXF_Int *sum2, JXF_Int *workspace)
{
#ifdef JXF_USING_OPENMP
   JXF_Int my_thread_num = jxf_GetThreadNum();
   JXF_Int num_threads = jxf_NumActiveThreads();
   jxf_assert(1 == num_threads || omp_in_parallel());

   workspace[(my_thread_num + 1)*2] = *in_out1;
   workspace[(my_thread_num + 1)*2 + 1] = *in_out2;

#pragma omp barrier
#pragma omp master
   {
      JXF_Int i;
      workspace[0] = 0;
      workspace[1] = 0;

      for (i = 1; i < num_threads; i++)
      {
         workspace[(i + 1)*2] += workspace[i*2];
         workspace[(i + 1)*2 + 1] += workspace[i*2 + 1];
      }
      *sum1 = workspace[num_threads*2];
      *sum2 = workspace[num_threads*2 + 1];
   }
#pragma omp barrier

   *in_out1 = workspace[my_thread_num*2];
   *in_out2 = workspace[my_thread_num*2 + 1];

#else

   *sum1 = *in_out1;
   *sum2 = *in_out2;
   *in_out1 = 0;
   *in_out2 = 0;

   workspace[0] = 0;
   workspace[1] = 0;
   workspace[2] = *sum1;
   workspace[3] = *sum2;
#endif
}
