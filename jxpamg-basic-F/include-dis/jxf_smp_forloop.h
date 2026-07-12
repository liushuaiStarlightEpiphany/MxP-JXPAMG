//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_smp_forloop.h -- Wrapper code for SMP compiler directives.  Translates 
 *  jxf SMP directives into the appropriate Open MP,
 *  IBM, SGI, or pgcc (Red) SMP compiler directives
 *  Date: 2011/09/08
 *  Modified Date: 2012/10/22
 *
 *  Created by peghoty
 *  Modified by Yue Xiaoqiang
 */

/* OpenMP */
#if JXF_USING_OPENMP
#ifndef JXF_SMP_MASTER
#ifndef JXF_SMP_BARRIER
#ifndef JXF_SMP_REDUCTION_OP
#ifndef JXF_SMP_PAR_REGION
#ifndef JXF_SMP_FOR
#ifndef JXF_SMP_PRIVATE
#pragma omp parallel for schedule(static)
#else
#pragma omp parallel for private(JXF_SMP_PRIVATE) schedule(static)
#endif
#endif
#endif
#endif
#endif
#endif

// /* OpenMP */
// #if JXF_USING_OPENMP
// #ifndef JXF_SMP_MASTER
// #ifndef JXF_SMP_BARRIER
// #ifndef JXF_SMP_REDUCTION_OP
// #ifndef JXF_SMP_PAR_REGION
// #ifndef JXF_SMP_FOR
// #pragma omp parallel for private(JXF_SMP_PRIVATE) schedule(static)
// #endif
// #endif
// #endif
// #endif
// #endif

#ifdef JXF_SMP_BARRIER
#pragma omp barrier
#endif

#ifdef JXF_SMP_MASTER
#pragma omp master
#endif

#ifdef JXF_SMP_PAR_REGION
#pragma omp parallel private(JXF_SMP_PRIVATE)
#endif

#ifdef JXF_SMP_FOR
#pragma omp for schedule(static)
#endif

#ifdef JXF_SMP_REDUCTION_OP
#pragma omp parallel for private(JXF_SMP_PRIVATE) \
reduction(JXF_SMP_REDUCTION_OP: JXF_SMP_REDUCTION_VARS) \
schedule(static)
#endif
#endif

/* SGI */
#ifdef JXF_USING_SGI_SMP
#pragma parallel
#pragma pfor
#pragma schedtype(gss)
#pragma chunksize(10)
#endif

/* IBM */
#ifdef JXF_USING_IBM_SMP
#pragma parallel_loop
#pragma schedule (guided,10)
#endif

/* PGCC */
#ifdef JXF_USING_PGCC_SMP
#ifndef JXF_SMP_REDUCTION_OP
#pragma parallel local(JXF_SMP_PRIVATE) pfor
#endif
#ifdef JXF_SMP_REDUCTION_OP
#endif
#endif

/* UnDef */
#undef JXF_SMP_PRIVATE
#undef JXF_SMP_REDUCTION_OP
#undef JXF_SMP_REDUCTION_VARS
#undef JXF_SMP_PAR_REGION
#undef JXF_SMP_FOR
#undef JXF_SMP_MASTER
#undef JXF_SMP_BARRIER
