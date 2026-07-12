//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_vector.c -- operations for parallel vectors.
 *  Date: 2014/07/05
 */

#include "jxf_mv.h"
#include "jxf_ilu.h"

/*!
 * \fn jxf_ParVector *jxf_BuildRhsParFromOneFile2
 * \brief Build a parallel vector from a given file
 *        including reorder by variables and along y-axis
 * \param *filename pointer to the file name.
 * \param *A pointer to the matrix which will provide the row partitioning.
 * \note The storage format of the vector is expected to be as follows:
 *         the 1st row:                   n                       // JXF_Int-type
 *         the 2 to (n+1)-th row:         f[i]    i=0(1)n-1       // JXF_Real-type
 *      where n is the total size of the vector, f[i] is the i-th entry of the vector.
 * \author Yue Xiaoqiang
 * \date 2014/04/10
 */
jxf_ParVector *
jxf_BuildRhsParFromOneFile2( char *filename, jxf_ParCSRMatrix *A, JXF_Int num_equns, JXF_Int nx, JXF_Int ny )
{
   jxf_ParVector *b = NULL;
   jxf_Vector    *b_CSR = NULL;
   jxf_Vector    *c_CSR = NULL;
   JXF_Int           my_id;
   JXF_Int          *partitioning;

   jxf_MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

   if (my_id == 0) 
   {
      b_CSR = jxf_SeqVectorRead(filename);
      c_CSR = jxf_VectorTDMGReorderByNodes(b_CSR, num_equns);
      jxf_SeqVectorDestroy(b_CSR);
      b_CSR = jxf_VectorTDMGReorderAlongY(c_CSR, num_equns, nx, ny);
      jxf_SeqVectorDestroy(c_CSR);
   }

   jxf_ParCSRMatrixGetRowPartitioning(A, &partitioning);

   b = jxf_VectorToParVector(MPI_COMM_WORLD, b_CSR, partitioning);

   if (my_id == 0) jxf_SeqVectorDestroy(b_CSR);

   return (b);
}

/*!
 * \fn JXF_Int jxf_GeneratePartitioning2
 * \brief Generates load balanced partitioning of a 1-d array, block-partitioned
 * \author Yue Xiaoqiang
 * \date 2014/04/13
 */
JXF_Int
jxf_GeneratePartitioning2(JXF_Int num_procs, JXF_Int **part_ptr, JXF_Int num_equns, JXF_Int length, JXF_Int ny)
{
   JXF_Int  ierr = 0;
   JXF_Int *part;
   JXF_Int  size, rest;
   JXF_Int  dist, eist;
   JXF_Int  i;

   part = jxf_CTAlloc(JXF_Int, num_procs+1);
   size = length / num_procs;
   rest = length - size * num_procs;
   dist = num_equns * ny;
   eist = dist * size;
   part[0] = 0;
   for (i = 0; i < num_procs; i ++)
   {
	part[i+1] = part[i] + eist;
	if (i < rest) part[i+1] += dist;
   }

  *part_ptr = part;
   return ierr;
}
