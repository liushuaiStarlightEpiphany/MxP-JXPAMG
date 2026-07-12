/*========================================================================*
 *  FSLS (Fast Solvers for Linear System)  (c) 2010-2012                  *
 *  School of Mathematics and Computational Science                       *
 *  Xiangtan University                                                   *
 *  Email: peghoty@163.com                                                *
 *========================================================================*/
 
/*!
 * jxf_multils.h
 *
 * Created by peghoty 2011/09/29
 * Xiangtan University
 * peghoty@163.com 
 *
 */

#ifndef FSLS_HEADER
#define FSLS_HEADER

#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#ifndef JXF_SEQUENTIAL
#include "mpi.h"
#endif

/*----------------------------------------------------------------*
 *                      Macro Definition                          *
 *----------------------------------------------------------------*/ 

struct fsls_Solver_struct;
typedef struct fsls_Solver_struct *FSLS_Solver;
 
/* time testing */  
#define GetTime(a) gettimeofday(&a,NULL)
#define mytime(a,b) ((b.tv_sec-a.tv_sec) + (float)(b.tv_usec-a.tv_usec)/1000000.0) 
#define DECIMALPlACE 2

/* ratio of the circumference of a circle to its diameter */
#define PI 3.141592653589793232

/* max and min  */ 
#define fsls_max(a,b)  (((a)<(b)) ? (b) : (a))
#define fsls_min(a,b)  (((a)<(b)) ? (a) : (b))

/* error message */
#define FSLS_ERROR_GENERIC      1
#define FSLS_ERROR_MEMORY       2
#define FSLS_ERROR_CONVERGENCE  3

/* print the position for program debugging */
#define fsls_CurrentFile()    jxf_printf("\n >>> Current File: \033[31m%s\033[00m\n", __FILE__)
#define fsls_CurrentFuction() jxf_printf("\n >>> Current Func: \033[31m%s\033[00m\n", __func__)
#define fsls_CurrentLine()    jxf_printf("\n >>> Current Line: \033[31m%d\033[00m\n", __LINE__)

/* system function to excute DOS commands */ 
#define fsls_CurrentDirTree()    system("tree")
#define fsls_CurrentDirectory()  system("dir")

/* debugging */
#define fsls_Debug(a) jxf_printf(" \033[4mFSLSDEBUG @ peghoty :\033[00m \033[31m%d\033[00m\n\n",(a))

/* allocate and free memory */
#define fsls_TFree(ptr) ( fsls_Free((char *)ptr), ptr = NULL )
#define fsls_CTAlloc(type, count) ( (type *)fsls_CAlloc((size_t)(count), (size_t)sizeof(type)) )

/* Approximate maximal number of DOF for GE solver */
#define GEMAXDOF 10000

/* check the zero diagonal elements during FEBS setup or not? */
#define NoZeroDiagCheck 1 // 0: with checking on; 1: without checking

/* small real number */
#define SMALLREAL 1.0e-15

/* pivot type */
#define NOPIVOT       0
#define COLUMNPIVOT   1
#define COMPLETEPIVOT 2 

/* point type for C/F splitting */
#define CPOINT         1
#define FPOINT        -1
#define SFPOINT       -3
#define COMMON_CPOINT  2
#define UNDECIDED      0

/* relaxation type */
#define Smooth_WJACOBI  0   // Jacobi 
#define Smooth_GS       1   // Gauss-Seidel
#define Smooth_SGS      2   // symmetric Gausss-Seidel
#define Smooth_SOR      3   // SOR(succesive over relaxation)
#define Smooth_KACZMARZ 4   // Kaczmarz
#define Smooth_POLY     5   // Polynomial relaxation
#define Smooth_MUMPS    6   // Direct Solver in MUMPS 
#define Smooth_ILUP     8   // ILU(p) 
#define Smooth_GE       9   // Gaussian Elimination
#define Smooth_GEP     10   // Gaussian Elimination with pivoting
 
/* relaxation order */
#define  CPFIRST       1   // C points first
#define  FPFIRST      -1   // F points first
#define  CP_ONLY       2   // C points only
#define  FP_ONLY      -2   // F points only
#define  ASCEND        3   // ascending order
#define  DESCEND      -3   // descending order 

/* cycle phase parameter */
#define  FINEST_LEFT   0   // finesest(L)
#define  DOWN_CYCLE    1   // down cycle
#define  COARSEST      2   // coarsest
#define  UP_CYCLE      3   // up cycle
#define  FINEST_RIGHT  4   // finesest(R) 

/* coarsening type */
#define  RS_ORIGINAL  -1   // original RS from hypre
#define  RS_MODIFIED   0   // modified RS by peghoty
#define  CLJP          1   // CLJP

/* interpolation type */
#define  CLASSIC       0   // modified classic interpolation
#define  DIRECT        1   // direct interpolation
#define  ENERGYMIN     2   // Energy-Minimizing based interpolation

/* GMRES type */
#define  CLASSICGMRES  0   // classic GMRES
#define  ADAPTIVEGMRES 1   // adaptive GMRES
#define  LOOSEGMRES    2   // Loose GMRES

#define  SUCCESS       0

#define FSLS_LIST_HEAD -1
#define FSLS_LIST_TAIL -2





/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/
 
struct fsls_double_linked_list
{
   JXF_Int                             data;
   struct fsls_double_linked_list *next_elt;
   struct fsls_double_linked_list *prev_elt;
   JXF_Int                             head;
   JXF_Int                             tail;
};

typedef struct fsls_double_linked_list fsls_ListElement;
typedef fsls_ListElement  *fsls_LinkList; 

typedef struct 
{
   JXF_Int      num_rows;      // number of rows
   JXF_Int      num_cols;      // number of columns 
   JXF_Real  *data;          // data part of the matrix stored row by row
	
} fsls_DenseMatrix;

#define fsls_DenseMatrixNumRows(full)   ((full) -> num_rows)
#define fsls_DenseMatrixNumCols(full)   ((full) -> num_cols)
#define fsls_DenseMatrixData(full)      ((full) -> data)

/*!
 * \struct fsls_CSRMatrix
 */
typedef struct
{
   JXF_Real  *data;         // a real array, contains the real value a_{ij}
                          // stored row by row, from row 0 to row n-1.
   JXF_Int     *i;            // an integer array, contains the pointers to the beginning
                          // of each row in the arrays 'data' or 'j'. Thus, the
                          // content of i[k],k = 0(1)n-1, is the position in arrays
                          // 'data' or 'j' where the k-th row starts, and
                          // i[n] = i[0] + num_nonzeros can be taken as the address in 'data'
                          // or 'j' of the beginning of a fictitious row number n.
   JXF_Int     *j;            // an integer array, contains the column indices of the
                          // elements a_{ij} as stored in the array 'data'.
   JXF_Int      num_rows;     // number of rows
   JXF_Int      num_cols;     // number of columns
   JXF_Int      num_nonzeros; // number of nonzero entries

} fsls_CSRMatrix;

#define fsls_CSRMatrixData(matrix)         ((matrix) -> data)
#define fsls_CSRMatrixI(matrix)            ((matrix) -> i)
#define fsls_CSRMatrixJ(matrix)            ((matrix) -> j)
#define fsls_CSRMatrixNumRows(matrix)      ((matrix) -> num_rows)
#define fsls_CSRMatrixNumCols(matrix)      ((matrix) -> num_cols)
#define fsls_CSRMatrixNumNonzeros(matrix)  ((matrix) -> num_nonzeros)

/*!
 * \struct fsls_CSRGroup
 */
typedef struct
{
   JXF_Int               n;       // number of CSR matrix
   fsls_CSRMatrix  **A_array; // A_array[i],i=0:n-1, the i-th CSR matrix in the group

} fsls_CSRGroup;

#define fsls_CSRGroupN(csrgroup)         ((csrgroup) -> n)
#define fsls_CSRGroupAArray(csrgroup)    ((csrgroup) -> A_array)

/*!
 * \struct fsls_ICSRMatrix
 */
typedef struct
{
   JXF_Int     *data;         // a integer array, contains the real value a_{ij}
                          // stored row by row, from row 0 to row n-1.
   JXF_Int     *i;            // an integer array, contains the pointers to the beginning
                          // of each row in the arrays 'data' or 'j'. Thus, the
                          // content of i[k],k = 0(1)n-1, is the position in arrays
                          // 'data' or 'j' where the k-th row starts, and
                          // i[n] = i[0] + num_nonzeros can be taken as the address in 'data'
                          // or 'j' of the beginning of a fictitious row number n.
   JXF_Int     *j;            // an integer array, contains the column indices of the
                          // elements a_{ij} as stored in the array 'data'.
   JXF_Int      num_rows;     // number of rows
   JXF_Int      num_cols;     // number of columns
   JXF_Int      num_nonzeros; // number of nonzero entries
 
} fsls_ICSRMatrix;

#define fsls_ICSRMatrixData(matrix)         ((matrix) -> data)
#define fsls_ICSRMatrixI(matrix)            ((matrix) -> i)
#define fsls_ICSRMatrixJ(matrix)            ((matrix) -> j)
#define fsls_ICSRMatrixNumRows(matrix)      ((matrix) -> num_rows)
#define fsls_ICSRMatrixNumCols(matrix)      ((matrix) -> num_cols)
#define fsls_ICSRMatrixNumNonzeros(matrix)  ((matrix) -> num_nonzeros)

/*!
 * \struct fsls_Vector
 */
typedef struct
{
   JXF_Real  *data;
   JXF_Int      size;

} fsls_Vector;

#define fsls_VectorData(vector) ((vector) -> data)
#define fsls_VectorSize(vector) ((vector) -> size)

/*!
 * \struct fsls_VecGroup
 */
typedef struct
{
   JXF_Int            n;       // number of Vector
   fsls_Vector  **x_array; // x_array[i],i=0:n-1, the i-th Vector in the group

} fsls_VecGroup;

#define fsls_VecGroupN(vecgroup)         ((vecgroup) -> n)
#define fsls_VecGroupXArray(vecgroup)    ((vecgroup) -> x_array)

/*!
 * \struct fsls_IVector
 */
typedef struct
{
   JXF_Int  *data;
   JXF_Int   size;

} fsls_IVector;

#define fsls_IVectorData(vector) ((vector) -> data)
#define fsls_IVectorSize(vector) ((vector) -> size)

/*!
 * \struct fsls_RelaxILUData
 */
typedef struct
{
   JXF_Int    **index_array;  // index part of the L/U in MSR format 
   JXF_Real **value_array;  // value part of the L/U in MSR format
   
   JXF_Real  *work;
} fsls_RelaxILUData;

#define fsls_RelaxILUDataIndexArray(relax_ilu_data)      ((relax_ilu_data) -> index_array)
#define fsls_RelaxILUDataValueArray(relax_ilu_data)      ((relax_ilu_data) -> value_array)
#define fsls_RelaxILUDataWork(relax_ilu_data)            ((relax_ilu_data) -> work)

/*!
 * \struct fsls_IterILUpData
 */
typedef struct 
{
   JXF_Int          level_fill_in;  // level of fill in in ILU(p) decomposition
   JXF_Int         *index;          // the index part of the L/U in MSR format
   JXF_Real      *value;          // the value part of the L/U in MSR format
   JXF_Real      *work;           // working array in the triangular system solving
   fsls_Vector *res;            // for the residual vector
   
} fsls_IterILUpData;

#define fsls_IterILUpDataLevelFillIn(iter_ilup_data)   ((iter_ilup_data) -> level_fill_in)
#define fsls_IterILUpDataIndex(iter_ilup_data)         ((iter_ilup_data) -> index)
#define fsls_IterILUpDataValue(iter_ilup_data)         ((iter_ilup_data) -> value)
#define fsls_IterILUpDataWork(iter_ilup_data)          ((iter_ilup_data) -> work)
#define fsls_IterILUpDataRes(iter_ilup_data)           ((iter_ilup_data) -> res)

/*!
 * \struct fsls_AMGData
 */
typedef struct
{
   /* setup params */
   JXF_Int      max_levels;        // maximal number of grid levels
   JXF_Real   strong_threshold;  // defines the strength matrix
   JXF_Real   max_row_sum;       // parameter to modify the definition of strength for diagonal dominant
                               // portions of the matrix. The default is 0.9. If max_row_sum = 1.0, no 
                               // checking for diagonally dominant rows is performed. 2010/11/11
   JXF_Real   A_trunc_factor;    // defines operator truncation factor
   JXF_Real   P_trunc_factor;    // defines interpolation truncation factor
   JXF_Int      A_max_elmts;       // defines max coeffs per row in operator
   JXF_Int      P_max_elmts;       // defines max coeffs per row in interpolation
   JXF_Int      coarsen_type;      // coarsening type
   JXF_Int      interp_type;       // interpolation type
   JXF_Int      coarse_threshold;  // number of dofs allowed on the coarsest level
   JXF_Int      S_mode;            // how to define the strength matrix S     

   /* solve params */
   JXF_Real   tol;               // tolerance    
   JXF_Int      max_iter;          // maximum number of iteration
   JXF_Int      min_iter;          // minimum number of iteration 
   JXF_Int      cycle_type;        // cycle type
   JXF_Int      corrective_type;   // corrective type or iterative type when cycling
   JXF_Int     *grid_relax_sweeps; // number of relaxation sweeps on each level
   JXF_Int     *grid_relax_type;   // type of relaxation on each level  
   JXF_Int     *grid_relax_order;  // order of relaxation on each level
   JXF_Real  *grid_relax_weight; // weight of relaxation on each level
   JXF_Real   relax_weight;      // relaxation weight
   JXF_Int      relax_type;        // relaxation type
   JXF_Int      relax_order;       // relaxation order    
   JXF_Int      relax_sweeps;      // relaxation sweeps
   JXF_Int      relax_sym;         // symmetricly relax or not? 
   JXF_Int      poly_degree;       // degree of polynomial for Polynomial relaxation
   JXF_Int      coarsest_solver;   // solver type on the coarsest level  
   JXF_Int      cycle_count;       // V-Cycle (W-Cycle etc.) iterations
   JXF_Real   last_rel_nrm;      // last relative residual norm
   JXF_Real   ave_conv_factor;   // average convergence factor
   JXF_Real  *full;              // full matrix for direct solver on the coarsest level
   JXF_Real   rhsnorm_threshold; // for the iteration control, the following strategy is used 
                               // when rhsnorm_threshold != 0:
                               // 1. if ||rhs|| > rhsnorm_threshold, use relative norm control
                               // 2. otherwise, use absolute norm control
                               // default of rhsnorm_threshold: 0.0

   /* data generated in the setup phase */
   fsls_CSRMatrix      *prematrix;        // preconditioner
   fsls_CSRMatrix     **A_array;          // grid operaters on each level
   fsls_Vector        **F_array;          // right hand side vectors on each level
   fsls_Vector        **U_array;          // approximate solution vectors on each level
   fsls_CSRMatrix     **P_array;          // interpolations on each level
   JXF_Int                **CF_marker_array;  // CF marker on each level
   JXF_Int                  num_levels;       // true number of grid levels
   fsls_Vector        **RNrm_array;       // row-norm data for Kaczmarz relaxation
   JXF_Int                  level_of_fill_in; // level of fill in for ILU(p) decomposition
   JXF_Int                  ilu_smooth_level; // number of levels use that employ ILU as smoother   
   fsls_RelaxILUData   *relax_ilu_data;   // data for ILU relaxation

   /* auxiliary data */
   fsls_Vector *Vtemp;  // temporary vector with the same size as 'F_array[0]'
   fsls_Vector *temp1;  // temporary vector for polynomial relaxation
   fsls_Vector *temp2;  // temporary vector for polynomial relaxation
   fsls_Vector *temp3;  // temporary vector for polynomial relaxation
   
   /* output params */
   JXF_Int  print_level;   // decide how much information to be output

} fsls_AMGData;

/* setup params */	  		      
#define fsls_AMGDataMaxLevels(amg_data)       ((amg_data) -> max_levels)
#define fsls_AMGDataStrongThreshold(amg_data) ((amg_data) -> strong_threshold)
#define fsls_AMGDataMaxRowSum(amg_data)       ((amg_data) -> max_row_sum)
#define fsls_AMGDataATruncFactor(amg_data)    ((amg_data) -> A_trunc_factor)
#define fsls_AMGDataPTruncFactor(amg_data)    ((amg_data) -> P_trunc_factor)
#define fsls_AMGDataAMaxElmts(amg_data)       ((amg_data) -> A_max_elmts)
#define fsls_AMGDataPMaxElmts(amg_data)       ((amg_data) -> P_max_elmts)
#define fsls_AMGDataCoarsenType(amg_data)     ((amg_data) -> coarsen_type)
#define fsls_AMGDataInterpType(amg_data)      ((amg_data) -> interp_type)
#define fsls_AMGDataCoarseThreshold(amg_data) ((amg_data) -> coarse_threshold)
#define fsls_AMGDataSMode(amg_data)           ((amg_data) -> S_mode)
/* solve params */
#define fsls_AMGDataTol(amg_data)             ((amg_data) -> tol)
#define fsls_AMGDataMaxIter(amg_data)         ((amg_data) -> max_iter)
#define fsls_AMGDataMinIter(amg_data)         ((amg_data) -> min_iter)
#define fsls_AMGDataCycleType(amg_data)       ((amg_data) -> cycle_type)
#define fsls_AMGDataCorrectiveType(amg_data)  ((amg_data) -> corrective_type)
#define fsls_AMGDataGridRelaxSweeps(amg_data) ((amg_data) -> grid_relax_sweeps)
#define fsls_AMGDataGridRelaxType(amg_data)   ((amg_data) -> grid_relax_type)
#define fsls_AMGDataGridRelaxOrder(amg_data)  ((amg_data) -> grid_relax_order)
#define fsls_AMGDataGridRelaxWeight(amg_data) ((amg_data) -> grid_relax_weight)
#define fsls_AMGDataRelaxWeight(amg_data)     ((amg_data) -> relax_weight)
#define fsls_AMGDataRelaxType(amg_data)       ((amg_data) -> relax_type)
#define fsls_AMGDataRelaxOrder(amg_data)      ((amg_data) -> relax_order)
#define fsls_AMGDataRelaxSweeps(amg_data)     ((amg_data) -> relax_sweeps)
#define fsls_AMGDataRelaxSym(amg_data)        ((amg_data) -> relax_sym)
#define fsls_AMGDataPolyDegree(amg_data)      ((amg_data) -> poly_degree)
#define fsls_AMGDataCoarsestSolver(amg_data)  ((amg_data) -> coarsest_solver)
#define fsls_AMGDataCycleCount(amg_data)      ((amg_data) -> cycle_count)
#define fsls_AMGDataLastRelNrm(amg_data)      ((amg_data) -> last_rel_nrm)
#define fsls_AMGDataAveConvFactor(amg_data)   ((amg_data) -> ave_conv_factor)
#define fsls_AMGDataFull(amg_data)            ((amg_data) -> full)
#define fsls_AMGDataRhsNrmThreshold(amg_data) ((amg_data) -> rhsnorm_threshold)
/* data generated by the setup phase */
#define fsls_AMGDataPreMatrix(amg_data)       ((amg_data) -> prematrix)
#define fsls_AMGDataAArray(amg_data)          ((amg_data) -> A_array)
#define fsls_AMGDataFArray(amg_data)          ((amg_data) -> F_array)
#define fsls_AMGDataUArray(amg_data)          ((amg_data) -> U_array)
#define fsls_AMGDataPArray(amg_data)          ((amg_data) -> P_array)
#define fsls_AMGDataCFMarkerArray(amg_data)   ((amg_data) -> CF_marker_array)
#define fsls_AMGDataNumLevels(amg_data)       ((amg_data) -> num_levels)
#define fsls_AMGDataRNrmArray(amg_data)       ((amg_data) -> RNrm_array)
#define fsls_AMGDataLevelOfFillIn(amg_data)   ((amg_data) -> level_of_fill_in)
#define fsls_AMGDataILUSmoothLevel(amg_data)  ((amg_data) -> ilu_smooth_level)   
#define fsls_AMGDataRelaxILUData(amg_data)    ((amg_data) -> relax_ilu_data)
/* auxiliary data */
#define fsls_AMGDataVtemp(amg_data)           ((amg_data) -> Vtemp)
#define fsls_AMGDataTemp1(amg_data)           ((amg_data) -> temp1)
#define fsls_AMGDataTemp2(amg_data)           ((amg_data) -> temp2)
#define fsls_AMGDataTemp3(amg_data)           ((amg_data) -> temp3)
/* output parameters */
#define fsls_AMGDataPrintLevel(amg_data)      ((amg_data) -> print_level)

typedef struct
{
   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();
   
} fsls_CSRPCGFunctions;

typedef struct
{
   JXF_Real   tol;
   JXF_Real   atolf;
   JXF_Real   cf_tol;
   JXF_Real   a_tol;
   JXF_Int      max_iter;
   JXF_Int      two_norm;
   JXF_Int      rel_change;
   JXF_Int      recompute_residual;
   JXF_Int      stop_crit;
   JXF_Int      converged;

   void    *p;
   void    *s;
   void    *r;

   void    *precond_data;

   fsls_CSRPCGFunctions *functions;

   JXF_Int      num_iterations;
   JXF_Real   rel_residual_norm;

   JXF_Int      print_level;
   JXF_Real  *norms;
   JXF_Real  *rel_norms;   
   
} fsls_CSRPCGData;

typedef struct
{
   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();

} fsls_CSRPGMRESFunctions;


typedef struct
{
   JXF_Int     k_dim;
   JXF_Int     min_iter;
   JXF_Int     max_iter;
   JXF_Int     rel_change;
   JXF_Int     stop_crit;
   JXF_Int     converged;
   JXF_Real  tol;
   JXF_Real  cf_tol;
   JXF_Real  a_tol;
   JXF_Real  rel_residual_norm;

   void   *r;
   void   *w;
   void   *w_2;
   void  **p;
   
   void   *precond_data;

   fsls_CSRPGMRESFunctions *functions;

   JXF_Int     num_iterations;
 
   JXF_Int     print_level;
   JXF_Real *norms;

} fsls_CSRPGMRESData;

typedef struct
{
   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();

} fsls_CSRPLGMRESFunctions;


typedef struct
{
   JXF_Int     k_dim;
   JXF_Int     min_iter;
   JXF_Int     max_iter;
   JXF_Int     stop_crit;
   JXF_Int     converged;
   JXF_Real  tol;
   JXF_Real  cf_tol;
   JXF_Real  a_tol;
   JXF_Real  rel_residual_norm;

   /*lgmres specific stuff */
   JXF_Int     aug_dim;
   JXF_Int     approx_constant;
   void  **aug_vecs;
   JXF_Int    *aug_order;
   void  **a_aug_vecs;

   void   *r;
   void   *w;
   void  **p;
   
   void   *precond_data;

   fsls_CSRPLGMRESFunctions *functions;

   JXF_Int     num_iterations;
 
   JXF_Int     print_level;
   JXF_Real *norms;

} fsls_CSRPLGMRESData;

typedef struct
{
   JXF_Int    (*precond)();
   JXF_Int    (*precond_setup)();

} fsls_CSRPAGMRESFunctions;


typedef struct
{
   JXF_Int     k_dim_max;  // upper bound for restart in each restart cycle
   JXF_Int     k_dim_min;  // lower bound for restart in each restart cycle (should be small)
   JXF_Int     k_dim_d;    // reduction for the restart parameter
   JXF_Int     min_iter;
   JXF_Int     max_iter;
   JXF_Int     stop_crit;
   JXF_Int     converged;
   JXF_Real  cr_max;    // maximum convergence rate, = cos(8^o)  (experimental) 
   JXF_Real  cr_min;    // minimum convergence rate, = cos(80^o) (experimental) 
   JXF_Real  tol;
   JXF_Real  cf_tol;
   JXF_Real  rel_residual_norm;

   void   *r;
   void   *w;
   void  **p;
   
   void   *precond_data;

   fsls_CSRPAGMRESFunctions *functions;

   JXF_Int     num_iterations;
 
   JXF_Int     print_level;
   JXF_Real *norms;

} fsls_CSRPAGMRESData;

typedef struct
{
   JXF_Int (*precond_setup)();
   JXF_Int (*precond)();

} fsls_CSRPBICGSTABFunctions;

typedef struct
{
   JXF_Int     min_iter;
   JXF_Int     max_iter;
   JXF_Int     stop_crit;
   JXF_Int     converged;
   JXF_Real  tol;
   JXF_Real  cf_tol;
   JXF_Real  rel_residual_norm;
   JXF_Real  a_tol;

   void   *r;
   void   *r0;
   void   *s;
   void   *v;
   void   *p;
   void   *q;
   
   void   *precond_data;

   fsls_CSRPBICGSTABFunctions *functions;

   JXF_Int     num_iterations;
   JXF_Int     print_level;
   JXF_Real *norms;
   
} fsls_CSRPBICGSTABData;

/*!
 * \struct fsls_PreILUpData
 */
typedef struct 
{
   JXF_Int             level_fill_in;  // level of fill in in ILU(p) decomposition
   JXF_Int            *index;          // the index part of the L/U in MSR format
   JXF_Real         *value;          // the value part of the L/U in MSR format
   JXF_Real         *work;           // working array in the triangular system solving
   fsls_CSRMatrix *prematrix;      // preconditioner
   
} fsls_PreILUpData;

#define fsls_PreILUpDataLevelFillIn(pre_ilup_data)   ((pre_ilup_data) -> level_fill_in)
#define fsls_PreILUpDataIndex(pre_ilup_data)         ((pre_ilup_data) -> index)
#define fsls_PreILUpDataValue(pre_ilup_data)         ((pre_ilup_data) -> value)
#define fsls_PreILUpDataWork(pre_ilup_data)          ((pre_ilup_data) -> work)
#define fsls_PreILUpDataPreMatrix(pre_ilup_data)     ((pre_ilup_data) -> prematrix)

/*!
 * \struct fsls_PreDSData
 */
typedef struct 
{
  fsls_CSRMatrix  *prematrix; // preconditioner
  	
} fsls_PreDSData;

#define fsls_PreDSDataPreMatrix(pre_ds_data)    ((pre_ds_data) -> prematrix)

/*!
 * \struct fsls_PreDiagData
 */
typedef struct
{
   JXF_Int                m1;             // number of DOFs of the first diagonal block
   JXF_Int                m2;             // number of DOFs of the second diagonal block
   JXF_Int                level_fill_in;  // level of fill-in for ILU(k)
   fsls_CSRMatrix    *prematrix;      // preconditioner
   fsls_CSRMatrix    *A1;             // the first diagonal block
   fsls_CSRMatrix    *A2;             // the second diagonal block 
   fsls_AMGData      *amg_solver1;    // amg_solver for the first diagonal block
   fsls_AMGData      *amg_solver2;    // amg_solver for the second diagonal block
   fsls_PreILUpData  *ilu_data1;      // ilu_data for the first diagonal block
   fsls_PreILUpData  *ilu_data2;      // ilu_data for the second diagonal block
   
   /* auxiliary vectors (only pointer, no data part) */
   fsls_Vector       *z1;
   fsls_Vector       *z2;
   fsls_Vector       *r1;
   fsls_Vector       *r2;

} fsls_PreDiagData;
	  		      
#define fsls_PreDiagDataM1(pre_diag_data)           ((pre_diag_data) -> m1)
#define fsls_PreDiagDataM2(pre_diag_data)           ((pre_diag_data) -> m2)
#define fsls_PreDiagDataLevelFillIn(pre_diag_data)  ((pre_diag_data) -> level_fill_in)
#define fsls_PreDiagDataPreMatrix(pre_diag_data)    ((pre_diag_data) -> prematrix)
#define fsls_PreDiagDataA1(pre_diag_data)           ((pre_diag_data) -> A1)
#define fsls_PreDiagDataA2(pre_diag_data)           ((pre_diag_data) -> A2)
#define fsls_PreDiagDataAMGSolver1(pre_diag_data)   ((pre_diag_data) -> amg_solver1)
#define fsls_PreDiagDataAMGSolver2(pre_diag_data)   ((pre_diag_data) -> amg_solver2)
#define fsls_PreDiagDataILUData1(pre_diag_data)     ((pre_diag_data) -> ilu_data1)
#define fsls_PreDiagDataILUData2(pre_diag_data)     ((pre_diag_data) -> ilu_data2)
#define fsls_PreDiagDataZ1(pre_diag_data)           ((pre_diag_data) -> z1)
#define fsls_PreDiagDataZ2(pre_diag_data)           ((pre_diag_data) -> z2)
#define fsls_PreDiagDataR1(pre_diag_data)           ((pre_diag_data) -> r1)
#define fsls_PreDiagDataR2(pre_diag_data)           ((pre_diag_data) -> r2)

typedef struct
{
   MPI_Comm          comm;          // communicator
   JXF_Int               local_num_ls;  // number of sub linear systems for the current processor
   JXF_Int               global_num_ls; // number of all the sub linear systems    
   JXF_Int               num_rows;      // number of rows of each coefficient matrix 
   JXF_Int               num_cols;      // number of cols of each coefficient matrix  
   JXF_Int               num_nonzeros;  // number of nonzeros of each coefficient matrix
   JXF_Int              *ls_partition;  // the partition for linear system data distributing
   JXF_Int              *num_ls_procs;  // num_ls_procs[i],i=0:nprocs-1, number of LS's in the i-th processor            
   fsls_CSRMatrix  **A_array;       // A_array[i],i=0:local_num_ls-1: the i-th coefficient matrix
   fsls_Vector     **b_array;       // b_array[i],i=0:local_num_ls-1: the i-th right hand side vector
   fsls_Vector     **x_array;       // x_array[i],i=0:local_num_ls-1: the i-th solution vector

} fsls_ParLSData;

#define fsls_ParLSDataComm(ls_data)         ((ls_data) -> comm)
#define fsls_ParLSDataLocalNumLS(ls_data)   ((ls_data) -> local_num_ls)
#define fsls_ParLSDataGlobalNumLS(ls_data)  ((ls_data) -> global_num_ls)
#define fsls_ParLSDataNumRows(ls_data)      ((ls_data) -> num_rows) 
#define fsls_ParLSDataNumCols(ls_data)      ((ls_data) -> num_cols) 
#define fsls_ParLSDataNumNonzeros(ls_data)  ((ls_data) -> num_nonzeros) 
#define fsls_ParLSDataLSPartition(ls_data)  ((ls_data) -> ls_partition)
#define fsls_ParLSDataNumLSProcs(ls_data)   ((ls_data) -> num_ls_procs) 
#define fsls_ParLSDataAArray(ls_data)       ((ls_data) -> A_array)
#define fsls_ParLSDataBArray(ls_data)       ((ls_data) -> b_array)
#define fsls_ParLSDataXArray(ls_data)       ((ls_data) -> x_array)

typedef struct
{
   MPI_Comm          comm;            // communicator
   JXF_Int               global_num_ls;   // number of all the sub linear systems    
   JXF_Int               global_num_dof;  // global number of dofs(grids) 
   JXF_Int               local_num_dof;   // local  number of dofs(grids)
   JXF_Int              *dof_partition;   // the partition for dof distributing
   JXF_Int              *num_dof_procs;   // num_dof_proc[i],i=0:nprocs-1, number of DOFs in the i-th processor           
   fsls_CSRMatrix  **subA_array;      // subA_array[i],i=0:global_num_ls-1: parts of the matrix of 
                                      // the i-th sub linear systems
   fsls_Vector     **subb_array;      // subb_array[i],i=0:global_num_ls-1: parts of the right hand side 
                                      // of the i-th sub linear systems
   fsls_Vector     **subx_array;      // subx_array[i],i=0:global_num_ls-1: parts of the approximation  
                                      // of the i-th sub linear systems                                      

} fsls_ParGLSData;

#define fsls_ParGLSDataComm(gls_data)            ((gls_data) -> comm)
#define fsls_ParGLSDataGlobalNumLS(gls_data)     ((gls_data) -> global_num_ls)
#define fsls_ParGLSDataGlobalNumDof(gls_data)    ((gls_data) -> global_num_dof)
#define fsls_ParGLSDataLocalNumDof(gls_data)     ((gls_data) -> local_num_dof)
#define fsls_ParGLSDataDofPartition(gls_data)    ((gls_data) -> dof_partition)
#define fsls_ParGLSDataNumDofProcs(gls_data)     ((gls_data) -> num_dof_procs)
#define fsls_ParGLSDataSubAArray(gls_data)       ((gls_data) -> subA_array)
#define fsls_ParGLSDataSubBArray(gls_data)       ((gls_data) -> subb_array)
#define fsls_ParGLSDataSubXArray(gls_data)       ((gls_data) -> subx_array)


typedef struct
{
   JXF_Int    solver_id;      // solver id
   JXF_Real tol;            // the tolerance   
   JXF_Int    max_iter;       // the maximal number of iterations
   JXF_Int    min_iter;       // the minimal number of iterations  
   JXF_Int    two_norm;       // L^2 norm(1) or preconditioner-norm(0)? used for PCG
   JXF_Int    k_dim;          // the dimension of krylov subspace, used for PGMRES 
   JXF_Int    print_level;    // how much information to be printed

} fsls_SolverParam;

#define fsls_SolverParamSolverID(solver)    ((solver) -> solver_id)
#define fsls_SolverParamTolerance(solver)   ((solver) -> tol)
#define fsls_SolverParamMaxIter(solver)     ((solver) -> max_iter)
#define fsls_SolverParamMinIter(solver)     ((solver) -> min_iter)
#define fsls_SolverParamTwoNorm(solver)     ((solver) -> two_norm)
#define fsls_SolverParamKDim(solver)        ((solver) -> k_dim)
#define fsls_SolverParamPrintLevel(solver)  ((solver) -> print_level)


/*----------------------------------------------------------------*
 *                  Functions Declaration                         *
 *----------------------------------------------------------------*/ 

/* tool.c */
void fsls_ISwap(JXF_Int *w, JXF_Int i, JXF_Int j);
void fsls_DSwap(JXF_Real *w, JXF_Int i, JXF_Int j);
void fsls_IISwap(JXF_Int *v, JXF_Int *w, JXF_Int i, JXF_Int j);
void fsls_IDSwap(JXF_Int *v, JXF_Real *w, JXF_Int i, JXF_Int j);
void fsls_IQuickSort(JXF_Int order, JXF_Int *a, JXF_Int left, JXF_Int right);
void fsls_DQuickSort(JXF_Int order, JXF_Real *a, JXF_Int left, JXF_Int right);
void fsls_IQuickSortIndex(JXF_Int order, JXF_Int *a, JXF_Int left, JXF_Int right, JXF_Int *index);
void fsls_DQuickSortIndex(JXF_Int order, JXF_Real *a, JXF_Int left, JXF_Int right, JXF_Int *index);
void fsls_IIQuickSort(JXF_Int order, JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right);
void fsls_IDQuickSort(JXF_Int order, JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right);
void fsls_DIQuickSort(JXF_Int order, JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right);
void fsls_IDfabsQuickSort(JXF_Int order, JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right);
void fsls_iQuickSort12( JXF_Int *data, JXF_Int left, JXF_Int right );
void fsls_iQuickSortIndex12( JXF_Int *data, JXF_Int left, JXF_Int right, JXF_Int *index );
void fsls_dQuickSort12( JXF_Real *data, JXF_Int left, JXF_Int right );
void fsls_dQuickSortIndex12( JXF_Real *data, JXF_Int left, JXF_Int right, JXF_Int *index );
void fsls_iiQuickSort12( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void fsls_ddQuickSort12( JXF_Real *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_idQuickSort12( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_idFabsQuickSort12( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_diQuickSort12( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void fsls_iQuickSort21( JXF_Int *data, JXF_Int left, JXF_Int right );
void fsls_iQuickSortIndex21( JXF_Int *data, JXF_Int left, JXF_Int right, JXF_Int *index );
void fsls_dQuickSort21( JXF_Real *data, JXF_Int left, JXF_Int right );
void fsls_dQuickSortIndex21( JXF_Real *data, JXF_Int left, JXF_Int right, JXF_Int *index );
void fsls_iiQuickSort21( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void fsls_ddQuickSort21( JXF_Real *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_idQuickSort21( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_idFabsQuickSort21( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void fsls_diQuickSort21( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void fsls_Free( char *ptr );
char *fsls_CAlloc( size_t count, size_t elt_size );
JXF_Int  fsls_OutOfMemory( size_t size );
void fsls_SeedRand( JXF_Int seed );
JXF_Real fsls_Rand();
void fsls_ArrayInitialize( JXF_Real *x, JXF_Int n );                     
JXF_Int fsls_ArrayPrint( JXF_Real *u, JXF_Int n, char *filename );
JXF_Int fsls_IntArrayPrint( JXF_Int *u, JXF_Int n, char *filename );
JXF_Int fsls_IntArrayIJPrint( JXF_Int *u, JXF_Int n, char *filename );
void fsls_ArrayCopy( JXF_Int size, JXF_Real *x, JXF_Real *y );
JXF_Real fsls_Arrayl2Norm( JXF_Real *x, JXF_Int n );
JXF_Real fsls_ArrayAveL2Norm( JXF_Real *x, JXF_Int n );
JXF_Int fsls_ArrayDoubleAbsMax( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr );
JXF_Int fsls_ArrayDoubleAbsMin( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr );
JXF_Int fsls_ArrayDoubleMax( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr );
JXF_Int fsls_ArrayDoubleMin( JXF_Real *x, JXF_Int n, JXF_Real *elm_ptr );
JXF_Int fsls_ArrayIntAbsMax( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr );
JXF_Int fsls_ArrayIntAbsMin( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr );
JXF_Int fsls_ArrayIntMax( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr );
JXF_Int fsls_ArrayIntMin( JXF_Int *x, JXF_Int n, JXF_Int *elm_ptr );
JXF_Int fsls_daxpy( JXF_Int n, JXF_Real alpha, JXF_Real *x, JXF_Real *y );
JXF_Int fsls_gselim(JXF_Real *A, JXF_Real *x, JXF_Int n);
JXF_Int fsls_gselim_piv(JXF_Real *A, JXF_Real *x, JXF_Int n);

/* matvec.c */
fsls_DenseMatrix *fsls_DenseMatrixCreate( JXF_Int num_rows, JXF_Int num_cols );
JXF_Int fsls_DenseMatrixInitialize( fsls_DenseMatrix *A );
JXF_Int fsls_DenseMatrixDestroy( fsls_DenseMatrix *A );
JXF_Int fsls_DenseMatrixPrint( JXF_Int row, JXF_Int col, JXF_Real *full, char *filename );
void fsls_DenseymAx( JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y );
void fsls_DenseaAxpby( JXF_Real alpha, JXF_Real beta, JXF_Int n, JXF_Real *A, JXF_Real *x, JXF_Real *y );
void fsls_DenseMatMul( JXF_Real *a, JXF_Real *b, JXF_Int n, JXF_Real *c );
void fsls_DenseMatvec( JXF_Real *a, JXF_Real *b, JXF_Int n, JXF_Real *c );
JXF_Int fsls_DenseInverse( JXF_Int pivot_type, JXF_Real *a, JXF_Int n, JXF_Int *is, JXF_Int *js );
JXF_Int fsls_BuildVecFromFile( char *filename, fsls_Vector **b_ptr );
fsls_Vector *fsls_SeqVectorRead( char *file_name );
fsls_Vector *fsls_IJVectorRead( char *file_name );
JXF_Int fsls_SeqVectorPrint( fsls_Vector *vector, char *file_name );
fsls_Vector *fsls_SeqVectorCreate( JXF_Int size );
JXF_Int fsls_SeqVectorInitialize( fsls_Vector *vector );
JXF_Int fsls_SeqVectorSetConstantValues( fsls_Vector *x, JXF_Real value );
JXF_Int fsls_SeqVectorDestroy( fsls_Vector *vector );
JXF_Real fsls_SeqVectorInnerProd( fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_SeqVectorCopy( fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_SeqVectorSetRandomValues( fsls_Vector *v, JXF_Int seed );
JXF_Int fsls_SeqVectorScale( JXF_Real alpha, fsls_Vector *y );
JXF_Int fsls_SeqVectorScaleBack( fsls_Vector *f, JXF_Real *v );
fsls_Vector *fsls_SeqVectorAllOne( JXF_Int n );
fsls_Vector *fsls_SeqVectorOne2N( JXF_Int n );
JXF_Int fsls_SeqVectorAxpy( JXF_Real alpha, fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_SeqVectorsDiff( fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_SeqVectorFilesDiff( char *file1, char *file2 );
JXF_Real fsls_SeqVectorAveL2Norm( fsls_Vector *x );
JXF_Real fsls_SeqVectorL2Norm( fsls_Vector *x );
JXF_Real fsls_SeqVectorDiagInvNorm( JXF_Real *diag, fsls_Vector *x );
fsls_Vector *fsls_SeqVectorReorderByDirection( fsls_Vector *x );
JXF_Int fsls_SeqVectorReorderByPoint( fsls_Vector *x );
JXF_Int fsls_SeqVectorRecoverByPoint( fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_SeqVectorPartitionByRow( fsls_Vector *x, JXF_Int np, JXF_Int *partition, fsls_Vector ***x_array_ptr );
fsls_Vector *fsls_SeqVectorPermute( fsls_Vector *x, JXF_Int *p );
fsls_Vector *fsls_SeqVectorPermuteT( fsls_Vector *x, JXF_Int *p );
fsls_Vector *fsls_VecfromREItoEIR( fsls_Vector *x );
JXF_Int fsls_IVectorSetConstantValues( fsls_IVector *x, JXF_Int value );
fsls_IVector *fsls_IVectorRead( char *file_name );
JXF_Int fsls_IVectorPrint( fsls_IVector *vector, char *file_name );
JXF_Int fsls_IVectorDestroy( fsls_IVector *vector );
fsls_IVector *fsls_IVectorCreate( JXF_Int size );
JXF_Int fsls_IVectorInitialize( fsls_IVector *vector );
JXF_Int fsls_IVectorSetNatureOrder( fsls_IVector *x );
fsls_VecGroup *fsls_VecGroupCreate( JXF_Int n );
JXF_Int fsls_VecGroupInitialize( fsls_VecGroup *vecgroup );
fsls_Vector *fsls_VecGroupCombine( fsls_VecGroup *vecgroup );
JXF_Int fsls_BuildCSRMatFromFile( char *filename, fsls_CSRMatrix **A_ptr );
fsls_CSRMatrix *fsls_CSRMatrixRead( char *file_name );
fsls_CSRMatrix *fsls_CSRMatrixCreate( JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros );
JXF_Int fsls_CSRMatrixInitialize( fsls_CSRMatrix *matrix );
JXF_Int fsls_CSRMatrixReorder( fsls_CSRMatrix *A );
JXF_Int fsls_CSRMatrixPrint( fsls_CSRMatrix *matrix, char *file_name );
JXF_Int fsls_CSRMatrixDestroy( fsls_CSRMatrix *matrix );
JXF_Int fsls_CSRMatrixTranspose( fsls_CSRMatrix *A, fsls_CSRMatrix **AT, JXF_Int data );
JXF_Int fsls_CSRMatrixScaledNorm( fsls_CSRMatrix *A, JXF_Real *scnorm );
JXF_Int fsls_CSRMatrixCopy( fsls_CSRMatrix *A, fsls_CSRMatrix *B, JXF_Int copy_data );
fsls_CSRMatrix *fsls_CSRMatrixAdd( fsls_CSRMatrix *A, fsls_CSRMatrix *B );
fsls_CSRMatrix *fsls_CSRMatrixMultiply( fsls_CSRMatrix *A, fsls_CSRMatrix *B );
fsls_CSRMatrix *fsls_CSRMatrixTriMultiply( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_opt( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_01( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fsls_CSRMatrixTriMultiply_02( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fasp_CSRMatrixTriMultiply( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_improve( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_improve2( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fasp_CSRMatrixTriMultiply_mixed( fsls_CSRMatrix *R, fsls_CSRMatrix *A, fsls_CSRMatrix *P );
fsls_CSRMatrix *fsls_CSRMatrixDeleteZeros( fsls_CSRMatrix *A, JXF_Real tol );
fsls_CSRMatrix *fsls_CSRMatrixSymmetrization( fsls_CSRMatrix *A );
JXF_Int fsls_CSRMatrixScale( fsls_CSRMatrix *A, JXF_Real alpha );
JXF_Real fsls_CSRMatrixFrobeniusNorm( fsls_CSRMatrix *A );
JXF_Real fsls_CSRMatrixInfiniteNorm( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_GetSubCSRMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc );
fsls_DenseMatrix *fsls_SubDenseMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc );
JXF_Int fsls_GetSubDenseMatrix( fsls_CSRMatrix *A, JXF_Int *ROW, JXF_Int *COL, JXF_Int mr, JXF_Int mc, JXF_Real *dense, JXF_Int *G2LCol );
JXF_Int fsls_BuildCSRMatFromArrays( JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n, JXF_Int nz, fsls_CSRMatrix **A_ptr );
fsls_CSRMatrix *fsls_CSRMatrixReorderByPoint( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_CSRMatrixReorderByDirection( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_CSRMatrixReorderByDirection2( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_CSRMatrixReorderByDirectionRect( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_CSRMatrixGet3Block( fsls_CSRMatrix *A );
void fsls_CSRMatrixReorderColumnNumber( fsls_CSRMatrix *A );
void fsls_CSRMatrixReorderColumnNumber2( fsls_CSRMatrix *A );
void fsls_CSRMatrixReorderColumnNumber3( fsls_CSRMatrix *A );
JXF_Real fsls_CSRMatrixFindMinElm( fsls_CSRMatrix *A );
JXF_Real fsls_CSRMatrixFindMaxElm( fsls_CSRMatrix *A );
JXF_Real fsls_CSRMatrixFindMinFabsElm( fsls_CSRMatrix *A );
JXF_Real fsls_CSRMatrixFindMaxFabsElm( fsls_CSRMatrix *A );
JXF_Int fsls_CSRMatrixSparsityPatternCheck( JXF_Int col_reordered, fsls_CSRMatrix *A, fsls_CSRMatrix *B );
JXF_Int fsls_CSRMatrixDiffFromFiles( char *MatFile1, char *MatFile2 );
JXF_Int fsls_CSRMatrixDiff( fsls_CSRMatrix *A, fsls_CSRMatrix *B );
JXF_Int fsls_CSRMatrixRelativeDiff( fsls_CSRMatrix *A, fsls_CSRMatrix *B, char *filename );
JXF_Int fsls_CSRMatrixPartitionByRow0( fsls_CSRMatrix *A, JXF_Int np, JXF_Int *partition, fsls_CSRMatrix ***recA_array_ptr );
JXF_Int fsls_CSRMatrixPartitionByRow( fsls_CSRMatrix *A, JXF_Int np, JXF_Int *partition, fsls_CSRMatrix ***recA_array_ptr );
fsls_CSRMatrix *fsls_CSRMatrixPermuteTest( fsls_CSRMatrix *A, JXF_Int *p ); 
fsls_CSRMatrix *fsls_CSRMatrixPermute( fsls_CSRMatrix *A, JXF_Int *p );
JXF_Int *fsls_CSRMatrixNumNZPerRow( fsls_CSRMatrix *A );
JXF_Int *fsls_CSRMatrixReOrderByNZPerRow( fsls_CSRMatrix *A, JXF_Int order );
void fsls_CSRMatrixIndexC2Fortran( JXF_Int n, JXF_Int nz, JXF_Int *ia, JXF_Int *ja ); 
void fsls_CSRMatrixIndexFortran2C( JXF_Int n, JXF_Int nz, JXF_Int *ia, JXF_Int *ja ); 
fsls_CSRMatrix *fsls_CSRMatrixFilter( fsls_CSRMatrix *A, JXF_Real strong_threshold );
fsls_CSRMatrix *fsls_CSRMatrixIdentity( JXF_Int n );
JXF_Int fsls_CSRMatrixLowerPart( fsls_CSRMatrix *A, JXF_Int reorder, fsls_CSRMatrix **L_ptr );
JXF_Int fsls_CSRMatrixUpperPart( fsls_CSRMatrix *A, JXF_Int reorder, fsls_CSRMatrix **U_ptr );
JXF_Real *fsls_CSRMatrixDiagonalEntries( fsls_CSRMatrix *A );
fsls_CSRMatrix *fsls_MatfromREItoEIR( fsls_CSRMatrix *A );
fsls_ICSRMatrix *fsls_ICSRMatrixCreate( JXF_Int num_rows, JXF_Int num_cols, JXF_Int num_nonzeros );              
JXF_Int fsls_ICSRMatrixInitialize( fsls_ICSRMatrix *matrix );
JXF_Int fsls_ICSRMatrixDestroy( fsls_ICSRMatrix *matrix );                
fsls_CSRGroup *fsls_CSRGroupCreate( JXF_Int n );
JXF_Int fsls_CSRGroupInitialize( fsls_CSRGroup *csrgroup );
fsls_CSRMatrix *fsls_CSRGroupDiagCombine( fsls_CSRGroup *csrgroup );
JXF_Int fsls_CSRMatrixMatvec00( fsls_CSRMatrix *A, fsls_Vector *x, fsls_Vector *y );
JXF_Int fsls_CSRMatrixMatvec01( JXF_Real alpha, fsls_CSRMatrix *A, fsls_Vector *x, JXF_Real beta, fsls_Vector *y );
JXF_Int fsls_CSRMatrixMatvec02( JXF_Real alpha, fsls_CSRMatrix *A, JXF_Real *x, JXF_Real beta, JXF_Real *y );
JXF_Int fsls_CSRMatrixMatvec( JXF_Real alpha, fsls_CSRMatrix *A, fsls_Vector *x, JXF_Real beta, fsls_Vector *y );
JXF_Int fsls_CSRMatrixMatvecT( JXF_Real alpha, fsls_CSRMatrix *A, fsls_Vector *x, JXF_Real beta, fsls_Vector *y );
JXF_Real fsls_CSRMatrixVecTMatVec( fsls_Vector *y, fsls_CSRMatrix *A, fsls_Vector *x );
JXF_Int fsls_ComputeResidual( JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n, JXF_Real *x, JXF_Real *b, JXF_Real *y );
JXF_Int fsls_MatVecMultiply( JXF_Int *ia, JXF_Int *ja, JXF_Real *a, JXF_Int n, JXF_Real *x, JXF_Real *y );
fsls_Vector *fsls_CreateRhsByMatrix( fsls_CSRMatrix *A );
JXF_Int 
fsls_CSRMatRhsRead( char    *MatFile,
                    char    *RhsFile,
                    JXF_Int    **ia_ptr,
                    JXF_Int    **ja_ptr,
                    JXF_Real **a_ptr,
                    JXF_Int     *n_ptr,
                    JXF_Int     *nz_ptr,
                    JXF_Real **f_ptr );
JXF_Int
fsls_Arrays2CSRMatVec( JXF_Int    *ia,
                       JXF_Int    *ja,
                       JXF_Real *a,
                       JXF_Int     n,
                       JXF_Int     nz,
                       JXF_Real *f,
                       JXF_Real *u,
                       fsls_CSRMatrix **A_ptr,
                       fsls_Vector    **b_ptr,
                       fsls_Vector    **x_ptr );
JXF_Int
fsls_Arrays2CSRMatVecNoAllocate( JXF_Int    *ia,
                                 JXF_Int    *ja,
                                 JXF_Real *a,
                                 JXF_Int     n,
                                 JXF_Int     nz,
                                 JXF_Real *f,
                                 JXF_Real *u,
                                 fsls_CSRMatrix **A_ptr,
                                 fsls_Vector    **b_ptr,
                                 fsls_Vector    **x_ptr );                             
JXF_Int
fsls_LinearSystemSymDiagScale0( fsls_CSRMatrix   *A,
                                fsls_Vector      *f,
                                fsls_CSRMatrix  **AA_ptr,
                                fsls_Vector     **ff_ptr,
                                JXF_Real          **diag_ptr );
JXF_Int
fsls_LinearSystemSymDiagScale( fsls_CSRMatrix   *A,
                               fsls_Vector      *f,
                               JXF_Real          **diag_ptr );
JXF_Int
fsls_LinearSystemUnSymDiagScale0( fsls_CSRMatrix   *A,
                                  fsls_Vector      *f,
                                  fsls_CSRMatrix  **AA_ptr,
                                  fsls_Vector     **ff_ptr );
JXF_Int
fsls_LinearSystemUnSymDiagScale( fsls_CSRMatrix   *A,
                                 fsls_Vector      *f  );                                          
JXF_Int fsls_CSR2FullMatrix( fsls_CSRMatrix *A, JXF_Real **full_ptr );
JXF_Int fsls_Full2CSRMatrix( JXF_Int n, JXF_Real *full, fsls_CSRMatrix **A_ptr );

/* ilu.c */ 
JXF_Int fsls_ilu1_FillInPredict( fsls_CSRMatrix *A );
JXF_Int fsls_ilu2_FillInPredict( fsls_CSRMatrix *A );
JXF_Int fsls_ilup_FillInPredict( fsls_CSRMatrix *A, JXF_Int p );
JXF_Int fsls_ILUp_FillInPredict( fsls_CSRMatrix *A, JXF_Int p );
JXF_Int
fsls_ILUp_Decomp( fsls_CSRMatrix  *A, 
                  JXF_Int              p,
                  JXF_Int            **indexLU_ptr,
                  JXF_Real         **valueLU_ptr );
JXF_Int
fsls_ILUp_DecompTest( fsls_CSRMatrix   *A, 
                      JXF_Int               p,
                      fsls_CSRMatrix  **L_ptr,
                      fsls_CSRMatrix  **U_ptr  );        

/* smoother.c */   
JXF_Int 
fsls_AMGRelaxJacobi( fsls_CSRMatrix  *A,
                     fsls_Vector     *f,
                     JXF_Int             *cf_marker,
                     JXF_Int              relax_order,
                     JXF_Real           relax_weight,
                     fsls_Vector     *u,
                     fsls_Vector     *Vtemp );                   
JXF_Int 
fsls_AMGRelaxGS( fsls_CSRMatrix  *A,
                 fsls_Vector     *f,
                 JXF_Int             *cf_marker,
                 JXF_Int              relax_order,
                 fsls_Vector     *u );
JXF_Int 
fsls_AMGRelaxGS_Partial( fsls_CSRMatrix  *A,
                         fsls_Vector     *f,
                         JXF_Int             *cf_marker,
                         JXF_Int              relax_order,
                         JXF_Int              begin,
                         JXF_Int              end,
                         fsls_Vector     *u );                                      
JXF_Int 
fsls_AMGRelaxSOR( fsls_CSRMatrix  *A,
                  fsls_Vector     *f,
                  JXF_Int             *cf_marker,
                  JXF_Int              relax_order,
                  JXF_Real           relax_weight,
                  fsls_Vector     *u  );
void fsls_RowNormCompute( fsls_CSRMatrix *A, fsls_Vector *x );
JXF_Int 
fsls_AMGRelaxKaczmarz( fsls_CSRMatrix  *A,
                       fsls_Vector     *f,
                       JXF_Int             *cf_marker,
                       JXF_Int              relax_order,
                       JXF_Real           relax_weight,
                       fsls_Vector     *u,
                       fsls_Vector     *Vtemp );
JXF_Int 
fsls_AMGRelaxPolynomial( JXF_Int              degree,
                         fsls_CSRMatrix  *A,
                         fsls_Vector     *f,
                         fsls_Vector     *u,
                         fsls_Vector     *R,
                         fsls_Vector     *V0,
                         fsls_Vector     *V,
                         fsls_Vector     *RAV );
JXF_Int  
fsls_AMGRelaxILUp( JXF_Int                 level,
                   fsls_RelaxILUData  *relax_ilu_data,
                   fsls_CSRMatrix     *A,
                   fsls_Vector        *f,
                   fsls_Vector        *u,
                   fsls_Vector        *r );            
JXF_Int fsls_AMGRelaxGE( JXF_Real *full, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_AMGRelaxGEP( JXF_Real *full, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u ); 

/* amg.c */
void 
fsls_dispose_elt(fsls_LinkList element_ptr);
void 
fsls_remove_point(fsls_LinkList  *LoL_head_ptr, 
                  fsls_LinkList  *LoL_tail_ptr, 
                  JXF_Int             measure,
                  JXF_Int             index, 
                  JXF_Int            *lists, 
                  JXF_Int            *where);
fsls_LinkList fsls_create_elt(JXF_Int Item);
void 
fsls_enter_on_lists(fsls_LinkList  *LoL_head_ptr, 
                    fsls_LinkList  *LoL_tail_ptr, 
                    JXF_Int             measure,
                    JXF_Int             index, 
                    JXF_Int            *lists, 
                    JXF_Int            *where);
void *fsls_AMGInitialize();
void fsls_AMGFinalize( void *data );
void fsls_AMGComplx( void *data, JXF_Real *grid, JXF_Real *operater );
JXF_Int fsls_AMGSetMaxLevels( void *data, JXF_Int max_levels );
JXF_Int fsls_AMGSetStrongThreshold( void *data, JXF_Real strong_threshold );
JXF_Int fsls_AMGSetMaxRowSum( void *data, JXF_Real max_row_sum );
JXF_Int fsls_AMGSetATruncFactor( void *data, JXF_Real A_trunc_factor );
JXF_Int fsls_AMGSetPTruncFactor( void *data, JXF_Real P_trunc_factor );
JXF_Int fsls_AMGSetAMaxElmts( void *data, JXF_Int A_max_elmts );
JXF_Int fsls_AMGSetPMaxElmts( void *data, JXF_Int P_max_elmts );
JXF_Int fsls_AMGSetCoarsenType( void *data, JXF_Int coarsen_type );
JXF_Int fsls_AMGSetInterpType( void *data, JXF_Int interp_type );
JXF_Int fsls_AMGSetCoarseThreshold( void *data, JXF_Int coarse_threshold );
JXF_Int fsls_AMGSetSMode( void *data, JXF_Int S_mode );
JXF_Int fsls_AMGSetTol( void *data, JXF_Real tol );
JXF_Int fsls_AMGSetMaxIter( void *data, JXF_Int max_iter );
JXF_Int fsls_AMGSetMinIter( void *data, JXF_Int min_iter );
JXF_Int fsls_AMGSetCycleType( void *data, JXF_Int cycle_type );
JXF_Int fsls_AMGSetCorrectiveType( void *data, JXF_Int corrective_type );
JXF_Int fsls_AMGSetRelaxWeight( void  *data, JXF_Real relax_weight );
JXF_Int fsls_AMGSetRelaxType( void  *data, JXF_Int relax_type );
JXF_Int fsls_AMGSetRelaxOrder( void  *data, JXF_Int relax_order );
JXF_Int fsls_AMGSetRelaxSweeps( void  *data, JXF_Int relax_sweeps );
JXF_Int fsls_AMGSetRelaxSym( void  *data, JXF_Int relax_sym );
JXF_Int fsls_AMGSetPolyDegree( void  *data, JXF_Int poly_degree );
JXF_Int fsls_AMGSetCoarsestSolver( void  *data, JXF_Int coarsest_solver );
JXF_Int fsls_AMGSetLevelOfFillIn( void  *data, JXF_Int level_of_fill_in );
JXF_Int fsls_AMGSetILUSmoothLevel( void  *data, JXF_Int ilu_smooth_level );
JXF_Int fsls_AMGSetRhsNrmThreshold( void  *data, JXF_Real rhsnorm_threshold );
JXF_Int fsls_AMGSetPrintLevel( void *data, JXF_Int print_level );
JXF_Int fsls_AMGSetCycleCount( void *data, JXF_Int cycle_count );
JXF_Int fsls_AMGGetCycleCount( void *data, JXF_Int *cycle_count );
JXF_Int fsls_AMGSetLastRelNrm( void *data, JXF_Real last_rel_nrm );
JXF_Int fsls_AMGGetLastRelNrm( void *data, JXF_Real *last_rel_nrm );
JXF_Int fsls_AMGSetAveConvFactor( void *data, JXF_Real ave_conv_factor );
JXF_Int fsls_AMGGetAveConvFactor( void *data, JXF_Real *ave_conv_factor );
JXF_Int fsls_AMGSetGridRelaxSweeps( void *data, JXF_Int *grid_relax_sweeps );
JXF_Int fsls_AMGSetGridRelaxType( void *data, JXF_Int *grid_relax_type );
JXF_Int fsls_AMGSetGridRelaxOrder( void *data, JXF_Int *grid_relax_order );
JXF_Int fsls_AMGSetGridRelaxWeight( void *data, JXF_Real *grid_relax_weight );
JXF_Int fsls_AMGCreateS( fsls_CSRMatrix *A, JXF_Real strength_threshold, JXF_Int mode, fsls_ICSRMatrix **S_ptr );
JXF_Int fsls_AMGCreateSNew( fsls_CSRMatrix *A, JXF_Real strength_threshold, 
                        JXF_Real max_row_sum, JXF_Int mode, fsls_ICSRMatrix **S_ptr );
JXF_Int
fsls_AMGCoarsenRS( fsls_CSRMatrix    *A,
                   fsls_ICSRMatrix   *S,
                   JXF_Int              **CF_marker_ptr,
                   JXF_Int               *coarse_size_ptr );
JXF_Int
fsls_AMGCoarsenRSNew( fsls_CSRMatrix    *A,
                      fsls_ICSRMatrix   *S,
                      JXF_Int              **CF_marker_ptr,
                      JXF_Int               *coarse_size_ptr );                   
JXF_Int
fsls_AMGCoarsenRugeLoL( fsls_CSRMatrix    *A,
                        fsls_ICSRMatrix   *S,
                        JXF_Int              **CF_marker_ptr,
                        JXF_Int               *coarse_size_ptr );
JXF_Int
fsls_AMGCoarsenCLJP( fsls_CSRMatrix   *A,
                     fsls_ICSRMatrix  *S,
                     JXF_Int             **CF_marker_ptr,
                     JXF_Int              *coarse_size_ptr );
JXF_Int
fsls_AMGCoarsenCLJPNew( fsls_CSRMatrix   *A,
                        fsls_ICSRMatrix  *S,
                        JXF_Int             **CF_marker_ptr,
                        JXF_Int              *coarse_size_ptr );                       
JXF_Int fsls_InitAMGIndepSet( fsls_ICSRMatrix *S, JXF_Real *measure_array, JXF_Real cconst );
JXF_Int
fsls_AMGIndepSet( fsls_ICSRMatrix *S,
                  JXF_Real          *measure_array,
                  JXF_Real           cconst,
                  JXF_Int             *graph_array,
                  JXF_Int              graph_array_size,
                  JXF_Int             *IS_marker );
JXF_Int fsls_AMGClassicalInterp( fsls_CSRMatrix *A, JXF_Int *CF_marker, fsls_ICSRMatrix *S, fsls_CSRMatrix **P_ptr );
JXF_Int fsls_AMGDirectInterp( fsls_CSRMatrix *A, JXF_Int *CF_marker, fsls_ICSRMatrix *S, fsls_CSRMatrix **P_ptr );
JXF_Int fsls_AMGPTruncation( fsls_CSRMatrix *A, JXF_Real trunc_factor, JXF_Int max_elmts );
JXF_Int 
fsls_AMGBuildCoarseOperator( fsls_CSRMatrix   *RT,
                             fsls_CSRMatrix   *A,
                             fsls_CSRMatrix   *P,
                             fsls_CSRMatrix  **RAP_ptr );
JXF_Int fsls_AMGAHTruncation( fsls_CSRMatrix *A, JXF_Real trunc_factor, JXF_Int max_elmts ); 
JXF_Int fsls_AMGSetup( void *amg_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_AMGSetupTT( void *amg_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_AMGSetupParams( void *amg_vdata );
JXF_Int fsls_AMGSetRelaxParams( void *data, JXF_Int use_Jacobi_on_coarsest );
JXF_Int fsls_AMGSolve( void *amg_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_AMGPrecond( void *amg_vdata, fsls_Vector *r, fsls_Vector *z ); // for preconditioner
JXF_Int fsls_AMGSolve_IterativeType( void *amg_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_AMGSolve_CorrectiveType( void *amg_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
void fsls_AMGSolveParams( void *data );
JXF_Int fsls_AMGCycle( void *amg_vdata, fsls_Vector **F_array, fsls_Vector **U_array ); 

/* krylov.c */
JXF_Int fsls_CSRPCGKrylovIdentitySetup( void *vdata, void *A );
JXF_Int fsls_CSRPCGKrylovIdentity( void *vdata, void *b, void *x );
fsls_CSRPCGFunctions *
fsls_CSRPCGFunctionsCreate( JXF_Int  (*PrecondSetup)  ( void *vdata, void *A ),
                            JXF_Int  (*Precond)       ( void *vdata, void *b, void *x ) );
fsls_CSRPCGData *fsls_CSRPCGCreate();
JXF_Int fsls_CSRPCGSetTol( void *pcg_vdata, JXF_Real tol );
JXF_Int fsls_CSRPCGSetAbsoluteTol( void *pcg_vdata, JXF_Real a_tol );
JXF_Int fsls_CSRPCGSetAbsoluteTolFactor( void *pcg_vdata, JXF_Real atolf );
JXF_Int fsls_CSRPCGSetConvergenceFactorTol( void *pcg_vdata, JXF_Real cf_tol );
JXF_Int fsls_CSRPCGSetMaxIter( void *pcg_vdata, JXF_Int max_iter );
JXF_Int fsls_CSRPCGSetTwoNorm( void *pcg_vdata, JXF_Int two_norm );
JXF_Int fsls_CSRPCGSetRelChange( void *pcg_vdata, JXF_Int rel_change );
JXF_Int fsls_CSRPCGSetRecomputeResidual( void *pcg_vdata, JXF_Int recompute_residual );
JXF_Int fsls_CSRPCGSetStopCrit( void *pcg_vdata, JXF_Int stop_crit );
JXF_Int fsls_CSRPCGSetConverged( void *pcg_vdata, JXF_Int converged );
JXF_Int fsls_CSRPCGSetPrintLevel( void *pcg_vdata, JXF_Int level );
JXF_Int fsls_CSRPCGGetConverged( void *pcg_vdata, JXF_Int *converged );
JXF_Int fsls_CSRPCGGetNumIterations( void *pcg_vdata, JXF_Int *num_iterations );
JXF_Int fsls_CSRPCGGetFinalRelativeResidualNorm( void *pcg_vdata, JXF_Real *relative_residual_norm );
JXF_Int fsls_CSRPCGSetPrecond( void *pcg_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int fsls_CSRPCGSetup( void *pcg_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPCGSolve( void *pcg_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPCGDestroy( void *pcg_vdata );
JXF_Int fsls_CSRPGMRESKrylovIdentitySetup( void *vdata, void *A );
JXF_Int fsls_CSRPGMRESKrylovIdentity( void *vdata, void *b, void *x );
fsls_CSRPGMRESFunctions *
fsls_CSRPGMRESFunctionsCreate( JXF_Int  (*PrecondSetup)  ( void *vdata, void *A ),
                               JXF_Int  (*Precond)       ( void *vdata, void *b, void *x ) );
fsls_CSRPGMRESData *fsls_CSRPGMRESCreate();
JXF_Int fsls_CSRPGMRESSetKDim( void *gmres_vdata, JXF_Int k_dim );
JXF_Int fsls_CSRPGMRESSetTol( void *gmres_vdata, JXF_Real tol );
JXF_Int fsls_CSRPGMRESSetAbsoluteTol( void *gmres_vdata, JXF_Real a_tol );
JXF_Int fsls_CSRPGMRESSetConvergenceFactorTol( void *gmres_vdata, JXF_Real cf_tol );
JXF_Int fsls_CSRPGMRESSetMinIter( void *gmres_vdata, JXF_Int min_iter );
JXF_Int fsls_CSRPGMRESSetMaxIter( void *gmres_vdata, JXF_Int max_iter );
JXF_Int fsls_CSRPGMRESSetRelChange( void *gmres_vdata, JXF_Int rel_change );
JXF_Int fsls_CSRPGMRESSetStopCrit( void *gmres_vdata, JXF_Int stop_crit );
JXF_Int fsls_CSRPGMRESSetPrintLevel( void *gmres_vdata, JXF_Int level );
JXF_Int fsls_CSRPGMRESGetNumIterations( void *gmres_vdata, JXF_Int *num_iterations );
JXF_Int fsls_CSRPGMRESGetConverged( void *gmres_vdata, JXF_Int  *converged );
JXF_Int fsls_CSRPGMRESGetFinalRelativeResidualNorm( void *gmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int fsls_CSRPGMRESSetPrecond( void *gmres_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int fsls_CSRPGMRESSetup( void *pgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPGMRESSolve( void *pgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPGMRESDestroy( void *pgmres_vdata );
JXF_Int fsls_CSRPLGMRESKrylovIdentitySetup( void *vdata, void *A );
JXF_Int fsls_CSRPLGMRESKrylovIdentity( void *vdata, void *b, void *x );
fsls_CSRPLGMRESFunctions *
fsls_CSRPLGMRESFunctionsCreate( JXF_Int  (*PrecondSetup) ( void *vdata, void *A ),
                                JXF_Int  (*Precond)      ( void *vdata, void *b, void *x ) );
fsls_CSRPLGMRESData *fsls_CSRPLGMRESCreate();
JXF_Int fsls_CSRPLGMRESSetKDim( void *plgmres_vdata, JXF_Int k_dim );
JXF_Int fsls_CSRPLGMRESSetAugDim( void *plgmres_vdata, JXF_Int aug_dim );
JXF_Int fsls_CSRPLGMRESSetTol( void *plgmres_vdata, JXF_Real tol );
JXF_Int fsls_CSRPLGMRESSetAbsoluteTol( void *plgmres_vdata, JXF_Real a_tol );
JXF_Int fsls_CSRPLGMRESSetConvergenceFactorTol( void *plgmres_vdata, JXF_Real cf_tol );
JXF_Int fsls_CSRPLGMRESSetMinIter( void *plgmres_vdata, JXF_Int min_iter );
JXF_Int fsls_CSRPLGMRESSetMaxIter( void *plgmres_vdata, JXF_Int max_iter );
JXF_Int fsls_CSRPLGMRESSetStopCrit( void *plgmres_vdata, JXF_Int stop_crit );
JXF_Int fsls_CSRPLGMRESSetPrintLevel( void *plgmres_vdata, JXF_Int level );
JXF_Int fsls_CSRPLGMRESGetNumIterations( void *plgmres_vdata, JXF_Int *num_iterations );
JXF_Int fsls_CSRPLGMRESGetConverged( void *plgmres_vdata, JXF_Int  *converged );
JXF_Int fsls_CSRPLGMRESGetFinalRelativeResidualNorm( void *plgmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int fsls_CSRPLGMRESSetPrecond( void *plgmres_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int fsls_CSRPLGMRESSetup( void *plgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPLGMRESSolve( void *plgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPLGMRESDestroy( void *plgmres_vdata );
JXF_Int fsls_CSRPAGMRESKrylovIdentitySetup( void *vdata, void *A );
JXF_Int fsls_CSRPAGMRESKrylovIdentity( void *vdata, void *b, void *x );
fsls_CSRPAGMRESFunctions *
fsls_CSRPAGMRESFunctionsCreate( JXF_Int  (*PrecondSetup)  ( void *vdata, void *A ),
                                JXF_Int  (*Precond)       ( void *vdata, void *b, void *x ) );
fsls_CSRPAGMRESData *fsls_CSRPAGMRESCreate();
JXF_Int fsls_CSRPAGMRESSetKDimMax( void *pagmres_vdata, JXF_Int k_dim_max );
JXF_Int fsls_CSRPAGMRESSetKDimMin( void *pagmres_vdata, JXF_Int k_dim_min );
JXF_Int fsls_CSRPAGMRESSetKDimD( void *pagmres_vdata, JXF_Int k_dim_d );
JXF_Int fsls_CSRPAGMRESSetCRMax( void *pagmres_vdata, JXF_Real cr_max );
JXF_Int fsls_CSRPAGMRESSetCRMin( void *pagmres_vdata, JXF_Real cr_min );
JXF_Int fsls_CSRPAGMRESSetTol( void *pagmres_vdata, JXF_Real tol );
JXF_Int fsls_CSRPAGMRESSetConvergenceFactorTol( void *pagmres_vdata, JXF_Real cf_tol );
JXF_Int fsls_CSRPAGMRESSetMinIter( void *pagmres_vdata, JXF_Int min_iter );
JXF_Int fsls_CSRPAGMRESSetMaxIter( void *pagmres_vdata, JXF_Int max_iter );
JXF_Int fsls_CSRPAGMRESSetStopCrit( void *pagmres_vdata, JXF_Int stop_crit );
JXF_Int fsls_CSRPAGMRESSetPrintLevel( void *pagmres_vdata, JXF_Int level );
JXF_Int fsls_CSRPAGMRESGetNumIterations( void *pagmres_vdata, JXF_Int *num_iterations );
JXF_Int fsls_CSRPAGMRESGetConverged( void *pagmres_vdata, JXF_Int  *converged );
JXF_Int fsls_CSRPAGMRESGetFinalRelativeResidualNorm( void *pagmres_vdata, JXF_Real *relative_residual_norm );
JXF_Int fsls_CSRPAGMRESSetPrecond( void *pagmres_vdata, JXF_Int (*precond)(), JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int fsls_CSRPAGMRESSetup( void *pgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPAGMRESSolve( void *pgmres_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPAGMRESDestroy( void *pgmres_vdata );
JXF_Int fsls_CSRPBICGSTABKrylovIdentitySetup( void *vdata, void *A );
JXF_Int fsls_CSRPBICGSTABKrylovIdentity( void *vdata, void *b, void *x );
fsls_CSRPBICGSTABFunctions *
fsls_CSRPBICGSTABFunctionsCreate( JXF_Int  (*PrecondSetup)  ( void *vdata, void *A ),
                                  JXF_Int  (*Precond)       ( void *vdata, void *b, void *x ) );
fsls_CSRPBICGSTABData *fsls_CSRPBICGSTABCreate();
JXF_Int fsls_CSRPBICGSTABSetTol( void *pbicgstab_vdata, JXF_Real tol );
JXF_Int fsls_CSRPBICGSTABSetAbsoluteTol( void *pbicgstab_vdata, JXF_Real a_tol );
JXF_Int fsls_CSRPBICGSTABSetConvergenceFactorTol( void *pbicgstab_vdata, JXF_Real cf_tol );
JXF_Int fsls_CSRPBICGSTABSetMinIter( void *pbicgstab_vdata, JXF_Int min_iter );
JXF_Int fsls_CSRPBICGSTABSetMaxIter( void *pbicgstab_vdata, JXF_Int max_iter );
JXF_Int fsls_CSRPBICGSTABSetStopCrit( void *pbicgstab_vdata, JXF_Int stop_crit );
JXF_Int fsls_CSRPBICGSTABSetPrintLevel( void *pbicgstab_vdata, JXF_Int print_level );
JXF_Int fsls_CSRPBICGSTABGetConverged( void *pbicgstab_vdata, JXF_Int *converged );
JXF_Int fsls_CSRPBICGSTABGetNumIterations( void *pbicgstab_vdata, JXF_Int *num_iterations );
JXF_Int fsls_CSRPBICGSTABGetFinalRelativeResidualNorm( void *pbicgstab_vdata, JXF_Real *relative_residual_norm );
JXF_Int fsls_CSRPBICGSTABSetPrecond( void *pbicgstab_vdata, JXF_Int (*precond)(), 
                                 JXF_Int (*precond_setup)(), void *precond_data );
JXF_Int fsls_CSRPBICGSTABSetup( void *pbicgstab_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPBICGSTABSolve( void *pbicgstab_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_CSRPBICGSTABDestroy( void *pbicgstab_vdata );

/* precond.c */
void *fsls_PreILUpDataInitialize( JXF_Int level_fill_in );
void fsls_PreILUpDataFinalize( void *ilup_vdata );
JXF_Int fsls_PreILUpSetup( void *pre_ilup_vdata, fsls_CSRMatrix *A, fsls_Vector *b, fsls_Vector *x );
JXF_Int fsls_PreILUpSolve( void *pre_ilup_vdata, fsls_Vector *r, fsls_Vector *z );
void *fsls_PreDSDataInitialize();
JXF_Int fsls_PreDSSetup( void *vdata, void *A, void *b, void *x );
JXF_Int fsls_PreDSSolve( void *vdata, void *b, void *x );
void *fsls_PreDiagDataInitialize( JXF_Int m1, JXF_Int m2, JXF_Int level_fill_in );
void fsls_PreDiagDataFinalize( void *pre_diag_vdata );
JXF_Int fsls_PreDiagSetup( void *pre_diag_vdata, fsls_CSRMatrix *A, fsls_Vector *f, fsls_Vector *u );
JXF_Int fsls_PreDiagSolve( void *pre_diag_vdata, fsls_Vector *r, fsls_Vector *z );

/* /csrc/multils/par_lsdata.c */
fsls_ParLSData *
fsls_ParLSDataCreate( MPI_Comm  comm,
                      JXF_Int       num_rows, 
                      JXF_Int       num_cols, 
                      JXF_Int       num_nonzeros,
                      JXF_Int       global_num_ls, 
                      JXF_Int      *ls_partition,
                      JXF_Int      *num_ls_procs );
JXF_Int fsls_ParLSDataInitialize( fsls_ParLSData *ls_data );
void fsls_ParLSDataDestroy( fsls_ParLSData *ls_data );
fsls_ParLSData *
fsls_GenerateLSData( MPI_Comm comm, JXF_Int global_num_ls, char **MatFile_array, char **RhsFile_array );

/* /csrc/multils/par_glsdata.c */
fsls_ParGLSData *
fsls_ParGLSDataCreate( MPI_Comm  comm, 
                       JXF_Int       global_num_ls, 
                       JXF_Int       global_num_dof, 
                       JXF_Int      *dof_partition,
                       JXF_Int      *num_dof_procs );
JXF_Int fsls_ParGLSDataInitialize( fsls_ParGLSData *gls_data );
void fsls_ParGLSDataDestroy( fsls_ParGLSData *gls_data );
fsls_ParGLSData *
fsls_ConstructGLSData( MPI_Comm  comm, 
                       JXF_Int       global_num_ls, 
                       JXF_Int       global_num_dof, 
                       char    **MatFile, 
                       char    **RhsFile );
fsls_ParGLSData *
fsls_GenerateGLSData0( MPI_Comm   comm, 
                       JXF_Int        global_num_ls, 
                       char      *MatFile, 
                       char      *RhsFile );
fsls_ParGLSData *
fsls_GenerateGLSData1( MPI_Comm   comm, 
                       JXF_Int        global_num_ls,
                       JXF_Int        global_num_dof, 
                       char     **MatFile_array, 
                       char     **RhsFile_array );
fsls_CSRMatrix  * 
fsls_GetSampleMatrix( MPI_Comm         comm, 
                      fsls_CSRMatrix  *A, 
                      JXF_Int             *partition );
fsls_Vector *
fsls_GetSampleVector( MPI_Comm      comm, 
                      fsls_Vector  *x, 
                      JXF_Int          *partition );

/* /csrc/multils/par_gls2ls.c */
fsls_ParLSData *fsls_ParGLS2LSData( fsls_ParGLSData *gls_data );

/* /csrc/multils/par_mllssolve.c */
fsls_SolverParam *fsls_SolverParamCreate();
JXF_Int fsls_SolverParamSetSolverID( void *vdata, JXF_Int solver_id );
JXF_Int fsls_SolverParamSetTolerance( void *vdata, JXF_Real tol );
JXF_Int fsls_SolverParamSetMaxIter( void *vdata, JXF_Int max_iter );
JXF_Int fsls_SolverParamSetMinIter( void *vdata, JXF_Int min_iter );
JXF_Int fsls_SolverParamSetTwoNorm( void *vdata, JXF_Int two_norm );
JXF_Int fsls_SolverParamSetKDim( void *vdata, JXF_Int k_dim );
JXF_Int fsls_SolverParamSetPrintLevel( void *vdata, JXF_Int print_level );
JXF_Int fsls_ParallelSolve( fsls_ParLSData *ls_data, fsls_SolverParam *solver );

/* /csrc/multils/par_util.c */
void fsls_GetWallTime( MPI_Comm comm, char *fctname, JXF_Real starttime, JXF_Real endtime, JXF_Int allproc, JXF_Int dp );
void fsls_MPICommInformation( MPI_Comm comm, JXF_Int *myid_ptr, JXF_Int *nprocs_ptr );
JXF_Int fsls_BalancedPartition( JXF_Int N, JXF_Int np, JXF_Int **partition_ptr, JXF_Int **each_ptr );

#endif                    
