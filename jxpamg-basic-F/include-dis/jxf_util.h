//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_util.h -- head files for basic tools
 *  Date: 2011/09/08
 *  Modified Date: 2012/10/22
 *
 *  Created by peghoty
 *  Modified by Yue Xiaoqiang
 */ 

#ifndef JXF_UTIL_HEADER
#define JXF_UTIL_HEADER

#if JXF_USING_MPI
#define TIMER_USE_MPI
#else
#define JXF_SEQUENTIAL
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/mman.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include <float.h>
#ifndef JXF_SEQUENTIAL
#include <mpi.h>
#endif
#if JXF_USING_OPENMP
#include <omp.h>
#endif

/*----------------------------------------------------------------*
 *                      Macro Definition                          *
 *----------------------------------------------------------------*/ 

/* BigInt */
typedef long int jxf_longint;
#if JXF_USING_BIG_INT
typedef long long int JXF_Int;
typedef unsigned long long int JXF_UInt;
#define JXF_MPI_INT MPI_LONG_LONG_INT
#else /* default */
typedef int JXF_Int;
typedef unsigned int JXF_UInt;
#define JXF_MPI_INT MPI_INT
#endif

/* Bool */
typedef JXF_Int JXF_Bool;
#define JXF_TRUE 1
#define JXF_FALSE 0

/* BigDouble */
#if JXF_USING_BIG_DOUBLE
typedef long double JXF_Real;
#define JXF_REAL_EPSILON LDBL_EPSILON
#define JXF_MPI_REAL MPI_LONG_DOUBLE
#else
typedef float JXF_Real;
#define JXF_REAL_EPSILON FLT_EPSILON
#define JXF_MPI_REAL MPI_FLOAT
#endif

/* Complex */
#if JXF_USING_COMPLEX
typedef double _Complex JXF_Complex;
#define JXF_MPI_COMPLEX MPI_C_DOUBLE_COMPLEX
JXF_Complex jxf_conj( JXF_Complex value );
JXF_Real jxf_cabs( JXF_Complex value );
JXF_Real jxf_creal( JXF_Complex value );
JXF_Real jxf_cimag( JXF_Complex value );
#else /* default */
typedef JXF_Real JXF_Complex;
#define JXF_MPI_COMPLEX JXF_MPI_REAL
#define jxf_conj(value)  value
#define jxf_cabs(value)  fabsf(value)
#define jxf_creal(value) value
#define jxf_cimag(value) 0.0f
#endif

/* Sequential MPI stuff */
#ifdef JXF_SEQUENTIAL
typedef JXF_Int MPI_Comm;
typedef JXF_Int MPI_Group;
typedef JXF_Int MPI_Request;
typedef JXF_Int MPI_Datatype;
typedef void (MPI_User_function) ();
typedef struct
{
   JXF_Int MPI_SOURCE;
   JXF_Int MPI_TAG;
} MPI_Status;
typedef JXF_Int MPI_Op;
typedef JXF_Int MPI_Aint;
#define MPI_COMM_SELF 1
#define MPI_COMM_WORLD 0
#define MPI_COMM_NULL -1
#define MPI_COMM_TYPE_SHARED 0
#define MPI_BOTTOM 0x0
#define MPI_FLOAT 0
#define MPI_DOUBLE 1
#define MPI_LONG_DOUBLE 2
#define MPI_INT 3
#define MPI_CHAR 4
#define MPI_LONG 5
#define MPI_BYTE 6
#define MPI_REAL 7
#define MPI_COMPLEX 8
#define MPI_SUM 0
#define MPI_MIN 1
#define MPI_MAX 2
#define MPI_LOR 3
#define MPI_LAND 4
#define MPI_SUCCESS 0
#define MPI_STATUSES_IGNORE 0
#define MPI_UNDEFINED -9999
#define MPI_REQUEST_NULL 0
#define MPI_INFO_NULL 0
#define MPI_ANY_SOURCE 1
#define MPI_ANY_TAG 1
#endif

/* small real number */
#define EPS 1.0e-6f
#define JXF_HPCSR 5555
#define JXF_UNITIALIZED -999
#define JXF_REAL_MIN_EXP FLT_MIN_EXP

/* max and min operation */
#define jxf_max(a,b)  (((a)<(b)) ? (b) : (a))
#define jxf_min(a,b)  (((a)<(b)) ? (a) : (b))
#define jxf_abs(a)    (((a)>0) ? (a) : -(a))
#define jxf_round(x)  (((x) < 0.0f) ? ((JXF_Int)(x - 0.5f)) : ((JXF_Int)(x + 0.5f)))
#define jxf_pow2(i)   (1 << (i))

/* Global variable used in jxf-error checking */
extern JXF_Int jxf__kdim_memory_error;
extern JXF_Int jxf__global_mvcpu_flag;
extern JXF_Int jxf__global_error;
#define jxf_error_flag     jxf__global_error
#define JXF_ERROR_ARG      4                /* argument error */
#define JXF_ERROR_GENERIC  1
#define JXF_ERROR_MEMORY   2
#define JXF_ERROR_CONV     256

/* Thread */
#if JXF_USING_OPENMP || defined (JXF_USING_PGCC_SMP)
JXF_Int jxf_NumThreads( void );
JXF_Int jxf_NumActiveThreads( void );
JXF_Int jxf_GetThreadNum( void );
#else
#define jxf_NumThreads() 1
#define jxf_NumActiveThreads() 1
#define jxf_GetThreadNum() 0
#endif

/* Memory */
#define jxf_TAlloc(type, count)        ( (type *)jxf_MAlloc((size_t)(sizeof(type) * (count))) )
#define jxf_CTAlloc(type, count)       ( (type *)jxf_CAlloc((size_t)(count), (size_t)sizeof(type)) )
#define jxf_TReAlloc(ptr, type, count) ( (type *)jxf_ReAlloc((char *)ptr, (size_t)(sizeof(type) * (count))) )
#define jxf_TReAlloc_v2(ptr, old_type, old_count, new_type, new_count) \
   ( (new_type *) jxf_ReAlloc_v2((char *)ptr, (size_t)(sizeof(old_type)*(old_count)), (size_t)(sizeof(new_type)*(new_count))) )
#define jxf_TFree(ptr)                 ( jxf_Free((char *)ptr), ptr = NULL )
#define jxf_TMemcpy(dst, src, type, count)   (memcpy((void *)(dst), (void *)(src), (size_t)(sizeof(type) * (count))))

/* Fortran, BLAS, LAPACK */
#define jxf_F90_NAME(name,NAME) name##_
#define jxf_F90_NAME_BLAS(name,NAME) name##_
#define jxf_F90_NAME_LAPACK(name,NAME) name##_

/* Partition */
#define JXF_OMP_GET_START_END(procid,nprocs,n,start,end) \
        if((procid)<(n)%(nprocs))                       \
        {                                               \
            (end)=(n)/(nprocs)+1;                       \
            (start)=(end)*(procid);                     \
        }                                               \
        else                                            \
        {                                               \
            (end)=(n)/(nprocs);                         \
            (start)=(end)*(procid)+(n)%(nprocs);        \
        }                                               \
        (end)=(end)+(start);

/* common struct */
struct jxf_Solver_struct;
typedef struct jxf_Solver_struct *JXF_Solver;

struct jxf_Matrix_struct;
typedef struct jxf_Matrix_struct *JXF_Matrix;

struct jxf_Vector_struct;
typedef struct jxf_Vector_struct *JXF_Vector;

struct jxf_ParCSRMatrix_struct;
typedef struct jxf_ParCSRMatrix_struct *JXF_ParCSRMatrix;

struct jxf_ParVector_struct;
typedef struct jxf_ParVector_struct *JXF_ParVector;

struct jxf_IJMatrix_struct;
typedef struct jxf_IJMatrix_struct *JXF_IJMatrix;

struct jxf_IJVector_struct;
typedef struct jxf_IJVector_struct *JXF_IJVector;

//hpcsr struct
struct jxf_hpCSRMatrix_struct;
typedef struct jxf_hpCSRMatrix_struct *JXF_hpCSRMatrix;

/*------------------------------------------------------------*/
/*--------------------- Extern Definition --------------------*/
/*------------------------------------------------------------*/

extern JXF_UInt jxf_long_size_length_rap;
extern JXF_UInt jxf_long_size_length_interp;
extern JXF_UInt jxf_total_calloc_global_rap;
extern JXF_UInt jxf_total_calloc_global_interp;
extern JXF_UInt jxf_real_long_size_length_rap;
extern JXF_UInt jxf_real_long_size_length_interp;
extern JXF_Real jxf_total_elapsed_time_matvec;
extern JXF_Int  jxf_total_iter_index;

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_DataExchangeResponse
 * \brief In the fill_response() function the user needs to set the recv__buf
 *        and the response_message_size.  Memory of size send_response_storage has been
 *        alllocated for the send_buf (in exchange_data) - if more is needed, then
 *        realloc and adjust the send_response_storage. 
 *        The realloc amount should be storage+overhead. 
 *        If the response is an empty "confirmation" message, then set
 *        response_message_size = 0 (and do not modify the send_buf).
 */
typedef struct
{
   JXF_Int  (*fill_response)(void* recv_buf, JXF_Int contact_size, 
                            JXF_Int contact_proc, void* response_obj, 
                            MPI_Comm comm, void** response_buf, 
                            JXF_Int* response_message_size);
   JXF_Int  send_response_overhead; /*set by exchange data */
   JXF_Int  send_response_storage;  /*storage allocated for send_response_buf*/
   void   *data1;                  /*data fields user may want to access in fill_response */
   void   *data2;
   
} jxf_DataExchangeResponse;


/*!
 * \struct jxf_BinaryTree
 */
typedef struct
{
   JXF_Int  parent_id;
   JXF_Int  num_child;
   JXF_Int *child_id;
   
} jxf_BinaryTree;

#define jxf_BinaryTreeParentId(tree)    (tree->parent_id)
#define jxf_BinaryTreeNumChild(tree)    (tree->num_child)
#define jxf_BinaryTreeChildIds(tree)    (tree->child_id)
#define jxf_BinaryTreeChildId(tree, i)  (tree->child_id[i])


/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/utilities/exchange.c */
JXF_Int 
jxf_DataExchangeList( JXF_Int                       num_contacts, 
		     JXF_Int                      *contact_proc_list, 
		     void                     *contact_send_buf, 
		     JXF_Int                      *contact_send_buf_starts, 
		     JXF_Int                       contact_obj_size, 
                     JXF_Int                       response_obj_size,
		     jxf_DataExchangeResponse  *response_obj, 
		     JXF_Int                       max_response_size, 
                     JXF_Int                       rnum, 
                     MPI_Comm                  comm,  
                     void                    **p_response_recv_buf, 
                     JXF_Int                     **p_response_recv_buf_starts );
JXF_Int jxf_CreateBinaryTree( JXF_Int myid, JXF_Int num_procs, jxf_BinaryTree *tree );
JXF_Int jxf_DestroyBinaryTree( jxf_BinaryTree *tree );

/* csrc/utilities/memory.c */
JXF_Int    jxf_OutOfMemory( size_t size );
char * jxf_MAlloc( size_t size );
char * jxf_CAlloc( size_t count, size_t elt_size );
char * jxf_ReAlloc( char *ptr, size_t size );
void * jxf_ReAlloc_v2( void *ptr, size_t old_size, size_t new_size );
void   jxf_Free( char *ptr );

/* csrc/utilities/error.c */
void jxf_error_handler( char *filename, JXF_Int line, JXF_Int ierr, const char *msg );
#define jxf_error(IERR)  jxf_error_handler(__FILE__, __LINE__, IERR, NULL)
#define jxf_error_w_msg(IERR, msg) jxf_error_handler(__FILE__, __LINE__, IERR, msg)
#define jxf_error_in_arg(IARG)  jxf_error(JXF_ERROR_ARG | IARG<<3)
#ifdef NDEBUG
#define jxf_assert(EX)
#else
#define jxf_assert(EX) \
 if(!(EX)) {jxf_fprintf(stderr,"jxf_assert failed: %s\n", #EX); jxf_error(1);}
#endif

/* /csrc/utilities/sort.c */
void jxf_swap( JXF_Int *v, JXF_Int i, JXF_Int j );
void jxf_swap1( JXF_Real *v, JXF_Int i, JXF_Int j );
void jxf_swap2( JXF_Int *v, JXF_Real *w, JXF_Int i, JXF_Int j );
void jxf_swap2i( JXF_Int *v, JXF_Int *w, JXF_Int i, JXF_Int j );
void jxf_swap3_d_perm( JXF_Int *v, JXF_Real *w, JXF_Int *z, JXF_Int i, JXF_Int j );
void jxf_swap3i( JXF_Int *v, JXF_Int *w, JXF_Int *z, JXF_Int i, JXF_Int j );
void jxf_qsortd( JXF_Real *v, JXF_Int left, JXF_Int right );
void jxf_qsort0( JXF_Int *v, JXF_Int left, JXF_Int right );
void jxf_qsort1( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void jxf_qsort2( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void jxf_qsort2i( JXF_Int *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void jxf_qsort2abs( JXF_Int *v, JXF_Real *w, JXF_Int left, JXF_Int right );
void jxf_qsort3i( JXF_Int *v, JXF_Int *w, JXF_Int *z, JXF_Int left, JXF_Int right );
void jxf_qsort3ir( JXF_Int *v, JXF_Real *w, JXF_Int *z, JXF_Int left, JXF_Int right );
void jxf_diQuickSort12( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void jxf_diQuickSort21( JXF_Real *v, JXF_Int *w, JXF_Int left, JXF_Int right );
void jxf_union2( JXF_Int n1, JXF_Int *arr1, JXF_Int n2, JXF_Int *arr2, JXF_Int *n3, JXF_Int *arr3, JXF_Int *map1, JXF_Int *map2 );
void jxf_swap3_d( JXF_Real *v, JXF_Int *w, JXF_Int *z, JXF_Int i, JXF_Int j );
void jxf_qsort3( JXF_Real *v, JXF_Int *w, JXF_Int *z, JXF_Int left, JXF_Int right );

/* /csrc/utilities/search.c */
JXF_Int jxf_BinarySearch( JXF_Int *list, JXF_Int value, JXF_Int list_length );
JXF_Int *jxf_LowerBound( JXF_Int *first, JXF_Int *last, JXF_Int value );
JXF_Int jxf_BinarySearch2( JXF_Int *list, JXF_Int value, JXF_Int low, JXF_Int high, JXF_Int *spot );

/* /csrc/utilities/array.c */
JXF_Int jxf_DoubleArrayPrint( JXF_Real *x, JXF_Int ndof, char *filename );
JXF_Int jxf_IntegerArrayGetInterp( JXF_Int *icor_interp, JXF_Int n_fine, JXF_Int nbl, JXF_Int nbr );
JXF_Int jxf_IntegerArraySetConstantValues( JXF_Int n, JXF_Int *x, JXF_Int val );
JXF_Int jxf_IntegerArrayModFine2Coarse( JXF_Int num_rows_A, JXF_Int *fine_to_coarse );
JXF_Int jxf_IntegerArrayGenerateIcor( JXF_Int num_rows_A,
                                    JXF_Int num_cols_P,
                                    JXF_Int *fine_to_coarse,
                                    JXF_Int nbl,
                                    JXF_Int nbr,
                                    JXF_Int *CF_marker,
                                    JXF_Int *icor );
JXF_Int jxf_IntegerArrayCopy( JXF_Int size, JXF_Int *x, JXF_Int *y );
JXF_Int jxf_DoubleArrayCopy( JXF_Int size, JXF_Real *x, JXF_Real *y );
JXF_Int jxf_DoubleArrayReciprocalMap( JXF_Real *u, JXF_Int n, JXF_Int *map );
JXF_Real **jxf_DyadicDoubleArrayCTAlloc( JXF_Int first, JXF_Int second );
JXF_Real ***jxf_TeracidicDoubleArrayCTAlloc( JXF_Int first, JXF_Int second, JXF_Int third );
JXF_Real ****jxf_TetravalentDoubleArrayCTAlloc( JXF_Int first, JXF_Int second, JXF_Int third, JXF_Int fourth );
void jxf_DyadicDoubleArrayFree( JXF_Real **ptr, JXF_Int first );
void jxf_TeracidicDoubleArrayFree( JXF_Real ***ptr, JXF_Int first, JXF_Int second );
void jxf_TetravalentDoubleArrayFree( JXF_Real ****ptr, JXF_Int first, JXF_Int second, JXF_Int third );
JXF_Real jxf_DoubleArrayAbsMaxElement( JXF_Real *x, JXF_Int n );
JXF_Real jxf_DoubleArrayAbsMinElement( JXF_Real *x, JXF_Int n );
JXF_Real jxf_DoubleArrayMaxElement( JXF_Real *x, JXF_Int n );

/* /csrc/utilities/intersection.c */
JXF_Int jxf_MultiSetsIntersection( JXF_Int numset, JXF_Int **set, JXF_Int *size, JXF_Int **c_ptr, JXF_Int *nc_ptr, JXF_Int *empty_ptr );
JXF_Int jxf_TwoSetsIntersection( JXF_Int *a, JXF_Int na, JXF_Int *b, JXF_Int nb, JXF_Int **c_ptr, JXF_Int *nc_ptr, JXF_Int *empty_ptr );

/* csrc/utilities/random.c */
void jxf_SeedRand( JXF_Int seed );
JXF_Real jxf_Rand();                   
                                         
/* csrc/utilities/timer.c */
JXF_Real 
jxf_time_getWallclockSeconds(void);
void
jxf_GetWallTime( MPI_Comm  comm, 
                char     *fctname, 
                JXF_Real    starttime, 
                JXF_Real    endtime, 
                JXF_Int       allproc,
                JXF_Int       dp );
void
jxf_GetWallTime2( MPI_Comm comm, char *fctname, JXF_Real totaltime, JXF_Int allproc, JXF_Int dp );
JXF_Real
jxf_GetWallTimeMax( MPI_Comm  comm, 
                   JXF_Real    starttime, 
                   JXF_Real    endtime );

/* csrc/utilities/printf.c */
JXF_Int
jxf_new_format( const char *format, char **newformat_ptr );
JXF_Int
jxf_free_format( char *newformat );
JXF_Int
jxf_printf( const char *format, ... );
JXF_Int
jxf_fprintf( FILE *stream, const char *format, ... );
JXF_Int
jxf_sprintf( char *s, const char *format, ... );
JXF_Int
jxf_scanf( const char *format, ... );
JXF_Int
jxf_fscanf( FILE *stream, const char *format, ... );
JXF_Int
jxf_sscanf( char *s, const char *format, ... );

/* csrc/utilities/set.c */
void jxf_set_mvcpu_handler( JXF_Int new_mvcpu_flag );

/* csrc/utilities/mpistubs.c */
JXF_Int
jxf_MPI_Init( int *argc, char ***argv );
JXF_Int
jxf_MPI_Finalize();
JXF_Int
jxf_MPI_Abort( MPI_Comm comm, JXF_Int errorcode );
JXF_Real
jxf_MPI_Wtime();
JXF_Real
jxf_MPI_Wtick();
JXF_Int
jxf_MPI_Barrier( MPI_Comm comm );
JXF_Int
jxf_MPI_Comm_create( MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm );
JXF_Int
jxf_MPI_Comm_dup( MPI_Comm comm, MPI_Comm *newcomm );
JXF_Int
jxf_MPI_Comm_size( MPI_Comm comm, JXF_Int *size );
JXF_Int
jxf_MPI_Comm_rank( MPI_Comm comm, JXF_Int *rank );
JXF_Int
jxf_MPI_Comm_free( MPI_Comm *comm );
JXF_Int
jxf_MPI_Comm_group( MPI_Comm comm, MPI_Group *group );
JXF_Int
jxf_MPI_Comm_split( MPI_Comm comm, JXF_Int n, JXF_Int m, MPI_Comm *comms );
JXF_Int
jxf_MPI_Group_incl( MPI_Group group, JXF_Int n, JXF_Int *ranks, MPI_Group *newgroup );
JXF_Int
jxf_MPI_Group_free( MPI_Group *group );
JXF_Int
jxf_MPI_Address( void *location, MPI_Aint *address );
JXF_Int
jxf_MPI_Get_count( MPI_Status *status, MPI_Datatype datatype, JXF_Int *count );
JXF_Int
jxf_MPI_Alltoall( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm );
JXF_Int
jxf_MPI_Allgather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm );
JXF_Int
jxf_MPI_Allgatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, MPI_Comm comm );
JXF_Int
jxf_MPI_Gather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Gatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Scatter( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Scatterv( void *sendbuf, JXF_Int *sendcounts, JXF_Int *displs, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Bcast( void *buffer, JXF_Int count, MPI_Datatype datatype, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Send( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm );
JXF_Int
jxf_MPI_Recv( void *buf, JXF_Int count, MPI_Datatype  datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status );
JXF_Int
jxf_MPI_Isend( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request );
JXF_Int
jxf_MPI_Irecv( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Request *request );
JXF_Int
jxf_MPI_Send_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request );
JXF_Int
jxf_MPI_Recv_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest,
                     JXF_Int tag, MPI_Comm comm, MPI_Request *request );
JXF_Int
jxf_MPI_Irsend( void *buf, JXF_Int count, MPI_Datatype datatype,
                  JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request );
JXF_Int
jxf_MPI_Startall( JXF_Int count, MPI_Request *array_of_requests );
JXF_Int
jxf_MPI_Probe( JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status );
JXF_Int
jxf_MPI_Iprobe( JXF_Int source, JXF_Int tag, MPI_Comm comm, JXF_Int *flag, MPI_Status *status );
JXF_Int
jxf_MPI_Test( MPI_Request *request, JXF_Int *flag, MPI_Status *status );
JXF_Int
jxf_MPI_Testall( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *flag, MPI_Status *array_of_statuses );
JXF_Int
jxf_MPI_Wait( MPI_Request *request, MPI_Status *status );
JXF_Int
jxf_MPI_Waitall( JXF_Int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses );
JXF_Int
jxf_MPI_Waitany( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *index, MPI_Status *status );
JXF_Int
jxf_MPI_Allreduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
JXF_Int
jxf_MPI_Reduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, JXF_Int root, MPI_Comm comm );
JXF_Int
jxf_MPI_Scan( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
JXF_Int
jxf_MPI_Request_free( MPI_Request *request );
JXF_Int
jxf_MPI_Type_contiguous( JXF_Int count, MPI_Datatype oldtype, MPI_Datatype *newtype );
JXF_Int
jxf_MPI_Type_vector( JXF_Int count, JXF_Int blocklength, JXF_Int stride, MPI_Datatype oldtype, MPI_Datatype *newtype );
JXF_Int
jxf_MPI_Type_hvector( JXF_Int count, JXF_Int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype );
JXF_Int
jxf_MPI_Type_struct( JXF_Int count, JXF_Int *array_of_blocklengths, MPI_Aint *array_of_displacements,
                       MPI_Datatype *array_of_types, MPI_Datatype *newtype );
JXF_Int
jxf_MPI_Type_commit( MPI_Datatype *datatype );
JXF_Int
jxf_MPI_Type_free( MPI_Datatype *datatype );
JXF_Int
jxf_MPI_Op_free( MPI_Op *op );
JXF_Int
jxf_MPI_Op_create( MPI_User_function *function, int commute, MPI_Op *op );

/* csrc/euclid/parcsr_euclid.c */
JXF_Int 
JXF_EuclidCreate( MPI_Comm comm, JXF_Solver *solver );
JXF_Int 
JXF_EuclidDestroy( JXF_Solver solver );
JXF_Int 
JXF_EuclidSetup( JXF_Solver solver, JXF_hpCSRMatrix A );
JXF_Int 
JXF_EuclidSolve( JXF_Solver solver,
                JXF_hpCSRMatrix A,
                JXF_ParVector bb,
                JXF_ParVector xx );
JXF_Int 
JXF_EuclidSetParams( JXF_Solver solver, JXF_Int argc, char *argv[] );
JXF_Int 
JXF_EuclidSetParamsFromFile( JXF_Solver solver, char *filename );
JXF_Int 
JXF_EuclidSetLevel( JXF_Solver solver, JXF_Int level );
JXF_Int 
JXF_EuclidSetBJ( JXF_Solver solver, JXF_Int bj );
JXF_Int 
JXF_EuclidSetStats( JXF_Solver solver, JXF_Int eu_stats );
JXF_Int 
JXF_EuclidSetMem( JXF_Solver solver, JXF_Int eu_mem );
JXF_Int 
JXF_EuclidSetILUT( JXF_Solver solver, JXF_Real ilut );
JXF_Int 
JXF_EuclidSetSparseA( JXF_Solver solver, JXF_Real sparse_A );
JXF_Int 
JXF_EuclidSetRowScale( JXF_Solver solver, JXF_Int row_scale );

/* csrc/utilities/thread.c */
void
jxf_GetSimpleThreadPartition( JXF_Int *begin, JXF_Int *end, JXF_Int n );

/* csrc/utilities/prefixsum.c */
void
jxf_prefix_sum_pair(JXF_Int *in_out1, JXF_Int *sum1, JXF_Int *in_out2, JXF_Int *sum2, JXF_Int *workspace);

#endif
