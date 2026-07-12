//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_euclid.h -- head files for Parallel ILU Solver
 *  Date: 2013/01/20
 *
 *  Created by Yue Xiaoqiang
 */

#ifndef JXF_EUCLID_HEADER
#define JXF_EUCLID_HEADER

#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

#define JXFPAMG_GET_ROW

#ifndef JXF_EUCLID_CONF_DH
#define JXF_EUCLID_CONF_DH

#define JXF_MAX_MPI_TASKS 50000
#define JXF_TRIPLES_FORMAT "%i %i %1.8e\n"

#define JXF_EUCLID_EXIT jxf_MPI_Abort(jxf_comm_dh, -1)

#define JXF_ERRCHKA \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("", __FUNC__, __FILE__, __LINE__); \
        if (jxf_logFile != NULL) { \
          jxf_printErrorMsg(jxf_logFile); \
          jxf_closeLogfile_dh(); \
        } \
        jxf_printErrorMsg(stderr); \
        if (jxf_myid_dh == 0) { \
          jxf_Mem_dhPrint(jxf_mem_dh, stderr, jxf_false); \
        } \
        JXF_EUCLID_EXIT; \
      }

#define JXF_PIVOT_FIX_DEFAULT 1e-3

#define JXF_MALLOC_DH(s) jxf_Mem_dhMalloc(jxf_mem_dh, (s))
#define JXF_FREE_DH(p) jxf_Mem_dhFree(jxf_mem_dh, p)

#define JXF_PRIVATE_MALLOC malloc
#define JXF_PRIVATE_FREE free

#endif

#ifndef JXF_MACROS_DH
#define JXF_MACROS_DH

#ifndef JXF_MAX
#define JXF_MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef JXF_MIN
#define JXF_MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define _JXF_MATLAB_ZERO_  1e-100

#define JXF_CHECK_MPI_V_ERROR(errCode) \
      { \
        if (errCode) { \
          jxf_setError_dh("MPI error!", __FUNC__, __FILE__, __LINE__); \
          return; \
        } \
      }

#define JXF_CHECK_MPI_ERROR(errCode) \
      { \
        if (errCode) { \
          jxf_setError_dh("MPI error!", __FUNC__, __FILE__, __LINE__); \
          return(errCode); \
        } \
      }

#define JXF_SET_V_ERROR(msg) \
      { jxf_setError_dh(msg, __FUNC__, __FILE__, __LINE__); \
        return; \
      }

#define JXF_SET_ERROR(retval, msg) \
      { jxf_setError_dh(msg, __FUNC__, __FILE__, __LINE__); \
        return (retval); \
      }

#define JXF_CHECK_V_ERROR \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("",  __FUNC__, __FILE__, __LINE__); \
        return; \
      }

#define JXF_CHECK_ERROR(retval) \
      if (jxf_errFlag_dh) { \
        jxf_setError_dh("",  __FUNC__, __FILE__, __LINE__); \
        return (retval); \
      }

#define JXF_SET_INFO(msg) jxf_setInfo_dh(msg, __FUNC__, __FILE__, __LINE__);

#define JXF_START_FUNC_DH \
      jxf_dh_StartFunc(__FUNC__, __FILE__, __LINE__, 1); \
      {

#define JXF_END_FUNC_DH \
      } \
      jxf_dh_EndFunc(__FUNC__, 1);

#define JXF_END_FUNC_VAL(a) \
      jxf_dh_EndFunc(__FUNC__, 1); \
      return a ; \
      }

#define JXF_START_FUNC_DH_2
#define JXF_END_FUNC_DH_2
#define JXF_END_FUNC_VAL_2(a) return a;

#endif

#ifndef JXF_COMMON_DH
#define JXF_COMMON_DH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#define JXF_REAL_DH JXF_Real

#include "jxf_mv.h"

typedef struct _jxf_matgenfd * jxf_MatGenFD;
typedef struct _jxf_subdomain_dh * jxf_SubdomainGraph_dh;
typedef struct _jxf_timer_dh * jxf_Timer_dh;
typedef struct _jxf_parser_dh * jxf_Parser_dh;
typedef struct _jxf_timeLog_dh * jxf_TimeLog_dh;
typedef struct _jxf_mem_dh * jxf_Mem_dh;
typedef struct _jxf_mat_dh * jxf_Mat_dh;
typedef struct _jxf_factor_dh * jxf_Factor_dh;
typedef struct _jxf_vec_dh * jxf_Vec_dh;
typedef struct _jxf_numbering_dh * jxf_Numbering_dh;
typedef struct _jxf_hash_dh * jxf_Hash_dh;
typedef struct _jxf_hash_i_dh * jxf_Hash_i_dh;
typedef struct _jxf_mpi_interface_dh * jxf_Euclid_dh;
typedef struct _jxf_sortedList_dh * jxf_SortedList_dh;
typedef struct _jxf_extrows_dh * jxf_ExternalRows_dh;
typedef struct _jxf_stack_dh * jxf_Stack_dh;
typedef struct _jxf_queue_dh * jxf_Queue_dh;
typedef struct _jxf_sortedset_dh * jxf_SortedSet_dh;

typedef JXF_Int jxf_bool;
#define jxf_true 1
#define jxf_false 0

extern jxf_Parser_dh jxf_parser_dh;
extern jxf_TimeLog_dh jxf_tlog_dh;
extern jxf_Mem_dh jxf_mem_dh;
extern FILE *jxf_logFile;
extern JXF_Int jxf_np_dh;
extern JXF_Int jxf_myid_dh;
extern MPI_Comm jxf_comm_dh;
extern jxf_bool jxf_ignoreMe;
extern JXF_Int jxf_ref_counter;
extern jxf_bool jxf_errFlag_dh;
extern void jxf_setInfo_dh( char *msg, char *function, char *file, JXF_Int line );
extern void jxf_setError_dh( char *msg, char *function, char *file, JXF_Int line );
extern void jxf_printErrorMsg( FILE *fp );

#ifndef JXF_MPI_MAX_ERROR_STRING
#define JXF_MPI_MAX_ERROR_STRING 256
#endif

#define JXF_MSG_BUF_SIZE_DH JXF_MAX(1024, JXF_MPI_MAX_ERROR_STRING)
extern char jxf_msgBuf_dh[JXF_MSG_BUF_SIZE_DH];

extern void jxf_openLogfile_dh( JXF_Int argc, char *argv[] );
extern void jxf_closeLogfile_dh();

extern jxf_bool jxf_logInfoToStderr;
extern jxf_bool jxf_logInfoToFile;
extern jxf_bool jxf_logFuncsToStderr;
extern jxf_bool jxf_logFuncsToFile;

extern void jxf_Error_dhStartFunc( char *function, char *file, JXF_Int line );
extern void jxf_Error_dhEndFunc( char *function );
extern void jxf_dh_StartFunc( char *function, char *file, JXF_Int line, JXF_Int priority );
extern void jxf_dh_EndFunc( char *function, JXF_Int priority );
extern void jxf_printFunctionStack( FILE *fp );
extern void jxf_EuclidInitialize( JXF_Int argc, char *argv[], char *help );
extern void jxf_EuclidFinalize();
extern jxf_bool jxf_EuclidIsInitialized();
extern void jxf_printf_dh( char *fmt, ... );
extern void jxf_fprintf_dh( FILE *fp, char *fmt, ... );
extern void jxf_echoInvocation_dh( MPI_Comm comm, char *prefix, JXF_Int argc, char *argv[] );

#endif

#ifndef JXF_EXTERNAL_ROWS_DH_H
#define JXF_EXTERNAL_ROWS_DH_H

extern void jxf_ExternalRows_dhCreate( jxf_ExternalRows_dh *er );
extern void jxf_ExternalRows_dhDestroy( jxf_ExternalRows_dh er );
extern void jxf_ExternalRows_dhInit( jxf_ExternalRows_dh er, jxf_Euclid_dh ctx );
extern void jxf_ExternalRows_dhRecvRows( jxf_ExternalRows_dh extRows );
extern void jxf_ExternalRows_dhSendRows( jxf_ExternalRows_dh extRows );
extern void jxf_ExternalRows_dhGetRow( jxf_ExternalRows_dh er,
               JXF_Int globalRow, JXF_Int *len, JXF_Int **cval, JXF_Int **fill, JXF_REAL_DH **aval );

struct _jxf_extrows_dh
{
    jxf_SubdomainGraph_dh sg;
    jxf_Factor_dh F;
    
    MPI_Status status[JXF_MAX_MPI_TASKS];
    MPI_Request req1[JXF_MAX_MPI_TASKS]; 
    MPI_Request req2[JXF_MAX_MPI_TASKS];
    MPI_Request req3[JXF_MAX_MPI_TASKS]; 
    MPI_Request req4[JXF_MAX_MPI_TASKS];
    MPI_Request cval_req[JXF_MAX_MPI_TASKS];
    MPI_Request fill_req[JXF_MAX_MPI_TASKS];
    MPI_Request aval_req[JXF_MAX_MPI_TASKS];
    
    JXF_Int rcv_row_counts[JXF_MAX_MPI_TASKS];
    JXF_Int rcv_nz_counts[JXF_MAX_MPI_TASKS];
    JXF_Int *rcv_row_lengths[JXF_MAX_MPI_TASKS];
    JXF_Int *rcv_row_numbers[JXF_MAX_MPI_TASKS];
    
    JXF_Int *cvalExt;
    JXF_Int *fillExt;
    JXF_REAL_DH *avalExt;
    
    jxf_Hash_dh rowLookup;
    
    JXF_Int *my_row_counts;
    JXF_Int *my_row_numbers;
    
    JXF_Int nzSend;
    JXF_Int *cvalSend;
    JXF_Int *fillSend;
    JXF_REAL_DH *avalSend;
    
    jxf_bool debug;
};

#endif

#ifndef JXF_FACTOR_DH
#define JXF_FACTOR_DH

struct _jxf_factor_dh
{
    JXF_Int m, n;
    
    JXF_Int id;
    JXF_Int beg_row;
    JXF_Int first_bdry;
    JXF_Int bdry_count;
    
    jxf_bool blockJacobi;
    
    JXF_Int *rp;
    JXF_Int *cval;
    JXF_REAL_DH *aval;
    JXF_Int *fill;
    JXF_Int *diag;
    JXF_Int alloc;
    
    JXF_Int num_recvLo, num_recvHi;
    JXF_Int num_sendLo, num_sendHi;
    JXF_Real *work_y_lo;
    JXF_Real *work_x_hi;
    JXF_Real *sendbufLo, *sendbufHi;
    JXF_Int *sendindLo, *sendindHi;
    JXF_Int sendlenLo, sendlenHi;
    jxf_bool solveIsSetup;
    jxf_Numbering_dh numbSolve;
    
    MPI_Request recv_reqLo[JXF_MAX_MPI_TASKS], recv_reqHi[JXF_MAX_MPI_TASKS];
    MPI_Request send_reqLo[JXF_MAX_MPI_TASKS], send_reqHi[JXF_MAX_MPI_TASKS];
    MPI_Request requests[JXF_MAX_MPI_TASKS];
    MPI_Status status[JXF_MAX_MPI_TASKS];
    
    jxf_bool debug;
};

extern void jxf_Factor_dhCreate( jxf_Factor_dh *mat );
extern void jxf_Factor_dhDestroy( jxf_Factor_dh mat );
extern void jxf_Factor_dhTranspose( jxf_Factor_dh matIN, jxf_Factor_dh *matOUT );
extern void jxf_Factor_dhInit( void *A, jxf_bool fillFlag, jxf_bool avalFlag,
                              JXF_Real rho, JXF_Int id, JXF_Int beg_rowP, jxf_Factor_dh *F );
extern void jxf_Factor_dhReallocate( jxf_Factor_dh F, JXF_Int used, JXF_Int additional );
extern void jxf_Factor_dhSolveSetup( jxf_Factor_dh mat, jxf_SubdomainGraph_dh sg );
extern void jxf_Factor_dhSolve( JXF_Real *rhs, JXF_Real *lhs, jxf_Euclid_dh ctx );
extern void jxf_Factor_dhSolveSeq( JXF_Real *rhs, JXF_Real *lhs, jxf_Euclid_dh ctx );
extern JXF_Real jxf_Factor_dhCondEst( jxf_Factor_dh mat, jxf_Euclid_dh ctx );
extern JXF_Real jxf_Factor_dhMaxValue( jxf_Factor_dh mat );
extern JXF_Real jxf_Factor_dhMaxPivotInverse( jxf_Factor_dh mat );
extern JXF_Int jxf_Factor_dhReadNz( jxf_Factor_dh mat );
extern void jxf_Factor_dhPrintTriples( jxf_Factor_dh mat, char *filename );
extern void jxf_Factor_dhPrintGraph( jxf_Factor_dh mat, char *filename );
extern void jxf_Factor_dhPrintDiags( jxf_Factor_dh mat, FILE *fp );
extern void jxf_Factor_dhPrintRows( jxf_Factor_dh mat, FILE *fp );

#endif

#ifndef JXF_VEC_DH_H
#define JXF_VEC_DH_H

struct _jxf_vec_dh
{
    JXF_Int n;
    JXF_Real *vals;
};

extern void jxf_Vec_dhCreate( jxf_Vec_dh *v );
extern void jxf_Vec_dhDestroy( jxf_Vec_dh v );
extern void jxf_Vec_dhInit( jxf_Vec_dh v, JXF_Int size );
extern void jxf_Vec_dhDuplicate( jxf_Vec_dh v, jxf_Vec_dh *out );
extern void jxf_Vec_dhCopy( jxf_Vec_dh x, jxf_Vec_dh y );
extern void jxf_Vec_dhSet( jxf_Vec_dh v, JXF_Real value );
extern void jxf_Vec_dhSetRand( jxf_Vec_dh v );
extern void jxf_Vec_dhRead( jxf_Vec_dh *v, JXF_Int ignore, char *filename );
extern void jxf_Vec_dhReadBIN( jxf_Vec_dh *v, char *filename );
extern void jxf_Vec_dhPrint( jxf_Vec_dh v, jxf_SubdomainGraph_dh sg, char *filename );
extern void jxf_Vec_dhPrintBIN( jxf_Vec_dh v, jxf_SubdomainGraph_dh sg, char *filename );

#endif

#ifndef JXF_MATGENFD_DH_DH
#define JXF_MATGENFD_DH_DH

struct _jxf_matgenfd
{
    jxf_bool allocateMem;
    JXF_Int px, py, pz;
    jxf_bool threeD;
    JXF_Int m;
    JXF_Int cc;
    JXF_Real hh;
    JXF_Int id;
    JXF_Int np;
    JXF_Real stencil[8];
    
    JXF_Real a, b, c, d, e, f, g, h;
    
    JXF_Int first;
    jxf_bool debug;
    
    JXF_Real bcX1, bcX2;
    JXF_Real bcY1, bcY2;
    JXF_Real bcZ1, bcZ2;
    
    JXF_Real (*A)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*B)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*C)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*D)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*E)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*F)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*G)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
    JXF_Real (*H)( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
};

extern void jxf_MatGenFD_Create( jxf_MatGenFD *mg );
extern void jxf_MatGenFD_Destroy( jxf_MatGenFD mg );
extern void jxf_MatGenFD_Run( jxf_MatGenFD mg, JXF_Int id, JXF_Int np, jxf_Mat_dh *A, jxf_Vec_dh *rhs );
extern JXF_Real jxf_konstant( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
extern JXF_Real jxf_e2_xy( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );

#define JXF_BOX1_X1 0.1
#define JXF_BOX1_X2 0.4
#define JXF_BOX1_Y1 0.1
#define JXF_BOX1_Y2 0.4
#define JXF_BOX2_X1 0.6
#define JXF_BOX2_X2 0.9
#define JXF_BOX2_Y1 0.1
#define JXF_BOX2_Y2 0.4
#define JXF_BOX3_X1 0.2
#define JXF_BOX3_X2 0.8
#define JXF_BOX3_Y1 0.6
#define JXF_BOX3_Y2 0.8
#define JXF_BOX1_DD 10
#define JXF_BOX2_DD 100
#define JXF_BOX3_DD 50

extern JXF_Real jxf_box_1( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );
extern JXF_Real jxf_box_2( JXF_Real coeff, JXF_Real x, JXF_Real y, JXF_Real z );

#endif

#ifndef JXF_MAT_DH_DH
#define JXF_MAT_DH_DH

#define JXF_MAT_DH_BINS 10
#define JXF_MATVEC_TIME 0
#define JXF_MATVEC_MPI_TIME 1
#define JXF_MATVEC_MPI_TIME2 5
#define JXF_MATVEC_TOTAL_TIME 2
#define JXF_MATVEC_RATIO 3
#define JXF_MATVEC_WORDS 4

struct _jxf_mat_dh
{
    JXF_Int m, n;
    JXF_Int beg_row;
    JXF_Int bs;
    
    JXF_Int *rp;
    JXF_Int *len;
    JXF_Int *cval;
    JXF_Int *fill;
    JXF_Int *diag;
    JXF_Real *aval;
    jxf_bool owner;
    
    JXF_Int len_private;
    JXF_Int rowCheckedOut;
    JXF_Int *cval_private;
    JXF_Real *aval_private;
    
    JXF_Int *row_perm;
    
    JXF_Real time[JXF_MAT_DH_BINS];
    JXF_Real time_max[JXF_MAT_DH_BINS];
    JXF_Real time_min[JXF_MAT_DH_BINS];
    jxf_bool matvec_timing;
    
    JXF_Int num_recv;
    JXF_Int num_send;
    MPI_Request *recv_req;
    MPI_Request *send_req;
    JXF_Real *recvbuf, *sendbuf;
    JXF_Int *sendind;
    JXF_Int sendlen;
    JXF_Int recvlen;
    jxf_bool matvecIsSetup;
    jxf_Numbering_dh numb;
    MPI_Status *status;
    
    jxf_bool debug;
};

extern void jxf_Mat_dhCreate( jxf_Mat_dh *mat );
extern void jxf_Mat_dhDestroy( jxf_Mat_dh mat );
extern void jxf_Mat_dhTranspose( jxf_Mat_dh matIN, jxf_Mat_dh *matOUT );
extern void jxf_Mat_dhMakeStructurallySymmetric( jxf_Mat_dh A );
extern void jxf_Mat_dhMatVecSetup( jxf_Mat_dh mat );
extern void jxf_Mat_dhMatVecSetdown( jxf_Mat_dh mat );
extern void jxf_Mat_dhMatVec( jxf_Mat_dh mat, JXF_Real *lhs, JXF_Real *rhs );
extern void jxf_Mat_dhMatVec_omp( jxf_Mat_dh mat, JXF_Real *lhs, JXF_Real *rhs );
extern void jxf_Mat_dhMatVec_uni( jxf_Mat_dh mat, JXF_Real *lhs, JXF_Real *rhs );
extern void jxf_Mat_dhMatVec_uni_omp( jxf_Mat_dh mat, JXF_Real *lhs, JXF_Real *rhs );
extern JXF_Int jxf_Mat_dhReadNz( jxf_Mat_dh mat );
extern void jxf_Mat_dhPrintGraph( jxf_Mat_dh mat, jxf_SubdomainGraph_dh sg, FILE *fp );
extern void jxf_Mat_dhPrintRows( jxf_Mat_dh mat, jxf_SubdomainGraph_dh sg, FILE *fp );
extern void jxf_Mat_dhPrintCSR( jxf_Mat_dh mat, jxf_SubdomainGraph_dh sg, char *filename );
extern void jxf_Mat_dhPrintTriples( jxf_Mat_dh mat, jxf_SubdomainGraph_dh sg, char *filename );
extern void jxf_Mat_dhPrintBIN( jxf_Mat_dh mat, jxf_SubdomainGraph_dh sg, char *filename );
extern void jxf_Mat_dhReadCSR( jxf_Mat_dh *mat, char *filename );
extern void jxf_Mat_dhReadTriples( jxf_Mat_dh *mat, JXF_Int ignore, char *filename );
extern void jxf_Mat_dhReadBIN( jxf_Mat_dh *mat, char *filename );
extern void jxf_Mat_dhPermute( jxf_Mat_dh Ain, JXF_Int *pIN, jxf_Mat_dh *Bout );
extern void jxf_Mat_dhFixDiags( jxf_Mat_dh A );
extern void jxf_Mat_dhPrintDiags( jxf_Mat_dh A, FILE *fp );
extern void jxf_Mat_dhGetRow( jxf_Mat_dh B, JXF_Int globalRow, JXF_Int *len, JXF_Int **ind, JXF_Real **val );
extern void jxf_Mat_dhRestoreRow( jxf_Mat_dh B, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val );
extern void jxf_Mat_dhPartition( jxf_Mat_dh mat, JXF_Int k, JXF_Int **beg_rowOUT,
                                JXF_Int **row_countOUT, JXF_Int **n2oOUT, JXF_Int **o2nOUT );
extern void jxf_Mat_dhZeroTiming( jxf_Mat_dh mat );
extern void jxf_Mat_dhReduceTiming( jxf_Mat_dh mat );
extern void jxf_Mat_dhRowPermute( jxf_Mat_dh );
extern void jxf_dldperm( JXF_Int job, JXF_Int n, JXF_Int nnz, JXF_Int colptr[],
              JXF_Int adjncy[], JXF_Real nzval[], JXF_Int *perm, JXF_Real u[], JXF_Real v[] );

#endif

#ifndef JXF_SUBDOMAIN_GRAPH_DH
#define JXF_SUBDOMAIN_GRAPH_DH

#define JXF_JXF_TIMING_BINS_SG 10

enum{ JXF_TOTAL_SGT, JXF_FIND_NABORS_SGT, JXF_ORDER_BDRY_SGT, JXF_FORM_GRAPH_SGT, JXF_EXCHANGE_PERMS_SGT };

struct _jxf_subdomain_dh
{
    JXF_Int blocks;
    JXF_Int *ptrs, *adj;
    JXF_Int *o2n_sub;
    JXF_Int *n2o_sub;
    JXF_Int colors;
    jxf_bool doNotColor;
    JXF_Int *colorVec;
    
    JXF_Int *beg_row;
    JXF_Int *beg_rowP;
    JXF_Int *row_count;
    JXF_Int *bdry_count;
    
    JXF_Int *loNabors, loCount;
    JXF_Int *hiNabors, hiCount;
    JXF_Int *allNabors, allCount;
    
    JXF_Int m;
    JXF_Int *n2o_row;
    JXF_Int *o2n_col;
    
    jxf_Hash_i_dh o2n_ext;
    jxf_Hash_i_dh n2o_ext;
    
    JXF_Real timing[JXF_JXF_TIMING_BINS_SG];
    jxf_bool debug;
};

extern void jxf_SubdomainGraph_dhCreate( jxf_SubdomainGraph_dh *s );
extern void jxf_SubdomainGraph_dhDestroy( jxf_SubdomainGraph_dh s );
extern void jxf_SubdomainGraph_dhInit( jxf_SubdomainGraph_dh s, JXF_Int blocks, jxf_bool bj, void *A );
extern void jxf_SubdomainGraph_dhColor( jxf_SubdomainGraph_dh s );
extern JXF_Int jxf_SubdomainGraph_dhFindOwner( jxf_SubdomainGraph_dh s, JXF_Int idx, jxf_bool permuted );
extern void jxf_SubdomainGraph_dhExchangePerms( jxf_SubdomainGraph_dh s );
extern void jxf_SubdomainGraph_dhPrintSubdomainGraph( jxf_SubdomainGraph_dh s, FILE *fp );
extern void jxf_SubdomainGraph_dhPrintStatsLong( jxf_SubdomainGraph_dh s, FILE *fp );
extern void jxf_SubdomainGraph_dhDump( jxf_SubdomainGraph_dh s, char *filename );
extern void jxf_SubdomainGraph_dhPrintRatios( jxf_SubdomainGraph_dh s, FILE *fp );
extern void jxf_SubdomainGraph_dhPrintStats( jxf_SubdomainGraph_dh sg, FILE *fp );

#endif

#ifndef JXF_TIMELOG_DH_DH
#define JXF_TIMELOG_DH_DH

extern void jxf_TimeLog_dhCreate( jxf_TimeLog_dh *t );
extern void jxf_TimeLog_dhDestroy( jxf_TimeLog_dh t );
extern void jxf_TimeLog_dhStart( jxf_TimeLog_dh t );
extern void jxf_TimeLog_dhStop( jxf_TimeLog_dh t );
extern void jxf_TimeLog_dhReset( jxf_TimeLog_dh t );
extern void jxf_TimeLog_dhMark( jxf_TimeLog_dh t, const char *description );
extern void jxf_TimeLog_dhPrint( jxf_TimeLog_dh t, FILE *fp, jxf_bool allPrint );

#endif

#ifndef JXF_SORTED_SET_DH
#define JXF_SORTED_SET_DH

struct _jxf_sortedset_dh
{
    JXF_Int n;
    JXF_Int *list;
    JXF_Int count;
};

extern void jxf_SortedSet_dhCreate( jxf_SortedSet_dh *ss, JXF_Int initialSize );
extern void jxf_SortedSet_dhDestroy( jxf_SortedSet_dh ss );
extern void jxf_SortedSet_dhInsert( jxf_SortedSet_dh ss, JXF_Int idx );
extern void jxf_SortedSet_dhGetList( jxf_SortedSet_dh ss, JXF_Int **list, JXF_Int *count );

#endif

#ifndef JXF_MEM_DH_DH
#define JXF_MEM_DH_DH

extern void jxf_Mem_dhCreate( jxf_Mem_dh *m );
extern void jxf_Mem_dhDestroy( jxf_Mem_dh m );
extern void *jxf_Mem_dhMalloc( jxf_Mem_dh m, size_t size );
extern void jxf_Mem_dhFree( jxf_Mem_dh m, void *ptr );
extern void jxf_Mem_dhPrint( jxf_Mem_dh m, FILE* fp, jxf_bool allPrint );

#endif

#ifndef JXF_SUPPORT_DH
#define JXF_SUPPORT_DH

extern void jxf_shellSort_int( const JXF_Int n, JXF_Int *x );
extern void jxf_shellSort_float( JXF_Int n, JXF_Real *v );

#endif

#ifndef JXF_NUMBERING_DH_H
#define JXF_NUMBERING_DH_H

struct _jxf_numbering_dh
{
    JXF_Int size;
    JXF_Int first;
    JXF_Int m;
    JXF_Int *idx_ext;
    JXF_Int *idx_extLo;
    JXF_Int *idx_extHi;
    JXF_Int num_ext;
    JXF_Int num_extLo;
    JXF_Int num_extHi;
    jxf_Hash_i_dh global_to_local;
    
    jxf_bool debug;
};

extern void jxf_Numbering_dhCreate( jxf_Numbering_dh *numb );
extern void jxf_Numbering_dhDestroy( jxf_Numbering_dh numb );
extern void jxf_Numbering_dhSetup( jxf_Numbering_dh numb, jxf_Mat_dh mat );
extern void jxf_Numbering_dhGlobalToLocal( jxf_Numbering_dh numb, JXF_Int len, JXF_Int *global_in, JXF_Int *local_out );

#endif

#ifndef JXF_HASH_I_DH
#define JXF_HASH_I_DH

extern void jxf_Hash_i_dhCreate( jxf_Hash_i_dh *h, JXF_Int size );
extern void jxf_Hash_i_dhDestroy( jxf_Hash_i_dh h );
extern void jxf_Hash_i_dhReset( jxf_Hash_i_dh h );
extern void jxf_Hash_i_dhInsert( jxf_Hash_i_dh h, JXF_Int key, JXF_Int data );
extern JXF_Int jxf_Hash_i_dhLookup( jxf_Hash_i_dh h, JXF_Int key );

#endif

#ifndef JXF_TIMER_DH_H
#define JXF_TIMER_DH_H

#include <time.h>
#include <unistd.h>

struct _jxf_timer_dh
{
    jxf_bool isRunning;
    jxf_longint sc_clk_tck;
    JXF_Real begin_wall;
    JXF_Real end_wall;
};

extern void jxf_Timer_dhCreate( jxf_Timer_dh *t );
extern void jxf_Timer_dhDestroy( jxf_Timer_dh t );
extern void jxf_Timer_dhStart( jxf_Timer_dh t );
extern void jxf_Timer_dhStop( jxf_Timer_dh t );
extern JXF_Real jxf_Timer_dhReadCPU( jxf_Timer_dh t );
extern JXF_Real jxf_Timer_dhReadWall( jxf_Timer_dh t );
extern JXF_Real jxf_Timer_dhReadUsage( jxf_Timer_dh t );

#endif

#ifndef JXF_PARSER_DH_DH
#define JXF_PARSER_DH_DH

extern void jxf_Parser_dhCreate( jxf_Parser_dh *p );
extern void jxf_Parser_dhDestroy( jxf_Parser_dh p );
extern jxf_bool jxf_Parser_dhHasSwitch( jxf_Parser_dh p, char *in );
extern jxf_bool jxf_Parser_dhReadString( jxf_Parser_dh p, char *in, char **out );
extern jxf_bool jxf_Parser_dhReadInt( jxf_Parser_dh p, char *in, JXF_Int *out );
extern jxf_bool jxf_Parser_dhReadDouble( jxf_Parser_dh p, char *in, JXF_Real *out );
extern void jxf_Parser_dhPrint( jxf_Parser_dh p, FILE *fp, jxf_bool allPrint );
extern void jxf_Parser_dhInsert( jxf_Parser_dh p, char *name, char *value );
extern void jxf_Parser_dhUpdateFromFile( jxf_Parser_dh p, char *name );
extern void jxf_Parser_dhInit( jxf_Parser_dh p, JXF_Int argc, char *argv[] );

#endif

#ifndef JXF_SORTEDLIST_DH_H
#define JXF_SORTEDLIST_DH_H

typedef struct _jxf_srecord
{
    JXF_Int col;
    JXF_Int level;
    JXF_Real val;
    JXF_Int next;
} jxf_SRecord;

extern void jxf_SortedList_dhCreate( jxf_SortedList_dh *sList );
extern void jxf_SortedList_dhDestroy( jxf_SortedList_dh sList );
extern void jxf_SortedList_dhInit( jxf_SortedList_dh sList, jxf_SubdomainGraph_dh sg );
extern void jxf_SortedList_dhEnforceConstraint( jxf_SortedList_dh sList, jxf_SubdomainGraph_dh sg );
extern void jxf_SortedList_dhReset( jxf_SortedList_dh sList, JXF_Int row );
extern JXF_Int jxf_SortedList_dhReadCount( jxf_SortedList_dh sList );
extern void jxf_SortedList_dhResetGetSmallest( jxf_SortedList_dh sList );
extern jxf_SRecord *jxf_SortedList_dhGetSmallest( jxf_SortedList_dh sList );
extern jxf_SRecord *jxf_SortedList_dhGetSmallestLowerTri( jxf_SortedList_dh sList );
extern void jxf_SortedList_dhInsert( jxf_SortedList_dh sList, jxf_SRecord *sr );
extern void jxf_SortedList_dhInsertOrUpdateVal( jxf_SortedList_dh sList, jxf_SRecord *sr );
extern jxf_bool jxf_SortedList_dhPermuteAndInsert( jxf_SortedList_dh sList, jxf_SRecord *sr, JXF_Real thresh );
extern void jxf_SortedList_dhInsertOrUpdate( jxf_SortedList_dh sList, jxf_SRecord *sr );
extern jxf_SRecord *jxf_SortedList_dhFind( jxf_SortedList_dh sList, jxf_SRecord *sr );
extern void jxf_SortedList_dhUpdateVal( jxf_SortedList_dh sList, jxf_SRecord *sr );

#endif

#ifndef JXF_HASH_D_DH
#define JXF_HASH_D_DH

typedef struct _jxf_hash_node
{
    JXF_Int iData;
    JXF_Real fData;
    JXF_Int *iDataPtr;
    JXF_Int *iDataPtr2;
    JXF_Real *fDataPtr;
} jxf_HashData;

typedef struct _jxf_hash_node_private jxf_HashRecord;

struct _jxf_hash_dh
{
    JXF_Int size;
    JXF_Int count;
    JXF_Int curMark;
    jxf_HashRecord *data;
};

extern void jxf_Hash_dhCreate( jxf_Hash_dh *h, JXF_Int size );
extern void jxf_Hash_dhDestroy( jxf_Hash_dh h );
extern void jxf_Hash_dhInsert( jxf_Hash_dh h, JXF_Int key, jxf_HashData *data );
extern jxf_HashData *jxf_Hash_dhLookup(jxf_Hash_dh h, JXF_Int key );
extern void jxf_Hash_dhReset( jxf_Hash_dh h );
extern void jxf_Hash_dhPrint( jxf_Hash_dh h, FILE *fp );

#define JXF_HASH_1(k,size,idxOut) \
      {  *idxOut = k % size;  }

#define JXF_HASH_2(k,size,idxOut) \
      { \
        JXF_Int r = k % (size-13); \
        r = (r % 2) ? r : r+1; \
       *idxOut = r; \
      }

#endif

#ifndef JXF_MAT_DH_PRIVATE
#define JXF_MAT_DH_PRIVATE

extern JXF_Int jxf_mat_jxf_find_owner( JXF_Int *beg_rows, JXF_Int *end_rows, JXF_Int index );
extern void jxf_mat_dh_transpose_private( JXF_Int m, JXF_Int *rpIN, JXF_Int **rpOUT,
                    JXF_Int *cvalIN, JXF_Int **cvalOUT, JXF_Real *avalIN, JXF_Real **avalOUT );
extern void jxf_mat_dh_transpose_reuse_private( JXF_Int m, JXF_Int *rpIN,
               JXF_Int *cvalIN, JXF_Real *avalIN, JXF_Int *rpOUT, JXF_Int *cvalOUT, JXF_Real *avalOUT );
extern void jxf_readMat( jxf_Mat_dh *Aout, char *fileType, char *fileName, JXF_Int ignore );
extern void jxf_readVec( jxf_Vec_dh *bout, char *fileType, char *fileName, JXF_Int ignore );
extern void jxf_writeMat( jxf_Mat_dh Ain, char *fileType, char *fileName );
extern void jxf_writeVec( jxf_Vec_dh b, char *fileType, char *fileName );
extern void jxf_readMat_par( jxf_Mat_dh *Aout, char *fileType, char *fileName, JXF_Int ignore );
extern void jxf_profileMat( jxf_Mat_dh A );
extern void jxf_mat_dh_print_graph_private( JXF_Int m, JXF_Int beg_row,
              JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, JXF_Int *n2o, JXF_Int *o2n, jxf_Hash_i_dh hash, FILE* fp );
extern void jxf_mat_dh_print_csr_private( JXF_Int m, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, FILE* fp );
extern void jxf_mat_dh_read_csr_private( JXF_Int *m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval, FILE* fp );
extern void jxf_mat_dh_read_triples_private( JXF_Int ignore, JXF_Int *m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval, FILE* fp );
extern void jxf_create_nat_ordering_private( JXF_Int m, JXF_Int **p );
extern void jxf_destroy_nat_ordering_private( JXF_Int *p );
extern void jxf_invert_perm( JXF_Int m, JXF_Int *pIN, JXF_Int *pOUT );
extern void jxf_make_full_private( JXF_Int m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval );
extern void jxf_make_symmetric_private( JXF_Int m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval );

#endif

#ifndef JXF_GET_ROW_DH
#define JXF_GET_ROW_DH

extern void jxf_EuclidGetDimensions( void *A, JXF_Int *beg_row, JXF_Int *rowsLocal, JXF_Int *rowsGlobal );
extern void jxf_EuclidGetRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val );
extern void jxf_EuclidRestoreRow( void *A, JXF_Int row, JXF_Int *len, JXF_Int **ind, JXF_Real **val );
extern JXF_Int jxf_EuclidReadLocalNz( void *A );
extern void jxf_PrintMatUsingGetRow( void *A, JXF_Int beg_row, JXF_Int m, JXF_Int *n2o_row, JXF_Int *n2o_col, char *filename );

#endif

#ifndef JXF_ILU_MPI_DH
#define JXF_ILU_MPI_DH

void jxf_reallocate_private( JXF_Int row, JXF_Int newEntries,
           JXF_Int *nzHave, JXF_Int **rp, JXF_Int **cval, float **aval, JXF_Real **avalD, JXF_Int **fill );
extern void jxf_ilu_mpi_pilu( jxf_Euclid_dh ctx );
extern void jxf_iluk_mpi_pilu( jxf_Euclid_dh ctx );
extern void jxf_compute_scaling_private( JXF_Int row, JXF_Int len, JXF_Real *AVAL, jxf_Euclid_dh ctx );
extern void jxf_iluk_mpi_bj( jxf_Euclid_dh ctx );
extern void jxf_iluk_seq( jxf_Euclid_dh ctx );
extern void jxf_iluk_seq_block( jxf_Euclid_dh ctx );
extern void jxf_ilut_seq( jxf_Euclid_dh ctx );

#endif

#ifndef JXF_EUCLID_MPI_INTERFACE_DH
#define JXF_EUCLID_MPI_INTERFACE_DH

#define JXF_DEFAULT_DROP_TOL 0.01

extern void jxf_Euclid_dhCreate( jxf_Euclid_dh *ctxOUT );
extern void jxf_Euclid_dhDestroy( jxf_Euclid_dh ctx );
extern void jxf_Euclid_dhSetup( jxf_Euclid_dh ctx );
extern void jxf_Euclid_dhSolve( jxf_Euclid_dh ctx, jxf_Vec_dh lhs, jxf_Vec_dh rhs, JXF_Int *its );
extern void jxf_Euclid_dhApply( jxf_Euclid_dh ctx, JXF_Real *lhs, JXF_Real *rhs );
extern void jxf_Euclid_dhPrintTestData( jxf_Euclid_dh ctx, FILE *fp );
extern void jxf_Euclid_dhPrintScaling( jxf_Euclid_dh ctx, FILE *fp );
extern void jxf_Euclid_dhPrintStatsShort( jxf_Euclid_dh ctx, JXF_Real setup, JXF_Real solve, FILE *fp );
extern void jxf_Euclid_dhPrintStatsShorter( jxf_Euclid_dh ctx, FILE *fp );
extern void jxf_Euclid_dhPrintJxpamgReport( jxf_Euclid_dh ctx, FILE *fp );
extern void jxf_Euclid_dhPrintStats( jxf_Euclid_dh ctx, FILE *fp );
extern void jxf_Euclid_dhInputJXFpamgMat( jxf_Euclid_dh ctx, JXF_ParCSRMatrix A );

#define JXF_MAX_OPT_LEN 20
#define JXF_TIMING_BINS 10

enum{ JXF_SOLVE_START_T, JXF_TRI_SOLVE_T, JXF_SETUP_T,
      JXF_SUB_GRAPH_T, JXF_FACTOR_T, JXF_SOLVE_SETUP_T,
      JXF_COMPUTE_RHO_T, JXF_TOTAL_SOLVE_TEMP_T, JXF_TOTAL_SOLVE_T };

#define JXF_STATS_BINS 10

enum{ JXF_NZA_STATS, JXF_NZF_STATS, JXF_NZA_USED_STATS, JXF_NZA_RATIO_STATS };

struct _jxf_mpi_interface_dh
{
    jxf_bool isSetup;
    
    JXF_Real rho_init;
    JXF_Real rho_final;
    
    JXF_Int m;
    JXF_Int n;
    JXF_Real *rhs;
    void *A;
    jxf_Factor_dh F;
    jxf_SubdomainGraph_dh sg;
    
    JXF_REAL_DH *scale;
    jxf_bool isScaled;
    
    JXF_Real *work;
    JXF_Real *work2;
    JXF_Int from, to;
    
    char algo_par[JXF_MAX_OPT_LEN];
    char algo_ilu[JXF_MAX_OPT_LEN];
    JXF_Int level;
    JXF_Real droptol;
    JXF_Real sparseTolA;
    JXF_Real sparseTolF;
    JXF_Real pivotMin;
    JXF_Real pivotFix;
    JXF_Real maxVal;
    
    jxf_SortedList_dh slist;
    jxf_ExternalRows_dh extRows;
    
    char krylovMethod[JXF_MAX_OPT_LEN];
    JXF_Int maxIts;
    JXF_Real rtol;
    JXF_Real atol;
    JXF_Int its;
    JXF_Int itsTotal;
    
    JXF_Int setupCount;
    JXF_Int logging;
    JXF_Real timing[JXF_TIMING_BINS];
    JXF_Real stats[JXF_STATS_BINS];
    jxf_bool timingsWereReduced;
    jxf_bool printStats;
};

#endif

#ifndef JXF_THREADED_KRYLOV_H
#define JXF_THREADED_KRYLOV_H

extern void jxf_bicgstab_euclid( jxf_Mat_dh A, jxf_Euclid_dh ctx, JXF_Real *x, JXF_Real *b, JXF_Int *itsOUT );
extern void jxf_cg_euclid( jxf_Mat_dh A, jxf_Euclid_dh ctx, JXF_Real *x, JXF_Real *b, JXF_Int *itsOUT );

#endif

#ifndef JXF_IO_DH
#define JXF_IO_DH

extern FILE *jxf_openFile_dh( const char *filenameIN, const char *modeIN );
extern void jxf_closeFile_dh( FILE *fpIN );
jxf_bool jxf_isSmallEndian();
extern void jxf_io_dh_print_ebin_mat_private( JXF_Int m, JXF_Int beg_row, JXF_Int *rp,
             JXF_Int *cval, JXF_Real *aval, JXF_Int *n2o, JXF_Int *o2n, jxf_Hash_i_dh hash, char *filename );
extern void jxf_io_dh_read_ebin_mat_private( JXF_Int *m, JXF_Int **rp, JXF_Int **cval, JXF_Real **aval, char *filename );
extern void jxf_io_dh_print_ebin_vec_private( JXF_Int n, JXF_Int beg_row, JXF_Real *vals,
                                JXF_Int *n2o, JXF_Int *o2n, jxf_Hash_i_dh hash, char *filename );
extern void jxf_io_dh_read_ebin_vec_private( JXF_Int *n, JXF_Real **vals, char *filename );

#endif

#ifndef JXF_THREADED_BLAS_DH
#define JXF_THREADED_BLAS_DH

extern void jxf_matvec_euclid_seq( JXF_Int n, JXF_Int *rp, JXF_Int *cval, JXF_Real *aval, JXF_Real *x, JXF_Real *y );
extern JXF_Real jxf_InnerProd( JXF_Int local_n, JXF_Real *x, JXF_Real *y );
extern void jxf_Axpy( JXF_Int n, JXF_Real alpha, JXF_Real *x, JXF_Real *y );
extern JXF_Real jxf_Norm2( JXF_Int n, JXF_Real *x );
extern void jxf_CopyVec( JXF_Int n, JXF_Real *xIN, JXF_Real *yOUT );
extern void jxf_ScaleVec( JXF_Int n, JXF_Real alpha, JXF_Real *x );

#endif

#endif
