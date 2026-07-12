#ifndef JX_HPCSRMV_HEADER 
#define JX_HPCSRMV_HEADER


#ifndef JX_MV_HEADER 
#include "jx_mv.h"
#endif


#include "jx_env_variable.h"


/*----------------------------------------------------------------*
 *                   hardware     Struct                          *
 *----------------------------------------------------------------*/
/*!
 * \struct jx_hardwareinfo
 */
typedef struct
{
   JX_Int         bind_type;  //0:bind to node, 1:bind to cpu, 2:bind to core
   JX_Int         num_core_in_node;
   JX_Int         num_core_in_cpu;
   JX_Int         num_cpu_in_node;
   JX_Int         proc_cpu_start;
   JX_Int         proc_cpu_end;
   JX_Int         proc_node_start;
   JX_Int         proc_node_end;
   MPI_Comm       Comm;
   MPI_Comm       Core_Comm;
   MPI_Comm       Intra_Cpu_Comm;
   MPI_Comm       Node_Comm;
   // 注释保留原因: Allreduce优化中可能会用到
   // MPI_Comm       Cpu_Comm;
   // MPI_Comm       Intra_Node_Comm;
} jx_HardwareInfo;

#define jx_hpBindType(Hardware)	            ((Hardware) -> bind_type)
#define jx_hpNumCoreInNode(Hardware)	      ((Hardware) -> num_core_in_node)
#define jx_hpNumCoreInCpu(Hardware)	         ((Hardware) -> num_core_in_cpu)
#define jx_hpNumCpuInNode(Hardware)	         ((Hardware) -> num_cpu_in_node)
#define jx_hpProcCpuStart(Hardware)	         ((Hardware) -> proc_cpu_start)
#define jx_hpProcCpuEnd(Hardware)	         ((Hardware) -> proc_cpu_end)
#define jx_hpProcNodeStart(Hardware)	      ((Hardware) -> proc_node_start)
#define jx_hpProcNodeEnd(Hardware)	         ((Hardware) -> proc_node_end)
#define jx_hpComm(Hardware)                  ((Hardware) -> Comm)
#define jx_hpCoreComm(Hardware)	            ((Hardware) -> Core_Comm)
#define jx_hpIntraCpuComm(Hardware)	         ((Hardware) -> Intra_Cpu_Comm)
#define jx_hpNodeComm(Hardware)	            ((Hardware) -> Node_Comm)
// 注释保留原因: Allreduce优化中可能会用到
// #define jx_hpCpuComm(Hardware)	            ((Hardware) -> Cpu_Comm)
// #define jx_hpIntraNodeComm(Hardware)	      ((Hardware) -> Intra_Node_Comm)

extern jx_HardwareInfo *hp_hardware;
/*!
 * \struct jx_hpCSRMatrix
 */
typedef struct
{
   // JX_Int           level_colstartnum;
   // JX_Int           level_colendnum;
   jx_CSRMatrix     *offd_level;
   JX_Int           *col_map_offd_level;
   jx_ParCSRCommPkg *level_comm_pkg;
} jx_hpCSRMatrixLevel;
// #define jx_hpCSRMatrixlevelColStartNum(matrix)  ((matrix) -> level_colstartnum)
// #define jx_hpCSRMatrixlevelColEndNum(matrix)    ((matrix) -> level_colendnum)
#define jx_hpCSRMatrixOffdlevel(matrix)         ((matrix) -> offd_level)
#define jx_hpCSRMatrixColMapOffdlevel(matrix)   ((matrix) -> col_map_offd_level)
#define jx_hpCSRMatrixlevelCommPkg(matrix)      ((matrix) -> level_comm_pkg)

/*!
 * \struct jx_hpCSRMatrix
 */
typedef struct
{
   jx_ParCSRMatrix      *par_matrix;
   jx_hpCSRMatrixLevel  *node_matrix;
   jx_hpCSRMatrixLevel  *cpu_matrix;
   jx_hpCSRMatrixLevel  *core_matrix;
   // jx_HardwareInfo      *hardware;
} jx_hpCSRMatrix;

#define jx_hpCSRMatrixComm(matrix)	           ((matrix) -> par_matrix -> comm)
#define jx_hpCSRMatrixGlobalNumRows(matrix)    ((matrix) -> par_matrix -> global_num_rows)
#define jx_hpCSRMatrixGlobalNumCols(matrix)    ((matrix) -> par_matrix -> global_num_cols)
#define jx_hpCSRMatrixFirstRowIndex(matrix)    ((matrix) -> par_matrix -> first_row_index)
#define jx_hpCSRMatrixFirstColDiag(matrix)     ((matrix) -> par_matrix -> first_col_diag)
#define jx_hpCSRMatrixLastRowIndex(matrix)     ((matrix) -> par_matrix -> last_row_index)
#define jx_hpCSRMatrixLastColDiag(matrix)      ((matrix) -> par_matrix -> last_col_diag)
#define jx_hpCSRMatrixDiag(matrix)  	        ((matrix) -> par_matrix -> diag)
#define jx_hpCSRMatrixOffd(matrix)  	        ((matrix) -> par_matrix -> offd)
#define jx_hpCSRMatrixDiagT(matrix)  	        ((matrix) -> par_matrix -> diagT)
#define jx_hpCSRMatrixOffdT(matrix)  	        ((matrix) -> par_matrix -> offdT)
#define jx_hpCSRMatrixColMapOffd(matrix)       ((matrix) -> par_matrix -> col_map_offd)
#define jx_hpCSRMatrixRowStarts(matrix)        ((matrix) -> par_matrix -> row_starts)
#define jx_hpCSRMatrixColStarts(matrix)        ((matrix) -> par_matrix -> col_starts)
#define jx_hpCSRMatrixCommPkg(matrix)          ((matrix) -> par_matrix -> comm_pkg)
#define jx_hpCSRMatrixCommPkgT(matrix)         ((matrix) -> par_matrix -> comm_pkgT)
#define jx_hpCSRMatrixOwnsData(matrix)         ((matrix) -> par_matrix -> owns_data)
#define jx_hpCSRMatrixOwnsRowStarts(matrix)    ((matrix) -> par_matrix -> owns_row_starts)
#define jx_hpCSRMatrixOwnsColStarts(matrix)    ((matrix) -> par_matrix -> owns_col_starts)
#define jx_hpCSRMatrixNumRows(matrix)          (jx_CSRMatrixNumRows((matrix) -> par_matrix -> diag))
#define jx_hpCSRMatrixNumCols(matrix)          (jx_CSRMatrixNumCols((matrix) -> par_matrix-> diag))
#define jx_hpCSRMatrixNumNonzeros(matrix)      ((matrix) -> par_matrix -> num_nonzeros)
#define jx_hpCSRMatrixDNumNonzeros(matrix)     ((matrix) -> par_matrix -> d_num_nonzeros)
#define jx_hpCSRMatrixRowindices(matrix)       ((matrix) -> par_matrix -> rowindices)
#define jx_hpCSRMatrixRowvalues(matrix)        ((matrix) -> par_matrix -> rowvalues)
#define jx_hpCSRMatrixGetrowactive(matrix)     ((matrix) -> par_matrix -> getrowactive)
#define jx_hpCSRMatrixAssumedPartition(matrix) ((matrix) -> par_matrix -> assumed_partition)
//hpcsr add
#define jx_hpCSRMatrixPar(matrix)	           ((matrix) -> par_matrix)
//cpu_matrix
#define jx_hpCSRMatrixCpu(matrix)	           ((matrix) -> cpu_matrix)
#define jx_hpCSRMatrixOffdCpu(matrix)         ((matrix) -> cpu_matrix -> offd_level)
#define jx_hpCSRMatrixColMapOffdCpu(matrix)    ((matrix) -> cpu_matrix  -> col_map_offd_level)
#define jx_hpCSRMatrixCpuCommPkg(matrix)       ((matrix) -> cpu_matrix  -> level_comm_pkg)
//node_matrix
#define jx_hpCSRMatrixNode(matrix)	           ((matrix) -> node_matrix)
#define jx_hpCSRMatrixOffdNode(matrix)         ((matrix) -> node_matrix -> offd_level)
#define jx_hpCSRMatrixColMapOffdNode(matrix)   ((matrix) -> node_matrix -> col_map_offd_level)
#define jx_hpCSRMatrixNodeCommPkg(matrix)      ((matrix) -> node_matrix -> level_comm_pkg)
//core_matrix
#define jx_hpCSRMatrixCore(matrix)	           ((matrix) -> core_matrix)
#define jx_hpCSRMatrixOffdCore(matrix)         ((matrix) -> core_matrix -> offd_level)
#define jx_hpCSRMatrixColMapOffdCore(matrix)   ((matrix) -> core_matrix -> col_map_offd_level)
#define jx_hpCSRMatrixCoreCommPkg(matrix)      ((matrix) -> core_matrix -> level_comm_pkg)


/*csrc/mat_vec_data/hpcsr_vector.c*/
jx_ParVector *jx_hpBuildRhsParFromOneFile( char *filename, jx_hpCSRMatrix *A, JX_Int file_type );

/*csrc/mat_vec_data/hpcsr_matrix.c*/
jx_hpCSRMatrix *jx_hpBuildMatParFromOneFile( char *filename, JX_Int num_functions, JX_Int file_base);
#ifdef USING_HWLOC
JX_Int jx_hpCreateHardwareInfo(MPI_Comm comm);
JX_Int jx_hpBindToType_MPI(jx_HardwareInfo *hardware);
JX_Int jx_hpBindToType_Hwloc(jx_HardwareInfo *hardware);
#endif
JX_Int jx_hpGetLevelBindId(JX_Int myid, JX_Int numprocess, JX_Int num_level, 
                           JX_Int *id_in_next_level, JX_Int *process_id_in_level, JX_Int *numprocess_in_level, JX_Int* start_process_id_in_level, JX_Int *end_process_id_in_level);
jx_hpCSRMatrix * jx_hpCSRMatrixCreate(MPI_Comm   comm,
                       JX_Int        global_num_rows,
                       JX_Int        global_num_cols,
                       JX_Int       *row_starts,
                       JX_Int       *col_starts,
                       JX_Int        num_cols_offd,
                       JX_Int        num_nonzeros_diag,
                       JX_Int        num_nonzeros_offd);
JX_Int jx_hpCSRMatrixInitialize(jx_hpCSRMatrix  *A);
jx_hpCSRMatrixLevel *jx_hpInitMatrixLevel();
JX_Int jx_hpCSRMatrixSetCpuLevel(jx_HardwareInfo *hardware);
JX_Int jx_hpCSRMatrixSetNodeLevel(jx_HardwareInfo *hardware);
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdNode(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware);
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdCpu(jx_hpCSRMatrix *par_matrix, jx_HardwareInfo * hardware);
jx_hpCSRMatrixLevel *jx_hpCreateMatrixNode(jx_hpCSRMatrix *A);
jx_hpCSRMatrixLevel *jx_hpCreateMatrixCpu(jx_hpCSRMatrix *A);
jx_hpCSRMatrixLevel *jx_hpCreateMatrixCore(jx_hpCSRMatrix *A );
JX_Int jx_hpCSRMatrixDestroy(jx_hpCSRMatrix *matrix);
JX_Int jx_hpCSRMatrixLevelDestroy(jx_hpCSRMatrixLevel *matrix);
JX_Int jx_hpCSRMatrixDestroyAssumedPartition( jx_hpCSRMatrix *matrix );
JX_Int jx_hpCSRMatrixPrint( jx_hpCSRMatrix *matrix, const char *file_name );
JX_Int jx_hpCSRMatrixSetRowStartsOwner( jx_hpCSRMatrix *matrix, JX_Int owns_row_starts );
jx_ParCSRMatrix * 
jx_hpMatrixLevelToPar( jx_hpCSRMatrix *hp_matrix,JX_Int Level);
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdCpuOutside(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware);
jx_hpCSRMatrixLevel *jx_hpCSRMatrixGenerateOffdNodeOutside(jx_hpCSRMatrix *hp_matrix, jx_HardwareInfo *hardware);
jx_hpCSRMatrix *jx_hpCreateMatrixLevelBlock(jx_hpCSRMatrix *A);
jx_hpCSRMatrix  *jx_hpInithpCSRMatrix();
JX_Int jx_hpCSRhardwareDestroy();

/*csrc/communication/hpcsr_comm.c*/
#if JX_REODER_SEND_RECV
JX_Int jx_hpCSRCommCommPkgSendReorder(JX_Int num_sends,jx_HardwareInfo  *hardware,
                                       JX_Int *SendProcs,  JX_Int **SendProcsReorder, JX_Int **send_reorder_map);
JX_Int jx_hpCSRCommCommPkgRecvReorder(JX_Int num_recvs,jx_HardwareInfo  *hardware, 
                                       JX_Int *RecvProcs, JX_Int **RecvProcsReorder, JX_Int **recv_reorder_map);
JX_Int jx_hpCSRCommPkgReorder(jx_ParCSRCommPkg *comm_pkg);
#endif

/*csrc/operation/mvops/hpcsr_matvec.c*/
JX_Int jx_hpCSRMatrixMatvec( JX_Real alpha, jx_hpCSRMatrix *A, jx_ParVector*x, JX_Real beta, jx_ParVector *y);
JX_Int jx_hpCSRMatrixMatvecT( JX_Real alpha, jx_hpCSRMatrix *A, jx_ParVector*x, JX_Real beta, jx_ParVector *y);
JX_Int
jx_hpCSRMatrixMatvecLevel( JX_Real           alpha,
              	       jx_hpCSRMatrix *hp_A,
                       jx_ParVector    *x,
                       JX_Real           beta,
                       jx_ParVector    *y,
                       JX_Int Level);

/*csrc/operation/matops/hpcsr_matop.c*/
JX_Int 
jx_hpCSRMatrixCopy( jx_hpCSRMatrix *A, 
                     jx_hpCSRMatrix *B, 
                     JX_Int              copy_data );

/* csrc/operation/relaxations/par_relax_0.c */
JX_Int
jx_hpPAMGRelax0( jx_hpCSRMatrix *par_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp);
JX_Int
jx_hpPAMGRelaxAI0( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_1.c */
JX_Int  
jx_hpPAMGRelax1( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI1( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_2.c */
JX_Int  
jx_hpPAMGRelax2( jx_hpCSRMatrix *par_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI2( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_3.c */
JX_Int  
jx_hpPAMGRelax3( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI3( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_4.c */
JX_Int  
jx_hpPAMGRelax4( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int
jx_hpPAMGRelaxAI4( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_5.c */
JX_Int  
jx_hpPAMGRelax5( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI5( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_6.c */
JX_Int  
jx_hpPAMGRelax6( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI6( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_7.c */
JX_Int  
jx_hpPAMGRelax7( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
JX_Int  
jx_hpPAMGRelaxAI7( jx_hpCSRMatrix *hp_matrix,
                 jx_ParVector    *par_rhs,
                 JX_Int             *cf_marker,
                 JX_Int             *relax_marker,
                 JX_Int              relax_points,
                 JX_Real           relax_weight,
                 JX_Real           omega,
                 jx_ParVector    *par_app,
                 jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_8.c */
JX_Int
jx_hpPAMGRelax8( jx_hpCSRMatrix *A,
               jx_ParVector *f,
               JX_Int *cf_marker,
               JX_Int relax_points,
               JX_Real relax_weight,
               JX_Real omega,
               JX_Real *l1_norms,
               jx_ParVector *u,
               jx_ParVector *Vtemp,
               jx_ParVector *Ztemp );
/* csrc/operation/relaxations/par_relax_9.c */
JX_Int  
jx_hpPAMGRelax9( jx_hpCSRMatrix *hp_matrix,
               jx_ParVector    *par_rhs,
               JX_Int             *cf_marker,
               JX_Int              relax_points,
               JX_Real           relax_weight,
               JX_Real           omega,
               jx_ParVector    *par_app,
               jx_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_13.c */
JX_Int
jx_hpPAMGRelax13( jx_hpCSRMatrix *A,
               jx_ParVector *f,
               JX_Int *cf_marker,
               JX_Int relax_points,
               JX_Real relax_weight,
               JX_Real omega,
               JX_Real *l1_norms,
               jx_ParVector *u,
               jx_ParVector *Vtemp,
               jx_ParVector *Ztemp );
/* csrc/operation/relaxations/par_relax_14.c */
JX_Int
jx_hpPAMGRelax14( jx_hpCSRMatrix *A,
               jx_ParVector *f,
               JX_Int *cf_marker,
               JX_Int relax_points,
               JX_Real relax_weight,
               JX_Real omega,
               JX_Real *l1_norms,
               jx_ParVector *u,
               jx_ParVector *Vtemp,
               jx_ParVector *Ztemp );
/* utilities/mpistubs.c */
JX_Int
jx_My_MPI_Allreduce( void *sendbuf, 
                     void *recvbuf, 
                     JX_Int count, 
                     MPI_Datatype datatype, 
                     MPI_Op op, 
                     MPI_Comm comm );
JX_Int
jx_my_bcast(void *sendbuf,JX_Int count, MPI_Datatype datatype, JX_Int root, MPI_Comm comm);
MPI_Request *
jx_My_MPI_Iallreduce_First( void *sendbuf, void *recvbuf, JX_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
JX_Int
jx_My_MPI_Iallreduce_Sencond(void *recvbuf, JX_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *req);
#endif

