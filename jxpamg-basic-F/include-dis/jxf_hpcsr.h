#ifndef JXF_HPCSRMV_HEADER 
#define JXF_HPCSRMV_HEADER


#ifndef JXF_MV_HEADER 
#include "jxf_mv.h"
#endif


#include "jxf_env_variable.h"


/*----------------------------------------------------------------*
 *                   hardware     Struct                          *
 *----------------------------------------------------------------*/
/*!
 * \struct jxf_hardwareinfo
 */
typedef struct
{
   JXF_Int         bind_type;  //0:bind to node, 1:bind to cpu, 2:bind to core
   JXF_Int         num_core_in_node;
   JXF_Int         num_core_in_cpu;
   JXF_Int         num_cpu_in_node;
   JXF_Int         proc_cpu_start;
   JXF_Int         proc_cpu_end;
   JXF_Int         proc_node_start;
   JXF_Int         proc_node_end;
   MPI_Comm       Comm;
   MPI_Comm       Core_Comm;
   MPI_Comm       Intra_Cpu_Comm;
   MPI_Comm       Node_Comm;
   // 注释保留原因: Allreduce优化中可能会用到
   // MPI_Comm       Cpu_Comm;
   // MPI_Comm       Intra_Node_Comm;
} jxf_HardwareInfo;

#define jxf_hpBindType(Hardware)	            ((Hardware) -> bind_type)
#define jxf_hpNumCoreInNode(Hardware)	      ((Hardware) -> num_core_in_node)
#define jxf_hpNumCoreInCpu(Hardware)	         ((Hardware) -> num_core_in_cpu)
#define jxf_hpNumCpuInNode(Hardware)	         ((Hardware) -> num_cpu_in_node)
#define jxf_hpProcCpuStart(Hardware)	         ((Hardware) -> proc_cpu_start)
#define jxf_hpProcCpuEnd(Hardware)	         ((Hardware) -> proc_cpu_end)
#define jxf_hpProcNodeStart(Hardware)	      ((Hardware) -> proc_node_start)
#define jxf_hpProcNodeEnd(Hardware)	         ((Hardware) -> proc_node_end)
#define jxf_hpComm(Hardware)                  ((Hardware) -> Comm)
#define jxf_hpCoreComm(Hardware)	            ((Hardware) -> Core_Comm)
#define jxf_hpIntraCpuComm(Hardware)	         ((Hardware) -> Intra_Cpu_Comm)
#define jxf_hpNodeComm(Hardware)	            ((Hardware) -> Node_Comm)
// 注释保留原因: Allreduce优化中可能会用到
// #define jxf_hpCpuComm(Hardware)	            ((Hardware) -> Cpu_Comm)
// #define jxf_hpIntraNodeComm(Hardware)	      ((Hardware) -> Intra_Node_Comm)

extern jxf_HardwareInfo *jxf_hp_hardware;
/*!
 * \struct jxf_hpCSRMatrix
 */
typedef struct
{
   // JXF_Int           level_colstartnum;
   // JXF_Int           level_colendnum;
   jxf_CSRMatrix     *offd_level;
   JXF_Int           *col_map_offd_level;
   jxf_ParCSRCommPkg *level_comm_pkg;
} jxf_hpCSRMatrixLevel;
// #define jxf_hpCSRMatrixlevelColStartNum(matrix)  ((matrix) -> level_colstartnum)
// #define jxf_hpCSRMatrixlevelColEndNum(matrix)    ((matrix) -> level_colendnum)
#define jxf_hpCSRMatrixOffdlevel(matrix)         ((matrix) -> offd_level)
#define jxf_hpCSRMatrixColMapOffdlevel(matrix)   ((matrix) -> col_map_offd_level)
#define jxf_hpCSRMatrixlevelCommPkg(matrix)      ((matrix) -> level_comm_pkg)

/*!
 * \struct jxf_hpCSRMatrix
 */
typedef struct
{
   jxf_ParCSRMatrix      *par_matrix;
   jxf_hpCSRMatrixLevel  *node_matrix;
   jxf_hpCSRMatrixLevel  *cpu_matrix;
   jxf_hpCSRMatrixLevel  *core_matrix;
   // jxf_HardwareInfo      *hardware;
} jxf_hpCSRMatrix;

#define jxf_hpCSRMatrixComm(matrix)	           ((matrix) -> par_matrix -> comm)
#define jxf_hpCSRMatrixGlobalNumRows(matrix)    ((matrix) -> par_matrix -> global_num_rows)
#define jxf_hpCSRMatrixGlobalNumCols(matrix)    ((matrix) -> par_matrix -> global_num_cols)
#define jxf_hpCSRMatrixFirstRowIndex(matrix)    ((matrix) -> par_matrix -> first_row_index)
#define jxf_hpCSRMatrixFirstColDiag(matrix)     ((matrix) -> par_matrix -> first_col_diag)
#define jxf_hpCSRMatrixLastRowIndex(matrix)     ((matrix) -> par_matrix -> last_row_index)
#define jxf_hpCSRMatrixLastColDiag(matrix)      ((matrix) -> par_matrix -> last_col_diag)
#define jxf_hpCSRMatrixDiag(matrix)  	        ((matrix) -> par_matrix -> diag)
#define jxf_hpCSRMatrixOffd(matrix)  	        ((matrix) -> par_matrix -> offd)
#define jxf_hpCSRMatrixDiagT(matrix)  	        ((matrix) -> par_matrix -> diagT)
#define jxf_hpCSRMatrixOffdT(matrix)  	        ((matrix) -> par_matrix -> offdT)
#define jxf_hpCSRMatrixColMapOffd(matrix)       ((matrix) -> par_matrix -> col_map_offd)
#define jxf_hpCSRMatrixRowStarts(matrix)        ((matrix) -> par_matrix -> row_starts)
#define jxf_hpCSRMatrixColStarts(matrix)        ((matrix) -> par_matrix -> col_starts)
#define jxf_hpCSRMatrixCommPkg(matrix)          ((matrix) -> par_matrix -> comm_pkg)
#define jxf_hpCSRMatrixCommPkgT(matrix)         ((matrix) -> par_matrix -> comm_pkgT)
#define jxf_hpCSRMatrixOwnsData(matrix)         ((matrix) -> par_matrix -> owns_data)
#define jxf_hpCSRMatrixOwnsRowStarts(matrix)    ((matrix) -> par_matrix -> owns_row_starts)
#define jxf_hpCSRMatrixOwnsColStarts(matrix)    ((matrix) -> par_matrix -> owns_col_starts)
#define jxf_hpCSRMatrixNumRows(matrix)          (jxf_CSRMatrixNumRows((matrix) -> par_matrix -> diag))
#define jxf_hpCSRMatrixNumCols(matrix)          (jxf_CSRMatrixNumCols((matrix) -> par_matrix-> diag))
#define jxf_hpCSRMatrixNumNonzeros(matrix)      ((matrix) -> par_matrix -> num_nonzeros)
#define jxf_hpCSRMatrixDNumNonzeros(matrix)     ((matrix) -> par_matrix -> d_num_nonzeros)
#define jxf_hpCSRMatrixRowindices(matrix)       ((matrix) -> par_matrix -> rowindices)
#define jxf_hpCSRMatrixRowvalues(matrix)        ((matrix) -> par_matrix -> rowvalues)
#define jxf_hpCSRMatrixGetrowactive(matrix)     ((matrix) -> par_matrix -> getrowactive)
#define jxf_hpCSRMatrixAssumedPartition(matrix) ((matrix) -> par_matrix -> assumed_partition)
//hpcsr add
#define jxf_hpCSRMatrixPar(matrix)	           ((matrix) -> par_matrix)
//cpu_matrix
#define jxf_hpCSRMatrixCpu(matrix)	           ((matrix) -> cpu_matrix)
#define jxf_hpCSRMatrixOffdCpu(matrix)         ((matrix) -> cpu_matrix -> offd_level)
#define jxf_hpCSRMatrixColMapOffdCpu(matrix)    ((matrix) -> cpu_matrix  -> col_map_offd_level)
#define jxf_hpCSRMatrixCpuCommPkg(matrix)       ((matrix) -> cpu_matrix  -> level_comm_pkg)
//node_matrix
#define jxf_hpCSRMatrixNode(matrix)	           ((matrix) -> node_matrix)
#define jxf_hpCSRMatrixOffdNode(matrix)         ((matrix) -> node_matrix -> offd_level)
#define jxf_hpCSRMatrixColMapOffdNode(matrix)   ((matrix) -> node_matrix -> col_map_offd_level)
#define jxf_hpCSRMatrixNodeCommPkg(matrix)      ((matrix) -> node_matrix -> level_comm_pkg)
//core_matrix
#define jxf_hpCSRMatrixCore(matrix)	           ((matrix) -> core_matrix)
#define jxf_hpCSRMatrixOffdCore(matrix)         ((matrix) -> core_matrix -> offd_level)
#define jxf_hpCSRMatrixColMapOffdCore(matrix)   ((matrix) -> core_matrix -> col_map_offd_level)
#define jxf_hpCSRMatrixCoreCommPkg(matrix)      ((matrix) -> core_matrix -> level_comm_pkg)


/*csrc/mat_vec_data/hpcsr_vector.c*/
jxf_ParVector *jxf_hpBuildRhsParFromOneFile( char *filename, jxf_hpCSRMatrix *A, JXF_Int file_type );

/*csrc/mat_vec_data/hpcsr_matrix.c*/
jxf_hpCSRMatrix *jxf_hpBuildMatParFromOneFile( char *filename, JXF_Int num_functions, JXF_Int file_base);
#ifdef USING_HWLOC
JXF_Int jxf_hpCreateHardwareInfo(MPI_Comm comm);
JXF_Int jxf_hpBindToType_MPI(jxf_HardwareInfo *hardware);
JXF_Int jxf_hpBindToType_Hwloc(jxf_HardwareInfo *hardware);
#endif
JXF_Int jxf_hpGetLevelBindId(JXF_Int myid, JXF_Int numprocess, JXF_Int num_level, 
                           JXF_Int *id_in_next_level, JXF_Int *process_id_in_level, JXF_Int *numprocess_in_level, JXF_Int* start_process_id_in_level, JXF_Int *end_process_id_in_level);
jxf_hpCSRMatrix * jxf_hpCSRMatrixCreate(MPI_Comm   comm,
                       JXF_Int        global_num_rows,
                       JXF_Int        global_num_cols,
                       JXF_Int       *row_starts,
                       JXF_Int       *col_starts,
                       JXF_Int        num_cols_offd,
                       JXF_Int        num_nonzeros_diag,
                       JXF_Int        num_nonzeros_offd);
JXF_Int jxf_hpCSRMatrixInitialize(jxf_hpCSRMatrix  *A);
jxf_hpCSRMatrixLevel *jxf_hpInitMatrixLevel();
JXF_Int jxf_hpCSRMatrixSetCpuLevel(jxf_HardwareInfo *hardware);
JXF_Int jxf_hpCSRMatrixSetNodeLevel(jxf_HardwareInfo *hardware);
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdNode(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware);
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdCpu(jxf_hpCSRMatrix *par_matrix, jxf_HardwareInfo * hardware);
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixNode(jxf_hpCSRMatrix *A);
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixCpu(jxf_hpCSRMatrix *A);
jxf_hpCSRMatrixLevel *jxf_hpCreateMatrixCore(jxf_hpCSRMatrix *A );
JXF_Int jxf_hpCSRMatrixDestroy(jxf_hpCSRMatrix *matrix);
JXF_Int jxf_hpCSRMatrixLevelDestroy(jxf_hpCSRMatrixLevel *matrix);
JXF_Int jxf_hpCSRMatrixDestroyAssumedPartition( jxf_hpCSRMatrix *matrix );
JXF_Int jxf_hpCSRMatrixPrint( jxf_hpCSRMatrix *matrix, const char *file_name );
JXF_Int jxf_hpCSRMatrixSetRowStartsOwner( jxf_hpCSRMatrix *matrix, JXF_Int owns_row_starts );
jxf_ParCSRMatrix * 
jxf_hpMatrixLevelToPar( jxf_hpCSRMatrix *hp_matrix,JXF_Int Level);
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdCpuOutside(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware);
jxf_hpCSRMatrixLevel *jxf_hpCSRMatrixGenerateOffdNodeOutside(jxf_hpCSRMatrix *hp_matrix, jxf_HardwareInfo *hardware);
jxf_hpCSRMatrix *jxf_hpCreateMatrixLevelBlock(jxf_hpCSRMatrix *A);
jxf_hpCSRMatrix  *jxf_hpInithpCSRMatrix();
JXF_Int jxf_hpCSRhardwareDestroy();

/*csrc/communication/hpcsr_comm.c*/
#if JXF_REODER_SEND_RECV
JXF_Int jxf_hpCSRCommCommPkgSendReorder(JXF_Int num_sends,jxf_HardwareInfo  *hardware,
                                       JXF_Int *SendProcs,  JXF_Int **SendProcsReorder, JXF_Int **send_reorder_map);
JXF_Int jxf_hpCSRCommCommPkgRecvReorder(JXF_Int num_recvs,jxf_HardwareInfo  *hardware, 
                                       JXF_Int *RecvProcs, JXF_Int **RecvProcsReorder, JXF_Int **recv_reorder_map);
JXF_Int jxf_hpCSRCommPkgReorder(jxf_ParCSRCommPkg *comm_pkg);
#endif

/*csrc/operation/mvops/hpcsr_matvec.c*/
JXF_Int jxf_hpCSRMatrixMatvec( JXF_Real alpha, jxf_hpCSRMatrix *A, jxf_ParVector*x, JXF_Real beta, jxf_ParVector *y);
JXF_Int jxf_hpCSRMatrixMatvecT( JXF_Real alpha, jxf_hpCSRMatrix *A, jxf_ParVector*x, JXF_Real beta, jxf_ParVector *y);
JXF_Int
jxf_hpCSRMatrixMatvecLevel( JXF_Real           alpha,
              	       jxf_hpCSRMatrix *hp_A,
                       jxf_ParVector    *x,
                       JXF_Real           beta,
                       jxf_ParVector    *y,
                       JXF_Int Level);

/*csrc/operation/matops/hpcsr_matop.c*/
JXF_Int 
jxf_hpCSRMatrixCopy( jxf_hpCSRMatrix *A, 
                     jxf_hpCSRMatrix *B, 
                     JXF_Int              copy_data );

/* csrc/operation/relaxations/par_relax_0.c */
JXF_Int
jxf_hpPAMGRelax0( jxf_hpCSRMatrix *par_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp);
JXF_Int
jxf_hpPAMGRelaxAI0( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_1.c */
JXF_Int  
jxf_hpPAMGRelax1( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI1( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_2.c */
JXF_Int  
jxf_hpPAMGRelax2( jxf_hpCSRMatrix *par_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI2( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_3.c */
JXF_Int  
jxf_hpPAMGRelax3( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI3( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_4.c */
JXF_Int  
jxf_hpPAMGRelax4( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int
jxf_hpPAMGRelaxAI4( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_5.c */
JXF_Int  
jxf_hpPAMGRelax5( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI5( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_6.c */
JXF_Int  
jxf_hpPAMGRelax6( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI6( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_7.c */
JXF_Int  
jxf_hpPAMGRelax7( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
JXF_Int  
jxf_hpPAMGRelaxAI7( jxf_hpCSRMatrix *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 JXF_Int             *cf_marker,
                 JXF_Int             *relax_marker,
                 JXF_Int              relax_points,
                 JXF_Real           relax_weight,
                 JXF_Real           omega,
                 jxf_ParVector    *par_app,
                 jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_8.c */
JXF_Int
jxf_hpPAMGRelax8( jxf_hpCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp );
/* csrc/operation/relaxations/par_relax_9.c */
JXF_Int  
jxf_hpPAMGRelax9( jxf_hpCSRMatrix *hp_matrix,
               jxf_ParVector    *par_rhs,
               JXF_Int             *cf_marker,
               JXF_Int              relax_points,
               JXF_Real           relax_weight,
               JXF_Real           omega,
               jxf_ParVector    *par_app,
               jxf_ParVector    *Vtemp );
/* csrc/operation/relaxations/par_relax_13.c */
JXF_Int
jxf_hpPAMGRelax13( jxf_hpCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp );
/* csrc/operation/relaxations/par_relax_14.c */
JXF_Int
jxf_hpPAMGRelax14( jxf_hpCSRMatrix *A,
               jxf_ParVector *f,
               JXF_Int *cf_marker,
               JXF_Int relax_points,
               JXF_Real relax_weight,
               JXF_Real omega,
               JXF_Real *l1_norms,
               jxf_ParVector *u,
               jxf_ParVector *Vtemp,
               jxf_ParVector *Ztemp );
/* utilities/mpistubs.c */
JXF_Int
jxf_My_MPI_Allreduce( void *sendbuf, 
                     void *recvbuf, 
                     JXF_Int count, 
                     MPI_Datatype datatype, 
                     MPI_Op op, 
                     MPI_Comm comm );
JXF_Int
jxf_my_bcast(void *sendbuf,JXF_Int count, MPI_Datatype datatype, JXF_Int root, MPI_Comm comm);
MPI_Request *
jxf_My_MPI_Iallreduce_First( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
JXF_Int
jxf_My_MPI_Iallreduce_Sencond(void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *req);
#endif

