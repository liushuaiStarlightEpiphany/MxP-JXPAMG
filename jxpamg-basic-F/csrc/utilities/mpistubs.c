//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

#ifdef JXF_SEQUENTIAL

JXF_Int
jxf_MPI_Init( int *argc, char ***argv )
{
   return(0);
}

JXF_Int
jxf_MPI_Finalize()
{
   return(0);
}

JXF_Int
jxf_MPI_Abort( MPI_Comm comm, JXF_Int errorcode )
{
   return(0);
}

JXF_Real
jxf_MPI_Wtime()
{
   return(0.0);
}

JXF_Real
jxf_MPI_Wtick()
{
   return(0.0);
}

JXF_Int
jxf_MPI_Barrier( MPI_Comm comm )
{
   return(0);
}

JXF_Int
jxf_MPI_Comm_create( MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm )
{
   return(0);
}

JXF_Int
jxf_MPI_Comm_dup( MPI_Comm comm, MPI_Comm *newcomm )
{
   return(0);
}

JXF_Int
jxf_MPI_Comm_size( MPI_Comm comm, JXF_Int *size )
{
   *size = 1;
   return(0);
}

JXF_Int
jxf_MPI_Comm_rank( MPI_Comm comm, JXF_Int *rank )
{
   *rank = 0;
   return(0);
}

JXF_Int
jxf_MPI_Comm_free( MPI_Comm *comm )
{
   return 0;
}

JXF_Int
jxf_MPI_Comm_group( MPI_Comm comm, MPI_Group *group )
{
   return(0);
}

JXF_Int
jxf_MPI_Comm_split( MPI_Comm comm, JXF_Int n, JXF_Int m, MPI_Comm *comms )
{
   return(0);
}

JXF_Int
jxf_MPI_Group_incl( MPI_Group group, JXF_Int n, JXF_Int *ranks, MPI_Group *newgroup )
{
   return(0);
}

JXF_Int
jxf_MPI_Group_free( MPI_Group *group )
{
   return 0;
}

JXF_Int
jxf_MPI_Address( void *location, MPI_Aint *address )
{
   return(0);
}

JXF_Int
jxf_MPI_Get_count( MPI_Status *status, MPI_Datatype datatype, JXF_Int *count )
{
   return(0);
}

JXF_Int
jxf_MPI_Alltoall( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm )
{
   return(0);
}

JXF_Int
jxf_MPI_Allgather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm )
{
   JXF_Int i;

   switch (sendtype)
   {
      case JXF_MPI_INT:
      {
         JXF_Int *crecvbuf = (JXF_Int *)recvbuf;
         JXF_Int *csendbuf = (JXF_Int *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      }
      break;

      case JXF_MPI_REAL:
      {
         JXF_Real *crecvbuf = (JXF_Real *)recvbuf;
         JXF_Real *csendbuf = (JXF_Real *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      }
      break;
   }

   return(0);
}

JXF_Int
jxf_MPI_Allgatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, MPI_Comm comm )
{
   return ( jxf_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, *recvcounts, recvtype, comm) );
}

JXF_Int
jxf_MPI_Gather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return ( jxf_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm) );
}

JXF_Int
jxf_MPI_Gatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return ( jxf_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, *recvcounts, recvtype, comm) );
}

JXF_Int
jxf_MPI_Scatter( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return ( jxf_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm) );
}

JXF_Int
jxf_MPI_Scatterv( void *sendbuf, JXF_Int *sendcounts, JXF_Int *displs, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return ( jxf_MPI_Allgather(sendbuf, *sendcounts, sendtype, recvbuf, recvcount, recvtype, comm) );
}

JXF_Int
jxf_MPI_Bcast( void *buffer, JXF_Int count, MPI_Datatype datatype, JXF_Int root, MPI_Comm comm )
{
   return(0);
}

JXF_Int
jxf_MPI_Send( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm )
{
   return(0);
}

JXF_Int
jxf_MPI_Recv( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status )
{
   return(0);
}

JXF_Int
jxf_MPI_Isend( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return(0);
}

JXF_Int
jxf_MPI_Irecv( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return(0);
}

JXF_Int
jxf_MPI_Send_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return 0;
}

JXF_Int
jxf_MPI_Recv_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return 0;
}

JXF_Int
jxf_MPI_Irsend( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return 0;
}

JXF_Int
jxf_MPI_Startall( JXF_Int count, MPI_Request *array_of_requests )
{
   return 0;
}

JXF_Int
jxf_MPI_Probe( JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status )
{
   return 0;
}

JXF_Int
jxf_MPI_Iprobe( JXF_Int source, JXF_Int tag, MPI_Comm comm, JXF_Int *flag, MPI_Status *status )
{
   return 0;
}

JXF_Int
jxf_MPI_Test( MPI_Request *request, JXF_Int *flag, MPI_Status *status )
{
   *flag = 1;
   return(0);
}

JXF_Int
jxf_MPI_Testall( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *flag, MPI_Status *array_of_statuses )
{
   *flag = 1;
   return(0);
}

JXF_Int
jxf_MPI_Wait( MPI_Request *request, MPI_Status *status )
{
   return(0);
}

JXF_Int
jxf_MPI_Waitall( JXF_Int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses )
{
   return(0);
}

JXF_Int
jxf_MPI_Waitany( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *index, MPI_Status *status )
{
   return(0);
}

JXF_Int
jxf_MPI_Allreduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
   JXF_Int i;

   switch (datatype)
   {
      case JXF_MPI_INT:
      {
         JXF_Int *crecvbuf = (JXF_Int *)recvbuf;
         JXF_Int *csendbuf = (JXF_Int *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      }
      break;

      case JXF_MPI_REAL:
      {
         JXF_Real *crecvbuf = (JXF_Real *)recvbuf;
         JXF_Real *csendbuf = (JXF_Real *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      }
      break;
   }

   return 0;
}

JXF_Int
jxf_MPI_Reduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, JXF_Int root, MPI_Comm comm )
{
   jxf_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
   return 0;
}

JXF_Int
jxf_MPI_Scan( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
   jxf_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
   return 0;
}

JXF_Int
jxf_MPI_Request_free( MPI_Request *request )
{
   return 0;
}

JXF_Int
jxf_MPI_Type_contiguous( JXF_Int count, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return(0);
}

JXF_Int
jxf_MPI_Type_vector( JXF_Int count, JXF_Int blocklength, JXF_Int stride, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return(0);
}

JXF_Int
jxf_MPI_Type_hvector( JXF_Int count, JXF_Int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return(0);
}

JXF_Int
jxf_MPI_Type_struct( JXF_Int count, JXF_Int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype )
{
   return(0);
}

JXF_Int
jxf_MPI_Type_commit( MPI_Datatype *datatype )
{
   return(0);
}

JXF_Int
jxf_MPI_Type_free( MPI_Datatype *datatype )
{
   return(0);
}

JXF_Int
jxf_MPI_Op_create( MPI_User_function *function, int commute, MPI_Op *op )
{
   return(0);
}

JXF_Int
jxf_MPI_Op_free( MPI_Op *op )
{
   return(0);
}

#else

JXF_Int
jxf_MPI_Init( int *argc, char ***argv )
{
   return (JXF_Int) MPI_Init(argc, argv);
}

JXF_Int
jxf_MPI_Finalize()
{
   return (JXF_Int) MPI_Finalize();
}

JXF_Int
jxf_MPI_Abort( MPI_Comm comm, JXF_Int errorcode )
{
   return (JXF_Int) MPI_Abort(comm, (int)errorcode);
}

JXF_Real
jxf_MPI_Wtime()
{
   return MPI_Wtime();
}

JXF_Real
jxf_MPI_Wtick()
{
   return MPI_Wtick();
}

JXF_Int
jxf_MPI_Barrier( MPI_Comm comm )
{
   return (JXF_Int) MPI_Barrier(comm);
}

JXF_Int
jxf_MPI_Comm_create( MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm )
{
   return (JXF_Int) MPI_Comm_create(comm, group, newcomm);
}

JXF_Int
jxf_MPI_Comm_dup( MPI_Comm comm, MPI_Comm *newcomm )
{
   return (JXF_Int) MPI_Comm_dup(comm, newcomm);
}

JXF_Int
jxf_MPI_Comm_size( MPI_Comm comm, JXF_Int *size )
{
   int mpi_size;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Comm_size(comm, &mpi_size);
   *size = (JXF_Int) mpi_size;
   return ierr;
}

JXF_Int
jxf_MPI_Comm_rank( MPI_Comm comm, JXF_Int *rank )
{ 
   int mpi_rank;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Comm_rank(comm, &mpi_rank);
   *rank = (JXF_Int) mpi_rank;
   return ierr;
}

JXF_Int
jxf_MPI_Comm_free( MPI_Comm *comm )
{
   return (JXF_Int) MPI_Comm_free(comm);
}

JXF_Int
jxf_MPI_Comm_group( MPI_Comm comm, MPI_Group *group )
{
   return (JXF_Int) MPI_Comm_group(comm, group);
}

JXF_Int
jxf_MPI_Comm_split( MPI_Comm comm, JXF_Int n, JXF_Int m, MPI_Comm *comms )
{
   return (JXF_Int) MPI_Comm_split(comm, (int)n, (int)m, comms);
}

JXF_Int
jxf_MPI_Group_incl( MPI_Group group, JXF_Int n, JXF_Int *ranks, MPI_Group *newgroup )
{
   int *mpi_ranks;
   JXF_Int  i;
   JXF_Int  ierr;

   mpi_ranks = jxf_TAlloc(int, n);
   for (i = 0; i < n; i++)
   {
      mpi_ranks[i] = (int) ranks[i];
   }
   ierr = (JXF_Int) MPI_Group_incl(group, (int)n, mpi_ranks, newgroup);
   jxf_TFree(mpi_ranks);

   return ierr;
}

JXF_Int
jxf_MPI_Group_free( MPI_Group *group )
{
   return (JXF_Int) MPI_Group_free(group);
}

JXF_Int
jxf_MPI_Address( void *location, MPI_Aint *address )
{
   return (JXF_Int) MPI_Get_address(location, address);
}

JXF_Int
jxf_MPI_Get_count( MPI_Status *status, MPI_Datatype datatype, JXF_Int *count )
{
   int mpi_count;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Get_count(status, datatype, &mpi_count);
   *count = (JXF_Int) mpi_count;
   return ierr;
}

JXF_Int
jxf_MPI_Alltoall( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm )
{
   return (JXF_Int) MPI_Alltoall(sendbuf, (int)sendcount, sendtype, recvbuf, (int)recvcount, recvtype, comm);
}

JXF_Int
jxf_MPI_Allgather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, MPI_Comm comm ) 
{
   return (JXF_Int) MPI_Allgather(sendbuf, (int)sendcount, sendtype, recvbuf, (int)recvcount, recvtype, comm);
}

JXF_Int
jxf_MPI_Allgatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, MPI_Comm comm ) 
{
   int *mpi_recvcounts, *mpi_displs, csize;
   JXF_Int  i;
   JXF_Int  ierr;

   MPI_Comm_size(comm, &csize);
   mpi_recvcounts = jxf_TAlloc(int, csize);
   mpi_displs = jxf_TAlloc(int, csize);
   for (i = 0; i < csize; i++)
   {
      mpi_recvcounts[i] = (int) recvcounts[i];
      mpi_displs[i] = (int) displs[i];
   }
   ierr = (JXF_Int) MPI_Allgatherv(sendbuf, (int)sendcount, sendtype, recvbuf, mpi_recvcounts, mpi_displs, recvtype, comm);
   jxf_TFree(mpi_recvcounts);
   jxf_TFree(mpi_displs);

   return ierr;
}

JXF_Int
jxf_MPI_Gather( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return (JXF_Int) MPI_Gather(sendbuf, (int) sendcount, sendtype, recvbuf, (int) recvcount, recvtype, (int)root, comm);
}

JXF_Int
jxf_MPI_Gatherv( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  JXF_Int *recvcounts, JXF_Int *displs, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   int *mpi_recvcounts = NULL;
   int *mpi_displs = NULL;
   int csize, croot;
   JXF_Int  i;
   JXF_Int  ierr;

   MPI_Comm_size(comm, &csize);
   MPI_Comm_rank(comm, &croot);
   if (croot == (int) root)
   {
      mpi_recvcounts = jxf_TAlloc(int, csize);
      mpi_displs = jxf_TAlloc(int, csize);
      for (i = 0; i < csize; i++)
      {
         mpi_recvcounts[i] = (int) recvcounts[i];
         mpi_displs[i] = (int) displs[i];
      }
   }
   ierr = (JXF_Int) MPI_Gatherv(sendbuf, (int)sendcount, sendtype, recvbuf, mpi_recvcounts, mpi_displs, recvtype, (int) root, comm);
   jxf_TFree(mpi_recvcounts);
   jxf_TFree(mpi_displs);

   return ierr;
}

JXF_Int
jxf_MPI_Scatter( void *sendbuf, JXF_Int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   return (JXF_Int) MPI_Scatter(sendbuf, (int)sendcount, sendtype, recvbuf, (int)recvcount, recvtype, (int)root, comm);
}

JXF_Int
jxf_MPI_Scatterv( void *sendbuf, JXF_Int *sendcounts, JXF_Int *displs, MPI_Datatype sendtype,
                   void *recvbuf, JXF_Int recvcount, MPI_Datatype recvtype, JXF_Int root, MPI_Comm comm )
{
   int *mpi_sendcounts = NULL;
   int *mpi_displs = NULL;
   int csize, croot;
   JXF_Int  i;
   JXF_Int  ierr;

   MPI_Comm_size(comm, &csize);
   MPI_Comm_rank(comm, &croot);
   if (croot == (int) root)
   {
      mpi_sendcounts = jxf_TAlloc(int, csize);
      mpi_displs = jxf_TAlloc(int, csize);
      for (i = 0; i < csize; i++)
      {
         mpi_sendcounts[i] = (int) sendcounts[i];
         mpi_displs[i] = (int) displs[i];
      }
   }
   ierr = (JXF_Int) MPI_Scatterv(sendbuf, mpi_sendcounts, mpi_displs, sendtype, recvbuf, (int) recvcount, recvtype, (int) root, comm);
   jxf_TFree(mpi_sendcounts);
   jxf_TFree(mpi_displs);

   return ierr;
}

JXF_Int
jxf_MPI_Bcast( void *buffer, JXF_Int count, MPI_Datatype datatype, JXF_Int root, MPI_Comm comm ) 
{ 
   return (JXF_Int) MPI_Bcast(buffer, (int)count, datatype, (int)root, comm);
}

JXF_Int
jxf_MPI_Send( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm ) 
{ 
   return (JXF_Int) MPI_Send(buf, (int)count, datatype, (int)dest, (int)tag, comm);
}

JXF_Int
jxf_MPI_Recv( void *buf, JXF_Int count, MPI_Datatype  datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status )
{ 
   return (JXF_Int) MPI_Recv(buf, (int)count, datatype, (int)source, (int)tag, comm, status);
}

JXF_Int
jxf_MPI_Isend( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{ 
   return (JXF_Int) MPI_Isend(buf, (int)count, datatype, (int)dest, (int)tag, comm, request);
}

JXF_Int
jxf_MPI_Irecv( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{ 
   return (JXF_Int) MPI_Irecv(buf, (int)count, datatype, (int)source, (int)tag, comm, request);
}

JXF_Int
jxf_MPI_Send_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return (JXF_Int) MPI_Send_init(buf, (int)count, datatype, (int)dest, (int)tag, comm, request);
}

JXF_Int
jxf_MPI_Recv_init( void *buf, JXF_Int count, MPI_Datatype datatype, JXF_Int dest,
                     JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return (JXF_Int) MPI_Recv_init(buf, (int)count, datatype, (int)dest, (int)tag, comm, request);
}

JXF_Int
jxf_MPI_Irsend( void *buf, JXF_Int count, MPI_Datatype datatype,
                  JXF_Int dest, JXF_Int tag, MPI_Comm comm, MPI_Request *request )
{
   return (JXF_Int) MPI_Irsend(buf, (int)count, datatype, (int)dest, (int)tag, comm, request);
}

JXF_Int
jxf_MPI_Startall( JXF_Int count, MPI_Request *array_of_requests )
{
   return (JXF_Int) MPI_Startall((int)count, array_of_requests);
}

JXF_Int
jxf_MPI_Probe( JXF_Int source, JXF_Int tag, MPI_Comm comm, MPI_Status *status )
{
   return (JXF_Int) MPI_Probe((int)source, (int)tag, comm, status);
}

JXF_Int
jxf_MPI_Iprobe( JXF_Int source, JXF_Int tag, MPI_Comm comm, JXF_Int *flag, MPI_Status *status )
{
   int mpi_flag;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Iprobe((int)source, (int)tag, comm, &mpi_flag, status);
  *flag = (JXF_Int) mpi_flag;
   return ierr;
}

JXF_Int
jxf_MPI_Test( MPI_Request *request, JXF_Int *flag, MPI_Status *status )
{
   int mpi_flag;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Test(request, &mpi_flag, status);
  *flag = (JXF_Int) mpi_flag;
   return ierr;
}

JXF_Int
jxf_MPI_Testall( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *flag, MPI_Status *array_of_statuses )
{
   int mpi_flag;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Testall((int)count, array_of_requests, &mpi_flag, array_of_statuses);
  *flag = (JXF_Int) mpi_flag;
   return ierr;
}

JXF_Int
jxf_MPI_Wait( MPI_Request *request, MPI_Status *status )
{
   return (JXF_Int) MPI_Wait(request, status);
}

JXF_Int
jxf_MPI_Waitall( JXF_Int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses )
{
   return (JXF_Int) MPI_Waitall((int)count, array_of_requests, array_of_statuses);
}

JXF_Int
jxf_MPI_Waitany( JXF_Int count, MPI_Request *array_of_requests, JXF_Int *index, MPI_Status *status )
{
   int mpi_index;
   JXF_Int ierr;
   ierr = (JXF_Int) MPI_Waitany((int)count, array_of_requests, &mpi_index, status);
  *index = (JXF_Int) mpi_index;
   return ierr;
}

JXF_Int
jxf_MPI_Allreduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
   return (JXF_Int) MPI_Allreduce(sendbuf, recvbuf, (int)count, datatype, op, comm);
}

JXF_Int
jxf_MPI_Reduce( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, JXF_Int root, MPI_Comm comm )
{
   return (JXF_Int) MPI_Reduce(sendbuf, recvbuf, (int)count, datatype, op, (int)root, comm);
}

JXF_Int
jxf_MPI_Scan( void *sendbuf, void *recvbuf, JXF_Int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{ 
   return (JXF_Int) MPI_Scan(sendbuf, recvbuf, (int)count, datatype, op, comm);
}

JXF_Int
jxf_MPI_Request_free( MPI_Request *request )
{
   return (JXF_Int) MPI_Request_free(request);
}

JXF_Int
jxf_MPI_Type_contiguous( JXF_Int count, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return (JXF_Int) MPI_Type_contiguous((int)count, oldtype, newtype);
}

JXF_Int
jxf_MPI_Type_vector( JXF_Int count, JXF_Int blocklength, JXF_Int stride, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return (JXF_Int) MPI_Type_vector((int)count, (int)blocklength, (int)stride, oldtype, newtype);
}

JXF_Int
jxf_MPI_Type_hvector( JXF_Int count, JXF_Int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype )
{
   return (JXF_Int) MPI_Type_create_hvector((int)count, (int)blocklength, stride, oldtype, newtype);
}

JXF_Int
jxf_MPI_Type_struct( JXF_Int count, JXF_Int *array_of_blocklengths, MPI_Aint *array_of_displacements,
                       MPI_Datatype *array_of_types, MPI_Datatype *newtype )
{
   int *mpi_array_of_blocklengths;
   JXF_Int  i;
   JXF_Int  ierr;

   mpi_array_of_blocklengths = jxf_TAlloc(int, count);
   for (i = 0; i < count; i++)
   {
      mpi_array_of_blocklengths[i] = (int) array_of_blocklengths[i];
   }
   ierr = (JXF_Int) MPI_Type_create_struct((int)count, mpi_array_of_blocklengths, array_of_displacements, array_of_types, newtype);
   jxf_TFree(mpi_array_of_blocklengths);

   return ierr;
}

JXF_Int
jxf_MPI_Type_commit( MPI_Datatype *datatype )
{
   return (JXF_Int) MPI_Type_commit(datatype);
}

JXF_Int
jxf_MPI_Type_free( MPI_Datatype *datatype )
{
   return (JXF_Int) MPI_Type_free(datatype);
}

JXF_Int
jxf_MPI_Op_free( MPI_Op *op )
{
   return (JXF_Int) MPI_Op_free(op);
}

JXF_Int
jxf_MPI_Op_create( MPI_User_function *function, int commute, MPI_Op *op )
{
   return (JXF_Int) MPI_Op_create(function, commute, op);
}

#endif
