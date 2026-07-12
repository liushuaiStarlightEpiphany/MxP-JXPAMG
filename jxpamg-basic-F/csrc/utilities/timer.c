//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"

/*!
 * \fn JXF_Real jxf_time_getWallclockSeconds
 * \brief Timing function.
 * \date 2011/09/02 
 */ 
JXF_Real 
jxf_time_getWallclockSeconds(void)
{
#ifdef TIMER_USE_MPI
   return(jxf_MPI_Wtime());
#else
#ifdef WIN32
   clock_t cl=clock();
   return(((JXF_Real) cl)/((JXF_Real) CLOCKS_PER_SEC));
#else
   struct tms usage;
   long wallclock = times(&usage);
   return(((JXF_Real) wallclock)/((JXF_Real) sysconf(_SC_CLK_TCK)));
#endif
#endif
}

/*!
 * \fn void jxf_GetWallTime
 * \brief Timing for jxf_MPI program
 * \param comm communicator
 * \param *fctname pointer to the name of the function that is timed
 * \param starttime starting time
 * \param endtime ending time
 * \param allproc whether to output the timing of all the processors
 * \param dp DecimalPlace of the output time
 * \author peghoty
 * \date 2010/11/22 
 */ 
void
jxf_GetWallTime( MPI_Comm  comm, 
                char     *fctname, 
                JXF_Real    starttime, 
                JXF_Real    endtime, 
                JXF_Int       allproc,
                JXF_Int       dp )
{
   JXF_Int myid, nprocs;

   JXF_Real totaltime    = 0.0;
   JXF_Real time_max     = 0.0;
   JXF_Real time_min     = 0.0;
   JXF_Real time_average = 0.0;
   
   jxf_MPI_Comm_rank(comm, &myid);      
   jxf_MPI_Comm_size(comm, &nprocs);
   
   totaltime = endtime - starttime;

   if (allproc)
   {   
      switch (dp)
      {
         case 1:
         jxf_printf(" >> %s: Proc %d took %.1f seconds!\n", fctname, myid, totaltime); 
         break;
         case 2:
         jxf_printf(" >> %s: Proc %d took %.2f seconds!\n", fctname, myid, totaltime); 
         break;
         case 3:
         jxf_printf(" >> %s: Proc %d took %.3f seconds!\n", fctname, myid, totaltime); 
         break;
         case 4:
         jxf_printf(" >> %s: Proc %d took %.4f seconds!\n", fctname, myid, totaltime); 
         break;
         case 5:
         jxf_printf(" >> %s: Proc %d took %.5f seconds!\n", fctname, myid, totaltime); 
         break; 
         case 6:
         jxf_printf(" >> %s: Proc %d took %.6f seconds!\n", fctname, myid, totaltime); 
         break;                                  
         default:
         jxf_printf(" >> %s: Proc %d took %.2f seconds!\n", fctname, myid, totaltime); 
         break;
      }    
   }
   jxf_MPI_Reduce(&totaltime, &time_min,     1, JXF_MPI_REAL, MPI_MIN, 0, comm);
   jxf_MPI_Reduce(&totaltime, &time_max,     1, JXF_MPI_REAL, MPI_MAX, 0, comm);
   jxf_MPI_Reduce(&totaltime, &time_average, 1, JXF_MPI_REAL, MPI_SUM, 0, comm);
   if (myid == 0)
   {
      switch (dp)
      {
         case 1:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.1f, %.1f, %.1f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 2:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.2f, %.2f, %.2f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 3:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.3f, %.3f, %.3f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 4:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.4f, %.4f, %.4f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 5:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.5f, %.5f, %.5f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 6:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.6f, %.6f, %.6f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;                           
         default:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.2f, %.2f, %.2f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
      }                
   } 
}

void
jxf_GetWallTime2( MPI_Comm  comm, 
                 char     *fctname, 
                 JXF_Real    totaltime,
                 JXF_Int       allproc,
                 JXF_Int       dp )
{
   JXF_Int myid, nprocs;

   JXF_Real time_max     = 0.0;
   JXF_Real time_min     = 0.0;
   JXF_Real time_average = 0.0;
   
   jxf_MPI_Comm_rank(comm, &myid);
   jxf_MPI_Comm_size(comm, &nprocs);
   
   if (allproc)
   {   
      switch (dp)
      {
         case 1:
         jxf_printf(" >> %s: Proc %d took %.1f seconds!\n", fctname, myid, totaltime); 
         break;
         case 2:
         jxf_printf(" >> %s: Proc %d took %.2f seconds!\n", fctname, myid, totaltime); 
         break;
         case 3:
         jxf_printf(" >> %s: Proc %d took %.3f seconds!\n", fctname, myid, totaltime); 
         break;
         case 4:
         jxf_printf(" >> %s: Proc %d took %.4f seconds!\n", fctname, myid, totaltime); 
         break;
         case 5:
         jxf_printf(" >> %s: Proc %d took %.5f seconds!\n", fctname, myid, totaltime); 
         break; 
         case 6:
         jxf_printf(" >> %s: Proc %d took %.6f seconds!\n", fctname, myid, totaltime); 
         break;                                  
         default:
         jxf_printf(" >> %s: Proc %d took %.2f seconds!\n", fctname, myid, totaltime); 
         break;
      }    
   }
   jxf_MPI_Reduce(&totaltime, &time_min,     1, JXF_MPI_REAL, MPI_MIN, 0, comm);
   jxf_MPI_Reduce(&totaltime, &time_max,     1, JXF_MPI_REAL, MPI_MAX, 0, comm);
   jxf_MPI_Reduce(&totaltime, &time_average, 1, JXF_MPI_REAL, MPI_SUM, 0, comm);
   if (myid == 0)
   {
      switch (dp)
      {
         case 1:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.1f, %.1f, %.1f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 2:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.2f, %.2f, %.2f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 3:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.3f, %.3f, %.3f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 4:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.4f, %.4f, %.4f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 5:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.5f, %.5f, %.5f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
         case 6:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.6f, %.6f, %.6f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;                           
         default:
         jxf_printf("\n >> %s: time(min,max,ave) = (%.2f, %.2f, %.2f) seconds\n\n", 
                fctname, time_min, time_max, time_average / nprocs);
         break;
      }                
   } 
}

/*!
 * \fn void jxf_GetWallTimeMax
 * \brief Timing for jxf_MPI program(get the maximal time of all processors)
 * \param comm communicator
 * \param starttime starting time
 * \param endtime ending time
 * \return the maximal of (endtime - starttime) among all the processors
 * \author peghoty
 * \date 2012/03/25 
 */ 
JXF_Real
jxf_GetWallTimeMax( MPI_Comm  comm, 
                   JXF_Real    starttime, 
                   JXF_Real    endtime )
{
   JXF_Int myid, nprocs;

   JXF_Real totaltime    = 0.0;
   JXF_Real time_max     = 0.0;
   
   jxf_MPI_Comm_rank(comm, &myid);      
   jxf_MPI_Comm_size(comm, &nprocs);
   
   totaltime = endtime - starttime;

   jxf_MPI_Reduce(&totaltime, &time_max, 1, JXF_MPI_REAL, MPI_MAX, 0, comm);

   return (time_max);
}
