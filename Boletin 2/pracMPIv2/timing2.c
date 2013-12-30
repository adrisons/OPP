#include <stdio.h>
#include <mpi.h>

/*NOTE: 
 The MPI_Wtime calls can be placed anywhere between the MPI_Init
 and MPI_Finalize calls.
*/

main(int argc, char **argv)
{
    int myrank,
        numprocs;
    double mytime,   /*variables used for gathering timing statistics*/
           maxtime,
           mintime,
           avgtime;


    MPI_Init(&argc,&argv);  
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Barrier(MPI_COMM_WORLD);  /*synchronize all processes*/
    mytime = MPI_Wtime();  /*get time just before work section */

    printf("Hello World from Node %d\n",myrank);

    mytime = MPI_Wtime() - mytime;  /*get time just after work section*/
 
    /*compute max, min, and average timing statistics*/
    MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    if (!myrank) {
      avgtime /= numprocs;
      printf("Min: %lf  Max: %lf  Avg:  %lf\n", mintime, maxtime,avgtime);
    }
   MPI_Finalize();

}
