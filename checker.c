#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

#include "check.h"
#include "util.h"
#include "problem.h"
#include "solution.h"

static problem_t  p;
static solution_t s;

#define DIFFTEMPS(a,b) (((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)

int main(int argc, char* argv[])
{
    double t1 = 0.,t2 =0.,t3 =0.;
    int score ,sum_score;
    int rang =0, size;
    struct timeval tv_begin, tv_end,inter;

    if (argc != 3) {
        fprintf(stderr, "usage: %s problem solution\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    CHECK(problem_read(argv[1], &p) == 0);
    CHECK(solution_read(argv[2], &s, &p) == 0);

    if( MPI_Init(NULL, NULL))
    {
    	fprintf(stderr, "Erreur MPI_Init\n" );
    	exit(1);
    }
    MPI_Comm_rank( MPI_COMM_WORLD, &rang );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    //printf("%d\n",size );
    MPI_Barrier(MPI_COMM_WORLD);

    if(rang == 0)
        t1 = MPI_Wtime();
    CHECK(solution_check(&s, &p) == 0);
    if(rang == 0)
       t2 = MPI_Wtime();
    score = solution_score(&s, &p);
    if(rang == 0)
       t3 = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);


    MPI_Reduce(&score, &sum_score, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //on a size score ....
    if(rang == 0)
    {
        printf("Temps de sol_check :  %lgs\n",t2-t1 );
        printf("Temps de sol_score:  %lgs\n",t3-t2 );
        printf("Temps de check :  %lgs\n",t3-t1 );
        fprintf(stderr, "Score %d\n", sum_score);

        // Write the score file
        util_write_score(argv[2], sum_score);
    }


    MPI_Finalize();
    return(0);
}
