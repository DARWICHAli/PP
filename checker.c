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
    int score;
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
    if(rang == 0)
        gettimeofday( &tv_begin, NULL);
    CHECK(solution_check(&s, &p) == 0);
    if(rang == 0)
        gettimeofday( &inter, NULL);
    score = solution_score(&s, &p);
    if(rang == 0)
        gettimeofday( &tv_end, NULL);


    //on a size score ....
    if(rang == 0)
    {
        printf("Temps de sol_check :  %lfs\n",DIFFTEMPS(tv_begin,inter));
        printf("Temps de sol_score:  %lfs\n",DIFFTEMPS(inter,tv_end));
        printf("Temps de check :  %lfs\n",DIFFTEMPS(tv_begin,tv_end));
        printf("\n");
        //fprintf(stderr, "Score %d\n", score);

        // Write the score file
        util_write_score(argv[2], score);
    }


    MPI_Finalize();
    return(0);
}
