#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <mpi.h>

#include "check.h"
#include "util.h"
#include "problem.h"
#include "solution.h"

#ifdef DEBUG
  #ifndef DEBUG_SCORE
    #define DEBUG_SCORE
  #endif
#endif
//#define DEBUG_SCORE

int solution_check(solution_t* const s, problem_t* const p )
{
    int rang =0, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rang );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    //printf("%d %d\n",size , rang );
  /* OK: errors = 0. */
  int errors = 0;
  int sum_errors = 0;

  //  const int nb_inter = p->NI;
  const int nb_streets = p->S;
  const int nb_inter_sol = s->A;
  int feu ;


  #pragma omp parallel for reduction(+:errors) private(feu)
  for(int i=rang*(nb_inter_sol/size); i<(rang +1) *(nb_inter_sol/size); i++)
  {

     // #pragma omp ordered
     // {
         // vérifie la solution pour l'intersection num i : s->schedule[i]
         if(s->schedule[i].nb < 1)
         {
           fprintf(stderr, "intersection has no light (%d)\n", i);
         }
         for(feu =0; feu<s->schedule[i].nb; feu++)
         {

             int rue;
             char* name;
             // s->schedule[i].t[feu] .rue et .duree sont valides
              #pragma omp critical
              {
                 rue = s->schedule[i].t[feu].rue;
                 name =(char *) street_table_find_name(p->table, rue);
             }

             //printf("%p %p %d \n",(void *)&rue ,(void*)&name  , omp_get_num_threads());
             if(rue >= nb_streets)
             {
                 fprintf(stderr, "invalid street number (%d -> \"%s\")\n", rue, name);
                 errors = rang;
             }
             errors = rang;
             int rid;
             // vérifie que cette rue (rue) arrive bien à cette intersection (i)
             for(rid=0; rid<nb_streets; rid++)
             {
                 //printf("%d\n" , omp_get_thread_num());
                 if(p->r[rid].street_id == rue)
                     break;
             }
             // p->r[rid] contient la rue, vérifie que la rue arrive bien à cette intersection
             if(p->r[rid].end != i)
             {
                 fprintf(stderr, "invalid street number (%d -> \"%s\"): not arriving to the intersection %d\n", rue, name, i);
                 errors = rang;
             }

             // durée > 0
             if(s->schedule[i].t[feu].duree <= 0)
             {
                 fprintf(stderr, "invalid schedule length (intersection %d light %d -> %d)\n", i, feu, s->schedule[i].t[feu].duree);
             }
          }
  }
  MPI_Reduce(&errors, &sum_errors, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  // if(rang == 0)
  //   printf("%d\n",sum_errors );
  //printf("%d\n",errors );

  /* OK */
  return errors;
}


typedef struct car_state {
  int street;     // Current street id,
  int distance;   // Remaining distance to end of street (0: end of street),
  int position;   // Position in the queue of the street (1: top position),
  int arrived;    // Arrived or not (Boolean),
  int nb_streets; // Number of streets already travelled
} car_state_t;

typedef struct street_state {
  int green;      // 1: green, 0: red
  int nb_cars;    // Number of cars in the street
  int max;        // Max number of cars in the street
  int out;        // A car just left the street (Boolean)
} street_state_t;


static car_state_t car_state[NB_CARS_MAX];
static street_state_t street_state[NB_STREETS_MAX];

void simulation_init(const problem_t* const p)
{
  memset(car_state, 0, NB_CARS_MAX * sizeof(car_state_t));
  memset(street_state, 0, NB_STREETS_MAX * sizeof(street_state_t));
  #pragma omp parallel for
  for (int i = 0; i < p->V; i++) {
    car_state[i].street = p->c[i].streets[0];
    car_state[i].distance = 0;
    // Queue the car
    street_state[car_state[i].street].nb_cars++;
    if (street_state[car_state[i].street].nb_cars > street_state[car_state[i].street].max)
      street_state[car_state[i].street].max = street_state[car_state[i].street].nb_cars;
    car_state[i].position = street_state[car_state[i].street].nb_cars;
    car_state[i].arrived = 0;
    car_state[i].nb_streets = 0;
  }
}

void simulation_update_intersection_lights(const solution_t* const s, int i, int T)
{
  int cycle = 0;
  int tick = 0;
  int no_green_light = 1;

  // Find the light cycle total time
  //#pragma omp parallel for reduction(+:cycle)
  for (int l = 0; l < s->schedule[i].nb; l++) {
    cycle += s->schedule[i].t[l].duree;
  }

  // Find at which time in the cycle we are
  tick = T % cycle;

  //printf("Inter %d, cycle %d, tick %d, T %d\n", i, cycle, tick, T);

  // Set the light state
  //#pragma omp parallel for reduction(-:tick)
  for (int l = 0; l < s->schedule[i].nb; l++) {

        // Remove duration, if we get below zero, this light is green and others are red
        tick -= s->schedule[i].t[l].duree;
        //printf("light %d, tick %d, duree %d\n", l, tick,  s->schedule[i].t[l].duree);
        if (tick < 0) {
          street_state[s->schedule[i].t[l].rue].green = 1;
          no_green_light = 0;
          //#pragma omp parallel for
          for (int next = l + 1; next < s->schedule[i].nb; next++) {
            street_state[s->schedule[i].t[next].rue].green = 0;
          }
          break;
        }
        street_state[s->schedule[i].t[l].rue].green = 0;
  }

  if (no_green_light) {
    printf("PROBLEM: NO GREEN LIGHT AT INTERSECTION %d (cycle %d)\n", i, cycle);
  }
}

int simulation_update_car(const problem_t* const p, int c, int T)
{
  // If already arrived, nothing to do
  if (car_state[c].arrived == 1)
    return 0;

  // If at the end of street, light green, queue 1 then move to next street
  if ((car_state[c].distance == 0) &&
      (street_state[car_state[c].street].green == 1) &&
      (car_state[c].position == 1)) {
    // Update number of street finished
    car_state[c].nb_streets++;
    // Signal a car left the street
    street_state[car_state[c].street].out = 1;
    // Set the new street where the car is
    car_state[c].street = p->c[c].streets[car_state[c].nb_streets];
    car_state[c].distance = p->r[car_state[c].street].len - 1;
    // Enqueue the car in the new street
    street_state[car_state[c].street].nb_cars++;
    if (street_state[car_state[c].street].nb_cars > street_state[car_state[c].street].max)
      street_state[car_state[c].street].max = street_state[car_state[c].street].nb_cars;
    car_state[c].position = street_state[car_state[c].street].nb_cars;
  } else if (car_state[c].distance > 0) {
    // If not at the end of street, advance
    car_state[c].distance--;
  }

  // If now at the last street AND at the end of the street: drive complete!
  if ((car_state[c].street == p->c[c].streets[p->c[c].P - 1]) &&
      (car_state[c].distance == 0)) {
    car_state[c].arrived = 1;
    // Remove the car immediately from the street
    street_state[car_state[c].street].nb_cars--;
    // If another car is in that street and was there before that car, dequeue it
    //#pragma omp parallel for
    for (int i = 0; i < p->V; i++) {
      if ((car_state[c].street == car_state[i].street) &&
          (car_state[c].position < car_state[i].position)) {
        car_state[i].position--;
      }
    }
    return p->F + (p->D - (T + 1));
  }

  return 0;
}

void simulation_print_state(const problem_t* const p, int T) {
  printf("Timestep: %d\n", T);
  #pragma omp parallel for
  for (int c = 0; c < p->V; c++) {
    printf("Car %d -> street %d, distance: %d, position: %d, "
           "arrived: %d, street#: %d\n",
        c,
        car_state[c].street,
        car_state[c].distance,
        car_state[c].position,
        car_state[c].arrived,
        car_state[c].nb_streets);
  }
  #pragma omp parallel for
  for (int s = 0; s < p->S; s++) {
    printf("Street %d -> green: %d, nb_cars: %d, out: %d\n",
        s,
        street_state[s].green,
        street_state[s].nb_cars,
        street_state[s].out);
  }
}

void simulation_dequeue(const problem_t* const p)
{
    int c = 0;
  //#pragma omp parallel for private(c)
  for (int street = 0; street < p->S; street++) {
    // If there is a street to dequeue
    if (street_state[street].out == 1) {
      // If a car is in that street, dequeue it
      for ( c = 0; c < p->V; c++) {
        if (car_state[c].street == street) {
          car_state[c].position--;
        }
      }
      street_state[street].nb_cars--;
      street_state[street].out = 0;
    }
  }
}

//#define DEBUG_SCORE

int simulation_run(const solution_t* const s, const problem_t* const p)
{
  int score = 0;

  #ifdef DEBUG_SCORE
  problem_write(stdout, p);
  solution_write(stdout, s, p);
  #endif

  // Init state
  simulation_init(p);

  // For each time step
  //omp_set_nested(true);
  int i= 0;
  #pragma omp parallel for private(i) schedule(dynamic)
  for (int T = 0; T < p->D; T++) {
    #ifdef DEBUG_SCORE
    printf("Score: %d\n", score);
    printf("- 1 Init:\n");
    simulation_print_state(p, T);
    #endif

    // Update light state for each intersection
    //#pragma omp parallel for
    for (i = 0; i < s->A; i++) {
      simulation_update_intersection_lights(s, i, T);
    }

    #ifdef DEBUG_SCORE
    printf("- 2 lights:\n");
    simulation_print_state(p, T);
    #endif

    // Update car state
    //#pragma omp parallel for reduction(+:score)
    for (int c = 0; c < p->V; c++) {

      score += simulation_update_car(p, c, T);
    }

    #ifdef DEBUG_SCORE
    printf("- 3 cars (score now = %d):\n", score);
    simulation_print_state(p, T);
    #endif

    simulation_dequeue(p);

    #ifdef DEBUG_SCORE
    printf("- 4 queues:\n");
    simulation_print_state(p, T);
    #endif
  }
  return score;
}

/*
static int score_descending_order(void const* r1, void const* r2) {
  int d1 = ((int*)r1)[1];
  int d2 = ((int*)r2)[1];
  return (d1 == d2) ? 0 : ((d1 > d2) ? -1 : 1);
}

static int score_ascending_order(void const* r1, void const* r2) {
  int d1 = ((int*)r1)[1];
  int d2 = ((int*)r2)[1];
  return (d1 == d2) ? 0 : ((d1 < d2) ? -1 : 1);
}
*/

int tab_street[NB_STREETS_MAX][2];
int tab_car[NB_CARS_MAX][2];

int solution_score(solution_t* s, const problem_t* const p)
{
  int score = 0;

  score = simulation_run(s, p);



#if 0
  printf("Score = %d\n", score);


  //#pragma omp parallel for
  #pragma omp parallel for
  for (int n = 0; n <= 10; n++) {

    memset(tab_street, 0, NB_STREETS_MAX * 2 * sizeof(int));
    memset(tab_car, 0, NB_CARS_MAX * 2 * sizeof(int));

    for (int street = 0; street < p->S; street++) {
      tab_street[street][0] = street;
      tab_street[street][1] = street_state[street].max;

      // Remove lights in streets without any car
      if (street_state[street].max == 0 && n >= 5) {
        for (int i = 0; i < s->A; i++) {
          if (s->schedule[i].nb > 1) {
            int rue;
            int nb_rues = s->schedule[i].nb;
            for (rue = 0; rue < nb_rues; rue++) {
              if (street == s->schedule[i].t[rue].rue) {
                s->schedule[i].nb--;
                break;
              }
          }
            rue++;
            for (; rue < nb_rues; rue++) {
              s->schedule[i].t[rue - 1].rue = s->schedule[i].t[rue].rue;
              s->schedule[i].t[rue - 1].duree = s->schedule[i].t[rue].duree;
            }
          }
        }
      }
    }

    // Add + 1 to green light time in jammed streets
    qsort(tab_street, p->S, sizeof(*tab_street), score_descending_order);
    //#pragma omp for
    for (int jam = 0; jam < MIN(p->S, 5); jam++) {
      //printf("[Street %d, max %d], ", tab_street[jam][0], street_state[tab_street[jam][0]].max);
      if (tab_street[jam][1] > p->r[tab_street[jam][0]].len) {
        //#pragma omp for
        for (int i = 0; i < s->A; i++) {
          //#pragma omp for
          for (int rue = 0; rue < s->schedule[i].nb; rue++) {
            if (tab_street[jam][0] == s->schedule[i].t[rue].rue) {
              (s->schedule[i].t[rue].duree)+=2;
            }
          }
        }
      }
    }
    #pragma omp singel
    score = simulation_run(s, p);
    printf("Score = %d\n", score);
  }
#endif

  return score;
}
