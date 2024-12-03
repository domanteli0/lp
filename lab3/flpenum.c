#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "wrappers.h"
#include "shared.h"
#include "matrix.h"

//===== Globalus kintamieji ===================================================

const int numDP = 5000; // Vietoviu skaicius (demand points, max 10000)
const int numPF = 5;	  // Esanciu objektu skaicius (preexisting facilities)
const int numCL = 50;	  // Kandidatu naujiems objektams skaicius (candidate locations)
const int numX = 3;	  // Nauju objektu skaicius

double **demandPoints;	 // Geografiniai duomenys
double *distanceMatrix; // Masyvas atstumu matricai saugoti

// int X[numX];        // Naujas sprendinys
// int bestX[numX];    // Geriausias rastas sprendinys
double u, bestU;			// Naujo sprendinio ir geriausio sprendinio naudingumas (utility)

//===== Funkciju prototipai ===================================================

double getTime();
void loadDemandPoints();
double HaversineDistance4(double lat1, double lon1, double lat2, double lon2);
double HaversineDistance2(int i, int j);
double evaluateSolution(int *X);
int increaseX(int *X, int index, int maxindex);

//=============================================================================

void printf_(char *_, ...) { }
// #define printf printf_

#define NOT(x) !(x)
#define AND &&

#define SIGNAL_DONE 10
#define SIGNAL_WORKER_DONE 11
#define DATA_X 1
#define DATA_U 2
#define DATA_BEST_U 3

void write_times(double t_start, double t_matrix, double t_finish, double t_iallgatherv);
void printX(int *X);

int world_size;
int world_rank;

int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);

   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // Get the rank of the process
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   loadDemandPoints();			// Duomenu nuskaitymas is failo
   double t_start = getTime(); // Algoritmo vykdymo pradzios laikas

   void *buffer;
   // int buffer_size = (2*world_size) * (sizeof(int) * numX * MPI_BSEND_OVERHEAD);
   int buffer_size = 1024 * (2*world_size) * (sizeof(int) * numX * MPI_BSEND_OVERHEAD);
   buffer = malloc(buffer_size);
   MPI_Buffer_attach(buffer, buffer_size);

   //----- Atstumu matricos skaiciavimas -------------------------------------
   int *lens = lengths(numDP, world_size);
   int *counts = calloc(sizeof(int), world_size);
   for (int ix = 0; ix < world_size; ++ix) {
      counts[ix] = (lens[ix + 1] - lens[ix]) * numDP;
   }

   int *disps = calloc(sizeof(int), world_size);
   for (int ix = 0; ix < world_size; ++ix) {
      disps[ix] = lens[ix] * numDP;
   }

   distanceMatrix = calloc(sizeof(double), numDP * numDP);
   for (int i = lens[world_rank]; i < lens[world_rank + 1]; i++) {
      for (int j = 0; j <= i; j++) {
         distanceMatrix[numDP * i + j] =
            HaversineDistance4(demandPoints[i][0], demandPoints[i][1], demandPoints[j][0], demandPoints[j][1]);
      }
   }

   MPI_Request req;
   MPI_Iallgatherv(
      distanceMatrix + disps[world_rank],
      counts[world_rank],
      MPI_DOUBLE,
      distanceMatrix,
      counts,
      disps,
      MPI_DOUBLE,
      MPI_COMM_WORLD
      ,&req
   );

   double t_iallgatherv = getTime();

   double t_matrix = 0;
   // printf("Matricos skaiciavimo trukme: %lf\n", t_matrix - t_start);

   //----- Pradines naujo ir geriausio sprendiniu reiksmes -------------------
   int *X = calloc(sizeof(int), numX);        // Naujas sprendinys
   int *bestX = calloc(sizeof(int), numX);    // Geriausias rastas sprendinys
   int *copy_X = calloc(sizeof(int), numX);

   u = 0; // Naujo sprendinio naudingumas (utility)
   bestU = u;               // Geriausio sprendinio sprendinio naudingumas
   double copy_u = u;

   //----- Visų galimų sprendinių perrinkimas --------------------------------
   bool increased = true;
   unsigned char dones[world_size];
   int dummy_load = 0;
   for (int ix = 0; ix < world_size; ++ix) { dones[ix] = dones[ix] & 0; }

   int sends = 0;
   int receives = 0;
   bool first_run = true;
   while (true) {
      if (world_rank == 0 && increased) {
         for(int ix = 1; ix < world_size; ++ix) {
            increased = increaseX(X, numX - 1, numCL);
            MPI_Bsend(X, numX, MPI_INT, ix, DATA_X, MPI_COMM_WORLD);
            sends += 1;

            if (!increased) {
               for (int ix = 1; ix < world_size; ++ix) {
                  MPI_Bsend(&dummy_load, 1, MPI_INT, ix, SIGNAL_DONE, MPI_COMM_WORLD);
               }
               break;
            }
         }
      }

      if (world_rank == 0 && increased) {
         increased = increaseX(X, numX - 1, numCL);
         if (first_run) { MPI_Wait(&req, MPI_STATUS_IGNORE); first_run = false; t_matrix = getTime(); }
         
         u = evaluateSolution(X);
         if (u > bestU) {
            bestU = u;
            memcpy(bestX, X, sizeof(int) * numX);
         }

         if (!increased) {
            for (int ix = 1; ix < world_size; ++ix) {
               MPI_Bsend(&dummy_load, 1, MPI_INT, ix, SIGNAL_DONE, MPI_COMM_WORLD);
            }
         }
      }

      if (world_rank != 0) {
         int master_done = WRP_Check_for(0, SIGNAL_DONE, MPI_COMM_WORLD);
         int master_sent_X = WRP_Check_for(0, DATA_X, MPI_COMM_WORLD);

         if (master_sent_X) {
            MPI_Recv(X, numX, MPI_INT, 0, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (first_run) { MPI_Wait(&req, MPI_STATUS_IGNORE); first_run = false; t_matrix = getTime(); }

            u = evaluateSolution(X);
            if (bestU < u) { bestU = u; bestX = memcpy(bestX, X, sizeof(int) * numX); }

            MPI_Bsend(&u, 1, MPI_DOUBLE, 0, DATA_U, MPI_COMM_WORLD);
            MPI_Bsend(X, numX, MPI_INT,  0, DATA_X, MPI_COMM_WORLD);
            sends += 1;
         }

         if (master_done && !master_sent_X) { break; }
      }

      if (world_rank == 0) {
         if (!increased && receives == sends) { break; }

         WRP_Check worker_sent_X = WRP_Check_for_(MPI_ANY_SOURCE, DATA_X, MPI_COMM_WORLD);
         while(worker_sent_X.flag) {
            MPI_Recv(&copy_u, 1, MPI_DOUBLE, worker_sent_X.status.MPI_SOURCE, DATA_U, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(copy_X, numX, MPI_INT, worker_sent_X.status.MPI_SOURCE, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            receives += 1;

            if (copy_u > bestU) {
               bestU = copy_u;
               memcpy(bestX, copy_X, sizeof(int) * numX);
            }

            worker_sent_X = WRP_Check_for_(MPI_ANY_SOURCE, DATA_X, MPI_COMM_WORLD);
         }
      }
   }

   double t_finish = getTime(); // Skaiciavimu pabaigos laikas
   //----- Rezultatu spausdinimas --------------------------------------------
   // if (world_rank == 0) {
   //    printf("Sprendinio paieskos trukme: %lf\n", t_finish - t_matrix);
   //    printf("Algoritmo vykdymo trukme: %lf\n", t_finish - t_start);

   //    display_results("stdout" , bestX, bestU);
   //    // display_results("new.dat", bestX, bestU);
   // }
   write_times(t_start, t_matrix, t_finish, t_iallgatherv);
	
   // printf("\n");

   MPI_Finalize();
}

//===== Funkciju implementacijos (siu funkciju LYGIAGRETINTI NEREIKIA) ========

void loadDemandPoints()
{
   FILE *f;
   f = fopen("demandPoints.dat", "r");
   demandPoints = calloc(sizeof(double *), numDP);
   for (int i = 0; i < numDP; i++) {
      demandPoints[i] = calloc(sizeof(double), 3);
      fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
   }
   fclose(f);
}

//=============================================================================

double HaversineDistance4(double lat1, double lon1, double lat2, double lon2) {
   double dlat = fabs(lat1 - lat2);
   double dlon = fabs(lon1 - lon2);
   double aa = pow((sin((double)dlat / (double)2 * 0.01745)), 2) + cos(lat1 * 0.01745) *
                                                      cos(lat2 * 0.01745) * pow((sin((double)dlon / (double)2 * 0.01745)), 2);
   double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
   double d = 6371 * c;
   return d;
}

double HaversineDistance2(int i, int j) {
   if (i >= j)
      return distanceMatrix[numDP * i + j];
   else
      return distanceMatrix[numDP * j + i];
}

//=============================================================================

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
   return rez;
}

//=============================================================================

double evaluateSolution(int *X) {
   double U = 0;
   double totalU = 0;
   int bestPF;
   int bestX;
   double d;
   for (int i = 0; i < numDP; i++) {
      totalU += demandPoints[i][2];

      bestPF = 1e5;
      for (int j = 0; j < numPF; j++) {
         d = HaversineDistance2(i, j);
         if (d < bestPF)
            bestPF = d;
      }

      bestX = 1e5;
      for (int j = 0; j < numX; j++) {
         d = HaversineDistance2(i, X[j]);
         if (d < bestX)
            bestX = d;
      }
      if (bestX < bestPF)
         U += demandPoints[i][2];
      else if (bestX == bestPF)
         U += 0.3 * demandPoints[i][2];
   }
   return U / totalU * 100;
}

//=============================================================================

int increaseX(int *X, int index, int maxindex) {
   if (X[index] + 1 < maxindex - (numX - index - 1)) {
      X[index]++;
   } else {
      if ((index == 0) && (X[index] + 1 == maxindex - (numX - index - 1))) {
         return 0;
      } else {
         if (increaseX(X, index - 1, maxindex))
            X[index] = X[index - 1] + 1;
         else
            return 0;
      }
   }
   return 1;
}

//=============================================================================

void write_times(double t_start, double t_matrix, double t_finish, double t_iallgatherv) {

   char *filename_buffer = (char *) calloc(sizeof(char), 1000);
   sprintf(filename_buffer, "parallel_%i_extra.tsv", world_size);
   FILE *fp = fopen(filename_buffer, "a+");

   // FILE *fp = stdout;

   fprintf(fp, "%i\t%i\t%f\t%f\t%f\t%lf\t%lf\n",
         world_size,
         world_rank,
         t_matrix - t_start,
         t_finish - t_matrix,
         t_finish - t_start,
         t_matrix - t_iallgatherv,
         t_iallgatherv - t_start);
   
   fclose(fp);
}
