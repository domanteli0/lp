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

//===== Globalus kintamieji ===================================================

const int numDP = 5000; // Vietoviu skaicius (demand points, max 10000)
const int numPF = 5;	  // Esanciu objektu skaicius (preexisting facilities)
const int numCL = 50;	  // Kandidatu naujiems objektams skaicius (candidate locations)
const int numX = 3;	  // Nauju objektu skaicius

double **demandPoints;	 // Geografiniai duomenys
double **distanceMatrix; // Masyvas atstumu matricai saugoti

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
#define printf printf_

#define NOT(x) !(x)
#define AND &&

#define SIGNAL_DONE 10
#define DATA_X 1
#define DATA_U 2
#define DATA_BEST_U 3

void display_results(char *filename, int *bestX, double bestU);
void write_times(double t_start, double t_matrix, double t_finish);
void printX(int *X);

int world_size;

struct UX {
   double u;
   int *X;
} ux;


int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);

   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   // MPI_Comm comm;
   // if (world_rank == 0) {
   //    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &comm);
   // } else {
   //    int color = 1;
   //    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &comm);
   // }

   loadDemandPoints();			// Duomenu nuskaitymas is failo
   double t_start = getTime(); // Algoritmo vykdymo pradzios laikas

   void *buffer;
   // int buffer_size = 1024 * 1024 * 1024 /*GiB*/ */; // malloc(sizeof(int) * world_size);
   // int buffer_size = 1024 * 1024 /*1 MiB*/; // malloc(sizeof(int) * world_size);
   int buffer_size = (2*world_size) * (sizeof(int) * numX * MPI_BSEND_OVERHEAD);
   buffer = malloc(buffer_size);
   MPI_Buffer_attach(buffer, buffer_size);

   //----- Atstumu matricos skaiciavimas -------------------------------------
   distanceMatrix = calloc(sizeof(double *), numDP);
   for (int i = 0; i < numDP; i++) {
      distanceMatrix[i] = calloc(sizeof(double), i + 1);
      for (int j = 0; j <= i; j++) {
         distanceMatrix[i][j] = HaversineDistance4(demandPoints[i][0], demandPoints[i][1], demandPoints[j][0], demandPoints[j][1]);
      }
   }

   printf("BARRIER A REACHED FROM %i\n", world_rank);
   MPI_Barrier(MPI_COMM_WORLD);

   double t_matrix = getTime();
   printf("Matricos skaiciavimo trukme: %lf\n", t_matrix - t_start);

   //----- Pradines naujo ir geriausio sprendiniu reiksmes -------------------
   int *X = calloc(sizeof(int), numX);        // Naujas sprendinys
   int *bestX = calloc(sizeof(int), numX);    // Geriausias rastas sprendinys
   int *copy_X = calloc(sizeof(int), numX);

   u = 0; // Naujo sprendinio naudingumas (utility)
   bestU = u;               // Geriausio sprendinio sprendinio naudingumas
   double copy_u = u;

   //----- Visų galimų sprendinių perrinkimas --------------------------------
   bool increased = true;
   unsigned long long times = 0;
   unsigned char dones[world_size];
   for (int ix = 0; ix < world_size; ++ix) { dones[ix] = dones[ix] & 0; }

   int counter = 0;
	// FILE *fp = fopen("c_main.tsv", "a+");
   while (true) {
      // counter += 1;
      // if (counter > 10 + 1) { break; }
      if (world_rank == 0 && increased) {
         printf("MAIN| about to send(X)\n");
         for(int ix = 1; ix < world_size; ++ix) {
            increased = increaseX(X, numX - 1, numCL);
            printf("MAIN| X: "); printX(X);

            // MPI_Request x_req;
            // IDEA: MPI_Scatter(..., 0, comm);
            MPI_Send(X, numX, MPI_INT, ix, DATA_X, MPI_COMM_WORLD);
            printf("MAIN| send(X)\n");

            if (!increased) {
               // MPI_Request done_req;
               // TODO: this could be acheived with MPI_Bcast
               for (int ix = 1; ix < world_size; ++ix) {
                  printf("MAIN| SENT DONE\n");
                  MPI_Send(X, 0, MPI_BYTE, ix, SIGNAL_DONE, MPI_COMM_WORLD);
               }
               break;
            }
         }
      }

      if (world_rank != 0) {
         bool master_done = WRP_Check_for(0, SIGNAL_DONE, MPI_COMM_WORLD);
         bool master_more_x = WRP_Check_for(0, DATA_X, MPI_COMM_WORLD);
         if (master_done && NOT(master_more_x)) {
            MPI_Send(X, 0, MPI_BYTE, 0, SIGNAL_DONE, MPI_COMM_WORLD);
            printf("WORKER %i| received DONE signal from MAIN\n", world_rank);
            break;
         }

         MPI_Recv(X, numX, MPI_INT, 0, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         printf("WORKER %i| x: ", world_rank); printX(X);

         u = evaluateSolution(X);

         printf("WORKER %i| u: %lf\n", world_rank, u);

         MPI_Bsend(&u, 1, MPI_DOUBLE, 0, DATA_U, MPI_COMM_WORLD);
         MPI_Bsend(X, numX, MPI_INT,  0, DATA_X, MPI_COMM_WORLD);
         printf("WORKER %i| send(X)\n", world_rank);
      }

      if (world_rank == 0) {
         WRP_Check worker_check = WRP_Check_for_(MPI_ANY_SOURCE, SIGNAL_DONE, MPI_COMM_WORLD);
         if (worker_check.flag) {
            dones[worker_check.status.MPI_SOURCE] = UCHAR_MAX;
            printf("MAIN: received DONE signal from worker %i\n", worker_check.status.MPI_SOURCE);
         }

         bool done = true;
         for (int ix = 1; ix < world_size; ++ix) { done = done && (dones[ix] == UCHAR_MAX); }
         if (done) { break; }

         WRP_Check worked_sent_X = WRP_Check_for_(MPI_ANY_SOURCE, DATA_X, MPI_COMM_WORLD);
         if (worked_sent_X.flag) {
            MPI_Recv(&copy_u, 1, MPI_DOUBLE, worked_sent_X.status.MPI_SOURCE, DATA_U, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(copy_X, numX, MPI_INT, worked_sent_X.status.MPI_SOURCE, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("MAIN: recv(X)\n");
            // fprint_X(fp, copy_u, copy_X);

            if (copy_u > bestU) {
               bestU = copy_u;
               memcpy(bestX, copy_X, sizeof(int) * numX);
            }
         }
      }
   }

   //----- Rezultatu spausdinimas --------------------------------------------
   if (world_rank == 0) {
      double t_finish = getTime(); // Skaiciavimu pabaigos laikas
      printf("Sprendinio paieskos trukme: %lf\n", t_finish - t_matrix);
      printf("Algoritmo vykdymo trukme: %lf\n", t_finish - t_start);
	   printf("TIMES: %lli\n", times);

      display_results("stdout" , bestX, bestU);
      display_results("new.dat", bestX, bestU);
      write_times(t_start, t_matrix, t_finish);
   }

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
      return distanceMatrix[i][j];
   else
      return distanceMatrix[j][i];
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

void display_results(char *filename, int *bestX, double bestU) {
   const char *cmp = "stdout";
   FILE *fp;

   if (strcmp(filename, cmp) == 0)
   {
      fp = stdout;
   }
   else
   {
      fp = fopen(filename, "w+");
   }

   for (int i = 0; i < numX; i++)
   {
      fprintf(fp, "%i\t", bestX[i]);
   }
   fprintf(fp, "\t%.3f\n", bestU);
}

void write_times(double t_start, double t_matrix, double t_finish) {
   char *filename_buffer = (char *) calloc(sizeof(char), 1000);
   sprintf(filename_buffer, "results/3_%i.tsv", world_size);
   FILE *fp = fopen(filename_buffer, "a+");

   // FILE *fp = stdout;

   fprintf(fp, "%i\t%f\t%f\t%f\n",
         world_size,
         t_matrix - t_start,
         t_finish - t_matrix,
         t_finish - t_start);
   
   fclose(fp);
}

void printX(int *X) {
   for(int ix = 0; ix < numX; ++ix) { printf("%i ", X[ix]); }
   printf("\n");
}
