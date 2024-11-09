#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdio.h>

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

void display_results(char *filename, int *bestX);
void write_times(double t_start, double t_matrix, double t_finish);

int world_size;

struct UX {
   double u;
   int *X;
} ux;

int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);

   int X[numX];        // Naujas sprendinys
   int bestX[numX];    // Geriausias rastas sprendinys

   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   loadDemandPoints();			// Duomenu nuskaitymas is failo
   double t_start = getTime(); // Algoritmo vykdymo pradzios laikas

   //----- Atstumu matricos skaiciavimas -------------------------------------
   distanceMatrix = calloc(sizeof(double *), numDP);
   for (int i = 0; i < numDP; i++) {
      distanceMatrix[i] = calloc(sizeof(double), i + 1);
      for (int j = 0; j <= i; j++) {
         distanceMatrix[i][j] = HaversineDistance4(demandPoints[i][0], demandPoints[i][1], demandPoints[j][0], demandPoints[j][1]);
      }
   }

   printf("BARRIER REACHED FROM %i\n", world_rank);
   MPI_Barrier(MPI_COMM_WORLD);

   double t_matrix = getTime();
   printf("Matricos skaiciavimo trukme: %lf\n", t_matrix - t_start);

   //----- Pradines naujo ir geriausio sprendiniu reiksmes -------------------
   for (int i = 0; i < numX; i++) { // Pradines naujo ir geriausio sprendiniu koordinates: [0,1,2,...]
      X[i] = i;
      bestX[i] = i;
   }
   u = evaluateSolution(X); // Naujo sprendinio naudingumas (utility)
   bestU = u;				 // Geriausio sprendinio sprendinio naudingumas

   //----- Visų galimų sprendinių perrinkimas --------------------------------
   bool increased = true;
   while (increased) {
      // TODO: calculate as many X values as there are processes (`world_size`).
      if (world_rank == 0) {
         // for(int ix = 1; ix < world_size; ++ix) {
            increased = increaseX(X, numX - 1, numCL);

            if (!increased) {
               // TODO: send the last X values
               // MPI_Bcast(...); // broadcast `DONE` to all other processes
            }
         // }

         // TODO: send the all X value
      }
      // TODO: send X values to other processes and self
      // MPI_Scatter(..., 0, comm);

      // TODO: receive the individual X value;
      // MPI_Receive(...);

      // if (message_type == DONE) { break; }

      u = evaluateSolution(X);
      // TODO: send u to master (world_rank == 0) process

      if (world_rank == 0) {
         // TODO receive asynchronously `u` and `X` from each process

         if (u > bestU) {
            bestU = u;
            memcpy(bestX, X, sizeof(int) * numX);
         }
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   //----- Rezultatu spausdinimas --------------------------------------------
   if (world_rank == 0) {
      double t_finish = getTime(); // Skaiciavimu pabaigos laikas
      printf("Sprendinio paieskos trukme: %lf\n", t_finish - t_matrix);
      printf("Algoritmo vykdymo trukme: %lf\n", t_finish - t_start);

      display_results("stdout" , bestX);
      display_results("new.dat", bestX);
      // write_times(t_start, t_matrix, t_finish);
   }

   MPI_Finalize();
}

//===== Funkciju implementacijos (siu funkciju LYGIAGRETINTI NEREIKIA) ========

void loadDemandPoints()
{
   FILE *f;
   f = fopen("demandPoints.dat", "r");
   demandPoints = calloc(sizeof(double *), numDP);
   for (int i = 0; i < numDP; i++)
   {
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

void display_results(char *filename, int *bestX) {
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
   sprintf(filename_buffer, "no_matrix_%i.tsv", world_size);
   FILE *fp = fopen(filename_buffer, "a+");

   // FILE *fp = stdout;

   fprintf(fp, "%i\t%f\t%f\t%f\n",
         world_size,
         t_matrix - t_start,
         t_finish - t_matrix,
         t_finish - t_start);
}