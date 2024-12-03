#include <stdio.h>

void fprint_X(FILE *fp, double u, int *X) {
    const int numX = 3;

    for (int ix = 0; ix < numX; ++ix) { fprintf(fp, "%i\t", X[ix]); }
    fprintf(fp, "%lf\n", u);
}

void printf_each_int(int *arr, int size) {
    for(int ix = 0; ix < size; ++ix) {
        printf("%i ", arr[ix]);
    }
}

void display_results(char *filename, int *bestX, double bestU) {
   const int numX = 3;
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