#include <stdio.h>

void fprint_X(FILE *fp, double u, int *X) {
    // const int numX = 3;

    // for (int ix = 0; ix < numX; ++ix) { fprintf(fp, "%i\t", X[ix]); }
    // fprintf(fp, "%lf\n", u);
}

void printf_each_int(int *arr, int size) {
    for(int ix = 0; ix < size; ++ix) {
        printf("%i ", arr[ix]);
    }
}