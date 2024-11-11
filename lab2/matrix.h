#include <stdlib.h>
#include <math.h>

int *lengths(int leg_length, int process_count) {
    int area = leg_length * leg_length / 2;

    int small_area = area / process_count;
    int *lenghts = calloc(sizeof(int), process_count + 1);
    for (int ix = 1; ix < process_count; ++ix) {
        lenghts[ix] = (int) sqrt(2 * ( ix ) * small_area);
    }
    lenghts[process_count] = leg_length;

    return lenghts;
}
