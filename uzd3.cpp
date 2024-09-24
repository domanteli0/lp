#include <iostream>
#include <random>
#include <omp.h>
#include <assert.h>
#include <chrono>

using namespace std;

//=============================================================================

void GenerateMatrix(int **&A, int rows, int cols)
{
    A = new int *[rows]; // A[ROW][COL]

    for (int i = 0; i < rows; i++)
    {
        A[i] = new int[cols];
        for (int j = 0; j < cols; j++)
        {
            A[i][j] = (double)rand() / RAND_MAX;
        }
    }
}

//-----------------------------------------------------------------------------

void SortRow(int *A, int cols)
{
    int t;
    for (int i = 0; i < cols - 1; i++)
        for (int j = 0; j < cols - 1; j++)
            if (A[j] > A[j + 1])
            {
                t = A[j];
                A[j] = A[j + 1];
                A[j + 1] = t;
            }
}

//=============================================================================

int main()
{
    srand(time(NULL));
    int rows = 128000; // Number of rows
    int cols = 200;    // Number of columns
    int **A;           // Matrix (MUTABLE)
    // double ts;

    int num_threads = 4;
    int chunk = rows / num_threads;

    assert(rows % num_threads == 0);

    omp_set_num_threads(num_threads);

    // Generate matrix
    auto t1 = chrono::high_resolution_clock::now();
    GenerateMatrix(A, rows, cols);
    auto t2 = chrono::high_resolution_clock::now();
    
    // Sort rows of the matrix
#pragma omp parallel shared(chunk, A, rows, cols)
    {
        // for (int i = 0; i < rows; i++)
        int id = omp_get_thread_num();
        for (int i = id * chunk; i < (id + 1) * chunk; i++)
        {
            // A[i] IS MUTATED HERE
            SortRow(A[i], cols);
        }
    }
    auto t3 = chrono::high_resolution_clock::now();

    printf("TIME: %lld ms\n", chrono::duration_cast<chrono::milliseconds>(t3 - t2).count());
}