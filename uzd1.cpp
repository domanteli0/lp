#include <iostream>
#include <random>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

using namespace std;

//=============================================================================

void GenerateMatrix(int **&A, int rows, int cols) {
   A = new int*[rows];
   for (int i=0; i<rows; i++) {
      A[i] = new int[cols];
      for (int j=0; j<cols; j++) {
         A[i][j] = (double)rand()/RAND_MAX;
      }
   }
}

//-----------------------------------------------------------------------------

void SortRow(int *A, int cols) {
   int t;
   for (int i=0; i<cols-1; i++)
      for (int j=0; j<cols-1; j++)
         if (A[j] > A[j+1]) { t = A[j]; A[j] = A[j+1]; A[j+1] = t; }
}

//=============================================================================

int main() {
   srand(time(NULL));
   int rows = 128000;         // Number of rows
   int cols = 50;            // Number of columns
   int **A;                   // Matrix
   double ts; 
   

   time_t ts1 = clock();
   auto t1 = Clock::now();
   
   // Generate matrix
   GenerateMatrix(A, rows, cols);

   time_t ts2 = clock();

   // Sort rows of the matrix
   for (int i=0; i<rows; i++) {
      SortRow(A[i], cols);
   }

   time_t ts3 = clock();

    cout << ts1 << endl;
    cout << ts2 << endl;
    cout << ts3 << endl;

    cout << "TIME" << endl;
    cout << (double) ts3 / (double) CLOCKS_PER_SEC << endl;

    // ---

    double alpha = (double) ts2 - (double) ts1;
    double beta = (double) ts3 - (double) ts2;
    double all = (double) ts3 - (double) ts1;

    cout << "ALPHA / BETTA " << endl;
    cout << alpha / all << endl;
    cout << beta / all << endl;
}