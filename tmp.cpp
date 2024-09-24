#include <omp.h> // <- čia yra visos funkcijos ir direktyvos

int main() {
  double A[1000];

  // jeigu nenurodoma kitaip automatiškai parenka thread'ų skaičių ir
  // tolygiai paskirsto skačiavimus, pvz:
  // 4 thread'ai:
  //   I-as thread'as: 0-254
  //   II-as thread'as: 255-449
  //   ...

  for(int i = 0; i < 1000; i++) {
   //  f(A[i]);
   // double a = A[i];
  }
}