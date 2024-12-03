#include <mpi.h>
#include <stdio.h>

typedef struct WRP_Check {
   MPI_Status status;
   int flag;
} WRP_Check;

int WRP_Check_for(int source, int tag, MPI_Comm comm) {
   MPI_Status status;
   int flag;
   MPI_Iprobe(source, tag, comm, &flag, &status);

   return flag;
}

WRP_Check WRP_Check_for_(int source, int tag, MPI_Comm comm) {
   MPI_Status status;
   int flag;
   MPI_Iprobe(source, tag, comm, &flag, &status);

   WRP_Check check = {
      status = status,
      flag = flag,
   };

   return check;
}

