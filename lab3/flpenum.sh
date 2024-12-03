#!/bin/bash
#SBATCH -p main
#SBATCH -n32

module load openmpi
mpicc -o flpenum flpenum.c -lm
seq 5 | xargs -I {} mpirun flpenum