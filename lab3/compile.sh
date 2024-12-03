#!/bin/bash
#SBATCH -p main
#SBATCH -n2
mpicc -W -Wall -Wextra flpenum.c -lm -o flpenum