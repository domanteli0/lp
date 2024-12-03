#!/bin/bash
#SBATCH -p main
#SBATCH -n1

g++ -W -Wall -Wextra flpenum.cpp -lm -o flpenum_cpp
seq 5 | xargs -I {} ./flpenum_cpp