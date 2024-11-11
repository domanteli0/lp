# mpirun --timeout 10 --report-state-on-timeout -n $1 ./flpenum 
mpirun -n $1 ./flpenum 
diff --color base.dat new.dat
# ./compile.sh && seq 5 | xargs -I{} ./run_compare.sh 9