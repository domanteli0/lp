g++ -W -Wall -Wextra flpenum.cpp -lm -o flpenum_cpp
seq 1 | xargs -I{} ./flpenum_cpp
diff --color base.dat new.dat