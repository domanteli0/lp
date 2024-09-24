run: compile
    ./uzd3

compile:
    /opt/homebrew/opt/llvm/bin/clang++ -Xpreprocessor -fopenmp uzd3.cpp -lomp -o uzd3 \
    -W -Wall -Wextra \
    -I /opt/homebrew/opt/llvm/include \
    -L /opt/homebrew/opt/llvm/lib

# compile:
#     /opt/homebrew/opt/llvm/bin/clang -v -Xpreprocessor -fopenmp uzd3.cpp -lomp -o uzd3 \
#     -I /opt/homebrew/include \
#     -I /opt/homebrew/opt/llvm/include \
#     -L /opt/homebrew/lib \
#     -L /opt/homebrew/opt/llvm/lib