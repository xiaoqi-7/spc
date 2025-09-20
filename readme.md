
1 compile:
g++ -O3 -fopenmp -std=c++11 stspc.cpp -o run
2 run
1) txt-to-bin
./run txt-to-bin fb/
@1 path
2) decompose_bt
./run decompose_bt fb/ 100 32
@1 path
@2 tree-width allowed
@3 number of threads
3) decompose_core
./run decompose_core fb/ 100 32
@1 path
@2 tree-width allowed
@3 number of threads
4) query
./run query-dis fb/ 100 100000
@1 path
@2 tree-width allowed
@3 number of threads

