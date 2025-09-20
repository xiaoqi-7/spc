### SPC Code: Shortest Path Counting Code

### compile:
g++ -O3 -fopenmp -std=c++11 stspc.cpp -o run

### run
* txt-to-bin
  ```
  ./run txt-to-bin facebook/
  @1 path
* decompose_bt
  ```
  ./run decompose_bt facebook/ 5 32
  @1 path
  @2 tree-width allowed
  @3 number of threads
* decompose_core
  ```
  ./run decompose_core facebook/ 5 32
  @1 path
  @2 tree-width allowed
  @3 number of threads
* make queries
  ```
  ./run make_queries facebook/ 10 facebook 5
  @1 path
  @2 number of pairs
  @3 query file name
  @4 tree-width allowed
* query
  ```
  ./run query_spc facebook/ 10 facebook 5
  @1 path
  @2 number of pairs
  @3 query file name
  @4 tree-width allowed
