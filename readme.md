### Shortest Path Counting System

### Shortest Path Counting System
This project involves the development, implementation, testing, and optimisation of a software system supporting shortest path search functionality based on the technical methods described in the following publication:
* Yiqi Wang, Long Yuan, Wenjie Zhang, Xuemin Lin, Zi Chen, Qing Liu, "Towards Efficient Shortest Path Counting on Billion-Scale Graphs", which published in ICDE2023.

### Technical Environment
The software was developed and evaluated in a Linux-based environment using an Intel Xeon CPU, 384 GB of memory, Debian GNU/Linux 12, and g++ 12.2.0.

### Implementation Details
The system was implemented in C++11 with compiler optimisation enabled using the -O3 flag. Development work included implementation of core software functions, performance-oriented optimisation, execution testing, and validation of system behaviour under large-scale data processing conditions.

### compile:
g++ -O3 -fopenmp -std=c++11 stspc.cpp -o run

### How to Test
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
Note: `qd/` folder is necessary, the generated query pairs are stored in `qd/`. Details are shown in the code.
