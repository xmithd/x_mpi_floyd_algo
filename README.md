# Parallel Floyd's all-pairs shortest path algorithm

This program measures the time it takes to find the distances between all pair of vertices in a weighted edge graph.

## Build
A Makefile is provided to ease building the programs. To build the applications,
the `mpicc` compiler is required. Use `make` to build.
* The serial version is called floyd_serial
* The parallel broadcast version is called floyd_parallel_bcast
* The parallel pipeline version is called floyd_parallel_pipeline

## Usage
An optional argument is a number. This will make the program generate an `input.txt`
with the specified number of vertices in the graph. The input file is stored as a list with MAX_INT representing infinity.
```
./floyd_serial [numelements]
```
```
mpirun -np [numprocessors] floyd_parallel_bcast [numelements]
```
```
mpirun -np [numprocessors] floyd_parallel_pipeline [numelements]
```
The computed distance matrix is then saved in a file called `output.txt`.  
When no argument is given, the program will attempt to read `input.txt`.
The programs will output the time it took to computed the shortest distance between all pairs of vertices. This does not include
the time it takes to generate/read the input file nor the time to scatter and gather the data among the processes.

# Verification
The serial and parallel versions output the same file given the same input.
You can view the distance matrix by enabling the PRINT_DEBUG flag when compiling the objects.  
Note: Only enable the flag for small inputs or else a buffer will occur and likely crash the program.
It will print the input distance matrix and the output distance matrix with the shortest paths.

# Limitations
The number of processes must be equally divisible by the square of the number of vertices so that each process has an equal portion of the distance matrix.

# Improvements
Performance was not the focus of this task. Enable the optimazation flags (such as -O2) to gain 4x the speed.  
Some optimizations are possible such as using less space for each process.
Another optimization in the pipelined version is to use asynchronous sending and managing the sent buffer memory separately (such as reusing it and releasing after all the computations are done). 
Some code re-factoring can be done in the functions that deal with row and column sending. These functions are very similar.
