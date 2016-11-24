#include "floyd_parallel.h"
#include "utils.h"
#include "file_ops.h"
#include "mpi.h"
#include "math.h"

int main(int argc, char *argv[]) {

  int rc; // return code
  matrix2d graph; // matrix representation of graph to use in Floyd's algorithm.
  int numprocs; // number of processes
  int sqrt_p; // square root of the number of processes
  int num_graph_nodes; // number of nodes in the graph
  int num_elements; // number of elements in the matrix. (num_graph_nodes ^ 2)

  proc_info info; // useful info for each process

  // timers
  double start, end;

  // information required for all the processes.
  int nodes; // number of nodes in the graph

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));

  sqrt_p = sqrt(numprocs);
  if (sqrt_p * sqrt_p != numprocs) {
    printf("Number of processes must be a perfect square.\n");
    return CODE_ERROR;
  }

  // generate graph if master and argv has a number argument (number of nodes)
  if (info.id == 0) {
    if (argc > 1) {
      nodes = atoi(argv[1]);
      if (nodes < 2 || nodes > MAX_ELEMENTS) {
        printf("Number of nodes in the graph must be greater than 2.\n");
        return CODE_ERROR;
      }
      rc = matrix2d_init(&d, nodes, nodes, sizeof(int));
      if (rc != CODE_SUCCESS) {
        printf("Error initializing matrix data structure.\n");
        return rc;
      }
      rc = generate_graph(&graph, 0.5, MAX_EDGE_LENGTH);
      if (rc != CODE_SUCCESS) {
          matrix2d_free(&graph);
          printf("Error generating graph!\n");
          return rc;
      }
      rc = write_graph(&graph, "input.txt");
      if (rc != CODE_SUCCESS) {
        matrix2d_free(&graph);
        printf("Error writing to input file!\n");
        return rc;
      }
      matrix2d_free(&graph);
    }

    //read input
    rc = read_graph(&graph, "input.txt");
    if (rc != CODE_SUCCESS) {
      printf("Input file could not be read...\n");
      return rc;
    }
    if (graph->rows < numprocs) {
      matrix2d_free(&graph);
      printf("There are too many processors (%d) for just %d elements.\n", numprocs, graph->rows);
      return CODE_ERROR;
    }

    if (graph->rows % numprocs != 0) {
      matrix2d_free(&graph);
      printf("Each processor must have an equal number of elements to process.\n");
      return CODE_ERROR;
    }
#ifdef PRINT_DEBUG
    printf("nodes in input: %zu\n", graphs.rows);
    {
      char buff[4096];
      printf("Input graph is: \n");
      matrix2d_print_int(&d, buff);
      printf("%s", buff);
    }
#endif
  }

  // broadcast the graph to all the processes.
  rc = MPI_Bcast(&nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS) {
    printf("Error broadcasting graph size.\n");
    return rc;
  }

  if (info.id != 0) {
    matrix2d_init(&graph, nodes, nodes, sizeof(int));
  }

  rc = MPI_Bcast(graph.list.data, graph.list.size, MPI_INT, 0, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS) {
    matrix2d_free(&graph);
    printf("Error broadcasting the initial data.\n");
    return rc;
  }
  
  if (info.id == 0) {
    start = MPI_Wtime();
  }

  rc = floyd_parallel_bcast(&graph, &info);

  if (rc != MPI_SUCCESS) {
    printf("Error %d while running parallel Floyd's algorithm.\n", rc);
    goto end;
  }

  // TODO put some kind of barrier here?

  if (info.id == 0) {
    end = MPI_Wtime();
  }
  // TODO gather data

  // TODO put gathered data into graph

  if (info.id == 0) {
    // Master stores the graph in a file.
    rc = write_graph(&graph, "output.txt");
    if (rc != CODE_SUCCESS) {
      printf("Failed to write to file!\n");
    }
  }

end:
  matrix2d_free(&graph);
  return rc;
}

