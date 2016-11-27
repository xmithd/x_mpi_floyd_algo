#include "floyd_parallel.h"
#include "utils.h"
#include "file_ops.h"
#include "mpi.h"
#include "math.h"
#include <stdio.h>

int main(int argc, char *argv[]) {

  int rc; // return code
  matrix2d graph; // matrix representation of graph to use in Floyd's algorithm.
  matrix2d extracted; // local matrix to each processor
  array_list gathered_list;
  int numprocs; // number of processes
  int sqrt_p; // square root of the number of processes

  proc_info info; // useful info for each process

  // timers
  double start, end;

  // information required for all the processes.
  int nodes; // number of nodes in the graph

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));

  sqrt_p = (int)sqrt(numprocs);
  if (sqrt_p * sqrt_p != numprocs) {
    printf("Number of processes must be a perfect square.\n");
    return CODE_ERROR;
  }
  info.sqrt_p = sqrt_p;

  // generate graph if master and argv has a number argument (number of nodes)
  if (info.id == 0) {
    if (argc > 1) {
      nodes = atoi(argv[1]);
      if (nodes < 2 || nodes > MAX_ELEMENTS) {
        printf("Number of nodes in the graph must be greater than 2.\n");
        return CODE_ERROR;
      }
      rc = matrix2d_init(&graph, nodes, nodes, sizeof(int));
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
    if (graph.list.size < numprocs) {
      matrix2d_free(&graph);
      printf("There are too many processors (%d) for just %zu elements.\n", numprocs, graph.rows);
      return CODE_ERROR;
    }

    if (graph.list.size % numprocs != 0) {
      matrix2d_free(&graph);
      printf("Each processor must have an equal number of elements to process.\n");
      return CODE_ERROR;
    }
    nodes = graph.rows;
#ifdef PRINT_DEBUG
    printf("nodes in input: %zu\n", graph.rows);
    {
      char buff[4096] = {0};
      printf("Input graph is: \n");
      matrix2d_print_int(&graph, buff);
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

  info.nodes = nodes;
  info.p_row = get_p_row(info.id, info.sqrt_p);
  info.p_column = get_p_column(info.id, info.sqrt_p);

#ifdef PRINT_DEBUG
    //printf("Pid: %d Info:\n\tnodes: %d\n\tsqrt_p: %d\n\tp_row: %d\n\tp_column: %d\n", info.id, info.nodes, info.sqrt_p, info.p_row, info.p_column);
#endif

  if (info.id != 0) {
    matrix2d_init(&graph, nodes, nodes, sizeof(int));
  }
#ifdef PRINT_DEBUG
  if (info.id == 0) {
    char buff[4096] = {0};
    printf("Broadcasting the data of size %zu\n", graph.list.size);
    array_list_print_int(&(graph.list), buff);
    //printf("%s\n", buff);
  }
#endif
  rc = MPI_Bcast(graph.list.data, graph.list.size, MPI_INT, 0, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS) {
    matrix2d_free(&graph);
    printf("Error broadcasting the initial data.\n");
    return rc;
  }
#ifdef PRINT_DEBUG
  {
    char buffer[4096] = {0};
    extracted = extract_local_matrix(&graph, &info);
    matrix2d_print_int(&extracted, buffer);
    //printf("Local matrix to process %d:\n%s", info.id, buffer);
  }
#endif
  
  if (info.id == 0) {
#ifdef PRINT_DEBUG
    printf("Initial broadcasat done!\n");
#endif
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
    // prepare to receive data from processes.
    array_list_init(&gathered_list, info.nodes * info.nodes, sizeof(int));
  }
  // extract each processe's part and gather data
  extracted = extract_local_matrix(&graph, &info);
#ifdef PRINT_DEBUG
  {
    char buffer[4096];
    sprintf(buffer, "Process %d's final matrix:\n", info.id);
    matrix2d_print_int(&extracted, buffer);
    //printf("%s", buffer);
  }
#endif
  rc = MPI_Gather(
      extracted.list.data,
      extracted.list.size,
      MPI_INT,
      gathered_list.data,
      extracted.list.size,
      MPI_INT,
      0,
      MPI_COMM_WORLD
      );
  if (rc != MPI_SUCCESS) {
    printf("Error gathering data after computation.\n");
  }

  // put gathered data into graph
  if (info.id == 0) {
    int i;
    gathered_list.size = gathered_list.capacity;
    for (i=0; i < numprocs; ++i) {
      array_list temp;
      size_t elements_per_proc = gathered_list.size / numprocs;
      size_t small_row_size = info.nodes / info.sqrt_p;
      array_list_init(&temp, elements_per_proc, sizeof(int));
      rc = array_list_from_buffer(&temp, array_list_get_at(&gathered_list, i*elements_per_proc), elements_per_proc, sizeof(int));
      if (rc != CODE_SUCCESS) {
        printf("Ugh, failed to split gathered list.\n");
      }
      // reuse extracted matrix
      matrix2d_free(&extracted);
      rc = matrix2d_init_from_list(&extracted, small_row_size, small_row_size, &temp);
      // extracted owns the temp list now.
      if (rc != CODE_SUCCESS) {
        printf("Ugh, failed to use the split list to create a matrix");
      }
      rc = put_local_matrix_in_global(&extracted, i, &info, &graph);
      if (rc != CODE_SUCCESS) {
        printf("Error putting local matrix into global one.\n");
      }
    }
#ifdef PRINT_DEBUG
    {
      char buff[4096] = {0};
      printf("Output graph is: \n");
      matrix2d_print_int(&graph, buff);
      printf("%s", buff);
    }
#endif
    printf("Computation time: %f seconds.\n", end-start);
    // Master stores the graph in a file.
    rc = write_graph(&graph, "output.txt");
    if (rc != CODE_SUCCESS) {
      printf("Failed to write to file!\n");
    }
  }

end:
  if (info.id == 0) {
    array_list_free(&gathered_list);
  }
  matrix2d_free(&extracted);
  matrix2d_free(&graph);
  return rc;
}

