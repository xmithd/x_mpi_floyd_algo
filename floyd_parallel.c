#include "floyd_parallel.h"
#include "utils.h"
#include "mpi.h"
#include "file_ops.h"
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/**
 * Main parallel Floyd algorithm with the broadcast
 */ 
int floyd_parallel_bcast(matrix2d *graph, proc_info *info)
{
  MPI_Comm ROW_COMM;
  MPI_Comm COL_COMM;
  int rc, k;

  // we need 2*sqrt_p communicators. One for every row and column
  MPI_Comm_split(MPI_COMM_WORLD,
      info->p_row, // row groups all the processes in this row
      info->p_column, // rank becomes the column
      &ROW_COMM);

  MPI_Comm_split(MPI_COMM_WORLD,
      info->p_column + info->sqrt_p, // groups all the processes in the column
      info->p_row, // rank becomes the row
      &COL_COMM);

  // Algorithm 10.4
  /* 
    for k 0 to nodes
      -get segment of row k belonging to pid
      -bcast to P*,j (processes in my row)
      -get segment of column k belonging to pid
      -bcast to Pi,* (processes in my column)
      --wait to receive segments of rows k belonging to other processes
      -- put these in my matrix
      compute floyd algo for my part.
  */
  for (k=0; k < info->nodes; ++k) {
    // broadcast the kth row to all processes in the column
    rc = broadcast_row(graph, k, info, COL_COMM);
    if (rc != CODE_SUCCESS) {
      printf("Something went wrong broadcasting rows for k = %d.\n", k);
    }
    // broadcast the kth column to all the processes in the row
    rc = broadcast_column(graph, k, info, ROW_COMM);
    if (rc != CODE_SUCCESS) {
      printf("Something went wrong broadcasting columns for k = %d.\n", k);
    }
    rc = compute_local_floyd(graph, info, k);
    if (rc != CODE_SUCCESS) {
      printf("Error for k = %d when computing local floyd for process %d.\n", k, info->id);
    }
  }

  MPI_Comm_free(&ROW_COMM);
  MPI_Comm_free(&COL_COMM);

  return rc;
}

/**
 * Maps local processor's row element to the global matrix based on the process id.
 */ 
int floyd_parallel_pipeline(matrix2d *graph, proc_info *info)
{
  //Same communicators as broadcast version
  MPI_Comm ROW_COMM;
  MPI_Comm COL_COMM;
  int rc, k;

  // we may need 2*sqrt_p communicators. One for every row and column
  MPI_Comm_split(MPI_COMM_WORLD,
      info->p_row, // row groups all the processes in this row
      info->p_column, // rank becomes the column
      &ROW_COMM);

  MPI_Comm_split(MPI_COMM_WORLD,
      info->p_column + info->sqrt_p, // groups all the processes in the column
      info->p_row, // rank becomes the row
      &COL_COMM);

  /*
    Algorithm: for k in 0..nodes
    -receive and propagate columnwise
    -receive and propagate rowwise
    -compute local with floyd algorithm
  */
  for (k=0; k < info->nodes; ++k) {
    rc = propagate_row(graph, k, info, COL_COMM);
    if (rc != CODE_SUCCESS) {
      printf("Something happened while propagating row for k = %d!\n", k);
    }

    rc = propagate_column(graph, k, info, ROW_COMM);
    if (rc != CODE_SUCCESS) {
      printf("Something went wrong while propagating column for k = %d!\n", k);
    }

    rc = compute_local_floyd(graph, info, k);
    if (rc != CODE_SUCCESS) {
      printf("Error for k = %d when computing local floyd for process %d.\n", k, info->id);
    }
  }

  MPI_Comm_free(&ROW_COMM);
  MPI_Comm_free(&COL_COMM);
  return rc;
}

/**
 * Maps local processor's row element to the global matrix based on the process id.
 */ 
int map_local_row_to_global(int local, int sqrt_p, int pid, int nodes)
{
  int Pi = pid / sqrt_p;
  return (Pi * (nodes / sqrt_p) + local);
}

/**
 * Maps local processor's column element to the global matrix based on the process id.
 */ 
int map_local_column_to_global(int local, int sqrt_p, int pid, int nodes)
{
  int Pj = pid % sqrt_p;
  return (Pj * (nodes / sqrt_p) + local);
}

/**
 * Get the process id given the row and column this process belongs to.
 */ 
int map_process_coord_to_pid(int i, int j, int sqrt_p)
{
  return i * sqrt_p + j;
}

/**
 * Get the process's row
 * (which row this process belongs to)
 */ 
int get_p_row(int pid, int sqrt_p)
{
  return pid / sqrt_p;
}

/**
 * Get the process's column
 * (which column this process belongs to)
 */ 
int get_p_column(int pid, int sqrt_p)
{
  return pid % sqrt_p;
}

/**
 * Compute the shortest paths after receiving the relevant data from other
 * processes for the @param k row/column.
 */
int compute_local_floyd(matrix2d *graph, proc_info const *info, int k)
{
  int i, j;
  int rc = CODE_SUCCESS;
  int elements = info->nodes / info->sqrt_p;
  int tr_i = map_local_row_to_global(0, info->sqrt_p, info->id, info->nodes);
  // optimization: no need to make a copy for the previous: D(k-1)
  //matrix2d prev = extract_local_matrix(graph, info);
#ifdef PRINT_DEBUG
{
  char buf[4096];
  sprintf(buf, "Computing Floyd for round %d with process %d\n", k, info->id);
  //matrix2d_print_int(&prev, buf);
  //printf("%s", buf);
}
#endif
  for (i = 0; i < elements; ++i) {
    int tr_j = map_local_column_to_global(0, info->sqrt_p, info->id, info->nodes);
    for (j = 0; j < elements; ++j) {
      //int current_shortest_path = *(int *)matrix2d_get_at(&prev, i, j);
      int current_shortest_path = *(int *)matrix2d_get_at(graph, tr_i, tr_j);
      int newly_computed_path, newVal;
      int p_i_k = *(int *)matrix2d_get_at(graph, tr_i, k);
      int p_k_j = *(int *)matrix2d_get_at(graph, k, tr_j);
      // Take care of 'infinity' path lengths to avoid overflow
      if (p_i_k == INT_MAX || p_k_j == INT_MAX) {
        newly_computed_path = INT_MAX;
      } else {
        newly_computed_path = p_i_k + p_k_j;
      }
      newVal = MIN(current_shortest_path, newly_computed_path);
      rc = matrix2d_set_at(graph, tr_i, tr_j, (char const *)(&newVal));
      if (rc != CODE_SUCCESS) {
        //matrix2d_free(&prev);
        return rc;
      }
      ++tr_j;
    }
    ++tr_i;
  }
  //matrix2d_free(&prev);
  return rc;
}

/**
 * Extracts the local matrix from the big one
 * @param graph the global matrix
 */ 
matrix2d extract_local_matrix(matrix2d const *graph, proc_info const *info)
{
  int i, j;
  int num_elements = info->nodes / info->sqrt_p; // n over sqrt(p);
  int tr_i = map_local_row_to_global(0, info->sqrt_p, info->id, info->nodes);
  matrix2d local;
#ifdef PRINT_DEBUG
/*
  int tr_jj = map_local_column_to_global(0, info->sqrt_p, info->id, info->nodes);
  printf("Extracting data... process: %d\n\tElements: %d tr_i: %d, tr_j: %d\n", info->id, num_elements, tr_i, tr_jj);
*/
#endif
  matrix2d_init(&local, num_elements, num_elements, sizeof(int));
  for (i=0; i < num_elements; ++i)
  {
    int tr_j = map_local_column_to_global(0, info->sqrt_p, info->id, info->nodes);
    for (j = 0; j < num_elements; ++j) { 
      int *element = (int *)matrix2d_get_at(graph, tr_i, tr_j);
#ifdef PRINT_DEBUG
      if (!element) {
        printf("(%d, %d) == NULL?\n", tr_i, tr_j);
      }
#endif
      matrix2d_set_at(&local, i, j, (const char *)element);
      ++tr_j;
    }
    ++tr_i;
  }
  return local;
}

/**
 * Sets the values in the global matrix from the smaller one.
 * @param pid is the process from which the local matrix was received
 */
int put_local_matrix_in_global(matrix2d *local, int pid, proc_info *info, matrix2d *global)
{
  int i,j, rc;
  int num_elements = info->nodes / info->sqrt_p;
  int tr_i = map_local_row_to_global(0, info->sqrt_p, pid, info->nodes);
  for (i = 0; i < num_elements; ++i) {
    int tr_j = map_local_column_to_global(0, info->sqrt_p, pid, info->nodes);
    for (j = 0; j < num_elements; ++j) {
      int *element = (int *)matrix2d_get_at(local, i, j);
      rc = matrix2d_set_at(global, tr_i, tr_j, (const char*)element);
      // TODO check return code?
      ++tr_j;
    }
    ++tr_i;
  }
  return rc;
}

/**
 * Builds a list of row elements to send.
 * Of course, list will be initialized and caller is
 * responsible for freeing the memory
 */
array_list build_row_elements_to_send(proc_info const *info, int row, matrix2d *graph)
{
  int i;
  array_list list;
  int num_elements = info->nodes / info->sqrt_p; // n over sqrt(p)
  array_list_init(&list, num_elements, sizeof(int));

  for (i=0; i < num_elements; ++i) {
    int * element;
    int rc;
    int tr_j = map_local_column_to_global(i, info->sqrt_p, info->id, info->nodes);
    element = (int *)matrix2d_get_at(graph, row, tr_j);
    rc = array_list_insert(&list, (char const *)element);
    if (rc != CODE_SUCCESS) {
      printf("Warning: something went wrong inserting row elements in list\n");
    }
  }

  return list;
}

/**
 * Build a list of columns elements to send.
 */
array_list build_column_elements_to_send(proc_info const *info, int column, matrix2d *graph)
{
  int i;
  array_list list;
  int num_elements = info->nodes / info->sqrt_p; // n over sqrt(p)
  array_list_init(&list, num_elements, sizeof(int));

  for (i=0; i < num_elements; ++i) {
    int * element;
    int rc;
    int tr_i = map_local_row_to_global(i, info->sqrt_p, info->id, info->nodes);
    element = (int *)matrix2d_get_at(graph, tr_i, column);
    rc = array_list_insert(&list, (char const *)element);
    if (rc != CODE_SUCCESS) {
      printf("Warning: something went wrong inserting column elements in list\n");
    }
  }

  return list;
}

/**
 * Inserts the elements from the list to the right row based
 * on the sender's process id.
 */
int insert_row_in_graph(array_list const *items, int senders_pid, int row,
    proc_info const *info, matrix2d  *graph)
{
  int i, rc;
  int num_elements = info->nodes / info->sqrt_p;
  if (items->size != num_elements) {
    printf("Warning: the list contains more items than what fits in this row.\n");
  }
  for (i=0; i < num_elements; ++i) {
    int tr_j = map_local_column_to_global(i, info->sqrt_p, senders_pid, info->nodes);
    int *element = (int *)array_list_get_at(items, i);
    rc = matrix2d_set_at(graph, row, tr_j, (char const *)element);
    if (rc != CODE_SUCCESS) {
      printf("Warning: cannot insert in matrix...\n");
    }
  }
  return rc;
}

/**
 * Inserts the elements from the list to the right column based
 * on the sender's process id.
 */
int insert_column_in_graph(array_list const *items, int senders_pid, int column,
    proc_info const *info, matrix2d  *graph)
{
  int i, rc;
  int num_elements = info->nodes / info->sqrt_p;
  if (items->size != num_elements) {
    printf("Warning: the list contains more items than what fits in this column.\n");
  }
  for (i=0; i < num_elements; ++i) {
    int tr_i = map_local_row_to_global(i, info->sqrt_p, senders_pid, info->nodes);
    int *element = (int *)array_list_get_at(items, i);
    rc = matrix2d_set_at(graph, tr_i, column, (char const *)element);
    if (rc != CODE_SUCCESS) {
      printf("Warning: cannot insert in matrix...\n");
    }
  }
  return rc;
}

/**
 * Broadcasts the elements of the given row from the current processor
 * and receives the data from other processors to populate the complete columns
 * of this process
 */ 
int broadcast_row(matrix2d *graph, int row, proc_info const *info, MPI_Comm communicator)
{
  int k, rc;
  int p_senders_row = row / (info->nodes / info->sqrt_p);
  for (k=0; k < info->sqrt_p; ++k) { // iterate through the processes of the row
    array_list data;
    // get the sender's process id for this row of processes.
    int process_id = map_process_coord_to_pid(p_senders_row, k, info->sqrt_p);
    // check if row belongs to this process.
    // build data to send
    if (p_senders_row == info->p_row) {
      data = build_row_elements_to_send(info, row, graph);
#ifdef PRINT_DEBUG
      {
        char buffer[4096] = {0};
        //printf("Row %d; Process: %d. Build row data to send:\n", row, info->id);
        array_list_print_int(&data, buffer);
        //printf("%s\n", buffer);
      }
#endif
    } else {
      // prepare a buffer to receive data
      array_list_init(&data, info->nodes / info->sqrt_p, sizeof(int));
      // set the size now
      data.size = data.capacity;
    }
    // root of the broadcast is the processor's column id that is sending this row
    rc = MPI_Bcast(data.data, data.size, MPI_INT, p_senders_row, communicator);
    if (rc != MPI_SUCCESS) {
      printf("Warning! MPI row %d Broadcast failed.\n", row);
    }
    if (p_senders_row != info->p_row) {
      // update local matrix with values received
      rc = insert_row_in_graph(&data, process_id , row, info, graph);
      if (rc != CODE_SUCCESS) {
        printf("T___T\n");
      }
    }
    array_list_free(&data);
  }
  return rc;
}

/**
 * Broadcasts the elements of the given column from the current processor
 * and receives the data from other processors to populate the complete column
 * of this process
 */ 
int broadcast_column(matrix2d *graph, int column, proc_info const *info, MPI_Comm communicator)
{
  int k, rc;
  int p_senders_column = column / (info->nodes / info->sqrt_p);
  for (k=0; k < info->sqrt_p; ++k) { // iterate through the processes in this column
    array_list data;
    // get the senders' process id for this row of processes.
    int process_id = map_process_coord_to_pid(k, p_senders_column, info->sqrt_p);
    if (p_senders_column == info->p_column) {
      data = build_column_elements_to_send(info, column, graph);
#ifdef PRINT_DEBUG
      {
        char buffer[4096] = {0};
        //printf("Column %d; Process: %d. Build column data to send:\n", column, info->id);
        array_list_print_int(&data, buffer);
        //printf("%s\n", buffer);
      }
#endif
    } else {
      array_list_init(&data, info->nodes / info->sqrt_p, sizeof(int));
      // set the size now
      data.size = data.capacity;
    }
    rc = MPI_Bcast(data.data, data.size, MPI_INT, p_senders_column, communicator);
    if (rc != MPI_SUCCESS) {
      printf("Warning! MPI column %d Broadcast failed.\n", column);
    }
    if (p_senders_column != info->p_column) {
      rc = insert_column_in_graph(&data, process_id, column, info, graph);
      if (rc != CODE_SUCCESS) {
        printf("T___T\n");
      }
    }
    array_list_free(&data);
  }
  return rc;
}

/**
 * Either:
 * -Receives the row needed from adjacent processor and propagate to the next processor
 * if not at boundary.
 * OR if current process owns the row:
 * -Extracts the row to send and sends it to adjacent nodes.
 */
int propagate_row(matrix2d *graph, int row, proc_info const *info, MPI_Comm communicator)
{
  int rc = CODE_SUCCESS;
  int p_senders_row = row / (info->nodes / info->sqrt_p);
  int p_senders_column = info->p_column; // senders column is same as mine for columnwise propagation
  array_list data;
  int sender = 0;
  // get the sender's process id for this row of processes.
  int process_id = map_process_coord_to_pid(p_senders_row, p_senders_column, info->sqrt_p);
    
  // check if row belongs to this process.
  // build data to send
  if (p_senders_row == info->p_row) {
    sender = 1;
    data = build_row_elements_to_send(info, row, graph);
#ifdef PRINT_DEBUG
    {
      char buffer[4096] = {0};
      //printf("Row %d; Process: %d. Build row data to send:\n", row, info->id);
      array_list_print_int(&data, buffer);
      //printf("%s\n", buffer);
    }
#endif
  } else {
    // prepare a buffer to receive data
    array_list_init(&data, info->nodes / info->sqrt_p, sizeof(int));
    // set the size now
    data.size = data.capacity;
  }
  if (!sender) {
    // receive
    int source = (p_senders_row > info->p_row) ? (info->p_row - 1) : (info->p_row + 1);
    int destination = (source > info->p_row) ? (info->p_row - 1) : (info->p_row + 1);
    rc = MPI_Recv(data.data, data.size, MPI_INT, source, 0, communicator, MPI_STATUS_IGNORE);
    if (rc != MPI_SUCCESS) {
      printf("Error receiving data from neighbour process!\n");
    }
    // propagate
    if (destination >= 0 || destination < info->sqrt_p)
    {
      rc = MPI_Send(data.data, data.size, MPI_INT, destination, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error propagating data to neighbour process!\n");
      }
    }
  }
  // propagate both ways if possible
  if (sender) {
    if (info->p_row > 0) { // not the first one
      // send to info->p_row - 1
      rc = MPI_Send(data.data, data.size, MPI_INT, info->p_row-1, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error sending data to process above!\n");
      }
    }
    if (info->p_row < (info->sqrt_p-1)) { // not the last one
      // send to info->p_row + 1
      rc = MPI_Send(data.data, data.size, MPI_INT, info->p_row+1, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error sending data to process below!\n");
      }
    }
  }
    
  if (!sender) {
    // update local matrix with values received
    rc = insert_row_in_graph(&data, process_id , row, info, graph);
    if (rc != CODE_SUCCESS) {
      printf("T___T\n");
    }
  }
  array_list_free(&data);
  return rc;
}

/**
 * Rowwise pipelining: 
 * Either:
 * -Receives the column needed from adjacent processor and propagate to the next processor
 * if not at boundary.
 * OR if current process owns the column
 * -Extracts the column to send and sends it to adjacent nodes.
 */
int propagate_column(matrix2d *graph, int column, proc_info const *info, MPI_Comm communicator)
{
  int rc = CODE_SUCCESS;
  int p_senders_row = info->p_row; // senders row is same as mine for rowwise propagation
  int p_senders_column = column / (info->nodes / info->sqrt_p); 
  array_list data;
  int sender = 0; //original sender (no need to receive if true)
  // get the sender's process id for this row of processes.
  int process_id = map_process_coord_to_pid(p_senders_row, p_senders_column, info->sqrt_p);
    
  // check if row belongs to this process.
  // build data to send
  if (p_senders_column == info->p_column) {
    sender = 1;
    data = build_column_elements_to_send(info, column, graph);
#ifdef PRINT_DEBUG
    {
      char buffer[4096] = {0};
      //printf("Column %d; Process: %d. Build column data to send:\n", column, info->id);
      array_list_print_int(&data, buffer);
      //printf("%s\n", buffer);
    }
#endif
  } else {
    // prepare a buffer to receive data
    array_list_init(&data, info->nodes / info->sqrt_p, sizeof(int));
    // set the size now
    data.size = data.capacity;
  }
  if (!sender) {
    // receive
    int source = (p_senders_column > info->p_column) ? (info->p_column - 1) : (info->p_column + 1);
    int destination = (source > info->p_column) ? (info->p_column - 1) : (info->p_column + 1);
    rc = MPI_Recv(data.data, data.size, MPI_INT, source, 0, communicator, MPI_STATUS_IGNORE);
    if (rc != MPI_SUCCESS) {
      printf("Error receiving data from neighbour process!\n");
    }
    // propagate
    if (destination >= 0 || destination < info->sqrt_p)
    {
      rc = MPI_Send(data.data, data.size, MPI_INT, destination, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error propagating data to neighbour process!\n");
      }
    }
  }
  // propagate both ways if possible
  if (sender) {
    if (info->p_column > 0) { // not the first one
      // send to info->p_column - 1
      rc = MPI_Send(data.data, data.size, MPI_INT, info->p_column-1, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error sending data to process to the left!\n");
      }
    }
    if (info->p_column < (info->sqrt_p-1)) { // not the last one
      // send to info->p_column + 1
      rc = MPI_Send(data.data, data.size, MPI_INT, info->p_column+1, 0, communicator);
      if (rc != MPI_SUCCESS) {
        printf("Error sending data to process to the right!\n");
      }
    }
  }
    
  if (!sender) {
    // update local matrix with values received
    rc = insert_column_in_graph(&data, process_id , column, info, graph);
    if (rc != CODE_SUCCESS) {
      printf("T___T\n");
    }
  }
  array_list_free(&data);
  return rc;
}

/**
 * Entry point.
 * Set pipeline to 0 (false) if use broadcast instead of pipeline
 */ 
int main_floyd_parallel(int argc, char *argv[], int pipeline)
{
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

  // init random seed
  srand(time(NULL));

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

  if (pipeline) {
    rc = floyd_parallel_pipeline(&graph, &info);
  } else {
    rc = floyd_parallel_bcast(&graph, &info);
  }

  if (rc != MPI_SUCCESS) {
    printf("Error %d while running parallel Floyd's algorithm.\n", rc);
    goto end;
  }


  if (info.id == 0) {
    // end timer - no barrier needed
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
