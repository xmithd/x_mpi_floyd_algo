#include "floyd_parallel.h"
#include "utils.h"
#include "mpi.h"
#include <limits.h>

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

  // TODO
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

  return CODE_SUCCESS;
}

/**
 * Maps local processor's row element to the global matrix based on the process id.
 */ 
int floyd_parallel_pipeline(matrix2d *graph, proc_info *info)
{
  // TODO
  return CODE_SUCCESS;
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
  matrix2d prev = extract_local_matrix(graph, info);
#ifdef PRINT_DEBUG
{
  char buf[4096];
  sprintf(buf, "Computing Floyd for round %d with process %d\n", k, info->id);
  matrix2d_print_int(&prev, buf);
  printf("%s", buf);
}
#endif
  for (i = 0; i < elements; ++i) {
    int tr_j = map_local_column_to_global(0, info->sqrt_p, info->id, info->nodes);
    for (j = 0; j < elements; ++j) {
      int current_shortest_path = *(int *)matrix2d_get_at(&prev, i, j);
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
        return rc;
      }
      ++tr_j;
    }
    ++tr_i;
  }
  matrix2d_free(&prev);
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
        printf("Row %d; Process: %d. Build row data to send:\n", row, info->id);
        array_list_print_int(&data, buffer);
        printf("%s\n", buffer);
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
        printf("Column %d; Process: %d. Build column data to send:\n", column, info->id);
        array_list_print_int(&data, buffer);
        printf("%s\n", buffer);
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
    if (column % info->sqrt_p != info->p_column) {
      rc = insert_column_in_graph(&data, process_id, column, info, graph);
      if (rc != CODE_SUCCESS) {
        printf("T___T\n");
      }
    }
    array_list_free(&data);
  }
  return rc;
}
