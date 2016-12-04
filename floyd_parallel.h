#ifndef FLOYD_PARALLEL_H
#define FLOYD_PARALLEL_H
#include "utils.h"
#include "mpi.h"

typedef struct _proc_info {
  int id; // the process id
  int sqrt_p; // p is the total number of processes. This is the square root of that.
  int nodes; // number of nodes in the graph
  int p_row; // process row
  int p_column; // process column
} proc_info;

/**
 * Main parallel Floyd algorithm with the broadcast
 */ 
int floyd_parallel_bcast(matrix2d *graph, proc_info *info);

/**
 * Main parallel Floyd algorithm with pipelining.
 */ 
int floyd_parallel_pipeline(matrix2d *graph, proc_info *info);

/**
 * Maps local processor's row element to the global matrix based on the process id.
 */ 
int map_local_row_to_global(int local, int sqrt_p, int pid, int nodes);

/**
 * Maps local processor's column element to the global matrix based on the process id.
 */ 
int map_local_column_to_global(int local, int sqrt_p, int pid, int nodes);

/**
 * Get the process id given the row and column this process belongs to.
 */ 
int map_process_coord_to_pid(int i, int j, int sqrt_p);

/**
 * Get the process's row
 * (which row this process belongs to)
 */ 
int get_p_row(int pid, int sqrt_p);

/**
 * Get the process's column
 * (which column this process belongs to)
 */ 
int get_p_column(int pid, int sqrt_p);

/**
 * Compute the shortest paths after receiving the relevant data from other
 * processes for the @param k row/column.
 */
int compute_local_floyd(matrix2d *graph, proc_info const *info, int k);

// utilities to build the data to send/receive.

/**
 * Extracts the local matrix from the big one
 * @param graph the global matrix
 */ 
matrix2d extract_local_matrix(matrix2d const *graph, proc_info const *info);

/**
 * Sets the values in the global matrix from the smaller one.
 * @param pid is the process from which the local matrix was received
 */
int put_local_matrix_in_global(matrix2d *local, int pid, proc_info *info, matrix2d *global);

/**
 * Builds a list of row elements to send.
 */
array_list build_row_elements_to_send(proc_info const *info, int row, matrix2d *graph);

/**
 * Build a list of columns elements to send.
 */
array_list build_column_elements_to_send(proc_info const *info, int column, matrix2d *graph);

/**
 * Inserts the elements from the list to the right row based
 * on the sender's process id.
 */
int insert_row_in_graph(array_list const *items, int senders_pid, int row,
    proc_info const *info, matrix2d  *graph);

/**
 * Inserts the elements from the list to the right column based
 * on the sender's process id.
 */
int insert_column_in_graph(array_list const *items, int senders_pid, int column,
    proc_info const *info, matrix2d  *graph);


// COMMUNICATION operations

/**
 * Broadcasts the elements of the given row from the current processor
 * and receives the data from other processors to populate the complete columns
 * of this process
 */ 
int broadcast_row(matrix2d *graph, int row, proc_info const *info, MPI_Comm communicator);

/**
 * Broadcasts the elements of the given column from the current processor
 * and receives the data from other processors to populate the complete rows
 * of this process
 */ 
int broadcast_column(matrix2d *graph, int column, proc_info const *info, MPI_Comm communicator);

/**
 * Either:
 * -Receives the row needed from adjacent processor and propagate to the next processor
 * if not at boundary.
 * OR if current process owns the row:
 * -Extracts the row to send and sends it to adjacent nodes.
 */
int propagate_row(matrix2d *graph, int row, proc_info const *info, MPI_Comm communicator);

/**
 * Rowwise pipelining: 
 * Either:
 * -Receives the column needed from adjacent processor and propagate to the next processor
 * if not at boundary.
 * OR if current process owns the column
 * -Extracts the column to send and sends it to adjacent nodes.
 */
int propagate_column(matrix2d *graph, int column, proc_info const *info, MPI_Comm communicator);

// MAIN PROGRAM

/**
 * Entry point.
 * Set pipeline to 0 (false) if use broadcast instead of pipeline
 */ 
int main_floyd_parallel(int argc, char *argv[], int pipeline);

#endif // FLOYD_PARALLEL_H
