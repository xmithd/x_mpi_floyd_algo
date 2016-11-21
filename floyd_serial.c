
#include "floyd_serial.h"
#include "utils.h"
#include <limits.h>

/**
 * Takes the adjacency matrix representing the graph
 * and modifies it to contain the path lengths between the points.
 */ 
int floyd_serial(matrix2d *d)
{
  size_t k, i, j;
  size_t nodes = d->rows;
  matrix2d prev;
  int rc = CODE_SUCCESS;
  if (nodes != d->columns)
    return CODE_ERROR;

  for (k=0; k < nodes; ++k) {
    // get a fresh copy of the previous iteration matrix
    if (matrix2d_copy(d, &prev) != CODE_SUCCESS) {
      return CODE_ERROR;
    }
    for (i = 0; i < d->rows; ++i) {
      for (j = 0; j < d->columns; ++j) {
        /*if (k == 0) {
          rc = matrix2d_set_at(d, i, j, matrix2d_get_at(&prev, i, j));
          if (rc != CODE_SUCCESS) {
            matrix2d_free(&prev);
            return rc;
          }
        } else {*/
          int current_shortest_path = *(int *)matrix2d_get_at(&prev, i, j); 
          int p_i_k =  *(int *)matrix2d_get_at(&prev, i, k);
          int p_k_j =  *(int *)matrix2d_get_at(&prev, k, j);
          int newly_computed_path, newVal;
          // Take care of 'infinity' path lengths to avoid overflow
          if (p_i_k == INT_MAX || p_k_j == INT_MAX) {
            newly_computed_path = INT_MAX;
          } else {
            newly_computed_path = p_i_k + p_k_j;
          }
          newVal = MIN(current_shortest_path, newly_computed_path);
          rc = matrix2d_set_at(d, i, j, (char const *)(&newVal));
          if (rc != CODE_SUCCESS) {
            matrix2d_free(&prev);
            return rc;
          }
        //}
      }
    }
    // release memory allocated for the previous iteration matrix
    matrix2d_free(&prev);
  }
  return rc;
}

