#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "floyd_serial.h"
#include "utils.h"
#include "file_ops.h"
#include <math.h>
#include <time.h>

// entry point
int main(int argc, char* argv[]) {
  matrix2d d;
  int rc;
  // timer values;
  double start, end;

  srand(time(NULL));

  // generate file
  if (argc > 1)
  {
    int elements = atoi(argv[1]);
    if (elements < 1 || elements > MAX_ELEMENTS)
    {
      printf("Invalid argument.\n");
      return CODE_ERROR;
    }
    rc = matrix2d_init(&d, elements, elements, sizeof(int));
    if (rc != CODE_SUCCESS) {
      printf("Error initializing d matrix.\n");
      return rc;
    }
    rc = generate_graph(&d, 0.5, MAX_EDGE_LENGTH);
    if (rc != CODE_SUCCESS) {
      printf("Error generating graph!\n");
      return rc;
    }
    rc = write_graph(&d, "input.txt");
    if (rc != CODE_SUCCESS) {
      printf("Error writing graph to file.\n");
      return rc;
    }
    matrix2d_free(&d);
  }
  // read input:
  rc = read_graph(&d, "input.txt");
  if (rc != CODE_SUCCESS) {
    printf("Input file could not be read...\n");
    return rc;
  }
#ifdef PRINT_DEBUG
  printf("nodes in input: %zu\n", d.rows);
  {
    char buff[4096];
    printf("Input graph is: \n");
    matrix2d_print_int(&d, buff);
    printf("%s", buff);
  }
#endif
  // start timer here
  start = MPI_Wtime();
  rc = floyd_serial(&d);
  if (rc != CODE_SUCCESS) {
    printf("ERROR %d\nThere was an error running floyd's algorithm\n", rc);
  }
  // stop timer
  end = MPI_Wtime();
#ifdef PRINT_DEBUG
  {
    char buff[4096];
    printf("Output graph is: \n");
    matrix2d_print_int(&d, buff);
    printf("%s", buff);
  }
#endif

  printf("Computed shortest all-pair path lengths for %zu nodes in %f seconds.\n", d.rows, end-start);
  rc = write_graph(&d, "output.txt");
  if (rc != CODE_SUCCESS) {
    printf("Writing to file failed!\n");
  }
  matrix2d_free(&d);

  return rc;
}
