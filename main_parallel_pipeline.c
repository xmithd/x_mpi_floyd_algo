#include "floyd_parallel.h"
#include "mpi.h"

int main(int argc, char *argv[])
{
  int rc; // return code
  rc = main_floyd_parallel(argc, argv, 1);
  MPI_Finalize();
  return rc;
}

