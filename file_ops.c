#include "file_ops.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "utils.h"

#define DEFAULT_LIST_SIZE 128

static const char* output_filename = "output.txt";

static int write_to_file(const char* filename, array_list const * content)
{
  FILE *file;
  int i;
  if (!content || !filename) {
    return CODE_ERROR;
  }
  file = fopen(filename, "w");
  if (!file)
    return CODE_ERROR;

  for (i=0; i < content->size; ++i) {
    int *val = (int *) array_list_get_at(content, i);
    fprintf(file, "%d ", *val);
  }
  fclose(file);

  return CODE_SUCCESS;
}

/**
 * Generates a file called input.txt with @param numelements
 * random elements between 0 and numelements.
 * If file exists, it will overwrite it.
 */
int generate_file(int numelements)
{
  const char* the_filename = "input.txt";
  int i, rc;
  array_list list;
  // use time as seed for pseudo random generator
  srand(time(NULL));

  array_list_init(&list, numelements, sizeof(unsigned int));

  for (i=0; i < numelements; ++i) {
    unsigned int num = (unsigned int) (rand()) % numelements;
    array_list_insert(&list, (char *)(&num));
  }

  rc = write_to_file(the_filename, &list);

  array_list_free(&list);

  return rc;
}

/**
 * Reads the input.txt file and stores the numbers in a given array_list.
 * The array list should not be initialized.
 * Caller is responsible for freeing the memory
 */ 
int read_file(array_list *array, const char* the_filename)
{
  int num;
  FILE* file = fopen(the_filename, "r");
  if (!file || !array) {
    fprintf(stderr, "Error reading file %s\n", the_filename);
    return CODE_ERROR;
  }
  if (!array) {
    fclose(file);
    fprintf(stderr, "Memory error!\n");
    return CODE_ERROR;
  }
  (void)array_list_init(array, DEFAULT_LIST_SIZE, sizeof(int));

  // Reading 
  while ( fscanf(file, "%d", &num) == 1) {
    if (ferror(file)) {
      array_list_free(array);
      fclose(file);
      fprintf(stderr, "Error reading file %s\n", the_filename);
      return CODE_ERROR;
    }
    (void)array_list_insert(array, (char const *)&num);
  }

  fclose(file);
  return CODE_SUCCESS;
}

/**
 * Writes to output.txt the contents of @param list
 * overwriting any existing file called output.txt
 */
int write_file(array_list const *list)
{
  return write_to_file(output_filename, list);
}

/**
 * Reads graph from the given file name.
 * graph must not be initialzed.
 */ 
int read_graph(matrix2d *graph, const char* filename)
{
  int rc;
  int nodes;

  if (!graph || !filename) {
    return CODE_ERROR;
  }
  rc = read_file(&(graph->list), filename);
  if (rc != CODE_SUCCESS)
    return rc;

  // check that number of elements is a perfect square (need a square matrix)
  nodes = (int) sqrt((double)(graph->list.size));
  if (nodes*nodes != (int)graph->list.size) {
    fprintf(stderr, "Number of elements while reading graph is not a perfect square.\n");
    matrix2d_free(graph);
    return CODE_ERROR;
  }
  // it is a perfect square
  graph->rows = nodes;
  graph->columns = nodes;

  return CODE_SUCCESS;

}

/**
 * Writes graph and stores it in data structure
 */ 
int write_graph(matrix2d *graph, const char* filename)
{
  if (!graph || !filename) {
    return CODE_ERROR;
  }

  return write_to_file(filename, &(graph->list));
}


