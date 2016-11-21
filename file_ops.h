#ifndef FILE_OPS_H
#define FILE_OPS_H

#include "utils.h"

/**
 * Generates a file called input.txt with @param numelements
 * random elements between 0 and numelements.
 * If file exists, it will overwrite it.
 */
int generate_file(int numelements);

/**
 * Reads the input.txt file and stores the numbers in a given array_list.
 * The array list should not be initialized.
 * Caller is responsible for freeing the memory
 */
int read_file(array_list *list, const char * filename);

/**
 * Writes to output.txt the contents of @param list
 * overwriting any existing file called output.txt
 */
int write_file(array_list const *list);

/**
 * Reads graph from the given file name.
 * Graph may not be initialised.
 */ 
int read_graph(matrix2d *graph, const char* filename);

/**
 * Writes graph and stores it in file.
 */ 
int write_graph(matrix2d *graph, const char* filename);

#endif // FILE_OPS_H
