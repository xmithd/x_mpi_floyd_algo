#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>

// some constants used in the programs

#ifdef PRINT_DEBUG
#define MAX_ELEMENTS (1 << 7)
#else
#define MAX_ELEMENTS (1 << 30)
#endif
#define CODE_SUCCESS 0
#define CODE_ERROR   1

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX_EDGE_LENGTH 10

/**
 * array list data structure
 * similar to std::vector in c++ (but much more basic)
 */
typedef struct _array_list {
  // data can be cast to any type but is stored as char* for convenience.
  char *data; 
  // the maximum number of elements that can be held by the data buffer.
  size_t capacity;
  // the number of elements contained in the data buffer.
  size_t size;
  // the size of a single element
  size_t element_size;
} array_list;

/**
 * Represents a 2D matrix using array_list as internal data structure.
 * Note: it has a fixed size for this project.
 */
typedef struct _matrix2d {
  size_t rows;
  size_t columns;
  array_list list;
} matrix2d;

/**
 * Initializes the array list data structure
 * @param list pointer to the array_list to initialize
 * @param capacity The number of elements to allocate for
 * @param element_size the size of one element.
 */ 
int array_list_init(array_list * list, size_t const capacity, size_t const element_size);

/**
 * Releases the memory allocated by the buffer.
 * @param list an initialized array_list
 */
int array_list_free(array_list * list);

/**
 * Inserts the an element at the end of the list.
 * Resizes the buffer (capacity is doubled) if necessary.
 * @param data: pointer to the data to insert.
 * The size of the element is given by list->element_size.
 */ 
int array_list_insert(array_list * list, char const *data);

/**
 * Inserts an element at the specified index.
 * Resizes the buffer is necessary.
 */ 
int array_list_insert_at(array_list *list, char const *data, size_t index);

/**
 * Copies the @param input array into @param output.
 */ 
int array_list_copy(array_list const * input, array_list *output);

/**
 * Performs a union operation on the sequence.
 * lhs U rhs
 * and stores the result in @param output.
 */ 
int array_list_union(array_list const *lhs, array_list const* rhs, array_list *output);

/**
 * Returns a char * pointer (that can be cast) to the 
 * data specified at the given index.
 * Returns NULL if index is invalid.
 */ 
char* array_list_get_at(array_list const * const list, size_t index);

/**
 * Copies the buffer into the given array_list (deletes any content that may be present)
 * Only copies the provided number of elements
 * with size @param element_size.
 */ 
int array_list_from_buffer(array_list * list, char * buffer, size_t size, size_t element_size);

#ifdef PRINT_DEBUG
#include <stdio.h>
/**
 * Helper function to print the contents of the @param list
 * into the @param buffer.
 * Note: unsafe (no checks on buffer size)
 */ 
void array_list_print_int(array_list const * list, char *buffer);

/**
 * Similar to above but puts a new line after every column.
 */ 
void matrix2d_print_int(matrix2d const * matrix, char *buffer);

#endif

/**
 * Simple function to return the log in base 2.
 * returns log2(number) (number must be a multiple of 2)
 * Otherwise, it returns -1.
 */
int log2int(int number);

/**
 * Initializes the matrix by allocating the necessary memory.
 */ 
int matrix2d_init(matrix2d *matrix, size_t rows, size_t columns, size_t element_size);

/**
 * Initializes the matrix by creating a shallow copy of the list structure
 * Matrix becomes the owner of the list.
 */ 
int matrix2d_init_from_list(matrix2d *matrix, size_t rows, size_t columns, array_list const * list);

/**
 * Releases the allocated buffer
 */ 
int matrix2d_free(matrix2d *matrix);

/**
 * Returns a pointer to the data at the specified row and column
 */ 
char* matrix2d_get_at(matrix2d const *matrix, size_t row, size_t column);

/**
 * Copies the specified value at the given location
 */ 
int matrix2d_set_at(matrix2d *matrix, size_t row, size_t column, const char * value);

/**
 * Gets a pointer to the 'raw' buffer used by the matrix
 */ 
char * matrix2d_get_buffer(matrix2d const *matrix);

/**
 * Size of the buffer (how many elements are stored)
 */ 
size_t matrix2d_get_buffer_size(matrix2d const *matrix);

/**
 * Generates a directed weighted graph (represented as adjacency matrix)
 * with given edge density.
 * If there is no edge, INT_MAX will be used to represent infinity.
 * The given matrix should be initialized and with equal number of rows and columns.
 */
int generate_graph(matrix2d *graph, double density, int max_length);

/**
 * Copies the matrix structure. @param out should not be initialized.
 */ 
int matrix2d_copy(matrix2d const *in, matrix2d *out);

#endif // UTILS_H
