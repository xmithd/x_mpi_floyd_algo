#include "utils.h"
#include <string.h>
#include <limits.h>

/**
 * 'zeroes' the list
 */ 
static int init_list(array_list * list)
{
  if (list) {
    memset(list, 0, sizeof(array_list));
    return CODE_SUCCESS;
  }
  return CODE_ERROR;
}

/**
 * Initializes the array list data structure (at least to fit one element)
 * @param list pointer to the array_list to initialize
 * @param capacity The number of elements to allocate for
 * @param element_size the size of one element.
 */ 
int array_list_init(array_list * list, size_t const capacity, size_t const element_size)
{
  if (init_list(list) != CODE_SUCCESS || element_size == 0) {
    return CODE_ERROR;
  }
  list->data = (char *) malloc((capacity == 0 ? 1 : capacity) *element_size);
  if (list->data) {
    list->capacity = capacity == 0 ? 1 : capacity;
    list->element_size = element_size;
    return CODE_SUCCESS;
  }
  return CODE_ERROR;
}

/**
 * Releases the memory allocated by the buffer.
 * @param list an initialized array_list
 */
int array_list_free(array_list * list)
{
  if (list->data) {
    free(list->data);
  }
  memset(list, 0, sizeof(array_list));
  return CODE_SUCCESS;
}

/**
 * Increases the capacity of the buffer by two,
 * copies the elements to the new buffer
 * and releases the old buffer
 */ 
static int resize(array_list * list)
{
  if (list->data && list->capacity > 0) {
    // double the size
    size_t newSize = list->capacity * 2 * list->element_size;
    char * newBuffer = (char *) malloc(newSize);
    if (!newBuffer)
      return CODE_ERROR;
    memset(newBuffer, 0, newSize);
    memmove(newBuffer, list->data, newSize / 2);
    free(list->data);
    list->data = newBuffer;
    list->capacity *= 2;
    return CODE_SUCCESS;
  }
  return CODE_ERROR;
}

/**
 * Clears the buffer
 */  
static int clear(array_list *list)
{
  if (!list)
    return CODE_ERROR;
  memset(list->data, 0, list->size*list->element_size);
  list->size = 0;
  return CODE_SUCCESS;
}

/**
 * Inserts the an element at the end of the list.
 * Resizes the buffer (capacity is doubled) if necessary.
 * @param data: pointer to the data to insert.
 * The size of the element is given by list->element_size.
 */  
int array_list_insert(array_list * list, char const *data)
{
  return array_list_insert_at(list, data, list->size);
}

/**
 * Inserts an element at the specified index.
 * Resizes the buffer is necessary.
 */ 
int array_list_insert_at(array_list *list, char const *data, size_t index)
{
  if (!data || index > list->size) {
    return CODE_ERROR;
  }
  if (list->capacity < (list->size + 1) ) {
    (void) resize(list);
  }
  // shift elements
  if (index < list->size) {
    memmove(&(list->data[(index+1)*(list->element_size)]), &(list->data[(index)*(list->element_size)]), (list->size - index)*list->element_size);
  }
  // insert element
  memmove(&(list->data[index*list->element_size]), data, list->element_size);
  // increment size:
  ++(list->size);
  return CODE_SUCCESS;
}

/**
 * Copies the @param input array into @param output.
 */ 
int array_list_copy(array_list const * input, array_list *output)
{
  int i;
  clear(output);
  for(i=0; i < input->size; ++i) {
    if (array_list_insert(output, array_list_get_at(input, i)) == CODE_ERROR)
      return CODE_ERROR;
  }
  return CODE_SUCCESS;
}

/**
 * Performs a union operation on the sequence.
 * lhs U rhs
 * and stores the result in @param output.
 */ 
int array_list_union(array_list const *lhs, array_list const* rhs, array_list *output)
{
  int i;
  clear(output);
  // first copy lhs
  for (i=0; i < lhs->size; ++i) {
    if (array_list_insert(output, array_list_get_at(lhs, i)) == CODE_ERROR)
      return CODE_ERROR;
  }
  // then just continue inserting rhs
  for (i=0; i < rhs->size; ++i) {
    if (array_list_insert(output, array_list_get_at(rhs, i)) == CODE_ERROR)
      return CODE_ERROR;
  }
  return CODE_SUCCESS;
}

/**
 * Returns a char * pointer (that can be cast) to the 
 * data specified at the given index.
 * Returns NULL if index is invalid.
 */ 
char* array_list_get_at(array_list const * const list, size_t index)
{
  return &(list->data[index * list->element_size]);
}

/**
 * Copies the buffer into the given array_list (deletes any content that may be present)
 * Only copies the provided number of elements
 * with size @param element_size.
 */ 
int array_list_from_buffer(array_list * list, char * buffer, size_t size, size_t element_size)
{
  if (list->data && list->size > 0) {
    free(list->data);
  }
  if (size > 0 && element_size > 0) {
    list->data = (char *) malloc(size*element_size);
    if (!list->data)
      return CODE_ERROR;
    memcpy(list->data, buffer, size*element_size);
    list->capacity = size;
    list->size = size;
    list->element_size = element_size;

    return CODE_SUCCESS;
  } else {
    return CODE_ERROR;
  }
}

/**
 * Simple function to return the log in base 2.
 * returns log2(number) (number must be a multiple of 2)
 * Otherwise, it returns -1.
 */
int log2int(int number)
{
  int count = 0;

  if (number == 1)
    return 0;
  // Invalid if negative, zero or greater than one and uneven
  else if (number < 1 || (number > 1 && number & 0x01)) 
    return -1;
  while ((number = number >> 1)) {
    ++count;
    if (number == 1)
      return count;
    if (number & 0x01) // invalid if first bit to the right is set and number is greater than 1 
      return -1;
  }
  return -1;
}

#ifdef PRINT_DEBUG
#include <string.h>
/**
 * Helper function to print the contents of the @param list
 * into the @param buffer.
 * Note: unsafe (no checks on buffer size)
 */ 
void array_list_print_int(array_list const * list, char *buffer)
{
  int i;
  if (list) {
    for (i=0; i < list->size; ++i) {
      sprintf(&buffer[strlen(buffer)], "%d ", *(int *)array_list_get_at(list, i));
    }
    sprintf(&buffer[strlen(buffer)], "\n");
  }
}

/**
 * Similar to above but puts a new line after every column.
 */ 
void matrix2d_print_int(matrix2d const * matrix, char *buffer)
{
  int i, j;
  if (matrix) {
    for (i=0; i < matrix->rows; ++i) {
      for (j=0; j < matrix->columns; ++j) {
        int val =  *(int *)matrix2d_get_at(matrix, i, j);
        if (val == INT_MAX) {
          sprintf(&buffer[strlen(buffer)], "inf ");
        } else {
          sprintf(&buffer[strlen(buffer)], "%d ", val);
        }
      }
      sprintf(&buffer[strlen(buffer)], "\n");
    }
  }
}

#endif


/**
 * Initializes the matrix by allocating the necessary memory.
 */ 
int matrix2d_init(matrix2d *matrix, size_t rows, size_t columns, size_t element_size)
{
  int rc;
  if (!matrix || rows == 0 || columns == 0)
    return CODE_ERROR;
  rc = array_list_init(&(matrix->list), rows*columns, element_size);
  if (rc != CODE_SUCCESS) {
    return rc;
  }
  matrix->rows = rows;
  matrix->columns = columns;
  // filled with zeroes.
  matrix->list.size = matrix->list.capacity;
  return rc;
}

/**
 * Initializes the matrix by creating a shallow copy of the list structure
 * Matrix becomes the owner of the list.
 */ 
int matrix2d_init_from_list(matrix2d *matrix, size_t rows, size_t columns, array_list const * list)
{
  if (!matrix || !list || rows == 0 || columns == 0)
    return CODE_ERROR;
  if (rows*columns != list->size)
    return CODE_ERROR;
  matrix->list = *list;
  matrix->rows = rows;
  matrix->columns = columns;
  return CODE_SUCCESS;
}



/**
 * Releases the allocated buffer
 */ 
int matrix2d_free(matrix2d *matrix)
{
  if (!matrix) {
    return CODE_ERROR;
  }
  matrix->rows = 0;
  matrix->columns = 0;
  return array_list_free(&(matrix->list));
}

/**
 * Returns a pointer to the data at the specified row and column
 */ 
char* matrix2d_get_at(matrix2d const *matrix, size_t row, size_t column)
{
  // no checks to get better performance.
  return array_list_get_at(&matrix->list, column + matrix->columns*row);
}

/**
 * Copies the specified value at the given location
 */ 
int matrix2d_set_at(matrix2d *matrix, size_t row, size_t column, const char * value)
{
  char *data = array_list_get_at(&(matrix->list), column + matrix->columns*row);
  if (!data)
    return CODE_ERROR;
  memcpy(data, value, matrix->list.element_size);
  return CODE_SUCCESS;
}

/**
 * Gets a pointer to the 'raw' buffer used by the matrix
 */ 
char * matrix2d_get_buffer(matrix2d const *matrix)
{
  if (!matrix)
    return NULL;
  return matrix->list.data;
}

/**
 * Size of the buffer (how many elements are stored)
 */ 
size_t matrix2d_get_buffer_size(matrix2d const *matrix)
{
  return matrix->list.size;
}

/**
 * Generates a directed weighted graph (represented as adjacency matrix)
 * with given edge density.
 * If there is no edge, INT_MAX will be used to represent infinity.
 */
int generate_graph(matrix2d *graph, double density, int max_length)
{
  int rc;
  size_t i,j;
  size_t num_max_edges;
  if (!graph || graph->rows != graph->columns)
    return CODE_ERROR;
  num_max_edges = matrix2d_get_buffer_size(graph);
  for (i=0; i < graph->rows; ++i) {
    for(j = 0; j < graph->columns; ++j) {
      int distance = INT_MAX;
      if ( i == j ) {
        distance = 0;
      } else {
        double val = (double)rand() / (double)RAND_MAX;
        val = val < 0 ? ( val * (-1) ) : val;
        if (val < density) { // there is an edge (it could be weight 0 or negative)
          distance = (int)rand() % (max_length + 1);
        }
      }
      rc = matrix2d_set_at(graph, i, j, (char *)(&distance));
      if (rc != CODE_SUCCESS)
        return rc;
    }
  }
  return rc;
}

/**
 * Copies the matrix structure
 */ 
int matrix2d_copy(matrix2d const *in, matrix2d *out)
{
  int rc;
  if (!in || !out)
    return CODE_ERROR;
  rc = array_list_init(&(out->list), in->list.size, in->list.element_size);
  if (rc != CODE_SUCCESS)
    return rc;
  rc = array_list_copy(&(in->list), &(out->list));
  if (rc != CODE_SUCCESS) {
    array_list_free(&(out->list));
    return rc;
  }
  out->rows = in->rows;
  out->columns = in->columns;
  return rc;
}

