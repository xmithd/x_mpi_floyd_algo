#include "utils.h"

int floyd_parallel_bcast(matrix2d *graph, proc_info *info);

int floyd_parallel_pipeline(matrix2d *graph, proc_info *info);

int map_local_row_to_global(int local, int sqrt_p, int pid);

int map_local_column_to_global(int local, int sqrt_p, int pid);

int map_process_coord_to_pid(int i, int j);
