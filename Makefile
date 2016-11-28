OBJS_SERIAL = floyd_serial.o main_serial.o utils.o file_ops.o
OBJS_PARALLEL_BCAST = utils.o floyd_parallel.o main_parallel_bcast.o file_ops.o
OBJS_PARALLEL_PIPELINE = utils.o floyd_parallel.o main_parallel_pipeline.o file_ops.o
CC = mpicc
DEBUG = #-DPRINT_DEBUG -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : floyd_serial floyd_parallel_bcast floyd_parallel_pipeline


floyd_serial : $(OBJS_SERIAL)
	$(CC) $(LFLAGS) $(OBJS_SERIAL) -o floyd_serial

floyd_parallel_bcast : $(OBJS_PARALLEL_BCAST)
	$(CC) $(LFLAGS) $(OBJS_PARALLEL_BCAST) -o floyd_parallel_bcast

floyd_parallel_pipeline : $(OBJS_PARALLEL_PIPELINE)
	$(CC) $(LFLAGS) $(OBJS_PARALLEL_BCAST) -o floyd_parallel_pipeline

main_parallel_bcast.o : main_parallel_bcast.c floyd_parallel.h floyd_parallel.o
	$(CC) $(CFLAGS) main_parallel_bcast.c

main_parallel_pipeline.o : main_parallel_pipeline.c floyd_parallel.h floyd_parallel.o
	$(CC) $(CFLAGS) main_parallel_pipeline.c

main_serial.o : main_serial.c floyd_serial.h utils.h utils.o floyd_serial.o
	$(CC) $(CFLAGS) main_serial.c

floyd_serial.o : floyd_serial.h utils.h utils.o floyd_serial.c
	$(CC) $(CFLAGS) floyd_serial.c

file_ops.o: utils.o file_ops.c file_ops.h
	$(CC) $(CFLAGS) file_ops.c

utils.o : utils.h utils.c
	$(CC) $(CFLAGS) utils.c

clean :
	rm *.o floyd_serial floyd_parallel_bcast floyd_parallel_pipeline
