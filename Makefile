OBJS_SERIAL = floyd_serial.o main_serial.o utils.o file_ops.o
#OBJS_PARALLEL = quicksort_serial.o utils.o quicksort_parallel.o main_parallel.o file_ops.o
CC = mpicc
DEBUG = -DPRINT_DEBUG
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

#all : floyd_serial floyd_parallel

all : floyd_serial

floyd_serial : $(OBJS_SERIAL)
	$(CC) $(LFLAGS) $(OBJS_SERIAL) -o floyd_serial

#floyd_parallel : $(OBJS_PARALLEL)
#	$(CC) $(LFLAGS) $(OBJS_PARALLEL) -o floyd_parallel

main_serial.o : main_serial.c floyd_serial.h utils.h utils.o floyd_serial.o
	$(CC) $(CFLAGS) main_serial.c

floyd_serial.o : floyd_serial.h utils.h utils.o floyd_serial.c
	$(CC) $(CFLAGS) floyd_serial.c

file_ops.o: utils.o file_ops.c file_ops.h
	$(CC) $(CFLAGS) file_ops.c

utils.o : utils.h utils.c
	$(CC) $(CFLAGS) utils.c

clean :
	rm *.o floyd_serial
