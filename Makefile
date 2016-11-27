OBJS_SERIAL = floyd_serial.o main_serial.o utils.o file_ops.o
OBJS_PARALLEL_BCAST = floyd_serial.o utils.o floyd_parallel.o main_parallel_bcast.o file_ops.o
CC = mpicc
DEBUG = -DPRINT_DEBUG
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : floyd_serial floyd_parallel_bcast


floyd_serial : $(OBJS_SERIAL)
	$(CC) $(LFLAGS) $(OBJS_SERIAL) -o floyd_serial

floyd_parallel_bcast : $(OBJS_PARALLEL_BCAST)
	$(CC) $(LFLAGS) $(OBJS_PARALLEL_BCAST) -o floyd_parallel_bcast

main_serial.o : main_serial.c floyd_serial.h utils.h utils.o floyd_serial.o
	$(CC) $(CFLAGS) main_serial.c

floyd_serial.o : floyd_serial.h utils.h utils.o floyd_serial.c
	$(CC) $(CFLAGS) floyd_serial.c

file_ops.o: utils.o file_ops.c file_ops.h
	$(CC) $(CFLAGS) file_ops.c

utils.o : utils.h utils.c
	$(CC) $(CFLAGS) utils.c

clean :
	rm *.o floyd_serial floyd_parallel_bcast
