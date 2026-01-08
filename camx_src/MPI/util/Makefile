SHELL=/bin/sh

#ifndef MPI_INST
  MPI_INST = /usr/local/mpich3
#endif

LIBNAME = libparlib.a

#ifndef CC
  CC = gcc
#endif

LIB = -L$(MPI_INST)/lib -lmpich

INC = $(MPI_INST)/include

CFLAGS = -mcmodel=medium -O2 -I $(INC)

OBJS= parlib.o

$(LIBNAME) : $(OBJS)
	$(CC) $(CFLAGS) -c $(?:.o=.c)
	ar -rv $@ $?
	chmod 644 $@
clean:
	rm -rf *.o $(LIBNAME)
