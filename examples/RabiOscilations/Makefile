# Compiler
GCC=g++

PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -L$(LIBDIR) -static 
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR) 


# Compilation rules
all: rabi.o main.o prog

rabi.o: rabi.cpp 
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) rabi.cpp 
main.o: main.cpp 
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) main.cpp 
prog: main.o 
	$(GCC) $(CFLAGS) $(INCCFLAGS) main.o rabi.o -lSQUIDS  -lgsl -lgslcblas   -o out.exe
clean:
	rm -f *.o  out.exe
