# Compiler
GCC=g++

#PATH_SQUIDS=/Users/carguelles/Workspace/SQuIDS/SQuIDS/
PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -L$(LIBDIR) -static
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR)


# Compilation rules
all: vacuum.o main.o out

vacuum.o: vacuum.cpp
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) vacuum.cpp 
main.o: main.cpp
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) main.cpp 
out: main.o
	$(GCC) $(CFLAGS) $(INCCFLAGS) main.o vacuum.o -lSQUIDS  -lgsl -lgslcblas   -o out
clean:
	rm -f *.o *.fo out
