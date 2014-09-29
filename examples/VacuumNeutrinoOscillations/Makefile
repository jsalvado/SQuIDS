# Compiler
GCC=g++

#PATH_SQUIDS=/Users/carguelles/Workspace/SQuIDS/SQuIDS/
PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -L$(LIBDIR) -std=c++11
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR)
LDFLAGS= -Wl,-rpath -Wl,$(LIBDIR) -lSQUIDS -lgsl -lgslcblas
EXE=vacuum.exe

# Compilation rules
all: $(EXE)

vacuum.o: vacuum.cpp vacuum.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) vacuum.cpp 
main.o: main.cpp vacuum.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) main.cpp 
$(EXE): main.o vacuum.o
	$(GCC) $(CFLAGS) $(INCCFLAGS) main.o vacuum.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.fo *.dat $(EXE)
