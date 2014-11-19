# Compiler
GCC=g++ 

PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -L$(LIBDIR) -std=c++11  
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR) 
LDFLAGS= -Wl,-rpath -Wl,$(LIBDIR) -lSQUIDS -lgsl -lgslcblas
EXE=collective.exe

# Compilation rules
all: $(EXE)

collective.o: collective.cpp collective.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) collective.cpp 
main.o: main.cpp collective.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) main.cpp 
$(EXE): main.o collective.o
	$(GCC) $(CFLAGS) $(INCCFLAGS) main.o collective.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
