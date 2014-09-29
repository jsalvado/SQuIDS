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
EXE=rabi.exe

# Compilation rules
all: $(EXE)

rabi.o: rabi.cpp rabi.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) rabi.cpp 
main.o: main.cpp rabi.h
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) main.cpp 
$(EXE): main.o rabi.o
	$(GCC) $(CFLAGS) $(INCCFLAGS) main.o rabi.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
