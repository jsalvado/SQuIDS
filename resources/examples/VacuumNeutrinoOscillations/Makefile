# Compiler
CXX ?= g++

PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CXXFLAGS= -O3 -std=c++11
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR)
LDFLAGS= -Wl,-rpath -Wl,$(LIBDIR) -L$(LIBDIR) -lSQuIDS -lgsl -lgslcblas
EXE=vacuum.exe

# Compilation rules
all: $(EXE)

vacuum.o: vacuum.cpp vacuum.h
	$(CXX) -c $(CXXFLAGS) $(INCCFLAGS) vacuum.cpp 
main.o: main.cpp vacuum.h
	$(CXX) -c $(CXXFLAGS) $(INCCFLAGS) main.cpp 
$(EXE): main.o vacuum.o
	$(CXX) $(CXXFLAGS) $(INCCFLAGS) main.o vacuum.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.fo *.dat $(EXE)
