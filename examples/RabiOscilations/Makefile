include ../../settings.mk

PATH_SQUIDS=../..
# Directories
LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CXXFLAGS= -O3 -std=c++11  
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR) $(GSL_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBDIR) -L$(LIBDIR) -lSQuIDS $(GSL_LDFLAGS)
EXE=rabi.exe

# Compilation rules
all: $(EXE)

rabi.o: rabi.cpp rabi.h
	$(CXX) -c $(CXXFLAGS) $(INCCFLAGS) rabi.cpp 
main.o: main.cpp rabi.h
	$(CXX) -c $(CXXFLAGS) $(INCCFLAGS) main.cpp 
$(EXE): main.o rabi.o
	$(CXX) $(CXXFLAGS) $(INCCFLAGS) main.o rabi.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
