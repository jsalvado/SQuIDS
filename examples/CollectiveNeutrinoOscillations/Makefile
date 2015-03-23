include ../../settings.mk

PATH_SQUIDS=../..
CXXFLAGS+=-I$(PATH_SQUIDS)/$(INCDIR)
LDFLAGS+= -Wl,-rpath -Wl,$(PATH_SQUIDS)/$(LIBDIR) -L$(PATH_SQUIDS)/$(LIBDIR) -lSQuIDS
EXE=collective.exe

# Compilation rules
all: $(EXE)

collective.o: collective.cpp collective.h
	$(CXX) -c $(CXXFLAGS) collective.cpp 
main.o: main.cpp collective.h
	$(CXX) -c $(CXXFLAGS) main.cpp 
$(EXE): main.o collective.o
	$(CXX) $(CXXFLAGS) main.o collective.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
