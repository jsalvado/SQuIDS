include ../../settings.mk

PATH_SQUIDS=../..
CXXFLAGS+=-I$(PATH_SQUIDS)/$(INCDIR)
LDFLAGS+= -Wl,-rpath -Wl,$(PATH_SQUIDS)/$(LIBDIR) -L$(PATH_SQUIDS)/$(LIBDIR) -lSQuIDS $(GSL_LDFLAGS)
EXE=collective.exe

# Compilation rules
all: $(EXE)

collective.o: collective.cpp collective.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) collective.cpp 
main.o: main.cpp collective.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) main.cpp 
$(EXE): main.o collective.o
	$(CXX) $(CXXFLAGS) $(CFLAGS) main.o collective.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
