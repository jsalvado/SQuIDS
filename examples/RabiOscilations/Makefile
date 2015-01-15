include ../../settings.mk

PATH_SQUIDS=../..
CXXFLAGS+=-I$(PATH_SQUIDS)/$(INCDIR)
LDFLAGS+= -Wl,-rpath -Wl,$(PATH_SQUIDS)/$(LIBDIR) -L$(PATH_SQUIDS)/$(LIBDIR) -lSQuIDS
EXE=rabi.exe

# Compilation rules
all: $(EXE)

rabi.o: rabi.cpp rabi.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) rabi.cpp 
main.o: main.cpp rabi.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) main.cpp 
$(EXE): main.o rabi.o
	$(CXX) $(CXXFLAGS) $(CFLAGS) main.o rabi.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.dat *.eps  $(EXE)
