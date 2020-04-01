include ../../settings.mk

PATH_SQUIDS=../..
CXXFLAGS+=-I$(PATH_SQUIDS)/$(INCDIR)
LDFLAGS+= -Wl,-rpath -Wl,$(PATH_SQUIDS)/$(LIBDIR) -L$(PATH_SQUIDS)/$(LIBDIR) -lSQuIDS $(GSL_LDFLAGS)
EXE=vacuum.exe

# Compilation rules
all: $(EXE)

vacuum.o: vacuum.cpp vacuum.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) vacuum.cpp 
main.o: main.cpp vacuum.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) main.cpp 
$(EXE): main.o vacuum.o
	$(CXX) $(CXXFLAGS) $(CFLAGS) main.o vacuum.o $(LDFLAGS)  -o $(EXE)
clean:
	rm -f *.o *.fo *.dat $(EXE)
