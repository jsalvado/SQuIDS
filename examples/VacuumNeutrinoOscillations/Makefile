include ../../settings.mk

PATH_SQUIDS=../..
INCCFLAGS+=-I$(PATH_SQUIDS)/$(INCDIR)
LDFLAGS+= -Wl,-rpath -Wl,$(PATH_SQUIDS)/$(LIBDIR) -L$(PATH_SQUIDS)/$(LIBDIR) -lSQuIDS
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
