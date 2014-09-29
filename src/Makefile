# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))

VERSION=1.0.0
PREFIX=/usr/local


CURRENT_DIR = $(shell pwd)
PATH_SQUIDS=$(CURRENT_DIR)/..

CXXFLAGS= -std=c++11

# Directories

GSL_CFLAGS= 
GSL_LDFLAGS=-lgsl -lgslcblas -lm 

LIBDIR=$(PATH_SQUIDS)/lib
INCDIR=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -fPIC -I$(INCDIR) -I$(SUINCDIR) $(GSL_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBDIR) -L$(LIBDIR) $(GSL_LDFLAGS)

# Project files
NAME=SQUIDS
STAT_PRODUCT=$(LIBDIR)/lib$(NAME).a
DYN_PRODUCT=$(LIBDIR)/lib$(NAME)$(DYN_SUFFIX)

OBJECTS= const.o SUNalg.o SQUIDS.o

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking `basename $(DYN_PRODUCT)`
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking `basename $(STAT_PRODUCT)`
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

const.o: const.cpp $(INCDIR)/const.h Makefile
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) const.cpp
SQUIDS.o: SQUIDS.cpp $(INCDIR)/SQUIDS.h Makefile
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) SQUIDS.cpp
SUNalg.o: SUNalg.cpp $(INCDIR)/SUNalg.h Makefile
	$(CXX) $(CXXFLAGS) -c -Warray-bounds $(CFLAGS) SUNalg.cpp

.PHONY: clean install
clean:
	rm -f *.o *.so *.dylib ../lib/*.a ../lib/*.so ../lib/*.dylib

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/SQuIDS
	@mkdir -p $(PREFIX)/include/SQuIDS
	@cp $(INCDIR)/*.h $(PREFIX)/include/SQuIDS
	@echo Installing libraries in $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp squids.pc $(PREFIX)/lib/pkgconfig

