#!/bin/bash

for var in "$@"
do
    GSLINC+=`expr "$var" : '\(-GSL_INC=.*\)'`
    GSLLIB+=`expr "$var" : '\(-GSL_LIB=.*\)'`
done

GSLINC_PATH=${GSLINC##-GSL_INC=}
GSLLIB_PATH=${GSLLIB##-GSL_LIB=}

echo '
# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

CURRENT_DIR = $(shell pwd)
PATH_SQUIDS=$(CURRENT_DIR)/..

# Directories

' > ./src/Makefile
if [ -n "$GSLINC" ]; then
echo "LIBDIR=-L${GSLLIB_PATH}" >> ./src/Makefile
echo "INCDIR=-L${GSLINC_PATH}" >> ./src/Makefile
fi
echo '

LIBDIR+=$(PATH_SQUIDS)/lib
INCDIR+=$(PATH_SQUIDS)/inc
SUINCDIR=$(PATH_SQUIDS)/inc/SU_inc
CFLAGS= -O3 -fPIC
LDFLAGS= -lgsl -lgslcblas
INCCFLAGS=-I$(INCDIR) -I$(SUINCDIR)

# Project files
NAME=SQUIDS
STAT_PRODUCT=lib$(NAME).a
DYN_PRODUCT=lib$(NAME)$(DYN_SUFFIX)

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
	DYN_SUFFIX=.so
	DYN_OPT=-shared -Wl,-soname,$(DYN_PRODUCT)
endif
ifeq ($(OS_NAME),Darwin)
  CC=clang
	CXX=clang++
	LD=clang++
	DYN_SUFFIX=.dylib
	DYN_OPT=-dynamiclib -install_name $(PATH_SQUIDS)/lib/$(DYN_PRODUCT)
endif

OBJECTS= const.o SUNalg.o SQUIDS.o

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)
	mv $(DYN_PRODUCT) $(PATH_SQUIDS)/lib/$(DYN_PRODUCT)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)
	mv $(STAT_PRODUCT) $(PATH_SQUIDS)/lib/$(STAT_PRODUCT)

const.o: const.cpp
	$(CXX) -c $(CFLAGS) $(INCCFLAGS) const.cpp
SQUIDS.o: SQUIDS.cpp
	$(CXX) -c $(CFLAGS) $(INCCFLAGS) SQUIDS.cpp
SUNalg.o: SUNalg.cpp
	$(CXX) -c -Warray-bounds $(CFLAGS) $(INCCFLAGS) SUNalg.cpp

.PHONY: clean
clean:
	rm -f *.o *.fo *.so
' >> ./src/Makefile
