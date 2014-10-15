#!/bin/bash

function check_pkgconfig(){
	if [ "$CHECKED_PKGCONFIG" ]; then return; fi
	echo "Looking for pkg-config..."
	which pkg-config 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then
		echo "Error: pkg-config not found; you will need to specify library locations manually" 1>&2
		exit 1
	fi
	CHECKED_PKGCONFIG=1
}

function find_package(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	check_pkgconfig
	echo "Looking for $PKG..."

	pkg-config --exists $PKG
	if [ "$?" -ne 0 ]; then
		echo "Error: $PKG not installed or not registered with pkg-config" 1>&2
		exit 1
	fi
	if [ $# -ge 2 ]; then
		MIN_VERSION=$2
		pkg-config --atleast-version $MIN_VERSION $PKG
		if [ "$?" -ne 0 ]; then
			echo "Error: installed $PKG verson ("`pkg-config --modversion $PKG`") is too old; version >=$MIN_VERSION is required" 1>&2
			exit 1
		fi
	fi
	eval ${VAR_PREFIX}_FOUND=1
	eval ${VAR_PREFIX}_CFLAGS=\"`pkg-config --cflags $PKG`\"
	eval ${VAR_PREFIX}_LDFLAGS=\"`pkg-config --libs $PKG`\"
}

PREFIX=/usr/local
VERSION=1.0.0

OS_NAME=`uname -s`

GUESS_CC=gcc
GUESS_CXX=g++
GUESS_AR=ar
GUESS_LD=ld
if [ "$OS_NAME" = Linux ]; then
	DYN_SUFFIX=.so
	DYN_OPT='-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))'
fi
if [ "$OS_NAME" = Darwin ]; then
	GUESS_CC=clang
	GUESS_CXX=clang++
	GUESS_LD=clang++
	DYN_SUFFIX=.dylib
	DYN_OPT='-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)'
fi

CC=${CC-$GUESS_CC}
CXX=${CXX-$GUESS_CXX}
AR=${AR-$GUESS_AR}
LD=${LD-$GUESS_LD}

HELP="Usage: ./config.sh [OPTION]... 

Installation directories:
  --prefix=PREFIX         install files in PREFIX
                          [$PREFIX]

By default, \`make install' will install all the files in
\`$PREFIX/bin', \`$PREFIX/lib' etc.  You can specify
an installation prefix other than \`$PREFIX' using \`--prefix',
for instance \`--prefix=\$HOME'.

The following options can be used to maunally specify the 
locations of dependencies:
  --with-gsl-incdir=DIR   use the copy of gsl in DIR
  --with-gsl-libdir=DIR   use the copy of gsl in DIR

Some influential environment variables:
CC          C compiler command
CXX         C++ compiler command
AR          Static linker command
LD          Dynamic linker command
" #`

for var in "$@"
do
	if [ "$var" = "--help" -o "$var" = "-h" ]; then
		echo "$HELP"
		exit 0
	fi

	TMP=`echo "$var" | sed -n 's/^--prefix=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then PREFIX="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_LIBDIR="$TMP"; continue; fi
done

if [ "$GSL_INCDIR" -a "$GSL_LIBDIR" ]; then
	echo "Checking manually specified GSL..."
	if [ -d "$GSL_INCDIR/gsl" \
         -a -f "$GSL_INCDIR/gsl/gsl_version.h" \
         -a -d "$GSL_LIBDIR" \
         -a -f "$GSL_LIBDIR/libgsl.a" ]; then
		GSL_FOUND=1
		GSL_CFLAGS="-I$GSL_INCDIR"
		GSL_LDFLAGS="-L$GSL_LIBDIR -lgsl -lgslcblas -lm"
	else
		echo "Warning: manually specifed GSL not found; will attempt auto detection"
	fi
fi

find_package gsl 1.14

if [ ! -d ./lib/ ]; then
    mkdir lib;
fi

echo "Generating config file..."
echo "prefix=$PREFIX" > src/squids.pc
echo '
libdir=${prefix}/lib
includedir=${prefix}/include

Name: SQuIDS
Description: Evolves quantum mechanical states
URL: https://github.com/jsalvado/SQuIDS' >> src/squids.pc
echo "Version: $VERSION" >> src/squids.pc
echo 'Requires: gsl >= 1.14
Libs: -L${libdir} -lSQUIDS
Cflags: -I${includedir}
' >> src/squids.pc

echo "Generating makefile..."
echo "# Compiler
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD

DYN_SUFFIX=$DYN_SUFFIX
DYN_OPT=$DYN_OPT

VERSION=$VERSION
PREFIX=$PREFIX
" > ./src/Makefile
echo '
CURRENT_DIR = $(shell pwd)
PATH_SQUIDS=$(CURRENT_DIR)/..

CXXFLAGS= -std=c++11

# Directories
' >> ./src/Makefile
echo "GSL_CFLAGS=$GSL_CFLAGS" >> ./src/Makefile
echo "GSL_LDFLAGS=$GSL_LDFLAGS" >> ./src/Makefile
echo '
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

doxygen:
	doxygen

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd ../test ; ./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/SQuIDS
	@mkdir -p $(PREFIX)/include/SQuIDS
	@cp $(INCDIR)/*.h $(PREFIX)/include/SQuIDS
	@mkdir -p $(PREFIX)/include/SQuIDS/detail
	@cp $(INCDIR)/detail/*.h $(PREFIX)/include/SQuIDS/detail
	@mkdir -p $(PREFIX)/include/SQuIDS/SU_inc
	@cp $(INCDIR)/SU_inc/*.txt $(PREFIX)/include/SQuIDS/SU_inc
	@echo Installing libraries in $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp squids.pc $(PREFIX)/lib/pkgconfig
' >> ./src/Makefile
echo "Done."
echo "To build, run the following: cd src; make"
