#!/bin/sh

check_pkgconfig(){
	if [ "$CHECKED_PKGCONFIG" ]; then return; fi
	echo "Looking for pkg-config..."
	which pkg-config 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then
		echo "Error: pkg-config not found; you will need to specify library locations manually" 1>&2
		exit 1
	fi
	CHECKED_PKGCONFIG=1
}

find_package(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	check_pkgconfig
	echo "Looking for $PKG..."

	pkg-config --exists $PKG
	if [ "$?" -ne 0 ]; then
		echo "Error: $PKG not installed or not registered with pkg-config" 1>&2
		lowername=`echo $PKG | tr [A-Z] [a-z]`
		echo "Please specify location using the --with-"$lowername"-incdir and --with-"$lowername"-libdir flags" 1>&2
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

VERSION_NUM=100200
VERSION=`echo $VERSION_NUM | awk '{
	major = int($1/100000);
	minor = ($1/100)%1000;
	patch = $1%100;
	print major"."minor"."patch;
}'`

OS_NAME=`uname -s`

GUESS_CC=gcc
GUESS_CXX=g++
GUESS_AR=ar
GUESS_LD=ld
if [ "$OS_NAME" = Darwin ]; then
	GUESS_CC=clang
	GUESS_CXX=clang++
	GUESS_LD=clang++
	DYN_SUFFIX=.dylib
	DYN_OPT='-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)'
elif [ "$OS_NAME" = FreeBSD ]; then
	DYN_SUFFIX=.so
	DYN_OPT='-shared'
else
	DYN_SUFFIX=.so
	DYN_OPT='-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))'
fi

CC=${CC-$GUESS_CC}
CXX=${CXX-$GUESS_CXX}
AR=${AR-$GUESS_AR}
LD=${LD-$GUESS_LD}

CXXFLAGS="$CXXFLAGS -std=c++11"

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
CFLAGS      C compiler flags
CXXFLAGS    C++ compiler flags
LDFLAGS     Linker flags
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

$CXX $CXXFLAGS resources/compiler_test.cpp -o lib/compiler_test.exe >/dev/null 2>&1
RESULT=$?

if [ "$RESULT" -ne 0 ];
then
	if $CXX --version | grep -q 'clang';
	then #maybe libc++ is available and will make compilation succeed
		CXXFLAGS="$CXXFLAGS -stdlib=libc++"
		LDFLAGS="$LDFLAGS -stdlib=libc++"
		$CXX $CXXFLAGS resources/compiler_test.cpp -o lib/compiler_test.exe >/dev/null 2>&1
		RESULT=$?
	fi
fi

if [ "$RESULT" -ne 0 ];
then
	rm -f lib/compiler_test.exe
	echo "Your C++ compiler ($CXX) is not able to compile this library." >&2
	echo "Plese set the CXX environment variable to point to a compiler which supports C++11." >&2
	echo "Version 4.8.1 or newer of gcc or version 3.3 or newer of clang are recommended." >&2
	exit 1
fi
rm -f lib/compiler_test.exe

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

find_package gsl 1.15

CFLAGS=$CFLAGS

if [ ! -d ./lib/ ]; then
    mkdir lib;
fi

echo "Generating config file..."
echo "prefix=$PREFIX" > lib/squids.pc
echo '
libdir=${prefix}/lib
includedir=${prefix}/include

Name: SQuIDS
Description: Evolves quantum mechanical states
URL: https://github.com/jsalvado/SQuIDS' >> lib/squids.pc
echo "Version: $VERSION" >> lib/squids.pc
echo 'Requires: gsl >= 1.15
Libs: -L${libdir} -lSQuIDS
Cflags: -I${includedir}
' >> lib/squids.pc

echo "Generating version header..."
sed -e "s|__SQUIDS_VERSION__|$VERSION_NUM|g" \
    -e "s|__SQUIDS_VERSION_STR__|$VERSION|g" \
    < resources/version.h.in > inc/version.h

echo "Generating makefile..."

echo "# Directories
LIBDIR:=lib
INCDIR:=inc
SUINCDIR:=inc/SU_inc
SRCDIR:=src

# Compiler
CC:=$CC
CXX:=$CXX
AR:=$AR
LD:=$LD

GSL_CFLAGS=$GSL_CFLAGS
GSL_LDFLAGS=$GSL_LDFLAGS

CFLAGS:=$CFLAGS "'-O3 -fPIC -I$(INCDIR) -I$(SUINCDIR) $(GSL_CFLAGS)'"
CXXFLAGS:=$CXXFLAGS
LDFLAGS:=$LDFLAGS "'-Wl,-rpath -Wl,$(LIBDIR) -L$(LIBDIR) $(GSL_LDFLAGS)' > settings.mk

echo "include settings.mk

VERSION:=$VERSION
PREFIX:=$PREFIX

DYN_SUFFIX:=$DYN_SUFFIX
DYN_OPT=$DYN_OPT
" > Makefile
echo '
# Project files
NAME:=SQuIDS
STAT_PRODUCT:=$(LIBDIR)/lib$(NAME).a
DYN_PRODUCT:=$(LIBDIR)/lib$(NAME)$(DYN_SUFFIX)

OBJECTS:= $(LIBDIR)/const.o $(LIBDIR)/SUNalg.o $(LIBDIR)/SQuIDS.o $(LIBDIR)/MatrixExp.o

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking `basename $(DYN_PRODUCT)`
	@$(CXX) $(DYN_OPT) $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking `basename $(STAT_PRODUCT)`
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

$(LIBDIR)/const.o: $(SRCDIR)/const.cpp $(INCDIR)/const.h Makefile
	@echo Compiling const.cpp to const.o
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(SRCDIR)/const.cpp -o $@
$(LIBDIR)/SQuIDS.o: $(SRCDIR)/SQuIDS.cpp $(INCDIR)/SQuIDS.h $(INCDIR)/SUNalg.h $(INCDIR)/const.h Makefile
	@echo Compiling SQuIDS.cpp to SQuIDS.o
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(SRCDIR)/SQuIDS.cpp -o $@
$(LIBDIR)/SUNalg.o: $(SRCDIR)/SUNalg.cpp $(INCDIR)/SUNalg.h $(INCDIR)/const.h Makefile
	@echo Compiling SUNalg.cpp to SUNalg.o
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(SRCDIR)/SUNalg.cpp -o $@
$(LIBDIR)/MatrixExp.o: $(SRCDIR)/MatrixExp.cpp $(INCDIR)/SUNalg.h  Makefile
	@echo Compiling MatrixExp.cpp to MatrixExp.o
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(SRCDIR)/MatrixExp.cpp -o $@

.PHONY: clean install doxygen docs test check
clean:
	@echo Erasing generated files
	@rm -f $(LIBDIR)/*.o $(LIBDIR)/*.a $(LIBDIR)/*.so $(LIBDIR)/*.dylib

doxygen:
	@doxygen src/doxyfile
docs:
	@doxygen src/doxyfile

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd test ; \
	./run_tests
check: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd test ; \
	./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/SQuIDS
	@mkdir -p $(PREFIX)/include/SQuIDS
	@cp $(INCDIR)/*.h $(PREFIX)/include/SQuIDS
	@mkdir -p $(PREFIX)/include/SQuIDS/detail
	@cp $(INCDIR)/detail/*.h $(PREFIX)/include/SQuIDS/detail
	@mkdir -p $(PREFIX)/include/SQuIDS/SU_inc
	@cp $(INCDIR)/SU_inc/*.h $(PREFIX)/include/SQuIDS/SU_inc
	@cp $(INCDIR)/SU_inc/*.txt $(PREFIX)/include/SQuIDS/SU_inc
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp $(LIBDIR)/squids.pc $(PREFIX)/lib/pkgconfig
' >> Makefile

echo "
export CXX=\"${CXX}\"
export CFLAGS=\"${CFLAGS} ${GSL_CFLAGS}\"
export CXXFLAGS=\"${CXXFLAGS}\"
export LDFLAGS=\"${LDFLAGS} ${GSL_LDFLAGS}\"
" > test/env_vars.sh

echo "Done."
echo "To build, run 'make'"
