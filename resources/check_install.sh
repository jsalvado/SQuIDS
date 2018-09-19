#!/bin/sh

PACKAGE_NAME=$1
INSTALL_PREFIX=$2

which pkg-config 2>&1 > /dev/null
if [ "$?" -ne 0 ]; then
	#don't have pkg-config, so none of this matters
	exit
fi

#See whether pkg-cofig finds the just-installed version
LOOKUP_PREFIX=`pkg-config --variable prefix $PACKAGE_NAME`

if [ "$LOOKUP_PREFIX" != "$INSTALL_PREFIX" ]; then
	#It does not, so see if it's fixable
	export PKG_CONFIG_PATH="$INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
	LOOKUP_PREFIX=`pkg-config --variable prefix $PACKAGE_NAME`
	if [ "$LOOKUP_PREFIX" = "$INSTALL_PREFIX" ]; then
		#That works, so tell the user about it
		echo "NOTE: This installation is not found by pkg-config by default; use"
		echo '      export PKG_CONFIG_PATH="'$INSTALL_PREFIX/lib/pkg-config':$PKG_CONFIG_PATH"'
		echo "      if you want it to be."
	fi
	
fi
