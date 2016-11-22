#!/bin/sh

libtoolize --force
aclocal -I config
autoheader
automake --force-missing --add-missing
autoconf
