#!/bin/sh

system=`uname -s`

if test $system = Darwin; then
  aclocal -I config && autoconf && autoheader && glibtoolize -f && automake -a -f
else
  aclocal -I config && autoconf && autoheader && libtoolize -f && automake -a -f
fi
