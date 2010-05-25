#!/bin/sh

aclocal -I config && autoconf && autoheader && libtoolize -f && automake -a -f
