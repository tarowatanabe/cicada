#!/bin/sh

echo "eigen_includedir = \$(pkgincludedir)"
echo
/bin/echo -n "nobase_eigen_include_HEADERS ="

for file in `find Eigen -type f`; do
  basename=`basename $file`
  if test $basename = "CMakeLists.txt"; then
    continue
  fi
  /bin/echo " \\"
  /bin/echo -n "$file"
done
echo
echo

echo "dist_noinst_SCRIPTS = \\"
echo "makefile.sh"
echo 

/bin/echo -n "EXTRA_DIST ="
for file in `ls COPYING.*`; do
  /bin/echo " \\"
  /bin/echo -n "$file"
done
echo


