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
