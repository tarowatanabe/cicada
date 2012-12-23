#!/bin/sh

echo "eigen_includedir = \$(pkgincludedir)"
echo
echo "nobase_eigen_include_HEADERS = \\"

for file in `find Eigen -type f`; do
  basename=`basename $file`
  if test $basename = "CMakeLists.txt"; then
    continue
  fi
  echo "$file \\"
done
echo
