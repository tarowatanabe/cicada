#!/bin/sh

echo "eigen_includedir = \$(pkgincludedir)"
echo
echo "nobase_eigen_include_HEADERS = \\"

for file in `find Eigen -type f`; do
  echo "$file \\"
done
echo
