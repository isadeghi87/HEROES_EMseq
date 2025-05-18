

set -ex



test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgfortran.so
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgfortran.a
test -f ${PREFIX}/bin/x86_64-conda_cos6-linux-gnu-gfortran
test -f ${PREFIX}/bin/x86_64-conda-linux-gnu-gfortran
find $PREFIX/lib -iname omp_lib.mod | grep '.'
find $PREFIX/lib -iname omp_lib.h | grep '.'
find $PREFIX/lib -iname ISO_Fortran_binding.h | grep '.'
echo 14
pushd tests/fortomp
sh test_fort.sh
exit 0
