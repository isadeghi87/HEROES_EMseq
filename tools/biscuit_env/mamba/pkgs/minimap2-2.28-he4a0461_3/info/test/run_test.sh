

set -ex



minimap2 2>&1 | grep 'Usage'
paftools.js version | grep 2.28
sdust  2>&1 | grep 'Usage'
test -e $PREFIX/lib/libminimap2.a
test -e $PREFIX/include/minimap.h
exit 0
