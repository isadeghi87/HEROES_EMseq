

set -ex



test -f $PREFIX/lib/libCgl${SHLIB_EXT}
test -f $PREFIX/include/coin/CglAllDifferent.hpp
exit 0
