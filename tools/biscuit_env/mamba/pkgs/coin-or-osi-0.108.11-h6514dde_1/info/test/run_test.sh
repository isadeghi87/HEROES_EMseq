

set -ex



test -f $PREFIX/lib/libOsi${SHLIB_EXT}
test -f $PREFIX/include/coin/OsiConfig.h
exit 0
