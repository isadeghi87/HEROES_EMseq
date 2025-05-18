

set -ex



test -f $PREFIX/lib/libCoinUtils${SHLIB_EXT}
test -f $PREFIX/include/coin/CoinSort.hpp
exit 0
