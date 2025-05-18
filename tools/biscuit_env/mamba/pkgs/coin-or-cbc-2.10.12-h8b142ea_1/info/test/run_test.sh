

set -ex



test -f $PREFIX/lib/libCbc${SHLIB_EXT}
test -f $PREFIX/include/coin/CbcConfig.h
test -f $PREFIX/lib/libCbcSolver${SHLIB_EXT}
test -f $PREFIX/lib/libOsiCbc${SHLIB_EXT}
cbc -import test.lp -solve -solution '$' | tr -s ' ' | grep '1 y 1 0'
cbc test.lp solve solution '$' | tr -s ' ' | grep '1 y 1 0'
exit 0
