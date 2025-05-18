

set -ex



test -f $PREFIX/lib/libClp${SHLIB_EXT}
test -f $PREFIX/lib/libOsiClp${SHLIB_EXT}
test -f $PREFIX/include/coin/ClpGubDynamicMatrix.hpp
test -f $PREFIX/include/coin/OsiClpSolverInterface.hpp
echo ? | clp
exit 0
