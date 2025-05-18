set -ex
export AS=nasm
# Build script need to be regenerated for cross compilation support
./autogen.sh
./configure --prefix=${PREFIX}
make -j${CPU_COUNT}
make -j${CPU_COUNT} install

rm -rf ${PREFIX}/share
rm ${PREFIX}/lib/libisal.a
