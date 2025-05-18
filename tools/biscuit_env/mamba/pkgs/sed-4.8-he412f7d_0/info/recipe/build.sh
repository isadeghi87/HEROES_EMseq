#!/bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* ./build-aux

set -ex

./configure --prefix="$PREFIX"

if [[ "$target_platform" == linux-* ]]; then
    # Remove test broken in sed 4.4
    sed -i.bak -e 's|testsuite/panic-tests.sh||' Makefile
fi

make $VERBOSE_AT

# These tests fail under emulation, still run them but ignore their result
if [[ ${target_platform} == linux-aarch64 ]]; then
    make check -j${NUM_CPUS} || true
elif [[ ${target_platform} == linux-ppc64le ]]; then
    make check -j${NUM_CPUS} || true
elif [[ "$CONDA_BUILD_CROSS_COMPILATION" != "1" ]]; then
    make check -j${NUM_CPUS} || (cat ./test-suite.log; false)
fi

if [[ -f ./test-suite.log ]]; then
    cat ./test-suite.log
fi

# remove installed sed from PATH so that it isn't called
export PATH=$(echo $PATH | sed "s@$PREFIX/bin:@@g")

make install
