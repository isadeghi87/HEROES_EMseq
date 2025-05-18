#!/bin/bash

# tolerate older osx sdks
sed -ie 's/PKG_CPPFLAGS=/PKG_CPPFLAGS=-D_LIBCPP_DISABLE_AVAILABILITY /g' src/Makevars

export DISABLE_AUTOBREW=1

# shellcheck disable=SC2086
${R} CMD INSTALL --build . ${R_ARGS}
