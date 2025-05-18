#!/bin/bash

export DISABLE_AUTOBREW=1

# shellcheck disable=SC2086
${R} CMD INSTALL --build . ${R_ARGS} --configure-vars="INCLUDE_DIR=$PREFIX/include/freetype2 LIB_DIR=$PREFIX/lib"
