#!/bin/bash
export DISABLE_AUTOBREW=1

# shellcheck disable=SC2086
if [ ${CONDA_BUILD_CROSS_COMPILATION} -eq 1 ]; then
    ${R} CMD INSTALL --build --no-demo --no-help --no-inst . ${R_ARGS}
else
    ${R} CMD INSTALL --build --install-tests . ${R_ARGS}
fi
