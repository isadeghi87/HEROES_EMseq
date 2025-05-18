#!/bin/bash -eux

${PYTHON} -m pip install . --no-deps --no-build-isolation -vv --global-option="--with-libyaml" --global-option="build_ext" --global-option="-I${PREFIX}/include" --global-option="-L${PREFIX}/lib"
