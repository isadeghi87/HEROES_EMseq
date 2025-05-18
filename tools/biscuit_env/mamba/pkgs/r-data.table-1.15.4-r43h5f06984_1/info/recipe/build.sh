#!/bin/bash

## manual patch
## NB: This should be removed in the next release!
bzip2 -d inst/tests/tests.Rraw.bz2
patch --verbose -V none inst/tests/tests.Rraw < ${RECIPE_DIR}/patches/0001-no-warn-old-dates.patch
bzip2 inst/tests/tests.Rraw

export DISABLE_AUTOBREW=1

# shellcheck disable=SC2086
${R} CMD INSTALL --build . ${R_ARGS}
