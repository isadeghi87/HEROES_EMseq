

set -ex



test ! -f "${PREFIX}/lib/libaccinj64.so"
test ! -f "${PREFIX}/lib/libcuinj64.so"
exit 0
