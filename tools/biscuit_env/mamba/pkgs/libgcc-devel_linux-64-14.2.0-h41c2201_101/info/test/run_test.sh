

set -ex



test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/crtbegin.o
test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/libgcc_eh.a
test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/libgcc.a
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgcc_s.so
exit 0
