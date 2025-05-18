

set -ex



test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/libstdc++.a
test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/libstdc++fs.a
test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/libsupc++.a
test -f ${PREFIX}/lib/gcc/x86_64-conda-linux-gnu/14.2.0/include/c++/cstdio
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libstdc++.so
exit 0
