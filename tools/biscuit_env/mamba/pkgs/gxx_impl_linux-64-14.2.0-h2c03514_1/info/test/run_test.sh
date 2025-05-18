

set -ex



test -f ${PREFIX}/bin/x86_64-conda_cos6-linux-gnu-g++
test -f ${PREFIX}/bin/x86_64-conda-linux-gnu-g++
CXX=$(${PREFIX}/bin/*-gcc -dumpmachine)-g++
${CXX} -Wall tests/aligned_alloc.cpp -c -o cpp_aligned.o --std=c++17
${CXX} -Wall cpp_aligned.o -o cpp_aligned --std=c++17 && ./cpp_aligned
${CXX} -Wall tests/hello_world.cpp -c -o hello_world.o --std=c++17
${CXX} -Wall hello_world.o -o hello_world --std=c++17
${CXX} -isystem ${PREFIX}/include --std=c++20 -Wall tests/tzdb-override.cpp -c -o tzdb-override.o
${CXX} -isystem ${PREFIX}/include --std=c++20 -Wall tests/tzdb.cpp -c -o tzdb.o
${CXX} tzdb-override.o -o tzdb-override && ./tzdb-override
${CXX} tzdb.o -o tzdb && ./tzdb
strace ./tzdb 2>&1 | grep ${CONDA_PREFIX}/lib/../share/zoneinfo/tzdata.zi
unset PREFIX
unset CONDA_PREFIX
./tzdb
exit 0
