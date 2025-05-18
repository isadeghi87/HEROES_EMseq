

set -ex



echo x86_64-conda_cos6-linux-gnu
echo sysroot
echo 2.12
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgomp.so
test `readlink ${PREFIX}/x86_64-conda-linux-gnu/lib/libgomp.so` == "../../lib/libgomp.so"
test -f ${PREFIX}/bin/x86_64-conda_cos6-linux-gnu-gcc
test -f ${PREFIX}/bin/x86_64-conda_cos6-linux-gnu-cpp
test -f ${PREFIX}/bin/x86_64-conda-linux-gnu-gcc
test -f ${PREFIX}/bin/x86_64-conda-linux-gnu-cpp
test ! -f ${PREFIX}/bin/gcc
test ! -f ${PREFIX}/bin/cpp
CC=$(${PREFIX}/bin/*-gcc -dumpmachine)-gcc
${CC} -Wall tests/aligned_alloc.c -c -o c_aligned.o -v -Wl,-v -march=native
${CC} -Wall tests/aligned_alloc.c -c -o c_aligned.o -v -fsanitize=address
${CC} -Wall tests/aligned_alloc.c -c -o c_aligned.o -v
${CC} -Wall c_aligned.o -o c_aligned -v && ./c_aligned
${CC} -Wall c_aligned.o -o c_aligned -Wl,-rpath,/foo && x86_64-conda-linux-gnu-readelf -d c_aligned | grep RPATH | grep "/foo:${PREFIX}/lib"
${CC} -Wall tests/hello_world.c -c -o hello_world.o -v
${CC} -Wall hello_world.o -o hello_world -v
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgcc_s.so
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libsanitizer.spec
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libasan_preinit.o
test -f ${PREFIX}/x86_64-conda-linux-gnu/lib/libgomp.spec
exit 0
