

set -ex



ldd $PREFIX/bin/bgzip | grep deflate
test -e test1.bed || printf 'chr1\t100\t200\n' > test1.bed
bgzip test1.bed
tabix test1.bed.gz
htsfile test1.bed.gz
exit 0
