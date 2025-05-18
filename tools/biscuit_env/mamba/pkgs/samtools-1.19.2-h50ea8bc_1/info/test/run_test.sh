

set -ex



samtools --help
samtools view 'https://example.com' 2>&1 | grep 'fail to read the header' -q
exit 0
