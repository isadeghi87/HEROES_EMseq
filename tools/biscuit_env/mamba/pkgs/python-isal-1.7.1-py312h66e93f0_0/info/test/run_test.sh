

set -ex



python -c "from isal import igzip"
python -c "from isal import isal_zlib; isal_zlib.adler32(b'test')"
exit 0
