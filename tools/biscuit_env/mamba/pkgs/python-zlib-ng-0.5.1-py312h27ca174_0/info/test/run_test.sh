

set -ex



python -c "from zlib_ng import gzip_ng"
python -c "from zlib_ng import zlib_ng; zlib_ng.adler32(b'test')"
exit 0
