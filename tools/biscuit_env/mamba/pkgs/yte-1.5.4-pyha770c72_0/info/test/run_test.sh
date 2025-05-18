

set -ex



pip check
yte --help
echo -e '?if True:
  foo: 1' | yte
exit 0
