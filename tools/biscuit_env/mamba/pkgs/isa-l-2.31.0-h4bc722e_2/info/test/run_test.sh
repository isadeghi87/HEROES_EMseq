

set -ex



igzip --help
igzip -k LICENSE
bash -c '[ $(stat -c %s LICENSE) -gt $(stat -c %s LICENSE.gz) ]'
bash -c '[ "$(igzip -cd LICENSE.gz)" == "$(cat LICENSE)" ]'
exit 0
