

set -ex



$R -e "library('data.table')"
$R -e "library('data.table');test.data.table()"
exit 0
