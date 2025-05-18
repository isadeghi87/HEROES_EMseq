

set -ex



$R -e "library('xml2')"
$R -e "testthat::test_package('xml2')"
exit 0
