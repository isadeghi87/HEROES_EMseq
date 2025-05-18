

set -ex



$R -e "library('ragg')"
$R -e "testthat::test_file('tests/testthat.R', stop_on_failure=TRUE)"
exit 0
