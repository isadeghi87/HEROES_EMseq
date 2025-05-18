

set -ex



$R -e "library('openssl')"
$R -e "testthat::test_dir('tests/testthat/', package='openssl', load_package='installed')"
exit 0
