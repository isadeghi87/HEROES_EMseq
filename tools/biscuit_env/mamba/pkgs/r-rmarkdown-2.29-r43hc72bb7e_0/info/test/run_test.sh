

set -ex



$R -e "library('rmarkdown')"
$R -e "capabilities()" || true
$R -e "grSoftVersion()" || true
$R -e "x <- tempfile(fileext = '.Rmd'); file.create(x); rmarkdown::render(x)"
$R -e "tinytex::install_tinytex();testthat::test_dir('tests/testthat/', package='rmarkdown', load_package='installed')"
exit 0
