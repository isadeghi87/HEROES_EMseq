

set -ex



R -e "library('curl');data.frame(curl_options())"
R -e "testthat::test_file('tests/testthat.R', stop_on_failure=TRUE)"
R -e "invisible(curl::curl_fetch_memory('http://rest.kegg.jp/conv/ncbi-geneid/dme'))"
exit 0
