

set -ex



R -h
R --version
R -e "print(.Platform\$dynlib.ext)"
R -e "library('tcltk')"
Rscript --version
Rscript -e  'cat("ok\\n")'
R -e "capabilities()"
R -e "grSoftVersion()"
Rscript test-svg.R
Rscript -e "stopifnot(capabilities('jpeg'), TRUE)"
Rscript -e "stopifnot(capabilities('png'), TRUE)"
R -e "options(warn=2);Sys.time()"
exit 0
