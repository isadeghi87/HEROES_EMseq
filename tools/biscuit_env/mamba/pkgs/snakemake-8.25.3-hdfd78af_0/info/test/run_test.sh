

set -ex



export GIT\_PYTHON_REFRESH=warn
snakemake --version
snakemake --version | grep "8.25.3"
printf %s\\n "rule empty:" > Snakefile && snakemake --cores 1 --report
rm -rf .snakemake Snakefile report.html
exit 0
