

set -ex



export GIT\_PYTHON_REFRESH=warn
snakemake --version
snakemake --version | grep "8.25.3"
exit 0
