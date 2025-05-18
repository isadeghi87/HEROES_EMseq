

set -ex



bowtie2 --help
bowtie2 --version
bowtie2 --version | grep -q "bowtie2-align-s version [0-9]"
bowtie2-align-l --help
bowtie2-align-s --help
bowtie2-build --help
bowtie2-build-l --help
bowtie2-build-s --help
bowtie2-inspect --help
bowtie2-inspect-l --help
bowtie2-inspect-s --help
exit 0
