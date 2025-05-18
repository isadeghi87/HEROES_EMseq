

set -ex



bam2nuc --help 2>&1 | grep SYNOPSIS > /dev/null
bismark --help > /dev/null
bismark2bedGraph --help 2>&1 | grep SYNOPSIS > /dev/null
bismark2report --help 2>&1 | grep SYNOPSIS > /dev/null
bismark2summary --help 2>&1 | grep SYNOPSIS > /dev/null
bismark_genome_preparation --help > /dev/null
bismark_methylation_extractor --help > /dev/null
coverage2cytosine --help 2>&1 | grep SYNOPSIS > /dev/null
exit 0
