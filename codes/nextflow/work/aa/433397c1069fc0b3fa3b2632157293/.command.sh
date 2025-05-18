#!/bin/bash -euo pipefail
printf "%s %s\n" AS-1383045-LR-78290_R1.fastq.gz 2LB-304_plasma-01-01_1.gz AS-1383045-LR-78290_R2.fastq.gz 2LB-304_plasma-01-01_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 6 2LB-304_plasma-01-01_1.gz 2LB-304_plasma-01-01_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
