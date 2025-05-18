#!/bin/bash -euo pipefail
printf "%s %s\n" AS-1271236-LR-72262_R1.fastq.gz OE0290-PED_2LB-217_plasma-01-01_1.gz AS-1271236-LR-72262_R2.fastq.gz OE0290-PED_2LB-217_plasma-01-01_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 6 OE0290-PED_2LB-217_plasma-01-01_1.gz OE0290-PED_2LB-217_plasma-01-01_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
