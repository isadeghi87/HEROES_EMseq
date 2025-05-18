#!/bin/bash -euo pipefail
printf "%s %s\n" AS-1052158-LR-69347_R1.fastq.gz OE0290-PED_I034-042_plasma-15-01_1.gz AS-1052158-LR-69347_R2.fastq.gz OE0290-PED_I034-042_plasma-15-01_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 6 OE0290-PED_I034-042_plasma-15-01_1.gz OE0290-PED_I034-042_plasma-15-01_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
