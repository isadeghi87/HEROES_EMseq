#!/bin/bash -euo pipefail
unset DISPLAY
mkdir -p tmp
export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
qualimap \
    --java-mem-size=29491M \
    bamqc \
    -gd HUMAN \
    -bam OE0290-PED_I034-042_plasma-15-01.sorted.bam \
     \
    -p non-strand-specific \
    --collect-overlap-pairs \
    -outdir OE0290-PED_I034-042_plasma-15-01 \
    -nt 6

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC":
    qualimap: $(echo $(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*$//')
END_VERSIONS
