#!/bin/bash -euo pipefail
preseq \
    lc_extrap \
     -verbose -bam \
    -pe \
    -output OE0290-PED_5LB-003_plasma-02-01.lc_extrap.txt \
    OE0290-PED_5LB-003_plasma-02-01.sorted.bam
cp .command.err OE0290-PED_5LB-003_plasma-02-01.command.log

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP":
    preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*$//')
END_VERSIONS
