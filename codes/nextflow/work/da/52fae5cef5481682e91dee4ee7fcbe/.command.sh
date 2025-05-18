#!/bin/bash -euo pipefail
bismark2summary 

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
