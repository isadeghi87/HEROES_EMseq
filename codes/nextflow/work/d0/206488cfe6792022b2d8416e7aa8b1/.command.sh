#!/bin/bash -euo pipefail
bismark2report 

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_REPORT":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
