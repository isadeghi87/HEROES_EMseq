bsub -R "rusage[mem=200G]" -q verylong -n 25  -J EMseq  <nextflow_run.sh
#qsub  nextflow_run.sh