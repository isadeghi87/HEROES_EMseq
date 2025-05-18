#bsub -R "rusage[mem=150G]" -q verylong -n 20 -J "biscuit_align" < biscuit_align.sh
bsub -R "rusage[mem=150G]" -q long -n 10 -J "biscuit" < biscuit.sh