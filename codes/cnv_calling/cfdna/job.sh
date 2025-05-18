#!/bin/bash
bsub -R "rusage[mem=150G]" -n 10 -q long < run_cfdna.sh