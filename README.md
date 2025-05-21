**EMSeq Data Analysis for Liquid Biopsy in Pediatric Tumors**

This repository contains code and results for analysis of EMSeq (Enzymatic Methyl-seq) data from liquid biopsy samples in pediatric tumors.

## Table of Contents

- [Project Overview](#project-overview)
- [Directory Structure](#directory-structure)
- [Code Modules](#code-modules)
- [Data Files](#data-files)
- [Results](#results)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Project Overview

The goal of this project is to perform comprehensive DNA methylation analysis using EMSeq data derived from liquid biopsy samples of pediatric tumors. Analyses include read alignment, quality control, differential methylation region (DMR) detection, and copy number variation (CNV) calling.

**Main directories:**

- `/omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/`
- `/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/`

## Directory Structure


├── codes/
│   ├── nextflow/           # nf-core/methylseq pipeline
│   ├── dmr/                # DMR analysis scripts
│   ├── qc_summary/         # QC summary scripts
│   ├── sample_prep/        # Input preparation for pipelines
│   └── cnv_calling/        # CNV calling (use cfdna/ only)
│       └── cfdna/
├── data/
│   └── methylation_calls/  # Normalized calls for DMR
├── datasets/
│   └── samples_sheet.tsv   # Sample metadata
├── results/
│   ├── nextflow/
│   │   ├── bismark/
│   │   │   ├── deduplicated/
│   │   │   └── methylation_calls/
│   ├── dmr/
│   │   └── figures/
│   └── biscuit_pipeline/
└── README.md


## Code Modules

### `codes/nextflow`
- nf-core/methylseq via Nextflow: alignment, deduplication, methylation calls.

### `codes/dmr`
- Identify & annotate differentially methylated regions.

### `codes/qc_summary`
- Aggregate QC metrics (coverage, duplication rate, etc.).

### `codes/sample_prep`
- Generate sample manifests & FASTQ lists.

### `codes/cnv_calling`
- CNV detection on cfDNA (use **cfdna/** only).

## Data Files

- `data/methylation_calls/`: normalized methylation-call tables.  
- `datasets/samples_sheet.tsv`: sample metadata for pipeline.

## Results

### `results/nextflow`
- `bismark/deduplicated/`: BAMs for CNV calling.  
- `bismark/methylation_calls/`: cytosine-call files for DMR.

### `results/dmr`
- `figures/`: DMR plots & heatmaps.

### `results/biscuit_pipeline`
- Prototype scripts for Biscuit benchmarking.

## Usage

1. Prepare samples: edit `datasets/samples_sheet.tsv`.  
2. Run Nextflow:
   ```bash
   cd codes/nextflow
   nextflow run nf-core/methylseq -profile odcf --input ../../datasets/samples_sheet.tsv
