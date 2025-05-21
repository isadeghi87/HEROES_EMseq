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




---

## Code Modules

### 1. `codes/nextflow`
- Implements the [nf-core/methylseq](https://github.com/nf-core/methylseq) pipeline using Nextflow.
- Produces alignment (Bismark), deduplication, and methylation call outputs.

### 2. `codes/dmr`
- Scripts for identifying and annotating differentially methylated regions across samples.

### 3. `codes/qc_summary`
- Aggregates QC metrics (e.g., coverage, duplication rate) across EMSeq samples into summary reports.

### 4. `codes/sample_prep`
- Generates input manifests and FASTQ file lists required by the Nextflow pipeline.

### 5. `codes/cnv_calling`
- Focuses on CNV detection from cfDNA data.
- **Use the `cfdna/` tool**—other scripts or tools in this directory are deprecated for liquid biopsy analyses.

---

## Data Files
- **Normalized methylation calls** are stored in `data/methylation_calls/`. These tables are used as input for downstream DMR analysis in the `codes/dmr/` module.
- **Sample sheets** in `datasets/` define sample metadata (IDs, paths, groups) for pipeline runs.

---

## Results

### Nextflow Pipeline Outputs (`results/nextflow`)
- **`bismark/deduplicated/`**: Deduplicated BAM files for CNV calling.
- **`bismark/methylation_calls/`**: Cytosine call files for DMR analysis.

### DMR Analysis (`results/dmr`)
- **Figures**: Plots and heatmaps summarizing differentially methylated regions.

### Biscuit Pipeline (`results/biscuit_pipeline`)
- Prototype scripts for running the [Biscuit](https://informatics.fas.harvard.edu/biscuit/) toolchain. Modify and re-run these to benchmark against Nextflow results.

---

## Usage
1. **Prepare samples**: Edit `datasets/samples_sheet.tsv` with sample metadata.  
2. **Run Nextflow**:
   ```bash
   cd codes/nextflow
   nextflow run nf-core/methylseq -profile odcf --input ../../datasets/samples_sheet.tsv```
3. **QC Summary**:

   ```bash
cd codes/qc_summary
Rscript summarize_qc.R --input ../../results/nextflow/bismark```

4. **DMR Analysis**: 
   ```bash
cd codes/dmr
python run_dmr.py --calls ../../data/methylation_calls/ --output ../../results/dmr/```

5. **CNV Calling**:

  ```bash
cd codes/cnv_calling/cfdna
bash run_cfdna.sh --bam ../../results/nextflow/bismark/deduplicated/*.bam```
