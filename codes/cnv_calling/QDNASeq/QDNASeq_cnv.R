setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("QDNAseq",force = T)
library(QDNAseq)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)

## parallel computing
future::plan("multisession", workers=6)

## bin size
binSize <- 50
 # create bins from the reference genome
bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=binSize)

# Path to your BAM file(s)
bamfiles <- c("results/nextflow/bismark/deduplicated/OE0290-PED_0LB-064_plasma-01-01_1_val_1_bismark_bt2_pe.deduplicated.bam.sorted.bam",
              "results/nextflow/bismark/deduplicated/OE0290-PED_I070-032_plasma-02-01_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam")

# Read count the BAM files using the custom bins
readCounts <- binReadCounts(bins, bamfiles = bamfiles)

# Apply quality filters to remove low-quality bins
readCountsFiltered <- applyFilters(readCounts)
# Estimate correction for the read counts
readCountsCorrected <- estimateCorrection(readCountsFiltered)

# Normalize the read counts to correct for biases
readCountsNormalized <- normalizeBins(readCountsFiltered)
print(class(readCountsFiltered))

# Segment the normalized bins
copyNumbersSegmented <- segmentBins(readCountsNormalized)

# Call CNVs based on segmented data
copyNumbersCalled <- callBins(copyNumbersSegmented)

# Plot CNVs
plot(copyNumbersCalled)

# Export CNV calls to a file
exportBins(copyNumbersCalled, file = "CNV_calls_hg38.txt")
