.libPaths("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/tools/R/4.0")
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/WGBS/')
# Main analysis package
library("methylKit")
# Annotation package
library("genomation")
library("GenomicRanges")

## dir 
library(data.table)
methyl= './results/nextflow/bismark/deduplicated/OE0290-PED_0LB-060_plasma-01-01_1_val_1_bismark_bt2_pe.deduplicated.sorted.chr1.bedgraph'
tumor='./results/nextflow/bismark/deduplicated/OE0290-PED_2LB-128_plasma-03-01_1_val_1_bismark_bt2_pe.deduplicated.sorted.chr1.bedgraph'
atac= './results/nextflow/bismark/deduplicated/pool3_atac.chr1.bedgraph'
# Load BEDGraph files
methylation_data <- fread(methyl, header = FALSE, col.names = c("chr", "start", "end", "methylation"))
tumor_data <- fread(tumor, header = FALSE, col.names = c("chr", "start", "end", "methylation"))
atac_seq_data <- fread(atac, header = FALSE, col.names = c("chr", "start", "end", "atac_signal"))

# Define Your Region of Interest
chromosome <- "chr1"
start_position <- 121500000  # Example start position
end_position <- 121900000    # Example end position

methylation_region <- methylation_data[chr == chromosome & start >= start_position & end <= end_position]
tumor_region <- tumor_data[chr == chromosome & start >= start_position & end <= end_position]
atac_seq_region <- atac_seq_data[chr == chromosome & start >= start_position & end <= end_position]

#Assuming methylation_region and atac_seq_region are already defined and filtered for your region of interest
# Aggregate data if necessary (e.g., calculate mean signal within smaller bins inside the region)
methylation_avg <- methylation_region[, .(avg_methylation = mean(methylation)), by = .(start, end)]
tumor_avg <- tumor_region[, .(avg_methylation = mean(methylation)), by = .(start, end)]
atac_seq_avg <- atac_seq_region[, .(avg_atac_signal = mean(atac_signal)), by = .(start, end)]

# Merge datasets by genomic coordinates
combined_data <- merge(methylation_region, atac_seq_region,by = c("start", "end"))
combined_tumor <- merge(tumor_region, atac_seq_region,by = c("start", "end"))


library(ggplot2)
ggplot(combined_data, aes(x = methylation, y = atac_signal)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Average Methylation", y = "Average ATAC-seq Signal", 
       title = "Correlation between Methylation in normal sample and ATAC-seq Signal") +
  geom_smooth(method = "lm", se = FALSE, color = "blue")  # Add a linear regression line


ggplot(combined_tumor, aes(x = methylation, y = atac_signal)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Average Methylation", y = "Average ATAC-seq Signal", 
       title = "Correlation between Methylation in tumor and ATAC-seq Signal") +
  geom_smooth(method = "lm", se = FALSE, color = "blue")  # Add a linear regression line


# Calculate Pearson correlation
correlation_result <- cor(methylation_region$methylation, atac_seq_region$atac_signal, method = "pearson")

# Print result
print(correlation_result)
