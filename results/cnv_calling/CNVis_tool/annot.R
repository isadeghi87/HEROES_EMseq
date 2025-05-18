library(GenomicRanges)   # For genomic data operations
library(rtracklayer)     # For importing genome annotations
library(ggplot2)         # For plotting
library(ggbio)           # For genomic visualizations

#1 Load a GTF file for hg38 annotation ####
gtf_file <- "../../../../pool_temp/datasets/references/human/gencode.v47.annotation.gtf.gz"  # Replace with the actual file path
genome_annotation <- rtracklayer::import(gtf_file)

# Filter for gene features
genome_annotation <- genome_annotation[genome_annotation$type == "gene"]


# Check the annotation
print(genome_annotation)

# Convert genomic bins to GRanges
cnv_granges <- GRanges(
  seqnames = gsub(":.*", "", cnv_data_filtered$genomic_bin),
  ranges = IRanges(
    start = as.numeric(gsub(".*:|-.*", "", cnv_data_filtered$genomic_bin)),
    end = as.numeric(gsub(".*-", "", cnv_data_filtered$genomic_bin))
  ),
  log2_ratio = cnv_data_filtered$log2_ratio_median,
  sample = cnv_data_filtered$sample
)
# Add metadata for visualization (e.g., color for gains and losses)
cnv_granges$cnv_type <- ifelse(cnv_granges$log2_ratio > 0.5, "Gain",
                               ifelse(cnv_granges$log2_ratio < -0.5, "Loss", "Neutral"))

# Check the GRanges object
print(cnv_granges)

# Plot the karyogram
p <- ggbio::autoplot(
  cnv_granges,
  aes(y = log2_ratio, fill = cnv_type),
  layout = "karyogram"
) +
  scale_fill_manual(
    values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
    name = "CNV Type"
  ) +
  theme_bw() +
  labs(
    title = "Genome-Wide CNV Karyogram with Gene Annotations",
    y = "Log2 Ratio",
    x = "Chromosome"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )

print(p)

# Filter for significant CNVs
significant_cnvs <- cnv_annotated[abs(cnv_annotated$log2_ratio) > 0.5, ]

# Add gene labels
p <- p + geom_text(
  data = as.data.frame(significant_cnvs),
  aes(
    x = (start + end) / 2,
    y = log2_ratio,
    label = gene
  ),
  size = 2,
  angle = 45,
  hjust = 0.5
)

print(p)

# Find overlaps between CNV regions and gene annotations
overlaps <- findOverlaps(cnv_granges, genome_annotation)
# Add gene names to CNV regions based on overlaps
cnv_granges$gene <- NA  # Initialize the gene column

# For overlapping regions, assign the gene name(s)
cnv_granges$gene[queryHits(overlaps)] <- genome_annotation$gene_name[subjectHits(overlaps)]

# Annotate CNVs with gene information
cnv_annotated <- as.data.frame(cnv_granges[queryHits(overlaps)])
cnv_annotated$gene <- genome_annotation$gene_name[subjectHits(overlaps)]  # Add gene names
cnv_annotated$gene <- genome_annotation$gene_name[subjectHits(overlaps)]  # Add gene names

# Check annotated CNVs
head(cnv_annotated)

library(data.table)

# Convert to data.table
cnv_annotated <- as.data.table(cnv_annotated)

# Summarize
cnv_summary <- cnv_annotated[, .(
  genes = paste(unique(gene), collapse = ","),
  log2_ratio = mean(log2_ratio),  # Or replace with another aggregation function
  sample = paste(unique(sample), collapse = ",")  # Combine all unique samples
), by = .(seqnames, start, end, width, strand)]



##Visualize Specific Genes Affected by CNVs: Plot the number of CNVs associated with each gene####
gene_counts <- cnv_annotated[, .N, by = gene]
gene_counts <- gene_counts[order(-N)]  # Sort by descending count

ggplot(gene_counts[1:20], aes(x = reorder(gene, -N), y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Top 20 Genes Affected by CNVs",
    x = "Gene",
    y = "Number of CNV Regions"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Check for Recurrent CNVs in Genes####
recurrent_genes <- cnv_annotated[, .(
  samples_with_cnv = length(unique(sample)),
  avg_log2_ratio = mean(log2_ratio, na.rm = TRUE)
), by = gene]

# Filter for recurrent genes (e.g., affected in at least 3 samples)
recurrent_genes <- recurrent_genes[samples_with_cnv >= 3]

print(recurrent_genes)

#Plot CNVs for a Specific Gene#
gene_of_interest <- "ADAM3A"

gene_cnvs <- cnv_annotated[gene == gene_of_interest]

ggplot(gene_cnvs, aes(x = sample, y = log2_ratio, fill = log2_ratio > 0)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = paste("CNV Profile for", gene_of_interest),
    x = "Sample",
    y = "Log2 Ratio"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Update the karyogram plot for better visibility
p <- ggbio::autoplot(
  cnv_granges_sig,
  aes(y = log2_ratio, fill = cnv_type),
  layout = "karyogram"
) +
  scale_fill_manual(
    values = c("Gain" = "red", "Loss" = "blue"),
    name = "CNV Type"
  ) +
  theme_bw() +
  labs(
    title = "Genome-Wide CNV Karyogram with Gene Annotations",
    y = "Log2 Ratio",
    x = "Chromosome"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  geom_text(
    data = as.data.frame(cnv_granges_sig[abs(cnv_granges_sig$log2_ratio) > 0.5, ]), # Filter significant CNVs
    aes(
      x = start + (end - start) / 2,  # Midpoint of the CNV
      y = log2_ratio + 0.5,           # Offset to avoid overlap
      label = gene                   # Gene names
    ),
    size = 3,                         # Adjust font size
    angle = 45,                       # Rotate text for better clarity
    hjust = 1,                        # Adjust horizontal justification
    vjust = 0.5,                      # Adjust vertical justification
    color = "black"                   # Label color
  )

print(p)

cnv_granges_sig = subset(cnv_granges,cnv_type != "Neutral")
# Update the karyogram plot for better visibility
p <- ggbio::autoplot(
  cnv_granges_sig,
  aes(y = log2_ratio, fill = cnv_type),
  layout = "karyogram"
) +
  scale_fill_manual(
    values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
    name = "CNV Type"
  ) +
  theme_bw() +
  labs(
    title = "Genome-Wide CNV Karyogram with Gene Annotations",
    y = "Log2 Ratio",
    x = "Chromosome"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  geom_text(
    data = as.data.frame(cnv_granges_sig[abs(cnv_granges_sig$log2_ratio) > 0.5, ]), # Filter significant CNVs
    aes(
      x = (start + end) / 2,  # Midpoint of the CNV
      y = log2_ratio,         # CNV log2 ratio
      label = gene            # Gene names
    ),
    size = 3,                 # Increase font size
    angle = 0,               # Rotate text for better visibility
    #hjust = 0.5,
    vjust = -1,
    color = "black"
  )
p

samples = unique(cnv_annotated$sample)

for (sample in samples) {
  sample_data = subset(cnv_annotated, sample == sample)
  
  # Add a CNV type classification column
  sample_data$cnv_type <- ifelse(
    sample_data$log2_ratio > 0.5, "Gain",
    ifelse(sample_data$log2_ratio < -0.5, "Loss", "Neutral")
  )
  
  # Plot CNV profile with genes and CNV type-based coloring
  p <- ggplot(sample_data, aes(x = (start + end) / 2, y = log2_ratio, color = cnv_type)) +
    geom_point(size = 1) +
    facet_wrap(~seqnames, scales = "free_x", ncol = 6) +
    theme_minimal() +
    labs(
      title = paste("CNV Profile for Sample", sample),
      x = "Genomic Position",
      y = "Normalized Log2 Ratio"
    ) +
    scale_color_manual(
      values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
      name = "CNV Type"
    ) +
    geom_text(
      data = subset(sample_data, abs(log2_ratio) > 0.5), # Filter significant CNVs
      aes(label = gene),
      size = 2,
      angle = 0,
      hjust = 0.5,
      vjust = -1,
      check_overlap = TRUE
    ) +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "right" # Enable legend for CNV types
    )
  print(p)
}
