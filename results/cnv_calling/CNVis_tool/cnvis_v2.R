setwd("~/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/results/cnv_calling/cfdna/")

# Load necessary libraries
library(data.table)      # For efficient data handling
library(ggplot2)         # For plotting
library(ComplexHeatmap)  # For advanced heatmaps
library(circlize)        # For color mapping in heatmaps
library(dendextend)      # For enhancing dendrograms
library(RColorBrewer)    # For color palettes
library(reshape2)        # For data reshaping
library(stats)           # For statistical functions
library(ggbio)
library(GenomicRanges)

# Custom functions
# Step 1: Read segmentation files ####
# Use the list of your actual segs files
# (Assuming seg_files is already defined with the list of file paths)
seg_files <- list.files(path = "./", pattern = "segs$", full.names = TRUE)
#seg_files = seg_files[grep('sorted',seg_files,invert = T)]

## functions ####
# Function to normalize CNV data
normalize_cnv <- function(cnv_data) {
  # Median centering normalization for each sample
  cnv_data$log2_ratio_normalized <- cnv_data$log2_ratio_median - median(cnv_data$log2_ratio_median, na.rm = TRUE)
  return(cnv_data)
}


# Function to create a CNV matrix for heatmap
create_cnv_matrix <- function(cnv_data) {
  cnv_matrix <- dcast(
    cnv_data, 
    genomic_bin ~ sample, 
    value.var = "log2_ratio_median", 
    fun.aggregate = mean
  )
  
  # Set genomic_bin as row names and remove the column
  rownames(cnv_matrix) <- cnv_matrix$genomic_bin
  cnv_matrix$genomic_bin <- NULL
  
  # Convert to matrix and handle missing values
  cnv_matrix <- as.matrix(cnv_matrix)
  cnv_matrix[is.na(cnv_matrix)] <- 0
  return(cnv_matrix)
}


# Function to perform hierarchical clustering
perform_clustering <- function(cnv_matrix) {
  distance_matrix <- dist(t(cnv_matrix), method = "euclidean")
  hc <- hclust(distance_matrix, method = "ward.D2")
  return(hc)
}

# Function to generate Manhattan plot
generate_manhattan_plot <- function(cnv_data) {
  cnv_data$chrom <- factor(cnv_data$chrom, levels = unique(cnv_data$chrom))
  cnv_data$position <- (cnv_data$start + cnv_data$end) / 2
  cnv_data$log10_pvalue <- -log10(cnv_data$p_value)
  
  ggplot(cnv_data, aes(x = position, y = log10_pvalue, color = chrom)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = rep(brewer.pal(12, "Paired"), length.out = length(unique(cnv_data$chrom)))) +
    theme_minimal() +
    labs(title = "Manhattan Plot of CNV Data", x = "Genomic Position", y = "-Log10(p-value)") +
    theme(legend.position = "none")
}

# Main script starts here

# # Function to extract sample ID from filename
# extract_sample_id <- function(filename) {
#   # Remove the directory path
#   base_name <- basename(filename)
#   # Remove the ".segs" extension
#   base_name <- sub("\\.segs$", "", base_name)
#   # Extract sample ID using a regex pattern
#   # For your filenames, extract everything up to the first "_plasma"
#   #sample_id <- sub("_plasma.*", "", base_name)
#   return(sample_id)
# }

# Initialize an empty data frame to hold all CNV data
cnv_data_all <- data.frame()

# Loop through each file and read the data
for (file in seg_files) {
  cnv_data <- fread(file)
  
  # Ensure the data has the required columns
  required_cols <- c("sample", "chrom", "start", "end", "log2_ratio_median")
  if (!all(required_cols %in% colnames(cnv_data))) {
    stop(paste("File", file, "is missing required columns."))
  }
  
  # Extract sample ID from the filename
  #sample_id <- extract_sample_id(file)
  #cnv_data$sample_id <- sample_id
  cnv_data$sample = gsub('_1_val_1_bismark_bt2_pe.deduplicated','',cnv_data$sample)
  
  # Append to the main data frame
  cnv_data_all <- rbind(cnv_data_all, cnv_data)
}

# Examine the distribution of raw log2_ratio_median values
hist(cnv_data_all$log2_ratio_median, breaks = 100, main = "Distribution of Raw Log2 Ratios", xlab = "Log2 Ratio")

# Step 2: Simulate annotations data frame for testing ####

# Extract unique sample IDs from the CNV data
cnv_data_all$sample = gsub('OE0290-PED_','',cnv_data_all$sample)
unique_sample_ids <- unique(cnv_data_all$sample)

# Check the number of unique samples
cat("Number of unique samples:", length(unique_sample_ids), "\n")


# Define possible categories
annotations = read.delim("~/Documents/emseq_annot.tsv")


# View the annotations data frame
print(annotations)
id = intersect(cnv_data_all$sample,annotations$sample)

cnv_data_all = cnv_data_all[cnv_data_all$sample %in% id, ]
annotations = annotations[annotations$sample %in% id,]

# Step 3: Preprocess and normalize data ####
# Create a unique identifier for genomic bins
cnv_data_all$genomic_bin <- paste0(cnv_data_all$chrom, ":", cnv_data_all$start, "-", cnv_data_all$end)
cnv_granges <- GRanges(
  seqnames = cnv_data_all$chrom,
  ranges = IRanges(start = cnv_data_all$start, end = cnv_data_all$end),
  log2_ratio_median = cnv_data_all$log2_ratio_median,
  sample = cnv_data_all$sample
)


tumor_samples = annotations$sample[annotations$diagnosis!='Control']

tumor_granges <- cnv_granges[cnv_granges$sample %in% tumor_samples]

cnv_filtered = cnv_data_all[cnv_data_all$sample %in%  tumor_samples]
cnv_filtered$genomic_bin = paste0(cnv_filtered$chrom, ":", cnv_filtered$start, "-", cnv_filtered$end)

# Normalize data for each sample
#cnv_data_normalized <- cnv_filtered[, normalize_cnv(.SD), by = sample]

# Ensure that only numeric columns are used
numeric_cols <- sapply(cnv_filtered, is.numeric)

# Check the range of log2_ratio_normalized
range(cnv_filtered$log2_ratio_median, na.rm = TRUE)

# Step 4: Merge CNV data with annotations
# Merge annotations with CNV data based on sample_id
cnv_merged <- merge(cnv_filtered, annotations, by = "sample", all.x = TRUE)

# Step 5: Create CNV matrix for visualization
main_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
cnv_merged <- cnv_merged[
  chrom %in% main_chromosomes & !is.na(chrom)
]

threshold <- 0.1  # Set your threshold value
cnv_merged_filtered <- cnv_merged[abs(log2_ratio_median) > threshold, ]
cnv_merged_filtered$chrom <- factor(cnv_merged_filtered$chrom, levels = unique(cnv_merged_filtered$chrom))

# Step 6: Perform clustering ####
cnv_matrix <- create_cnv_matrix(cnv_merged_filtered)
hc_samples <- perform_clustering(cnv_matrix)

# Step 7: Generate Heatmap ####
# Prepare annotation for heatmap columns (samples)
# Create a data frame with sample annotations
sample_annotations <- annotations[match(colnames(cnv_matrix), annotations$sample), ]


# Define colors for annotations, ensuring that each is a named vector
# Function to generate colors and handle cases where n < 3
generate_colors <- function(levels, palette_name) {
  n_levels <- length(levels)
  n_colors <- max(n_levels, 3)  # Ensure minimum n of 3 for brewer.pal
  colors <- brewer.pal(n = n_colors, name = palette_name)[1:n_levels]
  names(colors) <- levels
  return(colors)
}

# Get unique levels for each annotation
levels_diagnosis <- unique(sample_annotations$diagnosis)
levels_status <- unique(sample_annotations$status)
levels_entity <- unique(sample_annotations$entity)

# Generate colors for each annotation
colors_diagnosis <- generate_colors(levels_diagnosis, "Set3")
colors_status <- generate_colors(levels_status, "Set2")
colors_entity <- generate_colors(levels_entity, "Set1")

# Create the HeatmapAnnotation object
ha <- rowAnnotation(
  df = sample_annotations[, c("diagnosis", "status", "entity")],
  col = list(
    diagnosis = colors_diagnosis,
    status = colors_status,
    entity = colors_entity
  )
)

# Define colors for heatmap
heatmap_colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))



# Plot heatmap with annotations
Heatmap(
  t(cnv_matrix),
  name = "Log2 Ratio",
  col = heatmap_colors,
  right_annotation = ha,
  border = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_title = "Samples",
  column_title = "Genomic Bins",
  column_names_gp = gpar(fontsize = 2), # Smaller font
  column_names_rot = 90, # Rotate labels for readability
  heatmap_legend_param = list(title = "Log2 Ratio")
)

# Step 8: Generate Manhattan Plot
# Assuming p-values are available or need to be calculated
# For demonstration, simulate p-values
set.seed(123)
cnv_merged_filtered$p_value <- runif(nrow(cnv_merged_filtered), min = 0, max = 1)

print(manhattan_plot)

# Step 9: Analyze heterogeneity between samples
# Perform Principal Component Analysis (PCA)
pca_result <- prcomp(t(cnv_matrix))

# Merge PCA results with annotations
pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)
pca_data <- merge(pca_data, annotations, by.x = "Sample", by.y = "sample", all.x = TRUE)

# Plot PCA results colored by diagnosis and shaped by status
ggplot(pca_data, aes(x = PC1, y = PC2, color = entity, label = Sample)) +
  geom_point(size = 3) +  # Add the points
  geom_text(vjust = -1, size = 3) +  # Add sample names just above the points
  theme_minimal() +
  labs(title = "PCA of CNV Data", x = "PC1", y = "PC2") +
  scale_color_manual(values = brewer.pal(n = length(unique(pca_data$entity)), name = "Set1")) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Adjust shapes as needed
  theme(legend.position = "right")

# Step 10: Statistical analysis for heterogeneity
# Calculate pairwise distances between samples
distance_matrix <- as.matrix(dist(t(cnv_matrix), method = "euclidean"))

# Perform hierarchical clustering and plot dendrogram with annotations
hc <- hclust(as.dist(distance_matrix), method = "ward.D2")
dend <- as.dendrogram(hc)

# Color branches by entity
labels_colors <- annotations$entity[match(labels(dend), annotations$sample)]
labels_colors_factor <- as.numeric(factor(labels_colors))
labels_colors_palette <- brewer.pal(n = max(labels_colors_factor), name = "Set2")
labels_colors_assigned <- labels_colors_palette[labels_colors_factor]

dend <- dendextend::set(dend, "labels_col", labels_colors_assigned)
dend <- dendextend::set(dend, "branches_k_color", labels_colors_assigned)

# Plot dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram of Samples")

# Add legend for entity colors
legend("topright", legend = unique(labels_colors), fill = unique(labels_colors_assigned), border = NA, bty = "n")

# Step 11: Identify significant CNV regions across samples
# Aggregate CNV data by genomic bins and calculate mean log2 ratios
cnv_summary <- cnv_merged_filtered[, .(
  mean_log2_ratio = mean(log2_ratio_median, na.rm = TRUE),
  sd_log2_ratio = sd(log2_ratio_median, na.rm = TRUE)
), by = .(chrom, start, end, genomic_bin)]

# Identify regions with high variance (indicative of heterogeneity)
high_variance_threshold <- quantile(cnv_summary$sd_log2_ratio, 0.9, na.rm = TRUE)
high_variance_regions <- cnv_summary[sd_log2_ratio > high_variance_threshold]

# Step 12: Plot high variance regions
ggplot(high_variance_regions, aes(x = genomic_bin, y = sd_log2_ratio)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = "High Variance CNV Regions", x = "Genomic Bin", y = "Standard Deviation of Log2 Ratio") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save the high variance regions to a file
# fwrite(high_variance_regions, "high_variance_cnv_regions.csv")

# New Steps: Additional Plots Showing Loss or Gain Across Chromosomes

# Step 13: Aggregate CNV data per chromosome per sample
# Calculate the average log2_ratio_normalized per chromosome per sample
# Aggregate using median instead of mean
cnv_chromosome <- cnv_merged_filtered[, .(
  median_log2_ratio = median(log2_ratio_median, na.rm = TRUE)
), by = .(sample, chrom)]

main_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
cnv_chromosome_filtered <- cnv_chromosome[
  chrom %in% main_chromosomes & !is.na(chrom)
]

# Create a matrix with chromosomes as rows and samples as columns
cnv_chrom_matrix <- dcast(cnv_chromosome_filtered, chrom ~ sample, value.var = "median_log2_ratio",fun.aggregate = mean)
rownames(cnv_chrom_matrix) <- cnv_chrom_matrix$chrom
cnv_chrom_matrix$chrom <- NULL
cnv_chrom_matrix <- as.matrix(cnv_chrom_matrix)
# Handle missing values (e.g., if some chromosomes are missing in some samples)
cnv_chrom_matrix[is.na(cnv_chrom_matrix)] <- 0

# Step 14: Perform clustering on chromosome-level CNV data
hc_samples_chrom <- perform_clustering(cnv_chrom_matrix)

# Step 15: Create heatmap of chromosome-level CNV data
# Define colors for gains and losses
chrom_heatmap_colors <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

rownames = rownames(cnv_chrom_matrix)
rownames = rownames[order(as.numeric(sub("chr", "", rownames)))]
cnv_chrom_matrix = cnv_chrom_matrix[rownames,]

# Plot heatmap with annotations
Heatmap(
  t(cnv_chrom_matrix),
  name = "Avg Log2 Ratio",
  col = chrom_heatmap_colors,
  right_annotation = ha,
  cluster_rows = F,cluster_columns = F,
  #cluster_columns = as.dendrogram(hc_samples_chrom),
  show_row_names = TRUE,
  show_column_names = T,
  column_title = "Samples",
  row_title = "Chromosomes",
  heatmap_legend_param = list(title = "Avg Log2 Ratio")
)

# Step 16: Generate chromosome-wise CNV profiles for each sample
# Optionally, create plots for individual samples or groups

# For demonstration, plot CNV profiles for a subset of samples
# Adjust the number of samples as needed
samples_to_plot <- unique(cnv_merged_filtered$sample)  # Adjust as needed

for (sample in samples_to_plot) {
  sample_data <- cnv_merged_filtered[sample == sample]
  
  
  # Adjust chromosome bin as a factor to maintain order
  sample_data$chrom <- factor(sample_data$chrom, levels = unique(sample_data$chrom))
  
  # Plot the CNV profile for the specific sample
  samp_plot <- ggplot(sample_data, aes(x = (start + end) / 2, y = chrom, fill = log2_ratio_median)) +
    geom_tile() +  # Use tiles for better clarity and spacing
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = "Log2 Ratio"
    ) +
    theme_minimal() +
    labs(
      title = paste("Genome-wide CNV Profile for Sample", "OE0290-PED_0LB-066"),
      x = "Genomic Position (Mb)",
      y = "Chromosomes"
    ) +
    theme(
      axis.text.y = element_text(size = 8),  # Adjust font size for better readability
      legend.position = "right",  # Move legend to the side
      panel.grid.major = element_blank(),  # Remove gridlines
      panel.grid.minor = element_blank()
    )
  
  print(samp_plot)
}

# Step 17: Karyogram-like plots to show gains and losses across chromosomes
# Using the ggbio package to create idiograms

# Install and load ggbio if not already installed
if (!requireNamespace("ggbio", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ggbio")
}



# Prepare data for ggbio plotting
# Create a GRanges object from the CNV data
cnv_granges <- GRanges(
  seqnames = cnv_merged_filtered$chrom,
  ranges = IRanges(start = cnv_merged_filtered$start, end = cnv_merged_filtered$end),
  log2_ratio = cnv_merged_filtered$log2_ratio_median,
  sample = cnv_merged_filtered$sample
)

# Customize the karyogram plot
for (sample in samples_to_plot) {
  sample_granges <- cnv_granges[cnv_granges$sample == sample]
  
  # Adjust colors and themes for the karyogram
  p <- ggbio::autoplot(
    sample_granges,
    aes(fill = log2_ratio),  # Use fill for log2 ratio
    layout = "karyogram"
  ) +
    scale_fill_gradient2(
      low = "blue",  # Deletions
      mid = "white", # Neutral regions
      high = "red",  # Amplifications
      midpoint = 0,  # Center at 0 log2 ratio
      name = "Log2 Ratio"  # Legend title
    ) +
    theme_bw() +
    labs(
      title = paste("Genome-wide CNV Profile for Sample", sample)
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # Print the plot for each sample
  print(p)
}


# Step 18: Cluster samples based on chromosome-level CNV profiles
# The clustering was performed in Step 14, and the heatmap in Step 15 shows clustered samples.

# Optionally, you can perform hierarchical clustering and plot a dendrogram
distance_matrix_chrom <- as.matrix(dist(t(cnv_chrom_matrix), method = "euclidean"))
hc_chrom <- hclust(as.dist(distance_matrix_chrom), method = "ward.D2")
dend_chrom <- as.dendrogram(hc_chrom)

# Color branches by diagnosis
labels_colors_chrom <- annotations$diagnosis[match(labels(dend_chrom), annotations$sample)]
labels_colors_factor_chrom <- as.numeric(factor(labels_colors_chrom))
labels_colors_palette_chrom <- brewer.pal(n = max(labels_colors_factor_chrom), name = "Set2")
labels_colors_assigned_chrom <- labels_colors_palette_chrom[labels_colors_factor_chrom]

dend_chrom <- dendextend::set(dend_chrom, "labels_col", labels_colors_assigned_chrom)
dend_chrom <- dendextend::set(dend_chrom, "branches_k_color", labels_colors_assigned_chrom)

# Plot dendrogram
plot(dend_chrom, main = "Hierarchical Clustering Dendrogram (Chromosome-level CNV)")

# Add legend for diagnosis colors
legend("topright", legend = unique(labels_colors_chrom), fill = unique(labels_colors_assigned_chrom), border = NA, bty = "n")

hc <- hclust(dist(cnv_matrix_filtered), method = "ward.D2")
plot(hc, main = "Sample Clustering Based on CNV Patterns")

high_variance_regions <- cnv_merged_filtered[, .(
  mean_log2_ratio = mean(log2_ratio_median, na.rm = TRUE),
  variance = var(log2_ratio_median, na.rm = TRUE)
), by = genomic_bin]

ggplot(high_variance_regions, aes(x = genomic_bin, y = variance)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Variance of CNVs Across Genomic Regions")


##mutation burden ####
# Compute CNV burden for each chromosome
cnv_burden <- cnv_merged_filtered[, .(
  cnv_gain = sum(log2_ratio_median > 0.25, na.rm = TRUE),
  cnv_loss = sum(log2_ratio_median < -0.25, na.rm = TRUE)
), by = chrom]

# Reshape data for plotting gains and losses side by side
cnv_burden_melted <- melt(cnv_burden, id.vars = "chrom", variable.name = "type", value.name = "burden")

# Plot gains and losses
ggplot(cnv_burden_melted, aes(x = chrom, y = burden, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("cnv_gain" = "red", "cnv_loss" = "blue"),
    labels = c("CNV Gain", "CNV Loss")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    title = "CNV Burden Across Chromosomes",
    x = "Chromosome",
    y = "Mutation Burden"
  ) +
  geom_text(
    aes(label = burden),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  )

#. Comparative CNV Profiles Between Groups
# Example: t-test for each region
# Ensure diagnosis is a factor with two levels
cnv_merged_filtered$diagnosis <- factor(cnv_merged_filtered$diagnosis)

# Perform t-tests for each genomic bin
differential_cnv <- cnv_merged_filtered[, {
  if (length(unique(diagnosis[!is.na(diagnosis)])) == 2) {  # Ensure two levels in diagnosis
    test_result <- tryCatch(
      t.test(log2_ratio_median ~ diagnosis, na.action = na.omit),  # Perform t-test
      error = function(e) NULL  # Return NULL if t-test fails
    )
    if (!is.null(test_result)) {
      .(p_value = test_result$p.value, 
        mean_difference = diff(tapply(log2_ratio_median, diagnosis, mean, na.rm = TRUE)))
    } else {
      .(p_value = NA_real_, mean_difference = NA_real_)  # Explicitly return numeric NA
    }
  } else {
    .(p_value = NA_real_, mean_difference = NA_real_)  # Explicitly return numeric NA
  }
}, by = genomic_bin]

# Adjust p-values for multiple testing
differential_cnv[, adjusted_p_value := p.adjust(p_value, method = "BH")]

# Filter significant results
significant_cnv <- differential_cnv[!is.na(adjusted_p_value) & adjusted_p_value < 0.05]

# Inspect significant CNV regions
print(significant_cnv)



## umap
library(uwot)
umap_result <- umap(t(cnv_matrix), n_neighbors = 15)
# Define a color palette for diagnosis
diagnosis_colors <- c("Tumor" = "red", "Control" = "blue")

# Map the diagnosis labels to colors
diagnosis_colored <- diagnosis_colors[annotations$diagnosis]

# Ensure the `annotations` and `umap_result` rows are in the same order
if (!all(rownames(umap_result) == annotations$sample_id)) {
  stop("Mismatch between UMAP results and annotations")
}

# Assuming 'umap_result' contains the UMAP coordinates and 'annotations' has sample IDs
umap_data <- data.frame(
  Sample = annotations$sample, # Add sample names
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Diagnosis = annotations$diagnosis,
  Status = annotations$status
)

# Plot UMAP with sample labels
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Diagnosis, shape = Status)) +
  geom_point(size = 3, alpha = 0.8) +  # Add points
  geom_text(aes(label = Sample), hjust = -0.2, vjust = -0.2, size = 3) +  # Add labels
  theme_minimal() +
  labs(
    title = "UMAP of CNV Data",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors for diagnosis
  scale_shape_manual(values = c(16, 17, 18)) +  # Customize shapes for status
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )



## Time Series or Longitudinal Analysis ####
# Function to extract timepoint from sample names
extract_timepoint <- function(sample_name) {
  if (grepl("control", sample_name, ignore.case = TRUE)) {
    return(0)  # Assign 0 to control samples
  } else if (grepl("plasma-([0-9]+)-", sample_name)) {
    timepoint <- sub(".*plasma-([0-9]+)-.*", "\\1", sample_name)
    return(as.numeric(timepoint))
  } else {
    return(NA)  # Handle unexpected formats
  }
}

# Apply the function to extract timepoints
cnv_merged_filtered$timepoint <- sapply(cnv_merged_filtered$sample, extract_timepoint)

# View the updated annotations data frame
head(cnv_merged_filtered)


#For samples with repeated measures (e.g., relapse vs. primary), analyze CNV progression over time.
ggplot(cnv_merged_filtered, aes(timepoint, log2_ratio_median, color = genomic_bin)) +
  geom_line(show.legend = F) +
  theme_minimal() +
  labs(title = "Temporal Dynamics of CNVs")


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
  seqnames = gsub(":.*", "", cnv_merged_filtered$genomic_bin),
  ranges = IRanges(
    start = as.numeric(gsub(".*:|-.*", "", cnv_merged_filtered$genomic_bin)),
    end = as.numeric(gsub(".*-", "", cnv_merged_filtered$genomic_bin))
  ),
  log2_ratio = cnv_merged_filtered$log2_ratio_median,
  sample = cnv_merged_filtered$sample
)

# Add metadata for visualization (e.g., color for gains and losses)
cnv_granges$cnv_type <- ifelse(cnv_granges$log2_ratio > 0.1, "Gain",
                               ifelse(cnv_granges$log2_ratio < -0.1, "Loss", "Neutral"))

# Check the GRanges object
print(cnv_granges)


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

# Convert to data.table
cnv_annotated <- as.data.table(cnv_annotated)
write.table(cnv_annotated,"../CNVar/EMseq_CNV_annotated.tsv",quote = F,sep = '\t',row.names = F)

subset(cnv_annotated,seqnames=='chr7'& sample=='')
id = !duplicated(cnv_annotated[,1:n])
cnv_uniq = cnv_annotated[id,]

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


# Filter for significant CNVs
significant_cnvs <- cnv_uniq[abs(cnv_uniq$log2_ratio) > 0.1, ]

# Add gene labels
p <- p + geom_text(
  data = as.data.frame(significant_cnvs),
  aes(
    x = (start + end) / 2,
    y = log2_ratio,
    label = gene
  ),
  size = 2,
  angle = 0,
  vjust = -1,
  hjust = 0.5
)

print(p)

# Summarize
cnv_summary <- cnv_uniq[, .(
  genes = paste(unique(gene), collapse = ","),
  log2_ratio = mean(log2_ratio),  # Or replace with another aggregation function
  sample = paste(unique(sample), collapse = ",")  # Combine all unique samples
), by = .(seqnames, start, end, width, strand)]



##Visualize Specific Genes Affected by CNVs: Plot the number of CNVs associated with each gene####
gene_counts <- cnv_uniq[, .N, by = gene]
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
recurrent_genes <- cnv_uniq[, .(
  samples_with_cnv = length(unique(sample)),
  avg_log2_ratio = mean(log2_ratio, na.rm = TRUE)
), by = gene]

# Filter for recurrent genes (e.g., affected in at least 3 samples)
recurrent_genes <- recurrent_genes[samples_with_cnv >= 3]

print(recurrent_genes)

#Plot CNVs for a Specific Gene#
gene_of_interest <- "BBS9"

gene_cnvs <- cnv_uniq[gene == gene_of_interest]

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
cnv_granges_sig = subset(cnv_granges,cnv_type != "Neutral")

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
  ggrepel::geom_text_repel(
    data = as.data.frame(cnv_granges_sig[abs(cnv_granges_sig$log2_ratio) > 0.5, ]), # Filter significant CNVs
    aes(
      x = start + (end - start) / 2,  # Midpoint of the CNV
      y = log2_ratio + 0.5,           # Offset to avoid overlap
      label = gene                   # Gene names
    ),
    size = 3,                         # Adjust font size
    angle = 0,                       # Rotate text for better clarity
    hjust = 0.5,                        # Adjust horizontal justification
    vjust = -1,                      # Adjust vertical justification
    color = "black"                   # Label color
  )

print(p)




cnv_sig_df = as.data.frame(cnv_granges_sig)
samples = unique(cnv_sig_df$sample)

for (s in samples) {
  sample_data <- subset(cnv_sig_df, sample == s)
  
  # Plot CNV profile with genes and CNV type-based coloring
  p <- ggplot(sample_data, aes(x = (start + end) / 2, y = log2_ratio, color = cnv_type)) +
    geom_point(size = 1) +
    facet_wrap(~seqnames, scales = "free", ncol = 4) +
    theme_minimal() +
    labs(
      title = paste("CNV Profile for Sample", s),
      x = "Genomic Position",
      y = "Normalized Log2 Ratio"
    ) +
    scale_color_manual(
      values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
      name = "CNV Type"
    ) +
    ggrepel::geom_text_repel(
      data = subset(sample_data, abs(log2_ratio) > 0.1), # Filter significant CNVs
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
    )+
    theme_bw()
  print(p)
  
  # Update the karyogram plot for better visibility
  sample_granges = subset(cnv_granges_sig,sample ==s)
  text_data = subset(sample_data,abs(log2_ratio)>0.3)
  karyo_p <- ggbio::autoplot(
    sample_granges,
    aes(y = log2_ratio, fill = cnv_type),
    layout = "karyogram"
  ) +
    scale_fill_manual(
      values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
      name = "CNV Type"
    ) +
    theme_bw() +
    labs(
      title = paste0("CNV Karyogram for ",s),
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
    ggrepel::geom_text_repel(
      data = text_data, # Filter significant CNVs
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
  
  print(karyo_p)
  
}



