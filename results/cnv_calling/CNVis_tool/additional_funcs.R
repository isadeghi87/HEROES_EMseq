hc <- hclust(dist(cnv_matrix_filtered), method = "ward.D2")
plot(hc, main = "Sample Clustering Based on CNV Patterns")

high_variance_regions <- cnv_data_filtered[, .(
  mean_log2_ratio = mean(log2_ratio_median, na.rm = TRUE),
  variance = var(log2_ratio_median, na.rm = TRUE)
), by = genomic_bin]

ggplot(high_variance_regions, aes(x = genomic_bin, y = variance)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Variance of CNVs Across Genomic Regions")


##mutation burden ####
# Compute CNV burden for each chromosome
cnv_burden <- cnv_data_filtered[, .(
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
cnv_data_filtered$diagnosis <- factor(cnv_data_filtered$diagnosis)

# Perform t-tests for each genomic bin
differential_cnv <- cnv_data_filtered[, {
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
library(ggplot2)

# Assuming 'umap_result' contains the UMAP coordinates and 'annotations' has sample IDs
umap_data <- data.frame(
  Sample = annotations$sample, # Add sample names
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Diagnosis = annotations$diagnosis,
  Entity = annotations$entity
)

# Plot UMAP with sample labels
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Diagnosis, shape = Entity)) +
  geom_point(size = 3, alpha = 0.8) +  # Add points
  geom_text(aes(label = Sample), hjust = -0.2, vjust = -0.2, size = 3) +  # Add labels
  theme_minimal() +
  labs(
    title = "UMAP of CNV Data",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors for diagnosis
  #scale_shape_manual(values = c(16, 17, 18)) +  # Customize shapes for status
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )+
  theme_bw()



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
cnv_data_filtered$timepoint <- sapply(cnv_data_filtered$sample, extract_timepoint)

# View the updated annotations data frame
head(cnv_data_filtered)


#For samples with repeated measures (e.g., relapse vs. primary), analyze CNV progression over time.
ggplot(cnv_data_filtered, aes(timepoint, log2_ratio_median, color = genomic_bin)) +
  geom_line(show.legend = F) +
  theme_minimal() +
  labs(title = "Temporal Dynamics of CNVs")




