# Install and load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c( "RefFreeEWAS"))

library(TOAST)
library(RefFreeEWAS)
library(limma)
library(reshape2)
library(ggplot2)

# Load your methylation data
# Assuming pm_matrix is your methylation matrix and design_df is your sample information dataframe

pm_matrix <- as.matrix(read.csv("path_to_your_pm_matrix.csv", row.names = 1))
design_df <- read.csv("path_to_your_design_file.csv")

# Filter out CpG sites with low variance across samples
variance_threshold <- 0.02  # Adjust as necessary
pm_var <- apply(pm_matrix, 1, var)
pm_matrix_filtered <- pm_matrix[pm_var > variance_threshold, ]

# Check dimensions
dim(pm_matrix_filtered)

# Reference-Free Deconvolution
# Use RefFreeCellMix for reference-free deconvolution
K <- 6  # Number of cell types to estimate
outT_RF <- RefFreeCellMix(pm_matrix_filtered, K = K)
estProp_RF <- outT_RF$Omega

# Visualize the Results
# Combine results into a single data frame for visualization
cell_type_proportions <- data.frame(
  Sample = rownames(estProp_RF),
  estProp_RF
)

# Melt the data frame for ggplot2
cell_type_proportions_melted <- melt(cell_type_proportions, id.vars = "Sample")

# Plot the cell type proportions
ggplot(cell_type_proportions_melted, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Estimated Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot
ggsave("cell_type_proportions_plot.pdf")
