# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("EpiDISH")
BiocManager::install("methylKit")
BiocManager::install("reshape2")

library(EpiDISH)
library(methylKit)
library(dplyr)
library(reshape2)
library(ggplot2)

# Set working directory
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/')

# Load your methylation data (assuming `meth` is your methylation object)
meth <- readRDS('./data/processed/meth_unite.rds')

# Get the percent methylation matrix
pm <- percMethylation(meth)

# Extract chromosome and start position from the methylBase object
chr <- meth$chr
start <- meth$start

# Convert the row names to the format "chr:position"
rownames(pm) <- paste0(chr, ":", start)

# EpiDISH expects a matrix with CpG sites as rows and samples as columns
pm_matrix <- as.matrix(pm)

# Check for and remove columns with constant or nearly constant values
pm_matrix <- pm_matrix[, apply(pm_matrix, 2, var) > 1e-5]

# Remove rows with constant or nearly constant values
pm_matrix_filtered <- pm_matrix[apply(pm_matrix, 1, var) > 1e-5, ]

# Check the dimensions of the filtered matrix
print(dim(pm_matrix_filtered))

# Load reference data for RPC (Robust Partial Correlation)
data(centEpiFibIC.m)

# Run EpiDISH using RPC method
result_rpc <- tryCatch({
  epidish(beta.m = pm_matrix_filtered, ref.m = centEpiFibIC.m, method = "RPC")
}, error = function(e) {
  message("Failed to run EpiDISH with RPC method: ", e$message)
  NULL
})

if (!is.null(result_rpc)) {
  cell_type_proportions_rpc <- as.data.frame(result_rpc$estF)
  cell_type_proportions_rpc$Method <- "RPC"
}

# Combine results into one data frame if RPC was successful
if (!is.null(result_rpc)) {
  cell_type_proportions <- cell_type_proportions_rpc
  cell_type_proportions$Sample <- rownames(cell_type_proportions)
  
  # Melt the data frame for ggplot2
  melted_df <- melt(cell_type_proportions, id.vars = c("Sample", "Method"))
  
  # Plot the cell type proportions
  ggplot(melted_df, aes(x = Sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Method, scales = "free_x") +
    labs(title = "Estimated Cell Type Proportions", x = "Sample", y = "Proportion", fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  message("RPC method failed. Please check the data and try again.")
}

