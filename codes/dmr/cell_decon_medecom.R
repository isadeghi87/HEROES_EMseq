## here we do cell type deconvolution of EMseq using MeDeCom
# Step 1: Install and Load Necessary Libraries

BiocManager::install(c("limma", "quadprog", "nnls"))
BiocManager::install("FlowSorted.Blood.450k")
library(MeDeCom)
library(methylKit)
library(dplyr)
library(reshape2)
library(ggplot2)

## set wd
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/')

# Load your methylation data (assuming `meth` is your methylation object)
meth <- read_rds('./data/processed/meth_unite.rds')

#Step 2: Prepare the Data
# Get the percent methylation matrix
pm <- percMethylation(meth)

# Extract chromosome and start position from the methylBase object
chr <- meth$chr
start <- meth$start

# Convert the row names to the format "chr:position"
rownames(pm) <- paste0(chr, ":", start)

# EpiDISH expects a matrix with CpG sites as rows and samples as columns
pm_matrix <- as.matrix(pm)

#Step 3: Run EpiDISH for Deconvolution
# Load reference data for RPC (Robust Partial Correlation)
data(centEpiFibIC.m)

# Filter out CpG sites with low variance across samples
variance_threshold <- 0.01
pm_var <- apply(pm_matrix, 1, var)
pm_matrix_filtered <- pm_matrix[pm_var > variance_threshold, ]

# Check the dimensions of the filtered matrix
dim(pm_matrix_filtered)

# Filter out CpG sites with low variance across samples
variance_threshold <- 0.02  # Increased threshold for better variability
pm_var <- apply(pm_matrix, 1, var)
pm_matrix_filtered <- pm_matrix[pm_var > variance_threshold, ]

# Check the dimensions of the filtered matrix
dim(pm_matrix_filtered)

# Load reference data for TOAST method
data(centEpiFibIC.m)

# Run EpiDISH using TOAST method
result_toast <- tryCatch({
  epidish(beta.m = pm_matrix_filtered, ref.m = centEpiFibIC.m, method = "RPC")
}, error = function(e) {
  message("Failed to run EpiDISH with TOAST method: ", e$message)
  NULL
})

if (!is.null(result_toast)) {
  cell_type_proportions_toast <- as.data.frame(result_toast$estF)
  cell_type_proportions_toast$Method <- "TOAST"
}

# Combine results into one data frame if TOAST was successful
if (!is.null(result_toast)) {
  cell_type_proportions <- cell_type_proportions_toast
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
  message("TOAST method failed. Please check the data and try again.")
}