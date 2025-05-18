# Load necessary libraries
library(methylKit)
library(ggplot2)
library(dplyr)
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/')

# Load your methylation data (assuming `meth` is your methylation object)
meth <- readRDS('./data/processed/meth_unite.rds')

# Get the percent methylation matrix
pm <- percMethylation(meth)

# Transpose the pm matrix to have samples as rows and features as columns
pm_t <- t(pm)

# Calculate the Mahalanobis distance for each sample
mahalanobis_distance <- mahalanobis(pm_t, colMeans(pm_t), cov(pm_t))

# Calculate the p-value for each distance
p_values <- pchisq(mahalanobis_distance, df = ncol(pm_t), lower.tail = FALSE)

# Set a significance level (e.g., 0.01)
alpha <- 0.01

# Identify outliers
outliers <- names(p_values[p_values < alpha])

# Remove outlier samples from the methylation object
meth_filtered <- meth[!(meth@sample.ids %in% outliers)]

# Verify the removal
meth_filtered
