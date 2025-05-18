
# Install and load the rrcov package
install.packages("rrcov")
library(rrcov)

setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/')

# Load your methylation data (assuming `meth` is your methylation object)
meth <- readRDS('./data/processed/meth_unite.rds')

# Get the percent methylation matrix
pm <- percMethylation(meth)

# Transpose the pm matrix to have samples as rows and features as columns
pm_t <- t(pm)

# Perform robust PCA
rpca_result <- PcaHubert(pm_t, k = 2, alpha = 0.5)

# Get scores and distances
scores <- rpca_result@scores
distances <- rpca_result@od

# Identify outliers based on distances
outliers <- rownames(pm_t)[distances > quantile(distances, 0.95)]

# Remove outlier samples from the methylation object
meth_filtered <- meth[!(meth@sample.ids %in% outliers)]

# Verify the removal
meth_filtered


# Calculate correlation matrix
cor_matrix <- cor(pm)

# Plot the heatmap of the correlation matrix
heatmap(cor_matrix, symm = TRUE, main = "Correlation Matrix of Samples")
