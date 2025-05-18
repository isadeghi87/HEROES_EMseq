## pca analysis
setwd('/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/')

## read methyl processed data 
load('./data/processed/normalized_methyl_region.rdata')

# Load necessary libraries
library(ggplot2)
library(reshape2)

pm <- percMethylation(meth_filt_reg,rowids = T)
# Transpose the pm matrix to have samples as rows and features as columns
pm_t <- t(pm)

# Perform PCA
pca_result <- prcomp(pm_t)

# Extract PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
pca_data$diagnosis = methyldata$entity[match(rownames(pca_data),methyldata$samples)]

# Create a scree plot
(pca_plot <- ggplot(pca_data,aes(x = PC1,y=PC2))+
  geom_point(aes(fill = diagnosis),color = 'black',shape =21,size=3)+
  labs(title = 'PCA using methylation percentage')+
  # geom_text(aes(label = Sample), hjust = 0.5, vjust = 1.5, size = 3) +
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) +
  theme_bw())


# Save the PCA plot
ggsave(filename = "./results/figures/dmr/pca_plot.pdf", plot = pca_plot, width = 8, height = 5)
