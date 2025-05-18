library(GenomicRanges)
library(annotatr)
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(tibble)
library(ggrepel)

setwd('/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/')

## read methyl processed data 
load('./data/processed/normalized_methyl_region.rdata')

## subset for sample
#  Extract data from the methylBase object

# Subset the GRanges object to keep only rows for your sample of interest.
sample_of_interest <- "I034-042_plasma-01-02"
meth_sample <- meth_filt_reg[meth_filt_reg@sample.ids==sample_of_interest]

# Optionally, convert to a data frame for easier downstream processing:
meth_sample_df <- as.data.frame(meth_sample)

# Build a set of annotations; here we use basic genes for hg38.
annotations = c( "hg38_genes_1to5kb",               "hg38_genes_promoters",           
                 "hg38_genes_cds",                  "hg38_genes_5UTRs",               
                 "hg38_genes_exons",         
                 "hg38_genes_introns",              "hg38_genes_intronexonboundaries",
                 "hg38_genes_exonintronboundaries", "hg38_genes_3UTRs",               
                 "hg38_genes_intergenic",           "hg38_cpg_islands",               
                 "hg38_cpg_shores",                 "hg38_cpg_shelves",               
                 "hg38_cpg_inter",                  "hg38_enhancers_fantom", 
                 "hg38_basicgenes",                
                 "hg38_cpgs" )
annots <- build_annotations(genome = 'hg38', annotations = annotations)

# Annotate the methylation regions. This will add annotation columns to the GRanges object.
gr <- as(meth_sample_df, "GRanges")
methyl_annot <- annotate_regions(
  regions = gr,
  annotations = annots,
  ignore.strand = TRUE
)

# Check available annotation types
unique(mcols(methyl_annot)$annot$type)

## convert to df

meth_df <- as.data.frame(methyl_annot) %>%
  filter(!is.na(annot.symbol)) %>%
  mutate(gene = annot.symbol) %>%
  mutate(
    total_coverage = rowSums(dplyr::select(., starts_with("coverage")), na.rm = TRUE),
    total_numCs    = rowSums(dplyr::select(., starts_with("numCs")), na.rm = TRUE),
    avg_methyl     = total_numCs / total_coverage,
    methylation_scaled = ((total_numCs / total_coverage)*100 - 50) / 50
  )

### load scmultiome data 
sample_id = 'I034_042_0T1'
pool='pool27'

seurat_obj = readRDS("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/pool_temp/results/10x_multiomics/seurat/processed_obj/pools/pool27/I034_042_0T1_clusters_seurat.Rds")
seurat_obj <- ScaleData(seurat_obj, assay = "RNA")


# Get the unique cell types
seurat_obj$seurat_clusters <- paste0('cl',seurat_obj$seurat_clusters)
seurat_clusters <- unique(seurat_obj$seurat_clusters)

## calculate gene activity if not present
if (!"GeneActivity" %in% names(seurat_obj@assays)) {
  message("GeneActivity assay not found. Calculating gene activity...")
  # Calculate gene activity using Signac.
  # GeneActivity() takes the Seurat object (with ATAC data) and returns a matrix.
  DefaultAssay(seurat_obj)= 'ATAC'
  gene.activities <- GeneActivity(seurat_obj)
  
  # Add the gene activity data as a new assay called "GeneActivity"
  seurat_obj[["GeneActivity"]] <- CreateAssayObject(gene.activities, assay = "GeneActivity")
  
  # Normalize the GeneActivity assay.
  seurat_obj <- NormalizeData(seurat_obj, assay = "GeneActivity")
  seurat_obj <- ScaleData(seurat_obj, assay = "GeneActivity")
}

# -----------------------------------------------------------
# STEP 2: For each annotation type, compute cluster‐based correlations
# -----------------------------------------------------------
# For each annotation type in your methylation data, and for each Seurat cluster,
# compute the Pearson correlation across genes between the methylation values and:
#   - Average RNA expression (from the "RNA" assay)
#   - Average gene activity (from the "GeneActivity" assay, computed above)
results <- map_dfr(annot_types, function(a_type) {
  
  # Subset methylation data for the current annotation type.
  meth_sub <- meth_df %>% filter(annot.type == a_type)
  
  # Aggregate methylation values by gene for this annotation type.
  # (Each gene may have multiple regions; we take the average of the computed avg_methyl values.)
  # meth_agg <- meth_sub %>%
  #   group_by(gene) %>%
  #   summarise(avg_methyl = mean(avg_methyl, na.rm = TRUE)) %>%
  #   ungroup() %>%
  #   mutate(
  #     # Convert to a percentage (0–100)
  #     methylation_pct = avg_methyl * 100,
  #     # Scale to a -1 to +1 range:
  #     # 0%   -> (0 - 50) / 50  = -1
  #     # 50%  -> (50 - 50) / 50 =  0
  #     # 100% -> (100 - 50) / 50 = 1
  #     methylation_scaled = (methylation_pct - 50) / 50
  #   )
  # 
  # For each cell cluster, compute the correlations.
  cluster_res <- map_dfr(seurat_clusters, function(ct) {
    
    # --- Subset the Seurat object to the current cluster ---
    sub_obj <- subset(seurat_obj, seurat_clusters == ct)
    
    ## Compute average RNA expression per gene for this cluster
    DefaultAssay(sub_obj) <- "RNA"
    avg_rna_ct <- AverageExpression(sub_obj, assays = "RNA",slot = "scale.data")$RNA
    rna_df <- avg_rna_ct %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "dummy", values_to = "rna_expr") %>%
      dplyr::select(-dummy)
    
    ## Compute average gene activity (ATAC) per gene for this cluster
    DefaultAssay(sub_obj) <- "GeneActivity"
    avg_atac_ct <- AverageExpression(sub_obj, assays = "GeneActivity",slot = "scale.data")$GeneActivity
    atac_df <- avg_atac_ct %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "dummy", values_to = "atac_expr") %>%
      dplyr::select(-dummy)
    
    # --- Merge the aggregated methylation data with RNA data (by gene) ---
    merged_rna <- inner_join(meth_thub, rna_df, by = "gene")
    if (nrow(merged_rna) >= 2) {
      cor_test_rna <- cor.test(merged_rna$avg_methyl, merged_rna$rna_expr, method = "pearson")
      rna_cor <- cor_test_rna$estimate
      rna_p   <- cor_test_rna$p.value
    } else {
      rna_cor <- NA
      rna_p   <- NA
    }
    
    # --- Merge the aggregated methylation data with ATAC (gene activity) data and compute correlation ---
    merged_atac <- inner_join(meth_agg, atac_df, by = "gene")
    if (nrow(merged_atac) >= 2) {
      cor_test_atac <- cor.test(merged_atac$avg_methyl, merged_atac$atac_expr, method = "pearson")
      atac_cor <- cor_test_atac$estimate
      atac_p   <- cor_test_atac$p.value
    } else {
      atac_cor <- NA
      atac_p   <- NA
    }
    
    # Return the correlation results for this cluster as one row.
    tibble(
      seurat_clusters = ct,
      rna_cor = rna_cor,
      rna_p = rna_p,
      atac_cor = atac_cor,
      atac_p = atac_p
    )
    
    
  })  # End of cluster loop
  
  # Add the current annotation type to the cluster results and return.
  cluster_res %>% mutate(annot_type = a_type)
})

# Optionally, reorder the columns for clarity.
results <- results %>% dplyr::select(annot_type, seurat_clusters, everything())

## ---------------------------
## Step 7: Inspect and Save the Results
## ---------------------------
print(results)

#write.csv(results, "cluster_based_allAnnotation_methylation_rna_atac_correlations.csv", row.names = FALSE)

# Reshape the results from wide to long format.
# We have columns: rna_cor, rna_p, atac_cor, atac_p.
# We'll create a new column 'modality' (RNA vs. ATAC) and two new columns: 'cor' and 'p'.
results_long <- results %>%
  pivot_longer(
    cols = c(rna_cor, rna_p, atac_cor, atac_p),
    names_to = c("modality", ".value"),
    names_pattern = "(rna|atac)_(cor|p)"
  ) %>%
  mutate(modality = if_else(modality == "rna", "RNA", "ATAC"))

# Check the reshaped data.
print(results_long)

# Create a heatmap:
ggplot(results_long, aes(x = seurat_clusters, y = annot_type, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(cor, 2)), size = 3) +  # optional: add correlation value labels
  facet_wrap(~ modality) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  labs(
    x = "Cell Cluster",
    y = "Annotation Type",
    title = "Correlation of Aggregated Methylation with Expression/Activity"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )



### get results for point plot 
point_results <- map_dfr(annot_types, function(a_type) {
  
  # Subset methylation data for the current annotation type.
  meth_sub <- meth_df %>% filter(annot.type == a_type)
  
  # Aggregate methylation values by gene for this annotation type.
  # (Each gene may have multiple regions; we take the average of the computed avg_methyl values.)
  # meth_agg <- meth_sub %>%
  #   group_by(gene) %>%
  #   summarise(avg_methyl = mean(avg_methyl, na.rm = TRUE)) %>%
  #   ungroup() %>%
  #   mutate(
  #     # Convert to a percentage (0–100)
  #     methylation_pct = avg_methyl * 100,
  #     # Scale to a -1 to +1 range:
  #     # 0%   -> (0 - 50) / 50  = -1
  #     # 50%  -> (50 - 50) / 50 =  0
  #     # 100% -> (100 - 50) / 50 = 1
  #     methylation_scaled = (methylation_pct - 50) / 50
  #   )
  
  # For each cell cluster, compute the correlations.
  cluster_res <- map_dfr(seurat_clusters, function(ct) {
    
    # --- Subset the Seurat object to the current cluster ---
    sub_obj <- subset(seurat_obj, seurat_clusters == ct)
    
    ## Compute average RNA expression per gene for this cluster
    DefaultAssay(sub_obj) <- "RNA"
    avg_rna_ct <- AverageExpression(sub_obj, assays = "RNA",slot = "scale.data")$RNA
    rna_df <- avg_rna_ct %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "dummy", values_to = "rna_expr") %>%
      dplyr::select(-dummy)
    
    ## Compute average gene activity (ATAC) per gene for this cluster
    DefaultAssay(sub_obj) <- "GeneActivity"
    avg_atac_ct <- AverageExpression(sub_obj, assays = "GeneActivity",slot = "scale.data")$GeneActivity
    atac_df <- avg_atac_ct %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "dummy", values_to = "atac_expr") %>%
      dplyr::select(-dummy)
    
    merged_df <- meth_sub %>%
      inner_join(rna_df, by = "gene") %>%
      inner_join(atac_df, by = "gene") %>%
      dplyr::select(gene, methylation_scaled, rna_expr, atac_expr) %>%
      # Use mutate() to add the cluster column.
      mutate(cluster = ct)
    
    merged_df
    
  })  # End of cluster loop
  
  # Add the current annotation type to the cluster results and return.
  cluster_res %>% mutate(annot_type = a_type)
})


point_results = point_results[!duplicated(point_results),]

genes_desired <- subset(point_results,abs(methylation_scaled)>0.5)

print(genes_desired)

ggplot(genes_desired, aes(x = methylation_scaled, y = rna_expr)) +
  # Facet by annotation type (rows) and cluster (columns)
  facet_grid(annot_type ~ cluster, scales = "free_y") +
  
  # Shaded regions (with very high transparency)
  geom_rect(aes(xmin = -Inf, xmax = -0.5, ymin = 0, ymax = Inf),
            fill = "red", alpha = 0.005) +
  geom_rect(aes(xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = 0),
            fill = "blue", alpha = 0.005) +
  
  # Scatter points (point size reflects ATAC expression)
  geom_point(aes(size = atac_expr), color = "black", alpha = 0.7) +
  
  # Threshold lines at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  
  # Gene labels (only for genes_desired; remove or adjust if not needed)
  geom_text_repel(aes(label = gene),
                  size = 3,
                  box.padding = unit(0.25, "lines"),
                  point.padding = unit(0.2, "lines"),
                  segment.color = "grey50",
                  na.rm = TRUE,
                  max.overlaps = 15) +
  
  # Legends and labels
  scale_size_continuous(name = "ATAC Expression") +
  labs(x = "Scaled Methylation (-1 to +1)",
       y = "RNA Expression (scaled)",
       title = "Methylation vs. RNA Expression by Annotation and Cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



plots = list()
for(annot in unique(point_results$annot_type)){
  subdata = subset(point_results,annot_type == annot)
  genes_desired <- subset(subdata,abs(methylation_scaled)>0.5)
  
  
  ##plot
  p = ggplot(subdata, aes(x = methylation_scaled, y = rna_expr)) +
    # Facet by annotation type (rows) and cluster (columns)
    facet_grid(~ cluster, scales = "free_y") +
    
    # Shaded regions (with very high transparency)
    geom_rect(aes(xmin = -Inf, xmax = -0.5, ymin = -Inf, ymax = Inf),
              fill = "red", alpha = 0.005) +
    geom_rect(aes(xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "blue", alpha = 0.005) +
    
    # Scatter points (point size reflects ATAC expression)
    geom_point(aes(size = atac_expr), color = "black", alpha = 0.7) +
    
    # Threshold lines at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    
    # Gene labels (only for genes_desired; remove or adjust if not needed)
    geom_text_repel(data =genes_desired, aes(label = gene),
                    size = 3,
                    box.padding = unit(0.25, "lines"),
                    point.padding = unit(0.2, "lines"),
                    segment.color = "grey50",
                    na.rm = TRUE,
                    max.overlaps = 15) +
    
    # Legends and labels
    scale_size_continuous(name = "ATAC Expression") +
    labs(x = "Scaled Methylation (-1 to +1)",
         y = "RNA Expression (scaled)",
         title = annot) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  plots[[annot]] <- p
}

