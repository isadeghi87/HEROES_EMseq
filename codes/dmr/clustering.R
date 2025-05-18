


# Load required packages
library("methylKit")
library("genomation")
library("GenomicRanges")
library(bsseq)
library(magrittr)
library(readr)
library(dmrseq)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(Gviz)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(edgeR)
library(DSS)
library(ChAMP)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)
library(pheatmap)
library(ggpubr)


## read methyl processed data 
load('./data/processed/normalized_methyl_region.rdata')


# Get percent methylation matrix
pm <- percMethylation(meth_filt_reg,rowids = T)

# Hierarchical clustering
hc <- hclust(dist(pm))
plot(hc)

# Calculate standard deviation of CpGs
sds <- rowSds(pm)   
hist(sds, breaks = 100)
meth <- meth_filt_reg[sds > 2]

pdf("./results/figures/dmr/cpg_methyl_clustering.pdf", width = 10, height = 20)
clusterSamples(meth, dist = "correlation", method = 'ward', plot = TRUE)
dev.off()



# Calculate DMRs using multiple methods
# Method 1: methylKit
myDiff <- calculateDiffMeth(meth)
myDiff
myDiff25p <- getMethylDiff(myDiff,difference=10,qvalue=0.05)

# Annotate DMRs with annotatr
# Define the annotations to be used
annotations <- c(
  "hg38_basicgenes",
  "hg38_cpgs",
  "hg38_genes_intergenic",
  "hg38_genes_intronexonboundaries",
  "hg38_genes_5UTRs",
  "hg38_genes_3UTRs"
)

# Build the annotations
annots <- build_annotations(genome = 'hg38', annotations = annotations)

# Convert DMRs to GRanges
dmrs_gr <- as(myDiff25p, "GRanges")

# Annotate DMRs
dmrs_annotated <- annotate_regions(
  regions = dmrs_gr,
  annotations = annots,
  ignore.strand = TRUE,
  quiet = FALSE
)

# Visualize annotation summaries
pdf(file = './results/figures/dmr/dmrs_annotation_summary.pdf', width = 10, height = 7)
plot_annotation(
  annotated_regions = dmrs_annotated,
  annotation_order = annotations
)
dev.off()

# Detailed annotation
dmrs_annotated_detail <- summarize_annotations(dmrs_annotated)

annot_plot = ggplot(data = dmrs_annotated_detail,
       mapping = aes(y = reorder(annot.type,n),x = n, fill=annot.type))+
  geom_bar(stat = 'identity', position = 'dodge',width = 0.7,color = 'black')+
  scale_fill_brewer(palette = 'Set3')+
  scale_x_log10()+
  labs(title = 'Annotation of DMRs',
       x = 'n',
       y = '')+
  theme_bw()+
  theme(title = element_text(hjust = 0.5 ,face = 'bold'))

ggsave(filename = './results/figures/dmr/annotation_dmr_plot.pdf',plot = annot_plot,
       width = 8,height = 5)


# Save detailed annotation results
write.csv(as.data.frame(dmrs_annotated_detail), file = "./results/figures/dmr/dmrs_detailed_annotation.csv")

# Initialize the STRINGdb object
string_db <- STRINGdb$new(version = "11", species = 9606, score_threshold = 400)

# Map the gene symbols to STRING IDs
dmr_genes_string <- string_db$map(data.frame(gene = dmr_genes), "gene", removeUnmapped = FALSE)

# Remove unmapped entries
dmr_genes_string_mapped <- na.omit(dmr_genes_string[!is.na(dmr_genes_string$STRING_id), ])
# Ensure unique vertex names
dmr_genes_string_mapped$gene <- make.unique(dmr_genes_string_mapped$gene)

# Check the percentage of mapped genes
mapped_percentage <- nrow(dmr_genes_string_mapped) / nrow(dmr_genes_string) * 100
cat("Mapped percentage:", mapped_percentage, "%\n")

pdf(file = "./results/figures/dmr/string_network.pdf", width = 12, height = 12)
plot(network, main = "Protein-Protein Interaction Network of DMR-associated Genes")
dev.off()
# Get interactions for mapped genes
interactions <- string_db$get_interactions(dmr_genes_string_mapped$STRING_id)

# Inspect the interaction data
head(interactions)
# Ensure all vertices in the interaction data are present in the vertex data frame
# Ensure all vertices in the interaction data are present in the vertex data frame
valid_vertices <- union(interactions$from,interactions$to)

# Filter the vertex data frame to include only valid vertices
dmr_genes_string_mapped <- dmr_genes_string_mapped[dmr_genes_string_mapped$STRING_id %in% valid_vertices, ]

interactions_filtered <- interactions[grep(paste0(valid_vertices,collapse = '|'),interactions$from),]
interactions_filtered <- interactions[grep(paste0(valid_vertices,collapse = '|'),interactions$to),]

# Ensure there are no duplicate vertex names
dmr_genes_string_mapped <- dmr_genes_string_mapped[!duplicated(dmr_genes_string_mapped$STRING_id), ]


# Create an igraph object
g <- graph_from_data_frame(d = interactions_filtered, directed = FALSE, vertices = dmr_genes_string_mapped)

# Simplify the graph
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)


# Save results
write.csv(as.data.frame(myDiff25p), file = "./results/figures/dmr/significant_dmrs.csv")
write.csv(as.data.frame(ego), file = "./results/figures/dmr/go_enrichment.csv")
write.csv(as.data.frame(kegg), file = "./results/figures/dmr/kegg_enrichment.csv")

# Export significant DMRs as bed file for visualization in genome browsers
export.bed(as.data.frame(dmrs), con = "./results/figures/dmr/significant_dmrs.bed")


# Violin plots for selected DMRs
library(ggplot2)
selected_dmrs <- head(dmr.df, 10)
for (i in 1:nrow(selected_dmrs)) {
  plot_data <- data.frame(
    sample = rep(methyldata$samples, each = 2),
    treatment = rep(methyldata$diagnosis, each = 2),
    methylation = as.numeric(pm[selected_dmrs$chr[i], selected_dmrs$start[i]:selected_dmrs$end[i]])
  )
  p <- ggplot(plot_data, aes(x = as.factor(treatment), y = methylation, fill = as.factor(treatment))) +
    geom_violin() +
    theme_minimal() +
    ggtitle(paste("DMR in", selected_dmrs$chr[i], ":", selected_dmrs$start[i], "-", selected_dmrs$end[i]))
  #ggsave(filename = paste0("./results/figures/dmr/violin_plot_", i, ".pdf"), plot = p, width = 8, height = 6)
}



# Generate a report summarizing the analyses
library(knitr)
library(rmarkdown)
rmarkdown::render("report_template.Rmd", output_file = "results/report.html")
