# Manhattan Plot
library(EnhancedVolcano);library(pheatmap);library(ggpubr)
setwd('/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/')

dmr = readRDS("./results/tables/dmr/dmr_reg_obj.Rds")
load("./data/processed/normalized_methyl.rdata")

dmr.df = as.data.frame(dmr)

# Heatmap of Methylation
pm <- percMethylation(meth_filt_reg,rowids = T)
annot = data.frame(diagnosis = methyldata$entity)
rownames(annot) = methyldata$samples
heatmap <- pheatmap(pm, cluster_rows = TRUE, cluster_cols = TRUE, 
                    show_rownames = FALSE, show_colnames = TRUE, annotation_col = annot,
                    main = "Methylation Heatmap for Tumor vs Control")
## save heatmap
ggsave(filename = "./results/figures/allmeth_reg_heatmap.pdf", plot = heatmap$gtable,
       width = 12, height = 12)

#Heatmap of Top Differentially Methylated Regions
# Create unique identifiers for the DMRs
dmr_ids <- paste(dmr.df$chr, dmr.df$start, dmr.df$end, sep = ".")

rownames(pm) = paste(meth_filt_reg$chr, meth_filt_reg$start, meth_filt_reg$end, sep = ".")

# Filter the pm matrix to include only rows corresponding to DMRs
pm_filtered <- pm[rownames(pm) %in% dmr_ids, ]

# Ensure the filtered pm has the same number of rows as the unique DMRs
print(nrow(pm_filtered) == length(unique(dmr_ids)))

# Order the rows of pm_filtered to match the order in dmr.df
pm_filtered <- pm_filtered[match(dmr_ids, rownames(pm_filtered)), ]


# Plot the heatmap
ht_dmr = pheatmap(
  pm_filtered,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = F,
  clustering_distance_cols = "maximum",
  show_colnames = T,
  treeheight_col = 3,
  annotation_col = annot,annotation_names_col = F,
  main = "Heatmap of Percent Methylation in DMRs",
  color = colorRampPalette(c("DarkBlue", "white", "Orange"))(50)
)

ggsave(filename = "./results/figures/dmr/dmr_reg_heatmap.pdf",plot = ht_dmr$gtable, width = 12, height = 8)


# Volcano Plot
dmr.df = as.data.frame(dmrs_annotated)
topdata = subset(dmr.df,!is.na(annot.symbol)) %>% arrange(qvalue)

volcano.p = ggplot(topdata,aes(x = meth.diff,y = -log10(qvalue)))+
  geom_point(aes(fill = annot.type),shape=21,color='black',size=3)+
  geom_hline(yintercept = 1.3,lty =2)+
  geom_text_repel(data = topdata,aes(label = annot.symbol))+
  theme_bw()+
  labs(x  = 'Methylation difference',
       title = 'Top annotated genes with DMRs')
ggsave("./results/figures/dmr/volcano_topGenes_region.pdf",volcano.p,width = 10,height = 7)

##manhattan plot
manhattan_plot <- ggplot(topdata, aes(x = seqnames, y = -log10(qvalue),
                                  color = annot.type)) +
  geom_point(alpha = 0.5,show.legend = T) +
  geom_text_repel(aes(label = annot.symbol))+
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5))+
  labs(title = "Manhattan Plot for Tumor vs Control", x = "Chromosome", y = "-log10(qvalue)")

# Combine plots
combined_plot <- ggarrange(volcano.p,manhattan_plot, ncol = 1)

# Save combined plot
ggsave("./results/figures/dmr/Combined_Plot_volcano_manhattan.pdf", combined_plot, width = 12, height = 7)


## Boxplot of Methylation Differences by Annotation Type
annot.box = ggplot(dmr.df, aes(x = annot.type, y = meth.diff, fill = annot.type)) +
  geom_boxplot() +
  labs(title = "Boxplot of Methylation Differences by Annotation Type",
       x = "Annotation Type",
       y = "Methylation Difference") +
  scale_fill_brewer(palette = 'Set3')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "./results/figures/dmr/dmr_reg_boxplot_annotation_type.pdf", annot.box,width = 12, height = 8)

## genomic distribution of DMRs
library(Gviz)
library(GenomicRanges)

# Convert the data frame to GRanges
dmr_gr <- GRanges(
  seqnames = dmr.df$seqnames,
  ranges = IRanges(start = dmr.df$start, end = dmr.df$end),
  strand = dmr.df$strand,
  score = -log10(dmr.df$pvalue),
  annotation = dmr.df$annot.type
)

# Create a Genome Axis Track
gtrack <- GenomeAxisTrack()

# Create a Data Track for DMRs
dtrack <- DataTrack(
  range = dmr_gr,
  genome = "hg38",
  name = "DMRs",
  chromosome = seqlevels(dmr_gr),
  col.line = "blue",
  col = "blue",
  type = "p",
  baseline = 0,
  col.baseline = "black",
  lwd.baseline = 1,
  pch = 16,
  cex = 0.5
)

# Plot the DMRs with the Genome Axis
pdf(file = "./results/figures/dmr_reg_genomic_distribution.pdf", width = 12, height = 6)
plotTracks(list(gtrack, dtrack), from = min(dmr.df$start), to = max(dmr.df$end), 
           main = "Genomic Distribution of DMRs", 
           col.main = "black", 
           cex.main = 1.2)
dev.off()


## position of DMRs on circular plot
# Create a data frame for circlize
circlize_data <- data.frame(
  chr = dmr.df$seqnames,
  start = dmr.df$start,
  end = dmr.df$end,
  value = dmr.df$meth.diff
)

# Plot circular visualization
circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(circlize_data, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = ifelse(value > 0, "red", "blue"))
})

pdf(file = "./results/figures/dmr/dmr_reg_circular_plot.pdf", width = 18, height = 5)
circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(circlize_data, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = ifelse(value > 0, "red", "blue"))
})
dev.off()
