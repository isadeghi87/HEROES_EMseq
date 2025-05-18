
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

setwd('/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/')

## read methyl processed data 
load('./data/processed/normalized_methyl_region.rdata')

# Calculate DMRs using multiple methods
# Method 1: methylKit
myDiff <- calculateDiffMeth(meth_filt_reg)

# Get differentially methylated regions (DMRs)
dmr <- getMethylDiff(myDiff,difference=10,qvalue=0.05)

# Overview of percentage hyper and hypo CpGs per chromosome.
pdf("./results/figures/dmr/diffMethPerChr_region.pdf",width = 10,height = 8)
diffMethPerChr(myDiff,qvalue.cutoff = 0.05,meth.cutoff = 10)
dev.off()

# Save DMRs to a file
write.table(dmr, file = "./results/tables/dmr/DMRs_region.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Annotate DMRs with annotatr
# Define the annotations to be used
anns = builtin_annotations()
annotations =  anns[grep('hg38',anns)]
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
print(annotations)

# Build the annotations
annots <- build_annotations(genome = 'hg38', annotations = annotations)

# Convert DMRs to GRanges
dmrs_gr <- as(dmr, "GRanges")

# Annotate DMRs
dmrs_annotated <- annotate_regions(
  regions = dmrs_gr,
  annotations = annots,
  ignore.strand = TRUE,
  quiet = FALSE
)

saveRDS(dmrs_annotated,"./results/tables/dmr/dmr_region_annotated.Rds")
saveRDS(dmr,"./results/tables/dmr/dmr_reg_obj.Rds")


# Visualize annotation summaries
pdf(file = './results/figures/dmr/dmrs_annotation_summary.pdf', width = 10, height = 7)
plot_annotation(
  annotated_regions = dmrs_annotated,
  annotation_order = annotations
)
dev.off()

# Detailed annotation
dmrs_annotated_detail <- summarize_annotations(dmrs_annotated)

custom_colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", 
  "#ffd92f", "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", 
  "#7570b3", "#e7298a", "#66a61e", "#e6ab02","#ff5b22"
)
annot_plot = ggplot(data = dmrs_annotated_detail,
                    mapping = aes(y = reorder(annot.type,n),x = n, fill=annot.type))+
  geom_bar(stat = 'identity', position = 'dodge',width = 0.7,color = 'black',show.legend = F)+
  scale_fill_manual(values = custom_colors)+
  scale_x_log10()+
  labs(title = 'Annotation of DMRs',
       x = 'n',
       y = '')+
  theme_bw()+
  theme(title = element_text(hjust = 0.5 ,face = 'bold'))

ggsave(filename = './results/figures/dmr/annotation_dmr_reg_plot.pdf',plot = annot_plot,
       width = 8,height = 5)


# Save detailed annotation results
write.csv(as.data.frame(dmrs_annotated_detail), file = "./results/figures/dmr/dmrs_reg_detailed_annotation.csv")

dmr_genes  = unique(dmr_annot$annot$symbol)




