library(clusterProfiler);library(ReactomePA)
library(ggplot2);library(RColorBrewer)

## pathway analysis for DMRs 
dmrs_annotated = read_rds("./results/tables/dmr/dmr_region_annotated.rds")

# Extract gene symbols from annotated DMRs
dmr_genes <- unique(dmrs_annotated$annot$symbol)

# Convert to Entrez IDs
dmr_entrez <- bitr(dmr_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Enrichment analysis
go_enrichment <- enrichGO(gene = dmr_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP")

ego <- enrichGO(gene = dmr_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",pvalueCutoff = 0.01,readable = T)

ego_plot = dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment Analysis of DMR-associated Genes")

# Reactome Pathway Analysis
reactome_enrichment <- enrichPathway(gene = dmr_entrez$ENTREZID, pvalueCutoff = 0.05, readable = TRUE)

# Plot Reactome pathways
react_plot = dotplot(reactome_enrichment, showCategory = 20, title = "Reactome Pathways for Tumor vs Control")

go_plots = ego_plot+react_plot
ggsave("./results/figures/dmr/go_enrichment_dotplot.pdf",go_plots,width = 12,height = 8)

gene_id <- na.omit(unique(dmrs_annotated$annot$gene_id))
kegg <- enrichKEGG(gene =dmr_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
kegg_res = kegg@result[1:10,]

keggplot = ggplot(kegg_res,aes(y = reorder(Description,-log10(pvalue)),
                    x= -log10(pvalue),
                    fill = Description))+
  geom_bar(stat = 'identity',width = 0.5,color = 'black',show.legend = F)+
  labs(title = 'KEGG pathway enrichment for DMR genes',
       y = '')+
  scale_fill_brewer(palette = 'Set3')+
  geom_vline(xintercept = 1.3,lty=2)+
  theme_bw()
  
ggsave( "./results/figures/dmr/kegg_enrichment_dotplot.pdf",keggplot,width = 8, height = 5)
