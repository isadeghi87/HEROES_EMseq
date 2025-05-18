#step1: Identify top DMR genes based on qvalue####
top_dmr_genes <- dmr.df %>%
  filter(!is.na(annot.symbol)) %>%
  group_by(annot.symbol) %>%
  summarize(min_qvalue = min(qvalue)) %>%
  arrange(min_qvalue) %>%
  top_n(-5, min_qvalue)  # Select top 5 genes with lowest qvalue

# Identify top annotation types based on frequency
top_annotation_types <- dmr.df %>%
  group_by(annot.type) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(5, count)  # Select top 5 annotation types with highest frequency

# Print the top DMR genes and annotation types
print(top_dmr_genes)
print(top_annotation_types)

#step2: Create Tracks for the Genome, Gene Model, and Data####
# Function to create DataTrack for a given gene
create_dtrack_for_gene <- function(gene_name, dmr.df) {
  gene_dmr.df <- dmr.df %>%
    filter(annot.symbol == gene_name)
  
  gr <- GRanges(
    seqnames = gene_dmr.df$seqnames,
    ranges = IRanges(start = gene_dmr.df$start, end = gene_dmr.df$end),
    strand = gene_dmr.df$strand,
    score = gene_dmr.df$meth.diff,  # Using meth.diff for the score
    annotation = gene_dmr.df$annot.type
  )
  
  DataTrack(
    range = gr,
    genome = "hg38",
    name = paste(gene_name, "DMRs"),
    type = "h",
    col = "blue",
    baseline = 0,
    col.baseline = "black",
    lwd.baseline = 1,
    pch = 16,
    cex = 0.5
  )
}

# Create tracks for top DMR genes
top_dmr_gene_tracks <- lapply(top_dmr_genes$annot.symbol, function(gene) {
  create_dtrack_for_gene(gene, dmr.df)
})

# Create Genome Axis Track
gtrack <- GenomeAxisTrack()

# Attempt to create Ideogram Track with error handling
itrack <- tryCatch({
  IdeogramTrack(genome = "hg38", chromosome = "chr5")
}, error = function(e) {
  message("Failed to create IdeogramTrack: ", e$message)
  NULL
})

# Create Gene Region Track
ah <- AnnotationHub()
query(ah, "TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- ah[["AH52263"]]

grtrack <- GeneRegionTrack(txdb, genome = "hg38", chromosome = "chr5", name = "Genes")

# Step 3: Focus on the Region Around the Top DMR Genes or Annotation Types #### 

# Define function to plot gene region
plot_gene_region <- function(gene_name, gene_dtrack, grtrack, gtrack, itrack, output_file) {
  gene_range <- range(gene_dtrack@range)
  from <- min(start(gene_range))
  to <- max(end(gene_range))
  
  pdf(file = output_file, width = 14, height = 8)  # Increase the size of the output PDF
  if (!is.null(itrack)) {
    plotTracks(
      list(itrack, gtrack, grtrack, gene_dtrack),
      from = from - 10000,
      to = to + 10000,
      main = paste(gene_name, "DMRs"),
      col.main = "black",
      cex.main = 1.2
    )
  } else {
    plotTracks(
      list(gtrack, grtrack, gene_dtrack),
      from = from - 10000,
      to = to + 10000,
      main = paste(gene_name, "DMRs"),
      col.main = "black",
      cex.main = 1.2
    )
  }
  dev.off()
}

# Plot for top DMR genes
output_dir <- "./results/figures/dmr/top_dmr_genes/"
dir.create(output_dir, showWarnings = FALSE)
for (i in seq_along(top_dmr_gene_tracks)) {
  gene_name <- top_dmr_genes$annot.symbol[i]
  dtrack <- top_dmr_gene_tracks[[i]]
  output_file <- file.path(output_dir, paste0(gene_name, "_dmr_plot.pdf"))
  plot_gene_region(gene_name, dtrack, grtrack, gtrack, itrack, output_file)
}

# Function to create DataTrack for a given annotation type
create_dtrack_for_annotation <- function(annotation_type, dmr.df) {
  annot_dmr.df <- dmr.df %>%
    filter(annot.type == annotation_type) %>%
    top_n(100, -qvalue)  # Limit to top 100 DMRs for each annotation type
  
  gr <- GRanges(
    seqnames = annot_dmr.df$seqnames,
    ranges = IRanges(start = annot_dmr.df$start, end = annot_dmr.df$end),
    strand = annot_dmr.df$strand,
    score = annot_dmr.df$meth.diff,  # Using meth.diff for the score
    annotation = annot_dmr.df$annot.type
  )
  
  DataTrack(
    range = gr,
    genome = "hg38",
    name = paste(annotation_type, "DMRs"),
    type = "h",
    col = "red",
    baseline = 0,
    col.baseline = "black",
    lwd.baseline = 1,
    pch = 16,
    cex = 0.5
  )
}

# Create tracks for top annotation types
top_annotation_tracks <- lapply(top_annotation_types$annot.type, function(annot) {
  create_dtrack_for_annotation(annot, dmr.df)
})

# Plot for top annotation types
output_dir_annot <- "./results/figures/dmr/top_dmr_annotations/"
dir.create(output_dir_annot, showWarnings = FALSE)
for (i in seq_along(top_annotation_tracks)) {
  annotation_type <- top_annotation_types$annot.type[i]
  dtrack <- top_annotation_tracks[[i]]
  output_file <- file.path(output_dir_annot, paste0(annotation_type, "_dmr_plot.pdf"))
  plot_gene_region(annotation_type, dtrack, grtrack, gtrack, itrack, output_file)
}
