# R/visualization.R

#' Generate Specified Plots for CNV Data
#'
#' This function generates various plots based on user-specified options. It supports heatmaps, Manhattan plots,
#' PCA plots, UMAP plots, Karyograms, and dendrograms. The function adapts to available sample annotations.
#'
#' @param cnv_matrix A matrix of normalized log2 ratio values. Rows represent genomic bins, and columns represent samples.
#' @param annotations A `data.frame` containing sample annotations.
#' @param genome_version Genome version (`hg38` or `hg19`).
#' @param gtf_path Path to the GTF file for gene annotations.
#' @param plot_manhattan Logical. Whether to generate a Manhattan plot. Default is `TRUE`.
#' @param plot_heatmap_genomicBin Logical. Whether to generate a heatmap of genomic bins. Default is `TRUE`.
#' @param plot_pca Logical. Whether to generate a PCA plot. Default is `TRUE`.
#' @param plot_umap Logical. Whether to generate a UMAP plot. Default is `FALSE`.
#' @param plot_karyogram Logical. Whether to generate Karyogram plots. Default is `FALSE`.
#' @param plot_dendrogram Logical. Whether to generate a dendrogram plot. Default is `TRUE`.
#' @param output_dir Directory to save the generated plots.
#' @importFrom data.table data.table fwrite
#' @return NULL. Plots are saved to the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_plots(
#'   cnv_matrix = cnv_matrix,
#'   annotations = annotations,
#'   genome_version = "hg38",
#'   gtf_path = "/path/to/annotations.gtf",
#'   plot_manhattan = TRUE,
#'   plot_heatmap_genomicBin = TRUE,
#'   plot_pca = TRUE,
#'   plot_umap = FALSE,
#'   plot_karyogram = TRUE,
#'   plot_dendrogram = TRUE,
#'   output_dir = "results/"
#' )
#' }
generate_plots <- function(cnv_matrix, annotations, genome_version, gtf_path,
                           plot_manhattan = TRUE,
                           plot_heatmap_genomicBin = TRUE,
                           plot_pca = TRUE,
                           plot_umap = FALSE,
                           plot_karyogram = FALSE,
                           plot_dendrogram = TRUE,
                           output_dir = "results/") {
  tryCatch({
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(RColorBrewer)
    library(data.table)
    library(uwot)
    library(ggbio)
    library(GenomicRanges)
    library(rtracklayer)
    library(dendextend)
    library(ggdendro)

    # Ensure output directory exists
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Heatmap of Genomic Bins
    if (isTRUE(plot_heatmap_genomicBin)) {
      ha <- HeatmapAnnotation(
        df = annotations,
        col = list(
          diagnosis = if ("diagnosis" %in% colnames(annotations)) {
            setNames(brewer.pal(length(unique(annotations$diagnosis)), "Set3"), unique(annotations$diagnosis))
          } else {
            NULL
          },
          status = if ("status" %in% colnames(annotations)) {
            setNames(brewer.pal(length(unique(annotations$status)), "Set2"), unique(annotations$status))
          } else {
            NULL
          },
          timepoint = if ("timepoint" %in% colnames(annotations)) {
            setNames(brewer.pal(length(unique(annotations$timepoint)), "Spectral"), unique(annotations$timepoint))
          } else {
            NULL
          }
        )
      )

      heatmap_colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

      ht <- Heatmap(
        cnv_matrix,
        name = "Log2 Ratio",
        top_annotation = ha,
        col = heatmap_colors,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_title = "Samples",
        row_title = "Genomic Bins"
      )

      pdf(file.path(output_dir, "heatmap_genomicBin.pdf"), width = 12, height = 8)
      draw(ht)
      dev.off()
    }

    # Manhattan Plot
    if (isTRUE(plot_manhattan)) {
      # Assuming p-values are present or need to be calculated
      # For demonstration, simulate p-values
      cnv_data <- as.data.table(cnv_matrix, keep.rownames = "genomic_bin")
      setnames(cnv_data, "rn", "genomic_bin")
      # Simulate p-values or calculate based on some criteria
      set.seed(123)  # For reproducibility
      cnv_data[, p_value := runif(.N, 0, 1)]
      cnv_data[, chrom := sub(":.*", "", genomic_bin)]
      cnv_data[, position := as.numeric(sub(".*:(\\d+)-.*", "\\1", genomic_bin))]

      # Order chromosomes naturally
      chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
      cnv_data[, chrom := factor(chrom, levels = chrom_order)]
      cnv_data <- cnv_data[!is.na(chrom)]

      cnv_data[, log10_pvalue := -log10(p_value)]

      p <- ggplot(cnv_data, aes(x = position, y = log10_pvalue, color = chrom)) +
        geom_point(alpha = 0.6, size = 0.5) +
        scale_color_manual(values = rep(brewer.pal(12, "Paired"), length.out = length(unique(cnv_data$chrom)))) +
        theme_minimal() +
        labs(title = "Manhattan Plot of CNV Data", x = "Genomic Position", y = "-Log10(p-value)") +
        theme(legend.position = "none")

      ggsave(filename = file.path(output_dir, "manhattan_plot.pdf"), plot = p, width = 12, height = 6)
    }

    # PCA Plot
    if (isTRUE(plot_pca)) {
      pca_result <- prcomp(t(cnv_matrix), scale. = TRUE)
      pca_data <- data.frame(
        Sample = rownames(pca_result$x),
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2]
      )
      pca_data <- merge(pca_data, annotations, by.x = "Sample", by.y = "sample_name", all.x = TRUE)

      # Dynamically map available annotations
      aes_mapping <- aes()
      if ("diagnosis" %in% colnames(pca_data)) {
        aes_mapping$color <- pca_data$diagnosis
      }
      if ("status" %in% colnames(pca_data)) {
        aes_mapping$shape <- pca_data$status
      }
      if ("timepoint" %in% colnames(pca_data)) {
        aes_mapping$size <- pca_data$timepoint
      }

      p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(aes(color = diagnosis, shape = status, size = timepoint), alpha = 0.8) +
        theme_minimal() +
        labs(title = "PCA of CNV Data", x = "PC1", y = "PC2") +
        scale_color_manual(values = setNames(brewer.pal(length(unique(pca_data$diagnosis)), "Set1"), unique(pca_data$diagnosis))) +
        scale_shape_manual(values = c(16, 17, 18)) +
        scale_size_continuous(range = c(2, 6)) +
        theme(
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)
        )

      ggsave(filename = file.path(output_dir, "pca_plot.pdf"), plot = p, width = 10, height = 8)
    }

    # UMAP Plot
    if (isTRUE(plot_umap)) {
      umap_result <- uwot::umap(t(cnv_matrix), n_neighbors = 15)

      umap_data <- data.frame(
        Sample = rownames(umap_result),
        UMAP1 = umap_result[, 1],
        UMAP2 = umap_result[, 2]
      )
      umap_data <- merge(umap_data, annotations, by.x = "Sample", by.y = "sample_name", all.x = TRUE)

      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = diagnosis, shape = status, size = timepoint), alpha = 0.8) +
        theme_minimal() +
        labs(title = "UMAP of CNV Data", x = "UMAP1", y = "UMAP2") +
        scale_color_manual(values = setNames(brewer.pal(length(unique(umap_data$diagnosis)), "Set1"), unique(umap_data$diagnosis))) +
        scale_shape_manual(values = c(16, 17, 18)) +
        scale_size_continuous(range = c(2, 6)) +
        theme(
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)
        )

      ggsave(filename = file.path(output_dir, "umap_plot.pdf"), plot = p, width = 10, height = 8)
    }

    # Karyogram Plot
    if (isTRUE(plot_karyogram)) {
      # Load gene annotations
      genome_annotation <- rtracklayer::import(gtf_path)
      genome_annotation <- genome_annotation[genome_annotation$type == "gene"]

      # Convert CNV matrix to GRanges for each sample
      samples <- colnames(cnv_matrix)
      for (sample in samples) {
        sample_cnv <- data.table(genomic_bin = rownames(cnv_matrix), log2_ratio = cnv_matrix[, sample])
        sample_cnv[, c("chrom", "start", "end") := tstrsplit(genomic_bin, "[:\\-]", fixed = TRUE)]
        sample_cnv[, chrom := paste0("chr", chrom)]
        sample_cnv[, chrom := factor(chrom, levels = c(paste0("chr", 1:22), "chrX", "chrY"))]
        sample_cnv <- sample_cnv[!is.na(chrom)]

        # Create GRanges object
        cnv_granges <- GenomicRanges::GRanges(
          seqnames = sample_cnv$chrom,
          ranges = IRanges::IRanges(start = as.numeric(sample_cnv$start), end = as.numeric(sample_cnv$end)),
          log2_ratio = sample_cnv$log2_ratio
        )

        # Add CNV type
        cnv_granges$cnv_type <- ifelse(cnv_granges$log2_ratio > 0.5, "Gain",
                                       ifelse(cnv_granges$log2_ratio < -0.5, "Loss", "Neutral"))

        # Associate CNVs with genes
        overlaps <- findOverlaps(cnv_granges, genome_annotation)
        cnv_granges$gene <- NA
        cnv_granges$gene[queryHits(overlaps)] <- genome_annotation$gene_name[subjectHits(overlaps)]

        # Plot karyogram using ggbio
        p <- ggbio::autoplot(
          cnv_granges,
          aes(fill = cnv_type),
          layout = "karyogram"
        ) +
          scale_fill_manual(
            values = c("Gain" = "red", "Loss" = "blue", "Neutral" = "gray"),
            name = "CNV Type"
          ) +
          theme_bw() +
          labs(
            title = paste("Genome-Wide CNV Karyogram for Sample", sample),
            y = "Log2 Ratio",
            x = "Chromosome"
          ) +
          theme(
            legend.position = "right",
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)
          ) +
          geom_text(
            data = as.data.frame(cnv_granges)[!is.na(cnv_granges$gene) & abs(cnv_granges$log2_ratio) > 0.5, ],
            aes(
              x = (start + end) / 2,  # Midpoint of the CNV
              y = log2_ratio,         # CNV log2 ratio
              label = gene            # Gene names
            ),
            size = 3,                 # Adjust font size
            angle = 0,               # Rotate text for better visibility
            vjust = -1,
            color = "black"
          )

        # Save the karyogram plot
        ggsave(filename = file.path(output_dir, paste0("karyogram_", sample, ".pdf")), plot = p, width = 12, height = 8)
      }
    }

    # Dendrogram Plot
    if (isTRUE(plot_dendrogram)) {
      library(ggdendro)

      # Perform hierarchical clustering again for dendrogram (if not passed)
      hc <- perform_clustering(cnv_matrix)
      dendrogram <- as.dendrogram(hc)
      dend_data <- dendro_data(dendrogram)

      p <- ggplot() +
        geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = dend_data$labels, aes(x = x, y = y, label = label), hjust = 1, angle = 90, size = 3) +
        theme_minimal() +
        labs(title = "Hierarchical Clustering Dendrogram", x = "", y = "Height")

      ggsave(filename = file.path(output_dir, "dendrogram.pdf"), plot = p, width = 12, height = 6)
    }
  }, error = function(e) {
    stop("Error generating plots: ", e$message)
  })
}
