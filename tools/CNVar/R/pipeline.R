# R/pipeline.R

#' Run CNV Analysis Pipeline
#'
#' This function orchestrates the entire CNV analysis pipeline, including data loading, normalization,
#' matrix creation, clustering, statistical analysis, and plotting based on user-specified options.
#'
#' @param metadata_csv Path to the CSV file containing CNV metadata.
#' @param genome_version Genome version (`hg38` or `hg19`).
#' @param gtf_path Path to the GTF file for gene annotations.
#' @param plot_manhattan Logical. Whether to generate a Manhattan plot. Default is `TRUE`.
#' @param plot_heatmap_genomicBin Logical. Whether to generate a heatmap of genomic bins. Default is `TRUE`.
#' @param plot_pca Logical. Whether to generate a PCA plot. Default is `TRUE`.
#' @param plot_umap Logical. Whether to generate a UMAP plot. Default is `FALSE`.
#' @param plot_karyogram Logical. Whether to generate Karyogram plots. Default is `FALSE`.
#' @param plot_dendrogram Logical. Whether to generate a dendrogram plot. Default is `TRUE`.
#' @param output_dir Directory to save the results and plots. Default is `results/`.
#' @return NULL. Results and plots are saved to the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' run_pipeline(
#'   metadata_csv = "path/to/input_metadata.csv",
#'   genome_version = "hg38",
#'   gtf_path = "path/to/annotations.gtf",
#'   plot_manhattan = TRUE,
#'   plot_heatmap_genomicBin = TRUE,
#'   plot_pca = TRUE,
#'   plot_umap = FALSE,
#'   plot_karyogram = TRUE,
#'   plot_dendrogram = TRUE,
#'   output_dir = "results/"
#' )
#' }
run_pipeline <- function(metadata_csv,
                         genome_version,
                         gtf_path,
                         plot_manhattan = TRUE,
                         plot_heatmap_genomicBin = TRUE,
                         plot_pca = TRUE,
                         plot_umap = FALSE,
                         plot_karyogram = FALSE,
                         plot_dendrogram = TRUE,
                         output_dir = "results/") {
  tryCatch({
    library(data.table)

    # Validate genome_version
    if (!(genome_version %in% c("hg38", "hg19"))) {
      stop("Invalid genome_version. Choose 'hg38' or 'hg19'.")
    }

    # Load CNV data
    cnv_data <- load_cnv_data(metadata_csv)

    # Normalize CNV data
    cnv_normalized <- normalize_cnv(cnv_data)

    # Create CNV matrix
    cnv_matrix <- create_cnv_matrix(cnv_normalized)

    # Perform clustering
    hc_samples <- perform_clustering(cnv_matrix)

    # Prepare annotations
    # Extract all annotation columns except 'segs_file', 'sample_name', 'diagnosis', 'chrom', 'start', 'end', 'log2_ratio_median', 'log2_ratio_normalized', 'genomic_bin'
    annotation_cols <- setdiff(colnames(cnv_data), c("segs_file", "sample_name", "diagnosis", "chrom", "start", "end", "log2_ratio_median", "log2_ratio_normalized", "genomic_bin"))
    if (length(annotation_cols) > 0) {
      annotations <- unique(cnv_data[, ..annotation_cols])
    } else {
      annotations <- data.frame(sample_name = unique(cnv_data$sample_name), diagnosis = unique(cnv_data$diagnosis), stringsAsFactors = FALSE)
    }

    # Ensure all samples in cnv_matrix are present in annotations
    annotations <- annotations[sample_name %in% colnames(cnv_matrix), ]

    # Check if 'timepoint' is present and numeric if plot_umap or any plot needs it
    if ((isTRUE(plot_umap)) || any(c(plot_pca, plot_karyogram) & "timepoint" %in% colnames(annotations))) {
      if ("timepoint" %in% colnames(annotations)) {
        if (!is.numeric(annotations$timepoint)) {
          stop("The 'timepoint' column must be numeric for plotting purposes.")
        }
      } else {
        stop("Plotting timepoints requires a 'timepoint' column in the metadata CSV.")
      }
    }

    # Generate plots based on user-specified options
    generate_plots(
      cnv_matrix = cnv_matrix,
      annotations = annotations,
      genome_version = genome_version,
      gtf_path = gtf_path,
      plot_manhattan = plot_manhattan,
      plot_heatmap_genomicBin = plot_heatmap_genomicBin,
      plot_pca = plot_pca,
      plot_umap = plot_umap,
      plot_karyogram = plot_karyogram,
      plot_dendrogram = plot_dendrogram,
      output_dir = output_dir
    )

    # Perform statistical analysis
    diff_cnv_results <- differential_cnv_analysis(cnv_normalized)

    # Save statistical results
    fwrite(diff_cnv_results, file = file.path(output_dir, "differential_cnv_analysis.csv"))

    message("CNV analysis pipeline completed successfully. Results are saved in ", output_dir)

  }, error = function(e) {
    stop("Error running CNV analysis pipeline: ", e$message)
  })
}
