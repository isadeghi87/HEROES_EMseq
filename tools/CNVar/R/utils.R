# R/utils.R

#' Create CNV Matrix from CNV Data
#'
#' This function converts a `data.table` of normalized CNV data into a matrix suitable for clustering and visualization.
#'
#' @param cnv_data A `data.table` containing normalized CNV data with at least the following columns:
#'   - `genomic_bin`: Unique identifier for each genomic bin.
#'   - `sample_name`: Unique identifier for each sample.
#'   - `log2_ratio_normalized`: Normalized log2 ratio value for the CNV segment.
#' @return A matrix where rows are genomic bins and columns are samples.
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_matrix <- create_cnv_matrix(cnv_normalized)
#' }
create_cnv_matrix <- function(cnv_data) {
  tryCatch({
    library(data.table)

    # Ensure 'genomic_bin' exists
    if (!"genomic_bin" %in% colnames(cnv_data)) {
      cnv_data[, genomic_bin := paste0(chrom, ":", start, "-", end)]
    }

    # Cast to wide format
    cnv_matrix_dt <- dcast(
      cnv_data,
      genomic_bin ~ sample_name,
      value.var = "log2_ratio_normalized",
      fun.aggregate = mean
    )

    # Set genomic_bin as row names and remove the column
    cnv_matrix_df <- as.data.frame(cnv_matrix_dt)
    rownames(cnv_matrix_df) <- cnv_matrix_df$genomic_bin
    cnv_matrix_df$genomic_bin <- NULL

    # Convert to matrix and handle missing values
    cnv_matrix <- as.matrix(cnv_matrix_df)
    cnv_matrix[is.na(cnv_matrix)] <- 0

    return(cnv_matrix)
  }, error = function(e) {
    stop("Error creating CNV matrix: ", e$message)
  })
}

#' Extract Sample IDs from Filenames
#'
#' This function extracts sample IDs from filenames by removing the `.segs` extension.
#'
#' @param filenames A character vector of file paths.
#' @return A vector of sample IDs.
#' @export
#'
#' @examples
#' \dontrun{
#' sample_ids <- extract_sample_ids(c("/path/to/sample1.segs", "/path/to/sample2.segs"))
#' }
extract_sample_ids <- function(filenames) {
  tryCatch({
    sample_ids <- gsub("\\.segs$", "", basename(filenames))
    return(sample_ids)
  }, error = function(e) {
    stop("Error extracting sample IDs: ", e$message)
  })
}
