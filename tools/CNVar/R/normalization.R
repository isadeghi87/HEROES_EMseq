# R/normalization.R

#' Normalize CNV Data
#'
#' This function performs median centering normalization for each sample in the CNV data.
#'
#' @param cnv_data A `data.table` containing CNV data with at least the following columns:
#'   - `sample_name`: Unique identifier for each sample.
#'   - `log2_ratio_median`: The median log2 ratio value for the CNV segment.
#' @return A normalized CNV `data.table` with an additional column `log2_ratio_normalized`.
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_normalized <- normalize_cnv(cnv_data)
#' }
normalize_cnv <- function(cnv_data) {
  tryCatch({
    library(data.table)

    # Median centering normalization for each sample
    cnv_data[, log2_ratio_normalized := log2_ratio_median - median(log2_ratio_median, na.rm = TRUE), by = sample_name]

    return(cnv_data)
  }, error = function(e) {
    stop("Error normalizing CNV data: ", e$message)
  })
}
