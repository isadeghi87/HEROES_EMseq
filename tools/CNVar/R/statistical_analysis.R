# R/statistical_analysis.R

#' Perform Differential CNV Analysis
#'
#' This function performs differential CNV analysis between two groups based on the `diagnosis` annotation.
#' It calculates p-values and adjusted p-values for each genomic bin.
#'
#' @param cnv_data A `data.table` containing CNV data with at least the following columns:
#'   - `genomic_bin`: Unique identifier for each genomic bin.
#'   - `log2_ratio_median`: The median log2 ratio value for the CNV segment.
#'   - `diagnosis`: Diagnosis category (e.g., Tumor, Control).
#' @return A `data.table` with differential CNV results including p-values and adjusted p-values.
#' @export
#'
#' @examples
#' \dontrun{
#' diff_cnv_results <- differential_cnv_analysis(cnv_normalized)
#' }
differential_cnv_analysis <- function(cnv_data) {
  tryCatch({
    library(data.table)

    # Ensure 'diagnosis' is a factor
    cnv_data[, diagnosis := as.factor(diagnosis)]

    # Perform t-test for each genomic bin if there are exactly two groups
    results <- cnv_data[, {
      if (length(unique(diagnosis)) == 2) {
        test_result <- tryCatch(t.test(log2_ratio_median ~ diagnosis), error = function(e) NULL)
        if (!is.null(test_result)) {
          list(
            p_value = test_result$p.value,
            mean_diff = diff(tapply(log2_ratio_median, diagnosis, mean, na.rm = TRUE))
          )
        } else {
          list(
            p_value = NA_real_,
            mean_diff = NA_real_
          )
        }
      } else {
        list(
          p_value = NA_real_,
          mean_diff = NA_real_
        )
      }
    }, by = genomic_bin]

    # Adjust p-values for multiple testing using Benjamini-Hochberg
    results[, adjusted_p_value := p.adjust(p_value, method = "BH")]

    return(results)
  }, error = function(e) {
    stop("Error performing differential analysis: ", e$message)
  })
}
