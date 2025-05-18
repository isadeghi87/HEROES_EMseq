# R/clustering.R

#' Perform Hierarchical Clustering
#'
#' This function performs hierarchical clustering on the samples based on the normalized CNV matrix.
#'
#' @param cnv_matrix A matrix of normalized log2 ratio values. Rows represent genomic bins, and columns represent samples.
#' @return A hierarchical clustering object of class `hclust`.
#' @export
#'
#' @examples
#' \dontrun{
#' hc_samples <- perform_clustering(cnv_matrix)
#' }
perform_clustering <- function(cnv_matrix) {
  tryCatch({
    library(stats)

    # Compute distance matrix using Euclidean distance
    distance_matrix <- dist(t(cnv_matrix), method = "euclidean")

    # Perform hierarchical clustering using Ward's method
    hc <- hclust(distance_matrix, method = "ward.D2")

    return(hc)
  }, error = function(e) {
    stop("Error performing clustering: ", e$message)
  })
}
