# tests/testthat/test-clustering.R

context("Testing perform_clustering function")

test_that("perform_clustering correctly performs hierarchical clustering", {
  # Create dummy CNV matrix
  cnv_matrix <- matrix(
    c(0.5, 1.0,
      -0.5, -1.0),
    nrow = 2,
    dimnames = list(c("chr1:100-400", "chr1:500-900"), c("Sample1", "Sample2"))
  )

  # Perform clustering
  hc <- perform_clustering(cnv_matrix)

  # Check if hclust object is returned
  expect_s3_class(hc, "hclust")

  # Check if number of clusters is correct
  expect_equal(length(hc$labels), 2)
})
