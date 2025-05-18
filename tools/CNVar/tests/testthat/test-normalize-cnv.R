# tests/testthat/test-normalize-cnv.R

context("Testing normalize_cnv function")

test_that("normalize_cnv correctly normalizes log2_ratio_median", {
  # Create dummy CNV data
  cnv_data <- data.table(
    sample_name = c("Sample1", "Sample1", "Sample2", "Sample2"),
    chrom = c("chr1", "chr1", "chr2", "chr2"),
    start = c(100, 500, 200, 600),
    end = c(400, 900, 500, 1000),
    log2_ratio_median = c(0.5, 1.0, -0.5, -1.0)
  )

  # Normalize CNV data
  cnv_normalized <- normalize_cnv(cnv_data)

  # Check if normalization column exists
  expect_true("log2_ratio_normalized" %in% colnames(cnv_normalized))

  # Calculate expected normalized values
  sample1_median <- median(cnv_normalized[sample_name == "Sample1", log2_ratio_median], na.rm = TRUE)
  sample2_median <- median(cnv_normalized[sample_name == "Sample2", log2_ratio_median], na.rm = TRUE)

  # Check normalized values
  expect_equal(cnv_normalized[sample_name == "Sample1", log2_ratio_normalized], cnv_normalized[sample_name == "Sample1", log2_ratio_median] - sample1_median)
  expect_equal(cnv_normalized[sample_name == "Sample2", log2_ratio_normalized], cnv_normalized[sample_name == "Sample2", log2_ratio_median] - sample2_median)
})
