# tests/testthat/test-load-cnv-data.R

context("Testing load_cnv_data function")

test_that("load_cnv_data correctly loads and annotates CNV data", {
  # Create temporary metadata CSV
  temp_dir <- tempdir()
  seg_file1 <- file.path(temp_dir, "sample1.segs")
  seg_file2 <- file.path(temp_dir, "sample2.segs")

  # Create dummy .segs files
  data_table1 <- data.table::data.table(sample = "Sample1", chrom = "chr1", start = 100, end = 500, log2_ratio_median = 0.5)
  data_table2 <- data.table::data.table(sample = "Sample2", chrom = "chr2", start = 200, end = 600, log2_ratio_median = -0.5)
  data.table::fwrite(data_table1, seg_file1, sep = "\t")
  data.table::fwrite(data_table2, seg_file2, sep = "\t")

  # Create metadata CSV
  metadata_csv <- file.path(temp_dir, "metadata.csv")
  metadata <- data.table::data.table(
    segs_file = c(seg_file1, seg_file2),
    sample_name = c("Sample1", "Sample2"),
    diagnosis = c("Tumor", "Control")
  )
  data.table::fwrite(metadata, metadata_csv, sep = ",")

  # Load CNV data
  cnv_data <- load_cnv_data(metadata_csv)

  # Assertions
  testthat::expect_s3_class(cnv_data, "data.table")
  testthat::expect_equal(nrow(cnv_data), 2)
  testthat::expect_true(all(c("sample_name", "diagnosis", "chrom", "start", "end", "log2_ratio_median") %in% colnames(cnv_data)))

  # Check annotations
  testthat::expect_equal(unique(cnv_data$sample_name), c("Sample1", "Sample2"))
  testthat::expect_equal(unique(cnv_data$diagnosis), c("Tumor", "Control"))
})
