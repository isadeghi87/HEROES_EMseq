# R/data_loading.R

#' Load CNV Segmentation Data from Metadata CSV
#'
#' This function reads a metadata CSV file containing paths to `.segs` files and associated sample information.
#' It combines all CNV data into a single `data.table` with appropriate annotations.
#' @importFrom data.table data.table fwrite
#' @param metadata_csv Path to the CSV file containing CNV metadata.
#'   The CSV should contain at least the following columns:
#'   - `segs_file`: Path to each `.segs` file.
#'   - `sample_name`: Unique identifier for each sample.
#'   - `diagnosis`: Diagnosis category (e.g., Tumor, Control).
#'   - Optional columns: `timepoint`, `status`, etc.
#' @return A `data.table` containing combined and annotated CNV data.
#' @export
#'
#' @examples
#' \dontrun{
#' metadata_csv <- "path/to/input_metadata.csv"
#' cnv_data <- load_cnv_data(metadata_csv)
#' }
load_cnv_data <- function(metadata_csv) {
  tryCatch({
    library(data.table)

    # Read the metadata CSV
    metadata <- fread(metadata_csv)

    # Check for required columns
    required_cols <- c("segs_file", "sample_name", "diagnosis")
    if (!all(required_cols %in% colnames(metadata))) {
      stop("Metadata CSV is missing one or more required columns: ", paste(required_cols, collapse = ", "))
    }

    # Load each .segs file and annotate with metadata
    cnv_data_list <- lapply(1:nrow(metadata), function(i) {
      seg_file <- metadata$segs_file[i]
      sample_name <- metadata$sample_name[i]
      diagnosis <- metadata$diagnosis[i]

      # Read the .segs file
      seg_data <- fread(seg_file)

      # Check for required columns in .segs file
      seg_required_cols <- c("sample", "chrom", "start", "end", "log2_ratio_median")
      if (!all(seg_required_cols %in% colnames(seg_data))) {
        stop("The .segs file ", seg_file, " is missing required columns: ", paste(seg_required_cols, collapse = ", "))
      }

      # Rename 'sample' column to 'sample_name' to avoid confusion
      setnames(seg_data, "sample", "sample_name")

      # Assign metadata to each row
      seg_data[, sample_name := sample_name]
      seg_data[, diagnosis := diagnosis]

      return(seg_data)
    })

    # Combine all data.tables
    combined_cnv_data <- rbindlist(cnv_data_list, use.names = TRUE, fill = TRUE)

    return(combined_cnv_data)
  }, error = function(e) {
    stop("Error loading CNV data: ", e$message)
  })
}
