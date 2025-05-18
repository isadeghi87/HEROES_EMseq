setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/results/methylation_array/")

## packages 
library(rtracklayer)
library(dplyr)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Replace "path/to/your/idat/files" with the actual path where your .idat files are located
RGset_450k <- read.metharray.exp(base = "~/Downloads/idats/450/")
RGset_EPIC <- read.metharray.exp(base = "~/Downloads/idats/EPIC/",force = T)
RGset_unknown <- read.metharray.exp(base = "~/Downloads/idats/unknown/",force = T)

# Create a list of your RGset objects and label them accordingly.
RGsets <- list(
  "450k"    = RGset_450k,
  "EPIC"    = RGset_EPIC,
  "unknown" = RGset_unknown
)

# Optionally, define the corresponding annotation functions/packages.
# Note: For "unknown", you may not have a valid annotation.
getAnnot <- function(platform) {
  if(platform == "450k") {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    if(!"group_name" %in% colnames(annot)) annot$group_name <- annot$Name
    return(annot)
  } else if(platform == "EPIC") {
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    if(!"group_name" %in% colnames(annot)) annot$group_name <- annot$Name
    return(annot)
  } else {
    return(NULL)  # No annotation available for unknown platform.
  }
}

# Define the liftover chain file path
chain_file <- "/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/liftover/hg19ToHg38.over.chain"
chain <- import.chain(chain_file)

results_list <- list()

for(platform in names(RGsets)) {
  message("Processing platform: ", platform)
  RGset <- RGsets[[platform]]
  
  # Preprocess and extract beta values.
  Mset <- preprocessRaw(RGset)
  beta_values <- as.data.frame(getBeta(Mset))
  beta_values$group_name <- rownames(beta_values)
  
  # Get annotation if available.
  annot <- getAnnot(platform)
  
  if(!is.null(annot)) {
    # Liftover annotation from hg19 to hg38.
    gr_hg19 <- makeGRangesFromDataFrame(annot,
                                        seqnames.field = "chr",
                                        start.field = "pos",
                                        end.field = "pos",
                                        keep.extra.columns = TRUE)
    gr_hg38 <- unlist(liftOver(gr_hg19, chain))
    df_hg38 <- as.data.frame(gr_hg38)
    if(!"group_name" %in% colnames(df_hg38) && "Name" %in% colnames(df_hg38))
      df_hg38$group_name <- df_hg38$Name
    
    # Merge beta values with lifted annotation.
    merged_data <- merge(beta_values, df_hg38, by = "group_name")
    
    # For demonstration, use the first sample column for filtering.
    sample_col <- setdiff(colnames(beta_values), "group_name")[1]
    filtered_data <- merged_data %>%
      filter( (!!as.symbol(sample_col)) > 0.3,
              (!!as.symbol(sample_col)) < 0.7 )
    
    # Create BED coordinates (convert to 0-based start).
    filtered_data <- filtered_data %>%
      mutate(bed_start = start - 1,
             bed_end   = start)
    
    bed_data <- filtered_data %>%
      select(seqnames, bed_start, bed_end, group_name, !!as.symbol(sample_col), strand) %>%
      rename(score = !!as.symbol(sample_col))
    
    bed_data <- na.omit(bed_data)
    
    # Save output files with platform-specific filenames.
    bed_filename <- paste0("./", platform, "_filtered_bulk_methylation.bed")
    write.table(bed_data, file = bed_filename,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Save a BEDGRAPH file (chrom, start, end, score)
    bedgraph_data <- bed_data %>% select(seqnames, bed_start, bed_end, score)
    bedgraph_filename <- paste0("./", platform, "_filtered_bulk_methylation.bedgraph")
    write.table(bedgraph_data, file = bedgraph_filename,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Store results for later inspection.
    results_list[[platform]] <- list(beta_values = beta_values,
                                     merged_data = merged_data,
                                     bed_data = bed_data)
  } else {
    message("No annotation available for platform ", platform, ". Only beta values are saved.")
    results_list[[platform]] <- list(beta_values = beta_values)
    # Optionally, save the beta values to a CSV file.
    write.csv(beta_values,
              file = paste0("path/to/output/", platform, "_beta_values.csv"),
              row.names = TRUE)
  }
}
