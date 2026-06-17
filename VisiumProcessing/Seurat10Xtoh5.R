#!/usr/bin/env Rscript

# 10X Genomics H5 File Conversion Pipeline
# Convert sparse matrix counts to 10X HDF5 format for downstream analysis
# Usage: Rscript convert_to_10x_h5.R --input_dir <input_directory> --output_file <output_h5_file> [--species <species>]

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(DropletUtils)
  library(Seurat)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL,
              help="Input directory path containing matrix.mtx, genes.tsv, barcodes.tsv files (required)", metavar="directory"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="Output H5 file path (required), e.g., filtered_feature_bc_matrix.h5", metavar="file"),
  make_option(c("-s", "--species"), type="character", default="GRCh38",
              help="Species genome reference: GRCh38 (human) or mm10 (mouse), default: GRCh38", metavar="string"),
  make_option(c("-c", "--chemistry"), type="character", default="Single Cell 3' v3",
              help="10X chemistry version, default: Single Cell 3' v3", metavar="string"),
  make_option(c("-g", "--gene_column"), type="integer", default=2,
              help="Gene column to use in genes.tsv (1 for gene ID, 2 for gene symbol), default: 2", metavar="number"),
  make_option(c("--version"), type="character", default="3",
              help="H5 version format, default: 3", metavar="string"),
  make_option(c("--library_id"), type="character", default="slice1",
              help="Library ID for the dataset, default: slice1", metavar="string")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="Convert sparse matrix counts to 10X Genomics HDF5 format")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("ERROR: --input_dir parameter is required", call.=FALSE)
}

if (is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("ERROR: --output_file parameter is required", call.=FALSE)
}

#===============================================================================
# Function: convert_to_10x_h5
# Description: Read 10X sparse matrix files and convert to HDF5 format
# Parameters:
#   input_dir - Directory containing matrix.mtx, genes.tsv, barcodes.tsv
#   output_file - Output H5 file path
#   species - Species genome reference (GRCh38 or mm10)
#   chemistry - 10X chemistry version
#   gene_column - Column to use from genes.tsv (1=gene ID, 2=gene symbol)
#   version - H5 version format
#   library_id - Library identifier
#===============================================================================
convert_to_10x_h5 <- function(input_dir, 
                              output_file,
                              species = "GRCh38",
                              chemistry = "Single Cell 3' v3",
                              gene_column = 2,
                              version = "3",
                              library_id = "slice1") {
  
  # Print start message
  message("\n========================================")
  message("Starting 10X H5 File Conversion")
  message(paste("Input directory:", input_dir))
  message(paste("Output file:", output_file))
  message(paste("Species:", species))
  message(paste("Chemistry:", chemistry))
  message(paste("Gene column:", gene_column))
  message("========================================\n")
  
  #-----------------------------------------------------------------------------
  # Step 1: Check if input directory exists
  #-----------------------------------------------------------------------------
  message("[1/4] Checking input files...")
  
  if (!dir.exists(input_dir)) {
    stop(paste("Input directory does not exist:", input_dir))
  }
  
  # Check for required files
  required_files <- c("matrix.mtx.gz", "barcodes.tsv.gz")
  missing_files <- c()
  
  for (file in required_files) {
    if (!file.exists(file.path(input_dir, file))) {
      missing_files <- c(missing_files, file)
    }
  }
  
  # Also check for features.tsv (newer 10X format)
  if (file.exists(file.path(input_dir, "features.tsv.gz"))) {
    message("  - Found features.tsv.gz (using newer 10X format)")
    genes_file <- "features.tsv.gz"
  } else if (file.exists(file.path(input_dir, "genes.tsv.gz"))) {
    message("  - Found genes.tsv.gz (using older 10X format)")
    genes_file <- "genes.tsv.gz"
  } else {
    missing_files <- c(missing_files, "genes.tsv.gz or features.tsv.gz")
  }
  
  if (length(missing_files) > 0) {
    stop(paste("Missing required files:", paste(missing_files, collapse = ", ")))
  }
  
  message("  - All required files found")
  
  #-----------------------------------------------------------------------------
  # Step 2: Read 10X counts data
  #-----------------------------------------------------------------------------
  message("\n[2/4] Reading 10X counts data...")
  
  counts <- tryCatch({
    Read10X(input_dir, gene.column = gene_column)
  }, error = function(e) {
    stop(paste("Failed to read 10X data:", e$message))
  })
  
  #-----------------------------------------------------------------------------
  # Step 3: Extract gene_type (if available) and prepare data
  #-----------------------------------------------------------------------------
  message("\n[3/4] Preparing gene metadata...")
  
  # Initialize gene_type with a default value
  gene_type <- rep("Gene Expression", nrow(counts))
  
  # If using the newer features.tsv.gz format, try to extract the actual 3rd column
  if (genes_file == "features.tsv.gz") {
    # Use data.table or read.delim to read the gz file
    if (requireNamespace("data.table", quietly = TRUE)) {
      features_df <- data.table::fread(file.path(input_dir, genes_file), header = FALSE)
    } else {
      features_df <- read.delim(file.path(input_dir, genes_file), header = FALSE, stringsAsFactors = FALSE)
    }
    
    # Ensure we have at least 3 columns and match the rows with the counts matrix
    if (ncol(features_df) >= 3 && nrow(features_df) == nrow(counts)) {
      gene_type <- as.character(features_df[[3]])
      message("  - Successfully extracted gene types from features.tsv.gz")
    } else {
      message("  - Warning: features.tsv.gz column/row mismatch. Using default gene type 'Gene Expression'.")
    }
  } else {
    message("  - Using older genes.tsv format. Setting default gene type to 'Gene Expression'.")
  }

  #-----------------------------------------------------------------------------
  # Step 3: Write H5 file
  #-----------------------------------------------------------------------------
  message("\n[4/4] Writing H5 file...")
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("  - Created output directory:", output_dir))
  }
  
  # Write to HDF5 format
  tryCatch({
    write10xCounts(
      output_file,
      counts,
      barcodes = colnames(counts),
      gene.id = rownames(counts),
      gene.symbol = rownames(counts),
      gene.type = gene_type,
      overwrite = TRUE,
      type = "HDF5",
      genome = species,
      version = version,
      chemistry = chemistry,
      original.gem.groups = 1L,
      library.ids = library_id
    )
    message(paste("  - Successfully wrote H5 file to:", output_file))
  }, error = function(e) {
    stop(paste("Failed to write H5 file:", e$message))
  })
  
}

#===============================================================================
# Execute main function
#===============================================================================

# Run the main function with error handling
tryCatch({
  convert_to_10x_h5(
    input_dir = opt$input_dir,
    output_file = opt$output_file,
    species = opt$species,
    chemistry = opt$chemistry,
    gene_column = opt$gene_column,
    version = opt$version,
    library_id = opt$library_id
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)