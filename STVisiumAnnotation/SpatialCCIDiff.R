#!/usr/bin/env Rscript

# Stlearn CCI Enrichment Analysis Pipeline
# Usage: Rscript SpatialCCIDiff.R --outpath <output_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), containing SpaCET_tumor.rds and CCI files", metavar="directory"),
  make_option(c("--logfc_threshold"), type="double", default=0,
              help="Log fold change threshold for marker detection, default: 0", metavar="number"),
  make_option(c("--p_adj_threshold"), type="double", default=0.05,
              help="Adjusted p-value threshold for significance, default: 0.05", metavar="number"),
  make_option(c("--only_pos"), type="logical", default=TRUE,
              help="Only return positive markers, default: TRUE", metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="Stlearn cell-cell interaction enrichment analysis")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

#===============================================================================
# Function: pad_matrix
# Description: Pad matrix with missing columns to match reference column set
# Parameters:
#   mat - Input matrix
#   all_lr - Reference vector of all ligand-receptor pair names
# Returns:
#   Padded matrix with all columns from all_lr
#===============================================================================
pad_matrix <- function(mat, all_lr) {
  # Find missing ligand-receptor pairs
  missing_lr <- setdiff(all_lr, colnames(mat))
  
  # Add missing columns as zeros
  if (length(missing_lr) > 0) {
    pad <- matrix(0, nrow = nrow(mat), ncol = length(missing_lr))
    colnames(pad) <- missing_lr
    rownames(pad) <- rownames(mat)
    
    mat <- cbind(mat, pad)
  }
  
  # Reorder columns to match reference order
  mat <- mat[, all_lr, drop = FALSE]
  
  return(mat)
}

safe_read_cci <- function(file_path, all_lr) {

  if (!file.exists(file_path)) {
    message(paste("  - Missing file:", file_path, "-> creating empty matrix"))
    return(matrix(0, nrow = 0, ncol = length(all_lr),
                  dimnames = list(NULL, all_lr)))
  }

  mat <- tryCatch({
    read.csv(file_path, row.names = 1, check.names = FALSE)
  }, error = function(e) {
    message(paste("  - Failed reading:", file_path, "-> returning empty"))
    return(NULL)
  })

  if (is.null(mat) || nrow(mat) == 0) {
    message(paste("  - Empty file:", file_path, "-> creating empty matrix"))
    return(matrix(0, nrow = 0, ncol = length(all_lr),
                  dimnames = list(NULL, all_lr)))
  }

  return(mat)
}

#===============================================================================
# Function: SpatialCCIDiff
# Description: Perform analysis to identify enriched 
#              cell-cell interactions across tumor, interface, and stroma regions
# Parameters:
#   outpath - Path to directory containing SpaCET_tumor.rds and CCI files
#   logfc_threshold - Log fold change threshold for marker detection
#   p_adj_threshold - Adjusted p-value threshold for significance
#   only_pos - Only return positive markers
#===============================================================================
SpatialCCIDiff <- function(outpath, 
                              logfc_threshold = 0, 
                              p_adj_threshold = 0.05,
                              only_pos = TRUE) {
  
  # Print start message
  message("\n========================================")
  message("Starting SpaCET Spatial Trajectory Analysis")
  message(paste("Output path:", outpath))
  message("========================================\n")
  
  #-----------------------------------------------------------------------------
  # Step 1: Load SpaCET object and extract region annotations
  #-----------------------------------------------------------------------------
  message("[1/7] Loading SpaCET object...")
  
  input_rds <- file.path(outpath, "SpaCET_tumor.rds")
  if (!file.exists(input_rds)) {
          stop(paste("Input file does not exist. Expected:", 
                 file.path(outpath, "SpaCET_tumor.rds")))
  }
  
  SpaCET_obj <- tryCatch({
    readRDS(input_rds)
  }, error = function(e) {
    stop(paste("Failed to load SpaCET object:", e$message))
  })
  message(paste("  - Loaded from:", input_rds))
  
  # Extract interface regions
  if (is.null(SpaCET_obj@results$CCI$interface)) {
    stop("CCI interface data not found in SpaCET object. Please run SpaCET.identify.interface first.")
  }
  
  regions <- as.data.frame(t(SpaCET_obj@results$CCI$interface))
  if(!(is.null(SpaCET_obj@input$metaData))){
    rownames(regions) <- SpaCET_obj@input$metaData[rownames(regions),,drop=FALSE]$barcode
  }
  
  #-----------------------------------------------------------------------------
  # Step 2: Load and prepare Seurat object with spatial data
  #-----------------------------------------------------------------------------
  message("\n[2/7] Loading and preparing Seurat object...")
  
  # Load 10X Spatial data
  SeuratData <- CreateSeuratObject(counts = SpaCET_obj@input$counts,assay = "Spatial")
#  SeuratData <- tryCatch({
#    Load10X_Spatial(data.dir = rawpath)
#  }, error = function(e) {
#    stop(paste("Failed to load 10X Spatial data:", e$message))
#  })
#  message(paste("  - Loaded Seurat object with", ncol(SeuratData), "spots and", nrow(SeuratData), "genes"))
  
  # Subset to spots with region annotations
  SeuratData <- subset(SeuratData, cells = rownames(regions))
  message(paste("  - Subset to", ncol(SeuratData), "spots with region annotations"))
  
  # Add interface annotations to Seurat metadata
  SeuratData$Interface <- regions[rownames(SeuratData@meta.data), "Interface", drop = FALSE]
  message("  - Added Interface annotations to Seurat metadata")
  
  #-----------------------------------------------------------------------------
  # Step 3: Load cell-cell interaction scores
  #-----------------------------------------------------------------------------
  message("\n[3/7] Loading cell-cell interaction scores...")
  
  # Define CCI file paths
  tumor_file <- file.path(outpath, "lr_interaction_score.Tumor.csv")
  interface_file <- file.path(outpath, "lr_interaction_score.Interface.csv")
  stroma_file <- file.path(outpath, "lr_interaction_score.Stroma.csv")
  
  # Check if files exist
  if (!file.exists(tumor_file)) {
    stop(paste("Tumor CCI file not found:", tumor_file))
  }
  if (!file.exists(interface_file)) {
    stop(paste("Interface CCI file not found:", interface_file))
  }
  if (!file.exists(stroma_file)) {
    stop(paste("Stroma CCI file not found:", stroma_file))
  }
  
  # Load CCI scores
  TumorCCI <- safe_read_cci(tumor_file, NULL)
  InterfaceCCI <- safe_read_cci(interface_file, NULL)
  StromaCCI <- safe_read_cci(stroma_file, NULL)
  
  message(paste("  - Loaded Tumor CCI:", nrow(TumorCCI), "spots,", ncol(TumorCCI), "LR pairs"))
  message(paste("  - Loaded Interface CCI:", nrow(InterfaceCCI), "spots,", ncol(InterfaceCCI), "LR pairs"))
  message(paste("  - Loaded Stroma CCI:", nrow(StromaCCI), "spots,", ncol(StromaCCI), "LR pairs"))
  
  #-----------------------------------------------------------------------------
  # Step 4: Pad matrices to have consistent LR pairs
  #-----------------------------------------------------------------------------
  message("\n[4/7] Padding matrices to consistent LR pairs...")
  
  # Get union of all LR pairs
  all_lr <- union(colnames(TumorCCI), union(colnames(InterfaceCCI), colnames(StromaCCI)))
  message(paste("  - Total unique LR pairs:", length(all_lr)))
  
  # Pad matrices
  TumorCCI2 <- pad_matrix(TumorCCI, all_lr)
  InterfaceCCI2 <- pad_matrix(InterfaceCCI, all_lr)
  StromaCCI2 <- pad_matrix(StromaCCI, all_lr)
  
  #-----------------------------------------------------------------------------
  # Step 5: Combine all CCI scores and create Seurat assay
  #-----------------------------------------------------------------------------
  message("\n[5/7] Creating CCI assay in Seurat object...")
  
  # Combine all regions
  CCI_all <- rbind(TumorCCI2, InterfaceCCI2, StromaCCI2)
  
  # Transpose and reorder to match Seurat cells
  CCI_mat <- t(CCI_all)
  CCI_mat <- CCI_mat[, rownames(SeuratData@meta.data), drop = FALSE]
  
  # Create new assay
  SeuratData[["CCI"]] <- CreateAssayObject(data = CCI_mat)
  DefaultAssay(SeuratData) <- "CCI"
  message("  - Created CCI assay in Seurat object")
  
  #-----------------------------------------------------------------------------
  # Step 6: Find differentially enriched LR pairs across regions
  #-----------------------------------------------------------------------------
  message("\n[6/7] Finding differentially enriched LR pairs...")
  
  # Set identity to Interface regions
  Idents(SeuratData) <- "Interface"
  
  # Find markers for all regions
  message(paste("  - Finding markers with logFC threshold:", logfc_threshold))
  message(paste("  - Only positive markers:", only_pos))
  
  DEGs <- tryCatch({
    FindAllMarkers(SeuratData, 
                   assay = "CCI", 
                   logfc.threshold = logfc_threshold, 
                   only.pos = only_pos,
                   verbose = FALSE)
  }, error = function(e) {
    stop(paste("Failed to find markers:", e$message))
  })
  
  message(paste("  - Found", nrow(DEGs), "potential markers"))
  
  # Filter significant DEGs
  SigDEGs <- subset(DEGs,p_val_adj < p_adj_threshold & avg_log2FC > logfc_threshold & pct.1 >= pct.2)
  
  message(paste("  - Significant DEGs after filtering:", nrow(SigDEGs)))  
  #-----------------------------------------------------------------------------
  # Step 7: Save results
  #-----------------------------------------------------------------------------
  message("\n[7/7] Saving results...")
  
  # Save significant DEGs
  output_file <- file.path(outpath, "CCIEnriched.txt")
  write.table(SigDEGs, file = output_file, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Significant DEGs saved to:", output_file))
  
  # Save Seurat object
  seurat_file <- file.path(outpath, "seurat_CCI_obj.rds")
  saveRDS(SeuratData, file = seurat_file)
  message(paste("  - Seurat object saved to:", seurat_file))
  
  message("\nAnalysis completed successfully!")

}

#===============================================================================
# Execute main function
#===============================================================================

# Run the main function with error handling
tryCatch({
  SpatialCCIDiff(
    outpath = opt$outpath,
    logfc_threshold = opt$logfc_threshold,
    p_adj_threshold = opt$p_adj_threshold,
    only_pos = opt$only_pos
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)