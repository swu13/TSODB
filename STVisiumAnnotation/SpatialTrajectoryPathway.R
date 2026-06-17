#!/usr/bin/env Rscript

# SpaCET Spatial Trajectory Pathway Analysis Pipeline
# Usage: Rscript SpatialTrajectoryPathway.R --outpath <output_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(SpaCET)
  library(SPATA2)
  library(stringr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), containing SpaCET_tumor.rds, spata_obj.rds, and seurat_gene_obj.rds", metavar="directory"),
  make_option(c("--trajectory_fdr"), type="double", default=0.05,
              help="FDR threshold for trajectory screening, default: 0.05", metavar="number"),
  make_option(c("--logfc_threshold"), type="double", default=0,
              help="Log fold change threshold for pathway enrichment, default: 0", metavar="number"),
  make_option(c("--pathway_p_val_adj"), type="double", default=0.05,
              help="Adjusted p-value threshold for pathway enrichment, default: 0.05", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="SpaCET spatial trajectory pathway analysis using gene set scores")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

#===============================================================================
# Function: SpatialTrajectoryPathway
# Description: Perform pathway-level spatial trajectory analysis and identify
#              enriched pathways across tumor, interface, and stroma regions
# Parameters:
#   outpath - Path to directory containing required RDS files
#   trajectory_fdr - FDR threshold for trajectory screening
#   logfc_threshold - Log fold change threshold for pathway enrichment
#   pathway_p_val_adj - Adjusted p-value threshold for pathway enrichment
#===============================================================================
SpatialTrajectoryPathway <- function(outpath,
                                     trajectory_fdr = 0.05,
                                     logfc_threshold = 0,
                                     pathway_p_val_adj = 0.05) {
  
  # Print start message
  message("\n========================================")
  message("Starting SpaCET Spatial Trajectory Pathway Analysis")
  message(paste("Output path:", outpath))
  message(paste("Trajectory FDR threshold:", trajectory_fdr))
  message(paste("Pathway logFC threshold:", logfc_threshold))
  message(paste("Pathway adj. p-value threshold:", pathway_p_val_adj))
  message("========================================\n")
  
  #-----------------------------------------------------------------------------
  # Step 1: Load SpaCET object and calculate gene set scores
  #-----------------------------------------------------------------------------
  message("[1/6] Loading SpaCET object and calculating gene set scores...")
  
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
  
  # Calculate gene set scores for Hallmark, CancerCellState, and TLS
  message("  - Calculating Hallmark gene set scores...")
  SpaCET_obj <- tryCatch({
    SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "Hallmark")
  }, error = function(e) {
    message(paste("    WARNING: Hallmark calculation failed:", e$message))
  })
  
  message("  - Calculating CancerCellState gene set scores...")
  SpaCET_obj <- tryCatch({
    SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "CancerCellState")
  }, error = function(e) {
    message(paste("    WARNING: CancerCellState calculation failed:", e$message))
  })
  
  message("  - Calculating TLS gene set scores...")
  SpaCET_obj <- tryCatch({
    SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "TLS")
  }, error = function(e) {
    message(paste("    WARNING: TLS calculation failed:", e$message))
  })
  
  # Extract gene set scores
  if (is.null(SpaCET_obj@results$GeneSetScore)) {
    stop("Gene set scores not found in SpaCET object")
  }
  
  genesetscore <- as.data.frame(t(SpaCET_obj@results$GeneSetScore))
  message(paste("  - Extracted gene set scores:", nrow(genesetscore), "spots,", ncol(genesetscore), "pathways"))
  
  # Map barcodes
  if (!is.null(SpaCET_obj@input$metaData)) {
    rownames(genesetscore) <- SpaCET_obj@input$metaData[rownames(genesetscore),,drop=FALSE]$barcode
  }
  
  #-----------------------------------------------------------------------------
  # Step 2: Load SPATA2 object and add pathway scores
  #-----------------------------------------------------------------------------
  message("\n[2/6] Loading SPATA2 object and adding pathway scores...")
  
  spata_file <- file.path(outpath, "spata_obj.rds")
  if (!file.exists(spata_file)) {
    stop(paste("SPATA2 object not found:", spata_file))
  }
  
  spata2Data <- tryCatch({
    readRDS(spata_file)
  }, error = function(e) {
    stop(paste("Failed to load SPATA2 object:", e$message))
  })
  message(paste("  - Loaded SPATA2 object from:", spata_file))
  
  # Subset gene set scores to SPATA2 barcodes
  genesetscore_spata <- genesetscore[getBarcodes(spata2Data), , drop = FALSE]
  message(paste("  - Subset pathway scores to", nrow(genesetscore_spata), "spots"))
  
  # Convert to tibble for SPATA2
  genesetscore_spata <- tibble::rownames_to_column(genesetscore_spata, var = "barcodes") %>% 
    tibble::as_tibble()
  
  # Add features to SPATA2 object
  spata2Data <- addFeatures(object = spata2Data,
                            feature_df = genesetscore_spata,
                            feature_names = colnames(genesetscore_spata)[-1],
                            overwrite = TRUE)
  message("  - Added pathway scores as features to SPATA2 object")
  
  #-----------------------------------------------------------------------------
  # Step 3: Screen pathways along trajectory
  #-----------------------------------------------------------------------------
  message("\n[3/6] Screening pathways along tumor-stroma trajectory...")
  
  pathway_names <- colnames(genesetscore_spata)[-1]
  message(paste("  - Screening", length(pathway_names), "pathways"))
  
  # Run trajectory screening
  sts_out <- tryCatch({
    spatialTrajectoryScreening(object = spata2Data,
                               id = "Tumor-Stroma",
                               variables = pathway_names)
  }, error = function(e) {
    stop(paste("Trajectory screening failed:", e$message))
  })
  
  message(paste("  - Found", nrow(sts_out@results$significance), "pathways"))
  
  # Extract significant pathways
  sign_df <- sts_out@results$significance %>%
    filter(fdr < trajectory_fdr)
  
  message(paste("  - Found", nrow(sign_df), "pathways"))
  
  if (nrow(sign_df) > 0) {
    message(paste("  - Found", nrow(sign_df), "significant pathways (FDR <", trajectory_fdr, ")"))
    
    # Get best model fits
    best_fits <- sts_out@results$model_fits %>%
      filter(variables %in% sign_df[["variables"]]) %>%
      group_by(variables) %>%
      slice_min(mae, n = 1)
    
    # Create output table
    output <- left_join(best_fits, sign_df, by = "variables") %>%
      select(geneset = variables,
             models,
             variations = rel_var,
             p_value,
             fdr) %>%
      as.data.frame()
  }else {
    message(paste("  - No significant pathways found with FDR <", trajectory_fdr))
    output <- data.frame(geneset = character(),models = character(),variations = numeric(),p_value = numeric(),fdr = numeric(),stringsAsFactors = FALSE)
  }
  # Save trajectory pathway results
  trajectory_file <- file.path(outpath, "TrajectoryPathways.txt")
  write.table(output, file = trajectory_file,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Trajectory pathways saved to:", trajectory_file))
  
  #-----------------------------------------------------------------------------
  # Step 4: Extract and save line plot data
  #-----------------------------------------------------------------------------
  message("\n[4/6] Extracting trajectory line plot data...")
  
  sts_df <- tryCatch({
    getStsDf(spata2Data, id = "Tumor-Stroma", variables = pathway_names)
  }, error = function(e) {
    message(paste("  - WARNING: Failed to get STS data frame:", e$message))
    return(NULL)
  })
  
  if (!is.null(sts_df)) {
    sts_df <- as.data.frame(sts_df)
    
    # Remove unnecessary columns
    if ("expr_est_idx" %in% colnames(sts_df)) {
      sts_df$expr_est_idx <- NULL
    }
    if ("bins_order" %in% colnames(sts_df)) {
      sts_df$bins_order <- NULL
    }
  } else {
    message("  - No line plot data generated")
    sts_df <- setNames(data.frame(matrix(ncol = length(pathway_names)+2, nrow = 0)),c("dist","dist_unit",pathway_names))
  }
  # Save line plot data
  lineplot_file <- file.path(outpath, "TrajectoryGenesetsLinePlot.txt")
  write.table(sts_df, file = lineplot_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Line plot data saved to:", lineplot_file))
  
  # Save updated SPATA2 object
  spata_pathway_file <- file.path(outpath, "spata_obj.pathway.rds")
  saveRDS(spata2Data, file = spata_pathway_file)
  message(paste("  - Updated SPATA2 object saved to:", spata_pathway_file))
  
  #-----------------------------------------------------------------------------
  # Step 5: Perform pathway enrichment analysis across regions
  #-----------------------------------------------------------------------------
  message("\n[5/6] Performing pathway enrichment analysis across regions...")
  
  # Load Seurat object
  seurat_file <- file.path(outpath, "seurat_gene_obj.rds")
  if (!file.exists(seurat_file)) {
    stop(paste("Seurat object not found:", seurat_file))
  }
  
  SeuratData <- tryCatch({
    readRDS(seurat_file)
  }, error = function(e) {
    stop(paste("Failed to load Seurat object:", e$message))
  })
  message(paste("  - Loaded Seurat object from:", seurat_file))
  
  # Subset pathway scores to Seurat cells
  genesetscore <- genesetscore[rownames(SeuratData@meta.data), , drop = FALSE]
  message(paste("  - Subset pathway scores to", nrow(genesetscore), "spots"))
  
  # Create pathway assay
  SeuratData[["pathways"]] <- CreateAssayObject(data = as.matrix(t(genesetscore)))
  DefaultAssay(SeuratData) <- "pathways"
  message("  - Created pathways assay in Seurat object")
  
  # Set identity and find markers
  Idents(SeuratData) <- "Interface"
  
  # Increase future global size limit for parallel processing
  options(future.globals.maxSize = 6 * 1024^3)
  
  message(paste("  - Finding differentially enriched pathways (logFC threshold:", logfc_threshold, ")"))
  DEGs <- tryCatch({
    FindAllMarkers(SeuratData, 
                   assay = "pathways", 
                   logfc.threshold = 0, 
                   only.pos = TRUE,
                   verbose = FALSE)
  }, error = function(e) {
    message(paste("  - WARNING: FindAllMarkers failed:", e$message))
    return(data.frame())
  })
  
  if (nrow(DEGs) > 0) {
    # Filter significant pathways
    SigDEGs <- DEGs %>%
      filter(p_val_adj < pathway_p_val_adj) %>%
      filter(avg_log2FC > logfc_threshold) %>%
      filter(pct.1 >= pct.2)    
  } else {
    message("  - No DEGs found")
    SigDEGs <- data.frame(p_val = numeric(),avg_log2FC = numeric(),pct.1 = numeric(),pct.2 = numeric(),p_val_adj = numeric(),cluster=character(),gene=character(),stringsAsFactors = FALSE)
  }
  # Save pathway enrichment results
  pathway_file <- file.path(outpath, "PathwayEnriched.txt")
  write.table(SigDEGs, file = pathway_file,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Pathway enrichment results saved to:", pathway_file))
  
  #output deconvoluted results
  predictions <- SpaCET_obj@results$deconvolution$propMat
  if(!(is.null(SpaCET_obj@input$metaData))){
    colnames(predictions) <- SpaCET_obj@input$metaData[colnames(predictions),,drop=FALSE]$barcode
  }
  predictions <- predictions[,colnames(SeuratData)]
  SeuratData[["deconvoluted"]] <- CreateAssayObject(data = predictions)
  
  # Save updated Seurat object
  seurat_pathway_file <- file.path(outpath, "seurat_pathway_obj.rds")
  saveRDS(SeuratData, file = seurat_pathway_file)
  message(paste("  - Updated Seurat object saved to:", seurat_pathway_file))
  message("\nAnalysis completed successfully!")
}

# Run the main function with error handling
tryCatch({
  SpatialTrajectoryPathway(
    outpath = opt$outpath,
    trajectory_fdr = opt$trajectory_fdr,
    logfc_threshold = opt$logfc_threshold,
    pathway_p_val_adj = opt$pathway_p_val_adj
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)