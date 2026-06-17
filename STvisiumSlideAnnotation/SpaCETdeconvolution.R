#!/usr/bin/env Rscript

# SpaCET Deconvolution Pipeline
# Usage: Rscript SpaCET_deconvolution.R --outpath <output_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SpaCET)
  source("SpaCET.deconvolution.malignant.revised.R")
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), such as Output/BRCA.Brain.Visium.human.GSM5420753", metavar="directory"),
  make_option(c("-c", "--cancerType"), type="character", default=NULL,
              help="Cancer type (e.g., BRCA, LUCA, etc.) (required), such as the name in TCGA ", metavar="string"),
  make_option(c("-d", "--rawdata_dir"), type="character", default="NULL",
              help="Raw data directory (required), contained h5 and spatial folder", metavar="directory"),
  make_option(c("-n", "--cores"), type="integer", default=5,
              help="Number of CPU cores to use, default is 5", metavar="number"),
  make_option(c("-g", "--min_genes"), type="integer", default=100,
              help="Minimum number of genes per spot, default is 100", metavar="number"),
  make_option(c("-p", "--point_size"), type="double", default=2,
              help="Point size in main plot, default is 2", metavar="number"),
  make_option(c("-s", "--small_point_size"), type="double", default=0.1,
              help="Point size in individual cell type plot, default is 0.1", metavar="number"),
  make_option(c("--width"), type="integer", default=8,
              help="Main plot width (inches), default is 8", metavar="number"),
  make_option(c("--height"), type="integer", default=8,
              help="Main plot height (inches), default is 8", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="SpaCET spatial transcriptomics deconvolution analysis script")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

if (is.null(opt$cancerType)) {
  print_help(opt_parser)
  stop("ERROR: --cancerType parameter is required", call.=FALSE)
}

if (is.null(opt$rawdata_dir)) {
  print_help(opt_parser)
  stop("ERROR: --rawdata_dir parameter is required", call.=FALSE)
}

# Define main function
SpaCETdeconvolution <- function(outpath, cancerType, rawdata_dir,
                                coreNo = 5, min_genes = 100, 
                                pointSize = 2, small_pointSize = 0.1,
                                width = 8, height = 8) {
  
  
  # Check if raw data directory exists
  if (!dir.exists(rawdata_dir)) {
    stop(paste("Raw data directory does not exist:", rawdata_dir))
  }

  message(paste("Cancer type:", cancerType))
  message(paste("Raw data path:", rawdata_dir))
  
  # Create SpaCET object
  message("Creating SpaCET object...")
  SpaCET_obj <- tryCatch({
    create.SpaCET.object.10X(visiumPath = rawdata_dir)
  }, error = function(e1) {
    message( "initiateSpataObjectVisium failed: ", e1$message, "\nTrying initiateSpataObjectSlideseq...")
    tryCatch({
      expr_file <- list.files(rawdata_dir, pattern="(expression|counts).*\\.csv$", full.names=TRUE, recursive=TRUE)
      if(length(expr_file) == 0) stop("No expression matrix file found in cosmxPath.")
      expr_data <- read.csv(expr_file[1], row.names=1, check.names=FALSE)
      meta_file <- list.files(rawdata_dir, pattern="metadata_file\\.csv$|metadata\\.csv$", full.names=TRUE, recursive=TRUE)
      if(length(meta_file) == 0) stop("No metadata file found in cosmxPath.")
      spotCoordinates <- read.csv(meta_file[1], row.names=1, check.names=FALSE)
      common_cells <- intersect(colnames(expr_data), rownames(spotCoordinates))
      if(length(common_cells)==0){
        st.matrix.data <- Matrix::Matrix(t(as.matrix(expr_data)), sparse=TRUE)
      } else {
        st.matrix.data <- Matrix::Matrix(as.matrix(expr_data), sparse=TRUE)
      }
      common_cells <- intersect(colnames(st.matrix.data), rownames(spotCoordinates))
      if(length(common_cells) == 0) stop("No matching cell IDs between expression data and coordinates.")
      st.matrix.data <- st.matrix.data[, common_cells]
      spotCoordinates <- spotCoordinates[common_cells, ]
      create.SpaCET.object(counts=st.matrix.data,spotCoordinates=spotCoordinates,platform="CosMx")
    }, error = function(e2) {
      stop(paste0("Failed to create SPATA2 object.\n","Visium error: ", e1$message, "\n","VisiumHD error: ", e2$message))
    })
  })
  SpaCET_obj@input$platform="Visium"
  
  # Quality control
  message("Performing quality control...")
  SpaCET_obj <- tryCatch({
    SpaCET.quality.control(SpaCET_obj, min.genes = min_genes)
  }, error = function(e) {
    stop(paste("Quality control failed:", e$message))
  })
  
  # Deconvolution analysis
  message("Performing deconvolution analysis...")
  SpaCET_obj <- tryCatch({
    SpaCET.deconvolution(SpaCET_obj, cancerType = cancerType, coreNo = coreNo)
  }, error = function(e) {
    stop(paste("Deconvolution analysis failed:", e$message))
  })
  
  # Create output directory
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  rds_path <- file.path(outpath, "SpaCET_obj.rds")
  saveRDS(SpaCET_obj, file = rds_path)
  message(paste("Results saved to:", rds_path))
  
  # Define color vector
  colors_vector <- c(
    "Malignant" = "#f3c300", "CAF" = "#be0032", "Endothelial" = "#0067a5",
    "Unidentifiable" = "grey82", "B cell" = "#008856", "cDC" = "#782AB6",
    "pDC" = "#2b3d26", "Macrophage" = "#f6a600", "Neutrophil" = "#848482",
    "NK" = "#a1caf1", "Mast" = "#FA0087", "T CD4" = "#8db600", 
    "T CD8" = "#f99379", "Plasma" = "#875692", "B cell naive" = "#90AD1C",
    "B cell non-switched memory" = "#1CFFCE", "B cell switched memory" = "#2ED9FF",
    "B cell exhausted" = "#B10DA1", "cDC1 CLEC9A" = "#AA0DFE", 
    "cDC2 CD1C" = "#85660D", "cDC3 LAMP3" = "#f38400", 
    "Macrophage M1" = "#e68fac", "Macrophage M2" = "#882d17",
    "Macrophage other" = "#654522", "T CD4 naive" = "#1C8356", 
    "Th1" = "#c2b280", "Th2" = "#dcd300", "Th17" = "#F7E1A0", 
    "Tfh" = "#E2E2E2", "Treg" = "#604e97", "T CD8 naive" = "#1CBE4F",
    "T CD8 central memory" = "#C4451C", "T CD8 effector memory" = "#DEA0FD",
    "T CD8 effector" = "#FE00FA", "T CD8 exhausted" = "#325A9B"
  )
  
  # Generate main plot (cell type composition)
  svg_path <- file.path(outpath, "SpaCETdecon.pdf")
  message(paste("Generating main plot:", svg_path))
  pdf(svg_path)
  print(SpaCET.visualize.spatialFeature(
    SpaCET_obj, 
    spatialType = "CellTypeComposition",
    spatialFeatures = "MajorLineage",
    colors = colors_vector,
    pointSize = pointSize
  ))
  dev.off()
  
  # Generate individual cell type spatial distribution plots
  row_max <- apply(SpaCET_obj@results$deconvolution$propMat, 1, max)
  cell_types <- rownames(SpaCET_obj@results$deconvolution$propMat)[which(row_max > 0.1)]
  cell_types <- cell_types[cell_types != "Unidentifiable"]
  
  if (length(cell_types) > 0) {
    ncol <- min(length(cell_types), 3)
    nrow <- ceiling(length(cell_types) / 3)
    svg_path_cell <- file.path(outpath, "SpaCETdecon.eachcell.svg")
    message(paste("Generating individual cell type plots:", svg_path_cell))
    svg(svg_path_cell, width = 3 * ncol, height = 2 * nrow)
    print(SpaCET.visualize.spatialFeature(
      SpaCET_obj, 
      spatialType = "CellFraction",
      spatialFeatures = cell_types,
      sameScaleForFraction = FALSE,
      pointSize = small_pointSize,
      nrow = nrow
    ))
    dev.off()
  } else {
    message("Warning: No valid cell types found (proportion > 0.1)")
  }
  
  #deconvolution results output
  CellAnno <- SpaCET_obj@results$deconvolution$propMat
  colnames(CellAnno) <- SpaCET_obj@input$metaData[colnames(CellAnno),,drop=FALSE]$barcode
  CellAnno_filtered <- CellAnno[!grepl("Unidentifiable", rownames(CellAnno)), ]
  write.csv(t(CellAnno_filtered),paste0(outpath, "/CellAnno.txt"))
  
  message("Analysis completed successfully!")
  ##return(SpaCET_obj)
}

# Execute main function
tryCatch({
  SpaCETdeconvolution(
    outpath = opt$outpath,
    cancerType = opt$cancerType,
    rawdata_dir = opt$rawdata_dir,
    coreNo = opt$cores,
    min_genes = opt$min_genes,
    pointSize = opt$point_size,
    small_pointSize = opt$small_point_size,
    width = opt$width,
    height = opt$height
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)