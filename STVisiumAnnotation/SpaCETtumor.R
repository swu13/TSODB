#!/usr/bin/env Rscript

# SpaCET Tumor Interface and Malignant Cell Analysis Pipeline
# Usage: Rscript SpaCET_tumor_analysis.R --outpath <output_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SpaCET)
  source("SpaCET.deconvolution.malignant.revised.R")
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), containing SpaCET_obj.rds file", metavar="directory"),
  make_option(c("-n", "--cores"), type="integer", default=6,
              help="Number of CPU cores to use for deconvolution, default is 6", metavar="number"),
  make_option(c("--interface_width"), type="integer", default=8,
              help="Interface plot width (inches), default is 8", metavar="number"),
  make_option(c("--interface_height"), type="integer", default=8,
              help="Interface plot height (inches), default is 8", metavar="number"),
  make_option(c("--subcluster_width_factor"), type="double", default=4,
              help="Width factor per malignant subcluster plot (inches), default is 4", metavar="number"),
  make_option(c("--subcluster_height_factor"), type="double", default=4,
              help="Height factor per malignant subcluster plot (inches), default is 4", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="SpaCET tumor interface identification and malignant cell deconvolution analysis script")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

#===============================================================================
# Function: SpaCETtumor
# Description: Perform tumor interface identification and visualize malignant 
#              cell subclusters from a pre-processed SpaCET object
# Parameters:
#   outpath - Path to directory containing SpaCET_obj.rds file
#   coreNo - Number of CPU cores for malignant cell deconvolution
#   interface_width - Width of interface plot in inches
#   interface_height - Height of interface plot in inches
#   subcluster_width_factor - Width factor per malignant subcluster plot
#   subcluster_height_factor - Height factor per malignant subcluster plot
#===============================================================================
SpaCETtumor <- function(outpath, coreNo = 6, 
                        interface_width = 8, interface_height = 8,
                        subcluster_width_factor = 4, subcluster_height_factor = 4) {
  
  # Print start message
  message("\n========================================")
  message("Starting SpaCET Tumor Analysis")
  message(paste("Output path:", outpath))
  message(paste("CPU cores:", coreNo))
  message("========================================\n")
  
  # Check if input file exists
  input_rds <- file.path(outpath, "SpaCET_obj.rds")
  if (!file.exists(input_rds)) {
    stop(paste("Input file does not exist:", input_rds))
  }
  
  #-----------------------------------------------------------------------------
  # Step 1: Load SpaCET object
  #-----------------------------------------------------------------------------
  message("[1/5] Loading SpaCET object...")
  SpaCET_obj <- tryCatch({
    readRDS(input_rds)
  }, error = function(e) {
    stop(paste("Failed to load SpaCET object:", e$message))
  })
  message(paste("  - Loaded from:", input_rds))
  message("  - SpaCET object loaded successfully\n")
  
  #-----------------------------------------------------------------------------
  # Step 2: Identify tumor-stroma interface
  #-----------------------------------------------------------------------------
  
  message("[2/5] Identifying tumor-stroma interface...")
  SpaCET_obj <- tryCatch({
    SpaCET.identify.interface(SpaCET_obj)
  }, error = function(e) {
    stop(paste("Interface identification failed:", e$message))
  })
  message("  - Interface identification completed\n")
  
  interface <- as.data.frame(t(SpaCET_obj@results$CCI$interface))
  if(!is.null(SpaCET_obj@input$metaData)){
    rownames(interface) <- SpaCET_obj@input$metaData[rownames(interface),,drop=FALSE]$barcode
  }
  write.csv(interface,paste0(outpath, "/RegionAnno.txt"))
  
  #-----------------------------------------------------------------------------
  # Step 3: Perform malignant cell deconvolution
  #-----------------------------------------------------------------------------
  message("[3/5] Performing malignant cell deconvolution...")
  message(paste("  - Using", coreNo, "CPU cores"))
  SpaCET_obj <- tryCatch({
    SpaCET.deconvolution.malignant.revised(SpaCET_obj, coreNo = coreNo)
  }, error = function(e) {
    stop(paste("Malignant cell deconvolution failed:", e$message))
  })
  message("  - Malignant cell deconvolution completed\n")
  
  #-----------------------------------------------------------------------------
  # Step 4: Save updated SpaCET object
  #-----------------------------------------------------------------------------
  message("[4/5] Saving updated SpaCET object...")
  output_rds <- file.path(outpath, "SpaCET_tumor.rds")
  saveRDS(SpaCET_obj, file = output_rds)
  message(paste("  - Saved to:", output_rds))
  
  #-----------------------------------------------------------------------------
  # Step 5: Generate visualizations
  #-----------------------------------------------------------------------------
  message("[5/5] Generating visualizations...")
  
  # Visualization 1: Tumor-stroma interface
  message("  - Generating tumor-stroma interface plot...")
  interface_svg <- file.path(outpath, "TumorInterfere.svg")
  svg(interface_svg, width = interface_width, height = interface_height)
  tryCatch({
    print(SpaCET.visualize.spatialFeature(SpaCET_obj, 
                                          spatialType = "Interface", 
                                          spatialFeature = "Interface"))
    message(paste("    Saved to:", interface_svg))
  }, error = function(e) {
    message(paste("    WARNING: Interface plot generation failed:", e$message))
  })
  dev.off()
  
  # Visualization 2: Malignant cell subclusters
  message("  - Identifying malignant cell subclusters...")
  
  # Extract malignant cell types from deconvolution results
  cells <- rownames(SpaCET_obj@results$deconvolution$propMat)
  malignantcells <- cells[grepl("Malignant cell", cells, ignore.case = TRUE)]
  
  if (length(malignantcells) > 0) {
    message(paste("    Found", length(malignantcells), "malignant cell subclusters:"))
    # Calculate dynamic plot dimensions
    ncol <- min(length(malignantcells), 3)
    nrow <- ceiling(length(malignantcells) / 3)
    plot_width <- subcluster_width_factor * ncol
    plot_height <- subcluster_height_factor * nrow
    
    subcluster_svg <- file.path(outpath, "TumorSubcluster.svg")
    svg(subcluster_svg, width = plot_width, height = plot_height)
    tryCatch({
      print(SpaCET.visualize.spatialFeature(SpaCET_obj, 
                                            spatialType = "CellFraction",
                                            spatialFeatures = malignantcells, 
                                            nrow = nrow))
      message(paste("    Saved to:", subcluster_svg))
    }, error = function(e) {
      message(paste("    WARNING: Subcluster plot generation failed:", e$message))
    })
    dev.off()
  } else {
    message("    WARNING: No malignant cell subclusters found in deconvolution results")
    message("    Check if 'Malignant cell' pattern exists in cell type names")
  }
    
  message("\nAnalysis completed successfully!")
  ##return(SpaCET_obj)
}

#===============================================================================
# Execute main function
#===============================================================================

# Run the main function with error handling
tryCatch({
  SpaCETtumor(
    outpath = opt$outpath,
    coreNo = opt$cores,
    interface_width = opt$interface_width,
    interface_height = opt$interface_height,
    subcluster_width_factor = opt$subcluster_width_factor,
    subcluster_height_factor = opt$subcluster_height_factor
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)