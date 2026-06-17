#!/usr/bin/env Rscript

# SpaCET Cell-Cell Colocalization and Communication Analysis Pipeline
# Usage: Rscript SpaCET_colocalization.R --outpath <output_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SpaCET)
  library(mistyR)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), containing SpaCET_tumor.rds file", metavar="directory"),
  make_option(c("--min_cells"), type="integer", default=10,
              help="Minimum number of cells/spots required for region analysis, default: 10", metavar="number")
  )

opt_parser <- OptionParser(option_list=option_list, 
                           description="SpaCET cell-cell colocalization and communication analysis using mistyR")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

#===============================================================================
# Function: SpaCCcolocalization
# Description: Perform cell-cell colocalization analysis and communication 
#              network inference using mistyR across different spatial regions
# Parameters:
#   outpath - Path to directory containing SpaCET_tumor.rds file
#   min_cells - Minimum number of spots required for region analysis
#===============================================================================
SpaCCcolocalization <- function(outpath, min_cells = 10) {
  
  # Print start message
  message("\n========================================")
  message("Starting SpaCET Cell-Cell Colocalization Analysis")
  message(paste("Output path:", outpath))
  message(paste("Minimum cells per region:", min_cells))
  message("========================================\n")

  regions <- c("Tumor","Interface","Stroma")
  
  # Check if input file exists
  input_rds <- file.path(outpath, "SpaCET_tumor.rds")
  if (!file.exists(input_rds)) {
      stop(paste("Input file does not exist. Expected:", file.path(outpath, "SpaCET_tumor.rds")))
  } 
  
  #-----------------------------------------------------------------------------
  # Step 1: Load SpaCET object
  #-----------------------------------------------------------------------------
  message("[1/6] Loading SpaCET object...")
  SpaCET_obj <- tryCatch({
    readRDS(input_rds)
  }, error = function(e) {
    stop(paste("Failed to load SpaCET object:", e$message))
  })
  message(paste("  - Loaded from:", input_rds))
  
  # Check if deconvolution results exist
  if (is.null(SpaCET_obj@results$deconvolution$propMat)) {
    stop("Deconvolution results not found in SpaCET object. Please run SpaCET.deconvolution first.")
  }
  message("  - Deconvolution results found\n")
  
  #-----------------------------------------------------------------------------
  # Step 2: Extract and prepare predictions data
  #-----------------------------------------------------------------------------
  message("[2/6] Extracting cell type proportion data...")
  
  # Extract cell type proportion matrix
  predictions <- t(SpaCET_obj@results$deconvolution$propMat) %>% as.data.frame()
  message(paste("  - Cell types:", ncol(predictions)))
  message(paste("  - Spots/cells:", nrow(predictions)))
  
  #-----------------------------------------------------------------------------
  # Step 3: Initialize output dataframes
  #-----------------------------------------------------------------------------
  message("[3/6] Initializing output data structures...")
  
  medianoutput <- data.frame(cell = "", region = "", median = 0, pval = "")
  mistryoutput <- data.frame(region = "", Predictor = "", Target = "", Importance = 0)
  
  # Clean names for mistyR compatibility (remove spaces and special characters)
  original_names <- colnames(predictions)
  clean_names <- gsub("[ |-]", "_", original_names)  
  name_map <- setNames(original_names, clean_names)
  
  # Rename predictions columns for mistyR
  message("  - Output structures initialized\n")
  
  #-----------------------------------------------------------------------------
  # Step 4: Analyze each region
  #-----------------------------------------------------------------------------
  message("[4/6] Analyzing regions...")
  
  for (region in regions) {
    message(paste("\n  >>> Processing region:", region))
    
    # Identify spots belonging to current region
    spots <- colnames(SpaCET_obj@results$CCI$interface)[
      which(as.character(SpaCET_obj@results$CCI$interface) == region)
    ]
    
    # Filter spots that exist in predictions
    spots <- spots[spots %in% rownames(predictions)]
    
    if (length(spots) < min_cells) {
      message(paste("    WARNING: Only", length(spots), "spots found for region '", region, 
                    "'. Skipping (minimum required:", min_cells, ")"))
      next
    }
    
    # Extract region-specific predictions
    regionpreds <- predictions[spots, , drop = FALSE]
    regionothers <- predictions[!rownames(predictions) %in% spots, , drop = FALSE]
    
    # Calculate median proportions for this region
    message("    Calculating median proportions and Wilcoxon tests...")
    median_vals <- apply(regionpreds, 2, median, na.rm = TRUE)
    median_vals <- round(median_vals, digits = 4)
    
    # Perform Wilcoxon rank-sum tests for each cell type
    plist <- sapply(colnames(predictions), function(ct) {
      if (length(regionpreds[[ct]]) > 0 && length(regionothers[[ct]]) > 0) {
        test_result <- wilcox.test(regionpreds[[ct]], regionothers[[ct]])
        return(test_result$p.value)
      } else {
        return("")
      }
    })
    
    # Convert p-values to significance symbols
    plist_symbols <- case_when(
      plist < 0.001 ~ "***",
      plist < 0.01  ~ "**",
      plist < 0.05 ~ "*",
      TRUE ~ ""
    )
    
    # Add to median output
    temp_median <- data.frame(
      cell = names(median_vals),
      region = rep(region, length(median_vals)),
      median = as.numeric(median_vals),
      pval = plist_symbols,
      stringsAsFactors = FALSE
    )
    medianoutput <- rbind(medianoutput, temp_median)
    
    #-----------------------------------------------------------------------------
    # Step 5: Run mistyR analysis for current region
    #-----------------------------------------------------------------------------
    message("    Running mistyR analysis...")
    colnames(regionpreds) <- clean_names
    
    # Create mistyR view
    misty.views <- tryCatch({
      create_initial_view(regionpreds)
    }, error = function(e) {
      message(paste("      ERROR creating mistyR view:", e$message))
      return(NULL)
    })
    
    if (!is.null(misty.views)) {
      # Create region-specific output folder
      region_folder <- file.path(outpath, paste0("mistyR_", region))
      dir.create(region_folder, showWarnings = FALSE, recursive = TRUE)
      
      # Run mistyR
      run_misty(misty.views, results.folder = region_folder)
      
      # Collect results
      mistry_results <- tryCatch({
        collect_results(region_folder)
      }, error = function(e) {
        message(paste("      ERROR collecting mistyR results:", e$message))
        return(NULL)
      })
      
      if (!is.null(mistry_results)) {
        # Extract importances
        Importances <- as.data.frame(mistry_results$importances)
        
        # Handle NA and negative values
        Importances$Importance <- ifelse(Importances$Importance <= 0 | is.na(Importances$Importance), 
                                         0, Importances$Importance)
        Importances$Importance <- round(Importances$Importance, digits = 4)
        
        # Add region information
        mistry_importances <- cbind(
          region = rep(region, dim(Importances)[1]),
          Importances[, c("Predictor", "Target", "Importance")]
        )
        
        mistryoutput <- rbind(mistryoutput, mistry_importances)
        message(paste("      mistyR completed. Interactions found:", nrow(Importances)))
      } else {
        message("      WARNING: mistyR results collection failed")
      }
    } else {
      message("      WARNING: mistyR view creation failed")
    }
  }
  
  #-----------------------------------------------------------------------------
  # Step 6: Filter and save results
  #-----------------------------------------------------------------------------
  message("\n[5/6] Filtering and processing results...")
  
  # Remove the initial empty row
  medianoutput <- medianoutput[-1, ]
  mistryoutput <- mistryoutput[-1, ]
  
  # Filter median output
  message("  Filtering median proportion results...")
  medianoutput_filtered <- medianoutput %>%
    group_by(cell) %>%
    filter(any(median != 0)) %>%
    filter(cell != "Unidentifiable") %>%
    filter(!grepl("Malignant cell state", cell, ignore.case = TRUE)) %>%
    ungroup() %>%
    as.data.frame()
  
  # Filter mistyR output
  if (nrow(mistryoutput) > 0) {
    message("  Filtering mistyR interaction results...")
    
    # Map back to original names
    mistryoutput$Predictor <- name_map[as.character(mistryoutput$Predictor)]
    mistryoutput$Target <- name_map[as.character(mistryoutput$Target)]
    
    # Remove NA mappings and filter unwanted cell types
    mistryoutput_filtered <- mistryoutput %>%
      filter(Predictor != "Unidentifiable", Target != "Unidentifiable") %>%
      filter(!grepl("Malignant cell state", Predictor, ignore.case = TRUE)) %>%
      filter(!grepl("Malignant cell state", Target, ignore.case = TRUE))

  } else {
    message("  WARNING: No mistyR results to filter")
    mistryoutput_filtered <- mistryoutput
  }
  
  #-----------------------------------------------------------------------------
  # Step 7: Save results
  #-----------------------------------------------------------------------------
  message("\n[6/6] Saving results...")
  
  # Save median proportions
  median_file <- file.path(outpath, "median_CC.txt")
  write.table(medianoutput_filtered, file = median_file, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  Median proportions saved to:", median_file))
  
  # Save mistyR interactions
  if (nrow(mistryoutput_filtered) > 0) {
    misty_file <- file.path(outpath, "mistyR_CC.txt")
    write.table(mistryoutput_filtered, file = misty_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    message(paste("  mistyR interactions saved to:", misty_file))
  } else {
    message("  WARNING: No mistyR interactions to save")
  }
  
  # Generate summary statistics
  message("\nAnalysis completed successfully!")
}

#===============================================================================
# Execute main function
#===============================================================================


# Run the main function with error handling
tryCatch({
  SpaCCcolocalization(
    outpath = opt$outpath,
    min_cells = opt$min_cells
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)