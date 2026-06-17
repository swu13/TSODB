#!/usr/bin/env Rscript

# SpaCET Spatial Trajectory and Gene Enrichment Analysis Pipeline
# Usage: Rscript SpatialTrajectoryGene.R --outpath <output_directory> --rawpath <visium_raw_directory>

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(igraph)
  library(stringr)
  library(SPATA2)
  library(ggplot2)
  library(dbscan)
  library(ggpubr)
  source("plotSpatialTrajectories.revised.R")
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output directory path (required), containing seurat_CCI_obj.rds", metavar="directory"),
  make_option(c("-r", "--rawpath"), type="character", default=NULL,
              help="10X Visium raw data directory path (required)", metavar="directory"),
  make_option(c("-t", "--threshold_distance"), type="double", default=1.5,
              help="Distance thresholdfor spatial connectivity, default: 1.5", metavar="number"),
  make_option(c("--spark_pval"), type="double", default=0.05,
              help="P-value threshold for SPARKX analysis, default: 0.05", metavar="number"),
  make_option(c("--trajectory_fdr"), type="double", default=0.05,
              help="FDR threshold for trajectory screening, default: 0.05", metavar="number"),
  make_option(c("--logfc_threshold"), type="double", default=1,
              help="Log fold change threshold for gene enrichment, default: 1", metavar="number"),
  make_option(c("--gene_p_val_adj"), type="double", default=0.05,
              help="Adjusted p-value threshold for gene enrichment, default: 0.05", metavar="number"),
  make_option(c("--plot_width"), type="integer", default=8,
              help="Trajectory plot width in inches, default: 8", metavar="number"),
  make_option(c("--plot_height"), type="integer", default=8,
              help="Trajectory plot height in inches, default: 8", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list, 
                           description="SpaCET spatial trajectory analysis with SPATA2 and gene enrichment analysis")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$outpath)) {
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required", call.=FALSE)
}

if (is.null(opt$rawpath)) {
  print_help(opt_parser)
  stop("ERROR: --rawpath parameter is required", call.=FALSE)
}

#===============================================================================
# Function: SpatialTrajectoryGene
# Description: Perform spatial trajectory analysis for genes between tumor and stroma
#              regions using SPATA2 and identify enriched genes across regions
# Parameters:
#   outpath - Path to directory containing seurat_CCI_obj.rds
#   rawpath - Path to 10X Visium raw data directory
#   threshold_distance - Distance threshold for spatial connectivity
#   spark_pval - P-value threshold for SPARKX analysis
#   trajectory_fdr - FDR threshold for trajectory screening
#   logfc_threshold - Log fold change threshold for gene enrichment
#   gene_p_val_adj - Adjusted p-value threshold for gene enrichment
#   plot_width - Width of trajectory plot in inches
#   plot_height - Height of trajectory plot in inches
#===============================================================================
SpatialTrajectoryGene <- function(outpath, rawpath,
                              threshold_distance = 1.5,
                              spark_pval = 0.05,
                              trajectory_fdr = 0.05,
                              logfc_threshold = 1,
                              gene_p_val_adj = 0.05,
                              plot_width = 8,
                              plot_height = 8) {
  
  # Print start message
  message("\n========================================")
  message("Starting SpaCET Spatial Trajectory and Gene Enrichment Analysis")
  message(paste("Output path:", outpath))
  message(paste("Raw data path:", rawpath))
  message(paste("Distance threshold:", threshold_distance))
  message(paste("SPARKX p-value threshold:", spark_pval))
  message(paste("Trajectory FDR threshold:", trajectory_fdr))
  message(paste("Gene enrichment logFC threshold:", logfc_threshold))
  message(paste("Gene enrichment adj. p-value threshold:", gene_p_val_adj))
  message("========================================\n")
  
  #-----------------------------------------------------------------------------
  # Step 1: Load Seurat object and extract region annotations
  #-----------------------------------------------------------------------------
  message("[1/10] Loading Seurat object and region annotations...")
  
  seurat_file <- file.path(outpath, "seurat_CCI_obj.rds")
  if (!file.exists(seurat_file)) {
    stop(paste("Seurat object not found:", seurat_file))
  }
  
  SeuratData <- tryCatch({
    readRDS(seurat_file)
  }, error = function(e) {
    stop(paste("Failed to load Seurat object:", e$message))
  })
  message(paste("  - Loaded Seurat object from:", seurat_file))
  
  # Extract region annotations
  regions <- SeuratData@meta.data[, "Interface", drop = FALSE]
  
  #-----------------------------------------------------------------------------
  # Step 2: Create and subset SPATA2 object
  #-----------------------------------------------------------------------------
  message("\n[2/10] Creating SPATA2 object...")
  
  # Check if rawpath exists
  if (!dir.exists(rawpath)) {
    stop(paste("Raw data directory does not exist:", rawpath))
  }
  
  # Initialize SPATA2 object
  
  spata2Data <- tryCatch({
    initiateSpataObjectVisium(sample_name = "slide1",directory_visium = rawpath)
  }, error = function(e1) {
    message( "initiateSpataObjectVisium failed: ", e1$message, "\nTrying initiateSpataObjectVisiumHD...")
    tryCatch({
      initiateSpataObjectVisiumHD(sample_name = "slide1",directory_visium = rawpath)
    }, error = function(e2) {
      message(paste0("initiateSpataObjectVisiumHD failed: ", e2$message, "\n","Trying patched createSpatialDataVisium..."))
      source("/data_d/WSJ/SpatialMetsDB/Script/createSpatialDataVisium.revised.R")
      assignInNamespace(x = "createSpatialDataVisium",value = createSpatialDataVisium,ns = "SPATA2")
      tryCatch({
        initiateSpataObjectVisium(sample_name = "slide1",directory_visium = rawpath)
      }, error = function(e3) {
        message( "Patched Visium failed: ", e1$message, "\nTrying expression and metadata...")
        tryCatch({
          meta=read.csv(file.path(rawpath, "metadata.csv"),row.names=1)
          meta_df <- data.frame(barcodes = rownames(meta),x = meta$X,y = meta$Y)
          initiateSpataObject(sample_name = "slide1", modality = "gene",
                              count_mtr = GetAssayData(object = SeuratData, assay = "Spatial", layer = "counts"),
                              coords_df = tibble::as_tibble(meta_df))
        }, error = function(e4) {
          stop(paste0("\nAll methods failed.\n\n"))
        })
      })
    })
  })
  
  # Subset to spots with region annotations
  spata2Data <- subsetSpataObject(spata2Data, barcodes = rownames(regions))
  
  #-----------------------------------------------------------------------------
  # Step 3: Add region annotations as features to SPATA2
  #-----------------------------------------------------------------------------
  message("\n[3/10] Adding region annotations to SPATA2 object...")
  
  # Ensure regions match SPATA2 object
  regions <- regions[colnames(spata2Data@assays$gene@mtr_counts), , drop = FALSE]
  
  # Convert to tibble for SPATA2
  regions_tibble <- tibble::rownames_to_column(regions, var = "barcodes") %>% 
    tibble::as_tibble()
  
  # Add features
  spata2Data <- addFeatures(object = spata2Data, 
                            feature_df = regions_tibble, 
                            feature_names = c("Interface"), 
                            overwrite = TRUE)
  message("  - Added Interface annotations as features")


  #-----------------------------------------------------------------------------
  # Step 4: Identify tumor region centroid
  #-----------------------------------------------------------------------------
  message("\n[4/10] Identifying tumor region centroid...")
  
  # Subset tumor spots
  tumor_barcodes <- rownames(regions)[which(regions$Interface == "Tumor")]
  if (length(tumor_barcodes) == 0) {
    stop("No tumor spots found in region annotations")
  }
  
  spata2tumor <- subsetSpataObject(spata2Data, barcodes = tumor_barcodes)
  message(paste("  - Tumor spots:", ncol(spata2tumor)))
  
  if("x_orig" %in% names(spata2tumor@spatial@coordinates)) {
    names(spata2tumor@spatial@coordinates)[names(spata2tumor@spatial@coordinates) == "x_orig"] <- "row"
  }
  
  if("y_orig" %in% names(spata2tumor@spatial@coordinates)) {
    names(spata2tumor@spatial@coordinates)[names(spata2tumor@spatial@coordinates) == "y_orig"] <- "col"
  }
  
  # Extract coordinates
  tumor_coords <- spata2tumor@spatial@coordinates %>%
    select(barcodes, row, col) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    as.data.frame()
  
  # Create adjacency matrix based on distance threshold
#  dist_matrix <- as.matrix(dist(tumor_coords))
#  adj_matrix <- dist_matrix < threshold_distance
#  diag(adj_matrix) <- FALSE
#  
#  # Find connected components
#  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
#  clusters <- components(g) 
   coords <- as.matrix(tumor_coords)
   cl <- dbscan(coords, eps = threshold_distance, minPts = 1)
   tumor_cluster <- data.frame(cluster = cl$cluster,row.names = rownames(tumor_coords))  
  # Assign cluster IDs
#  tumor_cluster <- data.frame(cluster = as.numeric(clusters$membership),
#                              row.names = names(clusters$membership))
  tumor_cluster$row <- tumor_coords[rownames(tumor_cluster), "row"]
  tumor_cluster$col <- tumor_coords[rownames(tumor_cluster), "col"]
  
  # Find centroids of large clusters
  tumor_centroids <- tumor_cluster %>%
    group_by(cluster) %>%
    summarise(centroid_x = mean(row, na.rm = TRUE),
              centroid_y = mean(col, na.rm = TRUE),
              n_spots = n()) %>%
    arrange(desc(n_spots)) %>%
    slice(1) %>%
    ungroup()
  
  # Find spot closest to centroid
  tumor_centroids_spot <- tumor_cluster %>%
    inner_join(tumor_centroids, by = "cluster") %>%
    group_by(cluster) %>%
    mutate(dist = (row - centroid_x)^2 + (col - centroid_y)^2) %>%
    slice_min(dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster, centroid_row = row, centroid_col = col) %>%
    as.data.frame()
  
  message(paste("  - Selected tumor centroid with", tumor_centroids$n_spots[1], "spots"))
  message(paste("    Centroid spot: row", tumor_centroids_spot$centroid_row[1], 
                "col", tumor_centroids_spot$centroid_col[1]))
    
  #-----------------------------------------------------------------------------
  # Step 5: Identify stroma region centroid
  #-----------------------------------------------------------------------------
  message("\n[5/10] Identifying stroma region centroid...")
  
  # Subset stroma spots
  stroma_barcodes <- rownames(regions)[which(regions$Interface == "Stroma")]
  if (length(stroma_barcodes) == 0) {
    interface_type <- regions$Interface
    interface_type[interface_type == "Interface"] <- "Stroma"
    stroma_barcodes <- rownames(regions)[interface_type == "Stroma"]
    if (length(stroma_barcodes) == 0) {
      stop("No stroma and Interface spots found in region annotations")
    }
  }
  
  spata2Stroma <- subsetSpataObject(spata2Data, barcodes = stroma_barcodes)
  message(paste("  - Stroma spots:", ncol(spata2Stroma)))
  
  if("x_orig" %in% names(spata2Stroma@spatial@coordinates)) {
    names(spata2Stroma@spatial@coordinates)[names(spata2Stroma@spatial@coordinates) == "x_orig"] <- "row"
  }
  
  if("y_orig" %in% names(spata2Stroma@spatial@coordinates)) {
    names(spata2Stroma@spatial@coordinates)[names(spata2Stroma@spatial@coordinates) == "y_orig"] <- "col"
  }
  
  # Extract coordinates
  Stroma_coords <- spata2Stroma@spatial@coordinates %>%
    select(barcodes, row, col) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    as.data.frame()
  
  # Create adjacency matrix
#  dist_matrix <- as.matrix(dist(Stroma_coords))
#  adj_matrix <- dist_matrix < threshold_distance
#  diag(adj_matrix) <- FALSE
#  
#  # Find connected components
#  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
#  clusters <- components(g)
  
  coords <- as.matrix(Stroma_coords)
  cl <- dbscan(coords, eps = threshold_distance, minPts = 1)
  Stroma_cluster <- data.frame(cluster = cl$cluster,row.names = rownames(Stroma_coords)) 
  
  # Assign cluster IDs
#  Stroma_cluster <- data.frame(cluster = as.numeric(clusters$membership),
#                               row.names = names(clusters$membership))
  Stroma_cluster$row <- Stroma_coords[rownames(Stroma_cluster), "row"]
  Stroma_cluster$col <- Stroma_coords[rownames(Stroma_cluster), "col"]
  
  # Find centroids of large clusters
  Stroma_centroids <- Stroma_cluster %>%
    group_by(cluster) %>%
    summarise(centroid_x = round(mean(row, na.rm = TRUE)),
              centroid_y = round(mean(col, na.rm = TRUE)),
              n_spots = n()) %>%
    arrange(desc(n_spots)) %>%
    slice(1) %>%
    ungroup() #%>% filter(n_spots >= min_spots_per_cluster) 
  
  # Find spot closest to centroid
  Stroma_centroids_spot <- Stroma_cluster %>%
    inner_join(Stroma_centroids, by = "cluster") %>%
    group_by(cluster) %>%
    mutate(dist = (row - centroid_x)^2 + (col - centroid_y)^2) %>%
    slice_min(dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster, centroid_row = row, centroid_col = col) %>%
    as.data.frame()
  
  message(paste("  - Selected stroma centroid with", Stroma_centroids$n_spots[1], "spots"))
  message(paste("    Centroid spot: row", Stroma_centroids_spot$centroid_row[1], 
                "col", Stroma_centroids_spot$centroid_col[1]))
    
  #-----------------------------------------------------------------------------
  # Step 6: Convert coordinates to pixels and add trajectory
  #-----------------------------------------------------------------------------
  message("\n[6/10] Creating tumor-stroma trajectory...")
  
  
  # Get coordinates for centroid spots
  coords <- spata2Data@spatial@coordinates
  
  if(length(spata2Data@spatial@images)>0){
    # Tumor centroid pixel coordinates
    tumor_centroid_coords <- coords %>%
      filter(row == tumor_centroids_spot[1, "centroid_row"], 
             col == tumor_centroids_spot[1, "centroid_col"]) %>%
      select(x_orig, y_orig)
    tumor_pixels <- tumor_centroid_coords * spata2Data@spatial@images$lowres@scale_factors$image
    message(paste("  - Tumor centroid pixels: (", tumor_pixels$x_orig, ",", tumor_pixels$y_orig, ")"))
    # Stroma centroid pixel coordinates
    Stroma_centroid_coords <- coords %>%
      filter(row == Stroma_centroids_spot[1, "centroid_row"], 
             col == Stroma_centroids_spot[1, "centroid_col"]) %>%
      select(x_orig, y_orig)
    Stroma_pixels <- Stroma_centroid_coords * spata2Data@spatial@images$lowres@scale_factors$image
    message(paste("  - Stroma centroid pixels: (", Stroma_pixels$x_orig, ",", Stroma_pixels$y_orig, ")"))
    # Add trajectory
    spata2Data <- addSpatialTrajectory(object = spata2Data,
                                       id = "Tumor-Stroma",
                                       start = as.numeric(tumor_pixels),
                                       end = as.numeric(Stroma_pixels),
                                       overwrite = TRUE)
  }else{
    tumor_centroid_coords <- coords %>%
      filter(x_orig == tumor_centroids_spot[1, "centroid_row"], 
             y_orig == tumor_centroids_spot[1, "centroid_col"]) %>%
      select(x_orig, y_orig)
    Stroma_centroid_coords <- coords %>%
      filter(x_orig == Stroma_centroids_spot[1, "centroid_row"], 
             y_orig == Stroma_centroids_spot[1, "centroid_col"]) %>%
      select(x_orig, y_orig)
    spata2Data <- addSpatialTrajectory(object = spata2Data,
                                     id = "Tumor-Stroma",
                                     start = as.numeric(tumor_centroid_coords[1, ]),
                                     end = as.numeric(Stroma_centroid_coords[1, ]),
                                     overwrite = TRUE)
  }
  message("  - Added trajectory: Tumor -> Stroma")
    
  #-----------------------------------------------------------------------------
  # Step 7: Generate trajectory visualization
  #-----------------------------------------------------------------------------
  message("\n[7/10] Generating trajectory visualization...")
  
  svg_file <- file.path(outpath, "Trajectory.pdf")
  pdf(svg_file, width = plot_width, height = plot_height)
  
  tryCatch({
    print(plotSpatialTrajectories(object = spata2Data,
                                  color_by = "Interface",
                                  clrp_adjust = c(Tumor = "#F3C300",
                                                  Interface = "#000000",
                                                  Stroma = "#A9A9A9"),
                                  sgmt_clr = "#800000") + scale_y_reverse())
    message(paste("  - Trajectory plot saved to:", svg_file))
  }, error = function(e) {
    message(paste("  - WARNING: Plot generation failed:", e$message))
  })
  
  dev.off()
    
  #-----------------------------------------------------------------------------
  # Step 8: Identify genes along trajectory using SPARKX
  #-----------------------------------------------------------------------------
  message("\n[8/10] Identifying genes along trajectory...")
  
  # Run SPARKX analysis
  message("  - Running SPARKX analysis...")
  spata2Data <- tryCatch({
    runSPARKX(object = spata2Data)
  }, error = function(e) {
    message(paste("  - WARNING: SPARKX analysis failed:", e$message))
  })
  
  # Get significant SPARKX genes
  spark_df <- tryCatch({
    getSparkxGeneDf(object = spata2Data, threshold_pval = spark_pval)
  }, error = function(e) {
    message(paste("  - WARNING: Failed to get SPARKX results:", e$message))
    return(data.frame(genes = character()))
  })
  
  if (nrow(spark_df) > 0) {
    message(paste("  - Found", nrow(spark_df), "significant SPARKX genes"))
    
    # Perform trajectory screening
    message("  - Screening genes along trajectory...")
    sts_out <- tryCatch({
      spatialTrajectoryScreening(object = spata2Data,
                                 id = "Tumor-Stroma",
                                 variables = spark_df[["genes"]])
    }, error = function(e) {
      message(paste("  - WARNING: Trajectory screening failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(sts_out)) {
      # Extract significant genes
      sign_df <- sts_out@results$significance %>%
        filter(fdr < trajectory_fdr)
      
      if (nrow(sign_df) > 0) {
        message(paste("  - Found", nrow(sign_df), "significant genes along trajectory (FDR <", trajectory_fdr, ")"))
        
        # Get best model fits
        best_fits <- sts_out@results$model_fits %>%
          filter(variables %in% sign_df[["variables"]]) %>%
          group_by(variables) %>%
          slice_min(mae, n = 1)
        
        # Create output table
        output <- left_join(best_fits, sign_df, by = "variables") %>%
          select(gene = variables, 
                 models, 
                 variations = rel_var, 
                 p_value, 
                 fdr) %>%
          as.data.frame()
        
        sts_df <- tryCatch({
          getStsDf(spata2Data, id = "Tumor-Stroma", variables = spark_df[["genes"]])
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
          sts_df <- setNames(data.frame(matrix(ncol = length(spark_df[["genes"]])+2, nrow = 0)),c("dist","dist_unit",spark_df[["genes"]]))
        }
        
      } else {
        message(paste("  - No significant genes found with FDR <", trajectory_fdr))
        output <- data.frame(geneset = character(),models = character(),variations = numeric(),p_value = numeric(),fdr = numeric(),stringsAsFactors = FALSE)
        sts_df <- setNames(data.frame(matrix(ncol = length(spark_df[["genes"]])+2, nrow = 0)),c("dist","dist_unit",spark_df[["genes"]]))
      }
    } else {
      message("  - No significant SPARKX genes found")
      sts_df <- setNames(data.frame(matrix(ncol = length(spark_df[["genes"]])+2, nrow = 0)),c("dist","dist_unit",spark_df[["genes"]]))
      output <- data.frame(geneset = character(),models = character(),variations = numeric(),p_value = numeric(),fdr = numeric(),stringsAsFactors = FALSE)
    }
  }
  #Extract and save line plot data
  message("\n Extracting trajectory line plot data...") 
  # Save line plot data
  lineplot_file <- file.path(outpath, "TrajectoryGenesLinePlot.txt")
  write.table(sts_df, file = lineplot_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Line plot data saved to:", lineplot_file))
  # Save trajectory genes
  trajectory_file <- file.path(outpath, "TrajectoryGenes.txt")
  write.table(output, file = trajectory_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Trajectory genes saved to:", trajectory_file))
  # Save SPATA2 object
  spata_file <- file.path(outpath, "spata_obj.rds")
  saveRDS(spata2Data, file = spata_file)
  message(paste("  - SPATA2 object saved to:", spata_file))
  
  #-----------------------------------------------------------------------------
  # Step 9: Perform gene enrichment analysis across regions
  #-----------------------------------------------------------------------------
  message("\n[9/10] Performing gene enrichment analysis across regions...")
  
  # Set identity to Interface regions
  Idents(SeuratData) <- "Interface"
  DefaultAssay(SeuratData) <- "Spatial"
  
  # Run SCTransform normalization
  options(future.globals.maxSize = 8 * 1024^3)
  message("  - Running SCTransform normalization...")
  SeuratData <- tryCatch({
    SCTransform(SeuratData, assay = "Spatial", verbose = FALSE)
  }, error = function(e) {
    message(paste("  - WARNING: SCTransform failed:", e$message))
  })
  
  # Find markers for all regions
  message(paste("  - Finding markers with logFC threshold:", logfc_threshold))
  DEGs <- tryCatch({
    FindAllMarkers(SeuratData, assay = "SCT", logfc.threshold = 0, only.pos = TRUE, verbose = FALSE)
  }, error = function(e) {
    message(paste("  - WARNING: FindAllMarkers failed:", e$message))
    return(data.frame())
  })
  
  if (nrow(DEGs) > 0) {
    # Filter significant DEGs
    SigDEGs <- DEGs %>%
      filter(p_val_adj < gene_p_val_adj) %>%
      filter(avg_log2FC > logfc_threshold) %>%
      filter(pct.1 >= pct.2)
    message(paste("  - Found", nrow(SigDEGs), "significant enriched genes"))
  } else {
    message("  - No DEGs found")
    SigDEGs <- data.frame(p_val = numeric(),avg_log2FC = numeric(),pct.1 = numeric(),pct.2 = numeric(),p_val_adj = numeric(),cluster=character(),gene=character(),stringsAsFactors = FALSE)
  }
  # Save gene enrichment results
  gene_file <- file.path(outpath, "GeneEnriched.txt")
  write.table(SigDEGs, file = gene_file,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message(paste("  - Gene enrichment results saved to:", gene_file))
    
  # Save updated Seurat object
  seurat_gene_file <- file.path(outpath, "seurat_gene_obj.rds")
  saveRDS(SeuratData, file = seurat_gene_file)
  message(paste("  - Updated Seurat object saved to:", seurat_gene_file))
  
  message("\nAnalysis completed successfully!")
}

#===============================================================================
# Execute main function
#===============================================================================

# Run the main function with error handling
tryCatch({
  SpatialTrajectoryGene(
    outpath = opt$outpath,
    rawpath = opt$rawpath,
    threshold_distance = opt$threshold_distance,
    spark_pval = opt$spark_pval,
    trajectory_fdr = opt$trajectory_fdr,
    logfc_threshold = opt$logfc_threshold,
    gene_p_val_adj = opt$gene_p_val_adj,
    plot_width = opt$plot_width,
    plot_height = opt$plot_height
  )
}, error = function(e) {
  message(paste("\nERROR:", e$message))
  quit(status=1)
})

quit(status=0)