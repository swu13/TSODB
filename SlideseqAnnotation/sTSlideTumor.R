#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SpaCET)
  library(dbscan)
  library(ggplot2)
  source("SpaCET.deconvolution.malignant.revised.R")
})

#--------------------------------------------------
# Parse command line arguments
#--------------------------------------------------
option_list <- list(

  make_option(c("-o", "--outpath"),type = "character",default = NULL,help = "SpaCET result directory"),
  make_option(c("-r", "--radius"),type = "double",default = 50,help = "Neighbor radius [default %default]"),
  make_option(c("-m", "--malignant_cutoff"),type = "double",default = 0.5,help = "Malignant proportion cutoff [default %default]"),
  make_option(c("-s", "--stromal_tumor_cutoff"),type = "double",default = 0.3,help = "Malignant proportion cutoff around stromal cell [default %default]"),
  make_option(c("-n", "--cores"), type="integer", default=6,help="Number of CPU cores to use for deconvolution, default is 6", metavar="number")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "Tumor-Stroma Interface Analysis"
)

opt <- parse_args(opt_parser)

if(is.null(opt$outpath)){
  print_help(opt_parser)
  stop("ERROR: --outpath parameter is required")
}

#--------------------------------------------------
# Main function
#--------------------------------------------------
SpaCETInterface <- function(
    outpath,
    neighbor_radius = 50,
    malignant_cutoff = 0.5,
    stromal_tumor_cutoff = 0.3,
    coreNo = 6
){

  rds_file <- file.path(outpath, "SpaCET_obj.rds")

  if(!file.exists(rds_file)){
    stop(paste("Cannot find:", rds_file))
  }

  message("[1/4] Loading SpaCET object")
  
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
  SpaCET_obj <- readRDS(rds_file)
  propMat <- SpaCET_obj@results$deconvolution$propMat
  coords  <- SpaCET_obj@input$spotCoordinates
  coords <- coords[colnames(propMat), , drop = FALSE]  
  celltype <- rownames(propMat)[apply(propMat, 2, which.max)]
  xcol <- intersect(c("X","x"), colnames(coords))[1]
  ycol <- intersect(c("Y","y"), colnames(coords))[1]
  plot_df <- data.frame(X = coords[[xcol]],Y = coords[[ycol]],celltype = celltype,stringsAsFactors = FALSE)
  
  svg_path <- file.path(outpath, "SpaCETdecon.pdf")
  message(paste("Generating main plot:", svg_path))
  pdf(svg_path)
  print(ggplot(plot_df, aes(X, Y, color = celltype)) + geom_point(size = 0.5) + coord_equal() +
        scale_y_reverse() +scale_color_manual(values = colors_vector,na.value = "grey80") + 
        theme_bw() + theme(panel.grid = element_blank(),axis.title = element_blank()))
  dev.off()
  
  # Tumor / Stroma
  message("[2/4] Detecting tumor-stroma interface")
  if(!("Malignant" %in% rownames(propMat))){
    stop("Cannot find Malignant lineage.")
  }  
  Content <- propMat["Malignant", ]
  Content <- ifelse(Content >= malignant_cutoff, "Tumor", "Stroma")
  
  coords_mat <- as.matrix(coords[,c(xcol,ycol)])
  knn_res <- FNN::get.knn(coords_mat, k=2)
  nearest_dist <- knn_res$nn.dist[,1]
  neighbor_radius <- 3 * median(nearest_dist)
  message("radius = ", neighbor_radius)
  neighbor_list <- dbscan::frNN(coords_mat, eps = neighbor_radius)$id
  Content_new <- Content
  
  stroma_idx <- which(Content == "Stroma")
  for(i in stroma_idx){
    neigh <- setdiff(neighbor_list[[i]], i)
    if(length(neigh) == 0){next}
    frac <- mean(Content[neigh] == "Tumor")
    if(frac > stromal_tumor_cutoff){
      Content_new[i] <- "Interface"
    }
  }

  SpaCET_obj@results$CCI$interface <- matrix(Content_new,nrow = 1,byrow = TRUE,dimnames = list("Interface",colnames(propMat)))
  saveRDS(SpaCET_obj,file.path(outpath, "SpaCET_tumor.rds"))
  
  output <- data.frame(X = coords[[xcol]],Y = coords[[ycol]],celltype = celltype,Interface=Content_new,stringsAsFactors = FALSE) 
  rownames(output) <- rownames(coords)
  write.csv(output,file.path(outpath, "cell_annotation_spatial.csv"),quote = FALSE, row.names = TRUE)

  #--------------------------------------------------
  # Interface plot
  #--------------------------------------------------

  message("[3/4] Plotting interface")

  pdf(file.path(outpath, "TumorInterfere.pdf"))
  tryCatch({
    print(
      SpaCET.visualize.spatialFeature(
        SpaCET_obj,
        spatialType = "Interface",
        spatialFeature = "Interface"
      )
    )
  }, error = function(e){
    message(paste("WARNING: Interface plot failed:",e$message))
  })
  dev.off()

  #--------------------------------------------------
  # Malignant subclusters
  #--------------------------------------------------
  
  message("[3/5] Performing malignant cell deconvolution...")
  message(paste("  - Using", coreNo, "CPU cores"))
  SpaCET_obj <- tryCatch({
    SpaCET.deconvolution.malignant.revised(SpaCET_obj, coreNo = coreNo)
  }, error = function(e) {
    stop(paste("Malignant cell deconvolution failed:", e$message))
  })
  message("  - Malignant cell deconvolution completed\n")

  message("[4/4] Plotting malignant subclusters")
  
  propMat <- SpaCET_obj@results$deconvolution$propMat
  malignantcells <- grep("^Malignant",rownames(propMat),value = TRUE)
  malignantcells <- setdiff(malignantcells,"Malignant")
  if(length(malignantcells) > 0){
    mal_mat <- propMat[malignantcells, , drop=FALSE]
    
    tumor_score <- SpaCET_obj@results$deconvolution$propMat["Malignant", ]
    is_tumor <- tumor_score >= malignant_cutoff
    
    malignant_label <- rep("Stroma", ncol(propMat))
    names(malignant_label) <- colnames(propMat)
    
    sub_idx <- which(is_tumor)
    subtype_mat <- mal_mat[, sub_idx, drop = FALSE]
    malignant_label[sub_idx] <-  malignantcells[apply(subtype_mat, 2, which.max)]
    
    plot_df <- data.frame(X = coords[[xcol]],Y = coords[[ycol]],celltype = malignant_label,stringsAsFactors = FALSE)
    stroma_df <- plot_df[plot_df$celltype == "Stroma", ]
    tumor_df  <- plot_df[plot_df$celltype != "Stroma", ]
    svg_path <- file.path(outpath, "TumorSubcluster.pdf")
    message(paste("Generating main plot:", svg_path))
    pdf(svg_path)
    print(ggplot() + geom_point(data = stroma_df,aes(X, Y),color = "grey85",size = 0.4) +
          geom_point(data = tumor_df,aes(X, Y, color = celltype),size = 0.5) +
          coord_equal() + scale_y_reverse() + theme_bw() + theme(panel.grid = element_blank(),axis.title = element_blank()))
    dev.off()
  } else {
    message("No malignant subclusters detected.")
  }
  message("Analysis completed successfully!")
}

#--------------------------------------------------
# Run
#--------------------------------------------------
tryCatch({
  SpaCETInterface(
    outpath = opt$outpath,
    neighbor_radius = opt$radius,
    malignant_cutoff = opt$malignant_cutoff,
    stromal_tumor_cutoff = opt$stromal_tumor_cutoff,
    coreNo = opt$cores
  )
}, error = function(e){

  message(paste("\nERROR:",e$message))
  quit(status = 1)
})

quit(status = 0)