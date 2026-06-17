#!/usr/bin/env Rscript

# Seurat Differential Expression Analysis (DEG)
# Input: seurat.rds (must contain RNA assay + group metadata)
# Output: pairwise DEGs + merged table

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

#===============================================================================
# Command line arguments
#===============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type="character",help="Input seurat.rds file"),
  make_option(c("-o", "--outdir"), type="character",help="Output directory")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "Seurat DEG pairwise analysis"
)

opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

#===============================================================================
# Main Function
#===============================================================================

DEGAnalysis <- function(inputfile, outdir){

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: Load seuratdataect
  seuratdata <- tryCatch({
    readRDS(inputfile)
  }, error = function(e){
    stop(paste("Failed to read RDS:", e$message))
  })

  if (!"group" %in% colnames(seuratdata@meta.data)) {
    stop("metadata column 'group' not found")
  }

  DefaultAssay(seuratdata) <- "Spatial"
  Idents(seuratdata) <- "group"
  groups <- unique(seuratdata$group)
  seuratdata <- NormalizeData(seuratdata, assay = DefaultAssay(seuratdata), normalization.method = "LogNormalize", scale.factor = 1e4)
  
  # Step 2: DEG analysis
  all_results <- list()
  idx <- 1
  avg_exp <- AverageExpression(seuratdata,assays = DefaultAssay(seuratdata),slot = "data")[[1]]
  colnames(avg_exp) <- levels(factor(seuratdata$group))
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      markers <- tryCatch({FindMarkers(seuratdata,ident.1 = g1,ident.2 = g2,test.use = "wilcox",slot = "data",logfc.threshold = 0,min.pct = 0)},
                          error = function(e) {message("ERROR: ", g1, " vs ", g2, " | ", e$message);return(NULL)})
      if(is.null(markers)) next
      markers$group1 <- g1
      markers$group2 <- g2
      markers$mean1 <- avg_exp[rownames(markers), g1]
      markers$mean2 <- avg_exp[rownames(markers), g2]
      markers$gene <- rownames(markers)
      rownames(markers) <- seq_len(nrow(markers))
      all_results[[idx]] <- markers
      idx <- idx + 1
    }
  }
  # Step 3: Merge results
  all_results_df <- do.call(rbind, all_results)
  write.csv(all_results_df,file = file.path(outdir,"GeneEnriched.txt"),quote = FALSE)
  seuratdata@assays$Spatial@layers$counts=NULL
  saveRDS(seuratdata, file = file.path(outdir,"seurat.gene.rds"))
}

#===============================================================================
# Run
#===============================================================================

tryCatch({

  DEGAnalysis(
    inputfile = opt$input,
    outdir = opt$outdir
  )

}, error = function(e){

  message("\nERROR:")
  message(e$message)
  quit(status = 1)

})

quit(status = 0)