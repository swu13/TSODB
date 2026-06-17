#!/usr/bin/env Rscript

# SpatialDecon Cell Fraction Differential Analysis
# Usage:
# Rscript SpatialDecon_DE.R \
#   --input expression.txt \
#   --clinical clinical.txt \
#   --outdir result

suppressPackageStartupMessages({
  library(optparse)
  library(SpatialDecon)
  library(Seurat)
  library(CIBERSORT)
})

#===============================================================================
# Command line arguments
#===============================================================================

option_list <- list(
  make_option(c("-i", "--input"),type = "character",help = "Expression matrix file"),
  make_option(c("-c", "--clinical"),type = "character",help = "Clinical/group annotation file"),
  make_option(c("-o", "--outdir"),type = "character",help = "Output directory")
)

opt_parser <- OptionParser(option_list = option_list,description = "SpatialDecon differential cell fraction analysis")

opt <- parse_args(opt_parser)

if (is.null(opt$input) ||
    is.null(opt$clinical) ||
    is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

#===============================================================================
# Main Function
#===============================================================================

SpatialDeconAnalysis <- function(
    inputfile,
    clinicalfile,
    outdir
){
  dir.create(outdir,recursive = TRUE,showWarnings = FALSE)
  # Step 1. Read data
  data <- read.table(inputfile,sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  clinical <- read.table(clinicalfile,sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  colnames(clinical)[1] <- "group"
  colnames(clinical)[2] <- "slide"
  # Step 2. spatialdecon
  raw <- as.matrix(data)
  norm = sweep(raw, 2, colSums(raw), "/") * mean(colSums(raw))
  seuratdata <- CreateSeuratObject(counts = norm,assay = "Spatial")
  data("safeTME")
  data("safeTME.matches")
  res <- runspatialdecon(object = seuratdata,bg = 0.01,X = safeTME,align_genes = TRUE)
  celldeconvoluted <- as.data.frame(res$beta)
  # Step 3: cibersortx
  data(LM22)
  immunedeconvoluted <- cibersort(sig_matrix = LM22, mixture_file = norm)
  deconvoluted <- rbind(celldeconvoluted,t(immunedeconvoluted))
  deconvoluted <- deconvoluted[1:40,]
  write.table(deconvoluted,file = file.path(outdir,"CellFraction.txt"),sep = "\t",quote = FALSE)
  # Step 4. Differential analysis
  deconvoluted <- t(deconvoluted)
  common.samples <- intersect(rownames(clinical),rownames(deconvoluted))
  clinical <- clinical[common.samples,,drop=FALSE]
  deconvoluted <- deconvoluted[common.samples,,drop=FALSE]
  seuratdata[["SpatialDecon"]] <- CreateAssayObject(counts = as.matrix(t(deconvoluted)))      
  seuratdata$group <- clinical[colnames(seuratdata),"group"]
  seuratdata$slide <- clinical[colnames(seuratdata),"slide"]
  saveRDS(seuratdata, file = file.path(outdir,"seurat.rds"))
  DefaultAssay(seuratdata) <- "SpatialDecon"
  Idents(seuratdata) <- "group"
  avg_exp <- AverageExpression(seuratdata,assays = DefaultAssay(seuratdata),slot = "count")[[1]]
  colnames(avg_exp) <- levels(factor(seuratdata$group))
  groups <- unique(clinical$group)
  all_results <- list()
  idx <- 1
  for(i in 1:(length(groups)-1)){
    for(j in (i+1):length(groups)){
      g1 <- groups[i]
      g2 <- groups[j]
      markers <- tryCatch({FindMarkers(seuratdata,ident.1 = g1,ident.2 = g2,test.use = "wilcox")},
                          error = function(e) {message("ERROR: ", g1, " vs ", g2, " | ", e$message);return(NULL)})
      if(is.null(markers)) next
      markers$group1 <- g1
      markers$group2 <- g2
      markers$mean1 <- avg_exp[rownames(markers), g1]
      markers$mean2 <- avg_exp[rownames(markers), g2]
      markers$cell <- rownames(markers)
      rownames(markers) <- seq_len(nrow(markers))
      all_results[[idx]] <- markers
      idx <- idx + 1
    }
  }
  all_results_df <- do.call(rbind,all_results)
  write.csv(all_results_df,file = file.path(outdir,"CellEnriched.txt"),quote = FALSE)
  message("Analysis completed successfully!")
}

#===============================================================================
# Run
#===============================================================================

tryCatch({

  SpatialDeconAnalysis(
    inputfile = opt$input,
    clinicalfile = opt$clinical,
    outdir = opt$outdir
  )

}, error = function(e){

  message("")
  message("ERROR:")
  message(e$message)

  quit(status = 1)

})

quit(status = 0)