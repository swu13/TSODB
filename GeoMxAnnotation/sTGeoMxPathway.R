#!/usr/bin/env Rscript

# Gene Set Score Differential Analysis (Hallmark + AddModuleScore)
# Input: seurat.rds (RNA assay + group metadata)
# Output: geneset differential results

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(msigdbr)
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
  description = "Hallmark gene set scoring and differential analysis"
)

opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

safeAddModuleScore <- function(
    object,
    features,
    name = "HALLMARK",
    ctrl_candidates = c(100, 50, 25, 10, 5, 1),
    seed = 1234,
    outdir = "."
) {

    for(ctrl in ctrl_candidates){
        cat("Trying ctrl =", ctrl, "\n")
        res <- tryCatch({
            set.seed(seed)
            AddModuleScore(object = object,features = features,name = name,ctrl = ctrl)
        }, error = function(e){
            message("Failed ctrl=", ctrl, ": ", e$message)
            NULL
        })
        if(!is.null(res)){
            cat("Success with ctrl =", ctrl, "\n")
            return(res)
        }
    }
    all_results_df <- data.frame(p_val = numeric(0),avg_log2FC = numeric(0),pct.1 = numeric(0),pct.2 = numeric(0),p_val_adj = numeric(0),group1 = character(0),group2 = character(0),
        mean1 = numeric(0),mean2 = numeric(0),stringsAsFactors = FALSE)
    write.csv(all_results_df,file = file.path(outdir,"PathwayEnriched.txt"),quote = FALSE)
    saveRDS(object, file.path(outdir, "seurat.pathway.rds"))
    stop("AddModuleScore failed for all ctrl values.")
}

#===============================================================================
# Main Function
#===============================================================================

GenesetAnalysis <- function(inputfile, outdir){

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: Load Seurat seuratdataect

  seuratdata <- tryCatch({
    readRDS(inputfile)
  }, error = function(e){
    stop(paste("Failed to read RDS:", e$message))
  })

  if (!"group" %in% colnames(seuratdata@meta.data)) {
    stop("metadata column 'group' not found")
  }

  DefaultAssay(seuratdata) <- "Spatial"

  # Step 2: Build Hallmark gene sets
  hallmark <- msigdbr(species = "Homo sapiens",category = "H")
  hallmark_list <- split(hallmark$gene_symbol,hallmark$gs_name)

  # Step 3: Compute module scores
  seed=1234
  set.seed(seed)
  seuratdata <- safeAddModuleScore(seuratdata,features = hallmark_list,name = "HALLMARK",seed=seed,outdir=outdir)
  score_cols <- grep("^HALLMARK",colnames(seuratdata@meta.data),value = TRUE)
  geneset_mat <- t(seuratdata@meta.data[, score_cols, drop = FALSE])
  rownames(geneset_mat) <- names(hallmark_list)
  seuratdata[["GeneSet"]] <- CreateAssayObject(counts = geneset_mat)
  DefaultAssay(seuratdata) <- "GeneSet"
  seuratdata@meta.data[, score_cols] <- NULL
  avg_exp <- AverageExpression(seuratdata,assays = DefaultAssay(seuratdata),slot = "count")[[1]]
  colnames(avg_exp) <- levels(factor(seuratdata$group))
  # Step 4: Differential analysis
  Idents(seuratdata) <- "group"
  groups <- unique(seuratdata$group)
  all_results <- list()
  idx <- 1
  meta <- seuratdata@meta.data
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      markers <- tryCatch({FindMarkers(seuratdata,ident.1 = g1,ident.2 = g2,test.use = "wilcox")},
                          error = function(e) {message("ERROR: ", g1, " vs ", g2, " | ", e$message);return(NULL)})
      if(is.null(markers)) next
      markers$group1 <- g1
      markers$group2 <- g2
      markers$mean1 <- avg_exp[rownames(markers), g1]
      markers$mean2 <- avg_exp[rownames(markers), g2]
      markers$pathway <- rownames(markers)
      rownames(markers) <- seq_len(nrow(markers))
      all_results[[idx]] <- markers
      idx <- idx + 1
    }
  }
  
  # Step 5: Save merged results (with empty-safe handling)
  all_results_df <- do.call(rbind,all_results)
  write.csv(all_results_df,file = file.path(outdir,"PathwayEnriched.txt"),quote = FALSE)
  saveRDS(seuratdata, file.path(outdir, "seurat.pathway.rds"))
  message("Analysis completed successfully!")
}

#===============================================================================
# Run
#===============================================================================

tryCatch({

  GenesetAnalysis(
    inputfile = opt$input,
    outdir = opt$outdir
  )

}, error = function(e){

  message("\nERROR:")
  message(e$message)

  quit(status = 1)

})

quit(status = 0)