#!/bin/bash
set -euo pipefail


#1, Data Download
#https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1921
wget "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/921/S-BIAD1921/Files/GEO_oHSV_GBM_submission_20260109/supplementary/GBM_oHSV_seurat_CODEX_Full_preTX_20260114.rds" -P "../RawData/PMID41679299/"
wget "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/921/S-BIAD1921/Files/GEO_oHSV_GBM_submission_20260109/supplementary/GBM_oHSV_seurat_CODEX_Full_postTX_20260114.rds" -P "../RawData/PMID41679299/"

#2. Generate protein expression, spatial, samples and annotation files 
Rscript -e '
    library(Seurat)
    library(jsonlite)
    data1 = readRDS("../RawData/PMID41679299/GBM_oHSV_seurat_CODEX_Full_preTX_20260114.rds")
    data2 = readRDS("../RawData/PMID41679299/GBM_oHSV_seurat_CODEX_Full_postTX_20260114.rds")
    layer_names <- Layers(data1[["MFI"]])
    mat_list <- lapply(layer_names, function(layer){
        mat <- LayerData(object = data1,assay = "MFI",layer = layer)
        mat <- t(as.matrix(mat))    
        return(mat)
    })
    mat1 <- do.call(rbind, mat_list)
    mat1 <- mat1[, -which(colnames(mat1) == "DAPI-08")]
    mat2 <- GetAssayData(object = data2, assay = "MFI", layer = "counts")
    mat2 <- t(mat2)
    mat2 <- mat2[, -which(colnames(mat2) == "DAPI-03")]
    mat = rbind(mat1,mat2)
    protein_avg <- mat[, !grepl("Blank", colnames(mat))]
    
    protein_map <- fromJSON("proteinmap.json")
    alias_to_canonical <- c()
    for (canonical in names(protein_map)) {
      alias_to_canonical[canonical] <- canonical
      aliases <- protein_map[[canonical]]
      for (alias in aliases) {
        alias_to_canonical[alias] <- canonical
      }
    }
    colnames(protein_avg) <- ifelse(colnames(protein_avg) %in% names(alias_to_canonical), alias_to_canonical[colnames(protein_avg)], colnames(protein_avg))
    write.csv(as.matrix(protein_avg),"../RawData/PMID41679299/protein.csv",quote = FALSE)     
    samples <- data.frame(sample_id = c(data1@meta.data$orig.ident,data2@meta.data$orig.ident),row.names = c(rownames(data1@meta.data),rownames(data2@meta.data)))
    write.csv(samples,"../RawData/PMID41679299/samples.csv",quote = FALSE)  
    spatial <- data.frame(X = c(data1@meta.data$x,data2@meta.data$x),Y = c(data1@meta.data$y,data2@meta.data$y),row.names = c(rownames(data1@meta.data),rownames(data2@meta.data)))
    write.csv(spatial,"../RawData/PMID41679299/spatial.csv",quote = FALSE)  
    annotation <- data.frame(celltype = c(data1@meta.data$ann_cor,data2@meta.data$ann_cor),row.names = c(rownames(data1@meta.data),rownames(data2@meta.data)))
    write.csv(annotation,"../RawData/PMID41679299/annotation.csv",quote = FALSE)    
'

python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
protein_avg = pd.read_csv("../RawData/PMID41679299/protein.csv",index_col=0)
protein_avg = np.arcsinh(protein_avg)
protein_avg.to_csv("../RawData/PMID41679299/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID41679299/protein_avg.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID41679299/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "GBM\tBrain\tCODEX\tHuman\tTSOP%04d\tPMID41679299\t%s\n", FNR+num-1, $0}' >> SPsamples.txt

#3. Cell type annotaion
Rscript -e '
    library(ggplot2)
    library(RANN)
    library(Rcpp)
    sourceCpp("jaccard_coeff.cpp")
    library(igraph)
    source("functions_clustering.R")
    library(ggplot2)
    library(RANN)
    source("Rphenograph.R")
    xdata = read.csv("../RawData/PMID41679299/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","CD4","CD8A","PECAM1","CD68","CD163","HLA-DRA","LAMP3","CEACAM8","NCAM1","GFAP")]   
   
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID41679299/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(79L,68L),CD4T=c(31L,64L),CD8T=c(29L,22L,53L,44L,83L,42L,59L),Endothelial=c(56L,74L),Macrophage=c(17L),Dendritic=c(36L),Neutrophils=c(40L,9L,45L,34L),NK=c(10L),Tumor=c(35L,43L,89L))    
    good_clusters <- as.character(unlist(maps)) #define according to your running
    GC_cellIDs <- names(clusters[clusters %in% good_clusters])
    good_clusters_mean <- get_mean_intensity(expr = data[GC_cellIDs,],cell_clustering = clusters[clusters %in% good_clusters])
    good_clusters_mean <- as.matrix(good_clusters_mean)
    rownames(good_clusters_mean) <- good_clusters_mean[,'cell_clustering']
    good_clusters_mean <- good_clusters_mean[,-1]
    new_clust.km <- kmeans(x = data, centers = good_clusters_mean)
    mapped_cluster <- rownames(good_clusters_mean)[new_clust.km$cluster]
    named_clusters <- merge_clusters_alpha(cl = mapped_cluster, maplist = maps)
    dist <- get_nearest_center_with_dist(X = data, M = named_clusters,include.next.dist = TRUE)
    dist$cluster2 <- dist$cluster
    dist$cluster2[dist$dist.next.ratio > .97] <- 'unclassified'
    annotation <- data.frame(celltype=dist$cluster2)
    rownames(annotation) <- rownames(dist)    
    write.csv(annotation, file = "../RawData/PMID41679299/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID41679299/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID41679299/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID41679299/bubble.svg",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "B": ["MS4A1"],
    "T": ["CD3E"],
    "Macrophage": ["CD163"],
    "Dendritic": ["LAMP3"],
    "Neutrophils": ["CEACAM8"],
    "NK": ["NCAM1"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["GFAP"]
  }
}

import json
with open("../RawData/PMID41679299/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID41679299/protein_avg_nozscore.csv --astir_marker ../RawData/PMID41679299/astir_marker.json --output ../RawData/PMID41679299/ --anno_file ../RawData/PMID41679299/annotation.csv > ../RawData/PMID41679299/log.txt 2>&1


Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID41679299/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","CD163","HLA-DRA","LAMP3","CEACAM8","NCAM1","GFAP")]
    annotation = read.csv("../RawData/PMID41679299/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    write.csv(annotation,"../RawData/PMID41679299/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID41679299/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID41679299/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID41679299/bubble_astir.pdf ../RawData/PMID41679299/bubble_astir.svg
cp  ../RawData/PMID41679299/annotation_astir.svg ../RawData/PMID41679299/annotation_final.svg
cp  ../RawData/PMID41679299/bubble_astir.svg ../RawData/PMID41679299/bubble_final.svg


#4. Divide samples
head -n 1 ../RawData/PMID41679299/samples.csv  > ../RawData/PMID41679299/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID41679299/samples.csv  >> ../RawData/PMID41679299/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID41679299/ --output_dir ../SPOutput/