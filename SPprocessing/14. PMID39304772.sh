#!/bin/bash
set -euo pipefail


#1, Data Download
#https://springernature.figshare.com/ndownloader/files/48352324, normalized: an inverse hyperbolic sine transform (¡®base::asinh¡¯) on cell expression values for every marker, in every region of interest (ROI). Next, normalized values were z-scaled across both cells and markers.

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data = pd.read_csv('../RawData/PMID39304772/Celllevel_rawdata_8_12_24.csv',index_col=0)
protein_avg = data.iloc[:,4:]
data["sample_id"] = data.index.str.split(".").str[0]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
    
protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]  
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg.to_csv("../RawData/PMID39304772/protein_avg.csv")

protein_raw = protein_avg - protein_avg.min(axis=0)
protein_raw.to_csv("../RawData/PMID39304772/protein_avg_nonnegative.csv")

spatial = data[['x','y']]
spatial.columns=['X','Y']
spatial.to_csv("../RawData/PMID39304772/spatial.csv")
annotation = data[['Celltype']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID39304772/annotation.csv")
data["sample_id"].to_csv("../RawData/PMID39304772/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID39304772/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "LICA\tLive\tCODEX\tHuman\tTSOP%04d\tPMID39304772\t%s\n", FNR+num-1, $0}' >> SPsamples.txt


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
    xdata = read.csv("../RawData/PMID39304772/protein_avg.csv",row.names=1,check.names=FALSE)        
    #data <- xdata[,c("MS4A1","CD19","CD3E","CD4","FOXP3","CD8A","PECAM1","PDPN","CD68","CD163","ITGAX","CD1C","CD14","FUT4","NCAM1","KIT","SIGLEC8","ACTA2","PanCytoK")]   
    data <- xdata[,c("CD3E","CD4","CD8A","PECAM1","CD163","CD1C","CD14","FUT4","NCAM1","KIT","ACTA2","PanCytoK")]     
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 10)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID39304772/test.pdf", width = 12, height = 20)
    
    NK=c(11L),Mast=c(50L,44L,18L)
    
    maps=list(T=c(1L,29L,19L),Endothelial=c(5L,25L,35L),Macrophage=c(74L,4L),Dendritic=c(61L,31L),Monocyte=c(69L,2L),Neutrophils=c(70L,63L,54L,8L,57L,26L),Stromal=c(24L,28L),Tumor=c(48L,40L))   
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
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39304772/annotation.pdf")
    
    write.csv(annotation, file = "../RawData/PMID39304772/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID39304772/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39304772/annotation.pdf")
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID39304772/bubble.svg",width = 12,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "T": ["CD3E"],
    "Mac/Mono": ["CD68"],
    "Dendritic": ["CD1C"],
    "Neutrophils": ["FUT4"],
    "Endothelial": ["PECAM1"],
    "Mast": ["KIT"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/PMID39304772/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID39304772/protein_avg_nonnegative.csv --astir_marker ../RawData/PMID39304772/astir_marker.json --output ../RawData/PMID39304772/ --anno_file ../RawData/PMID39304772/annotation.csv > ../RawData/PMID39304772/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID39304772/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","CD163","CD14","CD1C","FUT4","KIT","ACTA2","PanCytoK")]  
    annotation = read.csv("../RawData/PMID39304772/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("Mast"))) 
    markers <- c(
      "T" = "CD3E",
      "Mac/Mono" = "CD68",
      "Dendritic" = "CD1C",
      "Neutrophils" = "FUT4",
      "Endothelial" = "PECAM1",
      "Mast" = "KIT",
      "Stromal" = "ACTA2",
      "Tumor" = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    MMcells = rownames(subset(annotation,celltype == "Mac/Mono"))
    annotation[MMcells, "celltype"] = ifelse(xdata[MMcells,"CD163"]>xdata[MMcells,"CD14"],"Macrophage","Monocyte") 
    
    write.csv(annotation,"../RawData/PMID39304772/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID39304772/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    #plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39304772/annotation_astir.pdf")
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID39304772/bubble_astir.pdf",width = 12,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID39304772/bubble_astir.pdf ../RawData/PMID39304772/bubble_astir.svg
cp  ../RawData/PMID39304772/annotation_astir.svg ../RawData/PMID39304772/annotation_final.svg
cp  ../RawData/PMID39304772/bubble_astir.svg ../RawData/PMID39304772/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID39304772/samples.csv  > ../RawData/PMID39304772/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID39304772/samples.csv  >> ../RawData/PMID39304772/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID39304772/ --output_dir ../SPOutput/