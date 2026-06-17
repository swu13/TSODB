#!/bin/bash
set -euo pipefail


#1, Data Download
#https://doi.org/10.5281/zenodo.5903179£¬ Clinical data (PMID35363543.txt) from Supplementary Tables

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data = pd.read_csv('../RawData/PMID35363543/ICI_meanIntensity.tsv',sep="\t")
protein_avg =  data.iloc[:,np.r_[4:33, 36:42]]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical

protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]  #S100 keep
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg.to_csv("../RawData/PMID35363543/protein_avg_nozscore.csv")
protein_avg = np.arcsinh(protein_avg)
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID35363543/protein_avg.csv")
spatial = data[['Location_Center_X','Location_Center_Y']]
spatial.columns=['X','Y']
spatial.to_csv("../RawData/PMID35363543/spatial.csv")
annotation = data[['Cluster.v2']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID35363543/annotation.csv")
samples = data[['Sample_ID']]
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID35363543/samples.csv")
EOF


#samples
cut -d"," -f2 ../RawData/PMID35363543/samples.csv | tail -n +2 |  sort | uniq > ../RawData/PMID35363543/samples.txt
num=`wc -l SPsamples.txt | cut -d" " -f1`
awk -F"\t" -v num="$num" 'ARGIND==1{samples[$5]=$1"\t"$2"\t"$3"\t"$4}ARGIND==2{printf samples[$0]"\tTSOP%04d\tPMID35363543\t%s\n", FNR+num-1, $0}' ../RawData/PMID35363543/PMID35363543.txt ../RawData/PMID35363543/samples.txt >> SPsamples.txt
rm -f ../RawData/PMID35363543/samples.txt


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
    xdata = read.csv("../RawData/PMID35363543/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","CD14","GZMB","ACTA2","SOX10","S100")]  
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID35363543/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(12L),CD4T=c(3L),Treg=c(6L),CD8T=c(7L),Endothelial=c(1L),Macrophage=c(5L),Monocyte=c(18L),NK=c(9L),Stromal=c(2L),Tumor=c(25L,14L))
    
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
    write.csv(annotation, file = "../RawData/PMID35363543/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID35363543/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID35363543/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID35363543/bubble.svg",width = 10,height = 4)
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
    "CD4T": ["CD4"],
    "CD8T": ["CD8A"],
    "Macrophage": ["CD68"],
    "Monocyte": ["CD14"],
    "NK": ["GZMB"],
    "Stromal": ["ACTA2"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["S100"]
  }
}
#SOX10
import json
with open("../RawData/PMID35363543/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID35363543/protein_avg_nozscore.csv --astir_marker ../RawData/PMID35363543/astir_marker.json --output ../RawData/PMID35363543/ --anno_file ../RawData/PMID35363543/annotation.csv > ../RawData/PMID35363543/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID35363543/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","CD14","ACTA2","GZMB","SOX10","S100")]  
    annotation = read.csv("../RawData/PMID35363543/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    NKcells = rownames(subset(annotation,celltype == "NK"))
    annotation[NKcells, "celltype"]="unclassified"
    Stomalcells = rownames(subset(annotation,celltype == "Stromal"))
    annotation[Stomalcells, "celltype"]="unclassified"
    
    unclassified <- rownames(subset(annotation,celltype == "unclassified")) 
    markers <- c(
      "B" = "MS4A1",
      "CD4T" = "CD4",
      "CD8T" = "CD8A",
      "Macrophage" = "CD68",
      "Monocyte" = "CD14",
      "NK" = "GZMB",
      "Stromal" = "ACTA2",
      "Endothelial" = "PECAM1",
      "Tumor" = "S100"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID35363543/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID35363543/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID35363543/bubble_astir.pdf",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID35363543/bubble_astir.pdf ../RawData/PMID35363543/bubble_astir.svg
cp  ../RawData/PMID35363543/annotation_astir.svg ../RawData/PMID35363543/annotation_final.svg
cp  ../RawData/PMID35363543/bubble_astir.svg ../RawData/PMID35363543/bubble_final.svg


#4. Divide samples
head -n 1 ../RawData/PMID35363543/samples.csv  > ../RawData/PMID35363543/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID35363543/samples.csv  >> ../RawData/PMID35363543/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID35363543/ --output_dir ../SPOutput/