#!/bin/bash
set -euo pipefail


#1, Data Download
#https://10.5281/zenodo.14531275

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data = pd.read_csv('../RawData/PMID35437329/SingleCells.csv')
data = data[data['is_tumour']==1]
protein_avg = data.iloc[:,11:48]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
    
protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]  
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg = np.arcsinh(protein_avg)
protein_avg.to_csv("../RawData/PMID35437329/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID35437329/protein_avg.csv")
spatial = data[['Location_Center_X','Location_Center_Y']]
spatial.columns=['X','Y']
spatial.to_csv("../RawData/PMID35437329/spatial.csv")
annotation = data[['cellPhenotype']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID35437329/annotation.csv")
samples = data[['metabric_id']]
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID35437329/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID35437329/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "BRCA\tBreast\tIMC\tHuman\tTSOP%04d\tPMID35437329\t%s\n", FNR+num-1, $0}' >> SPsamples.txt

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
    xdata = read.csv("../RawData/PMID35437329/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD38","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","CD163","HLA-DRA","ITGAX","FUT4","B3GAT1","ACTA2","PanCytoK")]      
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID35437329/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(24L,31L),Plasma=c(14L,16L),CD4T=c(4L),Treg=c(17L),CD8T=c(7L,27L),Endothelial=c(5L),Macrophage=c(3L),Dendritic=c(32L),Neutrophils=c(33L,28L),NK=c(37L),Stromal=c(10L),Tumor=c(42L))
    
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
    write.csv(annotation, file = "../RawData/PMID35437329/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID35437329/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID35437329/annotation.pdf")
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID35437329/bubble.svg",width = 10,height = 4)
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
    "Plasma": ["CD38"],
    "T": ["CD3E"],
    "Macrophage": ["CD163"],
    "Dendritic": ["ITGAX"],
    "NK": ["B3GAT1"],
    "Stromal": ["ACTA2"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/PMID35437329/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID35437329/protein_avg_nozscore.csv --astir_marker ../RawData/PMID35437329/astir_marker.json --output ../RawData/PMID35437329/ --anno_file ../RawData/PMID35437329/annotation.csv > ../RawData/PMID35437329/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID35437329/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD38","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","CD163","HLA-DRA","ITGAX","B3GAT1","ACTA2","PanCytoK")]   
    annotation = read.csv("../RawData/PMID35437329/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID35437329/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID35437329/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID35437329/bubble_astir.pdf",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID35437329/bubble_astir.pdf ../RawData/PMID35437329/bubble_astir.svg
cp  ../RawData/PMID35437329/annotation_astir.svg ../RawData/PMID35437329/annotation_final.svg
cp  ../RawData/PMID35437329/bubble_astir.svg ../RawData/PMID35437329/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID35437329/samples.csv  > ../RawData/PMID35437329/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID35437329/samples.csv  >> ../RawData/PMID35437329/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID35437329/ --output_dir ../SPOutput/