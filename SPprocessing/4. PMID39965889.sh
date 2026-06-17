#!/bin/bash
set -euo pipefail


#1, Data Download
#https://zenodo.org/records/13483663
unzip ../RawData/PMID39965889/CODEX.zip -d ../RawData/PMID39965889/
cp ../RawData/PMID39965889/CODEX/*/*.csv ../RawData/PMID39965889/
rm -f ../RawData/PMID39965889/CODEX.zip
cut -d"," -f3-4 "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/A7_Post_Tumor_Core_Analysis.csv" --complement > "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/A7_Post_Tumor_Core_Analysis_update.csv"
rm -f "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/A7_Post_Tumor_Core_Analysis.csv"
cut -d"," -f3 "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/B2_Post_Tumor_Core_Analysis.csv" --complement > "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/B2_Post_Tumor_Core_Analysis_update.csv"
rm -f "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889/B2_Post_Tumor_Core_Analysis.csv"
rm -f ../RawData/PMID39965889/*Non_Tumor*

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import re
import numpy as np
from scipy import stats
import glob
import os
import json
indir = "/data_d/WSJ/SpatialMetsDB/RawData/PMID39965889"
exclude_files = {"protein_avg.csv","samples.csv","spatial.csv","annotation.csv"}

files = [
    f for f in glob.glob(os.path.join(indir, "*.csv"))
    if os.path.basename(f) not in exclude_files
]
samples = [os.path.basename(f).replace("_Analysis.csv", "")for f in files]
#identical columns
ref = pd.read_csv(files[0], sep=",")
for f in files:
    df = pd.read_csv(f, sep=",")
    if not ref.columns.equals(df.columns):
        print(f)
#merge
all_data = []
for i, f in enumerate(files):
    data = pd.read_csv(f)
    cellid = (samples[i] + "_" + data.index.astype(str))
    data.index = cellid
    all_data.append(data)
merged = pd.concat(all_data, axis=0)
intensity_cols = [c for c in merged.columns if c.endswith("Cell Intensity")]
protein_avg = merged[intensity_cols].copy()
protein_avg.columns = [re.sub(r" Cell Intensity$", "", c) for c in intensity_cols]
protein_avg = protein_avg.filter(regex='^(?!DAPI)')

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
protein_avg.columns = [alias_to_canonical.get(c, c) for c in protein_avg.columns]
protein_avg = protein_avg.T.groupby(level=0).mean().T

protein_avg = np.arcsinh(protein_avg)
protein_avg.to_csv("../RawData/PMID39965889/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID39965889/protein_avg.csv")
merged["X"] = (merged["XMin"] + merged["XMax"]) / 2
merged["Y"] = (merged["YMin"] + merged["YMax"]) / 2
spatial = merged[['X','Y']]
spatial.to_csv("../RawData/PMID39965889/spatial.csv")
merged["sample_id"] = merged["Pt"].astype(str) + "_" + merged["Timepoint"].astype(str)
samples = merged["sample_id"]
samples.to_csv("../RawData/PMID39965889/samples.csv")
merged["celltype"] = merged.iloc[:, 7:23].idxmax(axis=1)
annotation = merged["celltype"]
annotation.to_csv("../RawData/PMID39965889/annotation.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID39965889/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "LICA\tLiver\tCODEX\tHuman\tTSOP%04d\tPMID39965889\t%s\n", FNR+num-1, $0}' >> SPsamples.txt


#3. Cell type annotaion  ??? More
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
    xdata = read.csv("../RawData/PMID39965889/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","CD34","CD68","CD163","LYVE1","HLA-DRA","ITGAX","CD86","MPO","FUT4","B3GAT1","ACTA2","PanCytoK","CPS1")]   
    
    pheno.out <- Rphenograph(data, k = 20)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID39965889/test.pdf", width = 12, height = 15)
    good_clusters <- c("18", "7","8", "12", "34", "4","17","49", "11", "33","10", "43", "13","9","28","41","52","44","31","47") #define according to your running
    GC_cellIDs <- names(clusters[clusters %in% good_clusters])
    good_clusters_mean <- get_mean_intensity(expr = data[GC_cellIDs,],cell_clustering = clusters[clusters %in% good_clusters])
    good_clusters_mean <- as.matrix(good_clusters_mean)
    rownames(good_clusters_mean) <- good_clusters_mean[,'cell_clustering']
    good_clusters_mean <- good_clusters_mean[,-1]
    new_clust.km <- kmeans(x = data, centers = good_clusters_mean)
    mapped_cluster <- rownames(good_clusters_mean)[new_clust.km$cluster]
    named_clusters <- merge_clusters_alpha(cl = mapped_cluster, maplist = list(B = c(18L),Treg = c(7L,8L),CD4T = 12L,CD8T = 34L,Endothelial = c(4L, 17L, 49L),Macrophage = c(11L,33L),Dendritic = c(10L, 43L),Neutrophils = c(13L, 9L),NK = 28L,Stromal = 41L,Tumor = c(52L,44L,31L,47L)))
    dist <- get_nearest_center_with_dist(X = data, M = named_clusters,include.next.dist = TRUE)
    dist$cluster2 <- dist$cluster
    dist$cluster2[dist$dist.next.ratio > .97] <- 'unclassified'
    annotation <- data.frame(celltype=dist$cluster2)
    rownames(annotation) <- rownames(dist)    
    write.csv(annotation, file = "../RawData/PMID39965889/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID39965889/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39965889/annotation.pdf", width = 15, height = 7)
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression),.groups = "drop")
    svg("../RawData/PMID39965889/bubble.svg",width = 8,height = 4)
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
    "Dendritic": ["ITGAX"],
    "Neutrophils": ["FUT4"],
    "NK": ["B3GAT1"],
    "Endothelial": ["CD34"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/PMID39965889/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID39965889/protein_avg_nozscore.csv --astir_marker ../RawData/PMID39965889/astir_marker.json --output ../RawData/PMID39965889/ --anno_file ../RawData/PMID39965889/annotation.csv > ../RawData/PMID39965889/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID39965889/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","CD34","CD68","CD163","HLA-DRA","ITGAX","FUT4","B3GAT1","ACTA2","PanCytoK")]   
    annotation = read.csv("../RawData/PMID39965889/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("Dendritic"))) 
    markers <- c(
      "B" = "MS4A1",
      "T" = "CD3E",
      "Macrophage" = "CD163",
      "Dendritic" = "ITGAX",
      "Neutrophils" = "FUT4",
      "NK" = "B3GAT1",
      "Endothelial" = "CD34",
      "Stromal" = "ACTA2",
      "Tumor" = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID39965889/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID39965889/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID39965889/bubble_astir.pdf",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID39965889/bubble_astir.pdf ../RawData/PMID39965889/bubble_astir.svg
cp  ../RawData/PMID39965889/annotation_astir.svg ../RawData/PMID39965889/annotation_final.svg
cp  ../RawData/PMID39965889/bubble_astir.svg ../RawData/PMID39965889/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID39965889/samples.csv  > ../RawData/PMID39965889/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID39965889/samples.csv  >> ../RawData/PMID39965889/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID39965889/ --output_dir ../SPOutput/