#!/bin/bash
set -euo pipefail


#1, Data Download
#https://humantumoratlas.org/
#CyCIF data BRCA (PMID40404617(6 samples, same as BioRxiv557464),PMID35243422(2 samples, same as BioRxiv557464),BioRxiv557464(32 samples)), rel7_imaging_level_4(67 samples)
unzip ../RawData/HTAN_OHSU/BioRxiv557464.zip -d ../RawData/HTAN_OHSU/
unzip ../RawData/HTAN_OHSU/rel7_imaging_level_4.zip -d ../RawData/HTAN_OHSU/
cp ../RawData/HTAN_OHSU/*/*csv ../RawData/HTAN_OHSU/

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import os
import glob
import re
import json
input_dir = "../RawData/HTAN_OHSU/"
exclude_files = {"protein_avg.csv","samples.csv","spatial.csv","annotation.csv"}
files = [
    f for f in glob.glob(os.path.join(input_dir, "*.csv"))
    if os.path.basename(f) not in exclude_files
]
samples = [re.sub(r'_(cellring|cellringmask)\.csv$', '', os.path.basename(f)) for f in files]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
 
#merge       
df_tables = []
xy_table = []
sample_table = []
for i, f in enumerate(files):
    df = pd.read_csv(f)
    mapped_aliases = {c: alias_to_canonical[c] for c in df.columns if c in alias_to_canonical}
    df_mapped = df.loc[:, df.columns.isin(alias_to_canonical)]
    df_mapped.columns = [alias_to_canonical[c] for c in df_mapped.columns]
    cells = [f"{samples[i]}_{idx}" for idx in df.index]
    df_mapped.index=cells
    df_mapped = df_mapped.T.groupby(level=0).mean().T
    df_tables.append(df_mapped)
    spatial = df[['X_centroid','Y_centroid']]
    spatial.index=cells
    spatial.columns = ['X','Y']
    xy_table.append(spatial)
    sample_df = pd.DataFrame({'sample_id': samples[i]}, index=cells)
    sample_table.append(sample_df)

merged_protein_tables = pd.concat(df_tables,axis=0).fillna(0)
merged_xy_table = pd.concat(xy_table, axis=0)
merged_sample_table = pd.concat(sample_table, axis=0)
 
merged_protein_tables = np.arcsinh(merged_protein_tables)
merged_protein_tables.to_csv("../RawData/HTAN_OHSU/protein_avg_nozscore.csv")

merged_protein_tables_zscore = stats.zscore(merged_protein_tables)
merged_protein_tables = pd.DataFrame(merged_protein_tables_zscore, columns=merged_protein_tables.columns, index=merged_protein_tables.index)
merged_protein_tables.to_csv("../RawData/HTAN_OHSU/protein_avg.csv")
merged_xy_table.to_csv("../RawData/HTAN_OHSU/spatial.csv")
merged_sample_table.to_csv("../RawData/HTAN_OHSU/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/HTAN_OHSU/samples.csv | tail -n +2 | sort | uniq | awk -F"\t" -v num="$num" '{printf "BRCA\tBreast\tCyCIF\tHuman\tTSOP%04d\tHTAN_OHSU\t%s\n", FNR+num-1, $0}' >> SPsamples.txt 

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
    xdata = read.csv("../RawData/HTAN_OHSU/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("PTPRC","ACTA2","PanCytoK")]   
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/HTAN_OHSU/test.pdf", width = 12, height = 20)
    
    maps=list(Immune=c(61L,11L,109L,110L,112L),Stromal=c(24L,31L,72L),Tumor=c(127L,125L,123L,103L,119L))
        
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
    write.csv(annotation, file = "../RawData/HTAN_OHSU/celltype_annotation.txt",quote = FALSE)
    
    
    svg("../RawData/HTAN_OHSU/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/HTAN_OHSU/bubble.svg",width = 12,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "Immune": ["PTPRC"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK","KRT19"]
  }
}

astir_marker_dict={
  "cell_type": {
    "B": ["PTPRC","MS4A1"],
    "Plasma": ["PTPRC","SDC1"],
    "CD4T": ["CD3E","CD4"],
    "Treg": ["CD4","FOXP3"],
    "CD8T": ["CD3E","CD8A"],
    "Macrophage": ["CD68","CD163"],
    "Monocyte": ["ITGAM","CD14"],
    "Dendritic": ["HLA-DRA","ITGAX"],
    "Neutrophils": ["FUT4"],
    "NK": ["B3GAT1","NCR1"],
    "Endothelial": ["PECAM1","PDPN"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK","KRT19"]
  }
}


import json
with open("../RawData/HTAN_OHSU/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/HTAN_OHSU/protein_avg_nozscore.csv --astir_marker ../RawData/HTAN_OHSU/astir_marker.json --output ../RawData/HTAN_OHSU/ --anno_file ../RawData/HTAN_OHSU/annotation.csv > ../RawData/HTAN_OHSU/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/HTAN_OHSU/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("PTPRC","ACTA2","PanCytoK")]   
    annotation = read.csv("../RawData/HTAN_OHSU/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    write.csv(annotation,"../RawData/HTAN_OHSU/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/HTAN_OHSU/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/HTAN_OHSU/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf ../RawData/HTAN_OHSU/bubble_astir.pdf ../RawData/HTAN_OHSU/bubble_astir.svg
cp  ../RawData/HTAN_OHSU/annotation.svg ../RawData/HTAN_OHSU/annotation_final.svg
cp  ../RawData/HTAN_OHSU/bubble.svg ../RawData/HTAN_OHSU/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/HTAN_OHSU/samples.csv  > ../RawData/HTAN_OHSU/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/HTAN_OHSU/samples.csv  >> ../RawData/HTAN_OHSU/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/HTAN_OHSU/ --output_dir ../SPOutput/