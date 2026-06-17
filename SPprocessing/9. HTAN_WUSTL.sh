#!/bin/bash
set -euo pipefail


#1, Data Download
#https://humantumoratlas.org/
unzip ../RawData/HTAN_WUSTL/HTAN_WUSTL.zip -d ../RawData/HTAN_WUSTL/

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import os
import glob
import json
input_dir = "../RawData/HTAN_WUSTL/HTAN_WUSTL/"
files = glob.glob(os.path.join(input_dir, "*.txt"))
samples = (pd.Series(files).apply(os.path.basename).str.extract(r'(u\d+)\.txt$')[0].tolist())
#merge
protein_tables = []
xy_table = []
annotation_table = []
sample_table = []
for i, f in enumerate(files):
    df = pd.read_csv(f, sep="\t")
    protein_cols = [c for c in df.columns if c.endswith("_intensity") and not c.endswith("_intensity_scaled")]
    protein_avg = df[protein_cols].copy()
    protein_avg.columns = [c.replace("_intensity", "") for c in protein_avg.columns]
    cells = [f"{samples[i]}_{idx}" for idx in protein_avg.index]
    protein_avg.index=cells
    protein_tables.append(protein_avg)
    spatial = df[['row','col']]
    spatial.index=cells
    spatial.columns = ['X','Y']
    xy_table.append(spatial)
    annotation = df[['cell_type']]
    annotation.index=cells
    annotation.columns = ['celltype']
    annotation_table.append(annotation)
    sample_df = pd.DataFrame({'sample_id': samples[i]}, index=cells)
    sample_table.append(sample_df)

merged_protein_tables = pd.concat(protein_tables,axis=0).fillna(0)
merged_xy_table = pd.concat(xy_table, axis=0)
merged_annotation_table = pd.concat(annotation_table, axis=0)
merged_sample_table = pd.concat(sample_table, axis=0)
merged_protein_tables = merged_protein_tables.drop(columns=["DAPI"])

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
merged_protein_tables.columns = [alias_to_canonical.get(c, c) for c in merged_protein_tables.columns]

merged_protein_tables = np.arcsinh(merged_protein_tables)
merged_protein_tables.to_csv("../RawData/HTAN_WUSTL/protein_avg_nozscore.csv")
merged_protein_tables_zscore = stats.zscore(merged_protein_tables)
merged_protein_tables = pd.DataFrame(merged_protein_tables_zscore, columns=merged_protein_tables.columns, index=merged_protein_tables.index)
merged_protein_tables.to_csv("../RawData/HTAN_WUSTL/protein_avg.csv")
merged_xy_table.to_csv("../RawData/HTAN_WUSTL/spatial.csv")
merged_annotation_table.to_csv("../RawData/HTAN_WUSTL/annotation.csv")
merged_sample_table.to_csv("../RawData/HTAN_WUSTL/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/HTAN_WUSTL/samples.csv | tail -n +2 | sort | uniq | awk -F"\t" -v num="$num" '{printf "BRCA\tBreast\tCODEX\tHuman\tTSOP%04d\tHTAN_WUSTL\t%s\n", FNR+num-1, $0}' >> SPsamples.txt  

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
    xdata = read.csv("../RawData/HTAN_WUSTL/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","KIT","ACTA2","PanCytoK","KRT19")]   
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/HTAN_WUSTL/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(47L),CD4T=c(15L,28L,16L),CD8T=c(10L,12L),Endothelial=c(48L,30L,9L,13L),Macrophage=c(31L),Mast=c(20L,34L),Stromal=c(35L,50L,19L),Tumor=c(21L,3L))
    
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
    write.csv(annotation, file = "../RawData/HTAN_WUSTL/celltype_annotation.txt",quote = FALSE)
    
    svg("../RawData/HTAN_WUSTL/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/HTAN_WUSTL/bubble.svg",width = 8,height = 4)
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
    "Macrophage": ["CD68"],
    "Stromal": ["ACTA2"],
    "Mast": ["KIT"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/HTAN_WUSTL/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/HTAN_WUSTL/protein_avg_nozscore.csv --astir_marker ../RawData/HTAN_WUSTL/astir_marker.json --output ../RawData/HTAN_WUSTL/ --anno_file ../RawData/HTAN_WUSTL/annotation.csv > ../RawData/HTAN_WUSTL/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/HTAN_WUSTL/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","KIT","ACTA2","PanCytoK","KRT19")] 
    annotation = read.csv("../RawData/HTAN_WUSTL/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    Maccells = rownames(subset(annotation,celltype == "Macrophage"))
    annotation[Maccells, "celltype"]="unclassified"
    Mascells = rownames(subset(annotation,celltype == "Mast"))
    annotation[Mascells, "celltype"]="unclassified"
    
    unclassified <- rownames(subset(annotation,celltype == "unclassified")) 
    markers <- c(
      "B" = "MS4A1",
      "T" = "CD3E",
      "Macrophage" = "CD68",
      "Stromal" = "ACTA2",
      "Mast" = "KIT",
      "Endothelial" = "PECAM1",
      "Tumor" = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD8A"]>0,"CD8T","CD4T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/HTAN_WUSTL/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/HTAN_WUSTL/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/HTAN_WUSTL/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/HTAN_WUSTL/bubble_astir.pdf ../RawData/HTAN_WUSTL/bubble_astir.svg
cp  ../RawData/HTAN_WUSTL/annotation_astir.svg ../RawData/HTAN_WUSTL/annotation_final.svg
cp  ../RawData/HTAN_WUSTL/bubble_astir.svg ../RawData/HTAN_WUSTL/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/HTAN_WUSTL/samples.csv  > ../RawData/HTAN_WUSTL/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/HTAN_WUSTL/samples.csv  >> ../RawData/HTAN_WUSTL/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/HTAN_WUSTL/ --output_dir ../SPOutput/