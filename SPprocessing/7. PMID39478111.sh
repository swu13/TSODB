#!/bin/bash
set -euo pipefail


#1, Data Download
#CODEX https://humantumoratlas.org/publications/hta1_2024_nature-medicine_johanna-klughammer?tab=codex (h5ad)
#https://singlecell.broadinstitute.org/single_cell/data/public/SCP2702/htapp-mbc?filename=codex.h5ad

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from scipy import stats
import json
adata = sc.read_h5ad("../RawData/PMID39478111/codex.h5ad")
protein_avg = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
protein_avg.columns = [alias_to_canonical.get(c, c) for c in protein_avg.columns]
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg.to_csv("../RawData/PMID39478111/protein_avg.csv")

protein_raw = protein_avg - protein_avg.min(axis=0)
protein_raw.to_csv("../RawData/PMID39478111/protein_avg_nonnegative.csv")

adata.var_names=[alias_to_canonical.get(c, c) for c in adata.var_names]
protein_avg = protein_avg.reindex(columns=adata.var_names)
mean = adata.var["mean"]
std = adata.var["std"]
protein_raw = protein_avg * std.values + mean.values

spatial = adata.obs[['x', 'y']].copy()
spatial.columns = ['X', 'Y']
spatial.to_csv("../RawData/PMID39478111/spatial.csv")
samples = adata.obs[['sample']].copy()
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID39478111/samples.csv")
annotation=adata.obs[['RCTD']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID39478111/annotation.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID39478111/samples.csv | tail -n +2 | sort | uniq | awk -F"\t" -v num="$num" '{printf "BRCA\tBreast\tCODEX\tHuman\tTSOP%04d\tPMID39478111\t%s\n", FNR+num-1, $0}' >> SPsamples.txt  


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
    xdata = read.csv("../RawData/PMID39478111/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("CD19","CD38","CD2","CD4","FOXP3","CD8A","PECAM1","CD1C","ITGAX","FUT4","B3GAT1","NCAM1","KIT","PDPN","KRT19")]   
    
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID39478111/test.pdf", width = 12, height = 15)
    maps=list(B=c(49L),Plasma=c(23L),CD4T=c(21L),CD8T=c(29L,3L),Endothelial=c(35L,33L),Dendritic=c(4L),Neutrophils=c(31L,37L),NK=c(50L),Mast=c(38L),Stromal=c(40L),Tumor=c(20L,46L,42L))
    
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
    write.csv(annotation, file = "../RawData/PMID39478111/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID39478111/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39478111/annotation.pdf", width = 15, height = 7)
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID39478111/bubble.svg",width = 20,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "B": ["CD19"],
    "Plasma": ["CD38"],
    "CD4T": ["CD4"],
    "CD8T": ["CD8A"],
    "Macrophage": ["CD68"],
    "Dendritic": ["ITGAX"],
    "Neutrophils": ["FUT4"],
    "NK": ["B3GAT1"],
    "Mast": ["KIT"],
    "Endothelial": ["PECAM1"],
    "Stromal": ["FAP"],
    "Tumor": ["KRT19"]
  }
}

import json
with open("../RawData/PMID39478111/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID39478111/protein_avg_nonnegative.csv --astir_marker ../RawData/PMID39478111/astir_marker.json --output ../RawData/PMID39478111/ --anno_file ../RawData/PMID39478111/annotation.csv > ../RawData/PMID39478111/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID39478111/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("CD19","CD38","CD2","CD8A","PECAM1","HLA-DRA","ITGAX","FUT4","B3GAT1","KIT","FAP","PDPN","KRT19")]   
    annotation = read.csv("../RawData/PMID39478111/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    write.csv(annotation,"../RawData/PMID39478111/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID39478111/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID39478111/bubble_astir.pdf",width = 12,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID39478111/bubble_astir.pdf ../RawData/PMID39478111/bubble_astir.svg
cp  ../RawData/PMID39478111/annotation_astir.svg ../RawData/PMID39478111/annotation_final.svg
cp  ../RawData/PMID39478111/bubble_astir.svg ../RawData/PMID39478111/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID39478111/samples.csv  > ../RawData/PMID39478111/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID39478111/samples.csv  >> ../RawData/PMID39478111/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID39478111/ --output_dir ../SPOutput/