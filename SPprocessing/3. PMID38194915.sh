#!/bin/bash
set -euo pipefail


#1, Data Download
#https://doi.org/10.5281/zenodo.10045066

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data = pd.read_csv('../RawData/PMID38194915/TNBC_immune_CODEX_counts.csv',index_col=0)
data = data.filter(regex='^(?!DAPI)')
parts = pd.Series(data.columns.to_list()).str.split('_', expand=True)
data.columns = pd.MultiIndex.from_arrays([parts[0], parts[1]], names=['protein', 'compartment'])
protein_avg = data.T.groupby(level='protein').mean().T

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
protein_avg.columns = [alias_to_canonical.get(c, c) for c in protein_avg.columns]

protein_avg = np.arcsinh(protein_avg)
protein_avg.to_csv("../RawData/PMID38194915/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID38194915/protein_avg.csv")
metadata = pd.read_csv('../RawData/PMID38194915/TNBC_immune_CODEX.csv',index_col=0)
spatial = metadata[['x','y']]
spatial.columns = ['X','Y']
spatial.to_csv("../RawData/PMID38194915/spatial.csv")
samples = metadata[['patient']]
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID38194915/samples.csv")
annotation = metadata[['celltype']]
annotation.to_csv("../RawData/PMID38194915/annotation.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID38194915/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "BRCA\tBreast\tCODEX\tHuman\tTSOP%04d\tPMID38194915\t%s\n", FNR+num-1, $0}' >> SPsamples.txt


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
    xdata = read.csv("../RawData/PMID38194915/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","CD8A","PECAM1","CD68","HLA-DRA","ITGAX","ACTA2","PanCytoK")]     
    
    set.seed(1234)  
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID38194915/test.pdf", width = 12, height = 15)

    maps=list(B=c(35L,65L,98L),Plasma=c(90L,107L),CD4T=c(15L,23L),CD8T=c(58L,2L),Endothelial=c(17L,7L,51L,105L),Macrophage=c(13L),Dendritic=c(27L,11L,74L),Stromal=c(75L,103L,29L),Tumor=c(3L,46L))
    
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
    write.csv(annotation, file = "../RawData/PMID38194915/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID38194915/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID38194915/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID38194915/bubble.svg",width = 8,height = 4)
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
    "Plasma": ["SDC1"],
    "T": ["CD3E"],
    "Macrophage": ["CD68"],
    "Dendritic": ["ITGAX"],
    "Endothelial": ["PECAM1"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/PMID38194915/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID38194915/protein_avg_nozscore.csv --astir_marker ../RawData/PMID38194915/astir_marker.json --output ../RawData/PMID38194915/ --anno_file ../RawData/PMID38194915/annotation.csv > ../RawData/PMID38194915/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID38194915/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","ITGAX","ACTA2","PanCytoK")] 
    annotation = read.csv("../RawData/PMID38194915/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    write.csv(annotation,"../RawData/PMID38194915/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID38194915/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID38194915/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID38194915/bubble_astir.pdf ../RawData/PMID38194915/bubble_astir.svg
cp  ../RawData/PMID38194915/annotation_astir.svg ../RawData/PMID38194915/annotation_final.svg
cp  ../RawData/PMID38194915/bubble_astir.svg ../RawData/PMID38194915/bubble_final.svg


#4. Divide samples
head -n 1 ../RawData/PMID38194915/samples.csv  > ../RawData/PMID38194915/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID38194915/samples.csv  >> ../RawData/PMID38194915/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID38194915/ --output_dir ../SPOutput/