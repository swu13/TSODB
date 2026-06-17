#!/bin/bash
set -euo pipefail


#1, Data Download
#MIBI https://www.synapse.org/Synapse:syn54117364 (csv)

#2. Generate protein expression, spatial, samples and annotation files 
python <<<EOF
import pandas as pd
import numpy as np
from scipy import stats
import json
data=pd.read_csv("../RawData/PMID39478111/MBC_all_Expression_ScaleSize.csv")
protein_avg = data.iloc[:,11:51]

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
protein_avg = pd.DataFrame(protein_avg, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID39478111/protein_avg_MIBI.csv")
spatial = data[['Cood_X', 'Cood_Y']].rename(columns={'Cood_X': 'X', 'Cood_Y': 'Y'})
spatial.to_csv("../RawData/PMID39478111/spatial_MIBI.csv")
samples = data[['pointNum']].rename(columns={'pointNum': 'sample_id'})
samples.to_csv("../RawData/PMID39478111/samples_MIBI.csv")
EOF


#samples
cut -d"," -f2 ../RawData/PMID39478111/samples_MIBI.csv | tail -n +2 | sort | uniq | awk -F"\t" '{print "BRCA\tBreast\tMIBI\tHuman\tPMID39478111_"$0}' >> SPsamples.txt 


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
    xdata = read.csv("../RawData/PMID39478111/protein_avg_MIBI.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","CD8A","CD34","CD68","CD163","CDH1","VIM")]   
    

    
    
    pheno.out <- Rphenograph(data, k = 20)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID39478111/test_MIBI.pdf", width = 12, height = 15)
    maps=list(B=c(45L,110L),Plasma=c(39L,124L,125L),CD4T=c(5L),CD8T=c(21L),Endothelial=c(100L,65L),Macrophage=c(14L,7L),Tumor=c(74L,129L,107L,11L,59L))

    
    
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
    write.csv(annotation, file = "../RawData/PMID39478111/celltype_annotation_MIBI.txt",quote = FALSE)
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID39478111/annotation_MIBI.pdf", width = 15, height = 7)
    
    svg("../RawData/PMID39478111/annotation_MIBI.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID39478111/bubble_MIBI.svg",width = 20,height = 4)
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
    "CD4T": ["CD4"],
    "Treg": ["FOXP3"],
    "CD8T": ["CD8A"],
    "Macrophage": ["CD163"],
    "Dendritic": ["ITGAX"],
    "Neutrophils": ["MPO"],
    "NK": ["NCAM1"],
    "Mast": ["KIT"],
    "Endothelial": ["PECAM1"],
    "Stromal": ["VIM"],
    "Tumor": ["CDH1"]
  }
}

import json
with open("../RawData/PMID39478111/astir_marker_MIBI.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID39478111/protein_avg_MIBI.csv --astir_marker ../RawData/PMID39478111/astir_marker_MIBI.json --output ./  > ../RawData/PMID39478111/log_MIBI.txt 2>&1

mv celltype_annotation_astir.txt ../RawData/PMID39478111/celltype_annotation_MIBI_astir.txt

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID39478111/protein_avg_MIBI.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD34","CD68","CD163","ITGAX","MPO","NCAM1","KIT","VIM","CDH1")]    
    annotation = read.csv("../RawData/PMID39478111/celltype_annotation_MIBI_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    #annotation$celltype[annotation$celltype == "Treg"] <- "CD4T"
    
    svg("../RawData/PMID39478111/annotation_MIBI_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID39478111/bubble_MIBI_astir.svg",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

cp ../RawData/PMID39478111/celltype_annotation.txt ../RawData/PMID39478111/celltype_annotation_final.txt 
cp  ../RawData/PMID39478111/annotation.svg ../RawData/PMID39478111/annotation_final.svg
cp  ../RawData/PMID39478111/bubble.svg ../RawData/PMID39478111/bubble_final.svg


#4. Divide samples
python spSampleDivision.py --input_dir ../RawData/PMID39478111/ --output_prefix ../Output/PMID39478111