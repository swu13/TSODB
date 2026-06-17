#!/bin/bash
set -euo pipefail


#1, Data Download
#https://zenodo.org/records/18434741

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data = pd.read_csv('../RawData/PMID38496566/Functional_data_FINAL.csv')
protein_avg =  data.iloc[:,np.r_[6:9,10:17,20:30]]
protein_avg.columns = protein_avg.columns.str.split('_').str[0]  #0-1
protein_avg.columns = np.where(protein_avg.columns == 'HLA.DR',protein_avg.columns, protein_avg.columns.str.split('.').str[0])

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
    
protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]   #latent membrane protein 1 (LMP1), ENV positive marker
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg.to_csv("../RawData/PMID38496566/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)

protein_avg.to_csv("../RawData/PMID38496566/protein_avg.csv")
spatial = data[['X_cent','Y_cent']]
spatial.columns=['X','Y']
spatial.to_csv("../RawData/PMID38496566/spatial.csv")
annotation = data[['Annotation']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID38496566/annotation.csv")
samples = data[['CoreID']]
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID38496566/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID38496566/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "HEME\tLymphNode\tCODEX\tHuman\tTSOP%04d\tPMID38496566\t%s\n", FNR+num-1, $0}' >> SPsamples.txt
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
    xdata = read.csv("../RawData/PMID38496566/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("PAX5","CD3E","FOXP3","CD8A","CD68","HLA-DRA","ITGAX","FUT4","GZMB","TNFRSF8")]         
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID38496566/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(27L,11L),CD4T=c(44L),Treg=c(35L),CD8T=c(17L),Macrophage=c(43L),Dendritic=c(42L),Neutrophils=c(32L),NK=c(33L),Tumor=c(46L))
    
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
    write.csv(annotation, file = "../RawData/PMID38496566/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID38496566/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID38496566/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID38496566/bubble.svg",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "B": ["PAX5"],
    "T": ["CD3E"],
    "Macrophage": ["CD68"],
    "Dendritic": ["ITGAX"],
    "Neutrophils": ["FUT4"],
    "NK": ["GZMB"],
    "Tumor": ["TNFRSF8"]
  }
}

import json
with open("../RawData/PMID38496566/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID38496566/protein_avg_nozscore.csv --astir_marker ../RawData/PMID38496566/astir_marker.json --output ../RawData/PMID38496566/ --anno_file ../RawData/PMID38496566/annotation.csv > ../RawData/PMID38496566/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID38496566/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("PAX5","CD3E","FOXP3","CD8A","CD68","HLA-DRA","ITGAX","FUT4","GZMB","TNFRSF8")]  
    annotation = read.csv("../RawData/PMID38496566/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL

    Neucells = rownames(subset(annotation,celltype == "Neutrophils"))
    annotation[Neucells, "celltype"] = ifelse(xdata[Neucells,"FUT4"]>0,"Neutrophils","unclassified")
    
    unclassified <- rownames(subset(annotation,celltype == "unclassified")) 
    markers <- c(
      "B" = "PAX5",
      "T" = "CD3E",
      "Macrophage" = "CD68",
      "Dendritic" = "ITGAX",
      "Neutrophils" = "FUT4",
      "NK" = "GZMB",
      "Tumor" = "TNFRSF8"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD8A"]>0,"CD8T","CD4T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID38496566/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID38496566/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()

    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID38496566/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID38496566/bubble_astir.pdf ../RawData/PMID38496566/bubble_astir.svg
cp  ../RawData/PMID38496566/annotation_astir.svg ../RawData/PMID38496566/annotation_final.svg
cp  ../RawData/PMID38496566/bubble_astir.svg ../RawData/PMID38496566/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID38496566/samples.csv  > ../RawData/PMID38496566/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID38496566/samples.csv  >> ../RawData/PMID38496566/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID38496566/ --output_dir ../SPOutput/