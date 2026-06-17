#!/bin/bash
set -euo pipefail


#1, Data Download
#https://zenodo.org/records/13901180

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
adata = sc.read_h5ad('../RawData/PMID40189697/20240215_COAD_cancer_control_umapped_annotated.h5ad')
adata_cancer = adata[adata.obs["batch"] == "Cancer"].copy()
data = pd.DataFrame(adata_cancer.X, columns=adata_cancer.var_names, index=adata_cancer.obs_names)
with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
data.columns=[alias_to_canonical.get(c, c) for c in data.columns]  
data = data.rename(columns={"keratin": "PanCytoK","Collagen": "COL1"})

data.to_csv("../RawData/PMID40189697/protein_avg.csv")
data_raw = data - data.min(axis=0)
data_raw.to_csv("../RawData/PMID40189697/protein_avg_nonnegative.csv")

spatial = adata_cancer.obs[['AreaShape_Center_X','AreaShape_Center_Y']]
spatial.columns=['X','Y']
spatial.to_csv("../RawData/PMID40189697/spatial.csv")
annotation = adata_cancer.obs[['cell_type']]
annotation.columns = ['celltype']
annotation.to_csv("../RawData/PMID40189697/annotation.csv")
samples = adata_cancer.obs[['alt_identifier']]
samples.columns = ['sample_id']
samples.to_csv("../RawData/PMID40189697/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID40189697/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "CRCA\tColorectum\tIMC\tHuman\tTSOP%04d\tPMID40189697\t%s\n", FNR+num-1, $0}' >> SPsamples.txt 

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
    xdata = read.csv("../RawData/PMID40189697/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","CD68","GZMB","ACTA2","PanCytoK")]     
    set.seed(1234)    
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID40189697/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(9L),CD4T=c(11L),Treg=c(18L),CD8T=c(6L),Macrophage=c(10L),NK=c(40L),Stromal=c(1L,21L),Tumor=c(14L,32L))
    
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
    write.csv(annotation, file = "../RawData/PMID40189697/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID40189697/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID40189697/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID40189697/bubble.svg",width = 8,height = 4)
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
    "NK": ["GZMB"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK"]
  }
}
import json
with open("../RawData/PMID40189697/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF


python spcellannotation_astir.py --input ../RawData/PMID40189697/protein_avg_nonnegative.csv --astir_marker ../RawData/PMID40189697/astir_marker.json --output ../RawData/PMID40189697/ --anno_file ../RawData/PMID40189697/annotation.csv > ../RawData/PMID40189697/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID40189697/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","CD68","GZMB","ACTA2","PanCytoK")]  
    annotation = read.csv("../RawData/PMID40189697/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    Maccells = rownames(subset(annotation,celltype == "Macrophage"))
    annotation[Maccells, "celltype"] = ifelse(xdata[Maccells,"CD68"]>0,"Macrophage","unclassified")
    
    unclassified <- rownames(subset(annotation,celltype == "unclassified")) 
    markers <- c(
      B = "MS4A1",
      T = "CD3E",
      Macrophage = "CD68",
      NK = "GZMB",
      Stromal = "ACTA2",
      Tumor = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID40189697/celltype_annotation_final.txt",quote=FALSE)

    svg("../RawData/PMID40189697/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID40189697/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID40189697/bubble_astir.pdf ../RawData/PMID40189697/bubble_astir.svg
cp  ../RawData/PMID40189697/annotation_astir.svg ../RawData/PMID40189697/annotation_final.svg
cp  ../RawData/PMID40189697/bubble_astir.svg ../RawData/PMID40189697/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID40189697/samples.csv  > ../RawData/PMID40189697/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID40189697/samples.csv  >> ../RawData/PMID40189697/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID40189697/ --output_dir ../SPOutput/