#!/bin/bash
set -euo pipefail


#1, Data Download
#https://humantumoratlas.org/
#HTAN Synapse and TableS1 for tumor location

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import json
data = pd.read_csv('../RawData/HTAN_HTAPP_NBL/neuroblastoma_all_clustering.csv')
protein_avg = data.iloc[:,8:47]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical

protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]
protein_avg = protein_avg.T.groupby(level=0).mean().T

#protein_avg = np.arcsinh(protein_avg)
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg.to_csv("../RawData/HTAN_HTAPP_NBL/protein_avg_nozscore.csv")
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/HTAN_HTAPP_NBL/protein_avg.csv")
spatial = data[['x_cent', 'y_cent']].rename(columns={'x_cent': 'X', 'y_cent': 'Y'})
spatial.to_csv("../RawData/HTAN_HTAPP_NBL/spatial.csv")
annolabels = data[['cell_type']].rename(columns={'cell_type': 'celltype'})
annolabels.to_csv("../RawData/HTAN_HTAPP_NBL/annotation.csv")
data["sample_id"] = data["HTAPP_ID"].astype(str) + "_" + data["FOV"].astype(str)
data["sample_id"].to_csv("../RawData/HTAN_HTAPP_NBL/samples.csv")
EOF


#samples
cut -d"," -f2 ../RawData/HTAN_HTAPP_NBL/samples.csv | tail -n +2 | sort | uniq > ../RawData/HTAN_HTAPP_NBL/samples.txt
num=`wc -l SPsamples.txt | cut -d" " -f1`
awk -F"\t" -v num="$num" 'ARGIND==1{split($5,m,"_");id=m[2]"_"m[5];pro=m[1]"_"m[2]"_"m[3];samples[id]=$1"\t"$2"\t"$3"\t"$4;project[id]=pro}ARGIND==2{split($1,m,"_");ids=m[1]"_"m[2];printf samples[ids]"\tTSOP%04d\t"project[ids]"\t%s\n", FNR+num-1, $0}' ../RawData/HTAN_HTAPP_NBL/HTAN_HTAPP_NBL.txt ../RawData/HTAN_HTAPP_NBL/samples.txt >> SPsamples.txt
rm -f ../RawData/HTAN_HTAPP_NBL/samples.txt  


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
    xdata = read.csv("../RawData/HTAN_HTAPP_NBL/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("PTPRC","MS4A1","CD38","SDC1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD34","CD68","CD163","ITGAX","KIT","MPO","NCAM1","CDH1","VIM","STAT5A")]   


    
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/HTAN_HTAPP_NBL/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(43L),CD4T=c(8L,14L,4L),CD8T=c(7L,10L),Endothelial=c(27L,22L,6L,23L),Macrophage=c(33L,34L),Mast=c(32L),Stromal=c(44L,31L,28L,16L),Tumor=c(19L,35L))
    
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
    write.csv(annotation, file = "../RawData/HTAN_HTAPP_NBL/celltype_annotation.txt",quote = FALSE)
    
    svg("../RawData/HTAN_HTAPP_NBL/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/HTAN_HTAPP_NBL/bubble.svg",width = 8,height = 4)
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
    "Macrophage": ["CD163"],
    "Dendritic": ["ITGAX"],
    "Neutrophils": ["MPO"],
    "NK": ["NCAM1"],
    "Stromal": ["VIM"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["CDH1"]
  }
}

import json
with open("../RawData/HTAN_HTAPP_NBL/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/HTAN_HTAPP_NBL/protein_avg_nozscore.csv --astir_marker ../RawData/HTAN_HTAPP_NBL/astir_marker.json --output ../RawData/HTAN_HTAPP_NBL/ --anno_file ../RawData/HTAN_HTAPP_NBL/annotation.csv > ../RawData/HTAN_HTAPP_NBL/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/HTAN_HTAPP_NBL/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD163","ITGAX","KIT","MPO","NCAM1","VIM","CDH1")]   
    annotation = read.csv("../RawData/HTAN_HTAPP_NBL/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("Macrophage","Endothelial"))) 
    markers <- c(
      "B" = "MS4A1",
      "Plasma" = "SDC1",
      "T" = "CD3E",
      "Macrophage" = "CD163",
      "Dendritic" = "ITGAX",
      "Neutrophils" = "MPO",
      "NK" = "NCAM1",
      "Stromal" = "VIM",
      "Endothelial" = "PECAM1",
      "Tumor" = "CDH1"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/HTAN_HTAPP_NBL/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/HTAN_HTAPP_NBL/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/HTAN_HTAPP_NBL/bubble_astir.pdf",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/HTAN_HTAPP_NBL/bubble_astir.pdf ../RawData/HTAN_HTAPP_NBL/bubble_astir.svg
cp  ../RawData/HTAN_HTAPP_NBL/annotation_astir.svg ../RawData/HTAN_HTAPP_NBL/annotation_final.svg
cp  ../RawData/HTAN_HTAPP_NBL/bubble_astir.svg ../RawData/HTAN_HTAPP_NBL/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/HTAN_HTAPP_NBL/samples.csv  > ../RawData/HTAN_HTAPP_NBL/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/HTAN_HTAPP_NBL/samples.csv  >> ../RawData/HTAN_HTAPP_NBL/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/HTAN_HTAPP_NBL/ --output_dir ../SPOutput/