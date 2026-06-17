#!/bin/bash
set -euo pipefail


#1, Data Download
#https://zenodo.org/records/10778429,brm_measurements.csv,pdg_measurements.csv

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import scanpy as sc
import json
data1 = pd.read_csv('../RawData/PMID38622132/brm_measurements.csv')
data1['sample_id'] = data1['Image'].str.split('.tif').str[0]
data1.index = (data1["sample_id"].astype(str) + "_" + data1.index.astype(str))
data2 = pd.read_csv('../RawData/PMID38622132/pdg_measurements.csv',encoding='latin1')
data2['sample_id'] = data2['Image'].str.split('.tif').str[0]
data2.index = (data2["sample_id"].astype(str) + "_" + data2.index.astype(str))
protein_avg1 =  data1.iloc[:,8:38]
protein_avg1.columns = protein_avg1.columns.str.split(':').str[0]
protein_avg2 =  data2.iloc[:,9:53]
protein_avg2.columns = protein_avg2.columns.str.split('\.\.').str[0]
protein_avg2 = protein_avg2.drop(columns=["Phalloidin","WGA","RedDot2","GFP","aTubulin"])
protein_avg2.columns = protein_avg2.columns.str.replace("NF.H", "NF-H")
protein_avg = pd.concat([protein_avg1, protein_avg2], axis=0).fillna(0)

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical

protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns]  
protein_avg = protein_avg.T.groupby(level=0).mean().T
protein_avg.to_csv("../RawData/PMID38622132/protein_avg_nozscore.csv")
protein_avg = np.arcsinh(protein_avg)
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID38622132/protein_avg.csv")
spatial1 = data1[['Centroid X ¦Ìm','Centroid Y ¦Ìm']]
spatial1.columns=['X','Y']
spatial2 = data2[['Centroid.X','Centroid.Y']]
spatial2.columns=['X','Y']
spatial = pd.concat([spatial1, spatial2], axis=0, ignore_index=True)
spatial.to_csv("../RawData/PMID38622132/spatial.csv")
samples = pd.concat([data1[['sample_id']],data2[['sample_id']]],axis=0, ignore_index=True)
samples.to_csv("../RawData/PMID38622132/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID38622132/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{if($0~".ome"){printf "GBM\tBrain\tHIFI\tMouse\tTSOP%04d\tPMID38622132\t%s\n", FNR+num-1, $0}else{printf "BRCA\tBrain\tHIFI\tMouse\tTSOP%04d\tPMID38622132\t%s\n", FNR+num-1, $0}}' >> SPsamples.txt
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
    xdata = read.csv("../RawData/PMID38622132/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("PTPRC","CD3E","CD8A","PECAM1","CD68","AIF1","P2RY12","Ly6b","S100A8","MPO","ACTA2","SOX9","SOX2","OLIG2")]     
    set.seed(1234)   

    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID38622132/test.pdf", width = 12, height = 15)
    
    maps=list(B=c(14L),CD4T=c(3L),Treg=c(5L),CD8T=c(6L),Endothelial=c(1L),Macrophage=c(8L),Monocyte=c(19L),Stromal=c(28L),Tumor=c(31L,13L,23L,38L))
    
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
    write.csv(annotation, file = "../RawData/PMID38622132/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID38622132/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID38622132/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID38622132/bubble.svg",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

#3. Cell type annotaion based on astir
python <<EOF
astir_marker_dict={
  "cell_type": {
    "T": ["CD3E"],
    "Macrophage": ["CD68"],
    "Neutrophils": ["MPO"],
    "Microglia": ["P2RY12"],
    "Stromal": ["ACTA2"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["SOX2"]
  }
}


import json
with open("../RawData/PMID38622132/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID38622132/protein_avg_nozscore.csv --astir_marker ../RawData/PMID38622132/astir_marker.json --output ../RawData/PMID38622132/ --anno_file ../RawData/PMID38622132/annotation.csv > ../RawData/PMID38622132/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID38622132/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("PTPRC","CD3E","CD8A","PECAM1","CD68","AIF1","P2RY12","Ly6b","S100A8","MPO","ACTA2","SOX9","SOX2","OLIG2")]        
    annotation = read.csv("../RawData/PMID38622132/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    write.csv(annotation,"../RawData/PMID38622132/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID38622132/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID38622132/bubble_astir.pdf",width = 8,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'
pdf2svg ../RawData/PMID38622132/bubble_astir.pdf ../RawData/PMID38622132/bubble_astir.svg
cp  ../RawData/PMID38622132/annotation_astir.svg ../RawData/PMID38622132/annotation_final.svg
cp  ../RawData/PMID38622132/bubble_astir.svg ../RawData/PMID38622132/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID38622132/samples.csv  > ../RawData/PMID38622132/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID38622132/samples.csv  >> ../RawData/PMID38622132/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID38622132/ --output_dir ../SPOutput/