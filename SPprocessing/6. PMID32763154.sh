#!/bin/bash
set -euo pipefail


#1, Data Download
#https://data.mendeley.com/datasets/mpjzbtfgfr/1; Supplementary File for sample information (PMID32763154.txt, only kept cancer)

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import json
data=pd.read_csv("../RawData/PMID32763154/CRC_clusters_neighborhoods_markers.csv",index_col=0)
clinical = pd.read_csv("../RawData/PMID32763154/PMID32763154.txt",sep="\t",index_col=0)
tumorsamples = clinical["Sample"].str.replace(r"^PMID\d+_", "", regex=True)
data_tumor = data[data["Region"].isin(set(tumorsamples))]
protein_avg=data_tumor.iloc[:, np.r_[11:60, 69:76]]
protein_avg.columns = [col.split(":")[0].split(" - ")[0] for col in protein_avg.columns]

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
protein_avg.columns = [alias_to_canonical.get(c, c) for c in protein_avg.columns]
protein_avg.columns = protein_avg.columns.str.replace("Cytokeratin", "PanCytoK")
protein_avg = protein_avg.T.groupby(level=0).mean().T

protein_avg = np.arcsinh(protein_avg)
protein_avg.to_csv("../RawData/PMID32763154/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID32763154/protein_avg.csv")
spatial = data_tumor.set_index("CellID")[['X:X', 'Y:Y']].rename(columns={'X:X': 'X', 'Y:Y': 'Y'})
spatial.to_csv("../RawData/PMID32763154/spatial.csv")
samples = data_tumor.set_index("CellID")[['Region']].rename(columns={'Region': 'sample_id'})
samples.to_csv("../RawData/PMID32763154/samples.csv")
annolabels = data_tumor.set_index("CellID")[['ClusterName']].rename(columns={'ClusterName': 'celltype'})
annolabels.to_csv("../RawData/PMID32763154/annotation.csv")
EOF 


#samples
cut -d"," -f2 ../RawData/PMID32763154/samples.csv | tail -n +2 | sort | uniq  >  ../RawData/PMID32763154/samples.txt
num=`wc -l SPsamples.txt | cut -d" " -f1`
awk -F"\t" -v num="$num" 'ARGIND==1{split($5,a,"_");sample[a[2]]=$1"\t"$2"\t"$3"\t"$4;project[a[2]]=a[1]"\t"a[2]}ARGIND==2{printf sample[$1]"\tTSOP%04d\t"project[$1]"\n", FNR+num-1}' ../RawData/PMID32763154/PMID32763154.txt ../RawData/PMID32763154/samples.txt  >> SPsamples.txt
rm -f ../RawData/PMID32763154/samples.txt

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
    xdata = read.csv("../RawData/PMID32763154/protein_avg.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","SDC1","CD38","CD3E","CD4","FOXP3","CD8A","PECAM1","CD34","CD68","CD163","HLA-DRA","ITGAX","FUT4","B3GAT1","NCAM1","ACTA2","PDPN","PanCytoK")]   
    
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID32763154/test.pdf", width = 12, height = 15)
    maps=list(B=c(40L,28L,19L,52L,13L),Plasma=c(54L,34L),CD4T=c(20L),CD8T=c(24L),Endothelial=c(51L,22L),Macrophage=c(17L),Neutrophils=c(11L,6L),NK=c(45L),Stromal=c(29L,43L,30L,31L),Tumor=c(49L,32L,55L))
    
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
    write.csv(annotation, file = "../RawData/PMID32763154/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID32763154/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID32763154/annotation.pdf", width = 15, height = 7)
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID32763154/bubble.svg",width = 20,height = 4)
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
    "Neutrophils": ["FUT4"],
    "NK": ["B3GAT1"],
    "Endothelial": ["PECAM1"],
    "Stromal": ["ACTA2"],
    "Tumor": ["PanCytoK"]
  }
}

import json
with open("../RawData/PMID32763154/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID32763154/protein_avg_nozscore.csv --astir_marker ../RawData/PMID32763154/astir_marker.json --output ../RawData/PMID32763154/ --anno_file ../RawData/PMID32763154/annotation.csv > ../RawData/PMID32763154/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID32763154/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","SDC1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","CD163","ITGAX","FUT4","B3GAT1","ACTA2","PanCytoK")]   
    annotation = read.csv("../RawData/PMID32763154/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf=NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("Plasma","Tumor"))) 
    markers <- c(
      "B" = "MS4A1",
      "Plasma" = "SDC1",
      "T" = "CD3E",
      "Macrophage" = "CD163",
      "Dendritic" = "ITGAX",
      "Neutrophils" = "FUT4",
      "NK" = "B3GAT1",
      "Endothelial" = "PECAM1",
      "Stromal" = "ACTA2",
      "Tumor" = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID32763154/celltype_annotation_final.txt",quote=FALSE)
    
    plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10,filename="../RawData/PMID32763154/annotation_astir.pdf")
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID32763154/bubble_astir.pdf",width = 13,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

pdf2svg ../RawData/PMID32763154/annotation_astir.pdf ../RawData/PMID32763154/annotation_astir.svg
pdf2svg ../RawData/PMID32763154/bubble_astir.pdf ../RawData/PMID32763154/bubble_astir.svg
cp  ../RawData/PMID32763154/annotation_astir.svg ../RawData/PMID32763154/annotation_final.svg
cp  ../RawData/PMID32763154/bubble_astir.svg ../RawData/PMID32763154/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID32763154/samples.csv  > ../RawData/PMID32763154/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID32763154/samples.csv  >> ../RawData/PMID32763154/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID32763154/ --output_dir ../SPOutput/