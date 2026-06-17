#!/bin/bash
set -euo pipefail


#1, Data Download
#https://figshare.com/s/4610a15363c8306dfa36, https://figshare.com/s/2005255a8b65de23109f and https://figshare.com/s/1d8c7ed76d4b3222ada4.

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import json
adata1 = sc.read_h5ad('../RawData/PMID37344596/TMA_1_processed.h5ad')
adata2 = sc.read_h5ad('../RawData/PMID37344596/TMA_2_processed.h5ad')
adata3 = sc.read_h5ad('../RawData/PMID37344596/TMA_3_processed.h5ad')
adata = sc.concat([adata1, adata2, adata3],axis=0,join="outer",label="batch",keys=["TMA1", "TMA2", "TMA3"])
mean_data = adata[:, adata.var_names.str.endswith('_mean')].X.toarray()
mean_vars = adata.var_names[adata.var_names.str.endswith('_mean')]
mean_vars_list = mean_vars.to_list()
df = pd.DataFrame(mean_data, columns=mean_vars, index=adata.obs_names)
parts = pd.Series(mean_vars_list).str.split('_', expand=True)
df.columns = pd.MultiIndex.from_arrays([parts[0], parts[1]], names=['protein', 'compartment'])
protein_avg = df.T.groupby(level='protein').mean().T
protein_avg = protein_avg.filter(regex='^(?!DAPI)')
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
protein_avg.to_csv("../RawData/PMID37344596/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID37344596/protein_avg.csv")
spatial = pd.DataFrame(adata.obsm["coordinates"], columns=['X','Y'], index=adata.obs_names)
spatial.to_csv("../RawData/PMID37344596/spatial.csv")
adata.obs["sample_id"].to_csv("../RawData/PMID37344596/samples.csv")
annotation = pd.DataFrame(adata.obs["celltype"],columns=["celltype"],index=adata.obs_names)
annotation.to_csv("../RawData/PMID37344596/annotation.csv")
EOF


#samples
echo -e "CancerType\tSite\tSequencing\tSpecies\tSample\tOrginalDataset\tOrginalSample" > SPsamples.txt
cut -d"," -f2 ../RawData/PMID37344596/samples.csv | tail -n +2 | sort | uniq | awk -F"\t" '{printf "BLCA\tUrothelium\tCODEX\tHuman\tTSOP%04d\tPMID37344596\t%s\n", NR, $0}' >> SPsamples.txt  


#3. Cell type annotaion based on cluster
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
    xdata = read.csv("../RawData/PMID37344596/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","ITGAX","ACTA2","PanCytoK")] 
    set.seed(1234)
    pheno.out <- Rphenograph(data, k = 15)
    
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID37344596/test.pdf", width = 12, height = 15)
    
    
    maps=list(B=c(22L),CD4T=c(9L),Treg=c(8L),CD8T=c(10L,16L),Endothelial=c(18L),Macrophage=c(13L),Dendritic=c(25L),Stromal=c(1L),Tumor=c(3L,30L))
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
    write.csv(annotation, file = "../RawData/PMID37344596/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID37344596/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID37344596/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID37344596/bubble.svg",width = 8,height = 4)
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
    "Dendritic": ["ITGAX"],
    "Stromal": ["ACTA2"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["PanCytoK"]
  }
}
import json
with open("../RawData/PMID37344596/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID37344596/protein_avg_nozscore.csv --astir_marker ../RawData/PMID37344596/astir_marker.json --output ../RawData/PMID37344596/ --anno_file ../RawData/PMID37344596/annotation.csv > ../RawData/PMID37344596/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID37344596/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","HLA-DRA","ITGAX","ACTA2","PanCytoK")] 
    annotation = read.csv("../RawData/PMID37344596/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("Dendritic"))) 
    markers <- c(
      "B" = "MS4A1",
      "T" = "CD3E",
      "Macrophage" = "CD68",
      "Dendritic" = "ITGAX",
      "Stromal" = "ACTA2",
      "Endothelial" = "PECAM1",
      "Tumor" = "PanCytoK"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID37344596/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID37344596/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID37344596/bubble_astir.pdf",width = 10,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

pdf2svg ../RawData/PMID37344596/bubble_astir.pdf ../RawData/PMID37344596/bubble_astir.svg
cp  ../RawData/PMID37344596/annotation_astir.svg ../RawData/PMID37344596/annotation_final.svg
cp  ../RawData/PMID37344596/bubble_astir.svg ../RawData/PMID37344596/bubble_final.svg

#4. Divide samples
mkdir ../SPOutput/
head -n 1 ../RawData/PMID37344596/samples.csv  > ../RawData/PMID37344596/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID37344596/samples.csv  >> ../RawData/PMID37344596/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID37344596/ --output_dir ../SPOutput/