#!/bin/bash
set -euo pipefail


#1, Data Download
#https://zenodo.org/records/7961844?preview_file=20210205_CompSlide_NSCLC_TMA_LC.zip
unzip ../RawData/PMID38242124/20210205_CompSlide_NSCLC_TMA_LC.zip -d ../RawData/PMID38242124

#2. Generate protein expression, spatial, samples and annotation files 
python <<EOF
import pandas as pd
import numpy as np
from scipy import stats
import os
import glob
import json
input_dir = "../RawData/PMID38242124/20210205_CompSlide_NSCLC_TMA_LC"
files = glob.glob(os.path.join(input_dir, "*.txt"))
samples = pd.Series(files).apply(os.path.basename).str.extract(r'_LC_(.*?)\.txt$')[0]
#columns idential
ref = pd.read_csv(files[0], sep="\t")
for f in files:
    df = pd.read_csv(f, sep="\t")
    if not ref.columns.equals(df.columns):
        print(f)

with open("proteinmap.json", "r", encoding="utf-8") as f:
    protein_map = json.load(f)
alias_to_canonical = {}
for canonical, aliases in protein_map.items():
    alias_to_canonical[canonical] = canonical
    for alias in aliases:
        alias_to_canonical[alias] = canonical
#merge
protein_tables = []
xy_table = []
sample_table = []
for i, f in enumerate(files):
    df = pd.read_csv(f, sep="\t")
    protein_cols = [c for c in df.columns if "_" in c]
    protein_avg = df[protein_cols].copy()
    protein_avg = protein_avg.iloc[:,3:]
    protein_avg = protein_avg.drop(columns=["Iridium_1033((1253))Ir191(Ir191Di)"])
    protein_avg = protein_avg.drop(columns=["Iridium_1033((1254))Ir193(Ir193Di)"])
    protein_avg = protein_avg.drop(columns=["anti-Hu_1950((2832))Yb173(Yb173Di)"])
    protein_avg.columns = [c.split("_")[0] for c in protein_avg.columns]
    protein_avg.columns = protein_avg.columns.str.replace(r"\(.*", "",regex=True)
    protein_avg.columns=[alias_to_canonical.get(c, c) for c in protein_avg.columns] 
    protein_avg = protein_avg.rename(columns={"Histone": "H3","Cadheri":"CDH11","Collage":"COL1"}) 
    cells = [f"{samples[i]}_{idx}" for idx in protein_avg.index]
    protein_avg.index=cells
    protein_tables.append(protein_avg)
    spatial = df[['X','Y']]
    spatial.index=cells
    xy_table.append(spatial)
    sample_df = pd.DataFrame({'sample_id': samples[i]}, index=cells)
    sample_table.append(sample_df)
merged_protein_tables = pd.concat(protein_tables, axis=0)
protein_avg = np.arcsinh(merged_protein_tables)
protein_avg.to_csv("../RawData/PMID38242124/protein_avg_nozscore.csv")
protein_avg_zscore = stats.zscore(protein_avg)
protein_avg = pd.DataFrame(protein_avg_zscore, columns=protein_avg.columns, index=protein_avg.index)
protein_avg.to_csv("../RawData/PMID38242124/protein_avg.csv")
merged_xy_table = pd.concat(xy_table, axis=0)
merged_xy_table.to_csv("../RawData/PMID38242124/spatial.csv")
merged_sample_table = pd.concat(sample_table, axis=0)
merged_sample_table.to_csv("../RawData/PMID38242124/samples.csv")
EOF


#samples
num=`wc -l SPsamples.txt | cut -d" " -f1`
cut -d"," -f2 ../RawData/PMID38242124/samples.csv | tail -n +2 | sort | uniq  | awk -F"\t" -v num="$num" '{printf "LUCA\tLung\tIMC\tHuman\tTSOP%04d\tPMID38242124\t%s\n", FNR+num-1, $0}' >> SPsamples.txt

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
    xdata = read.csv("../RawData/PMID38242124/protein_avg_nonzero.csv",row.names=1,check.names=FALSE)        
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD34","CD68","HLA-DRA","MPO","FUT4","ACTA2","FAP","PanCytoK")]       
    pheno.out <- Rphenograph(data, k = 100)
    selecteddata <- data[as.double(pheno.out$names),]
    clusters <- as.numeric(membership(pheno.out))
    names(clusters) <- rownames(selecteddata)
    plot_clustering_heatmap(expr = selecteddata, cell_clustering = clusters,fontsize = 10,filename="../RawData/PMID38242124/test.pdf", width = 12, height = 50)
    
    maps=list(B=c(195L,196L,197L),T=c(187L,188L,194L,204L,205L),Endothelial=c(182L,183L),Macrophage=c(191L,192L),Neutrophils=c(208L,145L),Stromal=c(218L,219L,220L,221L,175L,176L),Tumor=c(1L,2L,3L,4L,5L,6L,7L,8L,9L))
    
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
    write.csv(annotation, file = "../RawData/PMID38242124/celltype_annotation.txt",quote = FALSE)
    
    labels = read.csv("../RawData/PMID38242124/annotation.csv",row.names=1)
    print(table(labels$celltype, annotation$celltype))
    
    svg("../RawData/PMID38242124/annotation.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    svg("../RawData/PMID38242124/bubble.svg",width = 8,height = 4)
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
    "Neutrophils": ["MPO"],
    "Stromal": ["ACTA2"],
    "Endothelial": ["PECAM1"],
    "Tumor": ["PanCytoK"]
  }
}
import json
with open("../RawData/PMID38242124/astir_marker.json", "w", encoding="utf-8") as f:
    json.dump(astir_marker_dict, f, indent=2, ensure_ascii=False)
EOF

python spcellannotation_astir.py --input ../RawData/PMID38242124/protein_avg_nozscore.csv --astir_marker ../RawData/PMID38242124/astir_marker.json --output ../RawData/PMID38242124/  > ../RawData/PMID38242124/log.txt 2>&1

Rscript -e '
    library(ggplot2)
    source("functions_clustering.R")
    xdata = read.csv("../RawData/PMID38242124/protein_avg.csv",row.names=1,check.names=FALSE)    
    data <- xdata[,c("MS4A1","CD3E","CD4","FOXP3","CD8A","PECAM1","CD68","MPO","FUT4","ACTA2","PanCytoK")]     
     
    annotation = read.csv("../RawData/PMID38242124/celltype_annotation_astir.txt",row.names=1)
    annotation$astir_conf <- NULL
    
    unclassified = rownames(subset(annotation,celltype %in% c("B"))) 
    markers <- c(
      "B" = "MS4A1",
      "T" = "CD3E",
      "Macrophage" = "CD68"
    )
    expr <- xdata[unclassified, markers]    
    annotation[unclassified, "celltype"] <-names(markers)[max.col(expr, ties.method = "first")]
    
    Tcells = rownames(subset(annotation,celltype == "T"))
    annotation[Tcells, "celltype"] = ifelse(xdata[Tcells,"CD4"]>xdata[Tcells,"CD8A"],"CD4T","CD8T")
    
    CD4Tcells = rownames(subset(annotation,celltype == "CD4T"))
    annotation[CD4Tcells, "celltype"] = ifelse(xdata[CD4Tcells,"FOXP3"]>0,"Treg","CD4T") 
    
    write.csv(annotation,"../RawData/PMID38242124/celltype_annotation_final.txt",quote=FALSE)
    
    svg("../RawData/PMID38242124/annotation_astir.svg")
    print(plot_clustering_heatmap(expr = data, cell_clustering = setNames(annotation$celltype, rownames(annotation)),fontsize = 10))
    dev.off()
    
    library(dplyr)
    library(tidyr)
    df <- cbind(annotation, xdata)
    long_df <- df %>% pivot_longer(cols = -celltype,names_to = "marker",values_to = "expression")
    bubble_df <- long_df %>% group_by(celltype, marker) %>% summarise(avg_exp = mean(expression),pct_exp = mean(expression > 0) * 100,.groups = "drop")
    pdf("../RawData/PMID38242124/bubble_astir.pdf",width = 12,height = 4)
    print(ggplot(bubble_df,aes(x = marker,y = celltype,size = pct_exp,color = avg_exp)) + geom_point() + scale_size(range = c(1, 8)) +
        scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45,hjust = 1), panel.grid = element_blank()))
    dev.off()
'

pdf2svg ../RawData/PMID38242124/bubble_astir.pdf ../RawData/PMID38242124/bubble_astir.svg
cp  ../RawData/PMID38242124/annotation_astir.svg ../RawData/PMID38242124/annotation_final.svg
cp  ../RawData/PMID38242124/bubble_astir.svg ../RawData/PMID38242124/bubble_final.svg

#4. Divide samples
head -n 1 ../RawData/PMID38242124/samples.csv  > ../RawData/PMID38242124/samples.TSO.csv
awk -F"\t" 'ARGIND==1{ID[$NF]=$5}ARGIND==2 && (FNR>1){split($0,a,",");print a[1]","ID[a[2]]}' SPsamples.txt ../RawData/PMID38242124/samples.csv  >> ../RawData/PMID38242124/samples.TSO.csv
python spSampleDivision.py --input_dir ../RawData/PMID38242124/ --output_dir ../SPOutput/