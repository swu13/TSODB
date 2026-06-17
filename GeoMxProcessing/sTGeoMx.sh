#!/bin/bash
set -euo pipefail

############################################
# 1.GeoMx
############################################

#GSE200563
mkdir -p ../RawData/GSE200563
wget -O ../RawData/GSE200563/Expression.txt "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200563&format=file&file=GSE200563%5Fprocessed%5Fdata%2Etxt%2Egz"
head -n 1 ../RawData/GSE200563/Expression.txt | awk -F"\t" '{out="Gene"};for(i=2;i<=NF;i++){out=out"\t"$i};print out}' >  ../RawData/GSE200563/Expression1.txt
tail -n +2 ../RawData/GSE200563/Expression.txt >> ../RawData/GSE200563/Expression1.txt
mv ../RawData/GSE200563/Expression1.txt ../RawData/GSE200563/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE200563/clinical.txt
grep GSE200563 STsamples1.txt | awk -F"\t" '{print $10"\t"$2"\t"$6}' >> ../RawData/GSE200563/clinical.txt


#GSE254054
mkdir -p ../RawData/GSE254054
wget -O "../RawData/GSE254054/Hs_R_NGS_CTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE254054&format=file&file=GSE254054%5FHs%5FR%5FNGS%5FCTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE254054/GSE254054_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE254054&format=file"
tar -xvf ../RawData/GSE254054/GSE254054_RAW.tar -C ../RawData/GSE254054/
gunzip ../RawData/GSE254054/*gz
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE254054/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE254054/Hs_R_NGS_CTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE254054/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE254054/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE254054/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE254054/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'

echo -e "sample\tgroup" > ../RawData/GSE254054/clinical.txt
awk -F"\t" 'ARGIND==1{a[$1]=$2}\
            ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"a[$i]}}'\
            ../RawData/GSE254054/annotation.txt\
            ../RawData/GSE254054/Expression.txt |\
            awk -F"\t" '{if($2~"Non-epi"){print $1"\tStromal"}else{print $1"\tTumor"}}' >> ../RawData/GSE254054/clinical.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE254054/clinical.1.txt
awk -F"\t" 'ARGIND==1 && ($7=="GSE254054"){slide[$8]=$6}\
           ARGIND==2 && (FNR>1){split($1,a,"_");print $1"\t"$2"\t"slide[a[1]]}'\
           STsamples1.txt\
           ../RawData/GSE254054/clinical.txt >> ../RawData/GSE254054/clinical.1.txt
mv ../RawData/GSE254054/clinical.1.txt ../RawData/GSE254054/clinical.txt
