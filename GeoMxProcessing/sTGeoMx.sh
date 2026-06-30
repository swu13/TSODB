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


#GSE313748
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313748/soft/GSE313748_family.soft.gz
gunzip GSE313748_family.soft.gz
awk '
/^\^SAMPLE =/ {
    gsm=$3
}
/!Sample_characteristics_ch1 = biopsy site:/ {
    sub(/.*biopsy site: /,"")
    print gsm"\t"$0
}
' GSE313748_family.soft > GSM_biopsy_site.txt


#GSE290451,no expression
mkdir -p ../RawData/GSE290451
wget -O "../RawData/GSE290451/GSE290451_bi_bioqc_ntc1500_sa2500_v1.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE290451&format=file&file=GSE290451%5Fbi%5Fbioqc%5Fntc1500%5Fsa2500%5Fv1%2Exlsx"
pip install openpyxl

#GSE286410
mkdir -p ../RawData/GSE286410
wget -O "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE286410&format=file&file=GSE286410%5FGeoMx%5Fqnorm%2Etxt%2Egz"
gunzip "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt.gz"
sed -i 's/"//g' "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt"
sed -i 's/ /\t/g' "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt"
cp "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt" "../RawData/GSE286410/Expression.txt"
head -n 1 "../RawData/GSE286410/Expression.txt" | awk -F"\t" '{print "Gene\t"$0}' > "../RawData/GSE286410/Expression.1.txt"
tail -n +2 "../RawData/GSE286410/Expression.txt" >> "../RawData/GSE286410/Expression.1.txt"
mv "../RawData/GSE286410/Expression.1.txt" "../RawData/GSE286410/Expression.txt"

wget -O "../RawData/GSE286410/GSE286410_Annotation.txt.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE286410&format=file&file=GSE286410%5FAnnotation%2Etxt%2Egz"
gunzip "../RawData/GSE286410/GSE286410_Annotation.txt.gz"
sed -i 's/"//g' "../RawData/GSE286410/GSE286410_Annotation.txt"
echo -e "sample\tgroup\tslide" > ../RawData/GSE286410/clinical.txt
awk -F" " 'ARGIND==1 && (FNR>1){Sa[$2]=$13"\tSlide"$5}\
           ARGIND==2 && (FNR==1){split($0,a,",");for(i=1;i<=length(a);i++){print a[i]"\t"Sa[a[i]]}}'\
           "../RawData/GSE286410/GSE286410_Annotation.txt"\
           "../RawData/GSE286410/GSE286410_GeoMx_qnorm.txt" >> ../RawData/GSE286410/clinical.txt
#GSE287631, no clear group
#GSE288388
mkdir -p ../RawData/GSE288388
wget -O "../RawData/GSE288388/GSE288388_GeoMx_counts_from_GeoMx_DST_groups.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288388&format=file&file=GSE288388%5FGeoMx%5Fcounts%5Ffrom%5FGeoMx%5FDST%5Fgroups%2Exlsx"
wget -O "../RawData/GSE288388/GSE288388_GeoMX_experiment_ROI_information_and_primer_data.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288388&format=file&file=GSE288388%5FGeoMX%5Fexperiment%5FROI%5Finformation%5Fand%5Fprimer%5Fdata%2Exlsx"

python <<EOF
import pandas as pd
exp_raw = pd.read_excel("../RawData/GSE288388/GSE288388_GeoMx_counts_from_GeoMx_DST_groups.xlsx",engine="openpyxl")
metadata = pd.read_excel("../RawData/GSE288388/GSE288388_GeoMX_experiment_ROI_information_and_primer_data.xlsx",engine="openpyxl")
scan_name = exp_raw.iloc[3,3:]
roi_name  = exp_raw.iloc[4,3:]
treat     = exp_raw.iloc[5,3:]
scan_tag  = exp_raw.iloc[6,3:]
segment   = exp_raw.iloc[8,3:]
group = (treat + "_" + scan_tag + "_" + segment).astype(str).str.replace(",", "/")
sample_info = pd.DataFrame({"slide": scan_name.values,"roi": roi_name.values.astype(str),"group": group.values})
start_idx = exp_raw[exp_raw.iloc[:, 1] == "All Targets"].index[0]
exp = exp_raw.iloc[start_idx:].copy()
exp = exp.reset_index(drop=True)
gene_col = exp.columns[2] 
expr_mat = exp.iloc[:, 3:]
expr_mat.insert(0, "gene", exp[gene_col].values)
meta = metadata.copy()
meta = meta.rename(columns={meta.columns[0]: "sample", meta.columns[1]: "slide", meta.columns[3]: "Treatment", meta.columns[5]: "Scan Tag", meta.columns[7]: "roi", 
                            meta.columns[8]: "Segment Name", })
meta = meta[meta["sample"].astype(str).str.startswith("DSP-")]
meta["group"] = ( meta["Treatment"].astype(str) + "_" + meta["Scan Tag"].astype(str) + "_" + meta["Segment Name"].astype(str)).str.replace(",", "/")
meta = meta[["sample", "group", "slide","roi"]]
sample_info["key"] = (sample_info["slide"].astype(str) + "|" + sample_info["roi"].astype(str).str.zfill(3) + "|" + sample_info["group"].astype(str))
meta["key"] = (meta["slide"].astype(str) + "|" + meta["roi"].astype(str).str.zfill(3) + "|" + meta["group"].astype(str))
expr_mat.columns = ["gene"] + sample_info["key"].tolist()
key2sample = dict(zip(meta["key"], meta["sample"]))
expr_mat.columns = (["gene"] + [key2sample[k] for k in sample_info["key"]])
expr_mat.to_csv("../RawData/GSE288388/Expression.txt", sep=",", index=False)
meta = meta[["sample", "group", "slide"]]
sample_order = expr_mat.columns[1:]
meta = (meta.set_index("sample").loc[sample_order].reset_index())
meta.to_csv("../RawData/GSE288388/clinical.txt", sep="\t", index=False)
EOF

head -n 1 ../RawData/GSE288388/Expression.txt > ../RawData/GSE288388/Expression.1.txt
awk 'BEGIN{FS=OFS=","}\
   ARGIND==1{split($0,a,"\t");if((a[1]!="") && (a[2]!="")){Gene[a[1]]=a[2]}}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/GSE288388/Expression.txt >> ../RawData/GSE288388/Expression.1.txt
mv ../RawData/GSE288388/Expression.1.txt ../RawData/GSE288388/Expression.txt
sed -i 's/,/\t/g' ../RawData/GSE288388/Expression.txt
originalpath=`pwd`
rawsample="GSE288388"
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/GSE288388/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/GSE288388/Expression.1.txt
tail -n +2 ../RawData/GSE288388/Expression.txt >> ../RawData/GSE288388/Expression.1.txt
mv ../RawData/GSE288388/Expression.1.txt ../RawData/GSE288388/Expression.txt

wget -O "../RawData/GSE288388/GSE288388_family.soft.gz" https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288388/soft/GSE288388_family.soft.gz
gunzip "../RawData/GSE288388/GSE288388_family.soft.gz"
awk '
/^\^SAMPLE =/ {
    gsm=$3
}
/^!Sample_description = Library name:/ {
    sub(/^!Sample_description = Library name: /,"")
    print gsm"\t"$0
}' ../RawData/GSE288388/GSE288388_family.soft > ../RawData/GSE288388/GSM_biopsy_site.txt


#GSE288406
mkdir -p ../RawData/GSE288406
wget -O "../RawData/GSE288406/GSE288406_HNSCC_IPA_raw_counts.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288406&format=file&file=GSE288406%5FHNSCC%5FIPA%5Fraw%5Fcounts%2Ecsv%2Egz"
wget -O "../RawData/GSE288406/GSE288406_HNSCC_IPA_raw_metadata.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288406&format=file&file=GSE288406%5FHNSCC%5FIPA%5Fraw%5Fmetadata%2Ecsv%2Egz"
gunzip ../RawData/GSE288406/*
mv ../RawData/GSE288406/GSE288406_HNSCC_IPA_raw_counts.csv ../RawData/GSE288406/Expression.txt
dos2unix ../RawData/GSE288406/Expression.txt
sed -i 's/,/\t/g' ../RawData/GSE288406/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE288406/clinical.txt
awk -F"," 'ARGIND==1 && (FNR>1){split($1,a,"_");Sa[$1]=$4"_"$6"_"$10"_"$12"_"a[2]"\tHNSCC_"int(a[1])}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE288406/GSE288406_HNSCC_IPA_raw_metadata.csv\
           ../RawData/GSE288406/Expression.txt >> ../RawData/GSE288406/clinical.txt

           
#GSE289272
mkdir -p ../RawData/GSE289272
wget -O "../RawData/GSE289272/GSE289272_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE289272&format=file"
wget -O "../RawData/GSE289272/GSE289272_GeoMx_ROI_metadata.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE289272&format=file&file=GSE289272%5FGeoMx%5FROI%5Fmetadata%2Exlsx"
wget -O "../RawData/GSE289272/GSE289272_GeoMx_Hs_CTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE289272&format=file&file=GSE289272%5FGeoMx%5FHs%5FCTA%5Fv1%2E0%2Epkc%2Egz"
tar -xvf ../RawData/GSE289272/GSE289272_RAW.tar -C ../RawData/GSE289272/
gunzip ../RawData/GSE289272/*

wget -O "../RawData/GSE289272/GSE289272_family.soft.gz" https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289272/soft/GSE289272_family.soft.gz
gunzip "../RawData/GSE289272/GSE289272_family.soft.gz"
awk '
/^\^SAMPLE =/ {
    if (gsm != ""){print gsm "\t" title "\t" lib "\t" file}
    gsm=$3;title="";lib="";file="";next
}
/^!Sample_description = Library name:/ {sub(/^!Sample_description = Library name: /,"");lib=$0}
/^!Sample_title =/ {sub(/^!Sample_title = /,"");title=$0}
/^!Sample_supplementary_file_1 =/ {sub(/^!Sample_supplementary_file_1 = /,"");file=$0}
END {if (gsm != ""){print gsm "\t" title "\t" lib "\t" file}}
' ../RawData/GSE289272/GSE289272_family.soft |\
     awk -F"\t" '{split($4,a,"/");if($2~"nasopharyngeal carcinoma tissue"){group="NasopharyngealCarcinoma"}else{group="ChronicNasopharyngitis"};\
                  print $1"\t"group"\t"$3"\t"a[length(a)]}' > ../RawData/GSE289272/GSM_biopsy_site.txt
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE289272/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE289272/GSE289272_GeoMx_Hs_CTA_v1.0.pkc")       
    
    map <- read.table("../RawData/GSE289272/GSM_biopsy_site.txt",sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(map) <- c("GSM", "group", "ROI", "dcc")
    map$dcc <- sub("\\\\.gz$", "", basename(map$dcc))
    pheno <- data.frame(SampleID = map$dcc,GSM = map$GSM,ROI = map$ROI,stringsAsFactors = FALSE)

    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoData = pheno,phenoDataDccColName = "SampleID")
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE289272/Expression.txt",sep=",",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE289272/GSE289272_GeoMx_ROI_metadata.xlsx", sheet = 1)
    df2 <- merge(df, map[, c("ROI", "dcc", "group")], by.x = "SegmentDisplayName",by.y = "ROI",all.x = TRUE)
    metadata <- df2[,c("dcc","group","SlideName")]
    colnames(metadata) <- c("sample", "group", "slide")  
    
    write.table(metadata, "../RawData/GSE289272/clinical.txt", sep = ",", quote = FALSE, row.names = FALSE)
'
sed -i 's/ //g' "../RawData/GSE289272/Expression.txt"
sed -i 's/,/\t/g' "../RawData/GSE289272/Expression.txt"
sed -i 's/,/\t/g' "../RawData/GSE289272/clinical.txt"


#GSE284579 protein
#GSE284100 Visium withoug spatial location
#GSE281413
mkdir -p ../RawData/GSE281413
wget -O "../RawData/GSE281413/GSE281413_processed_data.txt.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE281413&format=file&file=GSE281413%5Fprocessed%5Fdata%2Etxt%2Egz"
gunzip ../RawData/GSE281413/*
cp "../RawData/GSE281413/GSE281413_processed_data.txt" "../RawData/GSE281413/Expression.txt"
dos2unix "../RawData/GSE281413/Expression.txt"
head -n 1 ../RawData/GSE281413/GSE281413_processed_data.txt | tr ',' '\n' | less
wget -O "../RawData/GSE281413/GSE281413_family.soft.gz" https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281413/soft/GSE281413_family.soft.gz
gunzip "../RawData/GSE281413/GSE281413_family.soft.gz"
awk '
/^\^SAMPLE =/ {
    if (gsm != ""){print gsm"\t"title"\t"cancer"\t"tissue"\t"treatment"\t"segment}
    gsm=$3;title="";tissue="";cancer="";treatment="";segment="";next
}
/^!Sample_title =/ {sub(/^!Sample_title = /,""); title=$0}
/^!Sample_characteristics_ch1 = tissue:/ {sub(/^!Sample_characteristics_ch1 = tissue: /,""); tissue=$0}
/^!Sample_characteristics_ch1 = treatment:/ {sub(/^!Sample_characteristics_ch1 = treatment: /,""); treatment=$0}
/^!Sample_characteristics_ch1 = cancer:/ {sub(/^!Sample_characteristics_ch1 = cancer: /,""); cancer=$0}
/^!Sample_characteristics_ch1 = segment:/ {sub(/^!Sample_characteristics_ch1 = segment: /,""); segment=$0}
END {if (gsm != ""){print gsm"\t"title"\t"cancer"\t"tissue"\t"treatment"\t"segment}}
' ../RawData/GSE281413/GSE281413_family.soft |\
     awk -F"\t" '{if($3=="No"){group="Normal_"$6}else{if(($4~"early onset") || ($4~"EO CRC")){group="EarlyOnset_"$6}else{group="AverageOnset_"$6}};print $1","$2","group","$5}' > ../RawData/GSE281413/GSM_biopsy_site.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE281413/clinical.txt
awk -F"," 'ARGIND==1{Sa[$2]=$3"\t"$4}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE281413/GSM_biopsy_site.txt\
           ../RawData/GSE281413/Expression.txt >> ../RawData/GSE281413/clinical.txt
awk -F"," 'ARGIND==1{Sa[$2]=$3"\t"$4"\t"$1}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE281413/GSM_biopsy_site.txt\
           ../RawData/GSE281413/Expression.txt > ../RawData/GSE281413/info.txt
           
#GSE281193
mkdir -p ../RawData/GSE281193
wget -O "../RawData/GSE281193/GSE281193_Q3_norm_all_batch.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE281193&format=file&file=GSE281193%5FQ3%5Fnorm%5Fall%5Fbatch%2Exlsx"

Rscript -e '
    library(readxl)
    df <- read_excel("../RawData/GSE281193/GSE281193_Q3_norm_all_batch.xlsx", sheet = 3)
    write.table(as.data.frame(df), "../RawData/GSE281193/Expression.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'

wget -O ../RawData/GSE281193/GSE281193_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281193/soft/GSE281193_family.soft.gz
gunzip ../RawData/GSE281193/GSE281193_family.soft.gz
awk '
/^\^SAMPLE =/ {
    gsm=$3
}
/!Sample_title = / {
    sub(/!Sample_title = /,"")
    print gsm","$0
}
' ../RawData/GSE281193/GSE281193_family.soft > ../RawData/GSE281193/GSM_biopsy_site.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE281193/clinical.txt
awk -F"," 'ARGIND==1{Sa[substr($6,2)]=$4"_"$5}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){split($i,a,"|");if($i in Sa){group=Sa[$i]}else{group="UnKnown"};print $i"\t"group"\t"a[1]}}'\
           ../RawData/GSE281193/GSM_biopsy_site.txt\
           ../RawData/GSE281193/Expression.txt | awk -F"\t" '{gsub(/ /,"",$2);print $1"\t"$2"\t"$3}' >> ../RawData/GSE281193/clinical.txt
awk -F"," 'ARGIND==1{split($6,a,"|");Sa[substr($6,2)]=$0}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE281193/GSM_biopsy_site.txt\
           ../RawData/GSE281193/Expression.txt > ../RawData/GSE281193/info.txt

#GSE279942
mkdir -p ../RawData/GSE279942
wget -O "../RawData/GSE279942/GSE279942_gene_signature_LARC_matrix.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE279942&format=file&file=GSE279942%5Fgene%5Fsignature%5FLARC%5Fmatrix%2Exlsx"
Rscript -e '
    library(readxl)
    df <- read_excel("../RawData/GSE279942/GSE279942_gene_signature_LARC_matrix.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE279942/Expression.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
wget -O ../RawData/GSE279942/GSE279942_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279942/soft/GSE279942_family.soft.gz
gunzip ../RawData/GSE279942/GSE279942_family.soft.gz
awk '
/^\^SAMPLE =/ {
    if (gsm != ""){print gsm "\t" title "\t" slide "\t" response}
    gsm=$3;title="";slide="";response="";next
}
/^!Sample_title = / {sub(/^!Sample_title = /,"");title=$0}
/^!Sample_characteristics_ch1 = slidename: / {sub(/^!Sample_characteristics_ch1 = slidename: /,"");slide=$0}
/^!Sample_characteristics_ch1 = response: / {sub(/^!Sample_characteristics_ch1 = response: /,"");response=$0}
END {if (gsm != ""){print gsm "\t" title "\t" slide "\t" response}}
' ../RawData/GSE279942/GSE279942_family.soft |\
     awk -F"\t" '{gsub(/ /,"",$4);print $1","$2","$4","$3}' > ../RawData/GSE279942/GSM_biopsy_site.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE279942/clinical.txt
awk -F"," 'ARGIND==1{split($2,a,"_");Sa[$2]=$3"_"a[length(a)]"\t"$4}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE279942/GSM_biopsy_site.txt\
           ../RawData/GSE279942/Expression.txt | awk -F"\t" '{gsub(/ /,"",$2);print $1"\t"$2"\t"$3}' >> ../RawData/GSE279942/clinical.txt
awk -F"," 'ARGIND==1{Sa[$2]=$0}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE279942/GSM_biopsy_site.txt\
           ../RawData/GSE279942/Expression.txt > ../RawData/GSE279942/info.txt


#GSE279576   Visium��ȡ

#GSE278670
mkdir -p ../RawData/GSE278670
wget -O "../RawData/GSE278670/GSE278670_Iyer_WTA_GeoMx_2022_GEO_Initial_Dataset.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278670&format=file&file=GSE278670%5FIyer%5FWTA%5FGeoMx%5F2022%5FGEO%5FInitial%5FDataset%2Exlsx"
Rscript -e '
    library(readxl)
    df <- read_excel("../RawData/GSE278670/GSE278670_Iyer_WTA_GeoMx_2022_GEO_Initial_Dataset.xlsx", sheet = 2)
    expr <- as.data.frame(df)[,c(4,13:ncol(df))]
    expr <- expr[which(!is.na(expr$HUGOSymbol)),]
    expr$HUGOSymbol <- make.unique(expr$HUGOSymbol)
    write.table(expr, "../RawData/GSE278670/Expression.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'

wget -O ../RawData/GSE278670/GSE278670_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278670/soft/GSE278670_family.soft.gz
gunzip ../RawData/GSE278670/GSE278670_family.soft.gz
awk '
/^\^SAMPLE =/ {
    if (gsm != ""){print gsm "\t" title "\t" stage}
    gsm=$3;title="";stage="";next
}
/^!Sample_title = / {sub(/^!Sample_title = /,"");title=$0}
/^!Sample_characteristics_ch1 = disease pathology: / {sub(/^!Sample_characteristics_ch1 = disease pathology: /,"");stage=$0}
END {if (gsm != ""){print gsm "\t" title "\t" stage}}
' ../RawData/GSE278670/GSE278670_family.soft |\
     awk -F"\t" '{split($2,a,"|");print $1","$2","$3"_"a[3]","a[1]}' | awk -F"," '{gsub(/ /,"",$3);print $1","$2","$3","$4}' > ../RawData/GSE278670/GSM_biopsy_site.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE278670/clinical.txt
awk -F"," 'ARGIND==1{Sa[$2]=$3}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){if($i in Sa){group=Sa[$i]}else{group="Unknown"};split($i,b,"|");print $i"\t"group"\t"b[1]}}'\
           ../RawData/GSE278670/GSM_biopsy_site.txt\
           ../RawData/GSE278670/Expression.txt >> ../RawData/GSE278670/clinical.txt
awk -F"," 'ARGIND==1{Sa[$2]=$0}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE278670/GSM_biopsy_site.txt\
           ../RawData/GSE278670/Expression.txt > ../RawData/GSE278670/info.txt

#GSE278436
mkdir -p ../RawData/GSE278436
wget -O "../RawData/GSE278436/matrix.mtx.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278436&format=file&file=GSE278436%5Fgene%5Fcounts%5Fmatrix%2Emtx%2Egz"
wget -O "../RawData/GSE278436/barcodes.tsv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278436&format=file&file=GSE278436%5Fbarcodes%2Etsv%2Egz"
wget -O "../RawData/GSE278436/features.tsv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278436&format=file&file=GSE278436%5Ffeatures%2Etsv%2Egz"
wget -O "../RawData/GSE278436/metadata.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278436&format=file&file=GSE278436%5Fmetadata%2Ecsv%2Egz"
gunzip ../RawData/GSE278436/metadata.csv.gz
Rscript -e '
    library(Seurat)
    seurat <- Read10X(data.dir = "../RawData/GSE278436",gene.column = 1)
    seu <- CreateSeuratObject(counts = seurat,project = "GSE278436")
    meta <- read.csv("../RawData/GSE278436/metadata.csv",row.names = 1,check.names = FALSE)
    common.cells <- intersect(colnames(seu), rownames(meta))
    seu <- subset(seu, cells = common.cells)
    seu <- AddMetaData(seu, meta[common.cells, , drop = FALSE])
    exp <- GetAssayData(seu, layer = "counts")
    write.table(as.data.frame(exp), "../RawData/GSE278436/Expression.txt", sep = "\t", quote = FALSE, row.names = TRUE)
    clinical <- data.frame(sample = rownames(seu@meta.data),group  = seu@meta.data$slide_neuLoc,slide  = seu@meta.data$Slide.Name,row.names = NULL)
    write.table(clinical, "../RawData/GSE278436/clinical.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'

rawsample="GSE278436"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

head -n 1 ../RawData/$rawsample/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/$rawsample/Expression.1.txt
tail -n +2 ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

head -n 1 ../RawData/$rawsample/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/$rawsample/Expression.1.txt
tail -n +2 ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

#GSE277104, seems no group??

#GSE276935
mkdir -p ../RawData/GSE276935
wget -O "../RawData/GSE276935/GSE276935_Processed_Data_File_All_Experiments.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE276935&format=file&file=GSE276935%5FProcessed%5FData%5FFile%5FAll%5FExperiments%2Exlsx"
Rscript -e '
    library(readxl)
    df <- read_excel("../RawData/GSE276935/GSE276935_Processed_Data_File_All_Experiments.xlsx", sheet = 3)
    write.table(as.data.frame(df), "../RawData/GSE276935/Expression.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
echo -e "sample\tgroup\tslide" > ../RawData/GSE276935/clinical.txt
head -n 1 ../RawData/GSE276935/Expression.txt | tr ',' '\n' | tail -n +2 | \
  awk -F"," '{match($0, /^([0-9]+[A-C]?|LN)([A-Za-z ]+)([0-9]{3})(.+)$/, a); print $0"\t"a[2]"_"a[4]"\t"a[1]}' |\
  awk -F"\t" '{gsub(/ /,"",$2);print $1"\t"$2"\t"$3}' >> ../RawData/GSE276935/clinical.txt

#GSE271689 no pkc files

#GSE271285
mkdir -p ../RawData/GSE271285
wget -O "../RawData/GSE271285/GSE271285_JeongLABAllDataWTA.xlsx" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE271285&format=file&file=GSE271285%5FJeongLABAllDataWTA%2Exlsx"
Rscript -e '
    library(readxl)
    df <- read_excel("../RawData/GSE271285/GSE271285_JeongLABAllDataWTA.xlsx", sheet = 4)
    write.table(as.data.frame(df), "../RawData/GSE271285/Expression.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
echo -e "sample\tgroup\tslide" > ../RawData/GSE271285/clinical.txt
head -n 1 ../RawData/GSE271285/Expression.txt | tr ',' '\n' | tail -n +2 |\
  awk -F"|" '{if((NR>=1 && NR<=8) || (NR>=21 && NR<=28)){dis="Normal"};if(NR>=9 && NR<=20){dis="ColonTumor"};if(NF>28){dis="LiverMets"}print $0"\t"dis"_"$3"\t"$1}' >> ../RawData/GSE271285/clinical.txt

#GSE271255
mkdir -p ../RawData/GSE271255
wget -O "../RawData/GSE271255/GSE271255_counts.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE271255&format=file&file=GSE271255%5Fcounts%2Ecsv%2Egz"
wget -O "../RawData/GSE271255/GSE271255_metadata.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE271255&format=file&file=GSE271255%5Fmetadata%2Ecsv%2Egz"
gunzip ../RawData/GSE271255/*
cp "../RawData/GSE271255/GSE271255_counts.csv" "../RawData/GSE271255/Expression.txt"
sed -i 's/,/\t/g' ../RawData/GSE271255/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE271255/clinical.txt
awk -F"," 'ARGIND==1 && (FNR>1){Sa[$1]=$6"\t"$2}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE271255/GSE271255_metadata.csv\
           ../RawData/GSE271255/Expression.txt >> ../RawData/GSE271255/clinical.txt

#GSE265899
mkdir -p ../RawData/GSE265899
wget -O "../RawData/GSE265899/GSE265899_Q3Norm.csv.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE265899&format=file&file=GSE265899%5FQ3Norm%2Ecsv%2Egz"
gunzip "../RawData/GSE265899/GSE265899_Q3Norm.csv.gz"
cp ../RawData/GSE265899/GSE265899_Q3Norm.csv ../RawData/GSE265899/Expression.txt
dos2unix ../RawData/GSE265899/Expression.txt
wget -O ../RawData/GSE265899/GSE265899_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE265nnn/GSE265899/soft/GSE265899_family.soft.gz
gunzip ../RawData/GSE265899/GSE265899_family.soft.gz
awk '
/^\^SAMPLE =/ {
    if (gsm != ""){print gsm "\t" title "\t" group "\t" slide}
    gsm=$3;title="";group="";slide="";next
}
/^!Sample_title = / {sub(/^!Sample_title = /,"");title=$0}
/^!Sample_characteristics_ch1 = group: / {sub(/^!Sample_characteristics_ch1 = group: /,"");group=$0}
/^!Sample_characteristics_ch1 = case: / {sub(/^!Sample_characteristics_ch1 = case: /,"");slide=$0}
END {if (gsm != ""){print gsm "\t" title "\t" group "\t" slide}}
' ../RawData/GSE265899/GSE265899_family.soft |\
     awk -F"\t" '{if($3~"PDL1H"){group=$3"igh"}else{group=$3"ow"};gsub(/ /,"",$4);split($2,a,",");print $1","group","$4","a[1]}' > ../RawData/GSE265899/GSM_biopsy_site.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE265899/clinical.txt     
awk -F"," 'ARGIND==1{match($4, /^[A-Za-z]+/);Sa[$4]=$2"_"substr($4, RSTART, RLENGTH)"\t"$3}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE265899/GSM_biopsy_site.txt\
           ../RawData/GSE265899/Expression.txt |\
           awk -F"\t" '{split($2,a,"_");if(a[1]=="ID"){immune="ImmuneDeficient"}else{immune="ImmuneInflamed"};print $1"\t"immune"_"a[2]"_"a[3]"\t"$3}'\
           >> ../RawData/GSE265899/clinical.txt  

awk -F"," 'ARGIND==1{match($4, /^[A-Za-z]+/);Sa[$4]=$1"\t"$2"_"substr($4, RSTART, RLENGTH)"\t"$3}\
           ARGIND==2 && (FNR==1){for(i=2;i<=NF;i++){print $i"\t"Sa[$i]}}'\
           ../RawData/GSE265899/GSM_biopsy_site.txt\
           ../RawData/GSE265899/Expression.txt > ../RawData/GSE265899/info.txt  
sed -i 's/,/\t/g' ../RawData/GSE265899/Expression.txt
           


#Revise
#GSE316098
dos2unix ../RawData/GSE316098/clinical.txt
dos2unix ../RawData/GSE316098/Expression.txt
head -n 186 ../RawData/GSE316098/clinical.txt > ../RawData/GSE316098/clinical.1.txt
mv ../RawData/GSE316098/clinical.1.txt ../RawData/GSE316098/clinical.txt

#GSE300750
head -n 1 ../RawData/GSE300750/Expression.txt > ../RawData/GSE300750/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/GSE300750/Expression.txt >> ../RawData/GSE300750/Expression.1.txt
mv ../RawData/GSE300750/Expression.1.txt ../RawData/GSE300750/Expression.txt

originalpath=`pwd`
rawsample="GSE300750"
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

#GSE303279
rawsample="GSE303279"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

#GSE302983
rawsample="GSE302983"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

#GSE296330
rawsample="GSE296330"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/GSE296330/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/GSE296330/Expression.1.txt
tail -n +2 ../RawData/GSE296330/Expression.txt >> ../RawData/GSE296330/Expression.1.txt
mv ../RawData/GSE296330/Expression.1.txt ../RawData/GSE296330/Expression.txt

# GSE330215
mkdir -p ../RawData/GSE330215
wget -O ../RawData/GSE330215/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE330215&format=file&file=GSE330215%5Fexpression%5Ffrom%5FDSP%5Fplatform%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE330215/Expression.xlsx > ../RawData/GSE330215/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE330215/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE330215/clinical.txt
paste <(head -n1 ../RawData/GSE330215/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE330215"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE330215/clinical.txt
sed -i '$d' ../RawData/GSE330215/clinical.txt


# GSE328332
mkdir -p ../RawData/GSE328332
wget -O ../RawData/GSE328332/Expression14.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE328332&format=file&file=GSE328332%5Fraw%5Fcount14%2Ecsv%2Egz"
wget -O ../RawData/GSE328332/Expression5.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE328332&format=file&file=GSE328332%5Fraw%5Fcount5%2Ecsv%2Egz"
gunzip ../RawData/GSE328332/Expression14.csv.gz
gunzip ../RawData/GSE328332/Expression5.csv.gz
n14=$(head -1 ../RawData/GSE328332/Expression14.csv | awk -F, '{print NF}')
awk -F, -v n="$n14" 'BEGIN{OFS="\t"} {
    gsub(/\r$/, "");
    if (NF == n+1 && $NF == "") NF=n;
    if (NF < n) for(i=NF+1; i<=n; i++) $i="";
    for(i=1; i<=NF; i++) gsub(/"/, "", $i);
    print
}' ../RawData/GSE328332/Expression14.csv > ../RawData/GSE328332/Expression14_fixed.tsv
n5=$(head -1 ../RawData/GSE328332/Expression5.csv | awk -F, '{print NF}')
awk -F, -v n="$n5" 'BEGIN{OFS="\t"} {
    gsub(/\r$/, "");
    if (NF == n+1 && $NF == "") NF=n;
    if (NF < n) for(i=NF+1; i<=n; i++) $i="";
    for(i=1; i<=NF; i++) gsub(/"/, "", $i);
    print
}' ../RawData/GSE328332/Expression5.csv > ../RawData/GSE328332/Expression5_fixed.tsv
paste -d "\t" <(head -n1 ../RawData/GSE328332/Expression14_fixed.tsv) <(head -n1 ../RawData/GSE328332/Expression5_fixed.tsv | cut -f2-) > ../RawData/GSE328332/Expression.txt
paste -d "\t" <(tail -n +2 ../RawData/GSE328332/Expression14_fixed.tsv) <(tail -n +2 ../RawData/GSE328332/Expression5_fixed.tsv | cut -f2-) >> ../RawData/GSE328332/Expression.txt
head -n1 ../RawData/GSE328332/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++) print $i}' > ../RawData/GSE328332/sample_list.tmp
awk -F"\t" '$7=="GSE328332" {print $2"\t"$6}' STsamples1.txt > ../RawData/GSE328332/group_slide.tmp
echo -e "sample\tgroup\tslide" > ../RawData/GSE328332/clinical.txt
paste -d "\t" ../RawData/GSE328332/sample_list.tmp ../RawData/GSE328332/group_slide.tmp >> ../RawData/GSE328332/clinical.txt
rm -f ../RawData/GSE328332/sample_list.tmp ../RawData/GSE328332/group_slide.tmp ../RawData/GSE328332/Expression14_fixed.tsv ../RawData/GSE328332/Expression5_fixed.tsv

# GSE327983
mkdir -p ../RawData/GSE327983
wget -O ../RawData/GSE327983/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE327983&format=file&file=GSE327983%5FWTAONLY%5FExp%5FNEW%2Ecsv%2Egz"
gunzip ../RawData/GSE327983/Expression.txt.gz
head -n 1 ../RawData/GSE327983/Expression.txt | awk -F"," '{out="Gene";for(i=2;i<=NF;i++){out=out"\t"$i};print out}' > ../RawData/GSE327983/Expression1.txt
tail -n +2 ../RawData/GSE327983/Expression.txt | sed 's/,/\t/g' >> ../RawData/GSE327983/Expression1.txt
mv ../RawData/GSE327983/Expression1.txt ../RawData/GSE327983/Expression.txt
head -n 3 ../RawData/GSE327983/Expression.txt | cut -f1-3
echo -e "sample\tgroup\tslide" > ../RawData/GSE327983/clinical.txt
paste <(head -n1 ../RawData/GSE327983/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++) print $i}' | tr -d '"') <(awk -F"\t" '$7=="GSE327983"{print $2"\t"$6}' STsamples1.txt | tr -d '"') >> ../RawData/GSE327983/clinical.txt

#GSE320059
mkdir -p ../RawData/GSE320059
wget -O "../RawData/GSE320059/Hs_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE320059&format=file&file=GSE320059%5FHs%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE320059/GSE320059_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE320059&format=file"
tar -xvf ../RawData/GSE320059/GSE320059_RAW.tar -C ../RawData/GSE320059/
gunzip ../RawData/GSE320059/*gz
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE320059/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE320059/Hs_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE320059/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE320059/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE320059/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE320059/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
echo -e "sample\tgroup\tslide" > ../RawData/GSE320059/clinical.txt
paste <(head -n1 ../RawData/GSE320059/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE320059"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE320059/clinical.txt

# GSE319968
mkdir -p ../RawData/GSE319968
wget -O ../RawData/GSE319968/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE319968&format=file&file=GSE319968%5F21019220628%5Frawdata%2Ecsv%2Egz"
gunzip ../RawData/GSE319968/Expression.txt.gz
head -n 1 ../RawData/GSE319968/Expression.txt | awk -F"," '{out="Gene";for(i=2;i<=NF;i++){out=out"\t"$i};print out}' > ../RawData/GSE319968/Expression1.txt
tail -n +2 ../RawData/GSE319968/Expression.txt | sed 's/,/\t/g' >> ../RawData/GSE319968/Expression1.txt
mv ../RawData/GSE319968/Expression1.txt ../RawData/GSE319968/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE319968/clinical.txt
paste <(head -n1 ../RawData/GSE319968/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++) print $i}' | tr -d '"') <(awk -F"\t" '$7=="GSE319968"{print $2"\t"$6}' STsamples1.txt | tr -d '"') >> ../RawData/GSE319968/clinical.txt

# GSE317400
mkdir -p ../RawData/GSE317400
wget -O ../RawData/GSE317400/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE317400&format=file&file=GSE317400%5FGEO%5Fprocessed%2Etxt%2Egz"
gunzip ../RawData/GSE317400/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE317400/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE317400/clinical.txt
paste <(head -n1 ../RawData/GSE317400/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE317400"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE317400/clinical.txt

# GSE310215
mkdir -p ../RawData/GSE310215
wget -O ../RawData/GSE310215/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE310215&format=file&file=GSE310215%5FTargetCountMatrix%5FallSamples%2Etxt%2Egz"
gunzip ../RawData/GSE310215/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE310215/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE310215/clinical.txt
paste <(head -n1 ../RawData/GSE310215/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE310215"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE310215/clinical.txt

# GSE316098
mkdir -p ../RawData/GSE316098
wget -O ../RawData/GSE316098/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE316098&format=file&file=GSE316098%5Fraw%2Ecounts%2Elinear%2Eformatted%2Ecsv%2Egz"
gunzip ../RawData/GSE316098/Expression.txt.gz
head -n 1 ../RawData/GSE316098/Expression.txt | awk -F"," '{out="Gene";for(i=2;i<=NF;i++){out=out"\t"$i};print out}' > ../RawData/GSE316098/Expression1.txt
tail -n +2 ../RawData/GSE316098/Expression.txt | sed 's/,/\t/g' >> ../RawData/GSE316098/Expression1.txt
mv ../RawData/GSE316098/Expression1.txt ../RawData/GSE316098/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE316098/clinical.txt
paste <(head -n1 ../RawData/GSE316098/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++) print $i}' | tr -d '"') <(awk -F"\t" '$7=="GSE316098"{print $2"\t"$6}' STsamples1.txt | tr -d '"') >> ../RawData/GSE316098/clinical.txt
dos2unix ../RawData/GSE316098/clinical.txt
dos2unix ../RawData/GSE316098/Expression.txt
head -n 186 ../RawData/GSE316098/clinical.txt > ../RawData/GSE316098/clinical.1.txt
mv ../RawData/GSE316098/clinical.1.txt ../RawData/GSE316098/clinical.txt

# GSE300750
mkdir -p ../RawData/GSE300750
wget -O ../RawData/GSE300750/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE300750&format=file&file=GSE300750%5F5%5Fpercent%5Ftarget%5FQ3%5FNorm%5FRead%5FCounts%5FDESeq%5F2%2Exlsx"
xlsx2csv -s 3 -d $'\t' ../RawData/GSE300750/Expression.xlsx > ../RawData/GSE300750/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE300750/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE300750/clinical.txt
paste <(head -n1 ../RawData/GSE300750/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE300750"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE300750/clinical.txt
head -n 1 ../RawData/GSE300750/Expression.txt > ../RawData/GSE300750/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/GSE300750/Expression.txt >> ../RawData/GSE300750/Expression.1.txt
mv ../RawData/GSE300750/Expression.1.txt ../RawData/GSE300750/Expression.txt

originalpath=`pwd`
rawsample="GSE300750"
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath



# GSE303279
mkdir -p ../RawData/GSE303279
wget -O ../RawData/GSE303279/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE303279&format=file&file=GSE303279%5Flog2%2Enormalized%2Ecounts%2Etxt%2Egz"
gunzip ../RawData/GSE303279/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE303279/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE303279/clinical.txt
paste <(head -n1 ../RawData/GSE303279/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE303279"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE303279/clinical.txt
rawsample="GSE303279"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

# GSE313748
mkdir -p ../RawData/GSE313748
wget -O "../RawData/GSE313748/Hs_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE313748&format=file&file=GSE313748%5FHs%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE313748/GSE313748_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE313748&format=file"
tar -xvf ../RawData/GSE313748/GSE313748_RAW.tar -C ../RawData/GSE313748/
gunzip ../RawData/GSE313748/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE313748/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE313748", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE313748/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE313748/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE313748/Hs_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE313748/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE313748/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE313748/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE313748/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM9373465_/ && $i !~ /^GSM9373560_/ && $i !~ /^GSM9373656_/ && 
           $i !~ /^GSM9373678_/ && $i !~ /^GSM9373774_/ && $i !~ /^GSM9373869_/ && 
           $i !~ /^GSM9373965_/ && $i !~ /^GSM9374050_/ && $i !~ /^GSM9374145_/ ) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE313748/Expression.txt > ../RawData/GSE313748/Expression_filtered.txt
mv ../RawData/GSE313748/Expression_filtered.txt ../RawData/GSE313748/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE313748/clinical.txt
paste <(head -n1 ../RawData/GSE313748/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE313748"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE313748/clinical.txt

# GSE300490
mkdir -p ../RawData/GSE300490
wget -O "../RawData/GSE300490/Hs_R_NGS_CTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE300490&format=file&file=GSE300490%5FHs%5FR%5FNGS%5FCTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE300490/GSE300490_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE300490&format=file"
tar -xvf ../RawData/GSE300490/GSE300490_RAW.tar -C ../RawData/GSE300490/
gunzip ../RawData/GSE300490/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE300490/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE300490", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE300490/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE300490/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE300490/Hs_R_NGS_CTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE300490/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE300490/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE300490/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE300490/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM9062619_/ && $i !~ /^GSM9062673_/ && $i !~ /^GSM9062741_/ && 
           $i !~ /^GSM9062828_/ && $i !~ /^GSM9062923_/ && $i !~ /^GSM9062931_/ && 
           $i !~ /^GSM9062995_/ && $i !~ /^GSM9063073_/ && $i !~ /^GSM9063116_/ && $i !~ /^GSM9063199_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE300490/Expression.txt > ../RawData/GSE300490/Expression_filtered.txt
mv ../RawData/GSE300490/Expression_filtered.txt ../RawData/GSE300490/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE300490/clinical.txt
paste <(head -n1 ../RawData/GSE300490/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE300490"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE300490/clinical.txt

# GSE302983
mkdir -p ../RawData/GSE302983
wget -O "../RawData/GSE302983/Mm_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE302983&format=file&file=GSE302983%5FMm%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE302983/GSE302983_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE302983&format=file"
tar -xvf ../RawData/GSE302983/GSE302983_RAW.tar -C ../RawData/GSE302983/
gunzip ../RawData/GSE302983/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE302983/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE302983", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE302983/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE302983/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE302983/Mm_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE302983/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE302983/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE302983/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE302983/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM9115138_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE302983/Expression.txt > ../RawData/GSE302983/Expression_filtered.txt
mv ../RawData/GSE302983/Expression_filtered.txt ../RawData/GSE302983/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE302983/clinical.txt
paste <(head -n1 ../RawData/GSE302983/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE302983"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE302983/clinical.txt
rawsample="GSE302983"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath

# GSE309894
mkdir -p ../RawData/GSE309894
wget -O "../RawData/GSE309894/Hs_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE309894&format=file&file=GSE309894%5FHs%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE309894/GSE309894_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE309894&format=file"
tar -xvf ../RawData/GSE309894/GSE309894_RAW.tar -C ../RawData/GSE309894/
gunzip ../RawData/GSE309894/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE309894/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE309894", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE309894/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE309894/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE309894/Hs_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE309894/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE309894/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE309894/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE309894/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'

echo -e "sample\tgroup\tslide" > ../RawData/GSE309894/clinical.txt
paste <(head -n1 ../RawData/GSE309894/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE309894"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE309894/clinical.txt

# GSE300414
mkdir -p ../RawData/GSE300414
wget -O "../RawData/GSE300414/Hs_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE300414&format=file&file=GSE300414%5FHs%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE300414/GSE300414_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE300414&format=file"
tar -xvf ../RawData/GSE300414/GSE300414_RAW.tar -C ../RawData/GSE300414/
gunzip ../RawData/GSE300414/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE300414/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE300414", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE300414/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE300414/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE300414/Hs_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE300414/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE300414/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE300414/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE300414/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM9059886_/ && $i !~ /^GSM9059961_/ && $i !~ /^GSM9060020_/ && 
           $i !~ /^GSM9060104_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE300414/Expression.txt > ../RawData/GSE300414/Expression_filtered.txt
mv ../RawData/GSE300414/Expression_filtered.txt ../RawData/GSE300414/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE300414/clinical.txt
paste <(head -n1 ../RawData/GSE300414/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE300414"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE300414/clinical.txt

# GSE299775
mkdir -p ../RawData/GSE299775
wget -O ../RawData/GSE299775/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE299775&format=file&file=GSE299775%5Fcount%5Fmatrix%2Etxt%2Egz"
gunzip ../RawData/GSE299775/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE299775/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE299775/clinical.txt
paste <(head -n1 ../RawData/GSE299775/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE299775"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE299775/clinical.txt

# GSE299880
mkdir -p ../RawData/GSE299880
wget -O ../RawData/GSE299880/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE299880&format=file&file=GSE299880%5FWTAONLY%5FExp%5FNEW%2Ecsv%2Egz"
gunzip ../RawData/GSE299880/Expression.txt.gz
head -n 1 ../RawData/GSE299880/Expression.txt | awk -F"," '{out="Gene";for(i=2;i<=NF;i++){out=out"\t"$i};print out}' > ../RawData/GSE299880/Expression1.txt
tail -n +2 ../RawData/GSE299880/Expression.txt | sed 's/,/\t/g' >> ../RawData/GSE299880/Expression1.txt
mv ../RawData/GSE299880/Expression1.txt ../RawData/GSE299880/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE299880/clinical.txt
paste <(head -n1 ../RawData/GSE299880/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++) print $i}' | tr -d '"') <(awk -F"\t" '$7=="GSE299880"{print $2"\t"$6}' STsamples1.txt | tr -d '"') >> ../RawData/GSE299880/clinical.txt

# GSE296330
mkdir -p ../RawData/GSE296330
wget -O ../RawData/GSE296330/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE296330&format=file&file=GSE296330%5FCCA%5Ffibroblasts%5Ffor%5Fanalysis%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE296330/Expression.xlsx | awk -F'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) if(i!=5) printf "%s%s",$i, (i==NF?"\t":"\t"); print ""}' | sed -e 's/|//g' -e 's/ \+/_/g' > ../RawData/GSE296330/Expression.txt
# 删除第2行
sed -i '2d' ../RawData/GSE296330/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE296330/clinical.txt
paste <(head -n1 ../RawData/GSE296330/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE296330"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE296330/clinical.txt
rawsample="GSE296330"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/GSE296330/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/GSE296330/Expression.1.txt
tail -n +2 ../RawData/GSE296330/Expression.txt >> ../RawData/GSE296330/Expression.1.txt
mv ../RawData/GSE296330/Expression.1.txt ../RawData/GSE296330/Expression.txt

# GSE296955
mkdir -p ../RawData/GSE296955
wget -O ../RawData/GSE296955/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE296955&format=file&file=GSE296955%5FProcessedDataGemCis%2Exlsx"
xlsx2csv -s 3 -d $'\t' ../RawData/GSE296955/Expression.xlsx > ../RawData/GSE296955/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE296955/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE296955/clinical.txt
paste <(head -n1 ../RawData/GSE296955/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE296955"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE296955/clinical.txt

# GSE292946
mkdir -p ../RawData/GSE292946
wget -O "../RawData/GSE292946/Hs_R_NGS_WTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292946&format=file&file=GSE292946%5FHs%5FR%5FNGS%5FWTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE292946/GSE292946_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292946&format=file"
tar -xvf ../RawData/GSE292946/GSE292946_RAW.tar -C ../RawData/GSE292946/
gunzip ../RawData/GSE292946/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE292946/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE292946", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE292946/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE292946/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE292946/Hs_R_NGS_WTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE292946/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE292946/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE292946/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE292946/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM8869943_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE292946/Expression.txt > ../RawData/GSE292946/Expression_filtered.txt
mv ../RawData/GSE292946/Expression_filtered.txt ../RawData/GSE292946/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE292946/clinical.txt
paste <(head -n1 ../RawData/GSE292946/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE292946"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE292946/clinical.txt

# GSE292226
mkdir -p ../RawData/GSE292226
wget -O "../RawData/GSE292226/Cf_R_NGS_CA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292226&format=file&file=GSE292226%5FCf%5FR%5FNGS%5FCA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE292226/GSE292226_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292226&format=file"
tar -xvf ../RawData/GSE292226/GSE292226_RAW.tar -C ../RawData/GSE292226/
gunzip ../RawData/GSE292226/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE292226/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE292226", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE292226/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE292226/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE292226/Cf_R_NGS_CA_v1.0.pkc")
    annotation_file <- "../RawData/GSE292226/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE292226/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE292226/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE292226/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM8853714_/ && $i !~ /^GSM8853618_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE292226/Expression.txt > ../RawData/GSE292226/Expression_filtered.txt
mv ../RawData/GSE292226/Expression_filtered.txt ../RawData/GSE292226/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE292226/clinical.txt
paste <(head -n1 ../RawData/GSE292226/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE292226"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE292226/clinical.txt

# GSE295028
mkdir -p ../RawData/GSE295028
wget -O ../RawData/GSE295028/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE295028&format=file&file=GSE295028%5Fexpression%5Ffrom%5FDSP%5Fplatform%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE295028/Expression.xlsx > ../RawData/GSE295028/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE295028/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE295028/clinical.txt
paste <(head -n1 ../RawData/GSE295028/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE295028"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE295028/clinical.txt
sed -i '$d' ../RawData/GSE295028/clinical.txt

# GSE261347
mkdir -p ../RawData/GSE261347
wget -O ../RawData/GSE261347/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261347&format=file&file=GSE261347%5FGEO%5Fprocessed%2Etxt%2Egz"
gunzip ../RawData/GSE261347/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE261347/Expression.txt
# 去除第1,2,4,5,6,7列（保留第3列作为基因名，第8列及之后作为表达量）
awk -F"\t" 'BEGIN{OFS="\t"} {line=$3; for(i=8;i<=NF;i++) line=line"\t"$i; print line}' ../RawData/GSE261347/Expression.txt > ../RawData/GSE261347/Expression_tmp.txt
mv ../RawData/GSE261347/Expression_tmp.txt ../RawData/GSE261347/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE261347/clinical.txt
paste <(head -n1 ../RawData/GSE261347/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE261347"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE261347/clinical.txt
sed -i 's/"//g' ../RawData/GSE261347/*.txt
rawsample="GSE261347"
originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/$rawsample/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/$rawsample/Expression.1.txt
tail -n +2 ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

# GSE261348
mkdir -p ../RawData/GSE261348
wget -O ../RawData/GSE261348/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261348&format=file&file=GSE261348%5FIMfirst%5FDSP%5Frawcounts%2Exlsx"
xlsx2csv -s 3 -d $'\t' ../RawData/GSE261348/Expression.xlsx > ../RawData/GSE261348/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE261348/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE261348/clinical.txt
paste <(head -n1 ../RawData/GSE261348/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE261348"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE261348/clinical.txt

# GSE247185
mkdir -p ../RawData/GSE247185
wget -O ../RawData/GSE247185/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE247185&format=file&file=GSE247185%5FBladderProject%5FRawData%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE247185/Expression.xlsx > ../RawData/GSE247185/Expression.txt
# 删除第一行
sed -i -e 's/|//g' -e 's/ \+/_/g' -e '1d' ../RawData/GSE247185/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE247185/clinical.txt
paste <(head -n1 ../RawData/GSE247185/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE247185"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE247185/clinical.txt
head -n 1 ../RawData/GSE247185/Expression.txt | awk -F"\t" '{print substr($0,2)}' > ../RawData/GSE247185/Expression.1.txt
tail -n +30 ../RawData/GSE247185/Expression.txt >> ../RawData/GSE247185/Expression.1.txt
mv ../RawData/GSE247185/Expression.1.txt ../RawData/GSE247185/Expression.txt

# GSE253937
mkdir -p ../RawData/GSE253937
wget -O ../RawData/GSE253937/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE253937&format=file&file=GSE253937%5FProcessed%5Fdata%5FMC38%5FTumor%5FGeoMX%5FNT%2BKD%5Ftumor%5Fdata%5Fanalysis%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE253937/Expression.xlsx > ../RawData/GSE253937/Expression.txt
# 去除第16列
cut -f1-15,17- ../RawData/GSE253937/Expression.txt > ../RawData/GSE253937/Expression.tmp
mv ../RawData/GSE253937/Expression.tmp ../RawData/GSE253937/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE253937/clinical.txt
paste <(head -n1 ../RawData/GSE253937/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE253937"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE253937/clinical.txt
rawsample="GSE253937"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/GSE253937/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/GSE253937/Expression.1.txt
tail -n +2 ../RawData/GSE253937/Expression.txt >> ../RawData/GSE253937/Expression.1.txt
mv ../RawData/GSE253937/Expression.1.txt ../RawData/GSE253937/Expression.txt

# GSE253938
mkdir -p ../RawData/GSE253938
wget -O ../RawData/GSE253938/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE253938&format=file&file=GSE253938%5FProcessed%5Fdata%5FMC38%5FTumor%5FGeoMX%5FNT%2BKD%5Ftumor%5Fdata%5Fanalysis%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE253938/Expression.xlsx > ../RawData/GSE253938/Expression.txt
# 去除第16列
cut -f1-15,17- ../RawData/GSE253938/Expression.txt > ../RawData/GSE253938/Expression.tmp
mv ../RawData/GSE253938/Expression.tmp ../RawData/GSE253938/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE253938/clinical.txt
paste <(head -n1 ../RawData/GSE253938/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE253938"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE253938/clinical.txt

# GSE253962
mkdir -p ../RawData/GSE253962
wget -O "../RawData/GSE253962/Hs_R_NGS_CTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE253962&format=file&file=GSE253962%5FHs%5FR%5FNGS%5FCTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE253962/GSE253962_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE253962&format=file"
tar -xvf ../RawData/GSE253962/GSE253962_RAW.tar -C ../RawData/GSE253962/
gunzip ../RawData/GSE253962/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE253962/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE253962", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE253962/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE253962/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE253962/Hs_R_NGS_CTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE253962/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE253962/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE253962/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE253962/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM8031057_/ && $i !~ /^GSM8031058_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE253962/Expression.txt > ../RawData/GSE253962/Expression_filtered.txt
mv ../RawData/GSE253962/Expression_filtered.txt ../RawData/GSE253962/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE253962/clinical.txt
paste <(head -n1 ../RawData/GSE253962/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE253962"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE253962/clinical.txt

# GSE255328
mkdir -p ../RawData/GSE255328
wget -O ../RawData/GSE255328/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE255328&format=file&file=GSE255328%5FProcessed%5Fdata%5Ffile%5F43%5FAOIs%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE255328/Expression.xlsx > ../RawData/GSE255328/Expression.txt
wget -O ../RawData/GSE255328/Annotation.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE255328&format=file&file=GSE255328%5FAnnotation%5F43%5FAOIs%2Exlsx"
xlsx2csv -d $'\t' ../RawData/GSE255328/Annotation.xlsx > ../RawData/GSE255328/Annotation.txt
exp_data=$(grep -A 100000 "Target Name" ../RawData/GSE255328/Expression.txt | tail -n +2)
{
    echo -ne "Target Name\t"
    tail -n +2 ../RawData/GSE255328/Annotation.txt | cut -f1 | tr '\n' '\t'
    echo ""
    echo "$exp_data"
} > ../RawData/GSE255328/Expression.txt
sed -i 's/[[:space:]]*$//' ../RawData/GSE255328/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE255328/clinical.txt
paste <(head -n1 ../RawData/GSE255328/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE255328"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE255328/clinical.txt
head -n 44 ../RawData/GSE255328/clinical.txt > ../RawData/GSE255328/clinical.1.txt
mv ../RawData/GSE255328/clinical.1.txt ../RawData/GSE255328/clinical.txt

# GSE260598
mkdir -p ../RawData/GSE260598
wget -O ../RawData/GSE260598/Expression.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE260598&format=file&file=GSE260598%5FGEO%5Fprocessed%2Etxt%2Egz"
gunzip ../RawData/GSE260598/Expression.txt.gz
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE260598/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE260598/clinical.txt
paste <(head -n1 ../RawData/GSE260598/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE260598"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE260598/clinical.txt

# GSE261345
mkdir -p ../RawData/GSE261345
wget -O ../RawData/GSE261345/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261345&format=file&file=GSE261345%5FCANTABRICO%5FDSP%5Frawcounts%2Exlsx"
xlsx2csv -s 3 -d $'\t' ../RawData/GSE261345/Expression.xlsx > ../RawData/GSE261345/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE261345/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE261345/clinical.txt
paste <(head -n1 ../RawData/GSE261345/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE261345"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE261345/clinical.txt


# GSE243408
mkdir -p ../RawData/GSE243408
wget -O "../RawData/GSE243408/Hs_R_NGS_CTA_v1.0.pkc.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243408&format=file&file=GSE243408%5FHs%5FR%5FNGS%5FCTA%5Fv1%2E0%2Epkc%2Egz"
wget -O "../RawData/GSE243408/GSE243408_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243408&format=file"
tar -xvf ../RawData/GSE243408/GSE243408_RAW.tar -C ../RawData/GSE243408/
gunzip ../RawData/GSE243408/*gz
Rscript -e '
    library(openxlsx)
    dcc_files <- dir("../RawData/GSE243408/", pattern = ".dcc$")
    gsm_ids <- gsub("_.*", "", dcc_files)
    st <- read.delim("STsamples1.txt", sep = "\t", stringsAsFactors = FALSE)
    st_sub <- st[st$OrginalDataset == "GSE243408", c("OrginalSample", "Site")]
    colnames(st_sub) <- c("GSM", "Group")
    annot <- data.frame(
        Sample_ID = dcc_files,
        GSM = gsm_ids,
        stringsAsFactors = FALSE
    )
    annot <- merge(annot, st_sub, by = "GSM", all.x = TRUE)
    annot$Group[is.na(annot$Group)] <- "Unknown"
    write.xlsx(annot[, c("Sample_ID", "Group")], 
                "../RawData/GSE243408/annotation.xlsx")
'
Rscript -e '
    library(GeomxTools)
    library(Biobase)
    library(dplyr)
    library(tibble)
    dcc_files <- dir("../RawData/GSE243408/", pattern = ".dcc$", full.names = TRUE)
    pkc_files <- c("../RawData/GSE243408/Hs_R_NGS_CTA_v1.0.pkc")
    annotation_file <- "../RawData/GSE243408/annotation.xlsx" #Sample_ID: DCC Filename, Sample_info: Group
    data <- readNanoStringGeoMxSet(dccFiles = dcc_files, pkcFiles = pkc_files, phenoDataFile=annotation_file)
    counts_matrix <- assayDataElement(data, "exprs")
    genes <- fData(data)
    genes <- genes[rownames(counts_matrix),]
    rownames(counts_matrix) <- genes$TargetName
    counts_df <- as.data.frame(counts_matrix)
    counts_df$Gene <- rownames(counts_matrix)
    counts_summarized <- counts_df %>% group_by(Gene) %>% summarise(across(where(is.numeric), mean))
    write.table(as.matrix(counts_summarized),"../RawData/GSE243408/Expression.txt",sep="\t",quote = FALSE, row.names = FALSE)
    
    library(readxl)
    df <- read_excel("../RawData/GSE243408/annotation.xlsx", sheet = 1)
    write.table(as.data.frame(df), "../RawData/GSE243408/annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
'
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i !~ /^GSM8031057_/ && $i !~ /^GSM8031058_/) {
            cols[++c]=i
        }
    }
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}
NR>1 {
    for(i=1;i<=c;i++) printf "%s%s", $cols[i], (i==c ? ORS : OFS)
}' OFS='\t' ../RawData/GSE243408/Expression.txt > ../RawData/GSE243408/Expression_filtered.txt
mv ../RawData/GSE243408/Expression_filtered.txt ../RawData/GSE243408/Expression.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE243408/clinical.txt
paste <(head -n1 ../RawData/GSE243408/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r|"/,""); print $i}}') <(awk -F"\t" '$7=="GSE243408"{gsub(/\r|"/,""); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE243408/clinical.txt


# GSE240531
mkdir -p ../RawData/GSE240531
wget -O ../RawData/GSE240531/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE240531&format=file&file=GSE240531%5F112AOIs%5Fn1222%5F5LOQ%5FFiltered%5Fall%5Flinear%5Fnormalized%5Fcounts%2Exlsx"
xlsx2csv -s 2 -d $'\t' ../RawData/GSE240531/Expression.xlsx > ../RawData/GSE240531/Expression1.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE240531/Expression1.txt
awk -F"\t" '
BEGIN { OFS = "\t" }
# 处理前5行元数据：逐列拼接非空内容、跳过空白单元格
NR <= 5 {
    for (i = 1; i <= NF; i++) {
        val = $i
        gsub(/^[ \t]+|[ \t]+$/, "", val) # 去除首尾空白
        if (val != "") {
            if (header[i] == "") {
                header[i] = val
            } else {
                header[i] = header[i] "_" val
            }
        }
    }
    next
}
NR == 6 {
    printf "%s", $1  # 第一列固定保留 Target Name
    for (i = 2; i <= NF; i++) {
        printf "%s%s", OFS, header[i]
    }
    printf "\n"
    next
}

{ print $0 }
' ../RawData/GSE240531/Expression1.txt > ../RawData/GSE240531/Expression.txt
rm ../RawData/GSE240531/Expression1.txt

echo -e "sample\tgroup\tslide" > ../RawData/GSE240531/clinical.txt
paste <(head -n1 ../RawData/GSE240531/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE240531"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE240531/clinical.txt


# GSE240138
mkdir -p ../RawData/GSE240138
wget -O ../RawData/GSE240138/Expression.xlsx "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE240138&format=file&file=GSE240138%5F112AOIs%5Fn1222%5F5LOQ%5FFiltered%5Fall%5Flinear%5Fnormalized%5Fcounts%2Exlsx"
xlsx2csv -s 1 -d $'\t' ../RawData/GSE240138/Expression.xlsx > ../RawData/GSE240138/Expression.txt
sed -i -e 's/|//g' -e 's/ \+/_/g' ../RawData/GSE240138/Expression.txt
echo -e "sample\tgroup\tslide" > ../RawData/GSE240138/clinical.txt
paste <(head -n1 ../RawData/GSE240138/Expression.txt | awk -F"\t" '{for(i=2;i<=NF;i++){gsub(/\r/,""); print $i}}') <(awk -F"\t" '$7=="GSE240138"{gsub(/\r/,"",$2); gsub(/\r/,"",$6); print $2"\t"$6}' STsamples1.txt) >> ../RawData/GSE240138/clinical.txt

#GSE208536
mkdir -p ../RawData/GSE208536
wget -O ../RawData/GSE208536/GSE208536_Processed_data.xlsx \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE208536&format=file&file=GSE208536%5FProcessed%5Fdata%2Exlsx"
python - <<'PY'
import pandas as pd
df = pd.read_excel(
    "../RawData/GSE208536/GSE208536_Processed_data.xlsx",
    sheet_name=0
)
df.to_csv(
    "../RawData/GSE208536/Expression.txt",
    sep="\t",
    index=False
)
PY
cat > ../RawData/GSE208536/roi_map.txt << 'EOF'
GSM6347941	TMA1_001
GSM6347942	TMA1_002
GSM6347943	TMA1_007
GSM6347944	TMA1_008
GSM6347945	TMA1_009
GSM6347946	TMA1_010
GSM6347947	TMA1_013
GSM6347948	TMA1_014
GSM6347949	TMA1_017
GSM6347950	TMA1_018
GSM6347951	TMA1_022
GSM6347952	TMA1_023
GSM6347953	TMA2_002
GSM6347954	TMA2_008
GSM6347955	TMA2_009
GSM6347956	TMA2_023
GSM6347957	TMA1_003
GSM6347958	TMA1_004
GSM6347959	TMA1_015
GSM6347960	TMA1_019
GSM6347961	TMA2_001
GSM6347962	TMA2_005
GSM6347963	TMA2_007
GSM6347964	TMA2_010
GSM6347965	TMA2_011
GSM6347966	TMA2_015
GSM6347967	TMA2_016
GSM6347968	TMA2_017
GSM6347969	TMA2_018
GSM6347970	TMA2_019
GSM6347971	TMA2_020
GSM6347972	TMA2_021
GSM6347973	TMA1_005
GSM6347974	TMA1_006
GSM6347975	TMA1_011
GSM6347976	TMA1_012
GSM6347977	TMA1_016
GSM6347978	TMA1_020
GSM6347979	TMA1_021
GSM6347980	TMA1_024
GSM6347981	TMA1_025
GSM6347982	TMA2_003
GSM6347983	TMA2_004
GSM6347984	TMA2_006
GSM6347985	TMA2_012
GSM6347986	TMA2_013
GSM6347987	TMA2_014
GSM6347988	TMA2_022
EOF
echo -e "sample\tgroup\tslide" > clinical.txt
awk -F'\t' '
BEGIN{
    while((getline < "roi_map.txt")>0){
        gsub(/\r/,"",$1)
        gsub(/\r/,"",$2)
        roi2gsm[$2]=$1
    }
    while((getline < "/data_d/WSJ/SpatialMetsDB/RawData/STsamples2.txt")>0){
        gsub(/\r/,"",$6)
        gsub(/\r/,"",$8)
        gsub(/\r/,"",$11)

        if($7=="GSE208536"){
            group[$8]=$11
            slide[$8]=$6
        }
    }
}
{
    gsm=roi2gsm[$2]

    gsub(/\r/,"",gsm)

    if(gsm!=""){
        print $2"\t"group[gsm]"\t"slide[gsm]
    }
}
' roi_map.txt >> clinical.txt

#GSE225456
mkdir -p ../RawData/GSE225456
wget -O ../RawData/GSE225456/Expression.txt.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE225456&format=file&file=GSE225456%5Fprocessed%5Fdata%2Etxt%2Egz"
gunzip -f ../RawData/GSE225456/Expression.txt.gz
#GSE225456
rawsample="GSE225456"
originalpath=`pwd`
head -n 1 ../RawData/$rawsample/Expression.txt > ../RawData/$rawsample/Expression.1.txt
awk 'BEGIN{FS=OFS="\t"}\
   ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt

originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/$rawsample/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/$rawsample/Expression.1.txt
tail -n +2 ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt
awk -F'\t' '
{
    for(i=2;i<=NF;i++){
        sample=$i
        gsub(/\r/,"",sample)
        if(sample=="") continue

        n=split(sample,a,/ *\| */)
        slide=a[2]
        slide_num=slide+0

        if(slide_num>=17 && slide_num<=20)
            group="Glut1_knockdown_MC38_CAFs"
        else if(slide_num>=21 && slide_num<=26)
            group="Control_MC38_CAFs"
        else
            group="Unknown"

        print sample "\t" group "\t" slide
    }
}' >> ../RawData/GSE225456/clinical.txt

head ../RawData/GSE225456/clinical.txt

#GSE208536
mkdir -p ../RawData/GSE226562
cd ../RawData/GSE226562
wget -O GSE226562_Wrenn_Apfelbaum_CHLA10_xenograft_DSP.xlsx \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226562&format=file&file=GSE226562%5FWrenn%5FApfelbaum%5FCHLA10%5Fxenograft%5FDSP%2Exlsx"
python - <<'PY'
import pandas as pd
xlsx="/data_d/WSJ/SpatialMetsDB/RawData/GSE226562/GSE226562_Wrenn_Apfelbaum_CHLA10_xenograft_DSP.xlsx"
sheet=pd.ExcelFile(xlsx).sheet_names[-1]
df=pd.read_excel(xlsx,sheet_name=sheet)
df=df.iloc[:,:-4]
df.to_csv(
    "/data_d/WSJ/SpatialMetsDB/RawData/GSE226562/Expression.txt",
    sep="\t",
    index=False
)
PY
echo -e "sample\tgroup\tslide" > /data_d/WSJ/SpatialMetsDB/RawData/GSE226562/clinical.txt
paste <(head -n1 /data_d/WSJ/SpatialMetsDB/RawData/GSE226562/Expression.txt | \
awk -F"\t" '{for(i=2;i<=NF;i++) print $i}') \
<(awk -F"\t" '$7=="GSE226562"{print $11"\t"$6}' \
/data_d/WSJ/SpatialMetsDB/RawData/STsamples2.txt) \
>> /data_d/WSJ/SpatialMetsDB/RawData/GSE226562/clinical.txt

# GSE227427
mkdir -p /data_d/WSJ/SpatialMetsDB/RawData/GSE227427
cd /data_d/WSJ/SpatialMetsDB/RawData/GSE227427
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227427/suppl/GSE227427_raw_data.csv.gz
gunzip GSE227427_raw_data.csv.gz
cd /data_d/WSJ/SpatialMetsDB/RawData/GSE227427
cp GSE227427_raw_data.csv Expression.txt
sed -i 's/,/\t/g' Expression.txt
echo -e "sample\tgroup\tslide" > clinical.txt
paste \
<(head -n1 Expression.txt | \
awk -F"\t" '{for(i=2;i<=NF;i++) print $i}') \
<(awk -F"\t" '
$7=="GSE227427"{
    print $11"\t"$6
}' /data_d/WSJ/SpatialMetsDB/RawData/STsamples2.txt) \
>> clinical.txt

# GSE231806
mkdir -p /data_d/WSJ/SpatialMetsDB/RawData/GSE231806
wget -c \
-O /data_d/WSJ/SpatialMetsDB/RawData/GSE231806/GSE231806_Data.tsv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE231806&format=file&file=GSE231806%5FData%2Etsv%2Egz"
gunzip -f /data_d/WSJ/SpatialMetsDB/RawData/GSE231806/GSE231806_Data.tsv.gz
rawsample="GSE231806"
originalpath=`pwd`
cd ../RawData/$rawsample/
Rscript -e '
  inputfile="Expression.txt"
  outputfile="Expression.1.txt"
  data <- read.table(inputfile,sep = "\t",header = TRUE,check.names = FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.table(data,file = outputfile,sep = "\t",quote = FALSE)
'
mv Expression.1.txt Expression.txt
cd $originalpath
head -n 1 ../RawData/$rawsample/Expression.txt | awk -F"\t" '{print "Gene\t"$0}' > ../RawData/$rawsample/Expression.1.txt
tail -n +2 ../RawData/$rawsample/Expression.txt >> ../RawData/$rawsample/Expression.1.txt
mv ../RawData/$rawsample/Expression.1.txt ../RawData/$rawsample/Expression.txt
awk -F'\t' '
BEGIN{
    OFS="\t"

    n=0
    while((getline < "../STsamples2.txt")>0){
        if($7=="GSE231806"){
            n++
            group[n]=$11
            slide[n]=$6
        }
    }

    print "sample\tgroup\tslide"
}

NR==1{
    for(i=12;i<=NF;i++){
        idx=i-11
        print $i"\t"group[idx]"\t"slide[idx]
    }
}
' GSE231806_Data.tsv > clinical.txt

# GSE235475
mkdir -p /data_d/WSJ/SpatialMetsDB/RawData/GSE235475 && \
cd /data_d/WSJ/SpatialMetsDB/RawData/GSE235475 && \
wget -c \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235475&format=file&file=GSE235475_raw_data_file.txt.gz" \
-O GSE235475_raw_data_file.txt.gz
gunzip GSE235475_raw_data_file.txt.gz
echo -e "sample\tgroup\tslide" > clinical.txt

awk -F'\t' '
$7=="GSE235475" {
    gsm=$10
    group=$11
    slide=$6

    gsub(/\r/,"",gsm)
    gsub(/\r/,"",group)
    gsub(/\r/,"",slide)

    print gsm "\t" group "\t" slide
}
' /data_d/WSJ/SpatialMetsDB/RawData/STsamples2.txt >> clinical.txt
cp /data_d/WSJ/SpatialMetsDB/RawData/GSE235475/GSE235475_raw_data_file.txt \
   /data_d/WSJ/SpatialMetsDB/RawData/GSE235475/Expression.txt