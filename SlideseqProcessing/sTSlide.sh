#!/bin/bash
set -euo pipefail

#GSE200278
mkdir ../RawData/GSE200278/
wget -O "../RawData/GSE200278/GSE200278_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200278&format=file"
tar -xvf "../RawData/GSE200278/GSE200278_RAW.tar" -C ../RawData/GSE200278/
gunzip ../RawData/GSE200278/*gz
sed -i 's/"//g' "../RawData/GSE200278/"*
ls ../RawData/GSE200278/*slide_raw_counts.csv | while read file
do
  rawsample=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]}' `
  mkdir ../RawData/$rawsample
  cp $file ../RawData/$rawsample/expression.csv
  echo -e ",X,Y" > ../RawData/$rawsample/metadata.csv
  cut -d"," -f1,3,4 ${file/"slide_raw_counts.csv"/"slide_spatial_info.csv"} | tail -n +2 >> ../RawData/$rawsample/metadata.csv
done

#GSE216055
mkdir ../RawData/GSE216055/
wget -O "../RawData/GSE216055/GSE216055_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE216055&format=file"
tar -xvf "../RawData/GSE216055/GSE216055_RAW.tar" -C ../RawData/GSE216055/
gunzip ../RawData/GSE216055/*gz
sed -i 's/"//g' "../RawData/GSE216055/"*
ls ../RawData/GSE216055/*raw_counts.csv | while read file
do
  rawsample=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]}' `
  mkdir ../RawData/$rawsample
  cp $file ../RawData/$rawsample/expression.csv
  echo -e ",X,Y" > ../RawData/$rawsample/metadata.csv
  cut -d"," -f2,3,4 ${file/"raw_counts.csv"/"spatial_coordinates.csv"} | tail -n +2 >> ../RawData/$rawsample/metadata.csv
done

#HTAPP
mkdir ../RawData/HTAPP-997-SMP-7789
wget -O "../RawData/HTAPP-997-SMP-7789/adata.h5ad" https://datasets.cellxgene.cziscience.com/55185bb7-deba-497f-bb8e-a2fe098bdcc7.h5ad
mkdir ../RawData/HTAPP-982-SMP-7629
wget -O "../RawData/HTAPP-982-SMP-7629/adata.h5ad" https://datasets.cellxgene.cziscience.com/e8ef7247-18e3-4d62-8357-bf446d3731cb.h5ad
mkdir ../RawData/HTAPP-944-SMP-7479
wget -O "../RawData/HTAPP-944-SMP-7479/adata.h5ad" https://datasets.cellxgene.cziscience.com/bffd3499-b106-48aa-a935-2a562f68e012.h5ad
mkdir ../RawData/HTAPP-917-SMP-4531
wget -O "../RawData/HTAPP-917-SMP-4531/adata.h5ad" https://datasets.cellxgene.cziscience.com/cecc7d77-0577-4600-b848-80742c472c3f.h5ad
mkdir ../RawData/HTAPP-895-SMP-7359
wget -O "../RawData/HTAPP-895-SMP-7359/adata.h5ad" https://datasets.cellxgene.cziscience.com/0c880b0d-ef83-411d-b6ca-7b2f543f47bb.h5ad
mkdir ../RawData/HTAPP-880-SMP-7179
wget -O "../RawData/HTAPP-880-SMP-7179/adata.h5ad" https://datasets.cellxgene.cziscience.com/f3b4d6ca-e7e6-43a4-8acb-272c210ac13a.h5ad
mkdir ../RawData/HTAPP-878-SMP-7149
wget -O "../RawData/HTAPP-878-SMP-7149/adata.h5ad" https://datasets.cellxgene.cziscience.com/ab8551c5-cafe-4a07-bf00-097167a87831.h5ad
mkdir ../RawData/HTAPP-853-SMP-4381
wget -O "../RawData/HTAPP-853-SMP-4381/adata.h5ad" https://datasets.cellxgene.cziscience.com/1130b63b-68f6-4d04-96ed-75cb7fabc1ee.h5ad
mkdir ../RawData/HTAPP-812-SMP-8239
wget -O "../RawData/HTAPP-812-SMP-8239/adata.h5ad" https://datasets.cellxgene.cziscience.com/f9743bd7-ac48-4fe5-ba47-516b4ecd1a58.h5ad
mkdir ../RawData/HTAPP-783-SMP-4081
wget -O "../RawData/HTAPP-783-SMP-4081/adata.h5ad" https://datasets.cellxgene.cziscience.com/0e4ba413-8673-4ac7-80e9-e977557753e7.h5ad
mkdir ../RawData/HTAPP-514-SMP-6760
wget -O "../RawData/HTAPP-514-SMP-6760/adata.h5ad" https://datasets.cellxgene.cziscience.com/ccc7ab02-1721-442e-858d-a77f1ce91c39.h5ad
mkdir ../RawData/HTAPP-364-SMP-1321
wget -O "../RawData/HTAPP-364-SMP-1321/adata.h5ad" https://datasets.cellxgene.cziscience.com/0dc10486-eca5-48af-8301-f61d9112c848.h5ad
mkdir ../RawData/HTAPP-330-SMP-1082
wget -O "../RawData/HTAPP-330-SMP-1082/adata.h5ad" https://datasets.cellxgene.cziscience.com/4b2fe7f7-e1a8-441b-a6e6-1693452f714b.h5ad
mkdir ../RawData/HTAPP-313-SMP-932
wget -O "../RawData/HTAPP-313-SMP-932/adata.h5ad" https://datasets.cellxgene.cziscience.com/b50ba134-9dac-49ee-974f-fc7c933a7acc.h5ad
mkdir ../RawData/HTAPP-213-SMP-6752
wget -O "../RawData/HTAPP-213-SMP-6752/adata.h5ad" https://datasets.cellxgene.cziscience.com/711fb5e6-32ec-4e2c-8fa7-004922ef6d9a.h5ad


ls ../RawData/HTAPP-*-SMP-*/adata.h5ad  |while read file
do
originalpath=`pwd`
cd ${file/"adata.h5ad"/}
python <<EOF
import scanpy as sc
import pandas as pd
adata = sc.read_h5ad("adata.h5ad")
expr = pd.DataFrame(adata.X.toarray(),index=adata.obs_names,columns=adata.var_names)
expr = expr.T
expr.to_csv("expression.csv")
coords = adata.obs[['x_orig','y_orig']]
coords.columns = ["X","Y"]
coords.to_csv("metadata.csv")
EOF
cd $originalpath
done

awk -F'\t' '$3=="gene"{match($9,/gene_id "([^"]+)"/,a);match($9,/gene_name "([^"]+)"/,b);if(a[1] && b[1]){sub(/\.[0-9]+$/,"",a[1]);print a[1]","b[1]}}' /data_d/WSJ/STADstudy/Reference/gencode.v32.annotation.gtf > /data_d/WSJ/STADstudy/Reference/ensg2symbol.csv
ls ../RawData/HTAPP-*-SMP-*/expression.csv | while read file
do
head -n 1 $file > $file.1
awk 'BEGIN{FS=OFS=","}\
     ARGIND==1{Gene[$1]=$2}\
     ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
     /data_d/WSJ/STADstudy/Reference/ensg2symbol.csv\
     $file >> $file.1
mv $file.1 $file 
originalpath=`pwd`
cd ${file/"expression.csv"/}
Rscript -e '
  inputfile="expression.csv"
  outputfile="expression.1.csv"
  data <- read.csv(inputfile, check.names=FALSE)
  rownames(data) <- make.unique(data[,1])
  data[,1] <- NULL
  write.csv(data,outputfile)
'
mv expression.1.csv expression.csv
cd $originalpath
done

#GSE223501
mkdir "../RawData/GSE223501/"
wget -O "../RawData/GSE223501/GSE223501_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223501&format=file"
tar -xvf "../RawData/GSE223501/GSE223501_RAW.tar" -C ../RawData/GSE223501/
gunzip ../RawData/GSE223501/*gz
sed -i 's/"//g' "../RawData/GSE223501/"*
ls ../RawData/GSE223501/*counts.csv| while read file
do
  rawsample=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]}' `
  mkdir ../RawData/$rawsample
  head -n 1 $file | awk '{gsub(/\./,"-");print $0}' > ../RawData/$rawsample/expression.csv
  tail -n +2 $file  >> ../RawData/$rawsample/expression.csv
  echo -e ",X,Y" > ../RawData/$rawsample/metadata.csv
  tail -n +2 ${file/"counts.csv"/"coords.csv"}  >> ../RawData/$rawsample/metadata.csv
done

#GSE260797
mkdir "../RawData/GSE260797/"
wget -O "../RawData/GSE260797/GSE260797_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE260797&format=file"
tar -xvf "../RawData/GSE260797/GSE260797_RAW.tar" -C ../RawData/GSE260797/
gunzip ../RawData/GSE260797/*gz
ls ../RawData/GSE260797/*.digital_expression.txt | while read file
do
  rawsample=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]}' `
  mkdir ../RawData/$rawsample
  tail -n +4 $file > ../RawData/$rawsample/expression.csv
  head -n 1 ../RawData/$rawsample/expression.csv > ../RawData/$rawsample/expression.1.csv
  awk 'BEGIN{FS=OFS="\t"}\
     ARGIND==1 && ($1!="") && ($2!=""){Gene[$1]=$2}\
     ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
     biomart_humna_mouse.txt\
     ../RawData/$rawsample/expression.csv >> ../RawData/$rawsample/expression.1.csv
  mv ../RawData/$rawsample/expression.1.csv ../RawData/$rawsample/expression.csv
  sed -i 's/\t/,/g' ../RawData/$rawsample/expression.csv
  
  originalpath=`pwd`
  cd ../RawData/$rawsample/
  Rscript -e '
    inputfile="expression.csv"
    outputfile="expression.1.csv"
    data <- read.csv(inputfile, check.names=FALSE)
    rownames(data) <- make.unique(data[,1])
    data[,1] <- NULL
    write.csv(data,outputfile)
  '
  mv expression.1.csv expression.csv
  cd $originalpath
  
  echo -e ",X,Y" > ../RawData/$rawsample/metadata.csv
  head -n 1 ../RawData/$rawsample/expression.csv | tr ',' '\n' | tail -n +2 > ../RawData/$rawsample/metadata.1.csv
  awk -F"\t" '{print $2","$3}' ${file/".digital_expression.txt"/"_matched_bead_locations.txt"}  > ../RawData/$rawsample/metadata.2.csv
  paste -d $',' ../RawData/$rawsample/metadata.1.csv ../RawData/$rawsample/metadata.2.csv >> ../RawData/$rawsample/metadata.csv
done

