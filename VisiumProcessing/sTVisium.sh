#!/bin/bash
set -euo pipefail

cd /data_d/WSJ/SpatialMetsDB/RawData
 mkdir -p zenodo14204217
cd zenodo14204217
wget -c https://zenodo.org/records/14204217/files/byArray.tar
tar -xvf byArray.tar 
wget -c https://zenodo.org/records/14204217/files/Clinical.tar
tar -xvf Clinical.tar
wget -c https://zenodo.org/records/14204217/files/rawCountsMatrices.tar
tar -xvf rawCountsMatrices.tar
 Rscript convert_zenodo_to_10x.R   --raw_dir /data_d/WSJ/SpatialMetsDB/RawData/zenodo14204217   --out_dir /data_d/WSJ/SpatialMetsDB/RawData   --ncores 4   --overwrite
for dir in zenodo14204217_*; do
    num=${dir##*_}
    mv "$dir" "$num"
done
Rscript /data_d/WSJ/SpatialMetsDB/Script/Batch_ENSG2Symbol.R

#GSM5420753
mkdir -p GSM5420753
cd GSM5420753
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM6177599
mkdir -p GSM6177599
cd GSM6177599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6177601
mkdir -p GSM6177601
cd GSM6177601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6177603
mkdir -p GSM6177603
cd GSM6177603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6592048
mkdir -p GSM6592048
cd GSM6592048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592049
mkdir -p GSM6592049
cd GSM6592049
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592050
mkdir -p GSM6592050
cd GSM6592050
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592051
mkdir -p GSM6592051
cd GSM6592051
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592052
mkdir -p GSM6592052
cd GSM6592052
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592053
mkdir -p GSM6592053
cd GSM6592053
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592054
mkdir -p GSM6592054
cd GSM6592054
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592055
mkdir -p GSM6592055
cd GSM6592055
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592056
mkdir -p GSM6592056
cd GSM6592056
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592057
mkdir -p GSM6592057
cd GSM6592057
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5“
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Fscalefactors%5Fjson%2Ejson%2Egz”
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz“
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz”
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592058
mkdir -p GSM6592058
cd GSM6592058
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592059
mkdir -p GSM6592059
cd GSM6592059
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592060
mkdir -p GSM6592060
cd GSM6592060
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592061
mkdir -p GSM6592061
cd GSM6592061
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592062
mkdir -p GSM6592062
cd GSM6592062
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7757971
mkdir -p GSM7757971
mkdir -p GSM7757971/filtered_feature_bc_matrix
cd GSM7757971/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757971/filtered_feature_bc_matrix --output_file GSM7757971/filtered_feature_bc_matrix.h5
mkdir -p GSM7757971/spatial
cd GSM7757971/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757972
mkdir -p GSM7757972
mkdir -p GSM7757972/filtered_feature_bc_matrix
cd GSM7757972/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757972/filtered_feature_bc_matrix --output_file GSM7757972/filtered_feature_bc_matrix.h5
mkdir -p GSM7757972/spatial
cd GSM7757972/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757973
mkdir -p GSM7757973
mkdir -p GSM7757973/filtered_feature_bc_matrix
cd GSM7757973/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757973/filtered_feature_bc_matrix --output_file GSM7757973/filtered_feature_bc_matrix.h5
mkdir -p GSM7757973/spatial
cd GSM7757973/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757974
mkdir -p GSM7757974
mkdir -p GSM7757974/filtered_feature_bc_matrix
cd GSM7757974/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757974/filtered_feature_bc_matrix --output_file GSM7757974/filtered_feature_bc_matrix.h5
mkdir -p GSM7757974/spatial
cd GSM7757974/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757975
mkdir -p GSM7757975
mkdir -p GSM7757975/filtered_feature_bc_matrix
cd GSM7757975/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757975/filtered_feature_bc_matrix --output_file GSM7757975/filtered_feature_bc_matrix.h5
mkdir -p GSM7757975/spatial
cd GSM7757975/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757976
mkdir -p GSM7757976
mkdir -p GSM7757976/filtered_feature_bc_matrix
cd GSM7757976/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757976/filtered_feature_bc_matrix --output_file GSM7757976/filtered_feature_bc_matrix.h5
mkdir -p GSM7757976/spatial
cd GSM7757976/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757977
mkdir -p GSM7757977
mkdir -p GSM7757977/filtered_feature_bc_matrix
cd GSM7757977/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757977/filtered_feature_bc_matrix --output_file GSM7757977/filtered_feature_bc_matrix.h5
mkdir -p GSM7757977/spatial
cd GSM7757977/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757978
mkdir -p GSM7757978
mkdir -p GSM7757978/filtered_feature_bc_matrix
cd GSM7757978/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757978/filtered_feature_bc_matrix --output_file GSM7757978/filtered_feature_bc_matrix.h5
mkdir -p GSM7757978/spatial
cd GSM7757978/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757979
mkdir -p GSM7757979
mkdir -p GSM7757979/filtered_feature_bc_matrix
cd GSM7757979/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757979/filtered_feature_bc_matrix --output_file GSM7757979/filtered_feature_bc_matrix.h5
mkdir -p GSM7757979/spatial
cd GSM7757979/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757980
mkdir -p GSM7757980
mkdir -p GSM7757980/filtered_feature_bc_matrix
cd GSM7757980/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757980/filtered_feature_bc_matrix --output_file GSM7757980/filtered_feature_bc_matrix.h5
mkdir -p GSM7757980/spatial
cd GSM7757980/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757981
mkdir -p GSM7757981
mkdir -p GSM7757981/filtered_feature_bc_matrix
cd GSM7757981/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757981/filtered_feature_bc_matrix --output_file GSM7757981/filtered_feature_bc_matrix.h5
mkdir -p GSM7757981/spatial
cd GSM7757981/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757982
mkdir -p GSM7757982
mkdir -p GSM7757982/filtered_feature_bc_matrix
cd GSM7757982/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757982/filtered_feature_bc_matrix --output_file GSM7757982/filtered_feature_bc_matrix.h5
mkdir -p GSM7757982/spatial
cd GSM7757982/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757983
mkdir -p GSM7757983
mkdir -p GSM7757983/filtered_feature_bc_matrix
cd GSM7757983/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757983/filtered_feature_bc_matrix --output_file GSM7757983/filtered_feature_bc_matrix.h5
mkdir -p GSM7757983/spatial
cd GSM7757983/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757984
mkdir -p GSM7757984
mkdir -p GSM7757984/filtered_feature_bc_matrix
cd GSM7757984/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757984/filtered_feature_bc_matrix --output_file GSM7757984/filtered_feature_bc_matrix.h5
mkdir -p GSM7757984/spatial
cd GSM7757984/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757985
mkdir -p GSM7757985
mkdir -p GSM7757985/filtered_feature_bc_matrix
cd GSM7757985/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757985/filtered_feature_bc_matrix --output_file GSM7757985/filtered_feature_bc_matrix.h5
mkdir -p GSM7757985/spatial
cd GSM7757985/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7782699
mkdir -p GSM7782699
cd GSM7782699
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
cd CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz
rm -f CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz
cd ../

#CytAssist_Fresh_Frozen_Human_Breast_Cancer
mkdir -p CytAssist_Fresh_Frozen_Human_Breast_Cancer
cd CytAssist_Fresh_Frozen_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz
rm -f CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Parent_Visium_Human_BreastCancer
mkdir -p Parent_Visium_Human_BreastCancer
cd Parent_Visium_Human_BreastCancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz"
tar -xzf Parent_Visium_Human_BreastCancer_spatial.tar.gz
rm -f Parent_Visium_Human_BreastCancer_spatial.tar.gz
cd ../

#V1_Breast_Cancer_Block_A_Section_1
mkdir -p V1_Breast_Cancer_Block_A_Section_1
cd V1_Breast_Cancer_Block_A_Section_1
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"
tar -xzf V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
rm -f V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
cd ../

#V1_Breast_Cancer_Block_A_Section_2
mkdir -p V1_Breast_Cancer_Block_A_Section_2
cd V1_Breast_Cancer_Block_A_Section_2
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz"
tar -xzf V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz
rm -f V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz
cd ../

#V1_Human_Invasive_Ductal_Carcinoma
mkdir -p V1_Human_Invasive_Ductal_Carcinoma
cd V1_Human_Invasive_Ductal_Carcinoma
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz"
tar -xzf V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz
rm -f V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz
cd ../

#Visium_FFPE_Human_Breast_Cancer
mkdir -p Visium_FFPE_Human_Breast_Cancer
cd Visium_FFPE_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz
rm -f Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Visium_Human_Breast_Cancer
mkdir -p Visium_Human_Breast_Cancer
cd Visium_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf Visium_Human_Breast_Cancer_spatial.tar.gz
rm -f Visium_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Targeted_Visium_Human_BreastCancer_Immunology
mkdir -p Targeted_Visium_Human_BreastCancer_Immunology
cd Targeted_Visium_Human_BreastCancer_Immunology
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz"
tar -xzf Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz
rm -f Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz
cd ../

#GSM7777520
mkdir -p GSM7777520
cd GSM7777520
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777520&format=file&file=GSM7777520%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777520&format=file&file=GSM7777520%5FA1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM7777521
mkdir -p GSM7777521
cd GSM7777521
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777521&format=file&file=GSM7777521%5FB1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777521&format=file&file=GSM7777521%5FB1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../


#GSM7777522
mkdir -p GSM7777522
cd GSM7777522
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777522&format=file&file=GSM7777522%5FC1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777522&format=file&file=GSM7777522%5FC1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM7777523
mkdir -p GSM7777523
cd GSM7777523
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777523&format=file&file=GSM7777523%5FD1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777523&format=file&file=GSM7777523%5FD1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM8291069
mkdir -p GSM8291069
cd GSM8291069
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
tail -n +2 tissue_positions_list_original.csv > tissue_positions_list.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM8291070
mkdir -p GSM8291070
cd GSM8291070
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
tail -n +2 tissue_positions_list_original.csv > tissue_positions_list.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM8291071
mkdir -p GSM8291071
cd GSM8291071
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
tail -n +2 tissue_positions_list_original.csv > tissue_positions_list.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
gunzip scalefactors_json.json.gz
cd ../../

#GSM5693665
mkdir -p GSM5693665
cd GSM5693665
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Efiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Etissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Etissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Escalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


#GSM6433585
mkdir -p GSM6433585
cd GSM6433585
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433586
mkdir -p GSM6433586
cd GSM6433586
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433587
mkdir -p GSM6433587
cd GSM6433587
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433588
mkdir -p GSM6433588
cd GSM6433588
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433589
mkdir -p GSM6433589
cd GSM6433589
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433590
mkdir -p GSM6433590
cd GSM6433590
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433591
mkdir -p GSM6433591
cd GSM6433591
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433592
mkdir -p GSM6433592
cd GSM6433592
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433593
mkdir -p GSM6433593
cd GSM6433593
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433594
mkdir -p GSM6433594
cd GSM6433594
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433595
mkdir -p GSM6433595
cd GSM6433595
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433596
mkdir -p GSM6433596
cd GSM6433596
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433597
cd GSM6433597
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433598
mkdir -p GSM6433598
cd GSM6433598
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433599
mkdir -p GSM6433599
cd GSM6433599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433600
mkdir -p GSM6433600
cd GSM6433600
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


#GSM6433601
mkdir -p GSM6433601
cd GSM6433601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433602
mkdir -p GSM6433602
cd GSM6433602
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433603
mkdir -p GSM6433603
cd GSM6433603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433604
mkdir -p GSM6433604
cd GSM6433604
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433605
mkdir -p GSM6433605
cd GSM6433605
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433606
mkdir -p GSM6433606
cd GSM6433606
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433607
mkdir -p GSM6433607
cd GSM6433607
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433608
mkdir -p GSM6433608
cd GSM6433608
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433609
mkdir -p GSM6433609
cd GSM6433609
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433610
mkdir -p GSM6433610
cd GSM6433610
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433611
mkdir -p GSM6433611
cd GSM6433611
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433612
mkdir -p GSM6433612
cd GSM6433612
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433613
mkdir -p GSM6433613
cd GSM6433613
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433614
mkdir -p GSM6433614
cd GSM6433614
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433615
mkdir -p GSM6433615
cd GSM6433615
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433616
mkdir -p GSM6433616
cd GSM6433616
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433617
mkdir -p GSM6433617
cd GSM6433617
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433618
mkdir -p GSM6433618
cd GSM6433618
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433619
mkdir -p GSM6433619
cd GSM6433619
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433620
mkdir -p GSM6433620
cd GSM6433620
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433621
mkdir -p GSM6433621
cd GSM6433621
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433622
mkdir -p GSM6433622
cd GSM6433622
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433623
mkdir -p GSM6433623
cd GSM6433623
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433624
mkdir -p GSM6433624
cd GSM6433624
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433625
mkdir -p GSM6433625
cd GSM6433625
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433626
mkdir -p GSM6433626
cd GSM6433626
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433627
mkdir -p GSM6433627
cd GSM6433627
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6177599
mkdir -p GSM6177599
cd GSM6177599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6177601
mkdir -p GSM6177601
cd GSM6177601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6177603
mkdir -p GSM6177603
cd GSM6177603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5420752
mkdir -p GSM5420752
cd GSM5420752
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415705
mkdir -p GSM6415705
cd GSM6415705
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415706
mkdir -p GSM6415706
cd GSM6415706
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM5924030
mkdir -p GSM5924030
cd GSM5924030
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924031
mkdir -p GSM5924031
cd GSM5924031
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924032
mkdir -p GSM5924032
cd GSM5924032
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924033
mkdir -p GSM5924033
cd GSM5924033
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924034
mkdir -p GSM5924034
cd GSM5924034
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924035
mkdir -p GSM5924035
cd GSM5924035
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924037
mkdir -p GSM5924037
cd GSM5924037
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924038
mkdir -p GSM5924038
cd GSM5924038
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924039
mkdir -p GSM5924039
cd GSM5924039
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924040
mkdir -p GSM5924040
cd GSM5924040
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924041
mkdir -p GSM5924041
cd GSM5924041
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924042
mkdir -p GSM5924042
cd GSM5924042
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924043
mkdir -p GSM5924043
cd GSM5924043
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924044
mkdir -p GSM5924044
cd GSM5924044
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924045
mkdir -p GSM5924045
cd GSM5924045
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924046
mkdir -p GSM5924046
cd GSM5924046
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924047
mkdir -p GSM5924047
cd GSM5924047
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924048
mkdir -p GSM5924048
cd GSM5924048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924049
mkdir -p GSM5924049
cd GSM5924049
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924050
mkdir -p GSM5924050
cd GSM5924050
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924051
mkdir -p GSM5924051
cd GSM5924051
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924052
mkdir -p GSM5924052
cd GSM5924052
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924053
mkdir -p GSM5924053
cd GSM5924053
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#HG_1.3
mkdir -p HG_1.3
cd HG_1.3
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_1.2
mkdir -p HG_1.2
cd HG_1.2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_3.2
mkdir -p HG_3.2
cd HG_3.2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_1
mkdir -p LG_1
cd LG_1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_2
mkdir -p LG_2
cd LG_2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_3
mkdir -p LG_3
cd LG_3
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_1.1
mkdir -p pHG_1.1
cd pHG_1.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_3.1
mkdir -p pHG_3.1
cd pHG_3.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


# GSM7974884
mkdir -p GSM7974884
cd GSM7974884
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../


# GSM7974885
mkdir -p GSM7974885
cd GSM7974885
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974886
mkdir -p GSM7974886
cd GSM7974886
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974887
mkdir -p GSM7974887
cd GSM7974887
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974888
mkdir -p GSM7974888
cd GSM7974888
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6716963
mkdir -p GSM6716963
cd GSM6716963
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716964
mkdir -p GSM6716964
cd GSM6716964
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716965
mkdir -p GSM6716965
cd GSM6716965
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716966
mkdir -p GSM6716966
cd GSM6716966
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_11mm_FFPE_Human_Colorectal_Cancer
mkdir -p CytAssist_11mm_FFPE_Human_Colorectal_Cancer
cd CytAssist_11mm_FFPE_Human_Colorectal_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Parent_Visium_Human_ColorectalCancer
mkdir -p Parent_Visium_Human_ColorectalCancer
cd Parent_Visium_Human_ColorectalCancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep1
mkdir -p CytAssist_FFPE_Human_Colon_Rep1
cd CytAssist_FFPE_Human_Colon_Rep1
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep2
mkdir -p CytAssist_FFPE_Human_Colon_Rep2
cd CytAssist_FFPE_Human_Colon_Rep2
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep2/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep2_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Targeted_Visium_Human_ColorectalCancer_GeneSignature
mkdir -p Targeted_Visium_Human_ColorectalCancer_GeneSignature
cd Targeted_Visium_Human_ColorectalCancer_GeneSignature
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM8041061 
mkdir -p GSM8041061
mkdir -p GSM8041061/filtered_feature_bc_matrix
cd GSM8041061/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041061/filtered_feature_bc_matrix --output_file GSM8041061/filtered_feature_bc_matrix.h5
mkdir -p GSM8041061/spatial
cd GSM8041061/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8041062
mkdir -p GSM8041062
mkdir -p GSM8041062/filtered_feature_bc_matrix
cd GSM8041062/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041062/filtered_feature_bc_matrix --output_file GSM8041062/filtered_feature_bc_matrix.h5
mkdir -p GSM8041062/spatial
cd GSM8041062/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6256810
mkdir -p GSM6256810
cd GSM6256810
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
# tissue positions
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256811
mkdir -p GSM6256811
cd GSM6256811
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256812
mkdir -p GSM6256812
cd GSM6256812
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256813
mkdir -p GSM6256813
cd GSM6256813
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7573205
mkdir -p GSM7573205
mkdir -p GSM7573205/filtered_feature_bc_matrix
cd GSM7573205/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM7573205/filtered_feature_bc_matrix \
--output_file GSM7573205/filtered_feature_bc_matrix.h5
mkdir -p GSM7573205/spatial
cd GSM7573205/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Escalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265211
mkdir -p GSM8265211
mkdir -p GSM8265211/filtered_feature_bc_matrix
cd GSM8265211/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265211/filtered_feature_bc_matrix \
--output_file GSM8265211/filtered_feature_bc_matrix.h5
mkdir -p GSM8265211/spatial
cd GSM8265211/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265212
mkdir -p GSM8265212
mkdir -p GSM8265212/filtered_feature_bc_matrix
cd GSM8265212/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265212/filtered_feature_bc_matrix \
--output_file GSM8265212/filtered_feature_bc_matrix.h5
mkdir -p GSM8265212/spatial
cd GSM8265212/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265213
mkdir -p GSM8265213
mkdir -p GSM8265213/filtered_feature_bc_matrix
cd GSM8265213/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265213/filtered_feature_bc_matrix \
--output_file GSM8265213/filtered_feature_bc_matrix.h5
mkdir -p GSM8265213/spatial
cd GSM8265213/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM5420750
mkdir RawData/
cd RawData/
mkdir -p GSM5420750
cd GSM5420750
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5420nnn/GSM5420750/suppl/GSM5420750%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420750&format=file&file=GSM5420750%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#CytAssist_FFPE_Human_Skin_Melanoma
mkdir -p CytAssist_FFPE_Human_Skin_Melanoma
cd CytAssist_FFPE_Human_Skin_Melanoma
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz"
tar -xzf CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz
rm -f CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz
tail -n +2 /data_d/WSJ/SpatialMetsDB/RawData/CytAssist_FFPE_Human_Skin_Melanoma/spatial/tissue_positions_list.csv > /data_d/WSJ/SpatialMetsDB/RawData/CytAssist_FFPE_Human_Skin_Melanoma/spatial/tissue_positions_list.1.csv
mv /data_d/WSJ/SpatialMetsDB/RawData/CytAssist_FFPE_Human_Skin_Melanoma/spatial/tissue_positions_list.1.csv /data_d/WSJ/SpatialMetsDB/RawData/CytAssist_FFPE_Human_Skin_Melanoma/spatial/tissue_positions_list.csv
cd ../

#GSE250636
mkdir -p GSE250636 && cd GSE250636
wget -c -O raw.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE250636&format=file"
tar -xf raw.tar
rm -f raw.tar
for h5 in *_filtered_feature_bc_matrix.h5; do
  GSM=${h5%%_*}
  mkdir -p $GSM/spatial
  cp $h5 $GSM/filtered_feature_bc_matrix.h5
  for f in scalefactors_json.json tissue_hires_image.png tissue_lowres_image.png tissue_positions_list.csv; do
    src=$(ls ${GSM}_*${f}.gz 2>/dev/null | head -n 1)
    [ -n "$src" ] && gunzip -c "$src" > $GSM/spatial/$f
  done
done

#GSM5420753
mkdir -p GSM5420753
cd GSM5420753
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM6177599
mkdir -p GSM6177599
cd GSM6177599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6177601
mkdir -p GSM6177601
cd GSM6177601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6177603
mkdir -p GSM6177603
cd GSM6177603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O aligned_fiducials.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Faligned%5Ffiducials%2Ejpg%2Egz"
gunzip aligned_fiducials.jpg.gz
wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip detected_tissue_image.jpg.gz
cd ../../

#GSM6592048
mkdir -p GSM6592048
cd GSM6592048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592049
mkdir -p GSM6592049
cd GSM6592049
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592049&format=file&file=GSM6592049%5FM2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592050
mkdir -p GSM6592050
cd GSM6592050
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592050&format=file&file=GSM6592050%5FM3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592051
mkdir -p GSM6592051
cd GSM6592051
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592051&format=file&file=GSM6592051%5FM4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592052
mkdir -p GSM6592052
cd GSM6592052
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592052&format=file&file=GSM6592052%5FM5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592053
mkdir -p GSM6592053
cd GSM6592053
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592053&format=file&file=GSM6592053%5FM6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592054
mkdir -p GSM6592054
cd GSM6592054
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592054&format=file&file=GSM6592054%5FM7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592055
mkdir -p GSM6592055
cd GSM6592055
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592055&format=file&file=GSM6592055%5FM8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592056
mkdir -p GSM6592056
cd GSM6592056
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592056&format=file&file=GSM6592056%5FM9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592057
mkdir -p GSM6592057
cd GSM6592057
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5鈥?mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Fscalefactors%5Fjson%2Ejson%2Egz鈥?gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz鈥?gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592057&format=file&file=GSM6592057%5FM10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz鈥?gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592058
mkdir -p GSM6592058
cd GSM6592058
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592058&format=file&file=GSM6592058%5FM11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592059
mkdir -p GSM6592059
cd GSM6592059
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592059&format=file&file=GSM6592059%5FM13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592060
mkdir -p GSM6592060
cd GSM6592060
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592060&format=file&file=GSM6592060%5FM14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592061
mkdir -p GSM6592061
cd GSM6592061
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592061&format=file&file=GSM6592061%5FM15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6592062
mkdir -p GSM6592062
cd GSM6592062
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592062&format=file&file=GSM6592062%5FM16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7757971
mkdir -p GSM7757971
mkdir -p GSM7757971/filtered_feature_bc_matrix
cd GSM7757971/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757971/filtered_feature_bc_matrix --output_file GSM7757971/filtered_feature_bc_matrix.h5
mkdir -p GSM7757971/spatial
cd GSM7757971/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757971&format=file&file=GSM7757971%5FR1%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757972
mkdir -p GSM7757972
mkdir -p GSM7757972/filtered_feature_bc_matrix
cd GSM7757972/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757972/filtered_feature_bc_matrix --output_file GSM7757972/filtered_feature_bc_matrix.h5
mkdir -p GSM7757972/spatial
cd GSM7757972/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757972&format=file&file=GSM7757972%5FR1%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757973
mkdir -p GSM7757973
mkdir -p GSM7757973/filtered_feature_bc_matrix
cd GSM7757973/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757973/filtered_feature_bc_matrix --output_file GSM7757973/filtered_feature_bc_matrix.h5
mkdir -p GSM7757973/spatial
cd GSM7757973/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757973&format=file&file=GSM7757973%5FR1%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757974
mkdir -p GSM7757974
mkdir -p GSM7757974/filtered_feature_bc_matrix
cd GSM7757974/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757974/filtered_feature_bc_matrix --output_file GSM7757974/filtered_feature_bc_matrix.h5
mkdir -p GSM7757974/spatial
cd GSM7757974/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757974&format=file&file=GSM7757974%5FR1%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757975
mkdir -p GSM7757975
mkdir -p GSM7757975/filtered_feature_bc_matrix
cd GSM7757975/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757975/filtered_feature_bc_matrix --output_file GSM7757975/filtered_feature_bc_matrix.h5
mkdir -p GSM7757975/spatial
cd GSM7757975/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757975&format=file&file=GSM7757975%5FR1%2DS6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757976
mkdir -p GSM7757976
mkdir -p GSM7757976/filtered_feature_bc_matrix
cd GSM7757976/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757976/filtered_feature_bc_matrix --output_file GSM7757976/filtered_feature_bc_matrix.h5
mkdir -p GSM7757976/spatial
cd GSM7757976/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757976&format=file&file=GSM7757976%5FR2%2DS1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757977
mkdir -p GSM7757977
mkdir -p GSM7757977/filtered_feature_bc_matrix
cd GSM7757977/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757977/filtered_feature_bc_matrix --output_file GSM7757977/filtered_feature_bc_matrix.h5
mkdir -p GSM7757977/spatial
cd GSM7757977/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757977&format=file&file=GSM7757977%5FR2%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757978
mkdir -p GSM7757978
mkdir -p GSM7757978/filtered_feature_bc_matrix
cd GSM7757978/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757978/filtered_feature_bc_matrix --output_file GSM7757978/filtered_feature_bc_matrix.h5
mkdir -p GSM7757978/spatial
cd GSM7757978/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757978&format=file&file=GSM7757978%5FR2%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757979
mkdir -p GSM7757979
mkdir -p GSM7757979/filtered_feature_bc_matrix
cd GSM7757979/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757979/filtered_feature_bc_matrix --output_file GSM7757979/filtered_feature_bc_matrix.h5
mkdir -p GSM7757979/spatial
cd GSM7757979/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757979&format=file&file=GSM7757979%5FR2%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757980
mkdir -p GSM7757980
mkdir -p GSM7757980/filtered_feature_bc_matrix
cd GSM7757980/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757980/filtered_feature_bc_matrix --output_file GSM7757980/filtered_feature_bc_matrix.h5
mkdir -p GSM7757980/spatial
cd GSM7757980/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757980&format=file&file=GSM7757980%5FR2%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757981
mkdir -p GSM7757981
mkdir -p GSM7757981/filtered_feature_bc_matrix
cd GSM7757981/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757981/filtered_feature_bc_matrix --output_file GSM7757981/filtered_feature_bc_matrix.h5
mkdir -p GSM7757981/spatial
cd GSM7757981/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757981&format=file&file=GSM7757981%5FR3%2DS1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757982
mkdir -p GSM7757982
mkdir -p GSM7757982/filtered_feature_bc_matrix
cd GSM7757982/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757982/filtered_feature_bc_matrix --output_file GSM7757982/filtered_feature_bc_matrix.h5
mkdir -p GSM7757982/spatial
cd GSM7757982/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757982&format=file&file=GSM7757982%5FR3%2DS2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757983
mkdir -p GSM7757983
mkdir -p GSM7757983/filtered_feature_bc_matrix
cd GSM7757983/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757983/filtered_feature_bc_matrix --output_file GSM7757983/filtered_feature_bc_matrix.h5
mkdir -p GSM7757983/spatial
cd GSM7757983/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757983&format=file&file=GSM7757983%5FR3%2DS3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757984
mkdir -p GSM7757984
mkdir -p GSM7757984/filtered_feature_bc_matrix
cd GSM7757984/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757984/filtered_feature_bc_matrix --output_file GSM7757984/filtered_feature_bc_matrix.h5
mkdir -p GSM7757984/spatial
cd GSM7757984/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757984&format=file&file=GSM7757984%5FR3%2DS4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7757985
mkdir -p GSM7757985
mkdir -p GSM7757985/filtered_feature_bc_matrix
cd GSM7757985/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7757985/filtered_feature_bc_matrix --output_file GSM7757985/filtered_feature_bc_matrix.h5
mkdir -p GSM7757985/spatial
cd GSM7757985/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757985&format=file&file=GSM7757985%5FR3%2DS5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM7782699
mkdir -p GSM7782699
cd GSM7782699
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
cd CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz
rm -f CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz
cd ../

#CytAssist_Fresh_Frozen_Human_Breast_Cancer
mkdir -p CytAssist_Fresh_Frozen_Human_Breast_Cancer
cd CytAssist_Fresh_Frozen_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz
rm -f CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Parent_Visium_Human_BreastCancer
mkdir -p Parent_Visium_Human_BreastCancer
cd Parent_Visium_Human_BreastCancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz"
tar -xzf Parent_Visium_Human_BreastCancer_spatial.tar.gz
rm -f Parent_Visium_Human_BreastCancer_spatial.tar.gz
cd ../

#V1_Breast_Cancer_Block_A_Section_1
mkdir -p V1_Breast_Cancer_Block_A_Section_1
cd V1_Breast_Cancer_Block_A_Section_1
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"
tar -xzf V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
rm -f V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
cd ../

#V1_Breast_Cancer_Block_A_Section_2
mkdir -p V1_Breast_Cancer_Block_A_Section_2
cd V1_Breast_Cancer_Block_A_Section_2
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz"
tar -xzf V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz
rm -f V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz
cd ../

#V1_Human_Invasive_Ductal_Carcinoma
mkdir -p V1_Human_Invasive_Ductal_Carcinoma
cd V1_Human_Invasive_Ductal_Carcinoma
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz"
tar -xzf V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz
rm -f V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz
cd ../

#Visium_FFPE_Human_Breast_Cancer
mkdir -p Visium_FFPE_Human_Breast_Cancer
cd Visium_FFPE_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz
rm -f Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Visium_Human_Breast_Cancer
mkdir -p Visium_Human_Breast_Cancer
cd Visium_Human_Breast_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf Visium_Human_Breast_Cancer_spatial.tar.gz
rm -f Visium_Human_Breast_Cancer_spatial.tar.gz
cd ../

#Targeted_Visium_Human_BreastCancer_Immunology
mkdir -p Targeted_Visium_Human_BreastCancer_Immunology
cd Targeted_Visium_Human_BreastCancer_Immunology
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_filtered_feature_bc_matrix.h5"
wget "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz"
tar -xzf Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz
rm -f Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz
cd ../

#GSM7777520
mkdir -p GSM7777520
cd GSM7777520
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777520&format=file&file=GSM7777520%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777520&format=file&file=GSM7777520%5FA1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM7777521
mkdir -p GSM7777521
cd GSM7777521
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777521&format=file&file=GSM7777521%5FB1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777521&format=file&file=GSM7777521%5FB1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../


#GSM7777522
mkdir -p GSM7777522
cd GSM7777522
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777522&format=file&file=GSM7777522%5FC1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777522&format=file&file=GSM7777522%5FC1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM7777523
mkdir -p GSM7777523
cd GSM7777523
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777523&format=file&file=GSM7777523%5FD1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7777523&format=file&file=GSM7777523%5FD1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

#GSM8291069
mkdir -p GSM8291069
cd GSM8291069
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291069&format=file&file=GSM8291069%5FV42Y04%2D372%5FA1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM8291070
mkdir -p GSM8291070
cd GSM8291070
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291070&format=file&file=GSM8291070%5FV42Y04%2D403%5FA1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM8291071
mkdir -p GSM8291071
cd GSM8291071
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
mv tissue_positions_list.csv tissue_positions_list_original.csv
tail -n +2 tissue_positions_list_original.csv > tissue_positions_list.csv
rm -f tissue_positions_list_original.csv
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8291071&format=file&file=GSM8291071%5FV42Y09%2D344%5FA1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
gunzip scalefactors_json.json.gz
cd ../../

#GSM5693665
mkdir -p GSM5693665
cd GSM5693665
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Efiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Etissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Etissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5693665&format=file&file=GSM5693665%5FG2021%2DST%2Escalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


#GSM6433585
mkdir -p GSM6433585
cd GSM6433585
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433585&format=file&file=GSM6433585%5F092A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433586
mkdir -p GSM6433586
cd GSM6433586
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433586&format=file&file=GSM6433586%5F092B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433587
mkdir -p GSM6433587
cd GSM6433587
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433587&format=file&file=GSM6433587%5F093A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433588
mkdir -p GSM6433588
cd GSM6433588
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433588&format=file&file=GSM6433588%5F093B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433589
mkdir -p GSM6433589
cd GSM6433589
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433589&format=file&file=GSM6433589%5F093C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433590
mkdir -p GSM6433590
cd GSM6433590
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433590&format=file&file=GSM6433590%5F093D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433591
mkdir -p GSM6433591
cd GSM6433591
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433591&format=file&file=GSM6433591%5F094A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433592
mkdir -p GSM6433592
cd GSM6433592
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433592&format=file&file=GSM6433592%5F094B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433593
mkdir -p GSM6433593
cd GSM6433593
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433593&format=file&file=GSM6433593%5F094C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433594
mkdir -p GSM6433594
cd GSM6433594
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433594&format=file&file=GSM6433594%5F094D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433595
mkdir -p GSM6433595
cd GSM6433595
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433595&format=file&file=GSM6433595%5F095A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433596
mkdir -p GSM6433596
cd GSM6433596
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433596&format=file&file=GSM6433596%5F095B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433597
cd GSM6433597
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433597&format=file&file=GSM6433597%5F117B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433598
mkdir -p GSM6433598
cd GSM6433598
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433598&format=file&file=GSM6433598%5F117C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433599
mkdir -p GSM6433599
cd GSM6433599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433599&format=file&file=GSM6433599%5F117D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433600
mkdir -p GSM6433600
cd GSM6433600
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433600&format=file&file=GSM6433600%5F117E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


#GSM6433601
mkdir -p GSM6433601
cd GSM6433601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433601&format=file&file=GSM6433601%5F118B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433602
mkdir -p GSM6433602
cd GSM6433602
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433602&format=file&file=GSM6433602%5F118C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433603
mkdir -p GSM6433603
cd GSM6433603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433603&format=file&file=GSM6433603%5F118D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433604
mkdir -p GSM6433604
cd GSM6433604
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433604&format=file&file=GSM6433604%5F118E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433605
mkdir -p GSM6433605
cd GSM6433605
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433605&format=file&file=GSM6433605%5F119B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433606
mkdir -p GSM6433606
cd GSM6433606
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433606&format=file&file=GSM6433606%5F119C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433607
mkdir -p GSM6433607
cd GSM6433607
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433607&format=file&file=GSM6433607%5F119D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433608
mkdir -p GSM6433608
cd GSM6433608
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433608&format=file&file=GSM6433608%5F119E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433609
mkdir -p GSM6433609
cd GSM6433609
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433609&format=file&file=GSM6433609%5F120B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433610
mkdir -p GSM6433610
cd GSM6433610
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433610&format=file&file=GSM6433610%5F120C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433611
mkdir -p GSM6433611
cd GSM6433611
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433611&format=file&file=GSM6433611%5F120D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433612
mkdir -p GSM6433612
cd GSM6433612
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433612&format=file&file=GSM6433612%5F120E%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433613
mkdir -p GSM6433613
cd GSM6433613
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433613&format=file&file=GSM6433613%5F395A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433614
mkdir -p GSM6433614
cd GSM6433614
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433614&format=file&file=GSM6433614%5F395B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433615
mkdir -p GSM6433615
cd GSM6433615
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433615&format=file&file=GSM6433615%5F395C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433616
mkdir -p GSM6433616
cd GSM6433616
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433616&format=file&file=GSM6433616%5F395D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433617
mkdir -p GSM6433617
cd GSM6433617
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433617&format=file&file=GSM6433617%5F396A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433618
mkdir -p GSM6433618
cd GSM6433618
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433618&format=file&file=GSM6433618%5F396C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433619
mkdir -p GSM6433619
cd GSM6433619
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433619&format=file&file=GSM6433619%5F396D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433620
mkdir -p GSM6433620
cd GSM6433620
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433620&format=file&file=GSM6433620%5F397A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433621
mkdir -p GSM6433621
cd GSM6433621
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433621&format=file&file=GSM6433621%5F397B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433622
mkdir -p GSM6433622
cd GSM6433622
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433622&format=file&file=GSM6433622%5F397C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433623
mkdir -p GSM6433623
cd GSM6433623
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433623&format=file&file=GSM6433623%5F397D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433624
mkdir -p GSM6433624
cd GSM6433624
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433624&format=file&file=GSM6433624%5F398A%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433625
mkdir -p GSM6433625
cd GSM6433625
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433625&format=file&file=GSM6433625%5F398B%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433626
mkdir -p GSM6433626
cd GSM6433626
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433626&format=file&file=GSM6433626%5F398C%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6433627
mkdir -p GSM6433627
cd GSM6433627
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz 
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6433627&format=file&file=GSM6433627%5F398D%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6177599
mkdir -p GSM6177599
cd GSM6177599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6177601
mkdir -p GSM6177601
cd GSM6177601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6177603
mkdir -p GSM6177603
cd GSM6177603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5420752
mkdir -p GSM5420752
cd GSM5420752
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415705
mkdir -p GSM6415705
cd GSM6415705
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415706
mkdir -p GSM6415706
cd GSM6415706
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM5924030
mkdir -p GSM5924030
cd GSM5924030
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924031
mkdir -p GSM5924031
cd GSM5924031
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924032
mkdir -p GSM5924032
cd GSM5924032
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924033
mkdir -p GSM5924033
cd GSM5924033
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924034
mkdir -p GSM5924034
cd GSM5924034
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924035
mkdir -p GSM5924035
cd GSM5924035
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924037
mkdir -p GSM5924037
cd GSM5924037
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924038
mkdir -p GSM5924038
cd GSM5924038
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924039
mkdir -p GSM5924039
cd GSM5924039
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924040
mkdir -p GSM5924040
cd GSM5924040
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924041
mkdir -p GSM5924041
cd GSM5924041
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924042
mkdir -p GSM5924042
cd GSM5924042
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924043
mkdir -p GSM5924043
cd GSM5924043
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924044
mkdir -p GSM5924044
cd GSM5924044
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924045
mkdir -p GSM5924045
cd GSM5924045
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924046
mkdir -p GSM5924046
cd GSM5924046
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924047
mkdir -p GSM5924047
cd GSM5924047
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924048
mkdir -p GSM5924048
cd GSM5924048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924049
mkdir -p GSM5924049
cd GSM5924049
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924050
mkdir -p GSM5924050
cd GSM5924050
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924051
mkdir -p GSM5924051
cd GSM5924051
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924052
mkdir -p GSM5924052
cd GSM5924052
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924053
mkdir -p GSM5924053
cd GSM5924053
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#HG_1.3
mkdir -p HG_1.3
cd HG_1.3
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_1.2
mkdir -p HG_1.2
cd HG_1.2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_3.2
mkdir -p HG_3.2
cd HG_3.2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_1
mkdir -p LG_1
cd LG_1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_2
mkdir -p LG_2
cd LG_2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_3
mkdir -p LG_3
cd LG_3
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_1.1
mkdir -p pHG_1.1
cd pHG_1.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_3.1
mkdir -p pHG_3.1
cd pHG_3.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


# GSM7974884
mkdir -p GSM7974884
cd GSM7974884
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../


# GSM7974885
mkdir -p GSM7974885
cd GSM7974885
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974886
mkdir -p GSM7974886
cd GSM7974886
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974887
mkdir -p GSM7974887
cd GSM7974887
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974888
mkdir -p GSM7974888
cd GSM7974888
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6716963
mkdir -p GSM6716963
cd GSM6716963
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716964
mkdir -p GSM6716964
cd GSM6716964
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716965
mkdir -p GSM6716965
cd GSM6716965
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716966
mkdir -p GSM6716966
cd GSM6716966
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_11mm_FFPE_Human_Colorectal_Cancer
mkdir -p CytAssist_11mm_FFPE_Human_Colorectal_Cancer
cd CytAssist_11mm_FFPE_Human_Colorectal_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Parent_Visium_Human_ColorectalCancer
mkdir -p Parent_Visium_Human_ColorectalCancer
cd Parent_Visium_Human_ColorectalCancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep1
mkdir -p CytAssist_FFPE_Human_Colon_Rep1
cd CytAssist_FFPE_Human_Colon_Rep1
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep2
mkdir -p CytAssist_FFPE_Human_Colon_Rep2
cd CytAssist_FFPE_Human_Colon_Rep2
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep2/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep2_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Targeted_Visium_Human_ColorectalCancer_GeneSignature
mkdir -p Targeted_Visium_Human_ColorectalCancer_GeneSignature
cd Targeted_Visium_Human_ColorectalCancer_GeneSignature
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM8041061 
mkdir -p GSM8041061
mkdir -p GSM8041061/filtered_feature_bc_matrix
cd GSM8041061/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041061/filtered_feature_bc_matrix --output_file GSM8041061/filtered_feature_bc_matrix.h5
mkdir -p GSM8041061/spatial
cd GSM8041061/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8041062
mkdir -p GSM8041062
mkdir -p GSM8041062/filtered_feature_bc_matrix
cd GSM8041062/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041062/filtered_feature_bc_matrix --output_file GSM8041062/filtered_feature_bc_matrix.h5
mkdir -p GSM8041062/spatial
cd GSM8041062/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6256810
mkdir -p GSM6256810
cd GSM6256810
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
# tissue positions
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256811
mkdir -p GSM6256811
cd GSM6256811
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256812
mkdir -p GSM6256812
cd GSM6256812
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256813
mkdir -p GSM6256813
cd GSM6256813
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7573205
mkdir -p GSM7573205
mkdir -p GSM7573205/filtered_feature_bc_matrix
cd GSM7573205/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM7573205/filtered_feature_bc_matrix \
--output_file GSM7573205/filtered_feature_bc_matrix.h5
mkdir -p GSM7573205/spatial
cd GSM7573205/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Escalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265211
mkdir -p GSM8265211
mkdir -p GSM8265211/filtered_feature_bc_matrix
cd GSM8265211/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265211/filtered_feature_bc_matrix \
--output_file GSM8265211/filtered_feature_bc_matrix.h5
mkdir -p GSM8265211/spatial
cd GSM8265211/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265212
mkdir -p GSM8265212
mkdir -p GSM8265212/filtered_feature_bc_matrix
cd GSM8265212/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265212/filtered_feature_bc_matrix \
--output_file GSM8265212/filtered_feature_bc_matrix.h5
mkdir -p GSM8265212/spatial
cd GSM8265212/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265213
mkdir -p GSM8265213
mkdir -p GSM8265213/filtered_feature_bc_matrix
cd GSM8265213/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265213/filtered_feature_bc_matrix \
--output_file GSM8265213/filtered_feature_bc_matrix.h5
mkdir -p GSM8265213/spatial
cd GSM8265213/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#zenodo14204217 ??? 
mkdir -p zenodo14204217
cd zenodo14204217
wget -c https://zenodo.org/records/14204217/files/byArray.tar
tar -xvf byArray.tar 
wget -c https://zenodo.org/records/14204217/files/Clinical.tar
tar -xvf Clinical.tar
wget -c https://zenodo.org/records/14204217/files/rawCountsMatrices.tar
tar -xvf rawCountsMatrices.tar
Rscript convert_zenodo_to_10x.R   --raw_dir /data_d/WSJ/SpatialMetsDB/RawData/zenodo14204217   --out_dir /data_d/WSJ/SpatialMetsDB/RawData   --ncores 4   --overwrite

#GSM6177599
mkdir -p GSM6177599
cd GSM6177599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6177601
mkdir -p GSM6177601
cd GSM6177601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177599&format=file&file=GSM6177599%5FNYU%5FBRCA0%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177601&format=file&file=GSM6177601%5FNYU%5FBRCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6177603
mkdir -p GSM6177603
cd GSM6177603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177603&format=file&file=GSM6177603%5FNYU%5FBRCA2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5420752
mkdir -p GSM5420752
cd GSM5420752
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420752&format=file&file=GSM5420752%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415705
mkdir -p GSM6415705
cd GSM6415705
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415705&format=file&file=GSM6415705%5FGLMF1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM6415706
mkdir -p GSM6415706
cd GSM6415706
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6415706&format=file&file=GSM6415706%5FGLMF2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ../

# GSM5924030
mkdir -p GSM5924030
cd GSM5924030
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924030&format=file&file=GSM5924030%5Fffpe%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924031
mkdir -p GSM5924031
cd GSM5924031
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924031&format=file&file=GSM5924031%5Fffpe%5Fc%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924032
mkdir -p GSM5924032
cd GSM5924032
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924032&format=file&file=GSM5924032%5Fffpe%5Fc%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924033
mkdir -p GSM5924033
cd GSM5924033
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924033&format=file&file=GSM5924033%5Fffpe%5Fc%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924034
mkdir -p GSM5924034
cd GSM5924034
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924034&format=file&file=GSM5924034%5Fffpe%5Fc%5F10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924035
mkdir -p GSM5924035
cd GSM5924035
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924035&format=file&file=GSM5924035%5Fffpe%5Fc%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924036
mkdir -p GSM5924036
cd GSM5924036
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924036&format=file&file=GSM5924036%5Fffpe%5Fc%5F21%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924037
mkdir -p GSM5924037
cd GSM5924037
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924037&format=file&file=GSM5924037%5Fffpe%5Fc%5F34%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924038
mkdir -p GSM5924038
cd GSM5924038
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924038&format=file&file=GSM5924038%5Fffpe%5Fc%5F36%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924039
mkdir -p GSM5924039
cd GSM5924039
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924039&format=file&file=GSM5924039%5Fffpe%5Fc%5F39%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924040
mkdir -p GSM5924040
cd GSM5924040
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924040&format=file&file=GSM5924040%5Fffpe%5Fc%5F45%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924041
mkdir -p GSM5924041
cd GSM5924041
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924041&format=file&file=GSM5924041%5Fffpe%5Fc%5F51%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924042
mkdir -p GSM5924042
cd GSM5924042
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924042&format=file&file=GSM5924042%5Ffrozen%5Fa%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924043
mkdir -p GSM5924043
cd GSM5924043
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924043&format=file&file=GSM5924043%5Ffrozen%5Fa%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924044
mkdir -p GSM5924044
cd GSM5924044
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924044&format=file&file=GSM5924044%5Ffrozen%5Fa%5F15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924045
mkdir -p GSM5924045
cd GSM5924045
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924045&format=file&file=GSM5924045%5Ffrozen%5Fa%5F17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924046
mkdir -p GSM5924046
cd GSM5924046
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924046&format=file&file=GSM5924046%5Ffrozen%5Fb%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924047
mkdir -p GSM5924047
cd GSM5924047
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924047&format=file&file=GSM5924047%5Ffrozen%5Fb%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924048
mkdir -p GSM5924048
cd GSM5924048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924048&format=file&file=GSM5924048%5Ffrozen%5Fb%5F13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924049
mkdir -p GSM5924049
cd GSM5924049
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924049&format=file&file=GSM5924049%5Ffrozen%5Fb%5F18%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924050
mkdir -p GSM5924050
cd GSM5924050
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924050&format=file&file=GSM5924050%5Ffrozen%5Fc%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924051
mkdir -p GSM5924051
cd GSM5924051
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924051&format=file&file=GSM5924051%5Ffrozen%5Fc%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924052
mkdir -p GSM5924052
cd GSM5924052
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924052&format=file&file=GSM5924052%5Ffrozen%5Fc%5F23%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5924053
mkdir -p GSM5924053
cd GSM5924053
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5924053&format=file&file=GSM5924053%5Ffrozen%5Fc%5F57%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#HG_1.3
mkdir -p HG_1.3
cd HG_1.3
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_1.2
mkdir -p HG_1.2
cd HG_1.2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_1.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_2
mkdir -p HG_2
cd HG_2
wget -O filtered_feature_bc_matrix.h5 \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_tissue_hires_image.png"
wget -O scalefactors_json.json \
"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#HG_3.2
mkdir -p HG_3.2
cd HG_3.2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/HG_3.2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_1
mkdir -p LG_1
cd LG_1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_2
mkdir -p LG_2
cd LG_2
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_2_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#LG_3
mkdir -p LG_3
cd LG_3
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/LG_3_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_1.1
mkdir -p pHG_1.1
cd pHG_1.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_1.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#pHG_3.1
mkdir -p pHG_3.1
cd pHG_3.1
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_filtered_feature_bc_matrix.h5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_positions_list.csv"
wget -O tissue_hires_image.png "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_tissue_hires_image.png"
wget -O scalefactors_json.json "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/pHG_3.1_scalefactors_json.json"
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


# GSM7974884
mkdir -p GSM7974884
cd GSM7974884
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974884&format=file&file=GSM7974884%5FS1%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../


# GSM7974885
mkdir -p GSM7974885
cd GSM7974885
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974885&format=file&file=GSM7974885%5FS1%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974886
mkdir -p GSM7974886
cd GSM7974886
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974886&format=file&file=GSM7974886%5FS2%2DB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974887
mkdir -p GSM7974887
cd GSM7974887
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974887&format=file&file=GSM7974887%5FS2%2DC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7974888
mkdir -p GSM7974888
cd GSM7974888
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7974888&format=file&file=GSM7974888%5FS2%2DD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

#GSM6716963
mkdir -p GSM6716963
cd GSM6716963
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716963&format=file&file=GSM6716963%5F19G081%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716964
mkdir -p GSM6716964
cd GSM6716964
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716964&format=file&file=GSM6716964%5F19G0619%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716965
mkdir -p GSM6716965
cd GSM6716965
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716965&format=file&file=GSM6716965%5F19G0635%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#GSM6716966
mkdir -p GSM6716966
cd GSM6716966
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6716966&format=file&file=GSM6716966%5F19G02977%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_11mm_FFPE_Human_Colorectal_Cancer
mkdir -p CytAssist_11mm_FFPE_Human_Colorectal_Cancer
cd CytAssist_11mm_FFPE_Human_Colorectal_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Parent_Visium_Human_ColorectalCancer
mkdir -p Parent_Visium_Human_ColorectalCancer
cd Parent_Visium_Human_ColorectalCancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep1
mkdir -p CytAssist_FFPE_Human_Colon_Rep1
cd CytAssist_FFPE_Human_Colon_Rep1
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep1_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#CytAssist_FFPE_Human_Colon_Rep2
mkdir -p CytAssist_FFPE_Human_Colon_Rep2
cd CytAssist_FFPE_Human_Colon_Rep2
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep2/CytAssist_FFPE_Human_Colon_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Rep1/CytAssist_FFPE_Human_Colon_Rep2_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

#Targeted_Visium_Human_ColorectalCancer_GeneSignature
mkdir -p Targeted_Visium_Human_ColorectalCancer_GeneSignature
cd Targeted_Visium_Human_ColorectalCancer_GeneSignature
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_ColorectalCancer_GeneSignature/Targeted_Visium_Human_ColorectalCancer_GeneSignature_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM8041061 
mkdir -p GSM8041061
mkdir -p GSM8041061/filtered_feature_bc_matrix
cd GSM8041061/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041061/filtered_feature_bc_matrix --output_file GSM8041061/filtered_feature_bc_matrix.h5
mkdir -p GSM8041061/spatial
cd GSM8041061/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041061&format=file&file=GSM8041061%5FCLM%5F01%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8041062
mkdir -p GSM8041062
mkdir -p GSM8041062/filtered_feature_bc_matrix
cd GSM8041062/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8041062/filtered_feature_bc_matrix --output_file GSM8041062/filtered_feature_bc_matrix.h5
mkdir -p GSM8041062/spatial
cd GSM8041062/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8041062&format=file&file=GSM8041062%5FCLM%5F02%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6256810
mkdir -p GSM6256810
cd GSM6256810
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
# tissue positions
wget -O tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256810&format=file&file=GSM6256810%5Fmeta1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256811
mkdir -p GSM6256811
cd GSM6256811
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256811&format=file&file=GSM6256811%5Fmeta2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256812
mkdir -p GSM6256812
cd GSM6256812
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256812&format=file&file=GSM6256812%5Fmeta3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6256813
mkdir -p GSM6256813
cd GSM6256813
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6256813&format=file&file=GSM6256813%5Fmeta4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7573205
mkdir -p GSM7573205
mkdir -p GSM7573205/filtered_feature_bc_matrix
cd GSM7573205/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM7573205/filtered_feature_bc_matrix \
--output_file GSM7573205/filtered_feature_bc_matrix.h5
mkdir -p GSM7573205/spatial
cd GSM7573205/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Etissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7573205&format=file&file=GSM7573205%5Fp1%2Escalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265211
mkdir -p GSM8265211
mkdir -p GSM8265211/filtered_feature_bc_matrix
cd GSM8265211/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265211/filtered_feature_bc_matrix \
--output_file GSM8265211/filtered_feature_bc_matrix.h5
mkdir -p GSM8265211/spatial
cd GSM8265211/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265211&format=file&file=GSM8265211%5FCTC21P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265212
mkdir -p GSM8265212
mkdir -p GSM8265212/filtered_feature_bc_matrix
cd GSM8265212/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265212/filtered_feature_bc_matrix \
--output_file GSM8265212/filtered_feature_bc_matrix.h5
mkdir -p GSM8265212/spatial
cd GSM8265212/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265212&format=file&file=GSM8265212%5FCTC21M%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265213
mkdir -p GSM8265213
mkdir -p GSM8265213/filtered_feature_bc_matrix
cd GSM8265213/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265213/filtered_feature_bc_matrix \
--output_file GSM8265213/filtered_feature_bc_matrix.h5
mkdir -p GSM8265213/spatial
cd GSM8265213/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265213&format=file&file=GSM8265213%5FCTC17P%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8265214
mkdir -p GSM8265214
mkdir -p GSM8265214/filtered_feature_bc_matrix
cd GSM8265214/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM8265214/filtered_feature_bc_matrix \
--output_file GSM8265214/filtered_feature_bc_matrix.h5
mkdir -p GSM8265214/spatial
cd GSM8265214/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8265214&format=file&file=GSM8265214%5FCTC17M%5Fscalefactors%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#SN048_A121573_Rep1
mkdir SN048_A121573_Rep1
cd SN048_A121573_Rep1
wget -O SN048_A121573_Rep1.zip "https://zenodo.org/records/7760264/files/SN048_A121573_Rep1.zip?download=1"
unzip SN048_A121573_Rep1.zip
mv SN048_A121573_Rep1/* . 2>/dev/null || true
rm -rf SN048_A121573_Rep1
rm -f SN048_A121573_Rep1.zip
cd ..

#SN048_A121573_Rep2
mkdir SN048_A121573_Rep2
cd SN048_A121573_Rep2
wget -O SN048_A121573_Rep2.zip "https://zenodo.org/records/7760264/files/SN048_A121573_Rep2.zip?download=1"
unzip SN048_A121573_Rep2.zip
mv SN048_A121573_Rep2/* . 2>/dev/null || true
rm -rf SN048_A121573_Rep2
rm -f SN048_A121573_Rep2.zip
cd ..

#SN048_A416371_Rep1
mkdir SN048_A416371_Rep1
cd SN048_A416371_Rep1
wget -O SN048_A416371_Rep1.zip "https://zenodo.org/records/7760264/files/SN048_A416371_Rep1.zip?download=1"
unzip SN048_A416371_Rep1.zip
mv SN048_A416371_Rep1/* . 2>/dev/null || true
rm -rf SN048_A416371_Rep1
rm -f SN048_A416371_Rep1.zip
cd ..

#SN048_A416371_Rep2
mkdir SN048_A416371_Rep2
cd SN048_A416371_Rep2
wget -O SN048_A416371_Rep2.zip "https://zenodo.org/records/7760264/files/SN048_A416371_Rep2.zip?download=1"
unzip SN048_A416371_Rep2.zip
mv SN048_A416371_Rep2/* . 2>/dev/null || true
rm -rf SN048_A416371_Rep2
rm -f SN048_A416371_Rep2.zip
cd ..

#SN123_A551763_Rep1
mkdir SN123_A551763_Rep1
cd SN123_A551763_Rep1
wget -O SN123_A551763_Rep1.zip "https://zenodo.org/records/7760264/files/SN123_A551763_Rep1.zip?download=1"
unzip SN123_A551763_Rep1.zip
mv SN123_A551763_Rep1/* . 2>/dev/null || true
rm -rf SN123_A551763_Rep1
rm -f SN123_A551763_Rep1.zip
cd ..

#SN123_A595688_Rep1
mkdir SN123_A595688_Rep1
cd SN123_A595688_Rep1
wget -O SN123_A595688_Rep1.zip "https://zenodo.org/records/7760264/files/SN123_A595688_Rep1.zip?download=1"
unzip SN123_A595688_Rep1.zip
mv SN123_A595688_Rep1/* . 2>/dev/null || true
rm -rf SN123_A595688_Rep1
rm -f SN123_A595688_Rep1.zip
cd ..

#SN123_A938797_Rep1_X
mkdir SN123_A938797_Rep1_X
cd SN123_A938797_Rep1_X
wget -O SN123_A938797_Rep1_X.zip "https://zenodo.org/records/7760264/files/SN123_A938797_Rep1_X.zip?download=1"
unzip SN123_A938797_Rep1_X.zip
mv SN123_A938797_Rep1_X/* . 2>/dev/null || true
rm -rf SN123_A938797_Rep1_X
rm -f SN123_A938797_Rep1_X.zip
cd ..

#SN124_A551763_Rep2
mkdir SN124_A551763_Rep2
cd SN124_A551763_Rep2
wget -O SN124_A551763_Rep2.zip "https://zenodo.org/records/7760264/files/SN124_A551763_Rep2.zip?download=1"
unzip SN124_A551763_Rep2.zip
mv SN124_A551763_Rep2/* . 2>/dev/null || true
rm -rf SN124_A551763_Rep2
rm -f SN124_A551763_Rep2.zip
cd ..

#SN124_A938797_Rep2
mkdir SN124_A938797_Rep2
cd SN124_A938797_Rep2
wget -O SN124_A938797_Rep2.zip "https://zenodo.org/records/7760264/files/SN124_A938797_Rep2.zip?download=1"
unzip SN124_A938797_Rep2.zip
mv SN124_A938797_Rep2/* . 2>/dev/null || true
rm -rf SN124_A938797_Rep2
rm -f SN124_A938797_Rep2.zip
cd ..

#SN84_A120838_Rep1
mkdir SN84_A120838_Rep1
cd SN84_A120838_Rep1
wget -O SN84_A120838_Rep1.zip "https://zenodo.org/records/7760264/files/SN84_A120838_Rep1.zip?download=1"
unzip SN84_A120838_Rep1.zip
mv SN84_A120838_Rep1/* . 2>/dev/null || true
rm -rf SN84_A120838_Rep1
rm -f SN84_A120838_Rep1.zip
cd ..

#SN84_A120838_Rep2
mkdir SN84_A120838_Rep2
cd SN84_A120838_Rep2
wget -O SN84_A120838_Rep2.zip "https://zenodo.org/records/7760264/files/SN84_A120838_Rep2.zip?download=1"
unzip SN84_A120838_Rep2.zip
mv SN84_A120838_Rep2/* . 2>/dev/null || true
rm -rf SN84_A120838_Rep2
rm -f SN84_A120838_Rep2.zip
cd ..

# GSM8594561
mkdir -p GSM8594561
cd GSM8594561
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8594561&format=file&file=GSM8594561%5FP2CRC%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8594561&format=file&file=GSM8594561%5FP2CRC%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8594561&format=file&file=GSM8594561%5FP2CRC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8594561&format=file&file=GSM8594561%5FP2CRC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7058756
mkdir -p GSM7058756
mkdir -p GSM7058756/filtered_feature_bc_matrix
cd GSM7058756/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R \
  --input_dir GSM7058756/filtered_feature_bc_matrix \
  --output_file GSM7058756/filtered_feature_bc_matrix.h5
mkdir -p GSM7058756/spatial
cd GSM7058756/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058756&format=file&file=GSM7058756%5FC1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7058757
mkdir -p GSM7058757
mkdir -p GSM7058757/filtered_feature_bc_matrix
cd GSM7058757/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7058757/filtered_feature_bc_matrix --output_file GSM7058757/filtered_feature_bc_matrix.h5
mkdir -p GSM7058757/spatial
cd GSM7058757/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058757&format=file&file=GSM7058757%5FC2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7058758
mkdir -p GSM7058758
mkdir -p GSM7058758/filtered_feature_bc_matrix
cd GSM7058758/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7058758/filtered_feature_bc_matrix --output_file GSM7058758/filtered_feature_bc_matrix.h5
mkdir -p GSM7058758/spatial
cd GSM7058758/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058758&format=file&file=GSM7058758%5FC3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7058759
mkdir -p GSM7058759
mkdir -p GSM7058759/filtered_feature_bc_matrix
cd GSM7058759/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7058759/filtered_feature_bc_matrix --output_file GSM7058759/filtered_feature_bc_matrix.h5
mkdir -p GSM7058759/spatial
cd GSM7058759/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058759&format=file&file=GSM7058759%5FC4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7058760
mkdir -p GSM7058760
mkdir -p GSM7058760/filtered_feature_bc_matrix
cd GSM7058760/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7058760/filtered_feature_bc_matrix --output_file GSM7058760/filtered_feature_bc_matrix.h5
mkdir -p GSM7058760/spatial
cd GSM7058760/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058760&format=file&file=GSM7058760%5FL1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7058761
mkdir -p GSM7058761
mkdir -p GSM7058761/filtered_feature_bc_matrix
cd GSM7058761/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%2Ematrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%2Efeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%2Ebarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7058761/filtered_feature_bc_matrix --output_file GSM7058761/filtered_feature_bc_matrix.h5
mkdir -p GSM7058761/spatial
cd GSM7058761/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7058761&format=file&file=GSM7058761%5FL2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM5420749
mkdir -p GSM5420749
cd GSM5420749
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420749&format=file&file=GSM5420749%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420749&format=file&file=GSM5420749%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM5420751
mkdir -p GSM5420751
cd GSM5420751
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420751&format=file&file=GSM5420751%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420751&format=file&file=GSM5420751%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM5420754
mkdir -p GSM5420754
cd GSM5420754
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420754&format=file&file=GSM5420754%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420754&format=file&file=GSM5420754%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# CytAssist_11mm_FFPE_Human_Lung_Cancer 
mkdir -p CytAssist_11mm_FFPE_Human_Lung_Cancer
cd CytAssist_11mm_FFPE_Human_Lung_Cancer
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz
rm -f spatial.tar.gz
cd ..

# GSM8284978
mkdir -p GSM8284978
mkdir -p GSM8284978/filtered_feature_bc_matrix
cd GSM8284978/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284978/filtered_feature_bc_matrix --output_file GSM8284978/filtered_feature_bc_matrix.h5
mkdir -p GSM8284978/spatial
cd GSM8284978/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284978&format=file&file=GSM8284978%5FLung%5F1%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284979 
mkdir -p GSM8284979
mkdir -p GSM8284979/filtered_feature_bc_matrix
cd GSM8284979/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284979/filtered_feature_bc_matrix --output_file GSM8284979/filtered_feature_bc_matrix.h5
mkdir -p GSM8284979/spatial
cd GSM8284979/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284979&format=file&file=GSM8284979%5FLung%5F2%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284980 
mkdir -p GSM8284980
mkdir -p GSM8284980/filtered_feature_bc_matrix
cd GSM8284980/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284980/filtered_feature_bc_matrix --output_file GSM8284980/filtered_feature_bc_matrix.h5
mkdir -p GSM8284980/spatial
cd GSM8284980/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284980&format=file&file=GSM8284980%5FLung%5F3%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284981 
mkdir -p GSM8284981
mkdir -p GSM8284981/filtered_feature_bc_matrix
cd GSM8284981/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284981/filtered_feature_bc_matrix --output_file GSM8284981/filtered_feature_bc_matrix.h5
mkdir -p GSM8284981/spatial
cd GSM8284981/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284981&format=file&file=GSM8284981%5FLung%5F4%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284982 
mkdir -p GSM8284982
mkdir -p GSM8284982/filtered_feature_bc_matrix
cd GSM8284982/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284982/filtered_feature_bc_matrix --output_file GSM8284982/filtered_feature_bc_matrix.h5
mkdir -p GSM8284982/spatial
cd GSM8284982/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284982&format=file&file=GSM8284982%5FLung%5F6%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284983
mkdir -p GSM8284983
mkdir -p GSM8284983/filtered_feature_bc_matrix
cd GSM8284983/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284983/filtered_feature_bc_matrix --output_file GSM8284983/filtered_feature_bc_matrix.h5
mkdir -p GSM8284983/spatial
cd GSM8284983/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284983&format=file&file=GSM8284983%5FLung%5F7%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284984
mkdir -p GSM8284984
mkdir -p GSM8284984/filtered_feature_bc_matrix
cd GSM8284984/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284984/filtered_feature_bc_matrix --output_file GSM8284984/filtered_feature_bc_matrix.h5
mkdir -p GSM8284984/spatial
cd GSM8284984/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284984&format=file&file=GSM8284984%5FLung%5F8%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284985 
mkdir -p GSM8284985
mkdir -p GSM8284985/filtered_feature_bc_matrix
cd GSM8284985/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284985/filtered_feature_bc_matrix --output_file GSM8284985/filtered_feature_bc_matrix.h5
mkdir -p GSM8284985/spatial
cd GSM8284985/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284985&format=file&file=GSM8284985%5FLung%5F9%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284986
mkdir -p GSM8284986
mkdir -p GSM8284986/filtered_feature_bc_matrix
cd GSM8284986/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284986/filtered_feature_bc_matrix --output_file GSM8284986/filtered_feature_bc_matrix.h5
mkdir -p GSM8284986/spatial
cd GSM8284986/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284986&format=file&file=GSM8284986%5FLung%5F10%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284987 
mkdir -p GSM8284987
mkdir -p GSM8284987/filtered_feature_bc_matrix
cd GSM8284987/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284987/filtered_feature_bc_matrix --output_file GSM8284987/filtered_feature_bc_matrix.h5
mkdir -p GSM8284987/spatial
cd GSM8284987/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284987&format=file&file=GSM8284987%5FLung%5F11%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284988 
mkdir -p GSM8284988
mkdir -p GSM8284988/filtered_feature_bc_matrix
cd GSM8284988/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284988/filtered_feature_bc_matrix --output_file GSM8284988/filtered_feature_bc_matrix.h5
mkdir -p GSM8284988/spatial
cd GSM8284988/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284988&format=file&file=GSM8284988%5FLung%5F12%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284989 
mkdir -p GSM8284989
mkdir -p GSM8284989/filtered_feature_bc_matrix
cd GSM8284989/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284989/filtered_feature_bc_matrix --output_file GSM8284989/filtered_feature_bc_matrix.h5
mkdir -p GSM8284989/spatial
cd GSM8284989/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284989&format=file&file=GSM8284989%5FLung%5F13%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284990
mkdir -p GSM8284990
mkdir -p GSM8284990/filtered_feature_bc_matrix
cd GSM8284990/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284990/filtered_feature_bc_matrix --output_file GSM8284990/filtered_feature_bc_matrix.h5
mkdir -p GSM8284990/spatial
cd GSM8284990/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284990&format=file&file=GSM8284990%5FLung%5F14%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284991 
mkdir -p GSM8284991
mkdir -p GSM8284991/filtered_feature_bc_matrix
cd GSM8284991/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284991/filtered_feature_bc_matrix --output_file GSM8284991/filtered_feature_bc_matrix.h5
mkdir -p GSM8284991/spatial
cd GSM8284991/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284991&format=file&file=GSM8284991%5FLung%5F15%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8284992
mkdir -p GSM8284992
mkdir -p GSM8284992/filtered_feature_bc_matrix
cd GSM8284992/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8284992/filtered_feature_bc_matrix --output_file GSM8284992/filtered_feature_bc_matrix.h5
mkdir -p GSM8284992/spatial
cd GSM8284992/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8284992&format=file&file=GSM8284992%5FLung%5F16%5Ftumor%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# P10_T1
mkdir -p P10_T1
cd P10_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P10_T2
mkdir -p P10_T2
cd P10_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P10_T3
mkdir -p P10_T3
cd P10_T3
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T3-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T3-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P10_T4
mkdir -p P10_T4
cd P10_T4
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T4-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P10_T4-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P11_T1
mkdir -p P11_T1
cd P11_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P11_T2
mkdir -p P11_T2
cd P11_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P11_T3
mkdir -p P11_T3
cd P11_T3
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T3-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T3-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P11_T4
mkdir -p P11_T4
cd P11_T4
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T4-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P11_T4-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P15_T1
mkdir -p P15_T1
cd P15_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P15_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P15_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P15_T2
mkdir -p P15_T2
cd P15_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P15_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P15_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P16_T1
mkdir -p P16_T1
cd P16_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P16_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P16_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P16_T2
mkdir -p P16_T2
cd P16_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P16_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P16_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P17_T1
mkdir -p P17_T1
cd P17_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P17_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P17_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P17_T2
mkdir -p P17_T2
cd P17_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P17_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P17_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P19_T1
mkdir -p P19_T1
cd P19_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P19_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P19_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P19_T2
mkdir -p P19_T2
cd P19_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P19_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P19_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P24_T1
mkdir -p P24_T1
cd P24_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P24_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P24_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P24_T2
mkdir -p P24_T2
cd P24_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P24_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P24_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P25_T1
mkdir -p P25_T1
cd P25_T1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P25_T1-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P25_T1-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# P25_T2
mkdir -p P25_T2
cd P25_T2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P25_T2-filtered_feature_bc_matrix.h5"
wget -O spatial.tar "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/P25_T2-spatial.tar"
tar -xf spatial.tar -C spatial
rm -f spatial.tar
cd ..

# GSM6177618
mkdir -p GSM6177618
cd GSM6177618
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177618&format=file&file=GSM6177618%5FNYU%5FPDAC1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177618&format=file&file=GSM6177618%5FNYU%5FPDAC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177618&format=file&file=GSM6177618%5FNYU%5FPDAC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177618&format=file&file=GSM6177618%5FNYU%5FPDAC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177618&format=file&file=GSM6177618%5FNYU%5FPDAC1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7021871
mkdir -p GSM7021871
cd GSM7021871
wget -O GSM7021871_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7021871&format=file&file=GSM7021871%5FPanIN%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM7021871_PanIN.tar.gz --strip-components=1
rm GSM7021871_PanIN.tar.gz
cd ..

# GSM7021872
mkdir -p GSM7021872
cd GSM7021872
wget -O GSM7021872_PDAClymphnode.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7021872&format=file&file=GSM7021872%5FPDAClymphnode%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM7021872_PDAClymphnode.tar.gz --strip-components=1
rm GSM7021872_PDAClymphnode.tar.gz
cd ..

# GSM7421780
mkdir -p GSM7421780
cd GSM7421780
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421780&format=file&file=GSM7421780%5FLG%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421780&format=file&file=GSM7421780%5FLG%5F1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421781
mkdir -p GSM7421781
cd GSM7421781
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421781&format=file&file=GSM7421781%5FLG%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421781&format=file&file=GSM7421781%5FLG%5F2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
rm -f spatial.tar.gz
cd ..

# GSM7421782
mkdir -p GSM7421782
cd GSM7421782
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421782&format=file&file=GSM7421782%5FLG%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421782&format=file&file=GSM7421782%5FLG%5F3%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421783
mkdir -p GSM7421783
cd GSM7421783
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421783&format=file&file=GSM7421783%5FLG%5F4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421783&format=file&file=GSM7421783%5FLG%5F4%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421784
mkdir -p GSM7421784
cd GSM7421784
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421784&format=file&file=GSM7421784%5FLG%5F5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421784&format=file&file=GSM7421784%5FLG%5F5%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421785
mkdir -p GSM7421785
cd GSM7421785
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421785&format=file&file=GSM7421785%5FLG%5F6%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421785&format=file&file=GSM7421785%5FLG%5F6%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421786
mkdir -p GSM7421786
cd GSM7421786
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421786&format=file&file=GSM7421786%5FLG%5F7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421786&format=file&file=GSM7421786%5FLG%5F7%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421787
mkdir -p GSM7421787
cd GSM7421787
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421787&format=file&file=GSM7421787%5FHG%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421787&format=file&file=GSM7421787%5FHG%5F1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421788
mkdir -p GSM7421788
cd GSM7421788
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421788&format=file&file=GSM7421788%5FHG%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421788&format=file&file=GSM7421788%5FHG%5F2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421789
mkdir -p GSM7421789
cd GSM7421789
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421789&format=file&file=GSM7421789%5FHG%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421789&format=file&file=GSM7421789%5FHG%5F3%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421790
mkdir -p GSM7421790
cd GSM7421790
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421790&format=file&file=GSM7421790%5FPDAC%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421790&format=file&file=GSM7421790%5FPDAC%5F1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7421791
mkdir -p GSM7421791
cd GSM7421791
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421791&format=file&file=GSM7421791%5FPDAC%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7421791&format=file&file=GSM7421791%5FPDAC%5F2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8058242
mkdir -p GSM8058242
cd GSM8058242
wget -O GSM8058242_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058242&format=file&file=GSM8058242%5FPanIN%2DLG1%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058242_PanIN.tar.gz --strip-components=1
rm GSM8058242_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058242/filtered_feature_bc_matrix --output_file GSM8058242/filtered_feature_bc_matrix.h5

# GSM8058243
mkdir -p GSM8058243
cd GSM8058243
wget -O GSM8058243_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058243&format=file&file=GSM8058243%5FPanIN%2DHG1%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058243_PanIN.tar.gz --strip-components=1
rm GSM8058243_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058243/filtered_feature_bc_matrix --output_file GSM8058243/filtered_feature_bc_matrix.h5

# GSM8058244
mkdir -p GSM8058244
cd GSM8058244
wget -O GSM8058244_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058244&format=file&file=GSM8058244%5FPanIN%2DLG2%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058244_PanIN.tar.gz --strip-components=1
rm GSM8058244_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058244/filtered_feature_bc_matrix --output_file GSM8058244/filtered_feature_bc_matrix.h5

# GSM8058245
mkdir -p GSM8058245
cd GSM8058245
wget -O GSM8058245_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058245&format=file&file=GSM8058245%5FPanIN%2DHG2%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058245_PanIN.tar.gz --strip-components=1
rm GSM8058245_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058245/filtered_feature_bc_matrix --output_file GSM8058245/filtered_feature_bc_matrix.h5

# GSM8058246
mkdir -p GSM8058246
cd GSM8058246
wget -O GSM8058246_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058246&format=file&file=GSM8058246%5FPanIN%2DLG3%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058246_PanIN.tar.gz --strip-components=1
rm GSM8058246_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058246/filtered_feature_bc_matrix --output_file GSM8058246/filtered_feature_bc_matrix.h5

# GSM8058247
mkdir -p GSM8058247
cd GSM8058247
wget -O GSM8058247_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058247&format=file&file=GSM8058247%5FPanIN%2DHG3%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058247_PanIN.tar.gz --strip-components=1
rm GSM8058247_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058247/filtered_feature_bc_matrix --output_file GSM8058247/filtered_feature_bc_matrix.h5

# GSM8058248
mkdir -p GSM8058248
cd GSM8058248
wget -O GSM8058248_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058248&format=file&file=GSM8058248%5FPanIN%2DHG4%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058248_PanIN.tar.gz --strip-components=1
rm GSM8058248_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058248/filtered_feature_bc_matrix --output_file GSM8058248/filtered_feature_bc_matrix.h5

# GSM8058249
mkdir -p GSM8058249
cd GSM8058249
wget -O GSM8058249_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058249&format=file&file=GSM8058249%5FPanIN%2DR%2DLG1%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058249_PanIN.tar.gz --strip-components=1
rm GSM8058249_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058249/filtered_feature_bc_matrix --output_file GSM8058249/filtered_feature_bc_matrix.h5

# GSM8058250
mkdir -p GSM8058250
cd GSM8058250
wget -O GSM8058250_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058250&format=file&file=GSM8058250%5FPanIN%2DR%2DLG2%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058250_PanIN.tar.gz --strip-components=1
rm GSM8058250_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058250/filtered_feature_bc_matrix --output_file GSM8058250/filtered_feature_bc_matrix.h5

# GSM8058251
mkdir -p GSM8058251
cd GSM8058251
wget -O GSM8058251_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058251&format=file&file=GSM8058251%5FPanIN%2DR%2DLG3%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058251_PanIN.tar.gz --strip-components=1
rm GSM8058251_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058251/filtered_feature_bc_matrix --output_file GSM8058251/filtered_feature_bc_matrix.h5

# GSM8058252
mkdir -p GSM8058252
cd GSM8058252
wget -O GSM8058252_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058252&format=file&file=GSM8058252%5FPanIN%2DR%2DLG4%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058252_PanIN.tar.gz --strip-components=1
rm GSM8058252_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058252/filtered_feature_bc_matrix --output_file GSM8058252/filtered_feature_bc_matrix.h5

# GSM8058253
mkdir -p GSM8058253
cd GSM8058253
wget -O GSM8058253_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058253&format=file&file=GSM8058253%5FPanIN%2DR%2DLG5%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058253_PanIN.tar.gz --strip-components=1
rm GSM8058253_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058253/filtered_feature_bc_matrix --output_file GSM8058253/filtered_feature_bc_matrix.h5

# GSM8058254
mkdir -p GSM8058254
cd GSM8058254
wget -O GSM8058254_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058254&format=file&file=GSM8058254%5FPanIN%2DR%2DLG6%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058254_PanIN.tar.gz --strip-components=1
rm GSM8058254_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058254/filtered_feature_bc_matrix --output_file GSM8058254/filtered_feature_bc_matrix.h5

# GSM8058255
mkdir -p GSM8058255
cd GSM8058255
wget -O GSM8058255_PanIN.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8058255&format=file&file=GSM8058255%5FPanIN%2DR%2DHG1%2Etar%2Egz"
tar --warning=no-unknown-keyword -xzf GSM8058255_PanIN.tar.gz --strip-components=1
rm GSM8058255_PanIN.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8058255/filtered_feature_bc_matrix --output_file GSM8058255/filtered_feature_bc_matrix.h5

# GSM8443449
mkdir -p GSM8443449
cd GSM8443449
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443449&format=file&file=GSM8443449%5FPDAC%2Dp1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443449&format=file&file=GSM8443449%5FPDAC%2Dp1%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8443450
mkdir -p GSM8443450
cd GSM8443450
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443450&format=file&file=GSM8443450%5FPDAC%2Dp2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443450&format=file&file=GSM8443450%5FPDAC%2Dp2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8443451
mkdir -p GSM8443451
cd GSM8443451
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443451&format=file&file=GSM8443451%5FPDAC%2Dp3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443451&format=file&file=GSM8443451%5FPDAC%2Dp3%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8443452
mkdir -p GSM8443452
cd GSM8443452
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443452&format=file&file=GSM8443452%5FPDAC%2Dp4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443452&format=file&file=GSM8443452%5FPDAC%2Dp4%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8443453
mkdir -p GSM8443453
cd GSM8443453
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443453&format=file&file=GSM8443453%5FPDAC%2Dp5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8443453&format=file&file=GSM8443453%5FPDAC%2Dp5%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452847
mkdir -p GSM8452847
cd GSM8452847
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452847&format=file&file=GSM8452847%5FPt%2D1A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452847&format=file&file=GSM8452847%5FPt%2D1A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452848
mkdir -p GSM8452848
cd GSM8452848
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452848&format=file&file=GSM8452848%5FPt%2D1B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452848&format=file&file=GSM8452848%5FPt%2D1B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452849
mkdir -p GSM8452849
cd GSM8452849
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452849&format=file&file=GSM8452849%5FPt%2D1C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452849&format=file&file=GSM8452849%5FPt%2D1C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452850
mkdir -p GSM8452850
cd GSM8452850
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452850&format=file&file=GSM8452850%5FPt%2D2A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452850&format=file&file=GSM8452850%5FPt%2D2A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452851
mkdir -p GSM8452851
cd GSM8452851
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452851&format=file&file=GSM8452851%5FPt%2D2B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452851&format=file&file=GSM8452851%5FPt%2D2B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452852
mkdir -p GSM8452852
cd GSM8452852
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452852&format=file&file=GSM8452852%5FPt%2D2C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452852&format=file&file=GSM8452852%5FPt%2D2C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452853
mkdir -p GSM8452853
cd GSM8452853
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452853&format=file&file=GSM8452853%5FPt%2D3A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452853&format=file&file=GSM8452853%5FPt%2D3A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452854
mkdir -p GSM8452854
cd GSM8452854
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452854&format=file&file=GSM8452854%5FPt%2D3B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452854&format=file&file=GSM8452854%5FPt%2D3B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452855
mkdir -p GSM8452855
cd GSM8452855
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452855&format=file&file=GSM8452855%5FPt%2D3C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452855&format=file&file=GSM8452855%5FPt%2D3C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452856
mkdir -p GSM8452856
cd GSM8452856
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452856&format=file&file=GSM8452856%5FPt%2D4A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452856&format=file&file=GSM8452856%5FPt%2D4A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452857
mkdir -p GSM8452857
cd GSM8452857
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452857&format=file&file=GSM8452857%5FPt%2D4B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452857&format=file&file=GSM8452857%5FPt%2D4B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452858
mkdir -p GSM8452858
cd GSM8452858
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452858&format=file&file=GSM8452858%5FPt%2D4C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452858&format=file&file=GSM8452858%5FPt%2D4C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452859
mkdir -p GSM8452859
cd GSM8452859
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452859&format=file&file=GSM8452859%5FPt%2D4D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452859&format=file&file=GSM8452859%5FPt%2D4D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452860
mkdir -p GSM8452860
cd GSM8452860
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452860&format=file&file=GSM8452860%5FPt%2D5A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452860&format=file&file=GSM8452860%5FPt%2D5A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452861
mkdir -p GSM8452861
cd GSM8452861
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452861&format=file&file=GSM8452861%5FPt%2D5B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452861&format=file&file=GSM8452861%5FPt%2D5B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452862
mkdir -p GSM8452862
cd GSM8452862
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452862&format=file&file=GSM8452862%5FPt%2D5C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452862&format=file&file=GSM8452862%5FPt%2D5C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452863
mkdir -p GSM8452863
cd GSM8452863
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452863&format=file&file=GSM8452863%5FPt%2D6A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452863&format=file&file=GSM8452863%5FPt%2D6A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452864
mkdir -p GSM8452864
cd GSM8452864
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452864&format=file&file=GSM8452864%5FPt%2D6B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452864&format=file&file=GSM8452864%5FPt%2D6B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452865
mkdir -p GSM8452865
cd GSM8452865
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452865&format=file&file=GSM8452865%5FPt%2D6C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452865&format=file&file=GSM8452865%5FPt%2D6C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452866
mkdir -p GSM8452866
cd GSM8452866
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452866&format=file&file=GSM8452866%5FPt%2D6D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452866&format=file&file=GSM8452866%5FPt%2D6D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452867
mkdir -p GSM8452867
cd GSM8452867
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452867&format=file&file=GSM8452867%5FPt%2D7A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452867&format=file&file=GSM8452867%5FPt%2D7A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452868
mkdir -p GSM8452868
cd GSM8452868
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452868&format=file&file=GSM8452868%5FPt%2D7B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452868&format=file&file=GSM8452868%5FPt%2D7B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452869
mkdir -p GSM8452869
cd GSM8452869
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452869&format=file&file=GSM8452869%5FPt%2D7C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452869&format=file&file=GSM8452869%5FPt%2D7C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452870
mkdir -p GSM8452870
cd GSM8452870
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452870&format=file&file=GSM8452870%5FPt%2D7D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452870&format=file&file=GSM8452870%5FPt%2D7D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452871
mkdir -p GSM8452871
cd GSM8452871
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452871&format=file&file=GSM8452871%5FPt%2D8A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452871&format=file&file=GSM8452871%5FPt%2D8A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452872
mkdir -p GSM8452872
cd GSM8452872
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452872&format=file&file=GSM8452872%5FPt%2D8B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452872&format=file&file=GSM8452872%5FPt%2D8B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452873
mkdir -p GSM8452873
cd GSM8452873
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452873&format=file&file=GSM8452873%5FPt%2D8C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452873&format=file&file=GSM8452873%5FPt%2D8C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452874
mkdir -p GSM8452874
cd GSM8452874
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452874&format=file&file=GSM8452874%5FPt%2D8D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452874&format=file&file=GSM8452874%5FPt%2D8D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452875
mkdir -p GSM8452875
cd GSM8452875
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452875&format=file&file=GSM8452875%5FPt%2D9A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452875&format=file&file=GSM8452875%5FPt%2D9A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452876
mkdir -p GSM8452876
cd GSM8452876
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452876&format=file&file=GSM8452876%5FPt%2D9B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452876&format=file&file=GSM8452876%5FPt%2D9B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452877
mkdir -p GSM8452877
cd GSM8452877
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452877&format=file&file=GSM8452877%5FPt%2D9C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452877&format=file&file=GSM8452877%5FPt%2D9C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452878
mkdir -p GSM8452878
cd GSM8452878
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452878&format=file&file=GSM8452878%5FPt%2D9D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452878&format=file&file=GSM8452878%5FPt%2D9D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452879
mkdir -p GSM8452879
cd GSM8452879
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452879&format=file&file=GSM8452879%5FPt%2D10A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452879&format=file&file=GSM8452879%5FPt%2D10A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452880
mkdir -p GSM8452880
cd GSM8452880
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452880&format=file&file=GSM8452880%5FPt%2D10B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452880&format=file&file=GSM8452880%5FPt%2D10B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452881
mkdir -p GSM8452881
cd GSM8452881
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452881&format=file&file=GSM8452881%5FPt%2D10C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452881&format=file&file=GSM8452881%5FPt%2D10C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452882
mkdir -p GSM8452882
cd GSM8452882
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452882&format=file&file=GSM8452882%5FPt%2D10D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452882&format=file&file=GSM8452882%5FPt%2D10D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452883
mkdir -p GSM8452883
cd GSM8452883
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452883&format=file&file=GSM8452883%5FPt%2D11A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452883&format=file&file=GSM8452883%5FPt%2D11A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452884
mkdir -p GSM8452884
cd GSM8452884
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452884&format=file&file=GSM8452884%5FPt%2D11B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452884&format=file&file=GSM8452884%5FPt%2D11B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452885
mkdir -p GSM8452885
cd GSM8452885
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452885&format=file&file=GSM8452885%5FPt%2D11C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452885&format=file&file=GSM8452885%5FPt%2D11C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452886
mkdir -p GSM8452886
cd GSM8452886
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452886&format=file&file=GSM8452886%5FPt%2D11D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452886&format=file&file=GSM8452886%5FPt%2D11D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452887
mkdir -p GSM8452887
cd GSM8452887
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452887&format=file&file=GSM8452887%5FPt%2D12A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452887&format=file&file=GSM8452887%5FPt%2D12A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452888
mkdir -p GSM8452888
cd GSM8452888
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452888&format=file&file=GSM8452888%5FPt%2D12B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452888&format=file&file=GSM8452888%5FPt%2D12B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452889
mkdir -p GSM8452889
cd GSM8452889
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452889&format=file&file=GSM8452889%5FPt%2D12C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452889&format=file&file=GSM8452889%5FPt%2D12C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452890
mkdir -p GSM8452890
cd GSM8452890
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452890&format=file&file=GSM8452890%5FPt%2D12D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452890&format=file&file=GSM8452890%5FPt%2D12D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452891
mkdir -p GSM8452891
cd GSM8452891
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452891&format=file&file=GSM8452891%5FPt%2D13A%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452891&format=file&file=GSM8452891%5FPt%2D13A%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452892
mkdir -p GSM8452892
cd GSM8452892
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452892&format=file&file=GSM8452892%5FPt%2D13B%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452892&format=file&file=GSM8452892%5FPt%2D13B%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452893
mkdir -p GSM8452893
cd GSM8452893
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452893&format=file&file=GSM8452893%5FPt%2D13C%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452893&format=file&file=GSM8452893%5FPt%2D13C%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452894
mkdir -p GSM8452894
cd GSM8452894
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452894&format=file&file=GSM8452894%5FPt%2D13D%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452894&format=file&file=GSM8452894%5FPt%2D13D%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8452895
mkdir -p GSM8452895
cd GSM8452895
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452895&format=file&file=GSM8452895%5FPt%2D13E%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8452895&format=file&file=GSM8452895%5FPt%2D13E%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714177
mkdir -p GSM8714177
cd GSM8714177
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714177&format=file&file=GSM8714177%5FPt3%5FPri2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714177&format=file&file=GSM8714177%5FPt3%5FPri2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714178
mkdir -p GSM8714178
cd GSM8714178
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714178&format=file&file=GSM8714178%5FPt6%5FPri2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714178&format=file&file=GSM8714178%5FPt6%5FPri2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714179
mkdir -p GSM8714179
cd GSM8714179
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714179&format=file&file=GSM8714179%5FPt6%5FPri3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714179&format=file&file=GSM8714179%5FPt6%5FPri3%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714180
mkdir -p GSM8714180
cd GSM8714180
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714180&format=file&file=GSM8714180%5FPt8%5FPri2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714180&format=file&file=GSM8714180%5FPt8%5FPri2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714181
mkdir -p GSM8714181
cd GSM8714181
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714181&format=file&file=GSM8714181%5FPt13%5FPri2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714181&format=file&file=GSM8714181%5FPt13%5FPri2%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714182
mkdir -p GSM8714182
cd GSM8714182
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714182&format=file&file=GSM8714182%5FPt13%5FPri3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714182&format=file&file=GSM8714182%5FPt13%5FPri3%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714183
mkdir -p GSM8714183
cd GSM8714183
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714183&format=file&file=GSM8714183%5FPDXPt10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714183&format=file&file=GSM8714183%5FPDXPt10%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8714184
mkdir -p GSM8714184
cd GSM8714184
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714184&format=file&file=GSM8714184%5FPDXPt13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8714184&format=file&file=GSM8714184%5FPDXPt13%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM8454387
mkdir -p GSM8454387
mkdir -p GSM8454387/filtered_feature_bc_matrix
cd GSM8454387/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8454387/filtered_feature_bc_matrix --output_file GSM8454387/filtered_feature_bc_matrix.h5
mkdir -p GSM8454387/spatial
cd GSM8454387/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454387&format=file&file=GSM8454387%5FPDAC1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8454388
mkdir -p GSM8454388
mkdir -p GSM8454388/filtered_feature_bc_matrix
cd GSM8454388/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8454388/filtered_feature_bc_matrix --output_file GSM8454388/filtered_feature_bc_matrix.h5
mkdir -p GSM8454388/spatial
cd GSM8454388/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454388&format=file&file=GSM8454388%5FPDAC2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8454389
mkdir -p GSM8454389
mkdir -p GSM8454389/filtered_feature_bc_matrix
cd GSM8454389/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8454389/filtered_feature_bc_matrix --output_file GSM8454389/filtered_feature_bc_matrix.h5
mkdir -p GSM8454389/spatial
cd GSM8454389/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454389&format=file&file=GSM8454389%5FPDAC3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8454390
mkdir -p GSM8454390
mkdir -p GSM8454390/filtered_feature_bc_matrix
cd GSM8454390/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8454390/filtered_feature_bc_matrix --output_file GSM8454390/filtered_feature_bc_matrix.h5
mkdir -p GSM8454390/spatial
cd GSM8454390/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8454390&format=file&file=GSM8454390%5FPDAC4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8651954
mkdir -p GSM8651954
mkdir -p GSM8651954/filtered_feature_bc_matrix
cd GSM8651954/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651954&format=file&file=GSM8651954%5FPC%5F3%5Fsp%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651954&format=file&file=GSM8651954%5FPC%5F3%5Fsp%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651954&format=file&file=GSM8651954%5FPC%5F3%5Fsp%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R \
--input_dir GSM8651954/filtered_feature_bc_matrix \
--output_file GSM8651954/filtered_feature_bc_matrix.h5
mkdir -p GSM8651954/spatial
cd GSM8651954/spatial
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651954&format=file&file=GSM8651954%5FPC%5F3%5Fspatial%2Etar%2Egz"
tar -xvzf spatial.tar.gz
for d in */ ; do
    [ -d "$d" ] || continue
    mv "$d"* .
    rm -rf "$d"
    break
done
rm -f spatial.tar.gz
cd ../../

# GSM8651955
mkdir -p GSM8651955
mkdir -p GSM8651955/filtered_feature_bc_matrix
cd GSM8651955/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651955&format=file&file=GSM8651955%5FLM%5F3%5Fsp%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651955&format=file&file=GSM8651955%5FLM%5F3%5Fsp%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651955&format=file&file=GSM8651955%5FLM%5F3%5Fsp%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8651955/filtered_feature_bc_matrix --output_file GSM8651955/filtered_feature_bc_matrix.h5
mkdir -p GSM8651955/spatial
cd GSM8651955/spatial
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8651955&format=file&file=GSM8651955%5FLM%5F3%5Fspatial%2Etar%2Egz"
tar -xvzf spatial.tar.gz
for d in */ ; do [ -d "$d" ] || continue; mv "$d"* .; rm -rf "$d"; break; done
rm -f spatial.tar.gz
cd ../../

# GSM7498811
mkdir -p GSM7498811
cd GSM7498811
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498811&format=file&file=GSM7498811%5FSS1905133%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS1905133_processed 2>/dev/null
cd ..

# GSM7498812
mkdir -p GSM7498812
cd GSM7498812
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498812&format=file&file=GSM7498812%5FSS1923404%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS1923404_processed 2>/dev/null
cd ..

# GSM7498813
mkdir -p GSM7498813
cd GSM7498813
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498813&format=file&file=GSM7498813%5FSS1945070%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS1945070_processed 2>/dev/null
cd ..

# GSM7498814
mkdir -p GSM7498814
cd GSM7498814
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498814&format=file&file=GSM7498814%5FSS1960050%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS1960050_processed 2>/dev/null
cd ..

# GSM7498815
mkdir -p GSM7498815
cd GSM7498815
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498815&format=file&file=GSM7498815%5FSS2002154%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS2002154_processed 2>/dev/null
cd ..

# GSM7498816
mkdir -p GSM7498816
cd GSM7498816
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498816&format=file&file=GSM7498816%5FSS2005094%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS2005094_processed 2>/dev/null
cd ..

# GSM7498817
mkdir -p GSM7498817
cd GSM7498817
wget -O raw.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7498817&format=file&file=GSM7498817%5FSS2021309%5Fprocessed%2Etar%2Egz"
tar -xvzf raw.tar.gz
rm -f raw.tar.gz
h5=$(find . -name "*.h5" | head -n 1)
if [ -n "$h5" ]; then mv "$h5" filtered_feature_bc_matrix.h5; fi
mkdir -p spatial
find . -type f -name "*tissue_lowres_image.png" -exec cp {} spatial/ \;
find . -type f -name "*tissue_hires_image.png" -exec cp {} spatial/ \;
find . -type f -name "*scalefactors_json.json" -exec cp {} spatial/ \;
find . -type f -name "*tissue_positions*.csv" -exec cp {} spatial/ \;
rm -rf SS2021309_processed 2>/dev/null
cd ..

# GSM8552941
mkdir -p GSM8552941
cd GSM8552941
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552941&format=file&file=GSM8552941%5FPA02%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552941&format=file&file=GSM8552941%5FPA02%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552942
mkdir -p GSM8552942
cd GSM8552942
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552942&format=file&file=GSM8552942%5FPA03%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552942&format=file&file=GSM8552942%5FPA03%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552943
mkdir -p GSM8552943
cd GSM8552943
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552943&format=file&file=GSM8552943%5FPA04%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552943&format=file&file=GSM8552943%5FPA04%2D1%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552944
mkdir -p GSM8552944
cd GSM8552944
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552944&format=file&file=GSM8552944%5FPA04%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552944&format=file&file=GSM8552944%5FPA04%2D2%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552945
mkdir -p GSM8552945
cd GSM8552945
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552945&format=file&file=GSM8552945%5FPA05%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552945&format=file&file=GSM8552945%5FPA05%2D1%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552946
mkdir -p GSM8552946
cd GSM8552946
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552946&format=file&file=GSM8552946%5FPA05%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552946&format=file&file=GSM8552946%5FPA05%2D2%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552947
mkdir -p GSM8552947
cd GSM8552947
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552947&format=file&file=GSM8552947%5FPA06%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552947&format=file&file=GSM8552947%5FPA06%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552948
mkdir -p GSM8552948
cd GSM8552948
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552948&format=file&file=GSM8552948%5FPA08%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552948&format=file&file=GSM8552948%5FPA08%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# GSM8552949
mkdir -p GSM8552949
cd GSM8552949
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552949&format=file&file=GSM8552949%5FPA11%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8552949&format=file&file=GSM8552949%5FPA11%2D1%5Fspatial%2Etar%2Egz"
gzip -dc spatial.tar.gz > layer1.gz
gzip -dc layer1.gz > spatial.tar
tar -xf spatial.tar -C spatial
rm -f spatial.tar.gz layer1.gz spatial.tar
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
cd ..

# Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma
mkdir -p Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma
cd Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma/Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma/Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Visium_FFPE_Human_Prostate_IF
mkdir -p Visium_FFPE_Human_Prostate_IF
cd Visium_FFPE_Human_Prostate_IF
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Visium_FFPE_Human_Prostate_Cancer
mkdir -p Visium_FFPE_Human_Prostate_Cancer
cd Visium_FFPE_Human_Prostate_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_Cancer/Visium_FFPE_Human_Prostate_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_Cancer/Visium_FFPE_Human_Prostate_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7211257
mkdir -p GSM7211257
mkdir -p GSM7211257/filtered_feature_bc_matrix
cd GSM7211257/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7211257/filtered_feature_bc_matrix --output_file GSM7211257/filtered_feature_bc_matrix.h5
mkdir -p GSM7211257/spatial
cd GSM7211257/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7211257
mkdir -p GSM7211257
mkdir -p GSM7211257/filtered_feature_bc_matrix
cd GSM7211257/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7211257/filtered_feature_bc_matrix --output_file GSM7211257/filtered_feature_bc_matrix.h5
mkdir -p GSM7211257/spatial
cd GSM7211257/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7211257&format=file&file=GSM7211257%5FEHU%2DW3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8162516
mkdir -p GSM8162516
mkdir -p GSM8162516/filtered_feature_bc_matrix
cd GSM8162516/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8162516/filtered_feature_bc_matrix --output_file GSM8162516/filtered_feature_bc_matrix.h5
mkdir -p GSM8162516/spatial
cd GSM8162516/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8162516&format=file&file=GSM8162516%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557976
mkdir -p GSM8557976
mkdir -p GSM8557976/filtered_feature_bc_matrix
cd GSM8557976/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557976/filtered_feature_bc_matrix --output_file GSM8557976/filtered_feature_bc_matrix.h5
mkdir -p GSM8557976/spatial
cd GSM8557976/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557976&format=file&file=GSM8557976%5FBPH%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557977
mkdir -p GSM8557977
mkdir -p GSM8557977/filtered_feature_bc_matrix
cd GSM8557977/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557977/filtered_feature_bc_matrix --output_file GSM8557977/filtered_feature_bc_matrix.h5
mkdir -p GSM8557977/spatial
cd GSM8557977/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557977&format=file&file=GSM8557977%5FBPH%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557978
mkdir -p GSM8557978
mkdir -p GSM8557978/filtered_feature_bc_matrix
cd GSM8557978/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557978/filtered_feature_bc_matrix --output_file GSM8557978/filtered_feature_bc_matrix.h5
mkdir -p GSM8557978/spatial
cd GSM8557978/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557978&format=file&file=GSM8557978%5FBPH%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557979
mkdir -p GSM8557979
mkdir -p GSM8557979/filtered_feature_bc_matrix
cd GSM8557979/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557979/filtered_feature_bc_matrix --output_file GSM8557979/filtered_feature_bc_matrix.h5
mkdir -p GSM8557979/spatial
cd GSM8557979/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557979&format=file&file=GSM8557979%5FBPH%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557980
mkdir -p GSM8557980
mkdir -p GSM8557980/filtered_feature_bc_matrix
cd GSM8557980/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557980/filtered_feature_bc_matrix --output_file GSM8557980/filtered_feature_bc_matrix.h5
mkdir -p GSM8557980/spatial
cd GSM8557980/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557980&format=file&file=GSM8557980%5FTRNA%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557981
mkdir -p GSM8557981
mkdir -p GSM8557981/filtered_feature_bc_matrix
cd GSM8557981/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557981/filtered_feature_bc_matrix --output_file GSM8557981/filtered_feature_bc_matrix.h5
mkdir -p GSM8557981/spatial
cd GSM8557981/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557981&format=file&file=GSM8557981%5FTRNA%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557982
mkdir -p GSM8557982
mkdir -p GSM8557982/filtered_feature_bc_matrix
cd GSM8557982/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557982/filtered_feature_bc_matrix --output_file GSM8557982/filtered_feature_bc_matrix.h5
mkdir -p GSM8557982/spatial
cd GSM8557982/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557982&format=file&file=GSM8557982%5FTRNA%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557983
mkdir -p GSM8557983
mkdir -p GSM8557983/filtered_feature_bc_matrix
cd GSM8557983/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557983/filtered_feature_bc_matrix --output_file GSM8557983/filtered_feature_bc_matrix.h5
mkdir -p GSM8557983/spatial
cd GSM8557983/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557983&format=file&file=GSM8557983%5FTRNA%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557984
mkdir -p GSM8557984
mkdir -p GSM8557984/filtered_feature_bc_matrix
cd GSM8557984/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557984/filtered_feature_bc_matrix --output_file GSM8557984/filtered_feature_bc_matrix.h5
mkdir -p GSM8557984/spatial
cd GSM8557984/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557984&format=file&file=GSM8557984%5FTRNA%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557985
mkdir -p GSM8557985
mkdir -p GSM8557985/filtered_feature_bc_matrix
cd GSM8557985/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557985/filtered_feature_bc_matrix --output_file GSM8557985/filtered_feature_bc_matrix.h5
mkdir -p GSM8557985/spatial
cd GSM8557985/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557985&format=file&file=GSM8557985%5FTRNA%5F6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557986
mkdir -p GSM8557986
mkdir -p GSM8557986/filtered_feature_bc_matrix
cd GSM8557986/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557986/filtered_feature_bc_matrix --output_file GSM8557986/filtered_feature_bc_matrix.h5
mkdir -p GSM8557986/spatial
cd GSM8557986/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557986&format=file&file=GSM8557986%5FTRNA%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557987
mkdir -p GSM8557987
mkdir -p GSM8557987/filtered_feature_bc_matrix
cd GSM8557987/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557987/filtered_feature_bc_matrix --output_file GSM8557987/filtered_feature_bc_matrix.h5
mkdir -p GSM8557987/spatial
cd GSM8557987/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557987&format=file&file=GSM8557987%5FTRNA%5F8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557988
mkdir -p GSM8557988
mkdir -p GSM8557988/filtered_feature_bc_matrix
cd GSM8557988/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557988/filtered_feature_bc_matrix --output_file GSM8557988/filtered_feature_bc_matrix.h5
mkdir -p GSM8557988/spatial
cd GSM8557988/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557988&format=file&file=GSM8557988%5FTRNA%5F9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557989
mkdir -p GSM8557989
mkdir -p GSM8557989/filtered_feature_bc_matrix
cd GSM8557989/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557989/filtered_feature_bc_matrix --output_file GSM8557989/filtered_feature_bc_matrix.h5
mkdir -p GSM8557989/spatial
cd GSM8557989/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557989&format=file&file=GSM8557989%5FTRNA%5F10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557990
mkdir -p GSM8557990
mkdir -p GSM8557990/filtered_feature_bc_matrix
cd GSM8557990/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557990/filtered_feature_bc_matrix --output_file GSM8557990/filtered_feature_bc_matrix.h5
mkdir -p GSM8557990/spatial
cd GSM8557990/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557990&format=file&file=GSM8557990%5FTRNA%5F11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557991
mkdir -p GSM8557991
mkdir -p GSM8557991/filtered_feature_bc_matrix
cd GSM8557991/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557991/filtered_feature_bc_matrix --output_file GSM8557991/filtered_feature_bc_matrix.h5
mkdir -p GSM8557991/spatial
cd GSM8557991/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557991&format=file&file=GSM8557991%5FTRNA%5F12%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557992
mkdir -p GSM8557992
mkdir -p GSM8557992/filtered_feature_bc_matrix
cd GSM8557992/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557992/filtered_feature_bc_matrix --output_file GSM8557992/filtered_feature_bc_matrix.h5
mkdir -p GSM8557992/spatial
cd GSM8557992/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557992&format=file&file=GSM8557992%5FTRNA%5F13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557993
mkdir -p GSM8557993
mkdir -p GSM8557993/filtered_feature_bc_matrix
cd GSM8557993/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557993/filtered_feature_bc_matrix --output_file GSM8557993/filtered_feature_bc_matrix.h5
mkdir -p GSM8557993/spatial
cd GSM8557993/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557993&format=file&file=GSM8557993%5FTRNA%5F14%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557994
mkdir -p GSM8557994
mkdir -p GSM8557994/filtered_feature_bc_matrix
cd GSM8557994/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557994/filtered_feature_bc_matrix --output_file GSM8557994/filtered_feature_bc_matrix.h5
mkdir -p GSM8557994/spatial
cd GSM8557994/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557994&format=file&file=GSM8557994%5FTRNA%5F15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557995
mkdir -p GSM8557995
mkdir -p GSM8557995/filtered_feature_bc_matrix
cd GSM8557995/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557995/filtered_feature_bc_matrix --output_file GSM8557995/filtered_feature_bc_matrix.h5
mkdir -p GSM8557995/spatial
cd GSM8557995/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557995&format=file&file=GSM8557995%5FTRNA%5F16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557996
mkdir -p GSM8557996
mkdir -p GSM8557996/filtered_feature_bc_matrix
cd GSM8557996/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557996/filtered_feature_bc_matrix --output_file GSM8557996/filtered_feature_bc_matrix.h5
mkdir -p GSM8557996/spatial
cd GSM8557996/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557996&format=file&file=GSM8557996%5FTRNA%5F17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557997
mkdir -p GSM8557997
mkdir -p GSM8557997/filtered_feature_bc_matrix
cd GSM8557997/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557997/filtered_feature_bc_matrix --output_file GSM8557997/filtered_feature_bc_matrix.h5
mkdir -p GSM8557997/spatial
cd GSM8557997/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557997&format=file&file=GSM8557997%5FNEADT%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557998
mkdir -p GSM8557998
mkdir -p GSM8557998/filtered_feature_bc_matrix
cd GSM8557998/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557998/filtered_feature_bc_matrix --output_file GSM8557998/filtered_feature_bc_matrix.h5
mkdir -p GSM8557998/spatial
cd GSM8557998/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557998&format=file&file=GSM8557998%5FNEADT%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8557999
mkdir -p GSM8557999
mkdir -p GSM8557999/filtered_feature_bc_matrix
cd GSM8557999/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8557999/filtered_feature_bc_matrix --output_file GSM8557999/filtered_feature_bc_matrix.h5
mkdir -p GSM8557999/spatial
cd GSM8557999/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8557999&format=file&file=GSM8557999%5FNEADT%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558000
mkdir -p GSM8558000
mkdir -p GSM8558000/filtered_feature_bc_matrix
cd GSM8558000/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558000/filtered_feature_bc_matrix --output_file GSM8558000/filtered_feature_bc_matrix.h5
mkdir -p GSM8558000/spatial
cd GSM8558000/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558000&format=file&file=GSM8558000%5FNEADT%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558001
mkdir -p GSM8558001
mkdir -p GSM8558001/filtered_feature_bc_matrix
cd GSM8558001/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558001/filtered_feature_bc_matrix --output_file GSM8558001/filtered_feature_bc_matrix.h5
mkdir -p GSM8558001/spatial
cd GSM8558001/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558001&format=file&file=GSM8558001%5FNEADT%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558002
mkdir -p GSM8558002
mkdir -p GSM8558002/filtered_feature_bc_matrix
cd GSM8558002/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558002/filtered_feature_bc_matrix --output_file GSM8558002/filtered_feature_bc_matrix.h5
mkdir -p GSM8558002/spatial
cd GSM8558002/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558002&format=file&file=GSM8558002%5FNEADT%5F6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558003
mkdir -p GSM8558003
mkdir -p GSM8558003/filtered_feature_bc_matrix
cd GSM8558003/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558003/filtered_feature_bc_matrix --output_file GSM8558003/filtered_feature_bc_matrix.h5
mkdir -p GSM8558003/spatial
cd GSM8558003/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558003&format=file&file=GSM8558003%5FNEADT%5F7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558004
mkdir -p GSM8558004
mkdir -p GSM8558004/filtered_feature_bc_matrix
cd GSM8558004/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558004/filtered_feature_bc_matrix --output_file GSM8558004/filtered_feature_bc_matrix.h5
mkdir -p GSM8558004/spatial
cd GSM8558004/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558004&format=file&file=GSM8558004%5FNEADT%5F8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558005
mkdir -p GSM8558005
mkdir -p GSM8558005/filtered_feature_bc_matrix
cd GSM8558005/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558005/filtered_feature_bc_matrix --output_file GSM8558005/filtered_feature_bc_matrix.h5
mkdir -p GSM8558005/spatial
cd GSM8558005/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558005&format=file&file=GSM8558005%5FNEADT%5F9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558006
mkdir -p GSM8558006
mkdir -p GSM8558006/filtered_feature_bc_matrix
cd GSM8558006/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558006/filtered_feature_bc_matrix --output_file GSM8558006/filtered_feature_bc_matrix.h5
mkdir -p GSM8558006/spatial
cd GSM8558006/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558006&format=file&file=GSM8558006%5FNEADT%5F10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558007
mkdir -p GSM8558007
mkdir -p GSM8558007/filtered_feature_bc_matrix
cd GSM8558007/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558007/filtered_feature_bc_matrix --output_file GSM8558007/filtered_feature_bc_matrix.h5
mkdir -p GSM8558007/spatial
cd GSM8558007/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558007&format=file&file=GSM8558007%5FNEADT%5F11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558008
mkdir -p GSM8558008
mkdir -p GSM8558008/filtered_feature_bc_matrix
cd GSM8558008/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558008/filtered_feature_bc_matrix --output_file GSM8558008/filtered_feature_bc_matrix.h5
mkdir -p GSM8558008/spatial
cd GSM8558008/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558008&format=file&file=GSM8558008%5FNEADT%5F12%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558009
mkdir -p GSM8558009
mkdir -p GSM8558009/filtered_feature_bc_matrix
cd GSM8558009/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558009/filtered_feature_bc_matrix --output_file GSM8558009/filtered_feature_bc_matrix.h5
mkdir -p GSM8558009/spatial
cd GSM8558009/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558009&format=file&file=GSM8558009%5FNEADT%5F13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558010
mkdir -p GSM8558010
mkdir -p GSM8558010/filtered_feature_bc_matrix
cd GSM8558010/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558010/filtered_feature_bc_matrix --output_file GSM8558010/filtered_feature_bc_matrix.h5
mkdir -p GSM8558010/spatial
cd GSM8558010/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558010&format=file&file=GSM8558010%5FNEADT%5F14%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558011
mkdir -p GSM8558011
mkdir -p GSM8558011/filtered_feature_bc_matrix
cd GSM8558011/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558011/filtered_feature_bc_matrix --output_file GSM8558011/filtered_feature_bc_matrix.h5
mkdir -p GSM8558011/spatial
cd GSM8558011/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558011&format=file&file=GSM8558011%5FNEADT%5F15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558012
mkdir -p GSM8558012
mkdir -p GSM8558012/filtered_feature_bc_matrix
cd GSM8558012/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558012/filtered_feature_bc_matrix --output_file GSM8558012/filtered_feature_bc_matrix.h5
mkdir -p GSM8558012/spatial
cd GSM8558012/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558012&format=file&file=GSM8558012%5FNEADT%5F16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558013
mkdir -p GSM8558013
mkdir -p GSM8558013/filtered_feature_bc_matrix
cd GSM8558013/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558013/filtered_feature_bc_matrix --output_file GSM8558013/filtered_feature_bc_matrix.h5
mkdir -p GSM8558013/spatial
cd GSM8558013/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558013&format=file&file=GSM8558013%5FNEADT%5F17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558014
mkdir -p GSM8558014
mkdir -p GSM8558014/filtered_feature_bc_matrix
cd GSM8558014/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558014/filtered_feature_bc_matrix --output_file GSM8558014/filtered_feature_bc_matrix.h5
mkdir -p GSM8558014/spatial
cd GSM8558014/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558014&format=file&file=GSM8558014%5FNEADT%5F18%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558015
mkdir -p GSM8558015
mkdir -p GSM8558015/filtered_feature_bc_matrix
cd GSM8558015/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558015/filtered_feature_bc_matrix --output_file GSM8558015/filtered_feature_bc_matrix.h5
mkdir -p GSM8558015/spatial
cd GSM8558015/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558015&format=file&file=GSM8558015%5FNEADT%5F19%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558016
mkdir -p GSM8558016
mkdir -p GSM8558016/filtered_feature_bc_matrix
cd GSM8558016/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558016/filtered_feature_bc_matrix --output_file GSM8558016/filtered_feature_bc_matrix.h5
mkdir -p GSM8558016/spatial
cd GSM8558016/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558016&format=file&file=GSM8558016%5FNEADT%5F20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558017
mkdir -p GSM8558017
mkdir -p GSM8558017/filtered_feature_bc_matrix
cd GSM8558017/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558017/filtered_feature_bc_matrix --output_file GSM8558017/filtered_feature_bc_matrix.h5
mkdir -p GSM8558017/spatial
cd GSM8558017/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558017&format=file&file=GSM8558017%5FNEADT%5F21%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558018
mkdir -p GSM8558018
mkdir -p GSM8558018/filtered_feature_bc_matrix
cd GSM8558018/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558018/filtered_feature_bc_matrix --output_file GSM8558018/filtered_feature_bc_matrix.h5
mkdir -p GSM8558018/spatial
cd GSM8558018/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558018&format=file&file=GSM8558018%5FNEADT%5F22%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558019
mkdir -p GSM8558019
mkdir -p GSM8558019/filtered_feature_bc_matrix
cd GSM8558019/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558019/filtered_feature_bc_matrix --output_file GSM8558019/filtered_feature_bc_matrix.h5
mkdir -p GSM8558019/spatial
cd GSM8558019/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558019&format=file&file=GSM8558019%5FCRPC%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558020
mkdir -p GSM8558020
mkdir -p GSM8558020/filtered_feature_bc_matrix
cd GSM8558020/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558020/filtered_feature_bc_matrix --output_file GSM8558020/filtered_feature_bc_matrix.h5
mkdir -p GSM8558020/spatial
cd GSM8558020/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558020&format=file&file=GSM8558020%5FCRPC%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558021
mkdir -p GSM8558021
mkdir -p GSM8558021/filtered_feature_bc_matrix
cd GSM8558021/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558021/filtered_feature_bc_matrix --output_file GSM8558021/filtered_feature_bc_matrix.h5
mkdir -p GSM8558021/spatial
cd GSM8558021/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558021&format=file&file=GSM8558021%5FCRPC%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558022
mkdir -p GSM8558022
mkdir -p GSM8558022/filtered_feature_bc_matrix
cd GSM8558022/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558022/filtered_feature_bc_matrix --output_file GSM8558022/filtered_feature_bc_matrix.h5
mkdir -p GSM8558022/spatial
cd GSM8558022/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558022&format=file&file=GSM8558022%5FCRPC%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558023
mkdir -p GSM8558023
mkdir -p GSM8558023/filtered_feature_bc_matrix
cd GSM8558023/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558023/filtered_feature_bc_matrix --output_file GSM8558023/filtered_feature_bc_matrix.h5
mkdir -p GSM8558023/spatial
cd GSM8558023/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558023&format=file&file=GSM8558023%5FCRPC%5F5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558024
mkdir -p GSM8558024
mkdir -p GSM8558024/filtered_feature_bc_matrix
cd GSM8558024/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558024/filtered_feature_bc_matrix --output_file GSM8558024/filtered_feature_bc_matrix.h5
mkdir -p GSM8558024/spatial
cd GSM8558024/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558024&format=file&file=GSM8558024%5FMET%5FA%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558025
mkdir -p GSM8558025
mkdir -p GSM8558025/filtered_feature_bc_matrix
cd GSM8558025/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558025/filtered_feature_bc_matrix --output_file GSM8558025/filtered_feature_bc_matrix.h5
mkdir -p GSM8558025/spatial
cd GSM8558025/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558025&format=file&file=GSM8558025%5FMET%5FB%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558026
mkdir -p GSM8558026
mkdir -p GSM8558026/filtered_feature_bc_matrix
cd GSM8558026/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558026/filtered_feature_bc_matrix --output_file GSM8558026/filtered_feature_bc_matrix.h5
mkdir -p GSM8558026/spatial
cd GSM8558026/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558026&format=file&file=GSM8558026%5FMET%5FC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8558027
mkdir -p GSM8558027
mkdir -p GSM8558027/filtered_feature_bc_matrix
cd GSM8558027/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8558027/filtered_feature_bc_matrix --output_file GSM8558027/filtered_feature_bc_matrix.h5
mkdir -p GSM8558027/spatial
cd GSM8558027/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8558027&format=file&file=GSM8558027%5FMET%5FD%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM5420753
mkdir -p GSM5420753
cd GSM5420753
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420753&format=file&file=GSM5420753%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6592048
mkdir -p GSM6592048
cd GSM6592048
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6592048&format=file&file=GSM6592048%5FM1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
mv tissue_hires_image.png.gz tissue_hires_image.png
mv scalefactors_json.json.gz scalefactors_json.json
mv tissue_positions_list.csv.gz tissue_positions_list.csv
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM7757970
mkdir -p GSM7757970
mkdir -p GSM7757970/filtered_feature_bc_matrix
cd GSM7757970/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM7757970/filtered_feature_bc_matrix \
--output_file GSM7757970/filtered_feature_bc_matrix.h5
mkdir -p GSM7757970/spatial
cd GSM7757970/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7757970&format=file&file=GSM7757970%5FR1-S1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7782699
mkdir -p GSM7782699
cd GSM7782699
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782699&format=file&file=GSM7782699%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
cd CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_Fresh_Frozen_Human_Breast_Cancer
mkdir -p CytAssist_Fresh_Frozen_Human_Breast_Cancer
cd CytAssist_Fresh_Frozen_Human_Breast_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# V1_Breast_Cancer_Block_A_Section_1
mkdir -p V1_Breast_Cancer_Block_A_Section_1
cd V1_Breast_Cancer_Block_A_Section_1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# V1_Breast_Cancer_Block_A_Section_2
mkdir -p V1_Breast_Cancer_Block_A_Section_2
cd V1_Breast_Cancer_Block_A_Section_2
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# V1_Human_Invasive_Ductal_Carcinoma
mkdir -p V1_Human_Invasive_Ductal_Carcinoma
cd V1_Human_Invasive_Ductal_Carcinoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/V1_Human_Invasive_Ductal_Carcinoma/V1_Human_Invasive_Ductal_Carcinoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Human_Skin_Melanoma
mkdir -p CytAssist_FFPE_Human_Skin_Melanoma
cd CytAssist_FFPE_Human_Skin_Melanoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Visium_FFPE_Human_Ovarian_Cancer
mkdir -p Visium_FFPE_Human_Ovarian_Cancer
cd Visium_FFPE_Human_Ovarian_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1
mkdir -p CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1
cd CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Protein_Expression_Human_Glioblastoma
mkdir -p CytAssist_FFPE_Protein_Expression_Human_Glioblastoma
cd CytAssist_FFPE_Protein_Expression_Human_Glioblastoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
cd CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/Script/extract_count_matrix.R

# CytAssist_11mm_FFPE_Human_Glioblastoma
mkdir -p CytAssist_11mm_FFPE_Human_Glioblastoma
cd CytAssist_11mm_FFPE_Human_Glioblastoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_11mm_FFPE_Human_Lung_Cancer
mkdir -p CytAssist_11mm_FFPE_Human_Lung_Cancer
cd CytAssist_11mm_FFPE_Human_Lung_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma
mkdir -p CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma
cd CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma
mkdir -p CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma
cd CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Visium_FFPE_Human_Cervical_Cancer
mkdir -p Visium_FFPE_Human_Cervical_Cancer
cd Visium_FFPE_Human_Cervical_Cancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Parent_Visium_Human_BreastCancer
mkdir -p Parent_Visium_Human_BreastCancer
cd Parent_Visium_Human_BreastCancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Targeted_Visium_Human_Glioblastoma_Pan_Cancer
mkdir -p Targeted_Visium_Human_Glioblastoma_Pan_Cancer
cd Targeted_Visium_Human_Glioblastoma_Pan_Cancer
mkdir -p spatial
mkdir -p filtered_feature_bc_matrix
wget -O filtered_feature_bc_matrix.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_Glioblastoma_Pan_Cancer/Targeted_Visium_Human_Glioblastoma_Pan_Cancer_filtered_feature_bc_matrix.tar.gz"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_Glioblastoma_Pan_Cancer/Targeted_Visium_Human_Glioblastoma_Pan_Cancer_spatial.tar.gz"
tar -xzf filtered_feature_bc_matrix.tar.gz -C filtered_feature_bc_matrix
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f filtered_feature_bc_matrix.tar.gz spatial.tar.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R \
  --input_dir /data_d/WSJ/SpatialMetsDB/RawData/Targeted_Visium_Human_Glioblastoma_Pan_Cancer/filtered_feature_bc_matrix \
  --output_file /data_d/WSJ/SpatialMetsDB/RawData/Targeted_Visium_Human_Glioblastoma_Pan_Cancer/filtered_feature_bc_matrix.h5

# Parent_Visium_Human_Glioblastoma
mkdir -p Parent_Visium_Human_Glioblastoma
cd Parent_Visium_Human_Glioblastoma
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Targeted_Visium_Human_OvarianCancer_Immunology
mkdir -p Targeted_Visium_Human_OvarianCancer_Immunology
cd Targeted_Visium_Human_OvarianCancer_Immunology
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_OvarianCancer_Immunology/Targeted_Visium_Human_OvarianCancer_Immunology_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_OvarianCancer_Immunology/Targeted_Visium_Human_OvarianCancer_Immunology_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Targeted_Visium_Human_OvarianCancer_Pan_Cancer
mkdir -p Targeted_Visium_Human_OvarianCancer_Pan_Cancer
cd Targeted_Visium_Human_OvarianCancer_Pan_Cancer
mkdir -p spatial
wget -O raw_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_OvarianCancer_Pan_Cancer/Targeted_Visium_Human_OvarianCancer_Pan_Cancer_raw_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_OvarianCancer_Pan_Cancer/Targeted_Visium_Human_OvarianCancer_Pan_Cancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# Parent_Visium_Human_OvarianCancer
mkdir -p Parent_Visium_Human_OvarianCancer
cd Parent_Visium_Human_OvarianCancer
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_OvarianCancer/Parent_Visium_Human_OvarianCancer_filtered_feature_bc_matrix.h5"
wget -O spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_OvarianCancer/Parent_Visium_Human_OvarianCancer_spatial.tar.gz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6612124
mkdir -p GSM6612124
mkdir -p GSM6612124/filtered_feature_bc_matrix
cd GSM6612124/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5F39738%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5F39738%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5F39738%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6612124/filtered_feature_bc_matrix \
--output_file GSM6612124/filtered_feature_bc_matrix.h5
mkdir -p GSM6612124/spatial
cd GSM6612124/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5FGS3%2DPlacebo%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5FGS3%2DPlacebo%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5FGS3%2DPlacebo%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612124&format=file&file=GSM6612124%5FGS3%2DPlacebo%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6612125
mkdir -p GSM6612125
mkdir -p GSM6612125/filtered_feature_bc_matrix
cd GSM6612125/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5F39739%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5F39739%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5F39739%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6612125/filtered_feature_bc_matrix \
--output_file GSM6612125/filtered_feature_bc_matrix.h5
mkdir -p GSM6612125/spatial
cd GSM6612125/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5FGS3%2DE2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5FGS3%2DE2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5FGS3%2DE2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612125&format=file&file=GSM6612125%5FGS3%2DE2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6612126
mkdir -p GSM6612126
mkdir -p GSM6612126/filtered_feature_bc_matrix
cd GSM6612126/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5F39740%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5F39740%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5F39740%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6612126/filtered_feature_bc_matrix \
--output_file GSM6612126/filtered_feature_bc_matrix.h5
mkdir -p GSM6612126/spatial
cd GSM6612126/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5FSC31%2DPlacebo%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5FSC31%2DPlacebo%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5FSC31%2DPlacebo%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612126&format=file&file=GSM6612126%5FSC31%2DPlacebo%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6612127
mkdir -p GSM6612127
mkdir -p GSM6612127/filtered_feature_bc_matrix
cd GSM6612127/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5F39741%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5F39741%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5F39741%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6612127/filtered_feature_bc_matrix \
--output_file GSM6612127/filtered_feature_bc_matrix.h5
mkdir -p GSM6612127/spatial
cd GSM6612127/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5FSC31%2DE2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5FSC31%2DE2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5FSC31%2DE2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6612127&format=file&file=GSM6612127%5FSC31%2DE2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6760693
mkdir -p GSM6760693
cd GSM6760693
wget -O GSM6760693_P1A_ER.tar.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760693&format=file&file=GSM6760693%5FP1A%5FER%2Etar%2Egz"
tar -xzf GSM6760693_P1A_ER.tar.gz
mv P1A_ER/filtered_feature_bc_matrix.h5 .
unzip -q P1A_ER/spatial.zip
rm -rf P1A_ER
rm -f GSM6760693_P1A_ER.tar.gz
cd ..

# GSM6760694
mkdir -p GSM6760694
cd GSM6760694
wget -O GSM6760694_P1B_ER.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760694&format=file&file=GSM6760694%5FP1B%5FER%2Etar%2Egz"
tar -xzf GSM6760694_P1B_ER.tar.gz
mv P1B_ER/filtered_feature_bc_matrix.h5 .
unzip -q P1B_ER/spatial.zip
rm -rf P1B_ER
rm -f GSM6760694_P1B_ER.tar.gz
cd ..

# GSM6760695
mkdir -p GSM6760695
cd GSM6760695
wget -O GSM6760695_P2A_TNBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760695&format=file&file=GSM6760695%5FP2A%5FTNBC%2Etar%2Egz"
tar -xzf GSM6760695_P2A_TNBC.tar.gz
mv P2A_TNBC/filtered_feature_bc_matrix.h5 .
unzip -q P2A_TNBC/spatial.zip
rm -rf P2A_TNBC
rm -f GSM6760695_P2A_TNBC.tar.gz
cd ..

# GSM6760696
mkdir -p GSM6760696
cd GSM6760696
wget -O GSM6760696_P2B_TNBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760696&format=file&file=GSM6760696%5FP2B%5FTNBC%2Etar%2Egz"
tar -xzf GSM6760696_P2B_TNBC.tar.gz
mv P2B_TNBC/filtered_feature_bc_matrix.h5 .
unzip -q P2B_TNBC/spatial.zip
rm -rf P2B_TNBC
rm -f GSM6760696_P2B_TNBC.tar.gz
cd ..

# GSM6760697
mkdir -p GSM6760697
cd GSM6760697
wget -O GSM6760697_P3A_MBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760697&format=file&file=GSM6760697%5FP3A%5FMBC%2Etar%2Egz"
tar -xzf GSM6760697_P3A_MBC.tar.gz
mv P3A_MBC/filtered_feature_bc_matrix.h5 .
unzip -q P3A_MBC/spatial.zip
rm -rf P3A_MBC
rm -f GSM6760697_P3A_MBC.tar.gz
cd ..

# GSM6760698
mkdir -p GSM6760698
cd GSM6760698
wget -O GSM6760698_P3B_MBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760698&format=file&file=GSM6760698%5FP3B%5FMBC%2Etar%2Egz"
tar -xzf GSM6760698_P3B_MBC.tar.gz
mv P3B_MBC/filtered_feature_bc_matrix.h5 .
unzip -q P3B_MBC/spatial.zip
rm -rf P3B_MBC
rm -f GSM6760698_P3B_MBC.tar.gz
cd ..

# GSM6760699
mkdir -p GSM6760699
cd GSM6760699
wget -O GSM6760699_P4A_MBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760699&format=file&file=GSM6760699%5FP4A%5FMBC%2Etar%2Egz"
tar -xzf GSM6760699_P4A_MBC.tar.gz
mv P4A_MBC/filtered_feature_bc_matrix.h5 .
unzip -q P4A_MBC/spatial.zip
rm -rf P4A_MBC
rm -f GSM6760699_P4A_MBC.tar.gz
cd ..

# GSM6760700
mkdir -p GSM6760700
cd GSM6760700
wget -O GSM6760700_P4B_MBC.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6760700&format=file&file=GSM6760700%5FP4B%5FMBC%2Etar%2Egz"
tar -xzf GSM6760700_P4B_MBC.tar.gz
mv P4B_MBC/filtered_feature_bc_matrix.h5 .
unzip -q P4B_MBC/spatial.zip
rm -rf P4B_MBC
rm -f GSM6760700_P4B_MBC.tar.gz
cd ..

# GSM7688224
mkdir -p GSM7688224
mkdir -p GSM7688224/filtered_feature_bc_matrix
cd GSM7688224/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7688224/filtered_feature_bc_matrix --output_file GSM7688224/filtered_feature_bc_matrix.h5
mkdir -p GSM7688224/spatial
cd GSM7688224/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688224&format=file&file=GSM7688224%5FTumor%2D897%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7688225
mkdir -p GSM7688225
mkdir -p GSM7688225/filtered_feature_bc_matrix
cd GSM7688225/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7688225/filtered_feature_bc_matrix --output_file GSM7688225/filtered_feature_bc_matrix.h5
mkdir -p GSM7688225/spatial
cd GSM7688225/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7688225&format=file&file=GSM7688225%5FTumor%2D899%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7996201
mkdir -p GSM7996201
mkdir -p GSM7996201/filtered_feature_bc_matrix
cd GSM7996201/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7996201/filtered_feature_bc_matrix --output_file GSM7996201/filtered_feature_bc_matrix.h5
mkdir -p GSM7996201/spatial
cd GSM7996201/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996201&format=file&file=GSM7996201%5FSXR%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7996202
mkdir -p GSM7996202
mkdir -p GSM7996202/filtered_feature_bc_matrix
cd GSM7996202/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7996202/filtered_feature_bc_matrix --output_file GSM7996202/filtered_feature_bc_matrix.h5
mkdir -p GSM7996202/spatial
cd GSM7996202/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996202&format=file&file=GSM7996202%5FSXR%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM7996203
mkdir -p GSM7996203
mkdir -p GSM7996203/filtered_feature_bc_matrix
cd GSM7996203/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7996203/filtered_feature_bc_matrix --output_file GSM7996203/filtered_feature_bc_matrix.h5
mkdir -p GSM7996203/spatial
cd GSM7996203/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996203&format=file&file=GSM7996203%5FYZL%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../


# GSM7996204
mkdir -p GSM7996204
mkdir -p GSM7996204/filtered_feature_bc_matrix
cd GSM7996204/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7996204/filtered_feature_bc_matrix --output_file GSM7996204/filtered_feature_bc_matrix.h5
mkdir -p GSM7996204/spatial
cd GSM7996204/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7996204&format=file&file=GSM7996204%5FYZL%5F2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8428408
mkdir -p GSM8428408
mkdir -p GSM8428408/filtered_feature_bc_matrix
cd GSM8428408/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F09%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F08%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F07%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8428408/filtered_feature_bc_matrix --output_file GSM8428408/filtered_feature_bc_matrix.h5
mkdir -p GSM8428408/spatial
cd GSM8428408/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F06%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F05%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F04%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428408&format=file&file=GSM8428408%5F03%2E9%5F627%5FTM%5FV11Y03%2D083%5FA%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8428410
mkdir -p GSM8428410
mkdir -p GSM8428410/filtered_feature_bc_matrix
cd GSM8428410/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F09%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F08%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F07%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8428410/filtered_feature_bc_matrix --output_file GSM8428410/filtered_feature_bc_matrix.h5
mkdir -p GSM8428410/spatial
cd GSM8428410/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F06%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F05%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F04%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8428410&format=file&file=GSM8428410%5F03%2E9%5F628%5FTM%5FV11Y03%2D083%5FC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../


############################################
# 1. HTAN Spatial data download
############################################

python DownloadSynapse.py



############################################
# 1. CROS Spatial data download
############################################

for /f %i in (synapseids.txt) do (

if not exist "%i_filtered_feature_bc_matrix.h5" (
curl -L ^
 -H "Referer: https://ngdc.cncb.ac.cn/crost/download" ^
 "https://ngdc.cncb.ac.cn/crost/api/download/files?name=download_page%%2F%i%%2Ffiltered_feature_bc_matrix.h5" ^
 -o "%i_filtered_feature_bc_matrix.h5"
)

if not exist "%i_spatial.zip" (
curl -L ^
 -H "Referer: https://ngdc.cncb.ac.cn/crost/download" ^
 "https://ngdc.cncb.ac.cn/crost/api/download/files?name=download_page%%2F%i%%2Fspatial.zip" ^
 -o "%i_spatial.zip"
)

)

############################################
# 1.10X Visium HD
############################################

#Visium_HD_6p5mm_Human_Colon_Cancer
wget https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Human_Colon_Cancer/Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz
tar -xzf Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz
mkdir Visium_HD_6p5mm_Human_Colon_Cancer
cp binned_outputs/square_016um/filtered_feature_bc_matrix.h5 Visium_HD_6p5mm_Human_Colon_Cancer
cp -r binned_outputs/square_016um/spatial/ Visium_HD_6p5mm_Human_Colon_Cancer
rm -r binned_outputs/
rm -f Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz

Rscript -e '
  output="Visium_HD_6p5mm_Human_Colon_Cancer/"
  library(arrow)
  tissue_pos <- read_parquet(paste0(output,"/spatial/tissue_positions.parquet"))
  write.table(tissue_pos,paste0(output,"/spatial/tissue_positions.csv"),col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
'

############################################
# 1.GeoMx
############################################

#GSE200563
mkdir -p ../RawData/GSE200563
wget -O ../RawData/GSE200563/Expression.txt "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE200563&format=file&file=GSE200563%5Fprocessed%5Fdata%2Etxt%2Egz"
head -n 1 ../RawData/GSE200563/Expression.txt | awk -F"\t" '{out="Gene"};for(i=2;i<=NF;i++){out=out"\t"$i};print out}' >  ../RawData/GSE200563/Expression1.txt
tail -n +2 ../RawData/GSE200563/Expression.txt >> ../RawData/GSE200563/Expression1.txt
mv ../RawData/GSE200563/Expression1.txt ../RawData/GSE200563/Expression.txt
echo -e "sample\tgroup" > ../RawData/GSE200563/clinical.txt
head -n 1 "../RawData/GSE200563/Expression.txt" | tr '\t' '\n' > ../RawData/GSE200563/clinical1.txt
dos2unix ../RawData/GSE200563/clinical1.txt
awk -F"\t" '(NR>1){if ($0 ~ /^LB/) {
        print $0 "\tBrainMets"
    } else if ($0 ~ /^L/) {
        print $0 "\tLungTumor"
    } else if ($0 ~ /^TIME-L/) {
        print $0 "\tLungImmuneTME"
    } else if ($0 ~ /^TIME-B/) {
        print $0 "\tBrainMetsImmuneTME"
    } else if ($0 ~ /^TBME/) {
        print $0 "\tBrainMetsTumorTME"
    } else if ($0 ~ /^mLN/) {
        print $0 "\tLymphNodeMets"
    } else if ($0 ~ /^BC/) {
        print $0 "\tBrainNormal"
    }}' ../RawData/GSE200563/clinical1.txt >> ../RawData/GSE200563/clinical.txt
rm -f ../RawData/GSE200563/clinical1.txt


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


# GSM6210024
mkdir -p GSM6210024
cd GSM6210024
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210024&format=file&file=GSM6210024%5FD15%5FIL7%5FVisium%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210024&format=file&file=GSM6210024%5FD15%5FIL7%5FVisium%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210024&format=file&file=GSM6210024%5FD15%5FIL7%5FVisium%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210024&format=file&file=GSM6210024%5FD15%5FIL7%5FVisium%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210024&format=file&file=GSM6210024%5FD15%5FIL7%5FVisium%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM6210023
mkdir -p GSM6210023
cd GSM6210023
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210023&format=file&file=GSM6210023%5FD15%5FPBS%5FVisium%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210023&format=file&file=GSM6210023%5FD15%5FPBS%5FVisium%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210023&format=file&file=GSM6210023%5FD15%5FPBS%5FVisium%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210023&format=file&file=GSM6210023%5FD15%5FPBS%5FVisium%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6210023&format=file&file=GSM6210023%5FD15%5FPBS%5FVisium%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM6505142
mkdir -p GSM6505142
cd GSM6505142
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505142&format=file&file=GSM6505142%5FMC38%5Fd10%5Fm5%5FC1%2DC2%5Fe1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505142&format=file&file=GSM6505142%5FMC38%5Fd10%5Fm5%5FC1%2DC2%5Fe1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505142&format=file&file=GSM6505142%5FMC38%5Fd10%5Fm5%5FC1%2DC2%5Fe1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505142&format=file&file=GSM6505142%5FMC38%5Fd10%5Fm5%5FC1%2DC2%5Fe1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505143
mkdir -p GSM6505143
cd GSM6505143
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505143&format=file&file=GSM6505143%5FMC38%5Fd10%5Fm5%5FD%5Fe2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505143&format=file&file=GSM6505143%5FMC38%5Fd10%5Fm5%5FD%5Fe2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505143&format=file&file=GSM6505143%5FMC38%5Fd10%5Fm5%5FD%5Fe2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505143&format=file&file=GSM6505143%5FMC38%5Fd10%5Fm5%5FD%5Fe2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505144
mkdir -p GSM6505144
cd GSM6505144
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505144&format=file&file=GSM6505144%5FMC38%5Fd16%5Fm6%5FD%5Fe2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505144&format=file&file=GSM6505144%5FMC38%5Fd16%5Fm6%5FD%5Fe2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505144&format=file&file=GSM6505144%5FMC38%5Fd16%5Fm6%5FD%5Fe2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505144&format=file&file=GSM6505144%5FMC38%5Fd16%5Fm6%5FD%5Fe2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505145
mkdir -p GSM6505145
cd GSM6505145
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505145&format=file&file=GSM6505145%5FMC38%5Fd16%5Fm7%5FD%5Fe3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505145&format=file&file=GSM6505145%5FMC38%5Fd16%5Fm7%5FD%5Fe3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505145&format=file&file=GSM6505145%5FMC38%5Fd16%5Fm7%5FD%5Fe3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505145&format=file&file=GSM6505145%5FMC38%5Fd16%5Fm7%5FD%5Fe3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505146
mkdir -p GSM6505146
cd GSM6505146
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505146&format=file&file=GSM6505146%5FMC38%5Fd16%5Fm8%5FD%5Fe3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505146&format=file&file=GSM6505146%5FMC38%5Fd16%5Fm8%5FD%5Fe3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505146&format=file&file=GSM6505146%5FMC38%5Fd16%5Fm8%5FD%5Fe3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505146&format=file&file=GSM6505146%5FMC38%5Fd16%5Fm8%5FD%5Fe3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6661400
mkdir -p GSM6661400
cd GSM6661400
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661400&format=file&file=GSM6661400%5FVisium4M1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6661400/filtered_feature_bc_matrix --output_file GSM6661400/filtered_feature_bc_matrix.h5

# GSM6661401
mkdir -p GSM6661401
cd GSM6661401
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6661401&format=file&file=GSM6661401%5FVisium6M1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6661401/filtered_feature_bc_matrix --output_file GSM6661401/filtered_feature_bc_matrix.h5

# GSM6859066
mkdir -p GSM6859066
cd GSM6859066
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859066&format=file&file=GSM6859066%5FSample3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859066&format=file&file=GSM6859066%5FSample3%5Fspatial%2Etar%2Egz"

tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859067
mkdir -p GSM6859067
cd GSM6859067
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859067&format=file&file=GSM6859067%5FSample4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859067&format=file&file=GSM6859067%5FSample4%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859068
mkdir -p GSM6859068
cd GSM6859068
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859068&format=file&file=GSM6859068%5FSample7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859068&format=file&file=GSM6859068%5FSample7%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859069
mkdir -p GSM6859069
cd GSM6859069
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859069&format=file&file=GSM6859069%5FSample10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859069&format=file&file=GSM6859069%5FSample10%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859070
mkdir -p GSM6859070
cd GSM6859070
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859070&format=file&file=GSM6859070%5FSample11%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859070&format=file&file=GSM6859070%5FSample11%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859071
mkdir -p GSM6859071
cd GSM6859071
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859071&format=file&file=GSM6859071%5FSample14%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859071&format=file&file=GSM6859071%5FSample14%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859072
mkdir -p GSM6859072
cd GSM6859072
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859072&format=file&file=GSM6859072%5FSample19%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859072&format=file&file=GSM6859072%5FSample19%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM6859073
mkdir -p GSM6859073
cd GSM6859073
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859073&format=file&file=GSM6859073%5FSample22%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6859073&format=file&file=GSM6859073%5FSample22%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

mkdir -p GSM7089855
cd GSM7089855
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7089855&format=file&file=GSM7089855%5FAjou%5FVisium%5FP1%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM7089856
mkdir -p GSM7089856
cd GSM7089856
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7089856&format=file&file=GSM7089856%5FAjou%5FVisium%5FP2%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM7089857
mkdir -p GSM7089857
cd GSM7089857
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7089857&format=file&file=GSM7089857%5FAjou%5FVisium%5FP3%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM7089858
mkdir -p GSM7089858
cd GSM7089858
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7089858&format=file&file=GSM7089858%5FAjou%5FVisium%5FP4%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# SN123_A798015_Rep1
mkdir SN123_A798015_Rep1
cd SN123_A798015_Rep1
wget -O SN123_A798015_Rep1.zip \
"https://zenodo.org/records/7760264/files/SN123_A798015_Rep1.zip?download=1"
unzip SN123_A798015_Rep1.zip
mv SN123_A798015_Rep1/* . 2>/dev/null || true
rm -rf SN123_A798015_Rep1
rm -f SN123_A798015_Rep1.zip
cd ..

# SN124_A798015_Rep2
mkdir SN124_A798015_Rep2
cd SN124_A798015_Rep2
wget -O SN124_A798015_Rep2.zip \
"https://zenodo.org/records/7760264/files/SN124_A798015_Rep2.zip?download=1"
unzip SN124_A798015_Rep2.zip
mv SN124_A798015_Rep2/* . 2>/dev/null || true
rm -rf SN124_A798015_Rep2
rm -f SN124_A798015_Rep2.zip
cd ..

# GSM8279108
mkdir -p GSM8279108
mkdir -p GSM8279108/filtered_feature_bc_matrix
cd GSM8279108/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8279108/filtered_feature_bc_matrix --output_file GSM8279108/filtered_feature_bc_matrix.h5
mkdir -p GSM8279108/spatial
cd GSM8279108/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279108&format=file&file=GSM8279108%5FGR4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8279109
mkdir -p GSM8279109
mkdir -p GSM8279109/filtered_feature_bc_matrix
cd GSM8279109/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8279109/filtered_feature_bc_matrix --output_file GSM8279109/filtered_feature_bc_matrix.h5
mkdir -p GSM8279109/spatial
cd GSM8279109/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279109&format=file&file=GSM8279109%5FGR4%5FPC%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8279110
mkdir -p GSM8279110
mkdir -p GSM8279110/filtered_feature_bc_matrix
cd GSM8279110/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8279110/filtered_feature_bc_matrix --output_file GSM8279110/filtered_feature_bc_matrix.h5
mkdir -p GSM8279110/spatial
cd GSM8279110/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279110&format=file&file=GSM8279110%5FGR7%5Fsite2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8279111
mkdir -p GSM8279111
mkdir -p GSM8279111/filtered_feature_bc_matrix
cd GSM8279111/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8279111/filtered_feature_bc_matrix --output_file GSM8279111/filtered_feature_bc_matrix.h5
mkdir -p GSM8279111/spatial
cd GSM8279111/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279111&format=file&file=GSM8279111%5FIC1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM8279112
mkdir -p GSM8279112
mkdir -p GSM8279112/filtered_feature_bc_matrix
cd GSM8279112/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM8279112/filtered_feature_bc_matrix --output_file GSM8279112/filtered_feature_bc_matrix.h5
mkdir -p GSM8279112/spatial
cd GSM8279112/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279112&format=file&file=GSM8279112%5FGR11%5Fsite1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../


# GSM5833528
mkdir -p GSM5833528
cd GSM5833528
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833528&format=file&file=GSM5833528%5FDMG1%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833529
mkdir -p GSM5833529
cd GSM5833529
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833529&format=file&file=GSM5833529%5FDMG2%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833530
mkdir -p GSM5833530
cd GSM5833530
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833530&format=file&file=GSM5833530%5FDMG3%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833531
mkdir -p GSM5833531
cd GSM5833531
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833531&format=file&file=GSM5833531%5FDMG4%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd .

# GSM5833532
mkdir -p GSM5833532
cd GSM5833532
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833532&format=file&file=GSM5833532%5FDMG5%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833533
mkdir -p GSM5833533
cd GSM5833533
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833533&format=file&file=GSM5833533%5FGBM1%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833534
mkdir -p GSM5833534
cd GSM5833534
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833534&format=file&file=GSM5833534%5FGBM2%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833535
mkdir -p GSM5833535
cd GSM5833535
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833535&format=file&file=GSM5833535%5FGBM3%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833536
mkdir -p GSM5833536
cd GSM5833536
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833536&format=file&file=GSM5833536%5FGBM4%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833537
mkdir -p GSM5833537
cd GSM5833537
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833537&format=file&file=GSM5833537%5FGBM5%5F1%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type  d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM5833538
mkdir -p GSM5833538
cd GSM5833538
wget -O data.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5833538&format=file&file=GSM5833538%5FGBM5%5F2%5Fspaceranger%5Fout%2Etar%2Egz"
tar -xzf data.tar.gz
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
rm -f data.tar.gz
cd ..

# GSM7507311
mkdir -p GSM7507311
cd GSM7507311
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507311&format=file&file=GSM7507311%5FST%5FPt1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507311&format=file&file=GSM7507311%5FPt1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507311&format=file&file=GSM7507311%5FPt1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507311&format=file&file=GSM7507311%5FPt1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507311&format=file&file=GSM7507311%5FPt1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
cd ../../

# GSM7507312
mkdir -p GSM7507312
cd GSM7507312
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507312&format=file&file=GSM7507312%5FST%5FPt2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507312&format=file&file=GSM7507312%5FPt2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507312&format=file&file=GSM7507312%5FPt2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507312&format=file&file=GSM7507312%5FPt2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507312&format=file&file=GSM7507312%5FPt2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507313
mkdir -p GSM7507313
cd GSM7507313
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507313&format=file&file=GSM7507313%5FST%5FPt3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507313&format=file&file=GSM7507313%5FPt3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507313&format=file&file=GSM7507313%5FPt3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507313&format=file&file=GSM7507313%5FPt3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507313&format=file&file=GSM7507313%5FPt3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507314
mkdir -p GSM7507314
cd GSM7507314
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507314&format=file&file=GSM7507314%5FST%5FPt4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507314&format=file&file=GSM7507314%5FPt4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507314&format=file&file=GSM7507314%5FPt4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507314&format=file&file=GSM7507314%5FPt4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507314&format=file&file=GSM7507314%5FPt4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507315
mkdir -p GSM7507315
cd GSM7507315
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507315&format=file&file=GSM7507315%5FST%5FPt5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507315&format=file&file=GSM7507315%5FPt5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507315&format=file&file=GSM7507315%5FPt5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507315&format=file&file=GSM7507315%5FPt5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507315&format=file&file=GSM7507315%5FPt5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507316
mkdir -p GSM7507316
cd GSM7507316
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507316&format=file&file=GSM7507316%5FST%5FPt5%5FPt7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507316&format=file&file=GSM7507316%5FPt5%5FPt7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507316&format=file&file=GSM7507316%5FPt5%5FPt7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507316&format=file&file=GSM7507316%5FPt5%5FPt7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507316&format=file&file=GSM7507316%5FPt5%5FPt7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507317
mkdir -p GSM7507317
cd GSM7507317
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507317&format=file&file=GSM7507317%5FST%5FPt6%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507317&format=file&file=GSM7507317%5FPt6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507317&format=file&file=GSM7507317%5FPt6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507317&format=file&file=GSM7507317%5FPt6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507317&format=file&file=GSM7507317%5FPt6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507318
mkdir -p GSM7507318
cd GSM7507318
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507318&format=file&file=GSM7507318%5FST%5FPt7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507318&format=file&file=GSM7507318%5FPt7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507318&format=file&file=GSM7507318%5FPt7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507318&format=file&file=GSM7507318%5FPt7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507318&format=file&file=GSM7507318%5FPt7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507319
mkdir -p GSM7507319
cd GSM7507319
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507319&format=file&file=GSM7507319%5FST%5FPt8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507319&format=file&file=GSM7507319%5FPt8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507319&format=file&file=GSM7507319%5FPt8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507319&format=file&file=GSM7507319%5FPt8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507319&format=file&file=GSM7507319%5FPt8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507320
mkdir -p GSM7507320
cd GSM7507320
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507320&format=file&file=GSM7507320%5FST%5FPt9%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507320&format=file&file=GSM7507320%5FPt9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507320&format=file&file=GSM7507320%5FPt9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507320&format=file&file=GSM7507320%5FPt9%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507320&format=file&file=GSM7507320%5FPt9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507321
mkdir -p GSM7507321
cd GSM7507321
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507321&format=file&file=GSM7507321%5FST%5FPt9%5FPt8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507321&format=file&file=GSM7507321%5FPt9%5FPt8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507321&format=file&file=GSM7507321%5FPt9%5FPt8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507321&format=file&file=GSM7507321%5FPt9%5FPt8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507321&format=file&file=GSM7507321%5FPt9%5FPt8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507322
mkdir -p GSM7507322
cd GSM7507322
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507322&format=file&file=GSM7507322%5FST%5FPt10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507322&format=file&file=GSM7507322%5FPt10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507322&format=file&file=GSM7507322%5FPt10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507322&format=file&file=GSM7507322%5FPt10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507322&format=file&file=GSM7507322%5FPt10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507323
mkdir -p GSM7507323
cd GSM7507323
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507323&format=file&file=GSM7507323%5FST%5FPt11%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507323&format=file&file=GSM7507323%5FPt11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507323&format=file&file=GSM7507323%5FPt11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507323&format=file&file=GSM7507323%5FPt11%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507323&format=file&file=GSM7507323%5FPt11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507324
mkdir -p GSM7507324
cd GSM7507324
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507324&format=file&file=GSM7507324%5FST%5FPt12%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507324&format=file&file=GSM7507324%5FPt12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507324&format=file&file=GSM7507324%5FPt12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507324&format=file&file=GSM7507324%5FPt12%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507324&format=file&file=GSM7507324%5FPt12%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507325
mkdir -p GSM7507325
cd GSM7507325
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507325&format=file&file=GSM7507325%5FST%5FPt13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507325&format=file&file=GSM7507325%5FPt13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507325&format=file&file=GSM7507325%5FPt13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507325&format=file&file=GSM7507325%5FPt13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507325&format=file&file=GSM7507325%5FPt13%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507326
mkdir -p GSM7507326
cd GSM7507326
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507326&format=file&file=GSM7507326%5FST%5FPt14%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507326&format=file&file=GSM7507326%5FPt14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507326&format=file&file=GSM7507326%5FPt14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507326&format=file&file=GSM7507326%5FPt14%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507326&format=file&file=GSM7507326%5FPt14%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7507327
mkdir -p GSM7507327
cd GSM7507327
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507327&format=file&file=GSM7507327%5FPt15%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7507327/filtered_feature_bc_matrix --output_file GSM7507327/filtered_feature_bc_matrix.h5

# GSM7507328
mkdir -p GSM7507328
cd GSM7507328
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507328&format=file&file=GSM7507328%5FPt16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7507328/filtered_feature_bc_matrix --output_file GSM7507328/filtered_feature_bc_matrix.h5

# GSM7507329
mkdir -p GSM7507329
cd GSM7507329
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507329&format=file&file=GSM7507329%5FPt17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7507329/filtered_feature_bc_matrix --output_file GSM7507329/filtered_feature_bc_matrix.h5

# GSM7507330
mkdir -p GSM7507330
cd GSM7507330
mkdir -p filtered_feature_bc_matrix
cd filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Fmatrix%2Emtx%2Egz"
cd ..
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7507330&format=file&file=GSM7507330%5FPt18%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM7507330/filtered_feature_bc_matrix --output_file GSM7507330/filtered_feature_bc_matrix.h5

#GSM7596587
mkdir -p GSM7596587
cd GSM7596587
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596587&format=file&file=GSM7596587%5Fmgh258%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596587&format=file&file=GSM7596587%5Fmgh258%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596587&format=file&file=GSM7596587%5Fmgh258%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596587&format=file&file=GSM7596587%5Fmgh258%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596588
mkdir -p GSM7596588
cd GSM7596588
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596588&format=file&file=GSM7596588%5Fzh881inf%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596588&format=file&file=GSM7596588%5Fzh881inf%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596588&format=file&file=GSM7596588%5Fzh881inf%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596588&format=file&file=GSM7596588%5Fzh881inf%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596589
mkdir -p GSM7596589
cd GSM7596589
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596589&format=file&file=GSM7596589%5Fzh881t1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596589&format=file&file=GSM7596589%5Fzh881t1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596589&format=file&file=GSM7596589%5Fzh881t1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596589&format=file&file=GSM7596589%5Fzh881t1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596590
mkdir -p GSM7596590
cd GSM7596590
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596590&format=file&file=GSM7596590%5Fzh916bulk%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596590&format=file&file=GSM7596590%5Fzh916bulk%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596590&format=file&file=GSM7596590%5Fzh916bulk%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596590&format=file&file=GSM7596590%5Fzh916bulk%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596591
mkdir -p GSM7596591
cd GSM7596591
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596591&format=file&file=GSM7596591%5Fzh916inf%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596591&format=file&file=GSM7596591%5Fzh916inf%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596591&format=file&file=GSM7596591%5Fzh916inf%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596591&format=file&file=GSM7596591%5Fzh916inf%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


#GSM7596592
mkdir -p GSM7596592
cd GSM7596592
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596592&format=file&file=GSM7596592%5Fzh916t1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596592&format=file&file=GSM7596592%5Fzh916t1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596592&format=file&file=GSM7596592%5Fzh916t1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596592&format=file&file=GSM7596592%5Fzh916t1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596593
mkdir -p GSM7596593
cd GSM7596593
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596593&format=file&file=GSM7596593%5Fzh1007inf%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596593&format=file&file=GSM7596593%5Fzh1007inf%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596593&format=file&file=GSM7596593%5Fzh1007inf%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596593&format=file&file=GSM7596593%5Fzh1007inf%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596594
mkdir -p GSM7596594
cd GSM7596594
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596594&format=file&file=GSM7596594%5Fzh1007nec%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596594&format=file&file=GSM7596594%5Fzh1007nec%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596594&format=file&file=GSM7596594%5Fzh1007nec%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596594&format=file&file=GSM7596594%5Fzh1007nec%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596595
mkdir -p GSM7596595
cd GSM7596595
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596595&format=file&file=GSM7596595%5Fzh1019inf%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596595&format=file&file=GSM7596595%5Fzh1019inf%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596595&format=file&file=GSM7596595%5Fzh1019inf%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596595&format=file&file=GSM7596595%5Fzh1019inf%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596596
mkdir -p GSM7596596
cd GSM7596596
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596596&format=file&file=GSM7596596%5Fzh1019t1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596596&format=file&file=GSM7596596%5Fzh1019t1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596596&format=file&file=GSM7596596%5Fzh1019t1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596596&format=file&file=GSM7596596%5Fzh1019t1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596597
mkdir -p GSM7596597
cd GSM7596597
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596597&format=file&file=GSM7596597%5Fzh8811a%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596597&format=file&file=GSM7596597%5Fzh8811a%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596597&format=file&file=GSM7596597%5Fzh8811a%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596597&format=file&file=GSM7596597%5Fzh8811a%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596598
mkdir -p GSM7596598
cd GSM7596598
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596598&format=file&file=GSM7596598%5Fzh8811b%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596598&format=file&file=GSM7596598%5Fzh8811b%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596598&format=file&file=GSM7596598%5Fzh8811b%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596598&format=file&file=GSM7596598%5Fzh8811b%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

mkdir -p GSM7596599
cd GSM7596599
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596599&format=file&file=GSM7596599%5Fzh8812%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596599&format=file&file=GSM7596599%5Fzh8812%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596599&format=file&file=GSM7596599%5Fzh8812%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596599&format=file&file=GSM7596599%5Fzh8812%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596600
mkdir -p GSM7596600
cd GSM7596600
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596600&format=file&file=GSM7596600%5Fbwh23%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596600&format=file&file=GSM7596600%5Fbwh23%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596600&format=file&file=GSM7596600%5Fbwh23%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596600&format=file&file=GSM7596600%5Fbwh23%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596601
mkdir -p GSM7596601
cd GSM7596601
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596601&format=file&file=GSM7596601%5Fbwh24%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596601&format=file&file=GSM7596601%5Fbwh24%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596601&format=file&file=GSM7596601%5Fbwh24%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596601&format=file&file=GSM7596601%5Fbwh24%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596602
mkdir -p GSM7596602
cd GSM7596602
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596602&format=file&file=GSM7596602%5Fbwh25%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596602&format=file&file=GSM7596602%5Fbwh25%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596602&format=file&file=GSM7596602%5Fbwh25%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596602&format=file&file=GSM7596602%5Fbwh25%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596603
mkdir -p GSM7596603
cd GSM7596603
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596603&format=file&file=GSM7596603%5Fbwh28%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596603&format=file&file=GSM7596603%5Fbwh28%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596603&format=file&file=GSM7596603%5Fbwh28%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596603&format=file&file=GSM7596603%5Fbwh28%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596604
mkdir -p GSM7596604
cd GSM7596604
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596604&format=file&file=GSM7596604%5Fbwh35%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596604&format=file&file=GSM7596604%5Fbwh35%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596604&format=file&file=GSM7596604%5Fbwh35%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596604&format=file&file=GSM7596604%5Fbwh35%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM7596605
mkdir -p GSM7596605
cd GSM7596605
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596605&format=file&file=GSM7596605%5Fmgh259%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596605&format=file&file=GSM7596605%5Fmgh259%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596605&format=file&file=GSM7596605%5Fmgh259%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7596605&format=file&file=GSM7596605%5Fmgh259%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_lowres_image.png tissue_hires_image.png
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

#GSM6176208
mkdir -p GSM6176208
cd GSM6176208
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176208&format=file&file=GSM6176208%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fp1%5Fvisium%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176208&format=file&file=GSM6176208%5Fscalefactors%5Fjson%5Fp1%5Fvisium%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176208&format=file&file=GSM6176208%5Ftissue%5Fhires%5Fimage%5Fp1%5Fvisium%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176208&format=file&file=GSM6176208%5Ftissue%5Fpositions%5Flist%5Fp1%5Fvisium%2Ecsv%2Egz"
gunzip -f scalefactors_json.json.gz
gunzip -f tissue_hires_image.png.gz
gunzip -f tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

#GSM6176213
mkdir -p GSM6176213
cd GSM6176213
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176213&format=file&file=GSM6176213%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fp2%5Fvisium%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176213&format=file&file=GSM6176213%5Fscalefactors%5Fjson%5Fp2%5Fvisium%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176213&format=file&file=GSM6176213%5Ftissue%5Fhires%5Fimage%5Fp2%5Fvisium%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176213&format=file&file=GSM6176213%5Ftissue%5Fpositions%5Flist%5Fp2%5Fvisium%2Ecsv%2Egz"
gunzip -f scalefactors_json.json.gz
gunzip -f tissue_hires_image.png.gz
gunzip -f tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

#GSM6176214
mkdir -p GSM6176214
cd GSM6176214
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176214&format=file&file=GSM6176214%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fp3%5Fvisium%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176214&format=file&file=GSM6176214%5Fscalefactors%5Fjson%5Fp3%5Fvisium%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176214&format=file&file=GSM6176214%5Ftissue%5Fhires%5Fimage%5Fp3%5Fvisium%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176214&format=file&file=GSM6176214%5Ftissue%5Fpositions%5Flist%5Fp3%5Fvisium%2Ecsv%2Egz"
gunzip -f scalefactors_json.json.gz
gunzip -f tissue_hires_image.png.gz
gunzip -f tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

#GSM6176215
mkdir -p GSM6176215
cd GSM6176215
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176215&format=file&file=GSM6176215%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%5Fp7%5Fvisium%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176215&format=file&file=GSM6176215%5Fscalefactors%5Fjson%5Fp7%5Fvisium%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176215&format=file&file=GSM6176215%5Ftissue%5Fhires%5Fimage%5Fp7%5Fvisium%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6176215&format=file&file=GSM6176215%5Ftissue%5Fpositions%5Flist%5Fp7%5Fvisium%2Ecsv%2Egz"
gunzip -f scalefactors_json.json.gz
gunzip -f tissue_hires_image.png.gz
gunzip -f tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM8324906
mkdir -p GSM8324906
cd GSM8324906
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324906&format=file&file=GSM8324906%5FHLS%5F78%5FS1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324906&format=file&file=GSM8324906%5FHLS%5F78%5FS1%5Fspaceranger%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ../..

# GSM8324907
mkdir -p GSM8324907
cd GSM8324907
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324907&format=file&file=GSM8324907%5FHLS%5F79%5FS2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324907&format=file&file=GSM8324907%5FHLS%5F79%5FS2%5Fspaceranger%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
find . -type f -path "*spatial/tissue_lowres_image.png" -exec mv {} spatial/ \;
find . -type f -path "*spatial/tissue_hires_image.png" -exec mv {} spatial/ \;
find . -type f -path "*spatial/scalefactors_json.json" -exec mv {} spatial/ \;
find . -type f -path "*spatial/tissue_positions_list.csv" -exec mv {} spatial/ \;
rm -rf HLS_* outs spatial/outs 2>/dev/null
rm -f spatial.tar.gz
cd ..

# GSM8324908
mkdir -p GSM8324908
cd GSM8324908
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324908&format=file&file=GSM8324908%5FHLS%5F80%5FS3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324908&format=file&file=GSM8324908%5FHLS%5F80%5FS3%5Fspaceranger%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz
find . -type f \( \
  -name "tissue_lowres_image.png" -o \
  -name "tissue_hires_image.png" -o \
  -name "scalefactors_json.json" -o \
  -name "tissue_positions_list.csv" \
\) -exec mv {} spatial/ \;
rm -rf HLS_* outs spatial/outs 2>/dev/null
rm -f spatial.tar.gz
cd ..

# GSM7422460
mkdir -p GSM7422460
cd GSM7422460
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7422460&format=file&file=GSM7422460%5F15wksDox%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7422460&format=file&file=GSM7422460%5F15wksDox%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
rm -f spatial.tar.gz
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
[ -f spatial/tissue_positions.csv ] && mv spatial/tissue_positions.csv spatial/tissue_positions_list.csv
cd ..

# GSM7422461
mkdir -p GSM7422461
cd GSM7422461
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7422461&format=file&file=GSM7422461%5F25wksDox%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7422461&format=file&file=GSM7422461%5F25wksDox%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
rm -f spatial.tar.gz
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
[ -f spatial/tissue_positions.csv ] && mv spatial/tissue_positions.csv spatial/tissue_positions_list.csv
cd ..



# GSM5808054
mkdir -p GSM5808054
cd GSM5808054
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808054&format=file&file=GSM5808054%5FKP%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808054&format=file&file=GSM5808054%5FKP%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808054&format=file&file=GSM5808054%5FKP%5F1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808054&format=file&file=GSM5808054%5FKP%5F1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808054&format=file&file=GSM5808054%5FKP%5F1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
cd ../../

# GSM5808056
mkdir -p GSM5808056
cd GSM5808056
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808056&format=file&file=GSM5808056%5FKP%5F3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808056&format=file&file=GSM5808056%5FKP%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808056&format=file&file=GSM5808056%5FKP%5F3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808056&format=file&file=GSM5808056%5FKP%5F3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808056&format=file&file=GSM5808056%5FKP%5F3%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
cd ../../

# GSM5808057
mkdir -p GSM5808057
cd GSM5808057
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808057&format=file&file=GSM5808057%5FKP%5F4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808057&format=file&file=GSM5808057%5FKP%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808057&format=file&file=GSM5808057%5FKP%5F4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808057&format=file&file=GSM5808057%5FKP%5F4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5808057&format=file&file=GSM5808057%5FKP%5F4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
cd ../../

# GSM6505136
mkdir -p GSM6505136
cd GSM6505136
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505136&format=file&file=GSM6505136%5FB16%5Fd10%5Fm1%5FC%5Fe2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505136&format=file&file=GSM6505136%5FB16%5Fd10%5Fm1%5FC%5Fe2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505136&format=file&file=GSM6505136%5FB16%5Fd10%5Fm1%5FC%5Fe2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505136&format=file&file=GSM6505136%5FB16%5Fd10%5Fm1%5FC%5Fe2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6505137
mkdir -p GSM6505137
cd GSM6505137
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505137&format=file&file=GSM6505137%5FB16%5Fd10%5Fm1%5FD1%2DD2%5Fe1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505137&format=file&file=GSM6505137%5FB16%5Fd10%5Fm1%5FD1%2DD2%5Fe1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505137&format=file&file=GSM6505137%5FB16%5Fd10%5Fm1%5FD1%2DD2%5Fe1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505137&format=file&file=GSM6505137%5FB16%5Fd10%5Fm1%5FD1%2DD2%5Fe1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../


# GSM6505138
mkdir -p GSM6505138
cd GSM6505138
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505138&format=file&file=GSM6505138%5FB16%5Fd16%5Fm2%5FC%5Fe2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505138&format=file&file=GSM6505138%5FB16%5Fd16%5Fm2%5FC%5Fe2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505138&format=file&file=GSM6505138%5FB16%5Fd16%5Fm2%5FC%5Fe2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505138&format=file&file=GSM6505138%5FB16%5Fd16%5Fm2%5FC%5Fe2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6505139
mkdir -p GSM6505139
cd GSM6505139
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505139&format=file&file=GSM6505139%5FB16%5Fd16%5Fm3%5FA%5Fe3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505139&format=file&file=GSM6505139%5FB16%5Fd16%5Fm3%5FA%5Fe3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505139&format=file&file=GSM6505139%5FB16%5Fd16%5Fm3%5FA%5Fe3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505139&format=file&file=GSM6505139%5FB16%5Fd16%5Fm3%5FA%5Fe3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6505140
mkdir -p GSM6505140
cd GSM6505140
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505140&format=file&file=GSM6505140%5FB16%5Fd16%5Fm3%5FC%5Fe3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505140&format=file&file=GSM6505140%5FB16%5Fd16%5Fm3%5FC%5Fe3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505140&format=file&file=GSM6505140%5FB16%5Fd16%5Fm3%5FC%5Fe3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505140&format=file&file=GSM6505140%5FB16%5Fd16%5Fm3%5FC%5Fe3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6505141
mkdir -p GSM6505141
cd GSM6505141
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505141&format=file&file=GSM6505141%5FB16%5Fd16%5Fm4%5FC%5Fe3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505141&format=file&file=GSM6505141%5FB16%5Fd16%5Fm4%5FC%5Fe3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505141&format=file&file=GSM6505141%5FB16%5Fd16%5Fm4%5FC%5Fe3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505141&format=file&file=GSM6505141%5FB16%5Fd16%5Fm4%5FC%5Fe3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7818789
mkdir -p GSM7818789
cd GSM7818789
wget -O GSM7818789.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7818789&format=file&file=GSM7818789%5FXV0001C1%2Etar%2Egz"
tar -xzf GSM7818789.tar.gz
SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | head -n 1 | sed 's|^\./||')
mv "$SAMPLE/filtered_feature_bc_matrix" .
mkdir -p spatial
find "$SAMPLE/spatial" -name "*tissue_positions*.csv" -exec cp {} spatial/tissue_positions_list.csv \;
find "$SAMPLE/spatial" -name "*scalefactors*.json" -exec cp {} spatial/scalefactors_json.json \;
find "$SAMPLE/spatial" -name "*hires*.png" -exec cp {} spatial/tissue_hires_image.png \;
find "$SAMPLE/spatial" -name "*lowres*.png" -exec cp {} spatial/tissue_lowres_image.png \;
cd filtered_feature_bc_matrix
mv *_barcodes.tsv.gz barcodes.tsv.gz
mv *_features.tsv.gz features.tsv.gz
mv *_matrix.mtx.gz matrix.mtx.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir filtered_feature_bc_matrix --output_file filtered_feature_bc_matrix.h5
rm -rf "$SAMPLE"
cd ..

# GSM7818790
mkdir -p GSM7818790
cd GSM7818790
wget -O GSM7818790.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7818790&format=file&file=GSM7818790%5FXV0001D1%2Etar%2Egz"
tar -xzf GSM7818790.tar.gz
SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | head -n 1 | sed 's|^\./||')
mv "$SAMPLE/filtered_feature_bc_matrix" .
mkdir -p spatial
find "$SAMPLE/spatial" -name "*tissue_positions*.csv" -exec cp {} spatial/tissue_positions_list.csv \;
find "$SAMPLE/spatial" -name "*scalefactors*.json" -exec cp {} spatial/scalefactors_json.json \;
find "$SAMPLE/spatial" -name "*hires*.png" -exec cp {} spatial/tissue_hires_image.png \;
find "$SAMPLE/spatial" -name "*lowres*.png" -exec cp {} spatial/tissue_lowres_image.png \;
cd filtered_feature_bc_matrix
[ -f barcodes.tsv.gz ] || mv *_barcodes.tsv.gz barcodes.tsv.gz
[ -f features.tsv.gz ] || mv *_features.tsv.gz features.tsv.gz
[ -f matrix.mtx.gz ] || mv *_matrix.mtx.gz matrix.mtx.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir filtered_feature_bc_matrix --output_file filtered_feature_bc_matrix.h5
rm -rf "$SAMPLE"
cd ..

# GSM7818791
mkdir -p GSM7818791
cd GSM7818791
wget -O GSM7818791.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7818791&format=file&file=GSM7818791%5FXV0008C1%2Etar%2Egz"
tar -xzf GSM7818791.tar.gz
SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | head -n 1 | sed 's|^\./||')
mv "$SAMPLE/filtered_feature_bc_matrix" .
mkdir -p spatial
find "$SAMPLE/spatial" -name "*tissue_positions*.csv" -exec cp {} spatial/tissue_positions_list.csv \;
find "$SAMPLE/spatial" -name "*scalefactors*.json" -exec cp {} spatial/scalefactors_json.json \;
find "$SAMPLE/spatial" -name "*hires*.png" -exec cp {} spatial/tissue_hires_image.png \;
find "$SAMPLE/spatial" -name "*lowres*.png" -exec cp {} spatial/tissue_lowres_image.png \;
cd filtered_feature_bc_matrix
[ -f barcodes.tsv.gz ] || mv *_barcodes.tsv.gz barcodes.tsv.gz
[ -f features.tsv.gz ] || mv *_features.tsv.gz features.tsv.gz
[ -f matrix.mtx.gz ] || mv *_matrix.mtx.gz matrix.mtx.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir filtered_feature_bc_matrix --output_file filtered_feature_bc_matrix.h5
rm -rf "$SAMPLE"
cd ..

# GSM7818792
mkdir -p GSM7818792
cd GSM7818792
wget -O GSM7818792.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7818792&format=file&file=GSM7818792%5FXV0008D1%2Etar%2Egz"
tar -xzf GSM7818792.tar.gz
SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | head -n 1 | sed 's|^\./||')
mv "$SAMPLE/filtered_feature_bc_matrix" .
mkdir -p spatial
find "$SAMPLE/spatial" -name "*tissue_positions*.csv" -exec cp {} spatial/tissue_positions_list.csv \;
find "$SAMPLE/spatial" -name "*scalefactors*.json" -exec cp {} spatial/scalefactors_json.json \;
find "$SAMPLE/spatial" -name "*hires*.png" -exec cp {} spatial/tissue_hires_image.png \;
find "$SAMPLE/spatial" -name "*lowres*.png" -exec cp {} spatial/tissue_lowres_image.png \;
cd filtered_feature_bc_matrix
[ -f barcodes.tsv.gz ] || mv *_barcodes.tsv.gz barcodes.tsv.gz
[ -f features.tsv.gz ] || mv *_features.tsv.gz features.tsv.gz
[ -f matrix.mtx.gz ] || mv *_matrix.mtx.gz matrix.mtx.gz
cd ..
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir filtered_feature_bc_matrix --output_file filtered_feature_bc_matrix.h5
rm -rf "$SAMPLE"
cd ..

# GSM7461403
mkdir -p GSM7461403
cd GSM7461403
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461403&format=file&file=GSM7461403%5FM28A26%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461403&format=file&file=GSM7461403%5FM28A26%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461403&format=file&file=GSM7461403%5FM28A26%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461403&format=file&file=GSM7461403%5FM28A26%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461403&format=file&file=GSM7461403%5FM28A26%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM7461404
mkdir -p GSM7461404
cd GSM7461404
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461404&format=file&file=GSM7461404%5FSFE65%5FNP137%5FT1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461404&format=file&file=GSM7461404%5FSFE65%5FNP137%5FT1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461404&format=file&file=GSM7461404%5FSFE65%5FNP137%5FT1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461404&format=file&file=GSM7461404%5FSFE65%5FNP137%5FT1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7461404&format=file&file=GSM7461404%5FSFE65%5FNP137%5FT1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM5643203
mkdir -p GSM5643203
cd GSM5643203
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5643203&format=file&file=GSM5643203%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5643203&format=file&file=GSM5643203%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5643203&format=file&file=GSM5643203%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5643203&format=file&file=GSM5643203%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6604700
mkdir -p GSM6604700
cd GSM6604700
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604700&format=file&file=GSM6604700%5FCtrl%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604700&format=file&file=GSM6604700%5FCtrl%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604700&format=file&file=GSM6604700%5FCtrl%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604700&format=file&file=GSM6604700%5FCtrl%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6604701
mkdir -p GSM6604701
cd GSM6604701
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604701&format=file&file=GSM6604701%5FRego%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604701&format=file&file=GSM6604701%5FRego%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604701&format=file&file=GSM6604701%5FRego%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6604701&format=file&file=GSM6604701%5FRego%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM6177612
mkdir -p GSM6177612
cd GSM6177612
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177612&format=file&file=GSM6177612%5FNYU%5FLIHC1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177612&format=file&file=GSM6177612%5FNYU%5FLIHC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177612&format=file&file=GSM6177612%5FNYU%5FLIHC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177612&format=file&file=GSM6177612%5FNYU%5FLIHC1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177612&format=file&file=GSM6177612%5FNYU%5FLIHC1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM7850822
mkdir -p GSM7850822
cd GSM7850822
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850822&format=file&file=GSM7850822%5FCHC20%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850822&format=file&file=GSM7850822%5FCHC20%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850822&format=file&file=GSM7850822%5FCHC20%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850822&format=file&file=GSM7850822%5FCHC20%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850822&format=file&file=GSM7850822%5FCHC20%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

#GSM7850823
mkdir -p GSM7850823
cd GSM7850823
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850823&format=file&file=GSM7850823%5FCHC23%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850823&format=file&file=GSM7850823%5FCHC23%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850823&format=file&file=GSM7850823%5FCHC23%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850823&format=file&file=GSM7850823%5FCHC23%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7850823&format=file&file=GSM7850823%5FCHC23%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

#GSM7968251
mkdir -p GSM7968251
mkdir -p GSM7968251/filtered_feature_bc_matrix
cd GSM7968251/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7968251/filtered_feature_bc_matrix --output_file GSM7968251/filtered_feature_bc_matrix.h5
mkdir -p GSM7968251/spatial
cd GSM7968251/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968251&format=file&file=GSM7968251%5FHB4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7968251
mkdir -p GSM7968252
mkdir -p GSM7968252/filtered_feature_bc_matrix
cd GSM7968252/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7968252/filtered_feature_bc_matrix --output_file GSM7968252/filtered_feature_bc_matrix.h5
mkdir -p GSM7968252/spatial
cd GSM7968252/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Ftissue%5Fpositions%2Ecsv%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7968252&format=file&file=GSM7968252%5FHB17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

mkdir -p GSM7476184
cd GSM7476184
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476184&format=file&file=GSM7476184%5FLAM1%5Ffilt%5Ffeat%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476184&format=file&file=GSM7476184%5FLAM1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476184&format=file&file=GSM7476184%5FLAM1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476184&format=file&file=GSM7476184%5FLAM1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476184&format=file&file=GSM7476184%5FLAM1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

mkdir -p GSM7476185
cd GSM7476185
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476185&format=file&file=GSM7476185%5FLAM2%5Ffilt%5Ffeat%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476185&format=file&file=GSM7476185%5FLAM2%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476185&format=file&file=GSM7476185%5FLAM2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476185&format=file&file=GSM7476185%5FLAM2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7476185&format=file&file=GSM7476185%5FLAM2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM5538852
mkdir -p GSM5538852
mkdir -p GSM5538852/filtered_feature_bc_matrix
cd GSM5538852/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538852&format=file&file=GSM5538852%5Fbarcodes%2Divuznh%5Fs8xyfv%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538852&format=file&file=GSM5538852%5Ffeatures%2Divuznh%5Fs8xyfv%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538852&format=file&file=GSM5538852%5Fmatrix%2Divuznh%5Fs8xyfv%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM5538852/filtered_feature_bc_matrix --output_file GSM5538852/filtered_feature_bc_matrix.h5
mkdir -p GSM5538852/spatial
cd GSM5538852/spatial
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538852&format=file&file=GSM5538852%5Fspatial%5Fivuznh%5Fs8xyfv%2Ezip"
unzip spatial.zip
rm -f spatial.zip
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
cd ../../

# GSM5538853
mkdir -p GSM5538853
mkdir -p GSM5538853/filtered_feature_bc_matrix
cd GSM5538853/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538853&format=file&file=GSM5538853%5Fbarcodes%2Dk886qc%5Frqix54%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538853&format=file&file=GSM5538853%5Ffeatures%2Dk886qc%5Frqix54%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538853&format=file&file=GSM5538853%5Fmatrix%2Dk886qc%5Frqix54%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM5538853/filtered_feature_bc_matrix --output_file GSM5538853/filtered_feature_bc_matrix.h5
mkdir -p GSM5538853/spatial
cd GSM5538853/spatial
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538853&format=file&file=GSM5538853%5Fspatial%2Dk886qc%5Frqix54%2Ezip"
unzip spatial.zip
rm -f spatial.zip
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
cd ../../

# GSM5538854
mkdir -p GSM5538854
mkdir -p GSM5538854/filtered_feature_bc_matrix
cd GSM5538854/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538854&format=file&file=GSM5538854%5Fbarcodes%2Ds2yq8o%5Ft5hk8u%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538854&format=file&file=GSM5538854%5Ffeatures%2Ds2yq8o%5Ft5hk8u%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538854&format=file&file=GSM5538854%5Fmatrix%2Ds2yq8o%5Ft5hk8u%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM5538854/filtered_feature_bc_matrix --output_file GSM5538854/filtered_feature_bc_matrix.h5
mkdir -p GSM5538854/spatial
cd GSM5538854/spatial
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538854&format=file&file=GSM5538854%5Fspatial%2Ds2yq8o%5Ft5hk8u%2Ezip"
unzip spatial.zip
rm -f spatial.zip
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
cd ../../

# GSM5538855
mkdir -p GSM5538855
mkdir -p GSM5538855/filtered_feature_bc_matrix
cd GSM5538855/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538855&format=file&file=GSM5538855%5Fbarcodes%2Duzntl1%5Ffkclxp%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538855&format=file&file=GSM5538855%5Ffeatures%2Duzntl1%5Ffkclxp%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538855&format=file&file=GSM5538855%5Fmatrix%2Duzntl1%5Ffkclxp%2Emtx%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM5538855/filtered_feature_bc_matrix --output_file GSM5538855/filtered_feature_bc_matrix.h5
mkdir -p GSM5538855/spatial
cd GSM5538855/spatial
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5538855&format=file&file=GSM5538855%5Fspatial%2Duzntl1%5Ffkclxp%2Ezip"
unzip spatial.zip
rm -f spatial.zip
if [ "$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find . -mindepth 1 -maxdepth 1 -type d); mv "$d"/* . && rmdir "$d"; fi
cd ../../

# GSM6177624
mkdir -p GSM6177624
cd GSM6177624
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177624&format=file&file=GSM6177624%5FPDX1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177624&format=file&file=GSM6177624%5FPDX1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177624&format=file&file=GSM6177624%5FPDX1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177624&format=file&file=GSM6177624%5FPDX1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177624&format=file&file=GSM6177624%5FPDX1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM6177625
mkdir -p GSM6177625
cd GSM6177625
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177625&format=file&file=GSM6177625%5FPDX2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177625&format=file&file=GSM6177625%5FPDX2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177625&format=file&file=GSM6177625%5FPDX2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177625&format=file&file=GSM6177625%5FPDX2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177625&format=file&file=GSM6177625%5FPDX2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../..

# GSM6248646
mkdir -p GSM6248646
mkdir -p GSM6248646/filtered_feature_bc_matrix
cd GSM6248646/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248646/filtered_feature_bc_matrix --output_file GSM6248646/filtered_feature_bc_matrix.h5
mkdir -p GSM6248646/spatial
cd GSM6248646/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248646&format=file&file=GSM6248646%5FNPC%5FST5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248647
mkdir -p GSM6248647
mkdir -p GSM6248647/filtered_feature_bc_matrix
cd GSM6248647/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248647/filtered_feature_bc_matrix --output_file GSM6248647/filtered_feature_bc_matrix.h5
mkdir -p GSM6248647/spatial
cd GSM6248647/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248647&format=file&file=GSM6248647%5FNPC%5FST6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248648
mkdir -p GSM6248648
mkdir -p GSM6248648/filtered_feature_bc_matrix
cd GSM6248648/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248648/filtered_feature_bc_matrix --output_file GSM6248648/filtered_feature_bc_matrix.h5
mkdir -p GSM6248648/spatial
cd GSM6248648/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248648&format=file&file=GSM6248648%5FNPC%5FST7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248649
mkdir -p GSM6248649
mkdir -p GSM6248649/filtered_feature_bc_matrix
cd GSM6248649/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248649/filtered_feature_bc_matrix --output_file GSM6248649/filtered_feature_bc_matrix.h5
mkdir -p GSM6248649/spatial
cd GSM6248649/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248649&format=file&file=GSM6248649%5FNPC%5FST8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248650
mkdir -p GSM6248650
mkdir -p GSM6248650/filtered_feature_bc_matrix
cd GSM6248650/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248650/filtered_feature_bc_matrix --output_file GSM6248650/filtered_feature_bc_matrix.h5
mkdir -p GSM6248650/spatial
cd GSM6248650/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248650&format=file&file=GSM6248650%5FNPC%5FST9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248651
mkdir -p GSM6248651
mkdir -p GSM6248651/filtered_feature_bc_matrix
cd GSM6248651/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248651/filtered_feature_bc_matrix --output_file GSM6248651/filtered_feature_bc_matrix.h5
mkdir -p GSM6248651/spatial
cd GSM6248651/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248651&format=file&file=GSM6248651%5FNPC%5FST10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248652
mkdir -p GSM6248652
mkdir -p GSM6248652/filtered_feature_bc_matrix
cd GSM6248652/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248652/filtered_feature_bc_matrix --output_file GSM6248652/filtered_feature_bc_matrix.h5
mkdir -p GSM6248652/spatial
cd GSM6248652/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248652&format=file&file=GSM6248652%5FNPC%5FST11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248653
mkdir -p GSM6248653
mkdir -p GSM6248653/filtered_feature_bc_matrix
cd GSM6248653/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248653/filtered_feature_bc_matrix --output_file GSM6248653/filtered_feature_bc_matrix.h5
mkdir -p GSM6248653/spatial
cd GSM6248653/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248653&format=file&file=GSM6248653%5FNPC%5FST12%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248654
mkdir -p GSM6248654
mkdir -p GSM6248654/filtered_feature_bc_matrix
cd GSM6248654/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248654/filtered_feature_bc_matrix --output_file GSM6248654/filtered_feature_bc_matrix.h5
mkdir -p GSM6248654/spatial
cd GSM6248654/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248654&format=file&file=GSM6248654%5FNPC%5FST16%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

# GSM6248655
mkdir -p GSM6248655
mkdir -p GSM6248655/filtered_feature_bc_matrix
cd GSM6248655/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript /data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R --input_dir GSM6248655/filtered_feature_bc_matrix --output_file GSM6248655/filtered_feature_bc_matrix.h5
mkdir -p GSM6248655/spatial
cd GSM6248655/spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248655&format=file&file=GSM6248655%5FNPC%5FST17%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#!/bin/bash
# Samples: cHC-1L cHC-1T HCC-1L HCC-1T HCC-2L HCC-2T HCC-3L HCC-3T HCC-4L HCC-4T HCC-5A HCC-5B HCC-5C HCC-5D ICC-1L
BASE_DIR="/data_d/WSJ/SpatialMetsDB/RawData"
URL_BASE="http://lifeome.net/supp/livercancer-st"
samples=(cHC-1L cHC-1T HCC-1L HCC-1T HCC-2L HCC-2T HCC-3L HCC-3T HCC-4L HCC-4T HCC-5A HCC-5B HCC-5C HCC-5D ICC-1L)
cd "$BASE_DIR"
for sample in "${samples[@]}"
do
    mkdir -p "$sample"
    cd "$sample" || exit
    wget -c "$URL_BASE/$sample.tar.gz"
    tar -xzf "$sample.tar.gz"
    mkdir -p spatial
    find . -type f -name "*.h5" ! -path "./spatial/*" | head -n1 | xargs -I{} cp "{}" filtered_feature_bc_matrix.h5
    find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
    find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
    find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
    find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
    find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} + 2>/dev/null
    cd "$BASE_DIR"
done

# GSM6339631
mkdir -p GSM6339631
cd GSM6339631
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339631&format=file&file=GSM6339631%5Fs1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339631&format=file&file=GSM6339631%5Fs1%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339631&format=file&file=GSM6339631%5Fs1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339631&format=file&file=GSM6339631%5Fs1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339632
mkdir -p GSM6339632
cd GSM6339632
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339632&format=file&file=GSM6339632%5Fs2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339632&format=file&file=GSM6339632%5Fs2%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339632&format=file&file=GSM6339632%5Fs2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339632&format=file&file=GSM6339632%5Fs2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339633
mkdir -p GSM6339633
cd GSM6339633
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339633&format=file&file=GSM6339633%5Fs3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339633&format=file&file=GSM6339633%5Fs3%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339633&format=file&file=GSM6339633%5Fs3%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339633&format=file&file=GSM6339633%5Fs3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339634
mkdir -p GSM6339634
cd GSM6339634
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339634&format=file&file=GSM6339634%5Fs4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339634&format=file&file=GSM6339634%5Fs4%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339634&format=file&file=GSM6339634%5Fs4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339634&format=file&file=GSM6339634%5Fs4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339635
mkdir -p GSM6339635
cd GSM6339635
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339635&format=file&file=GSM6339635%5Fs5%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339635&format=file&file=GSM6339635%5Fs5%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339635&format=file&file=GSM6339635%5Fs5%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339635&format=file&file=GSM6339635%5Fs5%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339636
mkdir -p GSM6339636
cd GSM6339636
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339636&format=file&file=GSM6339636%5Fs6%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339636&format=file&file=GSM6339636%5Fs6%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339636&format=file&file=GSM6339636%5Fs6%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339636&format=file&file=GSM6339636%5Fs6%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339637
mkdir -p GSM6339637
cd GSM6339637
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339637&format=file&file=GSM6339637%5Fs7%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339637&format=file&file=GSM6339637%5Fs7%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339637&format=file&file=GSM6339637%5Fs7%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339637&format=file&file=GSM6339637%5Fs7%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339638
mkdir -p GSM6339638
cd GSM6339638
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339638&format=file&file=GSM6339638%5Fs8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339638&format=file&file=GSM6339638%5Fs8%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339638&format=file&file=GSM6339638%5Fs8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339638&format=file&file=GSM6339638%5Fs8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339639
mkdir -p GSM6339639
cd GSM6339639
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339639&format=file&file=GSM6339639%5Fs9%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339639&format=file&file=GSM6339639%5Fs9%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339639&format=file&file=GSM6339639%5Fs9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339639&format=file&file=GSM6339639%5Fs9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339640
mkdir -p GSM6339640
cd GSM6339640
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339640&format=file&file=GSM6339640%5Fs10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339640&format=file&file=GSM6339640%5Fs10%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339640&format=file&file=GSM6339640%5Fs10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339640&format=file&file=GSM6339640%5Fs10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339641
mkdir -p GSM6339641
cd GSM6339641
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339641&format=file&file=GSM6339641%5Fs11%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339641&format=file&file=GSM6339641%5Fs11%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339641&format=file&file=GSM6339641%5Fs11%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339641&format=file&file=GSM6339641%5Fs11%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6339642
mkdir -p GSM6339642
cd GSM6339642
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339642&format=file&file=GSM6339642%5Fs12%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339642&format=file&file=GSM6339642%5Fs12%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339642&format=file&file=GSM6339642%5Fs12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6339642&format=file&file=GSM6339642%5Fs12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6833484
mkdir -p GSM6833484
mkdir -p GSM6833484/filtered_feature_bc_matrix
cd GSM6833484/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833484&format=file&file=GSM6833484%5FPatient1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833484&format=file&file=GSM6833484%5FPatient1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833484&format=file&file=GSM6833484%5FPatient1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6833484/filtered_feature_bc_matrix \
--output_file GSM6833484/filtered_feature_bc_matrix.h5
mkdir -p GSM6833484/spatial
cd GSM6833484
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833484&format=file&file=GSM6833484%5FPatient1%5Fimage%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then
    d=$(find spatial -mindepth 1 -maxdepth 1 -type d)
    mv "$d"/* spatial/
    rmdir "$d"
fi
rm -f spatial.tar.gz
cd spatial
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6833485
mkdir -p GSM6833485
mkdir -p GSM6833485/filtered_feature_bc_matrix
cd GSM6833485/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833485&format=file&file=GSM6833485%5FPatient2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833485&format=file&file=GSM6833485%5FPatient2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833485&format=file&file=GSM6833485%5FPatient2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6833485/filtered_feature_bc_matrix \
--output_file GSM6833485/filtered_feature_bc_matrix.h5
mkdir -p GSM6833485/spatial
cd GSM6833485
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833485&format=file&file=GSM6833485%5FPatient2%5Fimage%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then
    d=$(find spatial -mindepth 1 -maxdepth 1 -type d)
    mv "$d"/* spatial/
    rmdir "$d"
fi
rm -f spatial.tar.gz
cd spatial
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6833486
mkdir -p GSM6833486
mkdir -p GSM6833486/filtered_feature_bc_matrix
cd GSM6833486/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833486&format=file&file=GSM6833486%5FPatient3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833486&format=file&file=GSM6833486%5FPatient3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833486&format=file&file=GSM6833486%5FPatient3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6833486/filtered_feature_bc_matrix \
--output_file GSM6833486/filtered_feature_bc_matrix.h5
mkdir -p GSM6833486/spatial
cd GSM6833486
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833486&format=file&file=GSM6833486%5FPatient3%5Fimage%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then
    d=$(find spatial -mindepth 1 -maxdepth 1 -type d)
    mv "$d"/* spatial/
    rmdir "$d"
fi
rm -f spatial.tar.gz
cd spatial
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6833487
mkdir -p GSM6833487
mkdir -p GSM6833487/filtered_feature_bc_matrix
cd GSM6833487/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833487&format=file&file=GSM6833487%5FPatient4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833487&format=file&file=GSM6833487%5FPatient4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833487&format=file&file=GSM6833487%5FPatient4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" \
--input_dir GSM6833487/filtered_feature_bc_matrix \
--output_file GSM6833487/filtered_feature_bc_matrix.h5
mkdir -p GSM6833487/spatial
cd GSM6833487
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6833487&format=file&file=GSM6833487%5FPatient4%5Fimage%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then
    d=$(find spatial -mindepth 1 -maxdepth 1 -type d)
    mv "$d"/* spatial/
    rmdir "$d"
fi
rm -f spatial.tar.gz
cd spatial
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6177614
mkdir -p GSM6177614
cd GSM6177614
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177614&format=file&file=GSM6177614%5FNYU%5FOVCA1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177614&format=file&file=GSM6177614%5FNYU%5FOVCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177614&format=file&file=GSM6177614%5FNYU%5FOVCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177614&format=file&file=GSM6177614%5FNYU%5FOVCA1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177614&format=file&file=GSM6177614%5FNYU%5FOVCA1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6177617
mkdir -p GSM6177617
cd GSM6177617
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177617&format=file&file=GSM6177617%5FNYU%5FOVCA3%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177617&format=file&file=GSM6177617%5FNYU%5FOVCA3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177617&format=file&file=GSM6177617%5FNYU%5FOVCA3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_hires_image.png.gz
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177617&format=file&file=GSM6177617%5FNYU%5FOVCA3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177617&format=file&file=GSM6177617%5FNYU%5FOVCA3%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz
cd ../../

# GSM6506110
mkdir -p GSM6506110
mkdir -p GSM6506110/filtered_feature_bc_matrix
cd GSM6506110/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506110&format=file&file=GSM6506110%5FSP1%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506110&format=file&file=GSM6506110%5FSP1%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506110&format=file&file=GSM6506110%5FSP1%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506110/filtered_feature_bc_matrix --output_file GSM6506110/filtered_feature_bc_matrix.h5
mkdir -p GSM6506110/spatial
cd GSM6506110
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506110&format=file&file=GSM6506110%5FSP1%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP1 spatial.zip
find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} +
cd ..

# GSM6506111
mkdir -p GSM6506111
mkdir -p GSM6506111/filtered_feature_bc_matrix
cd GSM6506111/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506111&format=file&file=GSM6506111%5FSP2%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506111&format=file&file=GSM6506111%5FSP2%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506111&format=file&file=GSM6506111%5FSP2%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506111/filtered_feature_bc_matrix --output_file GSM6506111/filtered_feature_bc_matrix.h5
mkdir -p GSM6506111/spatial
cd GSM6506111
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506111&format=file&file=GSM6506111%5FSP2%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP2 spatial.zip
cd ..

# GSM6506112
mkdir -p GSM6506112
mkdir -p GSM6506112/filtered_feature_bc_matrix
cd GSM6506112/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506112&format=file&file=GSM6506112%5FSP3%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506112&format=file&file=GSM6506112%5FSP3%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506112&format=file&file=GSM6506112%5FSP3%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506112/filtered_feature_bc_matrix --output_file GSM6506112/filtered_feature_bc_matrix.h5
mkdir -p GSM6506112/spatial
cd GSM6506112
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506112&format=file&file=GSM6506112%5FSP3%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP3 spatial.zip
cd ..

# GSM6506113
mkdir -p GSM6506113
mkdir -p GSM6506113/filtered_feature_bc_matrix
cd GSM6506113/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506113&format=file&file=GSM6506113%5FSP4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506113&format=file&file=GSM6506113%5FSP4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506113&format=file&file=GSM6506113%5FSP4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506113/filtered_feature_bc_matrix --output_file GSM6506113/filtered_feature_bc_matrix.h5
mkdir -p GSM6506113/spatial
cd GSM6506113
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506113&format=file&file=GSM6506113%5FSP4%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP4 spatial.zip
cd ..

# GSM6506114
mkdir -p GSM6506114
mkdir -p GSM6506114/filtered_feature_bc_matrix
cd GSM6506114/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506114&format=file&file=GSM6506114%5FSP5%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506114&format=file&file=GSM6506114%5FSP5%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506114&format=file&file=GSM6506114%5FSP5%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506114/filtered_feature_bc_matrix --output_file GSM6506114/filtered_feature_bc_matrix.h5
mkdir -p GSM6506114/spatial
cd GSM6506114
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506114&format=file&file=GSM6506114%5FSP5%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP5 spatial.zip
find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} +
cd ..

# GSM6506115
mkdir -p GSM6506115
mkdir -p GSM6506115/filtered_feature_bc_matrix
cd GSM6506115/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506115&format=file&file=GSM6506115%5FSP6%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506115&format=file&file=GSM6506115%5FSP6%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506115&format=file&file=GSM6506115%5FSP6%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506115/filtered_feature_bc_matrix --output_file GSM6506115/filtered_feature_bc_matrix.h5
mkdir -p GSM6506115/spatial
cd GSM6506115
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506115&format=file&file=GSM6506115%5FSP6%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP6 spatial.zip
find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} +
cd ..

# GSM6506116
mkdir -p GSM6506116
mkdir -p GSM6506116/filtered_feature_bc_matrix
cd GSM6506116/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506116&format=file&file=GSM6506116%5FSP7%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506116&format=file&file=GSM6506116%5FSP7%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506116&format=file&file=GSM6506116%5FSP7%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506116/filtered_feature_bc_matrix --output_file GSM6506116/filtered_feature_bc_matrix.h5
mkdir -p GSM6506116/spatial
cd GSM6506116
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506116&format=file&file=GSM6506116%5FSP7%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP7 spatial.zip
find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} +
cd ..

# GSM6506117
mkdir -p GSM6506117
mkdir -p GSM6506117/filtered_feature_bc_matrix
cd GSM6506117/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506117&format=file&file=GSM6506117%5FSP8%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506117&format=file&file=GSM6506117%5FSP8%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506117&format=file&file=GSM6506117%5FSP8%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6506117/filtered_feature_bc_matrix --output_file GSM6506117/filtered_feature_bc_matrix.h5
mkdir -p GSM6506117/spatial
cd GSM6506117
wget -O spatial.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506117&format=file&file=GSM6506117%5FSP8%5Fspatial%2Ezip"
unzip -q spatial.zip
find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
rm -rf SP8 spatial.zip
find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} +
cd ..

# GSM7019835
mkdir -p GSM7019835
mkdir -p GSM7019835/filtered_feature_bc_matrix
cd GSM7019835/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7019835/filtered_feature_bc_matrix --output_file GSM7019835/filtered_feature_bc_matrix.h5
mkdir -p GSM7019835/spatial
cd GSM7019835/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019835&format=file&file=GSM7019835%5FOCCR2v%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
mv tissue_lowres_image.png tissue_hires_image.png 
jq '.tissue_hires_scalef = .tissue_lowres_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM7019836
mkdir -p GSM7019836
mkdir -p GSM7019836/filtered_feature_bc_matrix
cd GSM7019836/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7019836/filtered_feature_bc_matrix --output_file GSM7019836/filtered_feature_bc_matrix.h5
mkdir -p GSM7019836/spatial
cd GSM7019836/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7019836&format=file&file=GSM7019836%5FOCCS3v%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../../

# GSM7090083
mkdir -p GSM7090083
cd GSM7090083
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090083&format=file&file=GSM7090083%5FGenolabM%5FA1%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090083&format=file&file=GSM7090083%5FGenolabM%5FA1%2D1%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090083&format=file&file=GSM7090083%5FGenolabM%5FA1%2D1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090083&format=file&file=GSM7090083%5FGenolabM%5FA1%2D1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090083&format=file&file=GSM7090083%5FGenolabM%5FA1%2D1%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7090084
mkdir -p GSM7090084
cd GSM7090084
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090084&format=file&file=GSM7090084%5FGenolabM%5FA1%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090084&format=file&file=GSM7090084%5FGenolabM%5FA1%2D2%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090084&format=file&file=GSM7090084%5FGenolabM%5FA1%2D2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090084&format=file&file=GSM7090084%5FGenolabM%5FA1%2D2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090084&format=file&file=GSM7090084%5FGenolabM%5FA1%2D2%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7090085
mkdir -p GSM7090085
cd GSM7090085
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090085&format=file&file=GSM7090085%5FGenolabM%5FB1%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090085&format=file&file=GSM7090085%5FGenolabM%5FB1%2D2%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090085&format=file&file=GSM7090085%5FGenolabM%5FB1%2D2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090085&format=file&file=GSM7090085%5FGenolabM%5FB1%2D2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090085&format=file&file=GSM7090085%5FGenolabM%5FB1%2D2%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7090086
mkdir -p GSM7090086
cd GSM7090086
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090086&format=file&file=GSM7090086%5FNextSeq2000%5FA1%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090086&format=file&file=GSM7090086%5FNextSeq2000%5FA1%2D1%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090086&format=file&file=GSM7090086%5FNextSeq2000%5FA1%2D1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090086&format=file&file=GSM7090086%5FNextSeq2000%5FA1%2D1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090086&format=file&file=GSM7090086%5FNextSeq2000%5FA1%2D1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"gunzip scalefactors_json.json.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7090087
mkdir -p GSM7090087
cd GSM7090087
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090087&format=file&file=GSM7090087%5FNextSeq2000%5FA1%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090087&format=file&file=GSM7090087%5FNextSeq2000%5FA1%2D2%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090087&format=file&file=GSM7090087%5FNextSeq2000%5FA1%2D2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090087&format=file&file=GSM7090087%5FNextSeq2000%5FA1%2D2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090087&format=file&file=GSM7090087%5FNextSeq2000%5FA1%2D2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7090088
mkdir -p GSM7090088
cd GSM7090088
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090088&format=file&file=GSM7090088%5FNextSeq2000%5FB1%2D2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090088&format=file&file=GSM7090088%5FNextSeq2000%5FB1%2D2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090088&format=file&file=GSM7090088%5FNextSeq2000%5FB1%2D2%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090088&format=file&file=GSM7090088%5FNextSeq2000%5FB1%2D2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7090088&format=file&file=GSM7090088%5FNextSeq2000%5FB1%2D2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip scalefactors_json.json.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
cd ../..

# GSM6505133
mkdir -p GSM6505133
cd GSM6505133
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505133&format=file&file=GSM6505133%5FDonorA%5FFFPE%2Dprobes%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505133&format=file&file=GSM6505133%5FDonorA%5FFFPE%2Dprobes%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505133&format=file&file=GSM6505133%5FDonorA%5FFFPE%2Dprobes%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505133&format=file&file=GSM6505133%5FDonorA%5FFFPE%2Dprobes%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip -c tissue_hires_image.png.gz > tissue_hires_image.png
gunzip -c scalefactors_json.json.gz > scalefactors_json.json
gunzip -c tissue_positions_list.csv.gz > tissue_positions_list.csv
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505134
mkdir -p GSM6505134
cd GSM6505134
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505134&format=file&file=GSM6505134%5FDonorB%5FFFPE%2Dprobes%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505134&format=file&file=GSM6505134%5FDonorB%5FFFPE%2Dprobes%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505134&format=file&file=GSM6505134%5FDonorB%5FFFPE%2Dprobes%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505134&format=file&file=GSM6505134%5FDonorB%5FFFPE%2Dprobes%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip -c tissue_hires_image.png.gz > tissue_hires_image.png
gunzip -c scalefactors_json.json.gz > scalefactors_json.json
gunzip -c tissue_positions_list.csv.gz > tissue_positions_list.csv
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

# GSM6505135
mkdir -p GSM6505135
cd GSM6505135
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505135&format=file&file=GSM6505135%5FDonorC%5FFFPE%2Dprobes%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505135&format=file&file=GSM6505135%5FDonorC%5FFFPE%2Dprobes%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505135&format=file&file=GSM6505135%5FDonorC%5FFFPE%2Dprobes%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6505135&format=file&file=GSM6505135%5FDonorC%5FFFPE%2Dprobes%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip -c tissue_hires_image.png.gz > tissue_hires_image.png
gunzip -c scalefactors_json.json.gz > scalefactors_json.json
gunzip -c tissue_positions_list.csv.gz > tissue_positions_list.csv
cp tissue_hires_image.png tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' scalefactors_json.json > tmp.json && mv tmp.json scalefactors_json.json
cd ../..

#GSM7841727
mkdir -p GSM7841727
mkdir -p GSM7841727/filtered_feature_bc_matrix
cd GSM7841727/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841727/filtered_feature_bc_matrix --output_file GSM7841727/filtered_feature_bc_matrix.h5
mkdir -p GSM7841727/spatial
cd GSM7841727/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841727&format=file&file=GSM7841727%5FCDZ%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841728
mkdir -p GSM7841728
mkdir -p GSM7841728/filtered_feature_bc_matrix
cd GSM7841728/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841728/filtered_feature_bc_matrix --output_file GSM7841728/filtered_feature_bc_matrix.h5
mkdir -p GSM7841728/spatial
cd GSM7841728/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841728&format=file&file=GSM7841728%5FCFX%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841729
mkdir -p GSM7841729
mkdir -p GSM7841729/filtered_feature_bc_matrix
cd GSM7841729/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841729/filtered_feature_bc_matrix --output_file GSM7841729/filtered_feature_bc_matrix.h5
mkdir -p GSM7841729/spatial
cd GSM7841729/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841729&format=file&file=GSM7841729%5FHB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841730
mkdir -p GSM7841730
mkdir -p GSM7841730/filtered_feature_bc_matrix
cd GSM7841730/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841730/filtered_feature_bc_matrix --output_file GSM7841730/filtered_feature_bc_matrix.h5
mkdir -p GSM7841730/spatial
cd GSM7841730/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841730&format=file&file=GSM7841730%5FJLR%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841731
mkdir -p GSM7841731
mkdir -p GSM7841731/filtered_feature_bc_matrix
cd GSM7841731/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841731/filtered_feature_bc_matrix --output_file GSM7841731/filtered_feature_bc_matrix.h5
mkdir -p GSM7841731/spatial
cd GSM7841731/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841731&format=file&file=GSM7841731%5FLGL%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841732
mkdir -p GSM7841732
mkdir -p GSM7841732/filtered_feature_bc_matrix
cd GSM7841732/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841732/filtered_feature_bc_matrix --output_file GSM7841732/filtered_feature_bc_matrix.h5
mkdir -p GSM7841732/spatial
cd GSM7841732/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841732&format=file&file=GSM7841732%5FPRS%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841733
mkdir -p GSM7841733
mkdir -p GSM7841733/filtered_feature_bc_matrix
cd GSM7841733/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841733/filtered_feature_bc_matrix --output_file GSM7841733/filtered_feature_bc_matrix.h5
mkdir -p GSM7841733/spatial
cd GSM7841733/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841733&format=file&file=GSM7841733%5FRJZ%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM7841734
mkdir -p GSM7841734
mkdir -p GSM7841734/filtered_feature_bc_matrix
cd GSM7841734/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM7841734/filtered_feature_bc_matrix --output_file GSM7841734/filtered_feature_bc_matrix.h5
mkdir -p GSM7841734/spatial
cd GSM7841734/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7841734&format=file&file=GSM7841734%5FSLX%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM6534007
mkdir -p GSM6534007
mkdir -p GSM6534007/filtered_feature_bc_matrix
cd GSM6534007/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6534007/filtered_feature_bc_matrix --output_file GSM6534007/filtered_feature_bc_matrix.h5
mkdir -p GSM6534007/spatial
cd GSM6534007/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534007&format=file&file=GSM6534007%5FSRC106%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

#GSM6534008
mkdir -p GSM6534008
mkdir -p GSM6534008/filtered_feature_bc_matrix
cd GSM6534008/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6534008/filtered_feature_bc_matrix --output_file GSM6534008/filtered_feature_bc_matrix.h5
mkdir -p GSM6534008/spatial
cd GSM6534008/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534008&format=file&file=GSM6534008%5FSRC107%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

mkdir -p GSM6534010
mkdir -p GSM6534010/filtered_feature_bc_matrix
cd GSM6534010/filtered_feature_bc_matrix
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Fmatrix%2Emtx%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6534010/filtered_feature_bc_matrix --output_file GSM6534010/filtered_feature_bc_matrix.h5
mkdir -p GSM6534010/spatial
cd GSM6534010/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6534010&format=file&file=GSM6534010%5FSRC95%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
gunzip tissue_positions_list.csv.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_hires_image.png.gz
gunzip scalefactors_json.json.gz
cd ../../

mkdir -p GSM5420750
cd GSM5420750
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420750&format=file&file=GSM5420750%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5420750&format=file&file=GSM5420750%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..

# GSM7983359
mkdir -p GSM7983359
cd GSM7983359
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983359&format=file&file=GSM7983359%5Fsample4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983359&format=file&file=GSM7983359%5Fsample4%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983359&format=file&file=GSM7983359%5Fsample4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983359&format=file&file=GSM7983359%5Fsample4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983359&format=file&file=GSM7983359%5Fsample4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

# GSM7983366
mkdir -p GSM7983366
cd GSM7983366
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983366&format=file&file=GSM7983366%5Fsample16%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983366&format=file&file=GSM7983366%5Fsample16%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983366&format=file&file=GSM7983366%5Fsample16%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983366&format=file&file=GSM7983366%5Fsample16%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983366&format=file&file=GSM7983366%5Fsample16%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983358
mkdir -p GSM7983358
cd GSM7983358
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983358&format=file&file=GSM7983358%5Fsample2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983358&format=file&file=GSM7983358%5Fsample2%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983358&format=file&file=GSM7983358%5Fsample2%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983358&format=file&file=GSM7983358%5Fsample2%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983358&format=file&file=GSM7983358%5Fsample2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983361
mkdir -p GSM7983361
cd GSM7983361
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983361&format=file&file=GSM7983361%5Fsample8%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983361&format=file&file=GSM7983361%5Fsample8%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983361&format=file&file=GSM7983361%5Fsample8%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983361&format=file&file=GSM7983361%5Fsample8%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983361&format=file&file=GSM7983361%5Fsample8%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983362
mkdir -p GSM7983362
cd GSM7983362
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983362&format=file&file=GSM7983362%5Fsample12%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983362&format=file&file=GSM7983362%5Fsample12%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983362&format=file&file=GSM7983362%5Fsample12%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983362&format=file&file=GSM7983362%5Fsample12%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983362&format=file&file=GSM7983362%5Fsample12%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983363
mkdir -p GSM7983363
cd GSM7983363
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983363&format=file&file=GSM7983363%5Fsample13%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983363&format=file&file=GSM7983363%5Fsample13%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983363&format=file&file=GSM7983363%5Fsample13%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983363&format=file&file=GSM7983363%5Fsample13%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983363&format=file&file=GSM7983363%5Fsample13%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983364
mkdir -p GSM7983364
cd GSM7983364
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983364&format=file&file=GSM7983364%5Fsample14%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983364&format=file&file=GSM7983364%5Fsample14%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983364&format=file&file=GSM7983364%5Fsample14%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983364&format=file&file=GSM7983364%5Fsample14%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983364&format=file&file=GSM7983364%5Fsample14%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM7983365
mkdir -p GSM7983365
cd GSM7983365
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983365&format=file&file=GSM7983365%5Fsample15%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983365&format=file&file=GSM7983365%5Fsample15%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983365&format=file&file=GSM7983365%5Fsample15%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983365&format=file&file=GSM7983365%5Fsample15%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7983365&format=file&file=GSM7983365%5Fsample15%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM8376639
mkdir -p GSM8376639
cd GSM8376639
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376639&format=file&file=GSM8376639%5Fsample1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376639&format=file&file=GSM8376639%5Fsample1%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376639&format=file&file=GSM8376639%5Fsample1%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376639&format=file&file=GSM8376639%5Fsample1%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376639&format=file&file=GSM8376639%5Fsample1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM8376640
mkdir -p GSM8376640
cd GSM8376640
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376640&format=file&file=GSM8376640%5Fsample9%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376640&format=file&file=GSM8376640%5Fsample9%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376640&format=file&file=GSM8376640%5Fsample9%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376640&format=file&file=GSM8376640%5Fsample9%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376640&format=file&file=GSM8376640%5Fsample9%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM8376641
mkdir -p GSM8376641
cd GSM8376641
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376641&format=file&file=GSM8376641%5Fsample10%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376641&format=file&file=GSM8376641%5Fsample10%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376641&format=file&file=GSM8376641%5Fsample10%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376641&format=file&file=GSM8376641%5Fsample10%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8376641&format=file&file=GSM8376641%5Fsample10%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM6177607
mkdir -p GSM6177607
cd GSM6177607
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177607&format=file&file=GSM6177607%5FNYU%5FGIST1%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177607&format=file&file=GSM6177607%5FNYU%5FGIST1%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177607&format=file&file=GSM6177607%5FNYU%5FGIST1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177607&format=file&file=GSM6177607%5FNYU%5FGIST1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177607&format=file&file=GSM6177607%5FNYU%5FGIST1%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../..

#GSM6177609
mkdir -p GSM6177609
cd GSM6177609
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177609&format=file&file=GSM6177609%5FNYU%5FGIST2%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177609&format=file&file=GSM6177609%5FNYU%5FGIST2%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177609&format=file&file=GSM6177609%5FNYU%5FGIST2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177609&format=file&file=GSM6177609%5FNYU%5FGIST2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177609&format=file&file=GSM6177609%5FNYU%5FGIST2%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ..

#GSM7847146
mkdir -p GSM7847146
cd GSM7847146
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7847146&format=file&file=GSM7847146%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7847146&format=file&file=GSM7847146%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..


#GSM7990473
mkdir -p GSM7990473
cd GSM7990473
wget -O GSM7990473_20_00331_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990473&format=file&file=GSM7990473%5F20%5F00331%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990473_20_00331_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990473_20_00331_LI_SING.tar.gz
cd ..

#GSM7990474
mkdir -p GSM7990474
cd GSM7990474
wget -O GSM7990474_21_00731_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990474&format=file&file=GSM7990474%5F21%5F00731%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990474_21_00731_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990474_21_00731_LI_SING.tar.gz
cd ..

#GSM7990475
mkdir -p GSM7990475
cd GSM7990475
wget -O GSM7990475_21_00732_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990475&format=file&file=GSM7990475%5F21%5F00732%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990475_21_00732_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990475_21_00732_LI_SING.tar.gz
cd ..

#GSM7990476
mkdir -p GSM7990476
cd GSM7990476
wget -O GSM7990476_21_00733_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990476&format=file&file=GSM7990476%5F21%5F00733%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990476_21_00733_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990476_21_00733_LI_SING.tar.gz
cd ..

#GSM7990477
mkdir -p GSM7990477
cd GSM7990477
wget -O GSM7990477_21_00734_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990477&format=file&file=GSM7990477%5F21%5F00734%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990477_21_00734_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990477_21_00734_LI_SING.tar.gz
cd ..

#GSM7990478
mkdir -p GSM7990478
cd GSM7990478
wget -O GSM7990478_21_01251_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990478&format=file&file=GSM7990478%5F21%5F01251%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990478_21_01251_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990478_21_01251_LI_SING.tar.gz
cd ..

#GSM7990479
mkdir -p GSM7990479
cd GSM7990479
wget -O GSM7990479_21_01252_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990479&format=file&file=GSM7990479%5F21%5F01252%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990479_21_01252_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990479_21_01252_LI_SING.tar.gz
cd ..

#GSM7990480
mkdir -p GSM7990480
cd GSM7990480
wget -O GSM7990480_21_01253_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990480&format=file&file=GSM7990480%5F21%5F01253%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990480_21_01253_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990480_21_01253_LI_SING.tar.gz
cd ..

#GSM7990481
mkdir -p GSM7990481
cd GSM7990481
wget -O GSM7990481_21_01254_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990481&format=file&file=GSM7990481%5F21%5F01254%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990481_21_01254_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990481_21_01254_LI_SING.tar.gz
cd ..

#GSM7990482
mkdir -p GSM7990482
cd GSM7990482
wget -O GSM7990482_21_01675_LI_SING.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7990482&format=file&file=GSM7990482%5F21%5F01675%5FLI%5FSING%2Etar%2Egz"
tar -xzf GSM7990482_21_01675_LI_SING.tar.gz
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7990482_21_01675_LI_SING.tar.gz
cd ..

#GSM6177623
mkdir -p GSM6177623
cd GSM6177623
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177623&format=file&file=GSM6177623%5FNYU%5FUCEC3%5FVis%5Fprocessed%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
mkdir -p spatial
cd spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177623&format=file&file=GSM6177623%5FNYU%5FUCEC3%5FVis%5Fprocessed%5Fspatial%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177623&format=file&file=GSM6177623%5FNYU%5FUCEC3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177623&format=file&file=GSM6177623%5FNYU%5FUCEC3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6177623&format=file&file=GSM6177623%5FNYU%5FUCEC3%5FVis%5Fprocessed%5Fspatial%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ..

#GSM7054270
mkdir -p GSM7054270
cd GSM7054270
wget -O GSM7054270_01_034_C1d1.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7054270&format=file&file=GSM7054270%5F01%5F034%5FC1d1%2Ezip"
unzip GSM7054270_01_034_C1d1.zip
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7054270_01_034_C1d1.zip
cd ..

#GSM7054271
mkdir -p GSM7054271
cd GSM7054271
wget -O GSM7054271_01_034_C3d1.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7054271&format=file&file=GSM7054271%5F01%5F034%5FC3d1%2Ezip"
unzip GSM7054271_01_034_C3d1.zip
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7054271_01_034_C3d1.zip
cd ..

#GSM7054272
mkdir -p GSM7054272
cd GSM7054272
wget -O GSM7054272_01_039_C1d1.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7054272&format=file&file=GSM7054272%5F01%5F039%5FC1d1%2Ezip"
unzip GSM7054272_01_039_C1d1.zip
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7054272_01_039_C1d1.zip
cd ..

#GSM7054273
mkdir -p GSM7054273
cd GSM7054273
wget -O GSM7054273_01_039_C3d1.zip "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7054273&format=file&file=GSM7054273%5F01%5F039%5FC3d1%2Ezip"
unzip GSM7054273_01_039_C3d1.zip
find . -name "*filtered_feature_bc_matrix.h5" -exec mv {} ./filtered_feature_bc_matrix.h5 \;
mkdir -p spatial
find . -name "*tissue_hires_image.png" -exec mv {} spatial/tissue_hires_image.png \;
find . -name "*tissue_lowres_image.png" -exec mv {} spatial/tissue_lowres_image.png \;
find . -name "*tissue_positions_list.csv" -exec mv {} spatial/tissue_positions_list.csv \;
find . -name "*scalefactors_json.json*" -exec mv {} spatial/scalefactors_json.json \;
rm -rf $(ls | grep -v "filtered_feature_bc_matrix.h5" | grep -v spatial)
rm -f GSM7054273_01_039_C3d1.zip
cd ..




#!/bin/bash
# Samples: HCC-1L HCC-1T  HCC-5A  ICC-1L
BASE_DIR="/data_d/WSJ/SpatialMetsDB/RawData"
URL_BASE="http://lifeome.net/supp/livercancer-st"
samples=( HCC-1L HCC-1T HCC-5A HCC-5B HCC-5C )
cd "$BASE_DIR"
for sample in "${samples[@]}"
do
    mkdir -p "$sample"
    cd "$sample" || exit
    wget -c "$URL_BASE/$sample.tar.gz"
    tar -xzf "$sample.tar.gz"
    mkdir -p spatial
    find . -type f -name "*.h5" ! -path "./spatial/*" | head -n1 | xargs -I{} cp "{}" filtered_feature_bc_matrix.h5
    find . -type f -iname "*hires*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_hires_image.png
    find . -type f -iname "*lowres*.png" | head -n1 | xargs -I{} cp "{}" spatial/tissue_lowres_image.png
    find . -type f \( -iname "*tissue_positions*.csv" -o -iname "*positions*.csv" \) | head -n1 | xargs -I{} cp "{}" spatial/tissue_positions_list.csv
    find . -type f -iname "*scalefactors*.json" | head -n1 | xargs -I{} cp "{}" spatial/scalefactors_json.json
    find . -mindepth 1 ! -path "./filtered_feature_bc_matrix.h5" ! -path "./spatial" ! -path "./spatial/*" -exec rm -rf {} + 2>/dev/null
    cd "$BASE_DIR"
done



mkdir -p GSM8279107
mkdir -p GSM8279107/filtered_feature_bc_matrix
cd GSM8279107/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM8279107/filtered_feature_bc_matrix --output_file GSM8279107/filtered_feature_bc_matrix.h5
mkdir -p GSM8279107/spatial
cd GSM8279107/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8279107&format=file&file=GSM8279107%5FGR2%5Fsite4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../../

mkdir -p GSM6248656
mkdir -p GSM6248656/filtered_feature_bc_matrix
cd GSM6248656/filtered_feature_bc_matrix
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Fmatrix%2Emtx%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Ffeatures%2Etsv%2Egz"
wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Fbarcodes%2Etsv%2Egz"
cd ../../
Rscript "/data_d/WSJ/SpatialMetsDB/RawData/Seurat10Xtoh5.R" --input_dir GSM6248656/filtered_feature_bc_matrix --output_file GSM6248656/filtered_feature_bc_matrix.h5
mkdir -p GSM6248656/spatial
cd GSM6248656/spatial
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6248656&format=file&file=GSM6248656%5FNPC%5FST18%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip scalefactors_json.json.gz
gunzip tissue_hires_image.png.gz
gunzip tissue_lowres_image.png.gz
gunzip tissue_positions_list.csv.gz
cd ../../

mkdir -p GSM8324909
cd GSM8324909
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324909&format=file&file=GSM8324909%5FHLS%5F81%5FS4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8324909&format=file&file=GSM8324909%5FHLS%5F81%5FS4%5Fspaceranger%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then d=$(find spatial -mindepth 1 -maxdepth 1 -type d); mv "$d"/* spatial/ && rmdir "$d"; fi
rm -f spatial.tar.gz
cd ..



mkdir -p GSM5945509
cd GSM5945509
mkdir -p spatial
wget -O filtered_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE198353&format=file&file=GSE198353%5Fmmtv%5Fpymt%5FGEX%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O spatial.tar.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE198353&format=file&file=GSE198353%5Fmmtv%5Fpymt%5Fspatial%2Etar%2Egz"
tar -xzf spatial.tar.gz -C spatial
if [ "$(find spatial -mindepth 1 -maxdepth 1 -type d | wc -l)" -eq 1 ]; then
    d=$(find spatial -mindepth 1 -maxdepth 1 -type d)
    mv "$d"/* spatial/
    rmdir "$d"
fi
rm -f spatial.tar.gz
cd ..

# GSM6543813
mkdir -p GSM6543813
cd GSM6543813
wget -O sample.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6543813&format=file&file=GSM6543813%5FS%2D1%2Etar%2Egz"
tar -xzf sample.tar.gz
rm -f sample.tar.gz
d=$(find . -mindepth 1 -maxdepth 1 -type d | head -1)
mkdir -p spatial
find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" filtered_feature_bc_matrix.h5
find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" spatial/tissue_positions.csv
find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" spatial/
find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" spatial/tissue_hires_image.png
if [ ! -f spatial/tissue_lowres_image.png ]; then cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png; fi
rm -rf "$d"
cd ..

# GSM6543814
mkdir -p GSM6543814
cd GSM6543814
wget -O sample.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6543814&format=file&file=GSM6543814%5FS%2D2%2D2%2Etar%2Egz"
tar -xzf sample.tar.gz
rm -f sample.tar.gz
d=$(find . -mindepth 1 -maxdepth 1 -type d | head -1)
mkdir -p spatial
find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" filtered_feature_bc_matrix.h5
find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" spatial/tissue_positions.csv
find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" spatial/
find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" spatial/tissue_hires_image.png
if [ ! -f spatial/tissue_lowres_image.png ]; then cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png; fi
rm -rf "$d"
cd ..


# GSM6543815
mkdir -p GSM6543815
cd GSM6543815
wget -O sample.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6543815&format=file&file=GSM6543815%5FV19N18%2D027%2DA1%2Etar%2Egz"
tar -xzf sample.tar.gz
rm -f sample.tar.gz
d=$(find . -mindepth 1 -maxdepth 1 -type d | head -1)
mkdir -p spatial
find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" filtered_feature_bc_matrix.h5
find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" spatial/tissue_positions.csv
find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" spatial/
find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" spatial/tissue_hires_image.png
if [ ! -f spatial/tissue_lowres_image.png ]; then cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png; fi
rm -rf "$d"
cd ..

# GSM6543816
mkdir -p GSM6543816
cd GSM6543816
wget -O sample.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6543816&format=file&file=GSM6543816%5FL1101%2Etar%2Egz"
tar -xzf sample.tar.gz
rm -f sample.tar.gz
d=$(find . -mindepth 1 -maxdepth 1 -type d | head -1)
mkdir -p spatial
find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" filtered_feature_bc_matrix.h5
find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" spatial/tissue_positions.csv
find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" spatial/
find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" spatial/tissue_hires_image.png
if [ ! -f spatial/tissue_lowres_image.png ]; then cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png; fi
rm -rf "$d"
cd ..

# GSM6543817
mkdir -p GSM6543817
cd GSM6543817
wget -O sample.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6543817&format=file&file=GSM6543817%5FL1104%2Etar%2Egz"
tar -xzf sample.tar.gz
rm -f sample.tar.gz
d=$(find . -mindepth 1 -maxdepth 1 -type d | head -1)
mkdir -p spatial
find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" filtered_feature_bc_matrix.h5
find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" .
find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" spatial/tissue_positions.csv
find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" spatial/
find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" spatial/tissue_hires_image.png
if [ ! -f spatial/tissue_lowres_image.png ]; then cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png; fi
rm -rf "$d"
cd ..

# GSM5732357
mkdir -p GSM5732357
cd GSM5732357
wget -O raw_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732357&format=file&file=GSM5732357%5FA%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732357&format=file&file=GSM5732357%5FA%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732357&format=file&file=GSM5732357%5FA%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732357&format=file&file=GSM5732357%5FA%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip -f *.gz
mkdir -p spatial
mv raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
mv tissue_positions_list.csv spatial/tissue_positions.csv
mv scalefactors_json.json spatial/
mv tissue_hires_image.png spatial/
cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' spatial/scalefactors_json.json > tmp.json && mv tmp.json spatial/scalefactors_json.json
cd ..

# GSM5732358
mkdir -p GSM5732358
cd GSM5732358
wget -O raw_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732358&format=file&file=GSM5732358%5FB%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732358&format=file&file=GSM5732358%5FB%5Fscalefactors%5Fjson%2Ejson%2Egz"
wget -O tissue_hires_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732358&format=file&file=GSM5732358%5FB%5Ftissue%5Fhires%5Fimage%2Epng%2Egz"
wget -O tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5732358&format=file&file=GSM5732358%5FB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip -f *.gz
mkdir -p spatial
mv raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
mv tissue_positions_list.csv spatial/tissue_positions.csv
mv scalefactors_json.json spatial/
mv tissue_hires_image.png spatial/
cp spatial/tissue_hires_image.png spatial/tissue_lowres_image.png
jq '.tissue_lowres_scalef = .tissue_hires_scalef' spatial/scalefactors_json.json > tmp.json && mv tmp.json spatial/scalefactors_json.json
cd ..

#GSE195661
base=".RawData/GSE195661/Users/rf/Desktop/0222_clean/EPN_visium_012122/update/neuro_onc_spatial_files"
out="../RawData"
declare -A map
map[459]=GSM5844711
map[723]=GSM5844712
map[727]=GSM5844713
map[812]=GSM5844714
map[459_2]=GSM5844715
map[723_2]=GSM5844716
map[821]=GSM5844717
map[848]=GSM5844718
map[1101]=GSM5844719
map[1269]=GSM5844720
map[928]=GSM5844721
map[1239]=GSM5844722
map[1513]=GSM5844723
map[928_2]=GSM5844724
for k in "${!map[@]}"; do
    cd "$base"
    unzip -o "${k}.zip" -d "${k}"
    gsm=${map[$k]}
    mkdir -p "$out/$gsm"
    d=$(find "$k" -mindepth 1 -maxdepth 1 -type d | head -1)
    find "$d" -name "*feature_bc_matrix.h5" | head -1 | xargs -I{} mv "{}" "$out/$gsm/filtered_feature_bc_matrix.h5"
    find "$d" -name "barcodes.tsv.gz" | head -1 | xargs -I{} mv "{}" "$out/$gsm/"
    find "$d" -name "features.tsv.gz" | head -1 | xargs -I{} mv "{}" "$out/$gsm/"
    find "$d" -name "matrix.mtx.gz" | head -1 | xargs -I{} mv "{}" "$out/$gsm/"
    mkdir -p "$out/$gsm/spatial"
    find "$d" -name "tissue_positions*.csv" | head -1 | xargs -I{} mv "{}" "$out/$gsm/spatial/"
    find "$d" -name "scalefactors_json.json" | head -1 | xargs -I{} mv "{}" "$out/$gsm/spatial/"
    find "$d" -name "*hires*.png" | head -1 | xargs -I{} mv "{}" "$out/$gsm/spatial/tissue_hires_image.png"
    if [ ! -f "$out/$gsm/spatial/tissue_lowres_image.png" ]; then
        cp "$out/$gsm/spatial/tissue_hires_image.png" \
           "$out/$gsm/spatial/tissue_lowres_image.png"
    fi
    rm -rf "$k"
done

#GSM8207499
mkdir -p  GSM8207499
mkdir -p  GSM8207499/spatial
wget -O filtered_feature_bc_matrix.h5 
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207499&format=file&file=GSM8207499%5FCPS%5FOV1RtOV3%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O tissue_positions_list.csv.gz 
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207499&format=file&file=GSM8207499%5FCPS%5FOV1RtOV3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip tissue_positions_list.csv.gz
cd ..


#GSM8361988 GSM8361990 GSM8361991
mkdir "../RawData/GSM8361988/"
wget -O "../RawData/GSM8361988/GSM8361988_ST1.tar.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8361988&format=file&file=GSM8361988%5FST1%2Etar%2Egz"
tar -zxvf "../RawData/GSM8361988/GSM8361988_ST1.tar.gz" -C ../RawData/GSM8361988/
sed -i 's/"//g' "../RawData/GSM8361988/"*
mv ../RawData/GSM8361988/ST1_data.csv ../RawData/GSM8361988/expression.csv
echo -e ",X,Y" > ../RawData/GSM8361988/metadata.csv
tail -n +2 ../RawData/GSM8361988/ST1_loc.csv | cut -d"," -f1,3,4 >> ../RawData/GSM8361988/metadata.csv
mkdir "../RawData/GSM8361990/"
wget -O "../RawData/GSM8361990/GSM8361990_ST2.tar.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8361990&format=file&file=GSM8361990%5FST2%2Etar%2Egz"
tar -zxvf "../RawData/GSM8361990/GSM8361990_ST2.tar.gz" -C ../RawData/GSM8361990/
sed -i 's/"//g' "../RawData/GSM8361990/"*
mv ../RawData/GSM8361990/ST2_data.csv ../RawData/GSM8361990/expression.csv
echo -e ",X,Y" > ../RawData/GSM8361990/metadata.csv
tail -n +2 ../RawData/GSM8361990/ST2_loc.csv | cut -d"," -f1,3,4 >> ../RawData/GSM8361990/metadata.csv
mkdir "../RawData/GSM8361991/"
wget -O "../RawData/GSM8361991/GSM8361991_ST3.tar.gz" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8361991&format=file&file=GSM8361991%5FST3%2Etar%2Egz"
tar -zxvf "../RawData/GSM8361991/GSM8361991_ST3.tar.gz" -C ../RawData/GSM8361991/
sed -i 's/"//g' "../RawData/GSM8361991/"*
mv ../RawData/GSM8361991/ST3_data.csv ../RawData/GSM8361991/expression.csv
echo -e ",X,Y" > ../RawData/GSM8361991/metadata.csv
tail -n +2 ../RawData/GSM8361991/ST3_loc.csv | cut -d"," -f1,3,4 >> ../RawData/GSM8361991/metadata.csv

#GSM4838131 GSM4838132 GSM4838133 
mkdir -p GSM4838131 GSM4838132 GSM4838133 GSM8207500 
wget -O GSM4838131/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838131&format=file&file=GSM4838131%5FVisium%5FSample%5FA%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM4838131/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838131&format=file&file=GSM4838131%5FVisium%5FSample%5FA%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM4838132/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838132&format=file&file=GSM4838132%5FVisium%5FSample%5FB%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM4838132/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838132&format=file&file=GSM4838132%5FVisium%5FSample%5FB%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM4838133/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838133&format=file&file=GSM4838133%5FVisium%5FSample%5FC%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM4838133/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4838133&format=file&file=GSM4838133%5FVisium%5FSample%5FC%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
for gsm in GSM4838131 GSM4838132 GSM4838133 ; do gunzip -f ${gsm}/tissue_positions_list.csv.gz; done

#GSM8207500 GSM8207501 GSM8207502 GSM8207503 GSM8207504 GSM8207505 GSM8207506
wget -O GSM8207500/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207500&format=file&file=GSM8207500%5FCPS%5FOV34RtOV1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207500/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207500&format=file&file=GSM8207500%5FCPS%5FOV34RtOV1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207501/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207501&format=file&file=GSM8207501%5FCPS%5FOV71%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207501/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207501&format=file&file=GSM8207501%5FCPS%5FOV71%5F1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207502/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207502&format=file&file=GSM8207502%5FCPS%5FOV5LtOV4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207502/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207502&format=file&file=GSM8207502%5FCPS%5FOV5LtOV4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207503/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207503&format=file&file=GSM8207503%5FCPS%5FOV10RTOV4%2D1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207503/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207503&format=file&file=GSM8207503%5FCPS%5FOV10RTOV4%2D1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207504/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207504&format=file&file=GSM8207504%5FCPS%5FOV19%5FLtOV1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207504/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207504&format=file&file=GSM8207504%5FCPS%5FOV19%5FLtOV1%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207505/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207505&format=file&file=GSM8207505%5FCPS%5FOV20RtOV4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207505/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207505&format=file&file=GSM8207505%5FCPS%5FOV20RtOV4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
wget -O GSM8207506/filtered_feature_bc_matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207506&format=file&file=GSM8207506%5FCPS%5FOV24RTOV4%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM8207506/tissue_positions_list.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8207506&format=file&file=GSM8207506%5FCPS%5FOV24RTOV4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
for gsm in GSM8207500 GSM8207501 GSM8207502 GSM8207503 GSM8207504 GSM8207505 GSM8207506 ; do gunzip -f ${gsm}/tissue_positions_list.csv.gz; done


#GSM6319696
mkdir -p GSM6319696/spatial
wget -O GSM6319696/raw_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319696&format=file&file=GSM6319696%5FVisium%5F2%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM6319696/spatial/tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319696&format=file&file=GSM6319696%5FVisium%5F2%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip GSM6319696/spatial/tissue_positions_list.csv.gz

#GSM6319697
mkdir -p GSM6319697/spatial
wget -O GSM6319697/raw_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319697&format=file&file=GSM6319697%5FVisium%5F3%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM6319697/spatial/tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319697&format=file&file=GSM6319697%5FVisium%5F3%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip GSM6319697/spatial/tissue_positions_list.csv.gz

#GSM6319698
mkdir -p GSM6319698/spatial
wget -O GSM6319698/raw_feature_bc_matrix.h5 \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319698&format=file&file=GSM6319698%5FVisium%5F4%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget -O GSM6319698/spatial/tissue_positions_list.csv.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6319698&format=file&file=GSM6319698%5FVisium%5F4%5Ftissue%5Fpositions%5Flist%2Ecsv%2Egz"
gunzip GSM6319698/spatial/tissue_positions_list.csv.gz

#SCP2046
mkdir -p sham1/filtered_feature_bc_matrix sham1/spatial sham2/filtered_feature_bc_matrix sham2/spatial cancer1/filtered_feature_bc_matrix cancer1/spatial cancer2/filtered_feature_bc_matrix cancer2/spatial
wget -O sham1/filtered_feature_bc_matrix/barcodes.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham1_barcodes.tsv.gz"
wget -O sham1/filtered_feature_bc_matrix/features.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham1_features.tsv.gz"
wget -O sham1/filtered_feature_bc_matrix/matrix.mtx.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham1_matrix.mtx.gz"
wget -O sham2/filtered_feature_bc_matrix/barcodes.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham2_barcodes.tsv.gz"
wget -O sham2/filtered_feature_bc_matrix/features.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham2_features.tsv.gz"
wget -O sham2/filtered_feature_bc_matrix/matrix.mtx.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=sham2_matrix.mtx.gz"
wget -O cancer1/filtered_feature_bc_matrix/barcodes.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer1_barcodes.tsv.gz"
wget -O cancer1/filtered_feature_bc_matrix/features.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer1_features.tsv.gz"
wget -O cancer1/filtered_feature_bc_matrix/matrix.mtx.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer1_matrix.mtx.gz"
wget -O cancer2/filtered_feature_bc_matrix/barcodes.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer2_barcodes.tsv.gz"
wget -O cancer2/filtered_feature_bc_matrix/features.tsv.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer2_features.tsv.gz"
wget -O cancer2/filtered_feature_bc_matrix/matrix.mtx.gz "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=cancer2_matrix.mtx.gz"
wget -O sham1/spatial/tissue_positions_list.csv "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=spatial_s1.csv"
wget -O sham2/spatial/tissue_positions_list.csv "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=spatial_s2.csv"
wget -O cancer1/spatial/tissue_positions_list.csv "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=spatial_c1.csv"
wget -O cancer2/spatial/tissue_positions_list.csv "https://singlecell.broadinstitute.org/single_cell/data/public/SCP2046/liver-visium-spatial-transcriptomics?filename=spatial_c2.csv"


#SCP2046 Convert 10X directory to h5
cd ../RawData
for sample in sham1 sham2 cancer1 cancer2
do
Rscript Seurat10Xtoh5.R --input_dir ${sample}/filtered_feature_bc_matrix --output_file ${sample}/filtered_feature_bc_matrix.h5
done

#SCP2046 Mouse to Human
for id in sham1 sham2 cancer1 cancer2
do
echo "Processing $id ..."
python convert_mouse_to_human_10x.py --h5 "$id/filtered_feature_bc_matrix.h5" --biomart biomart_humna_mouse.txt --outdir "$id/filtered_feature_bc_matrix/"
gzip -f $id/filtered_feature_bc_matrix/*
mv $id/filtered_feature_bc_matrix.h5 $id/filtered_feature_bc_matrix_original.h5
Rscript Seurat10Xtoh5.R --input_dir "$id/filtered_feature_bc_matrix/" --output_file "$id/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5"
done

#GSM5538852 GSM5538853 GSM5538854 GSM5538855 GSM6612124 GSM6612125 GSM6612126 GSM6612127 Fix features.tsv.gz
samples=(GSM5538852 GSM5538853 GSM5538854 GSM5538855 GSM6612124 GSM6612125 GSM6612126 GSM6612127)
for gsm in "${samples[@]}"
do
cd ../RawData/${gsm}/filtered_feature_bc_matrix
cp features.tsv.gz features.tsv.gz.bak
zcat features.tsv.gz | awk 'BEGIN{FS=OFS="\t"}{gsub(/^GRCh38_+/,"",$2);print}' | gzip > features.tsv.gz.tmp
mv features.tsv.gz.tmp features.tsv.gz
echo "${gsm} fixed"
done

#GSM6177624 GSM6177625 Backup h5
for gsm in GSM6177624 GSM6177625
do
cp ../RawData/${gsm}/filtered_feature_bc_matrix.h5 ../RawData/${gsm}/filtered_feature_bc_matrix.h5.bak
done

#GSM6177624 GSM6177625 Fix gene names in h5
python fix_GSM6177624_6177625.py

#GSM5538852 GSM5538853 GSM5538854 GSM5538855 GSM6612124 GSM6612125 GSM6612126 GSM6612127 Rebuild h5
for gsm in GSM5538852 GSM5538853 GSM5538854 GSM5538855 GSM6612124 GSM6612125 GSM6612126 GSM6612127
do
Rscript Seurat10Xtoh5.R --input_dir ../RawData/${gsm}/filtered_feature_bc_matrix --output_file ../RawData/${gsm}/filtered_feature_bc_matrix.h5
done

#GSE195575 Download
mkdir -p GSE195575
cd GSE195575
wget -O GSE195575_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE195575&format=file"
tar -xvf GSE195575_RAW.tar

#GSE195575 Organize files
for gsm in GSM7016921 GSM7016922 GSM7016923 GSM7016924 GSM7016925 GSM7016926
do
mkdir -p ${gsm}/filtered_feature_bc_matrix
mkdir -p ${gsm}/spatial
mv ${gsm}*_barcodes.tsv.gz ${gsm}/filtered_feature_bc_matrix/barcodes.tsv.gz
mv ${gsm}*_features.tsv.gz ${gsm}/filtered_feature_bc_matrix/features.tsv.gz
mv ${gsm}*_matrix.mtx.gz ${gsm}/filtered_feature_bc_matrix/matrix.mtx.gz
gunzip -c ${gsm}*_tissue_positions_list.csv.gz > ${gsm}/spatial/tissue_positions_list.csv
done

rm -f GSE195575_RAW.tar

#GSE195575 Build h5
cd ../RawData/GSE195575
for gsm in GSM7016921 GSM7016922 GSM7016923 GSM7016924 GSM7016925 GSM7016926
do
Rscript ../Seurat10Xtoh5.R --input_dir ${gsm}/filtered_feature_bc_matrix --output_file ${gsm}/filtered_feature_bc_matrix.h5
done

#GSE195575 Generate expression.csv and metadata.csv
python GSE195575_expression_metadata.py

#GSE171351 Download
mkdir -p GSE171351
cd GSE171351
wget -O GSE171351_combined_visium.h5ad.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE171351&format=file&file=GSE171351%5Fcombined%5Fvisium%2Eh5ad%2Egz"
gunzip -f GSE171351_combined_visium.h5ad.gz

#GSE171351 Split samples and generate expression.csv metadata.csv
python GSE171351_expression_metadata.py

#Global tissue_positions fix
BASE="../RawData"
find "$BASE" -type d -path "*/spatial" | while read spatial_dir
do
file_wrong="$spatial_dir/tissue_positions.csv"
file_list="$spatial_dir/tissue_positions_list.csv"
if [ -f "$file_wrong" ]; then
mv "$file_wrong" "$file_list"
fi
if [ -f "$file_list" ] && head -n 1 "$file_list" | grep -qi "barcode"; then
tail -n +2 "$file_list" > "$file_list.tmp"
mv "$file_list.tmp" "$file_list"
fi
done

#convert_mouse_to_human_10x.py
samples=(
GSM5945509
GSM6543813
GSM6543814
GSM6543815
GSM6543816
GSM6543817
GSM7688224
GSM7688225
GSM8428408
GSM8428410
GSM6210024
GSM6505142
GSM6505143
GSM6505144
GSM6505145
GSM6505146
GSM6661400
GSM6661401
GSM6859066
GSM6859067
GSM6859068
GSM6859069
GSM6859070
GSM6859071
GSM6859072
GSM6859073
GSM6210023
GSM7422460
GSM7422461
GSM5808054
GSM5808055
GSM5808056
GSM5808057
GSM6505136
GSM6505137
GSM6505138
GSM6505139
GSM6505140
GSM6505141
GSM7818789
GSM7818790
GSM7818791
GSM7818792
GSM7461403
GSM7461404
GSM5643203
GSM6604700
GSM6604701
)
for sample in "${samples[@]}"
do
    echo "=============================="
    echo "Processing ${sample}"
    echo "=============================="

    python convert_mouse_to_human_10x.py \
        --h5 "/data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix.h5" \
        --biomart biomart_humna_mouse.txt \
        --outdir "/data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix/"
    gzip -f /data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix/*
    mv \
        /data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix.h5 \
        /data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix_original.h5
    Rscript Seurat10Xtoh5.R \
        --input_dir "/data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix/" \
        --output_file "/data_d/WSJ/SpatialMetsDB/RawData/${sample}/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5"
done

#Processing
echo -e "GSM6319696\nGSM6319697\tGSM6319698" | while read rawsample
do
  sed -i 's/"//g' ../RawData/$rawsample/expression.csv
  head -n 1 ../RawData/$rawsample/expression.csv > ../RawData/$rawsample/expression.1.csv
  awk 'BEGIN{FS=OFS=","}\
     ARGIND==1{split($0,a,"\t");if((a[1]!="") && (a[2]!="")){Gene[a[1]]=a[2]}}\
     ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
     biomart_humna_mouse.txt\
     ../RawData/$rawsample/expression.csv >> ../RawData/$rawsample/expression.1.csv
  mv ../RawData/$rawsample/expression.1.csv ../RawData/$rawsample/expression.csv
done

grep GSE195575 Retest.txt | tail -n +5 | cut -f8 | while read rawsample
do
  head -n 1 ../RawData/$rawsample/expression.csv > ../RawData/$rawsample/expression.1.csv
  awk 'BEGIN{FS=OFS=","}\
     ARGIND==1{split($0,a,"\t");if((a[1]!="") && (a[2]!="")){Gene[a[1]]=a[2]}}\
     ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
     biomart_humna_mouse.txt\
     ../RawData/$rawsample/expression.csv >> ../RawData/$rawsample/expression.1.csv
  mv ../RawData/$rawsample/expression.1.csv ../RawData/$rawsample/expression.csv
  
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
done


grep GSE263920 Retest.txt | cut -f8 | while read rawsample
do
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
done

grep GSE207592 Retest.txt | cut -f8 | while read rawsample
do
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
done

ls ../RawData/GSE260797/*.digital_expression.txt | while read file
do
  rawsample=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]}' `
  mkdir ../RawData/$rawsample
  grep -v "#" $file > ../RawData/$rawsample/expression.csv
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

rawsample="GSM6319698"
sed -i 's/"//g' ../RawData/$rawsample/expression.csv
head -n 1 ../RawData/$rawsample/expression.csv > ../RawData/$rawsample/expression.1.csv
awk 'BEGIN{FS=OFS=","}\
   ARGIND==1{split($0,a,"\t");if((a[1]!="") && (a[2]!="")){Gene[a[1]]=a[2]}}\
   ARGIND==2 && (FNR>1){if($1 in Gene){$1=Gene[$1]};print}'\
   biomart_humna_mouse.txt\
   ../RawData/$rawsample/expression.csv >> ../RawData/$rawsample/expression.1.csv
mv ../RawData/$rawsample/expression.1.csv ../RawData/$rawsample/expression.csv
GSM6319696

