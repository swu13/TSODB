#!/bin/bash
set -euo pipefail
############################################
# 1.10X Visium HD
############################################

#Visium_HD_6p5mm_Human_Colon_Cancer
mkdir Visium_HD_6p5mm_Human_Colon_Cancer
wget -P Visium_HD_6p5mm_Human_Colon_Cancer/ https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Human_Colon_Cancer/Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz 
tar -xzf Visium_HD_6p5mm_Human_Colon_Cancer/Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz
rm -r Visium_HD_6p5mm_Human_Colon_Cancer/binned_outputs/square_002um
rm -r Visium_HD_6p5mm_Human_Colon_Cancer/binned_outputs/square_008um



