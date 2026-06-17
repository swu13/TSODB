#!/bin/bash

# GeoMx
grep GeoMx STGeoMxsamples.txt | cut -f5,7 | sort | uniq | parallel -j 10 --colsep '\t' '
  rawpath="/data_d/WSJ/SpatialMetsDB/RawData/{2}"
  outputpath="/data_d/WSJ/SpatialMetsDB/STOutput/{1}"
  logFile=$outputpath"/log.txt"
  
  ExpressionFile=$rawpath"/Expression.txt"
  ClinicalFile=$rawpath"/clinical.txt"
  
  if [[ -f $ExpressionFile ]] && [[ -f $ClinicalFile ]]
  then
    mkdir -p $outputpath
    rm -f $logFile
    #1. Cell deconvolution
    if [[ ! -f $outputpath"/CellFraction.txt" ]] || [[ ! -f $outputpath"/CellEnriched.txt" ]]
    then
      Rscript sTGeoMxdeconvolution.R --input $ExpressionFile --clinical $ClinicalFile --outdir $outputpath >> $logFile 2>&1
    fi
    #2. Gene analysis among diverse ROIs
    if [[ ! -f $outputpath"/GeneEnriched.txt" ]]
    then
      Rscript sTGeoMxGene.R --input $outputpath"/seurat.rds" --outdir $outputpath >> $logFile 2>&1
    fi
    #3. pathways analysis among diverse ROIs
    if [[ ! -f $outputpath"/PathwayEnriched.txt" ]]
    then
      Rscript sTGeoMxPathway.R --input $outputpath"/seurat.gene.rds" --outdir $outputpath >> $logFile 2>&1
    fi
    #if [[ ! -f $outputpath"/seurat.qs2" ]]
    #then
    #4. Add the slide information and reduce the seize
    # Rscript sTGeoMxDataPre.R --outpath $outputpath >> $logFile 2>&1
    #fi
  fi
'