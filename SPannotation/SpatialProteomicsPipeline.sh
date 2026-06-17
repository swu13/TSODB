#!/bin/bash

tail -n +2 SPsamples.txt | grep -v PMID38194915 | grep -v PMID41679299 | grep -v HTAN_OHSU | grep -v PMID38622132 | parallel -j 5 --colsep '\t' '
  outputpath="/data_d/WSJ/SpatialMetsDB/SPOutput/{5}"
  logFile=$outputpath"/log.txt"
  rm -f $logFile
  #1. Spatial Map
  python spspatialplot.py --output $outputpath >> $logFile 2>&1
  #2. Tumor map
  python sptumor.py --output $outputpath --tumor_cutoff 0 --radius 80 >> $logFile 2>&1
  #3. Protein map
  python spprotein.py --output $outputpath >> $logFile 2>&1
  #4. Cell enrichment analysis
  python spcolocalization.py --output $outputpath >> $logFile 2>&1
'