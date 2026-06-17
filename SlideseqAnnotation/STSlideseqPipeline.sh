#!/bin/bash

#4. Slideseq
grep "Slide-seq" STsamples.txt | parallel -j 10 --colsep '\t' '
  rawpath="/data_d/WSJ/SpatialMetsDB/RawData/{8}"
  outputpath="/data_d/WSJ/SpatialMetsDB/STOutput/{5}"
  cancertype="{6}"
  species="{4}"
  species=${species,,}  #lower
  logFile=$outputpath"/log.txt"
  
  expfile=$rawpath"/expression.csv"
  clinicalfile=$rawpath"/metadata.csv"
  
  if [[ -f $expfile ]] && [[ -f $clinicalfile ]]
  then
    mkdir -p $outputpath
    rm -f $logFile
    #1. Cell deconvolution
    if [[ ! -f $outputpath"/CellAnno.txt" ]] || [[ ! -f $outputpath"/SpaCETdecon.eachcell.svg" ]] || [[ ! -f $outputpath"/SpaCETdecon.pdf" ]]
    then
      Rscript SpaCETdeconvolution.R --outpath $outputpath --cancerType $cancertype --rawdata_dir $rawpath >> $logFile 2>&1
    fi
    #2. Tumor edges and tumor subcluster
    if [[ ! -f $outputpath"/TumorInterfere.pdf" ]] || [[ ! -f $outputpath"/TumorSubcluster.pdf" ]]
    then
      Rscript sTSlideTumor.R --outpath $outputpath >> $logFile 2>&1
    fi
    #3. Tumor edge, tumor core, and strome spot colocalization
    if [[ ! -f $outputpath"/Cellcolocalization.txt" ]]
    then
      python sTSlideColocalization.py --rawpath $rawpath --output $outputpath --region_col "Interface" >> $logFile 2>&1
    fi
    #4. Tumor edge, tumor core, and strome spot interaction
    if [[ ! -f $outputpath"/CCIEnriched.txt" ]] || [[ ! -f $outputpath"/lr_interaction_score.Tumor.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Interface.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Stroma.csv" ]]
    then
      python sTcellinteraction.py --rawpath $rawpath --region_file $outputpath"/cell_annotation_spatial.csv" --output $outputpath --species "human" >> $logFile 2>&1
      Rscript SpatialCCIDiff.R --outpath $outputpath >> $logFile 2>&1
    fi
    #5. Genes analysis along with tumor-interface-stromal trajectory and gene enriched analysis
    if [[ ! -f $outputpath"/TrajectoryGenesLinePlot.txt" ]] || [[ ! -f $outputpath"/TrajectoryGenes.txt" ]] || [[ ! -f $outputpath"/GeneEnriched.txt" ]]
    then
      Rscript SpatialTrajectoryGene.R --outpath $outputpath --rawpath $rawpath >> $logFile 2>&1
    fi
    #6. pathways analysis along with tumor-interface-stromal trajectory and pathways enriched analysis
    if [[ ! -f $outputpath"/TrajectoryGenesetsLinePlot.txt" ]] || [[ ! -f $outputpath"/TrajectoryPathways.txt" ]] || [[ ! -f $outputpath"/PathwayEnriched.txt" ]]
    then
      Rscript SpatialTrajectoryPathway.R --outpath $outputpath >> $logFile 2>&1
    fi
    #7. reduce the rds size
    #Rscript ReduceSize.R --outpath $outputpath >> $logFile 2>&1  #seurat_pathway_obj.rds
  fi
'

