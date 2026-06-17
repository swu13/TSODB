#!/bin/bash

#1. Visium/10X Visuim/
tail -n +2 STsamples.txt | grep -v Mouse | parallel -j 10 --colsep '\t' '
  rawpath="/data_d/WSJ/SpatialMetsDB/RawData/{8}"
  outputpath="/data_d/WSJ/SpatialMetsDB/STOutput/{5}"
  cancertype="{6}"
  species="{4}"
  species=${species,,}  #lower
  logFile=$outputpath"/log.txt"
  
  #h5file=$rawpath"/filtered_feature_bc_matrix.h5"
  #spatialfile1=$rawpath"/spatial/tissue_lowres_image.png"
  #spatialfile2=$rawpath"/spatial/tissue_positions_list.csv"
  #spatialfile3=$rawpath"/spatial/scalefactors_json.json"
  
  #if [[ -f $h5file ]] && [[ -f $spatialfile1 ]] && [[ -f $spatialfile2 ]] && [[ -f $spatialfile3 ]]
  #then
    mkdir -p $outputpath
    #rm -f $logFile
    #1. Cell deconvolution
    if [[ ! -f $outputpath"/CellAnno.txt" ]] || [[ ! -f $outputpath"/SpaCETdecon.eachcell.svg" ]] || [[ ! -f $outputpath"/SpaCETdecon.svg" ]]
    then
      Rscript SpaCETdeconvolution.R --outpath $outputpath --cancerType $cancertype --rawdata_dir $rawpath >> $logFile 2>&1
    fi
    #2. Tumor edges and tumor subcluster
    if [[ ! -f $outputpath"/TumorInterfere.svg" ]] || [[ ! -f $outputpath"/TumorSubcluster.svg" ]]
    then
      Rscript SpaCETtumor.R --outpath $outputpath >> $logFile 2>&1
    fi
    #3. Tumor edge, tumor core, and strome spot colocalization
    if [[ ! -f $outputpath"/median_CC.txt" ]] || [[ ! -f $outputpath"/mistyR_CC.txt" ]]
    then
      Rscript SpaCCcolocalization.R --outpath $outputpath >> $logFile 2>&1
    fi
    #4. Tumor edge, tumor core, and strome spot interaction
    if [[ ! -f $outputpath"/CCIEnriched.txt" ]] || [[ ! -f $outputpath"/lr_interaction_score.Tumor.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Interface.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Stroma.csv" ]]
    then
      python sTcellinteraction.py --rawpath $rawpath --region_file $outputpath"/RegionAnno.txt" --output $outputpath --species "human" >> $logFile 2>&1
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
  #fi
'
#2. Visium HD
grep "VisiumHD" STsamples.txt | awk -F"\t" '($8!="")' | parallel -j 1 --colsep '\t' '
  rawpath="/data_d/WSJ/SpatialMetsDB/RawData/{8}"
  outputpath="/data_d/WSJ/SpatialMetsDB/STOutput/{5}"
  cancertype="{6}"
  species="{4}"
  species=${species,,}  #lower
  logFile=$outputpath"/log.txt"
  
  h5file=$rawpath"/binned_outputs/square_016um/filtered_feature_bc_matrix.h5"
  spatialfile1=$rawpath"/binned_outputs/square_016um/spatial/tissue_lowres_image.png"
  spatialfile2=$rawpath"/binned_outputs/square_016um/spatial/tissue_positions.parquet"
  spatialfile3=$rawpath"/binned_outputs/square_016um/spatial/scalefactors_json.json"
  
  if [[ -f $h5file ]] && [[ -f $spatialfile1 ]] && [[ -f $spatialfile2 ]] && [[ -f $spatialfile3 ]]
  then
    mkdir -p $outputpath
    #rm -f $logFile
    #1. Cell deconvolution
    if [[ ! -f $outputpath"/CellAnno.txt" ]] || [[ ! -f $outputpath"/SpaCETdecon.eachcell.svg" ]] || [[ ! -f $outputpath"/SpaCETdecon.svg" ]]
    then
      Rscript SpaCETdeconvolution.R --outpath $outputpath --cancerType $cancertype --rawdata_dir $rawpath"/binned_outputs/square_016um/" >> $logFile 2>&1
    fi
    #2. Tumor edges and tumor subcluster
    if [[ ! -f $outputpath"/TumorInterfere.svg" ]] || [[ ! -f $outputpath"/TumorSubcluster.svg" ]]
    then
      Rscript SpaCETtumor.R --outpath $outputpath --cores 1 >> $logFile 2>&1
    fi
    #3. Tumor edge, tumor core, and strome spot colocalization
    if [[ ! -f $outputpath"/median_CC.txt" ]] || [[ ! -f $outputpath"/mistyR_CC.txt" ]]
    then
      Rscript SpaCCcolocalization.R --outpath $outputpath >> $logFile 2>&1
    fi
    #4. Tumor edge, tumor core, and strome spot interaction
    if [[ ! -f $outputpath"/CCIEnriched.txt" ]] || [[ ! -f $outputpath"/lr_interaction_score.Tumor.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Interface.csv" ]] || [[ ! -f $outputpath"/lr_interaction_score.Stroma.csv" ]]
    then
      originalpath=`pwd`
      cd $rawpath"/binned_outputs/square_016um/"
      Rscript -e "
        library(arrow)
        pos <- read_parquet(\"spatial/tissue_positions.parquet\")
        write.table(pos,file = \"spatial/tissue_positions_list.csv\",sep = \",\",row.names = FALSE,col.names = FALSE,quote = FALSE)
      "
      cd $originalpath
      python sTcellinteraction.py --rawpath $rawpath"/binned_outputs/square_016um/" --region_file $outputpath"/RegionAnno.txt" --output $outputpath --species $species >> $logFile 2>&1
      rm -f $rawpath"/binned_outputs/square_016um/spatial/tissue_positions_list.csv"
      Rscript SpatialCCIDiff.R --outpath $outputpath --rawpath $rawpath"/binned_outputs/square_016um/" >> $logFile 2>&1
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

#3. GeoMx

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
    if [[ ! -f $outputpath"/seurat.qs2" ]]
    then
    #4. Add the slide information and reduce the seize
      Rscript sTGeoMxDataPre.R --outpath $outputpath >> $logFile 2>&1
    fi
  fi
'

#4. Slideseq
grep "Slide-seq" STSlidesamples.txt | parallel -j 10 --colsep '\t' '
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
      Rscript SpaCETdeconvolution_new.R --outpath $outputpath --cancerType $cancertype --rawdata_dir $rawpath >> $logFile 2>&1
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
      python sTcellinteraction_new.py --rawpath $rawpath --region_file $outputpath"/cell_annotation_spatial.csv" --output $outputpath --species "human" >> $logFile 2>&1
      Rscript SpatialCCIDiff_new.R --outpath $outputpath >> $logFile 2>&1
    fi
    #5. Genes analysis along with tumor-interface-stromal trajectory and gene enriched analysis
    if [[ ! -f $outputpath"/TrajectoryGenesLinePlot.txt" ]] || [[ ! -f $outputpath"/TrajectoryGenes.txt" ]] || [[ ! -f $outputpath"/GeneEnriched.txt" ]]
    then
      Rscript SpatialTrajectoryGene_new.R --outpath $outputpath --rawpath $rawpath >> $logFile 2>&1
    fi
    #6. pathways analysis along with tumor-interface-stromal trajectory and pathways enriched analysis
    if [[ ! -f $outputpath"/TrajectoryGenesetsLinePlot.txt" ]] || [[ ! -f $outputpath"/TrajectoryPathways.txt" ]] || [[ ! -f $outputpath"/PathwayEnriched.txt" ]]
    then
      Rscript SpatialTrajectoryPathway_new.R --outpath $outputpath >> $logFile 2>&1
    fi
    #7. reduce the rds size
    #Rscript ReduceSize.R --outpath $outputpath >> $logFile 2>&1  #seurat_pathway_obj.rds
  fi
'



