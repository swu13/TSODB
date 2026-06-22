### TSODB Data Processing and Annotation Pipeline



This repository contains all scripts used for data downloading, preprocessing, annotation, and downstream analyses in the **Tumor Spatial Omics Database (TSODB)**.



TSODB integrates multiple spatial omics technologies, including:



* Spatial Transcriptomics (ST)

* 10x Genomics Visium

* 10x Genomics Visium HD

* Slide-seq / Slide-seqV2

* GeoMx Digital Spatial Profiling (DSP)

* Spatial proteomics datasets (CODEX, IMC, MIBI, mIF, etc.)



##### Repository Structure



```text

TSODB/

├── VisiumProcessing/

├── SlideseqProcessing/

├── GeoMxProcessing/

├── SPprocessing/

├── STVisiumAnnotation/

├── SlideseqAnnotation/

├── GeoMxAnnotation/

└── SPannotation/

```



---



##### Part I. Data Download and Preprocessing



###### 1. VisiumProcessing



Data download and preprocessing pipeline for Spatial Transcriptomics (ST), 10x Visium and Visium HD datasets.


| Script | Description |
|---------|-------------|
| `sTVisium.sh` | Data download and preprocessing for ST and Visium datasets. |
| `sTVisiumHD.sh` | Data download and preprocessing for Visium HD datasets. |
| `DownloadSynapse.py` | Download datasets stored in Synapse. |
| `Seurat10Xtoh5.R` | Convert 10x barcode, features and matrix files into h5 format. |
| `convert_zenodo_to_10x.R` | Convert Zenodo datasets into standard 10x format and generate h5 files. |
| `GSE171351_expression_metadata.py` | Integrated gene expression matrix and sample metadata fromGSE171351 datasets. |
| `GSE195575_expression_metadata.py` | Integrated gene expression matrix and sample metadata from GSE195575 datasets. |
| `GSE159709_GSE207592_GSE263920_expression_metadata.py` | Integrated gene expression matrix and sample metadata from GSE159709, GSE207592, and GSE263920 datasets. |
| `convert_mouse_to_human_10x.py` | Convert mouse 10x gene expression matrices to human homolog gene symbols using BioMart mappings. |
| `GPMID39592577_Batch_ENSG2Symbol.R` | Convert Ensembl gene IDs (ENSG) to gene symbols for the PMID39592577 dataset. |
| `extract_count_matrix.R` | Extract gene expression count matrices from CytAssist FFPE Protein Expression datasets. |
---

###### 2. SlideseqProcessing
Data download and preprocessing pipeline for Slide-seq datasets.

| Script | Description |
|---------|-------------|
| `sTSlide.sh` | Data download and preprocessing for Slide-seq datasets.  |
| `biomart_humna_mouse.txt` | Human-mouse gene symbol conversion table downloaded from BioMart. |

---



###### 3. GeoMxProcessing

Data download and preprocessing pipeline for GeoMx DSP datasets.
| Script | Description |
|---------|-------------|
| `sTGeoMx.sh` | Data download and preprocessing for GeoMx DSP datasets.  |

---

###### 4. SPprocessing

Data download and preprocessing pipeline for spatial proteomics datasets.

| Script | Description |
|---------|-------------|
| `*.sh` | Download and preprocessing scripts for individual spatial proteomics platforms, including manual cell type annotation. |
| `sTVisiumHD.sh` | Data download and preprocessing for Visium HD datasets. |
| `DownloadSynapse.py` | Download datasets stored in Synapse. |
| `spcellannotation_astir.py` | Cell type annotation using ASTIR. |
| `spSampleDivision.py` | Divide spatial proteomics data into individual slides. |

---

##### Part II. Data Annotation and Downstream Analysis

###### 1. STVisiumAnnotation

Annotation pipeline for Spatial Transcriptomics (ST), 10x Visium and Visium HD datasets.

| Script | Description |
|---------|-------------|
| `STVisiumPipeline.sh` | Main analysis pipeline. |
| `SpaCETdeconvolution.R` | Cell type deconvolution using SpaCET. |
| `*revised.R` | Revised functions in SpaCET and SPATA2 to improve compatibility across datasets. |
| `SpaCETtumor.R` | Identification of tumor, interface and stromal regions, and tumor subclusters using SpaCET. |
| `SpaCCcolocalization.R` | Cell enrichment analysis using the median method and cell-cell interactions using MistyR. |
| `sTcellinteraction.py` | Ligand-receptor interaction analysis using stLearn. |
| `SpatialCCIDiff.R` | Differential ligand-receptor interactions among different regions. |
| `SpatialTrajectoryGene.R` | Differential gene expression analysis along tumor-interface-stromal directions and identification of enriched genes in each region. |
| `SpatialTrajectoryPathway.R` | Differential pathway activities along tumor-interface-stromal directions and identification of activated pathways in each region. |

---

###### 2. SlideseqAnnotation

Annotation pipeline for Slide-seq datasets.

| Script | Description |
|---------|-------------|
| `STSlideseqPipeline.sh` | Main analysis pipeline. |
| `SpaCETdeconvolution.R` | Cell type deconvolution using SpaCET. |
| `*revised.R` | Revised functions in SpaCET and SPATA2 to improve compatibility across datasets. |
| `sTSlideTumor.R` | Identification of tumor, interface and stromal regions using neighborhood analysis and tumor subclusters using SpaCET. |
| `sTSlideColocalization.py` | Cell enrichment analysis and cell-cell colocalization analysis using Squidpy. |
| `sTcellinteraction.py` | Ligand-receptor interaction analysis using stLearn. |
| `SpatialCCIDiff.R` | Differential ligand-receptor interactions among different regions. |
| `SpatialTrajectoryGene.R` | Differential gene expression analysis along tumor-interface-stromal directions and identification of enriched genes in each region. |
| `SpatialTrajectoryPathway.R` | Differential pathway activities along tumor-interface-stromal directions and identification of activated pathways in each region. |

---

###### 3. GeoMxAnnotation

Annotation pipeline for GeoMx DSP datasets.

| Script | Description |
|---------|-------------|
| `GeoMxPipeline.sh` | Main analysis pipeline. |
| `sTGeoMxdeconvolution.R` | Cell type deconvolution in each ROI using SpatialDecon and CIBERSORT and differential cell composition analysis among clinical groups. |
| `sTGeoMxGene.R` | Differential pathway activity analysis among clinical groups. |
| `sTGeoMxPathway.R` | Identification of tumor, interface and stromal regions using neighborhood analysis and tumor subclusters using SpaCET. |

---

###### 4. SPannotation

Annotation pipeline for spatial proteomics datasets.

| Script | Description |
|---------|-------------|
| `SpatialProteomicsPipeline.sh` | Main analysis pipeline. |
| `spspatialplot.py` | Visualization of cell type distributions across spatial locations. |
| `sptumor.py` | Identification of tumor core, tumor boundary, stromal boundary and stromal core regions using neighborhood analysis. |
| `spprotein.py` | Differential protein expression analysis between boundary and core regions. |
| `spcolocalization.py` | Identification of enriched cell-cell colocalizations using Squidpy. |

---

##### Citation

If you use this repository, please cite:

```text

Sijia Wu, et al. TumorSpatialOmics: A spatial omics database for cancer.Journal Name, 2026.

```
---

##### Contact

For questions and suggestions, please open an issue on GitHub or contact the corresponding authors of TSODB.

