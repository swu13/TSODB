suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

##################################################
# GTF
##################################################

gtf_file <- "/data_d/WSJ/STADstudy/Reference/gencode.v32.annotation.gtf"

cat("Loading GTF...\n")

gtf <- read.delim(
    gtf_file,
    comment.char="#",
    header=FALSE,
    sep="\t",
    stringsAsFactors=FALSE,
    quote=""
)

colnames(gtf) <- c(
    "chr","source","type",
    "start","end","score",
    "strand","phase","attr"
)

gtf_gene <- gtf[gtf$type=="gene",]

gene_id <- sub(
    '.*gene_id "([^"]+)".*',
    '\\1',
    gtf_gene$attr
)

gene_name <- sub(
    '.*gene_name "([^"]+)".*',
    '\\1',
    gtf_gene$attr
)

gene_map <- data.frame(
    ENSG=gsub("\\..*$","",gene_id),
    SYMBOL=gene_name,
    stringsAsFactors=FALSE
)

gene_map <- gene_map[
    !duplicated(gene_map$ENSG),
]

ensg2symbol <- setNames(
    gene_map$SYMBOL,
    gene_map$ENSG
)

cat("Loaded", nrow(gene_map), "genes\n")

##################################################
# samples
##################################################

samples <- setdiff(
    2:96,
    c(17,18)
)

##################################################
# loop
##################################################

for(id in samples){

    cat("\n=============================\n")
    cat("Sample:", id, "\n")
    cat("=============================\n")

    sample_dir <- paste0(
        "/data_d/WSJ/SpatialMetsDB/RawData/",
        id
    )

    old_h5 <- file.path(
        sample_dir,
        "filtered_feature_bc_matrix.h5"
    )

    backup_h5 <- file.path(
        sample_dir,
        "filtered_feature_bc_matrix_original.h5"
    )

    outdir <- file.path(
        sample_dir,
        "filtered_feature_bc_matrix"
    )

    ##################################################
    # rename original h5
    ##################################################

    if(!file.exists(backup_h5)){

        file.rename(
            old_h5,
            backup_h5
        )

        cat("Renamed original h5\n")
    }

    ##################################################
    # read original h5
    ##################################################

    counts <- Read10X_h5(
        backup_h5
    )

    ##################################################
    # ENSG -> Symbol
    ##################################################

    genes <- rownames(counts)

    genes2 <- gsub(
        "\\..*$",
        "",
        genes
    )

    symbol <- ensg2symbol[
        genes2
    ]

    keep <- !is.na(symbol)

    counts <- counts[keep,]

    symbol <- symbol[keep]

    rownames(counts) <- symbol

    ##################################################
    # merge duplicated genes
    ##################################################

    counts_dense <- as.matrix(counts)

    counts_merge <- rowsum(
        counts_dense,
        group=rownames(counts_dense)
    )

    rm(counts_dense)
    gc()

    counts <- as(
        counts_merge,
        "dgCMatrix"
    )

    rm(counts_merge)
    gc()

    ##################################################
    # write 10X directory
    ##################################################

    dir.create(
        outdir,
        recursive=TRUE,
        showWarnings=FALSE
    )

    writeMM(
        counts,
        file.path(outdir,"matrix.mtx")
    )

    system(
        paste(
            "gzip -f",
            file.path(outdir,"matrix.mtx")
        )
    )

    write.table(
        colnames(counts),
        file.path(outdir,"barcodes.tsv"),
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE
    )

    system(
        paste(
            "gzip -f",
            file.path(outdir,"barcodes.tsv")
        )
    )

    features <- data.frame(
        GeneID=rownames(counts),
        GeneName=rownames(counts),
        FeatureType="Gene Expression",
        stringsAsFactors=FALSE
    )

    write.table(
        features,
        file.path(outdir,"features.tsv"),
        sep="\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE
    )

    system(
        paste(
            "gzip -f",
            file.path(outdir,"features.tsv")
        )
    )

    ##################################################
    # build new h5
    ##################################################

    cmd <- paste(
        "Rscript",
        "/data_d/WSJ/SpatialMetsDB/Script/Seurat10Xtoh5.R",
        "-i", outdir,
        "-o", old_h5,
        "-s", "GRCh38"
    )

    cat(cmd,"\n")

    system(cmd)

    cat("Finished sample",id,"\n")
}
