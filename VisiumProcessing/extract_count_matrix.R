library(Seurat)
library(Matrix)

samples <- c(
    "CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer",
    "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma"
)

baseDir <- "/data_d/WSJ/SpatialMetsDB/RawData"

for(sample in samples){

    cat("Processing:", sample, "\n")

    visiumPath <- file.path(baseDir, sample)

    st.matrix.data <- Read10X_h5(
        file.path(visiumPath, "filtered_feature_bc_matrix.h5")
    )

    print(names(st.matrix.data))

    expr <- st.matrix.data[["Gene Expression"]]

    print(class(expr))
    print(dim(expr))

    outdir <- file.path(visiumPath, "rebuild_matrix")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    writeMM(
        expr,
        file.path(outdir, "matrix.mtx")
    )

    write.table(
        colnames(expr),
        file.path(outdir, "barcodes.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    write.table(
        data.frame(
            gene_id = rownames(expr),
            gene_name = rownames(expr)
        ),
        file.path(outdir, "features.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    cat("Finished:", sample, "\n\n")
}