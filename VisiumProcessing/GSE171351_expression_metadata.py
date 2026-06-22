
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import os

base_dir = "/data_d/WSJ/SpatialMetsDB/RawData/GSE171351"

sample2gsm = {
    "A1": "GSM5224027",
    "B1": "GSM5224028",
    "C1": "GSM5224029",
    "D1": "GSM5224030"
}

adata = sc.read_h5ad(
    f"{base_dir}/GSE171351_combined_visium.h5ad"
)

for sample in adata.obs["sampleID"].unique():

    gsm = sample2gsm[sample]

    print(f"Processing {sample} -> {gsm}")

    sub = adata[adata.obs["sampleID"] == sample].copy()

    outdir = f"{base_dir}/{gsm}"
    os.makedirs(outdir, exist_ok=True)

    #
    # expression.csv
    # Gene × Spot
    #
    if sp.issparse(sub.X):
        expr = pd.DataFrame(
            sub.X.toarray(),
            index=sub.obs_names,
            columns=sub.var_names
        ).T
    else:
        expr = pd.DataFrame(
            sub.X,
            index=sub.obs_names,
            columns=sub.var_names
        ).T

    expr.to_csv(
        f"{outdir}/expression.csv"
    )
    #
    # metadata.csv
    # Spot × Coordinate
    #
    meta = pd.DataFrame(
        {
            "X": sub.obsm["spatial"][:, 0],
            "Y": sub.obsm["spatial"][:, 1]
        },
        index=sub.obs_names
    )

    meta.to_csv(
        f"{outdir}/metadata.csv"
    )
    print(
        f"{gsm}: expression={expr.shape}, metadata={meta.shape}"
    )