import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import os

base_dir = "../RawData/GSE195575"

samples = [
    "GSM7016921",
    "GSM7016922",
    "GSM7016923",
    "GSM7016924",
    "GSM7016925",
    "GSM7016926"
]

for gsm in samples:
    print(f"\nProcessing {gsm}")

    sample_dir = os.path.join(base_dir, gsm)

    h5_file = os.path.join(
        sample_dir,
        "filtered_feature_bc_matrix.h5"
    )

    pos_file = os.path.join(
        sample_dir,
        "spatial",
        "tissue_positions_list.csv"
    )

    adata = sc.read_10x_h5(h5_file)

    if sp.issparse(adata.X):
        expr = pd.DataFrame(
            adata.X.toarray(),
            index=adata.obs_names,
            columns=adata.var_names
        ).T
    else:
        expr = pd.DataFrame(
            adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        ).T

    expr.to_csv(
        os.path.join(
            sample_dir,
            "expression.csv"
        )
    )

    pos = pd.read_csv(
        pos_file,
        header=None
    )

    pos.columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "Y",
        "X"
    ]

    pos = pos.set_index("barcode")

    pos = pos.loc[
        adata.obs_names
    ]

    meta = pos[
        ["X", "Y"]
    ]

    meta.to_csv(
        os.path.join(
            sample_dir,
            "metadata.csv"
        )
    )

    print(
        f"{gsm}: expression={expr.shape}, metadata={meta.shape}"
    )

    print(
        "Barcode match:",
        all(expr.columns == meta.index)
    )

print("\nAll samples finished.")