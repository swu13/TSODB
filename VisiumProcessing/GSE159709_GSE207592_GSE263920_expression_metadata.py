import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import os

base_dir = "../RawData"

samples = [
    "GSM6319696",
    "GSM6319697",
    "GSM6319698",
    "GSM4838131",
    "GSM4838132",
    "GSM4838133",
    "GSM8207499",
    "GSM8207500",
    "GSM8207501",
    "GSM8207502",
    "GSM8207503",
    "GSM8207504",
    "GSM8207505",
    "GSM8207506"
]

for gsm in samples:

    print(f"\nProcessing {gsm}")

    sample_dir = os.path.join(
        base_dir,
        gsm
    )

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

    # expression.csv
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

    # metadata.csv
    pos = pd.read_csv(
        pos_file,
        header=None
    )

    if pos.shape[1] == 6:

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

    elif pos.shape[1] == 3:

        pos.columns = [
            "barcode",
            "X",
            "Y"
        ]

        pos = pos.set_index("barcode")

        pos = pos.loc[
            adata.obs_names
        ]

        meta = pos[
            ["X", "Y"]
        ]

    else:

        raise ValueError(
            f"{gsm}: unexpected coordinate format ({pos.shape[1]} columns)"
        )

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