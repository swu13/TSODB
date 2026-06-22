#!/usr/bin/env python3
import os
import scanpy as sc
import pandas as pd
from scipy.io import mmwrite


def convert_mouse_to_human(
    h5_path,
    biomart_path,
    output_dir
):

    print("?? Reading 10X h5 file...")
    adata = sc.read_10x_h5(h5_path)

    print(f"Original genes: {adata.n_vars}")

    # =========================
    # 1. Load mouse-human mapping
    # =========================
    print("?? Loading gene mapping...")
    mapping_df = pd.read_csv(biomart_path, sep="\t")

    mouse_to_human = dict(zip(
        mapping_df["Mouse gene name"],
        mapping_df["Gene name"]
    ))

    # remove NaN
    mouse_to_human = {
        k: v for k, v in mouse_to_human.items()
        if pd.notna(k) and pd.notna(v)
    }

    # =========================
    # 2. Filter genes
    # =========================
    print("?? Filtering mouse genes with mapping...")
    mask = adata.var_names.isin(mouse_to_human.keys())
    adata = adata[:, mask].copy()

    print(f"After filtering: {adata.n_vars}")

    # =========================
    # 3. Convert gene names
    # =========================
    print("?? Converting mouse to human gene names...")

    adata.var_names = [
        mouse_to_human.get(g, g) for g in adata.var_names
    ]

    # Optional: ensure unique gene names
    adata.var_names_make_unique()

    # =========================
    # 4. Create output dir
    # =========================
    os.makedirs(output_dir, exist_ok=True)

    # =========================
    # 5. Write 10X format
    # =========================
    print("?? Writing 10X format files...")

    barcodes_path = os.path.join(output_dir, "barcodes.tsv")
    features_path = os.path.join(output_dir, "features.tsv")
    matrix_path = os.path.join(output_dir, "matrix.mtx")

    # barcodes
    pd.Series(adata.obs_names).to_csv(
        barcodes_path,
        index=False,
        header=False
    )

    # features
    pd.DataFrame({
        "gene_ids": adata.var_names,
        "gene_names": adata.var_names,
        "feature_types": ["Gene Expression"] * adata.n_vars
    }).to_csv(
        features_path,
        sep="\t",
        index=False,
        header=False
    )

    # matrix (IMPORTANT: transpose!)
    mmwrite(matrix_path, adata.X.T)

    print("? Done!")
    print(f"Output folder: {output_dir}")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Convert mouse 10X h5 to human gene 10X matrix"
    )

    parser.add_argument("--h5", required=True, help="Input 10X h5 file")
    parser.add_argument("--biomart", required=True, help="Mouse-human mapping file")
    parser.add_argument("--outdir", required=True, help="Output directory")

    args = parser.parse_args()

    convert_mouse_to_human(
        h5_path=args.h5,
        biomart_path=args.biomart,
        output_dir=args.outdir
    )
