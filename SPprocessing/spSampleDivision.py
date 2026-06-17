#!/usr/bin/env python3

"""
Split multi-sample dataset into per-sample subsets.

For each sample_id:
- Extract corresponding cells
- Subset all input tables
- Save results into:

Example:
Input:
    --output_dir ../Output/

Output:
    ../Output/TSOPxxx/
    ../Output/TSOPxxx/
"""

import pandas as pd
import os
import argparse


def split_by_sample(base_dir: str, output_dir: str):
    """
    Split dataset by sample_id and save into separate directories.

    Parameters
    ----------
    base_dir : str
        Directory containing input CSV files:
        - samples.csv
        - protein_log.csv
        - celltype_assignments.csv
        - celltype_annotation.csv

    output_dir : str
        Output directory prefix (e.g. ../Output/)
    """

    # -----------------------
    # Input file paths
    # -----------------------
    samples_file = os.path.join(base_dir, "samples.TSO.csv")
    protein_file = os.path.join(base_dir, "protein_avg.csv")
    annot_file   = os.path.join(base_dir, "celltype_annotation_final.txt")
    spatial_file   = os.path.join(base_dir, "spatial.csv")

    # -----------------------
    # Load data
    # -----------------------
    samples = pd.read_csv(samples_file, index_col=0)
    protein = pd.read_csv(protein_file, index_col=0)
    annot   = pd.read_csv(annot_file, index_col=0)
    spatial   = pd.read_csv(spatial_file, index_col=0)

    # -----------------------
    # Align indices across all tables
    # (critical to avoid mismatch)
    # -----------------------
    common_idx = samples.index.intersection(protein.index)\
                              .intersection(annot.index)\
                              .intersection(spatial.index)\

    samples = samples.loc[common_idx]
    protein = protein.loc[common_idx]
    annot   = annot.loc[common_idx]
    spatial   = spatial.loc[common_idx]

    print(f"Total cells after alignment: {len(common_idx)}")

    # -----------------------
    # Split by sample_id
    # -----------------------
    for sample_id, idx in samples.groupby("sample_id").groups.items():

        cells = idx

        # Construct output directory:
        out_dir = os.path.join(output_dir, sample_id)
        os.makedirs(out_dir, exist_ok=True)

        # Subset data
        samples_sub = samples.loc[cells]
        protein_sub = protein.loc[cells]
        annot_sub   = annot.loc[cells]
        spatial_sub   = spatial.loc[cells]

        # -----------------------
        # Save outputs
        # -----------------------
        samples_sub.to_csv(os.path.join(out_dir, "samples.csv"))
        protein_sub.to_csv(os.path.join(out_dir, "protein_avg.csv"))
        annot_sub.to_csv(os.path.join(out_dir, "celltype_annotation_final.txt"))
        spatial_sub.to_csv(os.path.join(out_dir, "spatial.csv"))

        print(f"Saved: {out_dir} ({len(cells)} cells)")


# =======================
# CLI
# =======================
def main():
    parser = argparse.ArgumentParser(
        description="Split dataset by sample_id"
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input CSV files"
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output prefix (e.g. ../Output/)"
    )

    args = parser.parse_args()

    split_by_sample(
        base_dir=args.input_dir,
        output_dir=args.output_dir
    )


if __name__ == "__main__":
    main()