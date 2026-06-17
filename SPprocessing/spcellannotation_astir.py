#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cell type annotation using Astir
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import argparse
# Astir
from astir.astir import Astir


def spcellannotation(
    csvfile,
    AstirMarkerFile,
    OutputFold,
    annotaionFile,
    filter_celltypes=None
):
    """
    Perform cell type annotation using Astir

    Parameters
    ----------
    csvfile : str
        Path to input protein expression CSV file.
        Rows = cells, Columns = proteins.

    AstirMarkerFile : str
        JSON marker file for Astir.

    OutputFold : str
        Output directory.

    annotaionFile : str
        Output annotation CSV filename.

    filter_celltypes : list
        Cell types to exclude from Astir annotation.

    Returns
    -------
    annotation : pd.DataFrame
        DataFrame containing predicted cell types
        and confidence scores.
    """

    if filter_celltypes is None:
        filter_celltypes = []

    # ==================================================
    # 1. Check input file
    # ==================================================
    if not os.path.exists(csvfile):
        print(f"ERROR: Input file does not exist:\n{csvfile}")
        sys.exit(1)

    if not os.path.exists(AstirMarkerFile):
        print(f"ERROR: Marker file does not exist:\n{AstirMarkerFile}")
        sys.exit(1)

    # ==================================================
    # 2. Create output folder
    # ==================================================
    os.makedirs(OutputFold, exist_ok=True)

    # ==================================================
    # 3. Read protein expression matrix
    # ==================================================
    try:
        protein_avg = pd.read_csv(csvfile, index_col=0)
    except Exception as e:
        print(f"Failed to read CSV file:\n{e}")
        sys.exit(1)

    # ==================================================
    # 4. Remove all-zero columns
    # ==================================================
    protein_avg = protein_avg.loc[:, (protein_avg != 0).any(axis=0)]
    protein_use = protein_avg.copy()

    # ==================================================
    # 5. Load Astir marker file
    # ==================================================
    with open(AstirMarkerFile, "r") as f:
        astir_marker_dict = json.load(f)

    # ==================================================
    # 6. Filter markers that exist in dataset
    # ==================================================
    protein_set = set(protein_use.columns)

    new_astir_marker_dict = {"cell_type": {}}

    for cell_type, markers in astir_marker_dict["cell_type"].items():

        # skip filtered cell types
        if cell_type in filter_celltypes:
            continue

        valid_markers = [
            marker for marker in markers
            if marker in protein_set
        ]

        if len(valid_markers) > 0:
            new_astir_marker_dict["cell_type"][cell_type] = valid_markers

    print("\nAstir cell types used:")
    print(list(new_astir_marker_dict["cell_type"].keys()))

    # ==================================================
    # 7. Check marker availability
    # ==================================================
    if len(new_astir_marker_dict["cell_type"]) == 0:
        raise ValueError("No valid Astir markers found.")

    # ==================================================
    # 8. Run Astir
    # ==================================================
    print("\nRunning Astir annotation...")

    Ast = Astir(
        protein_use,
        new_astir_marker_dict
    )

    Ast.fit_type()

    # ==================================================
    # 9. Get probabilities
    # ==================================================
    assignments = Ast.get_celltype_probabilities()

    print("\nAssignment matrix shape:")
    print(assignments.shape)

    # ==================================================
    # 10. Remove 'Other' category
    # ==================================================
    if 'Other' in assignments.columns:

        assignments_no_other = assignments.drop(
            columns=['Other']
        )

        print("\nRemoved 'Other' column.")
        print("Remaining cell types:")
        print(list(assignments_no_other.columns))

    else:
        assignments_no_other = assignments.copy()

    # ==================================================
    # 11. Predict labels
    # ==================================================
    pred_astir = assignments_no_other.idxmax(axis=1)

    astir_conf = assignments_no_other.max(axis=1)

    # ==================================================
    # 12. Save annotation
    # ==================================================
    annotation = pd.DataFrame({
        'celltype': pred_astir,
        'astir_conf': astir_conf
    })

    # keep original cell IDs
    annotation.index = protein_use.index

    # ==================================================
    # 13. Save outputs
    # ==================================================
    prob_path = os.path.join(
        OutputFold,
        "celltype_annotation_astir.txt"
    )

    annotation.to_csv(prob_path)
    
    # ==================================================
    # 14. Comparison
    # ==================================================
    if annotaionFile is not None and os.path.exists(annotaionFile):
        labels = pd.read_csv(annotaionFile, index_col=0)
        common_idx = labels.index.intersection(annotation.index)
        comparison_df = pd.DataFrame({'Ground_Truth': labels.loc[common_idx,'celltype'],'Predicted': annotation.loc[common_idx,'celltype']})
        ct = pd.crosstab(comparison_df['Ground_Truth'],  comparison_df['Predicted'], margins=True)
        print(ct.to_string())

# ======================================================
# Main
# ======================================================
def main():
    parser = argparse.ArgumentParser(
        description="Spatial Proteomics Cell Annotation Pipeline"
    )

    parser.add_argument("--input", required=True, help="Protein expression CSV")
    parser.add_argument("--astir_marker", required=True, help="Astir marker JSON")
    parser.add_argument("--output", required=True, help="Output folder")
    
    parser.add_argument("--anno_file", default=None, help="Optional: Ground truth annotation file for validation")
    parser.add_argument("--filter_celltypes", nargs="+", default=[], help="Cell types to exclude")


    args = parser.parse_args()
    try:
        print("\n[Step 1] Running cell type annotation...")

        annotation = spcellannotation(
            csvfile=args.input,
            AstirMarkerFile=args.astir_marker,
            OutputFold=args.output,
            annotaionFile=args.anno_file,
            filter_celltypes=args.filter_celltypes
        )

        print("\nPipeline finished successfully.")

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()