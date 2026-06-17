#!/usr/bin/env python3

"""
Differential protein analysis: Interface vs Stromal
"""

import pandas as pd
import numpy as np
import os
import argparse
import sys
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def sptumor(output):

    # =======================
    # Check files
    # =======================
    proteinfile   = os.path.join(output, "protein_avg.csv")
    region_file   = os.path.join(output, "region.csv")
    for f in [proteinfile, region_file]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    os.makedirs(output, exist_ok=True)

    # =======================
    # Load data
    # =======================
    protein_log = pd.read_csv(proteinfile, index_col=0)
    regions = pd.read_csv(region_file, index_col=0)

    # align index
    common_idx = protein_log.index.intersection(regions.index)
    protein_use = protein_log.loc[common_idx]
    regions = regions.loc[common_idx]
    
    columns = [
        "group", "celltype", "protein",
        "pvalue",
        "mean_boundary", "mean_core",
        "Diff",
        "n_boundary", "n_core"
    ]
    
    All_df = pd.DataFrame(columns=columns)
    
    results = []
    
    df = regions.groupby(["celltype", "region"]).size().reset_index(name="count")
    tumor_df = df[df["celltype"] == "Tumor"]
    if tumor_df.empty:
        valid_types = []
    else:
        tumor_pivot = tumor_df.pivot(
            index="celltype",
            columns="region",
            values="count"
        ).fillna(0)
        valid_types = tumor_pivot[
            (tumor_pivot.get("TumorBoundary", 0) > 20) &
            (tumor_pivot.get("TumorCore", 0) > 20)
        ].index.tolist()
        
        if len(valid_types)>0:
            tumor_I = regions[regions["region"] == "TumorBoundary"]
            tumor_S = regions[regions["region"] == "TumorCore"]
            
            idx_I = tumor_I.index
            idx_S = tumor_S.index
            
            for protein in protein_use.columns:
            
                x = protein_use.loc[idx_I, protein].dropna()
                y = protein_use.loc[idx_S, protein].dropna()
            
                if len(x) < 3 or len(y) < 3:
                    continue
            
                stat, p = mannwhitneyu(x, y, alternative="two-sided")
            
                results.append({
                    "group": "Tumor",
                    "celltype": "Tumor",
                    "protein": protein,
                    "pvalue": p,
                    "mean_boundary": x.mean(),   # boundary
                    "mean_core": y.mean(),   # core
                    "Diff": x.mean()-y.mean(),
                    "n_boundary": len(x),
                    "n_core": len(y)
                })
            res_df = pd.DataFrame(results)
            res_df["p_adj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
            All_df = pd.concat([All_df, res_df], ignore_index=True)
    
    df2 = df[df["celltype"] != "Tumor"]  
    if df2.empty:
        valid_types = [] 
    else:
        pivot = df2.pivot(index="celltype", columns="region", values="count").fillna(0)  
        valid_types = pivot[
            (pivot.get("StromalBoundary", 0) > 20) &
            (pivot.get("StromalCore", 0) > 20)
        ].index.tolist()
        
        for ct in valid_types:
            results = []
            stromal_I = regions[
                (regions["celltype"] == ct) &
                (regions["region"] == "StromalBoundary")
            ]
        
            stromal_S = regions[
                (regions["celltype"] == ct) &
                (regions["region"] == "StromalCore")
            ]
        
            idx_I = stromal_I.index
            idx_S = stromal_S.index
        
            for protein in protein_use.columns:
                x = protein_use.loc[idx_I, protein].dropna()
                y = protein_use.loc[idx_S, protein].dropna()
                if len(x) < 3 or len(y) < 3:
                    continue
                stat, p = mannwhitneyu(x, y, alternative="two-sided")
                results.append({               
                    "group": "Stromal",
                    "celltype": ct,
                    "protein": protein,
                    "pvalue": p,
                    "mean_boundary": x.mean(),   # boundary
                    "mean_core": y.mean(),   # core
                    "Diff": x.mean()-y.mean(),
                    "n_boundary": len(x),
                    "n_core": len(y)
                })
    
            res_df = pd.DataFrame(results)
            res_df["p_adj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
            All_df = pd.concat([All_df, res_df], ignore_index=True)
    
    All_df.to_csv(os.path.join(output, "protein_all_results.csv"),index=False)



# =======================
# CLI
# =======================
def main():
    parser = argparse.ArgumentParser(
        description="Differential protein markers in interface"
    )
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    try:
        sptumor(
            output=args.output,
        )

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()