#!/usr/bin/env python3

"""
Spatial Cell Type Visualization Pipeline
"""

# =======================
# Imports
# =======================
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import sys


# =======================
# Core function
# =======================
def spspatialplot(
    output
):

    # =======================
    # Check files
    # =======================
    celltype_file   = os.path.join(output, "celltype_annotation_final.txt")
    spatial_file   = os.path.join(output, "spatial.csv")
    for f in [celltype_file, spatial_file]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    # =======================
    # Load data
    # =======================
    celltypes = pd.read_csv(celltype_file, index_col=0)
    spatial = pd.read_csv(spatial_file, index_col=0)
    
    common_idx = celltypes.index.intersection(spatial.index)
    celltypes = celltypes.loc[common_idx]
    spatial = spatial.loc[common_idx]
    
    data = celltypes.join(spatial)

    print("Merged shape:", data.shape)

    # =======================
    # Colors
    # =======================
    celltype_colors = {
      "Tumor": "#e41a1c",          # red
      "Stromal": "#4daf4a",        # green
  
      "Macrophage": "#984ea3",     # purple
      "Monocyte": "#af7ac5",       # light purple
  
      "Dendritic": "#a65628",      # brown-orange
      "Neutrophils": "#f39c12",    # orange
      "Mast": "#ff7f00",
  
      "T": "#377eb8",              # blue (pan T)
      "CD4T": "#1f78b4",           # blue
      "CD8T": "#00bcd4",           # cyan-ish
      "NK": "#00cfd1",             # teal
      "B": "#f781bf",              # pink
      "Plasma": "#ff69b4",         # hot pink
  
      "Endothelial": "#6d4c41",    # brown
  
      "Immune": "#7f8c8d",         # gray (umbrella class)
    }


    # =======================
    # Plot per sample
    # =======================

    plt.figure(figsize=(6, 6))

    # plot according to celltype group
    for ct, group in data.groupby("celltype"):
        plt.scatter(
            group["X"],
            group["Y"],
            c=celltype_colors.get(ct, "black"),
            s=2,
            alpha=0.8,
            label=ct
        )
    plt.xlabel("X")
    plt.ylabel("Y")

    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().invert_yaxis()

    # legned
    plt.legend(
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        fontsize=8,
        frameon=False,
        markerscale=3
    )

    # =======================
    # Save
    # =======================
    svg_file = os.path.join(output, f"spatial.svg")

    plt.savefig(svg_file, bbox_inches='tight')

    plt.close()

    print("All done!")


# =======================
# CLI
# =======================
def main():
    parser = argparse.ArgumentParser(
        description="Spatial Cell Type Visualization Pipeline"
    )

    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    try:
        spspatialplot(
            output=args.output
        )

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()