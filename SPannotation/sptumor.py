#!/usr/bin/env python3

"""
Spatial Tumor-Stromal identification Pipeline
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import sys
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
import random

def set_seed(seed=42):
    os.environ["PYTHONHASHSEED"] = str(seed)

    random.seed(seed)
    np.random.seed(seed)
    
    matplotlib.use("Agg")


def sptumor(
    output,
    tumor_label="Tumor",
    tumor_cutoff=0,
    stromal_cutoff=0,
    radius=None,
    seed=42
):
    set_seed(seed)

    # =======================
    # Check files
    # =======================
    celltype_file   = os.path.join(output, "celltype_annotation_final.txt")
    spatial_file   = os.path.join(output, "spatial.csv")
    for f in [ celltype_file, spatial_file]:
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
        "TumorCore": "#d73027",        # deep red
        "TumorBoundary": "#fc8d59",    # orange
        "StromalBoundary": "#1a9850",  # green
        "StromalCore": "#2b2b2b"      # dark gray/black
    }

    # =======================
    # Process per sample
    # =======================

    coords = data[["X", "Y"]].values
    labels = data["celltype"].values

    # ===== radius neighbors =====
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    min_radius = 3 * np.median(distances[:, 1])
    
    if radius is None:
        radius_use = min_radius
    else:
        radius_use = radius
    
    print(f"min_radius={min_radius:.2f}, radius_use={radius_use:.2f}")
    
    nbrs = NearestNeighbors(radius=radius_use, algorithm='kd_tree').fit(coords)
    distances, indices = nbrs.radius_neighbors(coords)
    
    result = []
    
    for i in range(len(data)):

      cell_type = labels[i]
      neighbor_idx = indices[i]
      neighbor_idx = neighbor_idx[neighbor_idx != i]
  
      neighbor_types = labels[neighbor_idx]
      if len(neighbor_idx) == 0:
          tumor_ratio = 0
          stromal_ratio = 0
      else:
          tumor_ratio = np.mean(neighbor_types == tumor_label)
          stromal_ratio = np.mean(neighbor_types != tumor_label)
  
      if cell_type == tumor_label:
  
          if stromal_ratio > stromal_cutoff:
              result.append("TumorBoundary")
          else:
              result.append("TumorCore")
  
      else:
          if len(neighbor_idx) == 0:
              result.append("StromalCore")
              continue
          if tumor_ratio > tumor_cutoff:
              result.append("StromalBoundary")
          else:
              result.append("StromalCore")

    data["region"] = result
    
    csv_file = os.path.join(output, f"region.csv")
    data.to_csv(csv_file)

    # =======================
    # Plot dot 
    # =======================
    plt.figure(figsize=(6, 6))

    for ct, group in data.groupby("region"):
        plt.scatter(
            group["X"],
            group["Y"],
            c=celltype_colors.get(ct, "grey"),
            s=2,
            alpha=0.8,
            label=ct
        )

    plt.xlabel("X")
    plt.ylabel("Y")

    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().invert_yaxis()

    plt.legend(
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        fontsize=8,
        frameon=False,
        markerscale=3
    )
    svg_file = os.path.join(output, f"tumor.svg")
    plt.savefig(svg_file, bbox_inches='tight')
    plt.close()

    print("All done!")
    
# =======================
# CLI
# =======================
def main():
    parser = argparse.ArgumentParser(
        description="Spatial Tumor-Stromal identification Pipeline"
    )

    parser.add_argument("--output", required=True)
    parser.add_argument("--tumor_label", default="Tumor")
    parser.add_argument("--tumor_cutoff", type=float, default=0)
    parser.add_argument("--stromal_cutoff", type=float, default=0)
    parser.add_argument("--radius", type=int, default=None)
    
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    try:
        sptumor(
            output=args.output,
            tumor_label=args.tumor_label,
            tumor_cutoff=args.tumor_cutoff,
            stromal_cutoff=args.stromal_cutoff,
            radius=args.radius,
            seed=args.seed
        )

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()