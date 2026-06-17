#!/usr/bin/env python3

"""
Spatial cell colocalization analysis

Pipeline:
1. Load protein expression and spatial metadata
2. Construct AnnData object
3. Build radius-based spatial graph
4. Perform neighborhood enrichment (Squidpy)
5. Export results for:
   - All cells
   - Stromal region
   - Interface region
"""

import pandas as pd
import numpy as np
import os
import argparse
import sys
import scanpy as sc
import squidpy as sq
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
import random

def set_seed(seed=42):
    os.environ["PYTHONHASHSEED"] = str(seed)

    random.seed(seed)
    np.random.seed(seed)


# =======================
# Create AnnData
# =======================
def create_anndata_from_spatial(
    output: str,
    coord_columns: list,
    region_col: str,
    celltype_col: str,
    radius: float,
    seed,
    library_id: str = "spatial",
):
    """
    Create AnnData object and spatial graph.

    Parameters
    ----------
    output : str
        Output directory
    coord_columns : list
        Coordinate column names, e.g. ['X', 'Y']
    region_col : str
        Column name for region annotation
    celltype_col : str
        Column name for cell type annotation
    radius : float
        Radius for neighbor graph
    """
    set_seed(seed)

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

    # -----------------------
    # Load data
    # -----------------------
    
    counts = pd.read_csv(proteinfile, index_col=0)
    coords = pd.read_csv(region_file, index_col=0)

    # Align indices
    common_idx = counts.index.intersection(coords.index)
    counts = counts.loc[common_idx]
    coords = coords.loc[common_idx]

    # Save merged CSV (optional)
    merged = pd.concat(
        [counts, coords[[celltype_col, region_col] + coord_columns]],
        axis=1
    )
    merged.to_csv(os.path.join(output, "Data.csv"))
    # calculate region-cell ratio
    cross_tab = pd.crosstab(merged[celltype_col], merged[region_col])
    col_sums = cross_tab.sum(axis=0)
    prop_table = cross_tab.div(col_sums, axis=1)
    
    CellsReion = prop_table.reset_index().melt(id_vars=celltype_col,var_name=region_col, value_name='ratio')
    CellsReion.columns = ['cell', 'region', 'ratio']
    CellsReion.to_csv(os.path.join(output, "RegionCell.txt"), sep='\t', index=False)
    # calculate cell-protein ratio
    exclude_cols = [celltype_col, region_col] + coord_columns
    protein_cols = [col for col in merged.columns if col not in exclude_cols]
    protein_expression = merged.groupby(celltype_col)[protein_cols].mean()
    protein_expression_long = protein_expression.reset_index().melt(id_vars=celltype_col,var_name='protein',value_name='expression')
    protein_expression_long.columns = ['cell', 'protein', 'expression']
    protein_expression_long.to_csv(os.path.join(output, "CellProtein.txt"), sep='\t', index=False)

    # -----------------------
    # Build AnnData
    # -----------------------
    adata = sc.AnnData(
        X=counts.values,
        obs=pd.DataFrame(index=common_idx),
        var=pd.DataFrame(index=counts.columns)
    )

    # Add spatial coordinates
    adata.obsm["spatial"] = coords[coord_columns].values

    # Add annotations
    coords.index = coords.index.astype(str)
    adata.obs["region"] = coords[region_col].astype("category")
    adata.obs["celltype"] = coords[celltype_col].astype("category")

    # -----------------------
    # Build spatial graph
    # -----------------------
    coord_array = adata.obsm["spatial"]
    
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(coord_array)
    distances, indices = nbrs.kneighbors(coord_array)
    min_radius = 3 * np.median(distances[:, 1])
    
    nn = NearestNeighbors(radius=min_radius, algorithm="kd_tree")
    nn.fit(coord_array)

    distances, indices = nn.radius_neighbors(coord_array)

    # Connectivity graph
    connectivities = nn.radius_neighbors_graph(
        coord_array, mode="connectivity"
    )

    # Distance matrix (sparse)
    lengths = np.array([len(x) for x in indices])
    rows = np.repeat(np.arange(len(coord_array)), lengths)
    cols = np.concatenate(indices)
    data = np.concatenate(distances)
    
    mask = rows != cols

    rows = rows[mask]
    cols = cols[mask]
    data = data[mask]

    distance_sparse = csr_matrix(
        (data, (rows, cols)),
        shape=(len(coord_array), len(coord_array))
    )

    # Store graph
    adata.obsp["spatial_connectivities"] = connectivities
    adata.obsp["spatial_distances"] = distance_sparse

    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": "spatial_distances",
        "params": {"radius": radius}
    }
    
    adata.uns["spatial"] = {library_id: {}}

    # Save AnnData
    adata.write(os.path.join(output, "adata.h5ad"))

    return adata


# =======================
# Enrichment function
# =======================
def run_enrichment(adata, flag):
    """
    Run neighborhood enrichment and save result.
    """
    n_types = adata.obs["celltype"].dropna().nunique()

    if n_types < 2:
        print(
            f"Skip neighborhood enrichment for {flag}: "
            f"only {n_types} cell type(s) found."
        )
        return pd.DataFrame(
            columns=["Group", "cell1", "cell2", "enrichment_score"]
        )

    sq.gr.nhood_enrichment(adata, cluster_key="celltype")

    celltypes = adata.obs["celltype"].cat.categories
    zscore = adata.uns["celltype_nhood_enrichment"]["zscore"]

    df = pd.DataFrame(zscore, index=celltypes, columns=celltypes)

    long_df = df.stack().reset_index()
    long_df.columns = ["cell1", "cell2", "enrichment_score"]
    long_df.insert(0, "Group", flag)
    
    return long_df


# =======================
# Main pipeline
# =======================
def spcolocalization(
    output,
    coord_columns,
    region_col,
    celltype_col,
    radius,
    seed
):
    set_seed(seed)

    adata = create_anndata_from_spatial(
        output,
        coord_columns,
        region_col,
        celltype_col,
        radius,
        seed
    )

    # -------- All cells --------
    df = run_enrichment(adata,"All")

    # -------- Stromal boundary --------
    adata_SB = adata[adata.obs["region"] == "StromalBoundary"].copy()
    print(f"Stromal boundary cell number: {adata_SB.n_obs}")
    if adata_SB.n_obs > 0:
        df1 = run_enrichment(adata_SB,"StromalBoundary")
        df = pd.concat([df, df1])
        
    # -------- Stromal core --------
    adata_SB = adata[adata.obs["region"] == "StromalCore"].copy()
    print(f"Stromal boundary cell number: {adata_SB.n_obs}")
    if adata_SB.n_obs > 0:
        df1 = run_enrichment(adata_SB,"StromalCore")
        df = pd.concat([df, df1])
    
    df.to_csv(os.path.join(output, "Cellcolocalization.txt"), sep="\t", index=False)


# =======================
# CLI
# =======================
def main():
    parser = argparse.ArgumentParser(
        description="Spatial cell colocalization analysis"
    )

    parser.add_argument("--output", required=True,
                        help="Output directory")

    # Coordinate columns
    parser.add_argument("--x_col", default="X")
    parser.add_argument("--y_col", default="Y")

    # Annotation columns
    parser.add_argument("--region_col", default="region")
    parser.add_argument("--celltype_col", default="celltype")

    # Radius parameter
    parser.add_argument("--radius", type=float, default=80)
    
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    try:
        spcolocalization(
            output=args.output,
            coord_columns=[args.x_col, args.y_col],
            region_col=args.region_col,
            celltype_col=args.celltype_col,
            radius=args.radius,
            seed=args.seed
        )

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()