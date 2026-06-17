#!/usr/bin/env python3

"""
SpaCET T Cell Interaction Analysis Pipeline
Analyzes cell-cell interactions using stlearn with focus on T cell interactions

Usage: python sTcellinteraction.py --rawpath <visium_directory> --region_file <region_annotation_file> --output <output_directory> [--species <species>]
"""

import os
import sys
import argparse
import warnings
import random
import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st
import anndata as ad
from typing import List, Dict, Optional, Tuple

# Suppress warnings
warnings.filterwarnings("ignore")

#===============================================================================
# Function: sTcellinteraction
# Description: Perform cell-cell interaction analysis using stlearn with 
#              focus on T cell related ligand-receptor pairs
# Parameters:
#   rawpath - Path to 10X Visium raw data directory
#   RegionAnnoFile - Path to region annotation CSV file
#   OutputFold - Output directory path
#   species - Species ('human' or 'mouse'), default: 'human'
#   min_spots - Minimum number of spots for interaction analysis, default: 20
#   n_pairs - Number of ligand-receptor pairs to test, default: 2000
#   n_cpus - Number of CPUs for parallel processing, default: 5
#   seed - Random seed for reproducibility, default: 123
#===============================================================================

def sTcellinteraction(rawpath: str, 
                      RegionAnnoFile: str, 
                      OutputFold: str, 
                      species: str = 'human',
                      min_spots: int = 20,
                      n_pairs: int = 2000,
                      n_cpus: int = 5,
                      seed: int = 123) -> None:
    """
    Main function for T cell interaction analysis using stlearn
    """
    
    # Set random seeds for reproducibility
    np.random.seed(seed)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    
    # Print start message
    print("Starting SpaCET T Cell Interaction Analysis")
    
    #-----------------------------------------------------------------------------
    # Step 1: Load ligand-receptor pairs
    #-----------------------------------------------------------------------------
    print("[1/6] Loading ligand-receptor pairs...")
    try:
        lr_pairs = st.tl.cci.load_lrs(['connectomeDB2020_put'], species=species)
        lr_pairs = np.unique(lr_pairs)
        print(f"  - Loaded {len(lr_pairs)} unique ligand-receptor pairs")
    except Exception as e:
        print(f"  - ERROR loading LR pairs: {e}")
        raise
    
    #-----------------------------------------------------------------------------
    # Step 2: Load and preprocess Visium data
    #-----------------------------------------------------------------------------
    print("[2/6] Loading and preprocessing Visium data...")
    
    # Check if rawpath exists
    if not os.path.exists(rawpath):
        raise FileNotFoundError(f"Raw data directory not found: {rawpath}")
    
    # Load Visium data
    try:
        adata = sc.read_visium(rawpath,library_id="slice1")
        print(f"  - Loaded {adata.n_obs} spots and {adata.n_vars} genes")
    except Exception as e:
        try:
            adata = sc.read_h5ad(OutputFold + "/adata.h5ad")
        except Exception as e:
            try:
                expr = pd.read_csv(rawpath + "/expression.csv", index_col=0)
                meta = pd.read_csv(rawpath + "/metadata.csv", index_col=0)
                meta = meta.loc[expr.columns]
                adata = ad.AnnData(X=expr.T)
                adata.obs = meta
                adata.obsm["spatial"] = meta[["X", "Y"]].values
                adata.obs.columns = ["array_row","array_col"]
            except Exception as e:
                print(f"  - ERROR loading Visium data: {e}")
                raise
    
    adata.var_names = adata.var_names.str.replace("_", "-", regex=False)
    # Make variable names unique
    adata.var_names_make_unique()
    
    # Filter genes
    st.pp.filter_genes(adata, min_cells=3)
    print(f"  - After filtering: {adata.n_obs} spots and {adata.n_vars} genes")
    
    # Normalize total counts
    st.pp.normalize_total(adata)
    print("  - Normalization completed")
    
    # Rename columns for stlearn compatibility
    if 'array_row' in adata.obs.columns and 'array_col' in adata.obs.columns:
        adata.obs.rename(columns={'array_row': 'imagerow', 'array_col': 'imagecol'}, inplace=True)
        print("  - Column names renamed for stlearn compatibility")
    
    #-----------------------------------------------------------------------------
    # Step 3: Load region annotations
    #-----------------------------------------------------------------------------
    print("[3/6] Loading region annotations...")
    
    if not os.path.exists(RegionAnnoFile):
        raise FileNotFoundError(f"Region annotation file not found: {RegionAnnoFile}")
    
    try:
        RegionAnno = pd.read_csv(RegionAnnoFile, index_col=0)
        print(f"  - Loaded annotation for {len(RegionAnno)} spots")
        
        # Check if 'Interface' column exists
        if 'Interface' not in RegionAnno.columns:
            raise ValueError("Region annotation file must contain 'Interface' column")
        
        interface_types = RegionAnno['Interface'].unique()
        print(f"  - Found interface types: {list(interface_types)}")
        
    except Exception as e:
        print(f"  - ERROR loading region annotations: {e}")
        raise
   
    
    #-----------------------------------------------------------------------------
    # Step 4: Perform interaction analysis for each region
    #-----------------------------------------------------------------------------
    print("[5/6] Performing cell-cell interaction analysis...")
    miss_interface_types = list(set(["Tumor", "Interface", "Stroma"]) - set(interface_types))
    for interface_type in miss_interface_types:
        lr_interaction_score = pd.DataFrame()
        os.makedirs(OutputFold, exist_ok=True)
        # Save to CSV
        output_file = os.path.join(OutputFold, f'lr_interaction_score.{interface_type}.csv')
        lr_interaction_score.to_csv(output_file)
        print(f"    Saved to: {output_file}")

    for interface_type in interface_types:
        print(f"\n  >>> Processing region: {interface_type}")
        
        # Get cells for current region
        cells = RegionAnno[RegionAnno['Interface'] == interface_type].index.tolist()
        
        # Create subset for current region
        adata_subset = adata[cells, :].copy()
        
        # Run CCI analysis
        print(f"    Running CCI analysis with {n_pairs} LR pairs...")
        try:
            st.tl.cci.run(
                adata_subset, 
                lr_pairs,
                min_spots=min_spots, 
                distance=0, 
                n_pairs=n_pairs, 
                n_cpus=n_cpus
            )
            print("    CCI analysis completed")
            
            # Extract interaction scores
            lr_interaction_score = pd.DataFrame(
                adata_subset.obsm["lr_scores"],
                index=adata_subset.obs.index,
                columns=adata_subset.uns["lr_summary"].index
            )
        except Exception as e:
            print(f"    ERROR in CCI analysis: {e}")
            lr_interaction_score = pd.DataFrame(index=adata_subset.obs.index)
        os.makedirs(OutputFold, exist_ok=True)
        # Save to CSV
        output_file = os.path.join(OutputFold, f'lr_interaction_score.{interface_type}.csv')
        lr_interaction_score.to_csv(output_file)
        print(f"    Saved to: {output_file}")

#===============================================================================
# Command line interface
#===============================================================================

def main():
    """Main function to parse command line arguments and run analysis"""
    
    parser = argparse.ArgumentParser(
        description='SpaCET T Cell Interaction Analysis Pipeline using stlearn',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sTcellinteraction.py --rawpath /path/to/visium_data \\
                              --region_file /path/to/region_annotation.csv \\
                              --output /path/to/output \\
                              --species human
  
  python sTcellinteraction.py -r /data/visium -a annotations.csv -o results -s mouse -c 8
        """
    )
    
    parser.add_argument('-r', '--rawpath', type=str, required=True,
                        help='Path to 10X Visium raw data directory (required)')
    
    parser.add_argument('-a', '--region_file', type=str, required=True,
                        help='Path to region annotation CSV file with "Interface" column (required)')
    
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output directory path (required)')
    
    parser.add_argument('-s', '--species', type=str, default='human',
                        choices=['human', 'mouse'],
                        help='Species: human or mouse (default: human)')
    
    parser.add_argument('-m', '--min_spots', type=int, default=20,
                        help='Minimum number of spots for interaction analysis (default: 20)')
    
    parser.add_argument('-p', '--n_pairs', type=int, default=2000,
                        help='Number of ligand-receptor pairs to test (default: 2000)')
    
    parser.add_argument('-c', '--n_cpus', type=int, default=5,
                        help='Number of CPUs for parallel processing (default: 5)')
    
    parser.add_argument('--seed', type=int, default=123,
                        help='Random seed for reproducibility (default: 123)')
    
    args = parser.parse_args()
    
    # Run analysis
    try:
        sTcellinteraction(
            rawpath=args.rawpath,
            RegionAnnoFile=args.region_file,
            OutputFold=args.output,
            species=args.species,
            min_spots=args.min_spots,
            n_pairs=args.n_pairs,
            n_cpus=args.n_cpus,
            seed=args.seed
        )
    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()