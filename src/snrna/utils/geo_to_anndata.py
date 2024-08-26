import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gzip
import os
import glob
import argparse

def read_matrix(file_prefix, obs_names, var_names):
    chunk_files = sorted(glob.glob(f"{file_prefix}_chunk_*.txt.gz"))
    chunks = []
    for file in chunk_files:
        print(f"Processing {file}...")
        with gzip.open(file, 'rt') as f:
            chunk = pd.read_csv(f, sep='\t', index_col=0)
        chunks.append(chunk)
    matrix = pd.concat(chunks)
    return matrix

def reconstruct_anndata_from_geo(input_dir, output_file):
    input_prefix = os.path.join(input_dir, "geo_submission")
    
    print("Reading cell metadata...")
    with gzip.open(f"{input_prefix}_obs.txt.gz", 'rt') as f:
        obs = pd.read_csv(f, sep='\t', index_col=0)
    
    print("Reading gene metadata...")
    with gzip.open(f"{input_prefix}_var.txt.gz", 'rt') as f:
        var = pd.read_csv(f, sep='\t', index_col=0)
    
    print("Reading count matrix...")
    X = read_matrix(f"{input_prefix}_counts", obs.index, var.index)
    
    print("Creating AnnData object...")
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    
    # Read and add layers (if present)
    layer_files = glob.glob(f"{input_prefix}_*_chunk_0.txt.gz")
    layer_names = [os.path.basename(file).split('_chunk_')[0].replace(f"{os.path.basename(input_prefix)}_", "") 
                   for file in layer_files if "counts" not in file and "raw" not in file]
    
    for layer_name in layer_names:
        print(f"Reading layer: {layer_name}...")
        adata.layers[layer_name] = read_matrix(f"{input_prefix}_{layer_name}", obs.index, var.index)
    
    # Read raw data (if present)
    if os.path.exists(f"{input_prefix}_raw_var.txt.gz"):
        print("Reading raw data...")
        with gzip.open(f"{input_prefix}_raw_var.txt.gz", 'rt') as f:
            raw_var = pd.read_csv(f, sep='\t', index_col=0)
        
        print("Saving raw data to anndata object...")
        raw_X = read_matrix(f"{input_prefix}_raw_counts", obs.index, raw_var.index)
        adata.raw = anndata.AnnData(X=raw_X, obs = obs, var=raw_var)
    
    # Convert to sparse if more than 50% of values are 0
#    print("Converting matrices to sparse format if needed...")
#    for attr in ['X'] + list(adata.layers.keys()):
#        matrix = adata.X if attr == 'X' else adata.layers[attr]
#        if (matrix == 0).sum() / matrix.size > 0.5:
#            print(f"Converting {attr} to sparse format...")
#            if attr == 'X':
#                adata.X = sp.csr_matrix(matrix)
#            else:
#                adata.layers[attr] = sp.csr_matrix(matrix)
    
#    if adata.raw is not None and (adata.raw.X == 0).sum() / adata.raw.X.size > 0.5:
#        print("Converting raw matrix to sparse format...")
#        adata.raw.X = sp.csr_matrix(adata.raw.X)
    
    print(f"Saving reconstructed AnnData object to {output_file}...")
    adata.write_h5ad(output_file)
    
    print("Reconstruction completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reconstruct AnnData from GEO submission files")
    parser.add_argument("input_dir", help="Path to the directory containing GEO submission files")
    parser.add_argument("output_file", help="Path to save the reconstructed .h5ad file")
    
    args = parser.parse_args()
    
    reconstruct_anndata_from_geo(args.input_dir, args.output_file)
