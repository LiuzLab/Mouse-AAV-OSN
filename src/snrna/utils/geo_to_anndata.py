import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gzip
import os
import glob
import argparse

def reconstruct_anndata_from_geo(input_dir, output_file):
    input_prefix = os.path.join(input_dir, "geo_submission")
    
    # Read obs (cell metadata)
    print("Reading cell metadata...")
    with gzip.open(f"{input_prefix}_obs.txt.gz", 'rt') as f:
        obs = pd.read_csv(f, sep='\t', index_col=0)
    
    # Read var (gene metadata)
    print("Reading gene metadata...")
    with gzip.open(f"{input_prefix}_var.txt.gz", 'rt') as f:
        var = pd.read_csv(f, sep='\t', index_col=0)
    
    # Read X (count matrix) chunks
    print("Reading count matrix...")
    count_files = sorted(glob.glob(f"{input_prefix}_counts_chunk_*.txt.gz"))
    counts = []
    for file in count_files:
        print(f"Processing {file}...")
        with gzip.open(file, 'rt') as f:
            chunk = pd.read_csv(f, sep='\t', index_col=0)
        counts.append(chunk)
    X = pd.concat(counts)
    
    # Create AnnData object
    print("Creating AnnData object...")
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    
    # Read and add layers (if present)
    layer_files = glob.glob(f"{input_prefix}_*_chunk_0.txt.gz")
    layer_names = [os.path.basename(file).split('_chunk_')[0].replace(f"{os.path.basename(input_prefix)}_", "") 
                   for file in layer_files if "counts" not in file]
    
    for layer_name in layer_names:
        print(f"Reading layer: {layer_name}...")
        layer_chunks = sorted(glob.glob(f"{input_prefix}_{layer_name}_chunk_*.txt.gz"))
        layer_data = []
        for file in layer_chunks:
            print(f"Processing {file}...")
            with gzip.open(file, 'rt') as f:
                chunk = pd.read_csv(f, sep='\t', index_col=0)
            layer_data.append(chunk)
        adata.layers[layer_name] = pd.concat(layer_data)
    
    # Convert to sparse if the original was likely sparse
    if (adata.X == 0).sum() / adata.X.size > 0.5:  # If more than 50% of values are 0
        print("Converting main matrix to sparse format...")
        adata.X = sp.csr_matrix(adata.X)
    
    for layer_name, layer in adata.layers.items():
        if (layer == 0).sum() / layer.size > 0.5:  # If more than 50% of values are 0
            print(f"Converting {layer_name} layer to sparse format...")
            adata.layers[layer_name] = sp.csr_matrix(layer)
    
    # Save the reconstructed AnnData object
    print(f"Saving reconstructed AnnData object to {output_file}...")
    adata.write_h5ad(output_file)
    
    print("Reconstruction completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reconstruct AnnData from GEO submission files")
    parser.add_argument("input_dir", help="Path to the directory containing GEO submission files")
    parser.add_argument("output_file", help="Path to save the reconstructed .h5ad file")
    
    args = parser.parse_args()
    
    reconstruct_anndata_from_geo(args.input_dir, args.output_file)
