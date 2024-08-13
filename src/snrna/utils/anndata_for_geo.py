import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gzip
import os
import argparse

def prepare_anndata_for_geo(input_file, output_dir, chunk_size=1000000):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the AnnData file
    print(f"Reading {input_file}...")
    adata = anndata.read_h5ad(input_file)
    
    output_prefix = os.path.join(output_dir, "geo_submission")
    
    # Save obs (cell metadata)
    print("Saving cell metadata...")
    with gzip.open(f"{output_prefix}_obs.txt.gz", 'wt') as f:
        adata.obs.to_csv(f, sep='\t')
    
    # Save var (gene metadata)
    print("Saving gene metadata...")
    with gzip.open(f"{output_prefix}_var.txt.gz", 'wt') as f:
        adata.var.to_csv(f, sep='\t')
    
    # Save X (count matrix) in chunks
    print("Saving count matrix...")
    if sp.issparse(adata.X):
        X = adata.X.tocsr()
    else:
        X = adata.X
    
    n_cells, n_genes = X.shape
    for i in range(0, n_cells, chunk_size):
        print(f"Processing cells {i} to {min(i+chunk_size, n_cells)}...")
        chunk = X[i:i+chunk_size]
        if sp.issparse(chunk):
            chunk = chunk.toarray()
        df = pd.DataFrame(chunk, index=adata.obs_names[i:i+chunk_size], columns=adata.var_names)
        with gzip.open(f"{output_prefix}_counts_chunk_{i//chunk_size}.txt.gz", 'wt') as f:
            df.to_csv(f, sep='\t')
    
    # Save layers (if present)
    for layer_name, layer_data in adata.layers.items():
        print(f"Saving layer: {layer_name}...")
        if sp.issparse(layer_data):
            layer_data = layer_data.tocsr()
        for i in range(0, n_cells, chunk_size):
            print(f"Processing cells {i} to {min(i+chunk_size, n_cells)}...")
            chunk = layer_data[i:i+chunk_size]
            if sp.issparse(chunk):
                chunk = chunk.toarray()
            df = pd.DataFrame(chunk, index=adata.obs_names[i:i+chunk_size], columns=adata.var_names)
            with gzip.open(f"{output_prefix}_{layer_name}_chunk_{i//chunk_size}.txt.gz", 'wt') as f:
                df.to_csv(f, sep='\t')

    # Create a README file
    print("Creating README file...")
    with open(f"{output_prefix}_README.txt", 'w') as f:
        f.write("Dataset prepared for GEO submission\n\n")
        f.write(f"Original file: {input_file}\n")
        f.write(f"Number of cells: {n_cells}\n")
        f.write(f"Number of genes: {n_genes}\n")
        f.write(f"Files:\n")
        f.write(f"- geo_submission_obs.txt.gz: Cell metadata\n")
        f.write(f"- geo_submission_var.txt.gz: Gene metadata\n")
        f.write(f"- geo_submission_counts_chunk_*.txt.gz: Count matrix (split into chunks)\n")
        for layer_name in adata.layers.keys():
            f.write(f"- geo_submission_{layer_name}_chunk_*.txt.gz: {layer_name} layer data (split into chunks)\n")

    print("Conversion completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert AnnData to GEO submission format")
    parser.add_argument("input_file", help="Path to the input .h5ad file")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Chunk size for splitting large matrices (default: 1,000,000)")
    
    args = parser.parse_args()
    
    prepare_anndata_for_geo(args.input_file, args.output_dir, args.chunk_size)
