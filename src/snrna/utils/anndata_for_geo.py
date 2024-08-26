import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gzip
import os
import argparse

def save_matrix(matrix, obs_names, var_names, file_prefix, chunk_size):
    if sp.issparse(matrix):
        matrix = matrix.tocsr()
    
    n_cells, n_genes = matrix.shape
    for i in range(0, n_cells, chunk_size):
        print(f"Processing cells {i} to {min(i+chunk_size, n_cells)}...")
        chunk = matrix[i:i+chunk_size]
        if sp.issparse(chunk):
            chunk = chunk.toarray()
        df = pd.DataFrame(chunk, index=obs_names[i:i+chunk_size], columns=var_names)
        with gzip.open(f"{file_prefix}_chunk_{i//chunk_size}.txt.gz", 'wt') as f:
            df.to_csv(f, sep='\t')

def prepare_anndata_for_geo(input_file, output_dir, chunk_size=1000000):
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Reading {input_file}...")
    adata = anndata.read_h5ad(input_file)
    
    output_prefix = os.path.join(output_dir, "geo_submission")
    
    print("Saving cell metadata...")
    with gzip.open(f"{output_prefix}_obs.txt.gz", 'wt') as f:
        adata.obs.to_csv(f, sep='\t')
    
    print("Saving gene metadata...")
    with gzip.open(f"{output_prefix}_var.txt.gz", 'wt') as f:
        adata.var.to_csv(f, sep='\t')
    
    print("Saving count matrix...")
    save_matrix(adata.X, adata.obs_names, adata.var_names, f"{output_prefix}_counts", chunk_size)
    
    for layer_name, layer_data in adata.layers.items():
        print(f"Saving layer: {layer_name}...")
        save_matrix(layer_data, adata.obs_names, adata.var_names, f"{output_prefix}_{layer_name}", chunk_size)

    if adata.raw is not None:
        print("Saving raw data...")
        raw_prefix = f"{output_prefix}_raw"
        
        print("Saving raw var metadata...")
        with gzip.open(f"{raw_prefix}_var.txt.gz", 'wt') as f:
            adata.raw.var.to_csv(f, sep='\t')
        
        print("Saving raw count matrix...")
        save_matrix(adata.raw.X, adata.obs_names, adata.raw.var_names, f"{raw_prefix}_counts", chunk_size)

    print("Creating README file...")
    with open(f"{output_prefix}_README.txt", 'w') as f:
        f.write("Dataset prepared for GEO submission\n\n")
        f.write(f"Original file: {input_file}\n")
        f.write(f"Number of cells: {adata.n_obs}\n")
        f.write(f"Number of genes: {adata.n_vars}\n")
        f.write(f"Files:\n")
        f.write(f"- geo_submission_obs.txt.gz: Cell metadata\n")
        f.write(f"- geo_submission_var.txt.gz: Gene metadata\n")
        f.write(f"- geo_submission_counts_chunk_*.txt.gz: Count matrix (split into chunks)\n")
        for layer_name in adata.layers.keys():
            f.write(f"- geo_submission_{layer_name}_chunk_*.txt.gz: {layer_name} layer data (split into chunks)\n")
        if adata.raw is not None:
            f.write(f"- geo_submission_raw_var.txt.gz: Raw gene metadata\n")
            f.write(f"- geo_submission_raw_counts_chunk_*.txt.gz: Raw count matrix (split into chunks)\n")

    print("Conversion completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert AnnData to GEO submission format")
    parser.add_argument("input_file", help="Path to the input .h5ad file")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Chunk size for splitting large matrices (default: 1,000,000)")
    
    args = parser.parse_args()
    
    prepare_anndata_for_geo(args.input_file, args.output_dir, args.chunk_size)