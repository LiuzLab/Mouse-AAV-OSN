# Mouse-AAV-OSN [![DOI](https://zenodo.org/badge/836991156.svg)](https://zenodo.org/doi/10.5281/zenodo.13376699)
Computational analysis for "Comparative Analysis of AAV Serotypes for the Transduction of Olfactory Sensory Neurons" by Belfort and Jia et al. 2024

# GEO Submission Scripts Documentation: `src/snrna/utils`

This document describes three Python scripts designed to handle AnnData objects for Gene Expression Omnibus (GEO) submissions: a preprocessing script, an AnnData to GEO conversion script, and a GEO to AnnData reconstruction script.

## Requirements
`Python`
`pandas`
`anndata`
`scipy`
`numpy`

## 1. Preprocessing Script: `preprocess_geo_files.py`

### Purpose
Standardizes the file structure of GEO submission files by removing prefixes and organizing files into a consistent format.

### Input
- `input_dir`: Path to the directory containing all screen folders and final object folders.

### Output
- Standardized file structure in the specified output directory.

### Functionality
- Processes each subdirectory in the input directory.
- Removes prefixes from filenames.
- Copies files to a new directory structure with standardized names.

### Usage
```bash
python preprocess_geo_files.py /path/to/input/directory /path/to/output/directory
```

## 2. AnnData to GEO Conversion Script: `anndata_to_geo.py`

### Purpose
Converts an AnnData object to the format required for GEO submission.

### Input
- `input_file`: Path to the input .h5ad file.
- `output_dir`: Path to the output directory for GEO submission files.
- `--prefix` (optional): Prefix to add to output files.
- `--chunk_size` (optional): Chunk size for splitting large matrices (default: 1,000,000).

### Output
- GEO submission files in the specified output directory, including:
  - Cell metadata (obs)
  - Gene metadata (var)
  - Count matrix (in chunks)
  - Layer data (if present)
  - Raw data (if present)
  - README file

### Functionality
- Reads an AnnData object.
- Saves metadata, count matrix, layers, and raw data (if present) in GEO-compatible format.
- Splits large matrices into chunks.
- Generates a README file describing the dataset.

### Usage
```bash
python anndata_to_geo.py input.h5ad /path/to/output/directory --prefix optional_prefix_ --chunk_size 500000
```

## 3. GEO to AnnData Reconstruction Script: `geo_to_anndata.py`

### Purpose
Reconstructs an AnnData object from GEO submission files.

### Input
- `input_dir`: Path to the directory containing GEO submission files.
- `output_file`: Path to save the reconstructed .h5ad file.

### Output
- Reconstructed AnnData object saved as an .h5ad file.

### Functionality
- Automatically detects file prefixes.
- Reads cell metadata, gene metadata, and count matrix.
- Reconstructs layers and raw data if present.
- Converts matrices to sparse format if they contain more than 50% zero values.
- Creates and saves an AnnData object.

### Usage
```bash
python geo_to_anndata.py /path/to/geo/submission/files /path/to/output/reconstructed_file.h5ad
```

## Workflow Example

1. Standardize your file structure:
   ```bash
   python preprocess_geo_files.py /original/files /standardized/files
   ```

2. Convert AnnData to GEO format (if needed):
   ```bash
   python anndata_to_geo.py input.h5ad /geo/submission/files --prefix screen1_
   ```

3. Reconstruct AnnData from GEO files:
   ```bash
   python geo_to_anndata.py /geo/submission/files reconstructed_data.h5ad
   ```

Remember to adjust file paths and options according to your specific dataset and requirements.
