#!/usr/bin/python
# run this in cellbender environment
import os
import shutil
import subprocess

# Define the base path and the versions
base_path = '/home/johnathan/projects/arrenkiel_scrnaseq/align_count/output/starsolo'
versions = ['ASAP-screen-1', 'ASAP-screen-2', 'ASAP-screen-3', 'ASAP-screen-4']

# Loop through each version and perform the copying, renaming, and running CellBender
for version in versions:
    raw_path = os.path.join(base_path, version, 'Solo.out', 'Velocyto', 'raw')
    
    # Define the paths for the new spliced and unspliced directories
    spliced_path = os.path.join(raw_path, 'spliced')
    unspliced_path = os.path.join(raw_path, 'unspliced')
    
    # Create the spliced and unspliced directories if they don't exist
    os.makedirs(spliced_path, exist_ok=True)
    os.makedirs(unspliced_path, exist_ok=True)
    
    # Define the source files
    spliced_src = os.path.join(raw_path, 'spliced.mtx')
    unspliced_src = os.path.join(raw_path, 'unspliced.mtx')
    barcodes_src = os.path.join(raw_path, 'barcodes.tsv')
    features_src = os.path.join(raw_path, 'features.tsv')
    
    # Define the destination files
    spliced_dst = os.path.join(spliced_path, 'matrix.mtx')
    unspliced_dst = os.path.join(unspliced_path, 'matrix.mtx')
    barcodes_dst_spliced = os.path.join(spliced_path, 'barcodes.tsv')
    barcodes_dst_unspliced = os.path.join(unspliced_path, 'barcodes.tsv')
    features_dst_spliced = os.path.join(spliced_path, 'genes.tsv')
    features_dst_unspliced = os.path.join(unspliced_path, 'genes.tsv')
    
    # Copy and rename the files
    shutil.copyfile(spliced_src, spliced_dst)
    shutil.copyfile(unspliced_src, unspliced_dst)
    shutil.copyfile(barcodes_src, barcodes_dst_spliced)
    shutil.copyfile(barcodes_src, barcodes_dst_unspliced)
    shutil.copyfile(features_src, features_dst_spliced)
    shutil.copyfile(features_src, features_dst_unspliced)

    # Define the output paths for CellBender
    spliced_output = os.path.join(spliced_path, 'spliced_filt', 'output.lh5')
    unspliced_output = os.path.join(unspliced_path, 'unspliced_filt', 'output.lh5')
    
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(spliced_output), exist_ok=True)
    os.makedirs(os.path.dirname(unspliced_output), exist_ok=True)

    # Run CellBender on spliced matrix
    subprocess.run([
        'cellbender', 'remove-background',
        '--input', spliced_path,
        '--output', spliced_output,
        '--cuda'  # Use GPU if available
    ])

    # Run CellBender on unspliced matrix
    subprocess.run([
        'cellbender', 'remove-background',
        '--input', unspliced_path,
        '--output', unspliced_output,
        '--cuda'  # Use GPU if available
    ])

    print(f"Processed {version}")

print("All versions processed successfully.")
