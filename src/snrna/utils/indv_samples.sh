#!/bin/bash

# Define the base directories and script path
base_dir="/home/johnathan/projects/arrenkiel_scrnaseq/final_submission/git"
data_dir="${base_dir}/data"
output_dir="${data_dir}/geo_submission/indv_objects"
script_path="${base_dir}/src/snrna/utils/anndata_for_geo.py"

# Make sure the output directory exists
mkdir -p "${output_dir}"

# Define the function to process each screen
process_screen() {
    screen=$1
    file="${data_dir}/before_qc/${screen}_combined.h5ad"
    screen_output_dir="${output_dir}/${screen}"
    mkdir -p "${screen_output_dir}"
    python "${script_path}" "${file}" "${screen_output_dir}"
}

# Export the function so GNU Parallel can use it
export -f process_screen
export base_dir data_dir output_dir script_path

# Run the process in parallel
parallel -j 16 process_screen ::: screen-1 screen-2 screen-3 screen-4
