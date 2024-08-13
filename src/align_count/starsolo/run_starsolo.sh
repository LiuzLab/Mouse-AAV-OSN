#!/bin/bash

# running STAR v2.7.11b

#EmptyDrops_CR Params
nExpectedCells=10000        # Increase to allow more cells/nuclei; default 10k
maxPercentile=0.99       # Increase slightly to retain more droplets; default 0.99
maxMinRatio=10              # Reduce to capture droplets with lower UMI ratios; default 10
indMin=45000               # Lower to include droplets with fewer UMIs; default 45k
indMax=90000              # Raise to include droplets with more UMIs; default 90k
umiMin=500                 # Reduce to retain droplets with fewer UMIs; default 500
umiMinFracMedian=0.01     # Decrease fraction to include droplets below median UMI count; default 0.01
candMaxN=20000             # Increase to evaluate more candidate droplets; default 20k
FDR=0.01                   # Increase to be more lenient; default 0.01
simN=10000                  # Reduce if needed, but balance with accuracy; default 10k

# Set your directories
ref_dir="../../required_files/custom_ref/aav_ref/starsolo"
out_dir="../../output/starsolo"
fastq_dir="../../data/00_fastq" # Directory containing all FASTQ files

# Extract sample names, removing duplicates and empty names
samples=($(ls $fastq_dir | grep -oE '^[^_]+' | uniq))

# Process each sample
for sample_name in "${samples[@]}"; do
    if [ -z "$sample_name" ]; then
        continue # Skip processing if the sample name is empty
    fi

    echo "Processing $sample_name"

    # Construct the output directory path with the correct prefix
    sample_out_dir="${out_dir}/${sample_name}/"
    mkdir -p "$sample_out_dir"

    echo "Generating output in $sample_out_dir"

    # Run STAR Solo with the corrected output prefix
    STAR --genomeDir "$ref_dir" \
        --readFilesIn "${fastq_dir}/${sample_name}_S1_L001_R2_001.fastq.gz" "${fastq_dir}/${sample_name}_S1_L001_R1_001.fastq.gz" \
        --runThreadN 48 \
	--soloCBwhitelist None \
        --soloType CB_UMI_Simple \
        --soloFeatures Gene GeneFull Velocyto \
        --outFileNamePrefix "${sample_out_dir}" \
        --soloCBlen 16 \
	--soloBarcodeReadLength 150 \
        --outFilterScoreMin 30 \
        --soloUMIlen 12 \
        --soloMultiMappers EM \
	--soloCellFilter EmptyDrops_CR $nExpectedCells $maxPercentile $maxMinRatio $indMin $indMax $umiMin $umiMinFracMedian $candMaxN $FDR $simN \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --readFilesCommand zcat
done
