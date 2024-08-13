#!/usr/bin/bash

# set paths
req_files="../../required_files/custom_ref/aav_ref"
genome="${req_files}/starsolo"
fasta="${req_files}/fasta/genome.fa"
gtf="${req_files}/genes/genes.gtf"

mkdir -p "$genome"

#STAR generate reference

STAR --runMode genomeGenerate \
	--runThreadN 32 \
	--genomeDir "$genome" \
	--genomeFastaFiles "$fasta" \
	--sjdbGTFfile "$gtf" \
	--genomeSAsparseD 3

