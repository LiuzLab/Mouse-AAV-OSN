#!/usr/bin/bash

cellranger mkgtf \
    ../../required_files/reference/gencode.vM34.primary_assembly.annotation.gtf \
    ../../required_files/reference/gencode.vM34.primary_assembly.annotation_filtered.gtf \
    --attribute=gene_biotype:protein_coding
