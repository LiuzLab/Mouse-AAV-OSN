#!/usr/bin/bash

#make custom ref with aav sequences
cellranger mkref --genome=custom_ref \
		--fasta=GRCm39.primary_assembly.genome_aav.fa \
		--genes=gencode.vM34.primary_assembly.annotation_aav.gtf
