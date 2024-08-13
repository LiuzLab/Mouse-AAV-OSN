#!/bin/bash

# Check if an input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 fasta_file"
    exit 1
fi

input_file="$1"

# Use awk to count characters in each sequence, excluding headers
awk '/^>/ {if (seq) print length(seq); print $0; seq=""; next} {seq=seq""$0} END {if (seq) print length(seq)}' "$input_file"
