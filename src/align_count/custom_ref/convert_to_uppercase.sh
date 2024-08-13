#!/bin/bash

# Check if the input file was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file="$1"
output_file="$2"

# Use sed to convert lowercase to uppercase
sed 's/[a-z]/\U&/g' "$input_file" > "$output_file"

echo "Conversion complete. Uppercase sequences saved to '$output_file'."
