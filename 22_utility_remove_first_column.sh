#!/bin/bash

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_file [output_file]"
    echo "If output_file is not specified, will create input_file_minus_1st_col.txt"
    exit 1
fi

input_file="$1"

# Set output filename
if [ -z "$2" ]; then
    # Remove extension if exists
    base_name="${input_file%.*}"
    extension="${input_file##*.}"

    # Handle files without extension
    if [ "$base_name" = "$input_file" ]; then
        output_file="${input_file}_minus_1st_col"
    else
        output_file="${base_name}_minus_1st_col.${extension}"
    fi
else
    output_file="$2"
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found"
    exit 1
fi

# Detect delimiter (tab, comma, or space)
first_line=$(head -n 1 "$input_file")
delimiter=" "

if [[ "$first_line" == *$'\t'* ]]; then
    delimiter=$'\t'
elif [[ "$first_line" == *,* ]]; then
    delimiter=','
fi

# Remove first column and save to new file
if [ "$delimiter" == $'\t' ]; then
    cut --complement -f 1 -d $'\t' "$input_file" > "$output_file"
elif [ "$delimiter" == ',' ]; then
    cut --complement -f 1 -d ',' "$input_file" > "$output_file"
else
    awk '{$1=""; sub(/^ /, ""); print}' "$input_file" > "$output_file"
fi

echo "First column removed. Output saved to: $output_file"
                                                             
