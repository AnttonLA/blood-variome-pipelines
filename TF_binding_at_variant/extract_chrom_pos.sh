#!/bin/bash
# This script takes a .tsv with the columns  ID, Chrom, Pos, OA and EA outputs the string 'Chrom:Pos' for each row.

# Check if the input TSV file path is provided as an argument
if [ -z "$1" ]; then
    echo "Error: Please provide the path of the input TSV file as an argument."
    exit 1
fi

# Read the input TSV file, skip the first row, and extract Chrom and Pos columns
awk -F '\t' 'NR>1 {sub(/^chr/, "", $2); print $2 ":" $3}' "$1"
