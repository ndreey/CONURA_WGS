#!/bin/bash

# Define the input and output files
input_file=$1
output_file=$2

# Initialize the output file and add the header
echo "ACCESSION,SOURCE,HOST" > "$output_file"

# Use awk to parse and extract the desired information
awk '
BEGIN {
    FS = "\n"
    RS = ""
    OFS = ","
}

{
    accession = ""
    source = ""
    host = ""
    
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^ACCESSION/) {
            split($i, arr, " ")
            accession = arr[2]
        }
        if ($i ~ /^SOURCE/) {
            source = $i
            sub(/^SOURCE\s+/, "", source)
        }
        if ($i ~ /\/host=/) {
            match($i, /\/host="([^"]+)"/, arr)
            host = arr[1]
        }
    }
    
    if (accession && source && host) {
        print accession, source, host
    }
}
' "$input_file" >> "$output_file"

echo "CSV file has been created: $output_file"
