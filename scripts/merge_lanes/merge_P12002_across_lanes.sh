#!/bin/bash

#SBATCH --job-name P12002_merge
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 4
#SBATCH -t 01:30:00
#SBATCH --output=SLURM-%j-P12002_merge.out
#SBATCH --error=SLURM-%j-P12002_merge.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Each sample in the sequence project ($dir) are found in the .lst file.
# This script reads through each .lst and generates a concatenated fastq 
# at user specified destination ($out).
# OBS, each run over runs previous runs.


# Directories to find the project .lst files (dir) and where we will generate the
# new fastq.gz files (out).
dir="P12002"
out="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/00-RAW"

# Path from where .lst files path from.
data="/crex/proj/snic2020-6-222/Projects/Tconura/data/WGS/rawdata/${dir}"

# Finds the .lst files in dir that follows <sample_id>.lst
find "$data" -name "${dir}_*.lst" | while read lst; do

    sample=$(basename $lst .lst)
    echo "Processing files for sample: $sample"
    
    # Out files for merged R1 and R2
    R1_out="${out}/${sample}_R1.fastq.gz"
    R2_out="${out}/${sample}_R2.fastq.gz"
    
    # Generate empty fastq.gz files for R1 and R2 to avoid duplication by
    # generating and setting the files to 0 bytes. Thus overwriting.
    truncate -s 0 "$R1_out" "$R2_out"
    
    # Lets loop through the files in the .lst
    while read -r line; do
        # Only handle lines that end with .fastq.gz using pattern (=~)
        if [[ "$line" =~ ".fastq.gz" ]]; then
            fq_path="${data}/${line}"

            # Determines if its R1 or R2 fastq file
            if [[ "$line" =~ "_R1_" ]]; then
                cat "$fq_path" >> "$R1_out"
            
            elif [[ "$line" =~ "_R2_" ]]; then
                cat "$fq_path" >> "$R2_out"
            fi        
        fi
    done < "$lst"
done

