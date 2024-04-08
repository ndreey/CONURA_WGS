#!/bin/bash

#SBATCH --job-name fastp_P12002
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 8
#SBATCH -t 03:15:00
#SBATCH --output=SLURM-%j-P12002_R1_QC_RAW.out
#SBATCH --error=SLURM-%j-P12002_R1_QC_RAW.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load fastp/0.23.4

# Directory with raw reads
raw="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/00-RAW"

while read -r sample r1 r2; do

    R1="$raw/$r1"
    R2="$raw/$r2"

    fastp --in1 $R1 --in2 $R2 --out1 --out2


done < doc/raw_samples.txt






# End time and date
echo "$(date)       [End]"