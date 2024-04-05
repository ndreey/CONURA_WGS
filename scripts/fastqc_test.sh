#!/bin/bash

#SBATCH --job-name test_FastQC
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 8
#SBATCH -t 00:15:00
#SBATCH --output=SLURM-%j.out
#SBATCH --error=SLURM-%j.err


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load FastQC/0.11.9

# P12002 forward reads
fastqc 00-RAW/P12002_101_R1.fastq.gz 00-RAW/P12002_102_R1.fastq.gz 00-RAW/P12002_103_R1.fastq.gz \
    -t 8 \
    --outdir 01-QC/fastqc_raw_R1_P12002



# End time and date
echo "$(date)       [End]"