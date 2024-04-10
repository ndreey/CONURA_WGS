#!/bin/bash

#SBATCH --job-name P14052_R1_QC_raw
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 5
#SBATCH -t 01:15:00
#SBATCH --output=SLURM-%j-P14052_R1_QC_RAW.out
#SBATCH --error=SLURM-%j-P14052_R1_QC_RAW.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load FastQC/0.11.9 MultiQC/1.12

# Variables for R1/R2 and sequence project
proj="P14052"
R="R1"

# FastQC run with 8 cores
fastqc 02-TRIM/${proj}*${R}* \
    -t 5 \
    --outdir 01-QC/fastqc_trim/fastqc_trim_${R}_${proj}

# FastQC end-timestamp
echo "$(date)       [FastQC Complete]"

# We add --profile-runtime to see the runtime.
multiqc 01-QC/fastqc_trim/fastqc_trim_${R}_${proj} \
    --outdir 01-QC/multiqc_trim/multiqc_trim_${R}_${proj} \
    --profile-runtime
    
# End time and date
echo "$(date)       [End]"