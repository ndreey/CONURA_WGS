#!/bin/bash

#SBATCH --job-name P12002_R2_QC_trim-test
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 5
#SBATCH -t 0:35:00
#SBATCH --output=SLURM-%j-P12002_R2_QC_trim-test3.out
#SBATCH --error=SLURM-%j-P12002_R2_QC_trim-test3.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load FastQC/0.11.9 MultiQC/1.12

# Variables for R1/R2 and sequence project
proj="P12002"
R="R2"

# FastQC run with 8 cores
fastqc 02-TRIM/${proj}*${R}* \
    -t 5 \
    --outdir 01-QC/fastqc_trim/fastqc_trim_${R}_${proj}

# FastQC end timestamp
echo "$(date)       [FastQC Complete]"

# We add --profile-runtime to see the runtime.
multiqc 01-QC/fastqc_trim/fastqc_trim_${R}_${proj} \
    --outdir 01-QC/multiqc_trim/multiqc_trim_${R}_${proj} \
    --profile-runtime
    
# End time and date
echo "$(date)       [End]"