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

            ## Setting up directory
# Create directory for fastqc reports
fastqc_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/fastqc_trim/fastqc_trim_${R}_${proj}"
if [ ! -d "$fastqc_dir" ]; then
    mkdir -p "$fastqc_dir"
fi

# Create a directory for the multiqc reports
multiqc_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/multiqc_trim/multiqc_trim_${R}_${proj}"
if [ ! -d "$multiqc_dir" ]; then
    mkdir -p "$multiqc_dir"
fi


            ## Running FastQC and MultiQC
# FastQC run with 8 cores
fastqc 02-TRIM/${proj}*${R}* \
    -t 5 \
    --outdir ${fastqc_dir}

# FastQC end timestamp
echo "$(date)       [FastQC Complete]"

# MultiQC
multiqc ${fastqc_dir}
    --outdir ${multiqc_dir} \
    --title "${proj}_${R} FastQC report"
    
# End time and date
echo "$(date)       [End]"