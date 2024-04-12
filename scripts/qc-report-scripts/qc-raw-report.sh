#!/bin/bash

#SBATCH --job-name FastQC
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 10
#SBATCH -t 02:35:00
#SBATCH --output=SLURM-%j-raw_QC_reports.out
#SBATCH --error=SLURM-%j-raw_QC_reports.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load FastQC/0.11.9 MultiQC/1.12

# Variables for R1/R2 and sequence project
proj=$1
R=$2

# What sequence project and read is being analyzed
echo "Generating reports for: $proj $R"

# Create directory for fastqc reports
fastqc_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/fastqc_raw/fastqc_raw_${R}_${proj}"
if [ ! -d "$fastqc_dir" ]; then
    mkdir -p "$fastqc_dir"
fi

# Create a directory for the multiqc reports
multiqc_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/multiqc_raw/multiqc_raw_${R}_${proj}"
if [ ! -d "$multiqc_dir" ]; then
    mkdir -p "$multiqc_dir"
fi


# FastQC run with 8 cores
fastqc 00-RAW/${proj}*${R}* \
    -t 10 \
    --outdir ${fastqc_dir}

# FastQC end timestamp
echo "$(date)       [FastQC Complete]"

# MultiQC report
multiqc ${fastqc_dir} \
    --outdir ${multiqc_dir} \
    --title "RAW ${proj}_${R} FastQC report"
    
# End time and date
echo "$(date)       [End]"
echo "Finished generating reports for: $proj $R"