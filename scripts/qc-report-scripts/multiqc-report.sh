#!/bin/bash

#SBATCH --job-name multiqc
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 2
#SBATCH -t 00:35:00
#SBATCH --output=SLURM-%j-multiqc.out
#SBATCH --error=SLURM-%j-multiqc.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load MultiQC/1.12

# Create a directory for the fastp multiqc reports
MQC_DIR="01-QC/multiqc_fastp"
if [ ! -d "$MQC_DIR" ]; then
    mkdir -p "$MQC_DIR"
fi

# MultiQC report
multiqc 01-QC/fastp \
    --outdir ${MQC_DIR} \
    --title "fastp reports from ALL SAMPLES"

# End time and date
echo "$(date)       [End]"