#!/bin/bash

#SBATCH --job-name megahit-COGE
#SBATCH -A naiss2024-5-1
#SBATCH -p node -n 1
#SBATCH -t 00:35:00
#SBATCH -C mem256GB
#SBATCH --output=SLURM-%j-megahit-COGE.out
#SBATCH --error=SLURM-%j-megahit-COGE.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load megahit/1.2.9

# Set variables
POP="COGE"
SR_DIR="05-CLEAN-MERGED"
LR_DIR="04-CLEAN-FASTQ/hifi-pacbio"

# Load in modules
module load bioinfo-tools
module load megahit/1.2.9

megahit \
    --k-list 21,33,55 \
    --num-cpu-threads 20 \
    --memory 0.95 \
    -1 ${SR_DIR}/${POP}_R1-clean.fq.gz \
    -2 ${SR_DIR}/${POP}_R2-clean.fq.gz \
    -o 06-ASSEMBLY/COGE-mega \
    --out-prefix COGE_mega