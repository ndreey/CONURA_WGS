#!/bin/bash

#SBATCH --job-name metaSPAdes
#SBATCH -A naiss2024-5-1
#SBATCH -p node -n 1
#SBATCH -t 28:15:00
#SBATCH -C mem1TB
#SBATCH --output=SLURM-%j-metaSPAdes-CHST.out
#SBATCH --error=SLURM-%j-metaSPAdes-CHST.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load spades/3.15.5

# Set variables
POP="CHST"
SR_DIR="05-CLEAN-MERGED"
LR_DIR="04-CLEAN-FASTQ/hifi-pacbio"

# Create directory for trimmed reads if not existing
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/06-ASSEMBLY/${POP}"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Assembling the metagenome
spades.py \
    --meta \
    --only-assembler \
    -k auto \
    --threads 20 \
    --memory 1000 \
    -1 ${SR_DIR}/${POP}_R1-clean.fq.gz \
    -2 ${SR_DIR}/${POP}_R2-clean.fq.gz \
    --pacbio ${LR_DIR}/pt_042_001_cell1-clean.fastq.gz \
    --pacbio ${LR_DIR}/pt_042_001_cell2-clean.fastq.gz \
    --pacbio ${LR_DIR}/pt_042_001_cell3-clean.fastq.gz \
    -o $outdir

# Restarting from checkpoint
#spades.py --continue -o $outdir