#!/bin/bash

#SBATCH --job-name metaSPAdes-auto
#SBATCH -A naiss2023-22-412
#SBATCH -p node -n 1
#SBATCH -t 28:35:00
#SBATCH -C mem1TB
#SBATCH --output=SLURM-%j-spades-P12002_111-auto.out
#SBATCH --error=SLURM-%j-spades-P12002_111-auto.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load spades/3.15.5


# Create directory for trimmed reads if not existing
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/03-ASSEMBLY/spades-auto"
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
    -1 02-TRIM/P12002_111_R1-trim.fastq.gz \
    -2 02-TRIM/P12002_111_R2-trim.fastq.gz \
    -o $outdir

# Restarting from checkpoint
#spades.py --continue -o $outdir