#!/bin/bash

#SBATCH --job-name megahit-auto
#SBATCH -A naiss2023-22-412
#SBATCH -p node -n 1
#SBATCH -t 12:35:00
#SBATCH -C mem256GB
#SBATCH --output=SLURM-%j-megahit-P12002_111-auto.out
#SBATCH --error=SLURM-%j-megahit-P12002_111-auto.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load megahit/1.2.9


# Create directory for trimmed reads if not existing
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/03-ASSEMBLY/"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Load in modules
module load bioinfo-tools
module load megahit/1.2.9

megahit \
    --k-list 21,33,55 \
    --num-cpu-threads 20 \
    --memory 0.95 \
    -1 02-TRIM/P12002_111_R1-trim.fastq.gz \
    -2 02-TRIM/P12002_111_R2-trim.fastq.gz \
    -o $outdir/megahit-auto \
    --out-prefix P12002_111_auto