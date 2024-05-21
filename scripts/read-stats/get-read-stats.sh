#!/bin/bash

#SBATCH --job-name seqkit-stats
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 6
#SBATCH -t 00:45:00
#SBATCH --output=slurm-logs/decon-stats/SLURM-%j-seqkit-decon.out
#SBATCH --error=slurm-logs/decon-stats/SLURM-%j-seqkit-decon.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load SeqKit/2.4.0

# Start time and date
echo "$(date)  $POP     [Start]"

# Start slurm array jobs for raw reads
sbatch scripts/read-stats/raw-reads-stats.sh

# Start slurm array jobs for trimmed reads
sbatch scripts/read-stats/trim-reads-stats.sh

# Get the read stats for the decontaminated and merged reads
# Create directory
OUT_DIR="results/decon-reads-stats"
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

# Generate stats
seqkit stats \
    --basename \
    --tabular \
    --threads 6 \
    05-CLEAN-MERGED/* > ${OUT_DIR}/decontam-stats.tsv

# Start time and date
echo "$(date)  $POP     [End]"