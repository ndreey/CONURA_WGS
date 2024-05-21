#!/bin/bash

#SBATCH --job-name metaquast
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 8
#SBATCH -t 01:35:00
#SBATCH --output=SLURM-%j-metaquast.out
#SBATCH --error=SLURM-%j-metaquast.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load quast/5.0.2

quast.py -t 4 -l hybrid,sr,megahit --k-mer-stats --no-icarus -o 01-QC/quast-assembly/CHST 06-ASSEMBLY/CHST/contigs.fasta 06-ASSEMBLY/CHST-sr/contigs.fasta 06-ASSEMBLY/CHST-mega/CHST_mega.contigs.fa


