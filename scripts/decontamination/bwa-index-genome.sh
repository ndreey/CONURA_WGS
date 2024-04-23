#!/bin/bash

#SBATCH --job-name bwa-index
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 8
#SBATCH -t 08:35:00
#SBATCH --output=SLURM-%j-bwa-index.out
#SBATCH --error=SLURM-%j-bwa-index.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load bwa/0.7.17

cd data/Tconura_reference_genome/

# Index genome
bwa index -p Tconura_ref-filtered Tconura_ref-filtered.fasta
