#!/bin/bash

#SBATCH --job-name minimap2-index
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 4 
#SBATCH -t 01:35:00
#SBATCH --output=slurm-logs/decontamination/SLURM-%j-minimap2-index.out
#SBATCH --error=slurm-logs/decontamination/SLURM-%j-minimap2-index.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load minimap2/2.26-r1175

# Move to reference directory
cd data/Tconura_reference_genome/

# Index genome
minimap2 -t 4 -x map-hifi -d Tconura_ref-filtered.mmi Tconura_ref-filtered.fasta

# End time and date
echo "$(date)       [End]"
