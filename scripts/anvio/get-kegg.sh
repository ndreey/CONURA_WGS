#!/bin/bash

#SBATCH --job-name anvio-contigdb
#SBATCH -A naiss2024-22-580
#SBATCH -p core -n 6
#SBATCH -t 06:30:00
#SBATCH --output=slurm-logs/anvio/SLURM-%j-setup-kegg.out
#SBATCH --error=slurm-logs/anvio/SLURM-%j-setup-kegg.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)     [Start]"

# Activate the environment
mamba activate anvio-8

anvi-setup-kegg-data \
    --mode all \
    --kegg-data-dir ../databases/kegg \
    -T 6 \
    --reset

# End time and date
echo "$(date)     [End]"