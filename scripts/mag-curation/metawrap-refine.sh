#!/bin/bash

#SBATCH --job-name metaWRAP-refine
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1	
#SBATCH -t 03:00:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-metaWRAP-refine.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-metaWRAP-refine.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load metaWRAP/1.3.2

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
NUM_THREADS=20
BIN_DIR=05-metaWRAP/${POP}
COMP=50
CONT=10

metawrap bin_refinement \
    -o ${BIN_DIR}/bin_refinement_c${COMP}_x${CONT}/ \
    -t ${NUM_THREADS} \
    -m 128 \
    -A ${BIN_DIR}/concoct_bins/ \
    -B ${BIN_DIR}/maxbin2_bins/ \
    -C ${BIN_DIR}/metabat2_bins/ \
    -c ${COMP} \
    -x ${CONT}

# time and date
echo ">>THRESHOLD: ${COMP}% completeness and ${CONT}% contamination"  
echo "$(date)       ${POP}     [END]"