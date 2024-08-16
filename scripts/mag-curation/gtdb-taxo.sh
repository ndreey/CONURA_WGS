#!/bin/bash

#SBATCH --job-name GTDB-Tk
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1	
#SBATCH -t 06:00:00
#SBATCH --output=slurm-logs/gtdb-tk/SLURM-%j-make-mash.out
#SBATCH --error=slurm-logs/gtdb-tk/SLURM-%j-make-mash.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load GTDB-Tk/2.4.0

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
WRAP_DIR=05-metaWRAP/${POP}
NUM_THREADS=20


for BINNER in concoct_bins maxbin2_bins metabat2_bins metawrap_50_10_bins; do
    
    # Adjust path for bin-refinment directory
    if [ ${BINNER} == "metawrap_50_10_bins" ]; then
        BIN_DIR=${WRAP_DIR}/bin_refinement_c50_x10/${BINNER}
    else
        BIN_DIR=${WRAP_DIR}/${BINNER}
    fi

    # Create the checkm2 directory
    OUT_DIR=07-GTDB-Tk/${POP}/${BINNER}
    if [ ! -d "${OUT_DIR}" ]; then
        mkdir -p ${OUT_DIR}
    fi

    # Create the checkm2 directory
    MASH_DIR=07-GTDB-Tk/mash-db
    if [ ! -d "${MASH_DIR}" ]; then
        mkdir -p ${MASH_DIR}
    fi

    gtdbtk classify_wf \
        --genome_dir ${BIN_DIR} \
        --extension "fa" \
        --out_dir ${OUT_DIR} \
        --cpus ${NUM_THREADS} \
        --mash_db ${MASH_DIR}/R220_gtdb_ref_sketch.msh

done
