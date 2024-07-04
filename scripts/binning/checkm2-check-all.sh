#!/bin/bash

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
POP=all
WRAP_DIR=05-metaWRAP
DAS_DIR=06-DASTOOL/all/score-thresh-0.1_DASTool_bins
REFINED=05-metaWRAP/all/bin_refinement/metawrap_50_10_bins
NUM_THREADS=8
CHECKM2_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/checkm2-db/CheckM2_database/uniref100.KO.1.dmnd


for BINNER in concoct_bins maxbin2_bins metabat2_bins; do
    BIN_DIR=${WRAP_DIR}/${POP}/${BINNER}
    OUT_DIR=07-CHECKM2/${POP}/${BINNER}
    # Create the checkm2 directory
    if [ ! -d "${OUT_DIR}" ]; then
        mkdir -p ${OUT_DIR}
    fi    

    checkm2 predict \
        --input ${BIN_DIR} \
        --output-directory ${OUT_DIR} \
        --threads ${NUM_THREADS} \
        --allmodels \
        --extension fa \
        --database_path ${CHECKM2_DB}
done

# Separate run for DASTOOL bins
OUT_DIR=07-CHECKM2/${POP}/DASTOOL
# Create the checkm2 directory
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi   

checkm2 predict \
    --input ${DAS_DIR} \
    --output-directory ${OUT_DIR} \
    --threads ${NUM_THREADS} \
    --allmodels \
    --extension fa \
    --database_path ${CHECKM2_DB}

# Separate run for metaWRAP refinement
OUT_DIR=07-CHECKM2/${POP}/metaWRAP-refine
# Create the checkm2 directory
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi   

checkm2 predict \
    --input ${REFINED} \
    --output-directory ${OUT_DIR} \
    --threads ${NUM_THREADS} \
    --allmodels \
    --extension fa \
    --database_path ${CHECKM2_DB}


# End time and date
echo "$(date)       ${POP}      [End]"