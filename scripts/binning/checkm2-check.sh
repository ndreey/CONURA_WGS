#!/bin/bash

# Takes the population input
POP=$1

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
NF_DIR=04-BINNING
WRAP_DIR=05-metaWRAP
NUM_THREADS=12
CHECKM2_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/checkm2-db/CheckM2_database/uniref100.KO.1.dmnd


for BINNER in CONCOCT MaxBin2 MetaBAT2; do
    # Paths
    BIN_DIR=${NF_DIR}/${POP}/GenomeBinning/${BINNER}/bins
    OUT_DIR=07-CHECKM2/${POP}/nf-core/${BINNER}
    # Create the checkm2 directory
    if [ ! -d "${OUT_DIR}" ]; then
        mkdir -p ${OUT_DIR}
    fi

    checkm2 predict \
        --input ${BIN_DIR} \
        --output-directory ${OUT_DIR} \
        --threads ${NUM_THREADS} \
        --allmodels \
        --extension gz \
        --database_path ${CHECKM2_DB}

done

for BINNER in concoct_bins maxbin2_bins metabat2_bins; do
    BIN_DIR=${WRAP_DIR}/${POP}/${BINNER}
    OUT_DIR=07-CHECKM2/${POP}/metaWRAP/${BINNER}
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
# End time and date
echo "$(date)       ${POP}      [End]"