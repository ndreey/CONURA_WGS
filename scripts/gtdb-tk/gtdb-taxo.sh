#!/bin/bash

module load bioinfo-tools
module load GTDB-Tk/2.4.0

# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
NUM_THREADS=12
BINS1=05-metaWRAP/all/bin_refinement/metawrap_50_10_bins/
BINS2=06-DASTOOL/all/score-thresh-0.1_DASTool_bins/


for BINS in ${BINS1} ${BINS2}; do

    if [ "${BINS}" = "${BINS1}" ]; then
        OUT_DIR=08-GTDB-tk/metaWRAP-refine
    else
        OUT_DIR=08-GTDB-tk/DASTOOL
    fi

    mkdir -p ${OUT_DIR}

    gtdbtk classify_wf \
        --genome_dir ${BINS} \
        --extension "fa" \
        --out_dir ${OUT_DIR} \
        --cpus ${NUM_THREADS} \
        --skip_ani_screen

done
