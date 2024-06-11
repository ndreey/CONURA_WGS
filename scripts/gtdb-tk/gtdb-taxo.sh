#!/bin/bash

module load bioinfo-tools
module load GTDB-Tk/2.4.0

# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
BINS=06-DASTOOL/metaWRAP/CHST/score-thresh-0.1_DASTool_bins/
OUT_DIR=010-GTDBtk
NUM_THREADS=12

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi


gtdbtk classify_wf \
    --genome_dir ${BINS} \
    --extension "fa" \
    --out_dir ${OUT_DIR} \
    --cpus ${NUM_THREADS} \
    --skip_ani_screen