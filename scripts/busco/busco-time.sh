#!/bin/bash

module load bioinfo-tools
module load BUSCO/5.5.0


# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
BIN=06-DASTOOL/metaWRAP/CHST/score-thresh-0.1_DASTool_bins/maxbin2_bin.2.fa
BUS_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/busco-db/bacteria_odb10
OUT_DIR=09-BUSCO/maxbin2_2
NUM_THREADS=12


if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

busco \
    -i ${BIN} \
    --lineage_dataset ${BUS_DB} \
    --out maxbin2_2 \
    --out_path ${OUT_DIR} \
    --mode genome \
    --cpu ${NUM_THREADS}