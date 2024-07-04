#!/bin/bash

module load bioinfo-tools
module load BUSCO/5.5.0


# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
BIN=05-metaWRAP/all/bin_refinement/metawrap_50_10_bins/bin.7.fa
BUS_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/busco-db/bacteria_odb10
OUT_DIR=010-BUSCO/WOLBACHIA
NUM_THREADS=4

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

busco \
    -i ${BIN} \
    --lineage_dataset ${BUS_DB} \
    --out WOLBACHIA \
    --out_path ${OUT_DIR} \
    --mode genome \
    --cpu ${NUM_THREADS}