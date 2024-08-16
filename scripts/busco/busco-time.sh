#!/bin/bash

module load bioinfo-tools
module load BUSCO/5.5.0


# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
BIN=$1
TAXA=$2
BUS_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/busco-db/bacteria_odb10
NUM_THREADS=6
OUT_DIR=09-BUSCO/${TAXA}

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

busco \
    -i ${BIN} \
    --lineage_dataset ${BUS_DB} \
    --out ${TAXA} \
    --out_path ${OUT_DIR} \
    --mode genome \
    --cpu ${NUM_THREADS}