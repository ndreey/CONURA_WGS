#!/bin/bash

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
POP=all
WRAP_DIR=05-metaWRAP
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
CONT2BIN=contigs2bin/${POP}
NUM_THREADS=8

# Generate folder for contigs2bin .tsv files
if [ ! -d "${CONT2BIN}" ]; then
    mkdir -p ${CONT2BIN}
fi

# Loop through each bin outputs
for BINNER in concoct maxbin2 metabat2; do
    OUT_TSV=${CONT2BIN}/${POP}-${BINNER}-contigs2bin.tsv
    # Bin path
    BIN_PATH=${WRAP_DIR}/${POP}/${BINNER}_bins
    for BIN in ${BIN_PATH}/*.fa; do   
        # Get the basename of the bin without extension
        BIN_ID=$(basename ${BIN} .fa)
        # Get the contig-ids and concatenate to .tsv file with bin-id
        cat ${BIN} | grep "^>" | sed "s/>//g" | sed "s/$/\t${BIN_ID}/g" >> ${OUT_TSV}
    done
done

# End time and date
echo "$(date)       ${POP}      [End]"