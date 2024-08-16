#!/bin/bash

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
WRAP_DIR=05-metaWRAP/${POP}
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
CONT2BIN=contigs2bin/${POP}

# Generate folder for contigs2bin .tsv files
if [ ! -d "${CONT2BIN}" ]; then
    mkdir -p ${CONT2BIN}
fi

# Loop through each bin outputs
for BINNER in concoct maxbin2 metabat2 metawrap; do
    OUT_TSV=${CONT2BIN}/${POP}-${BINNER}-contigs2bin.tsv
    
     # Bin path
    if [ ${BINNER} == "metawrap" ]; then
        BIN_PATH=${WRAP_DIR}/bin_refinement_c50_x10/${BINNER}_50_10_bins
    else
        BIN_PATH=${WRAP_DIR}/${BINNER}_bins
    fi

    for BIN in ${BIN_PATH}/*.fa; do   
        # Get the basename of the bin without extension
        BIN_ID=$(basename ${BIN} .fa | tr "." "_")
        # Get the contig-ids and concatenate to .tsv file with bin-id
        cat ${BIN} | grep "^>" | sed "s/>//g" | sed "s/$/\t${BIN_ID}/g" >> ${OUT_TSV}
    done
done

# End time and date
echo "$(date)       ${POP}      [End]"