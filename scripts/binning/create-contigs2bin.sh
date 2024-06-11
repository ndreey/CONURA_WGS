#!/bin/bash

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO


# Paths and variables
POP=$1
NF_DIR=04-BINNING
WRAP_DIR=05-metaWRAP
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
NUM_THREADS=12



                ###### nf-core/mag BINS ##########

# Generate folder for nf-core contigs2bin .tsv files
NF_C2B=contigs2bin/nf-core
if [ ! -d "${NF_C2B}" ]; then
    mkdir -p ${NF_C2B}
fi

# Loop through each bin outputs
for BINNER in CONCOCT MaxBin2 MetaBAT2; do
    # Output name for .tsv
    NF_TSV=${NF_C2B}/${POP}-${BINNER}_contigs2bin.tsv
    # Bin path
    BIN_PATH=${NF_DIR}/${POP}/GenomeBinning/${BINNER}/bins
    for BIN in ${BIN_PATH}/*.fa.gz; do    
        # Get the basename of the bin without extension
        BIN_ID=$(basename ${BIN} .fa.gz)
        # Get the contig-ids and concatenate to .tsv file with bin-id
        zcat ${BIN} | grep "^>" | sed "s/>//g" | sed "s/$/\t${BIN_ID}/g" >> ${NF_TSV}        
    done
done



                ###### metaWRAP BINS ##########

# Generate folder for nf-core contigs2bin .tsv files
WRAP_C2B=contigs2bin/metaWRAP
if [ ! -d "${WRAP_C2B}" ]; then
    mkdir -p ${WRAP_C2B}
fi

# Loop through each bin outputs
for BINNER in concoct_bins maxbin2_bins metabat2_bins; do
    WRAP_TSV=${WRAP_C2B}/${POP}-${BINNER}_contigs2bin.tsv
    # Bin path
    BIN_PATH=${WRAP_DIR}/${POP}/${BINNER}
    for BIN in ${BIN_PATH}/*.fa; do
        # Check if the file name contains "unbinned"
        if [[ "${BIN}" == *unbinned* ]]; then
            echo "Skipping ${BIN} as it contains 'unbinned'"
            continue
        fi        
        # Get the basename of the bin without extension
        BIN_ID=$(basename ${BIN} .fa)
        # Get the contig-ids and concatenate to .tsv file with bin-id
        cat ${BIN} | grep "^>" | sed "s/>//g" | sed "s/$/\t${BIN_ID}/g" >> ${WRAP_TSV}
    done
done



# End time and date
echo "$(date)       ${POP}      [End]"