#!/bin/bash

# Load module
module load bioinfo-tools
module load blast/2.15.0+


# Paths and variables
POP=$1
TOOL=nf-core

# Move to 07-ANVIO
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Loop through each bin outputs
for BINNER in CONCOCT MaxBin2 MetaBAT2; do

    # BLAST directory path
    BLAST_DB_DIR=blastn-db/${POP}/${TOOL}/${BINNER}

    # Generate folder for nf-core blastdb
    if [ ! -d "${BLAST_DB_DIR}" ]; then
        mkdir -p ${BLAST_DB_DIR}
    fi

    # Bin path
    BIN_PATH=04-BINNING/${POP}/GenomeBinning/${BINNER}/bins

    # Bin BLASTDB fasta file
    BIN_FA=${BLAST_DB_DIR}/${BINNER}-contigs.fa

    # Create a .fasta file with all of $BINNER contigs and edit headers to include the filename
    for BIN in ${BIN_PATH}/*.fa.gz; do
        BIN_ID=$(basename ${BIN} .fa.gz)
        zcat ${BIN} | sed "s/>/>${BIN_ID}_/g" >> ${BIN_FA}
    done
    
    # Create a blastn-db
    makeblastdb \
        -in ${BIN_FA} \
        -input_type fasta \
        -dbtype nucl \
        -parse_seqids \
        -out ${BLAST_DB_DIR}/${BINNER}-db
done