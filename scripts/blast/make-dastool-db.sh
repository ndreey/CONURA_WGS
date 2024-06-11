#!/bin/bash

# Load module
module load bioinfo-tools
module load blast/2.15.0+


# Paths and variables
POP=$1
TOOL=DASTOOL

# Move to 07-ANVIO
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Loop through each bin outputs
for PROG in metaWRAP nf-core; do

    # BLAST directory path
    BLAST_DB_DIR=blastn-db/${POP}/${TOOL}/${PROG}

    # Generate folder for nf-core blastdb
    if [ ! -d "${BLAST_DB_DIR}" ]; then
        mkdir -p ${BLAST_DB_DIR}
    fi

    # Bin path
    BIN_PATH=06-DASTOOL/${PROG}/${POP}/score-thresh-0.1_DASTool_bins/

    # Bin BLASTDB fasta file
    BIN_FA=${BLAST_DB_DIR}/${PROG}-contigs.fa

    # Create a .fasta file with all of $BINNER contigs and edit headers to include the filename
    for BIN in ${BIN_PATH}/*.fa; do
        BIN_ID=$(basename ${BIN} .fa | sed "s/SPAdesHybrid-//g")
        cat ${BIN} | sed "s/>/>${BIN_ID}_/g" >> ${BIN_FA}
    done

    # Create a blastn-db
    makeblastdb \
        -in ${BIN_FA} \
        -input_type fasta \
        -dbtype nucl \
        -parse_seqids \
        -out ${BLAST_DB_DIR}/${PROG}-db
done