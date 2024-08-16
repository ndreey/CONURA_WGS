#!/bin/bash

# Load module
module load bioinfo-tools
module load blast/2.15.0+

# Move to 08-ANVIO
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
WRAP_DIR=05-metaWRAP/${POP}
BLAST_DB_DIR=blastn-db/${POP}
BIN_FA=${BLAST_DB_DIR}/merged-bins.fa

# Generate folder for nf-core blastdb
if [ ! -d "${BLAST_DB_DIR}" ]; then
    mkdir -p ${BLAST_DB_DIR}
fi

for BINNER in concoct_bins maxbin2_bins metabat2_bins metawrap_50_10_bins; do

    # Adjust path for bin-refinment directory
    if [ ${BINNER} == "metawrap_50_10_bins" ]; then
        BIN_DIR=${WRAP_DIR}/bin_refinement_c50_x10/${BINNER}
    else
        BIN_DIR=${WRAP_DIR}/${BINNER}
    fi

    # Create a .fasta file with all of $BINNER contigs and edit headers to include the filename
    for BIN in ${BIN_DIR}/*.fa; do
        BIN_ID=$(basename ${BIN} .fa)
        B_NAME="${BINNER%%_*}"
        cat ${BIN} | sed "s/>/>${B_NAME}_${BIN_ID}_/g" >> ${BIN_FA}
    done

done

# Create a blastn-db
makeblastdb \
    -in ${BIN_FA} \
    -input_type fasta \
    -dbtype nucl \
    -parse_seqids \
    -out ${BLAST_DB_DIR}/${POP}-db
