#!/bin/bash

# Load module
module load bioinfo-tools
module load blast/2.15.0+


# Paths and variables
POP=$1
BLAST_DB_DIR=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO/blastn-db/${POP}
QUERY=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/doc/Erwinia_Stammerula_Wolbachia-16SrRNA.fa
# Iterate over tools and binning methods using ls command
for TOOL in $(ls ${BLAST_DB_DIR}/); do
    for BINNER in $(ls ${BLAST_DB_DIR}/${TOOL}/); do

        DB_PATH=${BLAST_DB_DIR}/${TOOL}/${BINNER}/${BINNER}-db
        BLAST_OUT_DIR=08-BLAST-RESULTS
        BLAST_OUT=${BLAST_OUT_DIR}/blastn-${POP}-${TOOL}-${BINNER}.tsv

        if [ ! -d "${BLAST_OUT_DIR}" ]; then
            mkdir -p ${BLAST_OUT_DIR}
        fi

        # BLASTn search
        blastn \
            -query ${QUERY} \
            -db ${DB_PATH} \
            -outfmt "6 qseqid qlen sseqid slen nident pident length mismatch evalue bitscore score sstart ssend sseq" \
            -out ${BLAST_OUT} \
            -num_threads 4

    done
done
