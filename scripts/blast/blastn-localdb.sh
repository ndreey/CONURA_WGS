#!/bin/bash

# Load module
module load bioinfo-tools
module load blast/2.15.0+

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
NUM_THREADS=12
DB_PATH=blastn-db/${POP}/${POP}-db
QUERY=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/doc/Erwinia_Stammerula_Wolbachia-16SrRNA.fa
BLAST_OUT_DIR=08-BLAST-RESULTS
BLAST_OUT=${BLAST_OUT_DIR}/${POP}-blastn.tsv

# Make output dir
if [ ! -d "${BLAST_OUT_DIR}" ]; then
    mkdir -p ${BLAST_OUT_DIR}
fi

# BLASTn search
blastn \
    -query ${QUERY} \
    -db ${DB_PATH} \
    -outfmt "6 qseqid qlen sseqid slen nident pident length mismatch evalue score sstart send sseq" \
    -out ${BLAST_OUT} \
    -num_threads ${NUM_THREADS}