#!/bin/bash

# Start time and date
echo "$(date)     [Start]"

# Activate the environment
#mamba activate anvio-8

# Paths and variables
POP=all
FIX_DIR=00-FIXED-ASSEMBLY
DB_DIR=02-CONTIG-DB
CDB=${DB_DIR}/${POP}.db
NUM_THREADS=8
COG_DIR=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/ncbi-cogs
SCG_DIR=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/scg/

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Generate folder for contig database
if [ ! -d "${DB_DIR}" ]; then
    mkdir -p ${DB_DIR}
fi

# Generate the anvio contig database
anvi-gen-contigs-database \
    --contigs-fasta ${FIX_DIR}/${POP}-contigs.fa \
    --project-name Tconura_${POP} \
    --num-threads ${NUM_THREADS} \
    --output-db-path ${CDB} \
    --description ../doc/short-project-desc.txt

# Run hidden makarov models.
anvi-run-hmms \
    --contigs-db ${CDB} \
    --also-scan-trnas \
    --num-threads ${NUM_THREADS}

# Annotate the genes with functions from COG database
anvi-run-ncbi-cogs \
    --contigs-db ${CDB} \
    --cog-version COG20 \
    --cog-data-dir ${COG_DIR} \
    --num-threads ${NUM_THREADS}

# Classify the taxa for the database
anvi-run-scg-taxonomy \
    --contigs-db ${CDB} \
    --scgs-taxonomy-data-dir ${SCG_DIR} \
    --num-threads ${NUM_THREADS}

# Estimate the taxonomy
anvi-estimate-scg-taxonomy \
    --contigs-db ${CDB} \
    --metagenome-mode \
    --num-threads ${NUM_THREADS}

# End time and date
echo "$(date)  ${POP}    [End]"