#!/bin/bash

# Function to display help
usage() {
    echo "Creates an anvio contig database"
    echo "Usage: $0 [-p pop] [-m metagenome] [-t num_threads]"
    echo "  -p pop            Population name (required)"
    echo "  -m metagenome     Metagenome assembly name (required)"
    echo "  -t num_threads    Number of threads (required)"
    echo "  -s scg_name       SCG name to use for taxonomy (optional)"
    exit 1
}


while getopts ":p:m:t:s:h" opt; do
    case ${opt} in
        p )
            POP=$OPTARG
            ;;
        m )
            METAGENOME=$OPTARG
            ;;
        t )
            NUM_THREADS=$OPTARG
            ;;
        s )
            SCG=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check if all required options are set
if [ -z "${POP}" ] || [ -z "${METAGENOME}" ] || [ -z "${NUM_THREADS}" ]; then
    echo "Error: All options -p, -m, and -t are required." 1>&2
    usage
fi

# Start time and date
echo "$(date)     [Start]"

# Hard set variables.
DB_DIR="02-CONTIG-DB/${POP}"
CONTIG_DB_OUT="${DB_DIR}/${POP}.db"
PROJ_NAME="Tconura_${POP}"
DESC="../doc/short-project-desc.txt"
COG_DIR=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/ncbi-cogs
SCG_DIR=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/scg/

# Generate folder for contig database
if [ ! -d "${DB_DIR}" ]; then
    mkdir -p ${DB_DIR}
fi

# Generate the anvio contig database
anvi-gen-contigs-database \
    --contigs-fasta ${METAGENOME} \
    --project-name ${PROJ_NAME} \
    --num-threads ${NUM_THREADS} \
    --output-db-path ${CONTIG_DB_OUT} \
    --description ${DESC}

# Run hidden markov models.
anvi-run-hmms \
    --contigs-db ${CONTIG_DB_OUT} \
    --also-scan-trnas \
    --num-threads ${NUM_THREADS}

# Annotate the genes with functions from COG database
anvi-run-ncbi-cogs \
    --contigs-db ${CONTIG_DB_OUT} \
    --cog-version COG20 \
    --cog-data-dir ${COG_DIR} \
    --num-threads ${NUM_THREADS}

# Classify the taxa for the database
anvi-run-scg-taxonomy \
    --contigs-db ${CONTIG_DB_OUT} \
    --scgs-taxonomy-data-dir ${SCG_DIR} \
    --num-threads ${NUM_THREADS}

if [ -n ${SCG} ]; then
    # Estimate the taxonomy
    anvi-estimate-scg-taxonomy \
        --contigs-db ${CONTIG_DB_OUT} \
        --metagenome-mode \
        --num-threads ${NUM_THREADS} \
        --scg-name-for-metagenome-mode ${SCG}
else
    # Estimate the taxonomy
    anvi-estimate-scg-taxonomy \
        --contigs-db ${CONTIG_DB_OUT} \
        --metagenome-mode \
        --num-threads ${NUM_THREADS}
fi

# End time and date
echo "$(date)  ${POP}    [End]"