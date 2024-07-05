#!/bin/bash

# Function to display help
usage() {
    echo "Creates an anvio profile database"
    echo "Usage: $0 [-p pop] [-t num_threads] [-c contig_db] [-l min_len]"
    echo "  -p pop            Population name (required)"
    echo "  -t num_threads    Number of threads (required)"
    echo "  -c contig_db      Contig database (required)"
    echo "  -l min_len        Minimum length (optional, default: 2500)"
    exit 1
}

# Default set values
MIN_LEN=2500

# Parse command-line options
while getopts ":p:t:l:c:h" opt; do
    case ${opt} in
        p )
            POP=${OPTARG}
            ;;
        t )
            NUM_THREADS=${OPTARG}
            ;;
        l )
            MIN_LEN=${OPTARG}
            ;;
        c )
            CONTIG_DB=${OPTARG}
            ;;
        h )
            usage
            ;;
        \? )
            echo "Invalid option: -${OPTARG}" >&2
            usage
            ;;
        : )
            echo "Option -${OPTARG} requires an argument." >&2
            usage
            ;;
    esac
done

# Check if all required options are set
if [ -z "${POP}" ] || [ -z "${NUM_THREADS}" ] || [ -z "${CONTIG_DB}" ]; then
    echo "Error: Options -p, -t, and -c are required." >&2
    usage
fi

# Start time and date
echo "$(date)     [Start]"

# Hard set variables.
PROF_DB=03-PROFILES/${POP}
ALIGN_DIR=01-ALIGNMENT/${POP}
DESC="../doc/short-project-desc.txt"

# Create directory for profile files
if [ ! -d "${PROF_DB}" ]; then
    mkdir -p ${PROF_DB}
fi

# Specifies which populations to use
if [[ "${POP}" == "all" ]]; then
    # Generate profiles for each population
    while read SAMPLE; do

        # Dynamic variable
        BAM_IN=${ALIGN_DIR}/${SAMPLE}.bam

        anvi-profile \
            --input-file ${BAM_IN} \
            --contigs-db ${CONTIG_DB} \
            --sample-name ${SAMPLE} \
            --description ${DESC} \
            --min-contig-length ${MIN_LEN} \
            --output-dir ${PROF_DB}/${SAMPLE} \
            --num-threads ${NUM_THREADS}

    done < ../doc/populations.txt
else
    # Specific bam file for user input population
    BAM_IN=${ALIGN_DIR}/${POP}.bam

    anvi-profile \
        --input-file ${BAM_IN} \
        --contigs-db ${CONTIG_DB} \
        --sample-name ${POP} \
        --description ${DESC} \
        --min-contig-length ${MIN_LEN} \
        --output-dir ${PROF_DB}/${POP} \
        --num-threads ${NUM_THREADS} \
        --cluster-contigs
fi

# End time and date
echo "$(date)       ${POP}      [End]"