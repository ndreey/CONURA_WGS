#!/bin/bash

module load bioinfo-tools
module load minimap2/2.26-r1175

# Function to display help
usage() {
    echo "Creates PAF files using minimap2"
    echo "Usage: $0 [-T target] [-Q query] [-t num_threads] [-x preset] [-o output]"
    echo "  -T target         Target fasta (required)"
    echo "  -Q query          Query fasta (required)"
    echo "  -t num_threads    Number of threads (required)"
    echo "  -x preset         asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence"
    echo "  -o output         Output folder, it will create if folder does not exist"
    exit 1
}

# Parse command-line options
while getopts ":T:Q:t:x:o:h" opt; do
    case ${opt} in
        T )
            TARGET=${OPTARG}
            ;;
        Q )
            QUERY=${OPTARG}
            ;;
        t )
            NUM_THREADS=${OPTARG}
            ;;
        x )
            PRESET=${OPTARG}
            ;;
        o )
            OUTPUT=${OPTARG}
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
if [ -z "${TARGET}" ] || [ -z "${QUERY}" ] || [ -z "${NUM_THREADS}" ] || [ -z "${PRESET}" ] || [ -z "${OUTPUT}" ]; then
    echo "Error: Options -T, -Q, -t, -x, and -o are required." >&2
    usage
fi

# Create the output directory if it does not exist
if [ ! -d "${OUTPUT}" ]; then
    mkdir -p "${OUTPUT}"
fi


# Extract string to form output
# TARGET parts
POPt=$(echo ${TARGET} | awk -F'/' '{print $2}')
BINNERt=$(echo ${TARGET} | awk -F'/' '{print $3}' | cut -d'_' -f1)
BINt=$(basename ${TARGET} .fa | sed "s/\.//g")

# TQUERY parts
POPq=$(echo ${QUERY} | awk -F'/' '{print $2}')
BINNERq=$(echo ${QUERY} | awk -F'/' '{print $3}' | cut -d'_' -f1)
BINq=$(basename ${QUERY} .fa | sed "s/\.//g")

# Change BINNERt and BINNERq if mag is from bin refinement directory.
if [ ${BINNERt} == "bin" ]; then
    BINNERt="metawrap"
fi

if [ ${BINNERq} == "bin" ]; then
    BINNERq="metawrap"
fi


# Combine the parts
OUT_T="${POPt}-${BINNERt}-${BINt}"
OUT_Q="${POPq}-${BINNERq}-${BINq}"

# Output the result
echo "$result"

minimap2 \
    -x ${PRESET} \
    -t ${NUM_THREADS} \
    ${TARGET} \
    ${QUERY} \
    > ${OUTPUT}/${OUT_Q}_vs_${OUT_T}.paf