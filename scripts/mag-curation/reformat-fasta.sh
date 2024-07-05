#!/bin/bash

# Start time and date
echo "$(date)     [Start]"
echo "Reformat metagenome assembly to standardise contig headers and discard small contigs"

# Function to display help
usage() {
    echo "Usage: $0 [-p pop] [-m metagenome] [-t num_threads]"
    echo "  -p pop            Population name (required)"
    echo "  -m metagenome     Metagenome assembly name (required)"
    echo "  -t num_threads    Number of threads (required)"
    echo "  -l min_length     Minimum contig length (optional, default: 2500)"
    exit 1
}

# Default variable
MIN_LENGTH=2500

while getopts ":p:m:t:l:h" opt; do
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
        l )
            MIN_LENGTH=$OPTARG
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

# Hard set variables.
FIX_DIR="00-FIXED-ASSEMBLY/${POP}"
ALIGN_DIR="01-ALIGNMENT/${POP}"

# Generate folder for reformated fasta
if [ ! -d "${FIX_DIR}" ]; then
    mkdir -p ${FIX_DIR}
fi

# Reformat the contigs of the assembly, keeping only 
# contigs equal or above 2.5 Kbp.
anvi-script-reformat-fasta \
    --output-file ${FIX_DIR}/${POP}-contigs.fa \
    --prefix ${POP} \
    --min-len ${MIN_LENGTH} \
    --simplify-names \
    ../06-ASSEMBLY/${METAGENOME}/contigs.fasta

# Create directory for bam files
if [ ! -d "${ALIGN_DIR}" ]; then
    mkdir -p ${ALIGN_DIR}
fi

# Generate a temp file holding the clean reads for the chosen pop.
TMP_FILE=$(mktemp)

# Greps the reads to be used for alignment
if [[ "${POP}" == "all" ]]; then
    cat ../doc/all-clean-reads.txt > "${TMP_FILE}"
else 
    cat ../doc/all-clean-reads.txt | grep "${POP}" >  "${TMP_FILE}"
fi

# Build bowtie2 index with POP as prefix
bowtie2-build \
    --threads ${NUM_THREADS} \
    ${FIX_DIR}/${POP}-contigs.fa ${FIX_DIR}/${POP} 

while read R1 R2; do
    # Get the population sample name.
    SAMPLE=$(echo ${R1} | cut -d "/" -f2 | cut -d "_" -f1)
    echo "Mapping ${SAMPLE}"

    # Map the clean reads to the reformatted contigs file.
    # --no-unal suppresses SAM records for unaligned reads
    bowtie2 \
        --threads ${NUM_THREADS} \
        -x ${FIX_DIR}/${POP} \
        -1 ../${R1} \
        -2 ../${R2} \
        --no-unal \
        -S ${ALIGN_DIR}/${SAMPLE}.sam

    # Convert to bam and exclude all unmapped reads
    samtools view \
        -@ ${NUM_THREADS} \
        -F 4 \
        -b \
        ${ALIGN_DIR}/${SAMPLE}.sam > ${ALIGN_DIR}/${SAMPLE}-RAW.bam

    # Sort and creates index of BAM files
    anvi-init-bam \
        -o ${ALIGN_DIR}/${SAMPLE}.bam \
        --num-threads ${NUM_THREADS} \
        ${ALIGN_DIR}/${SAMPLE}-RAW.bam

    # Remove the .sam and raw bam file
    rm ${ALIGN_DIR}/${SAMPLE}.sam ${ALIGN_DIR}/${SAMPLE}-RAW.bam

done < ${TMP_FILE}

# Remove the temporary file
rm ${TMP_FILE}

# End time and date
echo "$(date)       [End]"
