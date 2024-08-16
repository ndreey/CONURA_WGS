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
if [ -z "${POP}" ] || [ -z "${NUM_THREADS}" ]; then
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


# Create directory for bam files
if [ ! -d "${ALIGN_DIR}" ]; then
    mkdir -p ${ALIGN_DIR}
fi

# Check if bowtie2 index files exist
INDEX_PREFIX="${FIX_DIR}/${POP}"
if [ ! -e "${INDEX_PREFIX}.1.bt2" ] || [ ! -e "${INDEX_PREFIX}.2.bt2" ] || [ ! -e "${INDEX_PREFIX}.3.bt2" ] || [ ! -e "${INDEX_PREFIX}.4.bt2" ] || [ ! -e "${INDEX_PREFIX}.rev.1.bt2" ] || [ ! -e "${INDEX_PREFIX}.rev.2.bt2" ]; then
    echo "Bowtie2 index files do not exist. Running bowtie2-build..."
    bowtie2-build \
        --threads ${NUM_THREADS} \
        ${FIX_DIR}/${POP}-contigs.fa ${INDEX_PREFIX}
else
    echo "Bowtie2 index files already exist. Skipping bowtie2-build."
fi

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

done < ../doc/all-clean-reads.txt

# End time and date
echo "$(date)       [End]"
