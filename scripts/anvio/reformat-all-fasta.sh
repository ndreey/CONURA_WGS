#!/bin/bash

# Start time and date
echo "$(date)     [Start]"

# Activate the environment
#mamba activate anvio-8

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and Variables
POP=all
METAGENOME=all-hybrid
NUM_THREADS=8
FIX_DIR=00-FIXED-ASSEMBLY
ALIGN_DIR=01-ALIGNMENT

# Generate folder for reformated fasta
if [ ! -d "${FIX_DIR}" ]; then
    mkdir -p ${FIX_DIR}
fi


# Reformat the contigs for CHST and COGE populations, keeping only 
# contigs equal or above 1kbp.
anvi-script-reformat-fasta \
    --output-file ${FIX_DIR}/${POP}-contigs.fa \
    --prefix ${POP} \
    --min-len 2500 \
    --simplify-names \
    ../06-ASSEMBLY/${METAGENOME}/contigs.fasta

# Create directory for bam files
if [ ! -d "${ALIGN_DIR}" ]; then
    mkdir -p ${ALIGN_DIR}
fi

# Build bowtie2 index with POP as prefix
bowtie2-build --threads ${NUM_THREADS} ${FIX_DIR}/${POP}-contigs.fa ${FIX_DIR}/${POP} 

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
