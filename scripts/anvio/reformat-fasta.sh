#!/bin/bash

#SBATCH --job-name anvio-alignment
#SBATCH --array=1-2
#SBATCH -A naiss2024-22-580
#SBATCH -p core -n 8
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/anvio/reformat-alignment/SLURM-%j-ref-align-%a.out
#SBATCH --error=slurm-logs/anvio/reformat-alignment/SLURM-%j-ref-align-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)     [Start]"

# Activate the environment
mamba activate anvio-8

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Assign POP based on JOBID
if [ "$JOBID" -eq 1 ]; then
    POP="CHST"
elif [ "$JOBID" -eq 2 ]; then
    POP="COGE"
else
    echo "Invalid JOBID: $JOBID"
    exit 1
fi

# Move to the anvio working directory
cd 07-ANVIO

# Generate folder for reformated fasta
FIX_DIR=00-FIXED-ASSEMBLY
if [ ! -d "${FIX_DIR}" ]; then
    mkdir -p ${FIX_DIR}
fi

# Reformat the contigs for CHST and COGE populations, keeping only 
# contigs equal or above 1kbp.
anvi-script-reformat-fasta \
    --output-file ${FIX_DIR}/${POP}-contigs.fa \
    --prefix ${POP}_ \
    --min-len 1000 \
    --simplify-names \
    ../06-ASSEMBLY/${POP}/contigs.fasta

# Create directory for bam files
ALIGN_DIR=01-ALIGNMENT
if [ ! -d "${ALIGN_DIR}" ]; then
    mkdir -p ${ALIGN_DIR}
fi

# Build bowtie2 index with POP as prefix
bowtie2-build --threads 8 ${FIX_DIR}/${POP}-contigs.fa ${FIX_DIR}/${POP} 

# R1 and R2 reads
READ_DIR="../05-CLEAN-MERGED"
R1=${READ_DIR}/${POP}_R1-clean.fq.gz
R2=${READ_DIR}/${POP}_R2-clean.fq.gz

# Map the clean reads to the reformatted contigs file.
# --no-unal suppresses SAM records for unaligned reads
bowtie2 \
    --threads 8 \
    -x ${FIX_DIR}/${POP} \
    -1 $R1 \
    -2 $R2 \
    --no-unal \
    -S ${ALIGN_DIR}/${POP}.sam

# Convert to bam and exclude all unmapped reads
samtools view \
    -@ 8 \
    -F 4 \
    -b \
    ${ALIGN_DIR}/${POP}.sam > ${ALIGN_DIR}/${POP}-RAW.bam

# Sort and creates index of BAM files
anvi-init-bam \
    -o ${ALIGN_DIR}/${POP}.bam \
    --num-threads 8 \
    ${ALIGN_DIR}/${POP}-RAW.bam

# Remove the .sam and raw bam file
rm ${ALIGN_DIR}/${POP}.sam ${ALIGN_DIR}/${POP}-RAW.bam

# End time and date
echo "$(date)       [End]"