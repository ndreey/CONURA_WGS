#!/bin/bash

#SBATCH --job-name bwa-decontamination
#SBATCH -A naiss2023-22-412
#SBATCH --array=1-304%40
#SBATCH -p core -n 8
#SBATCH -t 04:30:00
#SBATCH --output=slurm-logs/decontamination/SLURM-%j-bwa-decon-%a.out
#SBATCH --error=slurm-logs/decontamination/SLURM-%j-bwa-decon-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.19
module load BEDTools/2.31.1

# Path to trimmed reads and reference database of Tconura.
TRIM_DIR="02-TRIM"
REF="data/Tconura_reference_genome/Tconura_ref-filtered"

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Read the fastq files this task will work with.
FASTQ_FILES=$(sed -n "${JOBID}p" doc/trimmed_fastq.txt)

# Get the sample id, R1, R2 and lane id.
SAMPLE="$(echo "$FASTQ_FILES" | cut -d "_" -f 1,2)"
R1="$(echo "$FASTQ_FILES" | cut -f 1)"
R2="$(echo "$FASTQ_FILES" | cut -f 2)"
LANE_ID="${R1%%_R1*}"

# Get the population the sample belongs to
while IFS="," read SAMPLE_COMP POP HP REGION RANGE; do
    if [ "$SAMPLE" == "$SAMPLE_COMP" ]; then
        OUT_POP="$POP"
    fi
done < "doc/metadata-no-hybrids.csv"

# Info text
echo "$(date)   Processing: $LANE_ID"
echo "$(date)   Population: $OUT_POP"

# Create directory for sample specific bam files
BAM_DIR="03-HOST-ALIGNMENT-BAM/${OUT_POP}/${LANE_ID}"
if [ ! -d "${BAM_DIR}" ]; then
    mkdir -p "${BAM_DIR}"
fi

# Run bwa-mem algorithm and pipe it to .bam file
bwa mem -t 8 $REF ${TRIM_DIR}/${R1} ${TRIM_DIR}/${R2} | \
    samtools view - -b -@ 8 > ${BAM_DIR}/${LANE_ID}.bam
echo "$(date)   BWA Alignment finished"

# Generate a bam with all unmapped reads.
samtools view -@ 8 -b -f 4 ${BAM_DIR}/${LANE_ID}.bam \
    > ${BAM_DIR}/${LANE_ID}-unmapped-reads.bam

# Subset the pairs that did not align to reference
samtools view -@ 8 -b -f 12 ${BAM_DIR}/${LANE_ID}.bam \
    > ${BAM_DIR}/${LANE_ID}-unmapped-pairs.bam

# Sort by name so BEDtools can convert to fastq.
samtools sort -n -@ 8 \
    -o ${BAM_DIR}/${LANE_ID}-sorted-unmapped-pairs.bam \
    ${BAM_DIR}/${LANE_ID}-unmapped-pairs.bam
echo "$(date)   SAMtools complete"

# Create directory for decontaminated reads in population
CLEAN_DIR="04-CLEAN-FASTQ/${OUT_POP}"
if [ ! -d "$CLEAN_DIR" ]; then
    mkdir -p "$CLEAN_DIR"
fi

# Generate R1 and R2 clean fastq names.
R1_OUT=${R1%%_001*}
R2_OUT=${R2%%_001*}

# Convert the bam file to R1 and R2 fastq files
bedtools bamtofastq \
    -i ${BAM_DIR}/${LANE_ID}-sorted-unmapped-pairs.bam \
    -fq ${CLEAN_DIR}/${R1_OUT}-clean.fastq \
    -fq2 ${CLEAN_DIR}/${R2_OUT}-clean.fastq
echo "$(date)   bamtofastq complete"

# Compress the fastq files
gzip ${CLEAN_DIR}/${R1_OUT}-clean.fastq
gzip ${CLEAN_DIR}/${R2_OUT}-clean.fastq
echo "$(date)   clean-fastq compressed"

# End time and date
echo "$(date)       [End]"