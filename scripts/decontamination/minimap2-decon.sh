#!/bin/bash

#SBATCH --job-name minimap2-decontamination
#SBATCH -A naiss2023-22-412
#SBATCH --array=1-3
#SBATCH -p node -n 1
#SBATCH -t 06:00:00
#SBATCH --output=slurm-logs/decontamination/SLURM-%j-minimap2-decon-hifi.out
#SBATCH --error=slurm-logs/decontamination/SLURM-%j-minimap2-decon-hifi.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load samtools/1.19
module load minimap2/2.26-r1175
module load BEDTools/2.31.1

# Path to the hifi pacbio raw data
RAW="/crex/proj/snic2020-6-222/Projects/Tconura/data/reference/hifiasm_Assemb2020_pt_042/pt_042/ccsreads/pt_042_001"

# Path to trimmed reads and reference database of Tconura.
REF="data/Tconura_reference_genome/Tconura_ref-filtered"

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Read the fastq files this task will work with.
READ=$(sed -n "${JOBID}p" doc/pt_042_hifi-pacbio.txt)

# Get the sample id
SAMPLE="${READ%%.*}"

# Info text
echo "$(date)   Processing: $SAMPLE"

# Create directory for sample specific bam files
BAM_DIR="03-HOST-ALIGNMENT-BAM/hifi-pacbio/${SAMPLE}"
if [ ! -d "${BAM_DIR}" ]; then
    mkdir -p "${BAM_DIR}"
fi

# Run minimap2 alignment.
minimap2 -a -x map-hifi -t 16 Tconura_ref-filtered.mmi $READ | \
    samtools sort - -@ 16 -o ${BAM_DIR}/${SAMPLE}.bam
echo "$(date)   minimap2 alignment complete"

# Generate a bam with all unmapped reads.
samtools view -@ 16 -b -f 4 ${BAM_DIR}/${SAMPLE}.bam \
    > ${BAM_DIR}/${SAMPLE}-unmapped-reads.bam
echo "$(date)   samtools filtering complete"

# Create directory for decontaminated clean fastq files
CLEAN_DIR="04-CLEAN-FASTQ/hifi-pacbio"
if [ ! -d "$CLEAN_DIR" ]; then
    mkdir -p "$CLEAN_DIR"
fi

# Convert the bam file to fastq
bedtools bamtofastq \
    -i  ${BAM_DIR}/${SAMPLE}-unmapped-reads.bam\
    -fq ${CLEAN_DIR}/${SAMPLE}-clean.fastq 
echo "$(date)   bamtofastq complete"

# Compress the fastq files
gzip ${CLEAN_DIR}/${SAMPLE}-clean.fastq
echo "$(date)   clean-fastq compressed"

# End time and date
echo "$(date)       [End]"