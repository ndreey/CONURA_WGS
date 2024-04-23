#!/bin/bash

#SBATCH --job-name bwa-decontamination
#SBATCH --array=1-9
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 8
#SBATCH -t 05:35:00
#SBATCH --output=SLURM-%j-bwa-decon-%A_%a.out
#SBATCH --error=SLURM-%j-bwa-decon-%A_%a.err
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
trim_dir="02-TRIM"
REF="data/Tconura_reference_genome/Tconura_ref-filtered"

# SLURM array jobid
jobid=${SLURM_ARRAY_TASK_ID}

# Read the line corresponding to the SLURM array task ID from CHST_samples.txt
IFS="," read -r sample pop hp geo range < <(sed -n "${jobid}p" doc/CHST_samples.txt)

# Create directory for decontaminated reads
outdir="03-HOST-DECON/${pop}_${sample}"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# R1 and R2 input
R1="${trim_dir}/${sample}_R1-trim.fastq.gz"
R2="${trim_dir}/${sample}_R2-trim.fastq.gz"

# Run bwa-mem algorithm
bwa mem -t 8 -M $REF $R1 $R2 | \
    samtools sort -n --threads 8 -o $outdir/${sample}.bam

# Subset the pairs that did not align to reference
samtools view -@ 8 -b -F 2 $outdir/${sample}.bam \
    > $outdir/${sample}-unmapped-paired.bam

# Convert the bam file to R1 and R2 fastq files
bedtools bamtofastq \
    -i $outdir/${sample}-unmapped-paired.bam \
    -fq $outdir/${sample}_R1-clean.fastq \
    -fq2 $outdir/${sample}_R2-clean.fastq

# End time and date
echo "$(date)       [End]"