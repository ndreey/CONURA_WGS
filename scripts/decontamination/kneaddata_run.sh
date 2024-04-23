#!/bin/bash

#SBATCH --job-name KneadData
#SBATCH --array=1-9%3
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 16
#SBATCH -t 05:35:00
#SBATCH --output=SLURM-%j-kneaddata.out
#SBATCH --error=SLURM-%j-kneaddata.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load KneadData/0.12.0

# Create directory for decontaminated reads
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/03-DECON"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Path to trimmed reads and reference database of Tconura.
trim_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/02-TRIM"
DB_REF="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/data/Tconura_reference_genome/Tconura_ref"

# SLURM array jobid
jobid=${SLURM_ARRAY_TASK_ID}

# Read the line corresponding to the SLURM array task ID from CHST_samples.txt
IFS="," read -r sample pop hp geo range < <(sed -n "${jobid}p" doc/CHST_samples.txt)

# R1 and R2 input
R1="${trim_dir}/${sample}_R1-trim.fastq.gz"
R2="${trim_dir}/${sample}_R2-trim.fastq.gz"

# Create specific folder for sample.
out="$outdir/${pop}_${sample}"
mkdir -p $out


# Host decontamination without trimming or removal of tandem repeats.
kneaddata \
    --verbose \
    --bypass-trim \
    --bypass-trf \
    --threads 16 \
    --input1 $R1 \
    --input2 $R2 \
    --output $out \
    --reference-db $DB_REF \
    --output-prefix "${pop}_${sample}"

# End time and date
echo "$(date)       [End]"