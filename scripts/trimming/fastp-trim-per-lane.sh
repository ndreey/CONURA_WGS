#!/bin/bash

#SBATCH --job-name fastp
#SBATCH --array=1-114%20
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 16
#SBATCH -t 05:35:00
#SBATCH --output=SLURM-%j-fastp_task-%a.out
#SBATCH --error=SLURM-%j-fastp_task-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load fastp/0.23.4 

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Path to folder with .lst files
LST_DIR="00-RAW-LANES-LST"

# Get the .lst file for this task
SAMPLE_LST=$(sed -n "${JOBID}p" doc/sample_lst_lane_paths.txt)

# Get the project name and complete sample_id.
# %%_* removes everything starting from the first underscore
PROJ="${SAMPLE_LST%%_*}"
SAMPLE="${SAMPLE_LST%%.*}"

# Path to where .lst paths from.
RAWDIR="/crex/proj/snic2020-6-222/Projects/Tconura/data/WGS/rawdata/${PROJ}"

# Create directory for fastp reports
FASTP_DIR="01-QC/fastp"
if [ ! -d "$FASTP_DIR" ]; then
    mkdir -p "$FASTP_DIR"
fi

# Create directory for trimmed reads
TRIM_DIR="02-TRIM"
if [ ! -d "$TRIM_DIR" ]; then
    mkdir -p "$TRIM_DIR"
fi

# Loop through the raw files and trim.
while read R1 R2; do

    # Input for fastp
    R1_IN="${RAWDIR}/${R1}"
    R2_IN="${RAWDIR}/${R2}"

    # Basename for the lane fastq files
    R1_TRIMD=$(basename $R1 .fastq.gz)
    R2_TRIMD=$(basename $R2 .fastq.gz)

    # Output name for R1 and R2
    R1_OUT="${TRIM_DIR}/${R1_TRIMD}-trim.fastq.gz"
    R2_OUT="${TRIM_DIR}/${R2_TRIMD}-trim.fastq.gz"

    fastp \
        --in1 $R1_IN --in2 $R2_IN \
        --out1 $R1_OUT --out2 $R2_OUT \
        --html "${FASTP_DIR}/${R2_TRIMD}-fastp.html" \
        --json "${FASTP_DIR}/${R2_TRIMD}-fastp.json" \
        --average_qual 20 \
        --length_required 36 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --dedup \
        --thread 16 \
        --verbose 

done < ${LST_DIR}/${SAMPLE_LST}
