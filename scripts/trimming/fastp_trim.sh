#!/bin/bash

#SBATCH --job-name fastp_test
#SBATCH -A naiss2023-22-412
#SBATCH -p core -n 8
#SBATCH -t 01:35:00
#SBATCH --output=SLURM-%j-fastp_trim-test4.out
#SBATCH --error=SLURM-%j-fastp_trim-test4.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# Load in modules
module load bioinfo-tools
module load fastp/0.23.4 MultiQC/1.12

# Directory with raw reads
raw="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/00-RAW"

# Create directory for trimmed reads if not existing
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/02-TRIM"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Create directory for fastp reports
fastp_dir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/fastp"
if [ ! -d "$fastp_dir" ]; then
    mkdir -p "$fastp_dir"
fi

# Trim reads with fastp
while read sample r1 r2; do

    # Raw
    R1_in="$raw/$r1"
    R2_in="$raw/$r2"

    # Trimmed out
    R1_out="$outdir/${sample}_R1-trim.fastq.gz"
    R2_out="$outdir/${sample}_R2-trim.fastq.gz"

    fastp \
        --in1 $R1_in --in2 $R2_in --out1 $R1_out --out2 $R2_out \
        --html "${fastp_dir}/$sample-fastp.html" \
        --json "${fastp_dir}/$sample-fastp.json" \
        --report_title "$sample fastp report" \
        --average_qual 20 \
        --length_required 100 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --dedup \
        --thread 8 \
        --verbose 

done < test_trim_samples.txt

# MultiQC report of the fastp reports
multiqc 01-QC/fastp/ \
    --outdir 01-QC/multiqc_fastp/ \
    --profile-runtime

# End time and date
echo "$(date)       [End]"