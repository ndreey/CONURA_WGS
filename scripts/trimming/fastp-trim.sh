#!/bin/bash

#SBATCH --job-name fastp
#SBATCH -A naiss2023-22-412
#SBATCH -p node -n 1
#SBATCH -t 10:35:00
#SBATCH --output=SLURM-%j-fastp-trim.out
#SBATCH --error=SLURM-%j-fastp-trim.err
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

    # Store reads where either R1 or R2 does not meet thresholds
    # Discard reads below average quality of Q20
    # Discard reads that does not meet minimum length 36
    # De-duplicates and remove adapters, and polyG's
    fastp \
        --in1 $R1_in --in2 $R2_in --out1 $R1_out --out2 $R2_out \
        --html "${fastp_dir}/$sample-fastp.html" \
        --json "${fastp_dir}/$sample-fastp.json" \
        --report_title "$sample fastp report" \
        --unpaired1 ${fastp_dir}/00.unpaired.fasta.gz \
        --unpaired2 ${fastp_dir}/00.unpaired.fasta.gz \
        --average_qual 20 \
        --length_required 36 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --dedup \
        --thread 20 \
        --verbose 

done < doc/samples_to_use.txt


# Create directory for MultiQC report from fastp reports
fastp_mqc="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/01-QC/multiqc_fastp"

if [ ! -d "$fastp_mqc" ]; then
    mkdir -p "$fastp_mqc"
fi


# MultiQC report of the fastp reports
multiqc $fastp_dir \
    --outdir ${fastp_mqc} \
    --title "All fastp reports"

# End time and date
echo "$(date)       [End]"