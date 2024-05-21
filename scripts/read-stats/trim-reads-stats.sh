#!/bin/bash

#SBATCH --job-name seqkit-trim
#SBATCH --array=1-115%50
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 6
#SBATCH -t 00:35:00
#SBATCH --output=slurm-logs/trim-reads/SLURM-%j-seqkit-trim-%a.out
#SBATCH --error=slurm-logs/trim-reads/SLURM-%j-seqkit-trim-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load SeqKit/2.4.0

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Start time and date
echo "$(date)  $POP     [Start]"

# Path to trim reads
TRIM_DIR="02-TRIM"

# Get the sample
SAMPLE_LST=$(sed -n "${JOBID}p" doc/metadata-no-hybrids.csv)

# Get the project, sample and population. 
PROJ="${SAMPLE_LST%%_*}"
SAMPLE="${SAMPLE_LST%%,*}"
POP=$(echo $SAMPLE_LST | cut -d "," -f 2)

# Create directory
OUT_DIR="results/trim-reads-stats/${POP}"
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

seqkit stats \
    --basename \
    --tabular \
    --threads 6 \
    ${TRIM_DIR}/${SAMPLE}*.fastq.gz >> ${OUT_DIR}/${SAMPLE}-trim-stats.tsv

# End time and date
echo "$(date)  $POP     [End]"