#!/bin/bash

#SBATCH --job-name seqkit-raw
#SBATCH --array=1-115%50
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 6
#SBATCH -t 00:45:00
#SBATCH --output=slurm-logs/decon-stats/SLURM-%j-seqkit-raw-%a.out
#SBATCH --error=slurm-logs/decon-stats/SLURM-%j-seqkit-raw-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load SeqKit/2.4.0

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Start time and date
echo "$(date)  $POP     [Start]"

# Get the sample
SAMPLE_LST=$(sed -n "${JOBID}p" doc/metadata-no-hybrids.csv)

# Get the project, sample and population. 
PROJ="${SAMPLE_LST%%_*}"
SAMPLE="${SAMPLE_LST%%,*}"
POP=$(echo $SAMPLE_LST | cut -d "," -f 2)

# Path to where .lst paths from.
RAWDIR="/crex/proj/snic2020-6-222/Projects/Tconura/data/WGS/rawdata/${PROJ}"

# Path to folder with .lst files
LST_DIR="00-RAW-LANES-LST"

# Create directory
OUT_DIR="results/raw-reads-stats/${POP}"
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

# This will generate a file with all the stats from multiple lanes
# Hence we will get multiple headers, it will be fixed after while loop.

# Temporary file to store stats in
TMP_OUT=${OUT_DIR}/tmp-${SAMPLE}-raw-stats.tsv

while read R1 R2; do

    # R1 and R2 with absolute path
    R1_IN="${RAWDIR}/${R1}"
    R2_IN="${RAWDIR}/${R2}"

    # Get stats
    seqkit stats \
        --basename \
        --tabular \
        --threads 6 \
        $R1_IN $R2_IN >> ${TMP_OUT}

done < ${LST_DIR}/${SAMPLE}.txt

# Final file name
FINAL_OUT=${OUT_DIR}/${SAMPLE}-raw-stats.tsv

# Add header to empty final file
cat ${TMP_OUT} | head -n 1 > ${FINAL_OUT}

# Discard the extra headers and only add stats data
cat ${TMP_OUT} | grep -v "^file" >> ${FINAL_OUT}

# Removes the temporary file
rm $TMP_OUT

# Start time and date
echo "$(date)  $POP     [End]"