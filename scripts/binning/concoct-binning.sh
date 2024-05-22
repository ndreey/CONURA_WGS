#!/bin/bash

#SBATCH --job-name concoct
#SBATCH --array=1-2
#SBATCH -A naiss2024-22-580
#SBATCH -p core -n 12
#SBATCH -t 05:30:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-concoct-%a.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-concoct-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load CONCOCT/1.1.0

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

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Paths and variables
NUM_THREADS=20
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
BAM_FILE=01-ALIGNMENT/${POP}.bam
OUT_DIR=04-BINNING/CONCOCT/${POP}
BED_FILE=${OUT_DIR}/${POP}-contigs_10K.bed
FA_10K=${OUT_DIR}/${POP}-contigs_10K.fa
COV_TBL=${OUT_DIR}/${POP}-coverage_table.tsv

# Move to the anvio working directory
cd 07-ANVIO

# Generate folder for binning
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

# Cut contigs into 1kbp bits
cut_up_fasta.py \
    ${ASSEMBLY} \
    --chunk_size 10000 \
    --overlap_size 0 \
    --merge_last \
    --bedfile ${BED_FILE} > ${FA_10K}

# Generate table of coverage depth for each subcontig
concoct_coverage_table.py ${BED_FILE} ${BAM_FILE} > ${COV_TBL}

# Cluster
concoct \
    --coverage_file ${COV_TBL} \
    --composition_file ${FA_10K} \
    --threads ${NUM_THREADS} \
    --basename ${OUT_DIR}

# Merge the subcontigs
merge_cutup_clustering.py \
    ${OUT_DIR}/clustering_gt1000.csv > ${OUT_DIR}/clustering_merged.csv

# Create directory for bins
mkdir ${OUT_DIR}/bins

# Extract bins as individual FASTA files
extract_fasta_bins.py \
    ${ASSEMBLY} \
    ${OUT_DIR}/clustering_merged.csv \
    --output_path ${OUT_DIR}/bins

