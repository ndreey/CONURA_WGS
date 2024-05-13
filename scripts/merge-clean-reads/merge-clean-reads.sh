#!/bin/bash

#SBATCH --job-name merge-clean
#SBATCH -A naiss2023-22-412
#SBATCH --array=1-10
#SBATCH -p core -n 2
#SBATCH -t 04:35:00
#SBATCH --output=SLURM-%j-merge-clean-%a.out
#SBATCH --error=SLURM-%j-merge-clean-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date)       [Start]"

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}

# Variables
CLEAN_READS="04-CLEAN-FASTQ"

# Population folder
POP=$(ls -1 $CLEAN_READS/ | sed -n "${JOBID}p")

# Output directory
OUT_DIR="05-CLEAN-MERGED"
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p "${OUT_DIR}"
fi

# Write R1 and R2 to new merged file
cat ${CLEAN_READS}/${POP}/*R1* > ${OUT_DIR}/${POP}_R1-clean.fq.gz
cat ${CLEAN_READS}/${POP}/*R2* > ${OUT_DIR}/${POP}_R2-clean.fq.gz