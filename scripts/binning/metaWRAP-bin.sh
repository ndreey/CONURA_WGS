#!/bin/bash

#SBATCH --job-name metaWRAP
#SBATCH --array=1-2
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1
#SBATCH -C mem256GB	
#SBATCH -t 12:00:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-metaWRAP-binning-%a.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-metaWRAP-binning-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load metaWRAP/1.3.2

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
NUM_THREADS=16
OUT_DIR=04-metaWRAP/${POP}
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
READS=../05-CLEAN-MERGED/${POP}*

# Generate folder for metaWRAP
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

metawrap binning \
    -a ${ASSEMBLY} \
    -o ${OUT_DIR} \
    -t ${NUM_THREADS} \
    -m 256 \
    --metabat2 --maxbin2 --concoct \
    ${READS}

# End time and date
echo "$(date)  ${POP}    [End]"