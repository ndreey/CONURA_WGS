#!/bin/bash

#SBATCH --job-name anvio-profile
#SBATCH --array=1-2
#SBATCH -A naiss2024-22-580
#SBATCH -p core -n 12
#SBATCH -t 01:30:00
#SBATCH --output=slurm-logs/anvio/profile-db/SLURM-%j-profile-%a.out
#SBATCH --error=slurm-logs/anvio/profile-db/SLURM-%j-profile-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Activate the environment
mamba activate anvio-8

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

# Move to the anvio working directory
cd 07-ANVIO

# Paths and variables
PROF_DB=03-PROFILES/${POP}
BAM_IN=01-ALIGNMENT/${POP}.bam
CDB=02-CONTIG-DB/${POP}.db
DESC=../doc/short-project-desc.txt
MIN_LEN=1500
NUM_THREADS=12

anvi-profile \
    --input-file ${BAM_IN} \
    --contigs-db ${CDB} \
    --sample-name ${POP} \
    --description ${DESC} \
    --min-contig-length ${MIN_LEN} \
    --output-dir ${PROF_DB} \
    --num-threads ${NUM_THREADS}

# End time and date
echo "$(date)       ${POP}      [End]"