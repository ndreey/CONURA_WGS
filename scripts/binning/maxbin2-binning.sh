#!/bin/bash

#SBATCH --job-name MaxBin2
#SBATCH --array=1-2
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1
#SBATCH -C mem256GB	
#SBATCH -t 08:00:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-maxbin2-%a.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-maxbin2-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load MaxBin/2.2.7

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

run_MaxBin.pl \
    -contig ${ASSEMBLY} \
    