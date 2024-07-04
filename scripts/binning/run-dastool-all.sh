#!/bin/bash

# Activate the environment
#mamba activate dastool

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
POP=all
WRAP_DIR=05-metaWRAP
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
NUM_THREADS=8
CONT2BIN=contigs2bin/${POP}
SCORE_THRESH=0.1

                #### DASTOOL ####

# Run DASTOOL with these score thresholds.
# Set variable for output
DAS_DIR=06-DASTOOL/${POP}
DAS_OUT=${DAS_DIR}/score-thresh-${SCORE_THRESH}

# Create the DASTOOL directory
if [ ! -d "${DAS_DIR}" ]; then
    mkdir -p ${DAS_DIR}
fi

# Get the contig2bins.tsv files as tab delimited textform.
C2Bs=$(ls ${CONT2BIN}/* | paste - - - | tr "\t" ",")

# Run the program
DAS_Tool \
    --bins ${C2Bs} \
    --contigs ${ASSEMBLY} \
    --outputbasename ${DAS_OUT} \
    --labels concoct,maxbin2,metabat2 \
    --score_threshold ${SCORE_THRESH} \
    --write_bin_evals \
    --write_bins \
    --threads ${NUM_THREADS}


# End time and date
echo "$(date)       ${POP}      [End]"