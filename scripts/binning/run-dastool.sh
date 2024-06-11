#!/bin/bash

# Activate the environment
mamba activate dastool

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
POP=$1
NF_DIR=04-BINNING
WRAP_DIR=05-metaWRAP
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
NUM_THREADS=12
NF_C2B=contigs2bin/nf-core
WRAP_C2B=contigs2bin/metaWRAP


                #### DASTOOL ####

# Running the DASTOOL refinement for each tool
for C2B_DIR in ${NF_C2B} ${WRAP_C2B}; do

    # Run DASTOOL with these score thresholds.
    for SCORE in 0.1 0.2 0.3 0.4 0.5 0.6 0.7; do

        # Set variable for output
        DAS_OUT_PROG=$(basename ${C2B_DIR})
        DAS_DIR=06-DASTOOL/${DAS_OUT_PROG}
        DAS_OUT=${DAS_DIR}/${POP}/score-thresh-${SCORE}

        # Create the DASTOOL directory
        if [ ! -d "${DAS_DIR}" ]; then
            mkdir -p ${DAS_DIR}
        fi

        # Get the contig2bins.tsv files as tab delimited textform.
        C2Bs=$(ls ${C2B_DIR}/${POP}* | paste - - - | tr "\t" ",")

        # Run the program
        DAS_Tool \
            --bins ${C2Bs} \
            --contigs ${ASSEMBLY} \
            --outputbasename ${DAS_OUT} \
            --labels concoct,maxbin2,metabat2 \
            --score_threshold ${SCORE} \
            --write_bins \
            --write_bin_evals \
            --threads ${NUM_THREADS}
    done
done

# End time and date
echo "$(date)       ${POP}      [End]"