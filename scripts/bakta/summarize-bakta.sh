#!/bin/bash

# Start time and date
echo "$(date)        [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

		#### Combine quality reports into a single file

# Create the bakta directory
RES_DIR=010-BAKTA/results
if [ ! -d "${RES_DIR}" ]; then
    mkdir -p ${RES_DIR}
fi 

# Paths and variables
POP=$1
BAK_REP=${RES_DIR}/${POP}-bakta-report.csv
BAK_DIR=010-BAKTA/${POP}

# Build header
HEADER="Binner,Bin,23S,15S,5S,rRNA,tRNA,ori"

echo ${HEADER} > ${BAK_REP}

BINS=$(ls ${BAK_DIR}/)

for BIN in ${BINS}; do
    # Get the tsv and binner, bin variables
    BAK_TSV=${BAK_DIR}/${BIN}/${BIN}.tsv

    # Extract everything after the first "-"
    BIN_ID="${BIN#*-}"
    # Extract everything befoe the first "-"
    BINNER="${BIN%%-*}"

    # number tRNA
    NUM_tRNA=$(cat ${BAK_TSV} | grep tRNA- | grep -v cds | wc -l)

    # Number rRNA
    NUM_rRNA=$(cat ${BAK_TSV} | grep -v cds | grep rRNA | wc -l)

    # Number 16S
    NUM_16S=$(cat ${BAK_TSV}| grep -v cds | grep rRNA | grep 16S | wc -l)

    # Number 23S
    NUM_23S=$(cat ${BAK_TSV} | grep -v cds | grep rRNA | grep 23S | wc -l)

    # Number 5S
    NUM_5S=$(cat ${BAK_TSV} | grep -v cds | grep rRNA | grep 5S | wc -l)

    # Number of origin of replication
    NUM_ori=$(cat ${BAK_TSV}| grep "ori[CVT]" | wc -l)

    # Write to file
    echo "${BINNER},${BIN_ID},${NUM_23S},${NUM_16S},${NUM_5S},${NUM_rRNA},${NUM_tRNA},${NUM_ori}" \
        >> ${BAK_REP}
done

# End time and date
echo "$(date)       ${POP}      [End]"