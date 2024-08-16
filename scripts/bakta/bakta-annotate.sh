#!/bin/bash

### OBS, bakta has to be activated through mamba.

cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
WRAP_DIR=05-metaWRAP/${POP}
CHECKM2_DIR=06-CHECKM2/results
NUM_THREADS=16
OUT_DIR=010-BAKTA/${POP}
BAKTA_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/bakta-db/db
TMP_CHKM=${OUT_DIR}/tmp${POP}.tsv

# Ensure the output directory exists
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi  

# Prepare the temporary file
cat ${CHECKM2_DIR}/${POP}-bin-report.tsv | \
    cut -f 1,2,3 | grep -v Binner > ${TMP_CHKM}

# Loop through each line of the temporary file
while read -r BINNER BIN COMP; do

    # Check if completeness is less than 50.0
    if (( $(echo "${COMP} < 50.0" | bc -l) )); then
        continue
    fi

    # Determine the bin path based on the binner
    if [ "${BINNER}" == "metawrap" ]; then
        BIN_PATH=${WRAP_DIR}/bin_refinement_c50_x10/metawrap_50_10_bins/${BIN}.fa
    else
        BIN_PATH=${WRAP_DIR}/${BINNER}_bins/${BIN}.fa
    fi

    # Run bakta
    bakta \
        --db ${BAKTA_DB} \
        --prefix ${BINNER}-${BIN} \
        --output ${OUT_DIR}/${BINNER}-${BIN} \
        --keep-contig-headers \
        --threads ${NUM_THREADS} \
        ${BIN_PATH}

done < ${TMP_CHKM}

# Clean up
rm ${TMP_CHKM}