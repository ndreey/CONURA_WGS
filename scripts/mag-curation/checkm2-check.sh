#!/bin/bash

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
WRAP_DIR=05-metaWRAP/${POP}
NUM_THREADS=12
CHECKM2_DB=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/checkm2-db/CheckM2_database/uniref100.KO.1.dmnd


for BINNER in concoct_bins maxbin2_bins metabat2_bins metawrap_50_10_bins; do
    
    # Adjust path for bin-refinment directory
    if [ ${BINNER} == "metawrap_50_10_bins" ]; then
        BIN_DIR=${WRAP_DIR}/bin_refinement_c50_x10/${BINNER}
    else
        BIN_DIR=${WRAP_DIR}/${BINNER}
    fi

    # Create the checkm2 directory
    OUT_DIR=06-CHECKM2/${POP}/${BINNER}
    if [ ! -d "${OUT_DIR}" ]; then
        mkdir -p ${OUT_DIR}
    fi    

    # Run the MAG quality prediction
    checkm2 predict \
        --input ${BIN_DIR} \
        --output-directory ${OUT_DIR} \
        --threads ${NUM_THREADS} \
        --general \
        --extension fa \
        --ttable 11 \
        --database_path ${CHECKM2_DB}
done

		#### Combine quality reports into a single file

# Create the checkm2 directory
RES_DIR=06-CHECKM2/results
if [ ! -d "${RES_DIR}" ]; then
    mkdir -p ${RES_DIR}
fi 

# Paths for final report
TMP_REPORT=${RES_DIR}/${POP}-bin-report.tmp
FINAL_REPORT=${RES_DIR}/${POP}-bin-report.tsv
CHECKM_DIR=06-CHECKM2/${POP}

# Initialize the final report file
> "${TMP_REPORT}"

for BINNER in concoct_bins maxbin2_bins metabat2_bins metawrap_50_10_bins; do
	# Variables
	QUALITY_REPORT=${CHECKM_DIR}/${BINNER}/quality_report.tsv
	B_NAME="${BINNER%%_*}"

    if [ -f "${QUALITY_REPORT}" ]; then
        # Add BINNER column and append to the final report file
        awk -v binner="${B_NAME}" 'NR==1 && FNR==1 {print "Binner\t" $0} NR>1 {if(FNR>1) print binner "\t" $0}' "${QUALITY_REPORT}" >> "${TMP_REPORT}"
    else
        echo "Warning: ${QUALITY_REPORT} does not exist. Skipping..."
    fi
done

# Final touches of the combined quality and removal of temporary file
cat ${TMP_REPORT} | cut --complement -f5,12 | head -1 > ${FINAL_REPORT}
cat ${TMP_REPORT} | grep -v Binner | cut --complement -f5,12 >> ${FINAL_REPORT}
rm ${TMP_REPORT}

# End time and date
echo "$(date)       ${POP}      [End]"