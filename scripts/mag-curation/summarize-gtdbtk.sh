#!/bin/bash


# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

		#### Combine quality reports into a single file

# Create the checkm2 directory
RES_DIR=07-GTDB-Tk/results
if [ ! -d "${RES_DIR}" ]; then
    mkdir -p ${RES_DIR}
fi 

# Paths for final report
# Paths and variables
POP=$1
TMP_REPORT=${RES_DIR}/${POP}-gtdb-report.tmp
FINAL_REPORT=${RES_DIR}/${POP}-gtdb-report.tsv
GT_DIR=07-GTDB-Tk/${POP}

# Initialize the final report file
> "${TMP_REPORT}"

for BINNER in concoct_bins maxbin2_bins metabat2_bins metawrap_50_10_bins; do
	for DOMAIN in bac120 ar53; do
        # Variables
	    QUALITY_REPORT=${GT_DIR}/${BINNER}/gtdbtk.${DOMAIN}.summary.tsv
	    B_NAME="${BINNER%%_*}"

        if [ -f "${QUALITY_REPORT}" ]; then
            # Add BINNER column and append to the final report file
            awk -v binner="${B_NAME}" 'NR==1 && FNR==1 {print "Binner\t" $0} NR>1 {if(FNR>1) print binner "\t" $0}' "${QUALITY_REPORT}" >> "${TMP_REPORT}"
        else
            echo "Warning: ${QUALITY_REPORT} does not exist. Skipping..."
        fi
    done
done

# Final touches of the combined quality and removal of temporary file
cat ${TMP_REPORT} | cut --complement -f17 | head -1 > ${FINAL_REPORT}
cat ${TMP_REPORT} | grep -v Binner | cut --complement -f17 >> ${FINAL_REPORT}
rm ${TMP_REPORT}

# End time and date
echo "$(date)       ${POP}      [End]"
