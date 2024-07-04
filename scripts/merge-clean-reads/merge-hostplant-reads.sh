#!/bin/bash

# Start time and date
echo "$(date)       [Start]"

# Loops through CH and CO to generate .fq.gz files with all hostplant
# reads, respectively.
for HOSTP in CH CO; do
    echo "Merging reads for ${HOSTP} hostplant"
    READS_TXT=doc/all-${HOSTP}-clean-reads.txt
    R1_OUT=05-CLEAN-MERGED/${HOSTP}_R1-clean.fq.gz
    R2_OUT=05-CLEAN-MERGED/${HOSTP}_R2-clean.fq.gz

    while read R1 R2; do
        # Write R1 and R2 to new merged file
        cat ${R1} >> ${R1_OUT}
        cat ${R2} >> ${R2_OUT}
    done < ${READS_TXT}
done

# Start time and date
echo "$(date)       [End]"