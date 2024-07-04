#!/bin/bash

# Start time and date
echo "$(date)       [Start]"

# Loops through CH and CO to generate .fq.gz files with all hostplant
# reads, respectively.
READS_TXT=doc/all-clean-reads.txt
R1_OUT=05-CLEAN-MERGED/all_R1-clean.fq.gz
R2_OUT=05-CLEAN-MERGED/all_R2-clean.fq.gz
while read R1 R2; do
    # Write R1 and R2 to new merged file
    cat ${R1} >> ${R1_OUT}
    cat ${R2} >> ${R2_OUT}
done < ${READS_TXT}


# Start time and date
echo "$(date)       [End]"