#!/bin/bash



# Activate the environment
mamba activate anvio-8

POP=$1

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
PROF_DB=03-PROFILES/${POP}
BAM_IN=01-ALIGNMENT/${POP}.bam
CDB=02-CONTIG-DB/${POP}.db
DESC=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/doc/short-project-desc.txt
MIN_LEN=1500
NUM_THREADS=12

anvi-profile \
    --input-file ${BAM_IN} \
    --contigs-db ${CDB} \
    --sample-name ${POP} \
    --description ${DESC} \
    --min-contig-length ${MIN_LEN} \
    --output-dir ${PROF_DB} \
    --num-threads ${NUM_THREADS} \
    --cluster-contigs

# End time and date
echo "$(date)       ${POP}      [End]"