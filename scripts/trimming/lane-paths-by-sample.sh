#!/bin/bash

# Start time and date
echo "$(date)       [Start]"


IFS=","; while read SAMPLE POP HP REGION RANGE; do

    echo "Processing: $SAMPLE"
    # Get the project name.
    # %%_* removes everything after the first underscore
    PROJ="${SAMPLE%%_*}"

    # Path from where the .lst files path from
    RAWDIR="/crex/proj/snic2020-6-222/Projects/Tconura/data/WGS/rawdata/${PROJ}"

    # Path to .lst file
    LSTFILE="${RAWDIR}/${SAMPLE}.lst"

    # Create a temporary folder for sorted .lst files
    TMPLST="00-RAW-LANES-LST"
    if [ ! -d "$TMPLST" ]; then
        mkdir -p "$TMPLST"
    fi

    # Sorts and removes .md5 files from .lst
    cat $LSTFILE | grep -v ".md5" | sort | paste - - > ${TMPLST}/${SAMPLE}.txt


done < doc/metadata-no-hybrids.csv

# Store the sample.lst's into a text file so we can iterate in trimming script.
ls -1 $TMPLST > doc/sample_lst_lane_paths.txt

# Start time and date
echo "$(date)       [End]"