#!/bin/bash


####### OBS, you must activate the mamba env bbmap!

# Paths and variables
POP=$1
BINS=
OUT_DIR=results/MAG-assembly-stats
TMP_RES=${OUT_DIR}/tmp${POP}.tsv
FINAL_RES=${OUT_DIR}/${POP}-MAG-assembly-stats.csv

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

# bbmap stats.sh but for multiple files!
statswrapper.sh 05-metaWRAP/${POP}/*/*.fa format=5 \
    | grep -Ev "x10/metabat|x10/maxbin|x10/concoct|work_files" \
    > ${TMP_RES}

statswrapper.sh 05-metaWRAP/${POP}/*/*/*.fa format=5 \
    | grep -Ev "x10/metabat|x10/maxbin|x10/concoct|work_files|n_contigs" \
    >> ${TMP_RES}

# Converts tsv to csv, splits file name into binner and bin columns.
# inside awk an array is created from the file path. 
# the array is used to index out the binner and bin id.
cat ${TMP_RES} | awk '
BEGIN {
    FS = "\t"
    OFS = ","
    print "binner", "bin", "n_contigs", "contig_bp", "gap_pct", "ctg_N50", "ctg_L50", "ctg_N90", "ctg_L90", "ctg_max", "gc_avg", "gc_std"
}
NR > 1 {
    split($11, path_parts, "/")
    n = length(path_parts)
    bin = path_parts[n]
    binner_full = path_parts[n-1]
    split(binner_full, binner_parts, "_")
    binner = binner_parts[1]
    bin_name = substr(bin, 1, length(bin)-3)
    print binner, bin_name, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
}' > ${FINAL_RES}

# Removal of temp file.
rm ${TMP_RES}