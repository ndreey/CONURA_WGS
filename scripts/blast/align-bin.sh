#!/bin/bash

module load bioinfo-tools
# module load bowtie2/2.5.2
module load bwa/0.7.17
# module load samtools/1.19

# Set new work dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-ANVIO

# Paths and variables
BIN=05-metaWRAP/CHST/concoct_bins/bin.16.fa 
QUERY=../doc/stammerula1.tephritidis-16s.fa
OUT_DIR=08-BLAST-RESULTS/bam
REF_BIN=${OUT_DIR}/bin.16.fa
NUM_THREADS=12

cp ${BIN} ${OUT_DIR}/

# Build the BWA index
bwa index ${REF_BIN}

# Align the query sequence to the reference genome using BWA-MEM and sort the output by position
bwa mem -a -t ${NUM_THREADS} ${REF_BIN} ${QUERY} | \
    samtools view - -b -@ ${NUM_THREADS} | \
    samtools sort - -@ ${NUM_THREADS} -o ${OUT_DIR}/bin16.bam


#bowtie2-build --threads ${NUM_THREADS} ${BIN} ${OUT_DIR}/bin16
#
## Align the query sequence to the reference genome and sort the output by position
#bowtie2 -x ${OUT_DIR}/bin16 -f ${QUERY} | \
#    samtools view - -b -@ ${NUM_THREADS} | \
#    samtools sort - -@ ${NUM_THREADS} -o ${OUT_DIR}/bin16.bam
#
## Index the sorted BAM file (optional but recommended for further analysis)
#samtools index ${OUT_DIR}/bin16.bam