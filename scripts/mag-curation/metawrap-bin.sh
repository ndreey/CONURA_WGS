#!/bin/bash

#SBATCH --job-name metaWRAP-CHST
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1
#SBATCH -C mem256GB	
#SBATCH -t 05:00:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-metaWRAP-binning-CHST.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-metaWRAP-binning-CHST.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load metaWRAP/1.3.2

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

# Paths and variables
POP=$1
NUM_THREADS=16
OUT_DIR=05-metaWRAP/${POP}
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}/${POP}-contigs.fa
R1=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/05-CLEAN-MERGED/${POP}_R1-clean.fq.gz
R2=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/05-CLEAN-MERGED/${POP}_R2-clean.fq.gz
TMP_R1=${OUT_DIR}/${POP}-tmp-reads/${POP}_1.fastq
TMP_R2=${OUT_DIR}/${POP}-tmp-reads/${POP}_2.fastq

# Generate folder for metaWRAP
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

# Generate a temp directory inside metawrap folder.
mkdir -p ${OUT_DIR}/${POP}-tmp-reads/

# metawrap requires *_1.fastq as format. Hence, tmp versions of files
zcat ${R1} > ${TMP_R1}
zcat ${R2} > ${TMP_R2}

metawrap binning \
    -a ${ASSEMBLY} \
    -o ${OUT_DIR} \
    -t ${NUM_THREADS} \
    -m 256 \
    --metabat2 --maxbin2 --concoct \
    ${TMP_R1} ${TMP_R2}

# Remove the temporary files
rm -r ${OUT_DIR}/${POP}-tmp-reads/

# End time and date
echo "$(date)  ${POP}    [End]"