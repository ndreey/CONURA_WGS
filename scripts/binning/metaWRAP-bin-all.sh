#!/bin/bash

#SBATCH --job-name metaWRAP-all
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1
#SBATCH -C mem256GB	
#SBATCH -t 08:00:00
#SBATCH --output=slurm-logs/binning/SLURM-%j-metaWRAP-binning-all.out
#SBATCH --error=slurm-logs/binning/SLURM-%j-metaWRAP-binning-all.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in modules
module load bioinfo-tools
module load metaWRAP/1.3.2

# Start time and date
echo "$(date)       ${POP}     [Start]"

# Move to the anvio working directory
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/08-ANVIO-ALL

# Paths and variables
POP=all
NUM_THREADS=16
OUT_DIR=05-metaWRAP/${POP}
ASSEMBLY=00-FIXED-ASSEMBLY/${POP}-contigs.fa
R1=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/05-CLEAN-MERGED/${POP}_R1-clean.fq.gz
R2=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/05-CLEAN-MERGED/${POP}_R2-clean.fq.gz
TMP_R1=${OUT_DIR}/tmp-reads/${POP}_1.fastq
TMP_R2=${OUT_DIR}/tmp-reads/${POP}_2.fastq

# Generate folder for metaWRAP + tmp-reads
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi
mkdir -p ${OUT_DIR}/tmp-reads/

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
rm -r ${OUT_DIR}/tmp-reads/

# End time and date
echo "$(date)  ${POP}    [End]"