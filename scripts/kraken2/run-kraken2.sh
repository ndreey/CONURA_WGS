#!/bin/bash

#SBATCH --job-name Kraken2
#SBATCH -A naiss2024-22-580
#SBATCH -p node -n 1
#SBATCH -C mem1TB							
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/kraken2/SLURM-%j-kraken2.out
#SBATCH --error=slurm-logs/kraken2/SLURM-%j-kraken2.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load Kraken2/2.1.3-20231102-acc2248

# Move to anvio working dir
cd /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/07-MAG

OUT_DIR=012-KRAKEN2-k2_nt_20231129

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p ${OUT_DIR}
fi

POP=$1
KRAK_DB=/sw/data/Kraken2_data/prebuilt/k2_nt_20231129/
OUTPUT=${OUT_DIR}/${POP}-kraken2.output
CONF=$4
REPORT=${OUT_DIR}/${POP}-kraken2.report
NUM_THREADS=12
R1=$2
R2=$3

kraken2 \
    --db ${KRAK_DB} \
    --threads ${NUM_THREADS} \
    --output ${OUTPUT} \
    --confidence ${CONF} \
    --report ${REPORT} \
    --paired \
    --use-names \
    ${R1} ${R2}