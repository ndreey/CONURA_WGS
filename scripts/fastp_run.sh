#!/bin/bash

#SBATCH -job-name fastp
#SBATCH -A naiss2024-5-1
#SBATCH -p core -n 8
#SBATCH -t 02:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err


module load bioinf-tools
module load fastp/0.23.4