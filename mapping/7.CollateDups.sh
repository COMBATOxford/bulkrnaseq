#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N mergeduplogs
#$ -q short.qc

source config.sh
cd $QC_DIR
module load R/3.6.2-foss-2019b
Rscript --vanilla ../scripts/7.CollateDups.R
