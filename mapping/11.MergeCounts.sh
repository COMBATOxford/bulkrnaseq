#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N mergecounts
#$ -q short.qc

source config.sh
cd $COUNT_DIR
module load R/3.6.2-foss-2019b
Rscript --vanilla ../scripts/11.MergeCounts.R
