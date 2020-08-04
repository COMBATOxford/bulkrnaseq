#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N mergestarlogs
#$ -q short.qc

source config.sh
cd $MAPPING_DIR
module load R/3.6.2-foss-2019b
Rscript --vanilla ../scripts/6.CollateSTARLogs.R
