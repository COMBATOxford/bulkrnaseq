#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N multiqc
#$ -q short.qc

module load MultiQC/1.9-foss-2019b-Python-3.7.4
cd ..
multiqc .
