#!/bin/bash

######### Set parameters and paths
BASEDIR="/path/to/project/directory/"
echo "Base directory: "$BASEDIR
cd $BASEDIR
# contains fastq information file: list of fastq files and the corresponding sample ID
# contains mapping key: list of sample IDs and all the corresponding fastq files

##### make other subdirectories
if [[ ! -e "Mapping" ]]; then
mkdir ./Mapping
fi
MAPPING_DIR=$BASEDIR"Mapping/"

if [[ ! -e "QC" ]]; then
mkdir ./QC
fi
QC_DIR=$BASEDIR"QC/"

if [[ ! -e "Counts" ]]; then
mkdir ./Counts
fi
COUNT_DIR=$BASEDIR"Counts/"

# Fastq files (or softlinks to them)
FASTQDIR="/path/to/rawdata/"
echo "Fastq directory: "$FASTQDIR

# Genome references (build 38)
GTF="/path/to/annotation/file/genes.gtf"
GENOME="/path/to/genome/index/GRCh38/STARindex.101"
