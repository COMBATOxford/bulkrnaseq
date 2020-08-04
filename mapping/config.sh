#!/bin/bash

######### Set parameters and paths
BASEDIR="/well/combat/projects/rnaseq/CBD-RNA-00001/processing/"
echo "Base directory: "$BASEDIR
cd $BASEDIR

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

# Sample information file
KEY=${BASEDIR}P200133_fastq_key.txt
echo "Sample information file: "$KEY

# Fastq files (or softlinks to them)
FASTQDIR="/well/combat/projects/rnaseq/CBD-RNA-00001/rawdata/"
echo "Fastq directory: "$FASTQDIR

# Genome references (build 38)
GTF="/well/combat/shared/references/GRCh38/Annotation/Genes/genes.gtf"
GENOME="/well/combat/shared/references/GRCh38/STARindex.101"
