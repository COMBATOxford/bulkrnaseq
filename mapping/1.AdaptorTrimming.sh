#!/bin/bash

# 1. Adaptor trimming: for each fastq file remove adaptor sequences

# Set parameters and load required modules
source config.sh
module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

# get the fastq file name for this task
FASTQ=$(cat fastq_key.txt | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $FASTQ
FASTQ1=$FASTQDIR/$FASTQ"_1.fastq.gz"
FASTQ2=$FASTQDIR/$FASTQ"_2.fastq.gz"
READS=$FASTQ1" "$FASTQ2

trim_galore --paired --gzip --output_dir $MAPPING_DIR/AdapterTrimmed/ $READS
