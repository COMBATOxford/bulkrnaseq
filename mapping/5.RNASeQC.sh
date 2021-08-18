#!/bin/bash

# 5. RNA-SeQC
# task for each sample

# Set parameters
source config.sh

# Get sample ID for this task
SAMPLE_NAME=$(cat $BASEDIR"/mapping.info.txt" | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $SAMPLE_NAME

# Make sample directory
DIR_SAMPLE_NAME=$QC_DIR$SAMPLE_NAME
if [[ ! -e $DIR_SAMPLE_NAME ]]; then
              mkdir -p $DIR_SAMPLE_NAME
fi

rnaseqc.v2.3.5.linux $BASEDIR/"/Homo_sapiens.GRCh38.100.collapsed.gtf" $MAPPING_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam $DIR_SAMPLE_NAME/
