#!/bin/bash

# 4. Picard
# task for each sample

# Set parameters and load required modules
module load SAMtools/1.9-foss-2018b
MARKDUP=/path/to/gatk-4.0.3.0/gatk

# Get sample ID for this task
SAMPLE_NAME=$(cat $BASEDIR"/mapping.info.txt" | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $SAMPLE_NAME

# Make sample directory
DIR_SAMPLE_NAME=$QC_DIR$SAMPLE_NAME
if [[ ! -e $DIR_SAMPLE_NAME ]]; then
              mkdir -p $DIR_SAMPLE_NAME
fi

# Index bam fil
samtools index $MAPPING_DIR$SAMPLE_NAME/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam

# Quantify duplicates
$MARKDUP --java-options -Xmx8g MarkDuplicates -I $MAPPING_DIR$SAMPLE_NAME/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam -O $DIR_SAMPLE_NAME/$SAMPLE_NAME.dups.marked.bam -M $DIR_SAMPLE_NAME/${SAMPLE_NAME}.dup_metrics.txt

grep Unknown $DIR_SAMPLE_NAME/${SAMPLE_NAME}.dup_metrics.txt > $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup.counts.txt
rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.dups.marked.bam
rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup_metrics.txt
