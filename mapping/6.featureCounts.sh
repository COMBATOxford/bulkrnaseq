#!/bin/bash

# 6. Generate feature counts
# task for each sample

# Set parameters and load required modules
source config.sh
module load Subread/1.6.4-foss-2018b

# Get sample ID for this task
SAMPLE_NAME=$(cat $BASEDIR"/mapping.info.txt" | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $SAMPLE_NAME

# Make sample directory
DIR_SAMPLE_NAME=$COUNT_DIR$SAMPLE_NAME
if [[ ! -e $DIR_SAMPLE_NAME ]]; then
              mkdir -p $DIR_SAMPLE_NAME
fi

featureCounts -T 4 -a $GTF -g gene_id \
      -o ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt -p \
      -s 2 $MAPPING_DIR$SAMPLE_NAME/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam

cut -f 1,7 ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt > ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt   #  This
mv ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt             #  reduces the file size from ~ 30M to ~1M
