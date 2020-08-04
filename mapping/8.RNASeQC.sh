#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N rnaseqc
#$ -q short.qc

echo "------------------------------------------------"
echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}
echo SGE_TASK_FIRST=${SGE_TASK_FIRST}, SGE_TASK_LAST=${SGE_TASK_LAST}, SGE_TASK_STEPSIZE=${SGE_TASK_STEPSIZE}
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

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

~/rnaseqc.v2.3.5.linux $BASEDIR/"/Homo_sapiens.GRCh38.100.collapsed.gtf" $MAPPING_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam $DIR_SAMPLE_NAME/
