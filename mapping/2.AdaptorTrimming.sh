#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N trimming
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
module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

FASTQ=$(cat $KEY | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $FASTQ
FASTQ1=$FASTQDIR/$FASTQ"_1.fastq.gz"
FASTQ2=$FASTQDIR/$FASTQ"_2.fastq.gz"
READS=$FASTQ1" "$FASTQ2

trim_galore --paired --gzip --output_dir $MAPPING_DIR/AdapterTrimmed/ $READS
