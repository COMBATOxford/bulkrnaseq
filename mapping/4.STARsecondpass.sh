#!/bin/bash
#$ -cwd
#$ -pe shmem 8 -N star2ndpass
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
module load STAR/2.7.3a-GCC-8.3.0

# Get sample ID for this task
SAMPLE_NAME=$(cat $BASEDIR"/mapping.info.txt" | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo $SAMPLE_NAME

# Make sample directory
DIR_SAMPLE_NAME=$MAPPING_DIR$SAMPLE_NAME
if [[ ! -e $DIR_SAMPLE_NAME ]]; then
              mkdir -p $DIR_SAMPLE_NAME
fi

# Get fastq file names
FASTQ=$(cat $BASEDIR"/mapping.info.txt" | tail -n+${SGE_TASK_ID} | head -1 | cut -f2 )

FASTQ1="$(echo "$FASTQ" | sed -e "s|,|_1_val_1.fq.gz,$MAPPING_DIR/AdapterTrimmed/\/|g")"
FASTQ2="$(echo "$FASTQ" | sed -e "s|,|_2_val_2.fq.gz,$MAPPING_DIR/AdapterTrimmed/\/|g")"
FASTQ1=$FASTQ1"_1_val_1.fq.gz"
FASTQ2=$FASTQ2"_2_val_2.fq.gz"
FASTQ1=$MAPPING_DIR/AdapterTrimmed/$FASTQ1
FASTQ2=$MAPPING_DIR/AdapterTrimmed/$FASTQ2
echo $FASTQ1
echo $FASTQ2

STAR --genomeDir $GENOME \
      --runThreadN 6 \
      --readFilesIn $FASTQ1 $FASTQ2 \
      --readFilesCommand gunzip -c \
      --outFileNamePrefix $DIR_SAMPLE_NAME"/"${SAMPLE_NAME}. \
      --outSAMtype BAM SortedByCoordinate \
      --limitSjdbInsertNsj 10000000 \
      --sjdbGTFfile $GTF \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --sjdbFileChrStartEnd $MAPPING_DIR/*.SJ.out.tab \
      --outFilterType BySJout \
      --outReadsUnmapped Fastx

rm $DIR_SAMPLE_NAME"/"${SAMPLE_NAME}.Aligned.out.bam
