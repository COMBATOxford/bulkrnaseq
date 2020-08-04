#!/bin/bash

# sh ./RunPipeline.sh

##### make a key to link sample ID to multiple fastq files
# echo "Making sample key"
# module load R/3.6.2-foss-2019b
# Rscript --vanilla ./1.MakeSampleKey.R $KEY
# module purge

##### trim adaptor sequences
qsub -t 1-576 -P combat.prjc ./2.AdaptorTrimming.sh

##### alignment
echo "Submitting jobs for alignment"

qsub -hold_jid "trimming" -t 1-144 -P combat.prjc 3.STARfirstpass.sh
qsub -hold_jid "star1stpass" -t 1-144 -P combat.prjc 4.STARsecondpass.sh

##### QC
echo "Submitting jobs for QC"

qsub -hold_jid "star2ndpass" -t 1-144 -P combat.prjc 5.Picard.sh
qsub -hold_jid "star2ndpass" -P combat.prjc 6.CollateSTARLogs.sh
qsub -hold_jid "markdups" -P combat.prjc 7.CollateDups.sh
qsub -hold_jid "star2ndpass" -t 1-144 -P combat.prjc 8.RNASeQC.sh
qsub -hold_jid "rnaseqc" -P combat.prjc 9.CollateQCMetrics.sh

##### Counts
echo "Submitting jobs for counts"

qsub -hold_jid "star2ndpass" -t 1-144 -P combat.prjc 10.featureCounts.sh
qsub -hold_jid "counts" -P combat.prjc 11.MergeCounts.sh

qsub -hold_jid "counts" -P combat.prjc 12.multiqc.sh
