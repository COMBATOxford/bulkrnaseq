#!/usr/bin/env Rscript

library(data.table)
# Collate results
log.paths <- list.files(path="./", pattern="*.Aligned.sortedByCoord.out.bam.metrics.tsv",
                  recursive=TRUE, full.names=TRUE)
sample.names <- list.files(path="./", pattern="*.Aligned.sortedByCoord.out.bam.metrics.tsv",
                           recursive=TRUE, full.names=FALSE)
sample.names <- gsub(".Aligned.sortedByCoord.out.bam.metrics.tsv", "", sample.names)
logs <- lapply(log.paths, read.delim, header=TRUE, row.names=1)
logs <- lapply(logs, t)
logs <- lapply(logs, as.data.frame)

qc.logs <- rbindlist(logs)
rownames(qc.logs) <- sample.names

write.table(qc.logs, "RNA-SeqQC_metrics.txt", sep="\t", quote=F)
