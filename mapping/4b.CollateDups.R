#!/usr/bin/env Rscript

library(data.table)

# Collate Picard results
log.paths <- list.files(path=".", pattern="*.dup.counts.txt", full.names=TRUE, recursive=T)
sample.names <- list.files(path=".", pattern="*.dup.counts.txt", full.names=FALSE, recursive=T)
sample.names <- gsub(".dup.counts.txt", "", sample.names)

logs <- lapply(log.paths, read.delim, header=FALSE)
logs <- lapply(logs, as.data.frame)
dup.logs <- rbindlist(logs)
dup.logs <- as.data.frame(dup.logs)
rownames(dup.logs) <- sample.names

write.table(dup.logs, "Duplication_Metrics.txt", sep="\t", quote=FALSE)
