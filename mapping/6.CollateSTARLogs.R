#!/usr/bin/env Rscript

# Combine STAR logs
log.paths <- list.files(path=".", pattern="*Log.final.out", full.names=TRUE, recursive=TRUE)
sample.names <- list.files(path=".", pattern="*Log.final.out", full.names=FALSE, recursive=TRUE)
sample.names <- gsub(".Log.final.out", "", sample.names)
logs <- lapply(log.paths, read.delim, row.names=1, header=FALSE)
logs <- lapply(logs, as.data.frame)

cbindlist <- function(list) {
  n <- length(list)
  res <- list[[1]]
  for (i in 2:n) res <- cbind(res, list[[i]])
  return(res)
}

for (i in 1:length(logs)){
  colnames(logs[[i]]) <- sample.names[i]
}

mapping.stats <- cbindlist(logs)
mapping.stats <- t(mapping.stats)
colnames(mapping.stats) <- gsub(" ", "", colnames(mapping.stats))
colnames(mapping.stats) <- gsub("\\|", "", colnames(mapping.stats))
mapping.stats <- data.frame(mapping.stats)
head(mapping.stats)

mapping.stats$Numberofinputreads <- as.numeric(as.character(mapping.stats$Numberofinputreads))
mapping.stats$Uniquelymappedreadsnumber <- as.numeric(as.character(mapping.stats$Uniquelymappedreadsnumber))
mapping.stats$Numberofsplices.Total <- as.numeric(as.character(mapping.stats$Numberofsplices.Total))
mapping.stats$Numberofreadsmappedtomultipleloci <- as.numeric(as.character(mapping.stats$Numberofreadsmappedtomultipleloci))
head(mapping.stats)

write.table(mapping.stats, "mapping_stats.txt", sep="\t", quote=F)
