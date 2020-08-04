#!/usr/bin/env Rscript

# Make a sample-fastq key for use downstream

sample.key <- read.delim("../P200133_fastq_key.txt", stringsAsFactors=FALSE, header=F)
# This is provided with each lane of data received and looks like this:
# WTCHG_123456_123456	S1 
# WTCHG_123456_234567	S2
# WTCHG_123456_345678	S3

# All samples were pooled and sequenced across 4 lanes
# There are therefore 4 fastq files for each sample

mapping.key <- matrix(nrow=length(unique(sample.key[, 2])), ncol=2)
mapping.key[, 1] <- unique(sample.key[, 2])
for(i in 1:nrow(mapping.key)){
  files <- c(sample.key[, 1][sample.key[, 2] == mapping.key[i, 1]])
  mapping.key[i, 2] <- paste(files, collapse=",")
}

mapping.key <- data.frame(mapping.key)

write.table(mapping.key, "mapping.info.txt", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
