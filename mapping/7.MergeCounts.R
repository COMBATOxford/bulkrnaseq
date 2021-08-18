#!/usr/bin/env Rscript

# Make a list of the full paths to the final feature counts for all samples
basedir <- "./"
s.info <- read.delim("../mapping.info.txt", sep="", header=FALSE)
sample.names <- as.character(s.info$V1)
myfiles <- paste(basedir, sample.names, "/", sample.names, ".counts.txt", sep="")

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  tryCatch({
    DT[[myfiles[i]]] <- read.table(myfiles[i], header = T, stringsAsFactors = FALSE)
    colnames(DT[[myfiles[i]]]) <- c("ID", sample.names[i])
  }, error=function(e) print(myfiles[i]))
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# now add each additional table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[, -1]

# write data to file
write.table(data, file = "Full_count_data.txt", sep="\t")
