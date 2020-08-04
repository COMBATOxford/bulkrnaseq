# Read in data
counts <- read.delim("../processeddata/Counts_143_23063.txt")
s.info <- read.delim("../processeddata/Sample_info_143.txt")
logcpm <- read.delim("../processeddata/Logcpm_143_23063.txt")
gene.info <- read.delim("../processeddata/Gene_info_23063.txt")
load("SRSmodel.Rdata")

# Extract genes for SRS model
srs.genes <- c("DYRK2", "CCNB1IP1", "TDRD9", "ZAP70", "ARL14EP", "MDC1", "ADGRE3")
test.data <- data.frame(t(logcpm[gene.info$gene_name %in% srs.genes, ]))
colnames(test.data) <- gene.info$gene_name[match(colnames(test.data), gene.info$gene_id)]
test.data$SRS <- s.info$PCR_SRS

# Assign samples to SRS groups
new.samples <- predict(model, newdata=test.data, type="response")
test.data$SRS <- new.samples
s.info$SRSScore <- test.data$SRS
s.info$SRSModel <- as.factor((s.info$SRSScore > 0.5) + 1)

# Write out results
combat.sepsis.srs <- s.info[s.info$Source == "Sepsis",
                            c("RNASeq_sample_ID", "PCR_SRS", "SRSModel")]
write.table(combat.sepsis.srs, "SRSAssignmentsSepsis.txt", 
            sep="\t", quote=F, row.names=F)

