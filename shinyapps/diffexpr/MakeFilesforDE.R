# make files for DE shiny

library(limma)

counts <- read.delim("processeddata/Counts_143_23063.txt")
s.info <- read.delim("processeddata/Full_sample_info_143.txt")
logcpm <- read.delim("processeddata/Logcpm_143_23063.txt")
gene.info <- read.delim("processeddata/Gene_info_23063.txt")

# basic model
design <- model.matrix(~ 0 + Source, data = s.info)
colnames(design) <- gsub("Source", "", colnames(design))
corfit <- duplicateCorrelation(logcpm, design, block=s.info$baseID)
corcons <- corfit$consensus # 0.3701581
fit <- lmFit(logcpm, design, block=s.info$baseID, correlation=corcons)
save(logcpm, s.info, gene.info, fit, design, file="DEfit.Rdata")

# Age
designage <- model.matrix(~ 0 + Source + Age + AgeSquared, data = s.info)
colnames(designage) <- gsub("Source", "", colnames(designage))
corfit <- duplicateCorrelation(logcpm, designage, block=s.info$COMBAT_ID)
corcons <- corfit$consensus # 0.371435
fitage <- lmFit(logcpm, designage, block=s.info$baseID, correlation=corcons)
save(fitage, designage, file="DEfitage.Rdata")

# Age and sex
designagesex <- model.matrix(~ 0 + Source + Age + AgeSquared + Sex, data = s.info)
colnames(designagesex) <- gsub("Source", "", colnames(designagesex))
corfit <- duplicateCorrelation(logcpm, designagesex, block=s.info$baseID)
corcons <- corfit$consensus # 0.3761294
fitagesex <- lmFit(logcpm, designagesex, block=s.info$baseID, correlation=corcons)
save(fitagesex, designagesex, file="DEfitagesex.Rdata")

# Age, sex, cell counts
designagesexcells <- model.matrix(~ 0 + Source + Age + AgeSquared + Sex + 
                                    Neutr_prop + Mono_prop, data = s.info)
colnames(designagesexcells) <- gsub("Source", "", colnames(designagesexcells))
fitagesexcells <- lmFit(logcpm[, !is.na(s.info$Neutr_prop)], designagesexcells, 
                        block=s.info$baseID[!is.na(s.info$Neutr_prop)], correlation=corcons)
save(fitagesexcells, designagesexcells, file="DEfitagesexcells.Rdata")

# Charlson
designagesexcharlson <- model.matrix(~ 0 + Source + Age + AgeSquared + Sex + 
                                       Charlson, data = s.info)
colnames(designagesexcharlson) <- gsub("Source", "", colnames(designagesexcharlson))
corfit <- duplicateCorrelation(logcpm[, !(is.na(s.info$Charlson))], designagesexcharlson, block=s.info$baseID[!(is.na(s.info$Charlson))])
corcons <- corfit$consensus # 0.4167762
fitagesexcharlson <- lmFit(logcpm[, !(is.na(s.info$Charlson))], designagesexcharlson, 
                           block=s.info$baseID[!(is.na(s.info$Charlson))], correlation=corcons)
save(fitagesexcharlson, designagesexcharlson, file="DEfitagesexcharlson.Rdata")
