# Variance Partition
# Following vignette from the variancePartition package

# load library
library(variancePartition)

# load data:
geneExpr <- read.delim("data/Logcpm_143_23063.txt")
info <- read.delim("sampleinfo/Full_sample_info_143.txt")
info$Source[info$Source %in% c("Sepsis_SRS1", "Sepsis_SRS2")] <- "Sepsis"
rownames(info) <- gsub("-", ".", info$RNASeq_sample_ID)
all(rownames(info) == colnames(geneExpr))

# geneExpr <- geneExpr[, info$Source %in% c("COVID_MILD", "COVID_SEV", "COVID_CRIT")]
# info <- info[info$Source %in% c("COVID_MILD", "COVID_SEV", "COVID_CRIT"), ]
geneExpr <- geneExpr[, !duplicated(info$baseID)]
info <- info[!duplicated(info$baseID), ]

# Compute Canonical Correlation Analysis (CCA) between all pairs of variables
# returns absolute correlation value
form <- ~ Age + Sex + Source + SRSModel + Days_symptom_to_sample + 
  Neutr_prop + Mono_prop + Lymph_prop + SaO2_FiO2_ratio + ventilation_assistance + 
  Charlson + max_severity
C <- canCorPairs( form, info)
# Plot correlation matrix
plotCorrMatrix( C )

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual is categorical, so model as random effects
# form <- ~ Age + (1|Sex) + (1|Source) + SRSModel + Days_symptom_to_sample +
#   Neutr_prop + Mono_prop + SaO2_FiO2_ratio + # (1|ventilation_assistance) +
#   Charlson + (1|baseID)
form <- ~ Age + (1|Sex) + (1|Source) + SRSModel + Days_symptom_to_sample +
  Neutr_prop + Mono_prop + SaO2_FiO2_ratio + Charlson

# check colinearity
res <- fitVarPartModel( geneExpr[1:4,], form, info )
colinearityScore( res[[1]] )

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified, a linear mixed model is used
# If all variables are modeled as fixed effects, a linear model is used
# each entry in results is a regression model fit on a single gene

# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix, and EList output by voom() in the limma package,
# or an ExpressionSet
varPart <- fitExtractVarPartModel(geneExpr, form, info )
head(varPart)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 20 genes
plotPercentBars( vp[1:20, ] )

# violin plot of contribution of each variable to total variance
# pdf("variance_partition_COVID_all_samples.pdf", useDingbats=F)
plotVarPart( vp )
# dev.off()

# get gene with the highest variation across Source
# create data.frame with expression of gene i and Source
# type for each sample
i <- which.max( varPart$Source )
GE <- data.frame( "Expression" = t(geneExpr[i, ]), 
                  "Source" = as.factor(info$Source))
colnames(GE) <- c("Expression", "Source")

# plot expression stratified by source
plotStratify( Expression ~ Source, GE, main=rownames(geneExpr)[i],
              x.labels = TRUE, legend = FALSE)

# get gene with the highest variation across Individuals
# create data.frame with expression of gene i and individual type for each sample
i <- which.max( varPart$baseID )
GE <- data.frame( Expression = t(geneExpr[i,]),
                  Individual = info$baseID)
colnames(GE) <- c("Expression", "Individual")

# plot expression stratified by individual
label <- paste("Individual:", format(varPart$baseID[i]*100,
                                     digits=3), "%")
main <- rownames(geneExpr)[i]
plotStratify( Expression ~ Individual, GE, colorBy=NULL,
              text=label, main=main)


# get list of coagulation-related genes
coag <- read.delim("Coagulationgenesets.txt")
coag <- unique(coag$Reactome.term.Formation.of.Fibrin.Clot..Clotting.Cascade.)
gene.info <- read.delim("../processeddata/Gene_info_23063.txt")
coag <- geneExpr[gene.info$gene_name %in% coag, ]
rownames(coag) <- gene.info$gene_name[match(rownames(coag), gene.info$gene_id)]
varPart <- fitExtractVarPartModel(coag, form, info )
head(varPart)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 20 genes
plotPercentBars( vp[1:30,] )

# violin plot of contribution of each variable to total variance
pdf("variance_partition_coagulation_genes.pdf", useDingbats = F)
plotVarPart( vp )
dev.off()

