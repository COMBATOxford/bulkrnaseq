

suppressMessages(library(shiny))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(limma))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggrepel))
suppressMessages(library(DT))

load("DEfit.Rdata")
load("DEfitage.Rdata")
load("DEfitagesex.Rdata")
load("DEfitagesexcells.Rdata")
load("DEfitagesexcharlson.Rdata")

comparator.groups <- c("Healthy Volunteers" = "HV",
                          "COVID (Health Care Workers)" = "COVID_HCW",
                          "COVID (Mild)" = "COVID_MILD",
                          "COVID (Severe)" = "COVID_SEV",
                          "COVID (Critical)" = "COVID_CRIT",
                          "Sepsis" = "(Sepsis_SRS1 + Sepsis_SRS2)/2",
                       "SRS1" = "Sepsis_SRS1",
                       "SRS2" = "Sepsis_SRS2")

# Define server logic
shinyServer(function(input, output, session) {

  updateSelectizeInput(session, "geneSelector",
                       choices = gene.info,
                       selected="CD177",
                       options = list(labelField='gene_name', searchField='gene_name',
                                      valueField='gene_name',
                                      render = I("{ option: function(item, escape) {return '<div><strong>' + escape(item.gene_name) + '</span></div>';} }")),
                       server = TRUE)
  
  # DE output
  n.sig.table <- function(groupA, groupB){
    if(groupA != input$group2  &
       !(groupA == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS1') &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS2') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS1') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS2')){
	  
	  comp <- paste(groupA, "vs", groupB, sep="")
    my.contr <- paste(groupA, "-", groupB, sep="")
    
    if(input$covariates == "None"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=design)
      colnames(contrast.matrix) <- comp  
      fit2 <- contrasts.fit(fit, contrast.matrix)
    } else if(input$covariates == "Age"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designage)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitage, contrast.matrix)
    } else if(input$covariates == "AgeSex"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designagesex)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitagesex, contrast.matrix)
    } else if(input$covariates == "AgeSexCharlson"){
	    contrast.matrix <- makeContrasts(contrasts=my.contr,
                                             levels=designagesexcharlson)
	    colnames(contrast.matrix) <- comp
            fit2 <- contrasts.fit(fitagesexcharlson, contrast.matrix)
    } else if(input$covariates == "AgeSexCells"){
	    contrast.matrix <- makeContrasts(contrasts=my.contr,
					     levels=designagesexcells)
	    colnames(contrast.matrix) <- comp
	    fit2 <- contrasts.fit(fitagesexcells, contrast.matrix)
    }
    
    fit2 <- eBayes(fit2, trend=TRUE)
    
    de <- topTable(fit2, number = nrow(logcpm), coef=comp)
    de$Gene <- gene.info$gene_name[match(rownames(de), gene.info$gene_id)]
    de$Sig <- de$adj.P.Val < 0.05 & abs(de$logFC) >= log2(as.numeric(input$FCthreshold))
    results <- data.frame(table(de$Sig))
    colnames(results) <- c(paste("Significant (FDR<0.05 and FC>=", input$FCthreshold, ")", sep=""), "Number of genes")
    return(results)
    }
  }
  
  output$table1 <- renderTable(n.sig.table(input$group1, input$group2))
  
  volcano.plot <- function(groupA, groupB){
    if(groupA != input$group2  &
       !(groupA == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS1') &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS2') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS1') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS2')){

	namegp1 <- names(comparator.groups[comparator.groups == groupA])
    	namegp2 <- names(comparator.groups[comparator.groups == groupB])
    	comp <- paste(groupA, "vs", groupB, sep="")
    	my.contr <- paste(groupA, "-", groupB, sep="")
    
	if(input$covariates == "None"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=design)
      colnames(contrast.matrix) <- comp  
      fit2 <- contrasts.fit(fit, contrast.matrix)
    } else if(input$covariates == "Age"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designage)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitage, contrast.matrix)
    } else if(input$covariates == "AgeSex"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designagesex)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitagesex, contrast.matrix)
    } else if(input$covariates == "AgeSexCharlson"){
            contrast.matrix <- makeContrasts(contrasts=my.contr,
                                             levels=designagesexcharlson)
            colnames(contrast.matrix) <- comp
            fit2 <- contrasts.fit(fitagesexcharlson, contrast.matrix)
    }else if(input$covariates == "AgeSexCells"){
            contrast.matrix <- makeContrasts(contrasts=my.contr,
                                             levels=designagesexcells)
            colnames(contrast.matrix) <- comp
            fit2 <- contrasts.fit(fitagesexcells, contrast.matrix)
    }
    
    fit2 <- eBayes(fit2, trend=TRUE)
    
    de <- topTable(fit2, number = nrow(logcpm), coef=comp)
    de$Gene <- gene.info$gene_name[match(rownames(de), gene.info$gene_id)]
    de$Significant <- de$adj.P.Val < 0.05 & abs(de$logFC) >= log2(as.numeric(input$FCthreshold))
    
    ggplot(de, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(col=Significant)) +
      scale_color_manual(values=c("grey", "red")) +
      ylab("-log10(adjusted p value)") +
      xlab("log2(fold change)") +
      theme_bw() +
      geom_text_repel(data=subset(de, Significant == TRUE)[1:20, ], aes(label=Gene), fontface = "italic") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed") +
      geom_vline(xintercept=log2(as.numeric(input$FCthreshold)), linetype="dashed") +
      geom_vline(xintercept=-(log2(as.numeric(input$FCthreshold))), linetype="dashed") +
      ggtitle(paste(namegp1, "vs", namegp2), subtitle=paste("+ve FC indicates upregulation in", namegp1))
    }
  }
  
  output$figVP <- renderPlot({
    req(input$group1, input$group2)
    volcano.plot(groupA = input$group1, groupB = input$group2)
  })
  
  # Download the VP figure as a pdf
  output$download_VP <- downloadHandler(
    filename = function() {
      gsub(" ", "_", paste(names(comparator.groups[comparator.groups == input$group1]), 
            "vs", names(comparator.groups[comparator.groups == input$group2]),
	    "_VP", ".pdf", sep=""))
    },
    content = function(file){
      pdf(file, useDingbats=FALSE)
      print(volcano.plot(groupA = input$group1, groupB = input$group2))
      dev.off()
    })
  
  sig.genes <- function(groupA, groupB){
    if(groupA != input$group2  &
       !(groupA == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'HV' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupB == 'COVID_HCW' & input$covariates %in% c('AgeSexCells', 'AgeSexCharlson')) &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS1') &
       !(groupA == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupB == 'Sepsis_SRS2') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS1') &
       !(groupB == '(Sepsis_SRS1 + Sepsis_SRS2)/2' & groupA == 'Sepsis_SRS2')){

	 comp <- paste(groupA, "vs", groupB, sep="")
   	 my.contr <- paste(groupA, "-", groupB, sep="")
    if(input$covariates == "None"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=design)
      colnames(contrast.matrix) <- comp  
      fit2 <- contrasts.fit(fit, contrast.matrix)
    } else if(input$covariates == "Age"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designage)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitage, contrast.matrix)
    } else if(input$covariates == "AgeSex"){
      contrast.matrix <- makeContrasts(contrasts=my.contr, 
                                       levels=designagesex)
      colnames(contrast.matrix) <- comp
      fit2 <- contrasts.fit(fitagesex, contrast.matrix)
    } else if(input$covariates == "AgeSexCharlson"){
            contrast.matrix <- makeContrasts(contrasts=my.contr,
                                             levels=designagesexcharlson)
            colnames(contrast.matrix) <- comp
            fit2 <- contrasts.fit(fitagesexcharlson, contrast.matrix)
    } else if(input$covariates == "AgeSexCells"){
            contrast.matrix <- makeContrasts(contrasts=my.contr,
                                             levels=designagesexcells)
            colnames(contrast.matrix) <- comp
            fit2 <- contrasts.fit(fitagesexcells, contrast.matrix)
    }
    
    fit2 <- eBayes(fit2, trend=TRUE)
    
    de <- topTable(fit2, number = nrow(logcpm), coef=comp, sort.by = "P")
    de$Gene <- gene.info$gene_name[match(rownames(de), gene.info$gene_id)]
    de$Sig <- de$adj.P.Val < 0.05 & abs(de$logFC) >= log2(as.numeric(input$FCthreshold))
    de <- de[, c("Gene", "AveExpr", "logFC", "P.Value", "adj.P.Val", "Sig")]
    de$AveExpr <- round(de$AveExpr, digits=2)
    de$logFC <- round(de$logFC, digits=2)
    de$P.Value <- signif(de$P.Value, digits=3)
    de$adj.P.Val <- signif(de$adj.P.Val, digits=3)
    de$GeneType <- gene.info$gene_biotype[match(rownames(de), gene.info$gene_id)]
    
    return(de)
    }
  }
  
  # Table of DE results
  output$table2 <- renderDT(sig.genes(input$group1, input$group2),
                                   options=list(pageLength=10),
				   filter=list(position='bottom', plain=TRUE))

  # Download the DE results as a .txt file
  output$download_detable <- downloadHandler(
    filename = function() {
      gsub(" ", "_", paste(names(comparator.groups[comparator.groups == input$group1]), 
            "vs", names(comparator.groups[comparator.groups == input$group2]),
	   "_DEResults_", input$covariates, ".txt", sep=""))
    },
    content = function(file){
      write.table(sig.genes(input$group1, input$group2), file, quote=F, sep="\t")
    })
  
  # Function for plotting DE figure
  de.plot <- function(gene, variable){
    
    gene.id <- gene.info$gene_id[gene.info$gene_name == gene]
    exprn <- data.frame("Source"=s.info$Source,
                        "Expression"=as.numeric(logcpm[gene.id, ]),
                        "DaysofSymptoms"=s.info$Days_symptom_to_sample,
                        "Sex"=s.info$Sex,
                        "Age"=s.info$Age,
			"Neutr_prop"=s.info$Neutr_prop,
                        "Mono_prop"=s.info$Mono_prop,
                        "Lymph_prop"=s.info$Lymph_prop,
                        "SaO2_FiO2_ratio"=s.info$SaO2_FiO2_ratio,
                        "Ventilation_assistance"=s.info$ventilation_assistance)
    exprn$Ventilation_assistance <- factor(exprn$Ventilation_assistance,
                                           levels=c("nil", "NIV", "intubated"))
    if(variable == "Age"){
      
      ggplot(exprn, aes(Source, Expression, colour=Age)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitter(width=0.2)) +
        scale_colour_gradient(low="blue", high="red") +
        xlab(label="Source") +
        scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
        theme_bw() +
        ggtitle("Differential Gene Expression (Limma, all samples)")
      
    } else if(variable == "DaysofSymptoms"){
      
      ggplot(exprn, aes(Source, Expression, colour=DaysofSymptoms)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitter(width=0.2)) +
        scale_colour_gradient(low="darkblue", high="red") +
        xlab(label="Source") +
        scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
        theme_bw() +
        ggtitle("Differential Gene Expression (Limma, all samples)")
      
    } else if(variable == "Sex"){
      
      ggplot(exprn, aes(Source, Expression)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(aes(colour=Sex), position=position_jitter(width=0.2)) +
        xlab(label="Source") +
        scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
        theme_bw() +
        ggtitle("Differential Gene Expression (Limma, all samples)")
    } else if(variable == "Neutr_prop"){
          ggplot(exprn, aes(Source, Expression, colour=Neutr_prop)) +
            geom_boxplot(outlier.shape=NA) +
            geom_point(position=position_jitter(width=0.2)) +
            scale_colour_gradient2(low="darkblue", mid="yellow", high="red",
				   midpoint=median(exprn$Neutr_prop, na.rm=T)) +
            xlab(label="Source") +
            scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
            theme_bw() +
            ggtitle("Differential Gene Expression (Limma, all samples)")
        } else if(variable == "Mono_prop"){
          ggplot(exprn, aes(Source, Expression, colour=Mono_prop)) +
            geom_boxplot(outlier.shape=NA) +
            geom_point(position=position_jitter(width=0.2)) +
            scale_colour_gradient2(low="darkblue", mid="gold", high="red",
                                   midpoint=median(exprn$Mono_prop, na.rm=T)) +
            xlab(label="Source") +
            scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
            theme_bw() +
            ggtitle("Differential Gene Expression (Limma, all samples)")
        } else if(variable == "Lymph_prop"){
          ggplot(exprn, aes(Source, Expression, colour=Lymph_prop)) +
            geom_boxplot(outlier.shape=NA) +
            geom_point(position=position_jitter(width=0.2)) +
            scale_colour_gradient2(low="darkblue", mid="gold",high="red",
				   midpoint=median(exprn$Lymph_prop, na.rm=T)) +
            xlab(label="Source") +
            scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers", 
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
	    ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
            theme_bw() +
            ggtitle("Differential Gene Expression (Limma, all samples)")
        } else if(variable == "Ventilation_assistance"){
		ggplot(exprn, aes(Source, Expression)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(aes(colour=Ventilation_assistance), position=position_jitter(width=0.2)) +
        xlab(label="Source") +
        scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers",
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
            ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
        theme_bw() +
        ggtitle("Differential Gene Expression (Limma, all samples)")
  	} else if(variable == "SaO2_FiO2_ratio"){
	ggplot(exprn, aes(Source, Expression)) +
            geom_boxplot(outlier.shape=NA) +
            geom_point(aes(colour=SaO2_FiO2_ratio), position=position_jitter(width=0.2)) +
            scale_colour_gradient2(low="darkblue", mid="gold", high="red", 
				   midpoint=median(exprn$SaO2_FiO2_ratio, na.rm=T)) +
            xlab(label="Source") +
            scale_x_discrete(labels=c("HV" = "Healthy\n Volunteers",
                                      "COVID_HCW" = "COVID \n(Health Care Workers)",
                                      "COVID_MILD" = "COVID \n(Mild)",
                                      "COVID_SEV" = "COVID \n(Severe)",
                                      "COVID_CRIT" = "COVID \n(Critical)",
                                      "Sepsis_SRS1" = "Sepsis (SRS1)",
                                      "Sepsis_SRS2" = "Sepsis (SRS2)")) +
            ylab(label=paste("Log2(cpm + 1)", gene, "Expression")) +
            theme_bw() +
            ggtitle("Differential Gene Expression (Limma, all samples)")
    }
  }
  
  # Plot DE figure
  output$figBoxplot <- renderPlot({
    req(input$geneSelector, input$variable)
    de.plot(input$geneSelector, input$variable)
  })
  
  
  # Download the boxplot figure as a pdf
  output$download_BP <- downloadHandler(
    filename = function() {
      gsub(" ", "_", paste(input$geneSelector, "_boxplot", ".pdf", sep=""))
    },
    content = function(file){
      pdf(file, useDingbats=FALSE, width=9, height=6)
      print(de.plot(input$geneSelector, input$variable))
      dev.off()
    })

})
