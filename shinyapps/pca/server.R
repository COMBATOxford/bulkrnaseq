suppressMessages(library(shiny))
suppressMessages(library(ggplot2))

load("pcs.Rdata")
load("covidpcs.Rdata")

# Define server logic
shinyServer(function(input, output, session) {
  
pca.plot <- function(pcx, pcy, variable){
    if(pcx != pcy){
	if(input$data == "All"){
        	df <- data.frame("PCX"=as.numeric(pcs[, pcx]),
                         "PCY"=as.numeric(pcs[, pcy]),
                         "variable"=pcs[, variable])
      	pc_vc.n <- pc_vc.n.full
	} else if(input$data == "COVID"){
        	df <- data.frame("PCX"=as.numeric(covidpcs[, pcx]),
                         "PCY"=as.numeric(covidpcs[, pcy]),
                         "variable"=covidpcs[, variable])
      	pc_vc.n <- pc_vc.n.covid
	}

      if(class(df$variable) == "numeric"){
        ggplot(df, aes(PCX, PCY)) +
          geom_point(aes(col=variable)) +
	  xlab(paste0(pcx, "(", pc_vc.n[pcx], "%)")) +
          ylab(paste0(pcy, "(", pc_vc.n[pcy], "%)")) +
	  theme_bw() +
          labs(colour = paste(variable)) +
          scale_colour_gradient2(low="red", mid="yellow", high="blue",
                                 midpoint = median(df$variable, na.rm = T))
      } else {
        ggplot(df, aes(PCX, PCY)) +
          geom_point(aes(col=variable)) +
          xlab(paste0(pcx, "(", pc_vc.n[pcx], "%)")) +
          ylab(paste0(pcy, "(", pc_vc.n[pcy], "%)")) +
	  theme_bw() +
          labs(colour = paste(variable))
      }
    }
  }
  
  output$figPCA <- renderPlot({
    req(input$pcx, input$pcy, input$variable)
    pca.plot(input$pcx, input$pcy, input$variable)
  })
  
  # Download the PCA figure as a pdf
  output$download_PCA <- downloadHandler(
    filename = function() {
      paste(input$pcx, "_", input$pcy, "_", input$variable, ".pdf", sep="")
    },
    content = function(file){
      pdf(file, useDingbats=FALSE, width=9, height=6)
      print(pca.plot(input$pcx, input$pcy, input$variable))
      dev.off()
    })
  
  # Function for plotting correlation figure
  corr.plot <- function(pcx, variable){
    if(input$data == "All"){
      df <- data.frame("PCX"=as.numeric(pcs[, pcx]),
                       "Source"=pcs$Source,
                       "variable"=pcs[, variable])
    } else if(input$data == "COVID"){
      df <- data.frame("PCX"=as.numeric(covidpcs[, pcx]),
                       "Source"=covidpcs$Source,
                       "variable"=covidpcs[, variable])
    }

  if(class(df$variable) == "numeric"){
        ggplot(df, aes(PCX, variable)) +
          geom_point(aes(col=Source)) +
          xlab(pcx) +
          ylab(variable) +
          theme_bw() +
          geom_smooth(se=FALSE)

      } else {
        ggplot(df, aes(PCX, variable)) +
          geom_point(aes(col=Source)) +
          xlab(pcx) +
          ylab(variable) +
          theme_bw()
      }
    }
  
  # Plot correlation figure
  output$figCorrplot <- renderPlot({
    req(input$pcx, input$variable)
    corr.plot(input$pcx, input$variable)
  })
  
  # Download the correlation plot figure as a pdf
  output$download_CP <- downloadHandler(
    filename = function() {
      paste(input$pcx, "_", input$variable, "_corrplot", ".pdf", sep="")
    },
    content = function(file){
      pdf(file, useDingbats=FALSE, width=8, height=6)
      print(corr.plot(input$pcx, input$variable))
      dev.off()
    })
  
}
)
