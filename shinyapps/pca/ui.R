suppressMessages(library(shiny))
suppressMessages(library(ggplot2))

load("pcs.Rdata")
load("covidpcs.Rdata")

shinyUI(fluidPage(
  
  # App title
  titlePanel("PCA on the COMBAT bulk RNA-seq data"),
  
  br(),
  
  fluidRow(
	   column(8,
		  fluidRow(
			   column(4, align="center",
           			radioButtons("pcx",
                        		label = "X axis:",
                        		choiceNames = list("PC1", "PC2", "PC3", "PC4", "PC5",
                                        	 "PC6", "PC7", "PC8", "PC9", "PC10"),
                        		choiceValues = list("PC1", "PC2", "PC3", "PC4", "PC5",
                                	             "PC6", "PC7", "PC8", "PC9", "PC10"),
                        		inline=T)),
    			column(4, align="center",
          			 radioButtons("pcy",
                        		label = "Y axis:",
                        		choiceNames = list("PC1", "PC2", "PC3", "PC4", "PC5",
                                        	     "PC6", "PC7", "PC8", "PC9", "PC10"),
                        		choiceValues = list("PC1", "PC2", "PC3", "PC4", "PC5",
                                              "PC6", "PC7", "PC8", "PC9", "PC10"),
                        		selected = "PC2",
                        		inline=T))
			),
		  fluidRow(
			   column(8, align="center",
                    radioButtons("data",
                                 label="Data subset",
                                 choices = list("All samples" = "All",
                                                "Hospitalised COVID" = "COVID"),
                                 selected="All"))
			   )),
		
		column(4, align="center",
		       		  radioButtons("variable",
                        label="Please select a covariate to plot",
                        choices = list("Age", "Sex", "Source", 
                                       "Timepoint", 
                                       "max_severity",
				       "Days_symptoms_admission",
                                       "Days_symptoms_discharge",
                                       "days_symptoms_max",
                                       "On dexamethasone at time of sample" = "On_dexamethasone_at_time_of_sample",
                                       "Days_symptom_to_sample",
                                       "PCR_SRS", "SRS",
                                       "Neutrophil proportion" = "neutrophil_prop", 
                                       "Lymphocyte proportion" = "lymphoctye_prop",
                                       "Monocyte proportion" = "monocyte_prop",
                                       "SaO2:FiO2 ratio" = "SaO2_FiO2_ratio",
                                       "Ventilation" = "ventilation_assistance",
                                       "Maximum Charlson Index" = "maximum_charlson_comorbidity_2012",
				       "BMI", "Weight"),
                        selected="Source", inline=T))
  ),
  
  br(),
  
  hr(),
  
  fluidRow(
    conditionalPanel(
      condition = "pcx != pcy",
      column(6, align="center", plotOutput("figPCA"))),
    
    # correlation plot
    column(6, plotOutput("figCorrplot"))),
    
    br(),
     
    # download button
    fluidRow(
      column(6, align="center", downloadButton("download_PCA", label="Download the plot")),
      column(6, align="center", downloadButton("download_CP", label="Download the plot"))
      ),
  
  br(),
  
  hr(),
  
  br()
))
