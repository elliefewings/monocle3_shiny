#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(monocle3)
library(stringr)
library(SeuratWrappers)
library(mclust)
library(genesorteR)
library(pheatmap)
library(dqshiny)

# Change data load maximum
options(shiny.maxRequestSize = 800*1024^2)

#Source functions
source("./source.R")

# Define UI for application
ui <- shinyUI(fluidPage(
                        # Set style for vertical alignment
                        tags$head(tags$style(
                          HTML('
                               #vert { 
                               display: flex;
                               align-items: center;
                               margin-top: 50px;
                               }
                               .tooltip .tooltip-inner {
                               max-width: 100%;
                               }
                               '))),
                        shinyjs::useShinyjs(),
                        
                        # Application title
                        titlePanel(textOutput("title")),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Create Input
                        fluidRow(column(2, wellPanel(
                        fileInput("file", 'Choose Rdata/Rds to upload',
                                    accept = c('.Rdata', ".rds", ".Rds")),
                        #Set load button
                        actionButton("upload", "Upload"),
                        
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Set assay type
                        selectInput("assay",
                                    "Select assay",
                                    NULL),
                        
                        #Set data identity
                        selectInput("ident",
                                    "Select Identity for clustering",
                                    NULL),
                        
                        #Set progenitor gene
                        autocomplete_input("gene", "Progenitor gene:", NULL),
                        
                        #Set load button
                        actionButton("run", "Run Monocle3"),
                        
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        #Select plot to download
                        selectInput("plot",
                                    "Select plot to download",
                                    c("Dimensions plot", "Monocle3 results", "Heatmap")),
                        #Download plot
                        downloadButton("download", "Download Plot"),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        #Download dta
                        downloadButton("download2", "Download Results"))),
                        
                        #Plot dim and pseudotime
                        column(5, align="center", id="vert", plotOutput("dim", width="90%", height="800px")), 
                        column(5, id="vert", align="center", plotOutput("p2", width="90%", height="800px"))),
                        
                        # Heatmap
                        fluidRow(tags$hr(style="border-color: black;")),
                        fluidRow(column(10, align="center", offset=2, id="vert", h3("Heatmap of genes significantly associated with pseudotime"))),
                        fluidRow(column(10, align="center", offset=2, id="vert", plotOutput("p5", width="80%", height="700px"))),
                        fluidRow(tags$hr(style="border-color: black;"))
                        ))

# Define server logic
server <- shinyServer(function(input, output, session) {
  
  #Set text outputs
  output$title <- renderText("Pseudotime projection with Monocle3")
  
  # Load data
  load_Rdata <- function(){
    if(is.null(input$file)){return(NULL)} 
    rdata <- isolate({input$file})
    short <- input$file
    
    if(grepl(".Rdata", short$name)) {
    load(rdata$datapath, envir = .GlobalEnv)
    data <- get(names(which(unlist(eapply(.GlobalEnv, is.Seurat)))))
    } else {
    data <- readRDS(rdata$datapath)
    }
    return(data)
  }
  
  #Write function to find Seurat object
  is.Seurat <- function(x){
    any(class(x) == "Seurat")
  }
  
  #Disable buttons until Rdata is loaded
  shinyjs::disable("assay")
  shinyjs::disable("ident")
  shinyjs::disable("gene")
  shinyjs::disable("run")
  shinyjs::disable("plot")
  shinyjs::disable("download")
  
  observeEvent(input$upload,{
    data <- load_Rdata()
    
    #Update assay list
    updateSelectInput(session, "assay", "Select assay", choices = names(data@assays))
    
    #Update identity list
    updateSelectInput(session, "ident", "Select Identity for clustering", choices = colnames(data@meta.data))
    
    #Update autocomplete for gene searching
    genes <- unique(rownames(data))
    update_autocomplete_input(session, "gene",
                              options = genes)
    
    #Enable buttons
    shinyjs::enable("gene")
    shinyjs::enable("assay")
    shinyjs::enable("ident")
    shinyjs::enable("run")
  
  # Create event when report load button is activated
  observeEvent(input$run,{
    
    #Set assay
    DefaultAssay(data) <- input$assay
    
    #Create progress bar
    withProgress(message = 'Running Monocle3...please wait', {
    out <- mon.run(data=data, gene=input$gene, id=input$ident)
    })
    
    output$dim <- renderPlot({out$dim})
    
    output$p2 <- renderPlot({out$p2})
    
    output$p5 <- renderPlot({out$p5})
    
    #Enable buttons
    shinyjs::enable("plot")
    shinyjs::enable("download")

    #Create multiplot file needed for saving

      savePlot <- function(){
        pnames <- list("Dimensions plot"=out$dim, "Monocle3 results"=NA, "Heatmap"=out$p5)
        pnames[input$plot][[1]]
      }
      
      #Create action for download button
      output$download <- downloadHandler(
        filename = function() {
          fnames <- c("Dimensions plot"="clustering_dimensions.pdf", "Monocle3 results"=paste("monocle3_pseudotime_", input$gene, ".pdf", sep=""), "Heatmap"="monocle3_heatmap.pdf")
          fnames[input$plot][[1]]
        },
        content = function(file) {
          
          if (input$plot == "Dimensions plot") {
            pdf(file)
            print(out$dim)
            dev.off()  
          } else if (input$plot == "Monocle3 results"){
            pdf(file)
            print(out$p1)
            print(out$p2)
            print(out$p3)
            print(out$p4)
            dev.off()
          } else if (input$plot == "Heatmap"){
            pdf(file, width=15, height=9)
            print(out$p5)
            dev.off()
          }
        })
      
      #Create action for download of data
      output$download2 <- downloadHandler(
        filename = function(){"monocle3_pseudotime_results.csv"},
        content = function(file){
          write.csv(out$t1, file, row.names = FALSE)
          })
  
      })
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

