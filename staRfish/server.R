suppressPackageStartupMessages(library(cBioPortalData, quietly = T))
suppressPackageStartupMessages(library(shiny, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(cbioportalR, quietly = T))
suppressPackageStartupMessages(library(magrittr, quietly = T))

source("../R/grab.R")
source("../R/digest.R")
source("../R/squeeze.R")
source("../R/release_metadata.R")

function(input, output, session) {

  ## Studies table
  all_studies <- reactive({input$chosen_fields %>%
    strsplit(split = ",") %$%
    squeeze(n = input$n_samples, unlist(.))
    })

  output$study_table <- DT::renderDT({all_studies()})
  study_ids <- reactive({all_studies() %>% pull(studyId)})

  observe({
    updateSelectInput(session, "select_study", choices =   study_ids())
  })

  ## Metadata table
  output$metadata_table <- DT::renderDT({release_metadata(input$select_study)})
  # output$metadata_table <- metadata_table

  output$download_metadata <- downloadHandler(
    filename = "metadata.tsv",
    content = function(fname){
      write_tsv(release_metadata(input$select_study), fname)
    })


  ## Correlation plot
  data_plot <- reactive({
    if(input$gene_name %in% unique(iris$Species)){
      iris %>% filter(Species == input$gene_name)
    } else {
        iris
      }
    })

  output$plot_correlation <- renderPlot({plot(data_plot()$Sepal.Length, data_plot()$Sepal.Width)})



  ### ADD REAL PLOT CODE
  output$download_corr_plot <- downloadHandler(
    filename = "placeholder_filename.png",
    content = function(filename) {
      png(filename)
      plot(data_plot()$Sepal.Length, data_plot()$Sepal.Width)
      dev.off()
    })

  ## Placeholder plot
  output$plot_correlation2 <- renderPlot({plot(iris$Sepal.Length, iris$Sepal.Width)})

  ### ADD REAL PLOT CODE
  output$download_corr_plot2 <- downloadHandler(
    filename = "placeholder_filename2.png",
    content = function(filename) {
      png(filename)
      plot(iris$Sepal.Length, iris$Sepal.Width)
      dev.off()
    })
}

