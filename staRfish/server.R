library(cBioPortalData)
library(shiny)
library(data.table)
library(DT)
library(tidyverse)
library(cbioportalR)
library(magrittr)

source("../R/grab.R")
source("../R/digest.R")
source("../R/squeeze.R")
source("../R/release_metadata.R")

function(input, output, session) {

  all_studies <- reactive({input$chosen_fields %>%
    strsplit(split = ",") %$%
    squeeze(n = input$n_samples, unlist(.))
    })

  output$study_table <- DT::renderDT({all_studies()})
  study_ids <- reactive({all_studies() %>% pull(studyId)})

  observe({
    updateSelectInput(session, "select_study", choices =   study_ids())
  })

  metadata_table <- DT::renderDT({release_metadata(input$select_study)})
  output$metadata_table <- metadata_table

  output$download_metadata <- downloadHandler(
    filename = "metadata.tsv",
    content = function(fname){
      write_tsv(release_metadata(input$select_study), fname)
    })

}
