library(shiny)
library(data.table)
library(DT)
library(tidyverse)
library(cbioportalR)
library(cBioPortalData)
library(magrittr)

source("../R/grab.R")
source("../R/process.R")
source("../R/squeeze.R")
source("../R/release_data.R")
source("../R/release_metadata.R")


function(input, output, session) {

  all_studies <- reactive({input$chosen_fields %>%
      strsplit(split = ",") %$%
      squeeze(n = input$n_samples, unlist(.))
  })

  # output$study_ids <- NULL
  output$study_table <- renderTable({all_studies()})
  study_ids <- reactive({
    all_studies() %>% pull(studyId)
  })
  # print(study_ids)
  observe({
    updateSelectInput(session, "select_study", choices =   study_ids())
  })
  # %>% pull(studyId)
  # print(test %>% pull(studyId))
  # output$study_ids <- renderText(c("wait for results", "test"))#, renderTable(queried_table()) %>% pull(studyId))

 # active_study_id <- reactive({input$select_study})
  active_study_id <- reactive({
    input$select_study
  })

  molecular_types <- reactive({
    input$chosen_fields %>%
      strsplit(split = ",") %>%
      unlist()
  })

  rna_data <- reactive({
    types <- unlist(molecular_types())
    rna_type <- types[grepl("Rna", types)][1]
    rna_data <- release_data(study_id = active_study_id(),
                             molecular_type = rna_type)

    rna_data_filt <-
      rna_data |>
      select(patientId, hugoGeneSymbol, value)

    return(rna_data_filt)
  })

  protein_data <- reactive({
    protein_data <-  release_data(study_id = active_study_id(), molecular_type = "massSpectrometrySampleCount")
    return(protein_data)
  })

  output$rna_table <- renderTable({
    rna_data()
  })


  output$protein_table <- renderTable({
    protein_data()
  })




  # reactive({print(queried_table())})
  # print(queried_table)
  # isolate(queried_table)
  # print(input$chosen_fields)
  # print(observe(input$n_samples))
  # output$study_table <- DT::renderDataTable({DT::datatable((queried_data))})
  # output$study_table <- DT::renderDataTable({DT::datatable(renderDT(iris))})
  # output$study_table <- renderTable({iris})


  # rv <- reactiveValues()
  # rv$s <- NULL
  # rv$n_samples <- NULL
  # rv$chosen_fields <- NULL
  # rv$queried_data <- NULL
  #
  # observe({
  #   if(!is.null(input$study_table_rows_selected)){
  #     rv$s <- input$study_table_rows_selected
  #   } else {
  #     rv$s <- NULL
  #   }
  # })

  # observeEvent(input$n_samples, {
  #   rv$n_samples <- input$n_samples
  # })
  #
  # observeEvent(input$chosen_fields, {
  #   rv$chosen_fields <- input$chosen_fields
  # })


  # n_samples <- reactive({input$n_samples})
  # chosen_fields <- reactive({input$chosen_fields})

  # rv$queried_data <- reactive({
  #   input$chosen_fields %>%
  #     strsplit(split = ",") %$%
  #     grab(n = input$n_samples, unlist(.))
  # })


  # write_rds(queried_data(), "../data_test/test_cbioportal.RDS")
  # write.table(rv$big_table, "../data_test/cbioportal_export.csv")
  # readRDS("../data_test/test_cbioportal.RDS")

  # queried_dataqueried_data <- iris
  # queried_data <- fread(file = "../data_test/cbioportal_export.csv")
  # output$study_table <-renderDT(queried_data())
  # output$study_table <- DT::renderDataTable({
  #   DT::datatable(rv$queried_data())})
  # output$study_table <-renderDT(queried_data() %>% select(c(studyId, name, importDate,
                                                          # allSampleCount, cancerTypeId, referenceGenome, citation)))

  # output$filtered_study_table <- renderDT(
  #   # queried_data()[rv$s,]
  #   queried_data()[rv$s,] %>% select(c(studyId, name, importDate,
  #                       allSampleCount, cancerTypeId, referenceGenome, citation))
  # )
