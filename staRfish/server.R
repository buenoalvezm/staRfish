library(shiny)
library(data.table)
library(DT)
library(tidyverse)
library(cbioportalR)
library(magrittr)

source("../R/grab.R")
source("../R/process.R")
source("../R/squeeze.R")

function(input, output, session) {

  queried_table <- reactive({input$chosen_fields %>%
    strsplit(split = ",") %$%
    squeeze(n = input$n_samples, unlist(.))
    })

  # output$study_ids <- NULL
  output$study_table <- renderTable({queried_table()})
  output$study_ids <- renderText(c("wait for results", "test"))#, renderTable(queried_table()) %>% pull(studyId))








  }

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
