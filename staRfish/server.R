library(shiny)
# library(data.table)
library(DT)
library(tidyverse)
library(cbioportalR)

source("../R/grab.R")
source("../R/process.R")

function(input, output, session) {
  rv <- reactiveValues()
  rv$s <- NULL

  observe({
    if(!is.null(input$study_table_rows_selected)){
      rv$s <- input$study_table_rows_selected
    } else {
      rv$s <- NULL
    }
  })
  output$n_samples <- reactive({input$n_samples})
  output$chosen_fields <- reactive({input$chosen_fields})

  data <- reactive({
    grab(output$n_samples, output$chosen_fields)
  })
  # data <- iris

  output$study_table <- renderDT(data())
  output$filtered_study_table <- renderDT(
    data()[rv$s,]
  )

}
