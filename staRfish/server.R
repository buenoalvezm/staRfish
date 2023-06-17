library(shiny)
# library(data.table)
library(DT)

function(input, output, session) {
  rv <- reactiveValues()
  rv$s <- NULL

  observe({
    if(!is.null(input$study_table_rows_selected)){
      rv$s <- input$study_table_rows_selected
    }
  })
  data <- iris
  output$study_table <- DT::renderDataTable(DT::datatable(data))
  output$filtered_study_table <- DT::renderDataTable(
    datatable(data[rv$s,])
  )

}
