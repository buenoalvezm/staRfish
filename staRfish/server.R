library(cBioPortalData, quietly = T)
library(shiny, quietly = T)
library(cbioportalR, quietly = T)
library(magrittr, quietly = T)
library(org.Hs.eg.db)
library(reshape2)
library(viridis)
library(ggpubr)
library(ggrepel)
library(colorspace)
library(annotate)
library(clusterProfiler)
library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(ggbeeswarm)


source("../R/grab.R")
source("../R/digest.R")
source("../R/squeeze.R")
source("../R/release_data.R")
source("../R/release_metadata.R")

select <- dplyr::select

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





  metadata_in <- reactive({
    release_metadata(study_id = input$select_study)
    })

  observeEvent(input$start_analyses,{

    # Getting data
    active_study_id <- reactive({
      input$select_study
    })
    #
    molecular_types <- reactive({
      input$chosen_fields %>%
        strsplit(split = ",") %>%
        unlist()
    })
    # #
    rna <- reactive({
      types <- molecular_types()
      rna_type <- types[grepl("Rna", types)][[1]]
      rna_data <- release_data(study_id = active_study_id(),
                               molecular_type = rna_type)
    })

    protein_data <- reactive({
      release_data(study_id = active_study_id(), molecular_type = "massSpectrometrySampleCount")
    })

    protein <- reactive({
      protein_data() %>%
        na.omit()
    })

    gathered_data <- gather_rna_prot_data(rna = rna(), protein = protein())
    correlations_df <- create_rna_prot_correlation(rna = rna(), protein = protein())
    combined_pw_data <- create_pw_df(corr_df = correlations_df)
    plot_corr_matrix <- cor_matrix_samples(gathered_data = gathered_data)
    corr_plot_levels_df <- correlation_plot_levels(correlations_df)[1]
    corr_plot_levels_plot <- correlation_plot_levels(correlations_df)[[2]]
    KEGG_corr_plot <- KEGG_correlation_plot(combined_pw_data = combined_pw_data)


    ## Transcriptomics vs proteomics plot ###
    output$plot_correlation <- renderPlot({corr_plot_levels_plot})


    output$download_transcrp_prot <- downloadHandler(
      filename = "tran_prot_corr.png",
      content = function(filename) {
        png(filename)
        plot(corr_plot_levels_plot)
        dev.off()
      })
    ###


    ### Big correlation matrix ###
    output$plot_big_correlation <- renderPlot({cor_matrix_samples(gathered_data = gathered_data)})

    output$download_big_correlation <- downloadHandler(
      filename = "correlation_matrix.png",
      content = function(filename) {
        png(filename)
        plot(cor_matrix_samples(gathered_data = gathered_data))
        dev.off()
      })

    ###

    ## KEG plot ###
    output$plot_kegg <- renderPlot({KEGG_corr_plot})

    output$download_plot_kegg <- downloadHandler(
      filename = "kegg_corr.png",
      content = function(filename) {
        png(filename)
        plot(KEGG_corr_plot)
        dev.off()
      })

    ### Plot gene ###

    plot_gene_res <- reactive({
      plot_gene(gathered_data = gathered_data, gene_name = input$gene_names)
    })

    output$plot_gene <- renderPlot({plot_gene_res()})
    output$download_plot_gene <- downloadHandler(
      filename = "gene.png",
      content = function(filename) {
        png(filename)
        plot(plot_gene_res)
        dev.off()
      })
    all_genes <- reactive({gathered_data %>%
        dplyr::pull(gene) %>%
        unique()
    })
    #
    observe({
      updateSelectizeInput(session, "gene_names", choices =   all_genes())
    })
    #
    plot_gene_res <- reactive({
      if (length(renderText({input$gene_names})) >0){
        plot_gene(gathered_data = gathered_data, gene_name = {input$gene_names})
      } else {
        plot(0)
      }
    })

    all_metadata_cols <- reactive({
      colnames(metadata_in())
      })
    # Plot metadata
    output$plot_metadata <- renderPlot({metadata_correlations(metadata = metadata_in() ,
                                                              gathered_data = gathered_data,
                                                              clin_var = input$select_metadata,
                                                              gene_name = input$gene_names_meta,
                                                              regression = input$regression)})
    output$download_plot_metadata <- downloadHandler(
      filename = "plot_metadata.png",
      content = function(filename) {
        png(filename)
        plot(metadata_correlations(metadata = metadata_in() ,
                              gathered_data = gathered_data,
                              clin_var = input$select_metadata,
                              gene_name = input$gene_names_meta,
                              regression = T))
        dev.off()
      })

    observe({
      updateSelectizeInput(session, "gene_names_meta", choices =   all_genes())
    })
    observe({
      updateSelectizeInput(session, "select_metadata", choices =  all_metadata_cols()[2:length(all_metadata_cols())])
    })


  })



  ###
}


# plot_gene_res <-
## Correlation plot
# data_plot <- reactive({
#   if(input$gene_name %in% unique(iris$Species)){
#     iris %>% filter(Species == input$gene_name)
#   } else {
#       iris
#     }
#   })
