suppressPackageStartupMessages(library(cBioPortalData, quietly = T))
suppressPackageStartupMessages(library(shiny, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(cbioportalR, quietly = T))
suppressPackageStartupMessages(library(magrittr, quietly = T))
suppressPackageStartupMessages(library(org.Hs.eg.db))
library(reshape2)
library(viridis)
library(ggpubr)
library(ggrepel)
library(colorspace)
library(annotate)


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

  protein <- read_tsv("../data_test/data_protein_quantification.txt") %>%
    separate(Composite.Element.REF,into=c("gene","gene_2")) %>% select(!gene_2) %>% na.omit()
  rna <- read_tsv("../data_test/data_mrna_seq_fpkm.txt") %>% dplyr::rename("gene"=Hugo_Symbol)
  metadata <- read_tsv("../data_test/brca_cptac_2020_clinical_data.tsv") %>% select(-c(`Study ID`,`Patient ID`))

  observeEvent(input$start_analyses,{
    gathered_data <- gather_rna_prot_data(rna = rna, protein = protein)
    correlations_df <- create_rna_prot_correlation(rna = rna, protein = protein)
    combined_pw_data <- create_pw_df(corr_df = correlations_df)
    plot_corr_matrix <- cor_matrix_samples(gathered_data = gathered_data)
    corr_plot_levels_df <- correlation_plot_levels(correlations_df)[1]
    corr_plot_levels_plot <- correlation_plot_levels(correlations_df)[[2]]
    KEGG_corr_plot <- KEGG_correlation_plot(combined_pw_data = combined_pw_data)

    # plot_gene_res <-
    ## Correlation plot
    # data_plot <- reactive({
    #   if(input$gene_name %in% unique(iris$Species)){
    #     iris %>% filter(Species == input$gene_name)
    #   } else {
    #       iris
    #     }
    #   })

    ### Transcriptomics vs proteomics plot ###
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

    ### KEG plot ###
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
      print(head(gathered_data))
      print(input$gene_names)
      plot_gene(gathered_data = gathered_data, gene_name = input$gene_names)

    })
  }

  # output$plot_correlation <- renderPlot({plot(data_plot()$Sepal.Length, data_plot()$Sepal.Width)})



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
