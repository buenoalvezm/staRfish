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

  # active_study_id <- reactive({
  #   input$select_study
  # })
  #
  molecular_types <- reactive({
    input$chosen_fields %>%
      strsplit(split = ",") %>%
      unlist()
  })
  #
  rna_data <- reactive({
    types <- unlist(molecular_types())
    rna_type <- types[grepl("Rna", types)][1]
    rna_data <- release_data(study_id = study_ids(),
                             molecular_type = rna_type)

    rna_data_filt <-
      rna_data() |>
      select(patientId, hugoGeneSymbol, value)

    return(rna_data_filt)
  })
  #
  protein_data <- reactive({
    protein_data <-  release_data(study_id = study_ids(), molecular_type = "massSpectrometrySampleCount")
    return(protein_data)
  })
  #
  output$rna_table <- renderTable({
    rna_data()
  })
  #
  #
  output$protein_table <- renderTable({
    protein_data()
  })

  protein <- read_tsv("../data_test/data_protein_quantification.txt") %>%
    separate(Composite.Element.REF,into=c("gene","gene_2")) %>% dplyr::select(!gene_2) %>% na.omit()
  rna <- read_tsv("../data_test/data_mrna_seq_fpkm.txt") %>% dplyr::rename("gene"=Hugo_Symbol)
  metadata <- read_tsv("../data_test/brca_cptac_2020_clinical_data.tsv") %>% dplyr::select(-c(`Study ID`,`Patient ID`))

  observeEvent(input$start_analyses,{
    gathered_data <- gather_rna_prot_data(rna = rna, protein = protein)
    print(gathered_data)
    correlations_df <- create_rna_prot_correlation(rna = rna, protein = protein)
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
      print(head(gathered_data))
      print(input$gene_names)
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
      print(input$gene_names)
      if (length(renderText({input$gene_names})) >0){
        plot_gene(gathered_data = gathered_data, gene_name = {input$gene_names})
      } else {
        plot(0)
      }
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
