library(shiny)

fluidPage(

  # Application title
  titlePanel("StaRfish Data Explorer"),


  sidebarLayout(
    sidebarPanel(
      selectizeInput("chosen_fields", "Select desired data",
                     choices = c(CNA = "cnaSampleCount",
                                 mRNAseq = "mrnaRnaSeqSampleCount",
                                 mRNAseq2 = "mrnaRnaSeqV2SampleCount",
                                 mRNA_microarray = "mrnaMicroarraySampleCount",
                                 miRNA = "miRnaSampleCount",
                                 metHm27 = "methylationHm27SampleCount",
                                 rppa = "rppaSampleCount",
                                 massSpec ="massSpectrometrySampleCount"),
                     multiple = T,
                     selected = c("mrnaRnaSeqSampleCount",
                                  "massSpectrometrySampleCount")),
      numericInput("n_samples", "Min number of samples in study", value = 30, min = 0),
      selectInput("select_study",label="Select study", choices= "Please choose one", multiple = F),
      downloadButton('download_metadata',"Download metada"),
      actionButton("start_analyses", "Start analyses!")
      # textOutput("study_ids")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Study browser",
                 DT::DTOutput("study_table")
        ),
        tabPanel("Active study",
                 DT::DTOutput("metadata_table")
        ),
        tabPanel("Correlation",
                 downloadButton('download_big_correlation',"Download plot"),
                 plotOutput("plot_big_correlation")
        ),
        tabPanel("Transcriptomics vs Proteomics",
                 downloadButton('download_transcrp_prot',"Download plot"),
                 plotOutput("plot_correlation"),
                 selectizeInput("gene_names", "Select Gene",
                                choices = "Start selecting",
                                multiple =F),
                 plotOutput("plot_gene"),
                 downloadButton('download_corr_plot',"Download plot")

        ),
        # tabPanel("Test",
        # ),
        tabPanel("KEGG Correlation",
                 downloadButton('download_plot_kegg',"Download plot"),
                 plotOutput("plot_kegg")
        ),
        tabPanel("test",
                 tableOutput("protein_table")
                 )


      )
    )
  )
)
