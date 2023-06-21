library(shiny)

fluidPage(

    # Application title
    titlePanel("cBioPortal data selection"),


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
          # textOutput("study_ids")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Study browser", tableOutput("study_table")),
            tabPanel("Active study", tableOutput("metadata_table")),
            tabPanel("Correlation")

          )
          # , downloadButton('download_metadata',"Download the data")

        )
    )
)
