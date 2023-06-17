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
                      selected = c("mrnaRnaSeqSampleCount","mrnaRnaSeqV2SampleCount",
                                   "massSpectrometrySampleCount")),
          numericInput("n_samples", "Min number of samples in study", value = 30, min = 0),
          # fluidRow(
          #   column(8,
          #           numericInput("n_samples", "Min number of samples in study", value = 30, min = 0),
          #          ),
          #   column(4,
          #           checkboxInput("minimum_N_analysis", "Set per data type", value = F),
          #          ),
          # ),
          # conditionalPanel(condition = "input.minimum_N_analysis==1",
          #                   textInput("n_samples_diiff", "Min number of samples per data type (comma separated)",placeholder = "30,30,30", value = "30,30,30")
          #                  )
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Datasets",),
            tabPanel("Data Overview"),
            tabPanel("Analysis1")
          )

        )
    )
)
