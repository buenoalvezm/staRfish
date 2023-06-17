#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
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
        ),

        # Show a plot of the generated distribution
        mainPanel(

        )
    )
)
