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
shinyUI(fluidPage(

  # Application title
  titlePanel("Time profile in validated single cells"),

  wellPanel(
    fluidRow(
      shiny::column(
        width = 3,
        selectInput(
          "geneSymbol", label = "Gene symbol:",
          choices = annTable$GENENAME, "CD209")
      ),
      shiny::column(
        width = 3,
        selectInput(
          "geneId", label = "Gene identifier:",
          choices = character())
      ),
      shiny::column(
        width = 6,
        tags$h5("Summary"),
        uiOutput("countGeneIds")
      )
    )
  ),

  wellPanel(
    fluidRow(
      column(
        width = 3,
        selectInput(
          "facetRow", "Facet row",
          c("Infection", "Status", ".")
        )
      ),
      column(
        width = 3,
        selectInput(
          "facetCol", "Facet column",
          c("Status", "Infection", ".")
        )
      ),
      column(
        width = 3,
        selectInput(
          "colour", "Colour by",
          c("None", "Status", "Infection", "Time")
        )
      ),
      column(
        width = 3,
        selectInput(
          "shape", "Shape by",
          c("None", "Infection", "Status", "Time")
        )
      )
    ),
    fluidRow(
      column(
        width = 3,
        checkboxInput("yRangeFull", "Full Y range", TRUE)
      )
    ),
    fluidRow(
      column(
        width = 12,
        plotOutput("exprPlot")
      )
    )
  )
  # ,wellPanel(
  #   fluidRow(
  #     column(
  #       width = 12,
  #       tableOutput("normExprs")
  #     )
  #   )
  # )

))
