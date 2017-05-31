#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)

paletteNames <- rownames(brewer.pal.info)

# Define UI for application that draws a histogram
shinyUI(navbarPage(

  title = "VizExpress",
  id = "navbar",
  tabPanel(
    title = "Gene",

    wellPanel(
      fluidRow(
        shiny::column(1, uiOutput("ignoreGeneName")),
        shiny::column(6, uiOutput("geneName")),
        shiny::column(5, textOutput("countGeneNames"))
      )
    ),
    wellPanel(
      fluidRow(
        shiny::column(6, uiOutput("geneId"), offset = 1),
        shiny::column(5, textOutput("countGeneIds"))
      )
    ),
    tabsetPanel(
      tabPanel(
        title = "Violin",
        wellPanel(
          plotOutput("violin")
        )
      ),
      id = "genePlots"
    )
  ),

  tabPanel(
    title = "Pairs",

    wellPanel(
      fluidRow(
        shiny::column(1, uiOutput("ignoreGeneNamePairs")),
        shiny::column(6, uiOutput("geneNamePairs"))
      )
    ),
    wellPanel(
      fluidRow(
        shiny::column(6, uiOutput("geneIdPairs"), offset = 1)
      )
    ),
    tabsetPanel(
      tabPanel(
        title = "Pairs",
        wellPanel(
          plotOutput("pairs")
        )
      ),
      tabPanel(
        title = "Table",
        wellPanel(
          tableOutput("geneIdPairsTable")
        )
      ),
      id = "genePairs"
    ),
    wellPanel(
      fluidRow(
        shiny::column(4, selectInput(
          "pairsSubsetTime", "Time",
          levels(sce.norm$Time), levels(sce.norm$Time),
          multiple = TRUE
        )),
        shiny::column(4, selectInput(
          "pairsSubsetInfection", "Infection",
          levels(sce.norm$Infection), levels(sce.norm$Infection),
          multiple = TRUE
        )),
        shiny::column(4, selectInput(
          "pairsSubsetStatus", "Status",
          levels(sce.norm$Status), levels(sce.norm$Status),
          multiple = TRUE
        ))
      ),
      fluidRow(
        shiny::column(4, selectInput(
          "pairsColour", "Colour",
          choices = c("Time", "Infection", "Status"), selected = "Time"
        )),
        shiny::column(4, selectInput(
          "pairsShape", "Shape",
          choices = c("Time", "Infection", "Status"), "Infection"
        )),
        shiny::column(4, selectInput(
          "pairsPalette", "Colour palette",
          choices = paletteNames, selected = "Set1"
        ))
      )
    )
  )

))
