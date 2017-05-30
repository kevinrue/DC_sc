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
  )

))
