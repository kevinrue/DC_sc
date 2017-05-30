#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

sce.norm <- readRDS("data/sce.norm.tSNE.rds")
exprs_range <- range(norm_exprs(sce.norm))

library(shiny)
library(ggplot2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  geneNames <- reactive({

    geneNames <- sort(fData(sce.norm)[,"gene_name"])
    geneNames <- geneNames[!is.na(geneNames)]
    geneNames

  })

  output$geneName <- renderUI({

    geneNames <- geneNames()
    req(geneNames)

    selectInput(
      "geneName", "Gene name",
      choices = geneNames,
      selected = "MARCH1"
    )
  })

  output$countGeneNames <- renderText({
    geneNames <- geneNames()
    req(geneNames)

    sprintf("Choices: %i", length(geneNames))
  })

  output$ignoreGeneName <- renderUI({

    checkboxInput(
      "ignoreGeneName", "Ignore",
      value = FALSE
    )

  })

  geneIds <- reactive({
    geneName <- input$geneName
    ignoreGeneName <- input$ignoreGeneName

    req(geneName, is.logical(ignoreGeneName))

    fd <- fData(sce.norm)[,c("gene_id","gene_name")]

    if (!ignoreGeneName) {
      fd <- subset(fd, gene_name == geneName)
    }

    return(fd[,"gene_id"])
  })

  output$countGeneIds <- renderText({
    geneIds <- geneIds()
    req(geneIds)

    sprintf("Choices: %i", length(geneIds))
  })

  output$geneId <- renderUI({

    geneIds <- geneIds()
    req(geneIds)

    selectInput(
      "geneId", "Gene ID",
      choices = geneIds
    )
  })

  output$violin <- renderPlot({
    geneId <- input$geneId
    req(geneId)
    ignoreGeneName <- input$ignoreGeneName
    geneName <- input$geneName

    gg <- data.frame(
      exprs = norm_exprs(sce.norm)[geneId,],
      pData(sce.norm)[,c("Time","Infection","Status")]
    )

    ggTitle <- geneId
    if (!ignoreGeneName){
      ggTitle <- sprintf("%s: %s", geneId, geneName)
    }

    ggplot(gg, aes(Time, exprs)) +
      geom_violin(
        aes(fill = Time), alpha = 0.5,
        draw_quantiles = seq(0.25, 0.75, 0.25)) +
      geom_jitter(width = 0.15, alpha = 0.5) +
      facet_grid(~ Infection + Status) +
      labs(title = ggTitle, y = "norm_exprs") +
      theme_minimal()

  })

})
