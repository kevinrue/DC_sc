#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(scater)
library(GGally)

sce.norm <- readRDS("data/sce.norm.tSNE.rds")

exprs_range <- range(norm_exprs(sce.norm))

geneNames <- sort(fData(sce.norm)[,"gene_name"])
geneNames <- geneNames[!is.na(geneNames)]

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # Gene ----

  output$ignoreGeneName <- renderUI({

    checkboxInput(
      "ignoreGeneName", "Ignore",
      value = FALSE
    )

  })

  output$geneName <- renderUI({

    selectInput(
      "geneName", "Gene name",
      choices = geneNames,
      selected = "MARCH1"
    )
  })

  output$countGeneNames <- renderText({
    sprintf("Gene names: %i", length(geneNames))
  })

  geneIds <- reactive({
    geneName <- input$geneName
    ignoreGeneName <- input$ignoreGeneName

    validate(need(geneName, "Loading..."))
    validate(need(is.logical(ignoreGeneName), "Loading..."))

    fd <- fData(sce.norm)[,c("gene_id","gene_name")]

    if (!ignoreGeneName) {
      fd <- subset(fd, gene_name == geneName)
    }

    return(fd[,"gene_id"])
  })

  output$countGeneIds <- renderText({
    geneIds <- geneIds()
    validate(need(geneIds, "Loading..."))

    sprintf("Gene IDs: %i", length(geneIds))
  })

  output$geneId <- renderUI({

    geneIds <- geneIds()
    validate(need(geneIds, "Loading..."))

    selectInput(
      "geneId", "Gene ID",
      choices = geneIds
    )
  })

  output$violin <- renderPlot({
    geneId <- input$geneId
    ignoreGeneName <- input$ignoreGeneName
    geneName <- input$geneName

    validate(need(geneId, "Loading..."))
    validate(need(is.logical(ignoreGeneName), "Loading..."))
    validate(need(geneName, "Loading..."))

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
      scale_y_continuous(limits = exprs_range) +
      theme_minimal()

  })

  # Gene pairs ----

  output$ignoreGeneNamePairs <- renderUI({

    checkboxInput(
      "ignoreGeneNamePairs", "Ignore gene names",
      value = FALSE
    )

  })

  output$geneNamePairs <- renderUI({

    selectInput(
      "geneNamePairs", "Gene name for pairs",
      choices = geneNames,
      selected = c("TNF", "IL1A", "IL1B", "IFNB1", "CTSL"),
      multiple = TRUE
    )
  })

  geneIdPairs <- reactive({
    geneNames <- input$geneNamePairs
    print(geneNames)
    ignoreGeneNames <- input$ignoreGeneNamePairs

    validate(need(is.character(geneNames), "Loading..."))
    validate(need(is.logical(ignoreGeneNames), "Loading..."))

    fd <- fData(sce.norm)[,c("gene_id","gene_name")]

    if (!ignoreGeneNames) {
      # Subset to gene names
      fd <- subset(fd, gene_name %in% geneNames)
      # Order by selected gene names
      fd[,"gene_name"] <- factor(fd[,"gene_name"], geneNames)
      fd <- fd[order(fd[,"gene_name"]),]
      fd[,"gene_id"] <- factor(fd[,"gene_id"], unique(fd[,"gene_id"]))
    }

    return(fd)
  })

  output$geneIdPairs <- renderUI({

    geneIdPairs <- geneIdPairs()

    validate(need(is.data.frame(geneIdPairs), "Loading..."))

    selectedGeneIds <- geneIdPairs[,"gene_id"]

    if (length(selectedGeneIds) == nrow(fData(sce.norm))){
      selectedGeneIds <- c()
    }

    selectInput(
      "geneIdPairs", "Gene ID for pairs",
      choices = geneIdPairs[,"gene_id"],  selected = selectedGeneIds,
      multiple = TRUE
    )
  })

  output$pairs <- renderPlot({
    geneIds <- input$geneIdPairs
    ignoreGeneNamePairs <- input$ignoreGeneNamePairs
    geneNames <- input$geneNamePairs
    geneIdPairs <- geneIdPairs()

    validate(need(length(geneIds) > 1, "At least 2 genes required"))

    validate(need(geneIds, "Loading..."))
    validate(need(is.logical(ignoreGeneNamePairs), "Loading..."))
    validate(need(is.character(geneNames), "Loading..."))
    validate(need(geneIdPairs, "Loading..."))

    gg <- data.frame(
      exprs = t(norm_exprs(sce.norm)[geneIds,]),
      pData(sce.norm)[,c("Time","Infection","Status")]
    )

    colnames(gg)[1:length(geneIds)] <- as.character(
      geneIdPairs[match(geneIds, geneIdPairs[,"gene_id"]), "gene_name"]
    )

    ggpairs(
      gg,
      columns = 1:length(geneIds),
      lower = list(continuous = my_ggpairs)) +
      theme_minimal()

  })

  output$geneIdPairsTable <- renderTable({
    geneIdPairs <- geneIdPairs()
    geneIds <- input$geneIdPairs
    ignoreGeneNamePairs <- input$ignoreGeneNamePairs

    validate(need(length(geneIds) > 1, "At least 2 genes required"))

    validate(need(is.data.frame(geneIdPairs), "Loading..."))
    validate(need(is.character(geneIds), "Loading..."))
    validate(need(is.logical(ignoreGeneNamePairs), "Loading..."))

    if (ignoreGeneNamePairs){
      geneIdPairs <- geneIdPairs[,"gene_id", drop = FALSE]
    }

    geneIdPairs <- subset(geneIdPairs, gene_id %in% geneIds)

    geneIdPairs

  }, striped = TRUE)

})
