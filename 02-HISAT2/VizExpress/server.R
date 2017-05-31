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
library(RColorBrewer)

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
      selected = c("TNF","IL1B", "CTSL"),
      multiple = TRUE
    )
  })

  geneIdPairs <- reactive({
    geneNames <- input$geneNamePairs
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

    time <- input$pairsSubsetTime
    infection <- input$pairsSubsetInfection
    status <- input$pairsSubsetStatus

    colour_by <- input$pairsColour
    shape_by <- input$pairsShape
    paletteName <- input$pairsPalette

    validate(need(geneIds, "Loading..."))
    validate(need(is.logical(ignoreGeneNamePairs), "Loading..."))
    validate(need(is.character(geneNames), "Loading..."))
    validate(need(is.character(geneIds), "Loading..."))
    validate(need(geneIdPairs, "Loading..."))

    validate(need(colour_by, "Loading..."))
    validate(need(shape_by, "Loading..."))
    validate(need(paletteName, "Loading..."))

    validate(need(length(geneIds) > 1, "At least 2 genes required"))
    validate(need(time, "At least 1 time level required"))
    validate(need(infection, "At least 1 infection level required"))
    validate(need(status, "At least 1 status level required"))

    gg <- data.frame(
      exprs = t(norm_exprs(sce.norm)[geneIds,]),
      pData(sce.norm)[,c("Time","Infection","Status")]
    )

    if (length(time)){
      gg <- subset(gg, Time %in% time)
    }
    if (length(infection)){
      gg <- subset(gg, Infection %in% infection)
    }
    if (length(status)){
      gg <- subset(gg, Status %in% status)
    }

    geneIdPairs <- geneIdPairs[match(geneIds, geneIdPairs[,"gene_id"]),]

    colour_map <-
      brewer.pal(nlevels(pData(sce.norm)[,colour_by]), paletteName)
    names(colour_map) <- levels(pData(sce.norm)[,colour_by])

    ggpairs(
      gg,
      columns = 1:length(geneIds),
      lower = list(
        continuous = wrap(
          pairs_point, colour_map = colour_map,
          colour_by = colour_by, shape_by = shape_by, alpha = 0.5)
      ),
      diag = list(
        continuous = wrap(
          pairs_density, colour_map = colour_map,
          fill_by = colour_by, alpha = 0.5)
      ),
      columnLabels = with(geneIdPairs, paste(gene_name, gene_id, sep = "\n")),
      legend = c(2,1)
    ) +
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
