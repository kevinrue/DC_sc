#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  rv <- reactiveValues(
    matchedGeneIds = character(),
    normExprs = data.frame(),
    exprPlot = ggplot()
  )

  # New gene symbol selected, update choice of gene IDs
  observeEvent(
    input$geneSymbol, {
      message("Update choice of gene IDs")
      validate(need(
        input$geneSymbol %in% annTable$GENENAME,
        "Gene symbol not found in data set."
      ))

      rv[["matchedGeneIds"]] <-
        subset(
          annTable, GENENAME == input$geneSymbol, select = GENEID,
          drop = TRUE
        )

      updateSelectInput(session, "geneId", choices = rv[["matchedGeneIds"]])
    }
  )

  # Number of gene IDs available for selection
  output$countGeneIds <- renderUI({
    message("Update count of gene IDs")
    tagList(
      tags$code(length(rv[["matchedGeneIds"]])),
      "gene identifiers detected for gene symbol",
      tags$code(input$geneSymbol)
    )
  })

  # New gene ID selected, update expression data
  observeEvent(
    input$geneId,{
      req(input$geneId)
      message("New gene identifier")
      rv[["normExprs"]] <-
        cbind(
          data.frame(exprs = norm_exprs(sceset)[input$geneId,]),
          pData(sceset)[,c("Infection","Status","Time")]
        )
    }
  )

  # New expression data, update plot
  observeEvent(
    rv[["normExprs"]], {
      message("New expression data")
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # New row facet, update plot
  observeEvent(
    input$facetRow,{
      message("New facet row: ", input$facetRow)
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # New column facet, update plot
  observeEvent(
    input$facetCol,{
      message("New facet column: ", input$facetCol)
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # New colour, update plot
  observeEvent(
    input$colour,{
      message("New colour: ", input$colour)
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # New shape, update plot
  observeEvent(
    input$shape,{
      message("New shape: ", input$shape)
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # Toggle full Y range, update plot
  observeEvent(
    input$yRangeFull,{
      message("Full Y range: ", input$yRangeFull)
      rv[["exprPlot"]] <- drawExpProfile(
        rv[["normExprs"]], input$facetRow, input$facetCol,
        input$colour, input$shape, input$yRangeFull
      )
    }
  )

  # New plot, update output
  output$exprPlot <- renderPlot({
    message("Draw plot")
    validate(need(
      length(rv[["exprPlot"]]$data) > 0,
      "No data to plot."
    ))
    rv[["exprPlot"]]
  })

  # output$normExprs <- renderTable({
  #   message("Update sample table of expression data")
  #   validate(need(
  #     length(rv[["exprPlot"]]$data) > 0,
  #     "No data to show."
  #   ))
  #
  #   head(rv[["normExprs"]])
  # }, rownames = TRUE)

})
