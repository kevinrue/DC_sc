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
    normExprs = data.frame()
  )

  observeEvent(
    input$geneSymbol, {

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

  output$countGeneIds <- renderUI({
    tagList(
      tags$code(length(rv[["matchedGeneIds"]])),
      "gene identifiers detected for gene symbol",
      tags$code(input$geneSymbol)
    )
  })

  observeEvent(
    input$geneId,{
      req(input$geneId)

      rv[["normExprs"]] <-
        cbind(
          data.frame(exprs = norm_exprs(sceset)[input$geneId,]),
          pData(sceset)[,c("Infection","Status","Time")]
        )
    }
  )

  output$exprPlot <- renderPlot({
    validate(need(
      nrow(rv[["normExprs"]]) > 0,
      "No data to plot."
    ))
    message(input$geneId)
    ggplot(rv[["normExprs"]], aes(x = Time, y = exprs)) +
      geom_violin(draw_quantiles = seq(0.25, 0.75, 0.25)) + geom_jitter(width = 0.25) +
      facet_grid(Infection ~ Status)
  })

  output$normExprs <- renderTable({
    validate(need(
      nrow(rv[["normExprs"]]) > 0,
      "No data to show."
    ))

    head(rv[["normExprs"]])
  }, rownames = TRUE)

})
