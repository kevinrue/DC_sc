options(iSEE.maxlevels=50)
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(iSEE)
}))

# Load the (barely) preprocessed object ----

useHDF5 <- TRUE

if (useHDF5) {
    sce <- readRDS("sce.h5.rds")
} else {
    sce <- readRDS("sce.rds")
}

# Preconfigure the initial state of the app ----

colDataArgs <- colDataPlotDefaults(sce, 1)

colDataArgs$YAxis <- "Infection"

colDataArgs$XAxis <- "Column data"
colDataArgs$XAxisColData <- "Status"

colDataArgs$ColorBy <- "Column data"
colDataArgs$ColorByColData <- "Time"

rowDataArgs <- rowDataPlotDefaults(sce, 1)

rowDataArgs$YAxis <- "source"

rowDataArgs$XAxis <- "Row data"
rowDataArgs$XAxisRowData <- "gene_biotype"

rowDataArgs$BrushData[[1]] <- list(
      xmin = 13.5, xmax = 14.5, ymin = 1.5, ymax = 2.5, mapping = list(x = "X", y = "Y"),
      log = list(x = NULL, y = NULL), direction = "xy",
      brushId = "rowDataPlot1_Brush", outputId = "rowDataPlot1")

rowStatArgs <- rowStatTableDefaults(sce, 1)
rowStatArgs$SelectByPlot <- "Sample assay plot 1"

sampAssayArgs <- sampAssayPlotDefaults(sce, 1)

sampAssayArgs$YAxisSampName <- 1

sampAssayArgs$XAxis <- "Sample name"
sampAssayArgs$XAxisSampName <- 2

sampAssayArgs$SelectByPlot <- "Row data plot 1"
sampAssayArgs$SelectBoxOpen <- TRUE

sampAssayArgs$BrushData[[1]] <- list(
    xmin = 15393, xmax = 60034, ymin = 8759, ymax = 122569,
    mapping = list(x = "X", y = "Y"),
    log = list(x = NULL, y = NULL), direction = "xy", brushId = "sampAssayPlot1_Brush",
    outputId = "sampAssayPlot1")

heatMapArgs <- heatMapPlotDefaults(sce, 1)

heatMapArgs$FeatNameBoxOpen <- TRUE
heatMapArgs$FeatName[[1]] <- which(seqnames(rowRanges(sce)) == "MT")

heatMapArgs$ColDataBoxOpen <- TRUE
heatMapArgs$ColData[[1]] <- c("Time", "Infection", "Status")

initialPanels <- DataFrame(
    Name=c(
        "Column data plot 1",
        "Feature assay plot 1",
        "Row statistics table 1",
        "Row data plot 1",
        "Sample assay plot 1",
        "Column statistics table 1",
        "Heat map 1"),
    Width=c(5, 5, 12, 8, 4, 12, 12)
)

# gene_biotype has 46 unique values (including NA)
# however, it should be treated as a categorical covariate
options(iSEE.maxlevels=50)

app <- iSEE(
    se = sce,
    colDataArgs = colDataArgs, rowDataArgs = rowDataArgs, rowStatArgs = rowStatArgs, sampAssayArgs = sampAssayArgs, heatMapArgs = heatMapArgs,
    initialPanels = initialPanels, appTitle = "Aulicino & Rue-Albrecht et al., 2018, Nat. Comm.")

# Launch iSEE! ----

shiny::runApp(app, launch.browser = TRUE)
