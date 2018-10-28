
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(scater)
    require(org.Hs.eg.db)
    require(RColorBrewer)
    require(iSEE)
}))

# Load the (barely) preprocessed object ----

sce <- readRDS("sce.rds")

# Add log-transformed counts ----

sce <- normalize(sce)

# Compute scater QC metrics ----

sce <- calculateQCMetrics(object = sce,
    feature_controls = list(
        ERCC=which(rowData(sce)$source == "ERCC"),
        MT=which(seqnames(rowRanges(sce)) == "MT")),
    cell_controls = list(
        Blank=which(sce$Status == "Blank"),
        Bulk=which(sce$Status == "Bulk")
    ))

# Preconfigure the initial state of the app ----

colDataArgs <- colDataPlotDefaults(sce, 1)

colDataArgs$YAxis <- "pct_counts_MT"

colDataArgs$XAxis <- "Column data"
colDataArgs$XAxisColData <- "pct_counts_ERCC"

colDataArgs$ColorBy <- "Column data"
colDataArgs$ColorByColData <- "Status"

colDataArgs$PointSize <- 5
colDataArgs$PointAlpha <- 0.6

colDataArgs$BrushData[[1]] <- list(
      xmin = -1, xmax = 101, ymin = 10, ymax = 100, mapping = list(x = "X", y = "Y"),
      log = list(x = NULL, y = NULL), direction = "xy",
      brushId = "colDataPlot1_Brush", outputId = "colDataPlot1")

colStatArgs <- colStatTableDefaults(sce, 1)

colStatArgs$SelectByPlot <- "Column data plot 1"
colStatArgs$SelectBoxOpen <- TRUE

rowDataArgs <- rowDataPlotDefaults(sce, 1)

rowDataArgs$YAxis <- "log10_total_counts_Bulk"

rowDataArgs$XAxis <- "Row data"
rowDataArgs$XAxisRowData <- "log10_total_counts_Blank"

rowDataArgs$ColorBy <- "Row data"
rowDataArgs$ColorByRowData <- "source"

rowDataArgs$PointSize <- 3
rowDataArgs$PointAlpha <- 0.6

rowDataArgs$BrushData[[1]] <- list(
      xmin = 1, xmax = 5.2, ymin = -0.25, ymax = 6, mapping = list(x = "X", y = "Y"),
      log = list(x = NULL, y = NULL), direction = "xy",
      brushId = "rowDataPlot1_Brush", outputId = "rowDataPlot1")

rowStatArgs <- rowStatTableDefaults(sce, 1)
rowStatArgs$SelectByPlot <- "Row data plot 1"
# rowStatArgs$SelectByPlot <- "Sample assay plot 1"

sampAssayArgs <- sampAssayPlotDefaults(sce, 1)

sampAssayArgs$Assay <- "logcounts"

sampAssayArgs$YAxisSampName <- 1

sampAssayArgs$XAxis <- "Sample name"
sampAssayArgs$XAxisSampName <- 2

sampAssayArgs$BrushData[[1]] <- list(
      xmin = 14, xmax = 17, ymin = 14, ymax = 17, mapping = list(x = "X", y = "Y"),
      log = list(x = NULL, y = NULL), direction = "xy",
      brushId = "sampAssayPlot1_Brush", outputId = "sampAssayPlot1")

initialPanels <- DataFrame(
    Name=c(
        "Reduced dimension plot 1"
    ),
    Width=c()
)

# gene_biotype has 46 unique values (including NA)
# however, it should be treated as a categorical covariate
options(iSEE.maxlevels=50)

# Setting up the annotation function ----

annot.fun <- annotateEnsembl(sce, orgdb=org.Hs.eg.db, keytype="ENSEMBL", rowdata_col="gene_id")

# Setting up colormaps ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
# col.infection <- brewer.pal(12, "Paired")[c(9,4,2,5)]
# col.status <- brewer.pal(12, "Paired")[c(9,1,7,3,5)]
names(col.time) <- levels(sce$Time)[1:3]
# names(col.infection) <- levels(sce$Infection)
# names(col.status) <- levels(sce$Status)

colormapInfection <- function(n) {
    x <- brewer.pal(12, "Paired")[c(9, 4, 2, 5)]
    names(x) <- c("Mock", "STM-LT2", "STM-D23580", "Blank")
    x
}

colormapStatus <- function(n) {
    x <- brewer.pal(12, "Paired")[c(9, 1, 7, 3, 5)]
    names(x) <- c("Uninfected", "Exposed", "Infected", "Bulk", "Blank")
    x
}

colormapTime <- function(n) {
    x <- brewer.pal(9, "Set3")[c(2, 7, 8, 9)]
    names(x) <- c("2h", "4h", "6h", "Blank")
    x
}

ecm <- ExperimentColorMap(
    colData=list(
        Infection=colormapInfection,
        Status=colormapStatus,
        Time=colormapTime
    )
)

app <- iSEE(
    se = sce,
    colDataArgs = colDataArgs, colStatArgs = colStatArgs,
    rowDataArgs = rowDataArgs, rowStatArgs = rowStatArgs,
    sampAssayArgs = sampAssayArgs, annotFun = annot.fun,
    initialPanels = initialPanels, colormap = ecm,
    appTitle = "Aulicino & Rue-Albrecht et al., 2018, Nat. Comm.")

# Launch iSEE! ----

shiny::runApp(app, launch.browser = TRUE)
