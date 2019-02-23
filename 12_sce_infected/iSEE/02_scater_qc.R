options(iSEE.maxlevels=50)
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(scater)
    require(HDF5Array)
    require(org.Hs.eg.db)
    require(iSEE)
}))

doPreprocess <- FALSE
useHDF5 <- TRUE

if (doPreprocess) {
    # Load the (barely) preprocessed object ----

    if (useHDF5) {
        sce <- readRDS("sce.h5.rds")
    } else {
        sce <- readRDS("sce.rds")
    }

    # Add log-transformed counts ----

    sce <- normalize(sce)

    if (useHDF5) {
        h5file <- "sce.h5"
        assay(sce, "logcounts") <- writeHDF5Array(assay(sce, "logcounts"), h5file, "logcounts", chunkdim = c(100, 100), verbose=TRUE)
    }

    # Compute scater QC metrics ----

    sce <- calculateQCMetrics(object = sce,
        feature_controls = list(
            ERCC=which(rowData(sce)$source == "ERCC"),
            MT=which(seqnames(rowRanges(sce)) == "MT")),
        cell_controls = list(
            Blank=which(sce$Status == "Blank"),
            Bulk=which(sce$Status == "Bulk")
        ))

    if (useHDF5) {
        saveRDS(sce, "sce.02.h5.rds")
    } else {
        saveRDS(sce, "sce.02.rds")
    }

}

if (useHDF5) {
    sce <- readRDS("sce.02.h5.rds")
} else {
    sce <- readRDS("sce.02.rds")
}

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
        "Column data plot 1",
        "Column statistics table 1",
        "Row data plot 1",
        "Row statistics table 1",
        "Sample assay plot 1"),
    Width=c(4, 8, 4, 8, 8)
)

# gene_biotype has 46 unique values (including NA)
# however, it should be treated as a categorical covariate
options(iSEE.maxlevels=50)

# Setting up the annotation function.
annot.fun <- annotateEnsembl(sce, orgdb=org.Hs.eg.db, keytype="ENSEMBL", rowdata_col="gene_id")

app <- iSEE(
    se = sce,
    colDataArgs = colDataArgs, colStatArgs = colStatArgs, rowDataArgs = rowDataArgs, rowStatArgs = rowStatArgs, sampAssayArgs = sampAssayArgs,
    annotFun = annot.fun,
    initialPanels = initialPanels, appTitle = "Aulicino & Rue-Albrecht et al., 2018, Nat. Comm.")

# Launch iSEE! ----

shiny::runApp(app, launch.browser = TRUE)
