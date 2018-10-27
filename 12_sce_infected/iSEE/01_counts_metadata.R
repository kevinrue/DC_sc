
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(iSEE)
}))

# Load the (barely) preprocessed object ----

sce <- readRDS("sce.rds")

# Show the count matrix
assayNames(sce)
assay(sce, "counts")

# Show the gene metadata
rowData(sce)

# Show the cell metadata
colData(sce)

# Launch iSEE! ----

app <- iSEE(sce)

shiny::runApp(app)
